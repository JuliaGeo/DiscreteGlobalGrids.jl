# The per-task cache of derived node extents, shared by every tree cursor that
# computes an extent rather than storing one. One table per tree, a few tables
# per task, so a search that alternates between trees keeps every one of them.

"""
    EXTENT_MEMO_SLOTS

Slots in one [`ExtentTable`](@ref). A power of two, so a node's slot is a mask
rather than a modulo.
"""
const EXTENT_MEMO_SLOTS = 1024

"""
    EXTENT_MEMO_TABLES

Tables one task holds at once for a given node-key and tree-identity type. The
count bounds a task's memo memory; a handful covers a dual-tree join and a
per-chunk window tree read beside the holding it came from.
"""
const EXTENT_MEMO_TABLES = 4

# The key no node has, so a fresh or cleared slot never compares equal to one.
@inline no_extent_key(::Type{NTuple{N,Int}}) where {N} = ntuple(_ -> -1, Val(N))

# A node's slot: the whole key hashed down to the table width.
@inline extent_slot(key) = Int(hash(key) & UInt(EXTENT_MEMO_SLOTS - 1)) + 1

"""
    ExtentTable{K,I}(id)

One tree's [`EXTENT_MEMO_SLOTS`](@ref) direct-mapped extent slots: node key type
`K`, tree identity `id` of type `I`.

  - A node's key hashes to exactly one slot; a hit is a key compare and a load.
  - A miss derives the extent and overwrites whatever sat there. The slot holds
    the whole key and compares it, so a collision costs a re-derive, never a
    wrong extent.
  - `misses` counts the derivations since the table was keyed to `id`.
"""
mutable struct ExtentTable{K,I}
    id::I
    misses::Int
    stamp::Int              # last use, on the store's clock
    const keys::Vector{K}
    const vals::Vector{Cap}
end

ExtentTable{K,I}(id::I) where {K,I} = ExtentTable{K,I}(id, 0, 0,
    fill(no_extent_key(K), EXTENT_MEMO_SLOTS), Vector{Cap}(undef, EXTENT_MEMO_SLOTS))

Base.show(io::IO, t::ExtentTable{K}) where {K} =
    print(io, "ExtentTable{", K, "}(", t.misses, " derived)")

"""
    ExtentMemo{K,I}()

One task's [`ExtentTable`](@ref)s for node key type `K` and tree identity type
`I`, at most [`EXTENT_MEMO_TABLES`](@ref) of them.

  - A tree is found by an identity scan over the tables held, so turning to
    another tree costs a compare per table and leaves the tables already keyed
    intact: alternating between trees hits in every one of them.
  - A miss past the table count re-keys the least recently used table, so a tree
    read on every ask survives a stream of trees read once each.
  - Task-local and lock-free: tasks sharing a tree have separate tables.
"""
mutable struct ExtentMemo{K,I}
    const tables::Vector{ExtentTable{K,I}}
    clock::Int
end

ExtentMemo{K,I}() where {K,I} = ExtentMemo{K,I}(ExtentTable{K,I}[], 0)

"""
    extent_table(::Type{K}, id) -> ExtentTable{K,typeof(id)}

The calling task's table of node extents for the tree `id` stands for, keyed on
node keys of type `K`.

`id` fixes every extent the table holds: two nodes with equal keys under the
same `id` have the same extent. It is compared with `===` and nothing else, so
any value that pins the geometry will do — the tree object itself, or a tuple of
the grid, system and level a lazy cursor derives its boxes from.
"""
@inline function extent_table(::Type{K}, id::I) where {K,I}
    store = get!(ExtentMemo{K,I}, task_local_storage(),
        ExtentMemo{K,I})::ExtentMemo{K,I}
    tables = store.tables
    for n in eachindex(tables)
        table = @inbounds tables[n]
        if table.id === id
            table.stamp = (store.clock += 1)
            return table
        end
    end
    return _rekey_extent_table!(store, id)
end

# Take a table for `id`: a fresh one while the store is under its count, else
# the least recently used, cleared and re-keyed.
function _rekey_extent_table!(store::ExtentMemo{K,I}, id::I) where {K,I}
    tables = store.tables
    if length(tables) < EXTENT_MEMO_TABLES
        table = ExtentTable{K,I}(id)
        push!(tables, table)
    else
        table = @inbounds tables[1]
        for n in 2:length(tables)
            @inbounds tables[n].stamp < table.stamp && (table = @inbounds tables[n])
        end
        fill!(table.keys, no_extent_key(K))
        table.id = id
        table.misses = 0
    end
    table.stamp = (store.clock += 1)
    return table
end

"""
    memoized_extent(derive, table, key) -> Cap

`key`'s extent from `table`, deriving it with `derive()` and storing it on a
miss.

`derive` must answer for the tree the table is keyed to, and must answer the
same cap for the same key every time; the memo is then invisible, bit for bit,
whoever asks and in whatever order.
"""
@inline function memoized_extent(derive::F, table::ExtentTable{K}, key::K) where {F,K}
    s = extent_slot(key)
    @inbounds table.keys[s] == key && return @inbounds table.vals[s]
    extent = derive()::Cap
    table.misses += 1
    @inbounds table.keys[s] = key
    @inbounds table.vals[s] = extent
    return extent
end
