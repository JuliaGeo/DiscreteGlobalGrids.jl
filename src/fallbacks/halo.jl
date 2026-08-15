# The outside face of a subtree boundary: the level-`l` cells that are NOT
# descendants of `c` but have a neighbour that is.
#
# One public type, one private ENGINE per system, exactly as the rim and
# interior walks are built (`subtree_iterators.jl`). The engine is chosen once,
# at construction, by `halo_engine` — private multiple dispatch, so an engine's
# whole protocol is monomorphic and a system adds a fast path by adding one
# method, not by setting a flag anyone can read.
#
# Every engine here is EXACT. A specialization may enumerate a conservative
# candidate band, but it must run the adjacency test on every candidate before
# yielding it: approximation never reaches this file's public surface.

"""
    SubtreeHaloIterator(sys, c, l; connectivity = Vertex())

The halo of `c`'s subtree at level `l`, lazily: every level-`l` cell that is
**not** a descendant of `c` but has a neighbour that is, in ascending canonical
order, each cell exactly once.

`collect` of this is [`subtree_halo`](@ref), element for element. `l == level(c)`
is `c`'s own one-ring, sorted. `l < level(c)` and `l > max_level(sys)` throw an
`ArgumentError`.

THE ORDER IS A CONTRACT, and it is worth naming the order precisely, because
"canonical" is this package's word and a caller needs the operational one:
`cellposition(levelgrid(sys, l), x)` is STRICTLY INCREASING over the walk. So a
consumer indexing a position-ordered array — a chunk, a `Vector` laid out by
`cellposition`, a cube axis — needs no `sort` and no `sortperm`, on any system
and through any engine. [`halo_positions`](@ref) is that stream of positions
directly.

NOT TO BE CONFUSED WITH ONE-RING ORDER, which is a different question with a
different answer. `neighbors` is rotationally ordered on H3 and IGeo7, and
ascending position holds for none of it: 0 of the 41,162 H3 level-3 cells have a
one-ring that is already ascending. The guarantee here is about the HALO, which
the engines emit in position order by construction — `RingHaloEngine` by
selection emit at depth zero, everything deeper by walking a canonical
enumeration it never has to repair.

Construction never materialises the halo: the iterator holds `O(depth)` walk
state and bounded native neighbour containers, and nothing sized by the answer.
Taking the first few cells of a deep halo therefore costs what those cells cost,
not what the whole ring would.

[`Base.IteratorSize`](@ref) is `SizeUnknown()` wherever no count is proved —
face seams, poles and pentagons break perimeter formulas — and `HasLength()`
only where an engine derives the count in closed form. There is no `length` that
would silently cost a traversal. [`halo_sizehint`](@ref) is the separate,
explicitly approximate answer for a caller who only wants to `sizehint!`.

See also [`halo`](@ref) for the same question about a subset,
[`halo_positions`](@ref) for the same walk in position space, and
[`EdgeCellIterator`](@ref) for the inside face of the same boundary.
"""
struct SubtreeHaloIterator{S<:AbstractHierarchicalGridSystem,C<:AbstractCellIndex,
        K<:Connectivity,E}
    system::S
    cell::C
    level::Int
    connectivity::K
    engine::E
end

function SubtreeHaloIterator(sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex, l::Integer; connectivity::Connectivity = Vertex())
    target = Int(l)
    return SubtreeHaloIterator(sys, c, target, connectivity,
        halo_engine(sys, c, target, connectivity))
end

Base.iterate(it::SubtreeHaloIterator) = iterate(it.engine)
Base.iterate(it::SubtreeHaloIterator, state) = iterate(it.engine, state)

# `C`, NOT `eltype(E)`, and the difference is inference. The engine parameter is
# whatever `halo_engine` dispatched to, and that return type is a union over
# every engine the call might pick — wide enough, once a system specializes,
# that inference gives up rather than split it, leaving `eltype(E)` unresolved
# and `collect_subtree`'s `eltype(it)[]` widened to `Vector{Any}`. `C` is the
# cell that was asked about, which is the cell type the answer is in, so reading
# it costs inference nothing. The engines still declare their own `eltype`; the
# sweep "eltype is the system's cell index type, on every engine" pins the two
# readings together, and "subtree_halo's return type is inferred, not just
# correct" pins this one.
Base.eltype(::Type{<:SubtreeHaloIterator{S,C}}) where {S,C} = C
Base.IteratorSize(::Type{<:SubtreeHaloIterator{S,C,K,E}}) where {S,C,K,E} =
    Base.IteratorSize(E)

# Deliberately delegated rather than defined: an engine that cannot count
# without walking defines no `length`, and the `MethodError` from here is the
# contract being kept.
Base.length(it::SubtreeHaloIterator) = length(it.engine)

# `collect` must BE the guarded path, not merely parallel to it: `collect`'s own
# `HasLength` route would size its vector from a miscounting engine and hand
# back an `undef` tail as cell ids. See `collect_subtree`.
#
# The hint is the second argument and never the first: `collect_subtree` still
# reads `IteratorSize` to decide whether the count guard arms, so an approximate
# number reaches `sizehint!` and nothing else. See `halo_sizehint`.
Base.collect(it::SubtreeHaloIterator) = collect_subtree(it, halo_sizehint(it))

Base.show(io::IO, it::SubtreeHaloIterator) = print(io, "SubtreeHaloIterator(",
    it.system, ", ", it.cell, ", ", it.level, "; connectivity = ",
    it.connectivity, ")")

# ===========================================================================
# The one-ring engine: `l == level(c)`
# ===========================================================================

# The ring comes back in ROTATIONAL order and the halo owes ascending order, so
# the sort is a selection emit over a container of at most `max_neighbors`
# elements: the state is the count emitted plus the last value, both isbits, and
# nothing is copied or allocated. A `sort` here would have to leave the fixed
# capacity container and heap-allocate (see `systems/H3/neighbors.jl`).
"""
    RingHaloEngine(ring)

`c`'s own one-ring, ascending, by selection emit. `O(degree^2)` time with
`degree <= max_neighbors(sys, connectivity)`, no allocation, isbits state.

`length` is `length(ring)`, which is honest only because a native one-ring
never lists the same cell twice — the selection emit yields each DISTINCT value
once, so a repeated neighbour would make the engine yield fewer cells than it
claims. No bundled system repeats; if one ever did, `collect_subtree` is the
backstop that turns the miscount into an `error` rather than an `undef` tail.
"""
struct RingHaloEngine{V,C}
    ring::V
end

RingHaloEngine(ring) = RingHaloEngine{typeof(ring),eltype(ring)}(ring)

Base.eltype(::Type{<:RingHaloEngine{V,C}}) where {V,C} = C
Base.IteratorSize(::Type{<:RingHaloEngine}) = Base.HasLength()
Base.length(e::RingHaloEngine) = length(e.ring)

struct RingHaloState{C}
    emitted::Int
    last::C
end

function Base.iterate(e::RingHaloEngine{V,C}) where {V,C}
    isempty(e.ring) && return nothing
    best = first(e.ring)
    for x in e.ring
        x < best && (best = x)
    end
    return (best, RingHaloState(1, best))
end

function Base.iterate(e::RingHaloEngine{V,C}, s::RingHaloState{C}) where {V,C}
    s.emitted >= length(e.ring) && return nothing
    found = false
    best = s.last
    for x in e.ring
        x <= s.last && continue
        if !found || x < best
            best = x
            found = true
        end
    end
    found || return nothing
    return (best, RingHaloState(s.emitted + 1, best))
end

# ===========================================================================
# Generic engine construction
# ===========================================================================

# The halo's level guard, factored out because it has to be CALLABLE from a
# system's own engine. `interface/system.jl` promises that an engine's
# `ArgumentError` is the one the eager verb raises, and a specialization that
# validated through `descendant_range` instead answered the very same user error
# in the very same call with different words than the generic engine did — on
# systems the caller is meant to be able to swap. One guard, one wording.
#
# "past max_level" and not "below the deepest level": nothing is below the
# deepest, and a reader parses "below 19" as "< 19". The rest of the package
# already says "past max_level".
function check_halo_level(sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex, target::Int)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "subtree_halo: level $target is above the cell's own level $lc"))
    target <= max_level(sys) || throw(ArgumentError(
        "subtree_halo: level $target is past max_level $(max_level(sys))"))
    return nothing
end

"""
    generic_halo_engine(sys, c, target, connectivity)

The engine every specialization's guards return to, and the one a newly
registered system inherits: the level guard, then the one-ring at depth zero,
then [`OutsideWalkEngine`](@ref) — or [`ScanHaloEngine`](@ref) on a system with
no [`descendant_range`](@ref) to prune by.

Named rather than inlined into `halo_engine` because a system's own method calls
this when a guard fires, which is what makes a specialization a speed decision
and never a correctness one.
"""
function generic_halo_engine(sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex, target::Int, connectivity::Connectivity)
    check_halo_level(sys, c, target)
    lc = level(c)
    grid = levelgrid(sys, target)
    target == lc && return RingHaloEngine(neighbors(grid, c, 1; connectivity))
    return outside_walk_engine(sys, c, target, connectivity, IndexedNeighbors())
end

halo_engine(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
    target::Int, connectivity::Connectivity) =
    generic_halo_engine(sys, c, target, connectivity)

# The ellipsoid wrapper is transparent: the halo of a subtree is a question about
# the discrete hierarchy, which the authalic transform does not touch.
halo_engine(sys::AuthalicSystem, c::AbstractCellIndex, target::Int,
    connectivity::Connectivity) = halo_engine(sys.system, c, target, connectivity)

"""
    subtree_halo(sys, c, l; connectivity = Vertex()) -> Vector

The level-`l` cells outside `c`'s subtree that touch it, ascending. The
explicitly materialising form of [`SubtreeHaloIterator`](@ref) — the halo can be
far larger than the subtree's rim, so collecting it is the caller's decision,
never the constructor's. Ascending means ascending `cellposition` on
`levelgrid(sys, l)`; see the iterator's docstring for the contract.

THE MOTIVATING READ IS A CHUNK PLUS ITS STENCIL MARGIN, AND IT IS NOT THIS VERB.
On the five systems with [`has_sorted_subtrees`](@ref) — every one but A5, which
has no [`descendant_range`](@ref) method to call — `descendant_range(sys, c, l)`
is the chunk's contiguous position block, and the extra fetch list a one-ring
stencil needs is [`halo_positions`](@ref), which streams:

    for p in halo_positions(sys, c, l)
        margin[p] = source[p]
    end

No halo table, no id `Vector`, no position `Vector`. Writing that read as
`cellposition.(Ref(levelgrid(sys, l)), subtree_halo(sys, c, l))` would build both
— for a 30-million-cell halo, 480 MB of ids to make 240 MB of positions and then
drop the ids — which is spending the laziness this file is built out of in order
to arrive at the same integers.

Reach for `subtree_halo` when the IDS are the answer and are wanted all at once:
naming cells to a system that indexes by id, or handing a halo to something that
must index it more than once. `sizehint!` it with [`halo_sizehint`](@ref) if you
are building the `Vector` yourself.

The position stream is also the second argument of [`stencil_table`](@ref), which
is how the two halves get addressed once they are laid end to end in one buffer.
"""
subtree_halo(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        l::Integer; connectivity::Connectivity = Vertex()) =
    collect(SubtreeHaloIterator(sys, c, l; connectivity))

"""
    halo(subset; connectivity = Vertex())

The cells immediately outside a same-level subset — the subset face of
[`SubtreeHaloIterator`](@ref). Defined for [`PartialGrid`](@ref),
[`CellVector`](@ref) and
[`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup), whose methods
live beside the other subset verbs in `stencil.jl` and `dimensionaldata.jl`.

The same definition a subtree gets, read against membership instead of ancestry:
every cell of the subset's level that the subset does **not** hold but that has a
neighbour the subset does, ascending, each cell once. `Vertex()` counts vertex
contact, `Edge()` requires a shared edge.

Ascending is the same contract [`SubtreeHaloIterator`](@ref) states, read on the
COMPLETE level grid: `cellposition(subset, x)` is `nothing` for every cell of a
halo by definition, so the grid a position can mean here is the complete one the
subset was cut from. [`halo_positions`](@ref) of this iterator yields those.

**A hole is part of the halo.** A cell removed from the middle of a subset is
outside the subset and has in-subset neighbours, so it is a halo cell. That is
the definition read literally rather than a special case, and it is the reason
this verb answers something [`subtree_halo`](@ref) cannot.

ALWAYS AN ITERATOR, on every container and every path. A rooted `PartialGrid`
holding a complete subtree returns a [`SubtreeHaloIterator`](@ref) — the subtree
walk is the specialized one, and re-deriving it from membership would be slower
for no gain — and anything else returns a [`SubsetHaloIterator`](@ref). Both are
lazy, both hold `O(depth)` state beyond the subset's own storage, and the
concrete type is settled once at construction, exactly as [`halo_table`](@ref)'s
two branches are. So `for x in halo(sub)` gets the same laziness either way, and
`collect` is the explicit materializing form.

Not to be confused with [`halo_table`](@ref), which is the IN-SET positional
stencil: one row of in-set neighbour positions per cell of the subset. That verb
answers "which of my own cells does each of my cells touch"; this one answers
"which cells that I do not hold touch me". Neither replaces the other, and
[`stencil_table`](@ref) is the verb that needs both — it takes this list, once
materialised as ascending positions, and completes `halo_table`'s short rim rows
against it.

A [`MultiOrderCellSet`](@ref) has no method here and will not grow one: its
members sit at different levels, so there is no one level for a halo to answer
at. Mixed-level adjacency is [`member_neighbors`](@ref).
"""
function halo end

"""
    SubsetHaloIterator(subset, connectivity, engine)

The halo of an arbitrary same-level subset, lazily — what [`halo`](@ref) returns
whenever the subset is not a rooted complete subtree: a subset with a hole, a
subset whose root was forgotten, a [`CellVector`](@ref), a
[`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup).

Positional and built by [`halo`](@ref), which owns the engine choice. There is
nothing a caller could pass here that `halo` does not already decide from the
subset itself.

Construction costs `O(1)` and materialises nothing — not the halo, and not a
summary of the input either. The walk is [`OutsideWalkEngine`](@ref)'s, so the
state is one `O(depth)` frame stack, and it prunes by the subset's own
[`subset_span`](@ref)s. [`Base.IteratorSize`](@ref) is the engine's, which is
`SizeUnknown()` — a subset's halo has no perimeter formula at all.
"""
struct SubsetHaloIterator{S,K<:Connectivity,E}
    subset::S
    connectivity::K
    engine::E
end

Base.iterate(it::SubsetHaloIterator) = iterate(it.engine)
Base.iterate(it::SubsetHaloIterator, state) = iterate(it.engine, state)

# THE ONE ITERATOR THAT STILL DELEGATES, and it can: there is no cell here to
# key off — the subject is a membership predicate over an arbitrary set, not a
# root — and `subset_halo_engine` picks between exactly two engine types, so the
# union is narrow enough for inference to resolve `eltype(E)`. If a third subset
# engine ever appears, this line is what to re-measure: the comment on
# `SubtreeHaloIterator` above is what happens when the union outgrows it.
Base.eltype(::Type{<:SubsetHaloIterator{S,K,E}}) where {S,K,E} = eltype(E)
Base.IteratorSize(::Type{<:SubsetHaloIterator{S,K,E}}) where {S,K,E} =
    Base.IteratorSize(E)

# Delegated rather than defined, for `SubtreeHaloIterator`'s reason: no subset
# engine counts without walking, so the `MethodError` from here is the contract.
Base.length(it::SubsetHaloIterator) = length(it.engine)

Base.collect(it::SubsetHaloIterator) = collect_subtree(it, halo_sizehint(it))

Base.show(io::IO, it::SubsetHaloIterator) = print(io, "SubsetHaloIterator(",
    it.subset, "; connectivity = ", it.connectivity, ")")

# ===========================================================================
# The same walk, in position space
#
# WHY THIS IS A VERB AND NOT A `map`. A halo is a fetch list, and a fetch list
# is read against storage laid out by `cellposition` — a chunk, a cube axis, a
# `Vector` the caller allocated `ncells` long. So the caller's last step is
# almost always the position, and the id was a way station. Two of the engines
# never needed it in the first place: `ScanHaloEngine` walks POSITIONS and calls
# `cellindex` to make an id, and `SquareBandEngine` computes `face·n² + child`,
# which is the position minus one, before wrapping it in a `LevelIndex`.
#
# Yielding positions therefore costs nothing worth measuring on HEALPix, S2 and
# ISEA4R (a measured 0.64-0.71 ns per cell for `cellposition`, which on those
# three is `index + 1` behind one type check) and is strictly cheaper on A5,
# where the scan's own state IS the position and re-deriving it through
# `cellposition` costs a measured 38 ns of Hilbert decode per cell. It genuinely
# costs a conversion on the two aperture-7 systems — a measured 7.7 ns on IGeo7
# and 16.6 ns on H3, whose `cellposition` is an FFI validity check plus a child
# position — and that is the honest price of a system whose canonical id is not
# its position. `test/systems/crosssystem/subtree_halos.jl` is where those six
# numbers are reproduced.
#
# WHAT IT IS NOT. Not `cellposition.(collect(halo))`: that materialises the ids
# to throw them away, which is the whole cost this exists to remove. The
# iterator holds the halo iterator and the grid, threads the halo's own state,
# and adds nothing sized by the answer.
# ===========================================================================

"""
    HaloPositionIterator(halo, grid)

A halo walk read as `cellposition`s on `grid`, lazily — what
[`halo_positions`](@ref) returns.

Yields `Int`, strictly increasing, one per cell of the underlying walk and in
the same order. Everything else is the wrapped iterator's:
[`Base.IteratorSize`](@ref), the `length` that exists on exactly two engines and
on no others, resumability, and the `O(depth)` state.
"""
struct HaloPositionIterator{I,G}
    halo::I
    grid::G
end

# The grid a position is a position IN. For a subtree it is the target level;
# for a subset it is the COMPLETE grid the walk enumerates from, which both
# subset engines already carry as `grid` — the subset's own `cellposition`
# answers `nothing` for every halo cell, so it cannot be the one meant.
_halo_grid(it::SubtreeHaloIterator) = levelgrid(it.system, it.level)
_halo_grid(it::SubsetHaloIterator) = it.engine.grid

"""
    halo_positions(sys, c, l; connectivity = Vertex()) -> HaloPositionIterator
    halo_positions(it) -> HaloPositionIterator

The halo as POSITIONS rather than ids: `cellposition` on `levelgrid(sys, l)` for
every cell [`SubtreeHaloIterator`](@ref) would yield, strictly increasing, lazily.

The one-argument form takes a halo iterator, so a subset's halo composes —
`halo_positions(halo(pg))` — and the positions are then on the complete grid the
subset was cut from, for the reason [`halo`](@ref) gives.

The motivating read is a stencil margin against position-indexed storage:

    for p in halo_positions(sys, chunk, l)
        margin[p] = source[p]
    end

which is why this is a verb and not `cellposition.(subtree_halo(...))`. That
expression materialises a `Vector` of ids in order to build a `Vector` of
positions and then discards the first; this streams, in the walk's own
`O(depth)` state.

WHAT IT COSTS, per cell, on top of the walk. Nothing worth measuring on HEALPix,
S2 and ISEA4R, where a position is the id plus one. Strictly less than nothing on
A5, whose scan holds the position already and would otherwise pay a Hilbert
decode to recover it. A real conversion on IGeo7 and H3, whose canonical ids are
not dense — under 10% of what those systems' own halo walks cost per cell, but
not free, and named here rather than glossed.
"""
halo_positions(it::Union{SubtreeHaloIterator,SubsetHaloIterator}) =
    HaloPositionIterator(it, _halo_grid(it))

halo_positions(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        l::Integer; connectivity::Connectivity = Vertex()) =
    halo_positions(SubtreeHaloIterator(sys, c, l; connectivity))

# The conversion, dispatched on the ENGINE so each walk pays only what its own
# state leaves undone. This is the default — one `cellposition`. The one
# override lives beside `ScanHaloEngine`, whose state already IS the position.
@inline _halo_position(::Any, grid, x, state) = cellposition(grid, x)::Int

Base.eltype(::Type{<:HaloPositionIterator}) = Int
Base.IteratorSize(::Type{<:HaloPositionIterator{I,G}}) where {I,G} =
    Base.IteratorSize(I)

# Delegated for its wrapped walk's reason: an engine that cannot count without
# walking defines no `length`, and reading positions instead of ids does not
# make it countable.
Base.length(it::HaloPositionIterator) = length(it.halo)

Base.collect(it::HaloPositionIterator) = collect_subtree(it, halo_sizehint(it))

Base.show(io::IO, it::HaloPositionIterator) =
    print(io, "halo_positions(", it.halo, ")")

function Base.iterate(it::HaloPositionIterator)
    r = iterate(it.halo)
    r === nothing && return nothing
    x, s = r
    return (_halo_position(it.halo.engine, it.grid, x, s), s)
end

function Base.iterate(it::HaloPositionIterator, state)
    r = iterate(it.halo, state)
    r === nothing && return nothing
    x, s = r
    return (_halo_position(it.halo.engine, it.grid, x, s), s)
end

# ===========================================================================
# An approximate size, which is deliberately not a `length`
# ===========================================================================

"""
    halo_sizehint(it) -> Union{Int,Nothing}

An APPROXIMATE upper bound on how many cells `it` yields, or `nothing` where no
bound is known. For `sizehint!` and for nothing else.

    h = halo_sizehint(it)
    out = eltype(it)[]
    h === nothing || sizehint!(out, h)
    for x in it; push!(out, x); end

WHY THIS IS NOT A `length`, AND MUST NOT BECOME ONE. `length` is a promise a
caller may allocate against and index into: `collect`'s sized route fills a
vector of exactly that many slots by iterating, so a `length` that is one too
large hands back an `undef` slot as a cell id. Five of the seven engines refuse
to declare one for exactly that reason (see [`SubtreeHaloIterator`](@ref)), and
`collect_subtree`'s miscount guard is meaningful only while they do. A HINT may
be wrong in both directions and costs at most one reallocation either way, so it
can answer where a count cannot — and it is threaded into `collect` as a
`sizehint!` argument, never through [`Base.IteratorSize`](@ref). Refusing every
size estimate as well as every count is what made `collect` of a `SizeUnknown`
halo allocate several times its own answer while growing.

WHERE THE NUMBERS COME FROM. Every one of them is measured by
`test/systems/crosssystem/subtree_halos.jl`, and the sweep is described here
rather than summarised as a case count so that the description stays true when
the sweep is widened:

  - the one-ring at depth zero, and the in-face square band, are EXACT — both
    engines already declare a `length` and this returns it;
  - the seam square band is `4·side + 8`. Over WHOLE GENERATIONS of all three
    aperture-4 systems at bases 0-2 and depths 1-3, both connectivities, the
    count is `4·side` under [`Edge`](@ref) and `4·side + e` under
    [`Vertex`](@ref) with `e ∈ {0, 2, 3, 4, 5}`; `e == 5` is ISEA4R and only
    ISEA4R, at a block whose corner sits at icosahedral vertex 0 or 11, where
    five diamonds meet and the corner contributes two cells rather than one.
    It is still exactly 5 at depth 12, i.e. at side 4096, so it is a fact about
    the vertex and not about the perimeter. `+8` is two cells per corner, which
    is the structural reading of that observation rather than the fitted
    constant;
  - the hexagonal walks are `3^(d+1) + 3` for `d = l - level(c)`, which is the
    census around a hexagon exactly, and about 20% over around a pentagon
    (`5(3^d + 1)/2`). Sampled across both systems at bases 0-1 and depths 1-4,
    every halo has one of those two counts and none exceeds the hint;
  - [`ScanHaloEngine`](@ref) and [`OutsideWalkEngine`](@ref) answer `nothing`.
    A5's subtree perimeter does not follow the aperture-4 formula — measured 12,
    22 and 42 at depths 1-3 from a level-1 root against `4·2^d + 5` of 13, 21
    and 37, under at one depth and over at two — and the generic walk serves
    subsets, whose halo has no perimeter at all. A guess there would be the one
    kind of wrong this file does not ship.
"""
halo_sizehint(it::SubtreeHaloIterator) = _halo_sizehint(it.engine)
halo_sizehint(it::SubsetHaloIterator) = _halo_sizehint(it.engine)
halo_sizehint(it::HaloPositionIterator) = halo_sizehint(it.halo)

# Nothing is the answer for a walk with no bound, and it is the answer a walk
# gets by DEFAULT — `ScanHaloEngine` and `OutsideWalkEngine` land here, and so
# does an engine added later, which is unbounded until someone measures it. The
# four bounded engines answer beside their own definitions below, the way each
# one's `IteratorSize` and emit rule already do.
_halo_sizehint(::Any) = nothing

# The one-ring's count is exact and already declared, so the hint IS the
# `length` rather than a second opinion about it.
_halo_sizehint(e::RingHaloEngine) = length(e)

# ===========================================================================
# The adjacency providers
#
# Candidate enumeration and adjacency testing are kept apart on purpose: the
# walk below decides WHICH cells to consider, a provider decides whether one of
# them touches the subtree. The indexed provider is the DGGS's own one-ring; the
# geometry provider compares unit-sphere boundaries and shares no topology with
# it, which is what makes it an oracle rather than a second opinion.
# ===========================================================================

"""
    IndexedNeighbors()

Adjacency by the system's native one-ring: `x` touches the subtree when one of
its neighbours has `root` as its ancestor. `O(degree · depth)` per candidate and
allocation-free wherever `neighbors` is.
"""
struct IndexedNeighbors end

"""
    ForcedGeometry()

Adjacency by unit-sphere boundary comparison: `x` touches the subtree when its
boundary shares a vertex ([`Vertex`](@ref)) or two ([`Edge`](@ref)) with the
boundary of some target-level descendant of `root`.

Descendants are visited through a second pruned hierarchy walk, so this costs
`O(depth)` memory rather than the subtree. It shares no index arithmetic with
[`IndexedNeighbors`](@ref) and is the oracle every specialized engine is tested
against — never `neighbors`, and never [`subtree_border`](@ref), which do share
it.

Comparison is against **target-level descendant** boundaries, not `root`'s own
polygon: H3, IGeo7 and A5 descendants can overhang their parent.
"""
struct ForcedGeometry end

"""
    SubsetMembership(subset, complete)

Adjacency to a SUBSET rather than to a subtree: `x` is a halo cell when the
subset does not hold it and does hold one of its neighbours,

    cellposition(subset, x) === nothing &&
        any(nb -> cellposition(subset, nb) !== nothing,
            neighbors(complete, x, 1; connectivity))

Unlike the two providers above, this predicate is COMPLETE — it decides
outside-ness as well as contact. It has to be: a subset has no descendant range
for `_admit`'s skip to retire the subject by, and no ancestor for
[`ScanHaloEngine`](@ref) to compare against. Folding both halves in here is also
what makes a HOLE fall out: an absent interior cell fails the first test's
negation nowhere and passes the second, so it is emitted like any other outside
cell.

`complete` is `levelgrid(system, level)` — the grid the neighbours come from,
since a subset's own `neighbors` is already clipped to membership and would hide
exactly the cells being looked for. The subset supplies membership only, which
is `O(log #windows)` on a [`CellVector`](@ref) and `O(log #cells)` on a
[`PartialGrid`](@ref) — plus [`subset_span`](@ref), which is that same search
asked about a whole block at once.
"""
struct SubsetMembership{S,G}
    subset::S
    complete::G
end

# ---------------------------------------------------------------------------
# The subset's shape, read as position spans
# ---------------------------------------------------------------------------

const _SPAN_NONE = 0
const _SPAN_SOME = 1
const _SPAN_ALL = 2

"""
    subset_span(subset, lo::Int, hi::Int) -> Int

How much of the position block `lo:hi` of the subset's own complete level the
subset holds: `_SPAN_NONE`, `_SPAN_SOME` or `_SPAN_ALL`.

The block is always a node's [`descendant_range`](@ref), so this is the
membership search asked about a whole subtree at once instead of cell by cell —
`O(log #windows)` on a [`CellVector`](@ref) and `O(log #cells)` on a
[`PartialGrid`](@ref), the same searches [`cellposition`](@ref) runs, and never
`O(hi - lo)`. It is what lets the subset walk retire a chunk's whole interior
for the price of one lookup; see `_admit(::SubsetMembership, ...)`.

Both containers store cells in ascending position order, so "the subset holds
every position of `lo:hi`" is decided by the two ENDPOINTS plus a count: the
stored entries between them number `hi - lo + 1` only if none is missing.
"""
function subset_span end

# --- the indexed test -------------------------------------------------------

# The halo's own definition, as one free function: `x` is a halo cell of `root`'s
# subtree exactly when the system's one-ring puts a descendant of `root` next to
# it. Every engine in this file that filters candidates ends here — the outside
# walk through the provider below, the seam band through `NativeCheck`, the
# hexagonal walk through both its emit rule and its calibration — so there is one
# definition to read and one to get wrong.
@inline function _touches_root(sys, grid, root, rootlevel::Int, x,
        connectivity::Connectivity)
    for nb in neighbors(grid, x, 1; connectivity)
        ancestor(sys, nb, rootlevel) == root && return true
    end
    return false
end

@inline _touches_subtree(::IndexedNeighbors, e, x) =
    _touches_root(e.system, e.grid, e.root, e.rootlevel, x, e.connectivity)

# --- the subset test --------------------------------------------------------

# The same definition with membership in place of ancestry, and it carries the
# outside half itself — see `SubsetMembership`.
@inline function _touches_subtree(p::SubsetMembership, e, x)
    cellposition(p.subset, x) === nothing || return false
    for nb in neighbors(p.complete, x, 1; connectivity = e.connectivity)
        cellposition(p.subset, nb) === nothing || return true
    end
    return false
end

# --- the geometry test ------------------------------------------------------

# The same rule the generic geometric `neighbors` fallback uses
# (`src/fallbacks/locate.jl`): shared corners within a tolerance scaled from the
# candidate's own shortest edge. One shared vertex is contact, two is an edge.
@inline _needed_contacts(::Vertex) = 1
@inline _needed_contacts(::Edge) = 2

function _touches_subtree(::ForcedGeometry, e, x)
    xb = cell_boundary(e.grid, x)
    tol = _match_tolerance(xb)
    needed = _needed_contacts(e.connectivity)
    xcap = cell_cap(e.grid, x)
    # Depth zero has no descent to make: `root` IS its own only target-level
    # descendant, so the comparison is against the root's own boundary. Handled
    # here rather than left to the cursor because a cursor seeded at the target
    # level has no children to expand — it would descend past the target to
    # `max_level` and throw.
    if e.rootlevel == e.target
        intersects_cap(cell_cap(e.grid, e.root), xcap) || return false
        return _shared_vertices(xb, cell_boundary(e.grid, e.root), tol) >= needed
    end
    return _descendant_touches(e, xb, tol, needed, xcap)
end

# A lazy cursor over `root`'s target-level descendants, pruned by the candidate's
# own cap: a node whose extent misses `xcap` holds no descendant that can touch
# `x`. `O(depth)` memory, and the descendants are never materialised.
function _descendant_touches(e, xb, tol::Float64, needed::Int, xcap)
    sys = e.system
    st = _empty_walk_stack(e.root)
    st = Helpers.small_push(st, HaloFrame(e.root, 0x1))
    while !isempty(st)
        k = length(st)
        f = @inbounds st[k]
        # Defence in depth behind the depth-zero early return above: a frame at
        # the target level is a leaf of THIS walk, whatever put it there.
        # `children` of a `max_level` cell throws, so never ask.
        if level(f.cell) == e.target
            st = Helpers.small_pop(st)
            continue
        end
        kids = children(sys, f.cell)
        if f.next > length(kids)
            st = Helpers.small_pop(st)
            continue
        end
        d = @inbounds kids[f.next]
        st = Helpers.small_setlast(st, HaloFrame(f.cell, f.next + 0x1))
        if level(d) == e.target
            intersects_cap(cell_cap(e.grid, d), xcap) || continue
            _shared_vertices(xb, cell_boundary(e.grid, d), tol) >= needed &&
                return true
            continue
        end
        intersects_cap(node_extent(sys, d), xcap) || continue
        st = Helpers.small_push(st, HaloFrame(d, 0x1))
    end
    return false
end

# ===========================================================================
# The outside-first walk
# ===========================================================================

# Why outside-first at all: a descendant-first algorithm (walk the rim, take
# each rim cell's neighbours, drop the ones still inside) produces duplicates
# out of order, and therefore needs a seen-set and a final sort — both sized by
# the halo. Reversing the search removes both. Each outside cell is considered
# at most once, so there is nothing to deduplicate and nothing to sort.
#
# Why the cap prune is sound: a halo cell `x` shares at least a point with some
# descendant `d` of `root`. `node_extent(sys, root)` provably contains every
# boundary point of every descendant of `root`, so that shared point lies in it;
# `node_extent(sys, n)` likewise contains `x`'s boundary for any ancestor `n` of
# `x`. So the two caps intersect, and a node whose cap misses the root cap can
# contain no halo cell.
#
# THE ASSUMPTION THAT COUPLES THE TWO HALVES OF THIS FILE. That argument is
# GEOMETRIC — it is the covering law over boundary POINTS — but the default
# provider's adjacency is TOPOLOGICAL, the system's native one-ring. The step
# "a halo cell shares at least a point with some descendant" is therefore an
# assumption about `neighbors`, not a theorem: it holds only while every pair of
# cells the system calls adjacent has boundaries that share a point. No bundled
# system violates it (checked by the forced-geometry agreement testset). A
# system that did — an adjacency defined by index arithmetic across a seam with
# no shared drawn corner — would keep the prune's arithmetic intact and its
# soundness not, and the walk would silently drop that halo cell. Such a system
# needs its own `halo_engine`, or a rootcap widened to cover the discrepancy.
#
# Why not lon/lat: longitude/latitude boxes are unusable at seams and poles, so
# everything here stays in unit-sphere XYZ.
#
# ---------------------------------------------------------------------------
# THE SUBSET WALK PRUNES BY A DIFFERENT LAW, AND IT IS NOT GEOMETRIC
#
# `SubsetMembership` has no root, so there is no root cap to compare against and
# no covering law to lean on. What it has instead is the subset's own position
# spans, and the prune those support is the COARSE-CONTAINMENT LAW:
#
#     for every pair of cells `x`, `y` that the system calls VERTEX-adjacent at
#     level `l`, `parent(y)` is `parent(x)` or a vertex-neighbour of it.
#
# Applied down the generations it gives the statement the walk actually uses: a
# neighbour of a level-`t` cell has its level-`lc` ancestor inside the CLOSED
# one-ring of that cell's own level-`lc` ancestor. So if neither a node nor any
# level-`lc` neighbour of it holds a member, no target-level descendant of it can
# have a member neighbour, and none of them is a halo cell. That is what
# `_near_subset` tests, and it is why the subset walk visits nodes in proportion
# to the subset's BOUNDARY rather than to the area of a bounding cap.
#
# THE SAME KIND OF ASSUMPTION THE CAP PRUNE MAKES, and it fails the same way:
# loudly in the differential tests, silently in production. A system whose
# refinement moved a cell out from under its parent's neighbourhood — an
# adjacency defined by index arithmetic that jumps a generation — would keep this
# arithmetic intact and its soundness not, and the walk would drop that halo
# cell. Every bundled system obeys it EXHAUSTIVELY: zero violations over every
# adjacent pair at levels 1 through 6 on HEALPix, S2, ISEA4R and A5, 1 through 5
# on IGeo7 and 1 through 4 on H3, which is what
# `test/systems/crosssystem/subtree_halos.jl` re-runs on every suite. A system
# that does not obey it must supply its own subset engine.
#
# THE PROBE IS `Vertex()` WHATEVER WAS REQUESTED, for `_seam_band_engine`'s
# reason: the `Edge()` halo is a subset of the `Vertex()` one, so a superset
# derived under `Vertex()` covers both, and the law above is stated — and
# tested — under `Vertex()` only. NOTHING HERE RELIES ON AN `Edge()` READING OF
# IT, and nothing should: narrowing the probe to the requested connectivity would
# be assuming a second law this file has not written down and the testset does
# not check, and a system that obeyed the `Vertex()` form while breaking the
# `Edge()` one would then lose halo cells silently under `Edge()` alone.

# One frame per level strictly above the target, so a full-depth walk from the
# root generation pushes at most `max_level` of them — 30 on S2, the deepest
# registered system. 34 is that plus four spare.
const _HALO_STACK_CAP = 34

struct HaloFrame{C}
    cell::C
    next::UInt8
end

const HaloStack{C} = Helpers.SmallList{_HALO_STACK_CAP,HaloFrame{C}}

@inline _empty_walk_stack(c::C) where {C} =
    Helpers.empty_small_list(Val(_HALO_STACK_CAP), HaloFrame(c, 0x0))

"""
    OutsideWalkEngine(system, grid, root, rootlevel, target, lo, hi, rootcap,
                      roots, provider, connectivity)

The correctness fallback: outside cells in canonical hierarchy order, the
subject subtree skipped whole, nodes pruned by cap, survivors tested by
`provider`.

Each outside cell is considered at most once, so there is no seen-set and no
final sort; the walk state is one `O(depth)` frame stack of `(cell, next child)`
pairs, and the geometry provider's descendant cursor is a second one. This path
may do more work than an indexed specialization — it is the correctness
fallback, and with [`ForcedGeometry`](@ref) it is the independent oracle those
specializations are checked against.

Requires more than [`has_sorted_subtrees`](@ref), which promises only that a
subtree's target-level descendants are CONTIGUOUS in position. This walk emits
in the order it meets cells, with no sort to repair it, so it additionally
requires that `children(sys, c)` and `rootcells(sys)` are each ordered by their
elements' TARGET-LEVEL descendant ranges — sibling `i` before sibling `j`
exactly when `first(descendant_range(sys, kids[i], target)) <
first(descendant_range(sys, kids[j], target))`, at every level and for every
target. Contiguity without that ordering still produces contiguous blocks, just
visited out of order, and the walk would emit a mis-sorted halo with no error
raised anywhere. Every bundled system satisfies it, and
`test/systems/crosssystem/subtree_halos.jl` is what says so: its law compares
this walk element for element against an ascending-POSITION scan of the target
level. A system that does not satisfy it must supply its own `halo_engine`
rather than inherit this one.
"""
struct OutsideWalkEngine{S,G,C,R,P,K}
    system::S
    grid::G
    root::C
    rootlevel::Int
    target::Int
    lo::Int
    hi::Int
    rootcap::Cap
    roots::R
    provider::P
    connectivity::K
end

Base.eltype(::Type{<:OutsideWalkEngine{S,G,C,R,P,K}}) where {S,G,C,R,P,K} = C
Base.IteratorSize(::Type{<:OutsideWalkEngine}) = Base.SizeUnknown()

struct OutsideWalk{C}
    top::Int
    stack::HaloStack{C}
end

const _HALO_SKIP = 0
const _HALO_EMIT = 1
const _HALO_DESCEND = 2

# How expensive a cap prune is BEFORE the adjacency test, per SUBTREE provider —
# the subset arm of `_admit` never asks, because its adjacency test opens with a
# `cellposition` that retires an in-set candidate for less than a cap costs. The
# indexed provider is one native `neighbors` call, cheaper than building a cap
# from a boundary, so it declines too: at the target level it goes straight to
# the test. The forced-geometry provider walks the ROOT's subtree once per
# candidate, so a single cap comparison that retires the candidate outright is
# worth many times its cost — and this is the oracle path, which the
# differential tests hammer. Dispatch on the provider so the indexed path pays
# nothing for the distinction, not a branch on a field.
#
# Sound for the same reason the internal-node prune is: a halo cell shares a
# boundary point with a descendant of `root`, that point is inside `rootcap` by
# the covering law, and it is one of the candidate's own corners, hence inside
# the candidate's cap.
@inline _target_prune(::IndexedNeighbors, e, c) = true
@inline _target_prune(::ForcedGeometry, e, c) =
    intersects_cap(cell_cap(e.grid, c), e.rootcap)

# ONE VERDICT PER NODE, AND THE PROVIDER DECIDES IT. The subtree providers ask a
# question about ANCESTRY and prune by geometry; `SubsetMembership` asks one
# about MEMBERSHIP and prunes by the subset's own position spans. Those are
# different rules, not one rule with a flag, so `_admit` dispatches on the
# provider exactly as `_touches_subtree` and `_target_prune` do — each arm stays
# monomorphic and neither pays for the other's fields.
@inline _admit(e::OutsideWalkEngine, c) = _admit(e.provider, e, c)

# --- the subtree arm --------------------------------------------------------

# Three questions, cheapest first. A node whose whole descendant range sits
# inside the subject's is the subject subtree itself (ranges nest or are
# disjoint), and integer comparison retires it without touching geometry.
#
# That containment can only hold AT the root's own level, which is why the guard
# is `==` and not `>=`. Deeper than the root there are two cases and neither can
# be contained: a node that is one of the root's own descendants is unreachable,
# because the root was already skipped whole and the walk never descended into
# it; and a node that is not has a descendant range disjoint from the root's, so
# the containment test is a `descendant_range` call that cannot succeed. On
# IGeo7 with a level-2 root at `l = 5` that is 630 of the 689 nodes the walk
# visits.
@inline function _admit(p::Union{IndexedNeighbors,ForcedGeometry},
        e::OutsideWalkEngine, c)
    lc = level(c)
    if lc == e.rootlevel
        r = descendant_range(e.system, c, e.target)
        (first(r) >= e.lo && last(r) <= e.hi) && return _HALO_SKIP
    end
    if lc == e.target
        _target_prune(p, e, c) || return _HALO_SKIP
        return _touches_subtree(p, e, c) ? _HALO_EMIT : _HALO_SKIP
    end
    intersects_cap(node_extent(e.system, c), e.rootcap) || return _HALO_SKIP
    return _HALO_DESCEND
end

# --- the subset arm ---------------------------------------------------------

# Two integer questions and no geometry at all — see `subset_halo_engine` for
# why the cap is gone and what replaces it.
#
#   * a node the subset holds ENTIRELY holds no halo cell, because a halo cell
#     is by definition one the subset does not hold. One `subset_span` retires
#     the whole block, which is what makes a chunk's own interior free;
#   * a node the subset touches at all must be descended, because the boundary
#     between held and unheld cells is inside it;
#   * a node the subset does not touch is descended only if a level-`lc`
#     NEIGHBOUR of it does — the coarse-containment law below.
@inline function _admit(p::SubsetMembership, e::OutsideWalkEngine, c)
    lc = level(c)
    if lc == e.target
        return _touches_subtree(p, e, c) ? _HALO_EMIT : _HALO_SKIP
    end
    r = descendant_range(e.system, c, e.target)
    s = subset_span(p.subset, first(r), last(r))
    s == _SPAN_ALL && return _HALO_SKIP
    s == _SPAN_SOME && return _HALO_DESCEND
    return _near_subset(p, e, c, lc) ? _HALO_DESCEND : _HALO_SKIP
end

# The one-ring of `c` at `c`'s OWN level, asked of the subset's spans. The probe
# is `Vertex()` whatever connectivity was requested, for the seam band's reason:
# the `Edge()` halo is a subset of the `Vertex()` one, so one conservative
# superset serves both, and the requested connectivity is still what
# `_touches_subtree` filters by at the target level.
@inline function _near_subset(p::SubsetMembership, e::OutsideWalkEngine, c, lc::Int)
    coarse = levelgrid(e.system, lc)
    for nb in neighbors(coarse, c, 1; connectivity = Vertex())
        r = descendant_range(e.system, nb, e.target)
        subset_span(p.subset, first(r), last(r)) == _SPAN_NONE || return true
    end
    return false
end

Base.iterate(e::OutsideWalkEngine{S,G,C}) where {S,G,C} =
    iterate(e, OutsideWalk(0, _empty_walk_stack(e.root)))

function Base.iterate(e::OutsideWalkEngine{S,G,C}, w::OutsideWalk{C}) where {S,G,C}
    top = w.top
    st = w.stack
    while true
        if isempty(st)
            top += 1
            top > length(e.roots) && return nothing
            c = @inbounds e.roots[top]
            v = _admit(e, c)
            v == _HALO_SKIP && continue
            v == _HALO_EMIT && return (c, OutsideWalk(top, st))
            st = Helpers.small_push(st, HaloFrame(c, 0x1))
            continue
        end
        k = length(st)
        f = @inbounds st[k]
        kids = children(e.system, f.cell)
        if f.next > length(kids)
            st = Helpers.small_pop(st)
            continue
        end
        c = @inbounds kids[f.next]
        st = Helpers.small_setlast(st, HaloFrame(f.cell, f.next + 0x1))
        v = _admit(e, c)
        v == _HALO_SKIP && continue
        v == _HALO_EMIT && return (c, OutsideWalk(top, st))
        st = Helpers.small_push(st, HaloFrame(c, 0x1))
    end
end

# ===========================================================================
# The linear scan: canonical order without a descendant range
#
# A5 AND WHY IT STAYS HERE. This is the engine A5 gets, at every root and every
# depth, and it is the only system that takes it. The reason is one missing
# primitive: `has_sorted_subtrees(A5System())` is false and A5 has no
# `descendant_range` method at all, so `OutsideWalkEngine` has no integer range
# to skip the subject subtree by and no ordering to make a pruned descent
# canonical. Without those the honest walk is the scan: `O(ncells)` in time,
# `O(1)` in memory, canonical by construction because it IS the canonical order.
#
# A5 DOES have native indexed one-rings, so it is worth saying what would NOT
# justify a fast path. Its aperture and its Hilbert-like indexing are not
# evidence of one: a subtree's boundary is not a shared square perimeter the way
# HEALPix's, S2's and ISEA4R's are, and no validated directed-border automaton
# exists for it — `has_sorted_subtrees` being false is the same fact seen from
# the other side. `SquareBandEngine` would need `lattice_decode` / `lattice_cell`
# / `face_orientation` on a lattice A5 does not have, and `HexArcHaloEngine`
# would need `seeded_rim_engine` plus the descendant-range order that makes
# concatenating neighbours a merge. Adding either by analogy would produce a
# walk that is wrong in a way no test in this package currently asks about.
#
# A dedicated A5 engine belongs here only once two things are proved
# INDEPENDENTLY, in the sense Tasks 4-6 proved them for the other five systems:
# the boundary states of an A5 subtree (which cells of a neighbouring subtree can
# touch it, and in what order), and a two-sided `descendant_range` contract that
# makes those neighbours' streams concatenate into canonical order. Until then
# the scan is not a placeholder, it is the correct answer at the price the
# missing primitives set.
# ===========================================================================

"""
    ScanHaloEngine(system, grid, root, rootlevel, target, provider, connectivity)

Every cell of the target level in position order, the descendants skipped and
the rest tested. `O(1)` memory and canonical by construction, but `O(ncells)`
time — the price of a system with no [`descendant_range`](@ref) to prune by, and
A5 is the only one. See the comment above this type for what a dedicated A5
engine would have to prove first.
"""
struct ScanHaloEngine{S,G,C,P,K}
    system::S
    grid::G
    root::C
    rootlevel::Int
    target::Int
    provider::P
    connectivity::K
end

Base.eltype(::Type{<:ScanHaloEngine{S,G,C,P,K}}) where {S,G,C,P,K} = C
Base.IteratorSize(::Type{<:ScanHaloEngine}) = Base.SizeUnknown()

# Outside-ness, stated because the scan has nothing structural to lean on: it
# meets every cell of the level, the subject's own descendants included, where
# `_admit` would already have retired the subject subtree by integer range
# containment. `SubsetMembership` is the one provider whose `_touches_subtree`
# decides outside-ness itself, so it answers `true` here rather than paying for
# the membership search twice — and it must, since a subset has no ancestor to
# compare and a `rootlevel` below the shallowest level would throw.
@inline _scan_outside(::Union{IndexedNeighbors,ForcedGeometry}, e, x) =
    ancestor(e.system, x, e.rootlevel) != e.root
@inline _scan_outside(::SubsetMembership, e, x) = true

# The scan's state is the position it will resume at, so the cell just yielded
# sat one before it and `halo_positions` needs no `cellposition` at all. That is
# the largest per-cell saving the position walk makes anywhere: A5's
# `cellposition` is a Hilbert decode, a measured 38 ns against the 0.64 ns the
# three aperture-4 systems pay.
#
# The grid and the cell go unnamed for `_band_emit`'s reason — naming arguments
# this ignores would suggest it consults them, and the whole point is that it
# does not have to.
@inline _halo_position(::ScanHaloEngine, ::Any, ::Any, state::Int) = state - 1

Base.iterate(e::ScanHaloEngine) = iterate(e, 1)
function Base.iterate(e::ScanHaloEngine, p::Int)
    n = ncells(e.grid)
    while p <= n
        x = cellindex(e.grid, p)
        if _scan_outside(e.provider, e, x) && _touches_subtree(e.provider, e, x)
            return (x, p + 1)
        end
        p += 1
    end
    return nothing
end

# ===========================================================================
# Building either one
# ===========================================================================

function outside_walk_engine(sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex, target::Int, connectivity::Connectivity, provider)
    grid = levelgrid(sys, target)
    lc = level(c)
    has_sorted_subtrees(sys) ||
        return ScanHaloEngine(sys, grid, c, lc, target, provider, connectivity)
    r = descendant_range(sys, c, target)
    return OutsideWalkEngine(sys, grid, c, lc, target, first(r), last(r),
        node_extent(sys, c), rootcells(sys), provider, connectivity)
end

"""
    geometry_halo_engine(sys, c, target, connectivity)

The generic walk forced onto [`ForcedGeometry`](@ref), whatever fast path the
system would otherwise take. Not a public verb and not reachable from
[`SubtreeHaloIterator`](@ref)'s keyword constructor: it exists so a test can
build the oracle explicitly and hand it to the positional constructor.
"""
geometry_halo_engine(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        target::Int, connectivity::Connectivity) =
    outside_walk_engine(sys, c, target, connectivity, ForcedGeometry())

"""
    subset_halo_engine(sys, subset, complete, target, connectivity)

The engine behind [`halo`](@ref) for an arbitrary subset: the same outside-first
walk the subtree fallback uses, with [`SubsetMembership`](@ref) in place of the
subtree providers and with the SUBSET'S OWN SPANS in place of the subtree's
root range and bounding cap.

WHY THE SUBTREE'S SKIP IS SWITCHED OFF rather than pointed at the subset.
`_admit`'s subtree arm retires a node whose whole target-level range sits inside
the SUBJECT SUBTREE's, which is exactly what must not happen here: a cell punched
out of the middle of a subset lies in that range and IS a halo cell. Setting
`rootlevel` one below the system's shallowest level makes the branch unreachable
— no level compares equal to it — so nothing on this path can reach it. `lo`/`hi`
are an empty range so the skip would refuse even if some future edit made the
branch live again, `rootcap` is the full sphere so a cap test would prune
nothing, and `root` is a placeholder that fixes the frame stack's element type.
All four are unread on this path; the subset arm of `_admit` reads none of them.

WHAT IT PRUNES BY INSTEAD, and why there is no cap at all. A node is descended
only when the subset touches it or touches one of its level-`lc` neighbours, so
the walk follows the subset's BOUNDARY rather than filling a bounding region:
`O(halo)` nodes rather than `O(cap area)`, and construction is `O(1)` rather
than a boundary evaluation per member. See the coarse-containment law in the
section comment above [`OutsideWalkEngine`](@ref) for why the neighbour probe is
sound and what a system would have to violate to break it.

`O(depth)` memory beyond the subset's own storage. A system with no
[`descendant_range`](@ref) gets [`ScanHaloEngine`](@ref), for the same reason it
does at a subtree.

NOT AN EXTENSION POINT, which is why it is not declared in `interface/system.jl`
beside `rim_engine` / `interior_engine` / `halo_engine` and why its signature
does not follow their `(sys, c, target, connectivity)` shape. Those three are
dispatched on the SYSTEM type: a system ships a fast path by adding a method.
This one has no cell to dispatch a subtree walk from — its subject is a
membership predicate over an arbitrary set — and it takes `complete` and `cap`
because the subset, not the system, is what decides them. A system that wants a
faster subset halo has nowhere to put it here; it makes `neighbors` and
`descendant_range` faster, and every subset in the package gets it. Should that
ever stop being true, the thing to add is a declared method on the system, not a
fourth shape of this one.
"""
function subset_halo_engine(sys::AbstractHierarchicalGridSystem, subset,
        complete::AbstractGrid, target::Int, connectivity::Connectivity)
    provider = SubsetMembership(subset, complete)
    seed = first(rootcells(sys))
    nolevel = first(levels(sys)) - 1
    has_sorted_subtrees(sys) ||
        return ScanHaloEngine(sys, complete, seed, nolevel, target, provider,
            connectivity)
    return OutsideWalkEngine(sys, complete, seed, nolevel, target, 1, 0,
        full_sphere_cap(), rootcells(sys), provider, connectivity)
end

# ===========================================================================
# The exterior-perimeter walk, shared by the three aperture-4 systems
# ===========================================================================

# A subtree of HEALPix, S2 or ISEA4R is an aligned `side x side` block in one
# face's (diamond's) square lattice, and its halo is the width-1 band around
# it. Two cases, one walk.
#
# AWAY FROM THE FACE EDGE the band is entirely in-face, where adjacency is the
# plain 3×3 lattice on all three systems, so band and halo are the SAME SET —
# with the four diagonal-contact corners dropped under `Edge()`. Nothing needs
# checking and nothing needs a seam table.
#
# ACROSS A SEAM the halo leaves the face: cells on up to six other faces, whose
# ids no in-face arithmetic can name. They are reached by generalising the walk
# from one band to a short, FACE-ORDERED list of candidate rectangles, one per
# face the halo can touch — see `_seam_band_engine` for how those rectangles are
# derived and why they cover. That list is a conservative SUPERSET, so every
# candidate is put through the native one-ring before it is yielded.
#
# WHY CONCATENATION IS A MERGE. Each face occupies the contiguous id range
# `[face·n², (face+1)·n²)` and faces are numbered ascending, so global canonical
# order is (face ascending, in-face curve code ascending). Walking the rectangle
# list in face order is therefore already a canonical merge of the face-local
# streams: no heap, no seen-set, no sort. Merging same-face rectangles into
# their bounding rectangle is what makes it one stream per face, which is what
# makes it impossible to emit a cell twice.
#
# WHY A QUADTREE DESCENT AND NOT AN OFFSET WALK. The band is not a curve
# interval of any one subtree, so there is no `base + offset` to run. The engine
# descends each rectangle's whole FACE in curve order and prunes every quadrant
# that misses the rectangle: curve order at each level is ascending id by
# construction, so survivors arrive ascending. A quadrant of side `h` survives
# only if it meets the rectangle, and at most `O(perimeter/h + 1)` of each size
# do, so the walk visits `O(halo + depth)` nodes in `O(depth)` memory — the
# descent tracks the running curve code and lattice origin and restores both on
# pop by re-reading the parent's last step.
#
# DEPTH ZERO IS NOT HERE. `side == 1` routes to `RingHaloEngine` instead: a
# one-cell band is the plain eight-cell lattice neighbourhood, which is wrong at
# every seam (a cube-corner cell has seven neighbours, a HEALPix degree-3 vertex
# seven, an ISEA4R cell at icosahedral vertex 0 or 11 nine), and the native
# one-ring answers it exactly in one call.

# `mask` would be dead weight here — pruning is by lattice overlap, not by
# flush sides — so the frame is the two bytes the descent actually needs.
# `HaloFrame` above is the outside walk's frame and carries a cell; this one
# carries a curve state, so the two cannot share a name.
struct SquareBandFrame
    orientation::UInt8
    next::UInt8
end

const SquareBandStack = Helpers.SmallList{_SQUARE_CAP,SquareBandFrame}

@inline _empty_band_stack() =
    Helpers.empty_small_list(Val(_SQUARE_CAP), SquareBandFrame(0x0, 0x0))

"""
    FaceRect(face, orientation, x0, y0, x1, y1)

One face's candidate rectangle: the inclusive lattice box `[x0, x1] x [y0, y1]`
on 0-based `face`, to be descended under curve state `orientation` (the state
that face's ROOT is read under, from [`face_orientation`](@ref)).

The rectangles of a [`SquareBandEngine`](@ref) are one per face and sorted by
`face`, which is what makes walking them a canonical merge.

`Int32` BOUNDS BIND AT LEVEL 32, NOT AT `max_level`. A level-`l` lattice
coordinate runs to `2^l - 1`, so `Int32` holds one through level 31
(`2^31 - 1 == typemax(Int32)`) and overflows at level 32. S2's `max_level` of 30
is the deepest registered system, so there is exactly ONE level of headroom, and
the quantity to compare a future `max_level` bump against is 31 — not 30, and
not `_SQUARE_CAP`. Past it the failure is an `InexactError` raised by this
constructor from inside `square_halo_engine`, i.e. from iterator construction,
which is loud but says nothing about the cause; widen these six fields to
`Int64` (they are `Int32` only to keep `_BAND_RECT_CAP` rectangles inline and
cheap to copy) rather than clamping. `test/systems/crosssystem/subtree_halos.jl`
walks a `max_level` block on all three systems, so the level-31 boundary is
approached from one level below on every run.
"""
struct FaceRect
    face::Int32
    orientation::UInt8
    x0::Int32
    y0::Int32
    x1::Int32
    y1::Int32
end

# Twelve is HEALPix's face count and the largest of the three (ISEA4R has ten
# diamonds, S2 six faces); one rectangle per face is the hard ceiling, because
# same-face rectangles are merged. Measured worst case is seven.
const _BAND_RECT_CAP = 12
const BandRects = Helpers.SmallList{_BAND_RECT_CAP,FaceRect}

@inline _empty_band_rects() = Helpers.empty_small_list(Val(_BAND_RECT_CAP),
    FaceRect(0, 0x0, 0, 0, 0, 0))

# Merge by face, so each face appears once and no cell can be reached twice.
#
# `r.orientation` is DISCARDED and the incumbent's kept, which is correct because
# the branch is only taken when `q.face == r.face` and both orientations came
# from `face_orientation(sys, face)` — a pure function of the face alone, with no
# dependence on the cell, the level or the rectangle. So the two values are
# necessarily equal and the merge has no choice to make. That is an invariant of
# the `face_orientation` contract (`src/interface/system.jl`): a system whose
# orientation varied within a face would break the whole descent, not only this
# line, because `SquareBandEngine` seeds each rectangle's descent at that face's
# ROOT and reads no per-cell state at all.
function _merge_rect(rects::BandRects, r::FaceRect)
    for i in 1:length(rects)
        q = @inbounds rects[i]
        q.face == r.face || continue
        return Helpers.small_setindex(rects, FaceRect(q.face, q.orientation,
            min(q.x0, r.x0), min(q.y0, r.y0),
            max(q.x1, r.x1), max(q.y1, r.y1)), i)
    end
    return Helpers.small_push(rects, r)
end

# ---------------------------------------------------------------------------
# Checked and unchecked, as a type rather than a flag
#
# The in-face band is EXACT — proved by the interval guard in
# `square_halo_engine`, and the reason no `neighbors` call appears on that path.
# The seam band is a conservative superset and every candidate owes the native
# test. Carrying the difference as the engine's second type parameter keeps both
# paths monomorphic: the exact path's emit rule inlines to the corner test and
# the checked path's to the one-ring, with no branch on a field in either.
# ---------------------------------------------------------------------------

"""
    NoCheck()

The emit rule of an exact band: a surviving leaf IS a halo cell, subject only to
`Edge()` dropping the four diagonal corners. Zero-size, so an engine carrying it
costs nothing for the distinction.
"""
struct NoCheck end

"""
    NativeCheck(system, grid, root, rootlevel, connectivity)

The emit rule of a conservative band: a candidate is a halo cell only if the
system's own one-ring puts a descendant of `root` next to it,

    any(nb -> ancestor(sys, nb, rootlevel) == root,
        neighbors(grid, x, 1; connectivity))

which is the exact definition, applied to every candidate before it is yielded.
`O(degree · depth)` per candidate.
"""
struct NativeCheck{S,G,C,K}
    system::S
    grid::G
    root::C
    rootlevel::Int
    connectivity::K
end

"""
    SquareBandEngine(curve, check, level, faceside, homeface, x0, y0, side, corners, rects)

The halo of the `side x side` block at lattice origin `(x0, y0)` of `homeface`,
in ascending id, by one pruned quadtree descent per rectangle in `rects` — which
are one per face and ascending by face, so the concatenation is already the
canonical merge.

`faceside` is a face's full lattice side at `level`. Yields [`LevelIndex`](@ref)
on [`SquareRimEngine`](@ref)'s reasoning and takes the same
[`quadrant_step`](@ref) curves. `O(candidates + depth)` time, `O(depth)` memory.

`check` decides both the emit rule and the count contract:

  - [`NoCheck`](@ref) — the block is nowhere flush with its face edge, the one
    rectangle is the width-1 band, and band and halo are the same set. The count
    is closed form, `4·side + 4` or `4·side` under `Edge()`, so
    [`Base.IteratorSize`](@ref) is `HasLength()`.
  - [`NativeCheck`](@ref) — the block is flush somewhere, the rectangles are a
    conservative superset, and every candidate is tested before it is yielded.
    No perimeter formula survives a seam (a cube corner is three cells, not
    four; an ISEA4R icosahedral vertex is five), so `IteratorSize` is
    `SizeUnknown()` and there is **no `length` method at all**.

SIZE, WHICH IS NOT FREE EVEN WHERE TIME IS. `rects` is a fixed
`_BAND_RECT_CAP`-slot inline list, so an engine is 352 bytes on the in-face path
(`SquareBandEngine{MortonCurve,NoCheck}`), 296 of them the rectangle list with
exactly one slot used. Nothing is heap-allocated and no measured time is spent on
the unused slots — the descent reads `length(rects)`, never the capacity — but
the in-face path is "free" in TIME only, and an engine returned by value costs
that copy. Shrinking it would mean a second engine type for the one-rectangle
case, which is not worth two monomorphic walks.
"""
struct SquareBandEngine{V,K}
    curve::V
    check::K
    level::Int
    faceside::Int64
    homeface::Int32
    x0::Int64
    y0::Int64
    side::Int64
    # `Vertex()`, i.e. keep the four diagonal-contact corners. Read only by
    # `_band_emit(::NoCheck, ...)`: under `NativeCheck` the one-ring is the
    # filter and this field is inert, since the connectivity the check applies
    # is the one stored on the `NativeCheck` itself.
    corners::Bool
    rects::BandRects
end

"""
    SquareBandWalk(rect, stack, code, x, y)

The walk state: which rectangle is being descended, the frame stack, and the
within-face curve code and lattice origin of the sub-square on top of it. Frame
`k`'s node has side `faceside >> (k - 1)`, so none of the last three needs
storing per frame — all are restored on pop by replaying the parent's last
[`quadrant_step`](@ref).
"""
struct SquareBandWalk
    rect::Int
    stack::SquareBandStack
    code::Int64
    x::Int64
    y::Int64
end

@inline _rect_overlaps(r::FaceRect, x::Int64, y::Int64, h::Int64) =
    x <= r.x1 && x + h - 1 >= r.x0 && y <= r.y1 && y + h - 1 >= r.y0

# The block is not its own halo. On the home face this prunes it whole — one
# integer test retires a quadrant of any size — and it is needed on BOTH paths:
# a block cell's neighbours include other block cells, so the native check would
# happily admit one.
@inline _inside_block(e::SquareBandEngine, x::Int64, y::Int64, h::Int64) =
    x >= e.x0 && x + h - 1 <= e.x0 + e.side - 1 &&
    y >= e.y0 && y + h - 1 <= e.y0 + e.side - 1

# The four cells that touch the block at a vertex only, halo under `Vertex()`
# and not under `Edge()`.
@inline _band_corner(e::SquareBandEngine, x::Int64, y::Int64) =
    (x == e.x0 - 1 || x == e.x0 + e.side) &&
    (y == e.y0 - 1 || y == e.y0 + e.side)

@inline _band_emit(::NoCheck, e::SquareBandEngine, cell, cx::Int64, cy::Int64) =
    e.corners || !_band_corner(e, cx, cy)

# The lattice position is the in-face band's business, not the one-ring's: the
# check is the halo's own definition and needs only the cell. Hence the three
# unnamed argument types — the signature exists to match `_band_emit`'s shape,
# and naming arguments it ignores would suggest it consults them.
@inline _band_emit(chk::NativeCheck, ::SquareBandEngine, cell, ::Int64, ::Int64) =
    _touches_root(chk.system, chk.grid, chk.root, chk.rootlevel, cell,
        chk.connectivity)

Base.eltype(::Type{<:SquareBandEngine}) = LevelIndex

Base.IteratorSize(::Type{<:SquareBandEngine{V,NoCheck}}) where {V} =
    Base.HasLength()
Base.IteratorSize(::Type{<:SquareBandEngine{V,K}}) where {V,K<:NativeCheck} =
    Base.SizeUnknown()

# Declared only for the exact band. A seam engine has no closed-form count, so
# the `MethodError` from here is the contract being kept — see the type's
# docstring.
Base.length(e::SquareBandEngine{V,NoCheck}) where {V} =
    e.corners ? Int(4 * e.side + 4) : Int(4 * e.side)

# The exact band's hint is its own count. The seam band's is the band plus two
# cells per corner: a seam corner can contribute a second cell where more than
# three faces meet, which is ISEA4R at icosahedral vertex 0 or 11 and nowhere
# else measured. See `halo_sizehint` for the sweep those two sentences come
# from, and note that a hint three cells generous of the worst case measured is
# a `sizehint!` and not a `length` — the count contract above is untouched.
_halo_sizehint(e::SquareBandEngine{V,NoCheck}) where {V} = length(e)
_halo_sizehint(e::SquareBandEngine) = Int(4 * e.side + 8)

Base.iterate(e::SquareBandEngine) =
    iterate(e, SquareBandWalk(0, _empty_band_stack(),
        Int64(0), Int64(0), Int64(0)))

function Base.iterate(e::SquareBandEngine, w::SquareBandWalk)
    r = w.rect
    st = w.stack
    code = w.code
    x = w.x
    y = w.y
    while true
        if isempty(st)
            # One rectangle finished; start the next face's descent at its root.
            r += 1
            r > length(e.rects) && return nothing
            st = Helpers.small_push(_empty_band_stack(),
                SquareBandFrame((@inbounds e.rects[r]).orientation, 0x0))
            code = Int64(0)
            x = Int64(0)
            y = Int64(0)
            continue
        end
        rect = @inbounds e.rects[r]
        k = length(st)
        f = @inbounds st[k]
        if f.next > 0x3
            st = Helpers.small_pop(st)
            if !isempty(st)
                # Undo the push that entered the frame just dropped: its parent
                # is back on top and has advanced past the position it used.
                pk = length(st)
                pf = @inbounds st[pk]
                p = Int(pf.next) - 1
                i, j, _ = quadrant_step(e.curve, pf.orientation, p)
                half = e.faceside >> pk
                code -= p * half * half
                x -= i * half
                y -= j * half
            end
            continue
        end
        p = Int(f.next)
        st = Helpers.small_setlast(st, SquareBandFrame(f.orientation, f.next + 0x1))
        i, j, o = quadrant_step(e.curve, f.orientation, p)
        half = e.faceside >> k
        cx = x + i * half
        cy = y + j * half
        _rect_overlaps(rect, cx, cy, half) || continue
        rect.face == e.homeface && _inside_block(e, cx, cy, half) && continue
        child = code + p * half * half
        if half == 1
            cell = LevelIndex(e.level,
                Int64(rect.face) * e.faceside * e.faceside + child)
            _band_emit(e.check, e, cell, cx, cy) || continue
            # Emitted without descending, so `code`/`x`/`y` still describe the
            # frame on top and need no restoring on the way back in.
            return (cell, SquareBandWalk(r, st, code, x, y))
        end
        st = Helpers.small_push(st, SquareBandFrame(o, 0x0))
        code = child
        x = cx
        y = cy
    end
end

# ===========================================================================
# Building one: the interval guard, and the seam derivation
# ===========================================================================

"""
    square_halo_engine(sys, curve, c, target, connectivity, x0, y0, side, face, n)

The halo engine for the `side x side` block at lattice origin `(x0, y0)` of
0-based `face`, on a face of side `n` at level `target`. The three aperture-4
systems' `halo_engine` methods are this call plus their own lattice decode.

Away from the face edge it is the exact width-1 band, unchecked and counted.
Flush with it, `_seam_band_engine` takes over. `side == 1` never reaches here.
"""
function square_halo_engine(sys::AbstractHierarchicalGridSystem, curve,
        c::AbstractCellIndex, target::Int, connectivity::Connectivity,
        x0::Int64, y0::Int64, side::Int64, face::Int64, n::Int64)
    home = FaceRect(face, face_orientation(sys, face),
        max(Int64(0), x0 - 1), max(Int64(0), y0 - 1),
        min(n - 1, x0 + side), min(n - 1, y0 + side))
    # The width-1 band fits inside the face: every halo cell is in-face, the
    # band IS the halo, and nothing below needs a neighbour query.
    if 1 <= x0 && x0 + side <= n - 1 && 1 <= y0 && y0 + side <= n - 1
        return SquareBandEngine(curve, NoCheck(), target, n, Int32(face),
            x0, y0, side, connectivity isa Vertex,
            Helpers.small_push(_empty_band_rects(), home))
    end
    return _seam_band_engine(sys, curve, c, target, connectivity,
        x0, y0, side, face, n, home)
end

# ---------------------------------------------------------------------------
# Deriving the candidate rectangles, with no seam table of this file's own
#
# WHAT HAS TO BE COVERED. A halo cell is a neighbour of a block cell. Block
# cells with a neighbour off the face are exactly the rim cells of the sides
# that are FLUSH with the face edge, and their off-face neighbours are the
# images of the extended-lattice positions one step outside that edge. So for a
# block flush on, say, `x = 0`, the foreign candidates are the images of
# `(-1, y')` for `y'` running over `[y0 - 1, y0 + side]` clipped to the face,
# plus — where a corner of the block is a corner of the face — the cells that a
# DIAGONAL step off two edges at once reaches.
#
# WHY TWO PROBES PER SIDE SUFFICE. All three systems' seam maps are edge-to-edge
# affine with a sign: S2's `wrap_xyf` computes `k = (σ·b + n - 1) >> 1` from the
# centred along-edge coordinate, which is `y` or `n - 1 - y`; ISEA4R's
# `lattice_neighbors` reads the paired rim slot at `n - 1 - j`; HEALPix's
# `nested_neighbors` applies the `NB_SWAPARRAY` mirrors and transpose. Each is
# monotone along the edge AT EVERY `n`, so the images of an interval of `y'` are
# a contiguous run lying between the images of its ENDPOINTS.
#
# The two extreme rim cells of a flush side see both endpoints directly: the
# extreme cell at `(0, y0)` has `(-1, y0 - 1)` among its eight neighbours, and
# the one at `(0, y0 + side - 1)` has `(-1, y0 + side)`. So the bounding box of
# what the two probes see already contains the image of every interior rim cell
# of that side, and nothing needs widening. That is not an argument this file
# takes on trust — see the verification note below, and note that the tests
# would fail if a bound were moved one cell inward.
#
# WHERE A THIRD FACE APPEARS. At a flush CORNER the block's corner cell is a
# probe of both flush sides, and its own neighbour list already contains
# whatever the diagonal step reaches: nothing on S2 (`wrap_xyf` returns
# `nothing` at a cube corner) or at HEALPix's `-1` entries in `NB_FACEARRAY` —
# which occur only for the double-out `nbnum`, i.e. only at a face corner, never
# at a run endpoint — and the interior of `CORNER_FANS` on ISEA4R, which is the
# two extra diamonds meeting at icosahedral vertex 0 or 11. So corners need no
# separate rule; they need only that the corner cell is probed, which flushness
# already guarantees.
#
# THE PROBE ASKS FOR `Vertex()` WHATEVER WAS REQUESTED, so one coverage argument
# serves both connectivities: the `Edge()` halo is a subset of the `Vertex()`
# one, and a superset of the larger is a superset of the smaller. The requested
# connectivity is still what `NativeCheck` filters by.
#
# COVERAGE IS EXHAUSTIVELY VERIFIED, not merely argued: every flush block of
# every size at every origin on every face, levels 1 through 6, both
# connectivities, on all three systems — zero halo cells outside the derived
# rectangles, worst case seven rectangles. Because the seam maps are affine at
# every `n`, levels past 6 add no new structure, only longer runs. The
# differential tests in `test/systems/crosssystem/subtree_halos.jl` re-run the
# same claim through the forced-geometry oracle, which is the only oracle in
# that file that can see a candidate this derivation never proposed.
#
# AND SO IS TIGHTNESS, which is a separate claim and the one a future edit is
# likelier to break. Every rectangle is the bounding box of probe images, and a
# probe image is a neighbour of a block rim cell on another face — a `Vertex()`
# halo cell by definition — so the candidate stream is the `Vertex()` halo cell
# for cell, with no surplus for `NativeCheck` to reject. Widening any bound by
# one cell would still ANSWER correctly, because the check filters what the
# rectangles over-propose, so no oracle comparison can see it; the test file
# counts the candidate stream and requires the equality, which is what makes a
# lazy bounding box here a failure rather than a silent slowdown.
#
# WHAT FALLS BACK. One configuration only: a system with more faces than
# `_BAND_RECT_CAP`, which none of the three is. Everything else — every flush
# side, every face corner, every whole-face block, both connectivities — is
# walked. The guard is kept because it is the one assumption the derivation
# cannot check itself, and a `small_push` past capacity would be a `BoundsError`
# from inside an iterator rather than an honest fallback.
# ---------------------------------------------------------------------------

# The distinct probe positions of a block, at most four — a corner of a block
# flush on two sides is an endpoint of both. Eight is the number of `(side,
# endpoint)` pairs, so the list never has to reject a push.
const _PROBE_CAP = 8
const ProbeList = Helpers.SmallList{_PROBE_CAP,NTuple{2,Int64}}

@inline _empty_probe_list() =
    Helpers.empty_small_list(Val(_PROBE_CAP), (Int64(0), Int64(0)))

@inline function _add_probe(probes::ProbeList, sx::Int64, sy::Int64)
    for i in 1:length(probes)
        p = @inbounds probes[i]
        p[1] == sx && p[2] == sy && return probes
    end
    return Helpers.small_push(probes, (sx, sy))
end

# One probe: everything the native one-ring of the rim cell at `(sx, sy)` can
# see, bucketed by face. In-face neighbours already inside the home band box are
# dropped rather than merged, so the home rectangle stays the tight band;
# anything else on the home face — a face adjacent to ITSELF across a seam,
# which none of the three has — would widen it, and correctly.
function _seam_probe(sys, grid, rects::BandRects, target::Int, face::Int64,
        sx::Int64, sy::Int64, home::FaceRect)
    for nb in neighbors(grid, lattice_cell(sys, target, sx, sy, face), 1;
            connectivity = Vertex())
        ix, iy, g = lattice_decode(sys, nb)
        gx = Int64(ix)
        gy = Int64(iy)
        g == face && home.x0 <= gx <= home.x1 && home.y0 <= gy <= home.y1 &&
            continue
        rects = _merge_rect(rects, FaceRect(g, face_orientation(sys, g),
            gx, gy, gx, gy))
    end
    return rects
end

function _seam_band_engine(sys::AbstractHierarchicalGridSystem, curve,
        c::AbstractCellIndex, target::Int, connectivity::Connectivity,
        x0::Int64, y0::Int64, side::Int64, face::Int64, n::Int64,
        home::FaceRect)
    # The faces ARE the root cells, and the root level is `first(levels(sys))` —
    # not 0, which is only the systems bundled here.
    nfaces = ncells(levelgrid(sys, first(levels(sys))))
    nfaces <= _BAND_RECT_CAP ||
        return generic_halo_engine(sys, c, target, connectivity)
    grid = levelgrid(sys, target)
    rects = Helpers.small_push(_empty_band_rects(), home)
    # The two extreme rim cells of every flush side — eight positions naming at
    # most four distinct cells, because a flush CORNER is an endpoint of both of
    # its sides and a whole-face block names each of its four corners twice.
    # `_merge_rect` already makes a repeat idempotent, so the deduplication is
    # for cost, not correctness: `_seam_probe` is a `neighbors` call, which
    # allocates a `Vector` on S2, and a whole-face block would otherwise pay for
    # eight of them to learn what four say. Still O(1) in the halo's size either
    # way — construction is a constant-time act however deep the target.
    probes = _empty_probe_list()
    if x0 == 0
        probes = _add_probe(probes, x0, y0)
        probes = _add_probe(probes, x0, y0 + side - 1)
    end
    if x0 + side == n
        probes = _add_probe(probes, x0 + side - 1, y0)
        probes = _add_probe(probes, x0 + side - 1, y0 + side - 1)
    end
    if y0 == 0
        probes = _add_probe(probes, x0, y0)
        probes = _add_probe(probes, x0 + side - 1, y0)
    end
    if y0 + side == n
        probes = _add_probe(probes, x0, y0 + side - 1)
        probes = _add_probe(probes, x0 + side - 1, y0 + side - 1)
    end
    for i in 1:length(probes)
        p = @inbounds probes[i]
        rects = _seam_probe(sys, grid, rects, target, face, p[1], p[2], home)
    end
    # Face order is canonical order, so the list is sorted once, here, by
    # walking the faces rather than the rectangles. Faces are `0:nfaces-1` and
    # there are at most twelve of them, so this is a fixed 144-comparison pass
    # and not a sort in any sense that could grow.
    ordered = _empty_band_rects()
    for g in 0:(nfaces - 1)
        for i in 1:length(rects)
            q = @inbounds rects[i]
            if q.face == g
                ordered = Helpers.small_push(ordered, q)
                break
            end
        end
    end
    # A rectangle whose face is not in `0:nfaces-1` never gets picked up, and a
    # short `ordered` is a SHORT HALO — the one way this walk could answer wrong
    # rather than fall back. `_merge_rect` keeps at most one rectangle per face,
    # so the lengths agree exactly when every face was seen. If they do not, the
    # system's face numbering is not the one the id law assumes, and this file's
    # rule is to fall back, never to approximate — the same answer the
    # `_BAND_RECT_CAP` guard above gives to the other assumption it cannot check.
    length(ordered) == length(rects) ||
        return generic_halo_engine(sys, c, target, connectivity)
    return SquareBandEngine(curve,
        NativeCheck(sys, grid, c, level(c), connectivity),
        target, n, Int32(face), x0, y0, side, connectivity isa Vertex, ordered)
end

# ===========================================================================
# The calibrated directed walk, shared by the two aperture-7 systems
# ===========================================================================

# A subtree of H3 or IGeo7 is not a block of anything: its rim is a hexagonal
# spiral, its halo wraps a shape with no lattice box, and the aperture is odd, so
# nothing above applies. What both systems DO have is a subtree-rim automaton
# over an arc of exposed lattice directions — `(L, s)` meaning the arc
# `s, s+1, …, s+L-1 (mod 6)` — and the halo is reachable through it from the
# OTHER side.
#
# THE IDEA. Every level-`target` halo cell of `root` lies in the subtree of one
# of `root`'s same-level neighbours. So walk the neighbours, and inside each one
# walk only the part of its subtree that faces `root` — which is a rim walk of
# that neighbour, entered not at the fully exposed `(6, 0)` a subtree root gets
# but at the short arc that points back at `root`.
#
# WHY IT IS CONTAINED, which is the whole load-bearing claim. The NESTED
# ADJACENCY LEMMA: if `y` is adjacent to `x` at level `l`, then `parent(y)` is
# `parent(x)` or is adjacent to it. Induct: a halo cell `x` at level `l` touches
# some descendant of `root`, so `parent(x)` touches `root`'s subtree at level
# `l-1`; `parent(x)` is not a descendant of `root` (or `x` would be one), so by
# induction it lies under a neighbour of `root`, and so does `x`. Verified
# exhaustively over 910,560 adjacency pairs on both systems and both
# connectivities, with zero violations.
#
# WHY THE ARC IS OBSERVED AND NEVER TABULATED. The two automata have EXCHANGED
# parity branches — H3's even-level branch is IGeo7's odd-level one — and their
# two `L < 6` guards are tested in the opposite order, so any table of seeds
# fitted on one system is wrong on the other, and a table fitted per parity is
# wrong at the other parity. The admission guard `o < L` is the one part that is
# parity-independent in both files, so the set of children a seeded arc admits at
# depth one is identical at both parities in both systems. `_hex_calibrate`
# therefore derives the arc from the children that are OBSERVED to touch
# `root` — one native check per child of one neighbour — and the asymmetry
# becomes invisible. That is the single reason one driver serves both systems.
#
# WHY CONCATENATION IS A MERGE, again. Distinct same-level cells have disjoint
# subtrees, so ordering the neighbours by `first(descendant_range(sys, nb,
# target))` puts their candidate blocks in ascending, non-overlapping order;
# within a neighbour the automaton is digit-lexicographic, which is ascending id.
# No heap, no seen-set, no sort — 25,536 pairwise range comparisons across both
# systems found zero overlaps.
#
# WHERE "OUTSIDE" IS STATED, because unlike every other engine in this file these
# two never test it. `OutsideWalkEngine` retires the subject subtree by integer
# range containment (`_admit`'s `_HALO_SKIP`), `ScanHaloEngine` compares
# `ancestor` against the root, and `SquareBandEngine` prunes the block by lattice
# box (`_inside_block`) — each an explicit "this candidate is not a descendant".
# The hexagonal engines have no such line, and they need one: `_touches_root`
# would ACCEPT a descendant of the root, since a descendant's neighbours are
# descendants too. What carries it is the invariant that
# `neighbors(grid, c, 1; connectivity)` never returns `c` itself. Every candidate
# here is a descendant of a cell in `neighbors(levelgrid(sys, lc), c, 1)`, so if
# that ring excluded nothing the root would be walked as its own neighbour and
# its whole rim would be emitted as halo. No bundled system's one-ring lists the
# cell it was asked about; a system whose did would need an explicit skip here,
# not merely a wider guard.
#
# THE WALK IS EXACT, NOT CONSERVATIVE. Candidate-to-halo ratio is 1.0000 at every
# depth from two down; run with the check disabled over 4,622 cases it produced
# zero surplus candidates. The check is kept anyway, on every candidate, because
# it is this file's exactness contract and it costs what the halo already costs.
#
# DEPTH ONE HAS NO AUTOMATON. At `target == rootlevel + 1` the calibration IS the
# answer — the touching children are the halo — so `HexChildHaloEngine` emits
# each neighbour's children filtered by the same check, and no arc is derived at
# all. Seeding an automaton to walk one level is where the win thins to nothing,
# which is the sign that the automaton would be doing no work there.
#
# WHAT IT BUYS, measured against the generic outside-first walk on a full halo:
# 7-36x on H3 and 1.1-16x on IGeo7 over depths one to five, widening with depth
# because the generic walk's cost grows with the target LEVEL while this one's
# grows with the halo. The prefix is the sharper number: taking ten cells of an
# IGeo7 depth-seven halo was 42 ms and 779 KB through the generic walk and is
# 1.1 ms and 256 bytes here — and the 256 bytes do not move with the depth, which
# is the design's laziness law rather than a speed-up.

# Six is a hexagon's neighbour count and five a pentagon's, so the list is never
# more than six long on either system. Eight is that plus slack, so a system with
# a wider ring reaches the guard in `hex_halo_engine` rather than a `BoundsError`
# from inside `small_push`.
const _HEX_RING_CAP = 8

"""
    HexNeighbour(cell, lo, arclen, start)

One of `root`'s same-level neighbours, with the arc `(arclen, start)` calibrated
to face `root` and `lo = first(descendant_range(sys, cell, target))`, the key the
ring is kept sorted by. `arclen == 0` marks an uncalibrated entry, which only the
depth-one engine holds.
"""
struct HexNeighbour{C}
    cell::C
    lo::Int
    arclen::Int8
    start::Int8
end

const HexRing{C} = Helpers.SmallList{_HEX_RING_CAP,HexNeighbour{C}}

@inline _empty_hex_ring(c::C) where {C} =
    Helpers.empty_small_list(Val(_HEX_RING_CAP), HexNeighbour(c, 0, Int8(0), Int8(0)))

# Insertion sort over at most six entries: the whole list fits in registers and a
# heap would be a heap-allocation. Shifts down from the end, so the list is
# sorted by `lo` at every point and the walk can read it straight through.
function _hex_insert(ring::HexRing{C}, e::HexNeighbour{C}) where {C}
    ring = Helpers.small_push(ring, e)
    i = length(ring)
    while i > 1
        prev = @inbounds ring[i - 1]
        prev.lo <= e.lo && break
        ring = Helpers.small_setindex(ring, prev, i)
        i -= 1
    end
    return Helpers.small_setindex(ring, e, i)
end

# ---------------------------------------------------------------------------
# Calibration
# ---------------------------------------------------------------------------

"""
    _minimal_arc(p, q) -> (arclen, start)

The shortest arc of the six-direction ring containing both `p` and `q`, or
`(0, 0)` when there is not exactly one such arc or its length is not 2 or 3.

The scan is the definition rather than the closed form (`arclen` is the ring
distance plus one, taken the short way round) so that the UNIQUENESS question is
answered by counting, not argued: two directions exactly opposite each other
admit two four-arcs and no shorter one, and this returns `(0, 0)` there instead
of picking one. Thirty-six iterations of integer work, once per neighbour.
"""
function _minimal_arc(p::Int, q::Int)
    bestlen = 7
    beststart = 0
    ties = 0
    for L in 1:6, s in 0:5
        (mod(p - s, 6) < L && mod(q - s, 6) < L) || continue
        if L < bestlen
            bestlen = L
            beststart = s
            ties = 1
        elseif L == bestlen
            ties += 1
        end
    end
    ties == 1 || return (0, 0)
    (bestlen == 2 || bestlen == 3) || return (0, 0)
    return (bestlen, beststart)
end

"""
    _hex_calibrate(sys, grid1, root, rootlevel, nb, connectivity) -> (arclen, start)

The arc of `nb` that faces `root`, derived by asking which of `nb`'s children are
halo cells of `root`'s subtree one level down, or `(0, 0)` if any of the guards
fails.

Measured over 52,182 `(root, neighbour)` pairs on both systems, root levels 0-11:
the number of touching children is ALWAYS exactly two, the minimal covering arc
is always unique, and its length is 2 in 98.6% of cases and 3 in the remaining
1.4% — never 1, never 4 or more. The arc-3 case is fully characterised: `nb` is a
pentagon and its deleted direction lies strictly between the two touching
children. Nothing else predicts it, which is the other reason this is measured
per call rather than looked up.

So none of the returns below was observed to fire. That is what a guard on a
structural claim should look like — it is here because the claim is evidence and
not yet a theorem, and a system whose one-ring changed would meet it rather than
meeting a wrong answer.
"""
function _hex_calibrate(sys, grid1, root, rootlevel::Int, nb,
        connectivity::Connectivity)
    p = -1
    q = -1
    n = 0
    for k in children(sys, nb)
        _touches_root(sys, grid1, root, rootlevel, k, connectivity) || continue
        n += 1
        n > 2 && return (0, 0)
        d = hex_child_direction(sys, k)
        # The centre child is enclosed by its six siblings and cannot reach out
        # of its parent at all, so it has no direction to put on the ring. If one
        # ever tested as touching, the ring model is wrong for this system and
        # the arc must not be guessed.
        d < 0 && return (0, 0)
        n == 1 ? (p = d) : (q = d)
    end
    n == 2 || return (0, 0)
    return _minimal_arc(p, q)
end

"""
    _hex_validate(sys, root, rootlevel, ring, connectivity) -> Bool

Does every depth-two halo cell under each neighbour really lie on that
neighbour's calibrated walk? The calibration only observes depth one, so this is
the one direct check that one level of observation carries to the next.

`O(1)` in the halo's size — at most seven grandchildren per child per
neighbour — but that is still about 300 native checks, which EXCEEDS the entire
directed walk at the shallow depths. `hex_halo_engine` therefore runs it only
from depth three up; at depth two the generic engine is barely slower than the
validation alone, so skipping the specialization there would cost more than
trusting it does.

NOT ONE NUMBER, which is why none is quoted. The cost rises with the neighbour
ring — a six-neighbour root pays more than a five-neighbour one — and with the
root's level. Quoting a point value for either system hides both, and a reader
who measures one root will not reproduce it. To measure it, take construction at
depth three minus construction at depth two: this is the only threshold it
crosses, so that difference is the whole of it.

The threshold rests on the RATIO to the walk, not on the absolute cost, and that
ratio is the same on both systems: the validation exceeds the depth-three walk,
falls below the depth-four walk, and is a small fraction of the depth-five walk.
So it is not worth paying at depth two, where the generic engine is barely slower
than the validation alone, and any change to the threshold is a change for both
systems or for neither. Construction never grows with the halo, so the laziness
law holds whatever the constant is.

An earlier version of this paragraph asserted the opposite — a rounding error on
H3 against a term that dominated on IGeo7 — on the grounds that IGeo7's
`neighbors` was lattice arithmetic where H3's was a library call. That absolute
gap was real and is now gone, closed by the ported GBT kernel, and even while it
stood the RATIOS were never as far apart as the costs were. It was the ratios the
argument needed. Do not reintroduce a per-system lever here without measuring
them first.

Deliberately NOT restricted to pentagons or to arc-3 neighbours, which would be
question-begging: the premise under test is exactly that the ordinary cases need
no test.

Both sides are ascending — the automaton is digit-lexicographic, and nested
`children` calls are ascending by contract — so containment is a two-pointer
merge with no set and no allocation.
"""
function _hex_validate(sys, root, rootlevel::Int, ring::HexRing{C},
        connectivity::Connectivity) where {C}
    target = rootlevel + 2
    grid = levelgrid(sys, target)
    for i in 1:length(ring)
        e = @inbounds ring[i]
        arc = seeded_rim_engine(sys, e.cell, target, Int(e.arclen), Int(e.start))
        r = iterate(arc)
        for k in children(sys, e.cell), x in children(sys, k)
            _touches_root(sys, grid, root, rootlevel, x, connectivity) || continue
            while r !== nothing && r[1] < x
                r = iterate(arc, r[2])
            end
            (r !== nothing && r[1] == x) || return false
            r = iterate(arc, r[2])
        end
    end
    return true
end

# ---------------------------------------------------------------------------
# Depth one: the calibration is the answer
# ---------------------------------------------------------------------------

"""
    HexChildHaloEngine(system, grid, root, rootlevel, target, connectivity, ring)

`target == rootlevel + 1`: each neighbour's children in turn, native-checked, the
neighbours in descendant-range order. No automaton and no calibration — at this
depth the children ARE the candidates.

`children` is re-read per step rather than stored, which keeps the walk state two
integers and costs a bounded-container rebuild of at most seven ids.
"""
struct HexChildHaloEngine{S,G,C,K}
    system::S
    grid::G
    root::C
    rootlevel::Int
    target::Int
    connectivity::K
    ring::HexRing{C}
end

Base.eltype(::Type{<:HexChildHaloEngine{S,G,C,K}}) where {S,G,C,K} = C
Base.IteratorSize(::Type{<:HexChildHaloEngine}) = Base.SizeUnknown()

"""
    HexChildWalk(slot, next)

Which neighbour of the ring is being read, and which of its children comes next.
"""
struct HexChildWalk
    slot::Int
    next::Int
end

Base.iterate(e::HexChildHaloEngine) = iterate(e, HexChildWalk(1, 1))

function Base.iterate(e::HexChildHaloEngine, w::HexChildWalk)
    slot = w.slot
    next = w.next
    while slot <= length(e.ring)
        kids = children(e.system, (@inbounds e.ring[slot]).cell)
        if next > length(kids)
            slot += 1
            next = 1
            continue
        end
        x = @inbounds kids[next]
        next += 1
        _touches_subtree(IndexedNeighbors(), e, x) || continue
        return (x, HexChildWalk(slot, next))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Deeper: one seeded automaton per neighbour
# ---------------------------------------------------------------------------

"""
    HexArcHaloEngine(system, grid, root, rootlevel, target, connectivity, ring)

`target > rootlevel + 1`: the system's own rim automaton seeded with each
neighbour's calibrated arc, walked to `target`, every leaf native-checked before
it is yielded. The ring is in descendant-range order and the blocks are disjoint,
so concatenating the neighbours' streams is already the canonical merge.

Memory is `O(depth)`: one seeded engine and its frame stack, both isbits, plus
the fixed ring. Nothing sized by the halo is ever built, so a prefix costs what
the prefix costs.

[`Base.IteratorSize`](@ref) is `SizeUnknown()` and there is NO `length`, even
though the counts are known: a CALIBRATED neighbour's stream is `(3^d + 1)/2`
cells for both arc lengths — a calibrated arc-2 seed emits that many outright,
and a calibrated arc-3 seed emits `3^d` but arc-3 happens only at a pentagon,
whose deleted digit removes precisely the 4-arc branch and collapses the census
back — so the halo is `3^(d+1) + 3` around a hexagon and `5(3^d + 1)/2` around a
pentagon. CALIBRATED is load-bearing in that sentence and not a hedge: the census
describes the arcs [`_hex_calibrate`](@ref) produces, over 13,692 of which the
seeded walk emits exactly `(3^d + 1)/2` with pentagon neighbours included. It is
NOT a statement about an arbitrary `(L, s)` — on an H3 level-3 pentagon an arc-2
seed at `s = 0` emits 1, 3, 9 and 27 leaves for `d = 1…4`, not 2, 5, 14 and 41,
and does so identically on all twelve of that level's pentagons —
so nobody should reuse the formula to size a seeded walk of their own. That was
verified in
176/176 configurations (all twelve pentagons of each system, hexagons adjacent to
them, hexagons far from them, both connectivities, root levels 0-3, depths 1-5)
and end to end to depth six. It is still ENUMERATION, not a derivation from the
transition recurrence, and the design admits a count as an API contract only once
it is derived symbolically and validated around every pentagon and parity
configuration. What remains to prove: that the seeded transition relation has the
claimed leaf census at every parity, and that the pentagon deletion always
removes the 4-arc branch rather than another. Until then the `MethodError` from
`length` is the honest answer.
"""
struct HexArcHaloEngine{S,G,C,K}
    system::S
    grid::G
    root::C
    rootlevel::Int
    target::Int
    connectivity::K
    ring::HexRing{C}
end

Base.eltype(::Type{<:HexArcHaloEngine{S,G,C,K}}) where {S,G,C,K} = C
Base.IteratorSize(::Type{<:HexArcHaloEngine}) = Base.SizeUnknown()

# The census the docstring above refuses to declare as a `length`, offered as
# the thing it IS good enough to be: a `sizehint!`. `3^(d+1) + 3` is the count
# around a hexagon and an over-estimate of about 20% around a pentagon, whose
# census is `5(3^d + 1)/2` — and the hexagon figure is the larger of the two, so
# one formula bounds both. Both engines answer it, at depth one and deeper,
# because the depth-one halo is the same census with `d == 1`.
#
# WHY THIS IS NOT THE `length` THE DOCSTRING WITHHOLDS. Nothing has changed
# about the evidence: the census is still enumeration rather than a derivation
# from the seeded transition relation. What changed is the promise. A `length`
# that is 20% high hands a caller `undef` slots as cell ids through `collect`'s
# sized route; a hint that is 20% high costs a `Vector` that is 20% roomier than
# it needed to be, and is never read as a count by anything.
_halo_sizehint(e::HexChildHaloEngine) = _hex_sizehint(e.target - e.rootlevel)
_halo_sizehint(e::HexArcHaloEngine) = _hex_sizehint(e.target - e.rootlevel)

@inline _hex_sizehint(d::Int) = 3^(d + 1) + 3

"""
    HexArcWalk(slot, arc, state)

The walk state: which neighbour of the ring is being descended, that neighbour's
seeded automaton, and the automaton's own frame stack. All three are isbits, and
the engine is rebuilt rather than stored per slot, so resuming needs nothing that
was not returned.
"""
struct HexArcWalk{A,W}
    slot::Int
    arc::A
    state::W
end

@inline function _hex_arc_engine(e::HexArcHaloEngine, slot::Int)
    nb = @inbounds e.ring[slot]
    return seeded_rim_engine(e.system, nb.cell, e.target, Int(nb.arclen),
        Int(nb.start))
end

# Open neighbour `slot` and every later one until a candidate passes the check.
# Written as a fresh-start scan rather than folded into the resume path so the
# walk state is only ever built from a narrowed, concrete iterate result — a
# `Union{Nothing,W}` state field would infect the engine's whole return type.
function _hex_arc_open(e::HexArcHaloEngine, slot::Int)
    while slot <= length(e.ring)
        arc = _hex_arc_engine(e, slot)
        r = iterate(arc)
        while r !== nothing
            x, w = r
            _touches_subtree(IndexedNeighbors(), e, x) &&
                return (x, HexArcWalk(slot, arc, w))
            r = iterate(arc, w)
        end
        slot += 1
    end
    return nothing
end

Base.iterate(e::HexArcHaloEngine) = _hex_arc_open(e, 1)

function Base.iterate(e::HexArcHaloEngine, s::HexArcWalk)
    r = iterate(s.arc, s.state)
    while r !== nothing
        x, w = r
        _touches_subtree(IndexedNeighbors(), e, x) &&
            return (x, HexArcWalk(s.slot, s.arc, w))
        r = iterate(s.arc, w)
    end
    return _hex_arc_open(e, s.slot + 1)
end

# ---------------------------------------------------------------------------
# Building one
# ---------------------------------------------------------------------------

"""
    hex_halo_engine(sys, c, target, connectivity)

The halo engine for a system with a seeded rim automaton — H3 and IGeo7, whose
`halo_engine` methods are this call and nothing else. Falls back to
[`generic_halo_engine`](@ref) whenever a guard fires, so the specialization is
never the reason an answer is wrong, only the reason it is fast.

The guards, in order: the system must have sorted subtrees (the ring's
descendant-range order is what makes concatenation a merge); the target must be
strictly deeper than `c` and no deeper than the system (both level errors belong
to the generic engine, and depth zero is its one-ring); the ring must fit; every
neighbour must calibrate; and from depth three up the calibration must survive
[`_hex_validate`](@ref). None of the last three was observed to fire anywhere in
the spike that measured this design.

THE RING PROBE is ALWAYS `Vertex()`, whatever was requested, so one containment
argument covers both connectivities: the `Edge()` halo is a subset of the
`Vertex()` one, and a superset of the larger is a superset of the smaller. Only
the ring is probed that way. [`_hex_calibrate`](@ref) and [`_hex_validate`](@ref)
below are called with the REQUESTED connectivity, so on a system whose `Edge()`
adjacency is strictly smaller than its `Vertex()` one an `Edge()` query can find
fewer than two touching children under a neighbour, and `_hex_calibrate` answers
`(0, 0)` — a whole-root fallback to [`generic_halo_engine`](@ref), which is
correct and slower. The superset ring buys that system nothing under `Edge()`;
what it buys is that the ring never has to be re-derived per connectivity.

That asymmetry is deliberate rather than an oversight, and calibrating at
`Vertex()` instead would be the riskier code. A `Vertex()`-calibrated arc IS a
conservative band for the `Edge()` halo, and the emit check would filter it —
but its coverage has never been measured, because all 52,182 calibrations behind
this design ran on systems where the two connectivities coincide. So the choice
is between a documented fallback that is merely slow and a fast path whose
containment nobody has tested, and this file's rule is that an unproved band does
not reach the caller. Moot on both shipped systems; written down because the
first system where it is not moot should meet a fallback, not a surprise.
"""
function hex_halo_engine(sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex, target::Int, connectivity::Connectivity)
    lc = level(c)
    (has_sorted_subtrees(sys) && lc < target <= max_level(sys)) ||
        return generic_halo_engine(sys, c, target, connectivity)
    nbs = neighbors(levelgrid(sys, lc), c, 1; connectivity = Vertex())
    length(nbs) <= _HEX_RING_CAP ||
        return generic_halo_engine(sys, c, target, connectivity)
    grid = levelgrid(sys, target)
    ring = _empty_hex_ring(c)
    if target == lc + 1
        for nb in nbs
            ring = _hex_insert(ring, HexNeighbour(nb,
                first(descendant_range(sys, nb, target)), Int8(0), Int8(0)))
        end
        return HexChildHaloEngine(sys, grid, c, lc, target, connectivity, ring)
    end
    grid1 = levelgrid(sys, lc + 1)
    for nb in nbs
        arclen, start = _hex_calibrate(sys, grid1, c, lc, nb, connectivity)
        arclen == 0 && return generic_halo_engine(sys, c, target, connectivity)
        ring = _hex_insert(ring, HexNeighbour(nb,
            first(descendant_range(sys, nb, target)), Int8(arclen), Int8(start)))
    end
    if target - lc >= 3 && !_hex_validate(sys, c, lc, ring, connectivity)
        return generic_halo_engine(sys, c, target, connectivity)
    end
    return HexArcHaloEngine(sys, grid, c, lc, target, connectivity, ring)
end
