# The two-level hierarchy: 64 800 tiles at level 0, one tile's raster at level 1.
# Every relation here is arithmetic on the dense 0-based ordinals of `bands.jl`.

# `levelgrid(CopernicusDEMSystem(...), l)` is the package's `HierarchicalLevelGrid`:
# all 64 800 tiles, or every pixel on Earth, in ordinal order. This system's fast
# paths hang off the alias, and the five primitives it forwards to are the
# `(sys, ...)` methods below.
const LevelGrid{N} = DGG.HierarchicalLevelGrid{CopernicusDEMSystem{N}}

# ===========================================================================
# System interface
# ===========================================================================

DGG.cellindextype(::CopernicusDEMSystem) = DGG.LevelIndex
DGG.levels(::CopernicusDEMSystem) = 0:1
DGG.has_sorted_subtrees(::CopernicusDEMSystem) = true

"""
    max_neighbors(CopernicusDEMSystem(...), connectivity) -> Int

`8` under `Vertex()` and `4` under `Edge()`: the Moore and von Neumann bounds of a
raster lattice.

!!! warning "This system does not implement `neighbors`"
    The bound is the *interior lattice* bound, and it is stated because
    [`max_neighbors`](@ref) has no default. It is **not** a claim about the whole
    sphere: a pixel in the top row of a pole tile meets every other pixel of the pole
    ring at the pole itself, and a pixel just below a band boundary meets up to ten
    coarser pixels above it. Anyone adding [`neighbors`](@ref) to this system must
    revisit this number first.
"""
DGG.max_neighbors(::CopernicusDEMSystem, ::DGG.Vertex) = 8
DGG.max_neighbors(::CopernicusDEMSystem, ::DGG.Edge) = 4

"Lazy 0-based ids at one level; `rootcells` and `children` are windows over it."
struct IdRange <: AbstractVector{DGG.LevelIndex}
    level::Int32
    first::Int64      # 0-based index of element 1
    n::Int
end

Base.size(v::IdRange) = (v.n,)
Base.IndexStyle(::Type{IdRange}) = Base.IndexLinear()
Base.@propagate_inbounds function Base.getindex(v::IdRange, i::Int)
    @boundscheck checkbounds(v, i)
    return DGG.LevelIndex(v.level, v.first + i - 1)
end

"""
    rootcells(CopernicusDEMSystem(...))

All 64 800 level-0 tiles, `LevelIndex(0, 0:64799)`, ascending, as a **lazy** vector.

Lazy on purpose: `PartialGrid` reads `first(rootcells(sys))` on every construction
(`src/fallbacks/partial_grid.jl`), and a materialised 64 800-element vector would
allocate about a megabyte per chunk built. The count is far above the "small, cheap
collection" the contract has in mind (12 for HEALPix), which costs the generic tree
descent one cap evaluation per tile at the synthetic root; that is the price of a grid
whose base tessellation is the 1° graticule.
"""
DGG.rootcells(::CopernicusDEMSystem) = IdRange(Int32(0), Int64(0), NTILES)

# ===========================================================================
# The level grid: size, and positions <-> ids
# ===========================================================================

DGG.ncells(sys::CopernicusDEMSystem, l::Integer) =
    Int(l) == 0 ? NTILES :
    Int(l) == 1 ? Int(tables(sys).rowbase[end]) :
    throw(ArgumentError("level $l is outside $(DGG.levels(sys))"))

# The grid bounds-checks `i`, so this is the bijection and nothing else.
DGG.cellindex(::CopernicusDEMSystem, l::Integer, i::Int) = DGG.LevelIndex(l, i - 1)

"""
    cellposition(CopernicusDEMSystem(...), c) -> Union{Int,Nothing}

Closed form: `index + 1` for an in-range id, and `nothing` for one no cell has.
Never throws — a miss is an answer — so this is the one decoder that does not go
through the id guard.
"""
function DGG.cellposition(sys::CopernicusDEMSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    (l == 0 || l == 1) || return nothing
    0 <= c.index < DGG.ncells(sys, l) || return nothing
    return Int(c.index + 1)
end

# ===========================================================================
# The hierarchy
# ===========================================================================

"""
    parent(CopernicusDEMSystem(...), pixel) -> LevelIndex

The tile the pixel belongs to. Throws an `ArgumentError` on a level-0 cell, which is
a tile and is a root.
"""
function Base.parent(sys::CopernicusDEMSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    l == 1 || throw(ArgumentError(l == 0 ?
        "level-0 Copernicus DEM cell $c is a tile and has no parent" :
        "level $l is outside $(DGG.levels(sys))"))
    r, q, _, _ = decode(sys, c)
    return DGG.LevelIndex(0, tileordinal(r, q))
end

"""
    children(CopernicusDEMSystem(...), tile)

The tile's pixels in raster order — north row first, west to east — as a **lazy**
vector of `ncols * N` ids. Lazy because a GLO-30 tile in the 1x band has 12 960 000
of them and the tree cursor iterates without collecting.
"""
function DGG.children(sys::CopernicusDEMSystem{N}, c::DGG.LevelIndex) where {N}
    l = DGG.level(c)
    l == 0 || throw(ArgumentError(l == 1 ?
        "level-1 Copernicus DEM cell $c is a pixel, at max_level 1, and has no children" :
        "level $l is outside $(DGG.levels(sys))"))
    r, q, _, _ = decode(sys, c)
    return IdRange(Int32(1), tilebase(sys, r, q), Int(ncols(sys, r)) * N)
end

"""
    ancestor(CopernicusDEMSystem(...), c, l) -> LevelIndex

`c` itself at `l == level(c)`, and [`parent`](@ref) at `l == 0` for a pixel. Two
levels leave no third case.
"""
function DGG.ancestor(sys::CopernicusDEMSystem, c::DGG.LevelIndex, l::Integer)
    target = Int(l)
    lc = DGG.level(c)
    target <= lc || throw(ArgumentError(
        "ancestor level $target is deeper than the cell's own level $lc"))
    target >= 0 || throw(ArgumentError(
        "ancestor level $target is above the root level 0"))
    target == lc && return c
    return Base.parent(sys, c)
end

"""
    descendant_range(CopernicusDEMSystem(...), tile, 1) -> UnitRange{Int}

The tile's contiguous window of level-1 positions: `tilebase + 1 : tilebase + ncols*N`.

Exact and hole-free in both directions, which is what `has_sorted_subtrees == true`
asserts. It holds because the level-1 order is tile-major and raster-order within a
tile, so a tile's pixels are consecutive by construction and no id in the window
belongs to another tile.

`l == level(c)` is the cell's own one-element position range; `l < level(c)` throws an
`ArgumentError`.
"""
function DGG.descendant_range(sys::CopernicusDEMSystem{N}, c::DGG.LevelIndex,
        l::Integer) where {N}
    target = Int(l)
    lc = DGG.level(c)
    index = _checked_index(sys, c)          # also rejects a level outside 0:1
    target >= lc || throw(ArgumentError(
        "descendant level $target is above the cell's own level $lc"))
    target <= DGG.max_level(sys) || throw(ArgumentError(
        "descendant level $target is past max_level $(DGG.max_level(sys))"))
    if target == lc
        pos = Int(index + 1)
        return pos:pos
    end
    r, q, _, _ = decode(sys, c)
    base = tilebase(sys, r, q)
    return Int(base + 1):Int(base + ncols(sys, r) * Int64(N))
end

"""
    descendants(CopernicusDEMSystem(...), c, l)

Every level-`l` descendant of `c`, ascending, as the same **lazy** vector
[`children`](@ref) returns — the ids are [`descendant_range`](@ref) read off as
consecutive ordinals, and a GLO-30 tile's 12 960 000 of them are not worth
materialising.

!!! warning "This diverges from the interface"
    The interface docstring for [`descendants`](@ref) says the call *materializes*
    `O(subtree)` ids, and on every other system it hands back a freshly allocated
    `Vector` the caller owns. **This method does not.** It returns a **lazy, read-only
    `AbstractVector`** that computes each id on indexing. Reading is complete — `length`,
    `getindex`, iteration, `collect` — but nothing that writes works, because there is no
    array to write into: no `setindex!`, no `push!`, no `sort!`, and no passing it to an
    API that mutates its argument. **`collect` it first if you need any of those.** The
    divergence is deliberate: one GLO-30 tile has 12 960 000 level-1 descendants, and a
    `Vector{LevelIndex}` of them is 16 bytes apiece.
"""
function DGG.descendants(sys::CopernicusDEMSystem, c::DGG.LevelIndex, l::Integer)
    r = DGG.descendant_range(sys, c, l)     # validates `l` both ways
    return IdRange(Int32(l), Int64(first(r) - 1), length(r))
end
