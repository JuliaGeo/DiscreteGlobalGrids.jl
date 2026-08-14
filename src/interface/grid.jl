# ---------------------------------------------------------------------------
# Base-grid interface contracts. `Int` arguments are dense positions;
# `AbstractCellIndex` arguments are cell identities. Geometry uses
# `GO.UnitSphericalPoint` unless a converting wrapper states otherwise.
#
# A hierarchical system implements the four required primitives as system-level
# methods instead; see `src/interface/system.jl`. `levelgrid` returns a
# `HierarchicalLevelGrid`, which forwards to them.
# ---------------------------------------------------------------------------

# ===========================================================================
# Required primitives
# ===========================================================================

"""
    ncells(grid::AbstractGrid) -> Int

The number of cells in `grid`, which is the length of its canonical dense
order: positions run over `1:ncells(grid)`.

**Required** of every [`AbstractGrid`](@ref).

This is `ConservativeRegridding.Trees.ncells` — the same binding, extended here
— so any grid is a `Trees` source without an import or a wrapper.

Must be O(1) and must not change over the lifetime of the grid object.
"""
function ncells end

"""
    cellindex(grid::AbstractGrid, i::Int) -> AbstractCellIndex

The canonical typed id of the cell at **position** `i` in `grid`'s dense order.

**Required** of every [`AbstractGrid`](@ref).

The returned id is of type `cellindextype(system(grid))` for a grid that has a
system, and of the grid's own canonical id type otherwise. Together with
[`cellposition`](@ref) this is a bijection `1:ncells(grid)` ↔ the grid's cells:

    cellposition(grid, cellindex(grid, i)) == i   for all i in 1:ncells(grid)

`i` outside `1:ncells(grid)` throws a `BoundsError`.

See also [`cellindex(grid, i, T)`](@ref cellindex), which requests a specific
index scheme.
"""
function cellindex end

"""
    cell_boundary(grid::AbstractGrid, c::AbstractCellIndex) -> AbstractVector{<:GO.UnitSphericalPoint}

The exact boundary ring of cell `c`, as points on the unit sphere.

**Required** of every [`AbstractGrid`](@ref).

Contract:

  - The ring is **implicitly closed**: the first vertex is *not* repeated at the
    end. [`cell_polygon`](@ref) is what closes it.
  - Vertices are in **counter-clockwise order seen from outside the sphere**
    (right-hand rule about the outward normal), so the ring bounds the cell
    rather than its complement and spherical signed area comes out positive.
  - Consecutive vertices are joined by great-circle arcs. A boundary that is
    curved in the system's own chart must be densified here, because every
    consumer treats the result as a spherical polygon.
  - Every point is unit-norm to within a few `eps`.

The container may be any `AbstractVector` — a `Vector`, a static vector, or a
lazily computed one. Callers must not mutate it.
"""
function cell_boundary end

"""
    cell_centroid(grid::AbstractGrid, c::AbstractCellIndex) -> GO.UnitSphericalPoint

A representative point of cell `c` on the unit sphere.

**Required** of every [`AbstractGrid`](@ref).

The point must lie **strictly inside** the cell, never on its boundary — it is
what `cellat(grid, cell_centroid(grid, c)) == c` is tested against, and what
labelling, plotting and nearest-cell code uses as *the* location of the cell.
Systems that can compute a true area centroid should; systems that cannot may
return any interior point, and should say which in their own documentation.
"""
function cell_centroid end

# ===========================================================================
# Identity: levels, raw encodings, index schemes
# ===========================================================================

"""
    level(c::AbstractCellIndex) -> Int
    level(grid::AbstractGrid) -> Union{Int,Nothing}

The refinement level of a cell id or grid. `level(c)` is total on
[`AbstractCellIndex`](@ref): each id encodes its level without a system or grid.

`level(grid)` is the level of the system grid `grid` is drawn from, or
`nothing` for a standalone grid with no hierarchy. See [`system`](@ref).

Levels are non-negative, increase with refinement, and are compared with
`<`; the valid range for a system is [`levels(sys)`](@ref levels).
"""
function level end

"""
    rawid(c::AbstractCellIndex) -> Integer

The encoded integer behind a typed cell id: the `UInt64` of an `H3Cell` or
`Z7Cell`, the linear index of a [`LevelIndex`](@ref).

This is the value to write to disk, print in hex, or hand to a C library — not
a position, and not something to do arithmetic on unless the system documents
what the arithmetic means. Round-tripping it back into a typed id needs the id
type (and, for [`LevelIndex`](@ref)-style schemes, the level), so `rawid` is
lossy on its own; prefer passing typed ids.
"""
function rawid end

"""
    cellindex(grid::AbstractGrid, i::Int, T::Type{<:AbstractCellIndex}) -> T

The id of the cell at position `i`, in the requested index scheme `T` rather
than the system's canonical one — `cellindex(grid, i, Z7Cell)`, say, from a
grid whose canonical scheme is something else.

Provided generically as `reindex(T, system(grid), cellindex(grid, i))`;
overridable when a system can produce `T` directly. Throws if `T` is not in
[`cellindextypes`](@ref).
"""
cellindex(::AbstractGrid, ::Int, ::Type{<:AbstractCellIndex})

"""
    reindex(T::Type{<:AbstractCellIndex}, sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex) -> T

Convert cell id `c` to index scheme `T` within the same system: same cell, same
level, different encoding. `reindex(typeof(c), sys, c) === c`.

Conversion is exact and total between any two schemes the system supports —
they are alternative names for the same tessellation, not approximations of it.
Requesting an unsupported `T` throws an `ArgumentError` naming
[`cellindextypes(sys)`](@ref cellindextypes).
"""
function reindex end

"""
    cellindextypes(sys::AbstractHierarchicalGridSystem) -> Tuple{Vararg{Type{<:AbstractCellIndex}}}
    cellindextypes(grid::AbstractGrid)

The index schemes `sys` can name its cells in, canonical one first (so
`first(cellindextypes(sys)) === cellindextype(sys)`). Every listed type is
reachable through [`reindex`](@ref) in both directions.

Defaults to `(cellindextype(sys),)` — one canonical scheme, no conversions.

The grid form asks the grid's system; a standalone grid
(`system(grid) === nothing`) defaults to `(typeof(cellindex(grid, 1)),)` — the
one scheme it demonstrably names cells in (empty grids have no schemes to
report).
"""
function cellindextypes end

"""
    cellposition(grid::AbstractGrid, c::AbstractCellIndex) -> Union{Int,Nothing}

The **position** of cell `c` in `grid`'s dense order, or `nothing` when absent.
It is the inverse of [`cellindex`](@ref).

`c` may be given in any scheme in [`cellindextypes`](@ref); it is
[`reindex`](@ref)ed to canonical first. A `c` at a different level than the
grid is not an error either — it is simply not in the grid, so the answer is
`nothing`.
"""
function cellposition end

# ===========================================================================
# Geometry
# ===========================================================================

"""
    cell_polygon(grid::AbstractGrid, c::AbstractCellIndex) -> GI.Polygon

Cell `c` as a GeoInterface polygon on the unit sphere: the
[`cell_boundary`](@ref) ring, explicitly closed, wrapped in a
`GI.LinearRing` inside a `GI.Polygon`.

Coordinates are unit-sphere `(x, y, z)`, not longitude/latitude — this polygon
is meant for spherical predicates, spherical area, and
`ConservativeRegridding`, all of which work in that frame.
"""
function cell_polygon end

"""
    cell_area(grid::AbstractGrid, c::AbstractCellIndex) -> Float64

The true spherical area of cell `c`, in steradians. Multiply by `R^2` for area
on a sphere of radius `R`; using [`authalic_sphere`](@ref) gives ellipsoidal
area for an equal-area DGGS.

The generic fallback computes the spherical area of [`cell_polygon`](@ref).
Systems must override when the published boundary only approximates the true
cell, as for a densified curved boundary, and document the returned quantity.
The result is never a planar area.
"""
function cell_area end

"""
    cell_extent(grid::AbstractGrid, c::AbstractCellIndex) -> Extents.Extent{(:X, :Y)}

The longitude/latitude bounding box of cell `c`, in **degrees**: the
GeoInterface extent of the cell, for interoperating with extent-based
machinery that does not know about spheres.

This is a lon/lat rectangle, and it is conservative where a rectangle cannot be
exact: a cell containing a pole, or crossing the antimeridian, gets an `X` span
of `(-180, 180)`. It is therefore **not** the quantity tree pruning uses —
that is [`node_extent`](@ref), a `SphericalCap`, which has no such degeneracies.
Reach for `cell_extent` at an interoperability boundary, not in an algorithm.
"""
function cell_extent end

"""
    getcell(grid::AbstractGrid, i::Int) -> GI.Polygon

The cell at **position** `i` as a unit-sphere polygon. This extends
`ConservativeRegridding.Trees.getcell` and is implemented as

    cell_polygon(grid, cellindex(grid, i))

`getcell` takes a position; [`cell_polygon`](@ref) takes an id.
"""
function getcell end

# ===========================================================================
# Location
# ===========================================================================

"""
    cellat(grid::AbstractGrid, p::GO.UnitSphericalPoint) -> Union{AbstractCellIndex,Nothing}
    cellat(grid::AbstractGrid, lon::Real, lat::Real)

The cell of `grid` containing point `p`, or `nothing` if the point is outside
the grid's coverage.

The unit-sphere method is the primitive; the `(lon, lat)` method is a
converting wrapper and takes **degrees**.

`nothing` is a real answer, not an error: grids need not cover the sphere.

**Ties.** A shared-boundary point is assigned deterministically per platform to
one incident cell. The generic fallback selects the first candidate in canonical
id order. Floating-point boundary tests do not guarantee cross-platform
bit-identical ties; each system documents its rule. A tie must never select a
nonincident cell or return `nothing` inside coverage.

The generic implementation descends [`treeify(grid)`](@ref treeify) to a
candidate set and then tests point-in-cell; systems with a closed-form inverse
projection override it and should say what their tie rule is.
"""
function cellat end

# ===========================================================================
# Topology
# ===========================================================================

"""
    neighbors(grid::AbstractGrid, c::AbstractCellIndex, k::Int = 1; connectivity::Connectivity = Vertex())
    neighbors(grid::AbstractGrid, p::Int, k::Int = 1; connectivity::Connectivity = Vertex()) -> Vector{Int}

All cells of `grid` within `k` adjacency steps of `c`, excluding `c`.

# Connectivity

[`Vertex()`](@ref Vertex), which includes vertex contact, is the default.
[`Edge()`](@ref Edge) requires a shared edge. They coincide where exactly three
cells meet at every vertex — the icosahedral hexagon-with-pentagon family
(IGeo7, H3) — and differ wherever a vertex carries more, including one
*pentagonal* system: A5's Cairo-style tiling gives 11 vertex against 3 edge
neighbours at resolution 1. Quadrilateral grids add corner neighbours under
`Vertex()`: four cells to a vertex in the lattice interior, five where ISEA4R's
diamonds meet an icosahedral vertex.

# Order

Order is part of the contract. Results concatenate rings outward:

    neighbors(grid, c, k) == vcat(ring(grid, c, 1), ring(grid, c, 2), ..., ring(grid, c, k))

Within each ring, cells run counter-clockwise in the tangent plane at
`cell_centroid(grid, c)`, viewed from outside the sphere, from a system-defined
start. Results are not id-sorted. The geometric fallback uses the first ring-1
neighbour as zero azimuth and breaks exact azimuth ties by canonical id.

Cells with fewer neighbours yield shorter rings without padding.

# Container

Any ordered, indexable collection with the grid's cell-index `eltype`.

# Coverage, and what a subset means

Neighbours outside the grid's coverage are omitted, not padded.

On a **subset of a complete level** — [`PartialGrid`](@ref),
[`CellVector`](@ref), `CellLookup` — adjacency and distance are the complete
level's, clipped to membership:

    ring(sub, c, k) == filter(in(sub), ring(levelgrid(system(sub), level(sub)), c, k))

Distance is therefore measured in the *system*, never inside the subset. A hole
in the subset removes cells; it does not lengthen the path around itself, and a
cell reachable only by leaving the subset and coming back keeps the distance the
complete level gives it. The two readings coincide at `k == 1` and part company
from `k == 2`, which is why this is stated rather than left to the reader.

`c` outside the subset is an `ArgumentError`, not a complete-grid answer.

# Positions

Given a position, both verbs answer with **in-set positions in ascending order**
— the form a stencil table is indexed by. Ids keep the rotational order above;
positions do not, because an index list is read by membership rather than by
direction. [`halo_table`](@ref) is this form for a whole grid at once.

`k` must be ≥ 0; `k == 0` returns an empty collection. See [`ring`](@ref) for
the cells at *exactly* distance `k`.
"""
function neighbors end

"""
    ring(grid::AbstractGrid, c::AbstractCellIndex, k::Int; connectivity::Connectivity = Vertex())
    ring(grid::AbstractGrid, p::Int, k::Int; connectivity::Connectivity = Vertex()) -> Vector{Int}

The cells at adjacency distance exactly `k` from `c`. `ring(grid, c, 0)` is
`c` alone. The ordered result satisfies

    neighbors(grid, c, k) == vcat(ring(grid, c, 1), ..., ring(grid, c, k))

so ring `k` is the final ordered block of `neighbors(grid, c, k)`. Overrides
must preserve this equality.

`ring` carries the same order (counter-clockwise seen from outside the sphere,
from the system's documented start), container, coverage and subset-clipping
contracts as [`neighbors`](@ref), including the position form's ascending order.
"""
function ring end

"""
    halo_table(grid::AbstractGrid, k::Int = 1; connectivity::Connectivity = Vertex()) -> Vector{Vector{Int}}
    halo_table(cv::CellVector, k::Int = 1; connectivity::Connectivity = Vertex())
    halo_table(lk::CellLookup, k::Int = 1; connectivity::Connectivity = Vertex())

The whole stencil at once: entry `p` is `neighbors(grid, p, k)`, the in-set
positions within `k` adjacency steps of position `p`, ascending.

    halo_table(sub, k)[p] == neighbors(sub, p, k)

is the law, so this is never a second answer — only a faster route to the same
one. On a subset the clipping is [`neighbors`](@ref)'s: system adjacency
intersected with membership, omitted rather than padded, so rows have varying
length and a cell whose neighbours all lie outside gets an empty one.

A stencil pass is then one comprehension over the table.
"""
function halo_table end

# ===========================================================================
# Trees
# ===========================================================================

"""
    treeify(grid::AbstractGrid)
    treeify(manifold::GeometryOpsCore.Manifold, grid::AbstractGrid)

A spatial tree over `grid`, used for traversal pruning and
`ConservativeRegridding`. This extends `ConservativeRegridding.Trees.treeify`
and is total on [`AbstractGrid`](@ref).

The result implements `GeometryOps.SpatialTreeInterface`. **Node extents are
`GO.UnitSpherical.SphericalCap`s** at every level of every tree here, which is
the whole predicate vocabulary tree descent needs.

A hierarchical grid uses its [`node_extent`](@ref) hierarchy; other grids use a
fallback tree over position space.

The one-argument form picks the manifold with `best_manifold(grid)`.
"""
function treeify end

# ===========================================================================
# Queries
# ===========================================================================

"""
    query(grid::AbstractGrid, pred::DE9IM.DE9IMPredicate) -> Vector{<:AbstractCellIndex}
    query(sys::AbstractHierarchicalGridSystem, pred::DE9IM.DE9IMPredicate; level::Integer) -> Vector{<:AbstractCellIndex}

Every cell satisfying the spatial predicate `pred`, as a **sorted** `Vector` of
typed cell ids.

# Predicates

Predicates are re-exported DE9IM.jl wrappers such as `Intersects(target)`,
`Covers(target)`, and `Touches(target)`. `Base.parent(pred)` returns the target.
This package defines their spherical semantics.

The target may be a GeoInterface geometry, an `Extents.Extent`, or a
`GO.UnitSpherical.SphericalCap`. Longitude/latitude targets are lifted to the
unit sphere once, at the boundary of the call.

# Semantics

Tree pruning may over-select candidates; prepared spherical predicates determine
the exact result.

The `sys` method answers at the requested `level` without materialising the
level grid.
"""
function query end

# ===========================================================================
# Provenance
# ===========================================================================

"""
    system(grid::AbstractGrid) -> Union{AbstractHierarchicalGridSystem,Nothing}

The hierarchical system containing `grid`, or `nothing` for a standalone grid.
The default is `nothing`.

See also [`level`](@ref).
"""
function system end

# The two provenance defaults. A grid built by a system overrides both; a
# standalone grid inherits "no hierarchy" and is complete at the base
# interface.
system(::AbstractGrid) = nothing
level(::AbstractGrid) = nothing
