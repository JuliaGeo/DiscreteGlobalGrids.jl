# ---------------------------------------------------------------------------
# The base grid interface.
#
# Four required primitives, and every other generic declared with the contract
# it is implemented against. Implementations of the non-required generics live
# in `src/fallbacks/`; this file is the contract, and it is what a third-party
# implementor reads.
#
# Two conventions hold throughout, without exception:
#
#   * A bare `Int` argument is a POSITION in `1:ncells(grid)`. A typed
#     `AbstractCellIndex` argument is an IDENTITY. See `AbstractGrid`.
#   * All geometry is on the unit sphere, as `GO.UnitSphericalPoint`. Longitude
#     and latitude appear only in explicitly named converting wrappers.
#
# A hierarchical system does not implement the four on a grid type of its own:
# `levelgrid` hands back a `HierarchicalLevelGrid`, and the four arrive as the
# system-level methods documented in `src/interface/system.jl`.
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

The refinement level of a cell id, or of a grid.

`level(c)` is **total on [`AbstractCellIndex`](@ref)**: every canonical id is
self-describing about its level, either through in-band bits (`Z7Cell`,
`H3Cell`) or an explicit field ([`LevelIndex`](@ref)). No system, grid or table
is needed to answer it. This is precisely what lets the interface drop the old
`(system, level, id)` argument triple: an id already carries its level, so
passing one alongside it can only create a way to disagree.

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

The **position** of cell `c` in `grid`'s dense order, or `nothing` if `c` is
not in `grid`.

The inverse of [`cellindex`](@ref). `nothing` — rather than an error — is the
answer for a cell outside the grid, because asking whether a cell is present is
the normal way to intersect an id set with a partial grid, and because
neighbour and query results at a coverage edge are full of legitimate misses.

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

The area of cell `c` **in steradians**. Multiply by `R^2` for a physical area on
a sphere of radius `R`; for an equal-area DGGS read on its authalic sphere (see
[`authalic_sphere`](@ref)) that product is the true ellipsoidal area.

# What "the cell" means here

`cell_area` reports the area of the **true cell** — the region the system says
this cell *is*. For most systems the true cell is exactly the polygon
[`cell_boundary`](@ref) publishes, and then this is that polygon's spherical
area and the two are the same statement. But the two can come apart, in both
directions, and this method follows the cell rather than the ring:

  - A system may publish a boundary that is a **densification** of a curved true
    cell. HEALPix is this case: a pixel is an analytic equal-area diamond of
    area exactly `4π/(12·4^level)`, and the published ring is a polyline
    approximation of it that converges from below. `cell_area` returns the
    closed form; the ring is the approximation, not the definition.
  - A system may publish a boundary that tiles the sphere exactly while its
    *chart* cells do not have that boundary. IGeo7 is this case: the corner
    rings partition the sphere, so their areas are the honest answer to "how
    much sphere does this cell own", even though the equal-area chart assigns a
    slightly different figure (`IGeo7.equal_area_steradians`, which differs by
    +1.6% on hexagons and −9.9% on pentagons at level 1). There the ring *is*
    the true cell, and `cell_area` is ring-derived.

So the rule is: **the true cell, not the ring** — and for the systems where the
true cell is the ring, which is the common case, they coincide. A system whose
`cell_area` is not simply the area of its own `cell_boundary` must say so in its
`cell_area` docstring and say which quantity it returns.

The generic fallback computes the spherical area of the [`cell_polygon`](@ref)
ring, which is the right answer whenever the ring is the cell. A system
overrides only to *correct* that, never to speed it up at the cost of changing
the answer.

Never planar. Spherical throughout, so it is right for cells that span a large
solid angle, where a planar formula is not merely inaccurate but wrong in sign.
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

The cell at **position** `i` as a unit-sphere polygon.

This is `ConservativeRegridding.Trees.getcell`, extended here with its `Trees`
meaning intact: `Trees` addresses cells by dense position, and this is how a
regridder walks a grid. Implemented once, generically, as

    cell_polygon(grid, cellindex(grid, i))

Grid authors never write it. Note the argument convention — `getcell` takes a
position, `cell_polygon` takes an id; they are the two halves of the same
lookup, kept separate so neither can be confused for the other.
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

**Ties.** A point exactly on a shared boundary belongs to exactly one cell, and
which one is deterministic and documented per system. The generic fallback
resolves a tie by taking the first candidate in canonical id order.

"Deterministic" here means **per platform**: the same call, in the same process
or in another process on the same machine and Julia build, always gives the same
answer, and a system must document the rule it uses to get there. It does *not*
promise bit-identity across platforms. It cannot: a tie is decided by comparing
floating-point quantities that a different CPU, libm or `--math-mode` may round
differently in the last place, and the set of points where that changes the
answer is a measure-zero curve. Requiring cross-platform bit-identity would
force every system to carry exact arithmetic on its boundary test for a class of
input that no real workload lands on by accident — the wrong trade. What *is*
promised everywhere is the part that matters: the answer is always one of the
cells genuinely incident to the point, never a third cell and never `nothing`
inside coverage. Systems whose tie rule is an exact integer or lattice
comparison get cross-platform identity as a bonus and may say so.

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

All cells of `grid` within `k` adjacency steps of `c`, **excluding `c` itself**.

# Connectivity

[`Vertex()`](@ref Vertex) — Moore, share at least a vertex — is the **default**,
because it is the superset a caller cannot rebuild from an edge-only answer.
[`Edge()`](@ref Edge) is the opt-in restriction. The two coincide exactly where
**three cells meet at every vertex** — the icosahedral hexagons-with-twelve-
pentagons family (IGeo7, H3) and the dodecahedral resolution-0 shell — because
there a shared vertex always comes with a shared edge. Where a vertex carries
more than three cells the two differ, and that includes one *pentagonal*
system: A5's Cairo-style tiling has 11 vertex neighbours against 3 edge
neighbours at resolution 1. On quadrilateral grids `Vertex()` adds the corners
— four cells to a vertex in the lattice interior, five where ISEA4R's diamonds
meet an icosahedral vertex.

# Order

**The order is part of the contract, and it is rotational.** The result is the
rings 1, 2, …, `k` concatenated outward, shell by shell:

    neighbors(grid, c, k) == vcat(ring(grid, c, 1), ring(grid, c, 2), ..., ring(grid, c, k))

Within one ring the cells run **counter-clockwise seen from outside the
sphere** — counter-clockwise in the tangent plane at `cell_centroid(grid, c)`,
viewed from above that plane — beginning at a cell each system documents. The
rings are never interleaved and the result is never sorted by id.

Because ring `k` is the last block, [`ring`](@ref) is the *tail* of `neighbors`:

    ring(grid, c, k) == neighbors(grid, c, k)[end - length(ring(grid, c, k)) + 1 : end]

That is the law, and it is what the rotational order is for. Position `j` of a
ring always names the same direction relative to `c`, so a stencil weight vector
means the same thing at every cell and can be rotated as a unit; a coarser
neighbourhood extends a finer one by appending, without recomputing it. Sorting
by id destroys both properties, which is why no method here sorts.

Each system picks its own ring-1 start (a lattice direction, a compass point);
what a system must not do is leave it unstated. Where a system has a natural
direction only for the immediate ring, the recommended extension — and what the
geometric fallback does — is to measure the outer rings' azimuths **relative to
the first ring-1 neighbour**, so every ring starts on the same spoke.

The geometric fallback has no lattice to read a direction off, so it realises
the same contract by measurement: it orders each ring by azimuth about
`cell_centroid(grid, c)`, counter-clockwise seen from outside, with the
first ring-1 neighbour as the zero direction and exact ties broken by ascending
canonical id. No method in this package is ever allowed to return neighbours in
an unspecified order.

A cell with fewer neighbours than the lattice maximum — a pentagon, or a cell on
the edge of a partial grid — simply yields a shorter ring. The winding is
unchanged; there is no padding and no gap marker.

# Container

The grid's choice: any ordered, indexable collection whose `eltype` is the
grid's cell index type. A `SmallCollections.SmallVector` sized by
[`max_neighbors`](@ref) is recommended at small `k`, so a neighbour sweep does
not allocate.

# Coverage

Neighbours that fall outside a partial grid's coverage are **absent** from the
result — never zero-padded, never `nothing`-filled. Padding is a lookup-layer
convenience; it is not interface semantics, and a caller that wants it can get
it from [`cellposition`](@ref) returning `nothing`.

`k` must be ≥ 0; `k == 0` returns an empty collection. See [`ring`](@ref) for
the cells at *exactly* distance `k`.
"""
function neighbors end

"""
    ring(grid::AbstractGrid, c::AbstractCellIndex, k::Int; connectivity::Connectivity = Vertex())

The cells of `grid` at adjacency distance **exactly** `k` from `c` — the shell,
not the disc. `ring(grid, c, 0)` is `c` alone.

Derived from [`neighbors`](@ref) and overridable. The two are related exactly,
not merely as sets: the disc is the rings concatenated outward,

    neighbors(grid, c, k) == vcat(ring(grid, c, 1), ..., ring(grid, c, k))

so `ring(grid, c, k)` is the **tail block** of `neighbors(grid, c, k)`,
element for element and in order. A system that overrides one of these must
override the other consistently — computing the two with independent walks that
happen to agree as sets is the way this law is usually broken.

`ring` carries the same order (counter-clockwise seen from outside the sphere,
from the system's documented start), container and coverage contracts as
[`neighbors`](@ref).
"""
function ring end

# ===========================================================================
# Trees
# ===========================================================================

"""
    treeify(grid::AbstractGrid)
    treeify(manifold::GeometryOpsCore.Manifold, grid::AbstractGrid)

A spatial tree over `grid`, for pruning traversals and for
`ConservativeRegridding`.

This is `ConservativeRegridding.Trees.treeify`, extended here, and it is
**total on [`AbstractGrid`](@ref)** — every grid can be treeified, so the full
`Trees` surface ([`ncells`](@ref), [`getcell`](@ref),
`GeometryOpsCore.best_manifold`, `Trees.should_parallelize`) works for every
grid this package can describe.

The result implements `GeometryOps.SpatialTreeInterface`. **Node extents are
`GO.UnitSpherical.SphericalCap`s** at every level of every tree here, which is
the whole predicate vocabulary tree descent needs.

A grid from a hierarchical system treeifies to a `HierarchicalGridCursor` whose
node extents are the system's own [`node_extent`](@ref) — the hierarchy *is*
the tree, so there is nothing to build. Any other grid gets a fallback tree
built over position space from cell extents. Hierarchy is purely a fast path.

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

The predicate vocabulary is [DE9IM.jl](https://github.com/rafaqz/DE9IM.jl)'s
functor wrappers, re-exported here: `Intersects(target)`, `Covers(target)`,
`Touches(target)`, `Within(target)`, and the rest. `Base.parent(pred)` unwraps
the target; keywords ride along with it. DE9IM.jl supplies only the types —
every semantic is implemented in this package, on the sphere.

The target may be a GeoInterface geometry, an `Extents.Extent`, or a
`GO.UnitSpherical.SphericalCap`. Longitude/latitude targets are lifted to the
unit sphere once, at the boundary of the call.

# Semantics

**Exact.** The tree prunes with [`node_extent`](@ref) under the covering law,
which can only ever over-select; the surviving candidates are then decided by
prepared spherical predicates. Pruning is an optimisation and never appears in
the answer.

The `sys` method answers at the requested `level` without materialising the
level grid.
"""
function query end

# ===========================================================================
# Provenance
# ===========================================================================

"""
    system(grid::AbstractGrid) -> Union{AbstractHierarchicalGridSystem,Nothing}

The hierarchical system `grid` is a level (or a subset of a level) of, or
`nothing` for a standalone grid that has no hierarchy.

`nothing` is not a defect: a tripolar or gaussian grid is a perfectly good
[`AbstractGrid`](@ref) and simply stops at the base interface. Generic code
that wants a fast path tests for `nothing` and falls back; it must never
*require* a system.

The default is `nothing`, so a standalone grid implements nothing extra.

See also [`level`](@ref).
"""
function system end

# The two provenance defaults. A grid built by a system overrides both; a
# standalone grid inherits "no hierarchy" and is complete at the base
# interface.
system(::AbstractGrid) = nothing
level(::AbstractGrid) = nothing
