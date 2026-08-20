# ---------------------------------------------------------------------------
# Base-grid interface contracts. `Int` arguments are dense positions;
# `AbstractCellIndex` arguments are cell identities. Geometry uses
# `GO.UnitSphericalPoint` unless a converting wrapper states otherwise: `GO` is
# `GeometryOps`, and this package re-exports the point type, so an implementor
# writes `UnitSphericalPoint` with no module path.
#
# A hierarchical system implements the five level-grid primitives as
# system-level methods instead; see `src/interface/system.jl`. `levelgrid`
# returns a `HierarchicalLevelGrid`, which forwards to them.
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
`GI.LinearRing` inside a `GI.Polygon`. Neither wrapper allocates, so a system
whose `cell_boundary` uses inline storage gets an `isbits` polygon; read the
polygon through GeoInterface rather than depending on its container types.

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
    neighbors(grid::AbstractGrid, c::AbstractCellIndex, k::Integer = 1; connectivity::Connectivity = Vertex())
    neighbors(grid::AbstractGrid, p::Int, k::Integer = 1; connectivity::Connectivity = Vertex()) -> Vector{Int}

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

Order is part of the contract, and it is ONE order: every verb in this package
that hands back a neighbourhood hands it back counter-clockwise. Results
concatenate rings outward:

    neighbors(grid, c, k) == vcat(ring(grid, c, 1), ring(grid, c, 2), ..., ring(grid, c, k))

**Counter-clockwise, exactly.** Take `n = cell_centroid(grid, c)` as the outward
normal and project each ring member's centroid into the tangent plane at `n`.
Read in order, the azimuths of those projections in a right-handed frame
`(e₁, e₂ = n × e₁)` increase and wrap through `2π` **exactly once**. That is the
same rotational sense [`cell_boundary`](@ref)`(grid, c)` winds in: the package
has one handedness, fixed by the boundary contract, and every ring agrees with
it. A ring read backwards wraps `length - 1` times; an id-sorted one wraps some
arbitrary number of times. Neither is this order.

**Start.** The direction is guaranteed everywhere; the phase is guaranteed only
*within* a system. Ring 1 begins at the system's own documented direction — `+s`
for S2, `SW` for HEALPix, `NW` for CopernicusDEM, the development frame's `+1`
for IGeo7, the smallest-id neighbour for A5 and for the geometric fallback — and
rings `2:k` of the same call begin on the **same spoke** as ring 1, the azimuth
of `ring(grid, c, 1)[1]`. So a disc reads as concentric rings all starting in
one direction, but which cell that is differs by system and is not a portable
fact. Exact azimuth ties break by canonical id. The start is deterministic — a
property of the system and the cell alone, the same first member every time the
cell is asked, in every idiom, independent of any region or table the ring is
read through.

Cells with fewer neighbours yield shorter rings without padding, and an omitted
neighbour leaves no gap in the sequence.

# The idioms that carry this order

[`ring`](@ref), the position forms below, [`adjacency`](@ref),
[`member_neighbors`](@ref), the one-argument [`neighbors`](@ref) iterator, and
the rings `mapneighbors` and `foreachneighbors` pass to their callbacks — all of
them. Nothing in this package answers a neighbourhood question in ascending id
or position.

The verb that is ascending is [`halo`](@ref), and it is not a ring: it is a
fetch list, ordered so that a read is sequential. It says so where it is
documented.

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

The filter also pins the rotation: a clipped ring is the complete ring, read
from its canonical start, with non-members dropped in place. Its length is the
in-set degree, and which absolute slot a surviving member occupied is
deliberately not recoverable from the clipped ring — a consumer that needs slot
identity reads the complete level's ring.

`c` outside the subset is an `ArgumentError`, not a complete-grid answer.

`k` must be ≥ 0; `k == 0` returns an empty collection. See [`ring`](@ref) for
the cells at *exactly* distance `k`.

# Positions

Given a position, both verbs answer with **in-set positions in the rotational
order above**: `neighbors(grid, p, k)` is `neighbors(grid, cellindex(grid, p),
k)` mapped through [`cellposition`](@ref), element for element, with non-members
dropped. [`adjacency`](@ref) is this form for a whole region at once.

The result is therefore not sorted. An index list read only by membership does
not care, and one read by direction — a gradient, an upwind stencil — cannot be
written at all against a sorted list, which is why the position forms carry the
same order the id forms do.
"""
function neighbors end

"""
    ring(grid::AbstractGrid, c::AbstractCellIndex, k::Integer; connectivity::Connectivity = Vertex())
    ring(grid::AbstractGrid, p::Int, k::Integer; connectivity::Connectivity = Vertex()) -> Vector{Int}

The cells at adjacency distance exactly `k` from `c`. `ring(grid, c, 0)` is
`c` alone. The ordered result satisfies

    neighbors(grid, c, k) == vcat(ring(grid, c, 1), ..., ring(grid, c, k))

so ring `k` is the final ordered block of `neighbors(grid, c, k)`. Overrides
must preserve this equality.

`ring` carries [`neighbors`](@ref)' order, container, coverage and
subset-clipping contracts unchanged — including the position form's, which is
the same counter-clockwise order read through [`cellposition`](@ref).
"""
function ring end

"""
    one_ring(grid, c, connectivity) -> ordered neighbours of `c`

The `k == 1` primitive: `c`'s immediate neighbours in one counter-clockwise turn
from the system's own start direction, as [`neighbors`](@ref) documents. An
internal hook, not part of the public API.

A system implements this one method and inherits `neighbors`, `ring`, and the
shell walk behind both from `Fallbacks.adjacency_shells`. The generic method is
the geometric one — `Fallbacks.adjacent_cells` wound about the cell's centroid —
and a system with native adjacency overrides it.

The return may be any ordered, indexable collection with the grid's cell-index
`eltype`; a fixed-capacity `SmallVector` sized by [`maxneighbors`](@ref) is
what keeps the one-ring sweeps allocation-free.
"""
function one_ring end

"""
    halo(region; connectivity = Vertex(), cells = false)

The cells immediately OUTSIDE `region` that touch one of its members, lazily,
each exactly once, in ascending complete-level position.

A **region** is a subset of one complete level — [`PartialGrid`](@ref),
[`CellVector`](@ref), `CellLookup` — or a complete [`AbstractGrid`](@ref), which
has no outside and therefore an empty halo. `Vertex()` counts vertex contact,
`Edge()` requires a shared edge.

Yields **positions on the complete level grid**: a halo cell is by definition
absent from the region and has no position in it. `cells = true` yields cell ids
instead.

A cell punched out of the middle of a region is outside it and touches it, so it
joins the halo.

`collect` gives a `Vector` and `Set` gives a membership-queryable set — both are
Base's contracts over an iterator, and neither is overloaded to reach some other
product. [`sizehint`](@ref DiscreteGlobalGrids.sizehint) gives a cheap size
estimate where one exists.

The walk is serial; [`adjacency`](@ref) is the verb that threads.

A [`MultiOrderCellSet`](@ref) has no `halo`: its members may sit at different
levels, so there is no single level to answer at. Use
[`member_neighbors`](@ref).

See also [`border`](@ref) and [`interior`](@ref), the same boundary from inside.
"""
function halo end

"""
    border(region; connectivity = Vertex(), cells = false)
    interior(region; connectivity = Vertex(), cells = false)

The members of `region` that do (`border`) or do not (`interior`) have a
neighbour outside it, lazily, in ascending in-region position. The two are
disjoint and together are the region.

Yields **in-region positions** — `1:length(region)`, the index a data vector
laid out against the region is read by. `cells = true` yields cell ids instead.
([`halo`](@ref) yields complete-level positions instead, for the reason it
gives.) The region types are `halo`'s; a complete grid has no cell with an
absent neighbour, so its border is empty and its interior is all of it.

A region holding a whole rooted subtree walks its system's `O(border)`
automaton. Every other region scans its own cells and compares each clipped
one-ring against the complete one, `O(cells · degree)`.

Both walks are serial; [`adjacency`](@ref) is the verb that threads.
"""
function border end

@doc (@doc border)
function interior end

"""
    adjacency(region; halo = 0, connectivity = Vertex(), threaded = true) -> AdjacencyTable
    adjacency(region, hpos::AbstractVector{<:Integer}; connectivity = Vertex(), threaded = true)

The one-ring of every cell of `region` at once, cached: `adj[p]` is the ring of
in-region position `p` as a non-allocating view, in the counter-clockwise order
[`neighbors`](@ref) states. `length(adj)` is the region size, which is what an
entry is compared against to tell a region slot from a halo slot.

The region types are [`halo`](@ref)'s, the complete grid included —
`adjacency(levelgrid(sys, l))` is a whole level's adjacency.

# The three row shapes

  - `halo = 0` — rings CLIPPED to the region. Members outside it are dropped,
    the survivors keep their order, and the row length is the in-region degree.
  - `halo = 1` — rings COMPLETE, addressing a `[region; halo]` buffer: `1:n`
    names a cell of the region and `n + j` the `j`-th cell of `halo(region)`
    under the same connectivity. The halo is walked once, here, and the table
    keeps it — see [`halopositions`](@ref).
  - `halo = :mark` — rings COMPLETE, with `0` where a member is outside the
    region. Slot geometry with no halo to walk, materialise or fetch, which is
    what a direction codec reads.

`halo` above 1 throws. A row exists only for an in-region position, so a wider
receptive field is a wider region: `adjacency(grow(region, n); halo = 1)`.

# The anchor

Every row is a window onto the canonical ring the COMPLETE level answers —
`neighbors(levelgrid(system, level), cell)` — whose start is a deterministic
property of the system and the cell alone. The
complete-width shapes preserve SLOT INDICES — slot `k` of a row is ring member
`k`, in every table of every region containing that cell — so a direction code
is a property of the cell, persistable and stable across tables. The clipped
shape preserves ORDER but not slots: dropping a member shifts the ones after it,
and recovering slot identity from a clipped row is deliberately impossible.
Reach for `:mark` when the slot is the answer.

The second form takes a halo the caller already walked, as strictly ascending
complete-level positions, and builds the `halo = 1` table against it. A
neighbour in neither half is an `ArgumentError` naming the cell, not a short
row. The list need not be minimal, only ascending and covering.

`threaded` builds contiguous chunks in separate tasks and accepts `Bool` or
GeometryOps' `True()`/`False()`; the result is identical either way.
"""
function adjacency end

"""
    neighborcount(grid::AbstractGrid, c::AbstractCellIndex; connectivity::Connectivity = Vertex()) -> Int

Return `length(neighbors(grid, c))`. Systems may compute structural degrees
without constructing the ring. Subsets count only neighbours within the subset,
and an out-of-set cell throws as [`neighbors`](@ref) does.

A subset cell is interior exactly when
`length(neighbors(sub, c)) == neighborcount(complete, c)`, so a border scan needs
one ring and one count rather than two rings.
"""
function neighborcount end

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

"""
    subcursor(grid::AbstractGrid, inds::AbstractUnitRange) -> tree or `nothing`

The [`treeify`](@ref) tree restricted to the grid positions `inds`, with leaf
indices still in `grid`'s own position space, or `nothing` when this grid cannot
express that restriction.

The default is `nothing`. A grid whose tree is a lattice rectangle rather than a
cell hierarchy — a raster-backed system's — implements this so that a chunk of
it keeps a prunable tree instead of the regridder's bounding-cap fallback. A
hierarchical grid needs no method: regridding descends its own hierarchy to the
chunk's ancestor.

The result must cover exactly the cells at `inds`, no more.
"""
function subcursor end

subcursor(::AbstractGrid, ::AbstractUnitRange{<:Integer}) = nothing

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
