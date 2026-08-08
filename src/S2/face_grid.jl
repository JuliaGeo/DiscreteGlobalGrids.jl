# ---------------------------------------------------------------------------
# S2 dense face grids — the S2 instance of `src/core/face_grid.jl`
#
# `chart.jl` gives the six closed-form charts `[0, 1]² → S²` and the index maps
# over the `nside × nside` lattice. The package's shared face-grid layer
# (`src/core/face_grid.jl`) turns any such chart family into something a
# `ConservativeRegridding.Regridder` can consume: a *dense* S2 grid — all
# `6 nside²` cells of one resolution — as a `SpatialTreeInterface` tree. All
# that lives here is what is S2 about it: the system singleton, the two
# orderings, seven contract methods, the O(1) block-cap opt-in with its
# soundness proof, and the aliases this module's users type.
#
# Provenance: this file and its HEALPix and ISEA4R siblings were written as
# near-verbatim copies of one another, section for section — a deliberate
# rule-of-three duplication, since two instances is not enough to fix the shape
# of a shared abstraction. The third (ISEA4R) fixed it, and the shared layer is
# what the three files collapsed into; keeping the names identical down to
# `FaceChartGrid` is what made that hoist mechanical.
# ---------------------------------------------------------------------------

import ..DiscreteGlobalGrids as DGG
# Imported, not `using`: the first three are extended below with S2 methods, and
# the import is also what makes `using DiscreteGlobalGrids.S2: data_index`
# resolve. `num_cells` is the package generic the shared layer defines the
# `S2FaceSpace` / `S2FaceGrid` methods of.
import ..DiscreteGlobalGrids: data_index, lattice_index, validate_ordering, num_cells

# ---------------------------------------------------------------------------
# The system
# ---------------------------------------------------------------------------

"""
    S2FaceSystem()

S2 as a [`DiscreteGlobalGrids.FaceGridSystem`](@ref): six cube-face charts
`[0, 1]² → S²` (`stf_to_point`) over an `nside × nside` lattice per face.

This is the *chart* system, deliberately distinct from the `S2DGGS` registry
entry: the two disagree about what a resolution is (`max_nside` here bounds
lattice arithmetic; `max_level` there bounds the `s2_cellid` hierarchy, which
is not ported), and only this one is a legal parameter of `FaceGrid` and
friends.
"""
struct S2FaceSystem <: DGG.FaceGridSystem end

# ---------------------------------------------------------------------------
# The resolution
#
# Defined here, ahead of the orderings and the contract methods, because it is
# what the resolution-bound half of the shared contract dispatches on — see the
# "Resolution travels as a `FaceGridSpace`" section of `src/core/face_grid.jl`.
# ---------------------------------------------------------------------------

"""
    S2FaceSpace(nside)

The resolution of a dense S2 grid: `nside` cells along each edge of each of the
six cube faces, hence `6 * nside^2` cells in total.

Any `nside >= 1` is admissible — this is the *chart* resolution, not an
`s2_cellid` refinement level. The power-of-two restriction belongs to the
Hilbert index alone, and is enforced by [`HilbertOrder`](@ref) when (and only
when) a grid is built with it. `S2FaceSpace(3)` and `S2FaceSpace(5)` are
perfectly good S2 face grids that simply have no S2 id space.

```julia
space = S2FaceSpace(4)          # 96 cells
space.nside                     # 4
```

`num_cells(space)` is `6 * nside^2` — the length of the data vector any
[`AbstractS2Ordering`](@ref) indexes. (Mirrors HEALPix's `num_pixels`; S2 says
"cell" where HEALPix says "pixel".)

This is also the value the resolution-bound contract methods below dispatch on
(`face_cell_corners`, `data_index`, `lattice_index`, `validate_ordering`): it
carries a *checked* `nside`, so no path into the layer can supply one the chart
maps cannot be evaluated at. Each shim unwraps `space.nside` and calls the
raw-Int codec in `chart.jl`, which stays shared with `S2Kernel.jl`.

See also [`S2FaceGrid`](@ref).
"""
const S2FaceSpace = DGG.FaceGridSpace{S2FaceSystem}

# ---------------------------------------------------------------------------
# Orderings
#
# The lattice `(ix, iy, face)` says *where* a cell is; an ordering says *which
# slot of the data vector* it occupies. Keeping the two apart is what lets one
# grid type serve row-major layouts, Hilbert-ordered layouts, and (later) a
# custom permutation, without the tree code knowing which.
# ---------------------------------------------------------------------------

"""
    abstract type AbstractS2Ordering

How cells of a dense S2 grid map to positions in a data vector — the S2 branch
of [`DiscreteGlobalGrids.AbstractFaceOrdering`](@ref), which carries the full
extension contract — `data_index(o, space, ix, iy, face)`,
`lattice_index(o, space, j)` and the optional `validate_ordering(o, space)`,
where `space` is an [`S2FaceSpace`](@ref), `ix`, `iy` are 0-based in
`0:space.nside-1` and `face` is 0-based in `0:5`.

An [`S2FaceGrid`](@ref) is a lattice (`6 nside²` cells addressed by
`(ix, iy, face)`) plus one of these; the ordering is the *only* thing that
decides which column of a `ConservativeRegridding.Regridder` a cell lands in.
[`RowMajorOrder`](@ref) and [`HilbertOrder`](@ref) are the two shipped
instances; a Morton ordering is deliberately not shipped — S2 has no Morton
convention — but is about fifteen lines against the contract.
"""
abstract type AbstractS2Ordering <: DGG.AbstractFaceOrdering end

"""
    RowMajorOrder()

Row-major data ordering: position `j` holds the cell whose 0-based row-major id
is `j - 1`, i.e. `ix` fastest, then `iy`, then `face`.

Defined for **any** `nside >= 1` — the row-major closed forms
([`xyf_to_rowmajor`](@ref) / [`rowmajor_to_xyf`](@ref)) carry no power-of-two
restriction, which is what makes this the default and the ordering to reach for
at non-power-of-two resolutions.
"""
struct RowMajorOrder <: AbstractS2Ordering end

# `+ 1` because row-major ids are 0-based while data positions are 1-based; see
# the index-convention block at the top of `chart.jl`.
data_index(::RowMajorOrder, space::S2FaceSpace, ix::Integer, iy::Integer, face::Integer) =
    Int(xyf_to_rowmajor(ix, iy, face, space.nside)) + 1

lattice_index(::RowMajorOrder, space::S2FaceSpace, j::Integer) =
    rowmajor_to_xyf(j - 1, space.nside)

"""
    HilbertOrder()

Hilbert data ordering: position `j` holds the cell whose 0-based Hilbert id is
`j - 1`, where the id is `face * nside² + hilbert_position` — the scaffold
ordinal shape the `S2DGGS` registry entry records
(`face * 4^level + hilbert_position`, see `src/core/systems/s2.jl`). This is to
S2 what `NestedOrder` is to HEALPix.

Requires `nside = 2^k`: the Hilbert position is built two bits per level, which
only exists on a `2^k × 2^k` face. [`S2FaceGrid`](@ref) rejects any other
`nside` at construction rather than at first query.

!!! note "Alignment with native `s2_cellid` is intended, not yet verified"
    The Hilbert tables transcribe s2geometry's conventions, so a
    `HilbertOrder` grid at `nside = 2^level` is intended *by construction* to be
    column-for-column the native S2 cell order at `level`. That is not yet
    oracle-verified — there are no s2geometry fixtures in this repository, and
    the native 64-bit `s2_cellid` encoding is not ported, which is why
    `S2DGGS`'s hierarchy group still throws `NotPortedError`. What *is* pinned
    here: bijectivity, Hilbert locality, prefix nesting across resolutions, the
    per-face block structure, and — since `S2Kernel.jl` — that data position `j`
    of this ordering and scaffold ordinal `j - 1` of
    `cell_polygon(S2DGGS(), level, j - 1)` are bitwise the same polygon.
"""
struct HilbertOrder <: AbstractS2Ordering end

data_index(::HilbertOrder, space::S2FaceSpace, ix::Integer, iy::Integer, face::Integer) =
    Int(xyf_to_hilbert(ix, iy, face, space.nside)) + 1

lattice_index(::HilbertOrder, space::S2FaceSpace, j::Integer) =
    hilbert_to_xyf(j - 1, space.nside)

# The one thing that must be caught eagerly: a Hilbert-ordered grid at
# nside = 3 would construct fine and then throw from deep inside a dual-tree
# traversal, which is a miserable place to learn about it.
function validate_ordering(::HilbertOrder, space::S2FaceSpace)
    ispow2(space.nside) || throw(ArgumentError(
        "HilbertOrder requires nside = 2^k, got nside=$(space.nside); \
         use RowMajorOrder() for arbitrary nside"))
    return nothing
end

# ---------------------------------------------------------------------------
# The `FaceGridSystem` contract
# ---------------------------------------------------------------------------

# 6 cube faces. (Numerically `root_count(S2DGGS())` too — a coincidence of this
# system, not a contract, so the literal stays independent.)
DGG.nfaces(::S2FaceSystem) = 6

DGG.face_chart(::S2FaceSystem, s, t, face::Integer) = stf_to_point(s, t, face)

DGG.face_cell_corners(space::S2FaceSpace, ix::Integer, iy::Integer, face::Integer) =
    cell_corners(ix, iy, face, space.nside)

# Upper bound: `6 * nside^2` is `3 * 2^61 ≈ 6.9e18 < typemax(Int64)` at
# `nside = 2^30` and overflows silently at `2^31`. `2^30` is also the `nside` of
# the deepest native S2 level (`max_level(S2DGGS()) == 30`), so the arithmetic
# bound and the system bound coincide.
DGG.max_nside(::S2FaceSystem) = 2^30

DGG.default_ordering(::S2FaceSystem) = RowMajorOrder()

DGG.ordering_family(::S2FaceSystem) = AbstractS2Ordering

DGG.facegrid_prefix(::S2FaceSystem) = "S2"

# The cap of a block's four corners contains the whole block, so we can skip the
# generic method's inclusion of all perimeter vertices. For S2 this is provable
# rather than measured — here is the argument, which is the per-system soundness
# obligation this override carries:
#
#  1. Block edges are geodesic arcs: `u = const` on a cube face is a central
#     plane section, and so is `v = const`, so the block boundary is four
#     great-circle arcs. Any two points of one face are within `arccos(-1/3)`
#     ≈ 109.47° of each other (the cube-corner half-angle), and the endpoints of
#     a single within-face edge lie within a quarter turn: each edge subtends at
#     most π/2.
#  2. Along a geodesic arc, distance to a fixed point has an interior maximum
#     only where the arc passes through the antipode-side farthest point of its
#     full great circle — a point at distance >= π/2 from the reference point.
#  3. No block-boundary point is that far from the cap center `c ∝ p₁+p₂+p₃+p₄`.
#     For `p` on edge `p₁p₂` at arc parameter τ of an edge of length `L <= π/2`,
#     `p·p₁ + p·p₂ = cos τ + cos(L - τ) >= 1 + cos L >= 1`, while `p·p₃` and
#     `p·p₄` are each `>= -1/3` by the same-face bound in 1. Hence
#     `(Σᵢ pᵢ)·p >= 1 - 2/3 = 1/3 > 0`: every boundary point is strictly within
#     π/2 of `c`, so by 2. no interior maximum occurs and the farthest boundary
#     point from `c` is one of the four corners.
#  4. Corners in ⇒ block in: the block is the gnomonic image of a planar
#     rectangle, hence spherically convex, and the cap is geodesically convex
#     because its radius is below π/2 (`c·pᵢ >= (2/3)/‖Σp‖ >= 1/6`, so the
#     corner radius is at most `arccos(1/6)` ≈ 80.4°, and the `1.0001` slack
#     `circle_from_four_corners` applies keeps it there).
#
# So `circle_from_four_corners`' slerp midpoints and 1.0001 slack are pure
# insurance for S2 — for HEALPix, whose edges bulge off the geodesic, they are
# load-bearing.
DGG.cap_policy(::S2FaceSystem) = DGG.FourCornerCap()

# ---------------------------------------------------------------------------
# The S2 names for the shared types
# ---------------------------------------------------------------------------

"""
    S2FaceGrid(space::S2FaceSpace, ordering::AbstractS2Ordering)
    S2FaceGrid(nside::Integer; ordering = RowMajorOrder())

A complete S2 face grid at one resolution — all `6 nside²` cells — together
with the data ordering its cells are numbered by. `treeify(grid)` turns it into
a spatial tree ([`S2FaceRoot`](@ref)) that
`ConservativeRegridding.Regridder` consumes directly.

```julia
grid = S2FaceGrid(4; ordering = HilbertOrder())
R = ConservativeRegridding.Regridder(treeify(grid), treeify(other))
```

# Alignment

**Column `j` of a `Regridder` built on `treeify(grid)` is data-vector position
`j` under `ordering` — there is no permutation between the matrix and the
field.** This is structural, not a convention to remember: the per-face cursors
emit `data_index(ordering, ...)` as their leaf indices and `Trees.getcell(root,
j)` returns the polygon of `lattice_index(ordering, space, j)`, so the matrix is
*assembled* in data order. Concretely, with `RowMajorOrder` column `j` is
row-major id `j - 1`, and with `HilbertOrder` column `j` is scaffold ordinal
`j - 1` (`face * 4^level + hilbert_position` at `nside = 2^level`).

# Validation

`nside >= 1` is checked by [`S2FaceSpace`](@ref); the ordering gets a say too,
through `validate_ordering` — which is why
`S2FaceGrid(3; ordering = HilbertOrder())` throws an `ArgumentError` here rather
than failing later inside a traversal.
"""
const S2FaceGrid = DGG.FaceGrid{S2FaceSystem}

"""
    FaceChartGrid(manifold, space::S2FaceSpace, face, ordering)

One S2 cube face as an `nside × nside` `Trees.AbstractCurvilinearGrid`, so that
`Trees.TopDownQuadtreeCursor` can index it with no S2-specific descent logic.
`face` is 0-based (`0:5`). Internal: build an [`S2FaceGrid`](@ref) and
`treeify` it.

Cartesian cell index `(i, j)` (1-based) is lattice cell `(i-1, j-1)`, and
`Trees.getcell` returns its CCW polygon (see [`cell_polygon`](@ref) — CCW is
required by the convex-clip kernel, which clips a CW ring to EMPTY).
`Trees.getvertex(g, i, j)` is the *lattice* point `((i-1)/nside,
(j-1)/nside)`, `1:(nside+1)` in each direction, which is what drives the
bounding caps.

The two index maps are overridden away from the generic column-major default:
they target the *global* data layout of `ordering`, not a face-local one, so the
leaf indices this grid's cursor reports are already data-vector positions.
"""
const FaceChartGrid = DGG.FaceChartGrid{S2FaceSystem}

"""
    S2FaceRoot(manifold, space::S2FaceSpace, ordering)
    S2FaceRoot(nside::Integer, ordering = RowMajorOrder())

Root of the spatial tree over a dense S2 face grid: the whole sphere, with the
six cube faces as children. Each child is a `Trees.TopDownQuadtreeCursor` over
that face's [`FaceChartGrid`](@ref); leaf indices throughout are data-vector
positions under `ordering`, so `Trees.getcell(root, j)` is the cell of data
position `j` and a `Regridder`'s column `j` is that same position.

Build one with `treeify(::S2FaceGrid)` rather than by hand. Constructed
directly it re-runs the two ordering checks `S2FaceGrid` runs, so
`S2FaceRoot(3, HilbertOrder())` throws here rather than from inside a
traversal.
"""
const S2FaceRoot = DGG.FaceGridRoot{S2FaceSystem}

# ---------------------------------------------------------------------------
# Cell geometry
# ---------------------------------------------------------------------------

"""
    cell_polygon(ix, iy, face, nside) -> GI.Polygon

Cell `(ix, iy)` of `face` as a closed 4-gon on the unit sphere, from
[`cell_corners`](@ref).

Because S2 cell edges are great-circle arcs (see [`stf_to_point`](@ref)), this
4-gon is the cell *exactly* — not an approximation of it, as the HEALPix
`pixel_polygon` sibling necessarily is. Since every corner is a shared lattice
point evaluated by the same chart function, neighbouring polygons meet at
bit-identical vertices and the tessellation is exact.

The ring is **counter-clockwise as seen from outside the sphere**, which
`cell_corners` guarantees and which is a hard contract rather than a convention:
the convex-clip kernel that computes spherical intersections clips a clockwise
ring to EMPTY, so a reversed ring yields silent zero areas instead of an error.

`nside` is taken as a plain integer here — this is the chart-side vocabulary
wrapper, alongside [`cell_corners`](@ref) — but it is wrapped in an
[`S2FaceSpace`](@ref) before the shared layer sees it, so an inadmissible
`nside` is an `ArgumentError` rather than a silently wrong polygon.

This is `S2.cell_polygon` — a function in the `S2` namespace, not a method of
the top-level `cell_polygon(::AbstractDGGS, level, id)`.
"""
cell_polygon(ix::Integer, iy::Integer, face::Integer, nside::Integer) =
    DGG.face_cell_polygon(S2FaceSpace(nside), ix, iy, face)
