# ---------------------------------------------------------------------------
# ISEA4R dense diamond grids — the ISEA4R instance of `src/core/face_grid.jl`
#
# `chart.jl` gives the ten continuous charts `[0, 1]² → S²` and the index maps
# over the `nside × nside` lattice. The package's shared face-grid layer
# (`src/core/face_grid.jl`) turns any such chart family into something a
# `ConservativeRegridding.Regridder` can consume: a *dense* ISEA4R grid — all
# `10 nside²` cells of one resolution — as a `SpatialTreeInterface` tree. All
# that lives here is what is ISEA4R about it: the system singleton, the two
# orderings, seven contract methods, the O(1) block-cap opt-in with its
# pre-registered measurement, and the aliases this module's users type.
#
# Provenance: this file and its HEALPix and S2 siblings were written as
# near-verbatim copies of one another, section for section — a deliberate
# rule-of-three duplication, since two instances is not enough to fix the shape
# of a shared abstraction. This third instance is the one that fixed it (it is
# the one exercising chart piecewiseness, seam determinism, and a cap argument
# by measurement rather than proof), and the shared layer is what the three
# files collapsed into. The one vocabulary rename survives: the shared
# `FaceChartGrid` is aliased here as `DiamondChartGrid`, because in an
# icosahedral module "face" already means "icosahedron face", of which a diamond
# has two.
# ---------------------------------------------------------------------------

import ..DiscreteGlobalGrids as DGG
# Imported, not `using`: the first three are extended below with ISEA4R methods,
# and the import is also what makes `using DiscreteGlobalGrids.ISEA4R:
# data_index` resolve. `num_cells` is the package generic the shared layer
# defines the `Isea4rFaceSpace` / `Isea4rFaceGrid` methods of.
import ..DiscreteGlobalGrids: data_index, lattice_index, validate_ordering, num_cells

# ---------------------------------------------------------------------------
# The system
# ---------------------------------------------------------------------------

"""
    Isea4rFaceSystem()

ISEA4R as a [`DiscreteGlobalGrids.FaceGridSystem`](@ref): ten diamond charts
`[0, 1]² → S²` (`xyd_to_point`) over an `nside × nside` lattice per diamond.

This is the *chart* system, deliberately distinct from the `ISEA4RDGGS`
registry entry: the two disagree about what a resolution is — `max_nside` here
is `2^29`, the point where the lattice arithmetic overflows, while
`max_level(ISEA4RDGGS()) === nothing` because the id hierarchy is unbounded and
unported — and only this one is a legal parameter of `FaceGrid` and friends.
"""
struct Isea4rFaceSystem <: DGG.FaceGridSystem end

# ---------------------------------------------------------------------------
# The resolution
#
# Defined here, ahead of the orderings and the contract methods, because it is
# what the resolution-bound half of the shared contract dispatches on — see the
# "Resolution travels as a `FaceGridSpace`" section of `src/core/face_grid.jl`.
# ---------------------------------------------------------------------------

"""
    Isea4rFaceSpace(nside)

The resolution of a dense ISEA4R grid: `nside` cells along each edge of each of
the ten icosahedron diamonds, hence `10 * nside^2` cells in total.

Any `nside >= 1` is admissible — this is the *chart* resolution, not an
aperture-4 refinement level. The power-of-two restriction belongs to the Morton
index alone, and is enforced by [`MortonOrder`](@ref) when (and only when) a
grid is built with it. `Isea4rFaceSpace(3)` and `Isea4rFaceSpace(5)` are
perfectly good ISEA4R grids that simply have no aperture-4 id space.

```julia
space = Isea4rFaceSpace(4)      # 160 cells
space.nside                     # 4
```

`num_cells(space)` is `10 * nside^2` — the length of the data vector any
[`AbstractIsea4rOrdering`](@ref) indexes. (Mirrors HEALPix's `num_pixels` and
S2's `num_cells`.)

This is also the value the resolution-bound contract methods below dispatch on
(`face_cell_corners`, `data_index`, `lattice_index`, `validate_ordering`): it
carries a *checked* `nside`, so no path into the layer can supply one the chart
maps cannot be evaluated at. Each shim unwraps `space.nside` and calls the
raw-Int codec in `chart.jl`, which stays shared with `Isea4rKernel.jl`.

See also [`Isea4rFaceGrid`](@ref).
"""
const Isea4rFaceSpace = DGG.FaceGridSpace{Isea4rFaceSystem}

# ---------------------------------------------------------------------------
# Orderings
#
# The lattice `(ix, iy, diamond)` says *where* a cell is; an ordering says
# *which slot of the data vector* it occupies. Keeping the two apart is what
# lets one grid type serve row-major layouts, Morton-ordered layouts, and
# (later) a custom permutation, without the tree code knowing which.
# ---------------------------------------------------------------------------

"""
    abstract type AbstractIsea4rOrdering

How cells of a dense ISEA4R grid map to positions in a data vector — the ISEA4R
branch of [`DiscreteGlobalGrids.AbstractFaceOrdering`](@ref), which carries the
full extension contract — `data_index(o, space, ix, iy, diamond)`,
`lattice_index(o, space, j)` and the optional `validate_ordering(o, space)`,
where `space` is an [`Isea4rFaceSpace`](@ref), `ix`, `iy` are 0-based in
`0:space.nside-1` and the diamond is 0-based in `0:9`.

An [`Isea4rFaceGrid`](@ref) is a lattice (`10 nside²` cells addressed by
`(ix, iy, diamond)`) plus one of these; the ordering is the *only* thing that
decides which column of a `ConservativeRegridding.Regridder` a cell lands in.
[`RowMajorOrder`](@ref) and [`MortonOrder`](@ref) are the two shipped instances;
an SST-compatible permutation (see the recorded decision in `chart.jl`) or the
column order of an on-disk product can be added against the contract without
touching any grid or tree code.
"""
abstract type AbstractIsea4rOrdering <: DGG.AbstractFaceOrdering end

"""
    RowMajorOrder()

Row-major data ordering: position `j` holds the cell whose 0-based row-major id
is `j - 1`, i.e. `ix` fastest, then `iy`, then `diamond`.

Defined for **any** `nside >= 1` — the row-major closed forms
([`xyd_to_rowmajor`](@ref) / [`rowmajor_to_xyd`](@ref)) carry no power-of-two
restriction, which is what makes this the default and the ordering to reach for
at non-power-of-two resolutions.

This layout is already isomorphic to SphericalSpatialTrees.jl's
`LinearIndices((2^r, 2^r, 10))` order (first index fastest); what is *not*
pinned against SST is the diamond numbering and per-diamond axis orientation, so
no interoperability is claimed. See `diamonds.jl` and
`docs/design/isea4r_diamond_layout.md`.
"""
struct RowMajorOrder <: AbstractIsea4rOrdering end

# `+ 1` because row-major ids are 0-based while data positions are 1-based; see
# the index-convention block at the top of `chart.jl`.
data_index(::RowMajorOrder, space::Isea4rFaceSpace, ix::Integer, iy::Integer, diamond::Integer) =
    Int(xyd_to_rowmajor(ix, iy, diamond, space.nside)) + 1

lattice_index(::RowMajorOrder, space::Isea4rFaceSpace, j::Integer) =
    rowmajor_to_xyd(j - 1, space.nside)

"""
    MortonOrder()

Morton (Z-order) data ordering: position `j` holds the cell whose 0-based Morton
id is `j - 1`, where the id is `diamond * nside² + morton(ix, iy)` — exactly the
scaffold ordinal shape the `ISEA4RDGGS` registry entry records
(`diamond * 4^level + position`, see `src/core/systems/isea4r_isea9r.jl`). This
is to ISEA4R what `NestedOrder` is to HEALPix and `HilbertOrder` is to S2.

Requires `nside = 2^k`: the Morton code interleaves one bit of `ix` with one bit
of `iy` per level, which only exists on a `2^k × 2^k` diamond.
[`Isea4rFaceGrid`](@ref) rejects any other `nside` at construction rather than
at first query.

The radix-4 prefix arithmetic is exact across levels because the within-diamond
lattice nesting is bit-exact: `fl(ix/n) === fl(2ix/2n)`, so the cell at
`nside = 2n` really is a quarter of its parent at `nside = n`, corner for corner
(see [`xyd_to_point`](@ref)). What this ordering does *not* claim is
compatibility with any external ISEA4R product, DGGAL included: `position` is
hereby pinned to Morton in *this package's* diamond layout, and that layout has
no external oracle behind it.
"""
struct MortonOrder <: AbstractIsea4rOrdering end

data_index(::MortonOrder, space::Isea4rFaceSpace, ix::Integer, iy::Integer, diamond::Integer) =
    Int(xyd_to_morton(ix, iy, diamond, space.nside)) + 1

lattice_index(::MortonOrder, space::Isea4rFaceSpace, j::Integer) =
    morton_to_xyd(j - 1, space.nside)

# The one thing that must be caught eagerly: a Morton-ordered grid at nside = 3
# would construct fine and then throw from deep inside a dual-tree traversal,
# which is a miserable place to learn about it.
function validate_ordering(::MortonOrder, space::Isea4rFaceSpace)
    ispow2(space.nside) || throw(ArgumentError(
        "MortonOrder requires nside = 2^k, got nside=$(space.nside); \
         use RowMajorOrder() for arbitrary nside"))
    return nothing
end

# ---------------------------------------------------------------------------
# The `FaceGridSystem` contract
# ---------------------------------------------------------------------------

# 10 diamonds. (Numerically `root_count(ISEA4RDGGS())` too — a coincidence of
# this system, not a contract, so the literal stays independent.)
DGG.nfaces(::Isea4rFaceSystem) = 10

DGG.face_chart(::Isea4rFaceSystem, x, y, diamond::Integer) = xyd_to_point(x, y, diamond)

DGG.face_cell_corners(space::Isea4rFaceSpace, ix::Integer, iy::Integer, diamond::Integer) =
    cell_corners(ix, iy, diamond, space.nside)

# Upper bound: `10 * nside^2` is `10 * 2^58 ≈ 2.9e18 < typemax(Int64)` at
# `nside = 2^29` and wraps silently at `2^30`. The Morton id space needs the
# same room (`10 * 4^29 - 1 < 2^62`). Same bound and same rationale as
# `HealpixFaceSpace`, deliberately — it keeps the two layers diff-able.
DGG.max_nside(::Isea4rFaceSystem) = 2^29

DGG.default_ordering(::Isea4rFaceSystem) = RowMajorOrder()

DGG.ordering_family(::Isea4rFaceSystem) = AbstractIsea4rOrdering

DGG.facegrid_prefix(::Isea4rFaceSystem) = "Isea4r"

# The cap of a block's four corners contains the whole block, so we can skip the
# generic method's inclusion of all perimeter vertices (`cell_range_extent`,
# which is O(perimeter) per node). The O(1) form is the shared layer's
# `FourCornerCap`, and the soundness obligation selecting it carries was
# discharged by measurement rather than by proof — the decision rule was
# pre-registered before the numbers were in, and the measurement ships as a
# standing test:
#
#   Rule (pre-registered): over every block the cursor can create — every
#   range-bisection node from `(1:n, 1:n)` down, plus every 1x1 leaf range — on
#   every diamond, at `nside ∈ (1, 2, 3, 4, 5, 8, 16)`, compare a dense 17x17
#   sampling of the block's chart rectangle against the *pre-slack* radius of
#   this cap (its four corners and their slerp midpoints, without the `1.0001`
#   inflation `circle_from_four_corners` applies). Adopt the O(1) override if
#   the worst overhang is <= 0, or if the worst ratio is a stable <= 5%.
#
#   Measured: worst pre-slack overhang is EXACTLY `0.0` at every `nside` — the
#   farthest sampled point of every block is one of its own four corners,
#   including for seam-straddling blocks and whole-diamond roots. Adopted, with
#   the stock `1.0001` slack and no extra inflation.
#
# The geometric reading: the square → planar-rhombus map is globally affine (the
# two glued equilateral halves unfold to a single parallelogram, and the
# point-reflection identity between them holds exactly), so every block is the
# Snyder image of a planar parallelogram; Snyder's mild equal-area radial
# compression keeps the corners extremal relative to the corner-mean cap centre.
#
# If the standing assert in `test/ISEA4R/test_face_grid.jl` ever fires, the fix
# is to change this method to return `DGG.PerimeterWalkCap()`: the stock
# perimeter walk is sound and merely slower.
DGG.cap_policy(::Isea4rFaceSystem) = DGG.FourCornerCap()

# ---------------------------------------------------------------------------
# The ISEA4R names for the shared types
# ---------------------------------------------------------------------------

"""
    Isea4rFaceGrid(space::Isea4rFaceSpace, ordering::AbstractIsea4rOrdering)
    Isea4rFaceGrid(nside::Integer; ordering = RowMajorOrder())

A complete ISEA4R grid at one resolution — all `10 nside²` cells — together with
the data ordering its cells are numbered by. `treeify(grid)` turns it into a
spatial tree ([`Isea4rFaceRoot`](@ref)) that
`ConservativeRegridding.Regridder` consumes directly.

```julia
grid = Isea4rFaceGrid(4; ordering = MortonOrder())
R = ConservativeRegridding.Regridder(treeify(grid), treeify(other))
```

# Alignment

**Column `j` of a `Regridder` built on `treeify(grid)` is data-vector position
`j` under `ordering` — there is no permutation between the matrix and the
field.** This is structural, not a convention to remember: the per-diamond
cursors emit `data_index(ordering, ...)` as their leaf indices and
`Trees.getcell(root, j)` returns the polygon of `lattice_index(ordering, space,
j)`, so the matrix is *assembled* in data order. Concretely, with
`RowMajorOrder` column `j` is row-major id `j - 1`, and with `MortonOrder`
column `j` is scaffold ordinal `j - 1` (`diamond * 4^level + position` at
`nside = 2^level`).

# Validation

`nside >= 1` is checked by [`Isea4rFaceSpace`](@ref); the ordering gets a say
too, through `validate_ordering` — which is why
`Isea4rFaceGrid(3; ordering = MortonOrder())` throws an `ArgumentError` here
rather than failing later inside a traversal.
"""
const Isea4rFaceGrid = DGG.FaceGrid{Isea4rFaceSystem}

"""
    DiamondChartGrid(manifold, space::Isea4rFaceSpace, diamond, ordering)

One ISEA4R diamond as an `nside × nside` `Trees.AbstractCurvilinearGrid`, so
that `Trees.TopDownQuadtreeCursor` can index it with no ISEA4R-specific descent
logic. Internal: build an [`Isea4rFaceGrid`](@ref) and `treeify` it.

Named `DiamondChartGrid` rather than `FaceChartGrid` (the name both siblings
use, and the name of the shared type this aliases) for one reason only: in an
icosahedral module "face" already means "icosahedron face", of which a diamond
has two. The *field* keeps the shared layer's generic vocabulary — it is
`g.face`, 0-based in `0:9`, and it holds the diamond index.

Cartesian cell index `(i, j)` (1-based) is lattice cell `(i-1, j-1)`, and
`Trees.getcell` returns its CCW polygon (see [`cell_polygon`](@ref) — CCW is
required by the convex-clip kernel, which clips a CW ring to EMPTY).
`Trees.getvertex(g, i, j)` is the *lattice* point `((i-1)/nside,
(j-1)/nside)`, `1:(nside+1)` in each direction, which is what drives the
bounding caps.

The two index maps are overridden away from the generic column-major default:
they target the *global* data layout of `ordering`, not a diamond-local one, so
the leaf indices this grid's cursor reports are already data-vector positions.
"""
const DiamondChartGrid = DGG.FaceChartGrid{Isea4rFaceSystem}

"""
    Isea4rFaceRoot(manifold, space::Isea4rFaceSpace, ordering)
    Isea4rFaceRoot(nside::Integer, ordering = RowMajorOrder())

Root of the spatial tree over a dense ISEA4R grid: the whole sphere, with the
ten icosahedron diamonds as children. Each child is a
`Trees.TopDownQuadtreeCursor` over that diamond's [`DiamondChartGrid`](@ref);
leaf indices throughout are data-vector positions under `ordering`, so
`Trees.getcell(root, j)` is the cell of data position `j` and a `Regridder`'s
column `j` is that same position.

Build one with `treeify(::Isea4rFaceGrid)` rather than by hand. Constructed
directly it re-runs the two ordering checks `Isea4rFaceGrid` runs, so
`Isea4rFaceRoot(3, MortonOrder())` throws here rather than from inside a
traversal.
"""
const Isea4rFaceRoot = DGG.FaceGridRoot{Isea4rFaceSystem}

# ---------------------------------------------------------------------------
# Cell geometry
# ---------------------------------------------------------------------------

"""
    cell_polygon(ix, iy, diamond, nside) -> GI.Polygon

Cell `(ix, iy)` of `diamond` as a closed 4-gon on the unit sphere, from
[`cell_corners`](@ref).

Snyder cell edges are not great circles (unlike S2's, and like HEALPix's), so
this 4-gon describes the cell only approximately — densify through
[`xyd_to_point`](@ref) along the edges if more is needed. Since every corner is
a shared lattice point evaluated by the same chart function under the same
seam-ownership rule, neighbouring polygons meet at bit-identical vertices and
the tessellation is exact.

The ring is **counterclockwise as seen from outside the sphere**, which
`cell_corners` guarantees and which is a hard contract rather than a convention:
the convex-clip kernel that computes spherical intersections clips a clockwise
ring to EMPTY, so a reversed ring yields silent zero areas instead of an error.

`nside` is taken as a plain integer here — this is the chart-side vocabulary
wrapper, alongside [`cell_corners`](@ref) — but it is wrapped in an
[`Isea4rFaceSpace`](@ref) before the shared layer sees it, so an inadmissible
`nside` is an `ArgumentError` rather than a silently wrong polygon.

This is `ISEA4R.cell_polygon` — a function in the `ISEA4R` namespace, not a
method of the top-level `cell_polygon(::AbstractDGGS, level, id)`.
"""
cell_polygon(ix::Integer, iy::Integer, diamond::Integer, nside::Integer) =
    DGG.face_cell_polygon(Isea4rFaceSpace(nside), ix, iy, diamond)
