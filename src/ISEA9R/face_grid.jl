# ---------------------------------------------------------------------------
# ISEA9R dense diamond grids — the ISEA9R instance of `src/core/face_grid.jl`
#
# `chart.jl` records where the ten continuous charts `[0, 1]² → S²` come from
# (they are `ISEA4R`'s, imported unchanged — the rhombus chart is
# aperture-agnostic) and adds the base-9 index maps aperture 9 needs. The
# package's shared face-grid layer (`src/core/face_grid.jl`) turns any such
# chart family into something a `ConservativeRegridding.Regridder` can consume:
# a *dense* ISEA9R grid — all `10 nside²` cells of one resolution — as a
# `SpatialTreeInterface` tree. All that lives here is what is ISEA9R about it:
# the system singleton, the two orderings, seven contract methods, the O(1)
# block-cap opt-in with its own pre-registered measurement, and the aliases this
# module's users type.
#
# Provenance: this file is a near-verbatim sibling of `src/ISEA4R/face_grid.jl`,
# section for section, with `2^k → 3^k`, `radix 4 → radix 9` and the `max_nside`
# derivation redone from scratch (see below). That the two files diff to almost
# nothing is the intended reading: **ISEA9R is ISEA4R's chart at a different
# aperture**, which is what `docs/design/isea9r_layout.md` establishes from
# primary sources.
#
# One thing deliberately NOT inherited: the cap policy. `FourCornerCap` is a
# per-system soundness obligation under the shared layer's contract, and a new
# system has to earn it again rather than borrow the sibling's argument. It was
# earned by re-running the sibling's pre-registered measurement at the
# aperture-9 block shapes — see `cap_policy` below.
# ---------------------------------------------------------------------------

import ..DiscreteGlobalGrids as DGG
# Imported, not `using`: the first three are extended below with ISEA9R methods,
# and the import is also what makes `using DiscreteGlobalGrids.ISEA9R:
# data_index` resolve. `num_cells` is the package generic the shared layer
# defines the `Isea9rFaceSpace` / `Isea9rFaceGrid` methods of.
import ..DiscreteGlobalGrids: data_index, lattice_index, validate_ordering, num_cells

# ---------------------------------------------------------------------------
# The system
# ---------------------------------------------------------------------------

"""
    Isea9rFaceSystem()

ISEA9R as a [`DiscreteGlobalGrids.FaceGridSystem`](@ref): ten diamond charts
`[0, 1]² → S²` (`ISEA4R.xyd_to_point`, imported — see `chart.jl`) over an
`nside × nside` lattice per diamond.

This is the *chart* system, deliberately distinct from the `ISEA9RDGGS`
registry entry: the two disagree about what a resolution is — `max_nside` here
is `3^18`, the point where the lattice arithmetic runs out of the room the
borrowed chart's own argument covers, while `max_level(ISEA9RDGGS()) ===
nothing` because the id hierarchy is unbounded and unported — and only this one
is a legal parameter of `FaceGrid` and friends.

It is also deliberately distinct from [`ISEA4R.Isea4rFaceSystem`](@ref), even
though the two share every geometry method: they carry different orderings
(radix 9 against radix 4), different `max_nside` bounds, and a `FaceGrid`
parameterised by one must not silently accept the other's ordering. What they do
share is checked rather than assumed — `test/ISEA9R/test_delegation.jl` asserts
that `Isea4rFaceGrid(3^k)` and `Isea9rFaceGrid(3^k)` emit bitwise-identical
polygons.
"""
struct Isea9rFaceSystem <: DGG.FaceGridSystem end

# ---------------------------------------------------------------------------
# The resolution
#
# Defined here, ahead of the orderings and the contract methods, because it is
# what the resolution-bound half of the shared contract dispatches on — see the
# "Resolution travels as a `FaceGridSpace`" section of `src/core/face_grid.jl`.
# ---------------------------------------------------------------------------

"""
    Isea9rFaceSpace(nside)

The resolution of a dense ISEA9R grid: `nside` cells along each edge of each of
the ten icosahedron diamonds, hence `10 * nside^2` cells in total.

Any `1 <= nside <= 3^18` is admissible — this is the *chart* resolution, not an
aperture-9 refinement level. The power-of-three restriction belongs to the
base-9 Morton index alone, and is enforced by [`MortonOrder`](@ref) when (and
only when) a grid is built with it. `Isea9rFaceSpace(4)` and
`Isea9rFaceSpace(5)` are perfectly good ISEA9R grids that simply have no
aperture-9 id space.

```julia
space = Isea9rFaceSpace(3)      # 90 cells
space.nside                     # 3
```

`num_cells(space)` is `10 * nside^2` — the length of the data vector any
[`AbstractIsea9rOrdering`](@ref) indexes.

This is also the value the resolution-bound contract methods below dispatch on
(`face_cell_corners`, `data_index`, `lattice_index`, `validate_ordering`): it
carries a *checked* `nside`, so no path into the layer can supply one the chart
maps cannot be evaluated at. Each shim unwraps `space.nside` and calls the
raw-Int codec — `ISEA4R`'s for the geometry and the row-major index, this
module's for the base-9 Morton index (see `chart.jl`).

See also [`Isea9rFaceGrid`](@ref).
"""
const Isea9rFaceSpace = DGG.FaceGridSpace{Isea9rFaceSystem}

# ---------------------------------------------------------------------------
# Orderings
#
# The lattice `(ix, iy, diamond)` says *where* a cell is; an ordering says
# *which slot of the data vector* it occupies. Keeping the two apart is what
# lets one grid type serve row-major layouts, Morton-ordered layouts, and
# (later) a custom permutation, without the tree code knowing which.
#
# These are per-system types, distinct from `ISEA4R.RowMajorOrder` /
# `ISEA4R.MortonOrder` even where the arithmetic behind them is shared: the
# shared layer checks `ordering isa ordering_family(Sys())` at construction, so
# handing an `Isea9rFaceGrid` an ISEA4R ordering is an `ArgumentError` rather
# than a grid whose Morton index quietly means radix 4.
# ---------------------------------------------------------------------------

"""
    abstract type AbstractIsea9rOrdering

How cells of a dense ISEA9R grid map to positions in a data vector — the ISEA9R
branch of [`DiscreteGlobalGrids.AbstractFaceOrdering`](@ref), which carries the
full extension contract — `data_index(o, space, ix, iy, diamond)`,
`lattice_index(o, space, j)` and the optional `validate_ordering(o, space)`,
where `space` is an [`Isea9rFaceSpace`](@ref), `ix`, `iy` are 0-based in
`0:space.nside-1` and the diamond is 0-based in `0:9`.

An [`Isea9rFaceGrid`](@ref) is a lattice (`10 nside²` cells addressed by
`(ix, iy, diamond)`) plus one of these; the ordering is the *only* thing that
decides which column of a `ConservativeRegridding.Regridder` a cell lands in.
[`RowMajorOrder`](@ref) and [`MortonOrder`](@ref) are the two shipped instances.

**A DGGAL-parity ordering is deliberately not among them.** DGGAL's
within-rhombus index is row-major over a *transposed* square under a
non-identity root→diamond permutation, and both the permutation and the
transpose are derivations with no fixture behind them; shipping them would be
compatibility theatre. See `docs/design/isea9r_layout.md` §6-7 for the three
deltas and the fixture dump that would settle them.
"""
abstract type AbstractIsea9rOrdering <: DGG.AbstractFaceOrdering end

"""
    RowMajorOrder()

Row-major data ordering: position `j` holds the cell whose 0-based row-major id
is `j - 1`, i.e. `ix` fastest, then `iy`, then `diamond`.

Defined for **any** `nside >= 1` — the row-major closed forms
(`ISEA4R.xyd_to_rowmajor` / `ISEA4R.rowmajor_to_xyd`, imported by `chart.jl`)
carry no aperture at all, which is what makes this the default and the ordering
to reach for at non-power-of-three resolutions.

This is a *distinct type* from [`ISEA4R.RowMajorOrder`](@ref) — the shared layer
type-checks orderings against their system — but it computes the identical
index, since row-major over `10 nside²` cells is one arithmetic expression with
no aperture in it. At a common `nside = 3^k` the two systems' row-major grids
are therefore the same grid, cell for cell and bit for bit; that is asserted
rather than assumed in `test/ISEA9R/test_delegation.jl`.
"""
struct RowMajorOrder <: AbstractIsea9rOrdering end

# `+ 1` because row-major ids are 0-based while data positions are 1-based; see
# the index-convention block at the top of `chart.jl`.
data_index(::RowMajorOrder, space::Isea9rFaceSpace, ix::Integer, iy::Integer, diamond::Integer) =
    Int(xyd_to_rowmajor(ix, iy, diamond, space.nside)) + 1

lattice_index(::RowMajorOrder, space::Isea9rFaceSpace, j::Integer) =
    rowmajor_to_xyd(j - 1, space.nside)

"""
    MortonOrder()

Base-9 Morton (Z-order) data ordering: position `j` holds the cell whose 0-based
Morton id is `j - 1`, where the id is `diamond * nside² + morton9(ix, iy)` —
exactly the ordinal shape the `ISEA9RDGGS` registry entry records
(`diamond * 9^level + position`, see `src/core/systems/isea4r_isea9r.jl`). This
is to ISEA9R what `ISEA4R.MortonOrder` is to ISEA4R, `NestedOrder` is to HEALPix
and `HilbertOrder` is to S2.

Requires `nside = 3^k`: the base-9 Morton code interleaves one base-3 digit of
`ix` with one of `iy` per level, which only exists on a `3^k × 3^k` diamond.
[`Isea9rFaceGrid`](@ref) rejects any other `nside` at construction rather than
at first query.

The radix-9 prefix arithmetic is exact across levels because the within-diamond
lattice nesting is bit-exact: `fl(ix/n) === fl(3ix/3n)` — numerator and
denominator are integers exactly representable as `Float64` on both sides of
that equation and the real quotient is identical, so correctly-rounded division
returns the identical `Float64` — and hence the cell at `nside = 3n` really is a
ninth of its parent at `nside = n`, corner for corner (see
[`ISEA4R.xyd_to_point`](@ref)). It is this exactness, over this ordinal, that
`supports_prefix_ranges(ISEA9RDGGS()) == true` is a statement about.

What this ordering does *not* claim is compatibility with any DGGRID or DGGAL
ISEA9R product: `position` is hereby pinned to the base-9 Morton code in *this
package's* diamond layout, that layout has no external oracle behind it, and
DGGAL's own within-rhombus index is a different function entirely (row-major
over a transposed square, under a root permutation). See
`docs/design/isea9r_layout.md`.
"""
struct MortonOrder <: AbstractIsea9rOrdering end

data_index(::MortonOrder, space::Isea9rFaceSpace, ix::Integer, iy::Integer, diamond::Integer) =
    Int(xyd_to_morton(ix, iy, diamond, space.nside)) + 1

lattice_index(::MortonOrder, space::Isea9rFaceSpace, j::Integer) =
    morton_to_xyd(j - 1, space.nside)

# The one thing that must be caught eagerly: a Morton-ordered grid at nside = 4
# would construct fine and then throw from deep inside a dual-tree traversal,
# which is a miserable place to learn about it.
function validate_ordering(::MortonOrder, space::Isea9rFaceSpace)
    ispow3(space.nside) || throw(ArgumentError(
        "MortonOrder requires nside = 3^k, got nside=$(space.nside); \
         use RowMajorOrder() for arbitrary nside"))
    return nothing
end

# ---------------------------------------------------------------------------
# The `FaceGridSystem` contract
# ---------------------------------------------------------------------------

# 10 diamonds. Unlike the ISEA4R sibling, this one is *not* a coincidence with
# `root_count(ISEA9RDGGS()) == 10`: both are the ten root rhombuses OGC
# 21-038r1 Annex B.2 names ("The ten root rhombuses are formed by combining two
# icosahedron triangles at their base") and DGGAL's `countZones(0) = 10 * 9^0`
# computes. The literal still stays independent of the registry trait, because
# the shared layer's contract says these are chart-arithmetic facts.
DGG.nfaces(::Isea9rFaceSystem) = 10

# The chart, and the corner function that evaluates it on the lattice: ISEA4R's,
# imported by `chart.jl`. The seam-ownership rule (the half `y >= x` evaluates
# through the upper face) and the CCW corner order travel with them, so this
# system inherits *exactly* the bit-identity properties `ISEA4R` documents —
# there is no second evaluation path that could disagree.
DGG.face_chart(::Isea9rFaceSystem, x, y, diamond::Integer) = xyd_to_point(x, y, diamond)

DGG.face_cell_corners(space::Isea9rFaceSpace, ix::Integer, iy::Integer, diamond::Integer) =
    cell_corners(ix, iy, diamond, space.nside)

# Upper bound: `3^18 = 387_420_489`. Three constraints, and the derivation is
# which one binds — this is NOT the sibling's `2^29` carried over.
#
#  1. *Cell-count and id arithmetic must fit `Int64`.* `10 * nside^2` is the
#     grid size and also bounds both id spaces (the base-9 Morton code of a cell
#     is `< nside²`, so ids are `< 10 nside²`). `10 * n^2 <= typemax(Int64)`
#     holds up to `n = 960_383_883`-ish; at `n = 3^19 = 1_162_261_467` the
#     product is `1.35e19` and wraps silently. Not the binding constraint, but
#     it is the one that rules out `3^19`.
#
#  2. *The borrowed chart's own domain.* This system evaluates `ISEA4R`'s chart,
#     whose seam-branch exactness argument — `Float64(iy/nside) >=
#     Float64(ix/nside)` decides the branch identically to the integer predicate
#     `iy >= ix` — is stated for `0 <= ix, iy <= nside <= 2^29`. A delegating
#     system may not outrun the argument it delegates to, so `2^29 =
#     536_870_912` is a hard ceiling here whatever the count arithmetic allows.
#
#  3. *Cross-level nesting must stay bit-exact.* The radix-9 prefix property
#     rests on `fl(ix/3^k) === fl(3ix/3^(k+1))`, which holds whenever all four
#     integers are exactly representable as `Float64`, i.e. `3^(k+1) <= 2^53`,
#     i.e. `k <= 32`. Never binding at these sizes — `3^19 ≈ 1.2e9` is nine
#     orders under `2^53` — but it is the constraint that would bind first if
#     the other two were lifted, so it is recorded rather than assumed.
#
# The bound is then the largest power of three under (2): `3^18` (`3^19` fails
# both (1) and (2)). A power of three rather than the raw wrap point of (1)
# because `3^k` is the only resolution at which ISEA9R has an id space at all —
# `MortonOrder`, the ordinal, and every hierarchy statement live there — so
# admitting `nside = 700_000_000` would buy a row-major-only grid nobody can
# name a cell of, at the cost of a bound that no longer reads as a resolution.
# At `3^18` the grid has `10 * 3^36 = 1.50e18` cells, comfortably inside `Int64`.
DGG.max_nside(::Isea9rFaceSystem) = 3^18

DGG.default_ordering(::Isea9rFaceSystem) = RowMajorOrder()

DGG.ordering_family(::Isea9rFaceSystem) = AbstractIsea9rOrdering

DGG.facegrid_prefix(::Isea9rFaceSystem) = "Isea9r"

# The cap of a block's four corners contains the whole block, so we can skip the
# generic method's inclusion of all perimeter vertices (`cell_range_extent`,
# which is O(perimeter) per node). The O(1) form is the shared layer's
# `FourCornerCap`.
#
# The shared layer requires this override to carry its OWN soundness argument —
# proof or pre-registered measurement — per system, and that requirement is not
# satisfied by pointing at `cap_policy(::Isea4rFaceSystem)`: the chart is the
# same, but the blocks a cursor builds over a `3^k × 3^k` lattice are not the
# blocks it builds over a `2^k × 2^k` one, and it is the block shapes that the
# measurement is about. So the sibling's rule was re-run at this system's
# resolutions, before the numbers were in:
#
#   Rule (pre-registered, ISEA4R's verbatim with the resolutions changed): over
#   every block the cursor can create — every range-bisection node from
#   `(1:n, 1:n)` down, plus every 1x1 leaf range — on every diamond, at
#   `nside ∈ (3, 9, 27)`, compare a dense 17x17 sampling of the block's chart
#   rectangle against the *pre-slack* radius of this cap (its four corners and
#   their slerp midpoints, without the `1.0001` inflation
#   `circle_from_four_corners` applies). Adopt the O(1) override if the worst
#   overhang is <= 0, or if the worst ratio is a stable <= 5%.
#
#   Measured: worst pre-slack overhang is EXACTLY `0.0` at all three `nside`
#   (13, 117 and 1045 blocks per diamond respectively) — the farthest sampled
#   point of every block is one of its own four corners, including for
#   seam-straddling blocks and whole-diamond roots. Adopted, with the stock
#   `1.0001` slack and no extra inflation.
#
# The geometric reading is the sibling's, and it is aperture-free, which is why
# the aperture-9 numbers came out identical: the square → planar-rhombus map is
# globally affine (the two glued equilateral halves unfold to a single
# parallelogram), so every block is the Snyder image of a planar parallelogram;
# Snyder's mild equal-area radial compression keeps the corners extremal
# relative to the corner-mean cap centre.
#
# If the standing assert in `test/ISEA9R/test_face_grid.jl` ever fires, the fix
# is to change this method to return `DGG.PerimeterWalkCap()`: the stock
# perimeter walk is sound and merely slower.
DGG.cap_policy(::Isea9rFaceSystem) = DGG.FourCornerCap()

# ---------------------------------------------------------------------------
# The ISEA9R names for the shared types
# ---------------------------------------------------------------------------

"""
    Isea9rFaceGrid(space::Isea9rFaceSpace, ordering::AbstractIsea9rOrdering)
    Isea9rFaceGrid(nside::Integer; ordering = RowMajorOrder())

A complete ISEA9R grid at one resolution — all `10 nside²` cells — together with
the data ordering its cells are numbered by. `treeify(grid)` turns it into a
spatial tree ([`Isea9rFaceRoot`](@ref)) that
`ConservativeRegridding.Regridder` consumes directly.

```julia
grid = Isea9rFaceGrid(9; ordering = MortonOrder())
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
column `j` is `isea9r_ordinal` `j - 1` (`diamond * 9^level + position` at
`nside = 3^level`).

# Validation

`nside >= 1` is checked by [`Isea9rFaceSpace`](@ref); the ordering gets a say
too, through `validate_ordering` — which is why
`Isea9rFaceGrid(4; ordering = MortonOrder())` throws an `ArgumentError` here
rather than failing later inside a traversal.
"""
const Isea9rFaceGrid = DGG.FaceGrid{Isea9rFaceSystem}

"""
    DiamondChartGrid(manifold, space::Isea9rFaceSpace, diamond, ordering)

One ISEA9R diamond as an `nside × nside` `Trees.AbstractCurvilinearGrid`, so
that `Trees.TopDownQuadtreeCursor` can index it with no ISEA9R-specific descent
logic. Internal: build an [`Isea9rFaceGrid`](@ref) and `treeify` it.

Named `DiamondChartGrid` rather than `FaceChartGrid` for the reason the ISEA4R
sibling is: in an icosahedral module "face" already means "icosahedron face", of
which a diamond has two. The *field* keeps the shared layer's generic vocabulary
— it is `g.face`, 0-based in `0:9`, and it holds the diamond index. (This is
`ISEA9R.DiamondChartGrid`, a different `const` from `ISEA4R.DiamondChartGrid`:
they alias `FaceChartGrid` at different system parameters.)

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
const DiamondChartGrid = DGG.FaceChartGrid{Isea9rFaceSystem}

"""
    Isea9rFaceRoot(manifold, space::Isea9rFaceSpace, ordering)
    Isea9rFaceRoot(nside::Integer, ordering = RowMajorOrder())

Root of the spatial tree over a dense ISEA9R grid: the whole sphere, with the
ten icosahedron diamonds as children. Each child is a
`Trees.TopDownQuadtreeCursor` over that diamond's [`DiamondChartGrid`](@ref);
leaf indices throughout are data-vector positions under `ordering`, so
`Trees.getcell(root, j)` is the cell of data position `j` and a `Regridder`'s
column `j` is that same position.

Build one with `treeify(::Isea9rFaceGrid)` rather than by hand. Constructed
directly it re-runs the two ordering checks `Isea9rFaceGrid` runs, so
`Isea9rFaceRoot(4, MortonOrder())` throws here rather than from inside a
traversal.
"""
const Isea9rFaceRoot = DGG.FaceGridRoot{Isea9rFaceSystem}

# ---------------------------------------------------------------------------
# Cell geometry
# ---------------------------------------------------------------------------

"""
    cell_polygon(ix, iy, diamond, nside) -> GI.Polygon

Cell `(ix, iy)` of `diamond` as a closed 4-gon on the unit sphere, from
[`ISEA4R.cell_corners`](@ref) (imported — the chart is ISEA4R's, see
`chart.jl`).

Snyder cell edges are not great circles (unlike S2's, and like HEALPix's), so
this 4-gon describes the cell only approximately — densify through
[`ISEA4R.xyd_to_point`](@ref) along the edges if more is needed. Since every
corner is a shared lattice point evaluated by the same chart function under the
same seam-ownership rule, neighbouring polygons meet at bit-identical vertices
and the tessellation is exact.

The ring is **counterclockwise as seen from outside the sphere**, which
`cell_corners` guarantees and which is a hard contract rather than a convention:
the convex-clip kernel that computes spherical intersections clips a clockwise
ring to EMPTY, so a reversed ring yields silent zero areas instead of an error.

`nside` is taken as a plain integer here — this is the chart-side vocabulary
wrapper — but it is wrapped in an [`Isea9rFaceSpace`](@ref) before the shared
layer sees it, so an inadmissible `nside` is an `ArgumentError` rather than a
silently wrong polygon. That wrapping is also the one place this function
differs from `ISEA4R.cell_polygon`: the two emit identical geometry wherever
both are admissible, but they bound `nside` by their own systems' `max_nside`.

This is `ISEA9R.cell_polygon` — a function in the `ISEA9R` namespace, not a
method of the top-level `cell_polygon(::AbstractDGGS, level, id)`.
"""
cell_polygon(ix::Integer, iy::Integer, diamond::Integer, nside::Integer) =
    DGG.face_cell_polygon(Isea9rFaceSpace(nside), ix, iy, diamond)
