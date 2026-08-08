# ---------------------------------------------------------------------------
# HEALPix dense face grids — the HEALPix instance of `src/core/face_grid.jl`
#
# `chart.jl` gives the twelve closed-form charts `[0, 1]² → S²` and the index
# maps over the `nside × nside` lattice. The package's shared face-grid layer
# (`src/core/face_grid.jl`) turns any such chart family into something a
# `ConservativeRegridding.Regridder` can consume: a *dense* HEALPix grid — all
# `12 nside²` pixels of one resolution — as a `SpatialTreeInterface` tree. All
# that lives here is what is HEALPix about it: the system singleton, the two
# orderings, seven contract methods, and the aliases and vocabulary wrappers
# (`num_pixels`, `pixel_polygon`) this module's users type.
#
# Provenance: this layer is a port of ConservativeRegridding.jl's RingGrids
# extension (`ext/ConservativeRegriddingRingGridsExt/healpix.jl`, the
# `HEALPixFaceGrid` / `HEALPixRootNode` pair). The differences are deliberate
# and are the point of the port:
#
#   * the closed forms live in `chart.jl` and are shared with everything else
#     here, rather than being inlined into the tree file;
#   * the data ordering is a component instead of being hard-wired to RING, so
#     the same grid serves ring-ordered and nested-ordered fields (and, later,
#     a permutation/space-filling ordering) with no new tree code.
#
# The block-extent method for a block of pixels and the 12-child root were ported
# essentially verbatim, comments included, and have since been hoisted into the
# shared layer — the HEALPix-specific part of the first is the `cap_policy`
# override below.
# ---------------------------------------------------------------------------

import ..DiscreteGlobalGrids as DGG
# Imported, not `using`: these are extended below with HEALPix methods, and the
# import is also what makes `using DiscreteGlobalGrids.HEALPix: data_index`
# resolve.
import ..DiscreteGlobalGrids: data_index, lattice_index, validate_ordering

# ---------------------------------------------------------------------------
# The system
# ---------------------------------------------------------------------------

"""
    HealpixFaceSystem()

HEALPix as a [`DiscreteGlobalGrids.FaceGridSystem`](@ref): twelve charts
`[0, 1]² → S²` (`xyf_to_point`) over an `nside × nside` lattice per face.

This is the *chart* system, deliberately distinct from the `HEALPixDGGS`
registry entry: the two disagree about what a resolution is (`max_nside` here
bounds lattice arithmetic; `max_level` there bounds the nested id hierarchy),
and only this one is a legal parameter of `FaceGrid` and friends.
"""
struct HealpixFaceSystem <: DGG.FaceGridSystem end

# ---------------------------------------------------------------------------
# The resolution
#
# Defined here, ahead of the orderings and the contract methods, because it is
# what the resolution-bound half of the shared contract dispatches on — see the
# "Resolution travels as a `FaceGridSpace`" section of `src/core/face_grid.jl`.
# ---------------------------------------------------------------------------

"""
    HealpixFaceSpace(nside)

The resolution of a dense HEALPix grid: `nside` pixels along each edge of each
of the twelve base faces, hence `12 * nside^2` pixels in total.

Any `nside >= 1` is admissible — this is the *chart* resolution, not a nested
refinement level. The power-of-two restriction belongs to the NESTED index
alone, and is enforced by [`NestedOrder`](@ref) when (and only when) a grid is
built with it. `HealpixFaceSpace(3)` and `HealpixFaceSpace(5)` are perfectly
good HEALPix grids that simply have no nested id space.

```julia
space = HealpixFaceSpace(4)     # 192 pixels
space.nside                     # 4
```

This is also the value the resolution-bound contract methods below dispatch on
(`face_cell_corners`, `data_index`, `lattice_index`, `validate_ordering`): it
carries a *checked* `nside`, so no path into the layer can supply one the chart
maps cannot be evaluated at. Each shim unwraps `space.nside` and calls the
raw-Int codec in `chart.jl`, which stays shared with `HealpixKernel.jl`.

See also [`HealpixFaceGrid`](@ref).
"""
const HealpixFaceSpace = DGG.FaceGridSpace{HealpixFaceSystem}

# ---------------------------------------------------------------------------
# Orderings
#
# The lattice `(ix, iy, face)` says *where* a pixel is; an ordering says *which
# slot of the data vector* it occupies. Keeping the two apart is what lets one
# grid type serve ring-ordered files, nested-ordered files, and (later) a
# custom permutation, without the tree code knowing which.
# ---------------------------------------------------------------------------

"""
    abstract type AbstractHealpixOrdering

How pixels of a dense HEALPix grid map to positions in a data vector — the
HEALPix branch of [`DiscreteGlobalGrids.AbstractFaceOrdering`](@ref), which
carries the full extension contract — `data_index(o, space, ix, iy, face)`,
`lattice_index(o, space, j)` and the optional `validate_ordering(o, space)`,
where `space` is a [`HealpixFaceSpace`](@ref), `ix`, `iy` are 0-based in
`0:space.nside-1` and `face` is 0-based in `0:11`.

A [`HealpixFaceGrid`](@ref) is a lattice (`12 nside²` pixels addressed by
`(ix, iy, face)`) plus one of these; the ordering is the *only* thing that
decides which column of a `ConservativeRegridding.Regridder` a pixel lands in.
[`RingOrder`](@ref) and [`NestedOrder`](@ref) are the two shipped instances.
"""
abstract type AbstractHealpixOrdering <: DGG.AbstractFaceOrdering end

"""
    RingOrder()

RING data ordering: pixels numbered north→south along iso-latitude rings and
west→east within a ring, which is how HEALPix FITS products and
`SpeedyWeather`/`RingGrids` fields are laid out.

Defined for **any** `nside >= 1` — the ring closed forms
([`xyf_to_ring`](@ref) / [`ring_to_xyf`](@ref)) carry no power-of-two
restriction.
"""
struct RingOrder <: AbstractHealpixOrdering end

data_index(::RingOrder, space::HealpixFaceSpace, ix::Integer, iy::Integer, face::Integer) =
    xyf_to_ring(ix, iy, face, space.nside)

lattice_index(::RingOrder, space::HealpixFaceSpace, j::Integer) = ring_to_xyf(j, space.nside)

"""
    NestedOrder()

NESTED data ordering: position `j` holds the pixel whose 0-based EOPF nested id
is `j - 1`. This is the id space `HealpixKernel.jl` wires
(`DGGSGrid(HEALPixDGGS(), level)` has ordinal `p + 1` for nested id `p`), so a
`HealpixFaceGrid` built with `NestedOrder` is column-for-column interchangeable
with the id-hierarchy tree at `nside = 2^level`.

Requires `nside = 2^k`: the nested id is a Morton code, which only exists on a
`2^k × 2^k` face. [`HealpixFaceGrid`](@ref) rejects any other `nside` at
construction rather than at first query.
"""
struct NestedOrder <: AbstractHealpixOrdering end

# `+ 1` because nested ids are 0-based (EOPF) while data positions are 1-based;
# see the index-convention block at the top of `chart.jl`.
data_index(::NestedOrder, space::HealpixFaceSpace, ix::Integer, iy::Integer, face::Integer) =
    Int(xyf_to_nested(ix, iy, face, space.nside)) + 1

lattice_index(::NestedOrder, space::HealpixFaceSpace, j::Integer) =
    nested_to_xyf(j - 1, space.nside)

# The one thing that must be caught eagerly: a nested-ordered grid at
# nside = 3 would construct fine and then throw from deep inside a dual-tree
# traversal, which is a miserable place to learn about it.
function validate_ordering(::NestedOrder, space::HealpixFaceSpace)
    ispow2(space.nside) || throw(ArgumentError(
        "NestedOrder requires nside = 2^k, got nside=$(space.nside); \
         use RingOrder() for arbitrary nside"))
    return nothing
end

# ---------------------------------------------------------------------------
# The `FaceGridSystem` contract
# ---------------------------------------------------------------------------

# 12 base faces. (Numerically `root_count(HEALPixDGGS())` too — a coincidence of
# this system, not a contract, so the literal stays independent.)
DGG.nfaces(::HealpixFaceSystem) = 12

DGG.face_chart(::HealpixFaceSystem, x, y, face::Integer) = xyf_to_point(x, y, face)

DGG.face_cell_corners(space::HealpixFaceSpace, ix::Integer, iy::Integer, face::Integer) =
    pixel_corners(ix, iy, face, space.nside)

# Upper bound: the ring/pixel-count arithmetic (`12 * nside^2`, `ring_first`'s
# `12nside^2 - 2js * (js + 1)` in `chart.jl`) overflows `Int64` past
# `nside = 2^29` and wraps silently rather than throwing. Mirrors
# `HEALPixDGGS`'s `max_level == 29` (`nside = 2^level`).
DGG.max_nside(::HealpixFaceSystem) = 2^29

DGG.default_ordering(::HealpixFaceSystem) = RingOrder()

DGG.ordering_family(::HealpixFaceSystem) = AbstractHealpixOrdering

DGG.facegrid_prefix(::HealpixFaceSystem) = "Healpix"

DGG.facegrid_cell_noun(::HealpixFaceSystem) = "pixels"

# A HEALPix pixel never bulges outside the cap of its 4 corners + great-circle edge
# midpoints, so we can skip the generic method's inclusion of all perimeter vertices.
# (Ported verbatim from ConservativeRegridding's RingGrids ext.)
DGG.cap_policy(::HealpixFaceSystem) = DGG.FourCornerCap()

# ---------------------------------------------------------------------------
# The HEALPix names for the shared types
# ---------------------------------------------------------------------------

"""
    HealpixFaceGrid(space::HealpixFaceSpace, ordering::AbstractHealpixOrdering)
    HealpixFaceGrid(nside::Integer; ordering = RingOrder())

A complete HEALPix grid at one resolution — all `12 nside²` pixels — together
with the data ordering its pixels are numbered by. `treeify(grid)` turns it into
a spatial tree ([`HealpixFaceRoot`](@ref)) that
`ConservativeRegridding.Regridder` consumes directly.

```julia
grid = HealpixFaceGrid(4; ordering = NestedOrder())
R = ConservativeRegridding.Regridder(treeify(grid), treeify(other))
```

# Alignment

**Column `j` of a `Regridder` built on `treeify(grid)` is data-vector position
`j` under `ordering` — there is no permutation between the matrix and the
field.** This is structural, not a convention to remember: the per-face cursors
emit `data_index(ordering, ...)` as their leaf indices and `Trees.getcell(root,
j)` returns the polygon of `lattice_index(ordering, space, j)`, so the matrix
is *assembled* in data order. Concretely, with `RingOrder` a column indexes a
ring-ordered HEALPix FITS field as-is, and with `NestedOrder` column `j` is
nested id `j - 1` — the same identity the id-hierarchy tree
(`DGGSGrid(HEALPixDGGS(), level)`) uses, so the two are directly comparable at
`nside = 2^level`.

# Validation

`nside >= 1` is checked by [`HealpixFaceSpace`](@ref); the ordering gets a say
too, through `validate_ordering` — which is why
`HealpixFaceGrid(3; ordering = NestedOrder())` throws an `ArgumentError` here
rather than failing later inside a traversal.
"""
const HealpixFaceGrid = DGG.FaceGrid{HealpixFaceSystem}

"""
    FaceChartGrid(manifold, space::HealpixFaceSpace, face, ordering)

One HEALPix face as an `nside × nside` `Trees.AbstractCurvilinearGrid`, so that
`Trees.TopDownQuadtreeCursor` can index it with no HEALPix-specific descent
logic. `face` is 0-based (`0:11`). Internal: build a
[`HealpixFaceGrid`](@ref) and `treeify` it.

Cartesian cell index `(i, j)` (1-based) is lattice pixel `(i-1, j-1)`, and
`Trees.getcell` returns its CCW polygon (see [`pixel_polygon`](@ref) — CCW is
required by the convex-clip kernel, which clips a CW ring to EMPTY).
`Trees.getvertex(g, i, j)` is the *lattice* point `((i-1)/nside,
(j-1)/nside)`, `1:(nside+1)` in each direction, which is what drives the
bounding caps.

The two index maps are overridden away from the generic column-major default:
they target the *global* data layout of `ordering`, not a face-local one, so
the leaf indices this grid's cursor reports are already data-vector positions.
"""
const FaceChartGrid = DGG.FaceChartGrid{HealpixFaceSystem}

"""
    HealpixFaceRoot(manifold, space::HealpixFaceSpace, ordering)
    HealpixFaceRoot(nside::Integer, ordering = RingOrder())

Root of the spatial tree over a dense HEALPix grid: the whole sphere, with the
twelve base faces as children. Each child is a
`Trees.TopDownQuadtreeCursor` over that face's [`FaceChartGrid`](@ref); leaf
indices throughout are data-vector positions under `ordering`, so
`Trees.getcell(root, j)` is the pixel of data position `j` and a `Regridder`'s
column `j` is that same position.

Build one with `treeify(::HealpixFaceGrid)` rather than by hand. Constructed
directly it re-runs the two ordering checks `HealpixFaceGrid` runs, so
`HealpixFaceRoot(3, NestedOrder())` throws here rather than from inside a
traversal.
"""
const HealpixFaceRoot = DGG.FaceGridRoot{HealpixFaceSystem}

# ---------------------------------------------------------------------------
# HEALPix vocabulary over the shared layer
#
# HEALPix says "pixel" where the shared layer says "cell"; these three wrappers
# are that translation and nothing else.
# ---------------------------------------------------------------------------

"""
    num_pixels(space::HealpixFaceSpace) -> Int

`12 * nside^2`, the number of pixels of the grid — and the length of the data
vector any [`AbstractHealpixOrdering`](@ref) indexes.
"""
num_pixels(space::HealpixFaceSpace) = DGG.num_cells(space)

"""
    num_pixels(grid::HealpixFaceGrid) -> Int

`12 * nside^2`, the number of pixels of the grid — the length of the data
vector its ordering indexes, and the number of columns of a `Regridder` built
on it.
"""
num_pixels(grid::HealpixFaceGrid) = DGG.num_cells(grid)

"""
    pixel_polygon(ix, iy, face, nside) -> GI.Polygon

Pixel `(ix, iy)` of `face` as a closed 4-gon on the unit sphere, from
[`pixel_corners`](@ref).

The ring is **counter-clockwise as seen from outside the sphere**, which
`pixel_corners` guarantees and which is a hard contract rather than a
convention: the convex-clip kernel that computes spherical intersections clips
a clockwise ring to EMPTY, so a reversed ring yields silent zero areas instead
of an error. Since every corner is a shared lattice point evaluated by the same
chart function, neighbouring polygons meet at bit-identical vertices and the
tessellation is exact.

`nside` is taken as a plain integer here — this is the chart-side vocabulary
wrapper, alongside [`pixel_corners`](@ref) — but it is wrapped in a
[`HealpixFaceSpace`](@ref) before the shared layer sees it, so an inadmissible
`nside` is an `ArgumentError` rather than a silently wrong polygon.
"""
pixel_polygon(ix::Integer, iy::Integer, face::Integer, nside::Integer) =
    DGG.face_cell_polygon(HealpixFaceSpace(nside), ix, iy, face)
