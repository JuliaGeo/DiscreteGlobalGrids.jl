# ---------------------------------------------------------------------------
# Face-chart grids — the chart layer as a spatial tree
#
# Several systems in this package factor the same way: `nfaces` continuous
# charts `[0, 1]² → S²`, one per base face, plus an `nside × nside` square
# lattice on each face whose shared points make the tessellation exact. This
# file turns that factorization into something a
# `ConservativeRegridding.Regridder` can consume: a *dense* grid — all
# `nfaces * nside²` cells of one resolution — as a `SpatialTreeInterface` tree,
# with no per-system tree code at all. A system joins by defining a
# [`FaceGridSystem`](@ref) singleton plus seven small methods; the three
# shipped instances are `src/HEALPix/face_grid.jl` (12 faces),
# `src/S2/face_grid.jl` (6) and `src/ISEA4R/face_grid.jl` (10 diamonds), and
# each one is that singleton, its orderings, and nothing else.
#
# ## Why a per-face curvilinear grid rather than an id hierarchy
#
# The package also has a tree over each system's *id* hierarchy (`DGGSGrid`
# → `DGGSCursor`, see `src/core/generic_cursor.jl`). Those are radix quadtrees:
# they exist only for `nside = 2^k` (or `3^k`, ...) and they descend by id
# arithmetic. The *geometry* has no such restriction: a face is an
# `nside × nside` block of the chart lattice for any `nside >= 1`, and
# subdividing an index *range* needs no power of two at all. So a face is
# exactly a `Trees.AbstractCurvilinearGrid`, and the stock
# `Trees.TopDownQuadtreeCursor` (range bisection, no Morton/Hilbert arithmetic
# anywhere) indexes it for free; the faces are tied together by an
# `nfaces`-child root node. One code path covers every `nside`, and the
# ordering — which *data position* each cell occupies — is a separate,
# swappable component (see [`AbstractFaceOrdering`](@ref)).
#
# Provenance: the shape of this layer is a port of ConservativeRegridding.jl's
# RingGrids extension (`ext/ConservativeRegriddingRingGridsExt/healpix.jl`, the
# `HEALPixFaceGrid` / `HEALPixRootNode` pair). The differences are deliberate
# and are the point of the port: the closed forms live in each system's
# `chart.jl` rather than inlined into the tree file, and the data ordering is a
# component instead of being hard-wired to RING. `STI.node_extent` for a block
# of cells and the N-child root are ported essentially verbatim.
#
# ## The alignment rule
#
# Column `j` of a `Regridder` built on `treeify(grid)` is data-vector position
# `j` under the grid's ordering — no permutation anywhere. It is structural:
# the per-face cursors emit `data_index(ordering, ...)` as their leaf indices
# and `Trees.getcell(root, j)` returns the polygon of
# `lattice_index(ordering, space, j)`, so the matrix is *assembled* in data
# order.
#
# ## Resolution travels as a `FaceGridSpace`, never as a loose `Int`
#
# [`FaceGridSpace`](@ref) has an inner constructor and no other way in, so a
# value of that type is a *checked* resolution — `1 <= nside <= max_nside`,
# always. That guarantee is worth nothing if the layer below it accepts a bare
# `nside::Int` alongside the system singleton, because then every internal
# caller (and every user who reaches for the tree types directly) can fabricate
# a resolution the charts cannot be evaluated at. So the split is by
# *resolution-dependence*, not by convenience:
#
#   * the resolution-FREE contract methods dispatch on the system singleton —
#     `nfaces`, `face_chart` (continuous `[0, 1]²`, deliberately grid-free so
#     the id-hierarchy kernels can evaluate it at `2^level` lattice points with
#     no grid in hand), `max_nside`, `default_ordering`, `ordering_family`,
#     `facegrid_prefix`, `cap_policy`, `facegrid_cell_noun`;
#   * the resolution-BOUND ones take the space, which *is* system + nside —
#     `face_cell_corners`, `face_cell_polygon`, `data_index`, `lattice_index`,
#     `validate_ordering`.
#
# [`FaceChartGrid`](@ref) and [`FaceGridRoot`](@ref) store the space for the
# same reason: there is then no field for an unchecked `nside` to live in, and
# the tree internals cannot construct one even by accident. Each system's shim
# unwraps `space.nside` at the boundary and hands it to its own raw-Int codecs
# (`pixel_corners`, `ring_to_xyf`, ...), which stay resolution-agnostic
# building blocks shared with the id-hierarchy kernels.
#
# ## CCW is a hard contract
#
# Every polygon this layer emits winds counter-clockwise as seen from outside
# the sphere. That is not a convention: the convex-clip kernel that computes
# spherical intersections clips a clockwise ring to EMPTY, so a reversed ring
# yields silent zero areas instead of an error. The obligation lands on each
# system's `face_cell_corners`, which must emit its four corners in CCW order.
#
# ## Cap soundness
#
# Node extents are what make tree pruning sound: an extent that does not
# contain its subtree's geometry does not error, it silently drops candidate
# pairs in the dual DFS and produces a `Regridder` that quietly loses area.
# The generic perimeter walk (`Trees.cell_range_extent`, which samples every
# perimeter lattice vertex) is sound for any chart and is therefore the
# default. A system may opt into the O(1) four-corner cap by overriding
# [`cap_policy`](@ref) to [`FourCornerCap`](@ref) — but only together with a
# written soundness argument, proved or measured, at the override. See the
# three shipped overrides for what that looks like: S2 proves it from geodesic
# block edges, ISEA4R measures it under a pre-registered decision rule with the
# measurement kept as a standing test, HEALPix ports a measured property of the
# HEALPix chart.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# The system contract
# ---------------------------------------------------------------------------

"""
    abstract type FaceGridSystem

A grid system that factors into `nfaces` charts `[0, 1]² → S²` over a square
lattice — the contract this file's grid, tree and cap machinery is written
against. Instances are per-system singletons living in the system's submodule
(`HEALPix.HealpixFaceSystem`, `S2.S2FaceSystem`, `ISEA4R.Isea4rFaceSystem`).

# Extension contract

A new system `Sys <: FaceGridSystem` implements six methods on the singleton,
which are the facts that do not depend on a resolution:

```julia
DiscreteGlobalGrids.nfaces(::Sys)                 -> Int
DiscreteGlobalGrids.face_chart(::Sys, x, y, face) -> UnitSphericalPoint
DiscreteGlobalGrids.max_nside(::Sys)              -> Int
DiscreteGlobalGrids.default_ordering(::Sys)       -> AbstractFaceOrdering
DiscreteGlobalGrids.ordering_family(::Sys)        -> Type
DiscreteGlobalGrids.facegrid_prefix(::Sys)        -> String
```

and one on the *space*, which is the one fact that does:

```julia
DiscreteGlobalGrids.face_cell_corners(::FaceGridSpace{Sys}, ix, iy, face) -> NTuple{4}
```

The space is system + a checked `nside` in one value, so dispatching the
resolution-bound half of the contract on it is what keeps an unvalidated
`nside` out of the layer entirely — see the file header. Systems write that
method as a one-line unwrap onto their own raw-Int corner function
(`pixel_corners(ix, iy, face, space.nside)`).

Two more may be overridden, both of which have safe defaults:

```julia
DiscreteGlobalGrids.cap_policy(::Sys)          # default PerimeterWalkCap()
DiscreteGlobalGrids.facegrid_cell_noun(::Sys)  # default "cells"
```

plus its own orderings (`<: AbstractFaceOrdering`, see there). Everything else
— [`FaceGridSpace`](@ref), [`FaceGrid`](@ref), [`FaceChartGrid`](@ref),
[`FaceGridRoot`](@ref), the cursors, the extents and the `treeify` wiring — is
shared and needs no per-system code.

# Why not the registry singletons

These are deliberately *not* the `AbstractDGGS` registry types (`HEALPixDGGS`,
`S2DGGS`, `ISEA4RDGGS`). Face-grid facts are chart-arithmetic facts, and they
disagree with the id-hierarchy facts of the twin: `max_nside(Isea4rFaceSystem())
== 2^29` while `max_level(ISEA4RDGGS()) === nothing`, because one bounds the
lattice arithmetic and the other bounds a refinement hierarchy that is not
ported. Keeping the parameter separate also keeps `FaceGrid{H3DGGS}` from being
type-checkable nonsense: only systems that actually factor this way have a
singleton here.
"""
abstract type FaceGridSystem end

# ---------------------------------------------------------------------------
# Resolution
#
# The other half of the contract's dispatch: the system singleton above carries
# the facts that hold at every resolution, `FaceGridSpace` carries the
# resolution itself. It is defined here, ahead of the methods that take it,
# because the resolution-bound half of the contract — `face_cell_corners`,
# `face_cell_polygon`, `data_index`, `lattice_index`, `validate_ordering` —
# spells it in their signatures.
# ---------------------------------------------------------------------------

"""
    FaceGridSpace{Sys}(nside)

The resolution of a dense face grid: `nside` cells along each edge of each of
the system's `nfaces` base faces, hence `nfaces * nside^2` cells in total.

Any `1 <= nside <= max_nside(Sys())` is admissible — this is the *chart*
resolution, not a refinement level. Power-of-two (or power-of-three, ...)
restrictions belong to the individual orderings, and are enforced by
[`validate_ordering`](@ref) when (and only when) a grid is built with one:
`HealpixFaceSpace(3)` and `HealpixFaceSpace(5)` are perfectly good HEALPix
grids that simply have no nested id space.

Per-system aliases: `HealpixFaceSpace`, `S2FaceSpace`, `Isea4rFaceSpace`.

```julia
space = HealpixFaceSpace(4)     # 192 pixels
space.nside                     # 4
```
"""
struct FaceGridSpace{S<:FaceGridSystem}
    nside::Int
    # Inner constructor so no path — including the tree internals — can
    # fabricate an `nside` the chart maps cannot be evaluated at.
    function FaceGridSpace{S}(nside::Integer) where {S<:FaceGridSystem}
        n = Int(nside)
        n >= 1 || throw(ArgumentError("nside must be >= 1, got $n"))
        # Upper bound: past `max_nside` the system's own cell-count / index
        # arithmetic overflows `Int64` and wraps silently rather than throwing.
        n <= max_nside(S()) ||
            throw(ArgumentError("nside must be <= $(_nside_bound_string(max_nside(S()))), got $n"))
        return new{S}(n)
    end
end

# The bounds are powers of two and are documented as such (`2^29`, `2^30`);
# printing `536870912` in the error message would make them unrecognizable.
_nside_bound_string(n::Integer) = ispow2(n) ? "2^$(trailing_zeros(n))" : string(n)

"""
    num_cells(space::FaceGridSpace) -> Int

`nfaces * nside^2`, the number of cells of the grid — and the length of the
data vector any [`AbstractFaceOrdering`](@ref) indexes. HEALPix says
`num_pixels` for the same quantity.
"""
num_cells(space::FaceGridSpace{S}) where {S} = nfaces(S()) * space.nside^2

Base.show(io::IO, space::FaceGridSpace{S}) where {S} =
    print(io, facegrid_prefix(S()), "FaceSpace(", space.nside, ")")

# ---------------------------------------------------------------------------
# The contract methods
# ---------------------------------------------------------------------------

"""
    nfaces(sys::FaceGridSystem) -> Int

Number of base faces the system's charts cover — 12 for HEALPix, 6 for S2, 10
for ISEA4R. The root node of the tree has exactly this many children, and the
grid has `nfaces(sys) * nside^2` cells.

Numerically this equals `root_count` of the system's registry twin in all three
shipped instances; that is a coincidence of these systems, not a contract, and
the two are kept independent.
"""
function nfaces end

"""
    face_chart(sys::FaceGridSystem, x, y, face) -> UnitSphericalPoint

The system's chart: continuous `[0, 1]² → S²` for each 0-based `face`,
evaluated at *continuous* coordinates `(x, y)` (not lattice indices), so that
any `nside` and any refinement depth can be addressed. This is the single
geometry function of the contract — vertices, and through
[`face_cell_corners`](@ref) cell polygons, are all chart evaluations at lattice
points, which is what makes neighbouring cells share bit-identical corners and
the tessellation exact rather than consistent to rounding.

It dispatches on the singleton rather than on a [`FaceGridSpace`](@ref)
precisely because it carries no resolution: the id-hierarchy kernels
(`HealpixKernel.jl`, `S2Kernel.jl`, `Isea4rKernel.jl`) evaluate it at `2^level`
lattice points with no grid in hand at all, and a space parameter would make
them fabricate one to ask a question that does not depend on it.
"""
function face_chart end

"""
    face_cell_corners(space::FaceGridSpace{Sys}, ix, iy, face) -> NTuple{4}

The four corners of lattice cell `(ix, iy)` of `face` at the resolution of
`space`, **counter-clockwise as seen from outside the sphere**. `ix`, `iy` are
0-based in `0:space.nside-1` and `face` is 0-based in `0:nfaces-1`.

The resolution arrives as the space rather than as a trailing `nside::Int`
because the space is the only checked carrier of one (see the file header); the
system is recovered from its type parameter, so this is still one method per
system.

Systems delegate it to their own corner function (`HEALPix.pixel_corners`,
`S2.cell_corners`, `ISEA4R.cell_corners`), unwrapping `space.nside` at that
boundary, rather than have it re-derived here: those functions already own the
system's seam-ownership rules, which is what keeps shared lattice points
bit-identical between adjacent cells. Only the polygon assembly is hoisted,
into [`face_cell_polygon`](@ref).
"""
function face_cell_corners end

"""
    max_nside(sys::FaceGridSystem) -> Int

Largest admissible `nside` for the system's chart lattice, enforced by the
[`FaceGridSpace`](@ref) constructor. This is an arithmetic bound, not a
refinement-level bound: past it the system's own cell-count or index arithmetic
overflows `Int64` and wraps silently rather than throwing.
"""
function max_nside end

"""
    default_ordering(sys::FaceGridSystem) -> AbstractFaceOrdering

The ordering `FaceGrid{Sys}(nside)` uses when none is named. Always the one
defined for *every* `nside >= 1` (`RingOrder` for HEALPix, `RowMajorOrder` for
S2 and ISEA4R), never a power-of-two-only ordering.
"""
function default_ordering end

"""
    ordering_family(sys::FaceGridSystem) -> Type

The abstract ordering type of this system (e.g. `HEALPix.AbstractHealpixOrdering`).
[`FaceGrid`](@ref)'s constructor checks membership against it, so handing a
system one of *another* system's orderings is an `ArgumentError` at
construction rather than a `MethodError` from inside a traversal.
"""
function ordering_family end

"""
    facegrid_prefix(sys::FaceGridSystem) -> String

The system's name as it appears in the type aliases and `show` output —
`"Healpix"`, `"S2"`, `"Isea4r"` — spelled so that `prefix * "FaceGrid"` is the
alias a user actually types.
"""
function facegrid_prefix end

"""
    facegrid_cell_noun(sys::FaceGridSystem) -> String

What this system calls its cells in `show` output. Defaults to `"cells"`;
HEALPix overrides it to `"pixels"`.
"""
facegrid_cell_noun(::FaceGridSystem) = "cells"

"""
    face_cell_polygon(space::FaceGridSpace, ix, iy, face) -> GI.Polygon

Lattice cell `(ix, iy)` of `face` at the resolution of `space`, as a closed
4-gon on the unit sphere, from [`face_cell_corners`](@ref).

The ring is **counter-clockwise as seen from outside the sphere**, which
`face_cell_corners` guarantees and which is a hard contract rather than a
convention: the convex-clip kernel that computes spherical intersections clips
a clockwise ring to EMPTY, so a reversed ring yields silent zero areas instead
of an error. Since every corner is a shared lattice point evaluated by the same
chart function, neighbouring polygons meet at bit-identical vertices and the
tessellation is exact.

Whether the 4-gon *is* the cell or only approximates it is a per-system fact of
the chart (exact for S2, whose cell edges are geodesics; approximate for
HEALPix and ISEA4R, whose edges are not) and is documented on each system's
wrapper — `HEALPix.pixel_polygon`, `S2.cell_polygon`, `ISEA4R.cell_polygon`.
"""
function face_cell_polygon(space::FaceGridSpace, ix::Integer, iy::Integer,
                           face::Integer)
    a, b, c, d = face_cell_corners(space, ix, iy, face)
    return GI.Polygon([GI.LinearRing([a, b, c, d, a])])
end

# ---------------------------------------------------------------------------
# Orderings
#
# The lattice `(ix, iy, face)` says *where* a cell is; an ordering says *which
# slot of the data vector* it occupies. Keeping the two apart is what lets one
# grid type serve every on-disk layout a system has — ring-ordered files,
# nested/Hilbert/Morton-ordered files, and (later) arbitrary permutations —
# without the tree code knowing which.
#
# Orderings stay per-system: `S2.RowMajorOrder` and `ISEA4R.RowMajorOrder` are
# deliberately distinct types with per-system range checks, so the three
# contract functions below are shared but no concrete ordering is.
# ---------------------------------------------------------------------------

"""
    abstract type AbstractFaceOrdering

How cells of a dense face grid map to positions in a data vector.

A [`FaceGrid`](@ref) is a lattice (`nfaces * nside²` cells addressed by
`(ix, iy, face)`) plus one of these; the ordering is the *only* thing that
decides which column of a `ConservativeRegridding.Regridder` a cell lands in.
Each system subtypes this once (`HEALPix.AbstractHealpixOrdering`,
`S2.AbstractS2Ordering`, `ISEA4R.AbstractIsea4rOrdering`) and ships its
concrete orderings under that.

# Extension contract

A new ordering `O` implements two methods, which must be exact mutual inverses
over `1:num_cells(space)` for every space the ordering admits:

```julia
DiscreteGlobalGrids.data_index(::O, space, ix, iy, face) -> Int             # 1-based position
DiscreteGlobalGrids.lattice_index(::O, space, j)         -> (ix, iy, face)  # its inverse
```

and may implement a third, which is called once at grid construction and is the
place to reject resolutions the ordering cannot index (this is exactly how
`NestedOrder` / `HilbertOrder` / `MortonOrder` refuse a non-power-of-two
`nside`):

```julia
DiscreteGlobalGrids.validate_ordering(::O, space) -> nothing (or throw)
```

`space` is the [`FaceGridSpace`](@ref) the grid was built on — system and a
checked `nside` in one value, which is how a resolution reaches this layer at
all (see the file header). Orderings are per-system types, so the space's
system parameter is redundant with the ordering's own dispatch and each shim
simply reads `space.nside` and calls its system's raw-Int codec.

`ix`, `iy` are 0-based in `0:space.nside-1` and `face` is 0-based in
`0:nfaces-1`, matching each system's `chart.jl` throughout. Nothing else about
the tree changes: `treeify`, the per-face cursors and the node extents are all
written against these two functions, so a permutation ordering (a space-filling
curve, or the column order of an on-disk product) can be added later without
touching [`FaceGrid`](@ref), [`FaceChartGrid`](@ref) or
[`FaceGridRoot`](@ref).

Per-system aliases: `HealpixFaceGrid`, `S2FaceGrid`, `Isea4rFaceGrid`.
"""
abstract type AbstractFaceOrdering end

"""
    data_index(ordering, space::FaceGridSpace, ix, iy, face) -> Int

1-based position of lattice cell `(ix, iy)` of `face` — at the resolution of
`space` — in a data vector laid out according to `ordering`. Half of the
[`AbstractFaceOrdering`](@ref) contract; [`lattice_index`](@ref) is the other
half and its exact inverse.

There is deliberately **no generic fallback method here** — not even for
row-major. `S2.RowMajorOrder` and `ISEA4R.RowMajorOrder` are distinct types
whose codecs carry per-system range checks and face counts, and collapsing them
would silently accept out-of-range faces. Every ordering brings its own pair.
"""
function data_index end

"""
    lattice_index(ordering, space::FaceGridSpace, j) -> (ix, iy, face)

Lattice coordinates (0-based `ix`, `iy`, 0-based `face`) of the cell at 1-based
data position `j` under `ordering`, at the resolution of `space` — the inverse
of [`data_index`](@ref).
"""
function lattice_index end

"""
    validate_ordering(ordering, space::FaceGridSpace) -> nothing

Called once by the [`FaceGrid`](@ref) constructor, and again by
[`FaceGridRoot`](@ref)'s (which is independently constructible). Throw here if
`ordering` cannot index a grid of this resolution; the default accepts
everything.

The point is *when* it runs: a Morton-style ordering at `nside = 3` would
construct fine and then throw from deep inside a dual-tree traversal, which is
a miserable place to learn about it.
"""
validate_ordering(::AbstractFaceOrdering, ::FaceGridSpace) = nothing

# ---------------------------------------------------------------------------
# Cap policy
#
# The per-system soundness obligation of section "Cap soundness" above, made
# into a dispatch. `PerimeterWalkCap` is the default because it is sound for
# any chart; `FourCornerCap` is opt-in and must be accompanied by an argument.
# ---------------------------------------------------------------------------

"""
    abstract type CapPolicy

How a system's block (internal-node) extents are computed:
[`PerimeterWalkCap`](@ref), the sound default, or [`FourCornerCap`](@ref), the
O(1) shortcut. Selected per system by [`cap_policy`](@ref).
"""
abstract type CapPolicy end

"""
    FourCornerCap()

Bound an index block by `Trees.circle_from_four_corners` on its four *chart*
corners, skipping the generic method's inclusion of all perimeter vertices —
O(1) per node instead of O(perimeter).

**This is only sound if the block's geometry never bulges outside the cap of
its four corners plus the great-circle midpoints of the corner-to-corner
edges** (the `1.0001` slack `circle_from_four_corners` applies absorbs the
residual). That is a property of the specific chart, not of the face-grid
pattern, so a system selecting this policy must carry the argument — proof or
pre-registered measurement — at its [`cap_policy`](@ref) override. An unsound
cap does not error: it silently drops candidate pairs in the dual DFS and
produces a `Regridder` that quietly loses area.
"""
struct FourCornerCap <: CapPolicy end

"""
    PerimeterWalkCap()

Bound an index block with the stock `Trees.cell_range_extent` — the generic
spherical perimeter walk over every boundary lattice vertex of the block. Sound
for any chart (it needs only the much weaker property that a *single cell edge*
stays near the cap of its endpoints) and merely slower than
[`FourCornerCap`](@ref). This is the default, and it is the fallback any system
can return to by deleting its [`cap_policy`](@ref) override.
"""
struct PerimeterWalkCap <: CapPolicy end

"""
    cap_policy(sys::FaceGridSystem) -> CapPolicy

Which block-extent strategy this system's tree uses. Defaults to
[`PerimeterWalkCap`](@ref), which is sound for any chart; override to
[`FourCornerCap`](@ref) *only* with a written soundness argument at the
override site.
"""
cap_policy(::FaceGridSystem) = PerimeterWalkCap()

# ---------------------------------------------------------------------------
# User-facing grid
# ---------------------------------------------------------------------------

"""
    FaceGrid{Sys}(space::FaceGridSpace{Sys}, ordering::AbstractFaceOrdering)
    FaceGrid{Sys}(nside::Integer; ordering = default_ordering(Sys()))

A complete face grid at one resolution — all `nfaces * nside²` cells — together
with the data ordering its cells are numbered by. `treeify(grid)` turns it into
a spatial tree ([`FaceGridRoot`](@ref)) that
`ConservativeRegridding.Regridder` consumes directly.

Per-system aliases: `HealpixFaceGrid`, `S2FaceGrid`, `Isea4rFaceGrid`.

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
is *assembled* in data order. What data position `j` means concretely is the
ordering's business — see each system's ordering docstrings.

# Validation

`nside` is checked by [`FaceGridSpace`](@ref); the ordering gets two says of
its own — it must belong to this system's [`ordering_family`](@ref) (handing a
grid another system's ordering is an `ArgumentError`, not a `MethodError` from
inside a traversal), and it may reject the resolution through
[`validate_ordering`](@ref), which is why
`HealpixFaceGrid(3; ordering = NestedOrder())` throws here rather than failing
later during a traversal.

The space is built *first*, so `HealpixFaceGrid(0; ordering = NestedOrder())`
reports the `nside` problem rather than the ordering one.
"""
struct FaceGrid{S<:FaceGridSystem,O<:AbstractFaceOrdering}
    space::FaceGridSpace{S}
    ordering::O
    function FaceGrid{S,O}(space::FaceGridSpace{S}, ordering::O) where {S<:FaceGridSystem,O<:AbstractFaceOrdering}
        ordering isa ordering_family(S()) || throw(ArgumentError(
            "$(typeof(ordering)) is not an ordering of this system; \
             expected a subtype of $(ordering_family(S()))"))
        validate_ordering(ordering, space)
        return new{S,O}(space, ordering)
    end
end

(::Type{FaceGrid{S}})(space::FaceGridSpace{S}, ordering::O) where {S,O<:AbstractFaceOrdering} =
    FaceGrid{S,O}(space, ordering)

# The space is constructed before the ordering is looked at, so an invalid
# `nside` reports as an `nside` problem whatever the ordering is.
(::Type{FaceGrid{S}})(nside::Integer; ordering::AbstractFaceOrdering=default_ordering(S())) where {S} =
    FaceGrid{S}(FaceGridSpace{S}(nside), ordering)

"""
    num_cells(grid::FaceGrid) -> Int

`nfaces * nside^2`, the number of cells of the grid — and the length of the
data vector its ordering indexes.
"""
num_cells(grid::FaceGrid) = num_cells(grid.space)

Base.show(io::IO, grid::FaceGrid{S}) where {S} =
    print(io, facegrid_prefix(S()), "FaceGrid(nside=", grid.space.nside, ", ",
        grid.ordering, ", ", num_cells(grid), " ", facegrid_cell_noun(S()), ")")

# ---------------------------------------------------------------------------
# Per-face curvilinear grid
#
# A single face is an `nside × nside` block of cells on the chart lattice —
# exactly an `AbstractCurvilinearGrid`, so it plugs into the stock
# `TopDownQuadtreeCursor` for free. Cell/vertex geometry is analytic (no stored
# matrix), and the index maps are overridden so the cursor emits *global data*
# positions: `cartesian_to_linear_idx` returns the ordering's data index of a
# cell and `linear_to_cartesian_idx` inverts it.
# ---------------------------------------------------------------------------

"""
    FaceChartGrid{Sys}(manifold, space::FaceGridSpace{Sys}, face, ordering)

One face as an `nside × nside` `Trees.AbstractCurvilinearGrid`, so that
`Trees.TopDownQuadtreeCursor` can index it with no system-specific descent
logic. `face` is 0-based (`0:nfaces-1`). Internal: build a [`FaceGrid`](@ref)
and `treeify` it.

Cartesian cell index `(i, j)` (1-based) is lattice cell `(i-1, j-1)`, and
`Trees.getcell` returns its CCW polygon (see [`face_cell_polygon`](@ref) — CCW
is required by the convex-clip kernel, which clips a CW ring to EMPTY).
`Trees.getvertex(g, i, j)` is the *lattice* point `((i-1)/nside, (j-1)/nside)`,
`1:(nside+1)` in each direction, which is what drives the bounding caps.

The two index maps are overridden away from the generic column-major default:
they target the *global* data layout of `ordering`, not a face-local one, so
the leaf indices this grid's cursor reports are already data-vector positions.

The resolution is held as a [`FaceGridSpace`](@ref), not as a bare `nside`, and
that is the *only* validation story this type has: it deliberately checks
nothing itself, because it is internal and every field it receives arrives
already checked from [`FaceGridRoot`](@ref) — which is in turn built either by
`treeify(::FaceGrid)` or by its own validating constructor. Holding the space
rather than an `Int` is what makes "already checked" a property of the type
system instead of a convention the root has to keep.

Per-system aliases: `HEALPix.FaceChartGrid`, `S2.FaceChartGrid`,
`ISEA4R.DiamondChartGrid`.
"""
struct FaceChartGrid{S<:FaceGridSystem,M<:GOCore.Manifold,O<:AbstractFaceOrdering} <: Trees.AbstractCurvilinearGrid{M}
    manifold::M
    space::FaceGridSpace{S}
    face::Int                                   # 0-based
    ordering::O
end

(::Type{FaceChartGrid{S}})(manifold::M, space::FaceGridSpace{S}, face::Integer, ordering::O) where {S,M<:GOCore.Manifold,O<:AbstractFaceOrdering} =
    FaceChartGrid{S,M,O}(manifold, space, Int(face), ordering)

GOCore.manifold(g::FaceChartGrid) = g.manifold

Trees.ncells(g::FaceChartGrid, ::Int) = g.space.nside

Trees.getcell(g::FaceChartGrid, i::Int, j::Int) =
    face_cell_polygon(g.space, i - 1, j - 1, g.face)

# Lattice vertex (1-based point index 1:(nside+1)) → continuous coord (i-1)/nside.
# Drives the per-cell leaf caps (`cell_range_extent`) and the block cap
# (`node_extent` below). `face_chart` is the one contract method that takes no
# resolution, so the division by `nside` happens here rather than inside it.
Trees.getvertex(g::FaceChartGrid{S}, i::Int, j::Int) where {S} =
    face_chart(S(), (i - 1) / g.space.nside, (j - 1) / g.space.nside, g.face)

Trees.cartesian_to_linear_idx(g::FaceChartGrid, idx::CartesianIndex{2}) =
    data_index(g.ordering, g.space, idx[1] - 1, idx[2] - 1, g.face)

# The face is fixed per grid, so the inverse drops the face the ordering
# reports — as the reference extension does. Any data position handed to a face
# grid comes from that face's own cursor, so the dropped value always agrees
# with `g.face`.
function Trees.linear_to_cartesian_idx(g::FaceChartGrid, idx::Integer)
    ix, iy, _face = lattice_index(g.ordering, g.space, Int(idx))
    return CartesianIndex(ix + 1, iy + 1)
end

# Block extents, per the system's `cap_policy`. See the cap-soundness section
# of this file's header, and `FourCornerCap` / `PerimeterWalkCap`.
STI.node_extent(q::Trees.TopDownQuadtreeCursor{<:FaceChartGrid{S}}) where {S<:FaceGridSystem} =
    _node_extent(cap_policy(S()), q)

function _node_extent(::FourCornerCap, q::Trees.TopDownQuadtreeCursor)
    g = q.grid
    imin, imax = extrema(q.leafranges[1]); imax += 1
    jmin, jmax = extrema(q.leafranges[2]); jmax += 1
    bl = Trees.getvertex(g, imin, jmin)
    tl = Trees.getvertex(g, imin, jmax)
    br = Trees.getvertex(g, imax, jmin)
    tr = Trees.getvertex(g, imax, jmax)
    return Trees.circle_from_four_corners((bl, tl, br, tr), ())
end

# Exactly the stock `STI.node_extent(::TopDownQuadtreeCursor)` body: the O(1)
# override above is opt-in, so the generic path has to be spelled out here for
# the systems that do not take it.
_node_extent(::PerimeterWalkCap, q::Trees.TopDownQuadtreeCursor) =
    Trees.cell_range_extent(q.grid, q.leafranges[1], q.leafranges[2])

# Caching child extents costs a vector per visited node, which only the
# perimeter walk's O(block edge) charting earns back; `FourCornerCap` is the
# four-vertex shortcut it exists to avoid. Every system here wires the shortcut,
# so this is `false` throughout the package and `true` for a system that keeps
# the generic path.
STI.node_extent_is_expensive(::Type{<:Trees.TopDownQuadtreeCursor{<:FaceChartGrid{S}}}) where {S<:FaceGridSystem} =
    cap_policy(S()) isa PerimeterWalkCap

# ---------------------------------------------------------------------------
# Toplevel tree
# ---------------------------------------------------------------------------

"""
    FaceGridRoot{Sys}(manifold, space::FaceGridSpace{Sys}, ordering)
    FaceGridRoot{Sys}(nside::Integer, ordering = default_ordering(Sys()))

Root of the spatial tree over a dense face grid: the whole sphere, with the
system's base faces as children. Each child is a
`Trees.TopDownQuadtreeCursor` over that face's [`FaceChartGrid`](@ref); leaf
indices throughout are data-vector positions under `ordering`, so
`Trees.getcell(root, j)` is the cell of data position `j` and a `Regridder`'s
column `j` is that same position.

Build one with `treeify(::FaceGrid)` rather than by hand; the two-argument
form exists to default the manifold (and to build the space) for REPL and test
construction.

# Validation

The resolution is a [`FaceGridSpace`](@ref), so it is checked by construction
whichever way in you take: `treeify` hands over the grid's own space, and the
`nside` form builds one, which is where a `0` or an over-`max_nside` resolution
is rejected.

The two *ordering* checks — membership in this system's
[`ordering_family`](@ref), then [`validate_ordering`](@ref) against the space —
are deliberately run here as well as in [`FaceGrid`](@ref)'s constructor. That
redundancy is the point: this type is independently constructible, so without
it `FaceGridRoot{HealpixFaceSystem}(3, NestedOrder())` (or a root carrying
*another* system's ordering) would build a tree that only fails deep inside a
dual-tree traversal, going around the checks `FaceGrid` exists to run. Both are
two comparisons on a singleton, so paying for them twice costs nothing that
matters and closes the direct-construction hole.

Per-system aliases: `HealpixFaceRoot`, `S2FaceRoot`, `Isea4rFaceRoot`.
"""
struct FaceGridRoot{S<:FaceGridSystem,M<:GOCore.Manifold,O<:AbstractFaceOrdering}
    manifold::M
    space::FaceGridSpace{S}
    ordering::O
    function FaceGridRoot{S,M,O}(manifold::M, space::FaceGridSpace{S}, ordering::O) where {S<:FaceGridSystem,M<:GOCore.Manifold,O<:AbstractFaceOrdering}
        ordering isa ordering_family(S()) || throw(ArgumentError(
            "$(typeof(ordering)) is not an ordering of this system; \
             expected a subtype of $(ordering_family(S()))"))
        validate_ordering(ordering, space)
        return new{S,M,O}(manifold, space, ordering)
    end
end

(::Type{FaceGridRoot{S}})(manifold::M, space::FaceGridSpace{S}, ordering::O) where {S,M<:GOCore.Manifold,O<:AbstractFaceOrdering} =
    FaceGridRoot{S,M,O}(manifold, space, ordering)

# These charts are inherently spherical; default the manifold for REPL/test
# construction. The space is built *first*, so an invalid `nside` reports as an
# `nside` problem whatever the ordering is — same rule as `FaceGrid`.
(::Type{FaceGridRoot{S}})(nside::Integer, ordering::AbstractFaceOrdering=default_ordering(S())) where {S} =
    FaceGridRoot{S}(GO.Spherical(), FaceGridSpace{S}(nside), ordering)

Base.show(io::IO, root::FaceGridRoot{S}) where {S} =
    print(io, facegrid_prefix(S()), "FaceRoot(nside=", root.space.nside, ", ", root.ordering, ")")

STI.isspatialtree(::Type{<:FaceGridRoot}) = true
STI.isleaf(::FaceGridRoot) = false
STI.nchild(::FaceGridRoot{S}) where {S} = nfaces(S())

STI.getchild(root::FaceGridRoot{S,M,O}, i::Int) where {S,M,O} =
    Trees.TopDownQuadtreeCursor(
        FaceChartGrid{S,M,O}(root.manifold, root.space, i - 1, root.ordering))

STI.getchild(root::FaceGridRoot{S}) where {S} =
    (STI.getchild(root, i) for i in 1:nfaces(S()))

# The faces tile the sphere, so the root bounds nothing tighter than all of it.
# `full_sphere_extent` is the package's own (`src/core/interface.jl`), identical
# to the reference extension's `SphericalCap(+z, nextfloat(π))`.
STI.node_extent(::FaceGridRoot) = full_sphere_extent()

Trees.ncells(root::FaceGridRoot) = num_cells(root.space)

"""
    Trees.getcell(root::FaceGridRoot, j::Int) -> GI.Polygon

Polygon of the cell at 1-based data position `j` under the root's ordering —
the geometry side of the alignment rule documented on [`FaceGrid`](@ref). This
is what `ConservativeRegridding` calls to compute both cell areas and
intersection areas, which is why "column `j` is data position `j`" needs no
permutation anywhere.

The ring is **counter-clockwise as seen from outside the sphere**, the same
hard contract every other polygon emitter here carries: the convex-clip kernel
that computes spherical intersections clips a clockwise ring to EMPTY, so a
reversed ring yields silent zero areas instead of an error.
"""
function Trees.getcell(root::FaceGridRoot, j::Int)
    count = Trees.ncells(root)
    1 <= j <= count || throw(BoundsError(1:count, j))
    ix, iy, face = lattice_index(root.ordering, root.space, j)
    return face_cell_polygon(root.space, ix, iy, face)
end

Trees.getcell(root::FaceGridRoot) =
    (Trees.getcell(root, j) for j in 1:Trees.ncells(root))

# ---------------------------------------------------------------------------
# Introspection
# ---------------------------------------------------------------------------

"""
    facesystem(x) -> FaceGridSystem

The system singleton behind a [`FaceGridSpace`](@ref), [`FaceGrid`](@ref),
[`FaceChartGrid`](@ref) or [`FaceGridRoot`](@ref) — the value every contract
method dispatches on, recovered from the type parameter.
"""
facesystem(::FaceGridSpace{S}) where {S} = S()
facesystem(::FaceGrid{S}) where {S} = S()
facesystem(::FaceChartGrid{S}) where {S} = S()
facesystem(::FaceGridRoot{S}) where {S} = S()

# --------------------------------------------------------------------------
# ConservativeRegridding.Trees wiring
#
# `best_manifold` is what makes the one-argument `treeify(grid)` work — it is
# `Trees`' own generic `treeify(grid) = treeify(best_manifold(grid), grid)`
# resolving through these methods, exactly as `src/core/generic_cursor.jl` does
# for the DGGS grids. Face-chart geometry is on the unit sphere by
# construction, so naming the manifold at a call site carries no information.
# --------------------------------------------------------------------------

GOCore.best_manifold(::FaceGrid) = GO.Spherical()
GOCore.best_manifold(root::FaceGridRoot) = root.manifold

Trees.treeify(m::GO.Spherical, grid::FaceGrid{S}) where {S} =
    FaceGridRoot{S}(m, grid.space, grid.ordering)

Trees.treeify(::GO.Spherical, root::FaceGridRoot) = root
