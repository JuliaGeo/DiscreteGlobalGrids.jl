# ---------------------------------------------------------------------------
# The ISEA4R system on the grid interface
#
# ISEA4R is the second DENSE-ORDINAL system in this package, and it is the same
# arithmetic as nested HEALPix with ten diamonds in place of twelve faces: at
# level `l` the cells are exactly `0:10*4^l - 1`, cell `p`'s children are
# `4p:4p+3`, and its parent is `p ÷ 4`. Every hierarchy method below is that
# arithmetic and nothing here needs a lookup table.
#
# ## The id, and the base of its numbering
#
# The canonical id is the interface's own `LevelIndex(level, index)`, where
#
#     index  ==  diamond * 4^level + morton(ix, iy),  0-BASED
#     position in `levelgrid(sys, level)`  ==  index + 1
#
# 0-based because the radix-4 prefix arithmetic that is the whole point of the
# Morton code — `4p:4p+3` for the children, `p ÷ 4` for the parent, `p >> 2Δ`
# for an ancestor — is child/parent arithmetic only in 0-based numbering. The
# `+ 1` therefore lives in exactly one place, the position/index conversion in
# `cellindex`/`cellposition`, where the interface's "a bare `Int` is a position"
# rule makes it visible rather than ambient.
#
# ## Why Morton and not row-major
#
# `chart.jl` carries both codecs, and the row-major one is valid at every
# `nside` while Morton needs `nside = 2^k`. Morton is nonetheless the canonical
# one, and the reason is the trait: `diamond * 4^level + morton` is depth-first
# curve order, so subtrees are contiguous runs of positions and
# `has_sorted_subtrees` is `true`. Row-major order scatters a cell's four
# children across two rows, and every prefix fast path in the package would be
# unavailable. Row-major is not offered as an alternate id scheme (see the
# deferral note at `cellindextypes`).
#
# ## What is a fast path here
#
#   * `cellat` is closed-form (`point_to_morton`) — no tree descent.
#   * `descendant_range` is `[p*4^Δ, (p+1)*4^Δ)` shifted into position space.
#   * `node_extent` is the EXACT subtree cap, not the inflated default.
#   * `neighbors`/`ring` walk the lattice and the seam tables, not the geometry.
#   * `ancestor` drops `2Δ` bits in one shift.
#   * `cell_area` is the closed-form equal-area value.
# ---------------------------------------------------------------------------

# ===========================================================================
# Types
# ===========================================================================

"""
    ISEA4RSystem() <: AbstractHierarchicalGridSystem

The ISEA4R hierarchical grid system: ten equal-area rhombus charts over the
icosahedron's ten diamonds, each refined by aperture-4 quadrant subdivision, on
the unit sphere.

Level `l` has `10 * 4^l` cells of exactly equal solid angle `4π / (10 * 4^l)`
— the Snyder ISEA projection is equal-area per face and the chart's two affine
halves each carry exactly `2π/5`, so this is exact rather than nominal. Cells
are quadrilaterals in the diamond chart; their edges follow chart lines, not
great circles, so [`cell_boundary`](@ref) densifies them — see that method.

# Ids

The canonical cell index is [`LevelIndex`](@ref), whose `index` field holds
`diamond * 4^level + morton(ix, iy)`, **0-based**; a cell's position in
`levelgrid(sys, l)` is `index + 1`. The ten-diamond layout those ids are written
against is this package's own convention with no external oracle behind it —
read the [`ISEA4R`](@ref) module docstring before inferring DGGAL / SST or any
other identifier compatibility.

# Levels

`0:29`. The bound is the `Int64` codec's: at level 29 there are
`10 * 4^29 = 2882303761517117440` cells, and level 30 would overflow a signed
64-bit cell count.

# Traits

`has_sorted_subtrees` is `true` — the Morton ordering *is* depth-first curve
order, so a subtree occupies a contiguous run of positions and
[`descendant_range`](@ref) is exact and hole-free.
[`max_neighbors`](@ref) is **9** under `Vertex()` and 4 under `Edge()`; the 9 is
not a typo and not slack — see [`neighbors`](@ref) and `topology.jl` on the
corner cells at icosahedron vertices 0 and 11.
[`node_extent`](@ref) is overridden with the cell's own bounding cap: children
tile their parent's chart rectangle exactly, so nothing needs inflating and
[`cap_inflation`](@ref) is never consulted.
"""
struct ISEA4RSystem <: DGG.AbstractHierarchicalGridSystem end

"""
    ISEA4RGrid(level) <: AbstractGrid

The complete ISEA4R grid at refinement `level`: all `10 * 4^level` cells in
Morton order. Built by [`levelgrid`](@ref); a lightweight descriptor, not a
materialised cell list.
"""
struct ISEA4RGrid <: DGG.AbstractGrid
    level::Int
end

# The one place `nside` is derived from a level, and the guard that keeps a bad
# level from silently producing a shift of 64.
@inline _nside(level::Integer) = Int64(1) << Int(level)
@inline _ncells(level::Integer) = 10 * (Int64(1) << (2 * Int(level)))

"The deepest level at which `10 * 4^level` still fits a signed 64-bit integer."
const MAX_LEVEL = 29

"The number of diamonds — the ten rhombi of `diamonds.jl`, and the root count."
const NDIAMONDS = 10

# ===========================================================================
# System interface
# ===========================================================================

DGG.cellindextype(::ISEA4RSystem) = DGG.LevelIndex
DGG.levels(::ISEA4RSystem) = 0:MAX_LEVEL
DGG.has_sorted_subtrees(::ISEA4RSystem) = true

"""
    max_neighbors(ISEA4RSystem(), connectivity) -> Int

`9` under `Vertex()`, `4` under `Edge()`.

The nine is the icosahedron showing through a square lattice. A cell in a
diamond's interior has the usual eight, and a corner cell on one of the ten
valence-3 icosahedron vertices has seven — but vertices 0 and 11 each carry
FIVE diamond-corners, so the corner cell there meets four other cells at the
vertex, two of which are not reached by any axis offset. Eight becomes nine at
exactly twenty cells per level. See `topology.jl` for the derivation and
`test/systems/ISEA4R/runtests.jl`, which counts the whole level.
"""
DGG.max_neighbors(::ISEA4RSystem, ::DGG.Vertex) = 9
DGG.max_neighbors(::ISEA4RSystem, ::DGG.Edge) = 4

# Deferred, deliberately: `chart.jl` also carries a row-major codec, which would
# make a perfectly good `ISEA4RRowMajorIndex` alternate scheme. It is not
# offered, because unlike HEALPix's RING there is no external file layout that
# wants it — nothing reads or writes ISEA4R row-major on disk — so it would be
# an id space with no consumer. `xyd_to_rowmajor`/`rowmajor_to_xyd` stay
# available inside the submodule for the chart's own arbitrary-`nside` use.
DGG.cellindextypes(::ISEA4RSystem) = (DGG.LevelIndex,)

function DGG.levelgrid(sys::ISEA4RSystem, l::Integer)
    lvl = Int(l)
    lvl in DGG.levels(sys) || throw(ArgumentError(
        "level $lvl is outside $(DGG.levels(sys)) for $(nameof(typeof(sys)))"))
    return ISEA4RGrid(lvl)
end

"""
    rootcells(ISEA4RSystem())

The ten level-0 cells: one per diamond, `LevelIndex(0, 0:9)`. At level 0 the
Morton code is empty, so the id *is* the diamond number.
"""
DGG.rootcells(::ISEA4RSystem) = [DGG.LevelIndex(0, d) for d in 0:(NDIAMONDS - 1)]

"""
    parent(ISEA4RSystem(), c) -> LevelIndex

The Morton parent: `index ÷ 4`, one level up. Throws an `ArgumentError` on a
level-0 cell, which has no parent.
"""
function Base.parent(::ISEA4RSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    l > 0 || throw(ArgumentError(
        "level-0 ISEA4R cell $c is a root and has no parent"))
    return DGG.LevelIndex(l - 1, c.index >> 2)
end

"""
    children(ISEA4RSystem(), c)

The four Morton children `4*index .+ (0:3)`, one level down, ascending.

Always exactly four: ISEA4R refinement is a uniform quadtree on every diamond,
with no pentagons and no exceptional cells — the icosahedron's twelve vertices
distort the *neighbourhood*, never the *subdivision*, because every vertex is a
corner of the chart square and quartering a square keeps it one.

In `(ix, iy)` terms the four are `(2ix, 2iy)`, `(2ix+1, 2iy)`, `(2ix, 2iy+1)`,
`(2ix+1, 2iy+1)` in that order, since `morton(2ix + a, 2iy + b) == 4*morton(ix,
iy) + a + 2b`. Throws an `ArgumentError` at `max_level`.
"""
function DGG.children(sys::ISEA4RSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    l < DGG.max_level(sys) || throw(ArgumentError(
        "ISEA4R cell $c is at max_level $(DGG.max_level(sys)) and has no children"))
    base = c.index << 2
    return [DGG.LevelIndex(l + 1, base + k) for k in 0:3]
end

"""
    ancestor(ISEA4RSystem(), c, l) -> LevelIndex

The ancestor at level `l`, in one shift: `index >> 2Δ`.

Sound because the id is `diamond * 4^level + morton(ix, iy)` and the Morton code
is positional — dropping the low `2Δ` bits drops `Δ` bits from each of `ix` and
`iy`, which is exactly `Δ` steps up the quadtree, and the diamond term divides
through untouched.
"""
function DGG.ancestor(sys::ISEA4RSystem, c::DGG.LevelIndex, l::Integer)
    target = Int(l)
    lc = DGG.level(c)
    target <= lc || throw(ArgumentError(
        "ancestor level $target is deeper than the cell's own level $lc"))
    target >= 0 || throw(ArgumentError(
        "ancestor level $target is above the root level 0"))
    return DGG.LevelIndex(target, c.index >> (2 * (lc - target)))
end

"""
    descendant_range(ISEA4RSystem(), c, l) -> UnitRange{Int}

The contiguous **positions** in `levelgrid(sys, l)` occupied by `c`'s level-`l`
descendants: `index * 4^Δ` through `(index + 1) * 4^Δ - 1` in 0-based Morton
ids, shifted into 1-based positions.

Exact and hole-free in both directions, which is what
`has_sorted_subtrees(sys) == true` asserts: Morton order is depth-first curve
order by construction, every position in the range names a real cell (ISEA4R
has no id gaps — ten dense diamonds, no pentagons), and sibling ranges tile the
parent's exactly.
"""
function DGG.descendant_range(sys::ISEA4RSystem, c::DGG.LevelIndex, l::Integer)
    target = Int(l)
    lc = DGG.level(c)
    target >= lc || throw(ArgumentError(
        "descendant level $target is above the cell's own level $lc"))
    target <= DGG.max_level(sys) || throw(ArgumentError(
        "descendant level $target is past max_level $(DGG.max_level(sys))"))
    shift = 2 * (target - lc)
    lo = c.index << shift
    hi = ((c.index + 1) << shift) - 1
    return Int(lo + 1):Int(hi + 1)
end

"""
    descendants(ISEA4RSystem(), c, l)

Every level-`l` descendant of `c`, ascending.

Morton ids are dense and subtree-contiguous, so this is
[`descendant_range`](@ref) read off as consecutive ids, with no `children`
recursion and no sort.
"""
function DGG.descendants(sys::ISEA4RSystem, c::DGG.LevelIndex, l::Integer)
    r = DGG.descendant_range(sys, c, l)      # validates `l` both ways
    target = Int(l)
    return [DGG.LevelIndex(target, i - 1) for i in r]
end

# ===========================================================================
# Grid interface
# ===========================================================================

DGG.system(::ISEA4RGrid) = ISEA4RSystem()
DGG.level(g::ISEA4RGrid) = g.level
DGG.ncells(g::ISEA4RGrid) = Int(_ncells(g.level))

function DGG.cellindex(g::ISEA4RGrid, i::Int)
    1 <= i <= DGG.ncells(g) || throw(BoundsError(g, i))
    return DGG.LevelIndex(g.level, i - 1)
end

"""
    cellposition(grid, c) -> Union{Int,Nothing}

Closed form: `index + 1` for a cell at the grid's own level and in range, and
`nothing` otherwise (a different level, or an id no cell has). Replaces the
fallback's linear scan, and never throws — a miss is an answer.
"""
function DGG.cellposition(g::ISEA4RGrid, c::DGG.LevelIndex)
    DGG.level(c) == g.level || return nothing
    0 <= c.index < _ncells(g.level) || return nothing
    return Int(c.index + 1)
end

# The id guard every geometry entry point needs: `morton_to_xyd` would happily
# de-interleave an id no cell has, yielding the geometry of a cell that does not
# exist rather than an error. (`cellposition` deliberately does NOT use this —
# there, a miss is `nothing`.)
@inline function _checked_index(g::ISEA4RGrid, c::DGG.LevelIndex)
    DGG.level(c) == g.level || throw(ArgumentError(
        "cell $c is at level $(DGG.level(c)), not the grid's level $(g.level)"))
    0 <= c.index < _ncells(g.level) || throw(ArgumentError(
        "ISEA4R id $(c.index) is out of range 0:$(_ncells(g.level) - 1) at level $(g.level)"))
    return c.index
end

# ===========================================================================
# Geometry
# ===========================================================================

# How finely a cell edge is broken into great-circle segments.
#
# Snyder ISEA maps the chart's straight lines to CURVES on the sphere, so the
# four corners alone are a poor spherical polygon — the same situation as
# HEALPix's chart lines, and unlike S2, whose cell edges are great circles.
# `cell_boundary` is contractually a ring of great-circle arcs, so the edges are
# densified.
#
# THE COUNT DOES NOT DEPEND ON THE LEVEL, for the reason set out at HEALPix's
# `BOUNDARY_SEGMENTS`: refinement is self-similar, so a cell has essentially the
# same shape at every level and a k-gon approximation of it has a
# level-independent RELATIVE error. Measured max relative area error of the
# densified ring against the exact `4π/(10·4^level)`, over whole levels 0-4
# (`test/systems/ISEA4R/runtests.jl` keeps this honest):
#
#     segments/edge |    1    |    2    |    4    |    8    |   16
#     rel. area err | 1.5e-1  | 2.1e-2  | 5.8e-3  | 1.6e-3  | 4.3e-4
#
# — O(segments^-2) from 2 on, as the chord error predicts, and flat across
# levels (1.0e-3, 1.3e-3, 1.6e-3, 1.5e-3, 1.7e-3 at levels 1-5). Eight segments
# per edge, 32 vertices per cell, 0.16% on area; the same count HEALPix settled
# on. The bare 4-gon at the head of that row is why the pre-redesign kernel's
# rings, which shipped four corners at every level, were not good enough for a
# spherical predicate.
#
# LEVEL 0 IS EXACT and is not evidence about any other level: a diamond's four
# rim edges are icosahedron edges, which Snyder maps to great-circle arcs, so a
# level-0 cell IS its 4-gon and its ring area matches the closed form to 7e-16.
# Every deeper level has interior chart edges, which curve.
#
# Powers of two matter beyond cost. An interior densification point of one cell
# is a lattice point of its neighbour's edge at the same level, and with a
# power-of-two count the two `xyd_to_point` arguments are the same `Float64`
# (`ix + 1 - i/nseg` and `ix + (nseg-i)/nseg` are both exact and equal), so
# shared edges come out BIT-identical WITHIN A DIAMOND and the tessellation is
# exact there. Across a diamond rim the two sides are two developments of one
# icosahedron edge and agree to ~4e-15 rad but on NO coordinate exactly — that
# is a property of the icosahedron, not of this constant.
#
# So: an incidence test over these rings may use a hashed container only within
# a diamond, and must use a tolerance across a rim. Two further traps, both
# pinned by `signed zeros cannot break a hashed incidence test` in the suite:
# `Set`/`Dict` compare with `isequal`, under which `-0.0` and `0.0` are DISTINCT
# though `==` calls them equal, so a closed form that emits a negative zero on
# one side of a join and a positive one on the other loses the match silently
# (the sibling S2 port hit exactly this on its cube seams); and cross-rim points
# are not bit-equal at all, so a `Set` intersection there reports a seam
# neighbour as sharing nothing rather than sharing an edge.
const BOUNDARY_SEGMENTS = 8

"""
    _perimeter_points(ix, iy, diamond, nside, nseg) -> Vector{UnitSphericalPoint}

The cell's boundary walked in `nseg` equal chart steps per edge, starting at the
`(x+, y+)` corner and running in the [`cell_corners`](@ref) order
`(x+,y+) → (x-,y+) → (x-,y-) → (x+,y-)`. `nseg == 1` reproduces
[`cell_corners`](@ref) exactly.

Implicitly closed — each edge contributes its start vertex and its interior
points, never its end vertex, so the next edge's start is not duplicated.
"""
function _perimeter_points(ix::Integer, iy::Integer, diamond::Integer,
        nside::Integer, nseg::Integer)
    n = nside
    x0 = Int64(ix)
    y0 = Int64(iy)
    pts = Vector{GO.UnitSphericalPoint{Float64}}(undef, 4 * nseg)
    k = 0
    for i in 0:(nseg - 1)          # (x+,y+) -> (x-,y+), along y = (iy+1)/n
        t = i / nseg
        pts[k += 1] = xyd_to_point((x0 + 1 - t) / n, (y0 + 1) / n, diamond)
    end
    for i in 0:(nseg - 1)          # (x-,y+) -> (x-,y-), along x = ix/n
        t = i / nseg
        pts[k += 1] = xyd_to_point(x0 / n, (y0 + 1 - t) / n, diamond)
    end
    for i in 0:(nseg - 1)          # (x-,y-) -> (x+,y-), along y = iy/n
        t = i / nseg
        pts[k += 1] = xyd_to_point((x0 + t) / n, y0 / n, diamond)
    end
    for i in 0:(nseg - 1)          # (x+,y-) -> (x+,y+), along x = (ix+1)/n
        t = i / nseg
        pts[k += 1] = xyd_to_point((x0 + 1) / n, (y0 + t) / n, diamond)
    end
    return pts
end

"""
    cell_boundary(grid, c) -> Vector{UnitSphericalPoint}

The cell's boundary ring, **counter-clockwise seen from outside the sphere** and
implicitly closed (the first vertex is not repeated).

The ring starts at the cell's `(x+, y+)` chart corner and runs
`(x+,y+) → (x-,y+) → (x-,y-) → (x+,y-)`, the [`cell_corners`](@ref) order.
Because ISEA4R cell edges are chart lines rather than geodesics, each edge is
densified into `BOUNDARY_SEGMENTS` great-circle segments — at every level, for
the reason set out there — so the ring has 32 vertices and the four corners are
vertices 1, 9, 17 and 25.

The winding is structural, not measured: the chart's corner order is
counter-clockwise in the `(x, y)` plane, both affine halves are
orientation-preserving (`imag(conj(a)·b) == 2π/5 > 0`, asserted at load for all
twenty half-maps), and `snyder_inv_xyz` preserves orientation from a face's
right-handed `(u, w)` frame onto the sphere seen from outside.

For the cell's **area**, prefer [`cell_area`](@ref), which is the exact
equal-area value in closed form rather than this polygon's.
"""
function DGG.cell_boundary(g::ISEA4RGrid, c::DGG.LevelIndex)
    nside = _nside(g.level)
    ix, iy, d = morton_to_xyd(_checked_index(g, c), nside)
    return _perimeter_points(ix, iy, d, nside, BOUNDARY_SEGMENTS)
end

"""
    cell_area(grid, c) -> Float64

The cell's area in steradians: `4π / (10 * 4^level)`, **exactly**, for every
cell of every level.

This is an override of the generic polygon area, and it is a correction rather
than only a speedup. Equal-areaness is the defining property of the chart —
Snyder ISEA is exactly equal-area per face and each affine half of a diamond
carries exactly `2π/5`, so a chart rectangle of area `A` covers solid angle
`A · 4π/10` everywhere, on every diamond — while the generic answer is the area
of the densified boundary polygon, which approaches this value from below and
still differs by up to 0.16% at `BOUNDARY_SEGMENTS = 8`. The closed form is the
true semantic; the polygon is the approximation of it. (At level 0 the two
agree to 7e-16, because a diamond's rim edges are icosahedron edges and
therefore great circles — see `BOUNDARY_SEGMENTS`.)

O(1), and independent of the boundary densification, so tightening
`BOUNDARY_SEGMENTS` changes geometric predicates but never an area.
"""
DGG.cell_area(g::ISEA4RGrid, c::DGG.LevelIndex) =
    (_checked_index(g, c); 4 * Float64(π) / _ncells(g.level))

"""
    cell_centroid(grid, c) -> UnitSphericalPoint

The cell centre: the chart evaluated at the lattice cell's midpoint
`((ix + 0.5)/nside, (iy + 0.5)/nside)`.

This is the centre *by definition* — the chart is exactly equal-area, so the
midpoint of the chart rectangle is the canonical centre — and it is strictly
interior, as [`cell_centroid`](@ref) requires. It is not the spherical centroid
of the published 4-gon; no equal-area DGGS claims that of its cell centres.
"""
function DGG.cell_centroid(g::ISEA4RGrid, c::DGG.LevelIndex)
    nside = _nside(g.level)
    ix, iy, d = morton_to_xyd(_checked_index(g, c), nside)
    return cell_center(ix, iy, d, nside)
end

# ===========================================================================
# node_extent — the subtree cap
# ===========================================================================

# How many chart samples per edge the subtree cap is built from. See
# `_subtree_cap` for why this number, and not the boundary's, sets the bound.
const CAP_EDGE_SEGMENTS = 8

"""
    _subtree_cap(ix, iy, diamond, nside) -> SphericalCap

The bounding cap of an ISEA4R cell — *and therefore of its whole subtree*.

# Why the cell's own cap bounds the subtree

A cell's four children tile its chart rectangle exactly, at every depth: the
child rectangles quarter the parent's, and the lattice division is BIT-exact
across resolutions (`fl(ix/n) === fl(2ix/2n)`, see [`xyd_to_point`](@ref)), so
the shared chart coordinates are the same `Float64`s. The chart is a
homeomorphism of the closed square onto the diamond's spherical patch — Snyder
is a per-face homeomorphism and the two affine halves agree on the seam — so a
descendant's geometry at any depth lies in the closed chart rectangle of this
cell, and a cap that covers the rectangle covers the subtree.

(This is what the aperture-7 icosahedral systems cannot say, and why the generic
`node_extent` has to inflate. Here `cap_inflation` is never consulted.)

# Why the radius is what it is

The centre is the cell centre. The maximum of `d(centre, ·)` over the closed
chart rectangle is attained on the PERIMETER, which is sampled at
`CAP_EDGE_SEGMENTS` points per edge, and in fact at a CORNER: all four corners
are samples (each edge's sampling starts at its start vertex), and the
measurement — a dense 17×17 sampling of every cell's chart rectangle, on every
diamond, at `nside ∈ (1, 2, 3, 4, 5, 8, 16)` — finds a worst overhang past the
four-corner cap of exactly `0.0`, seam-straddling cells included: the farthest
point of a cell from its centre is always one of its own four corners. That
figure is inherited from the pre-redesign face-grid layer's `cap_policy`, and
`test/systems/ISEA4R/runtests.jl` re-runs it as a standing test rather than
trusting the record.

`gap/2` on top is measured insurance, not the proof; read the same caveat as
HEALPix's `_subtree_cap` carries — it is a first-order slack that absorbs the
third-order shortfall of a Lipschitz argument the corner case makes unnecessary
anyway. Over-covering costs only pruning time, while under-covering is a silent
correctness bug (see the covering law).
"""
function _subtree_cap(ix::Integer, iy::Integer, diamond::Integer, nside::Integer)
    center = cell_center(ix, iy, diamond, nside)
    pts = _perimeter_points(ix, iy, diamond, nside, CAP_EDGE_SEGMENTS)
    rmax = 0.0
    gap = 0.0
    prev = pts[end]
    for p in pts
        rmax = max(rmax, US.spherical_distance(center, p))
        gap = max(gap, US.spherical_distance(prev, p))
        prev = p
    end
    radius = min(Float64(π), rmax + gap / 2)
    return SphericalCap(center, nextfloat(radius))
end

"""
    node_extent(ISEA4RSystem(), c) -> SphericalCap

The subtree cap of `c` — an override of the generic inflated default, and the
reason [`cap_inflation`](@ref) is never consulted for this system.

What is **exact** is the nesting: children tile their parent's chart rectangle
bit-exactly and the chart is a homeomorphism, so the cell's own bounding cap
already bounds the whole subtree and there is nothing to inflate for.

What is **measured** is the cap's RADIUS — a sampled perimeter (which does
capture all four corners exactly, and the corners are where the maximum sits)
plus a slack term. See `_subtree_cap`.

The extent is geodesically convex at every level: the widest cap in the system
is a level-0 diamond's, whose radius is the median of a spherical face triangle
(58.3°) plus the sampling slack, **62.34°** as shipped — comfortably inside the
90° the conformance suite requires for vertex sampling to be a sound proxy for
the whole boundary. Level 1 is 32.7°, and each level roughly halves it.
"""
function DGG.node_extent(::ISEA4RSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    nside = _nside(l)
    ix, iy, d = morton_to_xyd(c.index, nside)
    return _subtree_cap(ix, iy, d, nside)
end

# ===========================================================================
# Location
# ===========================================================================

"""
    cellat(grid, p::UnitSphericalPoint) -> LevelIndex

The cell containing `p`, in closed form — [`point_to_morton`](@ref), the chart's
analytic inverse, with no tree descent and no point-in-polygon test.

Never `nothing`: a complete ISEA4R level grid covers the sphere, because the ten
diamonds do.

**Ties.** A point exactly on a cell boundary is legitimately contained by every
cell meeting there, so which one is returned is a tie. It is broken by the
arithmetic of [`point_to_xyd`](@ref) — `snyder_fwd`'s nearest-face-centre choice
(lowest face index on an exact tie) for the diamond, then `floor`, which puts
the point on the higher side of each chart cut line. Deterministic per platform
and self-consistent: the returned cell's own centroid maps back to it, which
`test/systems/ISEA4R/runtests.jl` asserts over whole levels.
"""
DGG.cellat(g::ISEA4RGrid, p::GO.UnitSphericalPoint) =
    DGG.LevelIndex(g.level, point_to_morton(p, _nside(g.level)))

# ===========================================================================
# Topology
# ===========================================================================

"""
    _one_ring(grid, c, connectivity) -> SmallVector{9,LevelIndex}

The immediate neighbours of `c` in **counter-clockwise rotational order seen
from outside the sphere, starting at the `(+1, 0)` chart direction**, from
[`lattice_neighbors`](@ref).
"""
function _one_ring(g::ISEA4RGrid, c::DGG.LevelIndex, connectivity::DGG.Connectivity)
    nside = _nside(g.level)
    ix, iy, d = morton_to_xyd(_checked_index(g, c), nside)
    out = SmallVector{9,DGG.LevelIndex}()
    for (jx, jy, jd) in lattice_neighbors(ix, iy, d, nside, connectivity)
        out = SmallCollections.push(out,
            DGG.LevelIndex(g.level, xyd_to_morton(jx, jy, jd, nside)))
    end
    return out
end

"""
    neighbors(grid, c, k = 1; connectivity = Vertex())

The cells within `k` lattice steps of `c`, excluding `c`, in **rotational
order**: the rings `1:k` concatenated outward, each ring counter-clockwise seen
from outside the sphere.

So `ring(grid, c, k)` is exactly the trailing block of `neighbors(grid, c, k)`,
and `neighbors(grid, c, k)` is `vcat(ring(grid, c, 1), ..., ring(grid, c, k))`.

# Connectivity

`Vertex()` (the default) is the 3×3 chart neighbourhood: 8 cells in a diamond's
interior, 7 at a corner cell on one of the ten valence-3 icosahedron vertices,
and 9 at the twenty corner cells on vertices 0 and 11, where five diamonds meet
and the diagonal offset yields two cells instead of none.

`Edge()` keeps only the four that share a whole cell edge. Because an ISEA4R
cell is an axis-aligned square in its chart, the edge-sharing neighbours are the
four AXIS offsets and the corner-only ones are the four diagonals — the opposite
pairing from HEALPix, whose pixel is a diamond rotated 45° against its lattice.

# Order

**Counter-clockwise seen from outside the sphere, starting at the `(+1, 0)`
chart direction** — the offset cycle `(+1,0), (+1,+1), (0,+1), (-1,+1), (-1,0),
(-1,-1), (0,-1), (+1,-1)` under `Vertex()` and its restriction to the four axis
offsets under `Edge()`. These are *chart* directions, not compass ones: a
diamond is not aligned with anything on the globe, and the same slot points a
different way on each of the ten.

Where an offset yields no cell (a valence-3 corner) it simply drops out of the
cycle; where it yields two (vertices 0 and 11) both are emitted in fan order
about the vertex, which is the same sweep. There is no padding and no gap
marker.

Rings beyond the first are ordered by azimuth about the cell centre from the
same starting spoke; see [`ring`](@ref).

`k == 0` returns an empty container; `k == 1` returns a
`SmallCollections.SmallVector` sized by [`max_neighbors`](@ref).
"""
function DGG.neighbors(g::ISEA4RGrid, c::DGG.LevelIndex, k::Integer = 1;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return SmallVector{9,DGG.LevelIndex}()
    steps == 1 && return _one_ring(g, c, connectivity)
    shells = _shells(g, c, steps, connectivity)
    isempty(shells) && return DGG.LevelIndex[]
    return reduce(vcat, shells)
end

"""
    ring(grid, c, k; connectivity = Vertex())

The cells at lattice distance **exactly** `k`, counter-clockwise seen from
outside the sphere. `ring(grid, c, 0)` is `[c]`.

`k == 1` is the chart offset cycle (see [`neighbors`](@ref)). For `k >= 2` there
is no lattice cycle to read off — an outer ring crosses diamond rims and
icosahedron vertices arbitrarily — so the shell is ordered **by azimuth about
the cell centre, measured counter-clockwise from the direction of the FIRST
ring-1 neighbour**. That is the extension the [`neighbors`](@ref) contract
recommends, and it makes every ring start on the same spoke by construction
rather than by a second convention. Ties in azimuth break by canonical id, so
the order is total and deterministic.
"""
function DGG.ring(g::ISEA4RGrid, c::DGG.LevelIndex, k::Integer;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return DGG.LevelIndex[c]
    shells = _shells(g, c, steps, connectivity)
    steps <= length(shells) || return DGG.LevelIndex[]
    return shells[steps]
end

# Breadth-first expansion over the lattice one-ring; shell `j` is the set at
# distance exactly `j`, each returned in CCW rotational order. Shared by
# `neighbors` and `ring` so the two cannot disagree about what a shell is or
# what order it is in.
function _shells(g::ISEA4RGrid, c::DGG.LevelIndex, steps::Int,
        connectivity::DGG.Connectivity)
    shells = Vector{DGG.LevelIndex}[]
    steps == 0 && return shells
    seen = Set{DGG.LevelIndex}((c,))
    frontier = DGG.LevelIndex[c]
    for j in 1:steps
        next = DGG.LevelIndex[]
        for x in frontier
            for y in _one_ring(g, x, connectivity)
                y in seen && continue
                push!(seen, y)
                push!(next, y)
            end
        end
        # Shell 1 is already the lattice cycle, in order. Outer shells come out
        # of the breadth-first walk in whatever order the frontier happened to
        # reach them, so they are wound geometrically, from the spoke shell 1
        # starts on.
        if j > 1 && !isempty(first(shells))
            _sort_ccw!(next, g, c, first(first(shells)))
        end
        push!(shells, next)
        isempty(next) && break
        frontier = next
    end
    return shells
end

# Order a shell counter-clockwise about `c`'s centre, from the azimuth of the
# first ring-1 neighbour. `e1 × e2 == centre` makes `(e1, e2)` right-handed SEEN
# FROM OUTSIDE, which is what puts increasing `atan(u·e2, u·e1)`
# counter-clockwise from outside rather than from inside.
function _sort_ccw!(cells::Vector{DGG.LevelIndex}, g::ISEA4RGrid,
        c::DGG.LevelIndex, reference::DGG.LevelIndex)
    length(cells) <= 1 && return cells
    centre = DGG.cell_centroid(g, c)
    e1, e2 = _tangent_basis(centre, DGG.cell_centroid(g, reference))
    key(x) = begin
        p = DGG.cell_centroid(g, x)
        (mod(_azimuth(centre, e1, e2, p), 2 * Float64(π)), x)
    end
    sort!(cells; by = key)
    return cells
end

# A right-handed-from-outside tangent basis at `centre`, with `e1` pointing at
# `toward` (projected into the tangent plane).
function _tangent_basis(centre, toward)
    u = (toward[1] - centre[1], toward[2] - centre[2], toward[3] - centre[3])
    dot = u[1] * centre[1] + u[2] * centre[2] + u[3] * centre[3]
    t = (u[1] - dot * centre[1], u[2] - dot * centre[2], u[3] - dot * centre[3])
    n = sqrt(t[1]^2 + t[2]^2 + t[3]^2)
    # A degenerate reference (a neighbour whose centroid projects onto the cell
    # centre) cannot happen for a real cell, but a fallback keeps this total.
    n <= eps(Float64) && (t = abs(centre[3]) < 0.9 ? (0.0, 0.0, 1.0) : (1.0, 0.0, 0.0);
                          n = 1.0)
    e1 = (t[1] / n, t[2] / n, t[3] / n)
    e2 = (centre[2] * e1[3] - centre[3] * e1[2],
          centre[3] * e1[1] - centre[1] * e1[3],
          centre[1] * e1[2] - centre[2] * e1[1])
    return e1, e2
end

function _azimuth(centre, e1, e2, p)
    u = (p[1] - centre[1], p[2] - centre[2], p[3] - centre[3])
    return atan(u[1] * e2[1] + u[2] * e2[2] + u[3] * e2[3],
                u[1] * e1[1] + u[2] * e1[2] + u[3] * e1[3])
end
