# ---------------------------------------------------------------------------
# The HEALPix system on the grid interface
#
# Nested HEALPix is the interface's reference DENSE-ORDINAL system: at level `l`
# the cells are exactly `0:12*4^l - 1`, pixel `p`'s children are `4p:4p+3`, and
# its parent is `p ÷ 4`. Every hierarchy method below is that arithmetic, and
# nothing here needs a lookup table.
#
# ## The id, and the base of its numbering
#
# The canonical id is the interface's own `LevelIndex(level, index)`, where
#
#     index  ==  the 0-BASED nested (Morton) pixel id at that level
#     position in `levelgrid(sys, level)`  ==  index + 1
#
# 0-based is not a coin flip. It is the EOPF/ESA convention the pre-redesign
# `HealpixLookup` stored on disk and the one `chart.jl`'s NESTED codec is
# written in, so choosing 1-based would have put a `+1` on every boundary
# between this package and a HEALPix file, and — worse — broken the arithmetic
# that is the entire reason nested ordering exists: `4p:4p+3` and `p ÷ 4` are
# child/parent maps only in 0-based numbering. (Healpix.jl's own `nest2ring` /
# `ring2nest` take 1-based pixel numbers, which is why the oracle tests here
# convert with an explicit `+ 1`.)
#
# The `+ 1` therefore lives in exactly one place — the position/index
# conversion in `cellindex` / `cellposition` — where the interface's "a bare
# `Int` is a position" rule makes it visible rather than ambient.
#
# The RING numbering is offered as the alternate scheme `HEALPixRingIndex`,
# reached with `reindex`. It is 1-based, because a ring index doubles as the
# position in a ring-ordered data vector (the convention Healpix.jl and
# SpeedyWeather's RingGrids both use), and making it 0-based would move an
# off-by-one onto every one of those call sites instead.
#
# ## What is a fast path here
#
#   * `cellat` is closed-form (`point_to_nested`) — no tree descent.
#   * `descendant_range` is `[p*4^Δ, (p+1)*4^Δ)` shifted into position space.
#   * `node_extent` is the EXACT subtree cap, not the inflated default.
#   * `neighbors` / `ring` walk the lattice, not the geometry.
#   * `ancestor` drops `2Δ` bits in one shift.
# ---------------------------------------------------------------------------

# ===========================================================================
# Types
# ===========================================================================

"""
    HEALPixSystem() <: AbstractHierarchicalGridSystem

The HEALPix hierarchical grid system in **nested** (Morton) ordering, on the
unit sphere.

Twelve equal-area base pixels, each refined by aperture 4: level `l` has
`12 * 4^l` pixels of exactly equal solid angle `4π / (12 * 4^l)`. Cells are
quadrilaterals in the HEALPix chart (`chart.jl`); their edges follow chart
lines, not great circles, so [`cell_boundary`](@ref) densifies them — see that
method.

# Ids

The canonical cell index is [`LevelIndex`](@ref), whose `index` field holds the
**0-based** nested pixel id; a cell's position in `levelgrid(sys, l)` is
`index + 1`. [`HEALPixRingIndex`](@ref) is the alternate scheme, reached with
[`reindex`](@ref).

# Levels

`0:29`. The bound is the Int64 codec's, not HEALPix's: at level 29 there are
`12 * 4^29 = 3458764513820540928` pixels, and level 30 would overflow a signed
64-bit cell count.

# Traits

`has_sorted_subtrees` is `true` — nested order *is* depth-first curve order, so
a subtree occupies a contiguous run of positions and
[`descendant_range`](@ref) is exact and hole-free.
[`max_neighbors`](@ref) is 8 under `Vertex()` and 4 under `Edge()`.
[`node_extent`](@ref) is overridden with the pixel's own bounding cap — nested
parents *are* the union of their children, so nothing needs inflating and
[`cap_inflation`](@ref) is never consulted.
"""
struct HEALPixSystem <: DGG.AbstractHierarchicalGridSystem end

"""
    HEALPixGrid(level) <: AbstractGrid

The complete HEALPix grid at refinement `level`: all `12 * 4^level` pixels in
nested order. Built by [`levelgrid`](@ref); a lightweight descriptor, not a
materialised cell list.
"""
struct HEALPixGrid <: DGG.AbstractGrid
    level::Int
end

"""
    HEALPixRingIndex(level, index) <: AbstractCellIndex

A HEALPix cell named in **RING** ordering: pixels numbered north-to-south along
iso-latitude rings, west-to-east within a ring.

`index` is **1-based**, unlike the 0-based nested `index` of the canonical
[`LevelIndex`](@ref), because a ring index doubles as the position of the pixel
in a ring-ordered data vector — the convention Healpix.jl and RingGrids both
index HEALPix fields with. Converting to Healpix.jl's `ring2nest`/`nest2ring`
therefore needs no shift on the ring side, and a `+ 1` on the nested side.

This is an alternate scheme, not the canonical one: it does **not** sort in
subtree order (ring-ordered siblings are scattered across rings), which is
exactly why nested is canonical. Reach it with
[`reindex`](@ref)`(HEALPixRingIndex, HEALPixSystem(), c)`.
"""
struct HEALPixRingIndex <: DGG.AbstractCellIndex
    level::Int32
    index::Int64
end

DGG.level(c::HEALPixRingIndex) = Int(c.level)
DGG.rawid(c::HEALPixRingIndex) = c.index
Base.isless(a::HEALPixRingIndex, b::HEALPixRingIndex) =
    isless((a.level, a.index), (b.level, b.index))
Base.show(io::IO, c::HEALPixRingIndex) =
    print(io, "HEALPixRingIndex(", c.level, ", ", c.index, ")")

# The one place `nside` is derived from a level, and the guard that keeps a
# bad level from silently producing a shift of 64.
@inline _nside(level::Integer) = Int64(1) << Int(level)
@inline _npix(level::Integer) = 12 * (Int64(1) << (2 * Int(level)))

const MAX_LEVEL = 29

# ===========================================================================
# System interface
# ===========================================================================

DGG.cellindextype(::HEALPixSystem) = DGG.LevelIndex
DGG.cellindextypes(::HEALPixSystem) = (DGG.LevelIndex, HEALPixRingIndex)
DGG.levels(::HEALPixSystem) = 0:MAX_LEVEL
DGG.has_sorted_subtrees(::HEALPixSystem) = true

DGG.max_neighbors(::HEALPixSystem, ::DGG.Vertex) = 8
DGG.max_neighbors(::HEALPixSystem, ::DGG.Edge) = 4

function DGG.levelgrid(sys::HEALPixSystem, l::Integer)
    lvl = Int(l)
    lvl in DGG.levels(sys) || throw(ArgumentError(
        "level $lvl is outside $(DGG.levels(sys)) for $(nameof(typeof(sys)))"))
    return HEALPixGrid(lvl)
end

DGG.rootcells(::HEALPixSystem) = [DGG.LevelIndex(0, i) for i in 0:11]

"""
    parent(HEALPixSystem(), c) -> LevelIndex

The nested parent: `index ÷ 4`, one level up. Throws an `ArgumentError` on a
level-0 cell, which has no parent.
"""
function Base.parent(::HEALPixSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    l > 0 || throw(ArgumentError(
        "level-0 HEALPix cell $c is a root and has no parent"))
    return DGG.LevelIndex(l - 1, c.index >> 2)
end

"""
    children(HEALPixSystem(), c)

The four nested children `4*index .+ (0:3)`, one level down, ascending.

Always exactly four: HEALPix refinement is a uniform quadtree with no
pentagons and no exceptional cells, so unlike the icosahedral systems this
count never varies. Throws an `ArgumentError` at `max_level`.
"""
function DGG.children(sys::HEALPixSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    l < DGG.max_level(sys) || throw(ArgumentError(
        "HEALPix cell $c is at max_level $(DGG.max_level(sys)) and has no children"))
    base = c.index << 2
    return [DGG.LevelIndex(l + 1, base + k) for k in 0:3]
end

"""
    ancestor(HEALPixSystem(), c, l) -> LevelIndex

The ancestor at level `l`, in one shift: `index >> 2Δ`.

Sound because the nested id is `face * 4^level + morton(ix, iy)` and the Morton
code is positional — dropping the low `2Δ` bits drops `Δ` bits from each of
`ix` and `iy`, which is exactly `Δ` steps up the quadtree, and the face term
divides through untouched.
"""
function DGG.ancestor(sys::HEALPixSystem, c::DGG.LevelIndex, l::Integer)
    target = Int(l)
    lc = DGG.level(c)
    target <= lc || throw(ArgumentError(
        "ancestor level $target is deeper than the cell's own level $lc"))
    target >= 0 || throw(ArgumentError(
        "ancestor level $target is above the root level 0"))
    return DGG.LevelIndex(target, c.index >> (2 * (lc - target)))
end

"""
    descendant_range(HEALPixSystem(), c, l) -> UnitRange{Int}

The contiguous **positions** in `levelgrid(sys, l)` occupied by `c`'s level-`l`
descendants: `index * 4^Δ` through `(index + 1) * 4^Δ - 1` in 0-based nested
ids, shifted into 1-based positions.

Exact and hole-free in both directions, which is what
`has_sorted_subtrees(sys) == true` asserts: nested order is depth-first curve
order by construction, every position in the range names a real pixel (HEALPix
has no id gaps), and sibling ranges tile the parent's exactly.
"""
function DGG.descendant_range(sys::HEALPixSystem, c::DGG.LevelIndex, l::Integer)
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
    descendants(HEALPixSystem(), c, l)

Every level-`l` descendant of `c`, ascending.

Nested ids are dense and subtree-contiguous, so this is
[`descendant_range`](@ref) read off as consecutive ids — one subtraction per
cell, with no `children` recursion and no sort. The generic fallback would
reach the same answer through `levelgrid`/`cellindex` per position; this skips
that indirection.
"""
function DGG.descendants(sys::HEALPixSystem, c::DGG.LevelIndex, l::Integer)
    r = DGG.descendant_range(sys, c, l)      # validates `l` both ways
    target = Int(l)
    return [DGG.LevelIndex(target, i - 1) for i in r]
end

# ===========================================================================
# Grid interface
# ===========================================================================

DGG.system(::HEALPixGrid) = HEALPixSystem()
DGG.level(g::HEALPixGrid) = g.level
DGG.ncells(g::HEALPixGrid) = Int(_npix(g.level))

function DGG.cellindex(g::HEALPixGrid, i::Int)
    1 <= i <= DGG.ncells(g) || throw(BoundsError(g, i))
    return DGG.LevelIndex(g.level, i - 1)
end

"""
    cellposition(grid, c) -> Union{Int,Nothing}

Closed form: `index + 1` for a cell at the grid's own level and in range, and
`nothing` otherwise (a different level, or an id no pixel has). Replaces the
fallback's linear scan.
"""
function DGG.cellposition(g::HEALPixGrid, c::DGG.LevelIndex)
    DGG.level(c) == g.level || return nothing
    0 <= c.index < _npix(g.level) || return nothing
    return Int(c.index + 1)
end

# The range guard runs BEFORE the conversion, and that ordering is the whole
# point: `ring_to_xyf` throws an `ArgumentError` on an out-of-range ring index,
# but a cell that is not in the grid is contractually `nothing`, not an error.
function DGG.cellposition(g::HEALPixGrid, c::HEALPixRingIndex)
    DGG.level(c) == g.level || return nothing
    1 <= c.index <= _npix(g.level) || return nothing      # RING ids are 1-based
    return DGG.cellposition(g, DGG.reindex(DGG.LevelIndex, HEALPixSystem(), c))
end

# ===========================================================================
# Alternate index scheme
# ===========================================================================

function DGG.reindex(::Type{HEALPixRingIndex}, ::HEALPixSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    return HEALPixRingIndex(l, nested_to_ring(c.index, _nside(l)))
end

function DGG.reindex(::Type{DGG.LevelIndex}, ::HEALPixSystem, c::HEALPixRingIndex)
    l = DGG.level(c)
    return DGG.LevelIndex(l, ring_to_nested(c.index, _nside(l)))
end

# ===========================================================================
# Geometry
# ===========================================================================

# How finely a cell edge is broken into great-circle segments.
#
# HEALPix pixel edges are chart lines, NOT great circles, so the four corners
# alone are a poor spherical polygon: on a level-0 pixel the true edge departs
# from the corner-to-corner geodesic by about 4.4 degrees. `cell_boundary` is
# contractually a ring of great-circle arcs, so the edges are densified.
#
# THE COUNT DOES NOT DEPEND ON THE LEVEL, and that is the whole design note.
# The tempting schedule is to densify coarse levels and stop at deep ones, on
# the theory that deep cells are small enough to be flat. That reasoning holds
# the ABSOLUTE error constant, which is the wrong invariant for a tessellation:
# refinement is self-similar, so a HEALPix cell has essentially the SAME SHAPE
# at every level and a k-gon approximation of it therefore has a
# level-independent RELATIVE error. Measured max relative area error over whole
# levels (`test/systems/HEALPix/runtests.jl` keeps this honest):
#
#     segments/edge |    1    |    2    |    4    |    8    |   16    |   32
#     rel. area err | 9.9e-2  | 2.8e-2  | 7.3e-3  | 1.8e-3  | 4.6e-4  | 1.1e-4
#
# — flat across levels 0-6, and O(segments^-2) as the chord error predicts. So
# a fixed count it is: 8 segments per edge, 32 vertices per cell, 0.18% on
# area. (The pre-redesign tree shipped the bare 4 corners at every level, i.e.
# the 9.9e-2 column, to ConservativeRegridding.)
#
# Powers of two matter beyond cost: an interior densification point of one
# pixel is a lattice point of its neighbour's edge at the same level, and with
# a power-of-two count the two `xyf_to_point` arguments are the same Float64,
# so shared edges come out BIT-identical and the tessellation stays exact —
# which is what makes conservative regridding conserve regardless of the
# residual per-cell area error above.
const BOUNDARY_SEGMENTS = 8

"""
    _perimeter_points(ix, iy, face, nside, nseg) -> Vector{UnitSphericalPoint}

The pixel's boundary walked in `nseg` equal chart steps per edge, starting at
the north corner and running north → west → south → east: the
[`pixel_corners`](@ref) order, densified. `nseg == 1` reproduces
`pixel_corners` exactly.

Implicitly closed — each edge contributes its start vertex and its interior
points, never its end vertex, so the next edge's start is not duplicated.
"""
function _perimeter_points(ix::Integer, iy::Integer, face::Integer,
        nside::Integer, nseg::Integer)
    n = nside
    x0 = Int64(ix)
    y0 = Int64(iy)
    pts = Vector{GO.UnitSphericalPoint{Float64}}(undef, 4 * nseg)
    k = 0
    for i in 0:(nseg - 1)          # north -> west, along y = (iy+1)/n
        t = i / nseg
        pts[k += 1] = xyf_to_point((x0 + 1 - t) / n, (y0 + 1) / n, face)
    end
    for i in 0:(nseg - 1)          # west -> south, along x = ix/n
        t = i / nseg
        pts[k += 1] = xyf_to_point(x0 / n, (y0 + 1 - t) / n, face)
    end
    for i in 0:(nseg - 1)          # south -> east, along y = iy/n
        t = i / nseg
        pts[k += 1] = xyf_to_point((x0 + t) / n, y0 / n, face)
    end
    for i in 0:(nseg - 1)          # east -> north, along x = (ix+1)/n
        t = i / nseg
        pts[k += 1] = xyf_to_point((x0 + 1) / n, (y0 + t) / n, face)
    end
    return pts
end

"""
    cell_boundary(grid, c) -> Vector{UnitSphericalPoint}

The pixel's boundary ring, **counter-clockwise seen from outside the sphere**
and implicitly closed (the first vertex is not repeated).

The ring starts at the pixel's north corner and runs north → west → south →
east, the `pixel_corners` order that `Healpix.boundariesRing` also emits.
Because HEALPix edges are chart lines rather than geodesics, each edge is
densified into `BOUNDARY_SEGMENTS` great-circle segments — at every level, for
the reason set out there — so the ring has 32 vertices and the four corners are
vertices 1, 9, 17 and 25.

For the cell's **area**, prefer [`cell_area`](@ref), which is the exact
equal-area value in closed form rather than this polygon's.
"""
function DGG.cell_boundary(g::HEALPixGrid, c::DGG.LevelIndex)
    nside = _nside(g.level)
    ix, iy, face = nested_to_xyf(_checked_index(g, c), nside)
    return _perimeter_points(ix, iy, face, nside, BOUNDARY_SEGMENTS)
end

"""
    cell_area(grid, c) -> Float64

The pixel's area in steradians: `4π / (12 * 4^level)`, **exactly**, for every
pixel of every level.

This is an override of the generic polygon area, and it is a correction rather
than only a speedup. Equal-areaness is HEALPix's defining property — the chart
is an equal-area map by construction, so every pixel at a level subtends
*identically* the same solid angle — while the generic answer is the area of
the densified boundary polygon, which approaches this value from below and
still differs from it by ~0.18% at `BOUNDARY_SEGMENTS = 8`. The closed form is
the true semantic; the polygon is the approximation of it.

O(1), and independent of the boundary densification, so tightening
`BOUNDARY_SEGMENTS` changes geometric predicates but never an area.
"""
DGG.cell_area(g::HEALPixGrid, c::DGG.LevelIndex) =
    (_checked_index(g, c); 4 * Float64(π) / _npix(g.level))

"""
    cell_centroid(grid, c) -> UnitSphericalPoint

The HEALPix pixel centre: the chart evaluated at the lattice cell's midpoint.

This is the pixel centre *by definition* — the chart is equal-area, so the
midpoint of the chart square is the canonical centre — and it is strictly
interior, as [`cell_centroid`](@ref) requires. Agrees with
`Healpix.pix2vecNest` to ~9e-16 per coordinate.
"""
function DGG.cell_centroid(g::HEALPixGrid, c::DGG.LevelIndex)
    nside = _nside(g.level)
    ix, iy, face = nested_to_xyf(_checked_index(g, c), nside)
    return pixel_center(ix, iy, face, nside)
end

# The id guard every geometry entry point needs: `nested_to_xyf` will happily
# un-Morton an id no pixel has, yielding the geometry of a cell that does not
# exist rather than an error.
@inline function _checked_index(g::HEALPixGrid, c::DGG.LevelIndex)
    DGG.level(c) == g.level || throw(ArgumentError(
        "cell $c is at level $(DGG.level(c)), not the grid's level $(g.level)"))
    0 <= c.index < _npix(g.level) || throw(ArgumentError(
        "nested id $(c.index) is out of range 0:$(_npix(g.level) - 1) at level $(g.level)"))
    return c.index
end

# ===========================================================================
# node_extent — the subtree cap
# ===========================================================================

# How many chart samples per edge the subtree cap is built from. See
# `_subtree_cap` for why this number, and not the boundary's, sets the bound.
const CAP_EDGE_SEGMENTS = 8

"""
    _subtree_cap(ix, iy, face, nside) -> SphericalCap

The bounding cap of a HEALPix pixel — *and therefore of its whole subtree*.

# Why the pixel's own cap bounds the subtree

A nested HEALPix parent is the exact geographic union of its four children, at
every depth — refinement subdivides the chart square and nothing pokes out.
(This is what the aperture-7 icosahedral systems cannot say, and why the
generic `node_extent` has to inflate.) So every descendant boundary point at
every depth lies in the closed chart square of this pixel, and a cap that
covers the square covers the subtree. That is the whole content of the
`node_extent` override, and it is O(1) at every level.

# Why the radius is what it is

The centre is the pixel centre. `d(centre, ·)` attains its maximum over the
closed square on the square's PERIMETER — the centre is interior and the
region is well inside a hemisphere, so there is no interior maximum — and the
perimeter is sampled at `CAP_EDGE_SEGMENTS` points per edge.

**What actually bounds it is the corners.** On a chart square the perimeter
maximum of `d(centre, ·)` is attained at a CORNER, and all four corners are
samples (the sampling starts each edge at its start vertex). So `rmax` is not
an approximation of the maximum over the perimeter — it *is* that maximum,
measured exactly. That, plus the nesting argument above, is the bound.

**`gap/2` is measured insurance on top of it, not the proof.** The tempting
justification — `d(centre, ·)` is 1-Lipschitz, so between two consecutive
samples the true distance cannot exceed the larger sample's by more than half
their separation — is not airtight as written. Lipschitz-ness bounds the excess
by half the ARC LENGTH along the perimeter between the samples, whereas `gap`
is the great-circle distance between them, a CHORD of that path. The arc
exceeds the chord by O(gap³), so for a hypothetical non-corner maximum `gap/2`
would under-cover by that third-order term. Do not read it as an exact
Lipschitz bound. It does not matter here, because the corner argument already
gives the true maximum, and the first-order slack `gap/2` adds on top absorbs
the third-order shortfall many times over: the tests re-measure the true
perimeter at 32x this sampling and confirm the cap strictly contains it, with
the slack never even approached.

The result is looser than a corner-only cap by a few per cent and tighter than
the generic inflated default; over-covering costs only pruning time, while
under-covering is a silent correctness bug (see the covering law).
"""
function _subtree_cap(ix::Integer, iy::Integer, face::Integer, nside::Integer)
    center = pixel_center(ix, iy, face, nside)
    pts = _perimeter_points(ix, iy, face, nside, CAP_EDGE_SEGMENTS)
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
    node_extent(HEALPixSystem(), c) -> SphericalCap

The subtree cap of `c` — an override of the generic inflated default, and the
reason `cap_inflation` is never consulted for this system.

What is **exact** is the nesting: a nested HEALPix parent *is* the geographic
union of its children, at every depth, so the pixel's own bounding cap already
bounds the whole subtree and there is nothing to inflate for.

What is **measured** is the cap's RADIUS — a sampled perimeter (which does
capture all four corners exactly, and the corners are where the maximum sits)
plus a slack term. See `_subtree_cap` for the construction and for how far that
argument does and does not go.
"""
function DGG.node_extent(::HEALPixSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    nside = _nside(l)
    ix, iy, face = nested_to_xyf(c.index, nside)
    return _subtree_cap(ix, iy, face, nside)
end

# ===========================================================================
# Location
# ===========================================================================

"""
    cellat(grid, p::UnitSphericalPoint) -> LevelIndex

The pixel containing `p`, in closed form — `point_to_nested`, the chart's
analytic inverse, with no tree descent and no point-in-polygon test.

Never `nothing`: a complete HEALPix level grid covers the sphere.

**Ties.** A point exactly on a pixel boundary is legitimately contained by
every cell meeting there, so which one is returned is a tie. It is broken by
the `floor` arithmetic of `point_to_xyf`, which puts the point on the higher
side of each cut line: **deterministic and self-consistent** (the returned
cell's own centroid maps back to it), which is what the interface requires.

On points **interior** to a cell this agrees exactly with Healpix.jl —
`vec2pixNest` fed the identical Cartesian point, at every level, is asserted in
`test/systems/HEALPix/runtests.jl`. On **boundary** points the two libraries do
not always agree, and neither is wrong: at the two known 4-way lattice corners
`(1, 0, 0)` and `(45°, -asin(2/3))` this port answers pixels 17 and 33 at
level 1 where `vec2pixNest` answers 19 and 35, each library being
self-consistent. Do not "fix" this by matching Healpix.jl; assert the
contractual properties instead.
"""
DGG.cellat(g::HEALPixGrid, p::GO.UnitSphericalPoint) =
    DGG.LevelIndex(g.level, point_to_nested(p, g.level))

# ===========================================================================
# Topology
# ===========================================================================

"""
    _one_ring(grid, c, connectivity) -> SmallVector{8,LevelIndex}

The immediate neighbours of `c` in **counter-clockwise rotational order seen
from outside the sphere, starting at the `SW` lattice direction**
(`_neighbor_cycle`), deduplicated, with non-existent entries dropped.

The dedup matters only at level 0, where `nside == 1` makes every lattice
offset wrap through the face tables and two offsets could in principle name the
same base pixel; keeping the FIRST occurrence is what preserves the cycle.
"""
function _one_ring(g::HEALPixGrid, c::DGG.LevelIndex, connectivity::DGG.Connectivity)
    _checked_index(g, c)
    raw = nested_neighbors(c.index, g.level)
    out = SmallVector{8,DGG.LevelIndex}()
    for m in _neighbor_cycle(connectivity)
        id = raw[m]
        id < 0 && continue
        nb = DGG.LevelIndex(g.level, id)
        nb in out && continue
        out = SmallCollections.push(out, nb)
    end
    return out
end

"""
    neighbors(grid, c, k = 1; connectivity = Vertex())

The cells within `k` lattice steps of `c`, excluding `c`, in **rotational
order**: the rings `1:k` concatenated outward, each ring counter-clockwise seen
from outside the sphere.

So `ring(grid, c, k)` is exactly the trailing block of
`neighbors(grid, c, k)`, and `neighbors(grid, c, k)` is
`vcat(ring(grid, c, 1), ..., ring(grid, c, k))`.

# Connectivity

`Vertex()` (the default) is the HEALPix 3x3 lattice neighbourhood — 8 cells,
and 7 at the 24 pixels sitting on a degree-3 vertex of the base tiling.
`Edge()` keeps only the four that share a whole pixel edge. Because a HEALPix
pixel is a diamond on the lattice, the edge-sharing neighbours are the compass
directions SW, NW, NE, SE and the four corner-only ones are W, N, E, S; see
`neighbors.jl`.

# Order

**Counter-clockwise seen from outside, starting at the `SW` lattice
direction** — the compass cycle `SW, S, SE, E, NE, N, NW, W` under `Vertex()`
and its restriction `SW, SE, NE, NW` under `Edge()`. (The raw compass tuple in
`nested_neighbors` runs the other way; see `_neighbor_cycle`.) Absent
neighbours drop out of the cycle rather than leaving a hole, so a ring of seven
is still in rotational order.

Rings beyond the first are ordered by azimuth about the cell centre from the
same starting reference; see [`ring`](@ref).

`k == 0` returns an empty container; `k == 1` returns a
`SmallCollections.SmallVector` sized by `max_neighbors` and allocates nothing.
"""
function DGG.neighbors(g::HEALPixGrid, c::DGG.LevelIndex, k::Integer = 1;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return SmallVector{8,DGG.LevelIndex}()
    steps == 1 && return _one_ring(g, c, connectivity)
    shells = _shells(g, c, steps, connectivity)
    isempty(shells) && return DGG.LevelIndex[]
    return reduce(vcat, shells)
end

"""
    ring(grid, c, k; connectivity = Vertex())

The cells at lattice distance **exactly** `k`, counter-clockwise seen from
outside the sphere. `ring(grid, c, 0)` is `[c]`.

`k == 1` is the lattice compass cycle (see [`neighbors`](@ref)). For `k >= 2`
there is no lattice cycle to read off — an outer ring crosses face seams
arbitrarily — so the shell is ordered **by azimuth about the cell centre**,
measured counter-clockwise from the direction of the cell's own west corner.
That reference is chosen because it reproduces the `k == 1` cycle exactly: the
west corner sits at 135° in the lattice plane, so the first neighbour
counter-clockwise from it is `SW` at 180°. Ties in azimuth break by canonical
id, so the order is total and deterministic.
"""
function DGG.ring(g::HEALPixGrid, c::DGG.LevelIndex, k::Integer;
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
function _shells(g::HEALPixGrid, c::DGG.LevelIndex, steps::Int,
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
        # reach them, so they are wound geometrically.
        j > 1 && _sort_ccw!(next, g, c)
        push!(shells, next)
        isempty(next) && break
        frontier = next
    end
    return shells
end

# Order a shell counter-clockwise about `c`'s centre, from the azimuth of `c`'s
# own west corner. `e1 x e2 == centre` makes `(e1, e2)` right-handed SEEN FROM
# OUTSIDE, which is what puts increasing `atan(u.e2, u.e1)` counter-clockwise
# from outside rather than from inside.
function _sort_ccw!(cells::Vector{DGG.LevelIndex}, g::HEALPixGrid, c::DGG.LevelIndex)
    length(cells) <= 1 && return cells
    centre = DGG.cell_centroid(g, c)
    ring = DGG.cell_boundary(g, c)
    west = ring[1 + BOUNDARY_SEGMENTS]          # the ring's second corner
    e1, e2 = _tangent_basis(centre, west)
    ref = _azimuth(centre, e1, e2, west)
    key(x) = begin
        p = DGG.cell_centroid(g, x)
        (mod(_azimuth(centre, e1, e2, p) - ref, 2 * Float64(π)), x)
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
    # A degenerate reference (a cell whose west corner is at its centre) cannot
    # happen for a real pixel, but a fallback keeps this total.
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
