# test/fallbacks/runtests.jl — the generic substrate (src/fallbacks/)
#
# Three self-contained mocks, chosen so that every path through the substrate
# is reachable without a real grid system:
#
#   * `SortedMock` — an equirectangular quadtree with `has_sorted_subtrees`,
#     so the cursor descends in WINDOW mode and `descendant_range` is live.
#     Children are geometrically nested in their parents, which makes the
#     covering law easy to assert against.
#   * `UnsortedMock` — the same geometry with the trait off and no
#     `descendant_range` method at all, so the cursor must take SELECTION
#     mode. Every traversal answer must match `SortedMock`'s exactly; that
#     equality verifies both cursor modes.
#   * `OctantGrid` — eight spherical triangles, no system at all: the
#     standalone-grid path (position-space tree, geometric neighbours) with
#     areas that are known in closed form (pi/2 steradians each).
#
# The mocks cover a band of the sphere (latitude -45..45), not the whole of it,
# because partial coverage is a first-class property of the base interface and
# `cellat` outside the coverage has to answer `nothing` rather than guess.

module TestFallbacks

using Test
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
const FB = DGG.Fallbacks
import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
import Extents
import ConservativeRegridding
import ConservativeRegridding: Trees
import GeometryOps: SpatialTreeInterface as STI

const US = GO.UnitSpherical
const P = GO.UnitSphericalPoint

# ===========================================================================
# Mock hierarchical system: an equirectangular quadtree
#
# Level 0 is 4 longitude sectors x 2 latitude bands (8 roots) covering
# latitude -45..45; each cell splits into 2x2, with child `4id + k` taking the
# `(k & 1, k >> 1)` half in (lon, lat). Ids are `LevelIndex(level, id)` with a
# 0-based dense index, so a position in a level grid is `id + 1`.
# ===========================================================================

struct SortedMock <: AbstractHierarchicalGridSystem end
struct UnsortedMock <: AbstractHierarchicalGridSystem end

# The same hierarchy again, with one geometric difference: every cell's box is
# grown by `0.01 * (level + 1)` degrees, so a child sticks out past its parent's
# edge by 0.01 degrees. That is the aperture-7 pathology in miniature — children
# are NOT contained in their parents — and it is what makes `node_extent` (which
# still covers, with room to spare: see the covering-law testset) the only sound
# thing to prune a subtree with.
struct OverhangMock <: AbstractHierarchicalGridSystem end

const MockSystem = Union{SortedMock,UnsortedMock,OverhangMock}

const ROOTS = 8
const RADIX = 4
const MAXLEVEL = 6

mock_ncells(l::Int) = ROOTS * RADIX^l

DGG.cellindextype(::MockSystem) = LevelIndex
DGG.levels(::MockSystem) = 0:MAXLEVEL
DGG.max_neighbors(::MockSystem, ::Connectivity) = 8
DGG.has_sorted_subtrees(::SortedMock) = true
DGG.has_sorted_subtrees(::OverhangMock) = true

struct MockGrid{S<:MockSystem} <: AbstractGrid
    system::S
    level::Int
end

function DGG.levelgrid(sys::MockSystem, l::Integer)
    Int(l) in 0:MAXLEVEL || throw(ArgumentError("level $l is outside 0:$MAXLEVEL"))
    return MockGrid(sys, Int(l))
end

DGG.system(g::MockGrid) = g.system
DGG.level(g::MockGrid) = g.level
DGG.ncells(g::MockGrid) = mock_ncells(g.level)

function DGG.cellindex(g::MockGrid, i::Int)
    1 <= i <= mock_ncells(g.level) || throw(BoundsError(g, i))
    return LevelIndex(g.level, i - 1)
end

DGG.cellposition(g::MockGrid, c::LevelIndex) =
    (level(c) == g.level && 0 <= rawid(c) < mock_ncells(g.level)) ? Int(rawid(c)) + 1 : nothing

DGG.rootcells(::MockSystem) = [LevelIndex(0, i) for i in 0:(ROOTS-1)]

function Base.parent(::MockSystem, c::LevelIndex)
    level(c) > 0 || throw(ArgumentError("the root cell $c has no parent"))
    return LevelIndex(level(c) - 1, rawid(c) ÷ RADIX)
end

function DGG.children(::MockSystem, c::LevelIndex)
    level(c) < MAXLEVEL || throw(ArgumentError("$c is at max_level"))
    return [LevelIndex(level(c) + 1, rawid(c) * RADIX + k) for k in 0:(RADIX-1)]
end

# `descendant_range` returns one-based positions, not raw ids.
function DGG.descendant_range(::Union{SortedMock,OverhangMock}, c::LevelIndex, l::Integer)
    Int(l) >= level(c) || throw(ArgumentError("level $l is above the cell's own"))
    span = RADIX^(Int(l) - level(c))
    lo = Int(rawid(c)) * span
    return (lo+1):(lo+span)
end

# A second naming scheme, wired onto `UnsortedMock` ALONE so that `SortedMock`
# stays the single-scheme mock the assertions above it are written against.
# `MockPairIndex` names a cell by (root, offset within that root's subtree)
# rather than by a dense id — and, crucially for the `_canonical` regression
# tests, its converter REJECTS an offset that names no cell. A system converter
# is allowed to throw for a value that is simply not a cell; `cellposition` is
# not allowed to let that throw escape.
struct MockPairIndex <: AbstractCellIndex
    level::Int
    root::Int
    within::Int
end
DGG.level(c::MockPairIndex) = c.level
DGG.rawid(c::MockPairIndex) = c.root * RADIX^c.level + c.within
Base.isless(a::MockPairIndex, b::MockPairIndex) = isless(rawid(a), rawid(b))

DGG.cellindextypes(::UnsortedMock) = (LevelIndex, MockPairIndex)

function DGG.reindex(::Type{LevelIndex}, ::UnsortedMock, c::MockPairIndex)
    span = RADIX^c.level
    0 <= c.root < ROOTS ||
        throw(ArgumentError("root $(c.root) is outside 0:$(ROOTS - 1)"))
    0 <= c.within < span ||
        throw(ArgumentError("offset $(c.within) is outside 0:$(span - 1)"))
    return LevelIndex(c.level, c.root * span + c.within)
end

function DGG.reindex(::Type{MockPairIndex}, ::UnsortedMock, c::LevelIndex)
    root, within = divrem(Int(rawid(c)), RADIX^level(c))
    return MockPairIndex(level(c), root, within)
end

"Longitude/latitude box (degrees) of a mock cell."
function mock_box(c::LevelIndex)
    l = level(c)
    span = RADIX^l
    root, within = divrem(Int(rawid(c)), span)
    x = 0
    y = 0
    for bit in 0:(l-1)
        x |= ((within >> (2bit)) & 1) << bit
        y |= ((within >> (2bit + 1)) & 1) << bit
    end
    w = 90.0 / (1 << l)
    h = 45.0 / (1 << l)
    lon = -180.0 + 90.0 * (root % 4) + x * w
    lat = -45.0 + 45.0 * (root ÷ 4) + y * h
    return (lon, lat, lon + w, lat + h)
end

sph(lon, lat) = P(cosd(lat) * cosd(lon), cosd(lat) * sind(lon), sind(lat))

function DGG.cell_boundary(::MockGrid, c::LevelIndex)
    lon0, lat0, lon1, lat1 = mock_box(c)
    return [sph(lon0, lat0), sph(lon1, lat0), sph(lon1, lat1), sph(lon0, lat1)]
end

function DGG.cell_centroid(::MockGrid, c::LevelIndex)
    lon0, lat0, lon1, lat1 = mock_box(c)
    return sph((lon0 + lon1) / 2, (lat0 + lat1) / 2)
end

# Every cell's box grows by `OVERHANG_STEP * (level + 1)` degrees, so a child's
# outer edge lands `OVERHANG_STEP` degrees beyond its parent's: the children are
# strictly not contained in the parent, and a target sitting in that sliver
# meets the child while missing the parent entirely.
#
# The growth is tiny next to the default `node_extent`'s 1.2x cap inflation
# (checked in the covering-law testset), so the extents still cover — which is
# exactly the situation the substrate has to survive.
const OVERHANG_STEP = 0.01

overhang_pad(l::Int) = OVERHANG_STEP * (l + 1)

function DGG.cell_boundary(g::MockGrid{OverhangMock}, c::LevelIndex)
    lon0, lat0, lon1, lat1 = mock_box(c)
    pad = overhang_pad(level(c))
    lon0 -= pad; lat0 -= pad; lon1 += pad; lat1 += pad
    return [sph(lon0, lat0), sph(lon1, lat0), sph(lon1, lat1), sph(lon0, lat1)]
end

const SORTED = SortedMock()
const UNSORTED = UnsortedMock()
const OVERHANG = OverhangMock()

# ===========================================================================
# Mock standalone grid: the eight octants of the sphere
# ===========================================================================

const AXES = (P(1.0, 0.0, 0.0), P(0.0, 1.0, 0.0), P(0.0, 0.0, 1.0))

struct OctantIndex <: AbstractCellIndex
    code::Int          # 0:7, bits = (sign x, sign y, sign z), 0 = positive
end
DGG.level(c::OctantIndex) = 0
DGG.rawid(c::OctantIndex) = c.code
Base.isless(a::OctantIndex, b::OctantIndex) = isless(a.code, b.code)

struct OctantGrid <: AbstractGrid end

DGG.ncells(::OctantGrid) = 8
DGG.cellindex(::OctantGrid, i::Int) = (1 <= i <= 8 || throw(BoundsError(OctantGrid(), i));
OctantIndex(i - 1))
DGG.cellposition(::OctantGrid, c::OctantIndex) = c.code + 1

octant_signs(c::OctantIndex) = ((c.code & 1) == 0 ? 1.0 : -1.0,
    (c.code & 2) == 0 ? 1.0 : -1.0,
    (c.code & 4) == 0 ? 1.0 : -1.0)

function DGG.cell_boundary(::OctantGrid, c::OctantIndex)
    sx, sy, sz = octant_signs(c)
    v = [P(sx, 0.0, 0.0), P(0.0, sy, 0.0), P(0.0, 0.0, sz)]
    # Counter-clockwise seen from outside: flip when the signs multiply to -1.
    return sx * sy * sz > 0 ? v : reverse(v)
end

function DGG.cell_centroid(::OctantGrid, c::OctantIndex)
    sx, sy, sz = octant_signs(c)
    return P(sx, sy, sz) / sqrt(3)
end

# ===========================================================================
# A two-cell standalone grid for the degenerate `cell_extent` branches
# ===========================================================================

struct ExtentIndex <: AbstractCellIndex
    i::Int
end
DGG.level(c::ExtentIndex) = 0
DGG.rawid(c::ExtentIndex) = c.i
Base.isless(a::ExtentIndex, b::ExtentIndex) = isless(a.i, b.i)

struct ExtentGrid <: AbstractGrid end
DGG.ncells(::ExtentGrid) = 2
DGG.cellindex(::ExtentGrid, i::Int) = ExtentIndex(i)

function DGG.cell_boundary(::ExtentGrid, c::ExtentIndex)
    if c.i == 1                    # a polar cap: a ring at 60N, pole enclosed
        return [sph(lon, 60.0) for lon in 0.0:30.0:330.0]
    end                            # a box straddling the antimeridian
    return [sph(170.0, 10.0), sph(-170.0, 10.0), sph(-170.0, 20.0), sph(170.0, 20.0)]
end

DGG.cell_centroid(::ExtentGrid, c::ExtentIndex) =
    c.i == 1 ? sph(0.0, 90.0) : sph(180.0, 15.0)

# ===========================================================================
# Helpers
# ===========================================================================

all_cells(grid) = [cellindex(grid, i) for i in 1:ncells(grid)]

lonlat_ring(points) = GI.Polygon([GI.LinearRing([points..., points[1]])])

"""
Every cell of `grid` whose polygon relates to `geom` under `predicate` — the
oracle the tree descent is checked against, with no pruning and no sandwich.

`predicate` is the `pred_*` **constructor**, not an instance: GeometryOps'
predicates are mutable accumulators, so reusing one across cells silently
corrupts every answer after the first.
"""
function brute_force(grid, geom, predicate=GO.pred_intersects)
    prepared = GO.prepare(GO.RelateNG(; manifold=GO.Spherical()), geom)
    return [c for c in all_cells(grid)
            if GO.relate_predicate(prepared, predicate(), cell_polygon(grid, c))]
end

"Walk every node of a spatial tree, applying `f(node)`."
function walk(f, node)
    f(node)
    STI.isleaf(node) && return nothing
    for child in STI.getchild(node)
        walk(f, child)
    end
    return nothing
end

leaf_positions(tree) = STI.query(tree, _ -> true)

"`q`'s offset from `p`, projected into the tangent plane at `p`."
function tangent_offset(p, q)
    u = (q[1] - p[1], q[2] - p[2], q[3] - p[3])
    r = u[1] * p[1] + u[2] * p[2] + u[3] * p[3]
    return (u[1] - r * p[1], u[2] - r * p[2], u[3] - r * p[3])
end

"""
The signed volume `dot(cross(a, b), p)` of the tangent-plane offsets `a`, `b` of
`qa` and `qb` about `p`. Positive exactly when the step from `qa` to `qb` is
counter-clockwise **seen from outside** the sphere, which is the winding
[`neighbors`](@ref) promises. Built from the offsets directly rather than from
any azimuth, so it is an independent check on the implementation.
"""
function turn_sign(p, qa, qb)
    a = tangent_offset(p, qa)
    b = tangent_offset(p, qb)
    c = (a[2] * b[3] - a[3] * b[2],
        a[3] * b[1] - a[1] * b[3],
        a[1] * b[2] - a[2] * b[1])
    return c[1] * p[1] + c[2] * p[2] + c[3] * p[3]
end

"""
The counter-clockwise-seen-from-outside angle in `[0, 2pi)` from `qa` to `qb`
about `p`. Same independence as [`turn_sign`](@ref): no tangent basis, just the
dot and signed cross products of the two offsets.
"""
function ccw_angle(p, qa, qb)
    a = tangent_offset(p, qa)
    b = tangent_offset(p, qb)
    c = (a[2] * b[3] - a[3] * b[2],
        a[3] * b[1] - a[1] * b[3],
        a[1] * b[2] - a[2] * b[1])
    sine = c[1] * p[1] + c[2] * p[2] + c[3] * p[3]
    cosine = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
    return mod(atan(sine, cosine), 2 * Float64(pi))
end

# ===========================================================================

@testset "geometry generics" begin
    grid = levelgrid(SORTED, 2)
    c = cellindex(grid, 17)

    ring = cell_boundary(grid, c)
    @test length(ring) == 4
    @test all(p -> isapprox(sum(abs2, p), 1.0; atol=1e-12), ring)

    poly = cell_polygon(grid, c)
    @test GI.trait(poly) isa GI.PolygonTrait
    closed = collect(GI.getpoint(GI.getexterior(poly)))
    @test length(closed) == 5
    @test closed[1] == closed[end]
    @test closed[1:4] == ring

    # `getcell` is `cell_polygon . cellindex`, by position.
    @test collect(GI.getpoint(GI.getexterior(getcell(grid, 17)))) == closed
    @test getcell === Trees.getcell

    # Areas: an octant is exactly pi/2 steradians, and eight of them are the
    # sphere. This is the closed-form check on the spherical area path.
    octants = OctantGrid()
    for c in all_cells(octants)
        @test cell_area(octants, c) ≈ pi / 2 rtol = 1e-12
    end
    @test sum(cell_area(octants, c) for c in all_cells(octants)) ≈ 4pi rtol = 1e-12

    # The mock covers the band |lat| <= 45, whose exact area is
    # 2*pi*(sin 45 - sin -45). Its cells' top and bottom edges are great-circle
    # arcs rather than parallels, and a long arc bulges further poleward than
    # the two short arcs that replace it one level down — so a level's total
    # area decreases monotonically towards the band as the level refines, and
    # is never below it. That is the closed-form check on `cell_area` for
    # cells this package's own geometry produces.
    total(l) = sum(cell_area(levelgrid(SORTED, l), c)
                   for c in all_cells(levelgrid(SORTED, l)))
    band = 2pi * 2 * sind(45)
    @test total(0) > total(2) > total(4) > band
    @test total(4) ≈ band rtol = 1e-2
    @test all(c -> 0 < cell_area(grid, c) < 4pi, all_cells(grid))
end

@testset "cell_extent" begin
    grid = levelgrid(SORTED, 0)
    # Root cell 0 is lon -180..-90, lat -45..0. Its top edge is a great-circle
    # arc between two points at lat 0, which is the equator, so no bulge; the
    # bottom edge bulges *south*, past -45.
    ext = cell_extent(grid, cellindex(grid, 1))
    @test ext.X == (-180.0, -90.0)
    @test ext.Y[2] ≈ 0.0 atol = 1e-12
    @test ext.Y[1] < -45.0          # the great-circle bottom edge dips below

    # Northern band: the top edge bulges north past 45.
    north = cell_extent(grid, cellindex(grid, 5))
    @test north.Y[2] > 45.0
    @test north.Y[1] ≈ 0.0 atol = 1e-12

    degenerate = ExtentGrid()
    polar = cell_extent(degenerate, cellindex(degenerate, 1))
    @test polar.X == (-180.0, 180.0)     # a pole has no longitude
    @test polar.Y[2] == 90.0
    @test polar.Y[1] ≈ 60.0 rtol = 1e-6

    dateline = cell_extent(degenerate, cellindex(degenerate, 2))
    @test dateline.X == (-180.0, 180.0)  # antimeridian crossing is conservative
    @test dateline.Y[1] ≈ 10.0 rtol = 1e-3
    @test dateline.Y[2] > 20.0           # great-circle bulge again
end

@testset "identity generics" begin
    grid = levelgrid(SORTED, 3)
    @test cellindextypes(grid) == (LevelIndex,)
    @test cellindextypes(SORTED) == (LevelIndex,)
    @test cellindex(grid, 5, LevelIndex) === cellindex(grid, 5)
    @test_throws ArgumentError reindex(OctantIndex, SORTED, cellindex(grid, 1))

    # Bijection over the whole level.
    for i in (1, 2, 100, ncells(grid))
        @test cellposition(grid, cellindex(grid, i)) == i
    end
    # A cell from another level is simply not in this grid.
    @test cellposition(grid, LevelIndex(2, 0)) === nothing

    # The generic linear-scan `cellposition` on a grid with no override.
    octants = OctantGrid()
    @test cellindextypes(octants) == (OctantIndex,)
    for i in 1:8
        @test cellposition(octants, cellindex(octants, i)) == i
    end

    # ancestor / descendants
    c = LevelIndex(3, 100)
    @test ancestor(SORTED, c, 3) === c
    @test ancestor(SORTED, c, 2) === parent(SORTED, c)
    @test ancestor(SORTED, c, 0) === LevelIndex(0, 100 ÷ 64)
    @test_throws ArgumentError ancestor(SORTED, c, 4)

    @test descendants(SORTED, LevelIndex(0, 3), 0) == [LevelIndex(0, 3)]
    @test descendants(SORTED, LevelIndex(0, 3), 2) ==
          [LevelIndex(2, i) for i in 48:63]
    # The unsorted mock has no `descendant_range`, so this is the
    # children-expansion path — and it must agree.
    @test descendants(UNSORTED, LevelIndex(0, 3), 2) ==
          descendants(SORTED, LevelIndex(0, 3), 2)
    @test_throws ArgumentError descendants(SORTED, LevelIndex(2, 0), 1)

    # `descendant_range` is a MethodError for a system without the trait — the
    # trait and the method are declared together or not at all.
    @test_throws MethodError descendant_range(UNSORTED, LevelIndex(0, 0), 2)
end

@testset "cellposition: `nothing`, never a throw" begin
    # `cellposition(grid, c)` returns `nothing` for unsupported ids even when
    # the underlying `reindex` call throws an `ArgumentError`.

    # --- (a) a FOREIGN id type ------------------------------------------------
    # Taken at the grid's OWN level, so the level guard is not what saves it:
    # `OctantIndex` is level 0, and this is the level-0 grid.
    roots = levelgrid(SORTED, 0)
    @test level(OctantIndex(3)) == level(roots) == 0
    @test !(OctantIndex in cellindextypes(SORTED))
    @test cellposition(roots, OctantIndex(3)) === nothing
    # ... and on a PartialGrid, which binary-searches but routes through the
    # same `_canonical`.
    part0 = PartialGrid(SORTED, 0, [LevelIndex(0, i) for i in (1, 3, 5)])
    @test cellposition(part0, OctantIndex(3)) === nothing
    # ... and on a standalone grid, which has no system to ask at all.
    @test cellposition(OctantGrid(), LevelIndex(0, 1)) === nothing

    # --- (b) an OUT-OF-RANGE id in a scheme the system DOES support ------------
    # `UnsortedMock` names cells two ways, and the second scheme's converter
    # rejects an offset that is not a cell. That rejection is correct; letting
    # it out of `cellposition` is not.
    @test cellindextypes(UNSORTED) == (LevelIndex, MockPairIndex)
    @test cellindextype(UNSORTED) === LevelIndex
    @test reindex(LevelIndex, UNSORTED, MockPairIndex(2, 0, 3)) === LevelIndex(2, 3)
    @test reindex(MockPairIndex, UNSORTED, LevelIndex(2, 40)) === MockPairIndex(2, 2, 8)

    ids = [LevelIndex(2, i) for i in (3, 5, 8, 40, 41, 100)]
    partial = PartialGrid(UNSORTED, 2, ids)
    @test cellposition(partial, MockPairIndex(2, 0, 3)) == 1
    @test cellposition(partial, MockPairIndex(2, 0, 5)) == 2
    @test cellposition(partial, MockPairIndex(2, 2, 8)) == 4      # 2*16 + 8 = 40
    # In range, correctly converted, simply not one of this grid's cells.
    @test cellposition(partial, MockPairIndex(2, 0, 4)) === nothing

    # The regression, both halves: the converter throws, `cellposition` does not.
    @test_throws ArgumentError reindex(LevelIndex, UNSORTED, MockPairIndex(2, 0, 99))
    @test_throws ArgumentError reindex(LevelIndex, UNSORTED, MockPairIndex(2, 99, 0))
    @test cellposition(partial, MockPairIndex(2, 0, 99)) === nothing
    @test cellposition(partial, MockPairIndex(2, 99, 0)) === nothing
    # ... on the complete level grid too, through the generic linear scan.
    @test cellposition(levelgrid(UNSORTED, 2), MockPairIndex(2, 0, 5)) == 6
    @test cellposition(levelgrid(UNSORTED, 2), MockPairIndex(2, 0, 99)) === nothing

    # --- and `reindex` itself still throws ------------------------------------
    # The fix belongs at the `cellposition` boundary; pushing it down into
    # `reindex` would turn "I cannot name this" into a silent `nothing`.
    @test_throws ArgumentError reindex(OctantIndex, UNSORTED, LevelIndex(2, 0))
    @test_throws ArgumentError reindex(LevelIndex, UNSORTED, OctantIndex(1))
end

@testset "default node_extent and the covering law" begin
    for l in 0:2, i in 1:mock_ncells(l)
        c = LevelIndex(l, i - 1)
        extent = node_extent(SORTED, c)
        @test extent isa US.SphericalCap
        # Every descendant's every vertex, two levels down.
        for deeper in l:min(l + 2, MAXLEVEL)
            grid = levelgrid(SORTED, deeper)
            for d in descendants(SORTED, c, deeper)
                @test all(p -> US._contains(extent, p), cell_boundary(grid, d))
            end
        end
    end
    # The default is the cell's own cap inflated by `cap_inflation`, so it is
    # strictly bigger than the tight cap the cursor uses for leaves.
    grid = levelgrid(SORTED, 2)
    c = cellindex(grid, 7)
    @test node_extent(SORTED, c).radius > FB.cell_cap(grid, c).radius
end

@testset "cursor: window mode" begin
    grid = levelgrid(SORTED, 3)
    tree = treeify(grid)
    @test tree isa HierarchicalGridCursor
    @test STI.isspatialtree(typeof(tree))
    @test STI.node_extent_is_expensive(typeof(tree))
    @test Trees.ncells(tree) == ncells(grid)
    @test GOCore.best_manifold(tree) == GO.Spherical(; radius=1.0)

    # The tree's leaves are exactly the grid's positions, once each.
    @test leaf_positions(tree) == collect(1:ncells(grid))

    # Constant type through descent, and windows that partition their parent.
    nodes = 0
    walk(tree) do node
        nodes += 1
        @test typeof(node) === typeof(tree)
        STI.isleaf(node) && return
        window = FB.node_indices(node)
        covered = Int[]
        for child in STI.getchild(node)
            append!(covered, FB.node_indices(child))
        end
        @test sort(covered) == collect(window)
        @test STI.nchild(node) == length(collect(STI.getchild(node)))
    end
    @test nodes > ncells(grid)   # internal nodes exist

    # Node extents cover every leaf under the node — the property the whole
    # traversal rests on.
    walk(tree) do node
        extent = STI.node_extent(node)
        for i in FB.node_indices(node)
            for p in cell_boundary(grid, cellindex(grid, i))
                @test US._contains(extent, p)
            end
        end
    end

    # Leaf entries are (position, tight cap) pairs.
    leaf = tree
    while !STI.isleaf(leaf)
        leaf = first(STI.getchild(leaf))
    end
    entries = STI.child_indices_extents(leaf)
    @test length(entries) == 1
    index, cap = only(entries)
    @test index == first(FB.node_indices(leaf))
    @test cap == FB.cell_cap(grid, cellindex(grid, index))
    @test_throws ArgumentError STI.child_indices_extents(tree)

    # `Trees.getcell` addresses the same index space at the root.
    @test collect(GI.getpoint(GI.getexterior(Trees.getcell(tree, 5)))) ==
          collect(GI.getpoint(GI.getexterior(getcell(grid, 5))))
    @test_throws BoundsError Trees.getcell(tree, ncells(grid) + 1)
    # Parallelising is a question about chunk size, so the answer has to differ
    # between the whole grid and a single leaf however many threads are around.
    @test !Trees.should_parallelize(tree, FB.full_sphere_cap())
    @test Trees.should_parallelize(leaf, FB.full_sphere_cap())
end

@testset "cursor: selection mode" begin
    sorted = treeify(levelgrid(SORTED, 2))
    unsorted = treeify(levelgrid(UNSORTED, 2))
    @test unsorted isa HierarchicalGridCursor
    # Without `descendant_range` the cursor has to carry an explicit selection,
    # and at the root that selection is the whole grid, in position order.
    @test sorted.selection === nothing
    @test unsorted.selection == collect(1:ncells(levelgrid(UNSORTED, 2)))

    # The two modes must agree on everything a traversal can observe.
    @test leaf_positions(unsorted) == leaf_positions(sorted)
    walk(unsorted) do node
        @test typeof(node) === typeof(unsorted)
    end

    # Same node set, same windows, level by level.
    function levels_of(tree)
        out = Dict{Int,Vector{Vector{Int}}}()
        walk(tree) do node
            push!(get!(out, node.level, Vector{Int}[]), sort(collect(FB.node_indices(node))))
        end
        return Dict(k => sort(v) for (k, v) in out)
    end
    @test levels_of(unsorted) == levels_of(sorted)

    # The two modes must also agree on the answers a caller actually sees, not
    # only on the shape of the tree they descend. Both grids are the same
    # geometry, so every id must match exactly.
    gs = levelgrid(SORTED, 2)
    gu = levelgrid(UNSORTED, 2)
    for i in 1:ncells(gs)
        c = cellindex(gs, i)
        @test neighbors(gu, c) == neighbors(gs, c)
        @test neighbors(gu, c; connectivity=Edge()) ==
              neighbors(gs, c; connectivity=Edge())
    end
    for lon in -175.0:12.5:175.0, lat in -44.0:8.0:44.0
        @test cellat(gu, lon, lat) === cellat(gs, lon, lat)
    end
    # ... including where the answer is "no cell here".
    @test cellat(gu, 0.0, 70.0) === cellat(gs, 0.0, 70.0) === nothing
end

@testset "PartialGrid" begin
    ids = [LevelIndex(2, i) for i in (3, 5, 8, 40, 41, 100)]
    grid = PartialGrid(SORTED, 2, ids)

    @test grid isa AbstractGrid
    @test ncells(grid) == 6
    @test system(grid) === SORTED
    @test level(grid) == 2
    @test cellindex(grid, 3) === LevelIndex(2, 8)
    @test cellposition(grid, LevelIndex(2, 8)) == 3
    @test cellposition(grid, LevelIndex(2, 9)) === nothing
    @test cellposition(grid, LevelIndex(3, 8)) === nothing
    @test cell_boundary(grid, ids[1]) == cell_boundary(levelgrid(SORTED, 2), ids[1])
    @test cell_centroid(grid, ids[1]) == cell_centroid(levelgrid(SORTED, 2), ids[1])
    @test occursin("PartialGrid", sprint(show, grid))

    # --- validation --------------------------------------------------------
    @test_throws ArgumentError PartialGrid(SORTED, 2, [LevelIndex(2, 5), LevelIndex(2, 3)])
    @test_throws ArgumentError PartialGrid(SORTED, 2, [LevelIndex(2, 5), LevelIndex(2, 5)])
    @test_throws ArgumentError PartialGrid(SORTED, 2, [OctantIndex(1)])
    @test_throws ArgumentError PartialGrid(SORTED, 2, [LevelIndex(3, 5)])
    @test_throws ArgumentError PartialGrid(SORTED, 9, ids)
    @test_throws ArgumentError PartialGrid(SORTED, 2, ids; bucket_size=-1)
    # ... and root membership, which is what keeps a chunk's descent windowed.
    @test_throws ArgumentError PartialGrid(SORTED, 2, ids; root=LevelIndex(0, 0))
    @test PartialGrid(SORTED, 2, [LevelIndex(2, i) for i in 0:3];
        root=LevelIndex(0, 0)) isa PartialGrid
    @test_throws ArgumentError PartialGrid(SORTED, 1, ids; root=LevelIndex(2, 0))

    # --- the subtree constructor: lazy where the system allows it ----------
    subtree = PartialGrid(SORTED, LevelIndex(1, 5), 4)
    @test ncells(subtree) == RADIX^3
    @test subtree.ids isa FB.SubtreeIds
    @test collect(subtree.ids) == descendants(SORTED, LevelIndex(1, 5), 4)
    @test cellindex(subtree, 1) === LevelIndex(4, 5 * 64)

    materialised = PartialGrid(UNSORTED, LevelIndex(1, 5), 4)
    @test materialised.ids isa Vector
    @test collect(materialised.ids) == collect(subtree.ids)
end

@testset "cursor over a PartialGrid" begin
    # A sparse selection: six cells scattered over a level-4 grid.
    ids = [LevelIndex(4, i) for i in (0, 1, 2, 300, 1000, 2047)]
    grid = PartialGrid(SORTED, 4, ids)
    tree = treeify(grid)
    @test tree isa HierarchicalGridCursor
    @test Trees.ncells(tree) == 6
    @test leaf_positions(tree) == collect(1:6)

    # Every node's window is a slice of the id vector, and children partition it.
    walk(tree) do node
        STI.isleaf(node) && return
        covered = Int[]
        for child in STI.getchild(node)
            append!(covered, FB.node_indices(child))
        end
        @test sort(covered) == collect(FB.node_indices(node))
    end

    # A sparse internal node tightens its extent past the system's own
    # `node_extent`, which is the one thing the cursor adds to the hierarchy.
    tightened = false
    walk(tree) do node
        STI.isleaf(node) && return
        node.level < 0 && return
        count = length(FB.node_indices(node))
        if 0 < count < RADIX^(4 - node.level)
            if STI.node_extent(node).radius < node_extent(SORTED, node.id).radius
                tightened = true
            end
        end
    end
    @test tightened

    # ... and it still covers what it claims to.
    walk(tree) do node
        extent = STI.node_extent(node)
        for i in FB.node_indices(node)
            for p in cell_boundary(grid, cellindex(grid, i))
                @test US._contains(extent, p)
            end
        end
    end

    # A rooted chunk starts descent at its own root, not at the sphere.
    chunk = PartialGrid(SORTED, LevelIndex(1, 5), 4)
    rooted = treeify(chunk)
    @test rooted.level == 1
    @test rooted.id === LevelIndex(1, 5)
    @test leaf_positions(rooted) == collect(1:ncells(chunk))

    # Bucketed descent stops early and scans.
    bucketed = treeify(PartialGrid(SORTED, LevelIndex(1, 5), 4; bucket_size=16))
    @test leaf_positions(bucketed) == collect(1:ncells(chunk))
    leaves = 0
    walk(bucketed) do node
        STI.isleaf(node) && (leaves += 1)
    end
    @test leaves == ncells(chunk) ÷ 16

    # The selection-mode cursor over the same subset agrees.
    unsorted_ids = [LevelIndex(4, i) for i in (0, 1, 2, 300, 1000, 2047)]
    unsorted = treeify(PartialGrid(UNSORTED, 4, unsorted_ids))
    @test unsorted.selection isa Vector{Int}
    @test leaf_positions(unsorted) == collect(1:6)
end

@testset "position-space fallback tree" begin
    grid = OctantGrid()
    tree = treeify(grid)
    @test tree isa FB.PositionTreeNode
    @test STI.isspatialtree(typeof(tree))
    @test !STI.node_extent_is_expensive(typeof(tree))
    @test Trees.ncells(tree) == 8
    @test leaf_positions(tree) == collect(1:8)
    @test system(grid) === nothing && level(grid) === nothing

    # Node extents cover their leaves.
    walk(tree) do node
        extent = STI.node_extent(node)
        if STI.isleaf(node)
            for (index, _) in STI.child_indices_extents(node)
                for p in cell_boundary(grid, cellindex(grid, index))
                    @test US._contains(extent, p)
                end
            end
        end
    end

    # A bigger grid actually branches (8 cells fit in one leaf block).
    big = treeify(levelgrid(SORTED, 2))
    @test STI.nchild(big) > 0

    # `treeify` is idempotent and manifold-agnostic.
    @test treeify(GO.Spherical(), grid) isa FB.PositionTreeNode
    @test treeify(tree) === tree
    @test treeify(GO.Spherical(), tree) === tree
end

@testset "cellat" begin
    grid = levelgrid(SORTED, 2)
    for i in 1:ncells(grid)
        c = cellindex(grid, i)
        @test cellat(grid, cell_centroid(grid, c)) === c
    end
    # The lon/lat wrapper takes degrees.
    @test cellat(grid, -179.0, -44.0) === cellat(grid, sph(-179.0, -44.0))
    # Outside the mock's -45..45 band: a real `nothing`, not an error.
    @test cellat(grid, 0.0, 80.0) === nothing
    @test cellat(grid, 0.0, -80.0) === nothing

    # A standalone grid goes through the position tree.
    octants = OctantGrid()
    for c in all_cells(octants)
        @test cellat(octants, cell_centroid(octants, c)) === c
    end
    # A boundary point belongs to exactly one cell, deterministically the first
    # in canonical order among the cells that legitimately claim it. The
    # candidate set is recomputed here by exhaustive scan, so this pins the
    # tie-break rule rather than restating the call.
    corner = P(1.0, 0.0, 0.0)
    verdicts = [(c, FB.point_in_cell(cell_boundary(octants, c), corner))
                for c in sort(all_cells(octants))]
    claimants = [c for (c, v) in verdicts if v === true]
    undecided = [c for (c, v) in verdicts if v === nothing]
    @test length(claimants) + length(undecided) > 1     # genuinely shared
    onedge = cellat(octants, corner)
    @test onedge === (isempty(claimants) ? first(undecided) : first(claimants))
    @test onedge === OctantIndex(0)
end

@testset "geometric neighbors and ring" begin
    octants = OctantGrid()
    # An octant shares an edge with the three octants one sign flip away, and
    # exactly one vertex with the three two flips away; the antipode shares
    # nothing.
    for c in all_cells(octants)
        moore = neighbors(octants, c)
        vonneumann = neighbors(octants, c; connectivity=Edge())
        @test length(moore) == 6
        @test length(vonneumann) == 3
        # NOT `issorted`: the order is rotational now, and only the membership
        # is asserted here. The winding itself is pinned in the
        # "rotational neighbour order" testset below.
        @test allunique(moore)
        @test allunique(vonneumann)
        @test vonneumann ⊆ moore
        @test !(c in moore)
        @test all(d -> count_ones(xor(c.code, d.code)) == 1, vonneumann)
        @test OctantIndex(xor(c.code, 7)) ∉ moore
    end

    # Symmetry, both connectivities.
    for c in all_cells(octants), conn in (Vertex(), Edge())
        for d in neighbors(octants, c; connectivity=conn)
            @test c in neighbors(octants, d; connectivity=conn)
        end
    end

    # k = 0 is empty, rings shell out, and the disc is the union of the shells.
    c = OctantIndex(0)
    @test isempty(neighbors(octants, c, 0))
    @test ring(octants, c, 0) == [c]
    @test ring(octants, c, 1) == neighbors(octants, c, 1)
    @test ring(octants, c, 2) == [OctantIndex(7)]
    # The disc is the rings concatenated OUTWARD — no sort, or the shells would
    # interleave and `ring` would stop being the tail block.
    @test vcat(ring(octants, c, 1), ring(octants, c, 2)) ==
          neighbors(octants, c, 2)
    @test length(neighbors(octants, c, 3)) == 7
    @test_throws ArgumentError neighbors(octants, c, -1)
    @test_throws ArgumentError ring(octants, c, -1)

    # The mock grid: an interior cell has 8 Moore and 4 edge neighbours...
    grid = levelgrid(SORTED, 1)
    interior = LevelIndex(1, 2)            # lon -180..-135, lat -22.5..0
    @test length(neighbors(grid, interior; connectivity=Edge())) == 4
    @test length(neighbors(grid, interior)) == 8
    @test all(d -> d != interior, neighbors(grid, interior))
    # ... and a cell on the coverage edge simply has fewer: neighbours outside
    # the grid are ABSENT from the result, never padded.
    edge = LevelIndex(1, 5)                # lon -45..0, lat -45..-22.5
    @test length(neighbors(grid, edge; connectivity=Edge())) == 3
    @test length(neighbors(grid, edge)) == 5
    # Longitude wraps across the antimeridian even though the grid's own
    # coordinates jump from 180 to -180 there.
    @test LevelIndex(1, 15) in neighbors(grid, interior)   # lon 135..180
end

@testset "rotational neighbour order" begin
    # The order is contract, not convenience: it is what makes position `j` of
    # a ring name a fixed direction. These are the three laws that pin it.
    octants = OctantGrid()
    g1 = levelgrid(SORTED, 1)
    g2 = levelgrid(SORTED, 2)

    subjects = vcat([(octants, c) for c in all_cells(octants)],
        [(g1, c) for c in all_cells(g1)],
        [(g2, LevelIndex(2, i)) for i in (0, 21, 42, 100)])

    for (g, c) in subjects, conn in (Vertex(), Edge())
        # 1. Shell concatenation: the disc IS the rings, outward, never sorted.
        for k in 0:2
            disc = neighbors(g, c, k; connectivity=conn)
            rings = [ring(g, c, j; connectivity=conn) for j in 1:k]
            @test disc == reduce(vcat, rings; init=typeof(c)[])
            @test allunique(disc)
            @test !(c in disc)
        end
        # 2. The tail-block law, element for element — and its corollary, that
        #    the head of a disc is the next smaller disc, so a coarser
        #    neighbourhood extends a finer one by appending.
        for k in 1:2
            shell = ring(g, c, k; connectivity=conn)
            disc = neighbors(g, c, k; connectivity=conn)
            isempty(shell) && continue
            @test disc[(end-length(shell)+1):end] == shell
            @test disc[1:(end-length(shell))] ==
                  neighbors(g, c, k - 1; connectivity=conn)
        end
        # 3. Every shell is ordered by azimuth about the subject's centroid,
        #    measured from ONE spoke — ring 1's first cell — so ring 2 does not
        #    get a zero direction of its own. Recomputed here from the offsets
        #    rather than from the fallback's frame.
        r1 = ring(g, c, 1; connectivity=conn)
        if !isempty(r1)
            p = cell_centroid(g, c)
            spoke = cell_centroid(g, first(r1))
            @test ccw_angle(p, spoke, spoke) == 0.0
            for k in 1:2
                shell = ring(g, c, k; connectivity=conn)
                @test issorted([ccw_angle(p, spoke, cell_centroid(g, d))
                                for d in shell])
            end
        end
    end

    # 4. Ring 1 winds counter-clockwise SEEN FROM OUTSIDE, all the way round:
    #    the signed volume of consecutive tangent offsets about the subject's
    #    centroid is positive at every step, the wrap from last back to first
    #    included. Only a ring that actually encircles its subject can close, so
    #    the closure check is on cells with a full complement of neighbours —
    #    every octant, and the mock's interior cells. (A coverage-edge cell's
    #    ring is a partial arc; it still winds, but it does not come round.)
    encircling = vcat([(octants, c, conn) for c in all_cells(octants)
                       for conn in (Vertex(), Edge())],
        [(g1, c, conn) for c in all_cells(g1) for conn in (Vertex(), Edge())
         if length(neighbors(g1, c; connectivity=conn)) == (conn isa Edge ? 4 : 8)])
    @test count(t -> t[1] === g1, encircling) > 0     # the filter is not empty
    for (g, c, conn) in encircling
        r = ring(g, c, 1; connectivity=conn)
        @test length(r) >= 3
        p = cell_centroid(g, c)
        qs = [cell_centroid(g, d) for d in r]
        for i in eachindex(qs)
            @test turn_sign(p, qs[i], qs[i == length(qs) ? 1 : i+1]) > 0
        end
    end

    # ... and a partial arc still winds, step by step, even though it does not
    # close: this is the coverage-edge cell the testset above measured at 5.
    let c = LevelIndex(1, 5), p = cell_centroid(g1, c)
        r = ring(g1, c, 1)
        @test length(r) == 5
        qs = [cell_centroid(g1, d) for d in r]
        # Azimuths ascend from the spoke — the arc's own consecutive steps are
        # counter-clockwise wherever they turn by less than half a circle.
        @test any(i -> turn_sign(p, qs[i], qs[i+1]) > 0, 1:(length(qs)-1))
    end
end

@testset "query: correctness against brute force" begin
    grid = levelgrid(SORTED, 3)

    # A box, an antimeridian-crossing box, a thin sliver, and a triangle.
    targets = Dict(
        "box" => lonlat_ring([sph(10.0, 10.0), sph(40.0, 10.0), sph(40.0, 30.0), sph(10.0, 30.0)]),
        "antimeridian" => lonlat_ring([sph(170.0, -10.0), sph(-170.0, -10.0),
            sph(-170.0, 10.0), sph(170.0, 10.0)]),
        "sliver" => lonlat_ring([sph(-100.0, -5.0), sph(100.0, -5.0),
            sph(100.0, -4.0), sph(-100.0, -4.0)]),
        "triangle" => lonlat_ring([sph(-30.0, -20.0), sph(0.0, 20.0), sph(30.0, -20.0)]),
    )

    for (name, target) in targets
        expected = brute_force(grid, target)
        got = query(grid, Intersects(target))
        @test issorted(got)
        @test got == expected

        # `Disjoint` against its own oracle — the engine computes it as the
        # complement of `Intersects`, so checking it against that complement
        # would only restate the implementation.
        @test query(grid, Disjoint(target)) == brute_force(grid, target, GO.pred_disjoint)

        # `Within` and `Covers` against their own oracles.
        @test query(grid, Within(target)) ==
              brute_force(grid, target, GO.pred_contains)
        @test query(grid, Covers(target)) ==
              brute_force(grid, target, GO.pred_coveredby)
        @test query(grid, Touches(target)) ==
              brute_force(grid, target, GO.pred_touches)
        @test query(grid, Within(target)) ⊆ got
    end

    # The rim sandwich is an optimisation, so switching it off must not change
    # a single answer: this is the same query with the arcs suppressed.
    target = targets["triangle"]
    prepared = GO.prepare(GO.RelateNG(; manifold=GO.Spherical()), target)
    plain = FB.GeometryTarget(prepared, target, FB._geometry_cap(prepared, target), nothing)
    @test FB._run_query(grid, Intersects(target), plain) == query(grid, Intersects(target))

    # The system-level form answers at the requested level.
    @test query(SORTED, Intersects(targets["box"]); level=3) == query(grid, Intersects(targets["box"]))
    @test query(SORTED, Intersects(targets["box"]); level=1) ==
          query(levelgrid(SORTED, 1), Intersects(targets["box"]))

    # Selection mode gives the same answers as window mode.
    for (_, target) in targets
        @test query(levelgrid(UNSORTED, 2), Intersects(target)) ==
              query(levelgrid(SORTED, 2), Intersects(target))
    end

    # A lon/lat target (2-D coordinates, degrees) is lifted to the sphere at
    # the boundary of the call and must give exactly the same answer as the
    # same ring already expressed in unit-sphere xyz.
    lonlat_box = GI.Polygon([GI.LinearRing([(10.0, 10.0), (40.0, 10.0),
        (40.0, 30.0), (10.0, 30.0), (10.0, 10.0)])])
    @test query(grid, Intersects(lonlat_box)) == query(grid, Intersects(targets["box"]))

    # A standalone grid: the position tree, same oracle.
    octants = OctantGrid()
    small = lonlat_ring([sph(10.0, 10.0), sph(20.0, 10.0), sph(20.0, 20.0), sph(10.0, 20.0)])
    @test query(octants, Intersects(small)) == brute_force(octants, small)
    @test length(query(octants, Intersects(small))) == 1
end

@testset "query: predicate direction" begin
    # The engine prepares the TARGET and asks GeometryOps for the converse
    # relation, so a swapped mapping would be invisible to an oracle built the
    # same way. These oracles instead call the unprepared, direction-explicit
    # form with the cell first, exactly as `Within(target)` reads: "cell within
    # target".
    grid = levelgrid(SORTED, 3)
    target = lonlat_ring([sph(-30.0, -20.0), sph(30.0, -20.0), sph(30.0, 20.0),
        sph(-30.0, 20.0)])
    alg = GO.RelateNG(; manifold=GO.Spherical())

    for (predicate, relation) in ((Within, GO.within), (Contains, GO.contains),
        (Covers, GO.covers), (CoveredBy, GO.coveredby),
        (Intersects, GO.intersects), (Touches, GO.touches),
        (Overlaps, GO.overlaps), (Disjoint, GO.disjoint))
        expected = [c for c in all_cells(grid)
                    if relation(alg, cell_polygon(grid, c), target)]
        @test query(grid, predicate(target)) == expected
    end

    # ... and the asymmetric pair really is asymmetric on this target, so the
    # loop above has something to catch: cells are far smaller than the target.
    @test !isempty(query(grid, Within(target)))
    @test isempty(query(grid, Contains(target)))
    @test !isempty(query(grid, Overlaps(target)))

    # `Equals` needs a target that is one cell exactly, and then it must find
    # that cell and no other.
    twin = cell_polygon(grid, LevelIndex(3, 100))
    @test query(grid, Equals(twin)) == [LevelIndex(3, 100)]
end

@testset "query: targets and predicate errors" begin
    grid = levelgrid(SORTED, 2)
    box = lonlat_ring([sph(10.0, 10.0), sph(40.0, 10.0), sph(40.0, 30.0), sph(10.0, 30.0)])

    # An `Extents.Extent` in lon/lat degrees densifies to the same answer as
    # the equivalent polygon, up to the sag of the densified parallels.
    from_extent = query(grid, Intersects(Extents.Extent(X=(10.0, 40.0), Y=(10.0, 30.0))))
    @test from_extent ⊇ query(grid, Within(box))
    @test !isempty(from_extent)

    # A spherical cap target, answered exactly and without polygonising it.
    cap = US.SphericalCap(sph(20.0, 20.0), deg2rad(10.0))
    hits = query(grid, Intersects(cap))
    @test !isempty(hits)
    missed = query(grid, Disjoint(cap))
    near_checked = 0
    far_checked = 0
    for c in all_cells(grid)
        # Two one-sided oracles, both proofs. A cell with a point inside the
        # cap meets it; a cell whose every point is outside — bounded below by
        # the centroid distance less the circumradius — does not.
        centroid = cell_centroid(grid, c)
        boundary = cell_boundary(grid, c)
        if US.spherical_distance(cap.point, centroid) <= cap.radius ||
           any(p -> US.spherical_distance(cap.point, p) <= cap.radius, boundary)
            near_checked += 1
            @test c in hits
            @test !(c in missed)
        elseif US.spherical_distance(cap.point, centroid) >
               cap.radius + maximum(p -> US.spherical_distance(centroid, p), boundary)
            far_checked += 1
            @test !(c in hits)
            @test c in missed
        end
    end
    @test near_checked > 0 && far_checked > 0
    @test query(grid, Within(cap)) ⊆ hits
    @test_throws ArgumentError query(grid, Touches(cap))

    # Unimplemented predicate, unusable target, band extent.
    @test_throws ArgumentError query(grid, Crosses(box))
    @test_throws ArgumentError query(grid, Intersects(nothing))
    @test_throws ArgumentError query(grid, Intersects(Extents.Extent(X=(-180.0, 180.0), Y=(-10.0, 10.0))))
    # ... but a full-longitude extent reaching a pole IS a cap: cells that
    # reach past 40N are in the answer, cells that stay well south of it not.
    polar = query(grid, Intersects(Extents.Extent(X=(-180.0, 180.0), Y=(40.0, 90.0))))
    high = 0
    low = 0
    for c in all_cells(grid)
        lats = [asind(p[3]) for p in cell_boundary(grid, c)]
        if minimum(lats) > 41.0 || maximum(lats) > 41.0
            high += 1
            @test c in polar
        elseif maximum(lats) < 39.0
            low += 1
            @test !(c in polar)
        end
    end
    @test high > 0 && low > 0

    # An empty grid answers empty rather than erroring.
    empty_grid = PartialGrid(SORTED, 2, LevelIndex[])
    @test isempty(query(empty_grid, Intersects(box)))
    @test isempty(query(empty_grid, Disjoint(box)))
end

@testset "ConservativeRegridding end to end" begin
    # The whole `Trees` surface, exercised by the consumer it exists for:
    # `best_manifold` picks the manifold, `treeify` builds both cursors,
    # the dual depth-first search walks them against each other with
    # `node_extent`/`child_indices_extents`, and `getcell` resolves the pairs.
    src = levelgrid(SORTED, 2)
    dst = levelgrid(SORTED, 3)
    regridder = ConservativeRegridding.Regridder(dst, src)

    # A constant field is the sharp check: it comes back unchanged exactly when
    # the intersection weights of every destination cell sum to its own area —
    # i.e. when the dual tree walk found every overlapping source cell and got
    # each intersection area right. A single missed pair shows up immediately.
    out = zeros(ncells(dst))
    ConservativeRegridding.regrid!(out, regridder, fill(3.5, ncells(src)))
    @test all(v -> isapprox(v, 3.5; rtol=1e-10), out)

    # A varying field must stay a convex combination of its sources.
    ramp = Float64.(1:ncells(src))
    ConservativeRegridding.regrid!(out, regridder, ramp)
    @test minimum(out) >= minimum(ramp) - 1e-8
    @test maximum(out) <= maximum(ramp) + 1e-8

    # Most destination cells sit strictly inside one parent and take its value
    # exactly. The rest do not, and that is the mock's geometry rather than the
    # regridder's doing: a cell's edges are great-circle arcs, and the long arc
    # bounding a level-2 cell bulges further poleward than the two short arcs
    # bounding its children, so the lens between them belongs to a different
    # parent than to the children that tile it.
    exact = count(1:ncells(dst)) do i
        p = parent(SORTED, cellindex(dst, i))
        isapprox(out[i], ramp[cellposition(src, p)]; rtol=1e-8)
    end
    @test exact > ncells(dst) ÷ 2
end

@testset "MultiOrderCoverage" begin
    target = lonlat_ring([sph(-20.0, -20.0), sph(20.0, -20.0), sph(20.0, 20.0), sph(-20.0, 20.0)])
    coverage = MultiOrderCoverage(target)
    @test parent(coverage) === target

    set = query(SORTED, coverage; level=4)
    @test set isa MultiOrderCellSet
    @test !isempty(set)
    @test eltype(set) === LevelIndex
    @test occursin("MultiOrderCellSet", sprint(show, set))
    @test DGG.system(set) === SORTED

    cells = collect(set)
    # Indexing and iteration are separate methods; they must agree.
    @test all(set[i] === cells[i] for i in 1:length(set))
    @test allunique(cells)
    # Mixed levels, coarsest possible: at least one cell above the max depth.
    @test minimum(level, cells) < 4
    @test maximum(level, cells) <= 4

    # Curve order, checked against the intervals themselves rather than against
    # the stored keys: the cells' level-4 position intervals ascend and are
    # pairwise disjoint, which is both the curve order and the guarantee that no
    # cell is an ancestor of another (that would double-count a region).
    intervals = [descendant_range(SORTED, c, 4) for c in cells]
    @test issorted(intervals; by=first)
    for k in 2:length(intervals)
        @test first(intervals[k]) > last(intervals[k-1])
    end
    for a in cells, b in cells
        a === b && continue
        level(a) < level(b) && @test ancestor(SORTED, b, level(a)) !== a
    end

    # Per-level expansion: sorted, disjoint, merged position ranges.
    ranges = level_ranges(set, 4)
    @test issorted(first.(ranges))
    for k in 2:length(ranges)
        @test first(ranges[k]) > last(ranges[k-1]) + 1   # merged, so never adjacent
    end
    positions = reduce(vcat, collect.(ranges))
    @test allunique(positions)
    @test issorted(positions)
    # Merging changed the shape of the answer and nothing else.
    @test positions == sort!(reduce(vcat, collect.(intervals)))

    # The coverage brackets the two single-level queries at that depth: every
    # cell wholly inside the target is in it, and nothing outside the target is.
    grid = levelgrid(SORTED, 4)
    inside = query(grid, Within(target))
    touching = query(grid, Intersects(target))
    expanded = FB.cellindices(set, 4)
    @test inside ⊆ expanded
    @test expanded ⊆ touching

    # Expansion to a coarser level than the set's own cells is refused.
    @test_throws ArgumentError level_ranges(set, 1)
    # ... and a system without descendant ranges has no position ranges at all,
    # though it must still find the same cells: the two mocks are the same
    # hierarchy over the same geometry, so their coverages have to agree.
    unsorted_set = query(UNSORTED, coverage; level=2)
    @test !isempty(unsorted_set)
    @test_throws ArgumentError level_ranges(unsorted_set, 2)
    @test sort!(reduce(vcat, [descendants(UNSORTED, c, 2) for c in unsorted_set])) ==
          FB.cellindices(query(SORTED, coverage; level=2), 2)

    @test_throws ArgumentError query(SORTED, coverage; level=MAXLEVEL + 1)

    # A cap target works the same way. Oracle: a cell whose centroid is inside
    # the cap certainly meets it and so must be covered, and a cell further off
    # than the cap's radius plus its own circumradius certainly does not.
    cap = US.SphericalCap(sph(0.0, 0.0), deg2rad(15.0))
    cap_set = query(SORTED, MultiOrderCoverage(cap); level=3)
    @test !isempty(cap_set)
    covered = FB.cellindices(cap_set, 3)
    grid3 = levelgrid(SORTED, 3)
    checked_near = 0
    checked_far = 0
    for c in all_cells(grid3)
        centroid = cell_centroid(grid3, c)
        d = US.spherical_distance(centroid, cap.point)
        circum = maximum(p -> US.spherical_distance(centroid, p), cell_boundary(grid3, c))
        if d < cap.radius
            checked_near += 1
            @test c in covered
        elseif d > cap.radius + circum
            checked_far += 1
            @test !(c in covered)
        end
    end
    @test checked_near > 0 && checked_far > 0
end

@testset "query: caps wider than a hemisphere" begin
    # A cap of radius > pi/2 is not convex, so containment cannot be read off
    # the vertices; it is decided through the complement cap instead. "All of
    # the world except one region" is an ordinary enough target for this to be
    # worth getting right rather than refusing.
    grid = levelgrid(SORTED, 2)
    wide = US.SphericalCap(sph(0.0, 0.0), deg2rad(100.0))
    anti = sph(180.0, 0.0)

    # Oracle: sample every edge densely and ask whether all of it is inside.
    # Equivalent to containment here because a cell this size cannot swallow
    # the 80-degree complement cap.
    function boundary_inside(grid, c, cap)
        ring, n = FB.open_ring(cell_boundary(grid, c))
        for i in 1:n
            a, b = ring[i], ring[i == n ? 1 : i+1]
            for t in range(0.0, 1.0; length=33)
                x = (1 - t) * a[1] + t * b[1]
                y = (1 - t) * a[2] + t * b[2]
                z = (1 - t) * a[3] + t * b[3]
                s = sqrt(x * x + y * y + z * z)
                US.spherical_distance(cap.point, P(x / s, y / s, z / s)) <= cap.radius ||
                    return false
            end
        end
        return true
    end

    inside = query(grid, Within(wide))
    @test inside == [c for c in all_cells(grid) if boundary_inside(grid, c, wide)]
    @test !isempty(inside)                       # ... and the oracle is not trivial
    @test length(inside) < ncells(grid)
    @test inside ⊆ query(grid, Intersects(wide))

    # A cap covering the whole sphere holds every cell; the query engine must
    # not trip over the degenerate complement.
    @test query(grid, Within(FB.full_sphere_cap())) == all_cells(grid)
    @test query(grid, Within(US.SphericalCap(sph(0.0, 0.0), Float64(pi)))) == all_cells(grid)

    # Exercise both coverage paths: coarse cells emitted whole and boundary
    # cells reached by recursive descent.
    set = query(SORTED, MultiOrderCoverage(wide); level=3)
    cells = collect(set)
    @test minimum(level, cells) < 3              # emitted coarse
    @test maximum(level, cells) == 3             # recursed to the rim
    @test allunique(cells)

    covered = FB.cellindices(set, 3)
    deep = levelgrid(SORTED, 3)
    met = 0
    missed = 0
    for c in all_cells(deep)
        boundary = cell_boundary(deep, c)
        centroid = cell_centroid(deep, c)
        if any(p -> US.spherical_distance(wide.point, p) <= wide.radius, boundary)
            met += 1                              # a vertex inside: the cell meets it
            @test c in covered
        elseif US.spherical_distance(anti, centroid) +
               maximum(p -> US.spherical_distance(centroid, p), boundary) <
               Float64(pi) - wide.radius
            missed += 1                           # wholly inside the complement
            @test !(c in covered)
        end
    end
    @test met > 0 && missed > 0
end

@testset "MultiOrderCoverage: children that overhang their parents" begin
    # The regression test for the substrate's one real soundness trap. Under a
    # non-congruent refinement a cell can be disjoint from the target while one
    # of its children meets it, so the exact intersection test may decide what
    # is EMITTED and never what is DESCENDED into. Pruning a subtree on the
    # parent's own geometry drops that child silently — no error, just a
    # coverage with a hole in it.

    # First: the mock is legitimate. Its extents still cover every descendant...
    for l in 0:2, c in all_cells(levelgrid(OVERHANG, l))
        extent = node_extent(OVERHANG, c)
        for deeper in l:min(l + 2, MAXLEVEL)
            grid = levelgrid(OVERHANG, deeper)
            for d in descendants(OVERHANG, c, deeper)
                @test all(p -> US._contains(extent, p), cell_boundary(grid, d))
            end
        end
    end
    # ... and its children really do escape their parents.
    let p = LevelIndex(1, 2), g1 = levelgrid(OVERHANG, 1), g2 = levelgrid(OVERHANG, 2)
        outside = 0
        for child in children(OVERHANG, p), v in cell_boundary(g2, child)
            FB.point_in_cell(cell_boundary(g1, p), v) === false && (outside += 1)
        end
        @test outside > 0
    end

    # A sliver in the gap between the level-1 cell's east edge (-135 + 0.02) and
    # its children's (-135 + 0.03): inside two of the children, outside the
    # parent, and inside the neighbouring level-1 cell to the east.
    target = lonlat_ring([sph(-134.977, -12.0), sph(-134.973, -12.0),
        sph(-134.973, -11.0), sph(-134.977, -11.0)])

    grid1 = levelgrid(OVERHANG, 1)
    grid2 = levelgrid(OVERHANG, 2)
    met1 = query(grid1, Intersects(target))
    touching = query(grid2, Intersects(target))
    @test !isempty(touching)

    # The witness: level-2 cells that meet the target although their parents do
    # not. Without these the testset would pass for the wrong reason.
    witness = [c for c in touching if !(parent(OVERHANG, c) in met1)]
    @test !isempty(witness)
    @test LevelIndex(2, 9) in witness
    @test !(LevelIndex(1, 2) in met1)

    set = query(OVERHANG, MultiOrderCoverage(target); level=2)
    expanded = FB.cellindices(set, 2)
    # The coverage covers every cell that meets the target, overhang included.
    @test touching ⊆ expanded
    @test witness ⊆ expanded
end

@testset "generic subtree_border / subtree_interior" begin
    # Every shipped system specializes `rim_engine`, so this testset and A5's
    # are the only exercises the generic scan in `src/fallbacks/subtree.jl`
    # gets. Without it, the correct-for-everyone implementation a new system
    # inherits before it writes an automaton is dead code.
    #
    # The mock's subtree at depth `d` is a `2^d x 2^d` block of the lon/lat
    # lattice, so its rim is contained in the block's boundary ring — at most
    # `4 * 2^d - 4` cells. Not exactly that many: the mock's roots tile a
    # bounded lon/lat domain rather than a closed surface, so a block on the
    # domain's outer edge has cells with nothing beyond them, and those are not
    # exposed. The bound is still an independent geometric fact (an interior
    # lattice cell cannot be on the rim), and the definitional loops below pin
    # the exact answer.
    for sys in (SORTED, UNSORTED)
        root = LevelIndex(0, 3)

        # A depth-0 subtree is the cell itself, and it is all rim.
        @test subtree_border(sys, root, 0) == [root]
        @test isempty(subtree_interior(sys, root, 0))

        for d in 1:2
            border = subtree_border(sys, root, d)
            interior = subtree_interior(sys, root, d)
            kids = descendants(sys, root, d)

            @test !isempty(border)
            @test length(border) <= 4 * 2^d - 4
            @test allunique(border)
            @test issorted(border)

            # Border and interior partition the subtree, disjointly.
            @test isempty(intersect(Set(border), Set(interior)))
            @test union(Set(border), Set(interior)) == Set(kids)
            @test length(border) + length(interior) == length(kids)

            # Every rim cell really does have a neighbour outside the subtree,
            # and no interior cell does — the definition, spelled out.
            inside = Set(kids)
            grid = levelgrid(sys, d)
            for c in border
                @test any(nb -> !(nb in inside), neighbors(grid, c, 1))
            end
            for c in interior
                @test all(nb -> nb in inside, neighbors(grid, c, 1))
            end
        end

        # Edge() is the narrower adjacency, so its rim can only be a subset of
        # Vertex()'s — a cell exposed only diagonally drops out.
        @test issubset(Set(subtree_border(sys, root, 2; connectivity=Edge())),
            Set(subtree_border(sys, root, 2)))

        # Asking below the cell's own level is an error, not an empty answer.
        @test_throws ArgumentError subtree_border(sys, LevelIndex(2, 0), 1)
        @test_throws ArgumentError subtree_interior(sys, LevelIndex(2, 0), 1)
    end
end

end # module TestFallbacks

# ---------------------------------------------------------------------------
# The allocation-free primitives the substrate is built out of
# (`src/Helpers/small_list.jl`, `src/Helpers/ids.jl`).
#
# Tests for the allocation-free helper primitives shared by every system. The
# authalic-math helpers are covered in `test/fallbacks/authalic.jl`.
# ---------------------------------------------------------------------------

module HelpersUtilTests

using Test
using DiscreteGlobalGrids.Helpers

@testset "shared helpers" begin
    empty = empty_small_list(Val(3), 0)
    @test isempty(empty)
    @test isbitstype(typeof(empty))

    full = small_push(small_push(small_push(empty, 2), 1), 2)
    @test collect(full) == [2, 1, 2]
    @test collect(small_sort(full)) == [1, 2, 2]
    @test_throws BoundsError small_push(full, 3)
    @test_throws BoundsError full[4]

    # The stack half of the list, which the subtree walks in
    # `src/fallbacks/subtree_iterators.jl` push, peek, advance and pop once per
    # tree node. Both stay allocation-free, and both refuse an empty list rather
    # than reading or writing past its end.
    @test collect(small_pop(full)) == [2, 1]
    @test collect(small_pop(small_pop(full))) == [2]
    @test isempty(small_pop(small_pop(small_pop(full))))
    @test_throws BoundsError small_pop(empty)

    @test collect(small_setlast(full, 9)) == [2, 1, 9]
    @test collect(small_setlast(small_pop(full), 9)) == [2, 9]
    @test_throws BoundsError small_setlast(empty, 1)

    # The accumulator half: a slot found by a linear scan rather than by being
    # the top, which is how the halo walk's per-face rectangle table merges.
    @test collect(small_setindex(full, 9, 1)) == [9, 1, 2]
    @test collect(small_setindex(full, 9, 2)) == [2, 9, 2]
    @test collect(small_setindex(full, 9, 3)) == [2, 1, 9]
    @test_throws BoundsError small_setindex(full, 9, 0)
    @test_throws BoundsError small_setindex(full, 9, 4)
    @test_throws BoundsError small_setindex(empty, 9, 1)
    let list = full
        small_setindex(list, 9, 2)               # warm up
        @test (@allocated small_setindex(list, 9, 2)) == 0
    end

    # A pop leaves the dropped slot's contents in place; `len` is what puts it
    # out of reach, and a later push must overwrite rather than resurrect it.
    @test collect(small_push(small_pop(full), 7)) == [2, 1, 7]

    @test isbitstype(typeof(small_pop(full)))
    @test isbitstype(typeof(small_setlast(full, 9)))
    let stack = full
        small_pop(small_setlast(stack, 9))            # warm up
        @test (@allocated small_pop(small_setlast(stack, 9))) == 0
    end

    @test strictly_increasing(Int[])
    @test strictly_increasing([1, 2, 3])
    @test !strictly_increasing([1, 1, 2])
    @test !strictly_increasing([2, 1])

    @test sorted_index([1, 3, 5], 3) == 2
    @test sorted_index([1, 3, 5], 2) == 0
    @test to_uint64_id("ff") == 0xff
    @test to_uint64_id("0xFF") == 0xff
end

end # module HelpersUtilTests
