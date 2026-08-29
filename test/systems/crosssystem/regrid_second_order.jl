# Second-order conservative regridding on the DGG spaces, rather than on the
# analytic toy spaces that `lib/GlobalRegridding/test/test_second_order.jl`
# pins the method itself down on. Three things belong to this side:
#
#   * the two hooks the gradient is recovered from — `cellneighbors`, which
#     reads the grid's own edge one-ring, and `celldiameter`, which bounds the
#     cell width the method declares as its support radius;
#   * that conservation and the second-order accuracy survive real cell
#     geometry, in both directions, eagerly and across source chunk seams;
#   * the shapes a caller actually writes: a global lon/lat raster onto a level
#     grid, a Copernicus DEM tile onto a `MultiOrderCoverage` region, and back.
#
# "The exact answer" here is always a cell mean, computed by a degree-8 fan
# quadrature over each cell's own polygon — the rule the toy suite integrates
# with — because a cell mean is what a conservative method reproduces.

module RegridSecondOrderTests

using Test
import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import DimensionalData as DD
import ConservativeRegridding as CR
import GeoInterface as GI
import GeometryOps as GO
import Extents

const CD = DGG.CopernicusDEM
const US = GO.UnitSpherical
const TO_LONLAT = US.GeographicFromUnitSphere()

const FIRST = GR.Conservative()
const SECOND = GR.ConservativeSecondOrder()

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

# Vectors, spelled out: `LinearAlgebra` is not a dependency of the test
# environment and three lines are cheaper than making it one.
_dot(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
_cross(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])
_norm(a) = sqrt(_dot(a, a))

const RULE = CR.TriangleQuadrature.reference_rule(8)

# `∫ f dA` over a spherical polygon, as a fan of geodesic triangles under the
# degree-8 reference rule.
function polygonintegral(f, pts)
    bary, w = RULE
    total = 0.0
    a = pts[1]
    for i in 2:(length(pts) - 1)
        b, c = pts[i], pts[i + 1]
        det = _dot(a, _cross(b, c))
        for k in eachindex(w)
            l1, l2, l3 = bary[k]
            p = (l1 * a[1] + l2 * b[1] + l3 * c[1],
                 l1 * a[2] + l2 * b[2] + l3 * c[2],
                 l1 * a[3] + l2 * b[3] + l3 * c[3])
            np = _norm(p)
            total += w[k] * f((p[1] / np, p[2] / np, p[3] / np)) * abs(det) / np^3
        end
    end
    return total
end

# A cell's exterior ring as plain 3-vectors, the closing repeat dropped.
function cellring(space, i)
    pts = NTuple{3,Float64}[]
    for p in GI.getpoint(GI.getexterior(GR.getcell(space, i)))
        push!(pts, (Float64(GI.x(p)), Float64(GI.y(p)), Float64(GI.z(p))))
    end
    pop!(pts)
    return pts
end

"Cell means of `f` over every cell of `space`, by quadrature over its polygon."
cellmeans(f, space) = map(1:GR.ncells(space)) do i
    ring = cellring(space, i)
    return polygonintegral(f, ring) / polygonintegral(_ -> 1.0, ring)
end

lonlat(p) = TO_LONLAT(GO.UnitSphericalPoint(p[1], p[2], p[3]))

# The measured angular diameter of one cell: the widest vertex-to-vertex arc.
# `celldiameter` is a bound on this, over the whole space.
function measureddiameter(space, i)
    ring = cellring(space, i)
    d = 0.0
    for a in ring, b in ring
        d = max(d, US.spherical_distance(GO.UnitSphericalPoint(a...),
            GO.UnitSphericalPoint(b...)))
    end
    return d
end

# A smooth global field with structure at every scale the tests resolve.
smooth(p) = 2 + sin(2 * asin(clamp(p[3], -1, 1)))^2 * cos(3 * atan(p[2], p[1])) +
            0.5 * p[1] * p[3]

"Synthetic terrain, in metres, as a function of an ambient unit-sphere point."
function terrain(p)
    lon, lat = lonlat(p)
    return 1200.0 + 800.0 * sinpi((lon - 10) * 1.5) * cospi((lat - 46) * 1.3) +
           300.0 * cospi((lon + lat) / 1.7)
end

_axis(D, centres, step) = D(DD.Sampled(collect(centres); span = DD.Regular(step),
    sampling = DD.Intervals(DD.Center()), order = DD.ForwardOrdered()))

# Abutting intervals, so the raster's cell edges tile its extent exactly and a
# conservative regrid off it can conserve.
function lonlatdims(step; lon = (-180.0, 180.0), lat = (-90.0, 90.0))
    return (_axis(DD.X, (lon[1] + step / 2):step:lon[2], step),
            _axis(DD.Y, (lat[1] + step / 2):step:lat[2], step))
end

"A `DimArray` of the quadrature cell means of `f` over the lattice `dims`."
function meanraster(f, dims)
    space = GR.RasterGrid(dims)
    return DD.DimArray(reshape(cellmeans(f, space), map(length, dims)), dims), space
end

eager(data, dst, src, method; policy = GR.Weighted()) =
    parent(DGG.regrid(data; to = dst, from = src, method,
        missingpolicy = policy, lazy = false))

mass(f, space) = sum(f[i] * GR.cellarea(space, i) for i in 1:GR.ncells(space))

rms(a, b) = sqrt(sum(abs2, vec(a) .- vec(b)) / length(b))

# Is every vertex of destination cell `j` inside the lon/lat box, with `margin`
# degrees to spare? Vertices only, so `margin` has to cover the bulge of the
# geodesic edges between them — 0.02° is an order of magnitude more than the
# widest cell here needs.
function ringinside(space, j, box, margin)
    for p in GI.getpoint(GI.getexterior(GR.getcell(space, j)))
        lon, lat = lonlat((GI.x(p), GI.y(p), GI.z(p)))
        (box.X[1] + margin <= lon <= box.X[2] - margin &&
         box.Y[1] + margin <= lat <= box.Y[2] - margin) || return false
    end
    return true
end

# --------------------------------------------------------------------------
# (1) The adjacency hooks
# --------------------------------------------------------------------------

@testset "the hooks a gradient is recovered from" begin
    hex = DGG.DGGSpace(DGG.levelgrid(DGG.IGeo7System(), 2))

    @testset "a complete hexagonal level" begin
        rings = [GR.cellneighbors(hex, i) for i in 1:GR.ncells(hex)]
        # Six neighbours everywhere but at the twelve pentagons, no cell is its
        # own neighbour, and every index names a cell of this space. Two
        # non-collinear neighbours is all a gradient fit needs, so every cell
        # here carries a full second-order stencil.
        @test all(r -> length(r) in (5, 6), rings)
        @test count(r -> length(r) == 5, rings) == 12
        @test all(i -> !(i in rings[i]), 1:GR.ncells(hex))
        @test all(r -> all(k -> 1 <= k <= GR.ncells(hex), r), rings)
        # Adjacency is symmetric, which is what makes a chunk's halo — the
        # cells that read it — the same set as the cells it reads.
        @test all(i in rings[k] for i in 1:GR.ncells(hex) for k in rings[i])
    end

    @testset "a partial grid clips its rim" begin
        level3 = DGG.levelgrid(DGG.IGeo7System(), 3)
        region = DGG.covering(DGG.CellVector(level3),
            Extents.Extent(X = (-140.0, 40.0), Y = (-20.0, 60.0)))
        part = DGG.DGGSpace(DGG.PartialGrid(region))
        rings = [GR.cellneighbors(part, i) for i in 1:GR.ncells(part)]
        # Local indices into the collection, not the complete level's global
        # ones: a stencil column is a source column of this space.
        @test all(r -> all(k -> 1 <= k <= GR.ncells(part), r), rings)
        @test all(i in rings[k] for i in 1:GR.ncells(part) for k in rings[i])
        # The rim is real — a one-sided fit is the documented behaviour there,
        # not an error — and the interior is untouched by the clipping.
        @test count(r -> length(r) < 5, rings) > 0
        @test count(r -> length(r) >= 5, rings) > GR.ncells(part) ÷ 2
    end

    @testset "a Copernicus tile is a quadrilateral lattice" begin
        # `Edge()`, so an interior pixel answers its four edge neighbours and
        # not the eight-cell ring: four is already a well-conditioned fit for
        # two parameters.
        twin = CD.CopernicusDEMSystem{30}()
        nrows = Int(CD.lat_intervals(twin))
        ncols = Int(CD.ncols_at(twin, 46))
        tile = DGG.DGGSpace(DGG.subtree(twin, CD.tilecell(twin, 46, 10), 1))
        @test GR.ncells(tile) == nrows * ncols
        interior = [j * ncols + i + 1 for j in 1:(nrows - 2) for i in 1:(ncols - 2)]
        @test all(length(GR.cellneighbors(tile, k)) == 4 for k in interior)
        # The vertex-contact fallback counts the diagonals as well, which is
        # exactly the reason this space answers from its own topology.
        @test length(invoke(GR.cellneighbors, Tuple{GR.RegridSpace,Int},
            tile, first(interior))) == 8

        # At a pole the tile's top row of pixels are triangles meeting at the
        # apex, so a vertex ring there is the whole row. The edge ring is not:
        # the westmost pixel of the polemost row sees its eastward neighbour
        # and the pixel below it, and nothing else.
        polar = DGG.DGGSpace(DGG.subtree(twin, CD.tilecell(twin, 89, 10), 1))
        polecols = Int(CD.ncols_at(twin, 89))
        @test polecols >= 3
        ring = GR.cellneighbors(polar, 1)
        @test sort(ring) == [2, polecols + 1]
        # The cell across the apex is in the vertex ring and not in this one.
        @test polecols in invoke(GR.cellneighbors, Tuple{GR.RegridSpace,Int}, polar, 1)
        @test !(polecols in ring)
    end

    @testset "celldiameter bounds the cells and is the support radius" begin
        twin = CD.CopernicusDEMSystem{30}()
        tile = DGG.DGGSpace(DGG.subtree(twin, CD.tilecell(twin, 46, 10), 1))
        for space in (hex, tile)
            bound = GR.celldiameter(space)
            # A bound, not a measurement: overestimating costs discovery work
            # and nothing else, so the only law is that nothing exceeds it.
            @test bound > 0
            @test all(measureddiameter(space, i) <= bound
                      for i in 1:max(1, GR.ncells(space) ÷ 24):GR.ncells(space))
            # The method's stencil reaches one cell past its source cell, so
            # this is the radius discovery has to dilate a destination cap by.
            @test GR.supportradius(SECOND, space) == bound
        end
    end

    @testset "the geometric fallback is a superset" begin
        # On hexagons the vertex ring and the edge ring are the same six cells,
        # so the grid's answer is not merely contained in the fallback's — it
        # is the fallback's. What matters for the method is the containment.
        for i in (1, 17, 200, GR.ncells(hex))
            grid = GR.cellneighbors(hex, i)
            generic = invoke(GR.cellneighbors, Tuple{GR.RegridSpace,Int}, hex, i)
            @test issubset(grid, generic)
            @test sort(collect(grid)) == sort(collect(generic))
        end
    end
end

# --------------------------------------------------------------------------
# (2) A global lon/lat field onto IGeo7, and back
# --------------------------------------------------------------------------

const GLOBAL_DIMS = lonlatdims(1.0)
const GLOBAL_DATA, GLOBAL_SPACE = meanraster(smooth, GLOBAL_DIMS)
const GLOBAL_MEANS = vec(parent(GLOBAL_DATA))
const HEX4 = DGG.DGGSpace(DGG.levelgrid(DGG.IGeo7System(), 4))
const HEX4_MEANS = cellmeans(smooth, HEX4)
# The `Weighted` answers both directions are judged by, built once: a weight
# build over 64800 source cells is the cost of this file and nothing here needs
# two of them.
const ONTO_HEX = map(m -> eager(GLOBAL_DATA, HEX4, GLOBAL_SPACE, m), (FIRST, SECOND))
const ONTO_RASTER = map(m -> eager(HEX4_MEANS, GLOBAL_SPACE, HEX4, m), (FIRST, SECOND))

@testset "a global 1° field onto IGeo7 level 4" begin
    total = mass(GLOBAL_MEANS, GLOBAL_SPACE)
    @test total ≈ 4pi * 2 rtol = 1e-6      # the field's mean is 2 by symmetry

    @testset "Extensive conserves the source integral" begin
        # The gradient correction has zero mean over each source cell, so it
        # moves mass between destinations and never creates or destroys it.
        # First order is the same claim with a zero gradient; both hold to
        # rounding, which is the point of "conservative".
        for method in (FIRST, SECOND)
            sums = eager(GLOBAL_DATA, HEX4, GLOBAL_SPACE, method; policy = GR.Extensive())
            @test sum(sums) ≈ total rtol = 1e-12
        end
    end

    @testset "the second-order error is the smaller one" begin
        first, second = ONTO_HEX
        @test all(isfinite, first) && all(isfinite, second)
        e1, e2 = rms(first, HEX4_MEANS), rms(second, HEX4_MEANS)
        # An order of magnitude, not a factor: the source cell is about a
        # third of a destination cell across, so the flattening first order
        # does is the whole of its error and the correction removes most of it.
        @test e2 < e1 / 10
        @test maximum(abs, second .- HEX4_MEANS) < maximum(abs, first .- HEX4_MEANS) / 10
    end

    @testset "a chunked source gives the eager answer" begin
        # The seam is what this is for: a destination straddling two source
        # chunks takes weights from both blocks, and each block's stencils read
        # one ring past its own cells. Sixteen source chunks put a seam through
        # most of the destination.
        chunked = GR.RasterGrid(GLOBAL_DIMS;
            chunks = ([1:90, 91:180, 181:270, 271:360],
                      [1:45, 46:90, 91:135, 136:180]))
        @test GR.nchunks(chunked) == 16
        reference = ONTO_HEX[2]
        lazy = Array(parent(DGG.regrid(GLOBAL_DATA; to = HEX4, from = chunked,
            method = SECOND, lazy = true)))
        @test maximum(abs, lazy .- reference) < 1e-12
    end
end

@testset "IGeo7 level 4 back onto the 1° lattice" begin
    total = mass(HEX4_MEANS, HEX4)

    @testset "Extensive conserves the source integral" begin
        for method in (FIRST, SECOND)
            sums = eager(HEX4_MEANS, GLOBAL_SPACE, HEX4, method; policy = GR.Extensive())
            @test sum(sums) ≈ total rtol = 1e-12
        end
    end

    @testset "the second-order error is the smaller one" begin
        first, second = ONTO_RASTER
        @test all(isfinite, first) && all(isfinite, second)
        # The coarse-to-fine direction: here the source cell is the wide one,
        # so both errors are larger and the ratio is what is under test.
        @test rms(second, GLOBAL_MEANS) < rms(first, GLOBAL_MEANS) / 10
    end

    @testset "a chunked plan gives the eager answer" begin
        reference = ONTO_RASTER[2]
        plan = GR.ChunkedPlan(SECOND, GR.Weighted(), GLOBAL_SPACE, HEX4)
        chunked = collect(GR.regrid(HEX4_MEANS, plan))
        @test maximum(abs, vec(parent(chunked)) .- vec(reference)) < 1e-12
    end
end

# --------------------------------------------------------------------------
# (3) A Copernicus DEM tile and an IGeo7 covering
# --------------------------------------------------------------------------

# The scaled twin of the Copernicus lattice: 30 rows to a 1° tile rather than
# 3600, so one tile's pixels are 900 cells and the whole section runs in
# seconds. Nothing is downloaded; the pixels are the system's own geometry and
# the elevations are synthetic.
const TWIN = CD.CopernicusDEMSystem{30}()
const TILE_LAT, TILE_LON = 46, 10
const TILE = DGG.DGGSpace(DGG.subtree(TWIN, CD.tilecell(TWIN, TILE_LAT, TILE_LON), 1))
const TILE_BOX = Extents.Extent(X = (Float64(TILE_LON), TILE_LON + 1.0),
    Y = (Float64(TILE_LAT), TILE_LAT + 1.0))
# Padded, so the covering's cells reach past the tile on every side and the
# tile's pixels are covered whole. A `MultiOrderCoverage` of the tile box
# exactly would leave the destination short of the pixel rows along the
# northern and southern edges — the covering's edges are geodesics and the
# pixel rows' are not — and no conservation claim would survive it.
const COVER_LEVEL = 7
const REGION = DGG.query(DGG.IGeo7System(),
    DGG.MultiOrderCoverage(Extents.Extent(X = (TILE_LON - 0.05, TILE_LON + 1.05),
        Y = (TILE_LAT - 0.05, TILE_LAT + 1.05))); level = COVER_LEVEL)
const COVER_CELLS = DGG.CellVector(REGION)
const COVER = DGG.DGGSpace(DGG.PartialGrid(COVER_CELLS))
# The covering cells that lie wholly inside the tile, with room for the bulge
# of their geodesic edges: the ones the tile can answer for, and the ones the
# tile covers when the direction is reversed.
const INNER = [j for j in 1:GR.ncells(COVER) if ringinside(COVER, j, TILE_BOX, 0.02)]

@testset "a Copernicus tile onto an IGeo7 covering" begin
    dem = [terrain(Tuple(GR.cellcentroid(TILE, i))) for i in 1:GR.ncells(TILE)]
    lo, hi = extrema(dem)

    @test GR.ncells(TILE) == 900
    @test 0 < length(INNER) < GR.ncells(COVER)

    first = eager(dem, COVER, TILE, FIRST)
    second = eager(dem, COVER, TILE, SECOND)

    @testset "values" begin
        # Inside the tile every destination has full coverage and a value; the
        # covering overhangs the tile, and those cells are missing on both
        # methods alike.
        @test all(isfinite, first[INNER])
        @test all(isfinite, second[INNER])
        @test count(isnan, first) == count(isnan, second) > 0
        # The gradient is doing something — a smooth DEM is nowhere flat — but
        # a smooth field gives it nothing to overshoot with. The margin is a
        # tenth of the field's range: the weights are signed, so this is a
        # sanity bound, not a maximum principle the method claims.
        margin = (hi - lo) / 10
        @test maximum(abs, second[INNER] .- first[INNER]) > 1.0
        @test all(lo - margin <= v <= hi + margin for v in second[INNER])
    end

    @testset "Extensive totals agree" begin
        e1 = eager(dem, COVER, TILE, FIRST; policy = GR.Extensive())
        e2 = eager(dem, COVER, TILE, SECOND; policy = GR.Extensive())
        # The covering takes every pixel whole, so first order recovers the
        # tile's integral to rounding and second order has to move mass around
        # without changing the total. The gradient terms cancel numerically
        # rather than identically — each source cell's first moment about its
        # own mean is a sum of overlap moments that cancels to about 1e-10 of
        # the total, not to an ulp — so this is the tolerance, not 1e-15.
        @test sum(e1) ≈ mass(dem, TILE) rtol = 1e-12
        @test sum(e2) ≈ sum(e1) rtol = 1e-9
    end
end

@testset "an IGeo7 patch back onto the Copernicus tile" begin
    # The reverse direction with the roles that make conservation testable:
    # the source is the covering cells the tile contains, so the tile covers
    # every one of them whole.
    patch = DGG.DGGSpace(DGG.PartialGrid(DGG.IGeo7System(), COVER_LEVEL,
        sort!([COVER_CELLS[j] for j in INNER])))
    @test GR.ncells(patch) == length(INNER)
    f = [terrain(Tuple(GR.cellcentroid(patch, i))) for i in 1:GR.ncells(patch)]
    total = mass(f, patch)
    for (method, tol) in ((FIRST, 1e-12), (SECOND, 1e-9))
        sums = eager(f, TILE, patch, method; policy = GR.Extensive())
        @test all(isfinite, sums)
        @test sum(sums) ≈ total rtol = tol
    end
    # A constant field has a zero gradient, so the second-order answer is the
    # first-order one — to rounding, not bit for bit: the stencil's self term
    # is the negated sum of its neighbour coefficients, and that sum cancels
    # numerically rather than identically.
    ones1 = eager(ones(GR.ncells(patch)), TILE, patch, FIRST; policy = GR.Extensive())
    ones2 = eager(ones(GR.ncells(patch)), TILE, patch, SECOND; policy = GR.Extensive())
    @test ones1 ≈ ones2 rtol = 1e-12
    @test sum(ones1) ≈ sum(GR.cellarea(patch, i) for i in 1:GR.ncells(patch)) rtol = 1e-12
end

@testset "the tutorial's path: a lon/lat DEM onto a covering" begin
    # `RasterGrid` source, `DGGSpace(PartialGrid)` destination — what
    # `docs/src/tutorials/hydrology.jl` does, at a size that runs in a second.
    dims = lonlatdims(1 / 30; lon = (Float64(TILE_LON), TILE_LON + 1.0),
        lat = (Float64(TILE_LAT), TILE_LAT + 1.0))
    data, space = meanraster(terrain, dims)
    @test GR.ncells(space) == 900
    total = mass(vec(parent(data)), space)
    exact = cellmeans(terrain, COVER)

    results = map((FIRST, SECOND)) do method
        sums = eager(data, COVER, space, method; policy = GR.Extensive())
        means = eager(data, COVER, space, method)
        @test all(isfinite, means[INNER])
        return sums, means
    end
    # Both totals are the raster's integral: the covering takes it whole.
    @test sum(results[1][1]) ≈ total rtol = 1e-12
    @test sum(results[2][1]) ≈ total rtol = 1e-9
    # And on the cells the covering can actually answer for, the correction is
    # an improvement rather than merely a difference.
    @test rms(results[2][2][INNER], exact[INNER]) <
          rms(results[1][2][INNER], exact[INNER])
end

# --------------------------------------------------------------------------
# (4) A masked source
# --------------------------------------------------------------------------

@testset "a hole in the source" begin
    # Coverage is kept apart from the signed weights, so `Weighted` thresholds
    # on the area a destination actually draws valid values from, not on the
    # sum of its weights — which the gradient terms make signed.
    step = 3.0
    dims = lonlatdims(step)
    data, space = meanraster(smooth, dims)
    holed = copy(data)
    xs = [x for x in DD.lookup(dims[1]) if -10 < x < 10]
    ys = [y for y in DD.lookup(dims[2]) if -10 < y < 10]
    for (i, x) in enumerate(DD.lookup(dims[1])), (j, y) in enumerate(DD.lookup(dims[2]))
        (x in xs && y in ys) && (holed[i, j] = NaN)
    end
    # The blanked area is the union of those cells, which is their centres
    # grown by half a step — not the box that selected them.
    hole = Extents.Extent(X = (minimum(xs) - step / 2, maximum(xs) + step / 2),
        Y = (minimum(ys) - step / 2, maximum(ys) + step / 2))
    @test count(isnan, parent(holed)) == length(xs) * length(ys) > 16

    dst = DGG.DGGSpace(DGG.levelgrid(DGG.IGeo7System(), 3))
    full = eager(data, dst, space, SECOND; policy = GR.Weighted(0.5))
    masked = eager(holed, dst, space, SECOND; policy = GR.Weighted(0.5))

    # A destination lying wholly inside the hole draws no valid area at all.
    inside = [j for j in 1:GR.ncells(dst) if ringinside(dst, j, hole, 0.02)]
    @test length(inside) > 4
    @test all(isnan, masked[inside])

    # Far from the hole nothing moves. "Far" is two rings out: a source cell
    # next to the hole has a NaN in its own gradient stencil, which the fixed
    # operator reads as zero and which therefore biases that cell — the
    # documented limit of the method, and the reason `Conservative` is the
    # advice for a source with missing values in it.
    centre = GO.UnitSphericalPoint(1.0, 0.0, 0.0)
    far = [j for j in 1:GR.ncells(dst)
           if US.spherical_distance(GR.cellcentroid(dst, j), centre) > 0.6]
    @test length(far) > GR.ncells(dst) ÷ 2
    @test masked[far] ≈ full[far]
    @test all(isfinite, full)

    # `Extensive` reports the covered integral, which the hole shrinks but
    # never makes missing.
    sums = eager(holed, dst, space, SECOND; policy = GR.Extensive())
    @test all(isfinite, sums)
    @test 0 < sum(sums) < mass(vec(parent(data)), space)
end

end # module RegridSecondOrderTests
