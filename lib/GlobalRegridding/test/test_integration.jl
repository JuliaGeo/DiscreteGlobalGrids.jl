# End-to-end regridding with `RasterGrid` and real weights.

import DimensionalData as DD

# Cell centres whose midpointed outer edges are `lo` and `hi`.
t6_centres(lo, hi, n) =
    collect(range(lo + (hi - lo) / (2n); step = (hi - lo) / n, length = n))

t6_space(xs, ys) = RasterGrid(DD.DimArray(zeros(length(xs), length(ys)),
    (DD.X(collect(xs)), DD.Y(collect(ys)))))

"""
    t6_raster(f, xs, ys; yfirst = false) -> DimArray

`f(lon, lat)` sampled on the cell centres `xs × ys`. The lookups are whatever
`xs` and `ys` are — pass a descending vector for a reverse-ordered lookup — and
`yfirst` stores the array as `(Y, X)` instead of `(X, Y)`.
"""
function t6_raster(f, xs, ys; yfirst = false)
    data = [Float64(f(x, y)) for x in xs, y in ys]
    xd, yd = DD.X(collect(xs)), DD.Y(collect(ys))
    yfirst && return DD.DimArray(permutedims(data), (yd, xd))
    return DD.DimArray(data, (xd, yd))
end

# Select destination cells by location, not assumed index arithmetic.
t6_lat(space, i) = GO.UnitSpherical.GeographicFromUnitSphere()(cellcentroid(space, i))[2]

t6_mass(space, values) =
    sum(values[i] * GR.cellarea(space, i) for i in 1:ncells(space))

# Counts the weight builds a regrid asks for, whatever the method underneath.
mutable struct T6CountingMethod{M<:AbstractRegriddingMethod} <: AbstractRegriddingMethod
    inner::M
    builds::Int
end

T6CountingMethod(inner) = T6CountingMethod(inner, 0)

function build_weights!(coo::WeightCOO, method::T6CountingMethod,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    countbuild!(method)
    return build_weights!(coo, method.inner, dst_space, dst_inds, src_space, src_inds)
end

# A lookup that names cells of its own, as a DGGS cell axis does.
struct T6Grid end
Base.show(io::IO, ::T6Grid) = print(io, "T6Grid()")

struct T6Cell
    id::Int
end

GR.dimsource(::DD.Lookups.Lookup{T6Cell}) = T6Grid()

@testset "Integration" begin

    @testset "conservation" begin
        f(lon, lat) = 2.0 + sind(2 * lat) + 0.25 * cosd(lon)
        src = t6_raster(f, t6_centres(-180, 180, 36), t6_centres(-90, 90, 18))
        dst = t6_space(t6_centres(-180, 180, 9), t6_centres(-90, 90, 6))
        srcspace = RasterGrid(src)
        mass = t6_mass(srcspace, vec(parent(src)))

        # Extensive regridding preserves the global source integral.
        sums = regrid(src; to = dst, method = Conservative(),
            missingpolicy = Extensive())
        @test sum(sums) ≈ mass rtol = 1e-10

        # Weighted means recover the same integral through destination areas.
        means = regrid(src; to = dst, method = Conservative())
        @test t6_mass(dst, means) ≈ mass rtol = 1e-10
    end

    @testset "NaN padding invariance" begin
        f(lon, lat) = 2.0 + sind(2 * lat) + 0.25 * cosd(lon)
        # The same field over ±60° of latitude, and over the whole sphere with
        # everything past ±60° NaN.
        small = t6_raster(f, t6_centres(-180, 180, 36), t6_centres(-60, 60, 12))
        padded = t6_raster((lon, lat) -> abs(lat) < 60 ? f(lon, lat) : NaN,
            t6_centres(-180, 180, 36), t6_centres(-90, 90, 18))
        dst = t6_space(t6_centres(-180, 180, 9), t6_centres(-90, 90, 9))

        a = regrid(small; to = dst, method = Conservative())
        b = regrid(padded; to = dst, method = Conservative())

        # Invalid padding does not affect fully valid interior cells.
        inner = findall(i -> abs(t6_lat(dst, i)) <= 45, 1:ncells(dst))
        @test !isempty(inner)
        @test !any(isnan, b[inner])
        @test b[inner] ≈ a[inner]

        # All-NaN destinations are blanked.
        polar = findall(i -> abs(t6_lat(dst, i)) >= 75, 1:ncells(dst))
        @test !isempty(polar)
        @test all(isnan, b[polar])

        # Bilinear stencils renormalize after one point becomes invalid.
        constant = t6_raster((lon, lat) -> 5.0,
            t6_centres(-180, 180, 36), t6_centres(-90, 90, 18))
        constant[10, 9] = NaN
        # Half-cell offsets produce four equal stencil weights.
        offset = t6_space(t6_centres(-175, 185, 36), t6_centres(-85, 95, 18))
        @test all(≈(5.0), regrid(constant; to = offset, method = BilinearPoint()))
    end

    @testset "orientation" begin
        # Asymmetry exposes flipped axes or transposed flattening.
        g(lon, lat) = sind(lon) + 2 * lat / 90
        xs = t6_centres(-180, 180, 24)
        ys = t6_centres(-90, 90, 12)
        dst = t6_space(t6_centres(-180, 180, 8), t6_centres(-90, 90, 4))

        forward = regrid(t6_raster(g, xs, ys); to = dst, method = Conservative())
        # Reversing latitude lookup order preserves the result.
        reversed = regrid(t6_raster(g, xs, reverse(ys)); to = dst,
            method = Conservative())
        # And so must the same raster stored `(Y, X)`.
        transposed = regrid(t6_raster(g, xs, ys; yfirst = true); to = dst,
            method = Conservative())

        @test forward ≈ reversed
        @test forward ≈ transposed
    end

    @testset "N-D in one plan" begin
        f(lon, lat) = 2.0 + sind(2 * lat) + 0.25 * cosd(lon)
        xs = t6_centres(-180, 180, 24)
        ys = t6_centres(-90, 90, 12)
        cube = DD.DimArray([f(x, y) + m for x in xs, y in ys, m in 1:12],
            (DD.X(xs), DD.Y(ys), DD.Ti(1:12)))
        dst = t6_space(t6_centres(-180, 180, 8), t6_centres(-90, 90, 4))
        method = T6CountingMethod(Conservative())

        out = regrid(cube; to = dst, method)

        # Twelve monthly fields cost one weight construction, not twelve.
        @test method.builds == 1
        # The destination's own axes first, the month dimension after them
        # unchanged.
        @test size(out) == (8, 4, 12)
        @test collect(DD.lookup(out, DD.Ti)) == 1:12
        # Each slice is that slice's own 2-D regrid: the slice loop neither
        # transposes nor reorders.
        for k in 1:12
            @test Array(out)[:, :, k] ≈
                  Array(regrid(cube[:, :, k]; to = dst, method = Conservative()))
        end
        @test method.builds == 1
    end

    @testset "destination labelling" begin
        f(lon, lat) = 2.0 + sind(2 * lat) + 0.25 * cosd(lon)
        src = t6_raster(f, t6_centres(-180, 180, 36), t6_centres(-90, 90, 18))
        dst = t6_space(t6_centres(-180, 180, 9), t6_centres(-90, 90, 6))

        # An area method returns the destination's own axes, sampled as the
        # cells it integrated over, with the space's edges as their bounds.
        area = regrid(src; to = dst, method = Conservative())
        @test DD.dims(area) isa Tuple{<:DD.X,<:DD.Y}
        @test collect(DD.lookup(area, DD.X)) == t6_centres(-180, 180, 9)
        @test DD.sampling(DD.lookup(area, DD.X)) isa DD.Lookups.Intervals
        @test DD.sampling(DD.lookup(area, DD.Y)) isa DD.Lookups.Intervals
        @test DD.intervalbounds(DD.lookup(area, DD.X))[1] == (-180.0, -140.0)
        @test DD.intervalbounds(DD.lookup(area, DD.Y))[1] == (-90.0, -60.0)

        # Labelling only reshapes: the values are the cell-index vector the
        # same regrid off a bare array returns.
        @test vec(parent(area)) == regrid(parent(src); to = dst,
            from = RasterGrid(src), method = Conservative())

        # A point sample says so instead, over the same cell centres.
        point = regrid(src; to = dst, method = NearestCell())
        @test DD.sampling(DD.lookup(point, DD.X)) isa DD.Lookups.Points
        @test collect(DD.lookup(point, DD.X)) == collect(DD.lookup(area, DD.X))

        # `sampling` overrides the method's rule in both directions, and
        # changes nothing but the label.
        forced = regrid(src; to = dst, method = Conservative(),
            sampling = DD.Lookups.Points())
        @test DD.sampling(DD.lookup(forced, DD.X)) isa DD.Lookups.Points
        @test parent(forced) == parent(area)
        @test DD.sampling(DD.lookup(regrid(src; to = dst, method = NearestCell(),
                sampling = DD.Lookups.Intervals(DD.Lookups.Center())), DD.X)) isa
              DD.Lookups.Intervals

        # A `(Y, X)` destination labels in its own array order.
        yfirst = RasterGrid(DD.DimArray(zeros(6, 9),
            (DD.Y(t6_centres(-90, 90, 6)), DD.X(t6_centres(-180, 180, 9)))))
        @test DD.dims(regrid(src; to = yfirst)) isa Tuple{<:DD.Y,<:DD.X}
    end

    @testset "a cell axis names its own source" begin
        data = DD.DimArray(zeros(32), (DD.Dim{:Cells}(DD.Lookups.Categorical(
            [T6Cell(i) for i in 1:32]; order = DD.Lookups.Unordered())),))
        dst = t6_space(t6_centres(-180, 180, 8), t6_centres(-90, 90, 4))

        # Without `from` the source space is derived from the data, and a cell
        # axis is not a raster lattice: name the grid it holds, not `xdim`.
        @test_throws "from = T6Grid()" regrid(data; to = dst)
    end

    @testset "cross-method agreement" begin
        # All methods closely reproduce a smooth, slowly varying field.
        h(lon, lat) = 10 + sind(lat) + 0.5 * cosd(lon) * cosd(lat)
        src = t6_raster(h, t6_centres(-180, 180, 72), t6_centres(-90, 90, 36))
        dst = t6_space(t6_centres(-180, 180, 18), t6_centres(-90, 90, 9))

        conservative = regrid(src; to = dst, method = Conservative())
        nearest = regrid(src; to = dst, method = NearestCell())
        bilinear = regrid(src; to = dst, method = BilinearPoint())

        @test maximum(abs, nearest .- conservative) < 0.15
        @test maximum(abs, bilinear .- conservative) < 0.15
    end

    @testset "bilinear across the longitude seam" begin
        # One destination centre lies on the ±180° source seam.
        src = t6_raster((lon, lat) -> sind(lon),
            t6_centres(-180, 180, 36), t6_centres(-90, 90, 18))
        dst = t6_space(t6_centres(-175, 185, 36), t6_centres(-90, 90, 18))
        seam = GR.localindex(dst, 36, 9)

        # Global bilinear interpolation wraps across the longitude seam.
        @test regrid(src; to = dst, method = BilinearPoint())[seam] ≈ 0 atol = 1e-12
        @test abs(regrid(src; to = dst, method = NearestCell())[seam]) > 0.08

        # Regional rasters do not report a longitude period.
        @test GR.chartperiod(t6_space(t6_centres(-40, 40, 8),
            t6_centres(-20, 20, 4))) == (nothing, nothing)
    end

    @testset "bilinear west of a regional raster" begin
        # A point just west of a regional raster stays just west on its chart;
        # folded a period east, the non-periodic axis clamps it to the east column.
        xs, ys = t6_centres(0, 120, 24), t6_centres(-30, 30, 12)
        region = t6_space(xs, ys)
        @test GR._onbranch(region.xedges, -1.0, 360.0) ≈ -1.0
        @test GR._onbranch(region.xedges, 61.0, 360.0) ≈ 61.0
        @test GR._onbranch(region.xedges, -1.0, nothing) == -1.0
        point = GO.UnitSpherical.UnitSphereFromGeographic()((-1.0, 10.0))
        @test GR.chartcoords(region, point)[1] ≈ -1.0

        # Pad cells clamp to the adjacent edge column, and destination tiling
        # changes nothing: the pad tile discovers the stencil's source chunk.
        f(lon, lat) = 2.0 + 0.01 * lon + 0.03 * lat
        data = t6_raster(f, xs, ys)
        src = RasterGrid(data; chunks = ([1:8, 9:16, 17:24], [1:12]))
        dxs = t6_centres(-5, 125, 26)
        dst = RasterGrid(DD.DimArray(zeros(12, 26), (DD.Y(ys), DD.X(dxs)));
            chunks = ([1:7, 8:14, 15:20, 21:26], [1:12]))
        untiled = regrid(data; to = dst, from = src, method = BilinearPoint(),
            lazy = false)
        tiled = regrid(data; to = dst, from = src, method = BilinearPoint(),
            lazy = true)
        @test untiled[GR.localindex(dst, 1, 6)] ≈ f(xs[1], ys[6])
        @test untiled[GR.localindex(dst, 26, 6)] ≈ f(xs[end], ys[6])
        # Lazy and eager output share the destination's axes and values.
        @test DD.dims(tiled) == DD.dims(untiled)
        @test all(isequal.(vec(Array(tiled)), vec(Array(untiled))))
    end

    @testset "destination dims echo construction order, lazy and eager" begin
        f(lon, lat) = 1.0 + 0.02 * lon - 0.01 * lat
        data = t6_raster(f, t6_centres(-180, 180, 24), t6_centres(-90, 90, 12))
        dxs, dys = t6_centres(-180, 180, 12), t6_centres(-90, 90, 6)
        xfirst = RasterGrid(DD.DimArray(zeros(12, 6), (DD.X(dxs), DD.Y(dys))))
        yfirst = RasterGrid(DD.DimArray(zeros(6, 12), (DD.Y(dys), DD.X(dxs))))

        # A (Y, X)-constructed destination comes back (Y, X)-shaped.
        eagery = regrid(data; to = yfirst, lazy = false)
        @test DD.dims(eagery, 1) isa DD.Y
        @test DD.dims(eagery, 2) isa DD.X
        @test size(eagery) == (6, 12)

        # The lazy result carries the same dims and values in either order.
        for dst in (xfirst, yfirst)
            eager = regrid(data; to = dst, lazy = false)
            lazy = regrid(data; to = dst, lazy = true)
            @test parent(lazy) isa GR.ShapedRegridArray
            @test DD.dims(lazy) == DD.dims(eager)
            @test size(lazy) == size(eager)
            @test all(isequal.(Array(parent(lazy)), parent(eager)))
        end
    end

    @testset "raster subtrees" begin
        space = RasterGrid(DD.DimArray(zeros(8, 6),
                (DD.X(t6_centres(-180, 180, 8)), DD.Y(t6_centres(-90, 90, 6))));
            chunks = ([1:4, 5:8], [1:3, 4:6]))

        # Chunk rectangles retain the restricted CR cursor in either orientation.
        @test all(GR.subtree(space, cellindices(space, c)) isa
                  CR.Trees.TopDownQuadtreeCursor
                  for c in 1:nchunks(space))
        @test GR.subtree(space, [1, 5, 30]) isa GR.CellSpaceRTree
    end
end
