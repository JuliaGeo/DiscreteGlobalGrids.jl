# End-to-end regrids and the acceptance laws. Owned by task T6.
#
# Nothing here uses `toyspaces.jl`: the point of this file is that the shipped
# pieces — `RasterGrid`, real weight construction, the eager executor and the
# API — agree with each other on geometry, orientation, flattening and missing
# data. Every destination cell here is a four-corner graticule quad and so
# convex, which keeps the GeometryOps Sutherland–Hodgman non-convex-clip defect
# (isolated in `test_conservative.jl`) out of the conservation numbers.

import DimensionalData as DD

# `n` cell centres spanning `lo`..`hi`, chosen so that midpointing them puts the
# outer edges exactly on `lo` and `hi`.
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

# Destination cells are picked out by where they are rather than by an index
# arithmetic these tests would then be assuming.
t6_lat(space, i) = GR.SphereToLonLat()(cellcentroid(space, i))[2]

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

@testset "Integration" begin

    @testset "conservation" begin
        f(lon, lat) = 2.0 + sind(2 * lat) + 0.25 * cosd(lon)
        src = t6_raster(f, t6_centres(-180, 180, 36), t6_centres(-90, 90, 18))
        dst = t6_space(t6_centres(-180, 180, 9), t6_centres(-90, 90, 6))
        srcspace = RasterGrid(src)
        mass = t6_mass(srcspace, vec(parent(src)))

        # `Extensive` sums the intersection areas, so the destination field's
        # global integral is the source's: every source cell's weights must sum
        # to exactly its own area. A weight matrix that is transposed, scaled,
        # or missing the polar rows fails here whatever its row sums are.
        sums = regrid(src; to = dst, method = Conservative(),
            missingpolicy = Extensive())
        @test sum(sums) ≈ mass rtol = 1e-10

        # The destination-side dual: `Weighted` divides by the covered area, so
        # multiplying back by the destination cells' own areas returns the same
        # integral. This is what says the denominators are areas and not row
        # sums of something else.
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

        # Destination cells wholly inside the valid region are untouched by the
        # padding: `Weighted` divides by the *valid* coverage, so NaNs the cell
        # never overlapped cannot reach it — and NaN times a zero weight cannot
        # poison it either.
        inner = findall(i -> abs(t6_lat(dst, i)) <= 45, 1:ncells(dst))
        @test !isempty(inner)
        @test !any(isnan, b[inner])
        @test b[inner] ≈ a[inner]

        # A destination cell every one of whose sources is NaN is blank, not
        # zero and not a division by zero.
        polar = findall(i -> abs(t6_lat(dst, i)) >= 75, 1:ncells(dst))
        @test !isempty(polar)
        @test all(isnan, b[polar])

        # A point method's support reaches past the cells a destination
        # overlaps, so one NaN source cell enters four bilinear stencils.
        # `Weighted` renormalizes each over the three points that survived; a
        # finalize that skipped the division because the block reported no
        # denominator would return three quarters of a constant field instead.
        constant = t6_raster((lon, lat) -> 5.0,
            t6_centres(-180, 180, 36), t6_centres(-90, 90, 18))
        constant[10, 9] = NaN
        # Centres offset half a source cell both ways, so every stencil is four
        # points of equal weight and none of them is blanked at the default
        # threshold.
        offset = t6_space(t6_centres(-175, 185, 36), t6_centres(-85, 95, 18))
        @test all(≈(5.0), regrid(constant; to = offset, method = BilinearPoint()))
    end

    @testset "orientation" begin
        # Asymmetric in both coordinates, so any flip of either axis — or a
        # transposed flattening — moves the answer.
        g(lon, lat) = sind(lon) + 2 * lat / 90
        xs = t6_centres(-180, 180, 24)
        ys = t6_centres(-90, 90, 12)
        dst = t6_space(t6_centres(-180, 180, 8), t6_centres(-90, 90, 4))

        forward = regrid(t6_raster(g, xs, ys); to = dst, method = Conservative())
        # A north-to-south latitude lookup describes the same cells carrying the
        # same values, and must regrid to the same destination field: the cell
        # rings and the data flattening derive their sense from the same
        # lookups, so they cannot disagree.
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
        # Destination cells first, the month dimension after it unchanged.
        @test size(out) == (ncells(dst), 12)
        @test collect(DD.lookup(out, DD.Ti)) == 1:12
        # Each slice is that slice's own 2-D regrid: the slice loop neither
        # transposes nor reorders.
        for k in 1:12
            @test Array(out)[:, k] ≈
                  regrid(cube[:, :, k]; to = dst, method = Conservative())
        end
        @test method.builds == 1
    end

    @testset "cross-method agreement" begin
        # Smooth and slowly varying, so a cell mean, a centre sample and a
        # bilinear sample should all be close; the field spans about 3, so a
        # tolerance of 0.15 is 5% of the signal while a source cell mis-indexed
        # by even one row is O(1) out.
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
        # Destination centres sit half a source cell east of the source centres,
        # so one of them lands exactly on ±180 — the seam between the source's
        # last cell centre and its first.
        src = t6_raster((lon, lat) -> sind(lon),
            t6_centres(-180, 180, 36), t6_centres(-90, 90, 18))
        dst = t6_space(t6_centres(-175, 185, 36), t6_centres(-90, 90, 18))
        seam = GR.cellposition(dst, 36, 9)

        # A global raster's longitude axis closes, so the seam is an interval
        # like any other and the stencil averages sin(175°) with sin(-175°). An
        # axis reported as non-periodic degrades to the nearest centre there and
        # returns one of them instead.
        @test regrid(src; to = dst, method = BilinearPoint())[seam] ≈ 0 atol = 1e-12
        @test abs(regrid(src; to = dst, method = NearestCell())[seam]) > 0.08

        # A raster that does not span the globe must not report a period: a
        # stencil that wrapped would reach across the whole domain.
        @test GR.chartperiod(t6_space(t6_centres(-40, 40, 8),
            t6_centres(-20, 20, 4))) == (nothing, nothing)
    end

    @testset "raster subtrees" begin
        space = RasterGrid(DD.DimArray(zeros(8, 6),
                (DD.X(t6_centres(-180, 180, 8)), DD.Y(t6_centres(-90, 90, 6))));
            chunks = ([1:4, 5:8], [1:3, 4:6]))

        # A chunk is a rectangle of the lattice however its cell positions are
        # numbered, so a chunk's subtree stays the O(1) recursive tree rather
        # than falling back to one cap per cell.
        @test all(GR.subtree(space, cellindices(space, c)) isa GR.RasterCellTree
                  for c in 1:nchunks(space))
        @test GR.subtree(space, [1, 5, 30]) isa GR.RasterFlatTree
    end
end
