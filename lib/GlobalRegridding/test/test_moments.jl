# Spherical polygon area and vector first moment.

import ConservativeRegridding as CR

"""
    moment_ring(n, r, N) -> Vector{USPoint}

An open ring of `N` vertices on the circle of angular radius `r` about `n`, wound
counter-clockwise seen from outside. Inscribed in the cap, so its area and moment
approach the cap's from below as `O(1/N^2)`.
"""
function moment_ring(n, r::Real, N::Int)
    axis = LinearAlgebra.normalize(USPoint(n[1], n[2], n[3]))
    # Any vector off the axis gives the frame; `u × v = axis` sets the winding.
    off = abs(axis[3]) < 0.9 ? USPoint(0.0, 0.0, 1.0) : USPoint(1.0, 0.0, 0.0)
    u = LinearAlgebra.normalize(LinearAlgebra.cross(off, axis))
    v = LinearAlgebra.cross(axis, u)
    return [USPoint((cos(r) .* axis .+
                     sin(r) .* (cos(2pi * k / N) .* u .+ sin(2pi * k / N) .* v))...)
            for k in 0:(N - 1)]
end

"""
    moment_graticule(lon0, lon1, lat0, lat1, n) -> Vector{USPoint}

An open graticule ring with each parallel densified into `n` chords. The meridians
are great circles already, so only the parallels converge, as `O(1/n^2)`.
"""
function moment_graticule(lon0, lon1, lat0, lat1, n::Int)
    ring = USPoint[]
    for k in 0:n                # south parallel, west to east
        push!(ring, toy_point(lon0 + (lon1 - lon0) * k / n, lat0))
    end
    for k in 0:n                # north parallel, east to west
        push!(ring, toy_point(lon1 - (lon1 - lon0) * k / n, lat1))
    end
    return ring
end

"""
    moment_graticule_exact(lon0, lon1, lat0, lat1) -> (area, moment)

Exact area and first moment on the unit sphere, from `∫ x cos φ dφ dλ` with
`∫cos²φ dφ = (φ + sin φ cos φ)/2` and `∫sin φ cos φ dφ = sin²φ/2`.
"""
function moment_graticule_exact(lon0, lon1, lat0, lat1)
    l0, l1 = deg2rad(Float64(lon0)), deg2rad(Float64(lon1))
    f0, f1 = deg2rad(Float64(lat0)), deg2rad(Float64(lat1))
    c = ((f1 + sin(f1) * cos(f1)) - (f0 + sin(f0) * cos(f0))) / 2
    z = (sin(f1)^2 - sin(f0)^2) / 2
    area = (l1 - l0) * (sin(f1) - sin(f0))
    return (area, (c * (sin(l1) - sin(l0)), c * (cos(l0) - cos(l1)), z * (l1 - l0)))
end

"""
    moment_rotation(axis, angle) -> Matrix{Float64}

Rodrigues' rotation matrix, for the equivariance check.
"""
function moment_rotation(axis, angle::Real)
    a = LinearAlgebra.normalize(collect(Float64, axis))
    K = [0.0 -a[3] a[2]; a[3] 0.0 -a[1]; -a[2] a[1] 0.0]
    return LinearAlgebra.I + sin(angle) * K + (1 - cos(angle)) * K^2
end

@testset "Polygon moments" begin
    m = GOCore.Spherical(; radius = 1.0)

    @testset "PolygonMoments algebra" begin
        a = GR.PolygonMoments(1.5, (0.1, 0.2, 0.3))
        b = GR.PolygonMoments(2.5, (0.4, 0.5, 0.6))
        @test zero(GR.PolygonMoments) == GR.PolygonMoments(0.0, (0.0, 0.0, 0.0))
        @test zero(a) === zero(GR.PolygonMoments)
        @test iszero(zero(GR.PolygonMoments))
        @test !iszero(a)
        @test !iszero(GR.PolygonMoments(0.0, (0.0, 1.0, 0.0)))
        @test a + b === GR.PolygonMoments(4.0, (0.5, 0.7, 0.8999999999999999))
        @test a + zero(GR.PolygonMoments) === a
        @test 2 * a === a * 2 === GR.PolygonMoments(3.0, (0.2, 0.4, 0.6))
        @test a * 0.0 === zero(GR.PolygonMoments)
        @test occursin("PolygonMoments(area = 1.5", sprint(show, a))
        @test occursin("moment = (0.1, 0.2, 0.3)", sprint(show, a))
    end

    @testset "a sparse matrix assembles over PolygonMoments" begin
        a = GR.PolygonMoments(1.5, (0.1, 0.2, 0.3))
        b = GR.PolygonMoments(2.5, (0.4, 0.5, 0.6))
        # A split overlap produces duplicate coordinates; `sparse` combines them with `+`.
        S = GR.SparseArrays.sparse([1, 1], [1, 1], [a, b], 1, 1)
        @test eltype(S) === GR.PolygonMoments
        @test S[1, 1] === a + b

        T = GR.SparseArrays.sparse([1, 2], [1, 2], [a, b], 2, 2)
        @test T[1, 1] === a
        @test T[2, 2] === b
        @test T[1, 2] === zero(GR.PolygonMoments)
    end

    @testset "a spherical cap, against its analytic moment" begin
        r = 0.4
        for n in ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0), (0.0, -1.0, 0.0),
            (0.3, -0.7, 0.2), (0.02, 0.01, 1.0))
            axis = LinearAlgebra.normalize(collect(Float64, n))
            pm = GR.polygonmoments(m, moment_ring(n, r, 2000); closed = false)
            @test pm.area ≈ 2pi * (1 - cos(r)) rtol = 1e-5
            @test collect(pm.moment) ≈ pi * sin(r)^2 .* axis rtol = 1e-5
        end

        # Caps from a sliver to most of a hemisphere, and the O(1/N^2) rate.
        axis = (0.3, -0.7, 0.2)
        for r in (0.01, 0.4, 1.4)
            exact = pi * sin(r)^2 .*
                    LinearAlgebra.normalize(collect(Float64, axis))
            errs = map((250, 500, 1000)) do N
                pm = GR.polygonmoments(m, moment_ring(axis, r, N); closed = false)
                LinearAlgebra.norm(collect(pm.moment) .- exact) / max(1e-12,
                    LinearAlgebra.norm(exact))
            end
            @test errs[1] / errs[2] > 3.5
            @test errs[2] / errs[3] > 3.5
        end
    end

    @testset "a graticule cell, against its analytic moment" begin
        boxes = (((-20.0, 40.0), (10.0, 55.0)), ((150.0, 190.0), (-70.0, -30.0)),
            ((0.0, 30.0), (-15.0, 15.0)))
        for (lon, lat) in boxes
            area, moment = moment_graticule_exact(lon[1], lon[2], lat[1], lat[2])
            ring(n) = moment_graticule(lon[1], lon[2], lat[1], lat[2], n)
            pm = GR.polygonmoments(m, ring(400); closed = false)
            @test pm.area ≈ area rtol = 1e-5
            @test collect(pm.moment) ≈ collect(moment) rtol = 1e-5

            # Only the chorded parallels are approximate, and they converge as O(1/n^2).
            errs = map((100, 200, 400)) do n
                q = GR.polygonmoments(m, ring(n); closed = false)
                LinearAlgebra.norm(collect(q.moment) .- collect(moment)) /
                LinearAlgebra.norm(collect(moment))
            end
            @test errs[1] / errs[2] > 3.5
            @test errs[2] / errs[3] > 3.5
        end
    end

    @testset "both integrals scale with radius squared" begin
        pts = moment_graticule(-20.0, 40.0, 10.0, 55.0, 64)
        unit = GR.polygonmoments(m, pts; closed = false)
        for R in (0.5, 6371.0)
            scaled = GR.polygonmoments(GOCore.Spherical(; radius = R), pts;
                closed = false)
            @test scaled.area ≈ unit.area * R^2 rtol = 1e-14
            @test collect(scaled.moment) ≈ collect(unit.moment) .* R^2 rtol = 1e-14
            # The scaling leaves the mean position alone: a direction, at any radius.
            @test collect(scaled.moment) ./ scaled.area ≈
                  collect(unit.moment) ./ unit.area rtol = 1e-12
        end
    end

    @testset "rotating the ring rotates the moment and keeps the area" begin
        R = moment_rotation((0.2, -0.5, 0.84), 0.9)
        for pts in (moment_graticule(-20.0, 40.0, 10.0, 55.0, 64),
            moment_ring((0.3, -0.7, 0.2), 0.6, 256))
            pm = GR.polygonmoments(m, pts; closed = false)
            turned = [USPoint((R * collect(p))...) for p in pts]
            qm = GR.polygonmoments(m, turned; closed = false)
            @test qm.area ≈ pm.area rtol = 1e-13
            @test collect(qm.moment) ≈ R * collect(pm.moment) rtol = 1e-12
        end
    end

    @testset "closed and open rings, and degenerate ones" begin
        open = moment_graticule(-20.0, 40.0, 10.0, 55.0, 64)
        closed = vcat(open, [open[1]])
        a = GR.polygonmoments(m, open; closed = false)
        b = GR.polygonmoments(m, closed; closed = true)
        @test a === b

        # A clockwise ring encloses the same region, so both fields match.
        rev = GR.polygonmoments(m, reverse(open); closed = false)
        @test rev.area ≈ a.area rtol = 1e-14
        @test collect(rev.moment) ≈ collect(a.moment) rtol = 1e-12

        @test GR.polygonmoments(m, open[1:2]; closed = false) ===
              zero(GR.PolygonMoments)
        @test GR.polygonmoments(m, [open[1], open[1], open[1]]; closed = false) ===
              zero(GR.PolygonMoments)
    end

    @testset "cell moments tile the sphere" begin
        space = ToyLonLatSpace(36, 18)
        total = sum(GR.cellmoments(space, i) for i in 1:ncells(space))
        @test total.area ≈ 4pi rtol = 1e-12
        # A closed surface has no first moment: the cells' moments cancel.
        @test LinearAlgebra.norm(collect(total.moment)) < 1e-10

        # Cell areas are what `Conservative()` measures, to the last bit.
        @test all(GR.cellmoments(space, i).area ===
                  GO.area(manifold(space), getcell(space, i))
                  for i in 1:ncells(space))

        for i in (1, 200, ncells(space) - 1)
            pm = GR.cellmoments(space, i)
            lon_lo, lon_hi, lat_lo, lat_hi = cellbounds(space, i)
            lon, lat = toy_lonlat(collect(pm.moment) ./ pm.area)
            @test lon_lo <= lon <= lon_hi
            @test lat_lo <= lat <= lat_hi
        end
    end

    @testset "IntersectionMomentOperator" begin
        op = GR.IntersectionMomentOperator(m)
        @test op.manifold === m
        @test op.cache isa GO.SutherlandHodgmanCache

        local_op = CR.task_local_operator(op)
        @test local_op isa GR.IntersectionMomentOperator
        @test local_op.manifold === m
        # Each assembly task clips into a cache of its own.
        @test local_op.cache !== op.cache
        @test CR.task_local_operator(local_op).cache !== local_op.cache

        src = ToyLonLatSpace(12, 6)
        dst = ToyLonLatSpace(20, 10)

        @test op(getcell(src, 1), getcell(src, ncells(src))) ===
              zero(GR.PolygonMoments)
        inner = op(getcell(src, 1), getcell(dst, 1))
        @test inner.area > 0

        # The hooks `intersection_areas` needs to assemble over the type.
        @test CR.output_eltype(op) === GR.PolygonMoments
        @test CR.should_store_result(op, inner)
        @test !CR.should_store_result(op, zero(GR.PolygonMoments))
    end

    @testset "intersection_areas assembles a moment matrix" begin
        src = ToyLonLatSpace(12, 6)
        dst = ToyLonLatSpace(20, 10)
        op = GR.IntersectionMomentOperator(m)
        M = CR.intersection_areas(m, GOCore.True(), celltree(dst), celltree(src);
            intersection_operator = op)
        @test M isa GR.SparseMatrixCSC{GR.PolygonMoments,<:Integer}
        @test size(M) == (ncells(dst), ncells(src))
        # Each column is one source cell, tiled by the destination cells.
        for i in (1, 7, 40, ncells(src))
            covered = sum(M[:, i])
            whole = GR.cellmoments(src, i)
            @test covered.area ≈ whole.area rtol = 1e-10
            @test collect(covered.moment) ≈ collect(whole.moment) rtol = 1e-10
        end
    end

    @testset "overlap moments tile the source cell" begin
        src = ToyLonLatSpace(12, 6)
        dst = ToyLonLatSpace(20, 10)
        op = GR.IntersectionMomentOperator(m)
        areaop = GR.IntersectionAreaOperator(m)

        # Offset lattices: no edges coincide, so every source cell is genuinely cut up.
        for i in (1, 7, 40, ncells(src))
            subject = getcell(src, i)
            pieces = GR.PolygonMoments[]
            for j in 1:ncells(dst)
                pm = op(subject, getcell(dst, j))
                pm.area > 0 || continue
                # The clipped area is the area-only operator's, to the last bit.
                @test pm.area === areaop(subject, getcell(dst, j))
                push!(pieces, pm)
            end
            @test length(pieces) > 1
            covered = sum(pieces)
            whole = GR.cellmoments(src, i)
            @test covered.area ≈ whole.area rtol = 1e-10
            @test collect(covered.moment) ≈ collect(whole.moment) rtol = 1e-10
        end
    end
end

@testset "tiny polygons keep their precision" begin
    # A 30 m Copernicus pixel is 5e-6 rad across, and needs its mean position to a
    # small fraction of that. The plain edge sum loses the cell to cancellation; the
    # shifted sum is limited only by the vertices, giving the mean to `eps / h²` of
    # `h`, which is `tol`. The comparison is against the square's planar centroid,
    # which the spherical one matches to O(h²).
    for h in (5e-5, 5e-6), n in (USPoint(1.0, 0.0, 0.0), toy_point(10.0, 46.5))
        e1, e2 = GR.tangentframe(n)
        corner(u, v) = LinearAlgebra.normalize(USPoint((n .+ u .* e1 .+ v .* e2)...))
        square = [corner(-h, -h), corner(h, -h), corner(h, h), corner(-h, h)]
        pm = GR.polygonmoments(GO.Spherical(; radius = 1.0), square; closed = false)
        @test pm.area ≈ 4 * h^2 rtol = 1e-6
        tol = 20 * eps(Float64) / h^2
        mean = USPoint((pm.moment ./ pm.area)...)
        @test all(abs.(GR.tangentcoords(e1, e2, mean - n)) .< tol * h)
        # Shifting the square by a cell moves its mean by exactly that.
        shifted = [corner(2h + s[1], s[2]) for s in ((-h, -h), (h, -h), (h, h), (-h, h))]
        qm = GR.polygonmoments(GO.Spherical(; radius = 1.0), shifted; closed = false)
        d = USPoint((qm.moment ./ qm.area .- pm.moment ./ pm.area)...)
        @test abs(GR.tangentcoords(e1, e2, d)[1] - 2h) < tol * h
        @test abs(GR.tangentcoords(e1, e2, d)[2]) < tol * h
    end
end
