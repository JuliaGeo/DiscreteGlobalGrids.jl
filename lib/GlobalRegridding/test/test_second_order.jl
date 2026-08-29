# Second-order conservative weights on the analytic toy spaces.

import ConservativeRegridding as CR
import SparseArrays

# Cell means of `a + b ⋅ x` over a space's own polygons: exact, from the
# first moments, for a field linear in the ambient coordinates.
so_linear_means(space, a, b) = map(1:ncells(space)) do i
    pm = GR.cellmoments(space, i)
    return a + (b[1] * pm.moment[1] + b[2] * pm.moment[2] + b[3] * pm.moment[3]) / pm.area
end

# Cell means of a smooth `f(p)` by degree-8 fan quadrature over each polygon,
# the rule ConservativeRegridding's own test helpers integrate with.
const SO_RULE = CR.TriangleQuadrature.reference_rule(8)

function so_polygon_integral(f, pts)
    bary, w = SO_RULE
    total = 0.0
    a = pts[1]
    for i in 2:(length(pts)-1)
        b, c = pts[i], pts[i+1]
        det = LinearAlgebra.dot(a, LinearAlgebra.cross(b, c))
        for k in eachindex(w)
            l1, l2, l3 = bary[k]
            p = l1 * a + l2 * b + l3 * c
            np = LinearAlgebra.norm(p)
            total += w[k] * f(p / np) * abs(det) / np^3
        end
    end
    return total
end

function so_means(f, space)
    return map(1:ncells(space)) do i
        pts = [USPoint(GI.x(p), GI.y(p), GI.z(p))
               for p in GI.getpoint(GI.getexterior(getcell(space, i)))]
        pop!(pts)
        return so_polygon_integral(f, pts) / so_polygon_integral(_ -> 1.0, pts)
    end
end

so_smooth(p) = 2 + sin(2 * asin(clamp(p[3], -1, 1)))^2 * cos(3 * atan(p[2], p[1])) +
               0.5 * p[1] * p[3]

so_lat(space, i) = toy_lonlat(cellcentroid(space, i))[2]

so_block(method, dst, src) = GR.pairblock(method, dst, 1:ncells(dst), src, 1:ncells(src))

so_eager(data, dst, src, method; policy = Weighted(0.5)) =
    regrid(data; to = dst, from = src, method, missingpolicy = policy, lazy = false)

# The whole-domain weights, assembled one chunk pair at a time.
function so_assembled(method, dst, src)
    W = zeros(ncells(dst), ncells(src))
    D = zeros(ncells(dst))
    for dc in 1:nchunks(dst), sc in 1:nchunks(src)
        di, si = ownedindices(dst, dc), ownedindices(src, sc)
        block = GR.pairblock(method, dst, di, src, si)
        W[di, si] .+= block.weights
        D[di] .+= block.denom
    end
    return W, D
end

function so_chunked(data, dst, src, method; policy = Weighted(0.5))
    plan = ChunkedPlan(method, policy, dst, src)
    return collect(regrid(data, plan))
end

@testset "ConservativeSecondOrder" begin
    src = ToyLonLatSpace(36, 18; chunks = (9, 6))
    dst = ToyLonLatSpace(40, 20; chunks = (10, 5))
    method = ConservativeSecondOrder()

    @testset "block structure" begin
        block = so_block(method, dst, src)
        first = so_block(Conservative(), dst, src)
        # Coverage is the first-order overlap matrix, entry for entry, and the
        # denominator is its row sums; the value weights are signed and reach
        # one cell further, but sum to the same coverage.
        @test block.coverage == first.weights
        @test block.denom ≈ first.denom
        @test all(>=(0), SparseArrays.nonzeros(block.coverage))
        @test vec(sum(block.weights; dims = 2)) ≈ block.denom atol = 1e-12
        @test SparseArrays.nnz(block.weights) > SparseArrays.nnz(first.weights)
        @test any(<(0), SparseArrays.nonzeros(block.weights))
        @test supportradius(method, src) == GR.celldiameter(src)
        @test GR.outputsampling(method) isa DD.Lookups.Intervals
    end

    @testset "linear fields" begin
        # A field linear in the ambient coordinates is linear in no cell's
        # tangent chart, so it is carried to second order, not exactly: the
        # remaining error is the chart's curvature term, an order below the
        # first-order flattening.
        a, b = 3.0, (0.7, -0.4, 1.1)
        f = so_linear_means(src, a, b)
        exact = so_linear_means(dst, a, b)
        second = maximum(abs.(so_eager(f, dst, src, method) .- exact))
        firstorder = maximum(abs.(so_eager(f, dst, src, Conservative()) .- exact))
        @test second < firstorder / 10
        @test firstorder > 1e-2
    end

    @testset "conserves the covered integral" begin
        f = so_means(so_smooth, src)
        mass = sum(f[i] * GR.cellarea(src, i) for i in 1:ncells(src))
        sums = so_eager(f, dst, src, method; policy = Extensive())
        @test sum(sums) ≈ mass rtol = 1e-12
        means = so_eager(f, dst, src, method)
        @test sum(means[j] * GR.cellarea(dst, j) for j in 1:ncells(dst)) ≈ mass rtol = 1e-10
    end

    @testset "second-order convergence" begin
        target = ToyLonLatSpace(20, 10)
        exact = so_means(so_smooth, target)
        errors = map(((18, 9), (36, 18), (72, 36))) do (nlon, nlat)
            s = ToyLonLatSpace(nlon, nlat)
            f = so_means(so_smooth, s)
            e2 = maximum(abs.(so_eager(f, target, s, method) .- exact))
            e1 = maximum(abs.(so_eager(f, target, s, Conservative()) .- exact))
            return (e1, e2)
        end
        @test all(e[2] < e[1] for e in errors)
        # Halving the source cells quarters the second-order error.
        @test errors[2][2] / errors[3][2] > 3
        @test errors[1][2] / errors[2][2] > 3
    end

    @testset "chunked blocks agree with the whole domain" begin
        # Every chunk pair's block, summed, is the whole-domain block: a chunk's
        # weights read its ring's overlaps, and only its own cells are emitted.
        whole = so_block(method, dst, src)
        W, D = so_assembled(method, dst, src)
        @test W ≈ Matrix(whole.weights) atol = 1e-14
        @test D ≈ whole.denom atol = 1e-14
        # Through a lazy plan, where discovery pairs chunks by cap and support
        # radius: a tile touching only a source chunk's ring still reads it.
        f = so_means(so_smooth, src)
        rows = ToyLonLatSpace(36, 18; chunks = (36, 3))
        tiles = ToyLonLatSpace(40, 20; chunks = (40, 4))
        @test maximum(abs.(so_chunked(f, tiles, rows, method) .-
                           so_eager(f, tiles, rows, method))) < 1e-12
        coarse = ToyLonLatSpace(12, 6; chunks = (12, 1))
        g = so_means(so_smooth, coarse)
        @test maximum(abs.(so_chunked(g, tiles, coarse, method) .-
                           so_eager(g, tiles, coarse, method))) < 1e-12
    end

    @testset "coverage counts valid area, not signed weight" begin
        f = so_means(so_smooth, src)
        hole = copy(f)
        hole[cellat(src, toy_point(0.5, 0.5))] = NaN
        means = so_eager(hole, dst, src, method)
        block = so_block(method, dst, src)
        cover = block.coverage * (1.0 .- isnan.(hole))
        for j in 1:ncells(dst)
            if cover[j] < 0.5 * block.denom[j]
                @test isnan(means[j])
            else
                @test isfinite(means[j])
            end
        end
        # Beyond the hole's own ring — its neighbours' gradients read it as
        # zero, the documented limit of a fixed operator — nothing changes.
        full = so_eager(f, dst, src, method)
        far = [j for j in 1:ncells(dst) if
               US.spherical_distance(cellcentroid(dst, j), toy_point(0.5, 0.5)) > 0.7]
        @test length(far) > ncells(dst) ÷ 2
        @test means[far] ≈ full[far]
        sums = so_eager(hole, dst, src, method; policy = Extensive())
        @test all(isfinite, sums)
    end

    @testset "partial sources and degenerate stencils" begin
        patch = ToyLonLatSpace(12, 6; lat = (0.0, 40.0))
        f = so_means(so_smooth, patch)
        means = so_eager(f, dst, patch, method)
        outside = [j for j in 1:ncells(dst) if so_lat(dst, j) < -10 || so_lat(dst, j) > 50]
        @test all(isnan, means[outside])
        sums = so_eager(f, dst, patch, method; policy = Extensive())
        @test sum(sums) ≈ sum(f[i] * GR.cellarea(patch, i) for i in 1:ncells(patch)) rtol = 1e-12

        # One column of cells: every stencil is collinear, so every gradient
        # is zero and the weights are first order exactly.
        column = ToyLonLatSpace(1, 6; lon = (0.0, 30.0), lat = (-60.0, 60.0))
        @test so_block(method, dst, column).weights == so_block(Conservative(), dst, column).weights
    end
end
