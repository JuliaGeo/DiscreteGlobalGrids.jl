# Cell adjacency, tangent frames, and the least-squares gradient stencil.

import DimensionalData as DD

# Plain 3-vector arithmetic, so the tests never borrow the package's own.
grad_dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
grad_cross3(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])

"""
    grad_spiral(n) -> Vector{UnitSphericalPoint}

`n` well-spread unit vectors, from the golden spiral. Deterministic, so a
failure is reproducible, and unaligned with any axis, so a frame that quietly
depends on one shows up.
"""
function grad_spiral(n::Int)
    out = USPoint[]
    for k in 1:n
        z = 1 - 2 * (k - 0.5) / n
        r = sqrt(max(0.0, 1 - z^2))
        phi = 2pi * k / Base.MathConstants.golden
        push!(out, USPoint(r * cos(phi), r * sin(phi), z))
    end
    return out
end

# Neighbour offsets in the tangent plane: a regular ring, and a scattered set
# that is neither symmetric nor collinear.
grad_ring(m::Int, r::Float64) =
    [(r * cos(2pi * (k - 1) / m), r * sin(2pi * (k - 1) / m)) for k in 1:m]

grad_scatter(m::Int, r::Float64) =
    [(r * cos(2.4k) * (0.4 + 0.6 * mod(0.37k, 1.0)),
        r * sin(2.4k) * (0.5 + 0.5 * mod(0.61k, 1.0))) for k in 1:m]

"""
    grad_neighbourpoints(c, e1, e2, offsets) -> Vector{NTuple{3,Float64}}

The neighbour mean positions an offset set stands for: `c` displaced in its own
tangent plane. Off the sphere by `O(offset²)`, exactly as polygon mean positions
are.
"""
grad_neighbourpoints(c, e1, e2, offsets) =
    [(c[1] + o[1] * e1[1] + o[2] * e2[1], c[2] + o[1] * e1[2] + o[2] * e2[2],
        c[3] + o[1] * e1[3] + o[2] * e2[3]) for o in offsets]

"""
    grad_recover(n, offsets, b; scale = 1.0, a = 0.0) -> (recovered, exact)

The tangent gradient a stencil recovers for the affine field `p -> a + b ⋅ p`,
and the gradient that field actually has. The recovery applies the cell's own
coefficient, `-sum(coeffs)`, so a nonzero `a` fails unless it cancels.
"""
function grad_recover(n, offsets, b; scale::Float64 = 1.0, a::Float64 = 0.0)
    e1, e2 = GR.tangentframe(n)
    c = (scale * n[1], scale * n[2], scale * n[3])
    neighbours = grad_neighbourpoints(c, e1, e2, offsets)
    coeffs = GR.gradientstencil(e1, e2, c, neighbours)
    coeffs === nothing && return nothing
    field(p) = a + grad_dot3(b, p)
    g1 = g2 = 0.0
    for k in eachindex(neighbours)
        g1 += coeffs[k][1] * field(neighbours[k])
        g2 += coeffs[k][2] * field(neighbours[k])
    end
    self = (-sum(w[1] for w in coeffs), -sum(w[2] for w in coeffs))
    g1 += self[1] * field(c)
    g2 += self[2] * field(c)
    return ((g1, g2), GR.tangentcoords(e1, e2, b))
end

# The cell's actual angular diameter: its ring vertices are geodesic corners, so
# the widest pair of them is the widest pair of the cell.
grad_ringpoints(space, i) = [Tuple(p) for p in GI.getpoint(getcell(space, i))]

function grad_celldiameter(space, i)
    ring = grad_ringpoints(space, i)
    d = 0.0
    for p in ring, q in ring
        d = max(d, US.spherical_distance(USPoint(p...), USPoint(q...)))
    end
    return d
end

function grad_sharesvertex(space, i, j; tol = 1e-9)
    a, b = grad_ringpoints(space, i), grad_ringpoints(space, j)
    return any(p -> any(q -> sqrt(sum((p .- q) .^ 2)) <= tol, b), a)
end

@testset "gradient" begin

    @testset "tangent frames" begin
        points = [USPoint(0.0, 0.0, 1.0), USPoint(0.0, 0.0, -1.0),
            USPoint(1.0, 0.0, 0.0), USPoint(-1.0, 0.0, 0.0),
            USPoint(0.0, 1.0, 0.0), USPoint(0.0, -1.0, 0.0),
            toy_point(37.0, -12.5), grad_spiral(64)...]
        for n in points
            e1, e2 = GR.tangentframe(n)
            # Orthonormal, and tangent at `n`.
            @test grad_dot3(e1, e1) ≈ 1.0 atol = 1e-14
            @test grad_dot3(e2, e2) ≈ 1.0 atol = 1e-14
            @test abs(grad_dot3(e1, e2)) < 1e-14
            @test abs(grad_dot3(e1, n)) < 1e-14
            @test abs(grad_dot3(e2, n)) < 1e-14
            # Right-handed about `n`, so a counter-clockwise turn stays one.
            @test all(abs.(grad_cross3(e1, e2) .- Tuple(n)) .< 1e-14)
            # The point itself is the origin of its own frame.
            @test all(abs.(GR.tangentcoords(e1, e2, n)) .< 1e-14)
        end

        # The frame is a property of the point, not of the call.
        @test GR.tangentframe(USPoint(0.0, 0.0, 1.0)) ==
              GR.tangentframe(USPoint(0.0, 0.0, 1.0))
    end

    @testset "gradient stencils" begin
        b = (0.6, -1.3, 0.45)
        for n in [USPoint(0.0, 0.0, 1.0), toy_point(-72.5, 41.0),
                grad_spiral(16)...]
            for offsets in [grad_ring(3, 0.01), grad_ring(4, 0.01),
                    grad_ring(6, 0.004), grad_ring(8, 0.01),
                    grad_scatter(5, 0.01), grad_scatter(9, 0.002)]
                got = grad_recover(n, offsets, b)
                @test got !== nothing
                recovered, exact = got
                @test recovered[1] ≈ exact[1] atol = 1e-12
                @test recovered[2] ≈ exact[2] atol = 1e-12
            end
        end

        # A mean position is not a point of the sphere, and a field has an
        # offset the self coefficient must cancel.
        got = grad_recover(toy_point(15.0, -60.0), grad_ring(4, 0.01), b;
            scale = 0.87, a = 7.5)
        @test got !== nothing
        @test all(abs.(got[1] .- got[2]) .< 1e-11)

        # A constant field has no gradient at all, exactly.
        e1, e2 = GR.tangentframe(USPoint(0.0, 0.0, 1.0))
        c = (0.0, 0.0, 1.0)
        neighbours = grad_neighbourpoints(c, e1, e2, grad_ring(4, 0.01))
        coeffs = GR.gradientstencil(e1, e2, c, neighbours)
        @test coeffs !== nothing
        self = (-sum(w[1] for w in coeffs), -sum(w[2] for w in coeffs))
        @test sum(w[1] for w in coeffs) + self[1] == 0.0
        @test sum(w[2] for w in coeffs) + self[2] == 0.0
        # A symmetric ring leans nowhere.
        @test all(abs.(self) .< 1e-12)

        # Rank-deficient stencils are refused rather than solved.
        buffer = NTuple{2,Float64}[]
        @test !GR.gradientstencil!(buffer, e1, e2, c, NTuple{3,Float64}[])
        @test isempty(buffer)
        one_neighbour = grad_neighbourpoints(c, e1, e2, [(0.01, 0.0)])
        @test !GR.gradientstencil!(buffer, e1, e2, c, one_neighbour)
        @test isempty(buffer)
        collinear = grad_neighbourpoints(c, e1, e2,
            [(0.01, 0.0), (-0.01, 0.0), (0.005, 0.0)])
        @test !GR.gradientstencil!(buffer, e1, e2, c, collinear)
        @test isempty(buffer)
        @test GR.gradientstencil(e1, e2, c, collinear) === nothing
        # Neighbours all at the cell's own position are collinear too.
        @test !GR.gradientstencil!(buffer, e1, e2, c, [c, c, c])

        # One buffer serves a sweep: it is resized, and a refusal empties it.
        @test GR.gradientstencil!(buffer, e1, e2, c, neighbours)
        @test length(buffer) == length(neighbours)
        @test buffer == coeffs
        @test GR.gradientstencil!(buffer, e1, e2, c,
            grad_neighbourpoints(c, e1, e2, grad_ring(6, 0.01)))
        @test length(buffer) == 6
        @test !GR.gradientstencil!(buffer, e1, e2, c, collinear)
        @test isempty(buffer)
    end

    @testset "lattice adjacency" begin
        space = ToyLonLatSpace(8, 4)

        # An interior cell has its four edge neighbours.
        i = localindex(space, 3, 2)
        @test sort(GR.cellneighbors(space, i)) == sort([localindex(space, 2, 2),
            localindex(space, 4, 2), localindex(space, 3, 1),
            localindex(space, 3, 3)])

        # Longitude wraps across the dateline in a global space.
        @test localindex(space, 8, 2) in GR.cellneighbors(space, localindex(space, 1, 2))
        @test localindex(space, 1, 2) in GR.cellneighbors(space, localindex(space, 8, 2))

        # Latitude never wraps: a polar row has three neighbours.
        @test length(GR.cellneighbors(space, localindex(space, 5, 1))) == 3
        @test length(GR.cellneighbors(space, localindex(space, 5, 4))) == 3
        @test localindex(space, 5, 4) ∉ GR.cellneighbors(space, localindex(space, 5, 1))

        # No cell is its own neighbour, and every answer is a cell.
        @test all(i -> i ∉ GR.cellneighbors(space, i), 1:ncells(space))
        @test all(j -> 1 <= j <= ncells(space),
            reduce(vcat, GR.cellneighbors(space, i) for i in 1:ncells(space)))

        # A regional space has a rim instead of a wrap.
        patch = ToyLonLatSpace(4, 3; lon = (-40.0, 40.0), lat = (-20.0, 20.0))
        @test length(GR.cellneighbors(patch, localindex(patch, 1, 1))) == 2
        @test length(GR.cellneighbors(patch, localindex(patch, 4, 3))) == 2
        @test length(GR.cellneighbors(patch, localindex(patch, 2, 2))) == 4
        @test localindex(patch, 4, 1) ∉ GR.cellneighbors(patch, localindex(patch, 1, 1))
    end

    @testset "geometric adjacency" begin
        space = ToyLonLatSpace(8, 4)
        # `CountingSpace` forwards geometry and not topology, so it reaches the
        # generic fallback the way an unstructured space does.
        counting = CountingSpace(space)

        for i in 1:ncells(space)
            generic = GR.cellneighbors(counting, i)
            @test generic == invoke(GR.cellneighbors, Tuple{RegridSpace,Int}, space, i)
            @test i ∉ generic
            # A superset of the lattice answer, and every cell of it touches.
            @test GR.cellneighbors(space, i) ⊆ generic
            @test all(j -> grad_sharesvertex(space, i, j), generic)
        end

        # In the interior that is the eight-cell ring, diagonals included.
        i = localindex(space, 3, 2)
        @test length(GR.cellneighbors(counting, i)) == 8
        @test localindex(space, 2, 1) in GR.cellneighbors(counting, i)

        # Every cell of a polar row meets every other one at the pole.
        top = [localindex(space, ix, 4) for ix in 1:8]
        polar = GR.cellneighbors(counting, localindex(space, 3, 4))
        @test setdiff(top, [localindex(space, 3, 4)]) ⊆ polar
    end

    @testset "raster adjacency" begin
        space = RasterGrid(DD.DimArray(zeros(720, 360),
            (DD.X(-179.75:0.5:179.75), DD.Y(-89.75:0.5:89.75))))

        i = GR.localindex(space, 400, 200)
        @test sort(GR.cellneighbors(space, i)) == sort([
            GR.localindex(space, 399, 200), GR.localindex(space, 401, 200),
            GR.localindex(space, 400, 199), GR.localindex(space, 400, 201)])

        # A raster spanning its full period wraps in X and never in Y.
        @test GR.localindex(space, 720, 17) in
              GR.cellneighbors(space, GR.localindex(space, 1, 17))
        @test GR.localindex(space, 1, 17) in
              GR.cellneighbors(space, GR.localindex(space, 720, 17))
        @test length(GR.cellneighbors(space, GR.localindex(space, 1, 1))) == 3
        @test length(GR.cellneighbors(space, GR.localindex(space, 360, 360))) == 3

        # A regional raster does not wrap.
        patch = RasterGrid(DD.DimArray(zeros(4, 3),
            (DD.X(5.0:10.0:35.0), DD.Y(5.0:10.0:25.0))))
        @test GR.chartperiod(patch)[1] === nothing
        @test length(GR.cellneighbors(patch, GR.localindex(patch, 1, 1))) == 2
        @test sort(GR.cellneighbors(patch, GR.localindex(patch, 2, 2))) ==
              sort([GR.localindex(patch, 1, 2), GR.localindex(patch, 3, 2),
            GR.localindex(patch, 2, 1), GR.localindex(patch, 2, 3)])

        # A transposed raster answers in its own index order.
        transposed = RasterGrid(DD.DimArray(zeros(18, 36),
            (DD.Y(-85.0:10.0:85.0), DD.X(-175.0:10.0:175.0))))
        @test sort(GR.cellneighbors(transposed, GR.localindex(transposed, 4, 9))) ==
              sort([GR.localindex(transposed, 3, 9), GR.localindex(transposed, 5, 9),
            GR.localindex(transposed, 4, 8), GR.localindex(transposed, 4, 10)])

        # The lattice answer agrees with the geometry the fallback reads.
        coarse = RasterGrid(DD.DimArray(zeros(36, 18),
            (DD.X(-175.0:10.0:175.0), DD.Y(-85.0:10.0:85.0))))
        for (ix, iy) in [(1, 1), (5, 9), (36, 9), (18, 18), (12, 2)]
            j = GR.localindex(coarse, ix, iy)
            generic = invoke(GR.cellneighbors, Tuple{RegridSpace,Int}, coarse, j)
            @test GR.cellneighbors(coarse, j) ⊆ generic
            @test all(k -> grad_sharesvertex(coarse, j, k), generic)
        end
    end

    @testset "cell diameters" begin
        space = RasterGrid(DD.DimArray(zeros(720, 360),
            (DD.X(-179.75:0.5:179.75), DD.Y(-89.75:0.5:89.75))))
        bound = GR.celldiameter(space)
        @test bound ≈ 2 * deg2rad(0.5)

        # Sampled equatorial, mid-latitude and polar cells: the bound covers
        # every one of them, and is the right size for the widest.
        sampled = [GR.localindex(space, ix, iy)
                   for ix in (1, 97, 360, 719) for iy in (1, 2, 90, 180, 181, 359, 360)]
        widest = maximum(grad_celldiameter(space, i) for i in sampled)
        @test all(grad_celldiameter(space, i) <= bound for i in sampled)
        @test bound <= 3 * widest

        # The toy space answers from its resolution; the generic fallback reads
        # the same cells' caps. Both bound every cell, and the cap bound is the
        # tighter of the two.
        toy = ToyLonLatSpace(8, 4)
        analytic = GR.celldiameter(toy)
        fallback = invoke(GR.celldiameter, Tuple{RegridSpace}, toy)
        @test fallback == GR.celldiameter(CountingSpace(toy))
        largest = maximum(grad_celldiameter(toy, i) for i in 1:ncells(toy))
        @test analytic ≈ 2 * deg2rad(45.0)
        @test analytic >= largest
        @test fallback >= largest
        @test fallback <= 2.5 * largest

        # A regional toy space, and one whose cells span more than a hemisphere.
        patch = ToyLonLatSpace(4, 3; lon = (-40.0, 40.0), lat = (-20.0, 20.0))
        @test GR.celldiameter(patch) >=
              maximum(grad_celldiameter(patch, i) for i in 1:ncells(patch))
        @test GR.celldiameter(ToyLonLatSpace(1, 1)) == Float64(pi)
    end
end
