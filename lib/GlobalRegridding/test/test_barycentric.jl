import DimensionalData as DD

# The point contract: one reusable weight row, the dual cell a point falls in,
# and the two coordinate kernels that weight its nodes. The kernels are tested
# on bare node coordinates, with no space and no chunk anywhere near them; the
# sampler and the chunk-pair builder are tested on a toy space that answers one
# hand-built dual cell.

# Source indices that are not node slot numbers, so a kernel that emits `k`
# instead of `indices[k]` is caught rather than silently right.
const P1_TRI_INDS = [7, 3, 11]
const P1_TRI_NODES = [(0.0, 0.0), (4.0, 0.0), (1.0, 3.0)]

const P1_RECT_INDS = [2, 9, 5, 8]
const P1_RECT_NODES = [(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (0.0, 2.0)]

const P1_QUAD_INDS = [12, 4, 6, 1]
const P1_QUAD_NODES = [(0.0, 0.0), (5.0, 0.0), (6.0, 3.0), (1.0, 4.0)]

const P1_HEX_INDS = collect(21:26)
const P1_HEX_NODES = [(cosd(60.0 * (k - 1)), sind(60.0 * (k - 1))) for k in 1:6]

const P1_PENT_INDS = [3, 14, 15, 9, 2]
const P1_PENT_NODES = [(0.0, 0.0), (4.0, 0.0), (5.0, 2.5), (2.0, 4.0), (-1.0, 2.0)]

# Fields the kernels promise: a constant, an affine field, and the bilinear one
# only tensor Q1 promises.
p1_constant(_) = 3.5
p1_affine(c) = 5.0 - 0.75 * c[1] + 2.5 * c[2]
p1_bilinear(c) = 1.0 + 0.5 * c[1] - 2.0 * c[2] + 0.25 * c[1] * c[2]

"""
    p1_reproduce(row, indices, nodes, f) -> Float64

The row's value for the field `f` sampled at the cell's nodes. Every entry must
name one of `indices`, so a row that emits a node's slot number, or a source
cell the cell never held, fails here rather than in a value comparison.
"""
function p1_reproduce(row, indices, nodes, f)
    total = 0.0
    for k in 1:length(row)
        j = findfirst(==(row.indices[k]), indices)
        j === nothing &&
            error("the row names source cell $(row.indices[k]), which is not a node")
        total += row.weights[k] * f(nodes[j])
    end
    return total
end

# Where the row says the value sits: the coordinate field, which every basis
# here reproduces.
p1_site(row, indices, nodes) = (p1_reproduce(row, indices, nodes, c -> c[1]),
    p1_reproduce(row, indices, nodes, c -> c[2]))

p1_rowsum(row) = sum(row.weights; init = 0.0)

p1_weightof(row, i) = sum((row.weights[k] for k in 1:length(row) if row.indices[k] == i);
    init = 0.0)

"""
    p1_barycentric(nodes, p) -> NTuple{3,Float64}

The barycentric coordinates of `p` in the triangle `nodes`, from the three
signed areas. This is the closed form the mean-value kernel must reproduce on a
triangle, computed here rather than taken from the kernel under test.
"""
function p1_barycentric(nodes, p)
    a, b, c = nodes
    twicearea(u, v, w) = (v[1] - u[1]) * (w[2] - u[2]) - (v[2] - u[2]) * (w[1] - u[1])
    total = twicearea(a, b, c)
    return (twicearea(p, b, c) / total, twicearea(p, c, a) / total,
        twicearea(p, a, b) / total)
end

# A convex ring: points in angular order on a circle, so the ring is convex, put
# through an affine map, which keeps it convex while giving it a shape. Both
# arguments are deterministic, so the sweep is the same on every run.
function p1_convexring(n::Int, k::Int)
    θ = [2π * (j - 1 + 0.4 * sinpi((j * k) / 7)) / n for j in 1:n]
    a, b = 1.0 + 0.3 * cospi(k / 5), 0.4 * sinpi(k / 3)
    c, d = 0.2 * cospi(k / 4), 0.8 + 0.5 * sinpi(k / 6)
    ox, oy = 2.0 * cospi(k / 9), 3.0 * sinpi(k / 11)
    return [(a * cos(t) + b * sin(t) + ox, c * cos(t) + d * sin(t) + oy) for t in θ]
end

# A point strictly inside a convex ring: a convex combination of its nodes.
function p1_inside(nodes, k::Int)
    w = [0.2 + abs(sinpi((j * k) / 13)) for j in eachindex(nodes)]
    w ./= sum(w)
    return (sum(w[j] * nodes[j][1] for j in eachindex(nodes)),
        sum(w[j] * nodes[j][2] for j in eachindex(nodes)))
end

# Warm the row, then measure the steady-state call.
p1_kernelbytes(kernel, row, indices, nodes, p) =
    (kernel(row, indices, nodes, p); @allocated kernel(row, indices, nodes, p))

p1_samplerbytes(row, smp, p) =
    (GR.weightsat!(row, smp, p); @allocated GR.weightsat!(row, smp, p))

# A source space whose sampler answers one hand-built dual cell: the
# quadrilateral of its own four sample sites, held in two chunks so every
# stencil crosses a source-chunk seam.
struct P1DualSpace <: RegridSpace
    space::ToyLonLatSpace
    cell::GR.DualCell{Vector{Int},Vector{NTuple{2,Float64}}}
end

function P1DualSpace(space::ToyLonLatSpace, order::Vector{Int}, kind::GR.BasisKind)
    nodes = [Tuple{Float64,Float64}(toy_lonlat(cellcentroid(space, i))) for i in order]
    return P1DualSpace(space, GR.DualCell(order, nodes, kind))
end

ncells(s::P1DualSpace) = ncells(s.space)
getcell(s::P1DualSpace, i::Int) = getcell(s.space, i)
manifold(s::P1DualSpace) = manifold(s.space)
cellcentroid(s::P1DualSpace, i::Int) = cellcentroid(s.space, i)
cellat(s::P1DualSpace, p) = cellat(s.space, p)
nchunks(s::P1DualSpace) = nchunks(s.space)
ownedindices(s::P1DualSpace, chunk::Int) = ownedindices(s.space, chunk)
celltree(s::P1DualSpace) = celltree(s.space)
chunkextents(s::P1DualSpace) = chunkextents(s.space)

# The chart the hand-built cell is written in.
hascellchart(::P1DualSpace) = true
GR.chartcoords(::P1DualSpace, p) = toy_lonlat(p)

GR.samplerstate(s::P1DualSpace) = s.cell
GR.dualcellat(smp::GR.Sampler{<:P1DualSpace}, p) = smp.state

# A source with no cell chart, so its sampler prepares no state and it takes the
# generic dual-cell path with no cells of its own to offer.
struct P1PlainSpace <: RegridSpace
    space::ToyLonLatSpace
end

ncells(s::P1PlainSpace) = ncells(s.space)
cellcentroid(s::P1PlainSpace, i::Int) = cellcentroid(s.space, i)

# One operator from every chunk pair, in the spaces' own local indices.
function p1_assemble(plan, dst, src)
    M = zeros(Float64, Int(ncells(dst)), Int(ncells(src)))
    for d in 1:nchunks(dst), s in 1:nchunks(src)
        block = Matrix(GR.buildblock(plan, d, s).weights)
        dinds, sinds = ownedindices(dst, d), ownedindices(src, s)
        for (r, i) in enumerate(dinds), (c, j) in enumerate(sinds)
            M[Int(i), Int(j)] += block[r, c]
        end
    end
    return M
end

@testset "barycentric point contracts" begin

    @testset "a triangle is a mean-value cell" begin
        # A triangle has no kernel of its own: mean-value coordinates on three
        # nodes are the triangle's barycentric coordinates, which the test
        # computes from signed areas.
        row = GR.WeightRow()
        inds, nodes = P1_TRI_INDS, P1_TRI_NODES
        λ(p) = p1_barycentric(nodes, p)

        # An interior point takes all three nodes, at the closed-form
        # coordinates and not merely at some partition of unity.
        @test GR.meanvalueweights!(row, inds, nodes, (1.5, 1.5)) === GR.WeightsMapped
        @test length(row) == 3
        @test sort(row.indices) == sort(inds)
        @test p1_rowsum(row) ≈ 1.0
        @test all(isapprox(p1_weightof(row, inds[k]), λ((1.5, 1.5))[k]; atol = 1e-12)
                  for k in 1:3)
        @test p1_weightof(row, 11) ≈ 0.5

        # The coordinates at several interior points, and the fields they carry.
        for p in ((1.0, 1.0), (2.0, 0.5), (0.5, 0.5), (1.2, 2.0))
            @test GR.meanvalueweights!(row, inds, nodes, p) === GR.WeightsMapped
            @test all(isapprox(p1_weightof(row, inds[k]), λ(p)[k]; atol = 1e-12)
                      for k in 1:3)
            @test p1_reproduce(row, inds, nodes, p1_affine) ≈ p1_affine(p)
            @test p1_rowsum(row) ≈ 1.0
        end

        # A point on an edge emits that edge's two nodes, at the two
        # coordinates the closed form gives them, and drops the zero.
        @test GR.meanvalueweights!(row, inds, nodes, (2.0, 0.0)) === GR.WeightsMapped
        @test length(row) == 2
        @test sort(row.indices) == [3, 7]
        @test all(isapprox(p1_weightof(row, inds[k]), λ((2.0, 0.0))[k]; atol = 1e-12)
                  for k in 1:3)

        # A point at a node reproduces that node exactly, with one entry.
        for k in 1:3
            @test GR.meanvalueweights!(row, inds, nodes, nodes[k]) === GR.WeightsMapped
            @test row.indices == [inds[k]]
            @test row.weights == [1.0]
            @test λ(nodes[k])[k] ≈ 1.0
        end

        # Outside is unmapped, never clamped to the boundary, and leaves no row;
        # the closed form is the reason, one coordinate there being negative.
        for p in ((5.0, 0.0), (-1.0, -1.0), (2.0, 3.0), (1.0, -0.5))
            @test GR.meanvalueweights!(row, inds, nodes, p) === GR.WeightsOutside
            @test isempty(row)
            @test minimum(λ(p)) < 0.0
        end

        # Degeneracies are rejected rather than answered with a division by zero.
        @test GR.meanvalueweights!(row, inds, [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)],
            (1.0, 1.0)) === GR.WeightsDegenerate
        @test GR.meanvalueweights!(row, inds, [(0.0, 0.0), (4.0, 0.0), (0.0, 0.0)],
            (1.0, 0.0)) === GR.WeightsDegenerate
        @test GR.meanvalueweights!(row, inds, [(0.0, 0.0), (4.0, 0.0)],
            (1.0, 0.0)) === GR.WeightsDegenerate
        @test isempty(row)
    end

    @testset "inverse bilinear coordinates" begin
        row = GR.WeightRow()

        # On an axis-aligned rectangle Q1 reproduces `a + bx + cy + dxy`, which
        # no affine basis on the same four nodes does.
        inds, nodes = P1_RECT_INDS, P1_RECT_NODES
        for p in ((1.0, 0.5), (3.0, 1.5), (2.0, 1.0), (0.5, 1.75))
            @test GR.bilinearweights!(row, inds, nodes, p) === GR.WeightsMapped
            @test length(row) == 4
            @test p1_rowsum(row) ≈ 1.0
            @test p1_reproduce(row, inds, nodes, p1_bilinear) ≈ p1_bilinear(p)
            @test p1_reproduce(row, inds, nodes, p1_constant) ≈ 3.5
        end

        # On a general convex quadrilateral the inverse map converges: the
        # recovered coordinates, read back off the weights, map forward to `p`.
        inds, nodes = P1_QUAD_INDS, P1_QUAD_NODES
        for p in ((2.0, 1.0), (3.0, 2.0), (1.0, 2.0), (4.0, 2.0))
            @test GR.bilinearweights!(row, inds, nodes, p) === GR.WeightsMapped
            u = p1_weightof(row, inds[2]) + p1_weightof(row, inds[3])
            v = p1_weightof(row, inds[3]) + p1_weightof(row, inds[4])
            image = GR._q1image(nodes[1], nodes[2], nodes[3], nodes[4], u, v)
            @test all(isapprox.(image, p; atol = 1e-10))
            @test p1_rowsum(row) ≈ 1.0
            @test p1_reproduce(row, inds, nodes, p1_affine) ≈ p1_affine(p)
            @test all(isapprox.(p1_site(row, inds, nodes), p; atol = 1e-10))
        end

        # Nodes and edges emit only the nodes that carry them.
        for k in 1:4
            @test GR.bilinearweights!(row, inds, nodes, nodes[k]) === GR.WeightsMapped
            @test row.indices == [inds[k]]
            @test row.weights == [1.0]
        end
        edge = ((nodes[1][1] + nodes[2][1]) / 2, (nodes[1][2] + nodes[2][2]) / 2)
        @test GR.bilinearweights!(row, inds, nodes, edge) === GR.WeightsMapped
        @test sort(row.indices) == sort(inds[1:2])
        @test p1_rowsum(row) ≈ 1.0

        # Outside the quadrilateral is unmapped, on either side of every edge.
        for p in ((-1.0, 0.0), (7.0, 3.0), (2.0, -0.5), (0.0, 5.0))
            @test GR.bilinearweights!(row, inds, nodes, p) === GR.WeightsOutside
            @test isempty(row)
        end

        # A folded, flat, repeated-node or wrong-sized cell is rejected.
        @test GR.bilinearweights!(row, inds,
            [(0.0, 0.0), (4.0, 0.0), (0.0, 3.0), (4.0, 3.0)],
            (2.0, 1.5)) === GR.WeightsDegenerate
        @test GR.bilinearweights!(row, inds,
            [(0.0, 0.0), (2.0, 0.0), (4.0, 0.0), (6.0, 0.0)],
            (3.0, 0.0)) === GR.WeightsDegenerate
        @test GR.bilinearweights!(row, inds,
            [(0.0, 0.0), (4.0, 0.0), (4.0, 0.0), (0.0, 3.0)],
            (1.0, 1.0)) === GR.WeightsDegenerate
        @test GR.bilinearweights!(row, inds, P1_TRI_NODES,
            (1.0, 1.0)) === GR.WeightsDegenerate
        @test isempty(row)
    end

    @testset "mean-value coordinates" begin
        row = GR.WeightRow()

        # A regular hexagon: every node carries weight, and constants, affine
        # fields and the coordinate field come back.
        inds, nodes = P1_HEX_INDS, P1_HEX_NODES
        for p in ((0.0, 0.0), (0.3, 0.2), (-0.4, 0.1), (0.1, -0.5))
            @test GR.meanvalueweights!(row, inds, nodes, p) === GR.WeightsMapped
            @test length(row) == 6
            @test all(>(0.0), row.weights)
            @test p1_rowsum(row) ≈ 1.0
            @test p1_reproduce(row, inds, nodes, p1_constant) ≈ 3.5
            @test p1_reproduce(row, inds, nodes, p1_affine) ≈ p1_affine(p)
            @test all(isapprox.(p1_site(row, inds, nodes), p; atol = 1e-12))
        end
        # The centre of a regular hexagon weights its nodes equally.
        @test GR.meanvalueweights!(row, inds, nodes, (0.0, 0.0)) === GR.WeightsMapped
        @test all(≈(1 / 6), row.weights)

        # An irregular convex pentagon.
        inds, nodes = P1_PENT_INDS, P1_PENT_NODES
        for p in ((2.0, 1.5), (1.0, 1.0), (3.0, 2.0))
            @test GR.meanvalueweights!(row, inds, nodes, p) === GR.WeightsMapped
            @test p1_rowsum(row) ≈ 1.0
            @test p1_reproduce(row, inds, nodes, p1_affine) ≈ p1_affine(p)
            @test all(isapprox.(p1_site(row, inds, nodes), p; atol = 1e-12))
        end

        # A node, and a point on an edge, emit only the nodes that carry them.
        for k in 1:5
            @test GR.meanvalueweights!(row, inds, nodes, nodes[k]) === GR.WeightsMapped
            @test row.indices == [inds[k]]
            @test row.weights == [1.0]
        end
        @test GR.meanvalueweights!(row, inds, nodes, (2.0, 0.0)) === GR.WeightsMapped
        @test sort(row.indices) == sort(inds[1:2])
        @test p1_reproduce(row, inds, nodes, p1_affine) ≈ p1_affine((2.0, 0.0))

        # Outside is unmapped, including a point on an edge's line beyond it.
        for p in ((6.0, 0.0), (-2.0, 0.0), (2.0, 5.0), (5.0, 0.0))
            @test GR.meanvalueweights!(row, inds, nodes, p) === GR.WeightsOutside
            @test isempty(row)
        end

        # A reflex corner, a straight corner, a repeated node and too few nodes
        # are all rejected: the basis holds only on a simple convex polygon.
        @test GR.meanvalueweights!(row, inds,
            [(0.0, 0.0), (4.0, 0.0), (2.0, 1.5), (4.0, 4.0), (0.0, 4.0)],
            (1.0, 2.0)) === GR.WeightsDegenerate
        @test GR.meanvalueweights!(row, inds,
            [(0.0, 0.0), (2.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
            (2.0, 2.0)) === GR.WeightsDegenerate
        @test GR.meanvalueweights!(row, inds,
            [(0.0, 0.0), (4.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)],
            (2.0, 2.0)) === GR.WeightsDegenerate
        @test GR.meanvalueweights!(row, inds, [(0.0, 0.0), (4.0, 0.0)],
            (2.0, 0.0)) === GR.WeightsDegenerate
        @test isempty(row)
    end

    @testset "mean value on convex cells of every node count" begin
        # A sweep of convex rings, three to nine nodes, each queried at three
        # interior points: the two laws that hold on every convex cell, and the
        # reading that does not depend on which way the ring turns.
        row, flipped = GR.WeightRow(), GR.WeightRow()
        mapped = positive = turns = true
        sumerr = siteerr = fliperr = 0.0
        for n in 3:9, k in 1:8
            nodes = p1_convexring(n, k)
            inds = collect(10:(9+n))
            for j in 1:3
                p = p1_inside(nodes, k + j)
                st = GR.meanvalueweights!(row, inds, nodes, p)
                mapped &= st === GR.WeightsMapped && length(row) == n
                st === GR.WeightsMapped && length(row) == n || continue
                positive &= all(>(0.0), row.weights)
                # Partition of unity, and linear reproduction: the weights put
                # the point back where it was asked for.
                sumerr = max(sumerr, abs(p1_rowsum(row) - 1.0))
                siteerr = max(siteerr, maximum(abs.(p1_site(row, inds, nodes) .- p)))
                # Reversing the ring reverses the row and changes nothing else.
                turns &= GR.meanvalueweights!(flipped, reverse(inds), reverse(nodes),
                    p) === GR.WeightsMapped && length(flipped) == n
                length(flipped) == n || continue
                fliperr = max(fliperr,
                    maximum(abs.(reverse(flipped.weights) .- row.weights)))
            end
        end
        @test mapped
        @test positive
        @test turns
        @test sumerr < 1e-12
        @test siteerr < 1e-12
        @test fliperr < 1e-12
    end

    @testset "the row is reused and costs nothing warm" begin
        row = GR.WeightRow()

        # The row is cleared on entry, so a mapped point never inherits the
        # previous point's entries and an unmapped one leaves none.
        GR.meanvalueweights!(row, P1_TRI_INDS, P1_TRI_NODES, (1.0, 1.0))
        @test length(row) == 3
        GR.meanvalueweights!(row, P1_TRI_INDS, P1_TRI_NODES, (2.0, 0.0))
        @test length(row) == 2
        GR.meanvalueweights!(row, P1_TRI_INDS, P1_TRI_NODES, (9.0, 9.0))
        @test isempty(row)

        # Once warm, a point costs no allocation at all.
        @test p1_kernelbytes(GR.meanvalueweights!, row, P1_TRI_INDS, P1_TRI_NODES,
            (1.0, 1.0)) == 0 skip = VERSION < v"1.12"
        @test p1_kernelbytes(GR.bilinearweights!, row, P1_QUAD_INDS, P1_QUAD_NODES,
            (2.0, 1.0)) == 0 skip = VERSION < v"1.12"
        @test p1_kernelbytes(GR.meanvalueweights!, row, P1_HEX_INDS, P1_HEX_NODES,
            (0.3, 0.2)) == 0 skip = VERSION < v"1.12"
    end

    @testset "the sampler prepares a source space once" begin
        space = P1PlainSpace(ToyLonLatSpace(4, 2; lon = (-40.0, 40.0),
            lat = (-20.0, 20.0)))

        # Sample sites are the space's centroids, read on demand.
        sites = GR.samplesites(space)
        @test length(sites) == Int(ncells(space))
        @test all(sites[i] == cellcentroid(space, i) for i in 1:Int(ncells(space)))

        # A sampler holds the space and its sites, and is built once. A source
        # with no cell chart prepares no state.
        smp = GR.sampler(BarycentricPoint(), space)
        @test smp.space === space
        @test smp.state === nothing

        # With no dual cells to offer, a generic source maps nothing rather than
        # inventing a stencil.
        row = GR.WeightRow()
        @test GR.nodecount(GR.dualcellat(smp, cellcentroid(space, 1))) == 0
        @test all(GR.weightsat!(row, smp, cellcentroid(space, i)) === GR.WeightsOutside
                  for i in 1:Int(ncells(space)))
        @test isempty(row)

        # The method reports point sampling, so the shared builders put it on the
        # point path, and only its name is public.
        @test GR.outputsampling(BarycentricPoint()) === DD.Lookups.Points()
        @test Base.ispublic(GR, :BarycentricPoint)
        @test !any(Base.ispublic(GR, n) for n in
                   (:WeightRow, :weightsat!, :sampler, :samplesites, :dualcellat,
            :DualCell, :bilinearweights!, :meanvalueweights!))
    end

    @testset "a source space's dual cell reaches the executor" begin
        # Four source cells in two chunks, whose sample sites are the corners of
        # the one dual cell, and destinations inside it.
        src = P1DualSpace(ToyLonLatSpace(2, 2; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
                chunks = (2, 1)), [1, 2, 4, 3], GR.Bilinear)
        dst = ToyLonLatSpace(3, 2; lon = (-15.0, 15.0), lat = (-6.0, 6.0),
            chunks = (3, 1))
        @test nchunks(src) == 2 && nchunks(dst) == 2

        # The prepared cell answers every point, and the row it fills is the
        # source space's local indices.
        smp = GR.sampler(BarycentricPoint(), src)
        row = GR.WeightRow()
        @test GR.nodecount(GR.dualcellat(smp, cellcentroid(dst, 1))) == 4
        @test all(GR.ismapped(GR.weightsat!(row, smp, cellcentroid(dst, i)))
                  for i in 1:Int(ncells(dst)))
        @test sort(row.indices) == [1, 2, 3, 4]
        @test p1_rowsum(row) ≈ 1.0
        @test p1_samplerbytes(row, smp, cellcentroid(dst, 1)) == 0 skip =
            VERSION < v"1.12"

        # Eager and chunked build one operator, with every stencil split across
        # both source chunks.
        whole = Matrix(GR.wholeblock(BarycentricPoint(), dst, src).weights)
        plan = ChunkedPlan(BarycentricPoint(), Weighted(0.5), dst, src;
            storage = PerChunk())
        @test p1_assemble(plan, dst, src) ≈ whole
        @test all(!iszero(Matrix(GR.buildblock(plan, d, s).weights))
                  for d in 1:nchunks(dst), s in 1:nchunks(src))
        @test all(≈(1.0), sum(whole; dims = 2))

        # The operator reproduces the bilinear field its Q1 cell promises.
        chart(i) = toy_lonlat(cellcentroid(src, i))
        field = [p1_bilinear(chart(i)) for i in 1:Int(ncells(src))]
        @test whole * field ≈
              [p1_bilinear(toy_lonlat(cellcentroid(dst, i))) for i in 1:Int(ncells(dst))]

        # And the values agree eagerly and by chunk.
        eager = regrid(field; to = dst, from = src, method = BarycentricPoint(),
            lazy = false)
        chunked = Vector{Float64}(undef, Int(ncells(dst)))
        regrid!(chunked, field, plan)
        @test chunked == eager
        @test collect(eager) ≈ whole * field

        # A destination the dual cell does not contain takes no weights, so the
        # missing policy decides it rather than an extrapolation.
        far = ToyLonLatSpace(2, 1; lon = (60.0, 100.0), lat = (40.0, 60.0))
        @test isempty(GR.wholeblock(BarycentricPoint(), far, src).weights.nzval)
    end

    @testset "chart Q1 matches BilinearPoint inside the lattice" begin
        # Interior fixture: source centres 20 degrees apart, every
        # destination centroid strictly inside the lattice hull.
        qsrc = ToyLonLatSpace(6, 5; lon = (-60.0, 60.0), lat = (-50.0, 50.0))
        qdst = ToyLonLatSpace(7, 4; lon = (-49.0, 49.0), lat = (-39.0, 39.0))

        # A source with a cell chart prepares its two axes and takes the fused
        # path, whatever kind of space it is.
        smp = GR.sampler(BarycentricPoint(), qsrc)
        @test smp.state isa GR.ChartState
        @test (smp.state.x.n, smp.state.y.n) == (6, 5)

        # Entry for entry the same operator as the bilinear point writes: on a
        # chart source the two are one tensor Q1 stencil.
        bary = GR.wholeblock(BarycentricPoint(), qdst, qsrc).weights
        bilin = GR.wholeblock(BilinearPoint(), qdst, qsrc).weights
        @test Matrix(bary) == Matrix(bilin)

        # Four entries per interior destination, summing to one.
        dense = Matrix(bary)
        @test all(count(!iszero, view(dense, r, :)) == 4 for r in axes(dense, 1))
        @test all(isapprox(1.0), sum(dense; dims = 2))

        # The whole Q1 space, not merely affine fields: a stencil that is not a
        # tensor product loses the cross term.
        q1(x, y) = 2.0 + 0.01x + 0.03y + 0.0007x * y
        chart(space, i) = GR.chartcoords(qsrc, cellcentroid(space, i))
        field = [q1(chart(qsrc, i)...) for i in 1:Int(ncells(qsrc))]
        @test dense * field ≈ [q1(chart(qdst, i)...) for i in 1:Int(ncells(qdst))]

        # Once warm, a point on the chart path costs no allocation at all.
        row = GR.WeightRow()
        @test p1_samplerbytes(row, smp, cellcentroid(qdst, 1)) == 0 skip =
            VERSION < v"1.12"
    end

    @testset "a periodic chart axis wraps at the seam" begin
        # Seam fixture: 35 degrees east of the last centre on a
        # 90-degree lattice, so the two longitude weights are 11/18 and 7/18,
        # halved again across the two latitudes.
        src = ToyLonLatSpace(4, 2)
        src_inds = ownedindices(src, 1)
        seam = ToyLonLatSpace(1, 1; lon = (169.0, 171.0), lat = (-1.0, 1.0))
        entries = t4_entries(t4_build(BarycentricPoint(), seam, [1], src, src_inds),
            [1], src_inds)
        @test length(entries) == 4
        @test all(entries[(1, localindex(src, 4, iy))] ≈ 11 / 36 for iy in (1, 2))
        @test all(entries[(1, localindex(src, 1, iy))] ≈ 7 / 36 for iy in (1, 2))
        @test t4_rowsum(entries, 1) ≈ 1.0
    end

    @testset "a nonperiodic rim takes no weights" begin
        patch = ToyLonLatSpace(4, 2; lon = (-40.0, 40.0), lat = (-20.0, 20.0))
        patch_inds = ownedindices(patch, 1)
        smp = GR.sampler(BarycentricPoint(), patch)
        row = GR.WeightRow()
        barybuild(dst) = t4_entries(
            t4_build(BarycentricPoint(), dst, [1], patch, patch_inds), [1], patch_inds)

        # Past the last latitude centre the bilinear point drops to the rim line
        # and still writes a row; the barycentric point stops at the outermost
        # sample sites and leaves the destination to the missing policy.
        latrim = ToyLonLatSpace(1, 1; lon = (-6.0, -4.0), lat = (14.0, 16.0))
        clamped = t4_entries(t4_build(BilinearPoint(), latrim, [1], patch, patch_inds),
            [1], patch_inds)
        @test clamped[(1, localindex(patch, 2, 2))] ≈ 0.75
        @test clamped[(1, localindex(patch, 3, 2))] ≈ 0.25
        @test GR.weightsat!(row, smp, cellcentroid(latrim, 1)) === GR.WeightsRim
        @test isempty(row)
        @test isempty(barybuild(latrim))

        # And the same across a longitude rim.
        lonrim = ToyLonLatSpace(1, 1; lon = (34.0, 36.0), lat = (2.0, 4.0))
        clamped = t4_entries(t4_build(BilinearPoint(), lonrim, [1], patch, patch_inds),
            [1], patch_inds)
        @test clamped[(1, localindex(patch, 4, 1))] ≈ 0.35
        @test clamped[(1, localindex(patch, 4, 2))] ≈ 0.65
        @test GR.weightsat!(row, smp, cellcentroid(lonrim, 1)) === GR.WeightsRim
        @test isempty(barybuild(lonrim))

        # Beyond the rim the point is outside the source altogether, and the two
        # are told apart by the source boundary half a cell past the last site.
        rim = ToyLonLatSpace(1, 1; lon = (-1.0, 1.0), lat = (18.0, 20.0))
        beyond = ToyLonLatSpace(1, 1; lon = (-1.0, 1.0), lat = (20.0, 22.0))
        @test GR.weightsat!(row, smp, cellcentroid(rim, 1)) === GR.WeightsRim
        @test GR.weightsat!(row, smp, cellcentroid(beyond, 1)) === GR.WeightsOutside

        # Outside fixture, which the bilinear point serves from the rim
        # line, and a longitude past the span on whatever branch the chart
        # reports.
        outside = ToyLonLatSpace(1, 1; lon = (-1.0, 1.0), lat = (59.0, 61.0))
        @test GR.weightsat!(row, smp, cellcentroid(outside, 1)) === GR.WeightsOutside
        @test isempty(barybuild(outside))
        west = ToyLonLatSpace(1, 1; lon = (-171.0, -169.0), lat = (-1.0, 1.0))
        @test GR.chartcoords(patch, cellcentroid(west, 1))[1] ≈ 190.0
        @test GR.weightsat!(row, smp, cellcentroid(west, 1)) === GR.WeightsOutside
        @test isempty(row)
    end

    @testset "a destination at a source sample site takes that source alone" begin
        # A destination lattice on the source's own sites, on a nonperiodic
        # source and on a periodic one: one entry of weight exactly one, with no
        # second entry left by roundoff at the site.
        for src in (ToyLonLatSpace(6, 5; lon = (-60.0, 60.0), lat = (-50.0, 50.0)),
            ToyLonLatSpace(4, 2))
            n = Int(ncells(src))
            dense = Matrix(GR.wholeblock(BarycentricPoint(), src, src).weights)
            @test dense == Matrix(LinearAlgebra.I, n, n)
            @test count(!iszero, dense) == n
        end
    end

    @testset "descending and one-cell chart axes" begin
        # A raster whose latitude lookup descends. The lattice index a stencil
        # names must be the one that lookup order gives, so the operator onto
        # the same cells in ascending order matches them site for site.
        down = RasterGrid(DD.DimArray(zeros(8, 4),
            (DD.X(raster_lon()), DD.Y(reverse(raster_lat())))))
        up = RasterGrid(DD.DimArray(zeros(8, 4),
            (DD.X(raster_lon()), DD.Y(raster_lat()))))
        sitekey(space, i) = round.(Tuple(cellcentroid(space, i)), digits = 10)
        matching = Dict(sitekey(down, j) => j for j in 1:Int(ncells(down)))
        dense = Matrix(GR.wholeblock(BarycentricPoint(), up, down).weights)
        @test count(!iszero, dense) == Int(ncells(up))
        @test all(dense[i, matching[sitekey(up, i)]] == 1.0 for i in 1:Int(ncells(up)))

        # An interior destination on the descending source reproduces Q1.
        mid = RasterGrid(DD.DimArray(zeros(7, 3),
            (DD.X(-135.0:45.0:135.0), DD.Y(-45.0:45.0:45.0))))
        q1(x, y) = 2.0 + 0.01x + 0.03y + 0.0007x * y
        chart(space, i) = GR.chartcoords(down, cellcentroid(space, i))
        field = [q1(chart(down, j)...) for j in 1:Int(ncells(down))]
        W = Matrix(GR.wholeblock(BarycentricPoint(), mid, down).weights)
        @test all(count(!iszero, view(W, r, :)) == 4 for r in axes(W, 1))
        @test W * field ≈ [q1(chart(mid, i)...) for i in 1:Int(ncells(mid))]

        # A one-cell axis has no interval, so its single site takes the whole
        # axis weight as one entry rather than a repeated index.
        strip = RasterGrid(DD.DimArray(zeros(8, 1),
            (DD.X(raster_lon()), DD.Y(0.0:45.0:0.0))))
        across = RasterGrid(DD.DimArray(zeros(7, 1),
            (DD.X(-135.0:45.0:135.0), DD.Y(0.0:45.0:0.0))))
        S = Matrix(GR.wholeblock(BarycentricPoint(), across, strip).weights)
        @test all(count(!iszero, view(S, r, :)) == 2 for r in axes(S, 1))
        @test all(isapprox(1.0), sum(S; dims = 2))
        @test Matrix(GR.wholeblock(BarycentricPoint(), strip, strip).weights) ==
              Matrix(LinearAlgebra.I, Int(ncells(strip)), Int(ncells(strip)))
    end

    @testset "chart stencils survive chunking" begin
        src = ToyLonLatSpace(8, 6; lon = (-80.0, 80.0), lat = (-60.0, 60.0),
            chunks = (2, 2))
        dst = ToyLonLatSpace(5, 4; lon = (-50.0, 50.0), lat = (-40.0, 40.0),
            chunks = (5, 2))
        @test nchunks(src) == 12 && nchunks(dst) == 2

        # The support radius is the larger chart spacing, and it belongs to the
        # state a source prepares rather than to its type.
        @test supportradius(BarycentricPoint(), src) ≈
              deg2rad(max(dlon(src), dlat(src)))
        @test supportradius(BarycentricPoint(), src) ==
              supportradius(BilinearPoint(), src)
        @test supportradius(BarycentricPoint(), P1PlainSpace(src)) == 0.0

        # No source-chunk boundary changes a stencil: every pair's share of the
        # operator adds back up to the eager one.
        whole = Matrix(GR.wholeblock(BarycentricPoint(), dst, src).weights)
        plan = ChunkedPlan(BarycentricPoint(), Weighted(0.5), dst, src;
            storage = PerChunk())
        @test p1_assemble(plan, dst, src) ≈ whole
        @test all(isapprox(1.0), sum(whole; dims = 2))

        # Eager and chunked values agree, and a differently chunked source gives
        # the same values again.
        chart(i) = GR.chartcoords(src, cellcentroid(src, i))
        field = [p1_bilinear(chart(i)) for i in 1:Int(ncells(src))]
        eager = regrid(field; to = dst, from = src, method = BarycentricPoint(),
            lazy = false)
        @test collect(eager) ≈ whole * field
        chunked = Vector{Float64}(undef, Int(ncells(dst)))
        regrid!(chunked, field,
            ChunkedPlan(BarycentricPoint(), Weighted(0.5), dst,
                ToyLonLatSpace(8, 6; lon = (-80.0, 80.0), lat = (-60.0, 60.0),
                    chunks = (8, 2)); storage = PerChunk()))
        @test chunked == eager
        differently = Vector{Float64}(undef, Int(ncells(dst)))
        regrid!(differently, field,
            ChunkedPlan(BarycentricPoint(), Weighted(0.5), dst,
                ToyLonLatSpace(8, 6; lon = (-80.0, 80.0), lat = (-60.0, 60.0),
                    chunks = (8, 3)); storage = PerChunk()))
        @test differently == chunked

        # Discovery has no false negative. Each source cell here is its own
        # chunk, and the destination sits in the corner of one of them, so the
        # far corner of its stencil is a chunk its own cell corners are nowhere
        # near. Only the support radius reaches it.
        cells = ToyLonLatSpace(8, 6; lon = (-80.0, 80.0), lat = (-60.0, 60.0),
            chunks = (1, 1))
        off = ToyLonLatSpace(2, 2; lon = (-36.0, -32.0), lat = (-36.0, -32.0),
            chunks = (1, 1))
        offwhole = Matrix(GR.wholeblock(BarycentricPoint(), off, cells).weights)
        @test all(isapprox(1.0), sum(offwhole; dims = 2))
        far = localindex(cells, 2, 1)
        @test offwhole[1, far] != 0
        @test all(!GR.US._contains(GR.chunkextents(cells)[GR.chunkat(cells, far)], p)
                  for p in cellcorners(off, 1))
        offplan = ChunkedPlan(BarycentricPoint(), Weighted(0.5), off, cells;
            storage = PerChunk())
        @test GR.chunkat(cells, far) in GR.sourcesof(GR.dependencies(offplan), 1)
        @test p1_assemble(offplan, off, cells) ≈ offwhole
    end

end
