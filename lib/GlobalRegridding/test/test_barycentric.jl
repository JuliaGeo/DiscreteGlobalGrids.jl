import DimensionalData as DD

# The point contract: one reusable weight row, the dual cell a point falls in,
# and the three coordinate kernels that weight its nodes. The kernels are tested
# on bare node coordinates, with no space and no chunk anywhere near them; the
# sampler and the chunk-pair builder are tested on a toy space that answers one
# hand-built dual cell.

# Source indices that are not node positions, so a kernel that emits `k` instead
# of `indices[k]` is caught rather than silently right.
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
name one of `indices`, so a row that emits a node position, or a source cell the
cell never held, fails here rather than in a value comparison.
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

    @testset "triangle barycentric coordinates" begin
        row = GR.WeightRow()
        inds, nodes = P1_TRI_INDS, P1_TRI_NODES

        # An interior point takes all three nodes, sums to one, and reproduces
        # the constant, affine and coordinate fields the basis promises.
        @test GR.linearweights!(row, inds, nodes, (1.5, 1.5)) === GR.WeightsMapped
        @test length(row) == 3
        @test sort(row.indices) == sort(inds)
        @test p1_rowsum(row) ≈ 1.0
        @test p1_reproduce(row, inds, nodes, p1_constant) ≈ 3.5
        @test p1_reproduce(row, inds, nodes, p1_affine) ≈ p1_affine((1.5, 1.5))
        @test all(isapprox.(p1_site(row, inds, nodes), (1.5, 1.5); atol = 1e-12))
        # The coordinates themselves, not merely a partition of unity.
        @test p1_weightof(row, 11) ≈ 0.5

        # Every affine field, at several interior points.
        for p in ((1.0, 1.0), (2.0, 0.5), (0.5, 0.5), (1.2, 2.0))
            @test GR.linearweights!(row, inds, nodes, p) === GR.WeightsMapped
            @test p1_reproduce(row, inds, nodes, p1_affine) ≈ p1_affine(p)
            @test p1_rowsum(row) ≈ 1.0
        end

        # A point on an edge emits that edge's two nodes and drops the zero.
        @test GR.linearweights!(row, inds, nodes, (2.0, 0.0)) === GR.WeightsMapped
        @test length(row) == 2
        @test sort(row.indices) == [3, 7]
        @test p1_reproduce(row, inds, nodes, p1_affine) ≈ p1_affine((2.0, 0.0))

        # A point at a node reproduces that node exactly, with one entry.
        for k in 1:3
            @test GR.linearweights!(row, inds, nodes, nodes[k]) === GR.WeightsMapped
            @test row.indices == [inds[k]]
            @test row.weights == [1.0]
        end

        # Outside is unmapped, never clamped to the boundary, and leaves no row.
        for p in ((5.0, 0.0), (-1.0, -1.0), (2.0, 3.0), (1.0, -0.5))
            @test GR.linearweights!(row, inds, nodes, p) === GR.WeightsOutside
            @test isempty(row)
        end

        # Degeneracies are rejected rather than answered with a division by zero.
        @test GR.linearweights!(row, inds, [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)],
            (1.0, 1.0)) === GR.WeightsDegenerate
        @test GR.linearweights!(row, inds, [(0.0, 0.0), (4.0, 0.0), (0.0, 0.0)],
            (1.0, 0.0)) === GR.WeightsDegenerate
        @test GR.linearweights!(row, inds, [(0.0, 0.0), (4.0, 0.0)],
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

        # On a triangle, mean-value coordinates are the barycentric ones.
        tri = GR.WeightRow()
        @test GR.meanvalueweights!(row, P1_TRI_INDS, P1_TRI_NODES,
            (1.5, 1.5)) === GR.WeightsMapped
        @test GR.linearweights!(tri, P1_TRI_INDS, P1_TRI_NODES,
            (1.5, 1.5)) === GR.WeightsMapped
        @test row.indices == tri.indices
        @test all(isapprox.(row.weights, tri.weights; atol = 1e-12))

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

    @testset "the row is reused and costs nothing warm" begin
        row = GR.WeightRow()

        # The row is cleared on entry, so a mapped point never inherits the
        # previous point's entries and an unmapped one leaves none.
        GR.linearweights!(row, P1_TRI_INDS, P1_TRI_NODES, (1.0, 1.0))
        @test length(row) == 3
        GR.linearweights!(row, P1_TRI_INDS, P1_TRI_NODES, (2.0, 0.0))
        @test length(row) == 2
        GR.linearweights!(row, P1_TRI_INDS, P1_TRI_NODES, (9.0, 9.0))
        @test isempty(row)

        # Once warm, a point costs no allocation at all.
        @test p1_kernelbytes(GR.linearweights!, row, P1_TRI_INDS, P1_TRI_NODES,
            (1.0, 1.0)) == 0 skip = VERSION < v"1.12"
        @test p1_kernelbytes(GR.bilinearweights!, row, P1_QUAD_INDS, P1_QUAD_NODES,
            (2.0, 1.0)) == 0 skip = VERSION < v"1.12"
        @test p1_kernelbytes(GR.meanvalueweights!, row, P1_HEX_INDS, P1_HEX_NODES,
            (0.3, 0.2)) == 0 skip = VERSION < v"1.12"
    end

    @testset "the sampler prepares a source space once" begin
        space = ToyLonLatSpace(4, 2; lon = (-40.0, 40.0), lat = (-20.0, 20.0))

        # Sample sites are the space's centroids, read on demand.
        sites = GR.samplesites(space)
        @test length(sites) == Int(ncells(space))
        @test all(sites[i] == cellcentroid(space, i) for i in 1:Int(ncells(space)))

        # A sampler holds the space and its sites, and is built once.
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
            :DualCell, :linearweights!, :bilinearweights!, :meanvalueweights!))
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
end
