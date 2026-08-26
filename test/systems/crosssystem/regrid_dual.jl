# Dual cells behind a `DGGSpace`: what a point method interpolates on when the
# source is a conforming DGGS. The polygon is the dual cell of a primal vertex —
# the sample sites of the cells meeting there — and it is found from the host
# cell's one-rings, so the tests below ask three things: that the cells it finds
# are the ones the geometry has, that the interpolant on them obeys the point
# laws, and that no chunking changes either.

module RegridDualTests

using Test
import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import GeometryOps as GO
using Random

# The relation oracles are defined once, in the GlobalRegridding suite, and
# shared with it rather than re-spelled here.
include(joinpath(@__DIR__, "..", "..", "..", "lib", "GlobalRegridding", "test",
    "graphoracles.jl"))
using .ChunkGraphOracles: graph_pairs

const IGEO7 = DGG.IGeo7System()
const S2 = DGG.S2System()
const HEALPIX = DGG.HEALPixSystem()
const ISEA4R = DGG.ISEA4RSystem()

unit(v) = (n = sqrt(sum(abs2, v)); GO.UnitSphericalPoint(v[1] / n, v[2] / n, v[3] / n))
midpoint(a, b) = unit((a[1] + b[1], a[2] + b[2], a[3] + b[3]))

# A smooth field to interpolate, read at a source cell's own sample site.
field(q) = 0.3 + 1.7 * q[1] - 0.9 * q[2] + 2.1 * q[3]
sourcefield(space) = [field(GR.cellcentroid(space, i)) for i in 1:Int(GR.ncells(space))]

# What the stencil at `p` makes of that field, or `nothing` where `p` is unmapped.
function interpolate(space, smp, row, p)
    GR.ismapped(GR.weightsat!(row, smp, p)) || return nothing
    return sum(row.weights[k] * field(GR.cellcentroid(space, row.indices[k]))
               for k in 1:length(row))
end

# Every distinct dual cell of a whole level, by probing each cell's own boundary:
# the midpoint between a cell's site and one of its corners lies in that
# corner's dual cell, so probing every boundary point of every cell reaches
# every dual cell there is.
function alldualcells(sys, lvl)
    grid = DGG.levelgrid(sys, lvl)
    space = DGG.DGGSpace(grid)
    smp = GR.sampler(GR.BarycentricPoint(), space)
    seen = Set{Set{Int}}()
    unmapped = 0
    kinds = Set{GR.BasisKind}()
    for i in 1:Int(GR.ncells(space))
        c = DGG.cellindex(grid, i)
        centre = DGG.cell_centroid(grid, c)
        for v in DGG.cell_boundary(grid, c)
            cell = GR.dualcellat(smp, midpoint(centre, v))
            if GR.nodecount(cell) == 0
                unmapped += 1
            else
                push!(seen, Set(collect(cell.indices)))
                push!(kinds, cell.kind)
            end
        end
    end
    tally = Dict{Int,Int}()
    for s in seen
        tally[length(s)] = get(tally, length(s), 0) + 1
    end
    return tally, unmapped, kinds
end

# The oracle, allowed here and nowhere near production: the cells of `grid`
# whose drawn boundary passes through the point `v`, found by matching
# coordinates among the host and its vertex ring.
function incidentcells(grid, host, v, tol)
    out = Set{Int}()
    for c in Iterators.flatten(((host,),
        DGG.neighbors(grid, host, 1; connectivity = DGG.Vertex())))
        for w in DGG.cell_boundary(grid, c)
            if (w[1] - v[1])^2 + (w[2] - v[2])^2 + (w[3] - v[3])^2 <= tol^2
                push!(out, DGG.localindex(grid, c))
                break
            end
        end
    end
    return out
end

# One operator assembled from every chunk pair, in the spaces' local indices.
function assemble(plan, dst, src)
    M = zeros(Float64, Int(GR.ncells(dst)), Int(GR.ncells(src)))
    for d in 1:GR.nchunks(dst), s in 1:GR.nchunks(src)
        block = Matrix(GR.buildblock(plan, d, s).weights)
        dinds, sinds = GR.ownedindices(dst, d), GR.ownedindices(src, s)
        for (r, i) in enumerate(dinds), (c, j) in enumerate(sinds)
            M[Int(i), Int(j)] += block[r, c]
        end
    end
    return M
end

# The chunk pairs the stencils themselves name: the exact relation a point
# method has to be a superset of.
function stencilpairs(dst, src)
    smp = GR.sampler(GR.BarycentricPoint(), src)
    row = GR.WeightRow()
    pairs = Set{Tuple{Int,Int}}()
    for d in 1:GR.nchunks(dst), i in GR.ownedindices(dst, d)
        GR.ismapped(GR.weightsat!(row, smp, GR.cellcentroid(dst, Int(i)))) || continue
        for k in 1:length(row)
            push!(pairs, (d, GR.chunkat(src, row.indices[k])))
        end
    end
    return pairs
end

querybytes(row, smp, p) =
    (GR.weightsat!(row, smp, p); @allocated GR.weightsat!(row, smp, p))

@testset "dual cells on a conforming DGGS" begin

    @testset "the space says it has them, and prepares them once" begin
        space = DGG.DGGSpace(DGG.levelgrid(IGEO7, 2))
        @test GR.hasdualcells(space)
        state = @inferred GR.samplerstate(space)
        @test state isa DGG.DualTopology
        smp = GR.sampler(GR.BarycentricPoint(), space)
        # A complete grid is its own topology; a subset reads its level's.
        @test smp.state.grid === space.grid

        # The chart is centred on the queried point, so the point is the origin
        # and the sampler's chart says exactly that.
        p = GR.cellcentroid(space, 7)
        @test GR.chartat(smp, p) === (0.0, 0.0)

        # Both entry points infer concretely, so nothing on the hot path boxes.
        @test (@inferred GR.dualcellat(smp, p)) isa DGG.GridDualCell
        row = GR.WeightRow()
        @test (@inferred GR.weightsat!(row, smp, p)) === GR.WeightsMapped
    end

    @testset "the dual complex is the one the geometry has" begin
        # Every dual cell of a whole level, counted by node count. A grid of `F`
        # cells and `E` edges has `E - F + 2` primal vertices and therefore that
        # many dual cells, so these tallies are the sphere's own arithmetic: a
        # construction that merged two fans, split one, or missed a defect
        # cannot reproduce them.
        #
        #   IGeo7  480 hexagons and 12 pentagons: every vertex carries three
        #          cells, so every dual cell is a triangle.
        #   S2     six square faces: four cells to a vertex, except the eight
        #          cube corners where three meet.
        #   HEALPix four cells to a vertex, except its eight special points.
        #   ISEA4R  ten diamonds: four cells to a vertex, five at the two poles
        #          of the diamond layout and three at the ten other icosahedral
        #          vertices.
        for (name, sys, lvl, want) in (
            ("IGeo7", IGEO7, 2, Dict(3 => 980)),
            ("S2", S2, 2, Dict(4 => 90, 3 => 8)),
            ("HEALPix", HEALPIX, 2, Dict(4 => 186, 3 => 8)),
            ("ISEA4R", ISEA4R, 2, Dict(5 => 2, 4 => 150, 3 => 10)))
            @testset "$name" begin
                tally, unmapped, kinds = alldualcells(sys, lvl)
                @test tally == want
                @test unmapped == 0
                # A cell of any node count says which coordinates weight it, and
                # says the same one: nothing is silently triangulated, and a
                # four-node cell is not read as a tensor-product one.
                @test kinds == Set((GR.MeanValue,))
            end
        end
    end

    @testset "the ring walk finds the cells the coordinates do" begin
        # The oracle: match the drawn cell boundaries at a genuine primal vertex
        # — one carrying three or more cells — and compare with the cell the
        # ring walk answers. Production never does this; a test may.
        for (name, sys, lvl) in (("IGeo7", IGEO7, 3), ("S2", S2, 3),
            ("HEALPix", HEALPIX, 3), ("ISEA4R", ISEA4R, 3))
            @testset "$name" begin
                grid = DGG.levelgrid(sys, lvl)
                space = DGG.DGGSpace(grid)
                smp = GR.sampler(GR.BarycentricPoint(), space)
                rng = MersenneTwister(11)
                picks = unique(rand(rng, 1:Int(GR.ncells(space)), 30))
                # Every defect the level has, alongside the sample: a pentagon
                # or a low-valence corner is where a fan-finding mistake hides.
                defects = [i for i in 1:Int(GR.ncells(space))
                           if DGG.neighborcount(grid, DGG.cellindex(grid, i)) <
                              DGG.maxneighbors(sys, DGG.Vertex())]
                checked = 0
                agreed = true
                for i in unique(vcat(picks, first(defects, 12)))
                    host = DGG.cellindex(grid, i)
                    centre = DGG.cell_centroid(grid, host)
                    for v in DGG.cell_boundary(grid, host)
                        want = incidentcells(grid, host, v, 1e-7)
                        length(want) >= 3 || continue
                        checked += 1
                        got = Set(collect(GR.dualcellat(smp,
                            midpoint(centre, v)).indices))
                        agreed &= got == want
                    end
                end
                @test checked > 0
                @test agreed
            end
        end
    end

    @testset "the point laws hold on every dual cell" begin
        for (name, sys, lvl) in (("IGeo7", IGEO7, 3), ("S2", S2, 3),
            ("HEALPix", HEALPIX, 3), ("ISEA4R", ISEA4R, 3))
            @testset "$name" begin
                grid = DGG.levelgrid(sys, lvl)
                space = DGG.DGGSpace(grid)
                smp = GR.sampler(GR.BarycentricPoint(), space)
                row = GR.WeightRow()
                rng = MersenneTwister(4)
                points = [unit((randn(rng), randn(rng), randn(rng))) for _ in 1:2000]

                mapped = true
                positive = true
                sumerr = 0.0
                charterr = 0.0
                for p in points
                    mapped &= GR.ismapped(GR.weightsat!(row, smp, p))
                    positive &= all(>=(0.0), row.weights)
                    sumerr = max(sumerr, abs(sum(row.weights) - 1.0))
                    # Affine reproduction in the chart. The chart is centred on
                    # the query, so the affine field it must reproduce is the
                    # coordinate field at the origin: the weighted node
                    # coordinates are the point itself, which is `(0, 0)`.
                    cell = GR.dualcellat(smp, p)
                    x = y = 0.0
                    for k in 1:length(row)
                        j = findfirst(==(row.indices[k]), cell.indices)
                        x += row.weights[k] * cell.nodes[j][1]
                        y += row.weights[k] * cell.nodes[j][2]
                    end
                    charterr = max(charterr, abs(x), abs(y))
                end
                @test mapped
                @test positive
                @test sumerr < 1e-14
                # The cells are a hundredth of a radian across at these levels,
                # so this is roundoff on the coordinates themselves.
                @test charterr < 1e-15

                # A destination exactly at a source sample site takes that
                # source alone: the site is the chart's own origin, so the row
                # is one entry of weight one.
                sites = true
                for i in 1:Int(GR.ncells(space))
                    GR.weightsat!(row, smp, GR.cellcentroid(space, i))
                    sites &= length(row) == 1 && row.indices[1] == i &&
                             row.weights[1] == 1.0
                end
                @test sites

                # Nothing per destination, once warm.
                @test querybytes(row, smp, points[1]) == 0 skip = VERSION < v"1.12"
            end
        end
    end

    @testset "the interpolant is continuous" begin
        # Step across a primal edge and the interpolant may change dual cell,
        # but not value: the two cells share the dual edge between the two
        # sample sites, and both read it the same way. A construction that
        # depended on which cell hosted the point would jump by a cell's worth
        # of the field here.
        for (name, sys, lvl) in (("IGeo7", IGEO7, 3), ("S2", S2, 3),
            ("ISEA4R", ISEA4R, 3), ("HEALPix", HEALPIX, 3))
            @testset "$name" begin
                grid = DGG.levelgrid(sys, lvl)
                space = DGG.DGGSpace(grid)
                smp = GR.sampler(GR.BarycentricPoint(), space)
                row = GR.WeightRow()
                rng = MersenneTwister(9)
                step = 1e-8
                worst = 0.0
                seams = 0
                worstseam = 0.0
                for _ in 1:300
                    i = rand(rng, 1:Int(GR.ncells(space)))
                    c = DGG.cellindex(grid, i)
                    nbrs = collect(DGG.neighbors(grid, c, 1; connectivity = DGG.Edge()))
                    isempty(nbrs) && continue
                    d = rand(rng, nbrs)
                    a1, a2 = DGG.cell_centroid(grid, c), DGG.cell_centroid(grid, d)
                    m = midpoint(a1, a2)
                    dir = unit((a2[1] - a1[1], a2[2] - a1[2], a2[3] - a1[3]))
                    lo = interpolate(space, smp, row,
                        unit((m[1] - step * dir[1], m[2] - step * dir[2],
                            m[3] - step * dir[3])))
                    hi = interpolate(space, smp, row,
                        unit((m[1] + step * dir[1], m[2] + step * dir[2],
                            m[3] + step * dir[3])))
                    (lo === nothing || hi === nothing) && continue
                    jump = abs(lo - hi)
                    worst = max(worst, jump)
                    # A face seam: the two cells hang under different roots of
                    # the system, which is where a chart or an ordering built
                    # per face would come apart.
                    if DGG.ancestor(sys, c, 0) != DGG.ancestor(sys, d, 0)
                        seams += 1
                        worstseam = max(worstseam, jump)
                    end
                end
                @test seams > 0
                # The field varies by about `1e-2` across one of these cells, so
                # a discontinuity would be orders of magnitude above this.
                @test worst < 1e-6
                @test worstseam < 1e-6
            end
        end
    end

    @testset "no source chunk boundary changes a stencil" begin
        src = DGG.levelgrid(IGEO7, 3)
        dst = DGG.DGGSpace(DGG.levelgrid(S2, 3); chunklevel = 2)
        fine = DGG.DGGSpace(src; chunklevel = 2)
        coarse = DGG.DGGSpace(src; chunklevel = 1)
        @test GR.nchunks(fine) != GR.nchunks(coarse)

        # The same query on two different chunkings of one grid is the same row,
        # entry for entry: chunking is not an input to a stencil.
        finesmp = GR.sampler(GR.BarycentricPoint(), fine)
        coarsesmp = GR.sampler(GR.BarycentricPoint(), coarse)
        a, b = GR.WeightRow(), GR.WeightRow()
        same = true
        for i in 1:Int(GR.ncells(dst))
            p = GR.cellcentroid(dst, i)
            GR.weightsat!(a, finesmp, p)
            GR.weightsat!(b, coarsesmp, p)
            same &= a.indices == b.indices && a.weights == b.weights
        end
        @test same

        # Weights are geometry: the whole operator is built from the two spaces
        # with no field anywhere near it, and the chunk pairs add back up to it.
        whole = Matrix(GR.wholeblock(GR.BarycentricPoint(), dst, fine).weights)
        @test all(isapprox(1.0), sum(whole; dims = 2))
        @test assemble(GR.ChunkedPlan(GR.BarycentricPoint(), GR.Weighted(0.5), dst,
            fine; storage = GR.PerChunk()), dst, fine) == whole
        @test assemble(GR.ChunkedPlan(GR.BarycentricPoint(), GR.Weighted(0.5), dst,
            coarse; storage = GR.PerChunk()), dst, coarse) == whole

        # Eager, chunked and lazy read the same values off it.
        data = sourcefield(fine)
        eager = collect(GR.regrid(data; to = dst, from = fine,
            method = GR.BarycentricPoint(), lazy = false))
        @test eager ≈ whole * data
        chunked = Vector{Float64}(undef, Int(GR.ncells(dst)))
        GR.regrid!(chunked, data,
            GR.ChunkedPlan(GR.BarycentricPoint(), GR.Weighted(0.5), dst, coarse;
                storage = GR.PerChunk()))
        @test chunked ≈ eager atol = 1e-13
        @test collect(GR.regrid(data; to = dst, from = fine,
            method = GR.BarycentricPoint(), lazy = true)) ≈ eager atol = 1e-13

        # Discovery has no false negative: every chunk pair a stencil names is
        # in the plan's relation, which is what the declared support radius is
        # for and what the lazy read above would have refused without.
        @test GR.supportradius(GR.BarycentricPoint(), fine) > 0
        @test stencilpairs(dst, fine) ⊆
              graph_pairs(GR.dependencies(GR.ChunkedPlan(GR.BarycentricPoint(),
            GR.Weighted(0.5), dst, fine; storage = GR.PerChunk())))
    end

    @testset "a node the collection does not hold is a rim" begin
        # A rooted subtree of one level. Its border cells' dual cells need cells
        # outside it, which have no local index here at all — so the query is a
        # rim rather than a stencil closed over whatever the subtree does hold.
        sub = DGG.subtree(IGEO7, DGG.cellindex(DGG.levelgrid(IGEO7, 1), 3), 3)
        space = DGG.DGGSpace(sub)
        smp = GR.sampler(GR.BarycentricPoint(), space)
        # Topology is read on the complete level, so a ring is never one cell
        # short where the subtree ends.
        @test smp.state.grid === sub.complete
        row = GR.WeightRow()
        border = Set(DGG.border(sub))
        rims = Int[]
        for i in 1:Int(GR.ncells(space))
            status = GR.weightsat!(row, smp, GR.cellcentroid(space, i))
            status === GR.WeightsRim && push!(rims, i)
            GR.ismapped(status) || @test isempty(row)
        end
        @test !isempty(rims)
        # A rim is always a border cell: an interior cell's whole vertex ring is
        # inside the subtree, so every dual cell around it can be written here.
        @test Set(rims) ⊆ border

        # Every mapped stencil names only cells of this collection, at its own
        # local indices, and never borrows one from outside.
        inside = true
        for i in 1:Int(GR.ncells(space))
            GR.ismapped(GR.weightsat!(row, smp, GR.cellcentroid(space, i))) || continue
            inside &= all(1 <= j <= Int(GR.ncells(space)) for j in row.indices)
        end
        @test inside

        # A point the level covers and this collection does not is a rim too:
        # the cell holding it has no local index, so no dual cell around it can
        # be written here.
        whole = DGG.levelgrid(IGEO7, 3)
        away = DGG.cell_centroid(whole, DGG.cellindex(whole, 1))
        @test GR.weightsat!(row, smp, away) === GR.WeightsRim
        @test isempty(row)
    end

end

end # module
