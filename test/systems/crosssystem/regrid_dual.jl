module RegridDualTests

using Test
import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GlobalRegridding as GR
import GeometryOps as GO
using Random

# Shared oracles keep dependency expectations identical across suites.
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

function interpolate(space, smp, row, p)
    GR.ismapped(GR.weightsat!(row, smp, p)) || return nothing
    return sum(row.weights[k] * field(GR.cellcentroid(space, row.indices[k]))
               for k in 1:length(row))
end

# Site-to-corner midpoints cover every dual cell of a conforming level.
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

# Coordinate matching supplies an independent topology oracle.
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

# Exact stencil pairs verify dependency-relation reach.
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

# Mixed-level representative selection, affine laws, and native interpolation.

function mixed(sys, top::Int)
    roots = collect(DGG.CellVector(DGG.levelgrid(sys, top + 1)))
    kids(c, l) = collect(DGG.CellVector(DGG.subtree(sys, c, l)))
    cells = vcat(kids(roots[1], top + 3), kids(roots[2], top + 2), roots[3:end])
    return DGG.MultiOrderVector(sys, cells; reference_level = top + 3)
end

# Find a shallow cell whose interpolation fan crosses the refinement boundary.
function coarsebesidefine(sys, mov, top::Int)
    lvl = DGG.levelgrid(sys, top + 1)
    fine = DGG.cellindex(lvl, 1)
    for i in 1:length(mov)
        c = mov[i]
        DGG.level(c) == top + 1 || continue
        any(n -> n == fine, DGG.neighbors(lvl, c, 1; connectivity = DGG.Vertex())) &&
            return c
    end
    return nothing
end

# Stay away from boundaries that make containment ambiguous.
function insidepoints(grid, c, t::Float64)
    centre = DGG.cell_centroid(grid, c)
    return [unit((centre[1] + t * (v[1] - centre[1]),
        centre[2] + t * (v[2] - centre[2]),
        centre[3] + t * (v[3] - centre[3])))
            for v in DGG.cell_boundary(grid, c)]
end

# Full-container scan provides an independent oracle for refined neighbors.
function representatives(mov, hostid, p)
    sys = DGG.system(mov)
    out = Set{Int}()
    for n in DGG.neighbors(DGG.levelgrid(sys, DGG.level(hostid)), hostid, 1;
        connectivity = DGG.Vertex())
        j = DGG.covering_index(mov, n)
        if j !== nothing
            push!(out, j)
            continue
        end
        best, bestd = 0, Inf
        for i in 1:length(mov)
            c = mov[i]
            DGG.level(c) > DGG.level(n) || continue
            DGG.ancestor(sys, c, DGG.level(n)) == n || continue
            q = DGG.cell_centroid(DGG.levelgrid(sys, DGG.level(c)), c)
            d = (q[1] - p[1])^2 + (q[2] - p[2])^2 + (q[3] - p[3])^2
            d < bestd && (bestd = d; best = i)
        end
        best == 0 || push!(out, best)
    end
    return out
end

# A typed wrapper keeps the allocation measurement representative.
function sweepstatus(row, smp, pts)
    n = 0
    for p in pts
        GR.ismapped(GR.weightsat!(row, smp, p)) && (n += 1)
    end
    return n
end

function interpolate_at(space, smp, row, vals, p)
    GR.ismapped(GR.weightsat!(row, smp, p)) || return nothing
    return sum(row.weights[k] * vals[row.indices[k]] for k in 1:length(row))
end

function maxadjacentjump(grid, v)
    ok(x) = !(ismissing(x) || (x isa AbstractFloat && isnan(x)))
    worst = 0.0
    for i in 1:DGG.ncells(grid)
        ok(v[i]) || continue
        c = DGG.cellindex(grid, i)
        for d in DGG.neighbors(grid, c, 1; connectivity = DGG.Edge())
            j = DGG.localindex(grid, d)
            (j === nothing || !ok(v[j])) && continue
            worst = max(worst, abs(v[i] - v[j]))
        end
    end
    return worst
end

# A smooth field isolates interpolation artifacts from source discontinuities.
smooth(q) = 0.3 + 1.7 * q[1] - 0.9 * q[2] + 2.1 * q[3] + 0.6 * q[1] * q[2]

storedvalues(mov, f) = [f(DGG.cell_centroid(
    DGG.levelgrid(DGG.system(mov), DGG.level(mov[i])), mov[i])) for i in 1:length(mov)]

@testset "dual cells over mixed levels on $(nameof(typeof(sys)))" for (sys, top) in
                                                                     ((HEALPIX, 0), (IGEO7, 0))

    mov = mixed(sys, top)
    ref = DGG.reference_level(mov)
    space = GR.sourcespacefor(mov, GR.BarycentricPoint())
    smp = GR.sampler(GR.BarycentricPoint(), space)
    row = GR.WeightRow()
    lref = DGG.levelgrid(sys, ref)
    refpoints = [DGG.cell_centroid(lref, DGG.cellindex(lref, i))
                 for i in 1:DGG.ncells(lref)]

    @testset "the space builds them, off the container" begin
        @test space.grid isa DGG.MultiOrderGrid
        @test GR.hasdualcells(space)
        st = GR.samplerstate(space)
        @test st isa DGG.MultiOrderDualTopology
        @test st.cells === mov
        @test smp.state === st
        @test occursin("MultiOrderDualTopology", sprint(show, st))
        @test GR.chartat(smp, refpoints[1]) === (0.0, 0.0)
        # Triangles remain well-defined across mixed-level T-junctions.
        @test all(p -> GR.nodecount(GR.dualcellat(smp, p)) == 3, refpoints)
        @test all(p -> GR.dualcellat(smp, p).kind === GR.MeanValue, refpoints)
    end

    @testset "a stored site reproduces itself" begin
        # Exact node hits expose misordered and misresolved fans.
        exact = 0
        for i in 1:length(mov)
            c = mov[i]
            p = DGG.cell_centroid(DGG.levelgrid(sys, DGG.level(c)), c)
            GR.ismapped(GR.weightsat!(row, smp, p)) || continue
            (length(row) == 1 && row.indices[1] == i && row.weights[1] == 1.0) &&
                (exact += 1)
        end
        @test exact == length(mov)
    end

    @testset "affine reproduction inside a coarse cell beside fine ones" begin
        coarse = coarsebesidefine(sys, mov, top)
        @test coarse !== nothing
        lvl = DGG.levelgrid(sys, top + 1)
        host = DGG.localindex(mov, coarse)

        # Three depths exercise different wedges and neighbor representatives.
        pts = vcat((insidepoints(lvl, coarse, t) for t in (0.2, 0.5, 0.85))...)
        @test length(pts) >= 8

        charterr = 0.0
        allhosted = true
        allresolved = true
        reached = Set{Int}()
        for p in pts
            @test GR.ismapped(GR.weightsat!(row, smp, p))
            allhosted &= DGG.localindex(mov, p) == host
            cell = GR.dualcellat(smp, p)
            # The query-centered chart makes the expected affine value the origin.
            x = y = 0.0
            for k in 1:length(row)
                j = findfirst(==(row.indices[k]), cell.indices)
                x += row.weights[k] * cell.nodes[j][1]
                y += row.weights[k] * cell.nodes[j][2]
            end
            charterr = max(charterr, abs(x), abs(y))
            # Affine laws alone cannot detect the wrong descendant representative.
            want = representatives(mov, coarse, p)
            for i in row.indices
                allresolved &= (i == host || i in want)
                push!(reached, i)
            end
        end
        @test allhosted
        @test charterr < 1e-10
        @test allresolved
        @test any(i -> DGG.level(mov[i]) > DGG.level(coarse), reached)
    end

    @testset "the point laws hold across the level boundary" begin
        mapped = 0
        positive = true
        finite = true
        sumerr = 0.0
        charterr = 0.0
        for p in refpoints
            status = GR.weightsat!(row, smp, p)
            if GR.ismapped(status)
                mapped += 1
                positive &= all(>=(0.0), row.weights)
                finite &= all(isfinite, row.weights)
                sumerr = max(sumerr, abs(sum(row.weights) - 1.0))
                cell = GR.dualcellat(smp, p)
                x = y = 0.0
                for k in 1:length(row)
                    j = findfirst(==(row.indices[k]), cell.indices)
                    x += row.weights[k] * cell.nodes[j][1]
                    y += row.weights[k] * cell.nodes[j][2]
                end
                charterr = max(charterr, abs(x), abs(y))
            else
                # Degeneracies must clear the row and return a status.
                @test isempty(row)
            end
        end
        # Full coverage exposes flipped sliver triangles as unmapped sites.
        @test mapped == length(refpoints)
        @test positive
        @test finite
        @test sumerr < 1e-14
        @test charterr < 1e-14

        # Fixed-capacity fans allocate nothing after warm-up.
        sweepstatus(row, smp, refpoints)
        @test (@allocated sweepstatus(row, smp, refpoints)) == 0 skip = VERSION < v"1.12"
    end

    @testset "a hole is a rim, not a closed fan" begin
        # Closing a fan across a hole would assign weight to uncovered ground.
        gapat = findfirst(i -> DGG.level(mov[i]) == top + 1, 1:length(mov))
        gap = mov[gapat]
        holed = DGG.MultiOrderVector(sys,
            [mov[i] for i in 1:length(mov) if i != gapat]; reference_level = ref)
        hspace = GR.sourcespacefor(holed, GR.BarycentricPoint())
        hsmp = GR.sampler(GR.BarycentricPoint(), hspace)
        lvl = DGG.levelgrid(sys, top + 1)

        inside = DGG.cell_centroid(lvl, gap)
        @test (@inferred GR.weightsat!(row, hsmp, inside)) === GR.WeightsOutside
        @test isempty(row)

        rims = 0
        mapped = 0
        for n in DGG.neighbors(lvl, gap, 1; connectivity = DGG.Vertex())
            DGG.covering_index(holed, n) === nothing && continue
            for p in insidepoints(lvl, n, 0.75)
                status = GR.weightsat!(row, hsmp, p)
                status === GR.WeightsRim && (rims += 1; @test isempty(row))
                GR.ismapped(status) && (mapped += 1)
            end
        end
        @test rims > 0
        @test mapped > 0
    end

    @testset "native interpolation reduces expansion artifacts" begin
        # Replicated leaf values create flat interiors and steps at leaf spacing.
        vals = storedvalues(mov, smooth)
        cube = DD.DimArray(vals, DGG.Cells(DGG.MultiOrderLookup(mov)))
        expanded = DGG.expand(cube, ref)
        aspace = GR.sourcespacefor(mov, GR.Conservative())
        asmp = GR.sampler(GR.BarycentricPoint(), aspace)
        avals = collect(parent(expanded))
        @test DGG.ncells(aspace) == length(avals) > length(mov)

        coarse = coarsebesidefine(sys, mov, top)
        lvl = DGG.levelgrid(sys, top + 1)
        pts = insidepoints(lvl, coarse, 0.35)
        a = [interpolate_at(aspace, asmp, row, avals, p) for p in pts]
        b = [interpolate_at(space, smp, row, vals, p) for p in pts]
        @test all(!isnothing, a) && all(!isnothing, b)
        af, bf = Float64[x for x in a], Float64[x for x in b]
        # Expansion is flat inside the coarse cell; native interpolation varies.
        @test maximum(af) - minimum(af) < 1e-12
        @test maximum(bf) - minimum(bf) > 1e-3

        # Native interpolation reduces the largest adjacent-cell step.
        dst = lref
        nat = parent(DGG.regrid(cube; to = dst, method = GR.BarycentricPoint()))
        exp_ = parent(DGG.regrid(expanded; to = dst, method = GR.BarycentricPoint()))
        jn, je = maxadjacentjump(dst, nat), maxadjacentjump(dst, exp_)
        @test jn < je
        @test jn < 0.5 * je
        truth = [smooth(DGG.cell_centroid(dst, DGG.cellindex(dst, i)))
                 for i in 1:DGG.ncells(dst)]
        @test maximum(abs.(nat .- truth)) < maximum(abs.(exp_ .- truth))
    end

    @testset "a uniform container keeps the conforming construction" begin
        # One stored cell per leaf retains the conforming grid's native dual cells.
        cells = collect(DGG.CellVector(DGG.levelgrid(sys, top + 2)))
        uniform = DGG.MultiOrderVector(sys, cells; reference_level = top + 2)
        uspace = GR.sourcespacefor(uniform, GR.BarycentricPoint())
        @test uspace.grid isa DGG.PartialGrid
        @test GR.samplerstate(uspace) isa DGG.DualTopology
        usmp = GR.sampler(GR.BarycentricPoint(), uspace)
        lvl = DGG.levelgrid(sys, top + 2)
        probes = [midpoint(DGG.cell_centroid(lvl, DGG.cellindex(lvl, i)),
            first(DGG.cell_boundary(lvl, DGG.cellindex(lvl, i))))
                  for i in 1:min(200, DGG.ncells(lvl))]
        counts = Set(GR.nodecount(GR.dualcellat(usmp, p)) for p in probes)
        @test sys === HEALPIX ? (4 in counts) : counts == Set((3,))

        # The uniform container must match the equivalent level grid exactly.
        dst = DGG.levelgrid(sys, top + 1)
        vals = storedvalues(uniform, smooth)
        cube = DD.DimArray(vals, DGG.Cells(DGG.MultiOrderLookup(uniform)))
        @test isequal(parent(DGG.regrid(cube; to = dst,
                method = GR.BarycentricPoint())),
            collect(GR.regrid(vals; to = DGG.DGGSpace(dst),
                from = DGG.DGGSpace(lvl), method = GR.BarycentricPoint(),
                lazy = false)))
    end

    @testset "DGGS point sampling ignores the raster pole policy" begin
        # DGGS cells reach the poles, so raster pole synthesis is irrelevant.
        pole = unit((0.0, 0.0, 1.0))
        polar = GR.sampler(GR.BarycentricPoint(poles = nothing), space)
        a, b = GR.WeightRow(), GR.WeightRow()
        same = true
        for p in vcat([pole, unit((0.0, 0.0, -1.0))],
            refpoints[1:min(50, length(refpoints))])
            GR.weightsat!(a, smp, p)
            GR.weightsat!(b, polar, p)
            same &= a.indices == b.indices && a.weights == b.weights
        end
        @test same
        @test GR.ismapped(GR.weightsat!(row, smp, pole))

        @test (@inferred GR.dualcellat(smp, pole)) isa DGG.GridDualCell
        @test (@inferred GR.weightsat!(row, smp, pole)) === GR.WeightsMapped
        @test (@inferred GR.chartat(smp, pole)) === (0.0, 0.0)
        # The stencil bound is twice the widest stored cell.
        r = @inferred GR.supportradius(GR.BarycentricPoint(), space)
        @test 0 < r < Float64(pi)
        widest = maximum(Float64(DGG.Fallbacks.cell_cap(
            DGG.levelgrid(sys, DGG.level(c)), c).radius) for c in mov)
        @test r ≈ 2 * widest
    end
end

end # module
