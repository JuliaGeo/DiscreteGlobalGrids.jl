# The regridding face: a DGG grid read by `GlobalRegridding` through
# GridInterface's verbs, the target spellings that resolve into one, and the
# cube a result comes back as. One coarse S2 destination stands in for every
# system — cell geometry is each system's own business and is tested elsewhere
# — so what is under test here is the grid contract, the `to` resolution, the
# `Cells` axis, and that mass survives the DGG side of a conservative regrid.

module RegridTests

using Test
import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import DimensionalData as DD
import Extents
import GeometryOps as GO
using GeometryOps: SpatialTreeInterface as STI
import ConservativeRegridding as CR
import GridInterface as GrI
import GeometryOpsCore as GOCore
using SparseArrays: spzeros

const SYS = DGG.S2System()
const LEVEL = 3
const GRID = DGG.levelgrid(SYS, LEVEL)
const GLOBE = Extents.Extent(X = (-180.0, 180.0), Y = (-90.0, 90.0))

# The four dependency-graph relations — truth, demand, cap join, and the graph's
# own rows — are defined once, in the GlobalRegridding suite's `graphoracles.jl`,
# and shared with that suite and the benchmark harness. Do not re-spell any of
# them here.
include(joinpath(@__DIR__, "..", "..", "..", "lib", "GlobalRegridding", "test",
    "graphoracles.jl"))
using .ChunkGraphOracles: contributing_pairs, graph_pairs, demanded_pairs,
    capjoin_pairs

# A 15° global raster whose cells are declared as abutting intervals, so its
# edges tile the sphere exactly and a conservative regrid off it can conserve.
_axis(D, centres, step) = D(DD.Sampled(collect(centres); span = DD.Regular(step),
    sampling = DD.Intervals(DD.Center()), order = DD.ForwardOrdered()))

field(x, y, t) = 10 + sind(x) + cosd(2y) + t

function globalraster(step)
    lon = (-180 + step / 2):step:180
    lat = (-90 + step / 2):step:90
    return DD.DimArray([field(x, y, t) for x in lon, y in lat, t in 1:2],
        (_axis(DD.X, lon, step), _axis(DD.Y, lat, step), DD.Dim{:month}(1:2)))
end

const RASTER = globalraster(15.0)
const SRC = GR.RasterGrid(RASTER)

# A grid chunked by the ancestor level `a`, or by the ancestor level holding
# about `n` cells per chunk (the sizing the default `chunks(grid)` uses).
bylevel(g, a) = GrI.mockchunks(g, GrI.LevelChunks(a))
bycells(g, n) = bylevel(g,
    DGG._chunklevel(DGG.system(g), DGG.level(g), DGG.ncells(g), n))
# The DGG grid under a chunking wrapper.
gridof(g) = g isa GrI.ChunkedGrid ? parent(g) : g

# The cells the grid assigns to each chunk, concatenated in chunk order.
chunkcells(space) = reduce(vcat,
    [collect(GR._ownedindices(space, c)) for c in 1:GR._nchunks(space)])

capscover(space) = all(DGG.Fallbacks.cap_contains(GR._chunkextents(space)[c], p)
    for c in 1:GR._nchunks(space)
    for i in GR._ownedindices(space, c)
    for p in DGG.cell_boundary(gridof(space), DGG.cellindex(gridof(space), i)))

const REGION = DGG.covering(DGG.CellVector(GRID),
    Extents.Extent(X = (-140.0, 40.0), Y = (-20.0, 60.0)))

@testset "the grid contract" begin
    # A complete level, a scattered subset of one, and a hexagonal system whose
    # aperture-7 children reach outside their parent's own boundary — the case a
    # chunk cap taken from the ancestor cell rather than from its subtree extent
    # would fail to cover.
    for space in (bycells(GRID, 32),
                  bycells(DGG.PartialGrid(REGION), 8),
                  bycells(DGG.levelgrid(DGG.IGeo7System(), 2), 8))
        grid = gridof(space)
        @test GR._nchunks(space) > 1
        # Chunks partition the cells, in ascending index order, and each one
        # is contiguous — which is what lets `_ownedindices` be a range and a
        # chunk be one read.
        @test chunkcells(space) == 1:DGG.ncells(space)
        @test all(GR._ownedindices(space, c) isa UnitRange for c in 1:GR._nchunks(space))
        @test all(GR.chunkranges(space, c, (DGG.ncells(space),)) ==
                  (GR._ownedindices(space, c),) for c in 1:GR._nchunks(space))
        # A chunk's cap covers every boundary vertex of every cell in it: the
        # covering law the lazy path's chunk discovery prunes on.
        @test capscover(space)
        # A chunk of a grid is the grid's own cursor restricted to the chunk,
        # rooted at the chunk's node, with the sweep leaf size.
        chunktree = GR._subtree(space, GR._ownedindices(space, 2))
        @test chunktree isa GrI.MaskedCursor{<:DGG.HierarchicalGridCursor}
        @test chunktree.node.bucket_size == GrI.leaf_capacity(grid, GrI.ForSweeps())
        @test chunktree.cells === GR._ownedindices(space, 2)
        @test GOCore.best_manifold(space) == GOCore.best_manifold(SRC)
        i = DGG.ncells(space) ÷ 2
        @test GrI.cellat(space, GrI.cell_centroid(space, i)) == i
        # The grid's local answer is its own membership read as an index, and
        # it resolves that membership once. What it has to agree with is the
        # two-step form — locate, then look the located cell up — at points
        # inside a cell, on a cell's boundary, and off the collection
        # altogether. `localindex` on the grid is the same question.
        let twostep = function (p)
                c = DGG.cellat(grid, p)
                return c === nothing ? nothing : DGG.localindex(grid, c)
            end,
            probes = vcat(
                [GrI.cell_centroid(space, j) for j in (1, i, DGG.ncells(space))],
                collect(DGG.cell_boundary(grid, DGG.cellindex(grid, i))),
                collect(DGG.cell_boundary(grid, DGG.cellindex(grid, 1))),
                [DGG.Fallbacks.unit_point(lon, lat)
                 for (lon, lat) in ((170.0, -85.0), (0.0, 0.0), (-179.0, 12.0))])

            @test all(GrI.cellat(space, p) == twostep(p) for p in probes)
            @test all(DGG.localindex(grid, p) == twostep(p) for p in probes)
            @test (@inferred Union{Nothing,Int} GrI.cellat(space, first(probes))) ==
                  twostep(first(probes))
            # A subset does not cover the sphere, so some of those points are
            # off it and the agreement above is an agreement about `nothing`
            # as well as about indices.
            if grid isa DGG.PartialGrid
                @test any(GrI.cellat(space, p) === nothing for p in probes)
            end
        end
        # The sites a point method interpolates between are the grid's own
        # centroids, computed on read: nothing is materialised to prepare a
        # sampler, however many cells the grid has.
        sites = GR.samplesites(space)
        @test sites isa GR.CentroidSites
        @test Base.summarysize(sites) - Base.summarysize(space) < 64
        @test all(sites[j] == GrI.cell_centroid(space, j)
                  for j in (1, i, DGG.ncells(space)))
        # `_chunkat` inverts `_ownedindices`: every cell is placed back in the
        # chunk it came from. A subset's windows are the ones that can
        # disagree, since they are the ancestor's descendant range intersected
        # with the grid.
        @test all(GR._chunkat(space, j) == c
                  for c in 1:GR._nchunks(space) for j in GR._ownedindices(space, c))
        # The chunk index is GlobalRegridding's packed R-tree over the chunk
        # extents. A whole-sphere query reaches every chunk exactly once, and
        # every chunk's own covering cap reaches itself.
        index = GR.chunkindex(space)
        @test index isa GR.FlexibleRTrees.RTree
        whole = DGG.Fallbacks.full_sphere_cap()
        @test GR.candidatechunks!(Int[], index, whole) == collect(1:GR._nchunks(space))
        @test all(c in GR.candidatechunks!(Int[], index, GR._chunkextents(space)[c])
                  for c in 1:GR._nchunks(space))
        # The point form is `cellat` composed with it, and answers nothing
        # outside the grid's coverage exactly as `cellat` does.
        @test GR._chunkat(space, GrI.cell_centroid(space, i)) == GR._chunkat(space, i)
        @test_throws BoundsError GR._chunkat(space, DGG.ncells(space) + 1)
    end
    # No sorted subtrees, so no ancestor level to chunk by, and one chunk holds
    # everything rather than the grid refusing to be a regrid side.
    a5 = DGG.levelgrid(DGG.A5System(), 2)
    @test GrI.chunks(a5) === nothing
    @test GR._nchunks(a5) == 1
    @test GR._ownedindices(a5, 1) == 1:DGG.ncells(a5)

    # A system that answers the level-grid contract with a grid type of its own
    # rather than with `ncells(sys, l)` — the escape hatch `AbstractGrid`
    # documents, and the one shipped system that takes it. Sizing a chunk level
    # and area-matching a level both have to go through `levelgrid` to see it.
    auth = DGG.AuthalicSystem(DGG.IGeo7System())
    authspace = bycells(DGG.levelgrid(auth, 3), 32)
    @test GR._nchunks(authspace) > 1
    @test chunkcells(authspace) == 1:DGG.ncells(authspace)
    @test DGG.levelfor(auth, SRC) == DGG.levelfor(DGG.IGeo7System(), SRC)
end

@testset "the level chunks are the existing hierarchy's frontier" begin
    function frontier!(ids, ranges, node, chunklevel)
        if node.level >= chunklevel
            push!(ids, node.id)
            push!(ranges, DGG.Engine.node_indices(node))
            return
        end
        for child in STI.getchild(node)
            frontier!(ids, ranges, child, chunklevel)
        end
    end

    # `chunks(grid, a)` is the level-`a` frontier of the grid's own hierarchy:
    # one rooted chunk per non-empty ancestor, rooted at that ancestor and
    # owning its descendant positions, in frontier order.
    for grid in (DGG.levelgrid(DGG.IGeo7System(), 3), DGG.PartialGrid(REGION))
        space = bylevel(grid, 2)
        ch = GrI.chunks(space)
        ids = typeof(first(ch).root)[]
        ranges = AbstractVector{Int}[]
        frontier!(ids, ranges, DGG.HierarchicalGridCursor(grid; bucket_size = 0), 2)
        @test ids == [c.root for c in ch]
        @test ranges == [GrI.parentpositions(c) for c in ch]
        @test GR._nchunks(space) == length(ids)
    end

    # CopernicusDEM's normal cell tree is a block cursor. Chunk discovery reads
    # only the chunks' extents, so it needs no adapter for that separate
    # cell-tree optimization.
    cop = DGG.CopernicusDEM.CopernicusDEMSystem{30}()
    root = first(DGG.rootcells(cop))
    copgrid = DGG.subtree(cop, root, 1)
    copspace = bylevel(copgrid, 1)
    @test !(DGG.treeify(copgrid) isa DGG.HierarchicalGridCursor)
    @test GR.candidatechunks!(Int[], GR.chunkindex(copspace),
              DGG.Fallbacks.full_sphere_cap()) == collect(1:GR._nchunks(copspace))
end

@testset "the dependency graph holds every pair the chunk index answers" begin
    # `consumersof` is a source chunk's refcount, so the graph has to hold every
    # pair `candidatechunks!` can answer with. Both relations being conservative
    # supersets of the true overlap is NOT enough: they can be conservative in
    # different directions. The chunk index is the packed R-tree over the very
    # caps `_chunkextents` reports, so here the graph, the demanded relation and
    # the cap join are one relation, and a Copernicus level-0 source — whose
    # block cursor once answered a node-rectangle cap half a pixel wider than
    # its tile's boundary-ring cap — must not cross it.
    #
    # WHAT THIS CAN AND CANNOT CATCH. The builder is itself one
    # `candidatechunks!` per destination cap on `chunkindex(src)`, so this
    # sweep compares `candidatechunks!` against a graph BUILT FROM
    # `candidatechunks!`.
    # It still catches a builder that mis-assembles, mis-sorts or loses rows —
    # which is most of what a CSR builder can get wrong — but it can no longer
    # catch a wrong *choice* of index, because both sides would move together.
    # The independent checks on that are `contributing_pairs` below, and the cap
    # join at the end of this testset.
    cop = DGG.CopernicusDEMSystem(90)
    sys7 = DGG.IGeo7System()
    src = bylevel(DGG.levelgrid(cop, 1), 0)
    dst = bylevel(DGG.levelgrid(sys7, 3), 2)
    @test GR._nchunks(src) == DGG.ncells(cop, 0)

    # The case includes both poles and all twelve pentagons — the places where
    # a cell's cap and its hierarchy's node extents diverge most.
    pentagons = [foldl((c, _) -> first(DGG.children(sys7, c)), 1:3; init = r)
                 for r in DGG.rootcells(sys7)]
    @test length(unique(GR._chunkat(dst, DGG.localindex(parent(dst), p))
                        for p in pentagons)) == 12
    for pole in (GO.UnitSphericalPoint(0.0, 0.0, 1.0),
                 GO.UnitSphericalPoint(0.0, 0.0, -1.0))
        @test GR._chunkat(dst, pole) in 1:GR._nchunks(dst)
    end

    graph = GR.chunk_dependency_graph(dst, src)
    index = GR.chunkindex(src)
    @test index isa GR.FlexibleRTrees.RTree
    caps = GR._chunkextents(dst)
    buf = Int[]
    demanded, unpredicted, reached = 0, 0, falses(GR._nchunks(src))
    for d in 1:GR._nchunks(dst)
        GR.candidatechunks!(buf, index, caps[d]; radius = 0.0)
        demanded += length(buf)
        row = GR.sourcesof(graph, d)
        unpredicted += count(s -> !insorted(Int32(s), row), buf)
        reached[buf] .= true
    end
    # A complete destination reaches every tile, polar rows included, so the
    # sweep really did exercise the whole source lattice.
    @test all(reached)
    @test demanded > 0
    @test unpredicted == 0

    # The independent brute-force cap join over the same chunk extents is the
    # graph exactly: the packed index neither prunes a pair whose caps meet nor
    # adds one whose caps do not.
    capjoin = capjoin_pairs(dst, src)
    native = graph_pairs(graph)
    @test native == capjoin
end

@testset "no geometrically contributing pair is dropped" begin
    # The test above compares the graph against another conservative relation;
    # both could be conservative about the wrong thing together. This one builds
    # REAL hexagon and pentagon geometry and asks the only question that
    # matters: is every chunk pair a weight could be nonzero on an edge?
    sys7 = DGG.IGeo7System()
    rasterspace = GR.RasterGrid(globalraster(10.0);
        chunks = ([lo:min(lo + 6, 36) for lo in 1:7:36],
            [lo:min(lo + 4, 18) for lo in 1:5:18]))
    rooted = bylevel(
        DGG.subtree(sys7, DGG.cellindex(DGG.levelgrid(sys7, 1), 20), 4), 2)
    sparse = bylevel(DGG.PartialGrid(REGION), 2)

    cases = (
        # A complete grid against a complete grid three levels finer: all twelve
        # pentagons, both poles, and hexagon children that reach outside their
        # own parent's boundary.
        ("complete IGeo7", bylevel(DGG.levelgrid(sys7, 2), 1),
            bylevel(DGG.levelgrid(sys7, 3), 2), 0.0),
        # Two systems whose cells share no edges at all.
        ("IGeo7 from S2", bylevel(DGG.levelgrid(sys7, 2), 1),
            bylevel(DGG.levelgrid(SYS, 3), 1), 0.0),
        # A rooted subtree destination: most source chunks reach nothing, so a
        # relation that lost a pair would leave no other trace.
        ("rooted subtree", rooted, bylevel(DGG.levelgrid(sys7, 2), 1), 0.0),
        # A scattered, non-rooted subset of a level.
        ("sparse subset", sparse, bylevel(DGG.levelgrid(SYS, 3), 1), 0.0),
        # The production shape: a chunked raster source into a DGG destination.
        ("raster source", bylevel(DGG.levelgrid(sys7, 2), 1), rasterspace, 0.0),
        ("raster source, support", bylevel(DGG.levelgrid(sys7, 2), 1),
            rasterspace, 0.05),
    )

    for (name, dst, src, radius) in cases
        @testset "$name" begin
            truth = contributing_pairs(dst, src; radius)
            graph = graph_pairs(GR.chunk_dependency_graph(dst, src; radius))
            @test !isempty(truth)
            @test truth ⊆ graph
            # `truth ⊆ graph` alone would pass on a builder that returned
            # every pair, so pin the relation exactly too: the rows ARE the
            # `candidatechunks!` answers, with no `refine`, so the graph and
            # the demanded relation are equal and not merely nested. That
            # equality is what a builder swap has to preserve.
            @test graph == demanded_pairs(dst, src; radius)
            # A chunk cap covers its own cells, so the cap join holds the same
            # pairs. This is the `_chunkextents` half of the obligation: it
            # fails if a chunk cap is too tight.
            @test truth ⊆ capjoin_pairs(dst, src; radius)
            # Every source's chunk index is the packed R-tree over those same
            # caps, whatever kind of grid it is, so the graph is the cap join
            # exactly — the raster's own quadtree no longer answers here.
            @test graph == capjoin_pairs(dst, src; radius)
        end
    end
end

@testset "a rooted subset chunks without scanning the level" begin
    # A rooted `PartialGrid` holds nothing outside its root's subtree, so its
    # level chunks are its root's own non-empty descendants — and the answer
    # has to be the one the same cells give unrooted, ancestor for ancestor and
    # range for range.
    sys = DGG.IGeo7System()
    root = DGG.cellindex(DGG.levelgrid(sys, 2), 40)
    grid = DGG.subtree(sys, root, 5)
    unrooted = DGG.PartialGrid(sys, 5, collect(DGG.CellVector(grid)))
    for a in 2:5
        narrow = GrI.chunks(bylevel(grid, a))
        wide = GrI.chunks(bylevel(unrooted, a))
        @test [c.root for c in narrow] == [c.root for c in wide]
        @test [GrI.parentpositions(c) for c in narrow] ==
              [GrI.parentpositions(c) for c in wide]
        @test GR._nchunks(bylevel(grid, a)) == 7^(a - 2)
    end
    # A root DEEPER than the chunk level: the whole grid sits under one
    # ancestor, which is the one chunk.
    deep = bylevel(DGG.subtree(sys, root, 4), 1)
    @test GR._nchunks(deep) == 1
    @test only(GrI.chunks(deep)).root == DGG.ancestor(sys, root, 1)
    @test GrI.parentpositions(only(GrI.chunks(deep))) == 1:(7^2)
    # An UNROOTED subset still partitions its cells.
    scattered = bylevel(DGG.PartialGrid(REGION), 2)
    @test GR._nchunks(scattered) > 1
    @test chunkcells(scattered) == 1:DGG.ncells(scattered)
end

@testset "every spelling of `to` names the same cells" begin
    set = DGG.query(SYS, DGG.MultiOrderCoverage(GLOBE); level = LEVEL)
    reference = DGG.regrid(RASTER; to = GRID)
    for target in (SYS, DGG.CellLookup(GRID), DGG.CellVector(GRID), set,
                   bylevel(GRID, 1))
        out = DGG.regrid(RASTER; to = target)
        @test parent(out) == parent(reference)
        @test collect(DD.lookup(out, 1)) == collect(DGG.CellVector(GRID))
    end

    # Each spelling is one `_asgrid` method answering the grid over exactly
    # those cells — the grid itself, not merely a matching cell count — and the
    # same method answers `from`, which is given no source grid to match.
    for (target, cells) in ((GRID, DGG.CellVector(GRID)),
                            (DGG.CellLookup(GRID), DGG.CellVector(GRID)),
                            (DGG.CellVector(GRID), DGG.CellVector(GRID)),
                            (set, DGG.CellVector(GRID)),
                            (DGG.PartialGrid(REGION), REGION),
                            (REGION, REGION))
        grid = GR._asgrid(target, "to", SRC)
        @test grid isa DGG.AbstractGrid
        target isa DGG.AbstractGrid && @test grid === target
        @test collect(DGG.CellVector(grid)) == collect(cells)
        @test collect(DGG.CellVector(GR._asgrid(target, "from"))) == collect(cells)
    end
    # The union and the open conversion generic that stood between a spelling
    # and its grid are gone, and so is the space wrapper; a spelling is a method.
    @test !isdefined(DGG, :regridgrid)
    @test !isdefined(DGG, :RegridTarget)
    @test !isdefined(DGG, :DGGSpace)

    # A bare system needs the source to choose a level, and any `from` that
    # names cells is a measurable source — the data itself need not carry them.
    bare = GR.plan_regrid(vec(parent(RASTER)[:, :, 1]); to = SYS, from = SRC)
    @test DGG.level(bare.dst_space) == DGG.levelfor(SYS, SRC) == LEVEL
    # And a DGG grid is one of those sources, not only a raster lattice.
    fromgrid = GR.plan_regrid(zeros(DGG.ncells(GRID)); to = SYS, from = GRID)
    @test DGG.level(fromgrid.dst_space) == LEVEL
    # The level follows the source whatever kind of grid the source is.
    @test DGG.level(GR._asgrid(SYS, "to", GRID)) == DGG.levelfor(SYS, GRID)
    # As a source a bare system has nothing to be matched against, and saying
    # so — naming the keyword and the spelling that fixes it — is the whole of
    # that spelling's behaviour.
    @test_throws ArgumentError GR._asgrid(SYS, "from")
    @test_throws "`from = S2System()` names no cells until a level is chosen" GR._asgrid(
        SYS, "from")
    @test_throws "you must name it with `levelgrid(sys, l)`" GR._asgrid(SYS, "from")
end

@testset "a bare system takes the size-matched level" begin
    # The rule, restated independently: the level whose cell size is closest in
    # ratio to the median source cell. `radius = 1` measures both in steradians.
    areas = sort!([GR.cellarea(SRC, i) for i in 1:GR.ncells(SRC)])
    median = (areas[length(areas) ÷ 2] + areas[length(areas) ÷ 2 + 1]) / 2
    closest = argmin(l -> abs(2 * log(DGG.cellsize(SYS, l; radius = 1.0)) -
                              log(median)), 0:8)
    @test DGG.levelfor(SYS, SRC) == closest == LEVEL
    # Ratio, not difference: a source four times as coarse drops a level.
    @test DGG.levelfor(SYS, GR.RasterGrid(globalraster(30.0))) == LEVEL - 1
end

@testset "the destination axis is the cells" begin
    out = DGG.regrid(RASTER; to = GRID)
    @test DD.dims(out, 1) isa DGG.Cells
    @test DD.lookup(out, 1) isa DGG.CellLookup
    @test collect(DD.lookup(out, 1)) == collect(DGG.CellVector(GRID))
    # Non-spatial dimensions pass through untouched, in order, after the cells.
    @test DD.dims(out, 2) == DD.dims(RASTER, :month)
    @test size(out) == (DGG.ncells(GRID), 2)
    sub = DGG.regrid(RASTER; to = REGION)
    @test collect(DD.lookup(sub, 1)) == collect(REGION)

    # That axis is the grid's `destinationdims` and nothing else. A plan whose
    # destination is a DGG grid carries no application of its own on either
    # route, so what labelled the result above is the generic output wrapping
    # every grid reaches.
    eagerplan = DGG.plan_regrid(RASTER; to = GRID)
    lazyplan = DGG.plan_regrid(RASTER; to = GRID, lazy = true)
    @test GR.destinationdims(eagerplan) == (DD.dims(out, 1),)
    @test GR.destinationdims(lazyplan) == GR.destinationdims(eagerplan)
    for plan in (eagerplan, lazyplan)
        applications = methods(GR.regrid, Tuple{Any,typeof(plan)})
        @test !isempty(applications) &&
              all(m -> parentmodule(m) === GR, applications)
    end
end

@testset "a raster spells the way back off the cells" begin
    onto = DGG.regrid(RASTER; to = GRID)

    # The return leg wants the raster's lattice, and the raster — or just the
    # dimensions it carries — is that lattice. Neither side has to reach past
    # it for the `RasterGrid` it stands for.
    reference = DGG.regrid(onto; to = SRC, from = GRID)
    for to in (RASTER, DD.dims(RASTER))
        back = DGG.regrid(onto; to, from = GRID)
        @test parent(back) == parent(reference)
        @test DD.dims(back) == DD.dims(reference)
    end

    # A `Cells` axis is not a lattice, so it is refused wherever it is given.
    # It is the same guard either way, and it names the keyword it was handed
    # and the grid the axis holds — not `from` and not `xdim`.
    named = DGG.cellset(DD.lookup(onto, 1))
    @test_throws "`to` was given a dimensional raster" DGG.regrid(RASTER; to = onto)
    @test_throws "to = $named" DGG.regrid(RASTER; to = onto)
    @test_throws "from = $named" DGG.regrid(onto; to = GRID)
end

@testset "Extensive conserves the global integral" begin
    out = DGG.regrid(RASTER; to = GRID, missingpolicy = GR.Extensive())
    for m in 1:2
        slice = view(parent(RASTER), :, :, m)
        total = sum(GR.cellarea(SRC, i) * slice[i] for i in 1:GR.ncells(SRC))
        @test sum(view(parent(out), :, m)) ≈ total rtol = 1e-10
    end
end

@testset "a source's declared sentinel reaches the plan" begin
    # This package resolves `to` and supplies grids; every keyword beyond that
    # is `plan_regrid`'s, defaults included. A sentinel the source declares of
    # itself is one of those defaults, so the same field with the sentinel
    # spelled NaN is the oracle for it arriving.
    holed = collect(parent(RASTER))
    holed[3, 4, :] .= -9999.0
    nanned = replace(holed, -9999.0 => NaN)
    ds = DD.dims(RASTER)
    declared = DD.DimArray(holed, ds; metadata = DD.Metadata(Dict("_FillValue" => -9999.0)))

    @test GR.sourcemissingval(declared) == -9999.0
    @test DGG.plan_regrid(declared; to = GRID).missingval == -9999.0
    @test isequal(parent(DGG.regrid(declared; to = GRID)),
        parent(DGG.regrid(DD.DimArray(nanned, ds); to = GRID)))
    # And the caller still overrides what the source says.
    @test DGG.plan_regrid(declared; to = GRID, missingval = nothing).missingval === nothing
end

@testset "IGeo7 caps stay analytical through the authalic wrapper" begin
    sys = DGG.IGeo7System()
    root = DGG.cellindex(DGG.levelgrid(sys, 2), 100)
    grids = (
        DGG.subtree(sys, root, 6),
        DGG.subtree(DGG.AuthalicSystem(sys), root, 6),
    )

    function firstleaf(tree)
        while !STI.isleaf(tree)
            tree = first(STI.getchild(tree))
        end
        return tree
    end

    for grid in grids
        space = bylevel(grid, 2)
        # The whole space's tree is the grid's own sweep cursor, with the sweep
        # leaf size; a leaf's entries are built once per visit, as a vector.
        tree = GR._subtree(space, 1:DGG.ncells(grid))
        @test tree isa DGG.HierarchicalGridCursor
        @test tree.grid === grid
        @test tree.bucket_size == GrI.leaf_capacity(grid, GrI.ForSweeps())

        leaf = firstleaf(tree)
        entries = STI.child_indices_extents(leaf)
        @test entries isa Vector
        @test 1 < length(entries) <= tree.bucket_size
        for (index, cap) in entries
            c = DGG.cellindex(grid, index)
            @test all(p -> DGG.Fallbacks.cap_contains(cap, p),
                DGG.cell_boundary(grid, c))
        end
    end

    # The authalic wrapper warps the centre and inflates the radius by its
    # Lipschitz bound; it must not reuse an untransformed spherical cap.
    base, auth = grids
    c = DGG.cellindex(base, 173)
    basecap = DGG.Fallbacks.cell_cap(base, c)
    authcap = DGG.Fallbacks.cell_cap(auth, c)
    @test authcap.point == DGG.cell_centroid(auth, c)
    @test authcap.radius == nextfloat(
        DGG.Fallbacks.authalic_stretch(auth.complete.transform) * basecap.radius)
end

@testset "chunk trees give the whole grid's weights" begin
    leafindices(t) = sort!(GrI.subtree_indices(t))

    m = GOCore.best_manifold(GRID)
    op = CR.DefaultIntersectionOperator(m)
    src = GR._subtree(SRC, 1:GR.ncells(SRC))
    areas(space, tree) = CR.intersection_areas(m, GOCore.False(), space, SRC;
        dst_tree = tree, src_tree = src, intersection_operator = op)

    for space in (bycells(GRID, 32), bycells(DGG.PartialGrid(REGION), 8))
        grid = gridof(space)
        n = DGG.ncells(space)
        # The whole grid's sweep tree is the cursor `ForSweeps` names, built in
        # O(1); the one-cell cursor is the finest tree. Neither keeps a cap.
        sweep = GR._subtree(space, 1:n)
        @test sweep isa DGG.HierarchicalGridCursor
        @test sweep.bucket_size == GrI.leaf_capacity(space, GrI.ForSweeps())
        @test GR._celltree(space) == sweep
        bare = DGG.treeify(space, 1)
        @test bare.bucket_size == 1
        # The leaf size changes no weight.
        reference = areas(space, bare)
        @test areas(space, sweep) == reference

        # A chunk's tree names the chunk's cells exactly, and the chunks' weight
        # matrices sum to the whole grid's, entry for entry: each entry is
        # measured in exactly one chunk.
        total = spzeros(size(reference)...)
        for c in 1:GR._nchunks(space)
            inds = GR._ownedindices(space, c)
            tree = GR._subtree(space, inds)
            @test leafindices(tree) == collect(inds)
            total += areas(space, tree)
        end
        @test total == reference
    end
end

@testset "the sweep leaf size changes no candidate pair" begin
    @test GrI.leaf_capacity(GRID, GrI.ForSweeps()) > 1

    leaves(n) = STI.isleaf(n) ? [collect(STI.child_indices_extents(n))] :
                reduce(vcat, (leaves(c) for c in STI.getchild(n)); init = Vector{Any}())

    n = DGG.ncells(GRID)
    sweep = GR._subtree(GRID, 1:n)
    @test sweep.bucket_size == GrI.leaf_capacity(GRID, GrI.ForSweeps())
    # The leaves really did grow, and they still name every cell exactly once.
    ls = leaves(sweep)
    @test maximum(length, ls) > 1
    @test all(length(l) <= sweep.bucket_size for l in ls)
    @test sort!([i for l in ls for (i, _) in l]) == 1:n

    # The grid's own `treeify` is the one-cell cursor.
    bare = DGG.treeify(GRID)
    @test bare.bucket_size == 0
    @test length(leaves(bare)) == n

    # A system without sorted subtrees has no window cursor; a chunk of it
    # restricts the grid's selection cursor, which still names its cells.
    a5 = DGG.levelgrid(DGG.A5System(), 1)
    a5inds = 1:(DGG.ncells(a5) ÷ 2)
    a5tree = GR._subtree(a5, a5inds)
    @test a5tree isa GrI.MaskedCursor{<:DGG.Engine.SelectionCursor}
    @test sort!(GrI.subtree_indices(a5tree)) == collect(a5inds)

    # A rooted chunk of a level grid is the grid's own cursor restricted to the
    # chunk and rooted at the chunk's node, with the sweep leaf size.
    space = bycells(GRID, 32)
    chunk = GR._subtree(space, GR._ownedindices(space, 2))
    @test chunk isa GrI.MaskedCursor{<:DGG.HierarchicalGridCursor}
    @test chunk.node.bucket_size == GrI.leaf_capacity(GRID, GrI.ForSweeps())
    @test DGG.Engine.node_indices(chunk.node) == GR._ownedindices(space, 2)

    # The candidate pairs a dual search collects are the same set either way —
    # a node's cap covers its descendants', so stopping early can neither add
    # nor drop a pair.
    src = GR._subtree(SRC, 1:GR.ncells(SRC))
    function pairs(dst)
        out = Tuple{Int,Int}[]
        STI.dual_depth_first_search(GO.Extents.intersects, dst, src) do i, j
            push!(out, (i, j))
        end
        return sort!(out)
    end
    @test pairs(sweep) == pairs(bare)
end

@testset "a plan, and the lazy array, give the bare answer" begin
    reference = DGG.regrid(RASTER; to = GRID)
    plan = DGG.plan_regrid(RASTER; to = GRID)
    applied = DGG.regrid(RASTER, plan)
    @test parent(applied) == parent(reference)
    @test DD.dims(applied, 1) isa DGG.Cells
    lazy = DGG.regrid(RASTER; to = GRID, lazy = true)
    # Cells and the pass-through month, with their lookups, on the lazy route too.
    @test DD.dims(lazy) == DD.dims(reference)
    # The cells are the whole of the destination's shape, so a lazy result needs
    # no reshaping view over the array that computes it.
    @test parent(lazy) isa GR.LazyRegridArray
    @test Array(parent(lazy)) == parent(reference)
end

@testset "a point method lands on the same axis, eagerly and lazily" begin
    # A point method evaluates at the destination's sample sites rather than
    # averaging over its cells, and asks the output for `Points` sampling where
    # `Conservative` asks for `Intervals`. A cell holds one value either way, so
    # the axis is the same axis, and every route fills it identically.
    #
    # The second source is the same lattice under a different chunking. A
    # stencil is found whole before it is partitioned, so no chunking may move
    # a value: `NearestCell` takes one source value entire and must reproduce
    # it bit for bit, while a stencil of several sites is summed in as many
    # partial sums as there are blocks and may reassociate by an ulp.
    #
    # A destination the source cannot map is missing on every route, so the
    # looser comparison still holds one missing equal to another.
    approxsame(a, b) = length(a) == length(b) &&
        all(isequal(x, y) || x ≈ y for (x, y) in zip(a, b))
    chunked = GR.RasterGrid(DD.dims(RASTER); chunks = ([1:12, 13:24], [1:6, 7:12]))
    @test GR._nchunks(chunked) > GR._nchunks(SRC)
    for (method, same) in ((GR.BarycentricPoint(), approxsame),
        (GR.NearestCell(), isequal))

        point(; kwargs...) = DGG.regrid(RASTER; to = GRID, method, kwargs...)
        eager = point()
        lazy = point(; lazy = true)
        @test DD.dims(eager, 1) == DD.dims(DGG.regrid(RASTER; to = GRID), 1)
        @test DD.dims(lazy) == DD.dims(eager)
        @test parent(lazy) isa GR.LazyRegridArray
        @test isequal(Array(parent(lazy)), parent(eager))

        rechunked = DGG.regrid(RASTER; to = GRID, from = chunked, method,
            lazy = true)
        @test same(Array(parent(rechunked)), parent(eager))
    end
end

@testset "a column adopts the whole covering's relation and reads the same" begin
    # The shape `scripts/copdem_production.jl` actually has: a covering plan
    # over many level-`ancestor` chunks, and a per-column regrid whose
    # destination is a ROOTED one-chunk subtree grid of the same hierarchy.
    #
    # A row view of the covering's relation is refused there, because the view
    # stamps the whole covering. `GR.subspace_dependencies` re-stamps it onto
    # the column's own grid, and this asserts the two things that make it worth
    # having: the plan accepts it, and the values do not move.
    sys = DGG.IGeo7System()
    level, ancestor = 4, 2
    covering = bylevel(DGG.levelgrid(sys, level), ancestor)
    graph = GR.chunk_dependency_graph(covering, SRC)
    @test GR.hasextents(graph)

    a2 = DGG.levelgrid(sys, ancestor)
    for d in (1, 23, GR._nchunks(covering))
        column = bylevel(DGG.subtree(sys, DGG.cellindex(a2, d), level), ancestor)
        @test GR._nchunks(column) == 1

        # The column's chunk cap IS the covering's cap for that chunk, which is
        # what makes the re-stamp sound and what it checks.
        @test only(GR._chunkextents(column)) == GR._chunkextents(covering)[d]
        view = GR.subspace_dependencies(graph, column, [d])
        @test collect(GR.sourcesof(view, 1)) == collect(GR.sourcesof(graph, d))
        @test GR.destinationchunk(view, 1) == d

        # It is the relation the column would have built for itself...
        own = GR.chunk_dependency_graph(column, SRC)
        @test graph_pairs(view) == graph_pairs(own)
        @test GR.dependency_identity(view) == GR.dependency_identity(own)
        # ...and the plan adopts it by reference, where a plain row view and the
        # whole covering's relation are both still refused.
        adopted = GR.plan_regrid(RASTER; to = column, from = SRC, lazy = true,
            dependencies = view)
        @test GR.dependencies(adopted) === view
        @test_throws ArgumentError GR.plan_regrid(RASTER; to = column, from = SRC,
            lazy = true, dependencies = GR.restrict(graph, [d]))
        @test_throws ArgumentError GR.plan_regrid(RASTER; to = column, from = SRC,
            lazy = true, dependencies = graph)

        # The values: adopted view, own relation, and the eager whole-domain
        # answer for the same column, all the same numbers.
        built = GR.plan_regrid(RASTER; to = column, from = SRC, lazy = true)
        fromview = Array(parent(DGG.regrid(RASTER, adopted)))
        frombuilt = Array(parent(DGG.regrid(RASTER, built)))
        @test fromview == frombuilt
        eager = Array(parent(DGG.regrid(RASTER; to = column, from = SRC, lazy = false)))
        @test fromview ≈ eager rtol = 1e-12
    end
end

end # module RegridTests
