# The regridding face: a cell collection as a `GlobalRegridding.RegridSpace`,
# the target spellings that resolve into one, and the cube a result comes back
# as. One coarse S2 destination stands in for every system — cell geometry is
# each system's own business and is tested elsewhere — so what is under test
# here is the space contract, the `to` resolution, the `Cells` axis, and that
# mass survives the DGG side of a conservative regrid.

module RegridTests

using Test
import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import DimensionalData as DD
import Extents
import GeometryOps as GO
using GeometryOps: SpatialTreeInterface as STI
import ConservativeRegridding as CR
import ConservativeRegridding: Trees
import GeometryOpsCore as GOCore

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

# The cells the space assigns to each chunk, concatenated in chunk order.
chunkcells(space) = reduce(vcat,
    [collect(GR.ownedindices(space, c)) for c in 1:GR.nchunks(space)])

capscover(space) = all(DGG.Fallbacks.cap_contains(GR.chunkextents(space)[c], p)
    for c in 1:GR.nchunks(space)
    for i in GR.ownedindices(space, c)
    for p in DGG.cell_boundary(space.grid, DGG.cellindex(space.grid, i)))

const REGION = DGG.covering(DGG.CellVector(GRID),
    Extents.Extent(X = (-140.0, 40.0), Y = (-20.0, 60.0)))

@testset "the space contract" begin
    # A complete level, a scattered subset of one, and a hexagonal system whose
    # aperture-7 children reach outside their parent's own boundary — the case a
    # chunk cap taken from the ancestor cell rather than from its subtree extent
    # would fail to cover.
    for space in (DGG.DGGSpace(GRID; chunkcells = 32),
                  DGG.DGGSpace(DGG.PartialGrid(REGION); chunkcells = 8),
                  DGG.DGGSpace(DGG.levelgrid(DGG.IGeo7System(), 2); chunkcells = 8))
        @test GR.nchunks(space) > 1
        # Chunks partition the cells, in ascending index order, and each one
        # is contiguous — which is what lets `ownedindices` be a range and a
        # chunk be one read.
        @test chunkcells(space) == 1:DGG.ncells(space)
        @test all(GR.ownedindices(space, c) isa UnitRange for c in 1:GR.nchunks(space))
        @test all(GR.chunkranges(space, c, (DGG.ncells(space),)) ==
                  (GR.ownedindices(space, c),) for c in 1:GR.nchunks(space))
        # A chunk's cap covers every boundary vertex of every cell in it: the
        # covering law the lazy path's chunk discovery prunes on.
        @test capscover(space)
        # A chunk keeps the hierarchy rather than falling back to a cap list,
        # and small enough chunks carry their leaf caps with them.
        chunktree = GR.subtree(space, GR.ownedindices(space, 2))
        @test chunktree isa DGG.CapCachedTree
        @test chunktree.node isa DGG.HierarchicalGridCursor
        @test GR.manifold(space) == GR.manifold(SRC)
        i = DGG.ncells(space) ÷ 2
        @test GR.cellat(space, GR.cellcentroid(space, i)) == i
        # The space's answer is the grid's own membership read as an index, and
        # it resolves that membership once. What it has to agree with is the
        # two-step form — locate, then look the located cell up — at points
        # inside a cell, on a cell's boundary, and off the collection
        # altogether. `localindex` on the grid is the same question, and the
        # method the space now takes.
        let twostep = function (p)
                c = DGG.cellat(space.grid, p)
                return c === nothing ? nothing : DGG.localindex(space.grid, c)
            end,
            probes = vcat(
                [GR.cellcentroid(space, j) for j in (1, i, DGG.ncells(space))],
                collect(DGG.cell_boundary(space.grid, DGG.cellindex(space.grid, i))),
                collect(DGG.cell_boundary(space.grid, DGG.cellindex(space.grid, 1))),
                [DGG.Fallbacks.unit_point(lon, lat)
                 for (lon, lat) in ((170.0, -85.0), (0.0, 0.0), (-179.0, 12.0))])

            @test all(GR.cellat(space, p) == twostep(p) for p in probes)
            @test all(DGG.localindex(space.grid, p) == twostep(p) for p in probes)
            @test (@inferred Union{Nothing,Int} GR.cellat(space, first(probes))) ==
                  twostep(first(probes))
            # A subset does not cover the sphere, so some of those points are
            # off it and the agreement above is an agreement about `nothing`
            # as well as about indices.
            if space.grid isa DGG.PartialGrid
                @test any(GR.cellat(space, p) === nothing for p in probes)
            end
        end
        # The sites a point method interpolates between are the space's own
        # centroids, computed on read: nothing is materialised to prepare a
        # sampler, however many cells the space has.
        sites = GR.samplesites(space)
        @test sites isa GR.CentroidSites
        @test Base.summarysize(sites) - Base.summarysize(space) < 64
        @test all(sites[j] == GR.cellcentroid(space, j)
                  for j in (1, i, DGG.ncells(space)))
        # `chunkat` inverts `ownedindices`: every cell is placed back in the
        # chunk it came from, by binary search over the windows rather than by
        # a scan. A subset's windows are the ones that can disagree, since they
        # are the ancestor's descendant range intersected with the grid.
        @test all(GR.chunkat(space, j) == c
                  for c in 1:GR.nchunks(space) for j in GR.ownedindices(space, c))
        # The DGG space itself is the private query index. A whole-sphere query
        # reaches every chunk exactly once through the original grid hierarchy,
        # and every chunk's own covering cap reaches itself.
        @test GR.chunkindex(space) === space
        whole = DGG.Fallbacks.full_sphere_cap()
        @test GR.candidatechunks!(Int[], space, whole) == collect(1:GR.nchunks(space))
        @test all(c in GR.candidatechunks!(Int[], space, GR.chunkextents(space)[c])
                  for c in 1:GR.nchunks(space))
        # The point form is `cellat` composed with it, and answers nothing
        # outside the space's coverage exactly as `cellat` does.
        @test GR.chunkat(space, GR.cellcentroid(space, i)) == GR.chunkat(space, i)
        @test_throws BoundsError GR.chunkat(space, DGG.ncells(space) + 1)
    end
    # No sorted subtrees, so no ancestor level to chunk by, and one chunk holds
    # everything rather than the space refusing to exist.
    a5 = DGG.DGGSpace(DGG.levelgrid(DGG.A5System(), 2))
    @test GR.nchunks(a5) == 1
    @test GR.ownedindices(a5, 1) == 1:DGG.ncells(a5)

    # A system that answers the level-grid contract with a grid type of its own
    # rather than with `ncells(sys, l)` — the escape hatch `AbstractGrid`
    # documents, and the one shipped system that takes it. Sizing a chunk level
    # and area-matching a level both have to go through `levelgrid` to see it.
    auth = DGG.AuthalicSystem(DGG.IGeo7System())
    authspace = DGG.DGGSpace(DGG.levelgrid(auth, 3); chunkcells = 32)
    @test GR.nchunks(authspace) > 1
    @test chunkcells(authspace) == 1:DGG.ncells(authspace)
    @test DGG.levelfor(auth, SRC) == DGG.levelfor(DGG.IGeo7System(), SRC)
end

@testset "the DGG chunk index is the existing hierarchy's frontier" begin
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

    for grid in (DGG.levelgrid(DGG.IGeo7System(), 3), DGG.PartialGrid(REGION))
        space = DGG.DGGSpace(grid; chunklevel = 2)
        ids = eltype(space.chunkids)[]
        ranges = AbstractVector{Int}[]
        frontier!(ids, ranges, DGG.HierarchicalGridCursor(grid; bucket_size = 0),
            space.chunklevel)
        @test ids == space.chunkids
        @test ranges == space.ranges
    end

    # CopernicusDEM's normal cell tree is a block cursor. Chunk discovery still
    # reaches the system hierarchy directly and therefore needs no adapter for
    # that separate cell-tree optimization.
    cop = DGG.CopernicusDEM.CopernicusDEMSystem{30}()
    root = first(DGG.rootcells(cop))
    copgrid = DGG.subtree(cop, root, 1)
    copspace = DGG.DGGSpace(copgrid; chunklevel = 1)
    @test !(DGG.treeify(copgrid) isa DGG.HierarchicalGridCursor)
    @test GR.candidatechunks!(Int[], copspace, DGG.Fallbacks.full_sphere_cap()) ==
          collect(1:GR.nchunks(copspace))
end

@testset "the dependency graph holds every pair the chunk index answers" begin
    # `consumersof` is a source chunk's refcount, so the graph has to hold every
    # pair `candidatechunks!` can answer with. Both relations being conservative
    # supersets of the true overlap is NOT enough: they can be conservative in
    # different directions, and here they were. The level-0 Copernicus frontier
    # tests a cap derived from the block cursor's node rectangle, whose east
    # edge carries the whole tile width; `chunkextents` reports the cap derived
    # from the tile's own boundary ring, half a pixel narrower. Neither cap
    # contains the other, so a graph built from `chunkextents` directly CROSSED
    # the executor's relation — and a refcount that reaches zero early retires a
    # tile the next read is still going to ask for.
    #
    # WHAT THIS CAN AND CANNOT CATCH. The builder is itself one
    # `candidatechunks!` per destination cap on `chunkindex(src)`, so this
    # sweep compares `candidatechunks!` against a graph BUILT FROM
    # `candidatechunks!`.
    # It still catches a builder that mis-assembles, mis-sorts or loses rows —
    # which is most of what a CSR builder can get wrong — but it can no longer
    # catch a wrong *choice* of index, because both sides would move together.
    # The independent check on that is `contributing_pairs` below, and the cap
    # crossing at the end of this testset.
    cop = DGG.CopernicusDEMSystem(90)
    sys7 = DGG.IGeo7System()
    src = DGG.DGGSpace(DGG.levelgrid(cop, 1); chunklevel = 0)
    dst = DGG.DGGSpace(DGG.levelgrid(sys7, 3); chunklevel = 2)
    @test GR.nchunks(src) == DGG.ncells(cop, 0)

    # The case includes both poles and all twelve pentagons — the places where
    # a cell's cap and its hierarchy's node extents diverge most.
    pentagons = [foldl((c, _) -> first(DGG.children(sys7, c)), 1:3; init = r)
                 for r in DGG.rootcells(sys7)]
    @test length(unique(GR.chunkat(dst, DGG.localindex(dst.grid, p))
                        for p in pentagons)) == 12
    for pole in (GO.UnitSphericalPoint(0.0, 0.0, 1.0),
                 GO.UnitSphericalPoint(0.0, 0.0, -1.0))
        @test GR.chunkat(dst, pole) in 1:GR.nchunks(dst)
    end

    graph = GR.chunk_dependency_graph(dst, src)
    index = GR.chunkindex(src)
    caps = GR.chunkextents(dst)
    buf = Int[]
    demanded, unpredicted, reached = 0, 0, falses(GR.nchunks(src))
    for d in 1:GR.nchunks(dst)
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

    # The crossing itself, on a window CI can run: the whole level-0 Copernicus
    # frontier, no tile list and no download. The block cursor's relation and the
    # cap join differ in BOTH directions, so the cap join is not a bound on the
    # graph on this hierarchy either — the same fact
    # `benchmark/chunk_graph_gates.jl` measures at production scale (72 pairs
    # the index holds and the cap join rejects, on the GLO-90 x IGeo7-L12
    # pair). That production figure is harness-only; these two assertions are
    # the tested form of the finding.
    capjoin = capjoin_pairs(dst, src)
    native = graph_pairs(graph)
    @test !isempty(setdiff(native, capjoin))
    @test !isempty(setdiff(capjoin, native))
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
    rooted = DGG.DGGSpace(
        DGG.subtree(sys7, DGG.cellindex(DGG.levelgrid(sys7, 1), 20), 4);
        chunklevel = 2)
    sparse = DGG.DGGSpace(DGG.PartialGrid(REGION); chunklevel = 2)

    cases = (
        # A complete grid against a complete grid three levels finer: all twelve
        # pentagons, both poles, and hexagon children that reach outside their
        # own parent's boundary.
        ("complete IGeo7", DGG.DGGSpace(DGG.levelgrid(sys7, 2); chunklevel = 1),
            DGG.DGGSpace(DGG.levelgrid(sys7, 3); chunklevel = 2), 0.0, true),
        # Two systems whose cells share no edges at all.
        ("IGeo7 from S2", DGG.DGGSpace(DGG.levelgrid(sys7, 2); chunklevel = 1),
            DGG.DGGSpace(DGG.levelgrid(SYS, 3); chunklevel = 1), 0.0, true),
        # A rooted subtree destination: most source chunks reach nothing, so a
        # relation that lost a pair would leave no other trace.
        ("rooted subtree", rooted,
            DGG.DGGSpace(DGG.levelgrid(sys7, 2); chunklevel = 1), 0.0, true),
        # A scattered, non-rooted subset of a level.
        ("sparse subset", sparse,
            DGG.DGGSpace(DGG.levelgrid(SYS, 3); chunklevel = 1), 0.0, true),
        # The production shape: a chunked raster source into a DGG destination.
        # Its quadtree is NOT answering the caps `chunkextents` reports, so the
        # cap join is not a bound on it in either direction.
        ("raster source", DGG.DGGSpace(DGG.levelgrid(sys7, 2); chunklevel = 1),
            rasterspace, 0.0, false),
        ("raster source, support", DGG.DGGSpace(DGG.levelgrid(sys7, 2); chunklevel = 1),
            rasterspace, 0.05, false),
    )

    for (name, dst, src, radius, capbounded) in cases
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
            # pairs. This is the `chunkextents` half of the obligation: it fails
            # if a chunk cap is too tight.
            @test truth ⊆ capjoin_pairs(dst, src; radius)
            # Where a space's index descends the very hierarchy
            # `chunkextents` is derived from, the cap join is an UPPER bound on
            # the graph: the descent can prune a chunk whose own cap intersects,
            # never add one whose cap does not. That holds for
            # `_dggcandidatechunks!` and NOT for the raster quadtree, which
            # answers whole chunks for a straddling leaf without testing each
            # chunk's own cap — so for those cases the containment is asserted
            # in the other direction below, never here.
            capbounded && @test graph ⊆ capjoin_pairs(dst, src; radius)
        end
    end

    # The raster quadtree's relation and the cap join cross in BOTH directions,
    # so the cap join bounds the graph in neither. The CopernicusDEM level-0
    # frontier crosses the same way; that one is pinned in "the dependency graph
    # holds every pair the chunk index answers" above.
    rdst = DGG.DGGSpace(DGG.levelgrid(sys7, 2); chunklevel = 1)
    native = graph_pairs(GR.chunk_dependency_graph(rdst, rasterspace))
    capjoin = capjoin_pairs(rdst, rasterspace)
    @test !isempty(setdiff(native, capjoin))
    @test !isempty(setdiff(capjoin, native))
end

@testset "a rooted subset chunks without scanning the level" begin
    # `_chunkwindows` visits every level-`a` ancestor to find the non-empty
    # ones, which at production sizes is a scan of the whole level per space
    # built. A rooted `PartialGrid` holds nothing outside its root's subtree, so
    # the visit narrows to that root's own descendants — and the answer has to
    # be the same one, ancestor for ancestor and range for range.
    sys = DGG.IGeo7System()
    root = DGG.cellindex(DGG.levelgrid(sys, 2), 40)
    grid = DGG.subtree(sys, root, 5)
    unrooted = DGG.PartialGrid(sys, 5, collect(DGG.CellVector(grid)))
    for a in 2:5
        narrow = DGG.DGGSpace(grid; chunklevel = a)
        wide = DGG.DGGSpace(unrooted; chunklevel = a)
        @test narrow.chunkids == wide.chunkids
        @test narrow.ranges == wide.ranges
        @test GR.nchunks(narrow) == 7^(a - 2)
    end
    # A root DEEPER than the chunk level: the whole grid sits under one
    # ancestor, which is the one chunk, and finding it is arithmetic rather
    # than a scan.
    deep = DGG.DGGSpace(DGG.subtree(sys, root, 4); chunklevel = 1)
    @test GR.nchunks(deep) == 1
    @test only(deep.chunkids) == DGG.ancestor(sys, root, 1)
    @test only(deep.ranges) == 1:(7^2)
    # An UNROOTED subset still scans, because nothing bounds it, and still
    # partitions its cells.
    scattered = DGG.DGGSpace(DGG.PartialGrid(REGION); chunklevel = 2)
    @test GR.nchunks(scattered) > 1
    @test chunkcells(scattered) == 1:DGG.ncells(scattered)
end

@testset "every spelling of `to` names the same cells" begin
    set = DGG.query(SYS, DGG.MultiOrderCoverage(GLOBE); level = LEVEL)
    reference = DGG.regrid(RASTER; to = GRID)
    for target in (SYS, DGG.CellLookup(GRID), DGG.CellVector(GRID), set,
                   DGG.DGGSpace(GRID))
        out = DGG.regrid(RASTER; to = target)
        @test parent(out) == parent(reference)
        @test collect(DD.lookup(out, 1)) == collect(DGG.CellVector(GRID))
    end

    # Each spelling is one `_asspace` method answering the space over exactly
    # those cells — the grid itself, not merely a matching cell count — and the
    # same method answers `from`, which is given no source space to match.
    for (target, cells) in ((GRID, DGG.CellVector(GRID)),
                            (DGG.CellLookup(GRID), DGG.CellVector(GRID)),
                            (DGG.CellVector(GRID), DGG.CellVector(GRID)),
                            (set, DGG.CellVector(GRID)),
                            (DGG.PartialGrid(REGION), REGION),
                            (REGION, REGION))
        space = GR._asspace(target, "to", SRC)
        @test space isa DGG.DGGSpace
        @test collect(DGG.CellVector(space.grid)) == collect(cells)
        @test collect(DGG.CellVector(GR._asspace(target, "from").grid)) ==
              collect(cells)
    end
    # The union and the open conversion generic that stood between a spelling
    # and its space are gone; a spelling is a method.
    @test !isdefined(DGG, :regridgrid)
    @test !isdefined(DGG, :RegridTarget)

    # A bare system needs the source to choose a level, and any `from` that
    # names cells is a measurable source — the data itself need not carry them.
    bare = GR.plan_regrid(vec(parent(RASTER)[:, :, 1]); to = SYS, from = SRC)
    @test DGG.level(bare.dst_space.grid) == DGG.levelfor(SYS, SRC) == LEVEL
    # And a grid is one of those spellings, not only a `RegridSpace`.
    fromgrid = GR.plan_regrid(zeros(DGG.ncells(GRID)); to = SYS, from = GRID)
    @test DGG.level(fromgrid.dst_space.grid) == LEVEL
    # The level follows the source whatever kind of space the source is.
    dggsrc = DGG.DGGSpace(GRID)
    @test DGG.level(GR._asspace(SYS, "to", dggsrc).grid) == DGG.levelfor(SYS, dggsrc)
    # As a source a bare system has nothing to be matched against, and saying
    # so — naming the keyword and the spelling that fixes it — is the whole of
    # that spelling's behaviour.
    @test_throws ArgumentError GR._asspace(SYS, "from")
    @test_throws "`from = S2System()` names no cells until a level is chosen" GR._asspace(
        SYS, "from")
    @test_throws "you must name it with `levelgrid(sys, l)`" GR._asspace(SYS, "from")
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

    # That axis is the space's `destinationdims` and nothing else. A plan whose
    # destination is a DGG space carries no application of its own on either
    # route, so what labelled the result above is the generic output wrapping
    # every space reaches.
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

@testset "Extensive conserves the global integral" begin
    out = DGG.regrid(RASTER; to = GRID, missingpolicy = GR.Extensive())
    for m in 1:2
        slice = view(parent(RASTER), :, :, m)
        total = sum(GR.cellarea(SRC, i) * slice[i] for i in 1:GR.ncells(SRC))
        @test sum(view(parent(out), :, m)) ≈ total rtol = 1e-10
    end
end

@testset "a source's declared sentinel reaches the plan" begin
    # This package resolves `to` and supplies spaces; every keyword beyond that
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

@testset "a shifted cap vector is addressed by local index" begin
    v = DGG._ShiftedCaps(collect(10:19), 100)
    @test v isa AbstractVector{Int}
    @test length(v) == 10
    @test size(v) == (10,)
    @test axes(v, 1) == 101:110
    @test parent(v) === v.data
    @test v[101] == 10
    @test v[110] == 19
    @test [v[i] for i in axes(v, 1)] == 10:19
    @test_throws BoundsError v[100]
    @test_throws BoundsError v[111]
    # Offset zero is the whole-space case: a plain 1-based vector.
    @test axes(DGG._ShiftedCaps(collect(1:3), 0), 1) == 1:3
end

@testset "the cached trees cache caps without changing them" begin
    samecap(a, b) = a.point == b.point && a.radius == b.radius

    # Extents, leaf entries and polygons must be the raw cursor's, bit for bit.
    function checktree(a, b)
        samecap(STI.node_extent(a), STI.node_extent(b)) || return false
        STI.isleaf(a) == STI.isleaf(b) || return false
        if STI.isleaf(a)
            va, vb = STI.child_indices_extents(a), collect(STI.child_indices_extents(b))
            length(va) == length(vb) || return false
            all(x[1] == y[1] && samecap(x[2], y[2]) for (x, y) in zip(va, vb)) || return false
            return all(Trees.getcell(a, i) == Trees.getcell(b, i) for (i, _) in va)
        end
        return all(checktree(x, y) for (x, y) in zip(STI.getchild(a), STI.getchild(b)))
    end

    # A cached tree also carries the seam's leaf size, so the shape to compare
    # its caps against is a bare cursor bucketed the same way. The field copy is
    # written out here so the test pins the constructor's field order itself.
    reshape_leaves(c) = typeof(c)(c.grid, c.system, c.top_level, c.leaf_level,
        DGG._CACHED_BUCKET_SIZE, c.level, c.id, c.first_index, c.last_index,
        c.selection)

    for space in (DGG.DGGSpace(DGG.PartialGrid(REGION)), DGG.DGGSpace(GRID))
        n = DGG.ncells(space.grid)
        cached = GR.subtree(space, 1:n)
        @test cached isa DGG.CapCachedTree
        cursor = GR.celltree(space)
        @test Trees.ncells(cached) == Trees.ncells(cursor) == n
        raw = DGG.treeify(DGG._decodedgrid(space.grid))
        @test checktree(cached, reshape_leaves(raw))
        # The weight matrix the two trees build is identical, entry for entry —
        # and `cursor` here still has the default one-cell leaf, so this is also
        # the statement that the seam's leaf size changes no weight.
        m = GR.manifold(space)
        op = CR.DefaultIntersectionOperator(m)
        src = GR.subtree(SRC, 1:GR.ncells(SRC))
        wa = CR.intersection_areas(m, GOCore.False(), cached, src; intersection_operator = op)
        wb = CR.intersection_areas(m, GOCore.False(), cursor, src; intersection_operator = op)
        @test wa == wb
    end

    # A chunk's tree caches only its own indices, addressed by global
    # index, and must answer the raw chunk cursor's caps and cells the same.
    for space in (DGG.DGGSpace(GRID; chunkcells = 32),
                  DGG.DGGSpace(DGG.PartialGrid(REGION); chunkcells = 8))
        for c in (1, GR.nchunks(space) ÷ 2, GR.nchunks(space))
            inds = GR.ownedindices(space, c)
            cached = GR.subtree(space, inds)
            @test cached isa DGG.CapCachedTree
            @test length(cached.caps) == length(inds)
            @test axes(cached.caps, 1) == inds
            # `checktree` covers extents, leaf entries and polygons; a chunk
            # tree's `Trees.ncells` is its own count while its leaf indices are
            # global, so weights only come out of a block build's index maps.
            @test checktree(cached, reshape_leaves(DGG._chunkcursor(space, inds)))
        end
    end

    # A partial grid's ids are decoded once: the tree's grid stores a plain
    # vector with the same cells.
    space = DGG.DGGSpace(DGG.PartialGrid(REGION))
    cached = GR.subtree(space, 1:DGG.ncells(space.grid))
    @test cached.node.grid.ids isa Vector
    @test cached.node.grid.ids == collect(space.grid.ids)
end

@testset "the bigger leaf rides on the cap cache, and nowhere else" begin
    # A leaf of `_CACHED_BUCKET_SIZE` cells is only cheap because `caps` already
    # holds their extents; a bare cursor re-derives them on every visit, where
    # the same leaf size is a large loss. So the size is attached to the two
    # sites that return a `CapCachedTree`, and every plain-cursor return path
    # keeps whatever the grid asked for.
    @test DGG._CACHED_BUCKET_SIZE > 1

    leaves(n) = STI.isleaf(n) ? [collect(STI.child_indices_extents(n))] :
                reduce(vcat, (leaves(c) for c in STI.getchild(n)); init = Vector{Any}())

    space = DGG.DGGSpace(GRID)
    n = DGG.ncells(space.grid)
    cached = GR.subtree(space, 1:n)
    @test cached isa DGG.CapCachedTree
    @test cached.node.bucket_size == DGG._CACHED_BUCKET_SIZE
    # The leaves really did grow, and they still name every cell exactly once.
    ls = leaves(cached)
    @test maximum(length, ls) > 1
    @test all(length(l) <= DGG._CACHED_BUCKET_SIZE for l in ls)
    @test sort!([i for l in ls for (i, _) in l]) == 1:n

    # Bare path 1: `celltree` hands back the plain cursor.
    bare = GR.celltree(space)
    @test !(bare isa DGG.CapCachedTree)
    @test bare.bucket_size == 0
    @test length(leaves(bare)) == n

    # Bare path 2: a chunk too large to pay for its cap vector is returned
    # untouched — same object, same leaf size.
    @test DGG._cachedchunktree(bare, 1:(DGG._CHUNK_CAP_CACHE_MAX + 1)) === bare

    # Bare path 3: a system without sorted subtrees uses a selection cursor,
    # which the cache cannot index, so `_cachedcelltree` falls back.
    a5 = DGG.DGGSpace(DGG.levelgrid(DGG.A5System(), 1))
    a5tree = GR.subtree(a5, 1:DGG.ncells(a5.grid))
    @test !(a5tree isa DGG.CapCachedTree)
    @test a5tree.bucket_size == 0

    # A caller that names a leaf size keeps it: `0` is the grid default, not a
    # request, and is the only value the seam fills in.
    explicit = DGG.DGGSpace(DGG.subtree(SYS, first(DGG.rootcells(SYS)), LEVEL;
        bucket_size = 7))
    etree = GR.subtree(explicit, 1:DGG.ncells(explicit.grid))
    @test etree isa DGG.CapCachedTree
    @test etree.node.bucket_size == 7

    # The candidate pairs a dual search collects are the same set either way —
    # a node's cap covers its descendants', so stopping early can neither add
    # nor drop a pair.
    src = GR.subtree(SRC, 1:GR.ncells(SRC))
    function pairs(dst)
        out = Tuple{Int,Int}[]
        STI.dual_depth_first_search(GO.Extents.intersects, dst, src) do i, j
            push!(out, (i, j))
        end
        return sort!(out)
    end
    @test pairs(cached) == pairs(bare)
    # (the weight matrix off each of the two is compared entry for entry in the
    # testset above, which is the same statement one level further on.)
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
    @test GR.nchunks(chunked) > GR.nchunks(SRC)
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
    # the column's own space, and this asserts the two things that make it worth
    # having: the plan accepts it, and the values do not move.
    sys = DGG.IGeo7System()
    level, ancestor = 4, 2
    covering = DGG.DGGSpace(DGG.levelgrid(sys, level); chunklevel = ancestor)
    graph = GR.chunk_dependency_graph(covering, SRC)
    @test GR.hasextents(graph)

    a2 = DGG.levelgrid(sys, ancestor)
    for d in (1, 23, GR.nchunks(covering))
        column = DGG.DGGSpace(DGG.subtree(sys, DGG.cellindex(a2, d), level);
            chunklevel = ancestor)
        @test GR.nchunks(column) == 1

        # The column's chunk cap IS the covering's cap for that chunk, which is
        # what makes the re-stamp sound and what it checks.
        @test only(GR.chunkextents(column)) == GR.chunkextents(covering)[d]
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
