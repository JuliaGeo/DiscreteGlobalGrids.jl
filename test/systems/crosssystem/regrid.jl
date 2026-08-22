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
    [collect(DGG.cellindices(space, c)) for c in 1:GR.nchunks(space)])

capscover(space) = all(DGG.Fallbacks.cap_contains(GR.chunkextents(space)[c], p)
    for c in 1:GR.nchunks(space)
    for i in DGG.cellindices(space, c)
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
        # Chunks partition the cells, in ascending position order, and each one
        # is contiguous — which is what lets `cellindices` be a range and a
        # chunk be one read.
        @test chunkcells(space) == 1:DGG.ncells(space)
        @test all(DGG.cellindices(space, c) isa UnitRange for c in 1:GR.nchunks(space))
        @test all(GR.chunkranges(space, c, (DGG.ncells(space),)) ==
                  (DGG.cellindices(space, c),) for c in 1:GR.nchunks(space))
        # A chunk's cap covers every boundary vertex of every cell in it: the
        # covering law the lazy path's chunk discovery prunes on.
        @test capscover(space)
        # A chunk keeps the hierarchy rather than falling back to a cap list,
        # and small enough chunks carry their leaf caps with them.
        chunktree = GR.subtree(space, DGG.cellindices(space, 2))
        @test chunktree isa DGG.CapCachedTree
        @test chunktree.node isa DGG.HierarchicalGridCursor
        @test GR.manifold(space) == GR.manifold(SRC)
        i = DGG.ncells(space) ÷ 2
        @test GR.cellat(space, GR.cellcentroid(space, i)) == i
        # `chunkat` inverts `cellindices`: every cell is placed back in the
        # chunk it came from, by binary search over the windows rather than by
        # a scan. A subset's windows are the ones that can disagree, since they
        # are the ancestor's descendant range intersected with the grid.
        @test all(GR.chunkat(space, j) == c
                  for c in 1:GR.nchunks(space) for j in DGG.cellindices(space, c))
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
    @test DGG.cellindices(a5, 1) == 1:DGG.ncells(a5)

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
    # A bare system needs the source to choose a level, and any `from` that
    # names cells is a measurable source — the data itself need not carry them.
    bare = GR.plan_regrid(vec(parent(RASTER)[:, :, 1]); to = SYS, from = SRC)
    @test DGG.level(bare.dst_space.grid) == DGG.levelfor(SYS, SRC) == LEVEL
    # And a grid is one of those spellings, not only a `RegridSpace`.
    fromgrid = GR.plan_regrid(zeros(DGG.ncells(GRID)); to = SYS, from = GRID)
    @test DGG.level(fromgrid.dst_space.grid) == LEVEL
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
    # `to = ` resolution goes through a `plan_regrid` method of this package's
    # own, which forwards the rest of the keywords. A sentinel the source
    # declares of itself is one of the defaults that forwarding must not
    # swallow, so the same field with the sentinel spelled NaN is the oracle.
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

@testset "a shifted cap vector is addressed by global position" begin
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

    # A chunk's tree caches only its own positions, addressed by global
    # position, and must answer the raw chunk cursor's caps and cells the same.
    for space in (DGG.DGGSpace(GRID; chunkcells = 32),
                  DGG.DGGSpace(DGG.PartialGrid(REGION); chunkcells = 8))
        for c in (1, GR.nchunks(space) ÷ 2, GR.nchunks(space))
            inds = DGG.cellindices(space, c)
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
        STI.dual_depth_first_search(GO.UnitSpherical._intersects, dst, src) do i, j
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
    @test DD.dims(lazy, 1) isa DGG.Cells
    @test Array(parent(lazy)) == parent(reference)
end

end # module RegridTests
