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
        # A chunk keeps the hierarchy rather than falling back to a cap list.
        @test GR.subtree(space, DGG.cellindices(space, 2)) isa DGG.HierarchicalGridCursor
        @test GR.manifold(space) == GR.manifold(SRC)
        i = DGG.ncells(space) ÷ 2
        @test GR.cellat(space, GR.cellcentroid(space, i)) == i
        # `chunkat` inverts `cellindices`: every cell is placed back in the
        # chunk it came from, by binary search over the windows rather than by
        # a scan. A subset's windows are the ones that can disagree, since they
        # are the ancestor's descendant range intersected with the grid.
        @test all(GR.chunkat(space, j) == c
                  for c in 1:GR.nchunks(space) for j in DGG.cellindices(space, c))
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
