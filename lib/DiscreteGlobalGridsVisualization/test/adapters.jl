module VisualizationAdapterTests

using Test
import DiscreteGlobalGrids as DGG
import DiscreteGlobalGridsVisualization as DGGV
import Makie

struct StandaloneGrid <: DGG.AbstractGrid
    n::Int
end
DGG.ncells(grid::StandaloneGrid) = grid.n
DGG.cellindex(grid::StandaloneGrid, i::Int) =
    1 <= i <= grid.n ? DGG.LevelIndex(0, i - 1) : throw(BoundsError(grid, i))
DGG.cell_boundary(::StandaloneGrid, cell::DGG.AbstractCellIndex) =
    DGG.cell_boundary(DGG.levelgrid(DGG.HEALPixSystem(), 0), cell)
DGG.cell_centroid(::StandaloneGrid, cell::DGG.AbstractCellIndex) =
    DGG.cell_centroid(DGG.levelgrid(DGG.HEALPixSystem(), 0), cell)

@testset "standalone grid adapters" begin
    grid = StandaloneGrid(2)
    @test DGG.system(grid) === nothing
    for adapter in (DGGV.cellset, DGGV.cellregion)
        cells = adapter(grid)
        @test cells.source === grid
        @test eltype(cells.cells) == DGG.LevelIndex
        @test collect(cells.cells) == [DGG.cellindex(grid, i) for i in 1:2]
        @test DGG.cell_boundary(cells.source, first(cells.cells)) ==
              DGG.cell_boundary(grid, DGG.cellindex(grid, 1))
        emptycells = adapter(StandaloneGrid(0))
        @test length(emptycells) == 0
        @test isempty(collect(emptycells.cells))
        @test eltype(emptycells.cells) == DGG.AbstractCellIndex
    end
end

@testset "stored cell axes retain the resampling heights argument" begin
    grid = DGG.levelgrid(DGG.HEALPixSystem(), 0)
    stored = DGG.ChunkedCellLookup(DGG.cellaxis(DGG.ImplicitEncoding(), grid,
        DGG.ncells(grid)))
    heights = Float64.(1:length(stored))
    for cells in (stored, parent(stored))
        converted, zs = Makie.convert_arguments(DGGV.DGGResample, cells, heights)
        @test converted.cells === parent(stored)
        @test collect(converted.cells) == collect(stored)
        @test zs === heights
    end
    converted, zs = Makie.convert_arguments(DGGV.DGGResample,
        DGG.system(grid), collect(stored))
    @test converted.cells == collect(stored)
    @test zs === nothing
end

end
