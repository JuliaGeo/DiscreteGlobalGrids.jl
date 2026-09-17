# `cell_corners` picks the chart corners out of `cell_boundary`: the two agree
# point for point where the boundary is not densified, and at every stride-th
# vertex where it is.

module CellCornersTests

using Test
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: systems, levelgrid, ncells, cellindex, cell_boundary,
    cell_corners, LevelIndex

include(joinpath(@__DIR__, "..", "..", "helpers.jl"))
using .DGGTestHelpers: syslabel

const LEVEL = 2

function sample_cells(grid, n::Int)
    total = ncells(grid)
    step = max(1, total ÷ n)
    return [cellindex(grid, i) for i in 1:step:total]
end

# The corner cells of every HEALPix base pixel at `LEVEL`.
function healpix_base_corner_cells()
    nside = 2^LEVEL
    return [LevelIndex(LEVEL, DGG.HEALPix.xyf_to_nested(ix, iy, face, nside))
            for face in 0:11 for ix in (0, nside - 1) for iy in (0, nside - 1)]
end

igeo7_pentagons(grid) =
    [c for c in (cellindex(grid, i) for i in 1:ncells(grid)) if DGG.IGeo7.is_pentagon(c)]

function check_corners(grid, cells)
    for c in cells
        boundary = collect(cell_boundary(grid, c))
        corners = collect(cell_corners(grid, c))
        @test length(boundary) % length(corners) == 0
        stride = length(boundary) ÷ length(corners)
        # A5's densified ring starts one edge in, so the first corner is the
        # stride-th vertex there and the first everywhere else.
        first = findfirst(==(corners[1]), boundary)
        @test first !== nothing
        first === nothing && continue
        @test first in (1, stride)
        @test corners == boundary[first:stride:end]
    end
end

@testset "cell_corners agrees with cell_boundary" begin
    for sys in systems()
        @testset "$(syslabel(sys))" begin
            grid = levelgrid(sys, LEVEL)
            cells = sample_cells(grid, 64)
            sys isa DGG.HEALPixSystem && append!(cells, healpix_base_corner_cells())
            if sys isa DGG.IGeo7System
                pentagons = igeo7_pentagons(grid)
                @test length(pentagons) == 12
                append!(cells, pentagons)
            end
            check_corners(grid, cells)
        end
    end

    @testset "CopernicusDEMSystem" begin
        sys = DGG.CopernicusDEMSystem(90)
        grid = levelgrid(sys, first(DGG.levels(sys)))
        check_corners(grid, sample_cells(grid, 64))
    end

    @testset "corner cells stay on the stack" begin
        for sys in (DGG.HEALPixSystem(), DGG.ISEA4RSystem())
            grid = levelgrid(sys, LEVEL)
            c = cellindex(grid, 1)
            cell_corners(grid, c)
            @test @allocated(cell_corners(grid, c)) == 0
        end
    end
end

end # module
