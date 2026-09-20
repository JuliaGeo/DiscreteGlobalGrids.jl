module NeighborhoodHookTests

using Test
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
import SmallCollections: SmallVector

struct HookGrid{G,W} <: AbstractGrid
    base::G
    order::W
    calls::Base.RefValue{Int}
    boundaries::Base.RefValue{Int}
end
HookGrid(base, order) = HookGrid(base, order, Ref(0), Ref(0))
DGG.ncells(g::HookGrid) = ncells(g.base)
DGG.cellindex(g::HookGrid, i::Int) = cellindex(g.base, i)
DGG.localindex(g::HookGrid, c::AbstractCellIndex) = DGG.localindex(g.base, c)
DGG.cell_centroid(g::HookGrid, c::AbstractCellIndex) = cell_centroid(g.base, c)
function DGG.cell_boundary(g::HookGrid, c::AbstractCellIndex)
    g.boundaries[] += 1
    return cell_boundary(g.base, c)
end
DGG.winding(g::HookGrid, ::Connectivity) = g.order
DGG.maxring(g::HookGrid, k::Integer, conn::Connectivity) = maxring(g.base, k, conn)
function DGG.one_ring(g::HookGrid, c::AbstractCellIndex, conn::Connectivity)
    g.calls[] += 1
    out = DGG.one_ring(g.base, c, conn)
    g.order isa Clockwise && length(out) > 1 && (out = reverse(out, 2, length(out)))
    return out
end

struct GeometricGrid{G} <: AbstractGrid
    base::G
    trees::Base.RefValue{Int}
end
DGG.ncells(g::GeometricGrid) = ncells(g.base)
DGG.cellindex(g::GeometricGrid, i::Int) = cellindex(g.base, i)
DGG.cell_centroid(g::GeometricGrid, c::AbstractCellIndex) = cell_centroid(g.base, c)
DGG.cell_boundary(g::GeometricGrid, c::AbstractCellIndex) = cell_boundary(g.base, c)
function DGG.treeify(g::GeometricGrid)
    g.trees[] += 1
    return invoke(DGG.treeify, Tuple{AbstractGrid}, g)
end

@testset "one_ring supplies every generic neighborhood" begin
    base = levelgrid(HEALPixSystem(), 3)
    for order in (CounterClockwise(), Clockwise(), Unordered()), conn in (Vertex(), Edge())
        grid = HookGrid(base, order)
        c = cellindex(base, 1)
        for k in 1:3
            @test neighbors(grid, c, k; connectivity = conn) == neighbors(base, c, k; connectivity = conn)
            @test ring(grid, c, k; connectivity = conn) == ring(base, c, k; connectivity = conn)
            @test ring(grid, c, Val(k); connectivity = conn) == ring(base, c, k; connectivity = conn)
            @test neighbors(grid, c, Val(k); connectivity = conn) == neighbors(base, c, k; connectivity = conn)
        end
        @test grid.calls[] > 0
        @test grid.boundaries[] == 0
        @test neighbors(grid, c, 3; connectivity = conn) ==
              reduce(vcat, [ring(grid, c, k; connectivity = conn) for k in 1:3]; init = LevelIndex[])
        @test DGG.shell_ring(grid, c, Val(3), conn) == ring(base, c, 3; connectivity = conn)
        @test DGG.shell_disc(grid, c, Val(3), conn) == neighbors(base, c, 3; connectivity = conn)
        @test DGG.Fallbacks.adjacency_shells(grid, c, 3, conn) ==
              [ring(base, c, k; connectivity = conn) for k in 1:3]
    end
end

function one_ring_bytes(verb::F, grid, cell) where {F}
    verb(grid, cell, 1)
    return @allocated verb(grid, cell, 1)
end

@testset "clockwise one-rings preserve fixed-capacity storage" begin
    grid = HookGrid(levelgrid(HEALPixSystem(), 3), Clockwise())
    baseline = HookGrid(grid.base, CounterClockwise())
    c = cellindex(grid, 1)
    for verb in (neighbors, ring)
        ids = verb(grid, c, 1)
        indices = verb(grid, 1, 1)
        @test ids isa SmallVector{8,LevelIndex}
        @test indices isa SmallVector{8,Int}
        @test cellindex.(grid, indices) == ids
        @test one_ring_bytes(verb, grid, c) <= one_ring_bytes(verb, baseline, c)
        @test one_ring_bytes(verb, grid, 1) <= one_ring_bytes(verb, baseline, 1)
    end
end

@testset "geometric neighborhood walks share one tree" begin
    grid = GeometricGrid(levelgrid(HEALPixSystem(), 1), Ref(0))
    c = cellindex(grid, 1)
    for verb in (neighbors, ring)
        grid.trees[] = 0
        @test !isempty(verb(grid, c, 2))
        @test grid.trees[] == 1
    end
end

end
