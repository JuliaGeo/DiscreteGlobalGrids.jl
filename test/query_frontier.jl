module QueryFrontierTests
using Test
import DiscreteGlobalGrids as DGG
import DiscreteGlobalGrids.DE9IM as DE9IM
import GeoInterface as GI
import GeometryOps as GO

# GeoInterface's wrapper refuses a one-point line; this one stands for it.
struct OnePointLine end
GI.geomtrait(::OnePointLine) = GI.LineStringTrait()
GI.ncoord(::GI.LineStringTrait, ::OnePointLine) = 2
GI.npoint(::GI.LineStringTrait, ::OnePointLine) = 1
GI.getpoint(::GI.LineStringTrait, ::OnePointLine) = ((12.0, 20.0),)
GI.getpoint(::GI.LineStringTrait, ::OnePointLine, i) = (12.0, 20.0)

box(x, y, r) = GI.LinearRing([(x-r, y-r), (x+r, y-r), (x+r, y+r), (x-r, y+r), (x-r, y-r)])

# Cells of mixed levels: a grid that reports a system and no level, which is the
# only shape the engine answers with an IndexTree.
struct MixedLevelGrid{S,V} <: DGG.AbstractGrid
    system::S
    ids::V
end
DGG.system(g::MixedLevelGrid) = g.system
DGG.ncells(g::MixedLevelGrid) = length(g.ids)
DGG.cellindex(g::MixedLevelGrid, i::Int) = g.ids[i]
for f in (:cell_boundary, :cell_centroid)
    @eval DGG.$f(g::MixedLevelGrid, c::DGG.AbstractCellIndex) =
        DGG.$f(DGG.levelgrid(g.system, DGG.level(c)), c)
end

# Brute force over every cell: the descent's bulk accept and prune must agree
# with the exact predicate asked of each cell polygon or centroid.
function oracle(grid, target, pred)
    ask(c) = pred isa DGG.Engine.CentroidCovered ?
        GO.relate_predicate(target.prepared, GO.pred_intersects(), DGG.cell_centroid(grid, c)) :
        GO.relate_predicate(target.prepared, DGG.Engine._converse_predicate(pred), DGG.cell_polygon(grid, c))
    return [i for i in 1:DGG.ncells(grid) if ask(DGG.cellindex(grid, i))]
end

@testset "query frontier" begin
    # Children overhang parents here, so a node's cap centre is not a witness.
    grid = DGG.levelgrid(DGG.H3System(), 2)
    holed = GI.Polygon([box(10.0, 20.0, 30.0), box(10.0, 20.0, 8.0)])
    target = DGG.Engine._query_target(holed)
    for pred in (DE9IM.Intersects(nothing), DE9IM.Within(nothing), DE9IM.CoveredBy(nothing),
            DGG.Engine.CentroidCovered())
        hits = DGG.Engine._query_indices(grid, pred, target)
        @test hits == oracle(grid, target, pred)
        @test !isempty(hits)
    end
    # Cells deep inside the box are boundary-free, and none of them contains it.
    @test isempty(DGG.Engine._query_indices(grid, DE9IM.Contains(nothing), target))
    @test DGG.Engine._query_indices(grid, DE9IM.Touches(nothing), target) ==
        oracle(grid, target, DE9IM.Touches(nothing))

    # A target wider than a hemisphere has no bounding cap; only the frontier
    # prunes.
    wide = GI.Polygon([box(0.0, 0.0, 80.0)])
    target = DGG.Engine._query_target(wide)
    @test target.cap.radius >= pi
    for pred in (DE9IM.Intersects(nothing), DE9IM.Within(nothing), DGG.Engine.CentroidCovered())
        @test DGG.Engine._query_indices(grid, pred, target) == oracle(grid, target, pred)
    end

    # A cap target never narrows the frontier; the centroid rule reads the cap.
    cap = GO.UnitSpherical.SphericalCap(DGG.Fallbacks.unit_point(10.0, 20.0), 0.3)
    target = DGG.Engine._query_target(cap)
    hits = DGG.Engine._query_indices(grid, DGG.Engine.CentroidCovered(), target)
    @test hits == [i for i in 1:DGG.ncells(grid)
        if GO.UnitSpherical.spherical_distance(DGG.cell_centroid(grid, DGG.cellindex(grid, i)),
            target.cap.point) <= 0.3]

    # Two parts with a gap between them: an edge invented from the last vertex
    # of one ring to the first of the next would accept cells in that gap.
    parts = GI.MultiPolygon([GI.Polygon([box(-40.0, 10.0, 6.0)]),
                             GI.Polygon([box(40.0, 10.0, 6.0)])])
    target = DGG.Engine._query_target(parts)
    for pred in (DE9IM.Intersects(nothing), DGG.Engine.CentroidCovered())
        hits = DGG.Engine._query_indices(grid, pred, target)
        @test hits == oracle(grid, target, pred)
        @test !isempty(hits)
    end

    # Copernicus tiles bring their own cursor, so the descent's node arms have to
    # answer for a tree that is neither the hierarchical cursor nor an IndexTree.
    copdem = DGG.levelgrid(DGG.CopernicusDEMSystem(90), 0)
    @test typeof(DGG.Engine._query_tree(copdem)) == typeof(DGG.treeify(copdem))
    target = DGG.Engine._query_target(GI.Polygon([box(10.0, 20.0, 10.0)]))
    for pred in (DE9IM.Intersects(nothing), DGG.Engine.CentroidCovered())
        hits = DGG.Engine._query_indices(copdem, pred, target)
        @test hits == oracle(copdem, target, pred)
        @test !isempty(hits)
    end

    # Stored cells of two levels are the IndexTree's case. Its nodes carry the
    # whole subtree window, which is what the bulk accept appends and the witness
    # reads; a window narrowed to the node's own leaf block would lose cells.
    centre = DGG.Fallbacks.unit_point(10.0, 20.0)
    near(g) = [DGG.cellindex(g, i) for i in 1:DGG.ncells(g)
               if GO.UnitSpherical.spherical_distance(
                   DGG.cell_centroid(g, DGG.cellindex(g, i)), centre) < deg2rad(10)]
    h3 = DGG.H3System()
    stored = MixedLevelGrid(h3, vcat(near(DGG.levelgrid(h3, 2)),
                                     near(DGG.levelgrid(h3, 3))))
    @test DGG.Engine._query_tree(stored) isa DGG.Engine.IndexTreeNode
    target = DGG.Engine._query_target(GI.Polygon([box(10.0, 20.0, 4.0)]))
    for pred in (DE9IM.Intersects(nothing), DE9IM.Within(nothing),
            DGG.Engine.CentroidCovered())
        hits = DGG.Engine._query_indices(stored, pred, target)
        @test hits == oracle(stored, target, pred)
        @test !isempty(hits)
    end

    # A one-point part is a boundary point of its own; the frontier keeps it.
    arcs = DGG.Engine._boundary_arcs(OnePointLine())
    @test length(arcs) == 1 && arcs[1].a == arcs[1].b
end

end
