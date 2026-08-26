# What one dual-cell point query costs on a conforming DGGS.
#
#     julia -t 1 --project=benchmark benchmark/dggs_dual_cells.jl
#
# Every arm is the same seam on the same source space in the same session, and
# every arm pays the same cell location, so what they compare is what the dual
# cell adds on top of the location a point method already pays for:
#
#   1. `nearest`      `weightsat!` through `NearestCell`'s sampler: one
#                     `cellat` and one entry. The point-location floor.
#   2. `barycentric`  `weightsat!` through `BarycentricPoint`'s sampler on a
#                     `DGGSpace`: the same location, then the host's two
#                     one-rings, the chart, the containing dual cell and its
#                     mean-value coordinates.
#
# Two more arms split that surcharge between the pieces it is made of:
#
#   3. `rings`        the location, then the host's `Vertex()` and `Edge()`
#                     one-rings: the topology a dual cell is cut out of.
#   4. `sites`        the location, then one `cell_centroid`: the sample site a
#                     dual cell reads once per node.
#
# Destinations are the cells of a finer grid of another system, so no query
# lands on a source sample site and every one runs the whole search.

using Printf
import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR

const LEVEL = 5
const DSTLEVEL = 4

function timequeries(smp, points, row)
    for p in points
        GR.weightsat!(row, smp, p)
    end
    t = @elapsed for p in points
        GR.weightsat!(row, smp, p)
    end
    b = @allocated for p in points
        GR.weightsat!(row, smp, p)
    end
    return t / length(points) * 1e9, b / length(points)
end

# The two one-rings a dual cell is cut out of, and nothing else.
function ringpass(grid, points)
    n = 0
    for p in points
        c = DGG.cellindex(grid, DGG.localindex(grid, p))
        n += length(DGG.neighbors(grid, c, 1; connectivity = DGG.Vertex()))
        n += length(DGG.neighbors(grid, c, 1; connectivity = DGG.Edge()))
    end
    return n
end

# One sample site, the read a dual cell repeats once per node.
function sitepass(grid, points)
    s = 0.0
    for p in points
        s += DGG.cell_centroid(grid, DGG.cellindex(grid, DGG.localindex(grid, p)))[1]
    end
    return s
end

function timepass(f, points)
    f(points)
    t = @elapsed f(points)
    b = @allocated f(points)
    return t / length(points) * 1e9, b / length(points)
end

function main()
    src = DGG.DGGSpace(DGG.levelgrid(DGG.IGeo7System(), LEVEL))
    dstgrid = DGG.levelgrid(DGG.S2System(), DSTLEVEL)
    points = [DGG.cell_centroid(dstgrid, DGG.cellindex(dstgrid, i))
              for i in 1:DGG.ncells(dstgrid)]
    row = GR.WeightRow()

    @printf("IGeo7 level %d source, %d cells; %d destination points\n",
        LEVEL, GR.ncells(src), length(points))
    nt, nb = timequeries(GR.sampler(GR.NearestCell(), src), points, row)
    bt, bb = timequeries(GR.sampler(GR.BarycentricPoint(), src), points, row)
    grid = DGG.levelgrid(DGG.IGeo7System(), LEVEL)
    rt, rb = timepass(pts -> ringpass(grid, pts), points)
    st, sb = timepass(pts -> sitepass(grid, pts), points)
    @printf("nearest      %8.1f ns/query  %6.1f bytes/query\n", nt, nb)
    @printf("barycentric  %8.1f ns/query  %6.1f bytes/query\n", bt, bb)
    @printf("rings        %8.1f ns/query  %6.1f bytes/query\n", rt, rb)
    @printf("sites        %8.1f ns/query  %6.1f bytes/query\n", st, sb)
    @printf("ratio        %8.2fx\n", bt / nt)
    return nothing
end

main()
