# Benchmark Geomorphometry kernels over rooted and multi-window subsets.
#
#   julia --project=<bench-env> benchmark/geomorphometry.jl
#
# The environment needs this package, `Geomorphometry`, and `Rasters`:
#
#   pkg> dev <path-to-DiscreteGlobalGrids.jl>
#   pkg> add Geomorphometry#clipped-neighbors   # https://github.com/asinghvi17/GeoArrayOps.jl.git
#   pkg> add Rasters
#
# The synthetic elevation field matches the integration test. Each timing uses
# `@timed` after a warmup with the same subset shape.
#
# Output on unmodified main (980e34d), Julia 1.12.6, M-series macOS, 2026-08-16;
# this baseline predates the fork rev and was measured against
# `Geomorphometry#feat/generic` from Deltares/Geomorphometry.jl:
#
#   one rooted subtree: 343 cells, 1 window(s)
#     topographic_position_index    0.0001 s   (0.0 s gc, 0.0 MB)
#     flowaccumulation(D8)          0.0004 s   (0.0 s gc, 0.0 MB)
#     adjacency                     0.0001 s   (0.0 s gc, 0.1 MB)
#   multi-order coverage: 2313802 cells, 3715 window(s)
#     topographic_position_index    1.3839 s   (0.0 s gc, 15.0 MB)
#     flowaccumulation(D8)          4.4336 s   (0.01 s gc, 112.5 MB)
#     adjacency                     0.9216 s   (0.11 s gc, 405.8 MB)
#
# Output with `clipped-neighbors` @ 9a4e053, Julia 1.12.6, 8 threads,
# M-series macOS, 2026-08-17. TPI uses `mapneighbors`; D8 uses an
# `AdjacencyTable`. Both table builds use contiguous chunks and match their
# sequential arrays:
#
#   one rooted subtree: 343 cells, 1 window(s)
#     topographic_position_index    0.0001 s   (0.0 s gc, 0.0 MB)
#     flowaccumulation(D8)          0.0002 s   (0.0 s gc, 0.1 MB)
#     adjacency                     0.0001 s   (0.0 s gc, 0.1 MB)
#   multi-order coverage: 2313802 cells, 3715 window(s)
#     topographic_position_index    0.0653 s   (0.0 s gc, 8.8 MB)
#     flowaccumulation(D8)          1.1874 s   (0.041 s gc, 307.2 MB)
#     adjacency                     0.1894 s   (0.104 s gc, 264.7 MB)

import DiscreteGlobalGrids as DGG
import Geomorphometry as GM
using Rasters

const Extents = Rasters.Extents          # Available through Rasters.

sys = DGG.IGeo7System()

tile = Extents.Extent(X=(10.0, 11.0), Y=(46.0, 47.0))
root = DGG.coarsest_contained(DGG.query(sys, DGG.MultiOrderCoverage(tile); level=10))

# Sample a dome and small ripple at cell centroids: the dome gives the routers
# one real drainage tree, the ripple keeps the D8 tie-break off the
# all-neighbours-equidistant path.
function make_dem(cells)
    complete = DGG.levelgrid(sys, DGG.level(cells))
    apex = DGG.cell_centroid(DGG.levelgrid(sys, DGG.level(root)), root)
    elevation = [
        begin
            p = DGG.cell_centroid(complete, c)
            d = acos(clamp(p[1] * apex[1] + p[2] * apex[2] + p[3] * apex[3], -1, 1))
            1000.0 * exp(-(d / 0.02)^2) + 40.0 * sin(120 * p[1]) * cos(97 * p[2])
        end
        for c in cells
    ]
    return Raster(elevation, (DGG.Cells(DGG.CellLookup(cells)),); name=:height)
end

nwindows(cells) = DGG.Engine.nwindows(DGG.Engine.windows(cells))

function bench(label, cells)
    dem = make_dem(cells)
    println("$label: $(length(cells)) cells, $(nwindows(cells)) window(s)")
    for (name, f) in (
        ("topographic_position_index", () -> GM.topographic_position_index(dem)),
        ("flowaccumulation(D8)", () -> GM.flowaccumulation(dem; method=GM.D8())),
        ("adjacency", () -> DGG.adjacency(cells)),
    )
        f()                                       # Warm up this shape and type.
        t = @timed f()
        println("  ", rpad(name, 30),
            round(t.time; digits=4), " s   (",
            round(t.gctime; digits=3), " s gc, ",
            round(t.bytes / 2^20; digits=1), " MB)")
    end
    return nothing
end

bench("one rooted subtree",
    DGG.CellVector(DGG.subtree(sys, root, DGG.level(root) + 3)))
bench("multi-order coverage",
    DGG.CellVector(DGG.query(sys, DGG.MultiOrderCoverage(tile); level=12)))
