# Baseline driver for the Geomorphometry kernels over the two subset shapes.
#
#   julia --project=<bench-env> benchmark/geomorphometry.jl
#
# The env holds this package plus `Geomorphometry` and `Rasters`:
#
#   pkg> dev <path-to-DiscreteGlobalGrids.jl>
#   pkg> add Geomorphometry#clipped-neighbors   # https://github.com/asinghvi17/GeoArrayOps.jl.git
#   pkg> add Rasters
#
# (The baseline below predates the fork rev and was measured against
# `Geomorphometry#feat/generic` from Deltares/Geomorphometry.jl.)
#
# No downloads, no plotting: the elevation field is synthetic (the same
# dome-and-ripple `test/integration/geomorphometry_synthetic.jl` uses), and the
# two shapes are the ones every stencil claim in this package is measured on —
# the 343-cell rooted subtree (one position window; `cellposition` is a
# subtraction) and the level-12 multi-order coverage of a 1°×1° Alpine tile
# (2,313,802 cells over 3,715 disjoint windows; `cellposition` is a binary
# search over them). Timings are `@timed` after one same-shape warmup run.
#
# Output on unmodified main (980e34d), Julia 1.12.6, M-series macOS, 2026-08-16:
#
#   one rooted subtree: 343 cells, 1 window(s)
#     topographic_position_index    0.0001 s   (0.0 s gc, 0.0 MB)
#     flowaccumulation(D8)          0.0004 s   (0.0 s gc, 0.0 MB)
#     halo_table                    0.0001 s   (0.0 s gc, 0.1 MB)
#   multi-order coverage: 2313802 cells, 3715 window(s)
#     topographic_position_index    1.3839 s   (0.0 s gc, 15.0 MB)
#     flowaccumulation(D8)          4.4336 s   (0.01 s gc, 112.5 MB)
#     halo_table                    0.9216 s   (0.11 s gc, 405.8 MB)
#
# Output on the branch pair — this tree paired with the fork's
# `clipped-neighbors` @ 9a4e053 — Julia 1.12.6, 8 threads, same machine,
# 2026-08-17. TPI rides the threaded `mapneighbors` sweep; D8 settles in
# position space over one `HaloTable`; `halo_table` is the same sweep
# materialized:
#
#   one rooted subtree: 343 cells, 1 window(s)
#     topographic_position_index    0.0001 s   (0.0 s gc, 0.0 MB)
#     flowaccumulation(D8)          0.0001 s   (0.0 s gc, 0.0 MB)
#     halo_table                    0.0001 s   (0.0 s gc, 0.1 MB)
#   multi-order coverage: 2313802 cells, 3715 window(s)
#     topographic_position_index    0.057  s   (0.0 s gc, 8.8 MB)
#     flowaccumulation(D8)          1.3891 s   (0.022 s gc, 192.5 MB)
#     halo_table                    0.5016 s   (0.022 s gc, 264.7 MB)

import DiscreteGlobalGrids as DGG
import Geomorphometry as GM
using Rasters

const Extents = Rasters.Extents          # a Rasters dep, so no extra import

sys = DGG.IGeo7System()

tile = Extents.Extent(X=(10.0, 11.0), Y=(46.0, 47.0))
root = DGG.coarsest_contained(DGG.query(sys, DGG.MultiOrderCoverage(tile); level=10))

# A single dome with a small ripple on it, sampled at cell centroids — one real
# drainage tree instead of a plateau of ties, with the tie-break kept off the
# degenerate path. Copied from `test/integration/geomorphometry_synthetic.jl`.
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

nwindows(cells) = DGG.Fallbacks.nwindows(DGG.Fallbacks.windows(cells))

function bench(label, cells)
    dem = make_dem(cells)
    println("$label: $(length(cells)) cells, $(nwindows(cells)) window(s)")
    for (name, f) in (
        ("topographic_position_index", () -> GM.topographic_position_index(dem)),
        ("flowaccumulation(D8)", () -> GM.flowaccumulation(dem; method=GM.D8())),
        ("halo_table", () -> DGG.halo_table(cells)),
    )
        f()                                       # warmup, same shape and types
        t = @timed f()
        println("  ", rpad(name, 30),
            round(t.time; digits=4), " s   (",
            round(t.gctime; digits=3), " s gc, ",
            round(t.bytes / 2^20; digits=1), " MB)")
    end
    return nothing
end

bench("one rooted subtree",
    DGG.CellVector(DGG.PartialGrid(sys, root, DGG.level(root) + 3)))
bench("multi-order coverage",
    DGG.CellVector(DGG.query(sys, DGG.MultiOrderCoverage(tile); level=12)))
