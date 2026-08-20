# Compare `dggpoly` against Makie's `poly` on IGEO7 cell sets of growing size.
#
#     xvfb-run -a julia --project=lib/DiscreteGlobalGridVisualizations/bench -t8 \
#         lib/DiscreteGlobalGridVisualizations/bench/bench.jl gl --maxlevel=12
#
# Pass `cairo` (the default) or `gl` to pick a backend, and `--maxlevel=N` to
# stop earlier.  Two numbers are reported for each size:
#
#   * **convert** — everything a plot does to its argument before a backend sees
#     it: for `poly` that is one triangulation per cell plus the outline it
#     builds whether or not it is drawn; for `dggpoly` it is `tessellate`.  This
#     is the part this package replaces, and it does not depend on the backend.
#   * **plot** — creating the figure and saving it, which adds rasterisation.
#     Under `xvfb` GLMakie rasterises in software, so treat the GPU-backend
#     numbers as an upper bound on that half.

import DiscreteGlobalGrids as DGG
import DiscreteGlobalGridVisualizations as DGGV
using DiscreteGlobalGridVisualizations
using Printf

const BACKEND = any(==("gl"), ARGS) ? :gl : :cairo
const MAXLEVEL = let i = findfirst(a -> startswith(a, "--maxlevel="), ARGS)
    isnothing(i) ? 11 : parse(Int, split(ARGS[i], "=")[2])
end

if BACKEND === :gl
    using GLMakie
    GLMakie.activate!()
else
    using CairoMakie
    CairoMakie.activate!(type = "png")
end

const SYS = DGG.IGeo7System()
const OUT = mktempdir()

"A patch of IGEO7 cells at `level` covering `extent`."
patch(level, extent) = DGG.CellVector(DGG.query(SYS, DGG.MultiOrderCoverage(extent); level))

"Run `f` once to compile, then time the second run."
timeit(f) = (f(); GC.gc(); @elapsed f())

function report(name, cells)
    n = length(cells)
    values = Float32.(1:n)

    # Conversion only: resolve exactly the attributes that hold the geometry.
    convert_poly = timeit() do
        plot = poly!(Scene(), cells; color = values, strokewidth = 0)
        plot.meshes[], plot.outline[], plot.computed_strokecolor[]
    end
    convert_dgg = timeit() do
        plot = dggpoly!(Scene(), cells; color = values, strokewidth = 0, primitive = :mesh)
        plot.mesh_positions[], plot.mesh_faces[], plot.mesh_color[]
    end

    # Figure to file.
    plot_poly = timeit() do
        figure, _, _ = poly(cells; color = values, strokewidth = 0)
        save(joinpath(OUT, "poly.png"), figure)
    end
    plot_dgg = timeit() do
        figure, _, _ = dggpoly(cells; color = values, strokewidth = 0)
        save(joinpath(OUT, "dgg.png"), figure)
    end

    @printf(
        "%-15s n=%8d   convert %8.3f -> %7.3f s (%5.1fx)   plot %8.3f -> %7.3f s (%5.1fx)\n",
        name, n, convert_poly, convert_dgg, convert_poly / convert_dgg,
        plot_poly, plot_dgg, plot_poly / plot_dgg
    )
    return nothing
end

alps = DGG.Extents.Extent(X = (10.0, 11.0), Y = (46.0, 47.0))

println("backend: ", BACKEND, "   threads: ", Threads.nthreads())
for level in 7:MAXLEVEL
    report("IGEO7 level $level", patch(level, alps))
end
