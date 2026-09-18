# An interactive plot over a pyramid store: zoom in, and a finer level loads.
#
#   julia --project=docs -i scripts/pyramid_stream_plot.jl
#
# Only Makie core is used here, so the backend is the caller's: running the file
# picks GLMakie, and `using WGLMakie; include(...); streamplot()` works the same
# way. The descent runs where the store is, and only the chosen cells cross to
# the frontend.
#
# Everything about which cells to draw is `pyramid_stream_poc.jl`'s; this file
# is the camera. On every settled camera change the axis limits become a view,
# the descent answers with a level and its cells, and the plot swaps to them.
# The terminal traces what that cost, level by level, in CHUNK reads.
#
# The store is built from the hydrology tutorial's tile the first time and
# cached in the temporary directory afterwards. Set `DGG_DEMO_STORE` to reuse
# one somewhere else.
#
ENV["RASTERDATASOURCES_PATH"] =
    mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import Extents
using Rasters, RasterDataSources
import ArchGDAL
using Zarr
using Makie
using DiscreteGlobalGridsVisualization: dggpoly!, cellset
using Printf

include(joinpath(@__DIR__, "pyramid_stream_poc.jl"))
using .PyramidStream

const SYS = DGG.IGeo7System()

"""
    demostore(; centre = (10.5, 46.5)) -> String

The hydrology tutorial's Copernicus DEM tile, regridded onto IGEO7 and written
as a pyramid. Built once and cached.
"""
function demostore(; centre=(10.5, 46.5))
    path = get(ENV, "DGG_DEMO_STORE", joinpath(tempdir(), "ortler_pyramid.zarr"))
    isdir(path) && return path
    @info "building the demo store (one 1x1 degree Copernicus DEM tile)" path
    tile = only(skipmissing(RasterDataSources.getraster(CopernicusDEM;
        extent=GI.extent(centre))))
    Sys.isapple() && Rasters.checkmem!(false)
    dem = Raster(tile; lazy=false)
    level = DGG.levelfor(SYS, dem)
    region = DGG.query(SYS, DGG.MultiOrderCoverage(Rasters.extent(dem)); level=level)
    elevation = rebuild(DGG.regrid(dem; to=region); name=:elevation)
    DGG.dggwrite(path, elevation; layout=:pyramid)
    return path
end

"The longitude/latitude box a Makie axis is currently looking at."
function viewextent(limits)
    o, w = limits.origin, limits.widths
    return Extents.Extent(X=(o[1], o[1] + w[1]), Y=(o[2], o[2] + w[2]))
end

# A camera moves continuously and the descent is not free, so a frame is only
# rebuilt when the view has really changed: a fifth of a span in either
# direction. Without this a single drag would re-read the store fifty times.
function moved(old, new)
    old === nothing && return true
    ow, nw = old.X[2] - old.X[1], new.X[2] - new.X[1]
    (nw < ow / 1.2 || nw > ow * 1.2) && return true
    shift = max(abs(sum(new.X) - sum(old.X)), abs(sum(new.Y) - sum(old.Y))) / 2
    return shift > nw / 5
end

"""
    streamplot(store; maxcells, colorrange) -> Figure

Draw a pyramid store, reloading a finer level whenever the camera settles on a
smaller view.
"""
function streamplot(store=demostore(); maxcells::Integer=200_000,
    colorrange=(200, 3800), colormap=:terrain,
    limits=((-180.0, 180.0), (-90.0, 90.0)))

    pyr = DGG.dggread(store)[:elevation]
    @info "opened" store levels = DGG.levels(pyr)

    # The whole earth to begin with: the store holds one tile, so the opening
    # frame is a handful of coarse cells over the Alps and nothing else, which
    # is the layout working rather than a bug. Scroll in and the levels load.
    fig = Figure(size=(1000, 620))
    ax = Axis(fig[1, 1]; xlabel="longitude", ylabel="latitude",
        aspect=DataAspect(), limits=limits)

    # How wide the axis really is, in pixels: the descent stops when a cell is
    # about three of them across, so it has to be asked rather than assumed.
    axpixels() = round(Int, widths(viewport(ax.scene)[])[1])

    e0 = viewextent(ax.finallimits[])
    cells0, level0, values0, trace0 = PyramidStream.streamframe(pyr, e0;
        pixels=axpixels(), maxcells)
    # `dggpoly` draws a cell SET: the cells plus what to ask for their
    # boundaries, which for a bare vector of ids is the system.
    plot = dggpoly!(ax, cellset(SYS, cells0); color=Float32.(values0),
        colormap=colormap, colorrange=colorrange, strokewidth=0)
    Colorbar(fig[1, 2], plot; label="elevation (m)")
    ax.title = _titlefor(level0, cells0, trace0, 0.0)

    shown = Ref{Any}(e0)
    on(ax.finallimits) do limits
        e = viewextent(limits)
        moved(shown[], e) || return nothing
        shown[] = e
        @printf("\nview %.4f deg across at (%.3f, %.3f)\n",
            e.X[2] - e.X[1], sum(e.X) / 2, sum(e.Y) / 2)
        t0 = time()
        new, level, values, trace = PyramidStream.streamframe(pyr, e;
            pixels=axpixels(), maxcells)
        isempty(new) && return nothing
        elapsed = time() - t0
        # Cells and colours in ONE update. Setting them as two observables lets
        # the compute pipeline rebuild the mesh before it has the matching
        # colours, and it refuses -- rightly -- to draw 9159 cells in 6944
        # colours.
        Makie.update!(plot, cellset(SYS, new); color=Float32.(values))
        ax.title = _titlefor(level, new, trace, elapsed)
        println("  => ", ax.title[])
        return nothing
    end

    return fig
end

_titlefor(level, cells, trace, elapsed) = @sprintf(
    "level %d, %d cells, %d chunk reads, %.0f ms",
    level, length(cells), sum(t.chunks for t in trace), elapsed * 1000)

if abspath(PROGRAM_FILE) == (@__FILE__) || isinteractive()
    # The backend is chosen here and nowhere else; swap this line for WGLMakie
    # and the rest of the file is unchanged.
    @eval using GLMakie
    fig = streamplot()
    display(fig)
    println("\nZoom with the scroll wheel; the terminal traces every reload.")
end
