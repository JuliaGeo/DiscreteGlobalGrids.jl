isdefined(@__MODULE__, :DGGSTalkFigures) ||
    include(joinpath(@__DIR__, "00-dggs-theme.jl"))

module DaggerTalkFigures

using ..DGGSTalkFigures
using DiscreteGlobalGrids
using DiscreteGlobalGridsVisualization
using FlyThroughPaths
using Makie
import GeoInterface as GI
import GeometryOps as GO
import GlobalRegridding as GR
import DimensionalData as DD

export FPS, FIELD_COLORMAP, FIELD_RANGE
export synthetic_raster, synthetic_cube, regrid_fixture
export html_path, video_path, clean_axis, draw_box!, draw_arrow!
export camera_orbit, camera_times, record_video

const FPS = 30
const FIELD_COLORMAP = [JG.purple, JG.purple_100, JG.paper,
    JG.green_100, colorant"#7dda71"]
const FIELD_RANGE = (-900.0, 1100.0)

html_path(name) = joinpath(@__DIR__, "html", name)
video_path(name) = joinpath(@__DIR__, "video", name)

function sampled_axis(D, values, step)
    D(DD.Sampled(collect(values); span = DD.Regular(step),
        sampling = DD.Intervals(DD.Center()), order = DD.ForwardOrdered()))
end

function synthetic_raster(; step = 15.0)
    lon = collect((-180 + step / 2):step:(180 - step / 2))
    lat = collect((-90 + step / 2):step:(90 - step / 2))
    xdim = sampled_axis(DD.X, lon, step)
    ydim = sampled_axis(DD.Y, lat, step)
    values = [1000sind(3x) * cosd(2y) + 100 for x in lon, y in lat]
    data = DD.DimArray(values, (xdim, ydim); name = :signal)
    return (; data, lon, lat, xdim, ydim)
end

function synthetic_cube(; nt = 24, step = 15.0)
    src = synthetic_raster(; step)
    values = [src.data[x, y] + 120sinpi(2(t - 1) / nt)
        for x in eachindex(src.lon), y in eachindex(src.lat), t in 1:nt]
    time = DD.Dim{:time}(1:nt)
    data = DD.DimArray(values, (src.xdim, src.ydim, time); name = :signal)
    return (; src..., data, time)
end

function regrid_fixture(; nt = 24)
    src = synthetic_raster()
    cube = synthetic_cube(; nt)
    grid = levelgrid(IGeo7System(), 2)
    plan = plan_regrid(src.data; to = grid, missingpolicy = Weighted(0.01))
    result = regrid(src.data, plan)
    cubeplan = plan_regrid(cube.data; to = grid, missingpolicy = Weighted(0.01))
    cuberesult = regrid(cube.data, cubeplan)
    return (; src..., cube = cube.data, grid, plan, result, cubeplan, cuberesult)
end

function clean_axis(slot; limits = (0, 1, 0, 1), background = JG.paper)
    axis = Axis(slot; backgroundcolor = background)
    hidedecorations!(axis)
    hidespines!(axis)
    xlims!(axis, limits[1], limits[2])
    ylims!(axis, limits[3], limits[4])
    return axis
end

function draw_box!(axis, x, y, w, h, label;
        fill = JG.paper, stroke = JG.ink, textcolor = JG.ink,
        linewidth = 1.5, fontsize = 12)
    poly!(axis, Rect2f(x, y, w, h); color = fill,
        strokecolor = stroke, strokewidth = linewidth)
    text!(axis, x + w / 2, y + h / 2; text = label,
        align = (:center, :center), font = DGGSTalkFigures.FONT_BODY,
        fontsize, color = textcolor)
    return axis
end

function draw_arrow!(axis, a, b; color = JG.ink, linewidth = 1.5)
    lines!(axis, Point2f[a, b]; color, linewidth)
    direction = Point2f(b[1] - a[1], b[2] - a[2])
    scatter!(axis, [Point2f(b)]; marker = :rtriangle, markersize = 11,
        rotation = atan(direction[2], direction[1]), color)
    return axis
end

function rotate_z(v, theta)
    c, s = cos(theta), sin(theta)
    return [c * v[1] - s * v[2], s * v[1] + c * v[2], v[3]]
end

function camera_orbit(axis; seconds = 6.0, angle = pi / 5)
    start = capture_view(axis)
    stop = ViewState(; eyeposition = rotate_z(start.eyeposition, angle),
        lookat = start.lookat, upvector = start.upvector, fov = start.fov)
    return Path(start) * Pause(0.5) *
        ConstrainedMove(seconds - 1.0, stop; constraint = :rotation,
            speed = :sinusoidal) * Pause(0.5)
end

function camera_times(path; fps = FPS)
    n = max(2, round(Int, FlyThroughPaths.duration(path) * fps))
    range(0, FlyThroughPaths.duration(path); length = n)
end

function record_video(path, fig, frames, update!; fps = FPS, backend = nothing)
    mkpath(dirname(path))
    kwargs = (; framerate = fps, compression = 18, profile = "high",
        pixel_format = "yuv420p", px_per_unit = 1, visible = false,
        loglevel = "warning")
    if backend === nothing
        record(fig, path, frames; kwargs...) do frame
            update!(frame)
        end
    else
        record(fig, path, frames; backend, kwargs...) do frame
            update!(frame)
        end
    end
    println("wrote ", path)
    return path
end

end
