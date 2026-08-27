isdefined(@__MODULE__, :DaggerTalkFigures) ||
    include(joinpath(@__DIR__, "11-dagger-theme.jl"))

using .DGGSTalkFigures, .DaggerTalkFigures
using GLMakie, Makie

const COMMON_SHAPE_OUTPUT = joinpath(@__DIR__, "video", "17-common-shape.mp4")

smoothpulse(t, a, b) = begin
    x = clamp((t - a) / (b - a), 0, 1)
    x^2 * (3 - 2x)
end

lerp_point(a, b, t) = Point2f((1 - t) * a[1] + t * b[1],
    (1 - t) * a[2] + t * b[2])

function render_common_shape(path = COMMON_SHAPE_OUTPUT)
    GLMakie.activate!()
    times = collect(range(0, 6.5; length = 156))

    with_theme(JG_THEME) do
        fig, body = slide_figure()
        axis = clean_axis(body[1, 1]; limits = (0, 1, 0, 1))

        domain_alpha = Observable(0.20)
        runtime_alpha = Observable(0.12)
        packet_a = Observable(Point2f(0.25, 0.70))
        packet_b = Observable(Point2f(0.25, 0.30))
        packet_out = [Observable(Point2f(0.63, 0.50)) for _ in 1:4]
        packet_visibility = Observable(0.0)
        worker_fill = [Observable(JG.paper_off) for _ in 1:4]

        # Connectors first, so they stay behind the entities.
        lines!(axis, [Point2f(0.24, 0.70), Point2f(0.43, 0.55)];
            color = @lift((JG.green_dark, $domain_alpha)), linewidth = 2.2)
        lines!(axis, [Point2f(0.24, 0.30), Point2f(0.43, 0.45)];
            color = @lift((JG.green_dark, $domain_alpha)), linewidth = 2.2)
        for (i, y) in enumerate((0.76, 0.58, 0.40, 0.22))
            lines!(axis, [Point2f(0.64, 0.50), Point2f(0.78, y)];
                color = @lift((JG.purple, $runtime_alpha)), linewidth = 2.0)
        end

        draw_box!(axis, 0.05, 0.62, 0.19, 0.16, "REGRID\nsource chunks → output chunk";
            fill = JG.green_50, stroke = JG.green_dark, fontsize = 11)
        draw_box!(axis, 0.05, 0.22, 0.19, 0.16, "STENCIL\nhalo chunks → owned chunk";
            fill = JG.green_50, stroke = JG.green_dark, fontsize = 11)
        draw_box!(axis, 0.43, 0.38, 0.21, 0.24,
            "TASK\nowned output\n+ required inputs\n+ operation";
            fill = JG.purple_100, stroke = JG.purple, fontsize = 12)
        for (i, y) in enumerate((0.69, 0.51, 0.33, 0.15))
            poly!(axis, Rect2f(0.78, y, 0.16, 0.14); color = worker_fill[i],
                strokecolor = JG.purple, strokewidth = 1.6)
            text!(axis, 0.86, y + 0.07; text = "WORKER $(i)",
                font = FONT_BOLD, fontsize = 11, color = JG.purple,
                align = (:center, :center))
        end

        scatter!(axis, packet_a; color = JG.green_dark, markersize = 13)
        scatter!(axis, packet_b; color = JG.green_dark, markersize = 13)
        for packet in packet_out
            scatter!(axis, packet; color = @lift((JG.purple, $packet_visibility)),
                markersize = 12)
        end
        mkpath(dirname(path))
        record(fig, path, times; framerate = 24, backend = GLMakie,
            compression = 18, profile = "high", pixel_format = "yuv420p",
            px_per_unit = 1, visible = false, loglevel = "warning") do t
            enter = smoothpulse(t, 0.5, 2.3)
            packet_a[] = lerp_point((0.24, 0.70), (0.43, 0.55), enter)
            packet_b[] = lerp_point((0.24, 0.30), (0.43, 0.45), enter)
            domain_alpha[] = 0.20 + 0.80enter

            fan = smoothpulse(t, 2.5, 4.6)
            runtime_alpha[] = 0.12 + 0.88fan
            packet_visibility[] = fan
            for (i, y) in enumerate((0.76, 0.58, 0.40, 0.22))
                lag = smoothpulse(t, 2.6 + 0.18(i - 1), 4.0 + 0.18(i - 1))
                packet_out[i][] = lerp_point((0.64, 0.50), (0.78, y), lag)
                worker_fill[i][] = lag > 0.82 ? JG.purple_100 : JG.paper_off
            end
        end
    end
    println("wrote ", path)
    return path
end

abspath(PROGRAM_FILE) == (@__FILE__) && render_common_shape()
