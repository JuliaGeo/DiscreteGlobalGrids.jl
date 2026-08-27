isdefined(@__MODULE__, :DaggerTalkFigures) ||
    include(joinpath(@__DIR__, "11-dagger-theme.jl"))

using .DGGSTalkFigures, .DaggerTalkFigures
using Makie

function boundary_options_app()
    with_theme(JG_THEME) do
        fig, body = slide_figure()
        axis = clean_axis(body[1, 1]; limits = (0, 1, 0, 1))

        # A single continuum, smallest contract on the left.
        lines!(axis, [Point2f(0.13, 0.24), Point2f(0.87, 0.24)];
            color = JG.hairline, linewidth = 3)
        scatter!(axis, [Point2f(0.13, 0.24), Point2f(0.87, 0.24)];
            color = JG.caption, markersize = 8)

        draw_box!(axis, 0.08, 0.45, 0.25, 0.24,
            "CHUNK GRAPH\n\nrequired input chunks\nfor each output chunk";
            fill = JG.green_100, stroke = JG.green_dark,
            textcolor = JG.ink, linewidth = 2.4, fontsize = 12)
        draw_box!(axis, 0.375, 0.45, 0.25, 0.24,
            "CHUNK + CELL MAP\n\nexact halo cells\nand supplying chunks";
            fill = JG.paper_off, stroke = JG.caption,
            textcolor = JG.ink, fontsize = 12)
        draw_box!(axis, 0.67, 0.45, 0.25, 0.24,
            "TASK PLAN\n\ndomain-shaped tasks\nand execution order";
            fill = JG.paper_off, stroke = JG.caption,
            textcolor = JG.ink, fontsize = 12)

        text!(axis, 0.13, 0.17; text = "less domain knowledge in Dagger",
            font = FONT_BODY, fontsize = 10, color = JG.caption,
            align = (:left, :center))
        text!(axis, 0.87, 0.17; text = "more domain knowledge in Dagger",
            font = FONT_BODY, fontsize = 10, color = JG.caption,
            align = (:right, :center))
        static_app(fig)
    end
end

export_boundary_options(dir = joinpath(@__DIR__, "html")) =
    export_html(joinpath(dir, "18-boundary-options.html"), boundary_options_app())

abspath(PROGRAM_FILE) == (@__FILE__) && export_boundary_options()
