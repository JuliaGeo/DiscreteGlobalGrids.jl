isdefined(@__MODULE__, :DGGSTalkFigures) ||
    include(joinpath(@__DIR__, "00-dggs-theme.jl"))

using .DGGSTalkFigures
using DiscreteGlobalGrids, DiscreteGlobalGridsVisualization, Makie

function igeo7_hierarchy_app()
    system = IGeo7System()

    with_theme(JG_THEME) do
        fig, body = slide_figure()
        globe = globe_axis(body[1, 1]; camera_longlat = (20, 18),
            camera_altitude = 1.55)

        # Draw fine-to-coarse so every parent boundary remains readable.  The
        # hierarchy is one geometry, not three small multiples: each darker
        # outline contains the successively fainter level below it.
        # Levels 1–3 are deliberately all visible: the fine mesh keeps enough
        # contrast to read as cells, while the coarser boundaries still anchor
        # the containment story.
        for (level, alpha, width) in ((3, 0.32, 0.48),
                (2, 0.60, 1.02), (1, 0.94, 2.05))
            dggpoly!(globe, levelgrid(system, level);
                color = (JG.paper, 0.0),
                strokecolor = (JG.green_dark, alpha),
                strokewidth = width, zlevel = 0.012 + 0.003 * (3 - level))
        end
        coastlines!(globe)
        static_app(fig)
    end
end

export_igeo7_hierarchy(dir = joinpath(@__DIR__, "html")) =
    export_html(joinpath(dir, "10b-igeo7-hierarchy.html"), igeo7_hierarchy_app())

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_igeo7_hierarchy()
