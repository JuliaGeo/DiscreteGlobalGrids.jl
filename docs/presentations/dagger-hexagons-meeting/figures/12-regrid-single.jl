isdefined(@__MODULE__, :DaggerTalkFigures) ||
    include(joinpath(@__DIR__, "11-dagger-theme.jl"))

using .DGGSTalkFigures, .DaggerTalkFigures
using DiscreteGlobalGridsVisualization, Makie
import SparseArrays

function regrid_single_app()
    data = regrid_fixture(; nt = 4)
    weights = data.plan.block.weights
    rows, cols, vals = SparseArrays.findnz(weights)

    with_theme(JG_THEME) do
        fig, body = slide_figure()

        source_panel = GridLayout(body[1, 1])
        Label(source_panel[1, 1], "SOURCE RASTER · 24 × 12"; fontsize = 12,
            font = FONT_BODY, color = JG.ink, tellwidth = false)
        # The input is still a regular longitude–latitude raster; putting it
        # on the sphere makes its geographic relationship to the IGeo7
        # destination explicit.
        # Both geographic views use the same viewport width and camera.
        source = globe_axis(source_panel[2, 1]; camera_longlat = (20, 18),
            camera_altitude = 2.50)
        plot_cells!(source, longlat_cells(length(data.lat));
            color = vec(parent(data.data)), colormap = FIELD_COLORMAP,
            colorrange = FIELD_RANGE, strokecolor = (JG.ink, 0.20),
            strokewidth = 0.28)
        coastlines!(source)
        rowgap!(source_panel, 2)

        matrix_panel = GridLayout(body[1, 2])
        Label(matrix_panel[1, 1], "REGRIDDER"; fontsize = 12,
            font = FONT_BODY, color = JG.ink, tellwidth = false)
        matrix = Axis(matrix_panel[2, 1];
            xlabel = "source cells", ylabel = "destination cells",
            # Keep the reusable sparse operator legible, but make it a compact
            # bridge between the two geographic representations.
            aspect = 0.65)
        scatter!(matrix, cols, rows; color = vals, colormap = FIELD_COLORMAP,
            colorrange = extrema(vals), markersize = 2.2)
        xlims!(matrix, 1, size(weights, 2)); ylims!(matrix, size(weights, 1), 1)
        matrix.xticks = ([1, size(weights, 2)], ["1", string(size(weights, 2))])
        matrix.yticks = ([1, size(weights, 1)], ["1", string(size(weights, 1))])
        rowsize!(matrix_panel, 2, Relative(0.56))
        rowgap!(matrix_panel, 2)

        target_panel = GridLayout(body[1, 3])
        Label(target_panel[1, 1], "DESTINATION · 492 CELLS"; fontsize = 12,
            font = FONT_BODY, color = JG.ink, tellwidth = false)
        target = globe_axis(target_panel[2, 1]; camera_longlat = (20, 18),
            camera_altitude = 2.50)
        dggpoly!(target, data.grid; color = vec(parent(data.result)),
            colormap = FIELD_COLORMAP, colorrange = FIELD_RANGE,
            strokecolor = (JG.ink, 0.35), strokewidth = 0.55)
        coastlines!(target)
        rowgap!(target_panel, 2)

        # Plot insertion normally fits each GlobeAxis to its own geometry.
        # Lock both cameras after all layers exist so the generic raster mesh
        # and the native DGGS mesh cannot acquire different final zooms.
        for axis in (source, target)
            axis.center[] = false
            cameracontrols(axis.scene).settings.center[] = false
            Makie.update_cam!(axis; longlat = (20, 18), altitude = 2.50,
                fov = 45.0)
        end

        # Mirror the geographic panels around a compact central operator.
        colsize!(body, 1, Relative(0.38))
        colsize!(body, 2, Relative(0.24))
        colsize!(body, 3, Relative(0.38))
        colgap!(body, 16)
        static_app(fig)
    end
end

export_regrid_single(dir = joinpath(@__DIR__, "html")) =
    export_html(joinpath(dir, "12-regrid-single.html"), regrid_single_app())

abspath(PROGRAM_FILE) == (@__FILE__) && export_regrid_single()
