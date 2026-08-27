module DGGSTalkFigures

using Bonito, GeoMakie, Makie, WGLMakie
using DiscreteGlobalGrids
using DiscreteGlobalGridsVisualization
import ConservativeRegridding as CR
import Geodesy
import GeometryOps as GO
import GeoInterface as GI
import Oceananigans

export JG, AREA_COLORMAP, AREA_RANGE, area_globe_app, cell_area_ratios
export export_html, globe_outline_app, longlat_cells, projection_app
export ROTATED_POLE, globe_panels_app, rotated_pole_cells, tripolar_cells
export JG_THEME, FONT_BODY, FONT_DISPLAY, FONT_BOLD
export slide_figure, globe_axis, coastlines!, plot_cells!, static_app

const JG = (
    green = colorant"#389826", green_dark = colorant"#2c7a1e",
    green_50 = colorant"#f0faea", green_100 = colorant"#dcf5d7",
    purple = colorant"#9558b2", purple_100 = colorant"#eadaf2",
    red = colorant"#cb3c33", blue = colorant"#4063d8",
    ink = colorant"#212529", muted = colorant"#495057",
    caption = colorant"#6c757d", hairline = colorant"#dee2e6",
    paper = colorant"#ffffff", paper_off = colorant"#f8f9fa",
)

const FONT_DIR = normpath(joinpath(@__DIR__, "..", "public", "fonts"))
const FONT_BODY = joinpath(FONT_DIR, "Inter-Variable.ttf")
const FONT_DISPLAY = joinpath(FONT_DIR, "RaleGrotesk-Medium.otf")
const FONT_BOLD = joinpath(FONT_DIR, "RaleGrotesk-Bold.otf")

const JG_THEME = Theme(
    fonts = (; regular = FONT_BODY, bold = FONT_BOLD, display = FONT_DISPLAY),
    fontsize = 16, textcolor = JG.ink, backgroundcolor = JG.paper,
    Axis = (; backgroundcolor = JG.paper, xgridcolor = JG.hairline,
        ygridcolor = JG.hairline, spinecolor = JG.hairline,
        xticklabelcolor = JG.caption, yticklabelcolor = JG.caption),
)

const AREA_COLORMAP = [JG.purple, JG.purple_100, JG.paper, JG.green_100, JG.green]
const AREA_RANGE = (0.0, 2.0)
const CARTESIAN_SPHERE = "+proj=cart +R=1"

WGLMakie.activate!(; framerate = 30)

function slide_figure()
    fig = Figure(; size = (960, 540), figure_padding = 14,
        backgroundcolor = JG.paper)
    body = GridLayout(fig[1, 1])
    colsize!(fig.layout, 1, Relative(1))
    rowsize!(fig.layout, 1, Relative(1))
    return fig, body
end

function globe_axis(slot; camera_longlat = (20, 18), camera_altitude = 1.9)
    axis = GlobeAxis(slot; source = "+proj=longlat +R=1",
        dest = Geodesy.Ellipsoid(; a = "1", b = "1"), camera_longlat,
        camera_altitude, backgroundvisible = false)
    meshimage!(axis, -180..180, -90..90, fill(JG.green_50, 1, 1);
        zlevel = -0.05, npoints = 300)
    return axis
end

coastlines!(axis) = lines!(axis, GeoMakie.coastlines(); color = (JG.ink, 0.72),
    linewidth = 1.0, zlevel = 0.01)

# Every Oceananigans grid reaches Makie the same way: ConservativeRegridding
# treeifies it into unit-sphere cell polygons, one per (i, j) column.  A folded
# grid (tripolar) carries degenerate ghost cells along the fold row so the field
# dimensions still line up — they have no area and nothing to draw, so drop them.
degenerate_cell(cell) = allequal(GI.getpoint(cell))

function oceananigans_cells(grid)
    cells = collect(CR.Trees.getcell(CR.Trees.treeify(GO.Spherical(), grid))) |> vec
    filter(!degenerate_cell, cells)
end

function longlat_cells(nlat)
    oceananigans_cells(Oceananigans.LatitudeLongitudeGrid(; size = (2nlat, nlat, 1),
        longitude = (-180, 180), latitude = (-90, 90), z = (0, 1)))
end

# The two singularities a tripolar grid pushes off the geographic pole, and the
# one a rotated-pole grid moves.  Shared so the cameras can aim at them.
const ROTATED_POLE = (longitude = 70, latitude = 55)

function tripolar_cells(nlon, nlat)
    oceananigans_cells(Oceananigans.TripolarGrid(; size = (nlon, nlat, 1), z = (0, 1),
        first_pole_longitude = ROTATED_POLE.longitude,
        north_poles_latitude = ROTATED_POLE.latitude))
end

function rotated_pole_cells(nlat; north_pole = values(ROTATED_POLE))
    oceananigans_cells(Oceananigans.RotatedLatitudeLongitudeGrid(;
        size = (2nlat, nlat, 1), north_pole,
        longitude = (-180, 180), latitude = (-90, 90), z = (0, 1)))
end

function cell_area_ratios(cells)
    geographic = GO.transform(GO.GeographicFromUnitSphere(), cells)
    areas = GO.area.(Ref(GO.Spherical(; radius = 1.0)), geographic)
    areas ./ (4pi / length(areas))
end

function cell_area_ratios(grid::AbstractGrid)
    cells = CellVector(grid)
    areas = [cell_area(grid, cell) for cell in cells]
    areas ./ (4pi / length(areas))
end

# Package grids stay as grids all the way into Makie.  The specialized recipe
# reads `cell_boundary` itself and builds one mesh for the full set.  The
# Oceananigans comparison grids are ordinary unit-sphere polygons, so they keep
# the generic Makie path.
plot_cells!(axis, grid::AbstractGrid; kwargs...) = dggpoly!(axis, grid; kwargs...)
plot_cells!(axis, cells; kwargs...) = poly!(axis, cells; source = CARTESIAN_SPHERE, kwargs...)

const WEB_STYLE = Styles(
    CSS("html, body", "width" => "100%", "height" => "100%", "margin" => "0",
        "display" => "grid", "place-items" => "center", "overflow" => "hidden",
        "background" => "#f8f9fa"),
    CSS(".dggs-root", "position" => "relative",
        "width" => "960px", "height" => "540px",
        "overflow" => "hidden", "background" => "#ffffff"),
)

static_app(fig) = App() do _
    DOM.div(WEB_STYLE, WGLMakie.WithConfig(fig; resize_to = :parent); class = "dggs-root")
end

function globe_outline_app(cells; camera...)
    with_theme(JG_THEME) do
        fig, body = slide_figure()
        axis = globe_axis(body[1, 1]; camera...)
        plot_cells!(axis, cells; color = JG.green_100,
            strokecolor = JG.green_dark, strokewidth = 0.65)
        coastlines!(axis)
        static_app(fig)
    end
end

globe_outline_app(::AbstractString, ::AbstractString, cells) = globe_outline_app(cells)

# A gallery of small globes, one labelled panel per grid, filled row-major.
# Each panel is `(label, cells)`; every globe faces the same way so the set
# reads as a comparison, and a panel that needs its own view can pass a third
# element, the named tuple `globe_axis` takes.  The Label needs
# `tellwidth = false`: on default settings it drives the column width instead,
# and a caption longer than the globe squeezes the GlobeAxis to a clipped sliver.
function globe_panels_app(panels; ncols = 2, camera = NamedTuple())
    with_theme(JG_THEME) do
        fig, body = slide_figure()
        for (k, panel_spec) in enumerate(panels)
            label, cells = panel_spec[1], panel_spec[2]
            row, col = fldmod1(k, ncols)
            panel = GridLayout(body[row, col])
            Label(panel[1, 1], label; fontsize = 11, color = JG.caption,
                tellwidth = false)
            axis = globe_axis(panel[2, 1];
                (length(panel_spec) > 2 ? panel_spec[3] : camera)...)
            plot_cells!(axis, cells; color = JG.green_100,
                strokecolor = JG.green_dark, strokewidth = 0.5)
            coastlines!(axis)
            rowgap!(panel, 2)
        end
        rowgap!(body, 4); colgap!(body, 4)
        static_app(fig)
    end
end

function area_globe_app(cells; colorrange = AREA_RANGE,
    ticks = ([0.5, 1.0, 1.5], ["0.5×", "1×", "1.5×"]), camera...)
    with_theme(JG_THEME) do
        fig, body = slide_figure()
        axis = globe_axis(body[1, 1]; camera...)
        area = plot_cells!(axis, cells;
            color = cell_area_ratios(cells), colormap = AREA_COLORMAP,
            colorrange, strokecolor = (JG.ink, 0.42), strokewidth = 0.6)
        coastlines!(axis)
        Colorbar(body[1, 2], area; label = "area / global mean", width = 12,
            ticks, labelsize = 14, ticklabelsize = 12)
        colsize!(body, 2, 58)
        static_app(fig)
    end
end

function tissot_ring(lon, lat, radius; n = 48)
    lambda0, phi0, delta = deg2rad.((lon, lat, radius))
    [begin
        theta = 2pi * k / n
        phi = asin(sin(phi0) * cos(delta) + cos(phi0) * sin(delta) * cos(theta))
        lambda = lambda0 + atan(sin(theta) * sin(delta) * cos(phi0),
            cos(delta) - sin(phi0) * sin(phi))
        rad2deg(lambda), rad2deg(phi)
    end for k in 0:n]
end

function projection_app()
    with_theme(JG_THEME) do
        fig, body = slide_figure()
        axis = Axis(body[1, 1]; aspect = 2, backgroundcolor = JG.paper)
        hidedecorations!(axis); hidespines!(axis)
        vlines!(axis, -150:30:150; color = JG.hairline, linewidth = 0.7)
        hlines!(axis, -60:30:60; color = JG.hairline, linewidth = 0.7)
        centers = vec(collect(Iterators.product(-150:50:150, -75:25:75)))
        poly!(axis, [GI.Polygon([tissot_ring(lon, lat, 7)]) for (lon, lat) in centers];
            color = (JG.green, 0.13), strokecolor = JG.green_dark,
            strokewidth = 1.0, transparency = true)
        lines!(axis, GeoMakie.coastlines(); color = (JG.ink, 0.78), linewidth = 0.8)
        xlims!(axis, -180, 180); ylims!(axis, -90, 90)
        static_app(fig)
    end
end

function export_html(path, app)
    mkpath(dirname(path))
    Bonito.export_static(path, app)
    println("wrote ", path)
    return path
end

end
