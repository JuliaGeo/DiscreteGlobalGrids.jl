isdefined(@__MODULE__, :DGGSTalkFigures) ||
    include(joinpath(@__DIR__, "00-dggs-theme.jl"))
using .DGGSTalkFigures

import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import GeometryOps as GO
using Bonito, GeoMakie, Makie
import WGLMakie

# The theme's own helpers, which it keeps unexported; this figure borrows the
# globe and the palette but builds its own page around them.
const TF = DGGSTalkFigures

# Copernicus DEM fights the meridian convergence by dropping longitude columns
# at 50/60/70/80/85 degrees, so a pixel narrows across a band and springs back
# at every band edge.  GLO-90 is 69 billion pixels, far past anything a globe
# can show, so this is the same lattice coarsened until one cell is a shape.
# 2.5-degree rows put every band edge on a row boundary, and 120 columns at the
# equator keep every reduced count (120, 80, 60, 40, 24, 12) whole.  Halving
# the row and doubling the columns together keeps both true, one quadtree level
# at a time.
const ROW_STEP = 2.5
const EQUATOR_COLS = 120

# `band_factor2` is the package's own table, doubled so the division stays
# exact; reading it here keeps the figure honest if the band rule ever moves.
# It indexes whole degrees, and a row never straddles an edge, so flooring the
# row's equatorward latitude lands in that row's own band.
row_columns(lat_equatorward) =
    2 * EQUATOR_COLS ÷ DGG.CopernicusDEM.band_factor2(floor(Int, lat_equatorward))

# The parallels bounding a cell are small circles, and the polar row is 30
# degrees wide, so the long edges are sampled: drawn corner to corner they would
# cut visibly inside the globe.
function lattice_cell(lon_w, width, lat_s, lat_n; step = 2.5)
    lons = range(lon_w, lon_w + width; length = max(2, ceil(Int, width / step) + 1))
    ring = [(lon, lat_s) for lon in lons]
    append!(ring, ((lon, lat_n) for lon in Iterators.reverse(lons)))
    push!(ring, (lon_w, lat_s))
    return GI.Polygon([GI.LinearRing(ring)])
end

function copdem_lattice_cells()
    cells = map(-90:ROW_STEP:(90 - ROW_STEP)) do lat_s
        lat_n = lat_s + ROW_STEP
        ncols = row_columns(min(abs(lat_s), abs(lat_n)))
        width = 360 / ncols
        [lattice_cell(-180 + k * width, width, lat_s, lat_n) for k in 0:(ncols - 1)]
    end
    return GO.transform(GO.UnitSphereFromGeographic(), reduce(vcat, cells))
end

# Nothing embeds this one, so it fills the window instead of sitting in the
# deck's 960x540 card.  WGLMakie bakes the scene at the figure's size and a
# static export has no Julia left to redo the layout, so a canvas stretched to
# the window just leaves the globe drawn in one corner of it: the page keeps the
# canvas at the size the figure was drawn and scales that box to fit.  Drawing
# at 1920x1080 means the usual window scales it down rather than up.
const STAGE_SIZE = (1920, 1080)

const STAGE_STYLE = Styles(
    CSS("html, body", "width" => "100%", "height" => "100%", "margin" => "0",
        "overflow" => "hidden", "background" => "#ffffff"),
    CSS(".dggs-stage", "position" => "fixed", "left" => "50%", "top" => "50%",
        "width" => "$(STAGE_SIZE[1])px", "height" => "$(STAGE_SIZE[2])px",
        "transform-origin" => "center center", "background" => "#ffffff"),
)

const FIT_SCRIPT = """
const fit = () => {
    const stage = document.querySelector(".dggs-stage");
    if (!stage) return;
    const k = Math.min(window.innerWidth / $(STAGE_SIZE[1]),
                       window.innerHeight / $(STAGE_SIZE[2]));
    stage.style.transform = "translate(-50%, -50%) scale(" + k + ")";
};
window.addEventListener("resize", fit);
window.addEventListener("DOMContentLoaded", fit);
fit();
"""

# Line weights are the card's, scaled with the figure so the lattice reads the
# same at twice the size.
function copdem_lattice_figure(cells; camera...)
    Makie.with_theme(TF.JG_THEME) do
        fig = Figure(; size = STAGE_SIZE, figure_padding = 28,
            backgroundcolor = JG.paper)
        axis = TF.globe_axis(fig[1, 1]; camera...)
        TF.plot_cells!(axis, cells; color = JG.green_100,
            strokecolor = JG.green_dark, strokewidth = 1.3)
        lines!(axis, GeoMakie.coastlines(); color = (JG.ink, 0.72),
            linewidth = 2.0, zlevel = 0.01)
        fig
    end
end

function copdem_lattice_app(; camera_longlat = (20, 42), camera...)
    fig = copdem_lattice_figure(copdem_lattice_cells(); camera_longlat, camera...)
    App() do _
        # No `resize_to`: it measures the parent after the transform and would
        # hand the scene a canvas smaller than it was drawn for.  The canvas
        # keeps the figure's own size and the transform does the fitting.
        DOM.div(STAGE_STYLE, WGLMakie.WithConfig(fig),
            DOM.script(FIT_SCRIPT); class = "dggs-stage")
    end
end

export_copdem_lattice(dir = joinpath(@__DIR__, "html")) =
    export_html(joinpath(dir, "11-copdem-lattice.html"), copdem_lattice_app())

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_copdem_lattice()
