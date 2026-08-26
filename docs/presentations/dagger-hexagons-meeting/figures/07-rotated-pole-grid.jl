isdefined(@__MODULE__, :DGGSTalkFigures) || include(joinpath(@__DIR__, "00-dggs-theme.jl"))
using .DGGSTalkFigures

# Same 72x36 grid as 01-longlat-grid, with the pole moved to central Asia — so
# the singularity is in frame instead of hiding at the top of the sphere.
const ROTATED_CAMERA = (camera_longlat = (ROTATED_POLE.longitude, 40),)

rotated_pole_app() = globe_outline_app(rotated_pole_cells(36); ROTATED_CAMERA...)

export_rotated_pole(dir = joinpath(@__DIR__, "html")) =
    export_html(joinpath(dir, "07-rotated-pole-grid.html"), rotated_pole_app())

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_rotated_pole()
