isdefined(@__MODULE__, :DGGSTalkFigures) || include(joinpath(@__DIR__, "00-dggs-theme.jl"))
using .DGGSTalkFigures

# Looking down on the Arctic: the geographic pole is an ordinary cell, and the
# two singularities sit at 55°N on opposite meridians, out over land.
const TRIPOLAR_CAMERA = (camera_longlat = (ROTATED_POLE.longitude, 78),
    camera_altitude = 2.1)

tripolar_app() = globe_outline_app(tripolar_cells(96, 56); TRIPOLAR_CAMERA...)

export_tripolar(dir = joinpath(@__DIR__, "html")) =
    export_html(joinpath(dir, "06-tripolar-grid.html"), tripolar_app())

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_tripolar()
