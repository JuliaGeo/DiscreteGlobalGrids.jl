isdefined(@__MODULE__, :DGGSTalkFigures) || include(joinpath(@__DIR__, "00-dggs-theme.jl"))
using .DGGSTalkFigures

longlat_app() = globe_outline_app("Long–latitude grid",
    "A familiar grid, refined from 15° to 5° cells",
    longlat_cells(36))

export_longlat(dir = joinpath(@__DIR__, "html")) =
    export_html(joinpath(dir, "01-longlat-grid.html"), longlat_app())

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_longlat()
