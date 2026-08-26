isdefined(@__MODULE__, :DGGSTalkFigures) || include(joinpath(@__DIR__, "00-dggs-theme.jl"))
using .DGGSTalkFigures
using DiscreteGlobalGrids

h3_app() = globe_outline_app(levelgrid(H3System(), 1))

export_h3(dir = joinpath(@__DIR__, "html")) =
    export_html(joinpath(dir, "02-h3-grid.html"), h3_app())

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_h3()
