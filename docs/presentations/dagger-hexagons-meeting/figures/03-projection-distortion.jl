isdefined(@__MODULE__, :DGGSTalkFigures) || include(joinpath(@__DIR__, "00-dggs-theme.jl"))
using .DGGSTalkFigures

export_projection(dir = joinpath(@__DIR__, "html")) =
    export_html(joinpath(dir, "03-projection-distortion.html"), projection_app())

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_projection()
