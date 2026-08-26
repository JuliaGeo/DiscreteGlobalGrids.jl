isdefined(@__MODULE__, :DGGSTalkFigures) || include(joinpath(@__DIR__, "00-dggs-theme.jl"))
using .DGGSTalkFigures
using DiscreteGlobalGrids

# One camera for all four, so the panels compare rather than each showing off.
# Coarser than the full-slide versions in 06 and 07: a quarter-canvas globe
# turns a 96x56 grid into a grey smudge.
more_dggs_app() = globe_panels_app([
    ("TRIPOLAR", tripolar_cells(48, 28)),
    ("ROTATED POLE", rotated_pole_cells(18)),
    ("A5", levelgrid(A5System(), 3)),
    ("ISEA4R", levelgrid(ISEA4RSystem(), 3)),
])

export_more_dggs(dir = joinpath(@__DIR__, "html")) =
    export_html(joinpath(dir, "10-more-dggs.html"), more_dggs_app())

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_more_dggs()
