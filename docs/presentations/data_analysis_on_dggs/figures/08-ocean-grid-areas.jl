isdefined(@__MODULE__, :DGGSTalkFigures) || include(joinpath(@__DIR__, "00-dggs-theme.jl"))
using .DGGSTalkFigures

include(joinpath(@__DIR__, "06-tripolar-grid.jl"))
include(joinpath(@__DIR__, "07-rotated-pole-grid.jl"))

# Same colour scale as 04-cell-areas, on the same cameras as 06 and 07 — moving
# the pole moves the distortion, it does not remove it.
const OCEAN_AREA_GRIDS = (
    ("tripolar", tripolar_cells(96, 56), TRIPOLAR_CAMERA),
    ("rotated-pole", rotated_pole_cells(36), ROTATED_CAMERA),
)

function export_ocean_grid_areas(dir = joinpath(@__DIR__, "html"))
    [export_html(joinpath(dir, "08-cell-area-$slug.html"),
            area_globe_app(cells; camera...))
        for (slug, cells, camera) in OCEAN_AREA_GRIDS]
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_ocean_grid_areas()
