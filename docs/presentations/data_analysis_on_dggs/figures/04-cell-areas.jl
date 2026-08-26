isdefined(@__MODULE__, :DGGSTalkFigures) || include(joinpath(@__DIR__, "00-dggs-theme.jl"))
using .DGGSTalkFigures
using DiscreteGlobalGrids

const AREA_GRIDS = (
    ("long–latitude", "longlat", longlat_cells(36)),
    ("H3", "h3", dggs_cells(H3DGGS(), 1)),
    ("S2", "s2", dggs_cells(S2DGGS(), 3)),
    ("IGEO7", "igeo7", dggs_cells(IGEO7DGGS(), 2)),
    ("ISEA4R", "isea4r", dggs_cells(ISEA4RDGGS(), 3)),
)

area_app(_, cells) = area_globe_app(cells)

function export_cell_areas(dir = joinpath(@__DIR__, "html"))
    [export_html(joinpath(dir, "04-cell-area-$slug.html"), area_app(name, cells))
        for (name, slug, cells) in AREA_GRIDS]
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_cell_areas()
