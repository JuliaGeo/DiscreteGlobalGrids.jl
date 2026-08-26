include(joinpath(@__DIR__, "01-longlat-grid.jl"))
include(joinpath(@__DIR__, "02-h3-grid.jl"))
include(joinpath(@__DIR__, "03-projection-distortion.jl"))
include(joinpath(@__DIR__, "04-cell-areas.jl"))

function export_series(dir = joinpath(@__DIR__, "html"))
    [export_longlat(dir), export_h3(dir), export_projection(dir), export_cell_areas(dir)...]
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_series()
