include(joinpath(@__DIR__, "06-tripolar-grid.jl"))
include(joinpath(@__DIR__, "07-rotated-pole-grid.jl"))
include(joinpath(@__DIR__, "08-ocean-grid-areas.jl"))

function export_ocean_series(dir = joinpath(@__DIR__, "html"))
    [export_tripolar(dir), export_rotated_pole(dir), export_ocean_grid_areas(dir)...]
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_ocean_series()
