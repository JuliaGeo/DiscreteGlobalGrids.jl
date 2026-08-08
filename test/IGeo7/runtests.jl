# runtests.jl — IGeo7 test entry point.
#
# Part of the DiscreteGlobalGrids.jl test suite; run it from the package's test
# environment, either through the top-level test/runtests.jl or directly:
#   julia --project=<DiscreteGlobalGrids.jl> test/IGeo7/runtests.jl
#
# The suite lives in its own module so the generic vocabulary the systems share
# (`cell_center`, `lonlat_to_cell`, ...) cannot collide with another system's.
# The bindings established below — `Helpers`, `ISEA` (the shared icosahedron /
# Snyder machinery) and `IGeo7` (Z7 + the composed grid) — are what the included
# test files reach for. Each suite below is picked up only if its file exists,
# so partially implemented checkouts still run green (spec/design.md section 9
# builds the modules in parallel).

module IGeo7TestSuite

using Test

using DiscreteGlobalGrids
using DiscreteGlobalGrids: Helpers, ISEA, IGeo7
using .IGeo7

const TEST_FILES = (
    "test_icosahedron.jl",
    "test_z7.jl",
    "test_engine.jl",
    "test_chart.jl",
    "test_grid.jl",
    "test_indexing.jl",
    "test_lookups.jl",
    "test_igeo7_kernel.jl",
)

@testset "IGeo7" begin
    for f in TEST_FILES
        path = joinpath(@__DIR__, f)
        if isfile(path)
            include(path)
        else
            @info "skipping $f (not present yet)"
        end
    end
end

end # module IGeo7TestSuite
