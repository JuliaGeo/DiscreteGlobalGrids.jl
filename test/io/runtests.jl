# The store-IO layer: grid id arithmetic, pluggable cell encodings, and the
# chunk manifest the lazy lookup is planned against.
#
# `src/io/` is not yet included by `src/DiscreteGlobalGrids.jl`. Until those
# include lines land, this suite loads the sources into its own namespace; the
# `isdefined` branch picks the package's own submodules once they are wired, so
# the suite does not change when they are.

module DGGIOTests

using Test
import DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG

if isdefined(DiscreteGlobalGrids, :Encodings)
    using DiscreteGlobalGrids: Encodings, ChunkedLookups, DGGSFormatError,
        store_context, with_store_context
else
    # `errors.jl` first: both submodules `import ..DGGSFormatError` from the
    # module that includes them, which is this one.
    include(joinpath(@__DIR__, "..", "..", "src", "io", "errors.jl"))
    include(joinpath(@__DIR__, "..", "..", "src", "io", "encodings.jl"))
    include(joinpath(@__DIR__, "..", "..", "src", "io", "chunked_lookup.jl"))
end

using .Encodings
using .ChunkedLookups

@testset "io" begin
    include("encodings.jl")
    include("chunked_lookup.jl")
    # The conventions suite carries its own namespace and its own source
    # loading, so it is only run from here, never wired into this one.
    isfile(joinpath(@__DIR__, "conventions.jl")) && include("conventions.jl")
    # The Zarr-extension suites need `using Zarr` and self-skip when absent.
    isfile(joinpath(@__DIR__, "read.jl")) && include("read.jl")
    isfile(joinpath(@__DIR__, "write.jl")) && include("write.jl")
end

end # module DGGIOTests
