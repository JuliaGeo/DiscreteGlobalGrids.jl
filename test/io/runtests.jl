# The store-IO layer: grid id arithmetic, pluggable cell encodings, and the
# chunk manifest the lazy lookup is planned against.
#
# The two layer suites share this namespace, because they are two halves of one
# thing: an encoding decides what the axis looks like and the lookup resolves
# against it. The convention suite and the two Zarr-extension suites carry
# namespaces of their own and are only included from here.

module DGGIOTests

using Test
import DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG

using DiscreteGlobalGrids: Encodings, ChunkedLookups, DGGSFormatError,
    store_context, with_store_context

using .Encodings
using .ChunkedLookups

@testset "io" begin
    include("encodings.jl")
    include("chunked_lookup.jl")
    include("conventions.jl")
    # The Zarr-extension suites need `using Zarr` and self-skip when absent.
    include("read.jl")
    include("write.jl")
    include("stores.jl")
end

end # module DGGIOTests
