module A5

import ..Helpers

include("A5Native.jl")
include("A5Lookups.jl")
# Operations-kernel wiring for `A5DGGS` (see `src/core/kernel.jl`).
include("A5Kernel.jl")

end
