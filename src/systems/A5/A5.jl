# ---------------------------------------------------------------------------
# T10 — A5, step 1: the arithmetic.
#
# The pure-Julia port of upstream a5 carried across from `src/A5/A5Native.jl`,
# so that the interface wiring in the next commit has something to be a wiring
# OF. Nothing here implements a contract yet.
# ---------------------------------------------------------------------------

"""
    DiscreteGlobalGrids.A5

The [A5](https://a5geo.org) discrete global grid system: an equal-area
pentagonal tiling of a dodecahedron.

At this point the module is the ported upstream arithmetic and nothing else —
see [`A5Native`](@ref). The grid-interface wiring lands on top of it.
"""
module A5

# The ported upstream arithmetic. Everything this module gains later is a
# wiring of it.
include("native.jl")

end # module A5
