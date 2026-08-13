# ---------------------------------------------------------------------------
# T11 — S2. Chart layer.
#
# The six cube-face charts, ported from the pre-redesign `src/S2/chart.jl`.
# Nothing here depends on the s2geometry library and no s2geometry code is
# vendored. The system layer (`S2System`, the hierarchy, neighbours) lands on
# top of this in the next commit.
# ---------------------------------------------------------------------------

module S2

import ..DiscreteGlobalGrids as DGG

import GeometryOps as GO
const US = GO.UnitSpherical

include("chart.jl")

end # module S2
