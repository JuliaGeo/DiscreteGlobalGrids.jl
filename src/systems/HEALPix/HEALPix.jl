# ---------------------------------------------------------------------------
# T6 — HEALPix.
#
# Nested HEALPix on the grid interface. The layers, bottom up:
#
#   chart.jl      the twelve equal-area face charts and the nested/ring/xyf
#                 codecs — pure closed forms, no Healpix.jl. Copied wholesale
#                 from the pre-redesign `src/HEALPix/chart.jl`, plus the
#                 chart's analytic inverse (`point_to_xyf`).
#   neighbors.jl  the 3x3 lattice neighbourhood on the nested id, lifted off
#                 the old lookup layer and off Healpix.jl.
#   system.jl     `HEALPixSystem`, `HEALPixGrid`, `HEALPixRingIndex`, and every
#                 interface method.
#   border.jl     the Morton rim walk — a subtree's border with no neighbour
#                 queries.
#
# Note the capitalisation: this submodule is `HEALPix`, distinct from the
# registered `Healpix.jl` package. Nothing here depends on that package; it is
# used only by this system's TESTS, as an independent oracle for the codecs and
# for point location.
#
# Package-level exports are T7's job. Until then, reach these names as
# `DiscreteGlobalGrids.HEALPix.<name>`.
# ---------------------------------------------------------------------------

module HEALPix

import ..DiscreteGlobalGrids as DGG

import GeometryOps as GO
const US = GO.UnitSpherical
# The UnionAll, never `SphericalCap{Float64}`: the two-argument
# `(point, radius)` constructor is a method of the UnionAll, and a parametrised
# alias would silently reach the three-field default constructor instead.
const SphericalCap = US.SphericalCap

import SmallCollections
using SmallCollections: SmallVector

include("chart.jl")
include("neighbors.jl")
include("system.jl")
include("border.jl")

export HEALPixSystem, HEALPixGrid, HEALPixRingIndex

end # module HEALPix
