# Nested HEALPix implemented directly on its twelve equal-area face charts.

"""
    DiscreteGlobalGrids.HEALPix

The [HEALPix](https://healpix.sourceforge.io) equal-area grid: twelve base
pixels with aperture-4 refinement at levels `0:29`, where `nside = 2^level`
and `npix = 12 * 4^level`.

[`HEALPixSystem`](@ref) is the system singleton; HEALPix ships no grid type of
its own, so `levelgrid(HEALPixSystem(), l)` returns the package's
[`HierarchicalLevelGrid`](@ref).

Canonical `LevelIndex` ids are 0-based NESTED pixel numbers; alternate
[`HEALPixRingIndex`](@ref) ids are 1-based RING numbers. Nested subtrees are
contiguous, parent footprints equal their child unions, and every cell has
exact area `4π / npix`. Location uses the analytic chart inverse; topology uses
the Morton lattice. This module does not depend on `Healpix.jl`.
"""
module HEALPix

import ..DiscreteGlobalGrids as DGG

import GeometryOps as GO

import SmallCollections
using SmallCollections: SmallVector

include("chart.jl")
include("neighbors.jl")
include("system.jl")
include("border.jl")

export HEALPixSystem, HEALPixRingIndex

end # module HEALPix
