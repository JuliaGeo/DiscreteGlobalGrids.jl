module Isea9rTestSuite

# Tests for the ISEA9R submodule of DiscreteGlobalGrids: the dense ten-diamond
# grids at aperture 9 (test_face_grid.jl), the geometry-only kernel wiring of
# `ISEA9RDGGS` (test_isea9r_kernel.jl), and the claim that makes both of those
# short — that the chart is `ISEA4R`'s, imported rather than reimplemented
# (test_delegation.jl).
#
# There is no `test_diamonds.jl` here, and that absence is the point: ISEA9R has
# no layout table and no chart of its own to pin. `test/ISEA4R/test_diamonds.jl`
# pins them once, and `test_delegation.jl` checks that ISEA9R really is using
# that one — `===` on the function objects, then bitwise equality of the two
# systems' grids and of the Regridders built on them.
#
# There is no `ISEA9RLookups` layer, and `Isea9rKernel.jl` wires geometry only:
# `cell_boundary` / `cell_center` / `cell_cap` / `cell_polygon` answer over the
# canonical `isea9r_ordinal` `diamond * 9^level + morton9_position`, while
# `cell_children`, `cell_parent`, `descendant_range`, `num_cells` and the rest
# of the hierarchy group still throw `NotPortedError` — deferred rather than
# blocked, on the same line the ISEA4R sibling holds.
# `test_isea9r_kernel.jl` pins both halves.
#
# The ten-ROOT layout these tests rest on is normative (OGC 21-038r1 Annex B.2,
# "The ten root rhombuses are formed by combining two icosahedron triangles at
# their base"; DGGAL `RI9R.ec` `countZones(level) = 10 * 9^level`). The
# NUMBERING of those ten, the in-diamond axes and the in-diamond index are
# package conventions with no external oracle; see
# `docs/design/isea9r_layout.md` before reading any DGGRID / DGGAL /
# SphericalSpatialTrees compatibility into them.

# Dense per-resolution diamond grids as spatial trees, the ordering contract,
# and the pre-registered cap measurement at nside 3/9/27
# (src/ISEA9R/face_grid.jl). The file wraps itself in a module of its own.
include("test_face_grid.jl")

# The operations-kernel wiring of `ISEA9RDGGS` (src/ISEA9R/Isea9rKernel.jl):
# geometry over the canonical ordinal, bitwise-identical to the diamond grid's,
# and the milestone boundary the hierarchy group still stops at. Also its own
# module.
include("test_isea9r_kernel.jl")

# The chart delegation to `ISEA4R` and the cross-system bitwise gate
# (src/ISEA9R/chart.jl). Also its own module.
include("test_delegation.jl")

println("All testsets finished.")

end # module Isea9rTestSuite
