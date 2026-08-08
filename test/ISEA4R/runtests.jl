module Isea4rTestSuite

# Tests for the ISEA4R submodule of DiscreteGlobalGrids: the ten-diamond layout
# table and the piecewise-affine rhombus chart built on `ISEA`'s Snyder
# machinery (test_diamonds.jl), the dense diamond grids as spatial trees
# (test_face_grid.jl), and the geometry-only kernel wiring of `ISEA4RDGGS`
# (test_isea4r_kernel.jl).
#
# There is no `ISEA4RLookups` layer, and `Isea4rKernel.jl` wires geometry only:
# `cell_boundary` / `cell_center` / `cell_cap` / `cell_polygon` answer over the
# canonical `isea4r_ordinal` `diamond * 4^level + morton_position`, while
# `cell_children`, `cell_parent`, `descendant_range`, `num_cells` and the rest
# of the hierarchy group still throw `NotPortedError` — deferred rather than
# blocked. `test_isea4r_kernel.jl` pins both halves.
#
# The layout these tests pin is a package convention with no external oracle
# behind it; see `docs/design/isea4r_diamond_layout.md` before reading any
# compatibility with an external ISEA4R product's identifiers into the
# numbering — no such compatibility is claimed or tested here.

# The ten-diamond table + the rhombus chart and its index maps
# (src/ISEA4R/diamonds.jl, src/ISEA4R/chart.jl); the file wraps itself in a
# module of its own.
include("test_diamonds.jl")

# Dense per-resolution diamond grids built on the charts, as spatial trees
# (src/ISEA4R/face_grid.jl); also its own module.
include("test_face_grid.jl")

# The operations-kernel wiring of `ISEA4RDGGS` (src/ISEA4R/Isea4rKernel.jl):
# geometry over the canonical ordinal, bitwise-identical to the diamond grid's,
# and the milestone boundary the hierarchy group still stops at. Also its own
# module.
include("test_isea4r_kernel.jl")

println("All testsets finished.")

end # module Isea4rTestSuite
