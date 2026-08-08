module S2TestSuite

# Tests for the S2 submodule of DiscreteGlobalGrids: the closed-form cube-face
# charts and the row-major / Hilbert index maps (test_chart.jl), the dense face
# grids as spatial trees (test_face_grid.jl), and the geometry-only kernel
# wiring of `S2DGGS` (test_s2_kernel.jl).
#
# There is no `S2Lookups` layer, and `S2Kernel.jl` wires geometry only: the
# `s2_cellid` id hierarchy is a later milestone, so `cell_boundary` /
# `cell_center` / `cell_cap` / `cell_polygon` answer over the registry's
# scaffold ordinal `face * 4^level + hilbert_position` while `cell_children`,
# `cell_parent`, `descendant_range`, `num_cells` and the rest of that group
# still throw `NotPortedError`. `test_s2_kernel.jl` pins both halves.

# Pure closed-form cube-face charts + row-major/Hilbert index maps
# (src/S2/chart.jl); the file wraps itself in a module of its own.
include("test_chart.jl")

# Dense per-resolution face grids built on the charts, as spatial trees
# (src/S2/face_grid.jl); also its own module.
include("test_face_grid.jl")

# The operations-kernel wiring of `S2DGGS` (src/S2/S2Kernel.jl): geometry over
# scaffold ordinals, bitwise-identical to the face grid's, and the milestone
# boundary the hierarchy group still stops at. Also its own module.
include("test_s2_kernel.jl")

println("All testsets finished.")

end # module S2TestSuite
