# A subtree is a square lattice block, and a Morton id refines by
# `id * 4^Δ + morton(dx, dy)` (see `descendant_range`), so a descendant's offset
# within the subtree's contiguous id range is its Morton code in the block —
# the quad-face family's premise, which supplies the border, interior, and halo
# engines. A boundary cell always has an axis neighbour outside the block —
# across a diamond border as readily as within one — so neither connectivity nor a
# block sitting against a seam needs a special case.
#
# The seam band crosses `EDGE_NEIGHBORS` onto other diamonds, and at an
# icosahedral vertex `CORNER_FANS` onto a third, but no seam table is read here:
# `square_halo_engine` derives those candidates by asking `neighbors` about a few
# border cells and filters every one of them with the native one-ring.

# Square-walk hooks. Morton orientation is `0x0` on every diamond, which is the
# family default, so only the codec pair is written here.
DGG.lattice_decode(sys::ISEA4RSystem, c::DGG.LevelIndex) =
    morton_to_xyd(c.index, DGG.nside(DGG.level(c)))

DGG.lattice_cell(sys::ISEA4RSystem, l::Int, ix::Integer, iy::Integer,
    diamond::Integer) = DGG.LevelIndex(l, xyd_to_morton(ix, iy, diamond, DGG.nside(l)))

DGG.face_orientation(sys::ISEA4RSystem, diamond::Integer) = 0x0
