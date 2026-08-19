# A nested id is
# `face * nside^2 + morton(ix, iy)` (`xyf_to_nested`), and refining Δ levels
# scales the lattice by `s = 2^Δ` on the SAME face: pixel `(ix, iy, face)` at
# `level` covers exactly `[ix*s, (ix+1)*s) x [iy*s, (iy+1)*s)` of the leaf
# lattice. Bit-interleaving is positional, so that scaling is a shift of the
# Morton code:
#
#     morton(ix*s + dx, iy*s + dy) = morton(ix, iy) * s^2 + morton(dx, dy)
#
# So the subtree is the contiguous id range `[id * 4^Δ, (id+1) * 4^Δ)` and a
# descendant's offset within it *is* `morton(dx, dy)` — the quad-face family's
# premise, which supplies the rim, interior, and halo engines. Seam crossings
# leave the subtree face (no non-centre row of `NB_FACEARRAY` maps a face to
# itself), and missing diagonals at degree-3 vertices do not affect rim
# membership, so both connectivities yield the same rim.
#
# The seam band is `NB_SWAPARRAY` territory — ids on other faces — but no seam
# table is read here: `square_halo_engine` derives the candidate rectangles by
# asking `neighbors` about a few rim cells and filtering every candidate with the
# native one-ring.

# Square-walk hooks. Morton orientation is `0x0` on every face, which is the
# family default, so only the codec pair is written here.
DGG.lattice_decode(sys::HEALPixSystem, c::DGG.LevelIndex) =
    nested_to_xyf(c.index, DGG.nside(DGG.level(c)))

DGG.lattice_cell(sys::HEALPixSystem, l::Int, ix::Integer, iy::Integer,
    face::Integer) = DGG.LevelIndex(l, xyf_to_nested(ix, iy, face, DGG.nside(l)))

DGG.face_orientation(sys::HEALPixSystem, face::Integer) = 0x0
