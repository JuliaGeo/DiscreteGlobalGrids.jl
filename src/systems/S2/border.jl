# The subtree engines are the quad-face family's; S2 supplies the curve they walk
# in. The scaffold ordinal refines by `index * 4^Δ + pos`, so a descendant's
# offset within the subtree's contiguous ordinal range IS its Hilbert position
# inside the block (`descendant_range`, `xyf_to_hilbert`'s nesting note), and
# ascending ordinal order over the subtree is Hilbert order over it.
#
# Hilbert's quadrant order is orientation-dependent, which `quadrant_step`
# carries: `POS_TO_IJ` reads the curve position back to a quadrant and
# `POS_TO_ORIENTATION` advances the state, the same two tables `hilbert_to_xyf`
# descends with. Pure index arithmetic — no geometry, and no neighbour query at
# any level.

"""
    HilbertCurve()

S2's per-face Hilbert curve, as the quadrant order [`DGG.SquareBorderEngine`](@ref)
descends by. Unlike Morton, its state advances with depth.
"""
struct HilbertCurve end

@inline function DGG.quadrant_step(::HilbertCurve, orientation::UInt8, p::Int)
    ij = @inbounds POS_TO_IJ[orientation+1][p+1]
    return ((ij >> 1) & 1, ij & 1,
        orientation ⊻ UInt8(@inbounds POS_TO_ORIENTATION[p+1]))
end

"""
    _hilbert_orientation(index, nside) -> UInt8

The curve orientation state at the cell with scaffold ordinal `index`, i.e. the
state `xyf_to_hilbert` has reached once it has consumed that cell's position
bits. It is the orientation its children are read under.
"""
function _hilbert_orientation(index::Int64, nside::Int64)
    npface = nside^2
    face, pos = divrem(index, npface)
    k = trailing_zeros(nside)
    orientation = isodd(face) ? SWAP_MASK : 0
    for b in (k-1):-1:0
        p = Int((pos >> (2b)) & 3)
        orientation ⊻= POS_TO_ORIENTATION[p+1]
    end
    return UInt8(orientation)
end

# The subtree block's own root is read under the state `xyf_to_hilbert` has
# reached at that cell. The halo descent instead starts at the FACE root, before
# any position bits are consumed, so `face_orientation` seeds it separately —
# on a face the block crossed a seam onto, `_hilbert_orientation` of the subtree
# cell would be the wrong state entirely.
DGG.subtree_curve(::S2System) = HilbertCurve()

DGG.subtree_orientation(::S2System, c::DGG.LevelIndex) =
    _hilbert_orientation(Int64(c.index), DGG.nside(DGG.level(c)))

# Square-walk hooks. At a cube corner `wrap_xyf` returns `nothing` and there is
# simply no diagonal candidate for the seam band to find.
DGG.lattice_decode(sys::S2System, c::DGG.LevelIndex) =
    hilbert_to_xyf(c.index, DGG.nside(DGG.level(c)))

DGG.lattice_cell(sys::S2System, l::Int, ix::Integer, iy::Integer,
    face::Integer) = DGG.LevelIndex(l, xyf_to_hilbert(ix, iy, face, DGG.nside(l)))

DGG.face_orientation(sys::S2System, face::Integer) =
    isodd(face) ? UInt8(SWAP_MASK) : 0x0
