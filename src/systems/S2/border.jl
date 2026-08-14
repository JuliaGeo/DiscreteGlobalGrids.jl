# A subtree at depth Δ is the aligned `2^Δ x 2^Δ` lattice block on one face, and
# the scaffold ordinal refines by `index * 4^Δ + pos`, so a descendant's offset
# within the subtree's contiguous ordinal range IS its Hilbert position inside
# the block (`descendant_range`, `xyf_to_hilbert`'s nesting note). Ascending
# ordinal order over the subtree is therefore Hilbert order over the block.
#
# The rim is the block's perimeter, `4*2^Δ - 4` cells: in the 3x3 lattice
# neighbourhood a descendant has every neighbour inside the block iff
# `0 < dx < s-1` and `0 < dy < s-1`. Both connectivities agree — a perimeter cell
# has an AXIS neighbour outside already — and a block flush against a face edge
# is no special case, because its neighbours across the seam are on another face
# and so outside the subtree regardless.
#
# So the walk is the package's `SquareRimEngine`, shared with HEALPix and ISEA4R;
# only the curve differs. Hilbert's quadrant order is orientation-dependent,
# which is exactly what `quadrant_step` carries: `POS_TO_IJ` reads the curve
# position back to a quadrant and `POS_TO_ORIENTATION` advances the state, the
# same two tables `hilbert_to_xyf` descends with. Pure index arithmetic — no
# geometry, and no neighbour query at any level.

import ..DiscreteGlobalGrids: subtree_border

"""
    HilbertCurve()

S2's per-face Hilbert curve, as the quadrant order [`DGG.SquareRimEngine`](@ref)
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

# `descendant_range` is the only level guard needed — it raises both
# `ArgumentError`s — and its `lo` is the subtree's first ordinal.
function _s2_square(sys::S2System, c::DGG.LevelIndex, target::Int)
    r = DGG.descendant_range(sys, c, target)
    lc = DGG.level(c)
    _checked_index(c)
    lo = Int64(first(r)) - 1
    orientation = _hilbert_orientation(Int64(c.index), _nside(lc))
    return lo, Int64(1) << (target - lc), orientation
end

function DGG.rim_engine(sys::S2System, c::DGG.LevelIndex, target::Int,
        connectivity::DGG.Connectivity)
    lo, side, orientation = _s2_square(sys, c, target)
    return DGG.SquareRimEngine(HilbertCurve(), lo, target, side, orientation)
end

function DGG.interior_engine(sys::S2System, c::DGG.LevelIndex, target::Int,
        connectivity::DGG.Connectivity)
    lo, side, orientation = _s2_square(sys, c, target)
    return DGG.SquareInteriorEngine(HilbertCurve(), lo, target, side, orientation)
end

"""
    subtree_border(S2System(), c, l; connectivity = Vertex()) -> Vector{LevelIndex}

The **rim** of `c`'s subtree at level `l`: the perimeter of the `2^Δ × 2^Δ`
lattice block, `4*2^Δ - 4` cells in ascending ordinal order, and `[c]` at
`Δ == 0`. `O(rim)` time in `O(Δ)` memory, by index arithmetic alone.
`collect` of [`EdgeCellIterator`](@ref), which is the same walk lazily.

`connectivity` does not change the result: a perimeter cell already has an axis
neighbour outside the block.
"""
subtree_border(sys::S2System, c::DGG.LevelIndex, l::Integer;
    connectivity::DGG.Connectivity = DGG.Vertex()) =
    DGG.collect_subtree(DGG.EdgeCellIterator(sys, c, l; connectivity))
