# A subtree is a square lattice block. Enumerate its boundary directly and
# sort canonically, giving `O(rim log rim)` cost under either connectivity.
# A boundary cell always has an AXIS neighbour outside the block — across a
# diamond rim as readily as within one — so neither connectivity nor a block
# sitting against a seam needs a special case.

"""
    subtree_border(ISEA4RSystem(), c, l; connectivity = Vertex()) -> Vector{LevelIndex}

Return the boundary of the `2^Δ × 2^Δ` descendant block in ascending order:
`4*2^Δ - 4` cells. Connectivity does not change this rim. Cost is
`O(rim log rim)`.

`subtree_border(sys, c, level(c))` is `[c]`; `l < level(c)` throws an
`ArgumentError`, uniformly with every other system.
"""
function DGG.subtree_border(sys::ISEA4RSystem, c::DGG.LevelIndex, l::Integer;
        connectivity::DGG.Connectivity = DGG.Vertex())
    lc = DGG.level(c)
    target = Int(l)
    target >= lc || throw(ArgumentError(
        "subtree_border: level $target is above the cell's own level $lc"))
    target <= DGG.max_level(sys) || throw(ArgumentError(
        "subtree_border: level $target is past max_level $(DGG.max_level(sys))"))
    target == lc && return DGG.LevelIndex[c]

    ix, iy, d = morton_to_xyd(c.index, _nside(lc))
    m = 1 << (target - lc)              # the block's side, in level-`target` cells
    nside = _nside(target)
    x0 = ix * m
    y0 = iy * m
    out = Vector{DGG.LevelIndex}(undef, 4m - 4)
    k = 0
    # Enumerate the two rows and intervening columns in `O(m)`, without scanning
    # the `O(m^2)` block interior.
    @inbounds begin
        for i in 0:(m - 1)                  # the j == 0 row, in full
            out[k += 1] = DGG.LevelIndex(target,
                xyd_to_morton_unchecked(x0 + i, y0, d, nside))
        end
        for j in 1:(m - 2)                  # the two columns between the rows
            out[k += 1] = DGG.LevelIndex(target,
                xyd_to_morton_unchecked(x0, y0 + j, d, nside))
            out[k += 1] = DGG.LevelIndex(target,
                xyd_to_morton_unchecked(x0 + m - 1, y0 + j, d, nside))
        end
        for i in 0:(m - 1)                  # the j == m-1 row, in full
            out[k += 1] = DGG.LevelIndex(target,
                xyd_to_morton_unchecked(x0 + i, y0 + m - 1, d, nside))
        end
    end
    return sort!(out)
end
