# ---------------------------------------------------------------------------
# The subtree rim, walked on the lattice
#
# `subtree_border(sys, c, l)` has an O(rim) closed form here, against the
# generic fallback's O(subtree · degree · depth), and the reason is the same one
# that makes `descendant_range` exact: a cell's level-`l` descendants are a
# SQUARE BLOCK of the level-`l` lattice, `2^Δ` on a side, inside one diamond.
#
#     ix, iy at level lc   ->   [ix·2^Δ, (ix+1)·2^Δ) × [iy·2^Δ, (iy+1)·2^Δ)
#
# A descendant has a neighbour outside the subtree exactly when it sits on that
# block's boundary ring — and this is where the seam topology pays off for a
# second time: a block that touches a diamond rim has neighbours in ANOTHER
# diamond, which are outside the subtree just as surely as an in-diamond
# neighbour would be, so the rule needs no special case for a root cell or for a
# block against a rim.
#
# `Vertex()` and `Edge()` give the SAME rim, which is why the argument is
# accepted and then ignored. A cell strictly inside the block has all eight
# lattice neighbours in the block; a cell on the boundary ring has at least one
# AXIS neighbour outside it, so it is on the rim under either connectivity. The
# corner-only offsets never decide a membership question here.
#
# Order is ascending canonical order, as the interface documents for the
# default. The block is enumerated in lattice order and then sorted, which costs
# `O(rim log rim)` — nothing next to the `O(4^Δ)` the fallback pays, and it
# keeps the emission order a fact about `isless` rather than about the walk.
# ---------------------------------------------------------------------------

"""
    subtree_border(ISEA4RSystem(), c, l; connectivity = Vertex()) -> Vector{LevelIndex}

The rim of `c`'s subtree at level `l`: the boundary ring of the `2^Δ × 2^Δ`
lattice block the subtree occupies, in ascending canonical order. `4·2^Δ - 4`
cells against the block's `4^Δ`.

`connectivity` does not change the answer — see the file header — and is
accepted so that generic code can pass it through.

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
    for j in 0:(m - 1), i in 0:(m - 1)
        (i == 0 || i == m - 1 || j == 0 || j == m - 1) || continue
        out[k += 1] = DGG.LevelIndex(target, xyd_to_morton(x0 + i, y0 + j, d, nside))
    end
    return sort!(out)
end
