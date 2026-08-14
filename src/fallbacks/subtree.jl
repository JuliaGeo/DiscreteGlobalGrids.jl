# Generic subtree border and interior implementations. Systems with direct
# border traversal should override `subtree_border`; `subtree_interior` uses it.

"""
    subtree_border(sys, c, l; connectivity = Vertex()) -> Vector

Return level-`l` descendants of `c` with a neighbor outside the subtree. The
generic implementation costs `O(subtree · degree · depth)`.
"""
function subtree_border(sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex, l::Integer;
        connectivity::Connectivity=Vertex())
    lc = level(c)
    target = Int(l)
    target >= lc || throw(ArgumentError(
        "subtree_border: level $target is above the cell's own level $lc"))
    T = typeof(c)
    target == lc && return T[c]

    grid = levelgrid(sys, target)
    out = T[]
    for d in descendants(sys, c, target)
        for nb in neighbors(grid, d, 1; connectivity)
            # A neighbour is outside the subtree exactly when its ancestor at
            # `c`'s level is not `c`. Cheaper than a membership set over the
            # whole subtree, and it never materialises one.
            if ancestor(sys, nb, lc) != c
                push!(out, d)
                break
            end
        end
    end
    return out
end

"""
    subtree_interior(sys, c, l; connectivity = Vertex()) -> Vector

Return level-`l` descendants of `c` excluding [`subtree_border`](@ref). The
result costs `Θ(subtree)` and uses any specialized border implementation.
"""
function subtree_interior(sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex, l::Integer;
        connectivity::Connectivity=Vertex())
    lc = level(c)
    target = Int(l)
    target >= lc || throw(ArgumentError(
        "subtree_interior: level $target is above the cell's own level $lc"))
    T = typeof(c)
    # A cell is its own rim, so a depth-0 subtree has no interior at all.
    target == lc && return T[]

    rim = Set{T}(subtree_border(sys, c, target; connectivity))
    return T[d for d in descendants(sys, c, target) if !(d in rim)]
end
