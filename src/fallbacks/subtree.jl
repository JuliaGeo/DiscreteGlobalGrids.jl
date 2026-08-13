# ---------------------------------------------------------------------------
# Subtree rim and bulk.
#
# The generic realisation of `subtree_border` / `subtree_interior`. Both are
# written against `descendants` + `neighbors` + `ancestor`, so they work for any
# hierarchical system the moment it has those — and both are slow enough that
# every system with a border automaton should override the border.
#
# The asymmetry is deliberate. `subtree_border` is the hook: its answer is
# O(3^d) against an O(7^d) subtree, so a system that can walk the rim from its
# digits alone wins by an unbounded factor, and that is worth a method.
# `subtree_interior` is *most* of the subtree, so there is no asymptotic win to
# be had; it is written once, here, in terms of the border, which means a system
# that overrides the border also gets the faster interior with no second walker.
# ---------------------------------------------------------------------------

"""
    subtree_border(sys, c, l; connectivity = Vertex()) -> Vector

The rim of `c`'s subtree at level `l` — see the interface docstring for the
contract.

The generic implementation enumerates [`descendants`](@ref) and keeps the ones
with a neighbour outside the subtree, testing membership by walking the
neighbour back up to `level(c)` with [`ancestor`](@ref). That is
`O(subtree · degree · depth)`, which is correct for every system and the right
answer for none of them: it exists so that the hook is total, and so that a
system's automaton has something to be differentially tested against.
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

The level-`l` descendants of `c` that are not on its rim — see the interface
docstring for the contract.

Computed as [`descendants`](@ref) minus [`subtree_border`](@ref), so the border
call dispatches to whatever fast path the system has. This is the whole
implementation for every system: the interior is Θ(subtree) however it is
found, so there is nothing for an automaton to save here.
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
