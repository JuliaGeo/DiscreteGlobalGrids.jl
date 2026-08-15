"""
    HaloOracle

Brute-force reference implementations of the halo, built the naive way a user
would build one if `subtree_halo` did not exist: walk every cell of the subset,
take its one-ring, keep whatever falls outside.

Two directions are computed separately, because they are only the same set if
the native adjacency relation is symmetric:

* `inside_out` — from each member, collect outside neighbours.
* `outside_in` — for each cell of the whole grid that is outside, ask whether
  any of ITS neighbours is a member. This is the definition `halo.jl`'s
  `IndexedNeighbors` provider actually implements.

Any disagreement between the two is itself a finding.
"""
module HaloOracle

import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: Vertex, Edge, Connectivity

export oracle_subtree_halo, oracle_subset_halo, subtree_cells, member_cells,
       adjacency_asymmetry

"The level-`l` descendants of `c`, via `descendant_range` + `cellindex`."
function subtree_cells(sys, c, l::Integer)
    g = DGG.levelgrid(sys, l)
    return [DGG.cellindex(g, p) for p in DGG.descendant_range(sys, c, l)]
end

"""
    oracle_subtree_halo(sys, root, l, connectivity) -> (inside_out, outside_in)

Both directions, each sorted ascending by position.
"""
function oracle_subtree_halo(sys, root, l::Integer, connectivity::Connectivity)
    g = DGG.levelgrid(sys, l)
    r = DGG.descendant_range(sys, root, l)
    lo, hi = first(r), last(r)

    io = Set{typeof(DGG.cellindex(g, 1))}()
    for p in r
        d = DGG.cellindex(g, p)
        for n in DGG.neighbors(g, d, 1; connectivity)
            pn = DGG.cellposition(g, n)
            (lo <= pn <= hi) || push!(io, n)
        end
    end

    oi = Set{typeof(DGG.cellindex(g, 1))}()
    for p in 1:DGG.ncells(g)
        (lo <= p <= hi) && continue
        x = DGG.cellindex(g, p)
        for n in DGG.neighbors(g, x, 1; connectivity)
            pn = DGG.cellposition(g, n)
            if lo <= pn <= hi
                push!(oi, x); break
            end
        end
    end

    key = c -> DGG.cellposition(g, c)
    return sort!(collect(io); by = key), sort!(collect(oi); by = key)
end

"""
    oracle_subset_halo(g, members, connectivity) -> (inside_out, outside_in)

Same, for an arbitrary same-level member set (a `Set` of cells).
"""
function oracle_subset_halo(g, members::AbstractSet, connectivity::Connectivity)
    C = eltype(members)
    io = Set{C}()
    for m in members
        for n in DGG.neighbors(g, m, 1; connectivity)
            n in members || push!(io, n)
        end
    end
    oi = Set{C}()
    for p in 1:DGG.ncells(g)
        x = DGG.cellindex(g, p)
        x in members && continue
        for n in DGG.neighbors(g, x, 1; connectivity)
            if n in members
                push!(oi, x); break
            end
        end
    end
    key = c -> DGG.cellposition(g, c)
    return sort!(collect(io); by = key), sort!(collect(oi); by = key)
end

"""
    adjacency_asymmetry(sys, l, connectivity) -> Vector{Tuple}

Every ordered pair `(a, b)` with `b ∈ neighbors(a)` but `a ∉ neighbors(b)`.
Exhaustive over the level grid, so keep `l` small.
"""
function adjacency_asymmetry(sys, l::Integer, connectivity::Connectivity)
    g = DGG.levelgrid(sys, l)
    bad = Tuple{Any,Any}[]
    for p in 1:DGG.ncells(g)
        a = DGG.cellindex(g, p)
        for b in DGG.neighbors(g, a, 1; connectivity)
            a in DGG.neighbors(g, b, 1; connectivity) || push!(bad, (a, b))
        end
    end
    return bad
end

end # module
