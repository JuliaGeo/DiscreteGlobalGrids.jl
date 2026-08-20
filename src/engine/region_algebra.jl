# Region algebra over the same-level region types: growth by rings, set union,
# bulk movement between levels, and compaction to a mixed-level set.
#
# Every verb here answers in leaf-grid POSITION space and hands the result to
# `_windows`/`_windows_from_intervals`, which own the `CellVector` invariants.
# Nothing assumes a cell's children are contiguous or ascending: sorted-subtree
# systems reach for `descendant_range`, and every other system goes through
# `descendants` and sorts. That branch is what makes A5 work.

"""
    grow(region, n; connectivity = Vertex()) -> CellVector

The region plus `n` rings of level-grid neighbours, at the region's own level.
`n == 0` returns the region as a [`CellVector`](@ref) and nothing else.

`region` is a [`PartialGrid`](@ref), a [`CellVector`](@ref), or a complete grid;
a `CellLookup` is grown through `parent(lk)`. Each ring is one [`halo`](@ref)
walk, so a rooted subtree grid uses the subtree engine on its first step and the
window walk afterwards — growth has no root.

Cost is `n` halo walks plus `O(m log m)` per step to merge that step's `m`
windows and halo cells. A [`MultiOrderCellSet`](@ref) has no method here: it has
no single level to walk at, so [`expand`](@ref) it first.

```julia
grow(subset, 2)                      # two rings of receptive field
grow(subset, 1; connectivity = Edge())
```
"""
function grow(region::Union{AbstractGrid,CellVector}, n::Integer;
        connectivity::Connectivity = Vertex())
    n >= 0 || throw(ArgumentError("grow needs a non-negative ring count, got $n"))
    cv = CellVector(region)
    n == 0 && return cv
    grown = _grow_once(_walkable(region), cv, connectivity)
    for _ in 2:n
        grown = _grow_once(grown, grown, connectivity)
    end
    return grown
end

# `halo` is defined for subsets, not for a complete level grid; the vector form
# of one walks nothing and answers empty.
_walkable(cv::CellVector) = cv
_walkable(pg::PartialGrid) = pg
_walkable(grid::AbstractGrid) = CellVector(grid)

function _grow_once(walkable, cv::CellVector, connectivity::Connectivity)
    ivs = intervals(cv.windows)
    for p in halo(walkable; connectivity)
        push!(ivs, (p, p))
    end
    sort!(ivs)
    return CellVector(_windows_from_intervals(ivs), cv.grid, nothing, cv.level)
end

"""
    expand(region, l::Integer) -> CellVector
    expand(set::MultiOrderCellSet, l::Integer) -> CellVector
    expand(mov::MultiOrderVector, l::Integer) -> CellVector

Every level-`l` descendant of the region's cells, as one [`CellVector`](@ref).
`l` equal to the region's own level returns it unchanged; `l` above it throws.

The expansion never assumes a cell's descendants are contiguous or ascending in
the deeper level. Where [`has_sorted_subtrees`](@ref) holds it merges one
[`descendant_range`](@ref) per cell; elsewhere (A5) it resolves
[`descendants`](@ref) to positions and sorts them. Both paths visit every cell of
the region, and the second visits every leaf it names.

The set and [`MultiOrderVector`](@ref) forms are [`CellVector`](@ref)`(x; level = l)`
— the same expansion, from the mixed-level side. [`expand`](@ref)`(A, l)` on a
`DimArray` carries the values along with it.

Descendants of a member need not lie inside the member's own footprint under
non-congruent refinement, so an expanded coverage over-covers exactly where the
refinement does; [`MultiOrderCoverage`](@ref) sizes that margin per system.
"""
function expand(cv::CellVector, l::Integer)
    target = Int(l)
    target >= cv.level || throw(ArgumentError(
        "cannot expand to level $target: the vector is already at level $(cv.level)"))
    target == cv.level && return cv
    sys = system(cv)
    grid = levelgrid(sys, target)
    return CellVector(_expand_windows(sys, cv, grid, target), grid, nothing, target)
end

expand(region::AbstractGrid, l::Integer) = expand(CellVector(region), l)

expand(set::MultiOrderCellSet, l::Integer) = CellVector(set; level = l)

expand(mov::MultiOrderVector, l::Integer) = CellVector(mov; level = l)

function _expand_windows(sys::AbstractHierarchicalGridSystem, cv::CellVector,
        grid::AbstractGrid, target::Int)
    if has_sorted_subtrees(sys)
        ivs = Vector{Tuple{Int,Int}}(undef, length(cv))
        for (k, c) in enumerate(cv)
            r = descendant_range(sys, c, target)
            @inbounds ivs[k] = (first(r), last(r))
        end
        issorted(ivs) || sort!(ivs)
        return _windows_from_intervals(ivs)
    end
    positions = Int[]
    for c in cv, d in descendants(sys, c, target)
        p = cellposition(grid, d)
        p === nothing && throw(ArgumentError(
            "descendant $d of $c is not a cell of levelgrid($sys, $target)"))
        push!(positions, p)
    end
    sort!(positions)
    return _windows(positions)
end

"""
    compact(region) -> MultiOrderCellSet

The region as a mixed-level set, with every complete sibling group replaced by
its parent, recursively. The result names the same leaves at the region's own
level, which becomes the set's reference level: `expand(compact(cv), level(cv))`
is `cv` again.

Ascent runs level by level and costs `O(k log k)` in the `k` cells still standing
at that level; a region with no complete sibling group stops after one pass. A
group counts as complete against `length(children(sys, parent))`, the parent's
own child count — pentagon parents are not assumed to have the hexagonal one.

[`iscontained`](@ref) is `false` on every member: compaction has no coverage
target, so nothing was proven to lie inside anything.
"""
function compact(cv::CellVector)
    sys = system(cv)
    ID = cellindextype(sys)
    rootlevel = first(levels(sys))
    kept = ID[]
    current = collect(cv)
    l = cv.level
    while l > rootlevel && !isempty(current)
        parents = ID[ancestor(sys, c, l - 1) for c in current]
        # Siblings are adjacent in id order only where subtrees are sorted, so
        # the grouping is done by the parent key rather than by adjacency.
        perm = sortperm(parents)
        promoted = ID[]
        i = 1
        while i <= length(perm)
            p = @inbounds parents[perm[i]]
            j = i
            while j < length(perm) && @inbounds(parents[perm[j+1]]) == p
                j += 1
            end
            if j - i + 1 == length(children(sys, p))
                push!(promoted, p)
            else
                for k in i:j
                    push!(kept, @inbounds current[perm[k]])
                end
            end
            i = j + 1
        end
        current = promoted
        l -= 1
    end
    append!(kept, current)
    return _sorted_cell_set(sys, kept, falses(length(kept)), cv.level)
end

compact(region::AbstractGrid) = compact(CellVector(region))

# --- Base set and concatenation verbs --------------------------------------
#
# `intersect` and `issubset` live in `cell_vector.jl` beside the windows they
# read. `union` and `vcat` are here because their cross-level answer is a
# `MultiOrderCellSet`, which is the algebra's own product.

"""
    union(a::CellVector, b::CellVector, rest::CellVector...)

The cells of every operand, once each.

At one level the answer is a [`CellVector`](@ref) again, merged over the stored
windows in `O(m log m)` for `m` windows and **ascending**: the container's
invariant is strict ascent, which Base's first-appearance order cannot honour
for two already-ascending operands.

Across levels the answer is a [`MultiOrderCellSet`](@ref) whose reference level
is the deepest operand's, holding each cell that no coarser operand already
covers — the set invariant that no member descends from another. Complete
sibling groups are left alone there; [`compact`](@ref) is the verb that merges
them.

Operands from different systems throw.
"""
function Base.union(a::CellVector, b::CellVector, rest::CellVector...)
    vs = (a, b, rest...)
    for v in vs
        system(v) == system(a) || throw(ArgumentError(
            "cannot union cell vectors from $(typeof(system(a))) and $(typeof(system(v)))"))
    end
    all(v -> v.level == a.level, vs) || return _mixed_union(vs)
    return _union_level(vs, a.level)
end

function _union_level(vs, l::Int)
    ivs = Tuple{Int,Int}[]
    for v in vs
        v.level == l || continue
        append!(ivs, intervals(v.windows))
    end
    sort!(ivs)
    grid = first(v for v in vs if v.level == l).grid
    return CellVector(_windows_from_intervals(ivs), grid, nothing, l)
end

function _mixed_union(vs)
    sys = system(first(vs))
    ID = cellindextype(sys)
    lvls = sort!(unique(Int[v.level for v in vs]))
    merged = [_union_level(vs, l) for l in lvls]
    cells = ID[]
    for (i, cv) in enumerate(merged), c in cv
        covered = false
        for j in 1:(i-1)
            if ancestor(sys, c, lvls[j]) in merged[j]
                covered = true
                break
            end
        end
        covered || push!(cells, c)
    end
    return _sorted_cell_set(sys, cells, falses(length(cells)), last(lvls))
end

"""
    vcat(a::CellVector, b::CellVector, rest::CellVector...)

Concatenation, as an ordered vector of ids.

A [`CellVector`](@ref) is returned when the concatenation is itself strictly
ascending at one level — the only case in which the windows can carry it — and
an ordinary `Vector` of ids otherwise, exactly as indexing a cell vector by a
non-ascending index vector does. `union` is the order-free verb, and the one to
reach for when the operands overlap.
"""
function Base.vcat(a::CellVector, b::CellVector, rest::CellVector...)
    vs = (a, b, rest...)
    if all(v -> system(v) == system(a) && v.level == a.level, vs) && _ascends(vs)
        ivs = Tuple{Int,Int}[]
        for v in vs
            append!(ivs, intervals(v.windows))
        end
        return CellVector(_windows_from_intervals(ivs), a.grid, nothing, a.level)
    end
    return reduce(vcat, map(collect, vs))
end

# Ascending concatenation needs every non-empty operand to start past the end of
# the one before it.
function _ascends(vs)
    previous = 0
    for v in vs
        n = length(v.windows)
        n == 0 && continue
        leafposition(v.windows, 1) > previous || return false
        previous = leafposition(v.windows, n)
    end
    return true
end
