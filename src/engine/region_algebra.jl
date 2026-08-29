"""
    grow(region, n; connectivity = Vertex()) -> CellVector

Add `n` neighbor rings to a same-level region and return a
[`CellVector`](@ref). Each step performs one [`halo`](@ref) walk and merges its
cells with the current windows. Expand mixed-level inputs before growing them.

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

# A complete level becomes an empty-halo subset through its `CellVector` form.
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

Return every level-`l` descendant as one [`CellVector`](@ref). Sorted-subtree
systems merge [`descendant_range`](@ref)s; other systems enumerate and sort
[`descendants`](@ref). Non-congruent refinement can over-cover the original
footprint; [`MultiOrderCoverage`](@ref) documents that margin.
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
    indices = Int[]
    for c in cv, d in descendants(sys, c, target)
        p = globalindex(grid, d)
        p === nothing && throw(ArgumentError(
            "descendant $d of $c is not a cell of levelgrid($sys, $target)"))
        push!(indices, p)
    end
    sort!(indices)
    return _windows(indices)
end

"""
    compact(region) -> MultiOrderCellSet

Replace complete sibling groups recursively with their parent and return a
[`MultiOrderCellSet`](@ref) naming the same leaves. The input level becomes the
reference level. Completeness uses each parent's actual child count.
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
        # Parent-key grouping remains correct when sibling ids are nonadjacent.
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

# --- set and concatenation --------------------------------------------------

"""
    union(a::CellVector, b::CellVector, rest::CellVector...)

Return every operand cell once. Same-level inputs produce an ascending
[`CellVector`](@ref). Cross-level inputs produce a
[`MultiOrderCellSet`](@ref) at the deepest reference level and omit cells
covered by coarser operands. [`compact`](@ref) merges complete sibling groups.
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

Concatenate operands in order. A same-level, strictly ascending result remains
a [`CellVector`](@ref); other results use an ordinary id vector. Use `union`
for order-independent overlapping inputs.
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

# Each nonempty operand must begin after its predecessor ends.
function _ascends(vs)
    previous = 0
    for v in vs
        n = length(v.windows)
        n == 0 && continue
        leafindex(v.windows, 1) > previous || return false
        previous = leafindex(v.windows, n)
    end
    return true
end
