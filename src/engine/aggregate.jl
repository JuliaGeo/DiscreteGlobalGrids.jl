# Keep the default mean independent of the optional `Statistics` package.
_mean(vs) = sum(vs) / length(vs)

function _needs_ranges(sys::AbstractHierarchicalGridSystem, verb::AbstractString)
    has_sorted_subtrees(sys) || throw(ArgumentError(
        "$verb requires contiguous descendant ranges, which $(typeof(sys)) does not provide"))
    return nothing
end

function _paired(cv::CellVector, values::AbstractVector)
    length(values) == length(cv) || throw(ArgumentError(
        "values must line up with the cells they belong to: got $(length(values)) " *
        "values for $(length(cv)) cells"))
    return nothing
end

function _within_levels(sys::AbstractHierarchicalGridSystem, l::Int)
    l in levels(sys) || throw(ArgumentError(
        "level $l is outside $(typeof(sys))'s levels $(levels(sys))"))
    return nothing
end

# --- window helpers ---------------------------------------------------------
function _next_index(w::RangeWindows, lo::Int)
    j = searchsortedfirst(w.stops, lo)
    j <= length(w.stops) || return nothing
    return max(lo, @inbounds w.starts[j])
end

function _next_index(w::IndexWindows, lo::Int)
    j = searchsortedfirst(w.indices, lo)
    j <= length(w.indices) || return nothing
    return @inbounds w.indices[j]
end

# A complete subtree maps to one contiguous range in the parallel values vector.
function _complete_segment(w::CellWindows, r::UnitRange{Int})
    a = windowindex(w, first(r))
    a === nothing && return nothing
    b = windowindex(w, last(r))
    b === nothing && return nothing
    b - a == length(r) - 1 || return nothing
    return a:b
end

# Iterate intervals lazily because `IndexWindows` has one interval per cell.
eachinterval(w::RangeWindows) =
    ((@inbounds(w.starts[j]), @inbounds(w.stops[j])) for j in eachindex(w.starts))

eachinterval(w::IndexWindows) = ((p, p) for p in w.indices)

"""
    aggregate(f, cv::CellVector, values::AbstractVector, l::Integer) -> (CellVector, Vector)

Reduce present descendant values to level `l`. Each represented level-`l`
ancestor contributes one ascending output cell and `f(view(values, range))`.
Partial groups reduce only their present values; absent groups produce no
output. The target must be coarser than `level(cv)`, and the system must support
[`has_sorted_subtrees`](@ref). [`coarsen`](@ref) selects levels adaptively.
"""
function aggregate(f, cv::CellVector, values::AbstractVector, l::Integer)
    sys = system(cv)
    _needs_ranges(sys, "aggregated to a coarser level")
    _paired(cv, values)
    target = Int(l)
    _within_levels(sys, target)
    target < level(cv) || throw(ArgumentError(
        "aggregation goes UP the hierarchy: level $target is not coarser than the " *
        "vector's level $(level(cv))"))
    coarse = levelgrid(sys, target)
    indices, segments = _aggregate_segments(cv, coarse, target)
    out = map(r -> f(view(values, r)), segments)
    # `_windows` verifies the ordering promised by `has_sorted_subtrees`.
    return CellVector(_windows(indices), coarse, nothing, target), out
end

# Specialize on the window shape and compute each output ancestor once.
_aggregate_segments(cv::CellVector, coarse::AbstractGrid, target::Int) =
    _aggregate_segments(cv, cv.windows, coarse, target)

function _aggregate_segments(cv::CellVector, w::CellWindows, coarse::AbstractGrid,
        target::Int)
    sys = system(cv)
    leaf = level(cv)
    indices = Int[]
    segments = UnitRange{Int}[]
    k = 0
    kstart = 0
    stop = 0
    for (lo, hi) in eachinterval(w)
        p = lo
        while p <= hi
            if p > stop
                stop == 0 || push!(segments, kstart:k)
                a = ancestor(sys, cellindex(cv.grid, p), target)
                stop = last(descendant_range(sys, a, leaf))
                push!(indices, _coarse_index(coarse, a))
                kstart = k + 1
            end
            n = min(hi, stop) - p + 1
            k += n
            p += n
        end
    end
    stop == 0 || push!(segments, kstart:k)
    return indices, segments
end

function _coarse_index(grid::AbstractGrid, c::AbstractCellIndex)
    p = localindex(grid, c)
    p === nothing && throw(ArgumentError(
        "$c is not a cell of levelgrid($(system(grid)), $(level(grid)))"))
    return p
end

"""
    coarsen(cv::CellVector, values::AbstractVector; atol, by = mean,
            minlevel = first(levels(system(cv)))) -> (MultiOrderVector, Vector)

Build an adaptive mesh by replacing each complete eligible subtree with its
coarsest ancestor. A subtree is eligible when its values are all `missing` or
satisfy `maximum(vs) - minimum(vs) <= atol`. The stored value is `by(vs)` over
leaf values; `by` defaults to the mean. The result uses `level(cv)` as its
reference level and expands to the original cell windows. See
[`aggregate`](@ref) and [`expand`](@ref).
"""
function coarsen(cv::CellVector, values::AbstractVector; atol, by=_mean,
        minlevel::Integer=first(levels(system(cv))))
    cells, vals = _coarsen(cv, values; atol, by, minlevel)
    # `_coarsen` preserves value alignment by emitting cells in interval order.
    return MultiOrderVector(system(cv), cells; reference_level=level(cv)), vals
end

"""
    _coarsen(cv, values; atol, by, minlevel) -> (Vector{ID}, Vector)

Return [`coarsen`](@ref)'s stored cells in descendant-range order and their
aligned values.
"""
function _coarsen(cv::CellVector, values::AbstractVector; atol, by=_mean,
        minlevel::Integer=first(levels(system(cv))))
    sys = system(cv)
    _needs_ranges(sys, "coarsened")
    _paired(cv, values)
    stop_level = Int(minlevel)
    _within_levels(sys, stop_level)
    stop_level <= level(cv) || throw(ArgumentError(
        "minlevel $stop_level is deeper than the vector's level $(level(cv)); " *
        "coarsening only ever climbs"))
    cells = cellindextype(sys)[]
    # The function barrier specializes the walk on the accumulator's runtime type.
    vals = _accumulator(values, by)
    for c in rootcells(sys)
        _coarsen_visit!(cells, vals, cv, values, c, stop_level, atol, by)
    end
    return cells, _narrow(vals, values)
end

# Cover input and summary values with a concrete accumulator element type.
function _accumulator(values::AbstractVector, by)
    T = Union{eltype(values),Base.promote_op(by, typeof(view(values, 1:0)))}
    return Vector{T}(undef, 0)
end

function _coarsen_visit!(cells, vals, cv::CellVector, values::AbstractVector,
        c::AbstractCellIndex, minlevel::Int, atol, by)
    sys = system(cv)
    w = cv.windows
    leaf = level(cv)
    r = descendant_range(sys, c, leaf)
    p = _next_index(w, first(r))
    (p === nothing || p > last(r)) && return nothing
    lc = level(c)
    if lc == leaf
        # `by` summarizes merges; unmerged leaves retain their original value.
        push!(cells, c)
        push!(vals, values[windowindex(w, p)])
        return nothing
    end
    if lc >= minlevel
        ks = _complete_segment(w, r)
        if ks !== nothing
            merges, allmissing = _merges(values, ks, atol)
            if merges
                push!(cells, c)
                # All-missing groups have no non-missing value to summarize.
                push!(vals, allmissing ? missing : by(view(values, ks)))
                return nothing
            end
        end
    end
    # Descendant ranges enumerate children without allocating a child vector.
    kids = levelgrid(sys, lc + 1)
    for q in descendant_range(sys, c, lc + 1)
        _coarsen_visit!(cells, vals, cv, values, cellindex(kids, q), minlevel, atol, by)
    end
    return nothing
end

# Two booleans keep the merge test concrete while distinguishing all-missing groups.
function _merges(values::AbstractVector, ks::UnitRange{Int}, atol)
    vs = view(values, ks)
    n = length(vs)
    v1 = @inbounds vs[1]
    if ismissing(v1)
        for i in 2:n
            ismissing(@inbounds vs[i]) || return false, false
        end
        return true, true
    end
    lo = hi = v1
    for i in 2:n
        v = @inbounds vs[i]
        ismissing(v) && return false, false
        # `min` and `max` propagate NaN and make the merge criterion fail.
        lo, hi = min(lo, v), max(hi, v)
        _within(lo, hi, atol) || return false, false
    end
    # Settles the one-element group, which never enters the loop.
    _within(lo, hi, atol) || return false, false
    return true, false
end

# The sign check rejects wrapped integer spans and NaN differences.
function _within(lo, hi, atol)
    d = hi - lo
    return d >= zero(d) && d <= atol
end

# `typeof(v) === T` avoids Julia 1.12 folding `isa` on a `Union{}` local.
function _narrow(vals::Vector{S}, values::AbstractVector) where {S}
    isempty(vals) && return similar(values, 0)
    isconcretetype(S) && return vals
    T = Union{}
    for v in vals
        U = typeof(v)
        U === T || (T = Union{T,U})
    end
    T === S && return vals
    out = Vector{T}(undef, length(vals))
    copyto!(out, vals)
    return out
end

"""
    expand(A, l)

Present a mixed-level array at level `l`, reindexing each stored value over the
cells it covers. The lazy result retains one value per multi-order cell.
Methods accept arrays carrying a [`Cells`](@ref) axis.
"""
function expand end
