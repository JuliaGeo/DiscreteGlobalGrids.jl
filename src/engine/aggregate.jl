# Aggregation over the hierarchy: `aggregate` reduces leaf data to one fixed
# coarser level; `coarsen` merges complete sibling groups within a tolerance
# into a `MultiOrderVector` + values; `expand` is the inverse presentation.
# Both verbs read and return a `(CellVector, values)` pair, usable without the
# DimensionalData layer. Both require `has_sorted_subtrees`: a subtree must be
# one interval of the leaf level, so that a group's values are a contiguous
# slice and completeness is two window lookups.

# `coarsen`'s default summary; equals `Statistics.mean` wherever both apply
# (`Statistics` is not a dependency).
_mean(vs) = sum(vs) / length(vs)

# Without descendant ranges a sibling group is not a slice of anything.
function _needs_ranges(sys::AbstractHierarchicalGridSystem, verb::AbstractString)
    has_sorted_subtrees(sys) || throw(ArgumentError(
        "$(typeof(sys)) has no descendant ranges, so leaf data cannot be $verb"))
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

# --- the two window questions the verbs below ask --------------------------

# The first position `cv` holds at or after `lo`, or `nothing`.
function _next_position(w::RangeWindows, lo::Int)
    j = searchsortedfirst(w.stops, lo)
    j <= length(w.stops) || return nothing
    return max(lo, @inbounds w.starts[j])
end

function _next_position(w::PositionWindows, lo::Int)
    j = searchsortedfirst(w.positions, lo)
    j <= length(w.positions) || return nothing
    return @inbounds w.positions[j]
end

# The values under a subtree as a range of data-array positions, or `nothing`
# when the subtree is not complete in `cv`: both endpoints present and exactly
# `length(r)` positions between them (the windows are strictly ascending).
function _complete_segment(w::CellWindows, r::UnitRange{Int})
    a = windowposition(w, first(r))
    a === nothing && return nothing
    b = windowposition(w, last(r))
    b === nothing && return nothing
    b - a == length(r) - 1 || return nothing
    return a:b
end

# Lazy counterpart of `intervals`, which materialises one tuple per window —
# one per cell on `PositionWindows`.
eachinterval(w::RangeWindows) =
    ((@inbounds(w.starts[j]), @inbounds(w.stops[j])) for j in eachindex(w.starts))

eachinterval(w::PositionWindows) = ((p, p) for p in w.positions)

# ===========================================================================
# Fixed level
# ===========================================================================

"""
    aggregate(f, cv::CellVector, values::AbstractVector, l::Integer) -> (CellVector, Vector)

Reduce leaf data to level `l`: one output cell per distinct level-`l` ancestor
of `cv`'s cells, ascending, each carrying `f` of the values of its **present**
descendants. `values` is indexed against `cv` — `values[k]` is the datum of
`cv[k]` — and the answer is another such pair:

```julia
pyramid = [aggregate(sum, cv, data, l) for l in level(cv)-1:-1:3]
```

Partial groups are reduced over what is there; a group `cv` misses entirely is
absent from the output. `f` is handed a contiguous `AbstractVector` view into
`values`, never a copy, and sees the values as they stand: `mean` over a group
holding a `missing` answers `missing`, and `v -> mean(skipmissing(v))` is the
other reading.

`l` must be strictly coarser than `level(cv)`, and the system must have
[`has_sorted_subtrees`](@ref). See also [`coarsen`](@ref), which picks the
level per cell instead of fixing one.
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
    positions, segments = _aggregate_segments(cv, coarse, target)
    out = map(r -> f(view(values, r)), segments)
    # `_windows` re-checks that the coarse positions ascend, which
    # `has_sorted_subtrees` guarantees; a system breaking that gets an error.
    return CellVector(_windows(positions), coarse, nothing, target), out
end

# One forward pass; `ancestor` runs once per output cell, not once per leaf.
# Segments are ranges even when a group's leaf positions have gaps: `k` counts
# positions `cv` holds, so a gap in the leaves is no gap in the concatenation.
# The windows are a separate argument so the loop specialises on one shape
# rather than the `Union` in `cv.windows`.
_aggregate_segments(cv::CellVector, coarse::AbstractGrid, target::Int) =
    _aggregate_segments(cv, cv.windows, coarse, target)

function _aggregate_segments(cv::CellVector, w::CellWindows, coarse::AbstractGrid,
        target::Int)
    sys = system(cv)
    leaf = level(cv)
    positions = Int[]
    segments = UnitRange{Int}[]
    k = 0        # values consumed
    kstart = 0   # first value of the open group
    stop = 0     # last leaf POSITION of the open group; 0 while none is open
    for (lo, hi) in eachinterval(w)
        p = lo
        while p <= hi
            if p > stop
                stop == 0 || push!(segments, kstart:k)
                a = ancestor(sys, cellindex(cv.grid, p), target)
                stop = last(descendant_range(sys, a, leaf))
                push!(positions, _coarse_position(coarse, a))
                kstart = k + 1
            end
            n = min(hi, stop) - p + 1
            k += n
            p += n
        end
    end
    stop == 0 || push!(segments, kstart:k)
    return positions, segments
end

function _coarse_position(grid::AbstractGrid, c::AbstractCellIndex)
    p = cellposition(grid, c)
    p === nothing && throw(ArgumentError(
        "$c is not a cell of levelgrid($(system(grid)), $(level(grid)))"))
    return p
end

# ===========================================================================
# Adaptive
# ===========================================================================

"""
    coarsen(cv::CellVector, values::AbstractVector; atol, by = mean,
            minlevel = first(levels(system(cv)))) -> (MultiOrderVector, Vector)

Merge each subtree whose leaf values agree to within `atol` into one coarse
cell: an adaptively refined mesh built from data at one level. Cells where the
field is flat are stored once at whatever level flatness reaches; the rest
stay at `level(cv)` with their original values.

A cell between `minlevel` and `level(cv)` replaces its subtree iff:

  - **complete** — `cv` holds *every* leaf descendant of it. Incomplete groups
    never merge, so `CellVector(mov; level = level(cv))` recovers `cv` window
    for window, whatever `atol`.
  - **within tolerance** — over the leaf values `vs` under it, either all are
    `missing` (the merged value is `missing`), or none is and
    `maximum(vs) - minimum(vs) <= atol`. A group mixing `missing` with data
    never merges.

The criterion is monotone, so each leaf's stored cell is its coarsest merging
ancestor. The stored value is `by(vs)` over the LEAF values, never over child
summaries — on an equal-area system the default mean is the exact
area-weighted mean of the region. With the default `by` every leaf value is
within `atol` of the value stored for it; any `by` bounded by `extrema(vs)`
keeps that bound. The default equals `Statistics.mean` (`Statistics` is not a
dependency).

`values` is indexed against `cv`; the system must have
[`has_sorted_subtrees`](@ref); `minlevel` must be no deeper than `level(cv)`.
The container comes back with `reference_level = level(cv)`. See also
[`aggregate`](@ref), which fixes one output level, and [`expand`](@ref), the
inverse presentation.
"""
function coarsen(cv::CellVector, values::AbstractVector; atol, by=_mean,
        minlevel::Integer=first(levels(system(cv))))
    cells, vals = _coarsen(cv, values; atol, by, minlevel)
    # `_coarsen` emits cells already ascending by descendant-range start, so
    # the constructor's sort is the identity and `vals` stays aligned.
    return MultiOrderVector(system(cv), cells; reference_level=level(cv)), vals
end

"""
    _coarsen(cv, values; atol, by, minlevel) -> (Vector{ID}, Vector)

[`coarsen`](@ref)'s core: the stored cells, ascending by descendant-range
start at `level(cv)`, and their values.
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
    # Function barrier: the accumulator's eltype is a runtime value, and the
    # walk below is specialised on the concrete vector.
    vals = _accumulator(values, by)
    for c in rootcells(sys)
        _coarsen_visit!(cells, vals, cv, values, c, stop_level, atol, by)
    end
    return cells, _narrow(vals, values)
end

# Eltype covering everything the walk can push — an element of `values` or
# `by`'s answer — so values land in a concrete vector instead of a boxing
# `Vector{Any}`. A widening only: `_narrow` reports the narrowest eltype
# actually emitted.
function _accumulator(values::AbstractVector, by)
    T = Union{eltype(values),Base.promote_op(by, typeof(view(values, 1:0)))}
    return Vector{T}(undef, 0)
end

# Top-down: by monotonicity the first merging cell on the way down is the
# coarsest merging ancestor of everything beneath it, and empty subtrees are
# pruned, so the walk is O(#nodes the data touches). The pre-order descent
# emits cells already sorted by descendant-range start.
function _coarsen_visit!(cells, vals, cv::CellVector, values::AbstractVector,
        c::AbstractCellIndex, minlevel::Int, atol, by)
    sys = system(cv)
    w = cv.windows
    leaf = level(cv)
    r = descendant_range(sys, c, leaf)
    p = _next_position(w, first(r))
    (p === nothing || p > last(r)) && return nothing
    lc = level(c)
    if lc == leaf
        # A leaf keeps its own value; `by` summarises merges only, so
        # `atol = 0` on distinct data is the identity.
        push!(cells, c)
        push!(vals, values[windowposition(w, p)])
        return nothing
    end
    if lc >= minlevel
        ks = _complete_segment(w, r)
        if ks !== nothing
            merges, allmissing = _merges(values, ks, atol)
            if merges
                push!(cells, c)
                # `by` runs only for committed merges; an all-`missing`
                # group stores `missing`.
                push!(vals, allmissing ? missing : by(view(values, ks)))
                return nothing
            end
        end
    end
    # A cell's next-level descendant range names exactly its children in
    # order, without the per-node vector `children` can allocate.
    kids = levelgrid(sys, lc + 1)
    for q in descendant_range(sys, c, lc + 1)
        _coarsen_visit!(cells, vals, cv, values, cellindex(kids, q), minlevel, atol, by)
    end
    return nothing
end

# The criterion on one subtree's leaf values: `(merges, allmissing)`, in one
# pass that stops at the first value settling the question. Two `Bool`s rather
# than the merged value, which would box a `Union{Missing,T}` per candidate.
# All-missing merges to `missing`; mixed never merges.
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
        # `min`/`max` so a `NaN` propagates into both ends and refuses the
        # merge, where `v < lo` would ignore it.
        lo, hi = min(lo, v), max(hi, v)
        _within(lo, hi, atol) || return false, false
    end
    # Settles the one-element group, which never enters the loop.
    _within(lo, hi, atol) || return false, false
    return true, false
end

# An integer span wider than its own type wraps negative, so the sign test
# refuses that merge instead of reading the wrapped difference as small. `NaN`
# fails it too, as it failed the bare comparison.
function _within(lo, hi, atol)
    d = hi - lo
    return d >= zero(d) && d <= atol
end

# Narrow the accumulator to the eltype actually emitted — an `Int` field that
# never merged comes back a `Vector{Int}`. A concrete accumulator is already
# narrow. The scan tests `typeof(v) === T` rather than `v isa T`: on Julia
# 1.12 `isa` against a `Union{}`-typed local folds to `true`, which would
# answer `Vector{Union{}}`.
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

Present a mixed-level array at one level: every stored value reindexed into
the level-`l` cells it covers. A presentation, not a materialisation — the
result still stores one value per multi-order cell. Methods live in
`src/dimensionaldata.jl`, on arrays carrying a `Cells` axis.
"""
function expand end
