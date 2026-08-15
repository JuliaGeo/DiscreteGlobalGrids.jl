# Aggregation over the hierarchy: `aggregate` reduces leaf data to one fixed
# coarser level; `coarsen` merges complete sibling groups within a tolerance,
# bottom-up, into a `MultiOrderVector` + values — the adaptive-mesh
# constructor. `expand` is the inverse presentation (its DimArray methods live
# in `src/dimensionaldata.jl`; the function is born here).
#
# Contract: docs/design/moc-storage.md §2.
#
# Both verbs read a `(CellVector, values)` PAIR and hand one back, with no cube
# in sight: the compression already lives in this layer and the aggregation
# belongs beside it, so a regridder or a plain-`Array` caller gets both.
#
# Both need `has_sorted_subtrees`. Every statement below rests on a subtree
# being ONE interval of the leaf level: that is what makes a group's values a
# contiguous SLICE of the array rather than a gather, and what turns "is this
# group complete" into two window lookups instead of a descendant walk.

# `Statistics` is not a dependency of this package, and taking one on for a
# three-token reduction would be the wrong trade, so `coarsen`'s default
# summary is spelled out here. It is `Statistics.mean` on everything that has
# one; `by` is the hook for a weighted or robust summary.
_mean(vs) = sum(vs) / length(vs)

# The same refusal `level_ranges` makes, for the same reason and in the same
# words: without descendant ranges a sibling group is not a slice of anything.
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
#
# Both are answered on either window shape, because a bare position is a
# one-cell interval and neither question can tell the shapes apart.

# The first position `cv` holds at or after `lo`, or `nothing`. The prune:
# a subtree with no present leaf is not visited.
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

# The values under a subtree, as a range of positions in the data array, or
# `nothing` when the subtree is not COMPLETE in `cv`.
#
# Both endpoints present and exactly `length(r)` positions between them is the
# whole test: the windows are strictly ascending, so a missing interior
# position would make the count short. That is one binary search per end and no
# dependence on which shape stored the windows.
function _complete_segment(w::CellWindows, r::UnitRange{Int})
    a = windowposition(w, first(r))
    a === nothing && return nothing
    b = windowposition(w, last(r))
    b === nothing && return nothing
    b - a == length(r) - 1 || return nothing
    return a:b
end

# ===========================================================================
# Fixed level
# ===========================================================================

"""
    aggregate(f, cv::CellVector, values::AbstractVector, l::Integer) -> (CellVector, Vector)

Reduce leaf data to level `l`: one output cell per distinct level-`l` ancestor
of `cv`'s cells, ascending, each carrying `f` of the values of its **present**
descendants.

`values` is indexed against `cv` — `values[k]` is the datum of `cv[k]` — and the
answer is another such pair, so a pyramid is this call once per level and no
type of its own:

```julia
cv, data = CellVector(set), read_something(set)
pyramid  = [aggregate(sum, cv, data, l) for l in level(cv)-1:-1:3]
```

# What is aggregated, and what is not

A level-`l` cell appears in the output exactly when `cv` holds at least one of
its level-`level(cv)` descendants, and `f` sees the values of those descendants
— not of the ones `cv` does not hold. PARTIAL GROUPS are therefore reduced over
what is there, the way a raster `aggregate` treats an edge tile, and a group
`cv` misses entirely is absent rather than reduced over nothing.

`f` is handed an `AbstractVector` view into `values`, never a copy, and the
group's values are contiguous in it — a subtree is one interval of the leaf
level, so the cells `cv` holds under one ancestor occupy consecutive positions
of the data array however many windows they are split across.

# Policy stays with the caller

`f` is applied to the values as they stand. `mean` over a group holding a
`missing` answers `missing` because `mean` does; `v -> mean(skipmissing(v))` is
the other reading, and this function will not choose between them.

# Requirements

`l` must be strictly coarser than `level(cv)` — aggregating to a vector's own
level is the identity and is an `ArgumentError` rather than a slow copy — and
the system must have [`has_sorted_subtrees`](@ref).

See also [`coarsen`](@ref), which picks the level per cell instead of fixing
one.
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
    # `_windows` re-checks that the coarse positions ascend. They do wherever a
    # level's canonical order agrees with its subtrees' interval order, which is
    # what `has_sorted_subtrees` asserts — so the check has nothing to find, and
    # a system that broke that agreement gets an error instead of a shuffled
    # answer.
    return CellVector(_windows(positions), coarse, nothing, target), out
end

# One forward pass over the windows. The ancestor at `target` changes only when
# a leaf position leaves the previous ancestor's descendant range, so `ancestor`
# is called once per OUTPUT cell rather than once per leaf, and the group
# boundaries come from `descendant_range` rather than from comparing ids.
#
# The data-array segments are ranges even though a group's leaf positions need
# not be: `k` counts positions `cv` holds, so a gap in the middle of a group is
# a gap in the leaf positions and no gap at all in the concatenation.
function _aggregate_segments(cv::CellVector, coarse::AbstractGrid, target::Int)
    sys = system(cv)
    leaf = level(cv)
    positions = Int[]
    segments = UnitRange{Int}[]
    k = 0        # values consumed
    kstart = 0   # first value of the open group
    stop = 0     # last leaf POSITION of the open group; 0 while none is open
    for (lo, hi) in intervals(cv.windows)
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

Merge each subtree whose leaf values agree to within `atol` into the single
coarse cell that stands for them — the adaptively refined mesh, built from
data at one level.

The result is a mixed-level container and a value per stored cell: cells where
the field is flat are stored once at whatever level flatness reaches, and cells
where it is not stay at `level(cv)` carrying their original values.

```julia
mov, vals = coarsen(cv, temperature; atol = 1.0)      # °C
length(mov) < length(cv)                              # what the tolerance bought
```

# When a cell merges

A cell at a level between `minlevel` and `level(cv)` replaces its subtree iff
both hold:

  - **complete** — `cv` holds *every* level-`level(cv)` descendant of it. An
    incomplete sibling group never merges, and that is what makes the leaf cell
    set exactly recoverable: `CellVector(mov; level = level(cv))` is `cv`
    again, window for window, whatever `atol` was.
  - **within tolerance** — over the LEAF values `vs` under it, either all are
    `missing`, and the merged value is `missing`; or none is, and
    `maximum(vs) - minimum(vs) <= atol`. A group mixing `missing` with data
    never merges, so a coastline stays sharp instead of averaging the ocean
    into the land.

The criterion is monotone — a cell that merges has no descendant that would
not — so the stored cell is the COARSEST merging ancestor of each leaf, and
`minlevel` is where the climb stops.

# The value it stores, and the error that buys

The merged value is `by(vs)` over the LEAF values, never over the children's
summaries: on an equal-area system with the default `by = mean` that is the
exact area-weighted mean of the region, where a mean of means would weight a
sparse branch like a dense one. `by` is handed an `AbstractVector` view of the
leaf values.

The default is this package's own arithmetic mean rather than
`Statistics.mean`, which it equals on everything that has one: `Statistics` is
not a dependency here, and one reduction is not worth taking one on. Pass
`Statistics.mean` explicitly, or any other summary, through `by`.

With the default `by`, the criterion bounds the error outright: every leaf
value is within `atol` of the value stored for the cell covering it, because
`by(vs)` lies between `minimum(vs)` and `maximum(vs)` and those differ by at
most `atol`. Any `by` with that property keeps the bound; `maximum` keeps it,
`sum` does not.

# Requirements

`values` is indexed against `cv`. The system must have
[`has_sorted_subtrees`](@ref), and `minlevel` must be no deeper than
`level(cv)`. The container comes back with `reference_level = level(cv)`, which
is the level its intervals — and the recovery law above — are stated at.

See also [`aggregate`](@ref), which fixes one output level instead of choosing
one per cell, and [`expand`](@ref), the inverse presentation.
"""
function coarsen(cv::CellVector, values::AbstractVector; atol, by=_mean,
        minlevel::Integer=first(levels(system(cv))))
    cells, vals = _coarsen(cv, values; atol, by, minlevel)
    # The only step that knows about the container. `_coarsen` emits cells in
    # ascending descendant-range order already, and the constructor sorts by
    # that same key over disjoint intervals, so the permutation is the identity
    # and `vals` still lines up.
    return MultiOrderVector(system(cv), cells; reference_level=level(cv)), vals
end

"""
    _coarsen(cv, values; atol, by, minlevel) -> (Vector{ID}, Vector)

[`coarsen`](@ref)'s core: the stored cells, ascending by descendant-range start
at `level(cv)`, and their values. Internal — `coarsen` is this plus wrapping
the cells in a [`MultiOrderVector`](@ref).
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
    vals = Any[]
    for c in rootcells(sys)
        _coarsen_visit!(cells, vals, cv, values, c, stop_level, atol, by)
    end
    return cells, _narrow(vals, values)
end

# Top-down, which the monotonicity of the criterion makes equivalent to the
# bottom-up statement and much cheaper: the first merging cell on the way down
# IS the coarsest merging ancestor of every leaf beneath it, so nothing below it
# is ever looked at. A subtree holding no data at all is pruned before its range
# is read, which is what keeps the walk O(#nodes the data touches) rather than
# O(#cells of the level).
#
# `rootcells` and `children` are ascending by contract, so the pre-order walk
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
        # `r` is one position and the prune above proved it present. A leaf
        # keeps its OWN value: `by` summarises a merge, and this is not one —
        # which is what lets `atol = 0` on distinct data be the identity down to
        # the eltype.
        push!(cells, c)
        push!(vals, values[windowposition(w, p)])
        return nothing
    end
    if lc >= minlevel
        ks = _complete_segment(w, r)
        if ks !== nothing
            merges, v = _merged_value(values, ks, atol, by)
            if merges
                push!(cells, c)
                push!(vals, v)
                return nothing
            end
        end
    end
    for child in children(sys, c)
        _coarsen_visit!(cells, vals, cv, values, child, minlevel, atol, by)
    end
    return nothing
end

# The criterion, on one subtree's leaf values. Returns whether it merges and,
# when it does, what to store.
#
# The three-way reading of `missing` is the pinned one: all-missing is a
# perfectly flat region and merges to `missing`, mixed never merges. Counting
# first also keeps `extrema` off a vector that would answer `missing` for a
# reason that has nothing to do with the spread.
function _merged_value(values::AbstractVector, ks::UnitRange{Int}, atol, by)
    vs = view(values, ks)
    nmissing = count(ismissing, vs)
    nmissing == length(vs) && return true, missing
    nmissing == 0 || return false, missing
    lo, hi = extrema(vs)
    hi - lo <= atol || return false, missing
    return true, by(vs)
end

# The values are accumulated untyped because their type is not known until `by`
# has run — a `mean` of integers is not an integer, and an all-`missing` group
# contributes `Missing` whatever the input eltype was. Widening over what was
# actually emitted gives the narrowest container that holds it, rather than
# `promote_type`'s answer to a question nobody asked.
function _narrow(vals::Vector{Any}, values::AbstractVector)
    isempty(vals) && return similar(values, 0)
    T = Union{}
    for v in vals
        v isa T || (T = Union{T,typeof(v)})
    end
    out = Vector{T}(undef, length(vals))
    copyto!(out, vals)
    return out
end

"""
    expand(A, l)

Present a mixed-level array at one level: every stored value reindexed into the
level-`l` cells it covers.

The inverse of [`coarsen`](@ref) as a *presentation* — the answer names
`CellVector(mov; level = l)`'s cells and still stores one value per
multi-order cell, so it is the compression rather than a materialisation of it.

Methods live in `src/dimensionaldata.jl`, on arrays carrying a `Cells` axis;
the function is declared here beside the verbs it inverts.
"""
function expand end
