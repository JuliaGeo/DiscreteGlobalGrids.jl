# --- the collection surface ------------------------------------------------

Base.length(set::MultiOrderCellSet) = length(set.cells)
Base.eltype(::Type{MultiOrderCellSet{S,ID}}) where {S,ID} = ID
Base.eltype(set::MultiOrderCellSet) = eltype(typeof(set))
Base.isempty(set::MultiOrderCellSet) = isempty(set.cells)
Base.iterate(set::MultiOrderCellSet, state...) = iterate(set.cells, state...)
Base.getindex(set::MultiOrderCellSet, i::Int) = set.cells[i]
Base.firstindex(::MultiOrderCellSet) = 1
Base.lastindex(set::MultiOrderCellSet) = length(set.cells)
Base.eachindex(set::MultiOrderCellSet) = Base.OneTo(length(set.cells))
Base.collect(set::MultiOrderCellSet) = copy(set.cells)

"""
    system(set::MultiOrderCellSet)

The system the set's cells are named in.
"""
system(set::MultiOrderCellSet) = set.system

# --- which cells fit inside the target -------------------------------------

"""
    iscontained(set::MultiOrderCellSet, i::Integer) -> Bool

Whether the set's `i`th cell was **proven** to lie inside the coverage target.
`true` means `Within` was asked of it and held. `false` means one of two things,
told apart by the cell's level:

  - above the traversal's maximum depth, the cell was asked and the target's
    boundary crosses it. The flag is exact there, in both directions.
  - at the maximum depth, the cell was never asked. The traversal ran out of
    depth and emitted it to cover; it may or may not fit.

In `level` mode the maximum depth *is* the reference level, so the blind spot
sits exactly there. In `maxcells` mode they part: the reference level is the
deepest level the budget reached, the maximum depth is the `maxlevel` cap, and a
budget that stopped short of the cap — the ordinary case — carries an exact flag
on every member.

The asymmetry is the contract. `Within` allocates ~48 KB per call against ~600 B
for `Intersects`, and a deep coverage ends on thousands of reference-level
cells: labelling them exactly costs hundreds of megabytes and changes no member
of the set. The direction that matters survives — `true` implies inside — and
code needing exact containment at the reference level asks `Within` of the few
cells it cares about.

`argmin(level, set)` is therefore not "the coarsest cell inside the target":
every emission can be unproven. [`coarsest_contained`](@ref) reads this flag.
"""
iscontained(set::MultiOrderCellSet, i::Integer) = set.contained[i]

"""
    coarsest_contained(set::MultiOrderCellSet) -> cell id or `nothing`

The shallowest cell of `set` **proven** inside the coverage target, or `nothing`
when no cell of it was — see [`iscontained`](@ref) for what "proven" leaves
out. Maximum-depth cells are never tested, so a set of nothing but those answers
`nothing` even where some fit; a target smaller than one cell is the clearest
way there, not the only one.

A budget set answers `nothing` for a second, more ordinary reason: at ten cells
over a state nothing has been refined far enough to fit inside it, and the
accessor says so rather than hand back the shallowest crossing cell. Raise the
budget and the answer appears.

```julia
set = query(sys, MultiOrderCoverage(tile); level = 12)
cell = coarsest_contained(set)          # `nothing`, or a cell that fits in `tile`
```

Ties go to the first such cell in the set's own order.
"""
function coarsest_contained(set::MultiOrderCellSet)
    best = nothing
    for i in eachindex(set)
        set.contained[i] || continue
        (best === nothing || level(set.cells[i]) < level(best)) && (best = set.cells[i])
    end
    return best
end

# --- geometry, without a level grid per cell -------------------------------

# A set is not a grid — no indices, and its cells are at different levels —
# but it does know which level grid each cell belongs to. `levelgrid` is O(1),
# so nothing here is worth caching.

"""
    cell_boundary(set::MultiOrderCellSet, c) -> Vector{UnitSphericalPoint}
    cell_centroid(set::MultiOrderCellSet, c) -> UnitSphericalPoint
    cell_polygon(set::MultiOrderCellSet, c) -> GI.Polygon

Geometry of one member of a mixed-level set, read from
`levelgrid(system(set), level(c))`. Same values as that grid gives, without the
caller resolving it.
"""
cell_boundary(set::MultiOrderCellSet, c::AbstractCellIndex) =
    cell_boundary(levelgrid(set.system, level(c)), c)

cell_centroid(set::MultiOrderCellSet, c::AbstractCellIndex) =
    cell_centroid(levelgrid(set.system, level(c)), c)

cell_polygon(set::MultiOrderCellSet, c::AbstractCellIndex) =
    cell_polygon(levelgrid(set.system, level(c)), c)

"""
    cell_polygons(set::MultiOrderCellSet) -> Vector{<:GI.Polygon}

Every cell of the set as a unit-sphere polygon, in the set's own order: what a
plot of a coverage needs, in one call.

```julia
poly(GO.transform(GO.GeographicFromUnitSphere(), cell_polygons(set)))
```
"""
cell_polygons(set::MultiOrderCellSet) =
    [cell_polygon(levelgrid(set.system, level(c)), c) for c in set.cells]

"""
    curve_keys(set::MultiOrderCellSet) -> Vector{Int}

Return stored-cell sort keys. For sorted-subtree systems, each key is the start
of the cell's reference-level descendant range. Otherwise it is the cell's
index within its own level and is not comparable across levels.
"""
curve_keys(set::MultiOrderCellSet) = set.keys

function Base.show(io::IO, set::MultiOrderCellSet)
    print(io, "MultiOrderCellSet(", typeof(set.system).name.name, ", ",
        length(set.cells), " cells")
    isempty(set.cells) || print(io, ", levels ",
        minimum(level, set.cells), ":", maximum(level, set.cells))
    print(io, ")")
end

Base.show(io::IO, ::MIME"text/plain", set::MultiOrderCellSet) = show(io, set)

"""
    level_ranges(set::MultiOrderCellSet, l::Integer) -> Vector{UnitRange{Int}}

Expand the set to sorted, disjoint index ranges in `levelgrid(sys, l)`,
merging adjacent ranges. Requires sorted subtrees and `l` no shallower than any
cell in the set.

!!! warning "Two things the expansion is not"
    Not universal: it throws where [`has_sorted_subtrees`](@ref) is `false`
    (A5). Branch on the trait, or expand with `descendants(sys, c, l)`.

    Not a covering. A cell is in the set because the *cell* is inside the
    target; under non-congruent refinement its descendants need not be, so the
    expansion can name leaves the target does not touch — most visibly inside a
    hole. See [`MultiOrderCoverage`](@ref).
"""
function level_ranges(set::MultiOrderCellSet, l::Integer)
    has_sorted_subtrees(set.system) || throw(ArgumentError(
        "$(typeof(set.system)) has no descendant ranges, so a multi-order set " *
        "cannot be expanded to index ranges"))
    target = Int(l)
    out = UnitRange{Int}[]
    for c in set.cells
        level(c) <= target || throw(ArgumentError(
            "cannot expand to level $target: the set contains a level-$(level(c)) cell"))
        r = descendant_range(set.system, c, target)
        if !isempty(out) && first(r) == last(out[end]) + 1
            out[end] = first(out[end]):last(r)
        else
            push!(out, r)
        end
    end
    return out
end

"""
    cellindices(set::MultiOrderCellSet, l::Integer) -> Vector{<:AbstractCellIndex}

The set expanded to level `l` as typed ids, ascending — [`level_ranges`](@ref)
resolved through `cellindex`. O(cells at `l`), so reach for the ranges instead
wherever the indices are what is wanted.
"""
function cellindices(set::MultiOrderCellSet, l::Integer)
    grid = levelgrid(set.system, Int(l))
    out = cellindextype(set.system)[]
    for r in level_ranges(set, l), i in r
        push!(out, cellindex(grid, i))
    end
    return out
end
