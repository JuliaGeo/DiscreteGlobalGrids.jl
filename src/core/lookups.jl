# ---------------------------------------------------------------------------
# The DGGS lookup supertype
#
# `H3Lookup`, `A5Lookup`, `IGeo7Lookup` and `HealpixLookup` each subtyped
# `Lookups.Lookup{ID,1}` directly, with no common supertype, so nothing generic
# could dispatch on "a DGGS lookup" — every shared behaviour had to be written
# once per system and kept in sync by hand, as the four one-line `treeify`
# delegations at the bottom of the `<X>Kernel.jl` files were. One abstract type
# here is what lets the globe-lookup machinery (`globe_ids.jl`), the `treeify`
# wiring below and the selectors built on it be written once.
#
# The type lives in core rather than in `Helpers` because it needs
# `DimensionalData`, and `Helpers` is a deliberately dependency-free leaf
# module. The lookup modules are grandchildren of this one, so they reach it
# the same way the kernel wirings do: `import ...DiscreteGlobalGrids as DGG`.
# ---------------------------------------------------------------------------

"""
    AbstractDGGSLookup{ID} <: DimensionalData.Lookups.Lookup{ID,1}

Supertype of this package's `DimensionalData` lookups — `H3Lookup`,
`A5Lookup`, `IGeo7Lookup`, `HealpixLookup`.

A DGGS lookup is a one-dimensional lookup over cell ids of type `ID`, held
strictly ascending in canonical-id order (which is what every selector's binary
search, and the leaf ordering of the trees in `grid_types.jl`, assume). The ids
themselves live in the lookup's `data` field as any `AbstractVector{ID}` — an
ordinary `Vector` for partial coverage, a [`DGGSGlobeIds`](@ref) for a
globe-complete dimension.

Systems still own their own constructors, metadata conventions and native
geometry; what this supertype exists for is the code that is the same for all
of them.
"""
abstract type AbstractDGGSLookup{ID} <: DD.Lookups.Lookup{ID,1} end

# --------------------------------------------------------------------------
# ConservativeRegridding.Trees
#
# Treeifying a lookup directly is the shortest path from a `DimensionalData`
# dimension to a `Regridder`. Which grid it goes through is decided by the
# *type* of the id vector: a `DGGSGlobeIds` is every cell of its level, so the
# tree is the dense `DGGSGrid` cursor, and anything else is a stored subset, so
# it is the `DGGSPartialGrid` cursor.
#
# By typing rather than by sniffing `length(l.data) == num_cells(...)`, because
# that test is not merely inelegant but wrong: a lookup built with
# `validate=false` can hold an invalid interior id and still have
# globe-complete length, so the sniff would silently produce a complete-globe
# tree instead of throwing at the bad id.
#
# Alignment — leaf index `i` is lookup position `i`, which is what makes a
# `Regridder` line up with a `DimArray` over this dimension without a
# permutation — survives both branches, but by two different arguments. The
# partial cursor takes the lookup's id vector by reference, so leaf `i` is
# `l.data[i]` because it is literally the same vector (`grid.ids === l.data`,
# see `DGGSPartialGrid`); the dense cursor numbers its leaves by
# ordinal, and a `DGGSGlobeIds` position *is* an ordinal by definition. Two
# arguments, hence two tests.
#
# One method here rather than the four per-system delegations it replaces
# (which is what `AbstractDGGSLookup` exists for): those were four copies of
# one line, and the globe branch is exactly the kind of thing that has to be
# present in every copy — a system that missed it would silently take the
# O(globe) partial path rather than fail.
# --------------------------------------------------------------------------

Trees.treeify(m::GO.Spherical, l::AbstractDGGSLookup) = _treeify_ids(m, l, DD.parent(l))

# The globe branch cannot be written as a signature on the lookup: the only
# type parameter `AbstractDGGSLookup` has is the *element* type, so
# `AbstractDGGSLookup{<:DGGSGlobeIds}` would name a lookup whose cell ids are
# themselves globe-id vectors, not a lookup over a globe. The distinction lives
# in each concrete lookup's own `data` type parameter, so dispatching on the id
# vector is what lets one method serve all four.
_treeify_ids(m::GO.Spherical, l::AbstractDGGSLookup, ::AbstractVector) =
    Trees.treeify(m, DGGSPartialGrid(l))

_treeify_ids(m::GO.Spherical, ::AbstractDGGSLookup, ids::DGGSGlobeIds) =
    Trees.treeify(m, DGGSGrid(ids.system, ids.level))

# --------------------------------------------------------------------------
# Point selectors
#
# `At(id)` and `Contains(point)` both end at one question — which position of
# this lookup holds this cell id — and every `<X>Lookup` answers it by binary
# searching its stored ids. On a globe that search is not merely O(log N) but
# O(log N) *kernel calls*: every probe of a `DGGSGlobeIds` is an
# `ordinal_to_cell`, some 49 of them at H3 res 15, spent locating a cell whose
# position the ordinal contract already names outright. `cell_to_ordinal` *is*
# that position — element `i` of a `DGGSGlobeIds` is `ordinal_to_cell(system,
# level, i)`, so the position of `id` is `cell_to_ordinal(system, level, id)` —
# and it is one call.
#
# Partial lookups keep the binary search. O(log N) over an arbitrary sorted
# subset is where "O(1) everywhere" stops, and stops deliberately: the contents
# of a stored subset are data, not arithmetic.
#
# Only the id -> position step is generic. Point -> id stays where it already
# is, in each system's own `Contains` (`H3Native.lonlat_to_cell` and its
# siblings): there is no generic point->cell kernel operation, and none is
# needed here, since the step that was O(log N) is this one.
#
# Dispatch is on the id vector rather than the lookup, for the reason spelled
# out above `_treeify_ids`: the globe case cannot be named on the supertype.
# --------------------------------------------------------------------------

"""
    cell_position(ids, id) -> Union{Nothing,Int}

Position of cell `id` in a lookup's ascending id vector, or `nothing` when the
vector does not hold it — the one question `At` and `Contains(point)` ask of
every `<X>Lookup`, and hence the whole of what those selectors need to
specialize on.

Over stored ids this is a binary search. Over a [`DGGSGlobeIds`](@ref) it is
[`cell_to_ordinal`](@ref): O(1), and *exact* — an id that names no cell of the
globe's level is rejected by the round trip back through
[`ordinal_to_cell`](@ref), with no `is_valid_cell` call and no O(N) pass. That
is a stronger guarantee than the stored branch can make, where an id sitting
between two stored ids is only as impossible as the lookup's `validate` flag
made it.
"""
function cell_position(ids::AbstractVector, id)
    i = Helpers.sorted_index(ids, id)
    return iszero(i) ? nothing : i
end

function cell_position(ids::DGGSGlobeIds, id)
    system = ids.system
    level = ids.level
    # `cell_to_ordinal` is contracted only over ids that *are* cells at
    # `level`, and off that domain each wiring answers differently: H3 reads
    # the base-cell field of whatever it is handed and returns a perfectly
    # ordinary ordinal — or a `BoundsError`, once those bits exceed 121 — while
    # A5 and HEALPix throw `ArgumentError` and IGEO7 `InvalidZ7Error`. Every one
    # of those means one thing at this call site, which is a membership test:
    # the globe does not hold this id. So the failure is answered rather than
    # propagated, exactly as the branch above answers `nothing` instead of
    # throwing, and the selector that called this reports the absence in its own
    # words — the same sentence a partial lookup gives for the same id.
    ordinal = try
        cell_to_ordinal(system, level, id)
    catch
        return nothing
    end
    # The two O(1) checks that make the answer exact. The range is the one
    # `ordinal_to_cell` reports as `OrdinalRangeError`, tested here rather than
    # relayed from there because out of range is not a caller error at this call
    # site: the caller named an id, never an ordinal, so the answer it is owed
    # is "no such cell". The round trip then rejects everything that decodes
    # *into* range and still is not the cell standing there — ids of the wrong
    # resolution, and structural ids whose unused slots carry garbage.
    1 <= ordinal <= length(ids) || return nothing
    return ordinal_to_cell(system, level, ordinal) == id ? ordinal : nothing
end
