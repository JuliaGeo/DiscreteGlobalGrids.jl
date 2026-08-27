# ---------------------------------------------------------------------------
# Deprecations.
#
# `cellposition` became two verbs. The old name always meant "where does this
# collection keep the cell", which is `localindex` — on a complete grid that is
# also the `globalindex`, so the forward is exact rather than a best guess.
#
# Kept working rather than removed: several branches were open against the old
# name when it changed, and a rename they can absorb on their own schedule costs
# one method each.
# ---------------------------------------------------------------------------

"""
    cellposition(collection, c) -> Union{Int,Nothing}

Deprecated. Use [`localindex`](@ref) for the index in a collection's own
storage, or [`globalindex`](@ref) for the index in the complete grid at that
level.

The old name did not say which of the two it meant — it answered the local one,
and on a complete grid that is also the global one, which is why the confusion
was survivable for as long as it was. This forwards to [`localindex`](@ref), so
existing calls keep their old behaviour exactly.
"""
function cellposition end

@deprecate cellposition(collection, c::AbstractCellIndex) localindex(collection, c) false
@deprecate cellposition(collection, p::GO.UnitSphericalPoint) localindex(collection, p) false
@deprecate cellposition(collection, lon::Real, lat::Real) localindex(collection, lon, lat) false

# The one-argument handle form: a handle only ever carried a local index.
@deprecate cellposition(h::Engine.SubsetIndexedCell) localindex(h) false
