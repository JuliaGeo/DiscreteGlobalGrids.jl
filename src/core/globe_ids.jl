# ---------------------------------------------------------------------------
# Lazy globe-complete id vectors
#
# Every `<X>Lookup` in this package stores an explicit sorted vector of cell
# ids, so a globe-complete `DimensionalData` dimension costs one id per cell —
# 790 MB at H3 res 7, and simply unconstructible at res 9 (4.8e9 ids, ~39 GB).
# The contents of that vector are, by construction, computable arithmetic: the
# ordinal contract (`kernel.jl`, "Dense ordinals") already reads "the 1-based
# position of a cell among all valid cells of its level, in ascending canonical
# id order", so
#
#     i -> ordinal_to_cell(system, level, i)
#
# *is* the sorted vector of every cell at `level`, for every kernel-wired
# system, with no new per-system proof obligation. `DGGSGlobeIds` is that
# function seen as an `AbstractVector`.
#
# It is deliberately not a new lookup family. The four lookup structs already
# type their id field `A<:AbstractVector{...}` rather than `Vector`, so
# `H3Lookup{<:DGGSGlobeIds}` is an ordinary `H3Lookup` — same selectors, same
# `show`, same `rebuild` — and the whole globe case costs one core type instead
# of four per-system clones of the `DimensionalData` interface. See
# `docs/design/full_globe_lookups.md`.
# ---------------------------------------------------------------------------

"""
    DGGSGlobeIds(system, level) -> AbstractVector{cell_id_type(system)}

Every cell of `system` at `level`, in ascending canonical-id order, computed on
demand instead of stored: element `i` is [`ordinal_to_cell(system, level,
i)`](@ref) and the whole vector is two words wide however many cells it names.

This is what makes a globe-complete `<X>Lookup` — and hence a globe-complete
`DimArray` dimension — constructible at all:

```julia
H3Lookup(DGGSGlobeIds(H3DGGS(), 9))     # 4.8e9 cells, O(1) to build
```

Each lookup's own keyword constructor picks the resolution/level up from the
ids object, so the id vector is the only thing named.

Indexing with anything but a single position falls through to `AbstractArray`'s
default and materializes a `Vector` — `g[1:10]` is ten real ids — which is
exactly how a globe lookup degrades to an ordinary partial one under slicing.

Construction validates only what it can without a wired kernel, as
[`DGGSGrid`](@ref) does: `level` against [`max_level`](@ref). The cell count
comes from `num_cells`, so a level whose id encoding exists but whose grid does
not (`DGGSGlobeIds(A5DGGS(), 30)`, A5's res-30 gap) constructs here and throws
at the first operation that asks for `size`.
"""
struct DGGSGlobeIds{S<:AbstractDGGS,ID<:Integer} <: AbstractVector{ID}
    system::S
    level::Int
    # `ID` is not carried by any field, so the only way in is a constructor
    # that reads it off the system — a hand-written `DGGSGlobeIds{H3DGGS,Int8}`
    # would otherwise claim an element type `ordinal_to_cell` never returns.
    function DGGSGlobeIds(system::AbstractDGGS, level::Integer)
        lvl = Int(level)
        lvl >= 0 || throw(ArgumentError("level must be non-negative"))
        limit = max_level(system)
        limit === nothing || lvl <= limit ||
            throw(ArgumentError("$(system_name(system)) level must be in 0:$limit"))
        return new{typeof(system),cell_id_type(system)}(system, lvl)
    end
end

Base.size(g::DGGSGlobeIds) = (Int(num_cells(g.system, g.level)),)
Base.IndexStyle(::Type{<:DGGSGlobeIds}) = Base.IndexLinear()
Base.getindex(g::DGGSGlobeIds, i::Int) = ordinal_to_cell(g.system, g.level, i)

# `AbstractArray`'s `show` methods print elements, which on a globe is precisely
# the materialization this type exists to avoid: `Base.print_matrix` allocates
# one row range for the whole array unless the caller set `:limit`, so
# `sprint(show, MIME"text/plain"(), l)` on a res-15 H3 globe lookup is an
# `OutOfMemoryError` rather than a slow print. `DimensionalData`'s
# non-compact `show(::Lookup)` takes exactly that path (`Lookups/show.jl`,
# "wrapping: "), so a globe lookup is unprintable without these. The system and
# the level are the entire content anyway — every id follows from them.
Base.show(io::IO, g::DGGSGlobeIds) =
    print(io, "DGGSGlobeIds(", g.system, ", ", g.level, ")")
Base.show(io::IO, ::MIME"text/plain", g::DGGSGlobeIds) = show(io, g)

# One of the two constructor fast paths every `<X>Lookup` needs (design §1.2).
# The lookups check their ids with `Helpers.strictly_increasing`, an O(N) pass —
# 4.8e9 `ordinal_to_cell` calls on a res-9 globe — to establish exactly what the
# ordinal contract already guarantees: `cell_to_ordinal` is *required* to be
# strictly monotone in the id (`kernel.jl`, "Dense ordinals"), so its inverse
# enumerated over `1:num_cells` is strictly ascending by definition. Stating
# that once here covers all four lookups, since they share this function.
Helpers.strictly_increasing(::DGGSGlobeIds) = true
