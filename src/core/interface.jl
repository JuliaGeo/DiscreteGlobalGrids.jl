abstract type AbstractDGGS end

struct NotPortedError <: Exception
    system::Symbol
    operation::Symbol
    detail::String
end

function Base.showerror(io::IO, err::NotPortedError)
    print(io, "DGGS operation `", err.operation, "` is not ported for ",
        err.system, ". ", err.detail)
end

"""
    OrdinalRangeError(system, level, ordinal, total) <: Exception

An ordinal named no cell of its level: [`ordinal_to_cell`](@ref) was handed a
position outside `1:num_cells(system, level)`. Following `NotPortedError`
above, the struct carries only the *facts* and the sentence is built in
[`Base.showerror`](@ref), i.e. when the error is printed, never when it is
thrown — the guard sits on the traversal path that resolves every dense leaf,
so its failure edge stores four immediates and calls `throw`.

| field     | meaning |
|:----------|:--------|
| `system`  | [`system_name`](@ref) of the system whose level was indexed |
| `level`   | the level whose ordinal space was indexed |
| `ordinal` | the offending 1-based position |
| `total`   | `num_cells(system, level)`, so the valid ordinals are `1:total` |

Every wiring of [`ordinal_to_cell`](@ref) raises this one type, so the error a
caller sees does not depend on which system it asked. The native layers
underneath keep their own contracts: `IGeo7.index_to_cell`'s `BoundsError` is
pinned as API shape and still reaches direct native callers, and the kernel
repeats the range check above it rather than relaying it.
"""
struct OrdinalRangeError <: Exception
    system::Symbol
    level::Int
    ordinal::Int
    total::Int
end

function Base.showerror(io::IO, err::OrdinalRangeError)
    print(io, "ordinal ", err.ordinal, " is out of range for ", err.system,
        " level ", err.level, ": an ordinal is a cell's 1-based position among ",
        "the cells of its level in ascending canonical-id order, and level ",
        err.level, " has ", err.total, " of them, so its ordinals are 1:",
        err.total)
end

# ---------------------------------------------------------------------------
# Trait functions
#
# One generic function per fact about a grid system, dispatched on the system
# struct itself. Methods are defined per system in `src/core/systems/*.jl`;
# prose (report section, storage model, sources, references, notes) lives in
# each system struct's docstring.
#
# Every registered system defines a method for every universally-known fact, so
# those functions have no `AbstractDGGS` fallback: an unregistered type is a
# `MethodError`, not a silent default. `root_count` and `radix` are the two
# facts that can be genuinely unverified, so they do have fallbacks that throw
# `NotPortedError`.
# ---------------------------------------------------------------------------

"""
    system_name(system::AbstractDGGS) -> Symbol

Canonical name of the grid system, e.g. `:HEALPix`, `:H3`, `:ISEA4R`.
"""
function system_name end

"""
    grid_family(system::AbstractDGGS) -> Symbol

Family the system belongs to, e.g. `:icosahedral_hex`, `:healpix`, `:ivea`.
"""
function grid_family end

"""
    base_solid(system::AbstractDGGS) -> Symbol

Polyhedron (or polyhedron-like base) the grid is built on, e.g.
`:icosahedron`, `:cube`, `:rhombic_triacontahedron`.
"""
function base_solid end

"""
    cell_shape(system::AbstractDGGS) -> Symbol

Shape of a cell, e.g. `:hexagon`, `:square`, `:pentagon`, `:rhomb`.
"""
function cell_shape end

"""
    is_equal_area(system::AbstractDGGS) -> Bool

Whether cells of a given level all have the same area.
"""
function is_equal_area end

"""
    aperture(system::AbstractDGGS)

Refinement aperture: the area ratio between a parent cell and a child cell.
Usually an `Int` (3, 4, 7, 9); a `Symbol` (`:implementation_defined`,
`:family`) when the system does not pin one, or a `Tuple` for mixed apertures.
"""
function aperture end

"""
    canonical_index_name(system::AbstractDGGS) -> Symbol

Name of the canonical cell-id encoding, e.g. `:h3_index`, `:s2_cellid`,
`:nested`.
"""
function canonical_index_name end

"""
    max_level(system::AbstractDGGS) -> Union{Int,Nothing}

Deepest refinement level the canonical index can represent, or `nothing` when
the system imposes no bound.
"""
function max_level end

"""
    supports_prefix_ranges(system::AbstractDGGS) -> Bool

`true` means a node's descendants at a fixed leaf level are a contiguous range
in the canonical ordinal id space. That is the property needed for sparse
partial trees backed by sorted ids; the kernel's operational counterpart is
[`has_descendant_ranges`](@ref), which a wired system must verify before the
generic cursor prunes with it.
"""
function supports_prefix_ranges end

"""
    root_count(system::AbstractDGGS) -> Int

Number of level-0 (root) zones. Throws [`NotPortedError`](@ref) for systems
whose root layout is not verified yet.
"""
function root_count(system::AbstractDGGS)
    throw(NotPortedError(system_name(system), :root_count,
        "The root zone count is not verified yet."))
end

"""
    radix(system::AbstractDGGS) -> Int

Number of children of a non-root node. Throws [`NotPortedError`](@ref) for
systems without a fixed verified integer radix.
"""
function radix(system::AbstractDGGS)
    throw(NotPortedError(system_name(system), :radix,
        "The child radix is not a fixed verified integer for this system."))
end

# ---------------------------------------------------------------------------
# Generic operations derived from the traits
#
# Everything below is `radix^level` arithmetic on `Int`, which is one
# multiplication away from wrapping silently — `leaf_count(HEALPixDGGS(), 30)`
# used to answer `-4611686018427387904` cells, and level 32 answered `0`. Two
# guards, in this order:
#
#   * a system that pins a `max_level` rejects anything deeper outright. That
#     bound *is* the wrap point (HEALPix's 29 is where `12 * 4^level` stops
#     fitting `Int64`, see `src/core/systems/healpix.jl`), and past it there is
#     no cell to count, so the answer is an error rather than a number.
#   * a system that pins none (`max_level === nothing` — RHEALPix, RTEA, the
#     other `supports_prefix_ranges` registry entries) has no such line to
#     draw, so the products themselves are checked and the caller gets an
#     `OverflowError` instead of a wrapped count.
#
# Both are a handful of integer operations against arithmetic these functions
# were doing anyway.
# ---------------------------------------------------------------------------

# The level bound as one reusable throw site. Returns the level as an `Int`, so
# call sites read `lvl = _checked_level(system, level)`.
function _checked_level(system::AbstractDGGS, level::Integer)
    level >= 0 || throw(ArgumentError("level must be non-negative"))
    limit = max_level(system)
    limit === nothing || level <= limit ||
        throw(ArgumentError("$(system_name(system)) level must be in 0:$limit"))
    return Int(level)
end

# `base^exponent` with every step checked: an unbounded system overflows loudly.
function _checked_pow(base::Integer, exponent::Integer)
    result = 1
    for _ in 1:exponent
        result = Base.checked_mul(result, Int(base))
    end
    return result
end

"""
    leaf_count(system::AbstractDGGS, level::Integer) -> Int

Number of logical id slots at `level`: `root_count * radix^level`. Distinct
from [`num_cells`](@ref), which counts *valid* cells — aperture-7 systems have
pentagon gaps, and only ordinal-id systems make the two coincide.

Levels past [`max_level`](@ref) are an `ArgumentError`; where the system pins
no maximum, an id space too deep for `Int` is an `OverflowError`. Neither ever
returns a wrapped count.
"""
function leaf_count(system::AbstractDGGS, level::Integer)
    lvl = _checked_level(system, level)
    return Base.checked_mul(root_count(system), _checked_pow(radix(system), lvl))
end

"""
    child_ids(system::AbstractDGGS, level::Integer, id::Integer) -> Vector{Int}

Ids of the immediate children of cell `(level, id)`, ascending — the radix
block `id * radix .+ (0:radix-1)`. Only valid for systems whose canonical
ordinal ids are prefix-ranges.

The children live at `level + 1`, so that is the level required to be
representable: asking a system that pins a [`max_level`](@ref) for the children
of a cell *at* its maximum is an `ArgumentError`, not a list of ids the
encoding cannot hold. Systems that pin no maximum get an `OverflowError` where
the block leaves `Int`.
"""
function child_ids(system::AbstractDGGS, level::Integer, id::Integer)
    supports_prefix_ranges(system) || throw(NotPortedError(system_name(system), :child_ids,
        "This system needs implementation-specific child enumeration."))
    lvl = _checked_level(system, level)
    limit = max_level(system)
    limit === nothing || lvl < limit || throw(ArgumentError(
        "children of a $(system_name(system)) cell at level $lvl are at level $(lvl + 1), past max_level $limit"))
    b = radix(system)
    base = Base.checked_mul(Int(id), b)
    return collect(base:Base.checked_add(base, b - 1))
end

function root_child_ids(system::AbstractDGGS)
    return collect(0:(root_count(system) - 1))
end

"""
    leaf_interval(system, level, id, leaf_level) -> UnitRange{Int}

Return the inclusive 0-based leaf ordinal interval owned by node `(level, id)`.
Only valid for systems whose canonical ordinal ids are prefix-ranges.

`leaf_level` past [`max_level`](@ref) is an `ArgumentError` — the span
`radix^(leaf_level - level)` is exactly the product that wraps `Int64` there.
Systems that pin no maximum get an `OverflowError` instead, so no caller ever
receives a wrapped (possibly negative, possibly overlapping a sibling's)
interval.
"""
function leaf_interval(system::AbstractDGGS, level::Integer, id::Integer, leaf_level::Integer)
    supports_prefix_ranges(system) || throw(NotPortedError(system_name(system), :leaf_interval,
        "No verified prefix-range id model for this system."))
    0 <= level <= leaf_level || throw(ArgumentError("expected 0 <= level <= leaf_level"))
    _checked_level(system, leaf_level)
    b = radix(system)
    span = _checked_pow(b, Int(leaf_level) - Int(level))
    lo = Base.checked_mul(Int(id), span)
    return lo:Base.checked_add(lo, span - 1)
end

function cell_polygon(system::AbstractDGGS, level::Integer, id::Integer)
    throw(NotPortedError(system_name(system), :cell_polygon,
        "Port the authoritative cell-boundary math or wrap an implementation."))
end

function cell_extent(system::AbstractDGGS, level::Integer, id::Integer)
    poly = cell_polygon(system, level, id)
    return GI.extent(poly)
end

function full_sphere_extent()
    return GO.UnitSpherical.SphericalCap(GO.UnitSphericalPoint(0.0, 0.0, 1.0), nextfloat(Float64(pi)))
end

# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

"""
    all_systems() -> Tuple

Every DGGS registered in this package, as instances, in report order. Report
section 1.8 splits into two systems (`ISEA4RDGGS` and `ISEA9RDGGS`), so the 13
report sections give 14 systems.
"""
function all_systems()
    return (
        H3DGGS(), S2DGGS(), A5DGGS(), IGEO7DGGS(),
        ISEA3HDGGS(), ISEA4HDGGS(), ISEA4TDGGS(),
        ISEA4RDGGS(), ISEA9RDGGS(), IVEADGGS(:IVEA_family),
        RTEADGGS(:RTEA_family), RHEALPixDGGS(), AusPIXDGGS(),
        HEALPixDGGS(),
    )
end
