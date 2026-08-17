# The plain-data layer the convention code operates on. Nothing here knows
# about Zarr: the Zarr extension builds a `StoreSnapshot` from a live store and
# the test suite builds one by hand, so convention logic stays testable without
# a store.
#
# Include order: `errors.jl` first, then this, then `conventions.jl`. The
# `CellEncoding` abstract type of `encodings.jl` must be bound in the including
# module before this file is read: `StoreDescription` types its encoding field
# with it.

# ===========================================================================
# Snapshots
# ===========================================================================

"""
    ArrayEntry(; name, attrs, shape, eltype = Any, dims = String[])

One array of a [`StoreSnapshot`](@ref), as metadata only: no chunks, no
compressor, no values.

  - `name`: the array's name within its group.
  - `attrs`: the array's attributes, mutable, as read (`.zattrs` in Zarr v2,
    `attributes` in v3). Conventions stamp write attrs into this dictionary.
  - `shape`: the DECLARED shape. Chunk extents may exceed it; every length
    check is against this.
  - `eltype`: the element type, or `Any` where it is not known.
  - `dims`: the dimension names, outermost first, empty for a scalar array.
    In Zarr v2 these come from `_ARRAY_DIMENSIONS`.
"""
Base.@kwdef struct ArrayEntry
    name::String
    attrs::Dict{String,Any} = Dict{String,Any}()
    shape::Tuple{Vararg{Int}} = ()
    eltype::Any = Any
    dims::Vector{String} = String[]
end

Base.ndims(a::ArrayEntry) = length(a.shape)

Base.copy(a::ArrayEntry) = ArrayEntry(name=a.name, attrs=deepcopy(a.attrs),
    shape=a.shape, eltype=a.eltype, dims=copy(a.dims))

"""
    StoreSnapshot(; identifier = "", attrs = Dict{String,Any}(), arrays = ArrayEntry[])

A metadata-only view of one store group: its group attributes and its arrays.

This is the only thing a [`DGGSConvention`](@ref) ever sees. `identifier` is
the store URL or path, used in [`DGGSFormatError`](@ref) messages.

Both `attrs` and each array's `attrs` are mutable; [`encode!`](@ref) writes
into them. `copy` deep-copies the attribute dictionaries so a snapshot can be
re-stamped without disturbing the one it came from.

Subgroups are out of scope: a snapshot describes one group, and a hierarchy is
a sequence of snapshots.
"""
Base.@kwdef struct StoreSnapshot
    identifier::String = ""
    attrs::Dict{String,Any} = Dict{String,Any}()
    arrays::Vector{ArrayEntry} = ArrayEntry[]
end

Base.copy(s::StoreSnapshot) = StoreSnapshot(identifier=s.identifier,
    attrs=deepcopy(s.attrs), arrays=[copy(a) for a in s.arrays])

"""
    getarray(snapshot, name) -> Union{ArrayEntry, Nothing}

The array of `snapshot` called `name`, or `nothing`.
"""
function getarray(s::StoreSnapshot, name::AbstractString)
    for a in s.arrays
        a.name == name && return a
    end
    return nothing
end

"""
    arrays_on(snapshot, dim; exclude = ()) -> Vector{String}

Names of the arrays whose dimensions include `dim`, in snapshot order, minus
`exclude`.
"""
function arrays_on(s::StoreSnapshot, dim::AbstractString; exclude=())
    return [a.name for a in s.arrays if dim in a.dims && !(a.name in exclude)]
end

# ===========================================================================
# Grid parameters
# ===========================================================================

"""
    Ellipsoid(; name, semi_major_axis, semi_minor_axis, inverse_flattening, radius)

The reference ellipsoid a store declares, as read. Every field may be
`nothing`: stores supply whichever subset their producer had.

Read accepts both spellings of the axis keys — the schema-valid
`semi_major_axis`/`semi_minor_axis` and the `semimajor_axis`/`semiminor_axis`
form that xdggs and every store in the wild actually use — and normalizes to
these field names. [`ellipsoid_attrs`](@ref) writes the schema-valid spelling.
"""
Base.@kwdef struct Ellipsoid
    name::Union{String,Nothing} = nothing
    semi_major_axis::Union{Float64,Nothing} = nothing
    semi_minor_axis::Union{Float64,Nothing} = nothing
    inverse_flattening::Union{Float64,Nothing} = nothing
    radius::Union{Float64,Nothing} = nothing
end

"""
    DEFAULT_ELLIPSOID

The sphere of radius 6 370 997 m that `zarr-conventions/dggs` prescribes when a
store declares no ellipsoid. [`ellipsoid`](@ref) applies it; a decoded
[`StoreDescription`](@ref) keeps `nothing` so that a store's silence is never
confused with a declaration.
"""
const DEFAULT_ELLIPSOID = Ellipsoid(name="sphere", radius=6370997.0)

"""
    GridOrientation(; vert0_lon, vert0_lat, vert0_azimuth, rotation_pattern)

The icosahedron placement of an ISEA-family grid, as read. Two stores of the
same grid name and level but different orientations hold incomparable cell ids
and nothing in any convention flags it, so this travels with the description.

Fields carry the DGGRID vertex-0 parameters; `nothing` means the store said
nothing. Convention A spells them `dggs_vert0_lon` and friends inside the
`dggs` object, convention B spells the longitude `igeo7_dggs_vert0_lon` on the
coordinate.
"""
Base.@kwdef struct GridOrientation
    vert0_lon::Union{Float64,Nothing} = nothing
    vert0_lat::Union{Float64,Nothing} = nothing
    vert0_azimuth::Union{Float64,Nothing} = nothing
    rotation_pattern::Union{String,Nothing} = nothing
end

# Field-wise `==`: the default for an immutable struct compares fields with
# `===`, which is false for two equal `String`s parsed from different stores.
for T in (:Ellipsoid, :GridOrientation)
    @eval function Base.:(==)(a::$T, b::$T)
        return all(f -> getfield(a, f) == getfield(b, f), fieldnames($T))
    end
end

# ===========================================================================
# The description
# ===========================================================================

"""
    StoreDescription(; gridname, kwargs...)

The pivot between store attributes and a cell axis. Reading is
attrs → description → lookup; writing is lookup → description → attrs.

  - `gridname`: canonical grid name, e.g. `"igeo7"` or `"healpix"`.
  - `system`: the grid system it resolves to, from [`GRID_REFERENCE`](@ref).
  - `idscheme`: how a cell id is packed — `:z7int`, `:nested`, `:ring`, …
  - `level`: the refinement level, or `nothing` for a variable-sized axis.
  - `encoding`: a `CellEncoding` instance, from `ENCODING_REGISTRY`.
  - `coordinate`: name of the array encoding the cell ids, `nothing` when the
    axis is implicit in position.
  - `spatial_dimension`: the dimension the data variables share. It is not the
    coordinate's own dimension: a ranges coordinate is shaped `(n, 2)` and
    names neither.
  - `variables`: the data variable names, or `nothing` when unknown.
  - `ellipsoid`, `orientation`, `geodetic_conversion`: grid parameters as read.
  - `provenance`: verbatim source attributes keyed by convention name, enough
    to regenerate a value-identical store.

Every field but `gridname` and `provenance` may be `nothing`, meaning the store
did not say. `==` compares the semantic fields and deliberately ignores
`provenance`, which differs between a store and its faithful copy.
"""
Base.@kwdef struct StoreDescription
    gridname::String
    system::Any = nothing
    idscheme::Union{Symbol,Nothing} = nothing
    level::Union{Int,Nothing} = nothing
    encoding::Union{CellEncoding,Nothing} = nothing
    coordinate::Union{String,Nothing} = nothing
    spatial_dimension::Union{String,Nothing} = nothing
    variables::Union{Vector{String},Nothing} = nothing
    ellipsoid::Union{Ellipsoid,Nothing} = nothing
    orientation::Union{GridOrientation,Nothing} = nothing
    geodetic_conversion::Union{Bool,Nothing} = nothing
    provenance::Dict{String,Any} = Dict{String,Any}()
end

const DESCRIPTION_FIELDS = filter(!=(:provenance), fieldnames(StoreDescription))

function Base.:(==)(a::StoreDescription, b::StoreDescription)
    return all(f -> getfield(a, f) == getfield(b, f), DESCRIPTION_FIELDS)
end

"""
    ellipsoid(desc) -> Ellipsoid

The description's ellipsoid, or [`DEFAULT_ELLIPSOID`](@ref) where it declared
none.
"""
ellipsoid(d::StoreDescription) = d.ellipsoid === nothing ? DEFAULT_ELLIPSOID : d.ellipsoid

function Base.show(io::IO, d::StoreDescription)
    print(io, "StoreDescription(", d.gridname)
    d.idscheme === nothing || print(io, "/", d.idscheme)
    d.level === nothing || print(io, ", level ", d.level)
    d.encoding === nothing || print(io, ", ", nameof(typeof(d.encoding)))
    d.coordinate === nothing || print(io, ", coordinate ", repr(d.coordinate))
    d.variables === nothing || print(io, ", ", length(d.variables), " variables")
    print(io, ")")
end

# ===========================================================================
# Detection
# ===========================================================================

"""
    Detection(convention, rank, evidence; payload = Dict{Symbol,Any}())

What [`detect`](@ref) returns when a convention recognizes a snapshot.

`rank` is `:declared` when the store names the convention by identifier — the
`zarr_conventions` UUID — and `:fingerprint` when it was recognized by the
shape of its attributes. Declared outranks fingerprint in registry order, so a
store that says what it is is decoded by that convention first.

`evidence` is a short human phrase naming what matched; `payload` carries
whatever the convention's own `decode` needs, typically the array it found.
"""
Base.@kwdef struct Detection
    convention::Any
    rank::Symbol = :fingerprint
    evidence::String = ""
    payload::Dict{Symbol,Any} = Dict{Symbol,Any}()
end

Detection(c, rank::Symbol, evidence::AbstractString; payload=Dict{Symbol,Any}()) =
    Detection(convention=c, rank=rank, evidence=String(evidence), payload=payload)

_rankorder(r::Symbol) = r === :declared ? 1 : 2

# `DGGSFormatError`, which the docstrings above refer to, is defined in
# `errors.jl`: every layer of `src/io/` throws it, so it is included before all
# of them.
