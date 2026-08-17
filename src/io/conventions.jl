# Pluggable storage conventions: attrs -> `StoreDescription` on read, and
# `StoreDescription` -> attrs on write. A convention never touches a lookup or
# an array; it sees only a `StoreSnapshot`.
#
# `IGeo7System`, `HEALPixSystem`, `ENCODING_REGISTRY` and `encodingname` are
# referenced unqualified. The package binds all four before including this
# file; the standalone suite supplies them in its own module.
#
# Include order: `errors.jl`, then `description.jl`, then this file.

"""
    DGGSConvention

A storage convention: a way a store says which grid, level and encoding its
cell axis has.

Three verbs, all operating on a [`StoreSnapshot`](@ref):

    detect(c, snapshot)             -> Union{Detection, Nothing}
    decode(c, snapshot, detection)  -> StoreDescription
    encode!(c, snapshot, desc)      -> snapshot

Downstream packages subtype this and `push!` onto [`CONVENTION_REGISTRY`](@ref);
[`gridname`](@ref) is the hook for folding a private grid vocabulary in.
"""
abstract type DGGSConvention end

"""
    detect(c::DGGSConvention, snapshot) -> Union{Detection, Nothing}

Whether `c` recognizes `snapshot`, and on what evidence. `nothing` means "not
mine"; a malformed declaration of `c`'s own convention is a
[`DGGSFormatError`](@ref) rather than a silent `nothing`.
"""
function detect end

"""
    decode(c::DGGSConvention, snapshot, detection) -> StoreDescription

The part of the store's identity `c` can read. Fields `c` has no vocabulary for
are left `nothing` so another convention can supply them; see
[`describe_store`](@ref) for how several decodings are reconciled.
"""
function decode end

"""
    encode!(c::DGGSConvention, snapshot, desc) -> snapshot

Stamp `desc` into `snapshot`'s attribute dictionaries in `c`'s spelling.
Read-only conventions do not implement this.
"""
function encode! end

"""
    conventionname(c::DGGSConvention) -> String

Short name used in [`DGGSFormatError`](@ref) messages and as the provenance key.
"""
conventionname(c::DGGSConvention) = string(nameof(typeof(c)))

encode!(c::DGGSConvention, ::StoreSnapshot, ::StoreDescription) =
    throw(ArgumentError("$(conventionname(c)) is a read-only convention: it has no encode! method"))

# ===========================================================================
# The grid reference table
# ===========================================================================

"""
    GridReference(name, system, idscheme, schemes)

What a canonical grid name resolves to: the grid system, the id scheme to
assume when a store names none, and the schemes that name accepts.

A grid name pins the id packing, not just the tessellation: `"igeo7"` and
`"isea7h"` are the same cells under different ids, and a reader that treated
one as the other would misplace every cell. Unknown names are therefore
rejected rather than guessed.
"""
struct GridReference
    name::String
    system::Any
    idscheme::Symbol
    schemes::Tuple{Vararg{Symbol}}
end

"""
    GRID_REFERENCE

Canonical grid name → [`GridReference`](@ref). Shared by every convention;
[`gridname`](@ref) is what maps a store's own spelling onto a key of this table,
and [`register_grid!`](@ref) is how a downstream grid system gets one.
"""
const GRID_REFERENCE = Dict{String,GridReference}(
    "igeo7" => GridReference("igeo7", IGeo7System(), :z7int, (:z7int,)),
    "healpix" => GridReference("healpix", HEALPixSystem(), :nested,
        (:nested, :ring, :nuniq, :zuniq)))

"""
    register_grid!(name::AbstractString, ref::GridReference) -> GRID_REFERENCE

Register the canonical store spelling `name` for a grid system.

Every convention resolves a store's grid name through [`GRID_REFERENCE`](@ref),
and an unregistered name is refused rather than guessed, so this is what makes a
downstream grid system readable and writeable. The name pins the id packing as
well as the tessellation: register `"isea7h"` separately from `"igeo7"` if the
ids differ, rather than aliasing one onto the other.
"""
register_grid!(name::AbstractString, ref::GridReference) =
    setindex!(GRID_REFERENCE, ref, String(name))

"""
    gridname(c::DGGSConvention, name, attrs) -> String

Canonical grid name for the raw `name` a store carries, given the attribute
dictionary it came from.

The default normalizes case and whitespace and otherwise trusts the store, so
that unregistered names fail loudly in [`gridreference`](@ref) rather than
being decoded as something else. A convention overrides this to fold its own
vocabulary in:

```julia
DGG.gridname(::MyConvention, name, attrs) =
    name == "ISEA7H_Z7" ? "igeo7" : invoke(DGG.gridname, Tuple{DGG.DGGSConvention,Any,Any}, MyConvention(), name, attrs)
```
"""
gridname(::DGGSConvention, name, attrs) = lowercase(strip(String(name)))

"""
    gridreference(canonical; store = "", conventions = String[]) -> GridReference

Look `canonical` up in [`GRID_REFERENCE`](@ref), or raise a
[`DGGSFormatError`](@ref) listing every registered name.
"""
function gridreference(canonical::AbstractString; store="", conventions=String[])
    ref = get(GRID_REFERENCE, canonical, nothing)
    ref === nothing || return ref
    throw(DGGSFormatError(check=:unknown_grid_name, store=store,
        conventions=conventions, declared=String(canonical),
        detail="registered canonical grid names: " *
               join(sort!(collect(keys(GRID_REFERENCE))), ", ") *
               ". Add one with `register_grid!(name, GridReference(...))`, or give " *
               "the convention a `gridname` method that maps this spelling onto a " *
               "registered name."))
end

# ===========================================================================
# Attribute plumbing
# ===========================================================================

_isdict(x) = x isa AbstractDict

function _require(attrs, key, store, convs, what)
    haskey(attrs, key) && attrs[key] !== nothing && return attrs[key]
    throw(DGGSFormatError(check=:missing_required_key, store=store,
        conventions=convs, declared=key,
        detail="$what requires the `$key` attribute."))
end

function _asint(x, key, store, convs)
    x isa Integer && return Int(x)
    if x isa Real && isinteger(x)
        return Int(x)
    end
    throw(DGGSFormatError(check=:invalid_attribute_type, store=store,
        conventions=convs, declared=key, observed=x,
        detail="`$key` must be an integer."))
end

_asfloat(x) = x isa Real ? Float64(x) : nothing

_asstring(::Nothing, key, store, convs) = nothing

function _asstring(x, key, store, convs)
    x isa AbstractString && return String(x)
    throw(DGGSFormatError(check=:invalid_attribute_type, store=store,
        conventions=convs, declared=key, observed=x,
        detail="`$key` must be a string."))
end

function _asbool(x, key, store, convs)
    x isa Bool && return x
    throw(DGGSFormatError(check=:invalid_attribute_type, store=store,
        conventions=convs, declared=key, observed=x,
        detail="`$key` must be a boolean."))
end

"""
    encoding_for(vocabulary; store = "", conventions = String[])

The registered `CellEncoding` instance for a convention's vocabulary string
(`"none"`, `"ranges"`, `"implicit"`, …).
"""
function encoding_for(vocab::AbstractString; store="", conventions=String[])
    enc = get(ENCODING_REGISTRY, String(vocab), nothing)
    enc === nothing || return enc
    throw(DGGSFormatError(check=:unknown_encoding, store=store,
        conventions=conventions, declared=String(vocab),
        detail="registered encodings: " *
               join(sort!(collect(keys(ENCODING_REGISTRY))), ", ") * "."))
end

# The four schemes CF and `zarr-conventions/dggs` name, plus the `nest` boolean
# spelling that the DKRZ corpus and the xdggs alias dialect use.
function _scheme(raw, ref::GridReference; store="", conventions=String[])
    s = raw isa Bool ? (raw ? :nested : :ring) : Symbol(lowercase(strip(String(raw))))
    s === :nest && (s = :nested)
    s in ref.schemes && return s
    throw(DGGSFormatError(check=:unknown_indexing_scheme, store=store,
        conventions=conventions, declared=String(s), observed=ref.name,
        detail="`$(ref.name)` accepts " * join(ref.schemes, ", ") * "."))
end

# nside is 2^level and the two are both bare integers, so reading one as the
# other is a silent 2^level error; only the key name disambiguates.
function _level_from_nside(nside, store, convs)
    n = _asint(nside, "nside", store, convs)
    (n > 0 && ispow2(n)) && return trailing_zeros(n)
    throw(DGGSFormatError(check=:invalid_nside, store=store, conventions=convs,
        declared=n, detail="nside must be a positive power of two."))
end

"""
    parse_ellipsoid(x; store = "", conventions = String[]) -> Union{Ellipsoid, Nothing}

Read an ellipsoid declaration. Accepts a bare name string, or an object in
either axis-key spelling: the schema-valid `semi_major_axis` and the
`semimajor_axis` form every store in the wild emits.

`nothing` when the store declares none, so that silence stays distinguishable
from a declaration and [`ellipsoid`](@ref) can apply
[`DEFAULT_ELLIPSOID`](@ref). An object carrying no ellipsoid key at all — the
DKRZ `crs` variable — is silence.
"""
function parse_ellipsoid(x; store="", conventions=String[])
    x === nothing && return nothing
    x isa AbstractString && return Ellipsoid(name=String(x))
    _isdict(x) || return nothing
    major = get(x, "semi_major_axis", get(x, "semimajor_axis", nothing))
    minor = get(x, "semi_minor_axis", get(x, "semiminor_axis", nothing))
    return _ellipsoid_or_nothing(Ellipsoid(
        name=_asstring(get(x, "name", get(x, "reference_ellipsoid_name", nothing)),
            "name", store, conventions),
        semi_major_axis=_asfloat(major),
        semi_minor_axis=_asfloat(minor),
        inverse_flattening=_asfloat(get(x, "inverse_flattening", nothing)),
        radius=_asfloat(get(x, "radius", get(x, "earth_radius", nothing)))))
end

"""
    ellipsoid_attrs(e::Ellipsoid) -> Dict{String,Any}

An ellipsoid as write attributes, in the schema-valid `semi_major_axis`
spelling. This deliberately diverges from the stores in the wild, which emit
`semimajor_axis` and therefore fail the schema they declare.
"""
function ellipsoid_attrs(e::Ellipsoid)
    out = Dict{String,Any}()
    e.name === nothing || (out["name"] = e.name)
    e.semi_major_axis === nothing || (out["semi_major_axis"] = e.semi_major_axis)
    e.semi_minor_axis === nothing || (out["semi_minor_axis"] = e.semi_minor_axis)
    e.inverse_flattening === nothing || (out["inverse_flattening"] = e.inverse_flattening)
    e.radius === nothing || (out["radius"] = e.radius)
    return out
end

# A parameter every field of which is `nothing` is a store saying nothing, not a
# store declaring emptiness: `nothing` is what lets the convention default
# through. Both are applied where the parameter is parsed, so every convention
# gets the distinction from one place.
_orientation_or_nothing(o::GridOrientation) =
    o == GridOrientation() ? nothing : o

_ellipsoid_or_nothing(e::Ellipsoid) = e == Ellipsoid() ? nothing : e

# ===========================================================================
# Convention A -- zarr-conventions/dggs
# ===========================================================================

"""
    ZARR_DGGS_UUID

The identifier of the `zarr-conventions/dggs` convention. Its meta-spec says
the declaration's `name` MUST NOT be used to identify a convention, so this is
what detection matches on.
"""
const ZARR_DGGS_UUID = "7b255807-140c-42ca-97f6-7a1cfecdbc38"

const ZARR_DGGS_DECLARATION = Dict{String,Any}(
    "description" => "Discrete Global Grid Systems convention for zarr",
    "name" => "dggs",
    "schema_url" => "https://raw.githubusercontent.com/zarr-conventions/dggs/refs/tags/v1/schema.json",
    "spec_url" => "https://github.com/zarr-conventions/dggs/blob/v1/README.md",
    "uuid" => ZARR_DGGS_UUID)

# Convention A's own compression vocabulary. `implicit` is this package's name
# for a missing `coordinate`, and is expressed by omitting the key.
const ZARR_DGGS_COMPRESSIONS = ("none", "compacted", "ranges")

"""
    ZarrDGGSConvention()

`zarr-conventions/dggs`: a `zarr_conventions` declaration on the group naming
[`ZARR_DGGS_UUID`](@ref), plus a `dggs` object carrying `name`,
`refinement_level`, `spatial_dimension` and optionally `coordinate`,
`compression`, `indexing_scheme` and `ellipsoid`.

This is the only convention that can express a non-dense encoding, so on a
ranges store it is the one that says so; the flat attrs of
[`XdggsConvention`](@ref) ride along on the same store and describe less.

Group-level declarations only: the spec's array-level override is not read.
"""
struct ZarrDGGSConvention <: DGGSConvention end

conventionname(::ZarrDGGSConvention) = "zarr-conventions/dggs"

function detect(c::ZarrDGGSConvention, s::StoreSnapshot)
    decls = get(s.attrs, "zarr_conventions", nothing)
    decls isa AbstractVector || return nothing
    for entry in decls
        _isdict(entry) || continue
        # Identity precedence is uuid > schema_url > spec_url; `name` is not an
        # identifier and is never consulted.
        matched = get(entry, "uuid", nothing) == ZARR_DGGS_UUID ||
                  _url_names_dggs(get(entry, "schema_url", nothing)) ||
                  _url_names_dggs(get(entry, "spec_url", nothing))
        matched || continue
        haskey(s.attrs, "dggs") || throw(DGGSFormatError(
            check=:missing_required_key, store=s.identifier,
            conventions=[conventionname(c)], declared="dggs",
            detail="the group declares zarr-conventions/dggs but carries no `dggs` object; array-level overrides are not read."))
        return Detection(c, :declared, "zarr_conventions uuid $ZARR_DGGS_UUID")
    end
    return nothing
end

_url_names_dggs(u) = u isa AbstractString && occursin("zarr-conventions/dggs", u)

function decode(c::ZarrDGGSConvention, s::StoreSnapshot, ::Detection)
    store, convs = s.identifier, [conventionname(c)]
    dggs = s.attrs["dggs"]
    _isdict(dggs) || throw(DGGSFormatError(check=:invalid_attribute_type,
        store=store, conventions=convs, declared="dggs", observed=dggs,
        detail="the `dggs` attribute must be an object."))

    canonical = gridname(c, _require(dggs, "name", store, convs, "the dggs object"), dggs)
    ref = gridreference(canonical; store, conventions=convs)

    # The key is required, but an explicit null is a declaration in its own
    # right: the axis is variable-sized. Absence is a store that forgot to say.
    haskey(dggs, "refinement_level") || _require(dggs, "refinement_level", store, convs, "the dggs object")
    rl = dggs["refinement_level"]
    level = rl === nothing ? nothing : _asint(rl, "refinement_level", store, convs)

    spatial = _asstring(_require(dggs, "spatial_dimension", store, convs, "the dggs object"),
        "spatial_dimension", store, convs)
    coordinate = _asstring(get(dggs, "coordinate", nothing), "coordinate", store, convs)

    # A missing `coordinate` means the whole domain is covered and position is
    # the id; `compression` may only be present when `coordinate` is.
    encoding = if coordinate === nothing
        encoding_for("implicit"; store, conventions=convs)
    else
        encoding_for(String(get(dggs, "compression", "none")); store, conventions=convs)
    end

    # A healpix name does not pin the id packing, and this convention need not
    # say which it is: left `nothing`, the merge or `applydefaults` supplies it,
    # and refusing here would reject a store another convention describes.
    scheme = if haskey(dggs, "indexing_scheme")
        _scheme(dggs["indexing_scheme"], ref; store, conventions=convs)
    elseif canonical == "healpix"
        nothing
    else
        ref.idscheme
    end

    orientation = GridOrientation(
        vert0_lon=_asfloat(get(dggs, "dggs_vert0_lon", nothing)),
        vert0_lat=_asfloat(get(dggs, "dggs_vert0_lat", nothing)),
        vert0_azimuth=_asfloat(get(dggs, "dggs_vert0_azimuth", nothing)),
        rotation_pattern=_asstring(get(dggs, "rotation_pattern", nothing),
            "rotation_pattern", store, convs))

    return StoreDescription(gridname=canonical, system=ref.system, idscheme=scheme,
        level=level, encoding=encoding, coordinate=coordinate,
        spatial_dimension=spatial,
        variables=spatial === nothing ? nothing :
                  arrays_on(s, spatial; exclude=(coordinate,)),
        ellipsoid=parse_ellipsoid(get(dggs, "ellipsoid", nothing); store, conventions=convs),
        orientation=_orientation_or_nothing(orientation),
        provenance=Dict{String,Any}(conventionname(c) => deepcopy(dggs)))
end

# ===========================================================================
# Convention B -- xdggs flat coordinate attrs
# ===========================================================================

# Every alias xdggs still accepts for the level, and the two spellings of the
# indexing scheme. Supplying two aliases for one parameter is an error there
# and here: `nside` and `level` are both bare integers and differ by `2^level`.
const XDGGS_LEVEL_ALIASES = ("level", "refinement_level", "nside", "order",
    "resolution", "depth")
const XDGGS_SCHEME_ALIASES = ("indexing_scheme", "nest")

"""
    XdggsConvention()

The xdggs dialect: flat grid attributes on the cell-id coordinate —
`grid_name`, `level` (or an alias), `indexing_scheme` (or `nest`), `ellipsoid`,
plus per-plugin extras such as `igeo7_dggs_vert0_lon`.

The only convention most producers emit, and the only one on stores whose group
attributes are empty. It describes a dense 1-D axis and nothing else: on the
`(n, 2)` coordinate of a ranges store it decodes the grid and leaves the
encoding to [`ZarrDGGSConvention`](@ref).

Writing puts ONLY grid keys on the coordinate. xdggs forwards every attribute
but `grid_name` as a dataclass keyword, so one stray `units` makes the store
unreadable there.
"""
struct XdggsConvention <: DGGSConvention end

conventionname(::XdggsConvention) = "xdggs"

function detect(c::XdggsConvention, s::StoreSnapshot)
    cands = [a for a in s.arrays if haskey(a.attrs, "grid_name")]
    isempty(cands) && return nothing
    length(cands) == 1 || throw(DGGSFormatError(check=:multiple_cell_coordinates,
        store=s.identifier, conventions=[conventionname(c)],
        observed=[a.name for a in cands],
        detail="one group may carry grid attributes on one coordinate only."))
    return Detection(c, :fingerprint, "grid_name on array $(cands[1].name)";
        payload=Dict{Symbol,Any}(:coordinate => cands[1].name))
end

decode(c::XdggsConvention, s::StoreSnapshot, det::Detection) =
    decode_flat_attrs(c, s, det)

"""
    decode_flat_attrs(c, snapshot, detection) -> StoreDescription

Decode flat grid attributes from the array named by `detection.payload[:coordinate]`.

Shared by every convention whose metadata is a flat dictionary on the cell-id
coordinate, and the reason a convention that only overrides [`gridname`](@ref)
needs no decoder of its own.
"""
function decode_flat_attrs(c::DGGSConvention, s::StoreSnapshot, det::Detection)
    store, convs = s.identifier, [conventionname(c)]
    arr = getarray(s, det.payload[:coordinate])
    attrs = arr.attrs

    canonical = gridname(c, _require(attrs, "grid_name", store, convs, "the coordinate"), attrs)
    ref = gridreference(canonical; store, conventions=convs)
    level = _alias_level(attrs, store, convs)
    scheme = _alias_scheme(attrs, ref, store, convs)

    # A 1-D coordinate is a dense axis. Anything else -- the ranges stores put
    # these same attrs on an `(n, 2)` array -- is a layout with no name here.
    onedim = ndims(arr) == 1
    spatial = onedim && !isempty(arr.dims) ? arr.dims[1] : nothing
    encoding = onedim ? encoding_for("none"; store, conventions=convs) : nothing

    orientation = GridOrientation(
        vert0_lon=_asfloat(get(attrs, "igeo7_dggs_vert0_lon", nothing)))
    conversion = haskey(attrs, "igeo7_wgs84_geodetic_conversion") ?
                 _asbool(attrs["igeo7_wgs84_geodetic_conversion"],
        "igeo7_wgs84_geodetic_conversion", store, convs) : nothing

    return StoreDescription(gridname=canonical, system=ref.system, idscheme=scheme,
        level=level, encoding=encoding, coordinate=arr.name,
        spatial_dimension=spatial,
        variables=spatial === nothing ? nothing :
                  arrays_on(s, spatial; exclude=(arr.name,)),
        ellipsoid=parse_ellipsoid(get(attrs, "ellipsoid", nothing); store, conventions=convs),
        orientation=_orientation_or_nothing(orientation),
        geodetic_conversion=conversion,
        provenance=Dict{String,Any}(conventionname(c) => deepcopy(attrs)))
end

function _alias_level(attrs, store, convs)
    present = [k for k in XDGGS_LEVEL_ALIASES
               if haskey(attrs, k) && attrs[k] !== nothing]
    isempty(present) && return nothing
    length(present) == 1 || throw(DGGSFormatError(check=:duplicate_level_alias,
        store=store, conventions=convs, observed=present,
        detail="the level may be given once, under one of " *
               join(XDGGS_LEVEL_ALIASES, ", ") * "."))
    key = only(present)
    key == "nside" && return _level_from_nside(attrs[key], store, convs)
    return _asint(attrs[key], key, store, convs)
end

function _alias_scheme(attrs, ref, store, convs)
    present = [k for k in XDGGS_SCHEME_ALIASES
               if haskey(attrs, k) && attrs[k] !== nothing]
    isempty(present) && return nothing
    length(present) == 1 ||
        throw(DGGSFormatError(check=:duplicate_indexing_scheme_alias, store=store,
            conventions=convs, observed=present,
            detail="the indexing scheme may be given once, as `indexing_scheme` or `nest`."))
    return _scheme(attrs[only(present)], ref; store, conventions=convs)
end

function encode!(c::XdggsConvention, s::StoreSnapshot, d::StoreDescription)
    d.coordinate === nothing && return s
    arr = getarray(s, d.coordinate)
    arr === nothing && throw(ArgumentError(
        "cannot stamp xdggs attrs: the snapshot has no array named $(repr(d.coordinate))"))
    attrs = arr.attrs
    # Only grid keys, and only the ones the reader's dataclass takes.
    attrs["grid_name"] = d.gridname
    d.level === nothing || (attrs["level"] = d.level)
    if d.gridname == "healpix"
        d.idscheme === nothing || (attrs["indexing_scheme"] = String(d.idscheme))
        d.ellipsoid === nothing || (attrs["ellipsoid"] = ellipsoid_attrs(d.ellipsoid))
    elseif d.gridname == "igeo7"
        d.geodetic_conversion === nothing ||
            (attrs["igeo7_wgs84_geodetic_conversion"] = d.geodetic_conversion)
        if d.orientation !== nothing && d.orientation.vert0_lon !== nothing
            attrs["igeo7_dggs_vert0_lon"] = d.orientation.vert0_lon
        end
    end
    return s
end

# ===========================================================================
# The legacy HEALPix trio
# ===========================================================================

"""
    LegacyHealpixConvention()

`grid_name: "healpix"` with the `nside`/`nest` pair rather than
`level`/`indexing_scheme` — the dialect the Ifremer and Grid4Earth demo stores
ship.

The trio is also within [`XdggsConvention`](@ref)'s alias set, so both fire on
such a store and merge to the same description. This convention is the named
extension point for the dialect: it is what provenance and error messages call
it, and it is where a healpix-specific reading of the trio dispatches.
"""
struct LegacyHealpixConvention <: DGGSConvention end

conventionname(::LegacyHealpixConvention) = "legacy-healpix"

_legacy_trio(attrs) = get(attrs, "grid_name", nothing) == "healpix" &&
                      haskey(attrs, "nside") && haskey(attrs, "nest")

function detect(c::LegacyHealpixConvention, s::StoreSnapshot)
    for a in s.arrays
        _legacy_trio(a.attrs) || continue
        return Detection(c, :fingerprint, "grid_name/nside/nest on array $(a.name)";
            payload=Dict{Symbol,Any}(:coordinate => a.name))
    end
    return nothing
end

decode(c::LegacyHealpixConvention, s::StoreSnapshot, det::Detection) =
    decode_flat_attrs(c, s, det)

# ===========================================================================
# The DKRZ / nextGEMS dialect
# ===========================================================================

# easygems' fallbacks, in order. The dimension is not named by any attribute.
const DKRZ_SPATIAL_DIMS = ("cell", "cells", "value", "values")

"""
    DKRZConvention()

The nextGEMS / DestinE / easy.gems dialect: a `crs` variable carrying
`grid_mapping_name: "healpix"`, `healpix_nside` and `healpix_order`.

Two things distinguish it from CF 1.13, which has the same shape.
`healpix_nside` is an nside where CF's `refinement_level` is an order, and
`healpix_order` takes the value `"nest"`. Most of these stores hold **no cell
array at all**: position is the nested index and the axis is implicit. Where a
cell array does exist the store is regional and holds global indices, so its
presence is the sparsity signal.

Data variables' `grid_mapping` back-references are known to dangle, so
detection scans every array for the grid mapping rather than following one.
"""
struct DKRZConvention <: DGGSConvention end

conventionname(::DKRZConvention) = "dkrz-healpix"

function detect(c::DKRZConvention, s::StoreSnapshot)
    for a in s.arrays
        at = a.attrs
        (get(at, "grid_mapping_name", nothing) == "healpix" &&
         haskey(at, "healpix_nside")) || continue
        return Detection(c, :fingerprint, "healpix_nside on array $(a.name)";
            payload=Dict{Symbol,Any}(:crs => a.name))
    end
    return nothing
end

function decode(c::DKRZConvention, s::StoreSnapshot, det::Detection)
    store, convs = s.identifier, [conventionname(c)]
    crs = det.payload[:crs]
    attrs = getarray(s, crs).attrs

    canonical = gridname(c, attrs["grid_mapping_name"], attrs)
    ref = gridreference(canonical; store, conventions=convs)
    level = _level_from_nside(attrs["healpix_nside"], store, convs)
    scheme = haskey(attrs, "healpix_order") ?
             _scheme(attrs["healpix_order"], ref; store, conventions=convs) : nothing

    spatial = _dkrz_spatial(s, crs)
    coordinate = spatial !== nothing && getarray(s, spatial) !== nothing ?
                 spatial : nothing
    encoding = encoding_for(coordinate === nothing ? "implicit" : "none";
        store, conventions=convs)

    return StoreDescription(gridname=canonical, system=ref.system, idscheme=scheme,
        level=level, encoding=encoding, coordinate=coordinate,
        spatial_dimension=spatial,
        variables=spatial === nothing ? nothing :
                  arrays_on(s, spatial; exclude=(crs, coordinate)),
        ellipsoid=parse_ellipsoid(attrs; store, conventions=convs),
        provenance=Dict{String,Any}(conventionname(c) => deepcopy(attrs)))
end

function _dkrz_spatial(s::StoreSnapshot, crs::AbstractString)
    seen = Set{String}()
    for a in s.arrays
        a.name == crs && continue
        union!(seen, a.dims)
    end
    for d in DKRZ_SPATIAL_DIMS
        d in seen && return d
    end
    return nothing
end

# ===========================================================================
# The registry, detection and merging
# ===========================================================================

"""
    CONVENTION_REGISTRY

The conventions [`describe_store`](@ref) tries, in order: UUID-declared before
fingerprinted, so a store that names its convention is decoded by that one
first. Downstream packages extend it with [`register_convention!`](@ref).
"""
const CONVENTION_REGISTRY = DGGSConvention[
    ZarrDGGSConvention(), XdggsConvention(), LegacyHealpixConvention(),
    DKRZConvention()]

"""
    DEFAULT_WRITE_CONVENTIONS

What [`dggwrite`](@ref) stamps by default: `zarr-conventions/dggs` for the
encoding vocabulary a flat coordinate cannot express, and xdggs so the store
opens in the ecosystem's own reader. The `zarr-conventions/dggs` half is
written schema-VALID, which the stores in the wild are not.
"""
const DEFAULT_WRITE_CONVENTIONS = (ZarrDGGSConvention(), XdggsConvention())

"""
    register_convention!(c::DGGSConvention; first = false) -> CONVENTION_REGISTRY

Add `c` to [`CONVENTION_REGISTRY`](@ref), at the end or ahead of everything.
"""
register_convention!(c::DGGSConvention; first::Bool=false) =
    first ? pushfirst!(CONVENTION_REGISTRY, c) : push!(CONVENTION_REGISTRY, c)

"""
    detections(snapshot; conventions = CONVENTION_REGISTRY) -> Vector{Detection}

Every convention that recognizes `snapshot`, declared ones first and otherwise
in registry order.
"""
function detections(s::StoreSnapshot; conventions=CONVENTION_REGISTRY)
    found = Detection[]
    for c in conventions
        det = detect(c, s)
        det === nothing || push!(found, det)
    end
    return sort!(found; by=det -> _rankorder(det.rank), alg=MergeSort)
end

"""
    describe_store(snapshot; conventions = CONVENTION_REGISTRY) -> StoreDescription

Detect, decode and reconcile: the whole read-side metadata path.

Where several conventions fire their descriptions must agree field for field —
one may be silent where another speaks, but two answers to one question is the
"attrs lie" failure and raises [`DGGSFormatError`](@ref) rather than a guess.
Convention defaults are applied last, so a store's silence is never mistaken
for a declaration during reconciliation.

What comes back is usable: a description still missing its encoding or its
spatial dimension after every convention and every default has spoken names a
store nothing can open, and raises rather than being returned.
"""
function describe_store(s::StoreSnapshot; conventions=CONVENTION_REGISTRY)
    found = detections(s; conventions)
    isempty(found) && throw(DGGSFormatError(check=:no_convention_detected,
        store=s.identifier,
        detail="no registered convention recognized this group; pass `description` to dggread to assert the grid yourself."))
    names = String[conventionname(det.convention) for det in found]
    descs = [decode(det.convention, s, det) for det in found]
    desc = applydefaults(merge_descriptions(descs; store=s.identifier,
        conventions=names))
    return requirecomplete(desc; store=s.identifier, conventions=names)
end

"""
    requirecomplete(desc; store = "", conventions = String[]) -> StoreDescription

`desc` if it names both a layout and a data dimension, and a
[`DGGSFormatError`](@ref) naming what is absent otherwise.

The two fields a reader cannot do without: with no `encoding` there is no way to
turn the coordinate into an axis, and with no `spatial_dimension` there is no
way to say which arrays the axis belongs to. This is the backstop that lets each
convention decode only what its vocabulary covers.
"""
function requirecomplete(d::StoreDescription; store="", conventions=String[])
    absent = String[]
    d.encoding === nothing && push!(absent, "encoding")
    d.spatial_dimension === nothing && push!(absent, "spatial_dimension")
    isempty(absent) && return d
    throw(DGGSFormatError(check=:incomplete_description, store=store,
        conventions=conventions, declared=d.gridname, observed=absent,
        detail="none of the conventions that fired could supply " *
               join(absent, " or ") *
               (d.coordinate === nothing ? "" :
                " for the coordinate $(repr(d.coordinate))") *
               ". A non-dense layout is described by zarr-conventions/dggs " *
               "only; pass `description` to dggread to assert it yourself."))
end

"""
    merge_descriptions(descs; store = "", conventions = String[]) -> StoreDescription

Reconcile the descriptions several conventions decoded from one store.
`nothing` is silence and yields to a value; two values must be equal;
[`Ellipsoid`](@ref) and [`GridOrientation`](@ref) reconcile field-wise, since
conventions carry different subsets of them.
"""
function merge_descriptions(descs; store="", conventions=String[])
    length(descs) == 1 && return only(descs)
    vals = Dict{Symbol,Any}()
    for f in DESCRIPTION_FIELDS
        acc = nothing
        for d in descs
            acc = _mergefield(acc, getfield(d, f), f, store, conventions)
        end
        vals[f] = acc
    end
    prov = Dict{String,Any}()
    for d in descs
        merge!(prov, d.provenance)
    end
    vals[:provenance] = prov
    return StoreDescription(; vals...)
end

_mergefield(::Nothing, ::Nothing, field, store, convs) = nothing
_mergefield(a, ::Nothing, field, store, convs) = a
_mergefield(::Nothing, b, field, store, convs) = b

function _mergefield(a, b, field, store, convs)
    a == b && return a
    for T in (Ellipsoid, GridOrientation)
        a isa T && b isa T || continue
        return T(; (f => _mergefield(getfield(a, f), getfield(b, f),
            Symbol(field, ".", f), store, convs) for f in fieldnames(T))...)
    end
    throw(DGGSFormatError(check=Symbol(field, "_disagreement"), store=store,
        conventions=convs, declared=a, observed=b,
        detail="the store's conventions disagree about `$field`; refusing to guess which is right."))
end

"""
    applydefaults(desc) -> StoreDescription

Fill in what a convention prescribes for what a store left unsaid. Currently
the id scheme, which defaults to its grid's canonical one — nested for HEALPix,
as both xdggs and the DKRZ corpus assume.

Separate from decoding so that two conventions' silence never reconciles into a
contradiction between two defaults.
"""
function applydefaults(d::StoreDescription)
    d.idscheme === nothing || return d
    ref = get(GRID_REFERENCE, d.gridname, nothing)
    ref === nothing && return d
    return _with(d; idscheme=ref.idscheme)
end

function _with(d::StoreDescription; kwargs...)
    vals = Dict{Symbol,Any}(f => getfield(d, f) for f in fieldnames(StoreDescription))
    for (k, v) in kwargs
        vals[k] = v
    end
    return StoreDescription(; vals...)
end

# ===========================================================================
# Stamping convention A
# ===========================================================================

function encode!(c::ZarrDGGSConvention, s::StoreSnapshot, d::StoreDescription)
    dggs = Dict{String,Any}(
        "name" => d.gridname,
        "refinement_level" => d.level,
        "spatial_dimension" => d.spatial_dimension)

    if d.coordinate !== nothing
        d.encoding === nothing && throw(ArgumentError(
            "cannot stamp $(conventionname(c)): the description names a coordinate but no encoding"))
        vocab = encodingname(d.encoding)
        vocab in ZARR_DGGS_COMPRESSIONS || throw(ArgumentError(
            "$(conventionname(c)) has no `compression` value for the $vocab encoding"))
        dggs["coordinate"] = d.coordinate
        dggs["compression"] = vocab
    end

    d.gridname == "healpix" && d.idscheme !== nothing &&
        (dggs["indexing_scheme"] = String(d.idscheme))
    d.ellipsoid === nothing || (dggs["ellipsoid"] = ellipsoid_attrs(d.ellipsoid))

    if d.orientation !== nothing
        o = d.orientation
        o.vert0_lon === nothing || (dggs["dggs_vert0_lon"] = o.vert0_lon)
        o.vert0_lat === nothing || (dggs["dggs_vert0_lat"] = o.vert0_lat)
        o.vert0_azimuth === nothing || (dggs["dggs_vert0_azimuth"] = o.vert0_azimuth)
        o.rotation_pattern === nothing ||
            (dggs["rotation_pattern"] = o.rotation_pattern)
    end

    s.attrs["dggs"] = dggs
    decls = get!(s.attrs, "zarr_conventions", Any[])
    any(e -> _isdict(e) && get(e, "uuid", nothing) == ZARR_DGGS_UUID, decls) ||
        push!(decls, deepcopy(ZARR_DGGS_DECLARATION))
    return s
end
