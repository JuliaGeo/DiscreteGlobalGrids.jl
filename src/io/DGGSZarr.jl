# DGGSZarr.jl — reading Zarr archives written against the DGGS convention.
#
# The convention (zarr-conventions/dggs v1) puts a `dggs` object in a group's
# attributes naming the grid, its refinement level, the spatial dimension, and
# — the part that matters here — how the cell ids along that dimension are
# stored:
#
#   compression = "none"     a 1-D array of N ids, one per cell
#   compression = "ranges"   an (R, 2) table of inclusive [start, end] pairs
#   compression = "compacted"  (not implemented; see `_decode_coordinate`)
#
# The convention itself is grid-agnostic, so the parsing below is too; only the
# step from a stored coordinate to an id vector is per-system, and it is one
# function (`_decode_coordinate`) with one method per grid. IGEO7/Z7 is wired;
# adding HEALPix or H3 is a method, not a restructuring.
#
# Everything is lazy on purpose. `open_dggs_dataset` reads metadata plus, for a
# ranges archive, the R × 2 table — never a data variable and never a length-N
# coordinate. Data variables stay as `Zarr.ZArray`s inside `YAXArray`s (Zarr
# arrays are `DiskArrays`, so slicing reads chunks); dense coordinates stay
# behind a `Z7CachedIds` that reads on first use. Materialization is a
# consumer's decision — plotting, modelling, a numeric reduction — not a side
# effect of opening.
#
# Two things about Zarr that a reader must get right, both pinned by tests:
#
#  1. **Zarr.jl reverses dimension order.** Zarr stores C-ordered, Julia is
#     column-major, so an on-disk `(R, 2)` array arrives as `(2, R)` and the
#     `_ARRAY_DIMENSIONS` attribute has to be reversed to name Julia's axes.
#     For the 1-D data variables this is invisible; for the range table it is
#     the difference between R ranges and 2.
#
#  2. **`YAXArrays.open_dataset` cannot open a ranges archive at all.** It
#     requires every dimension to have a backing coordinate array, and the
#     `ranges` / `bounds` dimensions of the range table have none — it fails
#     with `KeyError: key "bounds" not found`. Hence the `Dataset` here is
#     assembled from a raw `zopen` rather than delegated.

"""
    DGGSZarr

Read Zarr archives written against the [DGGS Zarr
convention](https://github.com/zarr-conventions/dggs/blob/v1/README.md) as lazy
[`YAXArrays.Dataset`](@ref)s indexed by cell id.

The entry point is [`open_dggs_dataset`](@ref). The dimension it builds is a
real DGGS lookup, so the ordinary `DimensionalData` selectors work on the
result — `At(cell_id)`, `Contains(point)`, `Touching(geometry)` — and, for a
`compression: "ranges"` archive, they are answered against the range table in
`O(log R)` without the length-N coordinate ever existing.

```julia
using DiscreteGlobalGrids.DGGSZarr

ds = open_dggs_dataset("pori_z7_r12_ranges.zarr")   # reads ~17 kB, not 158k cells
info = dggs_info(ds)                                # igeo7, level 12, "ranges"
elevation = ds.elevation                            # lazy YAXArray
elevation[cell_ids = At(0x0042404000ffffff)]        # one chunk read
```

Geometry ([`cell_centers`](@ref), [`cell_boundaries`](@ref),
[`sel_latlon`](@ref)) is served through this module rather than through
`IGeo7Lookups` because it needs two things only the archive knows: the
icosahedron placement (`dggs_vert0_lon`, which is `11.2` in the Python-written
archives against the package default of `11.25`) and whether latitudes are to
be reported as WGS84 geodetic rather than on the authalic sphere the grid
actually tiles. Getting either wrong displaces every coordinate — the datum
alone by ~0.115° of latitude, some 13 km — without any operation failing.
"""
module DGGSZarr

using Zarr
using YAXArrays
using DimensionalData
using DimensionalData: Lookups
using GeoInterface

import ..DiscreteGlobalGrids as DGG
import ..Helpers
import ..IGeo7
import ..ISEA

const DD = DimensionalData
const GI = GeoInterface

export DGGSZarrInfo,
    cell_boundaries,
    cell_centers,
    dggs_cell_at,
    dggs_cell_ids,
    dggs_info,
    dggs_lookup,
    open_dggs_dataset,
    sel_bbox,
    sel_latlon

# The attribute key holding the convention object, and the coordinate-level
# keys the IGEO7 writer sets alongside it.
const DGGS_ATTR = "dggs"
const ARRAY_DIMENSIONS_ATTR = "_ARRAY_DIMENSIONS"

"""
    DGGSZarrInfo

The parsed `dggs` attribute block of an archive: what grid it is, at what
level, and how its cell ids are stored.

| field | convention key | notes |
|:--|:--|:--|
| `name` | `name` | lower-cased grid name, e.g. `"igeo7"` |
| `level` | `refinement_level` | `nothing` only for a variable-sized coordinate |
| `spatial_dimension` | `spatial_dimension` | the dimension the data lives on |
| `coordinate` | `coordinate` | array holding the ids; `nothing` means "whole globe" |
| `compression` | `compression` | `"none"`, `"ranges"` or `"compacted"` |
| `ellipsoid` | `ellipsoid` | `nothing` means the convention's default sphere |
| `extra` | — | every other key, e.g. IGEO7's `dggs_vert0_lon` |

`extra` is not a leftovers bin: the convention explicitly permits grid-specific
parameters, and for IGEO7 the placement of icosahedron vertex 0 lives there and
is required to compute any coordinate at all.
"""
struct DGGSZarrInfo
    name::String
    level::Union{Nothing,Int}
    spatial_dimension::String
    coordinate::Union{Nothing,String}
    compression::String
    ellipsoid::Union{Nothing,Dict{String,Any}}
    extra::Dict{String,Any}
end

function Base.show(io::IO, info::DGGSZarrInfo)
    print(io, "DGGSZarrInfo(", info.name, ", level=", info.level,
        ", dim=\"", info.spatial_dimension, "\"",
        ", compression=\"", info.compression, "\"")
    info.coordinate === nothing || print(io, ", coordinate=\"", info.coordinate, "\"")
    print(io, ")")
end
Base.show(io::IO, ::MIME"text/plain", info::DGGSZarrInfo) = show(io, info)

"""
    DGGSZarrInfo(attrs::AbstractDict) -> DGGSZarrInfo

Parse the `dggs` object out of a group's attributes.

Throws `ArgumentError` when the attributes carry no `dggs` block (the archive
does not follow the convention) or when one of the three required keys — `name`,
`refinement_level`, `spatial_dimension` — is missing.
"""
function DGGSZarrInfo(attrs::AbstractDict)
    if !haskey(attrs, DGGS_ATTR)
        found = join(sort!(String[String(k) for k in keys(attrs)]), ", ")
        throw(ArgumentError(
            "not a DGGS-convention archive: no \"$DGGS_ATTR\" key in the group " *
            "attributes (found: $found)"))
    end
    block = attrs[DGGS_ATTR]
    block isa AbstractDict || throw(ArgumentError(
        "the \"$DGGS_ATTR\" attribute must be an object, got $(typeof(block))"))

    required(key) = haskey(block, key) ? block[key] : throw(ArgumentError(
        "the \"$DGGS_ATTR\" object is missing the required key \"$key\""))

    name = lowercase(String(required("name")))
    raw_level = required("refinement_level")
    level = raw_level === nothing ? nothing : Int(raw_level)
    spatial_dimension = String(required("spatial_dimension"))
    coordinate = get(block, "coordinate", nothing)
    coordinate === nothing || (coordinate = String(coordinate))
    # `compression` is required only when a `coordinate` was given; with no
    # coordinate the archive covers the whole globe and there is nothing to
    # compress.
    compression = String(get(block, "compression", coordinate === nothing ? "none" : "none"))
    ellipsoid = get(block, "ellipsoid", nothing)
    ellipsoid === nothing || (ellipsoid = Dict{String,Any}(String(k) => v for (k, v) in ellipsoid))

    known = ("name", "refinement_level", "spatial_dimension", "coordinate",
        "compression", "ellipsoid")
    extra = Dict{String,Any}(String(k) => v for (k, v) in block if !(String(k) in known))

    return DGGSZarrInfo(name, level, spatial_dimension, coordinate, compression,
        ellipsoid, extra)
end

# --- Zarr plumbing ---------------------------------------------------------

"""
    julia_dimensions(array) -> Vector{String}

Names of `array`'s dimensions **in Julia's axis order**, i.e. the archive's
`_ARRAY_DIMENSIONS` reversed.

Zarr writes C-ordered (slowest-varying dimension first) and Zarr.jl presents
the array column-major, so an on-disk `(R, 2)` `["ranges", "bounds"]` array is
a Julia `(2, R)` array whose axes are `["bounds", "ranges"]`. One-dimensional
variables are unaffected, which is exactly why the mistake survives a casual
test.
"""
function julia_dimensions(array::Zarr.ZArray)
    dims = get(array.attrs, ARRAY_DIMENSIONS_ATTR, nothing)
    dims === nothing && throw(ArgumentError(
        "Zarr array is missing its \"$ARRAY_DIMENSIONS_ATTR\" attribute; " *
        "it was not written by an xarray-compatible writer"))
    return reverse(String.(dims))
end

"Read a Zarr array whole, as a plain Julia array."
_read_all(array::Zarr.ZArray) = array[ntuple(_ -> Colon(), ndims(array))...]

# --- coordinate decoding ---------------------------------------------------
#
# The one per-grid step. Everything above and below is convention plumbing;
# this is where a stored coordinate becomes an `AbstractVector` of cell ids and
# a `DimensionalData` lookup.

"""
    _decode_coordinate(::Val{name}, group, info, n) -> Lookups.Lookup

Build the DGGS lookup for a `compression` form of grid `name`, given the opened
Zarr `group`, the parsed `info`, and the length `n` of the spatial dimension.

One method per wired grid. The fallback below is what an unwired grid hits, and
it names the grid rather than failing obscurely inside the reader.
"""
function _decode_coordinate(::Val{name}, group, info::DGGSZarrInfo, n::Int) where {name}
    throw(ArgumentError(
        "no DGGS-Zarr coordinate decoder is wired for grid \"$name\" (wired: igeo7). " *
        "Adding one is a `_decode_coordinate(::Val{:$name}, ...)` method returning a " *
        "lookup over the archive's cell ids."))
end

_decode_coordinate(info::DGGSZarrInfo, group, n::Int) =
    _decode_coordinate(Val(Symbol(info.name)), group, info, n)

"""
    igeo7_vert0_lon(info) -> Float64

Longitude at which the archive places icosahedron vertex 0 — the convention's
grid-specific `dggs_vert0_lon`, defaulting to the standard ISEA
`$(ISEA.ISEA_LON0)` when absent.

The Python IGEO7 tooling writes `11.2`, so this is *not* a formality: read it
and pass it to [`IGeo7.vert0_lon_orientation`](@ref), or every centroid comes
out rotated by the difference.
"""
igeo7_vert0_lon(info::DGGSZarrInfo) =
    Float64(get(info.extra, "dggs_vert0_lon", ISEA.ISEA_LON0))

function _igeo7_check_placement(info::DGGSZarrInfo)
    azimuth = get(info.extra, "dggs_vert0_azimuth", 0.0)
    # `dggs_vert0_lat` and `dggs_vert0_azimuth` are not free parameters of this
    # implementation — the azimuth is baked into `ISEA.REFERENCE_EDGE` and the
    # latitude into the icosahedron — so an archive that moved either is
    # describing a grid we would silently mis-decode. Refuse instead.
    iszero(azimuth) || throw(ArgumentError(
        "archive declares dggs_vert0_azimuth=$azimuth; only 0 is supported " *
        "(the reference edge azimuth is fixed in ISEA.REFERENCE_EDGE)"))
    return nothing
end

function _decode_coordinate(::Val{:igeo7}, group, info::DGGSZarrInfo, n::Int)
    level = info.level
    level === nothing && throw(ArgumentError(
        "an igeo7 archive must declare a refinement_level (a variable-sized " *
        "coordinate has no single-resolution lookup)"))
    _igeo7_check_placement(info)

    coordinate = info.coordinate
    coordinate === nothing && throw(ArgumentError(
        "an igeo7 archive with no `coordinate` would mean the whole globe at " *
        "level $level ($(IGeo7.num_cells(level)) cells); that form is not read here"))
    haskey(group.arrays, coordinate) || throw(ArgumentError(
        "the dggs block names coordinate \"$coordinate\", which the archive does not contain"))
    array = group[coordinate]

    ids = if info.compression == "ranges"
        _igeo7_range_ids(array, level)
    elseif info.compression == "none"
        # Deliberately not read here — see `Z7CachedIds`. Opening an archive
        # must stay O(1) in N, and a dense coordinate at Estonia res-12 is
        # 95 MB of ids nobody has asked for yet.
        size(array, 1) == n || throw(ArgumentError(
            "the dense coordinate \"$coordinate\" has $(size(array, 1)) ids but the " *
            "spatial dimension \"$(info.spatial_dimension)\" is $n long"))
        IGeo7.Z7CachedIds(array, level, n)
    elseif info.compression == "compacted"
        throw(ArgumentError(
            "compression=\"compacted\" (a mixed-resolution compacted id array) is " *
            "not implemented; this reader handles \"none\" and \"ranges\""))
    else
        throw(ArgumentError("unknown dggs compression \"$(info.compression)\""))
    end

    length(ids) == n || throw(ArgumentError(
        "the decoded coordinate holds $(length(ids)) cells but the spatial dimension " *
        "\"$(info.spatial_dimension)\" is $n long"))

    metadata = Dict{String,Any}(
        "grid_name" => "igeo7",
        "resolution" => level,
        "indexing_scheme" => "z7-u64",
        "external_indexing_schemes" => ["z7-string"],
        "compression" => info.compression,
        "dggs_vert0_lon" => igeo7_vert0_lon(info),
        # Whether the archive's coordinates are quoted as WGS84 geodetic
        # latitudes rather than on the authalic sphere the grid tiles. The
        # IGEO7 writer sets this on the coordinate array; default to `true`
        # because that is what every archive written so far declares.
        "wgs84_geodetic_conversion" =>
            Bool(get(array.attrs, "igeo7_wgs84_geodetic_conversion", true)),
    )
    return IGeo7.IGeo7Lookups.IGeo7Lookup(ids; resolution=level, metadata)
end

"Build a `Z7RangeIds` from the archive's `(R, 2)` range table, transposition and all."
function _igeo7_range_ids(array::Zarr.ZArray, level::Int)
    ndims(array) == 2 || throw(ArgumentError(
        "a compression=\"ranges\" coordinate must be 2-D (n_ranges, 2), " *
        "got $(ndims(array)) dimensions"))
    table = _read_all(array)
    # Julia's axis order is the reverse of the stored one, so the (R, 2) table
    # arrives as (2, R): axis 1 is `bounds`. Rather than assume, locate the
    # bounds axis by name and fall back to the length-2 axis.
    dims = julia_dimensions(array)
    bounds_axis = findfirst(==("bounds"), dims)
    bounds_axis === nothing && (bounds_axis = findfirst(==(2), collect(size(table))))
    bounds_axis === nothing && throw(ArgumentError(
        "cannot tell which axis of the $(size(table)) range table holds " *
        "(start, end); its dimensions are $dims"))
    size(table, bounds_axis) == 2 || throw(ArgumentError(
        "the \"bounds\" axis of the range table must be length 2, got $(size(table, bounds_axis))"))
    starts, ends = bounds_axis == 1 ?
                   (view(table, 1, :), view(table, 2, :)) :
                   (view(table, :, 1), view(table, :, 2))
    return IGeo7.Z7RangeIds(starts, ends, level)
end

# --- opening ---------------------------------------------------------------

"""
    open_dggs_dataset(path; consolidated=false) -> YAXArrays.Dataset

Open a DGGS-convention Zarr archive lazily, with its spatial dimension indexed
by a real DGGS lookup.

Nothing about the data is read. For a `compression: "ranges"` archive the R × 2
range table is read (tens of kB) and expanded arithmetically on demand; for
`compression: "none"` even the dense id array waits until something selects by
cell id. Data variables are the archive's `Zarr.ZArray`s wrapped in
`YAXArray`s, so a slice reads only the chunks it touches.

Arrays that *are* the index — the declared `coordinate`, and any array named
after the spatial dimension — are consumed into the lookup rather than exposed
as variables. Arrays on other dimensions are skipped with the reason recorded
in the dataset's properties under `"dggs_skipped_variables"`.

```julia
ds = open_dggs_dataset(joinpath(working, "pori_z7_r10_ranges.zarr"))
dggs_info(ds)                                  # DGGSZarrInfo(igeo7, level=10, ...)
ds.elevation[cell_ids = At(first(dggs_cell_ids(ds)))]
```

Throws `ArgumentError` for an archive without a `dggs` attribute block, for a
grid with no wired decoder, and for a `compression` this reader does not
implement.
"""
function open_dggs_dataset(path::AbstractString; consolidated::Bool=false)
    group = zopen(path, "r"; consolidated)
    group isa Zarr.ZGroup || throw(ArgumentError(
        "$path is a Zarr array, not a group; a DGGS archive is a group"))
    attrs = group.attrs
    info = DGGSZarrInfo(attrs)

    spatial = info.spatial_dimension
    # The spatial dimension's length is not stated anywhere in the convention;
    # it is the length of the data variables that live on it. (The declared
    # coordinate cannot supply it: under "ranges" that array is R long, not N.)
    n = _spatial_length(group, spatial, info)

    lookup = _decode_coordinate(info, group, n)
    dim = DD.Dim{Symbol(spatial)}(lookup)

    cubes = Dict{Symbol,YAXArray}()
    skipped = Dict{String,String}()
    reserved = Set{String}(filter(!isnothing, [info.coordinate, spatial]))
    for name in sort!(collect(keys(group.arrays)))
        if name in reserved
            continue                       # it is the index, not data
        end
        array = group[name]
        dims = julia_dimensions(array)
        if dims != [spatial]
            skipped[name] = "on dimensions $(dims), not [\"$spatial\"]"
            continue
        end
        properties = Dict{String,Any}(String(k) => v for (k, v) in array.attrs
                                      if String(k) != ARRAY_DIMENSIONS_ATTR)
        cubes[Symbol(name)] = YAXArray((dim,), array, properties)
    end

    properties = Dict{String,Any}(String(k) => v for (k, v) in attrs)
    properties["dggs_info"] = info
    isempty(skipped) || (properties["dggs_skipped_variables"] = skipped)
    return YAXArrays.Dataset(; properties, (k => v for (k, v) in cubes)...)
end

function _spatial_length(group, spatial::AbstractString, info::DGGSZarrInfo)
    for name in keys(group.arrays)
        name == info.coordinate && continue
        array = group[name]
        dims = try
            julia_dimensions(array)
        catch
            continue
        end
        dims == [spatial] && return size(array, 1)
    end
    # No data variable on the spatial dimension. A dense coordinate still
    # states the length; a range table does not, and the archive is unreadable.
    if info.compression == "none" && info.coordinate !== nothing &&
       haskey(group.arrays, info.coordinate)
        return size(group[info.coordinate], 1)
    end
    throw(ArgumentError(
        "cannot determine the length of the spatial dimension \"$spatial\": the archive " *
        "has no data variable on it, and a compression=\"$(info.compression)\" coordinate " *
        "does not state it"))
end

# --- accessors -------------------------------------------------------------

"""
    dggs_info(ds) -> DGGSZarrInfo

The parsed convention block of a dataset opened by [`open_dggs_dataset`](@ref).
"""
function dggs_info(ds::YAXArrays.Dataset)
    haskey(ds.properties, "dggs_info") || throw(ArgumentError(
        "this dataset carries no DGGS convention info; was it opened with open_dggs_dataset?"))
    return ds.properties["dggs_info"]::DGGSZarrInfo
end

"""
    dggs_lookup(x) -> Lookups.Lookup

The DGGS lookup indexing `x`, which may be a `YAXArray`, a `DimArray`, a
`Dimension`, or a `Dataset` opened by [`open_dggs_dataset`](@ref).

Throws `ArgumentError` when `x` has no DGGS-indexed dimension — including the
case of a dataset whose variables were all skipped.
"""
function dggs_lookup end

dggs_lookup(l::DGG.AbstractDGGSLookup) = l
dggs_lookup(d::DD.Dimension) = dggs_lookup(DD.lookup(d))

function dggs_lookup(x::Union{YAXArray,DD.AbstractDimArray})
    for d in DD.dims(x)
        DD.lookup(d) isa DGG.AbstractDGGSLookup && return DD.lookup(d)
    end
    throw(ArgumentError("no DGGS-indexed dimension found on $(typeof(x))"))
end

function dggs_lookup(ds::YAXArrays.Dataset)
    for d in values(ds.axes)
        DD.lookup(d) isa DGG.AbstractDGGSLookup && return DD.lookup(d)
    end
    throw(ArgumentError(
        "no DGGS-indexed dimension found in this dataset; it may hold no variables " *
        "on the spatial dimension (see properties[\"dggs_skipped_variables\"])"))
end

"""
    dggs_cell_ids(x) -> AbstractVector

The cell ids indexing `x`, in dimension order.

For a `compression: "ranges"` archive this is a [`Z7RangeIds`](@ref) — an
`AbstractVector` that computes each id in `O(log R)` rather than storing it, so
holding on to it costs `O(R)`. `collect` it to materialize.
"""
dggs_cell_ids(x) = DD.parent(dggs_lookup(x))

# --- geometry --------------------------------------------------------------
#
# These read the archive's placement and datum off the lookup metadata, which is
# what separates them from the same-named `IGeo7Lookups` functions. See the
# module docstring for why that separation exists rather than a keyword.

"The `ISEA.Orientation` an IGEO7 lookup's archive was written in."
_orientation(l) = IGeo7.vert0_lon_orientation(
    get(Lookups.metadata(l), "dggs_vert0_lon", ISEA.ISEA_LON0))

_wgs84(l) = Bool(get(Lookups.metadata(l), "wgs84_geodetic_conversion", true))

# Authalic (the sphere the grid tiles) -> the datum the archive quotes.
@inline _out_lat(lat, wgs84) =
    wgs84 ? Helpers.authalic_to_geodeticd(Helpers.WGS84_AUTHALIC, lat) : lat
# ... and back, for a coordinate arriving from the caller.
@inline _in_lat(lat, wgs84) =
    wgs84 ? Helpers.geodetic_to_authalicd(Helpers.WGS84_AUTHALIC, lat) : lat

"""
    cell_centers(x) -> Vector{Tuple{Float64,Float64}}

`(lon, lat)` center of every cell of `x`, in degrees and in dimension order,
honouring the archive's icosahedron placement and latitude datum.

Latitudes are WGS84 geodetic when the archive declares
`igeo7_wgs84_geodetic_conversion` (all archives written so far do); otherwise
they are left on the authalic sphere, which is what
`IGeo7Lookups.cell_centers` returns unconditionally. The difference is about
0.115° — roughly 13 km — so the two are not interchangeable.

This materializes: it is `N` centroid computations, and the point at which a
lazily-opened archive becomes real data.
"""
function cell_centers(x)
    lookup = dggs_lookup(x)
    orientation = _orientation(lookup)
    wgs84 = _wgs84(lookup)
    ids = DD.parent(lookup)
    out = Vector{Tuple{Float64,Float64}}(undef, length(ids))
    for (i, id) in enumerate(ids)
        lon, lat = IGeo7.cell_center(UInt64(id); orientation)
        @inbounds out[i] = (lon, _out_lat(lat, wgs84))
    end
    return out
end

"""
    cell_boundaries(x; closed_ring=true) -> Vector{<:GI.Polygon}

Cell boundary polygons for every cell of `x`, in `(lon, lat)` degrees under the
archive's placement and datum (see [`cell_centers`](@ref)).

Rings are hexagonal except at the twelve pentagons, and counter-clockwise seen
from outside the sphere. Like `cell_centers`, this materializes.
"""
function cell_boundaries(x; closed_ring::Bool=true)
    lookup = dggs_lookup(x)
    orientation = _orientation(lookup)
    wgs84 = _wgs84(lookup)
    ids = DD.parent(lookup)
    return map(ids) do id
        ring = IGeo7.cell_boundary(UInt64(id); closed_ring, orientation)
        GI.Polygon([GI.LinearRing([(lon, _out_lat(lat, wgs84)) for (lon, lat) in ring])])
    end
end

"""
    dggs_cell_at(x, lon, lat) -> UInt64

Cell id of `x`'s grid containing the point `(lon, lat)` in degrees, with `lat`
interpreted in the archive's datum.

The cell is the one the *grid* holds at that point, which need not be a cell
the archive stores; [`sel_latlon`](@ref) is the version that selects.
"""
function dggs_cell_at(x, lon::Real, lat::Real)
    lookup = dggs_lookup(x)
    level = lookup.resolution
    return IGeo7.lonlat_to_cell(Float64(lon), _in_lat(Float64(lat), _wgs84(lookup)), level;
        orientation=_orientation(lookup))
end

"""
    sel_latlon(x, lon, lat) -> selection

The element of `x` whose cell contains `(lon, lat)` in degrees — the equivalent
of the xdggs `sel_latlon`, and shorthand for `x[cell_ids = At(dggs_cell_at(...))]`.

Throws `ArgumentError` when the point falls outside the cells the archive holds,
which for a regional archive is most of the world.
"""
function sel_latlon(x, lon::Real, lat::Real)
    id = dggs_cell_at(x, lon, lat)
    lookup = dggs_lookup(x)
    DGG.cell_position(DD.parent(lookup), id) === nothing && throw(ArgumentError(
        "($lon, $lat) falls in cell 0x$(string(id; base=16, pad=16)), which this archive " *
        "does not hold"))
    return _select(x, DD.name(_dggs_dim(x)), Lookups.At(id))
end

"""
    sel_bbox(x, lon_min, lon_max, lat_min, lat_max) -> selection

Every cell of `x` whose *center* falls inside the given degree bounding box —
the equivalent of the xdggs bbox sample.

Selection is by center, matching `Contains` rather than `Touching`: a cell
straddling the boundary is included exactly when its centroid is inside. The
centers are computed (see [`cell_centers`](@ref)), so this materializes the
coordinate even though the data slice that follows stays lazy.
"""
function sel_bbox(x, lon_min::Real, lon_max::Real, lat_min::Real, lat_max::Real)
    centers = cell_centers(x)
    keep = findall(centers) do (lon, lat)
        lon_min <= lon <= lon_max && lat_min <= lat <= lat_max
    end
    return _select(x, DD.name(_dggs_dim(x)), keep)
end

function _dggs_dim(x::Union{YAXArray,DD.AbstractDimArray})
    for d in DD.dims(x)
        DD.lookup(d) isa DGG.AbstractDGGSLookup && return d
    end
    throw(ArgumentError("no DGGS-indexed dimension found on $(typeof(x))"))
end

function _dggs_dim(ds::YAXArrays.Dataset)
    for d in values(ds.axes)
        DD.lookup(d) isa DGG.AbstractDGGSLookup && return d
    end
    throw(ArgumentError("no DGGS-indexed dimension found in this dataset"))
end

# `Dataset` has no dimension-keyword `getindex`, so a dataset-wide selection is
# applied variable by variable and reassembled. Arrays take the keyword path.
_select(x::Union{YAXArray,DD.AbstractDimArray}, dimname::Symbol, selector) =
    getindex(x; (dimname => selector,)...)

function _select(ds::YAXArrays.Dataset, dimname::Symbol, selector)
    cubes = (name => _select(cube, dimname, selector) for (name, cube) in ds.cubes)
    return YAXArrays.Dataset(; properties=ds.properties, cubes...)
end

end # module DGGSZarr
