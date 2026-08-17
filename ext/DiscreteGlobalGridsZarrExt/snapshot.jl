# StoreSnapshot construction from a live Zarr store: the one place the
# convention layer's plain-data view is filled in from bytes.
#
# Nothing here interprets an attribute. The group's attrs and each array's
# attrs are copied verbatim, so a convention sees exactly what the producer
# wrote, and `_ARRAY_DIMENSIONS` is read but never invented.

using DiscreteGlobalGrids: StoreSnapshot, ArrayEntry, DGGSFormatError

"The `zarr-conventions/xarray` key both v2 and v3 stores carry dimension names in."
const ARRAY_DIMENSIONS = "_ARRAY_DIMENSIONS"

"The bucket-relative HTTPS endpoint every public Google Cloud Storage store has."
const GCS_BASE = "https://storage.googleapis.com/"

"""
    snapshot(g::ZGroup; identifier = storeidentifier(g)) -> StoreSnapshot

The metadata-only view of one Zarr group: its attributes verbatim, and one
[`ArrayEntry`](@ref) per array it holds, in name order.

This is the whole of what a [`DGGSConvention`](@ref) is shown, and the reason
convention logic needs no store to be tested against. Subgroups are not
followed: a snapshot describes one group.

Shapes are reported the way Zarr declares them, outermost dimension first, so
that `entry.shape[i]` is the extent of `entry.dims[i]`. Julia's own `size` on
the same array is the reverse of that.
"""
function snapshot(g::Zarr.ZGroup; identifier::AbstractString=storeidentifier(g))
    id = String(identifier)
    names = sort!(collect(keys(g.arrays)))
    entries = ArrayEntry[arrayentry(name, g.arrays[name], id) for name in names]
    return StoreSnapshot(identifier=id, attrs=Dict{String,Any}(g.attrs), arrays=entries)
end

"""
    arrayentry(name, z::ZArray, identifier) -> ArrayEntry

One array of a snapshot: its attributes copied, its declared shape, its element
type, and the dimension names it carries.
"""
function arrayentry(name::AbstractString, z::Zarr.ZArray, identifier::AbstractString)
    attrs = Dict{String,Any}(z.attrs)
    return ArrayEntry(name=String(name), attrs=attrs,
        shape=reverse(size(z)), eltype=eltype(z),
        dims=arraydimensions(name, attrs, ndims(z), identifier))
end

# A store that names no dimensions gets none: a synthesized name would look to
# a convention exactly like a declared one, and the DKRZ dialect's whole
# spatial-dimension search is over names that a producer chose.
function arraydimensions(name, attrs, n::Int, identifier)
    raw = get(attrs, ARRAY_DIMENSIONS, nothing)
    raw === nothing && return String[]
    raw isa AbstractVector || throw(DGGSFormatError(check=:invalid_array_dimensions,
        store=identifier, declared=name, observed=raw,
        detail="`$ARRAY_DIMENSIONS` on array `$name` is not a list of names."))
    dims = String[x isa AbstractString ? String(x) :
                  throw(DGGSFormatError(check=:invalid_array_dimensions,
        store=identifier, declared=name, observed=x,
        detail="`$ARRAY_DIMENSIONS` on array `$name` holds a name that is not a string."))
           for x in raw]
    length(dims) == n || throw(DGGSFormatError(check=:invalid_array_dimensions,
        store=identifier, declared=n, observed=dims,
        detail="array `$name` has $n dimensions but `$ARRAY_DIMENSIONS` names " *
               "$(length(dims)) of them."))
    return dims
end

"""
    storeidentifier(g::ZGroup) -> String

Where the group came from, as the URL or path a person would type — what a
[`DGGSFormatError`](@ref) reports and what the stack's provenance records.
"""
function storeidentifier(g::Zarr.ZGroup)
    base = rstrip(storeurl(g.storage), '/')
    return isempty(g.path) ? base : string(base, "/", g.path)
end

storeurl(s::Zarr.DirectoryStore) = s.folder
storeurl(s::Zarr.GCStore) = string(GCS_BASE, s.bucket)
storeurl(s::Zarr.HTTPStore) = s.url
storeurl(s::Zarr.ConsolidatedStore) = storeurl(s.parent)
storeurl(s::Zarr.AbstractStore) = string(s)
