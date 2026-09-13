import Tables
import Statistics

struct _RasterUnset end
const _RASTER_UNSET = _RasterUnset()

# --- geometry and fill inputs ---------------------------------------------

# A vector of NamedTuple points is rows to Tables but geometries here; a
# NamedTuple of columns carries a GeoInterface feature trait but is a table.
_raster_istable(data) = Tables.istable(data) && !(data isa AbstractVector{<:GI.NamedTuplePoint})
_raster_isfeature(data) = GI.trait(data) isa GI.AbstractFeatureTrait && !_raster_istable(data)
_raster_issingle(data) = _raster_isfeature(data) || GI.trait(data) isa GI.AbstractGeometryTrait ||
    data isa Union{Extents.Extent,GO.UnitSpherical.SphericalCap}

"""
    _raster_geometries(data; geometrycolumn=nothing) -> Vector

One geometry per feature, in input order: a geometry, feature, feature
collection, extent, spherical cap, `missing`/`nothing`, a Tables table
(`geometrycolumn` names the column or a coordinate-column tuple), or any
iterable of these.
"""
_raster_geometries(data; geometrycolumn=nothing) = _raster_column(data, nothing; geometrycolumn)

# One walk serves geometries (`name === nothing`) and fill columns, so they align.
function _raster_column(data, name; geometrycolumn=nothing)
    if _raster_istable(data)
        cols = Tables.columns(data)
        name === nothing || return collect(Tables.getcolumn(cols, name))
        names = Tables.columnnames(cols)
        gc = geometrycolumn === nothing ? (:geometry in names ? :geometry : :geom in names ? :geom : nothing) : geometrycolumn
        gc === nothing && throw(ArgumentError("Specify geometrycolumn for a table without a :geometry column"))
        gc isa Tuple && return map(tuple, (Tables.getcolumn(cols, c) for c in gc)...)
        return collect(Tables.getcolumn(cols, gc))
    end
    (_raster_issingle(data) || data === missing || data === nothing) && return [_raster_value(data, name)]
    GI.trait(data) isa GI.AbstractFeatureCollectionTrait && return [_raster_value(f, name) for f in GI.getfeature(data)]
    applicable(iterate, data) || throw(ArgumentError("Expected GeoInterface geometry/features or a Tables table, got $(typeof(data))"))
    parts = [_raster_column(item, name) for item in data]
    return isempty(parts) ? Union{}[] : reduce(vcat, parts)
end
_raster_value(f, ::Nothing) = _raster_isfeature(f) ? GI.geometry(f) : f
function _raster_value(f, name::Symbol)
    _raster_isfeature(f) || throw(ArgumentError("fill=:$name needs table columns or feature properties"))
    p = GI.properties(f)
    return p isa AbstractDict ? p[name] : getproperty(p, name)
end

# Fill values resolve to one value per geometry, or stay a function updating each cell.
function _raster_fills(data, fill, n)
    fill isa Function && return fill
    fill isa Symbol && return _raster_column(data, fill)
    _raster_issingle(data) && return [fill]
    if !(fill isa Union{Number,AbstractString,Missing,Nothing}) && applicable(iterate, fill)
        values = collect(fill)
        length(values) == 1 && return [only(values) for _ in 1:n]
        length(values) == n || throw(DimensionMismatch("fill must have one value per feature"))
        return values
    end
    return [fill for _ in 1:n]
end

_raster_layerkw(x, key) = x isa NamedTuple ? getproperty(x, key) : x
_raster_layermissing(A, missingval, key) = missingval isa _RasterUnset ? _raster_missingval(A) : _raster_layerkw(missingval, key)
_raster_copy(x) = isbits(x) ? x : deepcopy(x)

# --- keywords ---------------------------------------------------------------

# The selection keywords every verb shares; anything else is outside the in-memory API.
function _raster_options(; boundary=:center, shape=nothing, geometrycolumn=nothing, threaded=false,
        progress=true, verbose=true, crs=nothing, mappedcrs=nothing, kw...)
    isempty(kw) || throw(ArgumentError("Unsupported keyword(s): $(join(keys(kw), ", ")); file output, res, size, and chunked execution are outside the in-memory API"))
    (crs === nothing && mappedcrs === nothing) || throw(ArgumentError("DGGS cell axes have intrinsic spherical coordinates; transform geometries to longitude/latitude before rasterization instead of setting crs/mappedcrs"))
    return (; boundary=_raster_boundary(boundary), shape=_raster_shape(shape), geometrycolumn, threaded)
end

# --- targets ----------------------------------------------------------------

# Stored cells keep their own order and, for mixed levels, their own levels;
# `base` is the system those cells belong to, or a grid without a level.
struct _RasterStoredGrid{B,V} <: AbstractGrid
    base::B
    ids::V
end
system(g::_RasterStoredGrid) = g.base isa AbstractHierarchicalGridSystem ? g.base : system(g.base)
ncells(g::_RasterStoredGrid) = length(g.ids)
cellindex(g::_RasterStoredGrid, i::Int) = g.ids[i]
localindex(g::_RasterStoredGrid, c::AbstractCellIndex) = findfirst(==(c), g.ids)
_raster_base(g::_RasterStoredGrid, c) = g.base isa AbstractHierarchicalGridSystem ? levelgrid(g.base, level(c)) : g.base
for f in (:cell_boundary, :cell_centroid, :cell_area)
    @eval $f(g::_RasterStoredGrid, c::AbstractCellIndex) = $f(_raster_base(g, c), c)
end
treeify(g::_RasterStoredGrid) = Engine.IndexTreeNode(Engine.IndexTree(g), 1)

# Its cell axis keeps the grid so labelled views round-trip back to a target.
struct _RasterStoredLookup{T,G} <: DD.Lookups.Lookup{T,1}
    grid::G
    data::Vector{T}
end
_RasterStoredLookup(g::_RasterStoredGrid) = _RasterStoredLookup(g, collect(g.ids))
Base.parent(l::_RasterStoredLookup) = l.data
Base.size(l::_RasterStoredLookup) = size(l.data)
Base.getindex(l::_RasterStoredLookup, i::Int) = l.data[i]
Base.IndexStyle(::Type{<:_RasterStoredLookup}) = IndexLinear()
DD.Lookups.order(::_RasterStoredLookup) = DD.Lookups.Unordered()
DD.Lookups.sampling(::_RasterStoredLookup) = DD.Lookups.Points()
DD.Lookups.metadata(::_RasterStoredLookup) = DD.NoMetadata()
DD.Dimensions.format(l::_RasterStoredLookup, ::Type, values, axis::AbstractRange) = l
function DD.Lookups.rebuild(l::_RasterStoredLookup; data=parent(l), kw...)
    data === parent(l) && return l
    return _RasterStoredLookup(_RasterStoredGrid(l.grid.base, collect(data)))
end
for f in (:getindex, :view)
    @eval Base.$f(l::_RasterStoredLookup, ::Colon) = l
    @eval Base.$f(l::_RasterStoredLookup, i::AbstractVector{<:Integer}) = DD.Lookups.rebuild(l; data=parent(l)[i])
end
# SmallCollections indexes every AbstractVector with its own vectors; this breaks the tie.
Base.getindex(l::_RasterStoredLookup, i::SmallCollections.AbstractFixedOrSmallOrPackedVector{<:Integer}) = l[collect(i)]

const _RasterCellLookup = Union{AbstractCellLookup,_RasterStoredLookup}

_raster_grid(l::AbstractCellLookup) = PartialGrid(parent(l))
_raster_grid(l::_RasterStoredLookup) = l.grid
function _raster_grid(d::DD.Dimension)
    l = DD.lookup(d)
    l isa _RasterCellLookup || throw(ArgumentError("The spatial dimension must retain a DGGS cell lookup"))
    return _raster_grid(l)
end
_raster_grid(A::DD.AbstractDimArray) = _raster_grid(DD.dims(A)[_raster_dimnum(A)])
function _raster_dimnum(dims::Tuple)
    found = findall(d -> DD.lookup(d) isa _RasterCellLookup, dims)
    length(found) == 1 || throw(ArgumentError("Expected exactly one DGGS cell dimension, found $(length(found))"))
    return only(found)
end
_raster_dimnum(A::DD.AbstractDimArray) = _raster_dimnum(DD.dims(A))
_raster_celldim(g::_RasterStoredGrid) = Cells(_RasterStoredLookup(g))
_raster_celldim(g::AbstractGrid) = level(g) === nothing ?
    _raster_celldim(_RasterStoredGrid(g, [cellindex(g, i) for i in 1:ncells(g)])) : Cells(CellLookup(g))

"""
    _raster_target(to; level=nothing) -> (grid, dims, template)

The grid selections run on, the output dimensions, and the array to rebuild
from (`nothing` when `to` is not an array).
"""
function _raster_target(to; level=nothing)
    level === nothing && return _raster_destination(to)
    to isa AbstractHierarchicalGridSystem || throw(ArgumentError("level is only valid with a system destination"))
    return _raster_destination(levelgrid(to, level))
end
_raster_destination(::AbstractHierarchicalGridSystem) = throw(ArgumentError("A system destination needs level=..."))
_raster_destination(A::DD.AbstractDimArray) = (_raster_grid(A), DD.dims(A), A)
function _raster_destination(st::DD.AbstractDimStack)
    isempty(keys(st)) && throw(ArgumentError("An empty stack has no target grid"))
    return _raster_destination(st[first(keys(st))])
end
_raster_destination(g::AbstractGrid) = (g, (_raster_celldim(g),), nothing)
_raster_destination(set::MultiOrderCellSet) = _raster_destination(_RasterStoredGrid(set.system, set.cells))
_raster_destination(cv::AbstractCellVector) = _raster_destination(cv isa CellVector ? CellLookup(cv) : ChunkedCellLookup(cv))
_raster_destination(l::_RasterCellLookup) = _raster_destination((Cells(l),))
_raster_destination(d::DD.Dimension) = _raster_destination((d,))
_raster_destination(dims::Tuple{Vararg{DD.Dimension}}) = (_raster_grid(dims[_raster_dimnum(dims)]), dims, nothing)
_raster_destination(to) = throw(ArgumentError("to must name a DGGS grid, cell lookup, dimensions, or dimensional array"))

# --- selections -------------------------------------------------------------

_raster_selections(grid, geoms; boundary, shape) =
    Vector{Int}[_raster_indices(grid, g; boundary, shape) for g in geoms]

# Per-cell incidence in CSR form: feature ids in input order, no per-cell vectors.
struct _RasterIncidence
    offsets::Vector{Int}
    features::Vector{Int}
end
function _RasterIncidence(ncells::Int, selections::Vector{Vector{Int}})
    offsets = zeros(Int, ncells + 1)
    for sel in selections, j in sel
        offsets[j + 1] += 1
    end
    cumsum!(offsets, offsets)
    offsets .+= 1
    features = Vector{Int}(undef, offsets[end] - 1)
    cursor = copy(offsets)
    for (i, sel) in enumerate(selections), j in sel
        features[cursor[j]] = i
        cursor[j] += 1
    end
    return _RasterIncidence(offsets, features)
end
_raster_hits(inc::_RasterIncidence, j::Int) = @view inc.features[inc.offsets[j]:inc.offsets[j+1]-1]
_raster_touched(inc::_RasterIncidence, j::Int) = inc.offsets[j+1] > inc.offsets[j]

# Coverage along the cell axis, shaped to broadcast over the whole array.
function _raster_covered(A, selections, invert)
    axis = _raster_dimnum(A)
    covered = falses(size(A, axis))
    for sel in selections, j in sel
        covered[j] = true
    end
    return reshape(covered .⊻ invert, ntuple(i -> i == axis ? length(covered) : 1, ndims(A)))
end

# Stack layers sharing a cell axis select once; `f(key, layer, selections)` builds each layer.
function _raster_eachlayer(f, st::DD.AbstractDimStack, geoms; boundary, shape, keys=keys(st))
    cache = Dict{Any,Vector{Vector{Int}}}()
    layers = map(keys) do key
        layer = st[key]
        selections = get!(cache, DD.lookup(layer, _raster_dimnum(layer))) do
            _raster_selections(_raster_grid(layer), geoms; boundary, shape)
        end
        f(key, layer, selections)
    end
    return NamedTuple{keys}(layers)
end

# In-place work over disjoint indices, and typed results per chunk; the
# concatenation is the one runtime-typed step of the threaded `map`.
function _raster_foreach(f, indices, threaded)
    if threaded && Threads.nthreads() > 1
        Threads.@threads :dynamic for k in eachindex(indices)
            f(indices[k])
        end
    else
        foreach(f, indices)
    end
    return nothing
end
function _raster_map(f, xs, threaded)
    nchunks = threaded ? min(Threads.nthreads(), length(xs)) : 1
    nchunks > 1 || return map(f, xs)
    chunks = Iterators.partition(xs, cld(length(xs), nchunks))
    tasks = [Threads.@spawn(map(f, chunk)) for chunk in chunks]
    return reduce(vcat, map(t -> fetch(t)::AbstractVector, tasks))
end

# --- outputs ----------------------------------------------------------------

_raster_missingval(A) = missing
# The output missing value: the keyword, else the template's, else `missing`.
function _raster_missing(missingval, template)
    mv = missingval isa _RasterUnset ? (template === nothing ? missing : _raster_missingval(template)) : missingval
    return mv === nothing ? missing : mv
end
_raster_ismissing(x, sentinel) = ismissing(x) || (sentinel !== nothing && isequal(x, sentinel))
function _raster_rebuild(template, data, dims; name=nothing, metadata=nothing, missingval=_RASTER_UNSET)
    template === nothing && return DD.DimArray(data, dims; name=something(name, :layer), metadata=something(metadata, DD.NoMetadata()))
    return DD.rebuild(template; data, dims, name=something(name, DD.name(template)),
        metadata=something(metadata, DD.metadata(template)))
end
function _raster_stack(template, layers::NamedTuple; metadata=nothing)
    labelled = template isa Union{DD.AbstractDimArray,DD.AbstractDimStack}
    meta = metadata === nothing ? (labelled ? DD.metadata(template) : DD.NoMetadata()) : metadata
    return DD.DimStack(layers; metadata=meta, refdims=labelled ? DD.refdims(template) : ())
end
