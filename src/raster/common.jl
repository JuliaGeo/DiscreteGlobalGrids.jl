import Tables
import Statistics

struct _RasterUnset end
const _RASTER_UNSET = _RasterUnset()

# Keep feature properties beside the geometry so normalization happens only once.
struct _RasterFeature{G,P}
    geom::G
    properties::P
end

function _raster_features(data; geometrycolumn=nothing)
    out = _RasterFeature[]
    _raster_features!(out, data, geometrycolumn)
    return out
end
function _raster_features!(out, data, column)
    if data === missing || data === nothing
        push!(out, _RasterFeature(data, NamedTuple()))
        return out
    end
    trait = GI.trait(data)
    if trait isa GI.AbstractFeatureTrait && !(data isa NamedTuple && Tables.istable(data))
        push!(out, _RasterFeature(GI.geometry(data), GI.properties(data)))
    elseif trait isa GI.AbstractFeatureCollectionTrait
        for feature in GI.getfeature(data)
            _raster_features!(out, feature, column)
        end
    elseif trait isa GI.AbstractGeometryTrait
        push!(out, _RasterFeature(data, NamedTuple()))
    elseif data isa DD.Dimension
        _raster_features!(out, parent(DD.lookup(data)), column)
    elseif Tables.istable(data)
        rows = Tables.rows(data)
        for row in rows
            cols = Tables.columnnames(row)
            gc = column === nothing ? (:geometry in cols ? :geometry : :geom in cols ? :geom : nothing) : column
            gc === nothing && throw(ArgumentError("Specify geometrycolumn for a table without a :geometry column"))
            geom = gc isa Tuple ? Tuple(Tables.getcolumn(row, c) for c in gc) : Tables.getcolumn(row, gc)
            push!(out, _RasterFeature(geom, NamedTuple{Tuple(cols)}(Tuple(Tables.getcolumn(row, c) for c in cols))))
        end
    elseif data isa AbstractArray || data isa Tuple || applicable(iterate, data)
        for item in data
            _raster_features!(out, item, column)
        end
    else
        throw(ArgumentError("Expected GeoInterface geometry/features or a Tables table, got $(typeof(data))"))
    end
    return out
end

# This adapter retains explicit mixed-level membership and storage order.
struct _RasterStoredGrid{S,V} <: AbstractGrid
    system::S
    ids::V
end
system(g::_RasterStoredGrid) = g.system
ncells(g::_RasterStoredGrid) = length(g.ids)
cellindex(g::_RasterStoredGrid, i::Int) = g.ids[i]
localindex(g::_RasterStoredGrid, c::AbstractCellIndex) = findfirst(==(c), g.ids)
cell_boundary(g::_RasterStoredGrid, c::AbstractCellIndex) = cell_boundary(levelgrid(system(g), level(c)), c)
cell_centroid(g::_RasterStoredGrid, c::AbstractCellIndex) = cell_centroid(levelgrid(system(g), level(c)), c)
cell_area(g::_RasterStoredGrid, c::AbstractCellIndex) = cell_area(levelgrid(system(g), level(c)), c)

# Generic grids and mixed levels cannot use the single-level CellLookup contract.
# Keep their grid with the dimension, including when taking labelled views.
struct _RasterGridLookup{T,G,V<:AbstractVector{T}} <: DD.Lookups.Lookup{T,1}
    grid::G
    data::V
end
Base.parent(l::_RasterGridLookup) = l.data
Base.size(l::_RasterGridLookup) = size(l.data)
Base.getindex(l::_RasterGridLookup, i::Int) = l.data[i]
Base.IndexStyle(::Type{<:_RasterGridLookup}) = IndexLinear()
DD.Lookups.order(::_RasterGridLookup) = DD.Lookups.Unordered()
DD.Lookups.sampling(::_RasterGridLookup) = DD.Lookups.Points()
DD.Lookups.metadata(::_RasterGridLookup) = DD.NoMetadata()
DD.Dimensions.format(l::_RasterGridLookup, ::Type, values, axis::AbstractRange) = l
function DD.Lookups.rebuild(l::_RasterGridLookup; data=parent(l), kw...)
    data === parent(l) && return l
    ids = [findfirst(isequal(c), parent(l)) for c in data]
    any(isnothing, ids) && throw(ArgumentError("Unknown cell in raster lookup"))
    return _RasterGridLookup(_RasterSubsetGrid(l.grid, Int[ids...]), data)
end
Base.getindex(l::_RasterGridLookup, i::AbstractVector{<:Integer}) = DD.Lookups.rebuild(l; data=parent(l)[i])
Base.getindex(l::_RasterGridLookup, ::Colon) = l
Base.view(l::_RasterGridLookup, i::AbstractVector{<:Integer}) = l[i]
Base.view(l::_RasterGridLookup, ::Colon) = l
struct _RasterSubsetGrid{G,V} <: AbstractGrid
    grid::G
    indices::V
end
ncells(g::_RasterSubsetGrid) = length(g.indices)
cellindex(g::_RasterSubsetGrid, i::Int) = cellindex(g.grid, g.indices[i])
localindex(g::_RasterSubsetGrid, c::AbstractCellIndex) = findfirst(i -> cellindex(g.grid,i)==c, g.indices)
cell_boundary(g::_RasterSubsetGrid, c::AbstractCellIndex) = cell_boundary(g.grid,c)
cell_centroid(g::_RasterSubsetGrid, c::AbstractCellIndex) = cell_centroid(g.grid,c)
cell_area(g::_RasterSubsetGrid, c::AbstractCellIndex) = cell_area(g.grid,c)

function _raster_celldim(g::AbstractGrid)
    if system(g) !== nothing && level(g) !== nothing
        return Cells(CellLookup(g))
    end
    return Cells(_RasterGridLookup(g, [cellindex(g,i) for i in 1:ncells(g)]))
end
_raster_grid(l::AbstractCellLookup) = PartialGrid(parent(l))
_raster_grid(l::_RasterGridLookup) = l.grid
function _raster_grid(d::DD.Dimension)
    l = DD.lookup(d)
    (l isa AbstractCellLookup || l isa _RasterGridLookup) || throw(ArgumentError("The spatial dimension must retain a DGGS cell lookup"))
    return _raster_grid(l)
end
function _raster_dimnum(A)
    found = findall(d -> DD.lookup(d) isa Union{AbstractCellLookup,_RasterGridLookup}, DD.dims(A))
    length(found) == 1 || throw(ArgumentError("Expected exactly one DGGS cell dimension, found $(length(found))"))
    return only(found)
end
function _raster_target(to; level=nothing)
    level === nothing || to isa AbstractHierarchicalGridSystem || throw(ArgumentError("level is only valid with a system destination"))
    if to isa DD.AbstractDimArray
        d = _raster_dimnum(to)
        return (_raster_grid(DD.dims(to)[d]), DD.dims(to), to)
    elseif to isa DD.AbstractDimStack
        isempty(keys(to)) && throw(ArgumentError("An empty stack has no target grid"))
        return _raster_target(to[first(keys(to))]; level)
    elseif to isa AbstractHierarchicalGridSystem
        level === nothing && throw(ArgumentError("A system destination needs level=..."))
        return _raster_target(levelgrid(to, level))
    elseif to isa AbstractGrid
        level === nothing || throw(ArgumentError("level is only valid with a system destination"))
        return (to, (_raster_celldim(to),), nothing)
    elseif to isa MultiOrderCellSet
        level === nothing || throw(ArgumentError("Explicit multi-order cells retain their levels; expand with CellVector first"))
        return _raster_target(_RasterStoredGrid(system(to), to.cells))
    elseif to isa AbstractCellVector
        return (PartialGrid(to), (Cells(to isa CellVector ? CellLookup(to) : ChunkedCellLookup(to)),), nothing)
    elseif to isa Union{AbstractCellLookup,_RasterGridLookup}
        return (_raster_grid(to), (Cells(to),), nothing)
    elseif to isa DD.Dimension
        return (_raster_grid(to), (to,), nothing)
    elseif to isa Tuple && all(d -> d isa DD.Dimension, to)
        k = findall(d -> DD.lookup(d) isa Union{AbstractCellLookup,_RasterGridLookup}, to)
        length(k)==1 || throw(ArgumentError("Target dimensions need exactly one DGGS cell lookup"))
        return (_raster_grid(to[only(k)]), to, nothing)
    end
    throw(ArgumentError("to must name a DGGS grid, cell lookup, dimensions, or dimensional array"))
end

_raster_missingval(A) = missing
_raster_ismissing(x, sentinel) = ismissing(x) || (sentinel !== nothing && isequal(x, sentinel))
function _raster_rebuild(template, data, dims; name=nothing, metadata=nothing,
        crs=nothing, mappedcrs=nothing, missingval=_RASTER_UNSET)
    (crs === nothing && mappedcrs === nothing) || throw(ArgumentError("DGGS cell axes have intrinsic spherical coordinates; transform geometries to longitude/latitude before rasterization instead of setting crs/mappedcrs"))
    if template === nothing
        return DD.DimArray(data, dims; name=something(name, :layer), metadata=something(metadata, DD.NoMetadata()))
    end
    return DD.rebuild(template; data, dims, name=something(name, DD.name(template)),
        metadata=something(metadata, DD.metadata(template)))
end
function _raster_stack(template, layers::NamedTuple; metadata=nothing)
    labelled = template isa Union{DD.AbstractDimArray,DD.AbstractDimStack}
    meta = metadata === nothing ? (labelled ? DD.metadata(template) : DD.NoMetadata()) : metadata
    refs = labelled ? DD.refdims(template) : ()
    return DD.DimStack(layers;metadata=meta,refdims=refs)
end

function _raster_check_keywords(kw)
    isempty(kw) || throw(ArgumentError("Unsupported keyword(s): $(join(keys(kw), ", "))"))
end

function _raster_point_tuple(p)
    GI.trait(p) isa GI.AbstractPointTrait || return missing
    return GI.ncoord(p) == 3 ? (GI.x(p),GI.y(p),GI.z(p)) : (GI.x(p),GI.y(p))
end
function _raster_center(g,i)
    p = cell_centroid(g,cellindex(g,i))
    return (rad2deg(atan(p[2],p[1])),rad2deg(asin(clamp(p[3],-1,1))))
end

treeify(g::_RasterStoredGrid) = Engine.IndexTreeNode(Engine.IndexTree(g),1)

Base.getindex(l::_RasterGridLookup, i::SmallCollections.AbstractFixedOrSmallOrPackedVector{<:Integer}) = DD.Lookups.rebuild(l; data=parent(l)[i])

_raster_issingle(data) = GI.trait(data) isa Union{GI.AbstractGeometryTrait,GI.AbstractFeatureTrait} && !(data isa NamedTuple && Tables.istable(data))

function _raster_features!(out, region::Union{Extents.Extent,GO.UnitSpherical.SphericalCap}, column)
    push!(out,_RasterFeature(region,NamedTuple()))
    return out
end
_raster_issingle(::Union{Extents.Extent,GO.UnitSpherical.SphericalCap}) = true

function _raster_cached_selections!(cache,A,features;boundary,shape)
    lookup=DD.lookup(A,_raster_dimnum(A))
    return get!(cache,lookup) do
        grid,_,_=_raster_target(A)
        [_raster_indices(grid,f.geom;boundary,shape) for f in features]
    end
end
