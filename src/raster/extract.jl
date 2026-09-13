"""
    DiscreteGlobalGrids.extract(A, geometries; boundary=:center, kw...)

Extract cell values as NamedTuple rows.

- Points use containing-cell lookup; lines intersect cells; polygons follow
  `boundary` (`:center`, `:intersects`/`:touches`, `:inside`).
- Rows carry `:geometry` (input point, or the cell representative lon/lat),
  `id=true` the feature number, `index=true` the local cell-axis position.
- Unsampled dimensions stay labelled slices; `name` selects stack layers.
- `skipmissing=true` drops missing rows — a row whose slice holds any missing
  element counts as missing; `flatten=false` groups non-point rows by feature.
"""
function extract(A::Union{DD.AbstractDimArray,DD.AbstractDimStack}, data; names=nothing, name=names,
        skipmissing=false, flatten=true, id=false, geometry=true, index=false, kw...)
    geoms, opts = _raster_inputs(data; kw...)
    layers = _raster_extract_layers(A, name)
    grid = _raster_grid(first(values(layers)))
    axis = DD.lookup(first(values(layers)), _raster_dimnum(first(values(layers))))
    all(l -> DD.lookup(l, _raster_dimnum(l)) == axis, values(layers)) ||
        throw(ArgumentError("Extraction stack layers must share the cell axis"))
    reserved = ((id ? (:id,) : ())..., (geometry ? (:geometry,) : ())..., (index ? (:index,) : ())...)
    any(in(reserved), keys(layers)) && throw(ArgumentError("Layer names conflict with requested extraction fields"))
    missingvals = map(_raster_missingval, layers)
    groups = _raster_map(eachindex(geoms), opts.threaded) do j
        g = geoms[j]
        indices = _raster_indices(grid, g; opts.boundary, opts.shape)
        # A point or absent geometry always yields one row, missing-valued when it hits nothing.
        onerow = isempty(indices) && (_raster_ispoint(g) || _raster_absent(g))
        rows = [_raster_row(layers, grid, g, j, i; id, geometry, index) for i in (onerow ? Union{Nothing,Int}[nothing] : indices)]
        skipmissing || return rows
        return filter(row -> !any(k -> _raster_extract_missing(row[k], missingvals[k]), keys(layers)), rows)
    end
    allpoints = all(g -> _raster_absent(g) || _raster_ispoint(g), geoms)
    flatten || allpoints || return groups
    return isempty(groups) ? NamedTuple[] : reduce(vcat, groups)
end

_raster_ispoint(g) = GI.trait(g) isa GI.AbstractPointTrait
_raster_absent(g) = g === nothing || g === missing
function _raster_row(layers, grid, g, j, i; id, geometry, index)
    vals = map(layer -> i === nothing ? missing : _raster_extract_value(layer, i), layers)
    coord = _raster_absent(g) ? missing : _raster_ispoint(g) ? _raster_point_tuple(g) :
        Fallbacks.lonlat(cell_centroid(grid, cellindex(grid, i)))
    return merge(id ? (; id=j) : (;), geometry ? (; geometry=coord) : (;),
        index ? (; index=something(i, missing)) : (;), vals)
end
_raster_point_tuple(p) = GI.ncoord(p) == 3 ? (GI.x(p), GI.y(p), GI.z(p)) : (GI.x(p), GI.y(p))

function _raster_extract_layers(A::DD.AbstractDimArray, name)
    n = name === nothing ? DD.name(A) : name
    n = n isa Tuple ? only(n) : n
    n = n isa Symbol && !isempty(string(n)) ? n : :layer
    return NamedTuple{(n,)}((A,))
end
function _raster_extract_layers(A::DD.AbstractDimStack, name)
    names = name === nothing ? keys(A) : name isa Symbol ? (name,) : Tuple(name)
    isempty(names) && throw(ArgumentError("Select at least one stack layer"))
    return NamedTuple{names}(Tuple(A[k] for k in names))
end
function _raster_extract_value(A, i)
    d = _raster_dimnum(A)
    return ndims(A) == 1 ? A[i] : view(A, ntuple(k -> k == d ? i : Colon(), ndims(A))...)
end
_raster_extract_missing(v::DD.AbstractDimArray, mv) = any(x -> _raster_ismissing(x, mv), v)
_raster_extract_missing(v, mv) = _raster_ismissing(v, mv)
