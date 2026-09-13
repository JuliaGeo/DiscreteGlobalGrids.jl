# Raster membership: normalise the geometry, then hand the query engine a
# predicate — `boundary` is a choice of predicate, nothing more.

function _raster_boundary(boundary)
    boundary === :touches && return :intersects
    boundary in (:center, :inside, :intersects) || throw(ArgumentError(
        "boundary must be :center, :inside, :intersects, or :touches"))
    return boundary
end

_raster_predicate(boundary) =
    boundary === :center ? Engine.CentroidCovered() :
    boundary === :inside ? DE9IM.Within(nothing) : DE9IM.Intersects(nothing)

function _raster_shape(shape)
    shape === nothing && return nothing
    shape in (:point, :points) && return :point
    shape in (:line, :lines, :linestring, :linestrings) && return :line
    shape in (:polygon, :polygons) && return :polygon
    throw(ArgumentError("shape must be :point, :line, :polygon, or nothing"))
end

function _raster_ring_line(ring)
    points = collect(GI.getpoint(ring))
    isempty(points) || first(points) == last(points) || push!(points, first(points))
    return GI.LineString(points)
end

function _raster_shaped_geometry(geom, shape)
    shape = _raster_shape(shape)
    shape === nothing && return geom
    trait = GI.trait(geom)
    shape === :point && return trait isa GI.PointTrait ? geom : GI.MultiPoint(collect(GI.getpoint(geom)))
    if shape === :line
        trait isa GI.PolygonTrait && return GI.MultiLineString([
            _raster_ring_line(r) for r in GI.getring(geom)])
        trait isa GI.MultiPolygonTrait && return GI.MultiLineString([
            _raster_ring_line(r) for p in GI.getpolygon(geom)
            for r in GI.getring(p)])
        trait isa GI.MultiLineStringTrait && return geom
        return GI.LineString(collect(GI.getpoint(geom)))
    end
    trait isa Union{GI.PolygonTrait,GI.MultiPolygonTrait} && return geom
    return GI.Polygon([GI.LinearRing(collect(GI.getpoint(geom)))])
end

function _raster_empty_geometry(geom)
    (ismissing(geom) || geom === nothing) && return false
    trait = GI.trait(geom)
    trait === nothing && return false
    trait isa GI.FeatureTrait && return _raster_empty_geometry(GI.geometry(geom))
    if trait isa Union{GI.GeometryCollectionTrait,GI.FeatureCollectionTrait}
        parts = trait isa GI.FeatureCollectionTrait ? GI.getfeature(geom) : GI.getgeom(geom)
        return all(_raster_empty_geometry, parts)
    end
    GI.isempty(geom) && return true
    trait isa GI.PointTrait && return false
    # Flattening getpoint can fail for a polygon with no rings. Count its
    # points through the ring hierarchy before asking for coordinates.
    trait isa Union{GI.PolygonTrait,GI.MultiPolygonTrait} && return GI.npoint(geom) == 0
    return isempty(GI.getpoint(geom))
end
_raster_empty_geometry(::Union{Extents.Extent,GO.UnitSpherical.SphericalCap}) = false

function _raster_indices(grid::AbstractGrid, geom; boundary=:center, shape=nothing)
    boundary = _raster_boundary(boundary)
    shape = _raster_shape(shape)
    (ismissing(geom) || geom === nothing) && return Int[]
    trait = GI.trait(geom)
    trait isa GI.FeatureTrait && return _raster_indices(grid, GI.geometry(geom); boundary, shape)
    if trait isa Union{GI.GeometryCollectionTrait,GI.FeatureCollectionTrait}
        out = Int[]
        parts = trait isa GI.FeatureCollectionTrait ? GI.getfeature(geom) : GI.getgeom(geom)
        for part in parts
            append!(out, _raster_indices(grid, part; boundary, shape))
        end
        return sort!(unique!(out))
    end
    _raster_empty_geometry(geom) && return Int[]
    geom = _raster_shaped_geometry(geom, shape)
    trait = GI.trait(geom)
    if trait isa Union{GI.PointTrait,GI.MultiPointTrait}
        out = Int[]
        for point in (trait isa GI.PointTrait ? (geom,) : GI.getpoint(geom))
            i = localindex(grid, Fallbacks.query_point(point))
            i === nothing || push!(out, i)
        end
        return sort!(unique!(out))
    end
    trait isa GI.LinearRingTrait && (geom = _raster_ring_line(geom))
    # A line has no inside, so every boundary rule asks the same question of it.
    pred = trait isa Union{GI.PolygonTrait,GI.MultiPolygonTrait} ?
        _raster_predicate(boundary) : DE9IM.Intersects(nothing)
    return Engine._query_indices(grid, pred, Engine._query_target(geom))
end

# Extents use the query API's densified longitude/latitude outline. Polar
# full-longitude extents have a native spherical-cap representation instead.
function _raster_indices(grid::AbstractGrid, extent::Extents.Extent;
        boundary=:center, shape=nothing)
    shape === nothing || throw(ArgumentError("shape overrides require a GeoInterface geometry"))
    return _raster_indices(grid, Engine._extent_target(extent); boundary)
end

function _raster_indices(grid::AbstractGrid, cap::GO.UnitSpherical.SphericalCap;
        boundary=:center, shape=nothing)
    shape === nothing || throw(ArgumentError("shape overrides require a GeoInterface geometry"))
    return Engine._query_indices(grid, _raster_predicate(_raster_boundary(boundary)),
        Engine._query_target(cap))
end
