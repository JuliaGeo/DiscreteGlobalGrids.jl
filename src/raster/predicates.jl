# GeometryOps 0.1.45 can reach a planar XY crossing calculation from spherical
# RelateNG and divide by zero for a singular projection. Retry in the other
# coordinate planes until upstream uses a spherical crossing calculation.
# The cyclic permutation is an exact proper 3-D rotation: it preserves every
# coordinate, orientation, great-circle arc, and relation without perturbation.
function _raster_rotate_geometry(geom)
    return GO.apply(GI.PointTrait(), geom) do point
        GO.UnitSphericalPoint(GI.y(point), GI.z(point), GI.x(point))
    end
end

function _raster_relate(target, predicate, geom)
    try
        return GO.relate_predicate(target.prepared, predicate(), geom)
    catch exception
        exception isa DivideError || rethrow()
    end
    rotated_target = target.geom
    rotated_geom = Engine._to_unit_sphere(geom)
    for attempt in 1:2
        rotated_target = _raster_rotate_geometry(rotated_target)
        rotated_geom = _raster_rotate_geometry(rotated_geom)
        prepared = GO.prepare(GO.RelateNG(; manifold=GO.Spherical()), rotated_target)
        try
            # RelateNG mutates predicate state, so each retry needs a new one.
            return GO.relate_predicate(prepared, predicate(), rotated_geom)
        catch exception
            (exception isa DivideError && attempt < 2) || rethrow()
        end
    end
end
