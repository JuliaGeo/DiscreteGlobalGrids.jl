# A `GeoAxis`'s transform function is a `Proj.Transformation`, so this extension
# is what makes `dggpoly` work in one.  Two things are needed from Proj that a
# generic callable cannot give: where the map's cut meridian is, and a way to
# project a few million points without paying a `ccall` for each of them.

module ProjExt

import DiscreteGlobalGridsVisualization as DGGV
using DiscreteGlobalGridsVisualization: PlanarTarget, plot_target, project!
using GeometryBasics: Point2d
import Proj

"""
    plot_target(tf::Proj.Transformation) -> PlanarTarget

A `GeoAxis`'s projection, with its cut meridian read out of the projection
itself.
"""
DGGV.plot_target(tf::Proj.Transformation) = PlanarTarget(tf, cut_longitude(tf))

"""
    cut_longitude(tf::Proj.Transformation) -> Float64

The longitude of the seam in a projection: the meridian opposite its central
meridian.

The central meridian is `lon_0` in the projection's own definition, which PROJ
will print for us.  When it will not — a transformation between two full CRSs
describes itself as a pipeline whose definition names no single `lon_0` — the
value is recovered by asking where the origin of the projected plane came from,
which is the central meridian for every projection whose plane is centred on it.
Failing both, the antimeridian is assumed.
"""
function cut_longitude(tf::Proj.Transformation)
    lon0 = definition_lon_0(tf)
    lon0 === nothing && (lon0 = probe_lon_0(tf))
    lon0 === nothing && return 180.0
    return lon0 + 180.0
end

function definition_lon_0(tf::Proj.Transformation)
    definition = try
        unsafe_string(Proj.proj_pj_info(tf.pj).definition)
    catch
        return nothing
    end
    # A pipeline lists the steps in order, and the projection is the last one to
    # touch longitude, so the last `lon_0` is the one that sets the seam.
    match_last = nothing
    for m in eachmatch(r"\blon_0=([-+]?[0-9.]+)", definition)
        match_last = m
    end
    match_last === nothing && return nothing
    return tryparse(Float64, match_last.captures[1])
end

function probe_lon_0(tf::Proj.Transformation)
    origin = try
        Proj.proj_trans(tf.pj, Proj.PJ_INV, (0.0, 0.0))
    catch
        return nothing
    end
    lon = origin.x
    isfinite(lon) && -180.0 <= lon <= 180.0 || return nothing
    return lon
end

"""
    project!(target::PlanarTarget{<:Proj.Transformation}, positions)

Project the whole vertex buffer with a single strided `proj_trans_generic`.

`positions` holds `(lon, lat)` pairs contiguously, which is exactly the strided
layout PROJ's bulk entry point wants, so the buffer is transformed in place with
one call instead of one call per vertex.  A `PJ` object is not safe to share
between threads, so this stays serial — it is a tight loop inside PROJ rather
than a few million trips across the `ccall` boundary.
"""
function DGGV.project!(target::PlanarTarget{<:Proj.Transformation}, positions::Vector{Point2d})
    n = length(positions)
    n == 0 && return positions
    tf = target.projection
    stride = sizeof(Point2d)
    GC.@preserve positions begin
        base = Ptr{Float64}(pointer(positions))
        Proj.proj_trans_generic(
            tf.pj, tf.direction,
            base, stride, n,
            base + sizeof(Float64), stride, n,
            C_NULL, 0, 0,
            C_NULL, 0, 0,
        )
    end
    # PROJ reports an unprojectable point as ±HUGE_VAL.  Makie reads an infinite
    # vertex as a real one and lets it swallow the axis limits, so mark them the
    # way Makie already knows to ignore.
    Threads.@threads for i in eachindex(positions)
        @inbounds p = positions[i]
        if !(isfinite(p[1]) && isfinite(p[2]))
            @inbounds positions[i] = Point2d(NaN, NaN)
        end
    end
    return positions
end

end # module
