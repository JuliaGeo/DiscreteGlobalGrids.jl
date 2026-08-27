"""
    GlobalRegriddingProjExt

Task-local projected-coordinate charts for `GlobalRegridding.RasterGrid`.

Pass a `Proj.Transformation` from the raster's native CRS to geographic
longitude/latitude through the `native_to_unit_sphere` keyword. The template
must be constructed with `always_xy = true`, so its coordinate order agrees
with the raster's X/Y axes. The selected pipeline is cloned rather than rebuilt
from CRS metadata.
"""
module GlobalRegriddingProjExt

import GeometryOps as GO
import GlobalRegridding as GR
import Proj

const _GEOGRAPHIC_TO_UNIT_SPHERE = GO.UnitSpherical.UnitSphereFromGeographic()
const _UNIT_SPHERE_TO_GEOGRAPHIC = GO.UnitSpherical.GeographicFromUnitSphere()
const _TEMPLATE_CLONE_LOCK = ReentrantLock()

"""Shared description of one projected raster chart; no operational clone lives here."""
mutable struct _ProjChartState
    template::Proj.Transformation
    direction::Proj.PJ_DIRECTION
    clone_lock::ReentrantLock
    task_key::Symbol
end

function _ProjChartState(template::Proj.Transformation)
    template.pj == C_NULL && throw(ArgumentError(
        "cannot build a projected RasterGrid from a finalized Proj.Transformation"))
    return _ProjChartState(template, template.direction, _TEMPLATE_CLONE_LOCK,
        gensym(:GlobalRegriddingProjChart))
end

struct _NativeToUnitSphere
    state::_ProjChartState
end

struct _UnitSphereToNative
    state::_ProjChartState
end

struct _PreparedNativeToUnitSphere
    transformation::Proj.Transformation
end

"""Own both transformation clones and their context as one finalizable resource."""
mutable struct _ProjTaskPair{R,D}
    ctx::Ptr{Proj.PJ_CONTEXT}
    forward::Union{Nothing,Proj.Transformation}
    reverse::Union{Nothing,Proj.Transformation}
    closed::Bool
    release_transformation::R
    destroy_context::D
end

function _close_proj_pair!(pair::_ProjTaskPair)
    pair.closed && return pair

    pair.closed = true
    forward, reverse, ctx = pair.forward, pair.reverse, pair.ctx
    pair.forward = nothing
    pair.reverse = nothing
    pair.ctx = Ptr{Proj.PJ_CONTEXT}(C_NULL)

    # Both PJ objects belong to `ctx` and therefore must be gone before it is.
    try
        forward === nothing || pair.release_transformation(forward)
    finally
        try
            reverse === nothing || pair.release_transformation(reverse)
        finally
            ctx == C_NULL || pair.destroy_context(ctx)
        end
    end
    return pair
end

Base.close(pair::_ProjTaskPair) = _close_proj_pair!(pair)

function _new_task_pair(state::_ProjChartState;
    context_clone = Proj.proj_context_clone,
    clone_object = Proj.proj_clone,
    wrap_transformation = Proj.Transformation,
    release_transformation = finalize,
    destroy_context = Proj.proj_context_destroy,
    destroy_raw = Proj.proj_destroy)

    ctx = Ptr{Proj.PJ_CONTEXT}(C_NULL)
    raw = Ptr{Proj.PJ}(C_NULL)
    forward = reverse = nothing
    try
        ctx = context_clone()
        ctx == C_NULL && error("Proj failed to clone its default context")

        # The task will use only its clones. This brief lock protects cloning
        # from the mutable shared template; coordinate evaluation is lock-free.
        lock(state.clone_lock) do
            state.template.pj == C_NULL && throw(ArgumentError(
                "the projected RasterGrid's Proj.Transformation template was finalized"))

            raw = clone_object(state.template.pj, ctx)
            raw == C_NULL && error("Proj failed to clone the forward transformation")
            forward = wrap_transformation(raw, state.direction)
            raw = Ptr{Proj.PJ}(C_NULL) # ownership transferred to `forward`

            raw = clone_object(state.template.pj, ctx)
            raw == C_NULL && error("Proj failed to clone the reverse transformation")
            reverse = wrap_transformation(raw, inv(state.direction))
            raw = Ptr{Proj.PJ}(C_NULL) # ownership transferred to `reverse`
        end

        pair = _ProjTaskPair(
            ctx, forward, reverse, false, release_transformation, destroy_context)
        finalizer(_close_proj_pair!, pair)
        return pair
    catch
        # Release every PJ before its context. `destroy_raw` is needed only when
        # wrapping failed before a non-null clone acquired a Julia owner.
        try
            raw == C_NULL || destroy_raw(raw)
        finally
            try
                reverse === nothing || release_transformation(reverse)
            finally
                try
                    forward === nothing || release_transformation(forward)
                finally
                    ctx == C_NULL || destroy_context(ctx)
                end
            end
        end
        rethrow()
    end
end

function _task_pair(state::_ProjChartState; kwargs...)
    storage = task_local_storage()
    pair = get(storage, state.task_key, nothing)
    pair isa _ProjTaskPair && !pair.closed && return pair
    pair = _new_task_pair(state; kwargs...)
    storage[state.task_key] = pair
    return pair
end

function GR._task_prepared_raster_transform(chart::_NativeToUnitSphere)
    transformation = _task_pair(chart.state).forward
    transformation === nothing && error("the task-local forward Proj transformation is closed")
    return _PreparedNativeToUnitSphere(transformation)
end

@inline (chart::_PreparedNativeToUnitSphere)(xy) =
    _GEOGRAPHIC_TO_UNIT_SPHERE(chart.transformation(xy))

@inline function (chart::_NativeToUnitSphere)(xy)
    transformation = _task_pair(chart.state).forward
    transformation === nothing && error("the task-local forward Proj transformation is closed")
    return _GEOGRAPHIC_TO_UNIT_SPHERE(transformation(xy))
end

@inline function (chart::_UnitSphereToNative)(point)
    transformation = _task_pair(chart.state).reverse
    transformation === nothing && error("the task-local reverse Proj transformation is closed")
    return transformation(_UNIT_SPHERE_TO_GEOGRAPHIC(point))
end

function GR._prepare_raster_transform_pair(
    template::Proj.Transformation, raw_backward)
    state = _ProjChartState(template)
    forward = _NativeToUnitSphere(state)
    backward = raw_backward === GR._UNSET_RASTER_KEYWORD ?
        _UnitSphereToNative(state) :
        GR._one_coordinate_unit_sphere_to_native(raw_backward)
    return (forward, backward)
end

end # module GlobalRegriddingProjExt
