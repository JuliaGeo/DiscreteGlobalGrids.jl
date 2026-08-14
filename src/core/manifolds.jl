# ---------------------------------------------------------------------------
# GeometryOps manifold adapters. An ellipsoid supplies an `AuthalicTransform`
# for latitude conversion and a `Spherical` compute manifold at authalic radius
# `R_A`. Unit-sphere areas scaled by `R_A²` equal ellipsoidal areas.
#
# The compute manifold is never `Geodesic`: `ConservativeRegridding` and the
# `cell_range_extent` methods dispatch on `Planar`/`Spherical` only, and GO
# supports `Geodesic` for `area`/`arclength`/`segmentize` alone, so handing one
# down is a `MethodError` today.
# ---------------------------------------------------------------------------

"""
    AuthalicTransform(m::GeometryOpsCore.Manifold)
    AuthalicTransform{T}(m::GeometryOpsCore.Manifold)

Build an [`Helpers.AuthalicTransform`](@ref) from a GeometryOps manifold.

`Geodesic` supplies its semi-major axis and inverse flattening. `Spherical`
produces an `e² = 0` identity transform at its radius. `Planar`, `AutoManifold`,
and other manifolds that do not define an ellipsoid throw `ArgumentError`.

`T` defaults to the manifold's own element type, so a `Geodesic{Float32}` gives
an `AuthalicTransform{Float32}`.

See [`authalic_sphere`](@ref) for the companion direction, which produces the
compute manifold rather than the transform.
"""
function Helpers.AuthalicTransform{T}(m::GOCore.Geodesic) where {T<:AbstractFloat}
    return Helpers.AuthalicTransform{T}(;
        semimajor_axis=m.semimajor_axis,
        inverse_flattening=m.inv_flattening,
    )
end

function Helpers.AuthalicTransform{T}(m::GOCore.Spherical) where {T<:AbstractFloat}
    return Helpers.AuthalicTransform{T}(;
        semimajor_axis=m.radius,
        eccentricity_squared=0,
    )
end

Helpers.AuthalicTransform(m::GOCore.Geodesic{T}) where {T} =
    Helpers.AuthalicTransform{float(T)}(m)

Helpers.AuthalicTransform(m::GOCore.Spherical{T}) where {T} =
    Helpers.AuthalicTransform{float(T)}(m)

# `Planar` and `AutoManifold` are rejected with a reason rather than left to a
# `MethodError`, because the reason is the whole design argument and a caller
# who reaches here is one defaulted constructor away from a 14 km offset.
Helpers.AuthalicTransform{T}(m::GOCore.Manifold) where {T<:AbstractFloat} =
    throw(ArgumentError(lazy"""
    cannot derive an ellipsoid from $(typeof(m)).

    A DGGS declares its ellipsoid as part of its reference system, so it must \
    come from the manifold rather than be guessed. Use `Geodesic(; \
    semimajor_axis, inv_flattening)` for an ellipsoidal grid, or `Spherical(; \
    radius)` for a spherical one."""))

Helpers.AuthalicTransform(m::GOCore.Manifold) = Helpers.AuthalicTransform{Float64}(m)

"""
    authalic_sphere(x) -> GeometryOpsCore.Spherical

Return the spherical compute manifold whose cell areas equal the corresponding
ellipsoidal areas. `x` may be a [`GeometryOpsCore.Manifold`](@ref) or
[`Helpers.AuthalicTransform`](@ref).

A `Geodesic` resolves to `Spherical(; radius=R_A)`. A `Spherical` is returned
unchanged.

!!! warning "Not `Spherical()`"
    GO's default `Spherical()` uses the WGS84 mean radius, 6371008.8 m, not the
    area-preserving authalic radius, 6371007.180918474 m. Do not substitute a
    bare `Spherical()` for this result.

`Planar` and `AutoManifold` throw, for the reason given in
[`AuthalicTransform`](@ref).
"""
function authalic_sphere end

authalic_sphere(m::GOCore.Spherical) = m

authalic_sphere(t::Helpers.AuthalicTransform) =
    GOCore.Spherical(; radius=Helpers.authalic_radius(t))

authalic_sphere(m::GOCore.Geodesic) = authalic_sphere(Helpers.AuthalicTransform(m))

# Same reasoning as the `AuthalicTransform` fallback: name the problem.
authalic_sphere(m::GOCore.Manifold) = throw(ArgumentError(lazy"""
    cannot derive an authalic sphere from $(typeof(m)).

    A DGGS declares its ellipsoid as part of its reference system. Use \
    `Geodesic(; semimajor_axis, inv_flattening)` for an ellipsoidal grid, or \
    `Spherical(; radius)` for a spherical one."""))
