# ---------------------------------------------------------------------------
# GeometryOps manifolds ↔ the authalic transform
#
# A DGGS in this package is tessellated on a sphere: every chart, every cell
# boundary and every cap lives on `GO.UnitSpherical`. An ellipsoidal DGGS is
# that same tessellation read on the *authalic* sphere — the sphere of equal
# surface area — with latitudes warped back to geodetic only at the lon/lat
# boundary. So an ellipsoid enters this package as two derived objects:
#
#   * an `AuthalicTransform`, which does the latitude warp at the I/O edge, and
#   * a `Spherical` compute manifold at the authalic radius, which is what the
#     tree, the caps and `ConservativeRegridding` actually see.
#
# This file derives both from a `GeometryOpsCore.Manifold`, so the ellipsoid is
# declared once, in GO's vocabulary, and never restated. Two notes on why the
# split exists at all rather than just handing `Geodesic` to everything:
#
#   * `ConservativeRegridding` imports `Manifold, Planar, Spherical` and
#     nothing else; every `cell_range_extent` method dispatches on `Planar` or
#     `Spherical`. GO itself supports `Geodesic` only for `area`, `arclength`
#     and `segmentize`. A `Geodesic` compute manifold is a `MethodError` today.
#   * It is also the frame the cells are genuinely defined in. Areas computed
#     on the authalic sphere, scaled by `R_A²`, are exactly the ellipsoidal
#     areas; areas computed from geodetic-warped coordinates are not, because
#     the warp is what makes the grid equal-area in the first place.
# ---------------------------------------------------------------------------

"""
    AuthalicTransform(m::GeometryOpsCore.Manifold)
    AuthalicTransform{T}(m::GeometryOpsCore.Manifold)

Read the ellipsoid a GeometryOps manifold declares and build the corresponding
[`Helpers.AuthalicTransform`](@ref).

  - [`GeometryOpsCore.Geodesic`](@ref) supplies `semimajor_axis` and
    `inv_flattening` directly, which is the ellipsoid.
  - [`GeometryOpsCore.Spherical`](@ref) has no flattening, so it yields the
    degenerate `e² = 0` transform at that radius. Every series coefficient is
    then zero and both directions are *exactly* the identity, which is the
    point: a spherical grid and an ellipsoidal one run the same code path with
    no branch and no cost.
  - `Planar` and `AutoManifold` throw. Neither names an ellipsoid, and for a
    DGGS the manifold is part of the reference-system definition — guessing one
    is precisely the silent misregistration this adapter exists to prevent.

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

The compute manifold for `x`: the sphere on which its cell areas are the true
ellipsoidal areas. `x` may be a [`GeometryOpsCore.Manifold`](@ref) or a
[`Helpers.AuthalicTransform`](@ref).

A `Geodesic` resolves to `Spherical(; radius = R_A)`, the *authalic* radius. A
`Spherical` is already a compute manifold and is returned unchanged, so this is
idempotent and a spherical grid needs no special case.

!!! warning "Not `Spherical()`"
    GO's default `Spherical()` radius is `WGS84_EARTH_MEAN_RADIUS`
    (6371008.8 m), the *mean* radius `R₁ = (2a + b)/3`. That sphere is not
    area-preserving. For WGS84 the authalic radius is 6371007.180918474 m —
    1.6 m smaller, a ~5e-7 relative area error if the two are confused. This
    function exists so that no call site has to remember which is which; never
    reach for `Spherical()` bare.

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
