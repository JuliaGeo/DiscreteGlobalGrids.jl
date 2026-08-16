# Authalic latitude conversion and radius for an ellipsoid of revolution.
# Authalic latitude preserves area on the sphere of radius `R_A`.
#
# The order-6 forward and inverse series use Karney (2024), eqs. (A19–A20),
# with exact rational coefficients. `authalic_q` and `authalic_radius` use
# Snyder (1987), eqs. (3-12–3-13).
#
# Unsuffixed functions use radians; `d`-suffixed functions use degrees and
# preserve the poles and equator exactly.

# ---------------------------------------------------------------------------
# Series order and coefficient tables
# ---------------------------------------------------------------------------

"""
    AUTHALIC_SERIES_ORDER

Truncation order of the auxiliary-latitude series, whose expansion parameter is
the third flattening `n = (a − b)/(a + b) = f/(2 − f)`. For WGS84, neglected
terms are below `Float64` roundoff.
"""
const AUTHALIC_SERIES_ORDER = 6

"""
    AUTHALIC_FORWARD_TABLE

Exact rational coefficients from Karney (2024), eq. (A19), for geodetic to
authalic latitude. Row `j` is the ascending-power polynomial in third flattening
`n`, with `C_j = n · (row j)(n) = O(nʲ)`.
"""
const AUTHALIC_FORWARD_TABLE = (
    (-4//3, -4//45, 88//315, 538//4725, 20824//467775, -44732//2837835),
    (0//1, 34//45, 8//105, -2482//14175, -37192//467775, -12467764//212837625),
    (0//1, 0//1, -1532//2835, -898//14175, 54968//467775, 100320856//1915538625),
    (0//1, 0//1, 0//1, 6007//14175, 24496//467775, -5884124//70945875),
    (0//1, 0//1, 0//1, 0//1, -23356//66825, -839792//19348875),
    (0//1, 0//1, 0//1, 0//1, 0//1, 570284222//1915538625),
)

"""
    AUTHALIC_INVERSE_TABLE

Exact rational coefficients from Karney (2024), eq. (A20), for authalic to
geodetic latitude, laid out as [`AUTHALIC_FORWARD_TABLE`](@ref).
"""
const AUTHALIC_INVERSE_TABLE = (
    (4//3, 4//45, -16//35, -2582//14175, 60136//467775, 28112932//212837625),
    (0//1, 46//45, 152//945, -11966//14175, -21016//51975, 251310128//638512875),
    (0//1, 0//1, 3044//2835, 3802//14175, -94388//66825, -8797648//10945935),
    (0//1, 0//1, 0//1, 6059//4725, 41072//93555, -1472637812//638512875),
    (0//1, 0//1, 0//1, 0//1, 768272//467775, 455935736//638512875),
    (0//1, 0//1, 0//1, 0//1, 0//1, 4210684958//1915538625),
)

# ---------------------------------------------------------------------------
# Ellipsoid shape
# ---------------------------------------------------------------------------

"""
    WGS84_SEMIMAJOR_AXIS

WGS84 semi-major axis `a = 6378137.0` m (defining constant of the datum).
"""
const WGS84_SEMIMAJOR_AXIS = 6378137.0

"""
    WGS84_INVERSE_FLATTENING

WGS84 inverse flattening `1/f = 298.257223563` (defining constant of the
datum; `f` and `e²` are derived from it, never quoted independently).
"""
const WGS84_INVERSE_FLATTENING = 298.257223563

"""
    EllipsoidShapeError(count)

Thrown when an [`AuthalicTransform`](@ref) is asked for from anything other
than exactly one shape parameter. Carries only the count; the message is
built lazily in `showerror`.
"""
struct EllipsoidShapeError <: Exception
    count::Int
end

function Base.showerror(io::IO, err::EllipsoidShapeError)
    print(io, "AuthalicTransform needs exactly one of `flattening`, ",
        "`inverse_flattening` or `eccentricity_squared`; got ", err.count,
        ". They are different parameterizations of the same shape and ",
        "silently mixing them is the classic ellipsoid bug.")
end

# ---------------------------------------------------------------------------
# The transform
# ---------------------------------------------------------------------------

"""
    AuthalicTransform{T}

Precomputed forward and inverse authalic-latitude series and authalic radius for
one ellipsoid of revolution. Latitude conversions evaluate the stored
coefficients without recomputing ellipsoid constants.

# Constructors

    AuthalicTransform{T}(; semimajor_axis, flattening)
    AuthalicTransform{T}(; semimajor_axis, inverse_flattening)
    AuthalicTransform{T}(; semimajor_axis, eccentricity_squared)

`T` and `semimajor_axis` default to `Float64` and the WGS84 semi-major axis;
`semimajor_axis` affects only [`authalic_radius`](@ref), as the latitude
transforms are scale-free. Exactly one shape parameter is required; otherwise
construction throws [`EllipsoidShapeError`](@ref). Coefficients are evaluated in at least `Float64`
before conversion to `T`.

# Fields

  - `semimajor_axis` — `a`.
  - `eccentricity_squared` — `e² = f(2 − f)`, the canonical shape field.
  - `authalic_radius` — `R_A`, precomputed (see [`authalic_radius`](@ref)).
  - `fwd`, `inv` — the geodetic→authalic and authalic→geodetic sine-series
    coefficients `C_j`, in radians, `j = 1…AUTHALIC_SERIES_ORDER`.

```julia
julia> geodetic_to_authalicd(WGS84_AUTHALIC, 45.0)
44.871702873433939
```
"""
struct AuthalicTransform{T<:AbstractFloat}
    semimajor_axis::T
    eccentricity_squared::T
    authalic_radius::T
    fwd::NTuple{AUTHALIC_SERIES_ORDER,T}
    inv::NTuple{AUTHALIC_SERIES_ORDER,T}
end

"""
    _authalic_series(::Type{T}, n, table) -> NTuple{AUTHALIC_SERIES_ORDER,T}

Evaluate one coefficient table at third flattening `n` by Horner, returning
`C_j = n · (row j)(n)`. Worked in `promote_type(T, Float64)` and narrowed once,
so `Float32` transforms get correctly rounded coefficients.
"""
function _authalic_series(::Type{T}, n::Real, table) where {T<:AbstractFloat}
    W = promote_type(T, Float64)
    nw = W(n)
    return ntuple(Val(AUTHALIC_SERIES_ORDER)) do j
        row = table[j]
        value = W(row[AUTHALIC_SERIES_ORDER])
        for k in (AUTHALIC_SERIES_ORDER - 1):-1:1
            value = muladd(nw, value, W(row[k]))
        end
        return T(nw * value)
    end
end

"""
    _authalic_qp(e2) -> T

`q(90°)` of Snyder eq. (3-12), i.e. `q_p = 1 + (1 − e²)·atanh(e)/e`. Exactly
`2` on the sphere (`e² = 0`), where the general form would be `0/0`.
"""
@inline function _authalic_qp(e2::T) where {T<:AbstractFloat}
    iszero(e2) && return T(2)
    e = sqrt(e2)
    return one(T) + (one(T) - e2) * (atanh(e) / e)
end

function AuthalicTransform{T}(;
    semimajor_axis::Real=WGS84_SEMIMAJOR_AXIS,
    flattening::Union{Real,Nothing}=nothing,
    inverse_flattening::Union{Real,Nothing}=nothing,
    eccentricity_squared::Union{Real,Nothing}=nothing,
) where {T<:AbstractFloat}
    given = (flattening !== nothing) + (inverse_flattening !== nothing) +
            (eccentricity_squared !== nothing)
    given == 1 || throw(EllipsoidShapeError(given))

    W = promote_type(T, Float64)
    # Compute `n` from the supplied parameter without cancellation near a sphere.
    if eccentricity_squared === nothing
        f = flattening === nothing ? inv(W(inverse_flattening)) : W(flattening)
        # The closed form below supports oblate and spherical shapes only.
        (0 <= f < 1) || throw(DomainError(f,
            "flattening must lie in [0, 1) — an oblate ellipsoid of revolution"))
        e2 = f * (2 - f)
        n = f / (2 - f)
    else
        e2 = W(eccentricity_squared)
        (0 <= e2 < 1) || throw(DomainError(e2, "eccentricity_squared must lie in [0, 1)"))
        root = sqrt(one(W) - e2)
        n = e2 / ((one(W) + root) * (one(W) + root))
    end

    a = W(semimajor_axis)
    radius = a * sqrt(_authalic_qp(e2) / 2)
    return AuthalicTransform{T}(
        T(a), T(e2), T(radius),
        _authalic_series(T, n, AUTHALIC_FORWARD_TABLE),
        _authalic_series(T, n, AUTHALIC_INVERSE_TABLE),
    )
end

AuthalicTransform(; kwargs...) = AuthalicTransform{Float64}(; kwargs...)

# Reprecision uses `e²`; converting to the existing element type returns the
# original transform and preserves its stored coefficient rounding.
AuthalicTransform{T}(t::AuthalicTransform{T}) where {T<:AbstractFloat} = t

AuthalicTransform{T}(t::AuthalicTransform) where {T<:AbstractFloat} =
    AuthalicTransform{T}(; semimajor_axis=t.semimajor_axis,
        eccentricity_squared=t.eccentricity_squared)

function Base.show(io::IO, t::AuthalicTransform{T}) where {T}
    print(io, "AuthalicTransform{", T, "}(a = ", t.semimajor_axis,
        ", e² = ", t.eccentricity_squared, ", R_A = ", t.authalic_radius, ")")
end

"""
    WGS84_AUTHALIC

The [`AuthalicTransform`](@ref) for WGS84, with authalic radius
`6.371007180918474e6` m.
"""
const WGS84_AUTHALIC = AuthalicTransform{Float64}(;
    semimajor_axis=WGS84_SEMIMAJOR_AXIS,
    inverse_flattening=WGS84_INVERSE_FLATTENING,
)

# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

# Clenshaw sum of `Σ cⱼ sin(jθ)` from `2cos(θ)` and `sin(θ)`.
@inline function _clenshaw_sin(two_cos::T, sinθ::T, coefficients::NTuple{N,T}) where {N,T}
    b1 = zero(T)   # b_{j+1}
    b2 = zero(T)   # b_{j+2}
    @inbounds for j in N:-1:1
        b1, b2 = muladd(two_cos, b1, coefficients[j] - b2), b1
    end
    return sinθ * b1
end

"""
    geodetic_to_authalic(t::AuthalicTransform, latitude) -> T

Convert geodetic latitude `φ` to authalic latitude `ξ`, in radians:

    ξ = φ + Σ_{j=1}^{6} C_j sin(2jφ)          [Karney 2024, eq. (A19)]

For WGS84, `Float64` error is at most 0.5 ulp against a 300-bit closed-form
reference. The order-6 truncation remains below `Float64` roundoff for geodetic
datums; use [`authalic_q`](@ref) for much larger flattenings. `±π/2` and `0`
are exact fixed points. The result has the transform's element type.

See [`geodetic_to_authalicd`](@ref) for the degree form and
[`authalic_to_geodetic`](@ref) for the inverse.
"""
@inline function geodetic_to_authalic(t::AuthalicTransform{T}, latitude::Real) where {T}
    φ = T(latitude)
    s, c = sincos(2φ)
    return φ + _clenshaw_sin(2c, s, t.fwd)
end

"""
    authalic_to_geodetic(t::AuthalicTransform, latitude) -> T

Convert authalic latitude `ξ` to geodetic latitude `φ`, in radians:

    φ = ξ + Σ_{j=1}^{6} C'_j sin(2jξ)         [Karney 2024, eq. (A20)]

This evaluates Karney's order-6 reverted series. For WGS84, `Float64` error is
at most 0.5 ulp and forward/inverse round trips close within 1 ulp against a
300-bit reference. `±π/2` and `0` are exact fixed points.
"""
@inline function authalic_to_geodetic(t::AuthalicTransform{T}, latitude::Real) where {T}
    ξ = T(latitude)
    s, c = sincos(2ξ)
    return ξ + _clenshaw_sin(2c, s, t.inv)
end

"""
    geodetic_to_authalicd(t::AuthalicTransform, latitude) -> T

[`geodetic_to_authalic`](@ref) in degrees. It uses `sincospi(latitude/90)`
so `±90°` and `0°` are exact fixed points. Accuracy matches the radian form.
"""
@inline function geodetic_to_authalicd(t::AuthalicTransform{T}, latitude::Real) where {T}
    φ = T(latitude)
    s, c = sincospi(φ / 90)
    return φ + rad2deg(_clenshaw_sin(2c, s, t.fwd))
end

"""
    authalic_to_geodeticd(t::AuthalicTransform, latitude) -> T

[`authalic_to_geodetic`](@ref) in **degrees**; see
[`geodetic_to_authalicd`](@ref) for why the degree path is exact at `±90` and
`0` rather than merely accurate there.
"""
@inline function authalic_to_geodeticd(t::AuthalicTransform{T}, latitude::Real) where {T}
    ξ = T(latitude)
    s, c = sincospi(ξ / 90)
    return ξ + rad2deg(_clenshaw_sin(2c, s, t.inv))
end

# ---------------------------------------------------------------------------
# Closed form and the authalic radius
# ---------------------------------------------------------------------------

"""
    authalic_q(e2, sinlat) -> T

Snyder eq. (3-12), the authalic area function:

    q(φ) = (1 − e²) · [ sin φ / (1 − e² sin²φ) + atanh(e sin φ) / e ]

Its value is proportional to ellipsoidal area between the equator and `φ`.
Then `ξ = asin(q(φ)/q(90°))` and `R_A = a·√(q(90°)/2)`. This is the
closed-form reference and the fallback for flattenings beyond the order-6
series' range. Near the poles, prefer [`geodetic_to_authalic`](@ref) to avoid
the ill-conditioning of `asin`.

Returns `2 sin φ` exactly on the sphere (`e² = 0`), where the general form
would evaluate `0/0`.
"""
@inline function authalic_q(e2::T, sinlat::T) where {T<:AbstractFloat}
    iszero(e2) && return 2 * sinlat
    e = sqrt(e2)
    return (one(T) - e2) *
           (sinlat / (one(T) - e2 * sinlat * sinlat) + atanh(e * sinlat) / e)
end

authalic_q(e2::Real, sinlat::Real) = authalic_q(promote(float(e2), float(sinlat))...)

"""
    authalic_radius(t::AuthalicTransform) -> T
    authalic_radius(semimajor_axis, eccentricity_squared) -> T

Radius of the sphere with the same area as the ellipsoid, Snyder eq. (3-13):

    R_A = a · √(q_p / 2),    q_p = q(90°) = 1 + (1 − e²)·atanh(e)/e

For WGS84 this is `6.371007180918474e6` m. It equals `a` when `e² = 0`. The
one-argument form returns the precomputed field.
"""
@inline authalic_radius(t::AuthalicTransform) = t.authalic_radius

function authalic_radius(semimajor_axis::Real, eccentricity_squared::Real)
    a, e2 = promote(float(semimajor_axis), float(eccentricity_squared))
    (0 <= e2 < 1) || throw(DomainError(e2, "eccentricity_squared must lie in [0, 1)"))
    return a * sqrt(_authalic_qp(e2) / 2)
end
