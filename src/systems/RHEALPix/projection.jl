# Analytic HEALPix and rHEALPix projection on the authalic unit sphere.
# Formulas follow Calabretta & Roukema (2007) and Gibb, Raichev & Speth
# (2013/2016).  The rearrangement is only rotations and translations, so it
# preserves HEALPix's equal-area Jacobian.

const HALFPI = Float64(pi) / 2
const QUARTERPI = Float64(pi) / 4
const THREEQUARTERPI = 3 * Float64(pi) / 4
const TRANSITION_LATITUDE = asin(2 / 3)

@inline function _wrap_longitude(lon::Real)
    x = mod(Float64(lon) + Float64(pi), 2 * Float64(pi)) - Float64(pi)
    # Preserve the canonical half-open interval at the positive endpoint.
    return x == Float64(pi) ? -Float64(pi) : x
end

"""
    healpix_forward(longitude, latitude) -> (x, y)

Forward HEALPix projection of the unit sphere.  Angles and output coordinates
are radians/authalic radii.  Longitude is normalized to `[-π, π)`.
"""
function healpix_forward(longitude::Real, latitude::Real)
    lon = _wrap_longitude(longitude)
    lat = Float64(latitude)
    -HALFPI <= lat <= HALFPI || throw(DomainError(latitude,
        "latitude must lie in [-π/2, π/2]"))
    s = sin(lat)
    if abs(s) <= 2 / 3
        return lon, (3 * Float64(pi) / 8) * s
    end

    sigma = sqrt(3 * (1 - abs(s)))
    c = min(3, floor(Int, 2 * lon / Float64(pi) + 2))
    lonc = -THREEQUARTERPI + c * HALFPI
    x = lonc + (lon - lonc) * sigma
    y = copysign(QUARTERPI * (2 - sigma), lat)
    return x, y
end

"""
    healpix_inverse(x, y) -> (longitude, latitude)

Inverse unit-sphere HEALPix projection.  A pole has no intrinsic longitude and
is canonicalized to `-π`, matching the sealed reference vectors.
"""
function healpix_inverse(x::Real, y::Real)
    X, Y = Float64(x), Float64(y)
    abs(Y) <= HALFPI + 32eps(Float64) || throw(DomainError(y,
        "HEALPix y must lie in [-π/2, π/2]"))
    tol = 32eps(Float64)
    -Float64(pi)-tol <= X <= Float64(pi)+tol || throw(DomainError(x,
        "HEALPix x must lie in [-π, π]"))
    if abs(Y) > QUARTERPI
        tau = max(0.0, 2 - 4abs(Y) / Float64(pi))
        c = clamp(floor(Int, 2X / Float64(pi) + 2), 0, 3)
        lonc = -THREEQUARTERPI + c * HALFPI
        abs(X - lonc) <= tau * QUARTERPI + tol || throw(DomainError((x, y),
            "point lies outside the interrupted HEALPix image"))
    end
    if abs(Y) <= QUARTERPI
        return _wrap_longitude(X), asin(clamp(8 * Y / (3 * Float64(pi)), -1.0, 1.0))
    end

    tau = 2 - 4 * abs(Y) / Float64(pi)
    if abs(tau) <= 32eps(Float64)
        return -Float64(pi), copysign(HALFPI, Y)
    end
    c = clamp(floor(Int, 2 * X / Float64(pi) + 2), 0, 3)
    lonc = -THREEQUARTERPI + c * HALFPI
    lon = lonc + (X - lonc) / tau
    lat = copysign(asin(clamp(1 - tau * tau / 3, -1.0, 1.0)), Y)
    return _wrap_longitude(lon), lat
end

@inline function _rotate_quarters(x::Float64, y::Float64, turns::Integer)
    k = mod(Int(turns), 4)
    k == 0 && return (x, y)
    k == 1 && return (-y, x)
    k == 2 && return (-x, -y)
    return (y, -x)
end

@inline function _healpix_triangle(x::Float64)
    x < -HALFPI && return 0
    x < 0 && return 1
    x < HALFPI && return 2
    return 3
end

# The inverse diagonal ownership table is Appendix B of Gibb et al.  Keeping
# it in one function makes every seam decision explicit.  `tol` is used only
# to absorb roundoff after a rigid quarter-turn; it is not a grid-depth nudge.
function _inverse_triangle(x::Float64, y::Float64, square::Int, north::Bool)
    tol = 8eps(Float64)
    if north
        l1 = x - (-THREEQUARTERPI + (square - 1) * HALFPI)
        l2 = -x + (-THREEQUARTERPI + (square + 1) * HALFPI)
        y < l1 - tol && y >= l2 - tol && return mod(square + 1, 4)
        y >= l1 - tol && y > l2 + tol && return mod(square + 2, 4)
        y > l1 + tol && y <= l2 + tol && return mod(square + 3, 4)
    else
        l1 = x - (-THREEQUARTERPI + (square + 1) * HALFPI)
        l2 = -x + (-THREEQUARTERPI + (square - 1) * HALFPI)
        y <= l1 + tol && y > l2 + tol && return mod(square + 1, 4)
        y < l1 - tol && y <= l2 + tol && return mod(square + 2, 4)
        y >= l1 - tol && y < l2 - tol && return mod(square + 3, 4)
    end
    return square
end

function _combine_triangles(x::Float64, y::Float64, north_square::Int,
        south_square::Int)
    abs(y) <= QUARTERPI && return (x, y)
    c = _healpix_triangle(x)
    north = y > 0
    square = north ? north_square : south_square
    tipx = -THREEQUARTERPI + c * HALFPI
    tipy = copysign(HALFPI, y)
    ux = -THREEQUARTERPI + square * HALFPI
    uy = tipy
    turns = north ? c - square : -(c - square)
    rx, ry = _rotate_quarters(x - tipx, y - tipy, turns)
    return rx + ux, ry + uy
end

function _uncombine_triangles(x::Float64, y::Float64, north_square::Int,
        south_square::Int)
    abs(y) <= QUARTERPI && return (x, y)
    north = y > 0
    square = north ? north_square : south_square
    c = _inverse_triangle(x, y, square, north)
    ux = -THREEQUARTERPI + square * HALFPI
    uy = copysign(HALFPI, y)
    tipx = -THREEQUARTERPI + c * HALFPI
    tipy = uy
    turns = north ? -(c - square) : c - square
    rx, ry = _rotate_quarters(x - ux, y - uy, turns)
    return rx + tipx, ry + tipy
end

"""
    in_rhealpix_image(x, y, north_square=0, south_square=0) -> Bool

Whether `(x,y)` lies in the closed rHEALPix projection image.  The image is
the equatorial `4 × 1` strip plus one `1 × 1` polar square on either side.
"""
function in_rhealpix_image(x::Real, y::Real, north_square::Integer=0,
        south_square::Integer=0)
    X, Y = Float64(x), Float64(y)
    n, s = mod(Int(north_square), 4), mod(Int(south_square), 4)
    -Float64(pi) <= X <= Float64(pi) || return false
    abs(Y) <= QUARTERPI && return true
    if QUARTERPI <= Y <= THREEQUARTERPI
        lo = -Float64(pi) + n * HALFPI
        return lo <= X <= lo + HALFPI
    elseif -THREEQUARTERPI <= Y <= -QUARTERPI
        lo = -Float64(pi) + s * HALFPI
        return lo <= X <= lo + HALFPI
    end
    return false
end

"""
    rhealpix_forward(longitude, latitude; north_square=0, south_square=0,
                     longitude_origin=0) -> (x, y)

Forward rHEALPix projection on the unit authalic sphere.  All quantities are
in radians/authalic radii.
"""
function rhealpix_forward(longitude::Real, latitude::Real;
        north_square::Integer=0, south_square::Integer=0,
        longitude_origin::Real=0)
    n, s = mod(Int(north_square), 4), mod(Int(south_square), 4)
    x, y = healpix_forward(Float64(longitude) - Float64(longitude_origin), latitude)
    return _combine_triangles(x, y, n, s)
end

"""
    rhealpix_inverse(x, y; north_square=0, south_square=0,
                     longitude_origin=0) -> (longitude, latitude)

Inverse unit-authalic-sphere rHEALPix projection.  Throws `DomainError` outside
the projection image.
"""
function rhealpix_inverse(x::Real, y::Real;
        north_square::Integer=0, south_square::Integer=0,
        longitude_origin::Real=0)
    n, s = mod(Int(north_square), 4), mod(Int(south_square), 4)
    in_rhealpix_image(x, y, n, s) || throw(DomainError((x, y),
        "point lies outside the rHEALPix projection image"))
    X, Y = _uncombine_triangles(Float64(x), Float64(y), n, s)
    lon, lat = healpix_inverse(X, Y)
    return _wrap_longitude(lon + Float64(longitude_origin)), lat
end

"""
    auspix_forward(longitude_degrees, latitude_degrees) -> (x_metres, y_metres)

AusPIX's fixed WGS84, Greenwich, `(north_square,south_square)=(0,0)` profile.
Input angles are geodetic degrees; output coordinates are metres on the WGS84
authalic sphere.
"""
function auspix_forward(longitude_degrees::Real, latitude_degrees::Real)
    beta = DGG.Helpers.geodetic_to_authalicd(
        DGG.Helpers.WGS84_AUTHALIC, latitude_degrees)
    x, y = rhealpix_forward(deg2rad(longitude_degrees), deg2rad(beta))
    radius = DGG.Helpers.authalic_radius(DGG.Helpers.WGS84_AUTHALIC)
    return radius * x, radius * y
end

"""
    auspix_inverse(x_metres, y_metres) -> (longitude_degrees, latitude_degrees)

Inverse of [`auspix_forward`](@ref), returning WGS84 geodetic degrees.
"""
function auspix_inverse(x_metres::Real, y_metres::Real)
    radius = DGG.Helpers.authalic_radius(DGG.Helpers.WGS84_AUTHALIC)
    lon, beta = rhealpix_inverse(Float64(x_metres) / radius,
        Float64(y_metres) / radius)
    lat = DGG.Helpers.authalic_to_geodeticd(
        DGG.Helpers.WGS84_AUTHALIC, rad2deg(beta))
    return rad2deg(lon), lat
end
