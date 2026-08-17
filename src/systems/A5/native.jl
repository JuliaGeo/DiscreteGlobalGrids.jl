# Pure-Julia A5 arithmetic compatible with felixpalmer/a5. This module's
# `A5Cell` is a decoded record; `DiscreteGlobalGrids.A5.A5Cell` is the public
# encoded id.

"""
    DiscreteGlobalGrids.A5.A5Native

Low-level A5 projection, id, Hilbert-lattice, adjacency, and boundary
operations. Inputs and outputs use raw `UInt64` ids, decoded
`A5Native.A5Cell` records, and degree coordinates. Public grid contracts are
implemented by the enclosing [`A5`](@ref) module.
"""
module A5Native

using ...Helpers: SmallList, empty_small_list, small_push

export FIRST_HILBERT_RESOLUTION,
    MAX_GRID_RESOLUTION,
    MAX_RESOLUTION,
    WORLD_CELL,
    cell_area,
    cell_boundary,
    cell_boundary_cartesian,
    cell_to_children,
    cell_to_lonlat,
    cell_to_parent,
    deserialize,
    get_resolution,
    lonlat_to_cell,
    num_cells,
    res0_cells,
    serialize

const FIRST_HILBERT_RESOLUTION = 2
const MAX_RESOLUTION = 30
# Resolution 30 is available only for the 42 quintants that fit the 64-bit
# encoding. A uniform full-world grid therefore stops at resolution 29.
const MAX_GRID_RESOLUTION = 29
const HILBERT_START_BIT = 58
const WORLD_CELL = UInt64(0)
const LONGITUDE_OFFSET = 93.0

const PHI = (1 + sqrt(5.0)) / 2
const TWO_PI_OVER_5 = 2pi / 5
const PI_OVER_5 = pi / 5
const PI_OVER_10 = pi / 10
const DIHEDRAL_ANGLE = 2 * atan(PHI)
const INTERHEDRAL_ANGLE = pi - DIHEDRAL_ANGLE
const DISTANCE_TO_EDGE = (sqrt(5.0) - 1) / 2
const DISTANCE_TO_VERTEX = 3 - sqrt(5.0)
const AUTHALIC_RADIUS_EARTH = 6_371_007.2
const AUTHALIC_AREA_EARTH = 4pi * AUTHALIC_RADIUS_EARTH^2

const YES = -1
const NO = 1

const CLOCKWISE_FAN = (:vu, :uw, :vw, :vw, :vw)
const CLOCKWISE_STEP = (:wu, :uw, :vw, :vu, :uw)
const COUNTER_STEP = (:wu, :uv, :wv, :wu, :uw)
const COUNTER_JUMP = (:vu, :uv, :wv, :wu, :uw)
const QUINTANT_ORIENTATIONS = (
    CLOCKWISE_FAN,
    COUNTER_JUMP,
    COUNTER_STEP,
    CLOCKWISE_STEP,
    COUNTER_STEP,
    COUNTER_JUMP,
    COUNTER_STEP,
    CLOCKWISE_STEP,
    CLOCKWISE_STEP,
    CLOCKWISE_STEP,
    COUNTER_JUMP,
    COUNTER_JUMP,
)
const QUINTANT_FIRST = (4, 2, 3, 2, 0, 4, 3, 2, 2, 0, 3, 0)
const ORIGIN_ORDER = (0, 1, 2, 4, 3, 5, 7, 8, 6, 11, 10, 9)

struct Origin
    id::Int
    axis::NTuple{2,Float64}
    axis_cartesian::NTuple{3,Float64}
    quat::NTuple{4,Float64}
    inverse_quat::NTuple{4,Float64}
    angle::Float64
    orientation::NTuple{5,Symbol}
    first_quintant::Int
end

struct A5Cell
    origin::Origin
    segment::Int
    S::UInt64
    resolution::Int
end

struct Anchor
    q::Int
    offset::NTuple{2,Float64}
    flips::NTuple{2,Int}
end

_clamp1(x) = clamp(x, -1.0, 1.0)
_add2(a, b) = (a[1] + b[1], a[2] + b[2])
_sub2(a, b) = (a[1] - b[1], a[2] - b[2])
_scale2(a, s) = (a[1] * s, a[2] * s)
_neg2(a) = (-a[1], -a[2])
_lerp2(a, b, t) = (a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t)
_rotate2(a, theta) = (cos(theta) * a[1] - sin(theta) * a[2], sin(theta) * a[1] + cos(theta) * a[2])
_norm2(a) = hypot(a[1], a[2])

_add3(a, b) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
_sub3(a, b) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
_scale3(a, s) = (a[1] * s, a[2] * s, a[3] * s)
_dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
_cross3(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
_norm3(a) = sqrt(_dot3(a, a))
_lerp3(a, b, t) = (a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t, a[3] + (b[3] - a[3]) * t)
function _normalize3(a)
    n = _norm3(a)
    n == 0 && return (0.0, 0.0, 0.0)
    return _scale3(a, inv(n))
end

_quat_conjugate(q) = (-q[1], -q[2], -q[3], q[4])

function _quat_normalize(q)
    n = sqrt(q[1]^2 + q[2]^2 + q[3]^2 + q[4]^2)
    return (q[1] / n, q[2] / n, q[3] / n, q[4] / n)
end

function _quat_from_axis_angle(axis, angle)
    half = angle / 2
    s = sin(half)
    return (axis[1] * s, axis[2] * s, axis[3] * s, cos(half))
end

function _quat_rotation_to(a, b)
    dotab = _dot3(a, b)
    if dotab < -0.999999
        axis = _cross3((1.0, 0.0, 0.0), a)
        if _norm3(axis) < 0.000001
            axis = _cross3((0.0, 1.0, 0.0), a)
        end
        return _quat_from_axis_angle(_normalize3(axis), pi)
    elseif dotab > 0.999999
        return (0.0, 0.0, 0.0, 1.0)
    else
        crossab = _cross3(a, b)
        return _quat_normalize((crossab[1], crossab[2], crossab[3], 1 + dotab))
    end
end

function _quat_rotate(q, v)
    x, y, z, w = q
    vx, vy, vz = v
    uv = _cross3((x, y, z), (vx, vy, vz))
    uuv = _cross3((x, y, z), uv)
    return _add3(v, _add3(_scale3(uv, 2w), _scale3(uuv, 2.0)))
end

function _to_cartesian(spherical)
    theta, phi = spherical
    sinphi = sin(phi)
    return (sinphi * cos(theta), sinphi * sin(theta), cos(phi))
end

function _to_spherical(xyz)
    x, y, z = xyz
    r = _norm3(xyz)
    return (atan(y, x), acos(_clamp1(z / r)))
end

_to_polar(xy) = (_norm2(xy), atan(xy[2], xy[1]))
_to_face(polar) = (polar[1] * cos(polar[2]), polar[1] * sin(polar[2]))

const GEODETIC_TO_AUTHALIC = (
    -2.2392098386786394e-3,
    2.1308606513250217e-6,
    -2.5592576864212742e-9,
    3.3701965267802837e-12,
    -4.6675453126112487e-15,
    6.6749287038481596e-18,
)
const AUTHALIC_TO_GEODETIC = (
    2.2392089963541657e-3,
    2.8831978048607556e-6,
    5.0862207399726603e-9,
    1.02018123778161e-11,
    2.1912872306767718e-14,
    4.9284235482523806e-17,
)

function _apply_authalic(phi, C)
    sinphi = sin(phi)
    cosphi = cos(phi)
    X = 2 * (cosphi - sinphi) * (cosphi + sinphi)
    u0 = X * C[6] + C[5]
    u1 = X * u0 + C[4]
    u0 = X * u1 - u0 + C[3]
    u1 = X * u0 - u1 + C[2]
    u0 = X * u1 - u0 + C[1]
    return phi + 2 * sinphi * cosphi * u0
end

_authalic_forward(phi) = _apply_authalic(phi, GEODETIC_TO_AUTHALIC)
_authalic_inverse(phi) = _apply_authalic(phi, AUTHALIC_TO_GEODETIC)

function _from_lonlat(lonlat)
    longitude, latitude = lonlat
    theta = deg2rad(longitude + LONGITUDE_OFFSET)
    authalic_lat = _authalic_forward(deg2rad(latitude))
    return (theta, pi / 2 - authalic_lat)
end

function _normalize_longitude(lon)
    return mod(lon + 180, 360) - 180
end

function _to_lonlat(spherical)
    theta, phi = spherical
    longitude = _normalize_longitude(rad2deg(theta) - LONGITUDE_OFFSET)
    authalic_lat = pi / 2 - phi
    latitude = rad2deg(_authalic_inverse(authalic_lat))
    return (longitude, latitude)
end

function _make_quaternions()
    sqrt5 = sqrt(5.0)
    invsqrt5 = sqrt(0.2)
    sin_alpha = sqrt((1 - invsqrt5) / 2)
    cos_alpha = sqrt((1 + invsqrt5) / 2)
    A = 0.5
    B = sqrt((2.5 - sqrt5) / 10)
    C = sqrt((2.5 + sqrt5) / 10)
    D = sqrt((1 + invsqrt5) / 8)
    E = sqrt((1 - invsqrt5) / 8)
    F = sqrt((3 - sqrt5) / 8)
    G = sqrt((3 + sqrt5) / 8)
    centers = [
        (0.0, 0.0),
        (sin_alpha, 0.0),
        (B, A),
        (-D, F),
        (-D, -F),
        (B, -A),
        (-cos_alpha, 0.0),
        (-E, -G),
        (C, -A),
        (C, A),
        (-E, G),
        (0.0, 0.0),
    ]
    return ntuple(12) do i
        id = i - 1
        if id == 0
            (0.0, 0.0, 0.0, 1.0)
        elseif id == 11
            (0.0, -1.0, 0.0, 0.0)
        else
            x, y = centers[i]
            ax, ay = -y, x
            (ax, ay, 0.0, id < 6 ? cos_alpha : sin_alpha)
        end
    end
end

const QUATERNIONS_BY_SOURCE_ID = _make_quaternions()

function _haversine(point, axis)
    theta, phi = point
    theta2, phi2 = axis
    dtheta = theta2 - theta
    dphi = phi2 - phi
    A1 = sin(dphi / 2)
    A2 = sin(dtheta / 2)
    return A1 * A1 + A2 * A2 * sin(phi) * sin(phi2)
end

function _make_origins()
    tmp = Origin[]
    function add_origin(axis, angle, quat)
        source_id = length(tmp)
        push!(tmp, Origin(
            source_id,
            axis,
            _to_cartesian(axis),
            quat,
            _quat_conjugate(quat),
            angle,
            QUINTANT_ORIENTATIONS[source_id + 1],
            QUINTANT_FIRST[source_id + 1],
        ))
    end
    add_origin((0.0, 0.0), 0.0, QUATERNIONS_BY_SOURCE_ID[1])
    for i in 0:4
        alpha = i * TWO_PI_OVER_5
        alpha2 = alpha + PI_OVER_5
        add_origin((alpha, INTERHEDRAL_ANGLE), PI_OVER_5, QUATERNIONS_BY_SOURCE_ID[i + 2])
        add_origin((alpha2, pi - INTERHEDRAL_ANGLE), PI_OVER_5,
            QUATERNIONS_BY_SOURCE_ID[((i + 3) % 5) + 7])
    end
    add_origin((0.0, pi), 0.0, QUATERNIONS_BY_SOURCE_ID[12])
    order_rank = Dict(v => i for (i, v) in enumerate(ORIGIN_ORDER))
    sort!(tmp; by=o -> order_rank[o.id])
    return Tuple(Origin(
        i - 1,
        o.axis,
        o.axis_cartesian,
        o.quat,
        o.inverse_quat,
        o.angle,
        o.orientation,
        o.first_quintant,
    ) for (i, o) in enumerate(tmp))
end

const ORIGINS = _make_origins()

function _find_nearest_origin(spherical)
    nearest = ORIGINS[1]
    min_distance = Inf
    for origin in ORIGINS
        distance = _haversine(spherical, origin.axis)
        if distance < min_distance
            min_distance = distance
            nearest = origin
        end
    end
    return nearest
end

function _find_nearest_origin_cartesian(c)
    nearest = ORIGINS[1]
    min_distance = Inf
    for origin in ORIGINS
        ax = origin.axis_cartesian
        distance = 1 - _dot3(c, ax)
        if distance < min_distance
            min_distance = distance
            nearest = origin
        end
    end
    return nearest
end

function _quintant_to_segment(quintant::Integer, origin::Origin)
    layout = origin.orientation
    step = (layout == CLOCKWISE_FAN || layout == CLOCKWISE_STEP) ? -1 : 1
    delta = mod(Int(quintant) - origin.first_quintant, 5)
    face_relative = mod(step * delta, 5)
    orientation = layout[face_relative + 1]
    segment = mod(origin.first_quintant + face_relative, 5)
    return segment, orientation
end

function _segment_to_quintant(segment::Integer, origin::Origin)
    layout = origin.orientation
    step = (layout == CLOCKWISE_FAN || layout == CLOCKWISE_STEP) ? -1 : 1
    face_relative = mod(Int(segment) - origin.first_quintant, 5)
    orientation = layout[face_relative + 1]
    quintant = mod(origin.first_quintant + step * face_relative, 5)
    return quintant, orientation
end

function get_resolution(index::UInt64)
    index == 0 && return -1
    ((index & 0x01) != 0 || (index & 0x07) == 0x04 || (index & 0x1f) == 0x10) && return MAX_RESOLUTION
    resolution = MAX_RESOLUTION - 1
    shifted = index >> 1
    shifted == 0 && return -1
    while resolution > -1 && (shifted & 0x01) == 0
        resolution -= 1
        shifted >>= resolution < FIRST_HILBERT_RESOLUTION ? 1 : 2
    end
    return resolution
end

function deserialize(index::UInt64)
    resolution = get_resolution(index)
    if resolution == -1
        return A5Cell(ORIGINS[1], 0, 0x00, resolution)
    end

    quintant_shift = HILBERT_START_BIT
    quintant_offset = 0
    if resolution == MAX_RESOLUTION
        marker_bits = (index & 0x01) != 0 ? 1 :
                      (index & 0x04) != 0 ? 3 : 5
        quintant_shift = HILBERT_START_BIT + marker_bits
        quintant_offset = marker_bits == 1 ? 0 : marker_bits == 3 ? 32 : 40
    end

    top_bits = Int(index >> quintant_shift) + quintant_offset
    if resolution == 0
        origin = ORIGINS[top_bits + 1]
        segment = 0
    else
        origin = ORIGINS[fld(top_bits, 5) + 1]
        segment = mod(top_bits + origin.first_quintant, 5)
    end

    if resolution < FIRST_HILBERT_RESOLUTION
        return A5Cell(origin, segment, 0x00, resolution)
    end

    hilbert_levels = resolution - FIRST_HILBERT_RESOLUTION + 1
    hilbert_bits = 2 * hilbert_levels
    removal_mask = (UInt64(1) << quintant_shift) - UInt64(1)
    S = (index & removal_mask) >> (quintant_shift - hilbert_bits)
    return A5Cell(origin, segment, S, resolution)
end

function serialize(cell::A5Cell)
    origin, segment, S, resolution = cell.origin, cell.segment, cell.S, cell.resolution
    resolution > MAX_RESOLUTION && throw(ArgumentError("Resolution ($resolution) is too large"))
    resolution == -1 && return WORLD_CELL

    quintant_shift = HILBERT_START_BIT
    R = resolution < FIRST_HILBERT_RESOLUTION ? resolution + 1 :
        2 * (1 + resolution - FIRST_HILBERT_RESOLUTION) + 1
    segment_n = mod(segment - origin.first_quintant, 5)

    index = UInt64(0)
    if resolution == 0
        index = UInt64(origin.id) << quintant_shift
    else
        quintant = 5 * origin.id + segment_n
        if resolution == MAX_RESOLUTION
            quintant_value = 0
            if quintant <= 31
                quintant_shift = HILBERT_START_BIT + 1
                quintant_value = quintant
            elseif quintant <= 39
                quintant_shift = HILBERT_START_BIT + 3
                quintant_value = quintant - 32
            elseif quintant <= 41
                quintant_shift = HILBERT_START_BIT + 5
                quintant_value = quintant - 40
            else
                return serialize(A5Cell(origin, segment, S >> 2, MAX_RESOLUTION - 1))
            end
            index = UInt64(quintant_value) << quintant_shift
        else
            index = UInt64(quintant) << quintant_shift
        end
    end

    if resolution >= FIRST_HILBERT_RESOLUTION
        hilbert_levels = resolution - FIRST_HILBERT_RESOLUTION + 1
        hilbert_bits = 2 * hilbert_levels
        S >= (UInt64(1) << hilbert_bits) &&
            throw(ArgumentError("S ($S) is too large for resolution level $resolution"))
        index += S << (quintant_shift - hilbert_bits)
    end

    index |= UInt64(1) << (quintant_shift - R)
    return index
end

function _write_children!(
        children::Vector{UInt64},
        child_index::Int,
        origin::Origin,
        first_segment::Int,
        last_segment::Int,
        shifted_S::UInt64,
        children_count::UInt64,
        resolution::Int,
    )
    for segment in first_segment:last_segment
        for i in UInt64(0):(children_count - UInt64(1))
            children[child_index] = serialize(A5Cell(origin, segment, shifted_S + i, resolution))
            child_index += 1
        end
    end
    return child_index
end

function cell_to_children(index::UInt64, child_resolution::Union{Nothing,Integer}=nothing)
    cell = deserialize(index)
    current_resolution = cell.resolution
    new_resolution = child_resolution === nothing ? current_resolution + 1 : Int(child_resolution)
    new_resolution < current_resolution &&
        throw(ArgumentError("Target resolution ($new_resolution) must be >= current resolution ($current_resolution)"))
    new_resolution > MAX_RESOLUTION &&
        throw(ArgumentError("Target resolution ($new_resolution) exceeds maximum resolution ($MAX_RESOLUTION)"))
    if new_resolution == current_resolution
        return small_push(empty_small_list(Val(1), UInt64(0)), index)
    end
    index == WORLD_CELL && new_resolution == 0 && return RES0_CELLS

    if new_resolution == current_resolution + 1
        if current_resolution == 0
            children = empty_small_list(Val(5), UInt64(0))
            for segment in 0:4
                child = serialize(A5Cell(cell.origin, segment, cell.S, new_resolution))
                children = small_push(children, child)
            end
            return children
        elseif current_resolution >= 1
            children = empty_small_list(Val(4), UInt64(0))
            shifted_S = cell.S << 2
            for i in UInt64(0):UInt64(3)
                child = serialize(A5Cell(cell.origin, cell.segment, shifted_S + i, new_resolution))
                get_resolution(child) == new_resolution ||
                    throw(ArgumentError(
                        "A5 resolution $new_resolution is not representable for this cell",
                    ))
                children = small_push(children, child)
            end
            return children
        end
    end

    if new_resolution == MAX_RESOLUTION
        current_resolution >= 1 ||
            throw(ArgumentError(
                "A5 resolution $MAX_RESOLUTION is not a uniform full-grid resolution",
            ))
        segment_n = mod(cell.segment - cell.origin.first_quintant, 5)
        quintant = 5 * cell.origin.id + segment_n
        quintant <= 41 ||
            throw(ArgumentError(
                "A5 resolution $MAX_RESOLUTION is not representable for this cell",
            ))
    end

    all_origins = current_resolution == -1
    all_segments = (all_origins && new_resolution > 0) || current_resolution == 0
    origin_count = all_origins ? length(ORIGINS) : 1
    segment_count = all_segments ? 5 : 1
    resolution_diff = max(0, new_resolution - max(current_resolution, FIRST_HILBERT_RESOLUTION - 1))
    children_count = UInt64(4)^UInt64(resolution_diff)
    shifted_S = cell.S << (2 * resolution_diff)
    total = children_count * UInt64(origin_count) * UInt64(segment_count)
    total <= UInt64(typemax(Int)) ||
        throw(OverflowError("A5 child count $total does not fit in Int"))
    children = Vector{UInt64}(undef, Int(total))
    child_index = 1
    first_segment = all_segments ? 0 : cell.segment
    last_segment = all_segments ? 4 : cell.segment
    if all_origins
        for origin in ORIGINS
            child_index = _write_children!(
                children, child_index, origin, first_segment, last_segment,
                shifted_S, children_count, new_resolution,
            )
        end
    else
        _write_children!(
            children, child_index, cell.origin, first_segment, last_segment,
            shifted_S, children_count, new_resolution,
        )
    end
    return children
end

function _is_max_resolution(index::UInt64)
    return (index & 0x01) != 0 || (index & 0x07) == 0x04 || (index & 0x1f) == 0x10
end

function _normalize_res30(index::UInt64)
    if (index & 0x01) != 0
        qshift, qoffset, marker_bits = 59, UInt64(0), 1
    elseif (index & 0x04) != 0
        qshift, qoffset, marker_bits = 61, UInt64(32), 3
    else
        qshift, qoffset, marker_bits = 63, UInt64(40), 5
    end
    quintant = (index >> qshift) + qoffset
    s58 = (index >> marker_bits) & ((UInt64(1) << 58) - UInt64(1))
    return (quintant << 58) | ((s58 >> 2) << 2) | (UInt64(1) << 1)
end

function cell_to_parent(index::UInt64, parent_resolution::Union{Nothing,Integer}=nothing)
    current_resolution = get_resolution(index)
    pres = parent_resolution === nothing ? current_resolution - 1 : Int(parent_resolution)
    pres <= current_resolution ||
        throw(ArgumentError(
            "Target resolution ($pres) must be <= current resolution ($current_resolution)",
        ))
    pres == -1 && return WORLD_CELL
    (pres < -1 || pres > MAX_RESOLUTION) && throw(ArgumentError("Target resolution ($pres) is out of range"))
    index == WORLD_CELL &&
        throw(ArgumentError("Target resolution ($pres) must be <= current resolution (-1)"))

    c = index
    if _is_max_resolution(index)
        pres == MAX_RESOLUTION && return index
        c = _normalize_res30(index)
        pres == MAX_RESOLUTION - 1 && return c
    end

    if pres >= FIRST_HILBERT_RESOLUTION
        keep_shift = 60 - 2 * pres
        return ((c >> keep_shift) << keep_shift) | (UInt64(1) << (59 - 2 * pres))
    elseif pres == 1
        return ((c >> 58) << 58) | (UInt64(1) << 56)
    else
        ((c & ((UInt64(1) << 57) - UInt64(1))) == 0) && return c
        return ((div(c >> 58, UInt64(5))) << 58) | (UInt64(1) << 57)
    end
end

const RES0_CELLS = SmallList{12,UInt64}(
    12,
    ntuple(i -> serialize(A5Cell(ORIGINS[i], 0, 0x00, 0)), 12),
)
res0_cells() = RES0_CELLS

function num_cells(resolution::Integer)
    res = Int(resolution)
    res == -1 && return UInt64(1)
    0 <= res <= MAX_GRID_RESOLUTION ||
        throw(ArgumentError("A5 full-grid resolution must be in 0:$MAX_GRID_RESOLUTION"))
    res == 0 && return UInt64(12)
    return UInt64(60) * UInt64(4)^UInt64(res - 1)
end

function cell_area(resolution::Integer)
    res = Int(resolution)
    -1 <= res <= MAX_RESOLUTION ||
        throw(ArgumentError("A5 resolution must be in -1:$MAX_RESOLUTION"))
    count = res == -1 ? UInt64(1) :
        res == 0 ? UInt64(12) :
        UInt64(60) * UInt64(4)^UInt64(res - 1)
    return AUTHALIC_AREA_EARTH / Float64(count)
end

function _shape_area(vertices)
    signed = 0.0
    n = length(vertices)
    for i in 1:n
        j = i == n ? 1 : i + 1
        signed += (vertices[j][1] - vertices[i][1]) * (vertices[j][2] + vertices[i][2])
    end
    return signed
end

_wind(vertices) = _shape_area(vertices) >= 0 ? vertices : reverse(vertices)

function _reflect_y(vertices)
    return reverse(map(v -> (v[1], -v[2]), vertices))
end

_translate(vertices, d) = map(v -> _add2(v, d), vertices)
_scale(vertices, s) = map(v -> _scale2(v, s), vertices)
_rotate(vertices, theta) = map(v -> _rotate2(v, theta), vertices)
_center(vertices) = (sum(v[1] for v in vertices) / length(vertices), sum(v[2] for v in vertices) / length(vertices))

function _contains_point(vertices, point)
    _shape_area(vertices) < 0 && error("Pentagon is not counter-clockwise")
    dmax = 1.0
    n = length(vertices)
    for i in 1:n
        v1 = vertices[i]
        v2 = vertices[i == n ? 1 : i + 1]
        dx = v1[1] - v2[1]
        dy = v1[2] - v2[2]
        px = point[1] - v1[1]
        py = point[2] - v1[2]
        cross_product = dx * py - dy * px
        if cross_product < 0
            plen = hypot(px, py)
            dmax = min(dmax, cross_product / plen)
        end
    end
    return dmax
end

function _split_edges(vertices, segments::Integer)
    segments <= 1 && return vertices
    new_vertices = NTuple{2,Float64}[]
    n = length(vertices)
    for i in 1:n
        v1 = vertices[i]
        v2 = vertices[i == n ? 1 : i + 1]
        push!(new_vertices, v1)
        for j in 1:(segments - 1)
            push!(new_vertices, _lerp2(v1, v2, j / segments))
        end
    end
    return _wind(new_vertices)
end

function _make_pentagon_basis()
    a = (0.0, 0.0)
    b = (0.0, 1.0)
    c = (0.7885966681787006, 1.6149108024237764)
    d = (1.6171013659387945, 1.054928690397459)
    e = (cos(PI_OVER_10), sin(PI_OVER_10))
    edge_midpoint_d = 2 * _norm2(c) * cos(PI_OVER_5)
    basis_rotation = PI_OVER_5 - atan(c[2], c[1])
    scale = (2 * DISTANCE_TO_EDGE) / edge_midpoint_d
    raw = (a, b, c, d, e)
    vertices = map(v -> _rotate2(_scale2(v, scale), basis_rotation), raw)
    pentagon = _wind(vertices)
    bisector_angle = atan(vertices[3][2], vertices[3][1]) - PI_OVER_5
    u = (0.0, 0.0)
    L = DISTANCE_TO_EDGE / cos(PI_OVER_5)
    Vangle = bisector_angle + PI_OVER_5
    v = (L * cos(Vangle), L * sin(Vangle))
    Wangle = bisector_angle - PI_OVER_5
    w = (L * cos(Wangle), L * sin(Wangle))
    triangle = _wind((u, v, w))
    return pentagon, triangle, u, v, w
end

const PENTAGON, TRIANGLE, U_VEC, V_VEC, W_VEC = _make_pentagon_basis()
const BASIS_DET = V_VEC[1] * W_VEC[2] - W_VEC[1] * V_VEC[2]
const BASIS_INVERSE = (W_VEC[2] / BASIS_DET, -V_VEC[2] / BASIS_DET, -W_VEC[1] / BASIS_DET, V_VEC[1] / BASIS_DET)

_ij_to_face(ij) = _add2(_scale2(V_VEC, ij[1]), _scale2(W_VEC, ij[2]))
_face_to_ij(face) = (
    BASIS_INVERSE[1] * face[1] + BASIS_INVERSE[3] * face[2],
    BASIS_INVERSE[2] * face[1] + BASIS_INVERSE[4] * face[2],
)
_kj_to_ij(kj) = (kj[1] - kj[2], kj[2])

function _get_pentagon_vertices(resolution::Integer, quintant::Integer, anchor::Anchor)
    vertices = PENTAGON
    if anchor.flips == (NO, YES)
        vertices = map(_neg2, vertices)
    end
    q = anchor.q
    F = anchor.flips[1] + anchor.flips[2]
    if (((F == -2 || F == 2) && q > 1) || (F == 0 && (q == 0 || q == 3)))
        vertices = _reflect_y(vertices)
    end
    if anchor.flips == (YES, YES)
        vertices = map(_neg2, vertices)
    elseif anchor.flips[1] == YES
        vertices = _translate(vertices, _neg2(W_VEC))
    elseif anchor.flips[2] == YES
        vertices = _translate(vertices, W_VEC)
    end
    vertices = _translate(vertices, _ij_to_face(anchor.offset))
    vertices = _scale(vertices, inv(2.0^resolution))
    return _rotate(vertices, TWO_PI_OVER_5 * quintant)
end

function _get_pentagon_flavor(anchor::Anchor)
    f = 0
    anchor.flips[2] == YES && (f += 2)
    q = anchor.q
    F = anchor.flips[1] + anchor.flips[2]
    if (((F == -2 || F == 2) && q > 1) || (F == 0 && (q == 0 || q == 3)))
        f += 1
    end
    (F == -2 || F == 2) && (f += 4)
    return f
end

function _get_quintant_vertices(quintant::Integer)
    return _rotate(TRIANGLE, TWO_PI_OVER_5 * quintant)
end

function _get_face_vertices()
    vertices = ntuple(i -> _rotate2(V_VEC, TWO_PI_OVER_5 * (i - 1)), 5)
    return _wind(reverse(vertices))
end

_js_round_int(x) = floor(Int, x + 0.5)
_get_quintant_polar(polar) = mod(_js_round_int(polar[2] / TWO_PI_OVER_5), 5)

const PATTERN = (0, 1, 3, 4, 5, 6, 7, 2)
const PATTERN_FLIPPED = (0, 1, 2, 7, 3, 4, 5, 6)
_reverse_pattern(pattern::NTuple{N,Int}) where {N} =
    ntuple(i -> findfirst(==(i - 1), pattern) - 1, Val(N))
const PATTERN_REVERSED = _reverse_pattern(PATTERN)
const PATTERN_FLIPPED_REVERSED = _reverse_pattern(PATTERN_FLIPPED)

_quaternary_digit(digits::UInt64, idx::Int) =
    Int((digits >> (2 * (idx - 1))) & UInt64(3))

function _set_quaternary_digit(digits::UInt64, idx::Int, digit::Integer)
    shift = 2 * (idx - 1)
    mask = UInt64(3) << shift
    return (digits & ~mask) | (UInt64(digit) << shift)
end

function _shift_digits(digits::UInt64, idx::Int, flips, invert_j::Bool, pattern)
    idx <= 1 && return digits
    parent_k = _quaternary_digit(digits, idx)
    child_k = _quaternary_digit(digits, idx - 1)
    F = flips[1] + flips[2]
    if invert_j != (F == 0)
        needs_shift = parent_k == 1 || parent_k == 2
        first = parent_k == 1
    else
        needs_shift = parent_k < 2
        first = parent_k == 0
    end
    needs_shift || return digits
    src = first ? child_k : child_k + 4
    dst = pattern[src + 1]
    digits = _set_quaternary_digit(digits, idx - 1, dst % 4)
    return _set_quaternary_digit(
        digits,
        idx,
        mod(parent_k + 4 + fld(dst, 4) - fld(src, 4), 4),
    )
end

function _quaternary_to_kj(n::Integer, flips)
    flipx, flipy = flips
    if flipx == NO && flipy == NO
        p, q = (1.0, 0.0), (0.0, 1.0)
    elseif flipx == YES && flipy == NO
        p, q = (0.0, -1.0), (-1.0, 0.0)
    elseif flipx == NO && flipy == YES
        p, q = (0.0, 1.0), (1.0, 0.0)
    else
        p, q = (-1.0, 0.0), (0.0, -1.0)
    end
    n == 0 && return (0.0, 0.0)
    n == 1 && return p
    n == 2 && return _add2(q, p)
    n == 3 && return _add2(q, _scale2(p, 2.0))
    throw(ArgumentError("Invalid quaternary value: $n"))
end

_quaternary_to_flips(n::Integer) =
    n == 0 ? (NO, NO) : n == 1 ? (NO, YES) : n == 2 ? (NO, NO) : (YES, NO)

function _ij_to_quaternary(ij, flips)
    i, j = ij
    a = flips[1] == YES ? -(i + j) : i + j
    b = flips[2] == YES ? -i : i
    c = flips[1] == YES ? -j : j
    if flips[1] + flips[2] == 0
        c < 1 ? 0 : b > 1 ? 3 : a > 1 ? 2 : 1
    else
        a < 1 ? 0 : b > 1 ? 3 : c > 1 ? 2 : 1
    end
end

function _s_to_anchor(s::UInt64, resolution::Integer, orientation::Symbol, do_shift_digits::Bool=true)
    input = s
    reverse = orientation in (:vu, :wu, :vw)
    invert_j = orientation in (:wv, :vw)
    flip_ij = orientation in (:wu, :uw)
    if reverse
        input = (UInt64(1) << (2 * resolution)) - input - UInt64(1)
    end
    anchor = __s_to_anchor(input, Int(resolution), invert_j, flip_ij, do_shift_digits)
    if flip_ij
        offset = (anchor.offset[2], anchor.offset[1])
        if anchor.flips[1] == YES
            offset = _add2(offset, (-1.0, 1.0))
        end
        if anchor.flips[2] == YES
            offset = _sub2(offset, (-1.0, 1.0))
        end
        anchor = Anchor(anchor.q, offset, anchor.flips)
    end
    if invert_j
        i, j = anchor.offset
        flips = (-anchor.flips[1], anchor.flips[2])
        anchor = Anchor(anchor.q, (i, (1 << resolution) - (i + j)), flips)
    end
    return anchor
end

function __s_to_anchor(s::UInt64, resolution::Int, invert_j::Bool, flip_ij::Bool, do_shift_digits::Bool)
    offset = (0.0, 0.0)
    flips = (NO, NO)
    digits = s
    pattern = flip_ij ? PATTERN_FLIPPED : PATTERN
    for idx in resolution:-1:1
        do_shift_digits && (digits = _shift_digits(digits, idx, flips, invert_j, pattern))
        qflips = _quaternary_to_flips(_quaternary_digit(digits, idx))
        flips = (flips[1] * qflips[1], flips[2] * qflips[2])
    end
    flips = (NO, NO)
    for idx in resolution:-1:1
        offset = _scale2(offset, 2.0)
        digit = _quaternary_digit(digits, idx)
        child_offset = _quaternary_to_kj(digit, flips)
        offset = _add2(offset, child_offset)
        qflips = _quaternary_to_flips(digit)
        flips = (flips[1] * qflips[1], flips[2] * qflips[2])
    end
    q = resolution == 0 ? 0 : _quaternary_digit(digits, 1)
    return Anchor(q, _kj_to_ij(offset), flips)
end

function _ij_to_s(input, resolution::Integer, orientation::Symbol=:uv, do_shift_digits::Bool=true)
    reverse = orientation in (:vu, :wu, :vw)
    invert_j = orientation in (:wv, :vw)
    flip_ij = orientation in (:wu, :uw)
    ij = input
    if flip_ij
        ij = (input[2], input[1])
    end
    if invert_j
        i, j = ij
        ij = (i, (1 << resolution) - (i + j))
    end
    S = __ij_to_s(ij, invert_j, flip_ij, Int(resolution), do_shift_digits)
    if reverse
        S = (UInt64(1) << (2 * resolution)) - S - UInt64(1)
    end
    return S
end

function __ij_to_s(input, invert_j::Bool, flip_ij::Bool, resolution::Int, do_shift_digits::Bool)
    digits = UInt64(0)
    flips = (NO, NO)
    pivot = (0.0, 0.0)
    for i0 in (resolution - 1):-1:0
        idx = i0 + 1
        relative = _sub2(input, pivot)
        scale = 1 << i0
        scaled = _scale2(relative, inv(scale))
        digit = _ij_to_quaternary(scaled, flips)
        digits = _set_quaternary_digit(digits, idx, digit)
        child_offset = _kj_to_ij(_quaternary_to_kj(digit, flips))
        pivot = _add2(pivot, _scale2(child_offset, scale))
        qflips = _quaternary_to_flips(digit)
        flips = (flips[1] * qflips[1], flips[2] * qflips[2])
    end
    pattern = flip_ij ? PATTERN_FLIPPED_REVERSED : PATTERN_REVERSED
    for idx in 1:resolution
        qflips = _quaternary_to_flips(_quaternary_digit(digits, idx))
        flips = (flips[1] * qflips[1], flips[2] * qflips[2])
        do_shift_digits && (digits = _shift_digits(digits, idx, flips, invert_j, pattern))
    end
    return digits
end

function _ij_to_flips(input, resolution::Integer)
    flips = (NO, NO)
    pivot = (0.0, 0.0)
    for i0 in (Int(resolution) - 1):-1:0
        relative = _sub2(input, pivot)
        scale = 1 << i0
        scaled = _scale2(relative, inv(scale))
        digit = _ij_to_quaternary(scaled, flips)
        child_offset = _kj_to_ij(_quaternary_to_kj(digit, flips))
        pivot = _add2(pivot, _scale2(child_offset, scale))
        qflips = _quaternary_to_flips(digit)
        flips = (flips[1] * qflips[1], flips[2] * qflips[2])
    end
    return flips
end

const PROBE_R = 0.1
const PROBE_OFFSETS = (
    (PROBE_R * cos(deg2rad(45.0)), PROBE_R * sin(deg2rad(45.0))),
    (PROBE_R * cos(deg2rad(113.0)), PROBE_R * sin(deg2rad(113.0))),
    (PROBE_R * cos(deg2rad(293.0)), PROBE_R * sin(deg2rad(293.0))),
    (PROBE_R * cos(deg2rad(225.0)), PROBE_R * sin(deg2rad(225.0))),
)

function _anchor_to_s(anchor::Anchor, resolution::Integer, orientation::Symbol=:uv)
    i, j = anchor.offset
    offset_index = 1 - anchor.flips[1] + fld(1 - anchor.flips[2], 2)
    probe = PROBE_OFFSETS[offset_index + 1]
    return _ij_to_s((i + probe[1], j + probe[2]), resolution, orientation)
end

function _anchor_to_triple(anchor::Anchor)
    shift_i = 0.25
    shift_j = 0.25
    flip0, flip1 = anchor.flips
    if flip0 == NO && flip1 == YES
        shift_i = -shift_i
        shift_j = -shift_j
    end
    if flip0 == YES && flip1 == YES
        shift_i = -shift_i
        shift_j = -shift_j
    elseif flip0 == YES
        shift_j -= 1
    elseif flip1 == YES
        shift_j += 1
    end
    i = anchor.offset[1] + shift_i
    j = anchor.offset[2] + shift_j
    r = i + j - 0.5
    c = i - j + r
    return (x=floor(Int, (c + 1) / 2 - r), y=_js_round_int(r), z=floor(Int, (1 - c) / 2))
end

_triple_parity(t) = t.x + t.y + t.z

function _triple_in_bounds(t, max_row::Integer)
    sumv = t.x + t.y + t.z
    (sumv == 0 || sumv == 1) || return false
    limit = t.y - sumv
    return t.x <= 0 && t.z <= 0 && t.y >= 0 && t.y <= max_row &&
           t.x >= -limit && t.z >= -limit
end

function _compute_q(offset, flips, orientation::Symbol=:uv)
    i = Int(offset[1])
    j = Int(offset[2])
    flip0, flip1 = flips
    imod2 = i & 1
    jmod2 = j & 1
    f0idx = (flip0 + 1) >> 1
    f1idx = (flip1 + 1) >> 1
    if orientation == :uw || orientation == :wu
        group2 = (
            (((0, 3), (3, 0)), ((3, 2), (2, 3))),
            (((2, 1), (1, 2)), ((1, 0), (0, 1))),
        )
        return group2[imod2 + 1][jmod2 + 1][f0idx + 1][f1idx + 1]
    end
    if imod2 == 0
        return jmod2 == 0 ? 0 : 2
    end
    odd_i = (
        ((3, 1), (1, 3)),
        ((1, 3), (3, 1)),
    )
    return odd_i[jmod2 + 1][f0idx + 1][f1idx + 1]
end

function _offset_flips_to_anchor(offset, flips, orientation::Symbol=:uv)
    return Anchor(_compute_q(offset, flips, orientation), (Float64(offset[1]), Float64(offset[2])), flips)
end

function _triple_to_anchor(t, resolution::Integer, orientation::Symbol=:uv)
    sumv = t.x + t.y + t.z
    (sumv == 0 || sumv == 1) || return nothing
    r = t.y
    cmin = max(2 * t.x + 2 * r - 1, -2 * t.z - 1 + 0.0001)
    cmax = min(2 * t.x + 2 * r + 1 - 0.0001, 1 - 2 * t.z)
    c = _js_round_int((cmin + cmax) / 2)
    center_i = (c + 0.5) / 2
    center_j = r - c / 2 + 0.25
    if orientation == :uv || orientation == :vu
        flips = _ij_to_flips((center_i, center_j), resolution)
        shift_i = 0.25
        shift_j = 0.25
        if flips[1] == NO && flips[2] == YES
            shift_i = -shift_i
            shift_j = -shift_j
        end
        if flips[1] == YES && flips[2] == YES
            shift_i = -shift_i
            shift_j = -shift_j
        elseif flips[1] == YES
            shift_j -= 1
        elseif flips[2] == YES
            shift_j += 1
        end
        offset = (_js_round_int(center_i - shift_i), _js_round_int(center_j - shift_j))
        return _offset_flips_to_anchor(offset, flips, orientation)
    end
    s = _ij_to_s((center_i, center_j), resolution, orientation)
    return _s_to_anchor(s, resolution, orientation)
end

function _triple_to_s(t, resolution::Integer, orientation::Symbol=:uv)
    anchor = _triple_to_anchor(t, resolution, orientation)
    anchor === nothing && return nothing
    return _anchor_to_s(anchor, resolution, orientation)
end

function _face_to_barycentric(p, tri)
    p1, p2, p3 = tri
    d31 = _sub2(p1, p3)
    d23 = _sub2(p3, p2)
    d3p = _sub2(p, p3)
    det = d23[1] * d31[2] - d23[2] * d31[1]
    b0 = (d23[1] * d3p[2] - d23[2] * d3p[1]) / det
    b1 = (d31[1] * d3p[2] - d31[2] * d3p[1]) / det
    return (b0, b1, 1 - b0 - b1)
end

function _barycentric_to_face(b, tri)
    p1, p2, p3 = tri
    return (
        b[1] * p1[1] + b[2] * p2[1] + b[3] * p3[1],
        b[1] * p1[2] + b[2] * p2[2] + b[3] * p3[2],
    )
end

function _triple_product(A, B, C)
    return _dot3(A, _cross3(B, C))
end

function _spherical_triangle_area(v1, v2, v3)
    midA = _normalize3(_lerp3(v2, v3, 0.5))
    midB = _normalize3(_lerp3(v3, v1, 0.5))
    midC = _normalize3(_lerp3(v1, v2, 0.5))
    S = _triple_product(midA, midB, midC)
    clamped = _clamp1(S)
    abs(clamped) < 1e-8 && return 2 * clamped
    return 2 * asin(clamped)
end

function _vector_difference(A, B)
    midpoint = _normalize3(_lerp3(A, B, 0.5))
    Dvec = _cross3(A, midpoint)
    D = _norm3(Dvec)
    if D < 1e-8
        return 0.5 * _norm3(_sub3(A, B))
    end
    return D
end

function _quadruple_product(A, B, C, D)
    crossCD = _cross3(C, D)
    tripleACD = _dot3(A, crossCD)
    tripleBCD = _dot3(B, crossCD)
    return _sub3(_scale3(B, tripleACD), _scale3(A, tripleBCD))
end

function _slerp(A, B, t)
    gamma = acos(_clamp1(_dot3(A, B) / max(_norm3(A) * _norm3(B), eps(Float64))))
    gamma < 1e-12 && return _lerp3(A, B, t)
    singamma = sin(gamma)
    weightA = sin((1 - t) * gamma) / singamma
    weightB = sin(t * gamma) / singamma
    return _add3(_scale3(A, weightA), _scale3(B, weightB))
end

function _compute_equal_area_constants(tri)
    A, B, C = tri
    c1 = _cross3(B, C)
    c12 = _dot3(B, C)
    return (
        V=_dot3(A, c1),
        c12=c12,
        s12=_norm3(c1),
        kQ=2 / acos(_clamp1(c12)),
        areaABC=_spherical_triangle_area(A, B, C),
    )
end

function _make_crs_vertices()
    vertices = NTuple{3,Float64}[]
    function add_vertex(v)
        normalized = _normalize3(v)
        for existing in vertices
            _norm3(_sub3(normalized, existing)) < 1e-5 && return false
        end
        push!(vertices, normalized)
        return true
    end
    for origin in ORIGINS
        add_vertex(_to_cartesian(origin.axis))
    end
    phi_vertex = atan(DISTANCE_TO_VERTEX)
    for origin in ORIGINS
        for i in 0:4
            theta = ((2 * i + 1) * pi) / 5 + origin.angle
            add_vertex(_quat_rotate(origin.quat, _to_cartesian((theta, phi_vertex))))
        end
    end
    phi_midpoint = atan(DISTANCE_TO_EDGE)
    for origin in ORIGINS
        for i in 0:4
            theta = (2 * i * pi) / 5 + origin.angle
            add_vertex(_quat_rotate(origin.quat, _to_cartesian((theta, phi_midpoint))))
        end
    end
    length(vertices) == 62 || error("Failed to construct CRS: vertices length is $(length(vertices))")
    return Tuple(vertices)
end

const CRS_VERTICES = _make_crs_vertices()
const CANONICAL_TRIANGLE = (CRS_VERTICES[1], CRS_VERTICES[33], CRS_VERTICES[13])
const EQUAL_AREA_CONSTANTS = _compute_equal_area_constants(CANONICAL_TRIANGLE)

function _crs_get_vertex(point)
    normalized = _normalize3(point)
    for vertex in CRS_VERTICES
        _norm3(_sub3(normalized, vertex)) < 1e-5 && return vertex
    end
    error("Failed to find vertex in CRS")
end

function _safe_acos(x)
    x < 1e-3 ? 2 * x + (x^3) / 3 : acos(1 - 2 * x^2)
end

function _equal_area_forward(v, spherical_triangle, face_triangle)
    A, B, C = spherical_triangle
    Z = _normalize3(_sub3(v, A))
    p = _normalize3(_quadruple_product(A, Z, B, C))
    h = _vector_difference(A, v) / _vector_difference(A, p)
    scaled_area = h / EQUAL_AREA_CONSTANTS.areaABC
    b = (
        1 - h,
        scaled_area * _spherical_triangle_area(A, p, C),
        scaled_area * _spherical_triangle_area(A, B, p),
    )
    return _barycentric_to_face(b, face_triangle)
end

function _equal_area_inverse(face_point, face_triangle, spherical_triangle)
    A, B, C = spherical_triangle
    b = _face_to_barycentric(face_point, face_triangle)
    threshold = 1 - 1e-14
    b[1] > threshold && return A
    b[2] > threshold && return B
    b[3] > threshold && return C
    h = 1 - b[1]
    R = b[3] / h
    alpha = R * EQUAL_AREA_CONSTANTS.areaABC
    S = sin(alpha)
    halfC = sin(alpha / 2)
    CC = 2 * halfC * halfC
    c01 = _dot3(A, B)
    c20 = _dot3(C, A)
    c12 = EQUAL_AREA_CONSTANTS.c12
    s12 = EQUAL_AREA_CONSTANTS.s12
    f = S * EQUAL_AREA_CONSTANTS.V + CC * (c01 * c12 - c20)
    g = CC * s12 * (1 + c01)
    q = EQUAL_AREA_CONSTANTS.kQ * atan(g, f)
    P = _slerp(B, C, q)
    K = _vector_difference(A, P)
    t = _safe_acos(h * K) / _safe_acos(K)
    return _slerp(A, P, t)
end

function _normalize_gamma(gamma)
    segment = gamma / TWO_PI_OVER_5
    s_center = round(segment)
    s_offset = segment - s_center
    return s_offset * TWO_PI_OVER_5
end

function _face_triangle_index(polar)
    return mod(floor(Int, polar[2] / PI_OVER_5), 10)
end

function _should_reflect(polar)
    D = _to_face((polar[1], _normalize_gamma(polar[2])))[1]
    return D > DISTANCE_TO_EDGE
end

function _get_face_triangle(face_triangle_index::Integer, reflected::Bool=false, squashed::Bool=false)
    if !reflected
        quintant = mod(fld(face_triangle_index + 1, 2), 5)
        verts = _get_quintant_vertices(quintant)
        vcenter, vcorner1, vcorner2 = verts
        vedge_midpoint = _lerp2(vcorner1, vcorner2, 0.5)
        even = iseven(face_triangle_index)
        return even ? (vcenter, vedge_midpoint, vcorner1) : (vcenter, vcorner2, vedge_midpoint)
    end
    A, B, C = _get_face_triangle(face_triangle_index, false, false)
    even = iseven(face_triangle_index)
    A = _neg2(A)
    midpoint = even ? B : C
    A = _add2(A, _scale2(midpoint, squashed ? 1 + 1 / cos(INTERHEDRAL_ANGLE) : 2))
    return (A, C, B)
end

function _compute_spherical_triangle(face_triangle_index::Integer, origin_id::Integer, reflected::Bool)
    origin = ORIGINS[origin_id + 1]
    face_triangle = _get_face_triangle(face_triangle_index, reflected, true)
    return map(face_triangle) do face
        rho, gamma = _to_polar(face)
        rotated = _to_cartesian((gamma + origin.angle, atan(rho)))
        _crs_get_vertex(_quat_rotate(origin.quat, rotated))
    end
end

function _make_spherical_triangles()
    triangle_type = NTuple{3,NTuple{3,Float64}}
    triangles = Array{triangle_type}(undef, 10, 12, 2)
    for reflected in (false, true), origin_id in 0:11, face_triangle_index in 0:9
        triangles[face_triangle_index + 1, origin_id + 1, Int(reflected) + 1] =
            _compute_spherical_triangle(face_triangle_index, origin_id, reflected)
    end
    return triangles
end

const SPHERICAL_TRIANGLES = _make_spherical_triangles()

@inline function _get_spherical_triangle(
        face_triangle_index::Integer,
        origin_id::Integer,
        reflected::Bool=false,
    )
    return @inbounds SPHERICAL_TRIANGLES[
        face_triangle_index + 1,
        origin_id + 1,
        Int(reflected) + 1,
    ]
end

function _dodeca_forward_cartesian(unprojected, origin_id::Integer)
    origin = ORIGINS[origin_id + 1]
    out = _quat_rotate(origin.inverse_quat, unprojected)
    projected_spherical = _to_spherical(out)
    polar = (tan(projected_spherical[2]), projected_spherical[1] - origin.angle)
    fti = _face_triangle_index(polar)
    reflect = _should_reflect(polar)
    face_triangle = _get_face_triangle(fti, reflect, false)
    spherical_triangle = _get_spherical_triangle(fti, origin_id, reflect)
    return _equal_area_forward(unprojected, spherical_triangle, face_triangle)
end

_dodeca_forward(spherical, origin_id::Integer) = _dodeca_forward_cartesian(_to_cartesian(spherical), origin_id)

function _dodeca_inverse(face, origin_id::Integer)
    polar = _to_polar(face)
    fti = _face_triangle_index(polar)
    reflect = _should_reflect(polar)
    face_triangle = _get_face_triangle(fti, reflect, false)
    spherical_triangle = _get_spherical_triangle(fti, origin_id, reflect)
    return _to_spherical(_equal_area_inverse(face, face_triangle, spherical_triangle))
end

function _face_to_estimate(dodec_point, origin::Origin, resolution::Integer)
    polar = _to_polar(dodec_point)
    quintant = _get_quintant_polar(polar)
    segment, orientation = _quintant_to_segment(quintant, origin)
    if resolution < FIRST_HILBERT_RESOLUTION
        return A5Cell(origin, segment, 0x00, resolution)
    end
    point = dodec_point
    if quintant != 0
        point = _rotate2(point, -TWO_PI_OVER_5 * quintant)
    end
    hilbert_resolution = 1 + resolution - FIRST_HILBERT_RESOLUTION
    point = _scale2(point, 2.0^hilbert_resolution)
    ij = _face_to_ij(point)
    S = _ij_to_s(ij, hilbert_resolution, orientation)
    return A5Cell(origin, segment, S, resolution)
end

function _spherical_to_estimate(spherical, resolution::Integer)
    origin = _find_nearest_origin(spherical)
    dodec_point = _dodeca_forward(spherical, origin.id)
    return _face_to_estimate(dodec_point, origin, resolution)
end

function _cartesian_to_estimate(cartesian, resolution::Integer)
    origin = _find_nearest_origin_cartesian(cartesian)
    dodec_point = _dodeca_forward_cartesian(cartesian, origin.id)
    return _face_to_estimate(dodec_point, origin, resolution)
end

function _get_pentagon(cell::A5Cell)
    quintant, orientation = _segment_to_quintant(cell.segment, cell.origin)
    if cell.resolution == FIRST_HILBERT_RESOLUTION - 1
        return _get_quintant_vertices(quintant)
    elseif cell.resolution == FIRST_HILBERT_RESOLUTION - 2
        return _get_face_vertices()
    end
    hilbert_resolution = cell.resolution - FIRST_HILBERT_RESOLUTION + 1
    anchor = _s_to_anchor(cell.S, hilbert_resolution, orientation)
    return _get_pentagon_vertices(hilbert_resolution, quintant, anchor)
end

function _a5cell_contains_point(cell::A5Cell, spherical)
    pentagon = _get_pentagon(cell)
    projected = _dodeca_forward(spherical, cell.origin.id)
    return _contains_point(pentagon, projected)
end

const FACE_ADJACENCY = (
    ((1, 2), (4, 3), (5, 4), (6, 0), (11, 1)),
    ((2, 3), (4, 4), (0, 0), (11, 0), (10, 1)),
    ((9, 2), (3, 0), (4, 0), (1, 0), (10, 0)),
    ((2, 1), (9, 1), (8, 1), (5, 1), (4, 1)),
    ((2, 2), (3, 4), (5, 0), (0, 1), (1, 1)),
    ((4, 2), (3, 3), (8, 0), (6, 1), (0, 2)),
    ((0, 3), (5, 3), (8, 4), (7, 1), (11, 2)),
    ((11, 3), (6, 3), (8, 3), (9, 4), (10, 3)),
    ((5, 2), (3, 2), (9, 0), (7, 2), (6, 2)),
    ((8, 2), (3, 1), (2, 0), (10, 4), (7, 3)),
    ((2, 4), (1, 4), (11, 4), (7, 4), (9, 3)),
    ((1, 3), (0, 4), (6, 4), (7, 0), (10, 2)),
)

const NEIGHBORS = (
    ((0, -2, -1, 1), (0, -2, -1, -1), (0, -1, 1, -1), (0, -1, -1, -1), (0, -1, 1, 1),
     (1, -2, -1, -1), (1, -1, -1, 1), (1, -1, 1, -1), (1, 0, 1, -1), (2, -1, 1, -1),
     (2, -2, -1, -1)),
    ((-1, -1, -1, 1), (0, -2, -1, -1), (0, -1, -1, -1), (0, -1, 1, -1), (0, 0, -1, 1),
     (0, 0, -1, -1), (0, 1, 1, -1), (0, 1, 1, 1), (1, -2, -1, -1), (1, -1, 1, -1),
     (1, -1, -1, -1), (1, 0, 1, -1)),
    ((-2, 2, -1, -1), (-2, 1, 1, -1), (-1, 0, 1, -1), (-1, 1, 1, -1), (-1, 1, -1, 1),
     (-1, 2, -1, -1), (0, 1, -1, -1), (0, 1, 1, -1), (0, 1, 1, 1), (0, 2, -1, -1),
     (0, 2, -1, 1)),
    ((-1, 0, 1, -1), (-1, 1, 1, -1), (-1, 1, -1, -1), (-1, 2, -1, -1), (0, -1, 1, -1),
     (0, -1, 1, 1), (0, 0, -1, -1), (0, 0, -1, 1), (0, 1, -1, -1), (0, 1, 1, -1),
     (0, 2, -1, -1), (1, 1, -1, 1)),
    ((0, -1, 1, -1), (0, -1, 1, 1), (0, 0, -1, -1), (0, 0, -1, 1), (0, 1, -1, -1),
     (1, 0, -1, -1), (1, 0, 1, -1), (1, -1, 1, -1), (1, 1, -1, 1), (2, -1, 1, -1),
     (2, 0, -1, -1)),
    ((-1, 1, -1, 1), (0, -1, 1, -1), (0, 0, -1, -1), (0, 1, -1, -1), (0, 1, 1, -1),
     (0, 1, 1, 1), (0, 2, -1, -1), (0, 2, -1, 1), (1, -1, 1, -1), (1, 0, -1, -1),
     (1, 0, 1, -1), (1, 1, -1, -1)),
    ((-2, 0, -1, -1), (-2, 1, 1, -1), (-1, -1, -1, 1), (-1, 0, -1, -1), (-1, 0, 1, -1),
     (-1, 1, 1, -1), (0, -1, -1, -1), (0, 0, -1, -1), (0, 0, -1, 1), (0, 1, 1, -1),
     (0, 1, 1, 1)),
    ((-1, -1, -1, -1), (-1, 0, -1, -1), (-1, 0, 1, -1), (-1, 1, 1, -1), (0, -2, -1, -1),
     (0, -2, -1, 1), (0, -1, -1, -1), (0, -1, 1, -1), (0, -1, 1, 1), (0, 0, -1, -1),
     (0, 1, 1, -1), (1, -1, -1, 1)),
)

function _is_neighbor(origin::Anchor, candidate::Anchor)
    origin_flavor = _get_pentagon_flavor(origin)
    candidate_flavor = _get_pentagon_flavor(candidate)
    origin_flavor == candidate_flavor && return false
    relative = (
        Int(candidate.offset[1] - origin.offset[1]),
        Int(candidate.offset[2] - origin.offset[2]),
        candidate.flips[1] * origin.flips[1],
        candidate.flips[2] * origin.flips[2],
    )
    return any(==(relative), NEIGHBORS[origin_flavor + 1])
end

function _find_quintant_neighbor_s(source_triple, uv_source_anchor, source_s::UInt64, resolution::Integer,
        orientation::Symbol, edge_only::Bool)
    max_s = UInt64(4)^UInt64(resolution)
    max_row = (1 << resolution) - 1
    out = UInt64[]
    for dx in -1:1, dy in -1:1, dz in -1:1
        dx == 0 && dy == 0 && dz == 0 && continue
        abs(dx) + abs(dy) + abs(dz) > 3 && continue
        edge_only && abs(dx) + abs(dy) + abs(dz) > 2 && continue
        neighbor_triple = (x=source_triple.x + dx, y=source_triple.y + dy, z=source_triple.z + dz)
        _triple_in_bounds(neighbor_triple, max_row) || continue
        uv_neighbor_anchor = _triple_to_anchor(neighbor_triple, resolution, :uv)
        (uv_neighbor_anchor === nothing || uv_source_anchor === nothing) && continue
        _is_neighbor(uv_source_anchor, uv_neighbor_anchor) || continue
        neighbor_s = _triple_to_s(neighbor_triple, resolution, orientation)
        if neighbor_s !== nothing && neighbor_s < max_s && neighbor_s != source_s
            push!(out, neighbor_s)
        end
    end
    return out
end

function _serialize_res1(origin::Origin, quintant::Integer)
    segment, _ = _quintant_to_segment(quintant, origin)
    return serialize(A5Cell(origin, segment, 0x00, 1))
end

function _get_res0_neighbors(origin::Origin)
    out = Set{UInt64}()
    for q in 0:4
        adjacent_face_id, _ = FACE_ADJACENCY[origin.id + 1][q + 1]
        push!(out, serialize(A5Cell(ORIGINS[adjacent_face_id + 1], 0, 0x00, 0)))
    end
    return sort!(collect(out))
end

function _get_res1_neighbors(origin::Origin, segment::Integer, edge_only::Bool)
    quintant, _ = _segment_to_quintant(segment, origin)
    out = Set{UInt64}()
    left_q = mod(quintant - 1, 5)
    right_q = mod(quintant + 1, 5)
    push!(out, _serialize_res1(origin, left_q))
    push!(out, _serialize_res1(origin, right_q))
    adjacent_face_id, adjacent_quintant = FACE_ADJACENCY[origin.id + 1][quintant + 1]
    adjacent_origin = ORIGINS[adjacent_face_id + 1]
    push!(out, _serialize_res1(adjacent_origin, adjacent_quintant))
    edge_only && return sort!(collect(out))
    push!(out, _serialize_res1(origin, mod(quintant - 2, 5)))
    push!(out, _serialize_res1(origin, mod(quintant + 2, 5)))
    push!(out, _serialize_res1(adjacent_origin, mod(adjacent_quintant - 1, 5)))
    push!(out, _serialize_res1(adjacent_origin, mod(adjacent_quintant + 1, 5)))
    left_adjacent_face_id, left_adjacent_quintant = FACE_ADJACENCY[origin.id + 1][left_q + 1]
    left_adjacent_origin = ORIGINS[left_adjacent_face_id + 1]
    push!(out, _serialize_res1(left_adjacent_origin, left_adjacent_quintant))
    push!(out, _serialize_res1(left_adjacent_origin, mod(left_adjacent_quintant - 1, 5)))
    right_adjacent_face_id, right_adjacent_quintant = FACE_ADJACENCY[origin.id + 1][right_q + 1]
    right_adjacent_origin = ORIGINS[right_adjacent_face_id + 1]
    push!(out, _serialize_res1(right_adjacent_origin, right_adjacent_quintant))
    push!(out, _serialize_res1(right_adjacent_origin, mod(right_adjacent_quintant + 1, 5)))
    return sort!(collect(out))
end

const LEFT_EDGE_DELTAS = (
    ((0, 0, 0, true), (0, 0, 1, false)),
    ((0, 0, 0, true), (0, 1, 0, true), (0, -1, 1, false), (0, 1, -1, false)),
    (),
    ((0, -1, 0, true), (0, 0, -1, false)),
)
const RIGHT_EDGE_DELTAS = (
    ((0, 0, 0, true), (0, 1, 0, true), (-1, 1, 0, false), (1, -1, 0, false)),
    ((0, 0, 0, true), (1, 0, 0, false)),
    ((0, -1, 0, true), (-1, 0, 0, false)),
    (),
)
const CROSS_FACE_DELTAS = (
    ((0, 0, 0, true), (1, 0, 0, true), (1, 0, -1, false)),
    ((0, 0, -1, true), (0, 0, 0, false)),
)

function _push_triple!(out, triple, orientation::Symbol, origin::Origin, segment::Integer, ctx)
    _triple_in_bounds(triple, ctx.max_row) || return nothing
    s = _triple_to_s(triple, ctx.hilbert_resolution, orientation)
    if s !== nothing && s < ctx.max_s
        push!(out, serialize(A5Cell(origin, segment, s, ctx.resolution)))
    end
    return nothing
end

function _push_deltas!(out, base, deltas, edge_only::Bool, orientation::Symbol, origin::Origin, segment::Integer, ctx)
    for (dx, dy, dz, is_edge) in deltas
        edge_only && !is_edge && continue
        _push_triple!(out, (x=base.x + dx, y=base.y + dy, z=base.z + dz), orientation, origin, segment, ctx)
    end
    return nothing
end

function _get_boundary_neighbors(ctx, edge_only::Bool, skip_corners::Bool=false)
    out = UInt64[]
    triple = ctx.triple
    parity = ctx.parity
    source_quintant = ctx.source_quintant
    origin = ctx.origin
    max_row = ctx.max_row
    y_odd = isodd(triple.y)
    delta_index = parity * 2 + (y_odd ? 1 : 0) + 1

    if triple.z == 0
        target_quintant = mod(source_quintant - 1, 5)
        segment, orientation = _quintant_to_segment(target_quintant, origin)
        _push_deltas!(out, (x=0, y=triple.y, z=triple.x), LEFT_EDGE_DELTAS[delta_index],
            edge_only, orientation, origin, segment, ctx)
    end
    if triple.x == 0
        target_quintant = mod(source_quintant + 1, 5)
        segment, orientation = _quintant_to_segment(target_quintant, origin)
        _push_deltas!(out, (x=triple.z, y=triple.y, z=0), RIGHT_EDGE_DELTAS[delta_index],
            edge_only, orientation, origin, segment, ctx)
    end
    if triple.y == max_row
        adjacent_face_id, adjacent_quintant = FACE_ADJACENCY[origin.id + 1][source_quintant + 1]
        adjacent_origin = ORIGINS[adjacent_face_id + 1]
        segment, orientation = _quintant_to_segment(adjacent_quintant, adjacent_origin)
        _push_deltas!(out, (x=triple.z, y=max_row, z=triple.x), CROSS_FACE_DELTAS[parity + 1],
            edge_only, orientation, adjacent_origin, segment, ctx)
    end
    if triple.x == 0 && triple.y == 0 && triple.z == 0
        for q in 0:4
            q == source_quintant && continue
            distance = min(mod(q - source_quintant, 5), mod(source_quintant - q, 5))
            edge_only && distance != 1 && continue
            segment, orientation = _quintant_to_segment(q, origin)
            _push_triple!(out, triple, orientation, origin, segment, ctx)
        end
    end
    if !skip_corners && triple.x == -max_row && triple.y == max_row && triple.z == 0
        prev_quintant = mod(source_quintant - 1, 5)
        prev_adj_face_id, prev_adj_quintant = FACE_ADJACENCY[origin.id + 1][prev_quintant + 1]
        prev_adj_origin = ORIGINS[prev_adj_face_id + 1]
        prev_segment, prev_orientation = _quintant_to_segment(prev_adj_quintant, prev_adj_origin)
        _push_triple!(out, triple, prev_orientation, prev_adj_origin, prev_segment, ctx)

        cross_face_id, cross_quintant = FACE_ADJACENCY[origin.id + 1][source_quintant + 1]
        cross_origin = ORIGINS[cross_face_id + 1]
        next_cross_quintant = mod(cross_quintant + 1, 5)
        cross_segment, cross_orientation = _quintant_to_segment(next_cross_quintant, cross_origin)
        _push_triple!(out, triple, cross_orientation, cross_origin, cross_segment, ctx)
    end
    return out
end

function _get_global_cell_neighbors(cell_id::UInt64; edge_only::Bool=false)
    cell = deserialize(cell_id)
    cell.resolution == 0 && return _get_res0_neighbors(cell.origin)
    cell.resolution == 1 && return _get_res1_neighbors(cell.origin, cell.segment, edge_only)
    hilbert_resolution = cell.resolution - FIRST_HILBERT_RESOLUTION + 1
    source_quintant, source_orientation = _segment_to_quintant(cell.segment, cell.origin)
    anchor = _s_to_anchor(cell.S, hilbert_resolution, source_orientation)
    triple = _anchor_to_triple(anchor)
    uv_source_anchor = _triple_to_anchor(triple, hilbert_resolution, :uv)
    out = Set{UInt64}()
    for neighbor_s in _find_quintant_neighbor_s(triple, uv_source_anchor, cell.S, hilbert_resolution,
            source_orientation, edge_only)
        push!(out, serialize(A5Cell(cell.origin, cell.segment, neighbor_s, cell.resolution)))
    end
    ctx = (
        triple=triple,
        parity=_triple_parity(triple),
        source_quintant=source_quintant,
        origin=cell.origin,
        hilbert_resolution=hilbert_resolution,
        max_s=UInt64(4)^UInt64(hilbert_resolution),
        max_row=(1 << hilbert_resolution) - 1,
        resolution=cell.resolution,
    )
    for neighbor in _get_boundary_neighbors(ctx, edge_only)
        push!(out, neighbor)
    end
    return sort!(collect(out))
end

const SPIRAL_SAMPLE_COUNT = 24
const SPIRAL_SCALE_RAD = 70pi / 180
const ANGLE_STEP_RAD = 1.4

function _spiral_sample(center_spherical, scale_rad, i::Integer)
    c = _to_cartesian(center_spherical)
    q = _quat_rotation_to((0.0, 0.0, 1.0), c)
    angle = (i + 1) * ANGLE_STEP_RAD
    tangent = _quat_rotate(q, (cos(angle), sin(angle), 0.0))
    R = ((i + 1) / (SPIRAL_SAMPLE_COUNT + 1)) * scale_rad
    return _add3(c, _scale3(tangent, R))
end

function _spherical_to_cell(spherical, resolution::Integer)
    resolution == -1 && return WORLD_CELL
    if resolution < FIRST_HILBERT_RESOLUTION
        return serialize(_spherical_to_estimate(spherical, resolution))
    end
    first_estimate = _spherical_to_estimate(spherical, resolution)
    first_key = serialize(first_estimate)
    first_distance = _a5cell_contains_point(first_estimate, spherical)
    first_distance > 0 && return first_key

    hilbert_resolution = 1 + resolution - FIRST_HILBERT_RESOLUTION
    scale = SPIRAL_SCALE_RAD / 2.0^hilbert_resolution
    seen = Set{UInt64}([first_key])
    cells = [(cell_id=first_key, distance=first_distance)]
    for i in 0:(SPIRAL_SAMPLE_COUNT - 1)
        estimate = _cartesian_to_estimate(_spiral_sample(spherical, scale, i), resolution)
        estimate_key = serialize(estimate)
        estimate_key in seen && continue
        push!(seen, estimate_key)
        distance = _a5cell_contains_point(estimate, spherical)
        distance > 0 && return estimate_key
        push!(cells, (cell_id=estimate_key, distance=distance))
    end

    sort!(cells; by=x -> x.distance, rev=true)
    for k in 1:min(3, length(cells))
        for neighbor_key in _get_global_cell_neighbors(cells[k].cell_id)
            neighbor_key in seen && continue
            push!(seen, neighbor_key)
            neighbor_cell = deserialize(neighbor_key)
            distance = _a5cell_contains_point(neighbor_cell, spherical)
            distance > 0 && return neighbor_key
            push!(cells, (cell_id=neighbor_key, distance=distance))
        end
    end

    sort!(cells; by=x -> x.distance, rev=true)
    return first(cells).cell_id
end

function lonlat_to_cell(lon::Real, lat::Real, resolution::Integer)
    res = Int(resolution)
    -1 <= res <= MAX_RESOLUTION ||
        throw(ArgumentError("A5 resolution must be in -1:$MAX_RESOLUTION"))
    id = _spherical_to_cell(_from_lonlat((Float64(lon), Float64(lat))), res)
    get_resolution(id) == res ||
        throw(ArgumentError("A5 resolution $res is not representable at this location"))
    return id
end

function cell_to_lonlat(cell_id::UInt64)
    cell_id == WORLD_CELL && return (0.0, 0.0)
    cell = deserialize(cell_id)
    pentagon = _get_pentagon(cell)
    return _to_lonlat(_dodeca_inverse(_center(pentagon), cell.origin.id))
end

function cell_boundary_cartesian(
        cell_id::UInt64;
        closed_ring::Bool=true,
        segments::Union{Integer,Symbol,String}=:auto,
    )
    cell_id == WORLD_CELL && return NTuple{3,Float64}[]
    cell = deserialize(cell_id)
    nsegments = segments == :auto || segments == "auto" ?
                (cell.resolution <= 6 ? 2^(6 - cell.resolution) : 1) :
                Int(segments)
    pentagon = _get_pentagon(cell)
    split = _split_edges(pentagon, nsegments)
    nvertices = length(split)
    boundary = Vector{NTuple{3,Float64}}(undef, nvertices + Int(closed_ring))
    for i in eachindex(split)
        boundary[i] = _to_cartesian(_dodeca_inverse(split[i], cell.origin.id))
    end
    closed_ring && (boundary[end] = boundary[1])
    reverse!(boundary)
    return boundary
end

function _normalize_longitudes!(contour, npoints::Int=length(contour))
    center = (0.0, 0.0, 0.0)
    for i in 1:npoints
        center = _add3(center, _to_cartesian(_from_lonlat(contour[i])))
    end
    center = _normalize3(center)
    center_lon, center_lat = _to_lonlat(_to_spherical(center))
    if center_lat > 89.99 || center_lat < -89.99
        center_lon = contour[1][1]
    end
    center_lon = _normalize_longitude(center_lon)
    for i in 1:npoints
        lon, lat = contour[i]
        while lon - center_lon > 180
            lon -= 360
        end
        while lon - center_lon < -180
            lon += 360
        end
        contour[i] = (lon, lat)
    end
    return contour
end

function cell_boundary(cell_id::UInt64; closed_ring::Bool=true, segments::Union{Integer,Symbol,String}=:auto)
    cell_id == WORLD_CELL && return NTuple{2,Float64}[]
    cell = deserialize(cell_id)
    nsegments = segments == :auto || segments == "auto" ?
                (cell.resolution <= 6 ? 2^(6 - cell.resolution) : 1) :
                Int(segments)
    pentagon = _get_pentagon(cell)
    split = _split_edges(pentagon, nsegments)
    nvertices = length(split)
    boundary = Vector{NTuple{2,Float64}}(undef, nvertices + Int(closed_ring))
    for i in eachindex(split)
        boundary[i] = _to_lonlat(_dodeca_inverse(split[i], cell.origin.id))
    end
    _normalize_longitudes!(boundary, nvertices)
    closed_ring && (boundary[end] = boundary[1])
    reverse!(boundary)
    return boundary
end

end
