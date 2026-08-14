# Six dense roots followed by base-9, row-major child digits.  The integer
# ordinal is merely the compact representation of the published prefix SUID:
# at level l, `root * 9^l + digits_base9`.

const ROOT_NAMES = ('N', 'O', 'P', 'Q', 'R', 'S')
const MAX_LEVEL = 19
const N_SIDE = 3
const N_CHILDREN = 9

"""
    RHEALPixCell(suid)
    RHEALPixCell(level, ordinal0)

An isbits rHEALPix cell identity.  `ordinal0` is the dense zero-based ordinal
at `level`; [`suid`](@ref) exposes the equivalent six-root prefix identifier.
The two representations are bijective for supported levels `0:19`.
"""
struct RHEALPixCell <: DGG.AbstractCellIndex
    level::Int32
    ordinal::Int64
end

RHEALPixCell(level::Integer, ordinal::Integer) =
    RHEALPixCell(Int32(level), Int64(ordinal))
RHEALPixCell(text::AbstractString) = parse_suid(text)

DGG.level(c::RHEALPixCell) = Int(c.level)
DGG.rawid(c::RHEALPixCell) = c.ordinal
Base.isless(a::RHEALPixCell, b::RHEALPixCell) =
    isless((a.level, a.ordinal), (b.level, b.ordinal))

@inline _pow9(level::Integer) = Int64(9)^Int(level)
@inline _ncells(level::Integer) = 6 * _pow9(level)

"""Return the published root-plus-base-9 string identifier for `cell`."""
function suid(cell::RHEALPixCell)
    l = DGG.level(cell)
    0 <= l <= MAX_LEVEL || throw(ArgumentError("rHEALPix level $l is outside 0:$MAX_LEVEL"))
    scale = _pow9(l)
    root, rest = divrem(cell.ordinal, scale)
    0 <= root < 6 || throw(ArgumentError(
        "rHEALPix ordinal $(cell.ordinal) is invalid at level $l"))
    bytes = Vector{UInt8}(undef, l + 1)
    bytes[1] = UInt8(ROOT_NAMES[root + 1])
    for i in (l + 1):-1:2
        rest, digit = divrem(rest, 9)
        bytes[i] = UInt8('0') + UInt8(digit)
    end
    return String(bytes)
end

"""Parse a published rHEALPix SUID such as `"N"`, `"R7"`, or `"S08"`."""
function parse_suid(text::AbstractString)
    bytes = codeunits(text)
    isempty(bytes) && throw(ArgumentError("an rHEALPix SUID cannot be empty"))
    root = findfirst(==(Char(bytes[1])), ROOT_NAMES)
    root === nothing && throw(ArgumentError(
        "rHEALPix root must be one of N, O, P, Q, R, S"))
    l = length(bytes) - 1
    l <= MAX_LEVEL || throw(ArgumentError(
        "rHEALPix SUID level $l exceeds the Int64 limit $MAX_LEVEL"))
    ordinal = Int64(root - 1)
    for byte in @view bytes[2:end]
        UInt8('0') <= byte <= UInt8('8') || throw(ArgumentError(
            "rHEALPix child digits must lie in 0:8"))
        ordinal = 9 * ordinal + (byte - UInt8('0'))
    end
    return RHEALPixCell(l, ordinal)
end

Base.string(c::RHEALPixCell) = suid(c)
Base.show(io::IO, c::RHEALPixCell) = print(io, "RHEALPixCell(\"", suid(c), "\")")

"""
    RHEALPixSystem(; north_square=0, south_square=0, longitude_origin=0)

rHEALPix on the unit authalic sphere with aperture 9.  The two square
parameters choose positions `0:3` for the assembled polar squares; longitude
origin is in radians.  Parameters are part of system identity.
"""
struct RHEALPixSystem <: DGG.AbstractHierarchicalGridSystem
    north_square::UInt8
    south_square::UInt8
    longitude_origin::Float64

    function RHEALPixSystem(north_square::UInt8, south_square::UInt8,
            longitude_origin::Float64)
        north_square <= 0x03 || throw(ArgumentError("north_square must lie in 0:3"))
        south_square <= 0x03 || throw(ArgumentError("south_square must lie in 0:3"))
        return new(north_square, south_square, _wrap_longitude(longitude_origin))
    end
end

function RHEALPixSystem(; north_square::Integer=0, south_square::Integer=0,
        longitude_origin::Real=0)
    north_square in 0:3 || throw(ArgumentError("north_square must lie in 0:3"))
    south_square in 0:3 || throw(ArgumentError("south_square must lie in 0:3"))
    return RHEALPixSystem(UInt8(north_square), UInt8(south_square),
        _wrap_longitude(longitude_origin))
end

RHEALPixSystem(north_square::Integer, south_square::Integer;
    longitude_origin::Real=0) = RHEALPixSystem(;
        north_square, south_square, longitude_origin)

"""
    AusPIXSystem()

The AusPIX profile: WGS84 geodetic geometry over the Greenwich `(0,0)`
rHEALPix kernel, with `N_side=3`.  This returns the package's standard
[`AuthalicSystem`](@ref), so hierarchy and SUIDs are exactly those of
`RHEALPixSystem()` while latitude I/O is correctly warped to WGS84 geodetic
latitude.
"""
AusPIXSystem() = DGG.AuthalicSystem(RHEALPixSystem(), DGG.Helpers.WGS84_AUTHALIC)

const LevelGrid = DGG.HierarchicalLevelGrid{RHEALPixSystem}

DGG.cellindextype(::RHEALPixSystem) = RHEALPixCell
DGG.cellindextypes(::RHEALPixSystem) = (RHEALPixCell,)
DGG.levels(::RHEALPixSystem) = 0:MAX_LEVEL
DGG.has_sorted_subtrees(::RHEALPixSystem) = true
DGG.max_neighbors(::RHEALPixSystem, ::DGG.Edge) = 4
DGG.max_neighbors(::RHEALPixSystem, ::DGG.Vertex) = 8
DGG.rootcells(::RHEALPixSystem) = [RHEALPixCell(0, i) for i in 0:5]

function Base.parent(::RHEALPixSystem, c::RHEALPixCell)
    l = DGG.level(c)
    l > 0 || throw(ArgumentError("level-0 rHEALPix cell $(suid(c)) has no parent"))
    return RHEALPixCell(l - 1, c.ordinal ÷ 9)
end

function DGG.children(sys::RHEALPixSystem, c::RHEALPixCell)
    l = DGG.level(c)
    l < DGG.max_level(sys) || throw(ArgumentError(
        "rHEALPix cell $(suid(c)) is at max_level $(DGG.max_level(sys))"))
    base = 9 * c.ordinal
    return [RHEALPixCell(l + 1, base + digit) for digit in 0:8]
end

function DGG.ancestor(::RHEALPixSystem, c::RHEALPixCell, level::Integer)
    target, own = Int(level), DGG.level(c)
    0 <= target <= own || throw(ArgumentError(
        "ancestor level $target must lie in 0:$own"))
    return RHEALPixCell(target, c.ordinal ÷ _pow9(own - target))
end

function DGG.descendant_range(sys::RHEALPixSystem, c::RHEALPixCell,
        level::Integer)
    target, own = Int(level), DGG.level(c)
    target >= own || throw(ArgumentError(
        "descendant level $target is above cell level $own"))
    target <= DGG.max_level(sys) || throw(ArgumentError(
        "descendant level $target exceeds max_level $(DGG.max_level(sys))"))
    scale = _pow9(target - own)
    lo = c.ordinal * scale
    hi = (c.ordinal + 1) * scale - 1
    return Int(lo + 1):Int(hi + 1)
end

function DGG.descendants(sys::RHEALPixSystem, c::RHEALPixCell, level::Integer)
    target = Int(level)
    r = DGG.descendant_range(sys, c, target)
    return [RHEALPixCell(target, i - 1) for i in r]
end

DGG.ncells(::RHEALPixSystem, level::Integer) = Int(_ncells(level))
DGG.cellindex(::RHEALPixSystem, level::Integer, i::Int) = RHEALPixCell(level, i - 1)

function DGG.cellposition(::RHEALPixSystem, c::RHEALPixCell)
    l = DGG.level(c)
    l in 0:MAX_LEVEL || return nothing
    0 <= c.ordinal < _ncells(l) || return nothing
    return Int(c.ordinal + 1)
end

@inline function _checked_index(c::RHEALPixCell)
    l = DGG.level(c)
    l in 0:MAX_LEVEL || throw(ArgumentError(
        "rHEALPix level $l is outside 0:$MAX_LEVEL"))
    0 <= c.ordinal < _ncells(l) || throw(ArgumentError(
        "rHEALPix ordinal $(c.ordinal) is invalid at level $l"))
    return c.ordinal
end

@inline function _checked_index(g::LevelGrid, c::RHEALPixCell)
    DGG.level(c) == g.level || throw(ArgumentError(
        "cell $(suid(c)) is at level $(DGG.level(c)), not grid level $(g.level)"))
    return _checked_index(c)
end

# Root upper-left corner.  The equatorial roots are consecutive; the polar
# roots move with their placement parameters.
@inline function _root_ul(sys::RHEALPixSystem, root::Int)
    root == 0 && return (-Float64(pi) + Int(sys.north_square) * HALFPI,
                         THREEQUARTERPI)
    root == 5 && return (-Float64(pi) + Int(sys.south_square) * HALFPI,
                         -QUARTERPI)
    return (-Float64(pi) + (root - 1) * HALFPI, QUARTERPI)
end

"""Return `(upper_left_x, upper_left_y, width)` for a cell's planar square."""
function cell_rectangle(sys::RHEALPixSystem, c::RHEALPixCell)
    _checked_index(c)
    l = DGG.level(c)
    scale = _pow9(l)
    root, rest = divrem(c.ordinal, scale)
    ulx, uly = _root_ul(sys, Int(root))
    width = HALFPI
    divisor = scale
    for _ in 1:l
        divisor ÷= 9
        digit, rest = divrem(rest, divisor)
        width /= 3
        ulx += (digit % 3) * width
        uly -= (digit ÷ 3) * width
    end
    return ulx, uly, width
end

@inline _project(sys::RHEALPixSystem, lon, lat) = rhealpix_forward(lon, lat;
    north_square=Int(sys.north_square), south_square=Int(sys.south_square),
    longitude_origin=sys.longitude_origin)

@inline _unproject(sys::RHEALPixSystem, x, y) = rhealpix_inverse(x, y;
    north_square=Int(sys.north_square), south_square=Int(sys.south_square),
    longitude_origin=sys.longitude_origin)

@inline function _sphere_point(sys::RHEALPixSystem, x::Real, y::Real)
    lon, lat = _unproject(sys, x, y)
    coslat = cos(lat)
    return GO.UnitSphericalPoint(coslat * cos(lon), coslat * sin(lon), sin(lat))
end

# As in the package's HEALPix implementation, eight great-circle segments per
# inverse-projected chart edge keep the finite polygon useful for predicates
# and rendering while `cell_area` remains the exact analytic area.
const BOUNDARY_SEGMENTS = 8

function _perimeter_points(sys::RHEALPixSystem, c::RHEALPixCell,
        segments::Integer=BOUNDARY_SEGMENTS)
    ulx, uly, width = cell_rectangle(sys, c)
    n = Int(segments)
    n > 0 || throw(ArgumentError("segments must be positive"))
    points = Vector{GO.UnitSphericalPoint{Float64}}(undef, 4n)
    k = 0
    # Upper-left -> lower-left -> lower-right -> upper-right -> upper-left is
    # CCW viewed from outside after the orientation-preserving inverse map.
    for i in 0:(n - 1)
        points[k += 1] = _sphere_point(sys, ulx, uly - width * i / n)
    end
    for i in 0:(n - 1)
        points[k += 1] = _sphere_point(sys, ulx + width * i / n, uly - width)
    end
    for i in 0:(n - 1)
        points[k += 1] = _sphere_point(sys, ulx + width, uly - width + width * i / n)
    end
    for i in 0:(n - 1)
        points[k += 1] = _sphere_point(sys, ulx + width - width * i / n, uly)
    end
    return points
end

"""Inverse image of the four chart edges, densified to 8 geodesic chords each."""
DGG.cell_boundary(sys::RHEALPixSystem, c::RHEALPixCell) =
    _perimeter_points(sys, c)

"""The inverse image of the planar square centre (the rHEALPix nucleus)."""
function DGG.cell_centroid(sys::RHEALPixSystem, c::RHEALPixCell)
    ulx, uly, width = cell_rectangle(sys, c)
    return _sphere_point(sys, ulx + width / 2, uly - width / 2)
end

function DGG.cell_area(g::LevelGrid, c::RHEALPixCell)
    _checked_index(g, c)
    return 2 * Float64(pi) / (3 * _pow9(g.level))
end

function DGG.node_extent(sys::RHEALPixSystem, c::RHEALPixCell)
    center = DGG.cell_centroid(sys, c)
    points = _perimeter_points(sys, c, 8)
    radius = 0.0
    gap = 0.0
    previous = points[end]
    for p in points
        radius = max(radius, US.spherical_distance(center, p))
        gap = max(gap, US.spherical_distance(previous, p))
        previous = p
    end
    return SphericalCap(center, nextfloat(min(Float64(pi), radius + gap / 2)))
end

@inline function _root_for_plane(x::Float64, y::Float64)
    y > QUARTERPI && return 0
    y < -QUARTERPI && return 5
    x < -HALFPI && return 1
    x < 0 && return 2
    x < HALFPI && return 3
    return 4
end

function _plane_to_cell(sys::RHEALPixSystem, x::Real, y::Real, level::Integer)
    l = Int(level)
    X, Y = Float64(x), Float64(y)
    root = _root_for_plane(X, Y)
    ulx, uly = _root_ul(sys, root)
    dx = clamp((X - ulx) / HALFPI, 0.0, 1.0)
    dy = clamp((uly - Y) / HALFPI, 0.0, 1.0)
    ordinal = Int64(root)
    for _ in 1:l
        dx *= 3
        dy *= 3
        col = min(floor(Int, dx), 2)
        row = min(floor(Int, dy), 2)
        ordinal = 9 * ordinal + 3row + col
        dx -= col
        dy -= row
    end
    return RHEALPixCell(l, ordinal)
end

function DGG.cellat(g::LevelGrid, p::GO.UnitSphericalPoint)
    lon = atan(p[2], p[1])
    lat = asin(clamp(Float64(p[3]), -1.0, 1.0))
    x, y = _project(g.system, lon, lat)
    return _plane_to_cell(g.system, x, y, g.level)
end
