# Fixed-level Eisenstein-lattice displacements. These retain the hydrology
# integration's current arithmetic contract.
# TODO: Replace these with a seam-safe origin/destination neighbor-step token.

"""
    HexIndex(i, j, k)

Signed cube coordinates for an IGeo7 lattice displacement. Valid coordinates
satisfy `i + j + k == 0`.
"""
struct HexIndex
    i::Int64
    j::Int64
    k::Int64

    function HexIndex(i::Integer, j::Integer, k::Integer)
        ii, jj, kk = Int64(i), Int64(j), Int64(k)
        Int128(ii) + Int128(jj) + Int128(kk) == 0 ||
            throw(ArgumentError("hex coordinates must satisfy i + j + k == 0"))
        return new(ii, jj, kk)
    end
end

"""
    RelativeIGEO7Index(a, b, level)
    RelativeIGEO7Index(hex, level)
    RelativeIGEO7Index(level)

Signed displacement `a + b*omega` in one fixed-level IGeo7 lattice frame.
"""
struct RelativeIGEO7Index
    a::Int64
    b::Int64
    level::UInt8

    function RelativeIGEO7Index(a::Integer, b::Integer, level::Integer)
        0 <= level <= MAX_RESOLUTION ||
            throw(ArgumentError("IGeo7 level must be in 0:$MAX_RESOLUTION"))
        return new(Int64(a), Int64(b), UInt8(level))
    end
end

RelativeIGEO7Index(level::Integer) = RelativeIGEO7Index(0, 0, level)
Base.zero(d::RelativeIGEO7Index) = RelativeIGEO7Index(d.level)
Base.iszero(d::RelativeIGEO7Index) = iszero(d.a) && iszero(d.b)

function RelativeIGEO7Index(hex::HexIndex, level::Integer)
    b = Int128(hex.i) + Int128(hex.j)
    return RelativeIGEO7Index(hex.i, Int64(b), level)
end

function Base.convert(::Type{HexIndex}, d::RelativeIGEO7Index)
    return HexIndex(d.a, Int128(d.b) - Int128(d.a), -Int128(d.b))
end

Base.show(io::IO, d::RelativeIGEO7Index) =
    print(io, "RelativeIGEO7Index(", repr(convert(HexIndex, d)), ", ",
        Int(d.level), ")")

@inline function _same_level(a::RelativeIGEO7Index, b::RelativeIGEO7Index)
    a.level == b.level ||
        throw(DimensionMismatch("IGeo7 displacements have different levels"))
    return Int(a.level)
end

Base.:-(d::RelativeIGEO7Index) = RelativeIGEO7Index(
    Base.Checked.checked_neg(d.a), Base.Checked.checked_neg(d.b), d.level)

function Base.:+(a::RelativeIGEO7Index, b::RelativeIGEO7Index)
    l = _same_level(a, b)
    return RelativeIGEO7Index(
        Base.Checked.checked_add(a.a, b.a),
        Base.Checked.checked_add(a.b, b.b),
        l,
    )
end

function Base.:-(a::RelativeIGEO7Index, b::RelativeIGEO7Index)
    l = _same_level(a, b)
    return RelativeIGEO7Index(
        Base.Checked.checked_sub(a.a, b.a),
        Base.Checked.checked_sub(a.b, b.b),
        l,
    )
end

Base.:*(n::Integer, d::RelativeIGEO7Index) = RelativeIGEO7Index(
    Base.Checked.checked_mul(Int64(n), d.a),
    Base.Checked.checked_mul(Int64(n), d.b),
    d.level,
)
Base.:*(d::RelativeIGEO7Index, n::Integer) = n * d

const DIRECTION_CODES = (
    0x05, 0x04, 0xff, 0xff,
    0x06, 0x00, 0x03, 0xff,
    0xff, 0x01, 0x02, 0xff,
    0xff, 0xff, 0xff, 0xff,
)

"""
    directioncode(displacement) -> UInt8

Code `1:6` for the six counterclockwise IGeo7 unit directions and `0` for the
zero displacement.
"""
@inline function directioncode(d::RelativeIGEO7Index)
    (-1 <= d.a <= 1 && -1 <= d.b <= 1) ||
        throw(ArgumentError("displacement does not point to an immediate neighbor"))
    key = ((Int(d.a) + 1) << 2) | (Int(d.b) + 1)
    code = @inbounds DIRECTION_CODES[key+1]
    code != 0xff ||
        throw(ArgumentError("displacement does not point to an immediate neighbor"))
    return code
end

@inline function _raw_lattice(c::Z7Cell)
    l = level(c)
    a = Int64(0)
    b = Int64(0)
    for k in 1:l
        a, b = mulchi(a, b, k)
        digit = _z7_digit(c.id, k)
        if digit != 0
            da, db = @inbounds SIGMA_AB[digit]
            a += da
            b += db
        end
    end
    return (a, b)
end

function _shared_nonpentagon_frame(a::Z7Cell, b::Z7Cell)
    l = level(a)
    l == level(b) || return false
    z7_base_cell(a.id) == z7_base_cell(b.id) || return false
    prefix_level = 0
    for k in 1:l
        _z7_digit(a.id, k) == _z7_digit(b.id, k) || break
        prefix_level = k
    end
    return !z7_is_pentagon(z7_parent(a.id, prefix_level))
end

function Base.:-(a::Z7Cell, b::Z7Cell)
    l = level(a)
    l == level(b) ||
        throw(DimensionMismatch("IGeo7 cells have different levels"))
    a == b && return RelativeIGEO7Index(l)
    _shared_nonpentagon_frame(a, b) || throw(ArgumentError(
        "IGeo7 cells do not share a supported non-pentagon lattice frame"))
    aa, ab = _raw_lattice(a)
    ba, bb = _raw_lattice(b)
    return RelativeIGEO7Index(
        Base.Checked.checked_sub(aa, ba),
        Base.Checked.checked_sub(ab, bb),
        l,
    )
end

function _decode_raw(base::Int, a::Int64, b::Int64, l::Int)
    z = @inbounds res0_cells()[base+1]
    for k in l:-1:1
        digit, a, b = decode_step(a, b, k)
        z = _z7_set_digit(z, k, digit)
    end
    (a == 0 && b == 0 && is_valid_cell(z)) || return nothing
    return Z7Cell(z)
end

"""
    trytranslate(cell, displacement) -> Union{Z7Cell,Nothing}

Translate a cell in its local non-pentagon lattice frame. Returns `nothing`
when the displacement leaves that frame.
"""
function trytranslate(c::Z7Cell, d::RelativeIGEO7Index)
    l = level(c)
    l == d.level ||
        throw(DimensionMismatch("IGeo7 cell and displacement have different levels"))
    a, b = _raw_lattice(c)
    ta, tb = Int128(a) + Int128(d.a), Int128(b) + Int128(d.b)
    typemin(Int64) <= ta <= typemax(Int64) || return nothing
    typemin(Int64) <= tb <= typemax(Int64) || return nothing
    target = _decode_raw(z7_base_cell(c.id), Int64(ta), Int64(tb), l)
    target === nothing && return nothing
    target == c && return iszero(d) ? target : nothing
    _shared_nonpentagon_frame(c, target) || return nothing
    aa, ab = _raw_lattice(target)
    aa - a == d.a && ab - b == d.b || return nothing
    return target
end

function Base.:+(c::Z7Cell, d::RelativeIGEO7Index)
    target = trytranslate(c, d)
    target === nothing && throw(DomainError((c, d),
        "IGeo7 translation leaves the supported non-pentagon lattice frame"))
    return target
end

Base.:+(d::RelativeIGEO7Index, c::Z7Cell) = c + d
Base.:-(c::Z7Cell, d::RelativeIGEO7Index) = c + (-d)
