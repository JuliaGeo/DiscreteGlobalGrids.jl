# Fixed-resolution cell identities and local Eisenstein-lattice displacement.

"""
    IGEO7Index(id)
    IGEO7Index(z7_string)

Canonical IGEO7 cell identity. Integer construction validates the Z7 encoding
and rejects resolution 20, which has no cell geometry.
"""
struct IGEO7Index
    id::UInt64

    function IGEO7Index(id::Integer)
        z = UInt64(id)
        get_resolution(z)
        return new(z)
    end
end

IGEO7Index(text::AbstractString) = IGEO7Index(z7_from_string(text))

Base.convert(::Type{UInt64}, index::IGEO7Index) = index.id
Base.convert(::Type{IGEO7Index}, id::Integer) = IGEO7Index(id)
Base.UInt64(index::IGEO7Index) = index.id
Base.isless(a::IGEO7Index, b::IGEO7Index) = isless(a.id, b.id)
Base.show(io::IO, index::IGEO7Index) =
    print(io, "IGEO7Index(", repr(z7_to_string(index.id)), ")")
get_resolution(index::IGEO7Index) = get_resolution(index.id)

"""
    HexIndex(i, j, k)

Signed cube coordinates for an IGEO7 lattice displacement. Valid coordinates
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
    RelativeIGEO7Index(a, b, resolution)
    RelativeIGEO7Index(hex, resolution)
    RelativeIGEO7Index(resolution)

Signed displacement `a + b*omega` in the raw GBT/Eisenstein lattice at one
IGEO7 resolution. The one-argument form constructs the zero displacement.
"""
struct RelativeIGEO7Index
    a::Int64
    b::Int64
    resolution::UInt8

    function RelativeIGEO7Index(a::Integer, b::Integer, resolution::Integer)
        0 <= resolution <= MAX_RESOLUTION ||
            throw(ArgumentError("IGEO7 resolution must be in 0:$MAX_RESOLUTION"))
        return new(Int64(a), Int64(b), UInt8(resolution))
    end
end

RelativeIGEO7Index(resolution::Integer) = RelativeIGEO7Index(0, 0, resolution)
Base.zero(d::RelativeIGEO7Index) = RelativeIGEO7Index(d.resolution)
Base.iszero(d::RelativeIGEO7Index) = iszero(d.a) && iszero(d.b)

function RelativeIGEO7Index(hex::HexIndex, resolution::Integer)
    b = Int128(hex.i) + Int128(hex.j)
    return RelativeIGEO7Index(hex.i, Int64(b), resolution)
end

function Base.convert(::Type{HexIndex}, d::RelativeIGEO7Index)
    return HexIndex(d.a, Int128(d.b) - Int128(d.a), -Int128(d.b))
end
Base.show(io::IO, d::RelativeIGEO7Index) =
    print(io, "RelativeIGEO7Index(", repr(convert(HexIndex, d)), ", ",
        Int(d.resolution), ")")

@inline function _same_resolution(a::RelativeIGEO7Index, b::RelativeIGEO7Index)
    a.resolution == b.resolution ||
        throw(DimensionMismatch("IGEO7 displacements have different resolutions"))
    return Int(a.resolution)
end

Base.:-(d::RelativeIGEO7Index) = RelativeIGEO7Index(
    Base.Checked.checked_neg(d.a), Base.Checked.checked_neg(d.b), d.resolution)

function Base.:+(a::RelativeIGEO7Index, b::RelativeIGEO7Index)
    resolution = _same_resolution(a, b)
    return RelativeIGEO7Index(
        Base.Checked.checked_add(a.a, b.a),
        Base.Checked.checked_add(a.b, b.b),
        resolution,
    )
end

function Base.:-(a::RelativeIGEO7Index, b::RelativeIGEO7Index)
    resolution = _same_resolution(a, b)
    return RelativeIGEO7Index(
        Base.Checked.checked_sub(a.a, b.a),
        Base.Checked.checked_sub(a.b, b.b),
        resolution,
    )
end

Base.:*(n::Integer, d::RelativeIGEO7Index) = RelativeIGEO7Index(
    Base.Checked.checked_mul(Int64(n), d.a),
    Base.Checked.checked_mul(Int64(n), d.b),
    d.resolution,
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

Direction code for a unit IGEO7 displacement: `1:6` follow the six
counterclockwise Eisenstein units and `0` means no displacement. Throws
`ArgumentError` for a displacement longer than one cell.

The implementation packs `(a + 1, b + 1)` into a four-bit key and performs one
constant-table lookup; it does not construct a `HexIndex` or search.
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

@inline function _raw_lattice(index::IGEO7Index)
    resolution = get_resolution(index)
    a = Int64(0)
    b = Int64(0)
    for level in 1:resolution
        a, b = mulchi(a, b, level)
        digit = _z7_digit(index.id, level)
        if digit != 0
            da, db = @inbounds SIGMA_AB[digit]
            a += da
            b += db
        end
    end
    return (a, b)
end

function _shared_nonpentagon_frame(a::IGEO7Index, b::IGEO7Index)
    resolution = get_resolution(a)
    resolution == get_resolution(b) || return false
    z7_base_cell(a.id) == z7_base_cell(b.id) || return false
    prefix_resolution = 0
    for level in 1:resolution
        _z7_digit(a.id, level) == _z7_digit(b.id, level) || break
        prefix_resolution = level
    end
    ancestor = z7_parent(a.id, prefix_resolution)
    return !z7_is_pentagon(ancestor)
end

function Base.:-(a::IGEO7Index, b::IGEO7Index)
    resolution = get_resolution(a)
    resolution == get_resolution(b) ||
        throw(DimensionMismatch("IGEO7 cells have different resolutions"))
    a == b && return RelativeIGEO7Index(resolution)
    _shared_nonpentagon_frame(a, b) || throw(ArgumentError(
        "IGEO7 cells do not share a supported non-pentagon lattice frame"))
    aa, ab = _raw_lattice(a)
    ba, bb = _raw_lattice(b)
    return RelativeIGEO7Index(
        Base.Checked.checked_sub(aa, ba),
        Base.Checked.checked_sub(ab, bb),
        resolution,
    )
end

function _decode_raw(base::Int, a::Int64, b::Int64, resolution::Int)
    z = @inbounds res0_cells()[base+1]
    for level in resolution:-1:1
        digit, a, b = decode_step(a, b, level)
        z = _z7_set_digit(z, level, digit)
    end
    (a == 0 && b == 0 && is_valid_cell(z)) || return nothing
    return IGEO7Index(z)
end

"""
    trytranslate(index, displacement) -> Union{IGEO7Index,Nothing}

Translate `index` in its local non-pentagon lattice frame. Returns `nothing`
when the displacement leaves the supported frame.
"""
function trytranslate(index::IGEO7Index, d::RelativeIGEO7Index)
    resolution = get_resolution(index)
    resolution == d.resolution ||
        throw(DimensionMismatch("IGEO7 cell and displacement have different resolutions"))
    a, b = _raw_lattice(index)
    ta, tb = Int128(a) + Int128(d.a), Int128(b) + Int128(d.b)
    typemin(Int64) <= ta <= typemax(Int64) || return nothing
    typemin(Int64) <= tb <= typemax(Int64) || return nothing
    target = _decode_raw(z7_base_cell(index.id), Int64(ta), Int64(tb), resolution)
    target === nothing && return nothing
    target == index && return iszero(d.a) && iszero(d.b) ? target : nothing
    _shared_nonpentagon_frame(index, target) || return nothing
    aa, ab = _raw_lattice(target)
    aa - a == d.a && ab - b == d.b || return nothing
    return target
end

function Base.:+(index::IGEO7Index, d::RelativeIGEO7Index)
    target = trytranslate(index, d)
    target === nothing && throw(DomainError((index, d),
        "IGEO7 translation leaves the supported non-pentagon lattice frame"))
    return target
end

Base.:+(d::RelativeIGEO7Index, index::IGEO7Index) = index + d
Base.:-(index::IGEO7Index, d::RelativeIGEO7Index) = index + (-d)

"""
    neighbors(index)
    neighbors(source, index)

Edge neighbors of `index`. The one-argument form returns its global
neighborhood. `source` may instead be `IGEO7DGGS()` for the same result or a
one-dimensional IGEO7 `DimArray`/`Raster` for its stored subset.
"""
function neighbors end

"""
    celldistance(A, from, to) -> Float64

Great-circle distance in metres between the centers of two cells stored in a
one-dimensional IGEO7 `DimArray`/`Raster`, measured on the WGS84 authalic
sphere.
"""
function celldistance end

"""
    cellarea(A, index) -> Float64

Area in square metres of a cell stored in a one-dimensional IGEO7
`DimArray`/`Raster`.
"""
function cellarea end

"""
    edges(A) -> Vector{IGEO7Index}

Cells on the stored coverage boundary: every returned cell has at least one
global edge neighbor not stored in `A`.
"""
function edges end

"""
    cell_to_position(A, index) -> Int

Position of `index` in a one-dimensional IGEO7 `DimArray`/`Raster`. Throws
`BoundsError` when the cell is not stored.
"""
function cell_to_position end

"""
    position_to_cell(A, position) -> IGEO7Index

Canonical cell identity at `position` in a one-dimensional IGEO7
`DimArray`/`Raster`.
"""
function position_to_cell end

"""
    IGEO7Indices

Lazy index domain returned by `eachindex` for a one-dimensional IGEO7 array.
It wraps the lookup IDs without materializing `IGEO7Index` objects and caches
whether those IDs occupy one contiguous global-ordinal interval.
"""
struct IGEO7Indices{V<:AbstractVector{UInt64}} <: AbstractVector{IGEO7Index}
    ids::V
    resolution::Int
    first_ordinal::Int
    contiguous::Bool

    function IGEO7Indices(ids::V, resolution::Integer) where {V<:AbstractVector{UInt64}}
        res = Int(resolution)
        0 <= res <= MAX_RESOLUTION ||
            throw(ArgumentError("IGEO7 resolution must be in 0:$MAX_RESOLUTION"))
        isempty(ids) && return new{V}(ids, res, 0, true)
        first_ordinal = cell_to_index(first(ids))
        last_ordinal = cell_to_index(last(ids))
        contiguous = last_ordinal - first_ordinal + 1 == length(ids)
        return new{V}(ids, res, first_ordinal, contiguous)
    end
end

Base.size(indices::IGEO7Indices) = size(indices.ids)
Base.IndexStyle(::Type{<:IGEO7Indices}) = Base.IndexLinear()
Base.getindex(indices::IGEO7Indices, i::Int) = IGEO7Index(indices.ids[i])

@inline function Base.in(index::IGEO7Index, indices::IGEO7Indices)
    isempty(indices) && return false
    get_resolution(index) == indices.resolution || return false
    if indices.contiguous
        ordinal = cell_to_index(index.id)
        return indices.first_ordinal <= ordinal <
               indices.first_ordinal + length(indices)
    end
    return !iszero(Helpers.sorted_index(indices.ids, index.id))
end

function Base.show(io::IO, indices::IGEO7Indices)
    print(io, "IGEO7Indices(")
    if isempty(indices)
        print(io, "resolution=", indices.resolution, ", 0 cells")
    else
        show(io, first(indices))
        print(io, ':')
        show(io, last(indices))
        print(io, ", ", length(indices), " cells")
    end
    return print(io, ')')
end

Base.show(io::IO, ::MIME"text/plain", indices::IGEO7Indices) = show(io, indices)
