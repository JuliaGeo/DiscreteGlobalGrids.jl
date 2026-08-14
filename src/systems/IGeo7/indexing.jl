"""
    RelativeZ7Cell(cell, offset)

A displacement from `cell` by `offset` positions in the canonical Z7
space-filling-curve order. It is produced by subtracting two [`Z7Cell`](@ref)s:
`target - origin == RelativeZ7Cell(origin, offset)`.

The origin is part of the displacement because a space-filling-curve offset is
not translation invariant. Consequently, `cell + displacement` requires
`cell == displacement.cell`.
"""
struct RelativeZ7Cell
    cell::Z7Cell
    offset::Int
end

Base.zero(d::RelativeZ7Cell) = RelativeZ7Cell(d.cell, 0)
Base.iszero(d::RelativeZ7Cell) = iszero(d.offset)

function Base.:-(target::Z7Cell, origin::Z7Cell)
    l = level(target)
    l == level(origin) ||
        throw(DimensionMismatch("IGeo7 cells have different levels"))
    target_position = DGG.cellposition(IGeo7System(), target)
    target_position === nothing &&
        throw(ArgumentError("$target is not a valid IGeo7 cell"))
    origin_position = DGG.cellposition(IGeo7System(), origin)
    origin_position === nothing &&
        throw(ArgumentError("$origin is not a valid IGeo7 cell"))
    return RelativeZ7Cell(
        origin,
        Base.Checked.checked_sub(target_position, origin_position),
    )
end

function Base.:+(c::Z7Cell, d::RelativeZ7Cell)
    c == d.cell || throw(DomainError(
        (c, d),
        "a RelativeZ7Cell can only be applied to its origin",
    ))
    position = DGG.cellposition(IGeo7System(), c)
    position === nothing && throw(ArgumentError("$c is not a valid IGeo7 cell"))
    target_position = Int128(position) + Int128(d.offset)
    1 <= target_position <= num_cells(level(c)) || throw(BoundsError(
        DGG.levelgrid(IGeo7System(), level(c)),
        target_position,
    ))
    return Z7Cell(index_to_cell(Int(target_position), level(c)))
end

Base.:+(d::RelativeZ7Cell, c::Z7Cell) = c + d

function Base.:-(d::RelativeZ7Cell)
    target = d.cell + d
    return RelativeZ7Cell(target, Base.Checked.checked_neg(d.offset))
end

Base.:-(c::Z7Cell, d::RelativeZ7Cell) = c + (-d)

const UNIT_CODES = (
    0x05, 0x04, 0xff, 0xff,
    0x06, 0x00, 0x03, 0xff,
    0xff, 0x01, 0x02, 0xff,
    0xff, 0xff, 0xff, 0xff,
)

@inline function _unitcode(a::Int64, b::Int64)
    (-1 <= a <= 1 && -1 <= b <= 1) || return 0xff
    key = ((Int(a) + 1) << 2) | (Int(b) + 1)
    return @inbounds UNIT_CODES[key+1]
end

"""
    directioncode(displacement) -> UInt8

Code `1:6` for an immediate neighbor, following [`neighbors`](@ref)'s
counterclockwise order, and `0` for the zero displacement.
"""
function directioncode(d::RelativeZ7Cell)
    iszero(d) && return UInt8(0)
    target = d.cell + d
    res = level(d.cell)
    base = z7_base_cell(d.cell.id)
    if z7_base_cell(target.id) == base
        a, b = _encode_lattice(d.cell.id, base, res)
        ta, tb = _encode_lattice(target.id, base, res)
        code = _unitcode(ta - a, tb - b)
        code != 0xff && return code

        # The two sides of a base's development-cone cut differ by one unit
        # rotation. Trying both exact rotations avoids geometry at the cut.
        for rotation in (-1, 1)
            ra, rb = unitmul(ta, tb, rotation)
            code = _unitcode(ra - a, rb - b)
            code != 0xff && return code
        end
    end

    # Crossing between base-cell development frames is rare and needs the
    # icosahedral topology decoder. Keep it as the complete fallback rather
    # than duplicating a second face-transition table here.
    for (code, neighbor) in enumerate(_cell_neighbors_ccw(d.cell.id))
        neighbor == target.id && return UInt8(code)
    end
    throw(ArgumentError("displacement does not point to an immediate neighbor"))
end
