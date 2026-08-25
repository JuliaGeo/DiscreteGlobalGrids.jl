"""
    RelativeZ7Cell(cell, offset)

A displacement from `cell` by `offset` indices in the canonical Z7
space-filling-curve order. It is produced by subtracting two [`Z7Cell`](@ref)s:
`target - origin == RelativeZ7Cell(origin, offset)`.

The origin is part of the displacement because a space-filling-curve offset is
not translation invariant. Consequently, `cell + displacement` requires
`cell == displacement.cell`.

Every rejection is a [`RelativeZ7Error`](@ref), except an id that names no cell
at all, which is Z7's own [`InvalidZ7Error`](@ref).
"""
struct RelativeZ7Cell
    cell::Z7Cell
    offset::Int
end

Base.zero(d::RelativeZ7Cell) = RelativeZ7Cell(d.cell, 0)
Base.iszero(d::RelativeZ7Cell) = iszero(d.offset)

# Like `InvalidZ7Error`, this stores structured facts and formats them in
# `showerror`, so the throw sites below carry no string machinery at all.

"""
    RelativeZ7Error(reason, cell, other, offset) <: Exception

A rejected [`RelativeZ7Cell`](@ref) operation. The fields record the reason,
the cell the operation started from, the second cell involved (the cell itself
when only one is), and the offset. [`Base.showerror`](@ref) constructs the
message lazily.

Reasons: `:level_mismatch` (subtracting cells from different levels),
`:foreign_origin` (applying a displacement to a cell that is not its origin),
`:out_of_range` (the offset leaves the level), `:not_a_neighbor`
([`directioncode`](@ref) of a displacement longer than one step).
"""
struct RelativeZ7Error <: Exception
    reason::Symbol
    cell::Z7Cell
    other::Z7Cell
    offset::Int
end

# Out of line so that no operator below carries a throw path in its own body:
# these are one-cache-line functions whose whole point is to be inlined into
# hot stencil loops.
@noinline _rel_error(reason::Symbol, cell::Z7Cell, other::Z7Cell, offset::Integer) =
    throw(RelativeZ7Error(reason, cell, other, offset % Int))

function Base.showerror(io::IO, e::RelativeZ7Error)
    r = e.reason
    if r === :level_mismatch
        print(io, "IGeo7 cells are at different levels: ", e.cell, " is at level ",
            level(e.cell), " and ", e.other, " at level ", level(e.other))
    elseif r === :foreign_origin
        print(io, "a RelativeZ7Cell displaces its own origin ", e.other,
            ", so it cannot be applied to ", e.cell,
            "; rebase it with `target - ", e.cell, "` first")
    elseif r === :out_of_range
        index = cell_to_index(e.cell.id)
        n = num_cells(level(e.cell))
        print(io, "offset ", e.offset, " leaves level ", level(e.cell), ": from ",
            e.cell, " (index ", index, " of ", n, ") the offsets in range are ",
            1 - index, ":", n - index)
    elseif r === :not_a_neighbor
        print(io, "offset ", e.offset, " from ", e.cell, " reaches ", e.other,
            ", which is not one of its immediate neighbors")
    else
        print(io, "invalid RelativeZ7Cell operation (", r, ") on ", e.cell,
            " with offset ", e.offset)
    end
    return nothing
end

# `globalindex`'s walk, but throwing Z7's own validation error instead of
# answering `nothing`: a displacement between ids that name no cell is a caller
# bug, exactly as it is on the geometry entry points.
@inline function _checked_index(c::Z7Cell)
    _geometry_checked(c.id)
    return cell_to_index(c.id)
end

function Base.:-(target::Z7Cell, origin::Z7Cell)
    level(target) == level(origin) || _rel_error(:level_mismatch, target, origin, 0)
    # Indices live in `1:num_cells(19)`, i.e. below 1.2e17, so the difference
    # is two decimal orders away from overflowing an `Int`. No widening needed.
    return RelativeZ7Cell(
        origin,
        _checked_index(target) - _checked_index(origin),
    )
end

function Base.:+(c::Z7Cell, d::RelativeZ7Cell)
    c == d.cell || _rel_error(:foreign_origin, c, d.cell, d.offset)
    index = _checked_index(c)
    n = num_cells(level(c))
    # Bound the OFFSET rather than the target index: `index ∈ 1:n` makes
    # both ends of this window exact in `Int`, so the check itself cannot
    # overflow and the addition it guards cannot either.
    (1 - index) <= d.offset <= (n - index) ||
        _rel_error(:out_of_range, c, c, d.offset)
    return Z7Cell(index_to_cell(index + d.offset, level(c)))
end

Base.:+(d::RelativeZ7Cell, c::Z7Cell) = c + d

function Base.:-(d::RelativeZ7Cell)
    # `+` has already confined `d.offset` to the level, so negating it is safe.
    target = d.cell + d
    return RelativeZ7Cell(target, -d.offset)
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
counterclockwise order, and `0` for the zero displacement. Throws
[`RelativeZ7Error`](@ref) (`:not_a_neighbor`) for anything longer.
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
    _rel_error(:not_a_neighbor, d.cell, target, d.offset)
end
