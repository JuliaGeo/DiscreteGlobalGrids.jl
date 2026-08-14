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

"""
    trytranslate(cell, displacement) -> Union{Z7Cell,Nothing}

Apply a [`RelativeZ7Cell`](@ref) to its origin. Returns `nothing` when `cell`
is not that origin or the offset leaves the origin's level.
"""
function trytranslate(c::Z7Cell, d::RelativeZ7Cell)
    c == d.cell || return nothing
    position = DGG.cellposition(IGeo7System(), c)
    position === nothing && return nothing
    target_position = Int128(position) + Int128(d.offset)
    1 <= target_position <= num_cells(level(c)) || return nothing
    return Z7Cell(index_to_cell(Int(target_position), level(c)))
end

function Base.:+(c::Z7Cell, d::RelativeZ7Cell)
    target = trytranslate(c, d)
    target === nothing && throw(DomainError(
        (c, d),
        "a RelativeZ7Cell can only be applied to its origin without leaving its level",
    ))
    return target
end

Base.:+(d::RelativeZ7Cell, c::Z7Cell) = c + d

function Base.:-(d::RelativeZ7Cell)
    target = d.cell + d
    return RelativeZ7Cell(target, Base.Checked.checked_neg(d.offset))
end

Base.:-(c::Z7Cell, d::RelativeZ7Cell) = c + (-d)

"""
    directioncode(displacement) -> UInt8

Code `1:6` for an immediate neighbor, following [`neighbors`](@ref)'s
counterclockwise order, and `0` for the zero displacement.
"""
function directioncode(d::RelativeZ7Cell)
    iszero(d) && return UInt8(0)
    target = d.cell + d
    code = findfirst(
        ==(target),
        DGG.neighbors(DGG.levelgrid(IGeo7System(), level(d.cell)), d.cell),
    )
    code === nothing &&
        throw(ArgumentError("displacement does not point to an immediate neighbor"))
    return UInt8(code)
end
