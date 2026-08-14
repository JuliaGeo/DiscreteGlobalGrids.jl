# `A5Cell` uses upstream a5's self-describing `UInt64` serialization:
#
#   bits 58-63  QUINTANT, `5·origin + segment`, 0:59 (only `origin`, 0:11, at
#               resolution 0)
#   bits below  the Hilbert state `S`, two bits per level from resolution 2
#   one bit     the RESOLUTION MARKER: the highest set bit below `S`, at bit
#               `58 - R` with `R = level + 1` below resolution 2 and
#               `R = 2·level - 1` from resolution 2 down
#   the rest    zero
#
# `get_resolution` finds the marker bit. Above it, `(quintant, S)` gives the
# canonical lexicographic order within a level. Two encodings lie outside the
# system's levels:
#
#   * `WORLD_CELL == 0` is upstream's whole-sphere cell at "resolution -1". It
#     is representable as an `A5Cell` and `level` reports `-1` for it, but it is
#     not a cell of any grid.
#   * resolution 30 exists, with the quintant pushed into bits 59-63 by a
#     shorter marker, but only for the 42 of 60 quintants that still fit 64
#     bits. There is therefore no complete res-30 grid, `levels` stops at 29,
#     and `cellposition` answers `nothing` for a res-30 id at every level.

# The decoded record `A5Native` works in — origin, segment, Hilbert state,
# resolution — under a name that cannot be confused with the id type below.
const NativeCell = A5Native.A5Cell

# The deepest level with a COMPLETE grid, which is a fact about the encoding
# rather than about the system: `A5Native.MAX_RESOLUTION` is 30, but only 42 of
# the 60 quintants fit the res-30 layout. `levels(A5System())` is `0:MAX_LEVEL`.
const MAX_LEVEL = A5Native.MAX_GRID_RESOLUTION      # 29

"""
    A5Cell(id::UInt64) <: AbstractCellIndex

Upstream a5's `UInt64` cell encoding. [`level`](@ref) reads its marker bit and
[`rawid`](@ref) returns the interoperable value. Within a level, raw order is
`(quintant, Hilbert state)` and matches grid positions; across levels it is not
level-major. Encoded levels include the world cell at `-1` and partial level
30, while [`isvalid`](@ref) accepts only complete system levels.
"""
struct A5Cell <: AbstractCellIndex
    id::UInt64
end

A5Cell(id::Integer) = A5Cell(UInt64(id))
A5Cell(c::A5Cell) = c

# Total on every `A5Cell`, valid or not — including `-1` for the world cell and
# `30` for the deep encoding.
level(c::A5Cell) = A5Native.get_resolution(c.id)

rawid(c::A5Cell) = c.id

Base.isless(a::A5Cell, b::A5Cell) = isless(a.id, b.id)

Base.show(io::IO, c::A5Cell) =
    print(io, "A5Cell(0x", string(c.id; base=16, pad=16), ", res ", level(c), ")")

# ---------------------------------------------------------------------------
# The one decoder every entry point that can be handed an arbitrary id goes
# through.
# ---------------------------------------------------------------------------

"""
    _decode(id::UInt64) -> Union{NativeCell,Nothing}

Return the decoded record, or `nothing` unless all encoding invariants hold:

 1. resolution is in `0:29`;
 2. the quintant is in `0:11` at level 0 or `0:59` otherwise — checked before
    `deserialize`, which would otherwise index `ORIGINS` out of bounds;
 3. the id round-trips through `deserialize`/`serialize`, rejecting padding bits.
"""
function _decode(id::UInt64)
    r = A5Native.get_resolution(id)
    0 <= r <= MAX_LEVEL || return nothing
    top = id >> A5Native.HILBERT_START_BIT
    (r == 0 ? top <= 0x0b : top <= 0x3b) || return nothing
    cell = A5Native.deserialize(id)
    A5Native.serialize(cell) == id || return nothing
    return cell
end

"""
    isvalid(c::A5Cell) -> Bool

Whether `c` is a well-formed cell on a complete [`A5System`](@ref) level.
This rejects the world cell, resolution-30 encodings, invalid quintants, and
nonzero padding bits.
"""
Base.isvalid(c::A5Cell) = _decode(c.id) !== nothing
