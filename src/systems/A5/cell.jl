# ---------------------------------------------------------------------------
# `A5Cell` — the canonical id
#
# One `UInt64`, in upstream a5's own serialization, which is already
# self-describing about its level. The layout, from the top:
#
#   bits 58-63  QUINTANT, `5·origin + segment`, 0:59 (only `origin`, 0:11, at
#               resolution 0)
#   bits below  the Hilbert state `S`, two bits per level from resolution 2
#   one bit     the RESOLUTION MARKER: the highest set bit below `S`, at bit
#               `58 - R` with `R = level + 1` below resolution 2 and
#               `R = 2·level - 1` from resolution 2 down
#   the rest    zero
#
# So the level is read by *finding the marker bit* rather than by masking a
# field — `A5Native.get_resolution` walks down from bit 58 two bits at a time —
# and everything above the marker is `(quintant, S)` in the high bits. That is
# what makes unsigned comparison of the raw integer, at a fixed level, exactly
# `(quintant, S)` lexicographic, which is the canonical dense order.
#
# THE TWO EDGES OF THE ENCODING, both outside `levels(A5System())`:
#
#   * `WORLD_CELL == 0` is upstream's whole-sphere cell at "resolution -1". It
#     is representable as an `A5Cell` and `level` reports `-1` for it, but it is
#     not a cell of any grid.
#   * resolution 30 exists, with the quintant pushed into bits 59-63 by a
#     shorter marker, but only for the 42 of 60 quintants that still fit 64
#     bits. There is therefore no complete res-30 grid, `levels` stops at 29,
#     and `cellposition` answers `nothing` for a res-30 id at every level.
# ---------------------------------------------------------------------------

# The decoded record `A5Native` works in — origin, segment, Hilbert state,
# resolution — under a name that cannot be confused with the id type below.
const NativeCell = A5Native.A5Cell

# The deepest level with a COMPLETE grid, which is a fact about the encoding
# rather than about the system: `A5Native.MAX_RESOLUTION` is 30, but only 42 of
# the 60 quintants fit the res-30 layout. `levels(A5System())` is `0:MAX_LEVEL`.
const MAX_LEVEL = A5Native.MAX_GRID_RESOLUTION      # 29

"""
    A5Cell(id::UInt64) <: AbstractCellIndex

The canonical cell id of [`A5System`](@ref): one `UInt64` in upstream a5's own
index encoding, carrying its resolution in-band as the position of its
lowest-but-one marker bit.

`level(c)` finds that marker, so it needs no system and no table, and
[`rawid`](@ref) hands back the `UInt64` to print in hex or pass to another a5
implementation.

# Order

`isless` is unsigned comparison of the raw index, which **at one level** is
`(quintant, Hilbert state)` lexicographic — the same order the level grid's
positions run in, so a level grid's cells come out of [`cellindex`](@ref)
sorted, as the base interface requires.

Across levels it is *not* level-major: the quintant sits in the top bits and the
resolution marker below the Hilbert state, so a coarse cell and a fine one in
the same quintant interleave. Nothing in the package needs level-major raw
order — the one place that wants levels grouped, `MultiOrderCellSet`, sorts by
`(level(c), c)` explicitly — but it is worth knowing before comparing ids from
two levels and expecting the coarser one first.

# Levels outside the system

`level` is a **read of the encoding**, not an assertion about it: it answers
`-1` for the world cell (`A5Cell(0)`) and `30` for the deep encoding neither of
which [`A5System`](@ref) offers a grid at. [`isvalid`](@ref) is the predicate
that says whether an id names a cell of a real level.

```jldoctest
julia> c = DiscreteGlobalGrids.A5.A5Cell(0x8708000000000000)
A5Cell(0x8708000000000000, res 4)

julia> DiscreteGlobalGrids.level(c)
4
```
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

The decoded record behind `id`, or `nothing` when `id` names no cell of any
level [`A5System`](@ref) offers.

Three things have to hold, and each rules out a different way an arbitrary
`UInt64` can look plausible:

 1. the marker bit puts it at a **full-grid resolution**, `0:29` — which rejects
    the world cell (`-1`) and the res-30 encoding (see the header of this file);
 2. the top bits name a **real quintant**: `0:11` at resolution 0, `0:59` below
    it. Without this check `A5Native.deserialize` indexes `ORIGINS` out of
    bounds and raises where the caller was promised an answer;
 3. the id **round-trips** through `deserialize`/`serialize`. `deserialize`
    discards every bit below the Hilbert state, so an id carrying junk there
    would otherwise decode to a perfectly good cell that it is not equal to,
    and would be given that cell's position.
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

Whether `c` names a real cell of a real level of [`A5System`](@ref).

Worth spelling out because the a5 arithmetic validates almost nothing on its
own: `cell_to_children` of a malformed index cheerfully returns a subtree of
malformed indices, and `deserialize` raises a `BoundsError` from inside
`ORIGINS` for a top-bit pattern that names no dodecahedron face.

Malformed does not mean "random bits". The world cell `A5Cell(0)`, a res-30 id,
an id with a quintant of 60 or more, and an id carrying junk in the padding
below its Hilbert state all look plausible and are not cells of any grid this
system offers. [`cellposition`](@ref) checks this and answers `nothing`;
[`neighbors`](@ref) checks it and throws, because there is no `nothing` in its
contract to return.
"""
Base.isvalid(c::A5Cell) = _decode(c.id) !== nothing
