# ---------------------------------------------------------------------------
# `H3Cell` — the canonical id
#
# One `UInt64`, in Uber's own H3 index encoding, which is already
# self-describing about its level. The bit layout, from the top:
#
#   bit  63     reserved, 0
#   bits 59-62  mode (1 for a cell index)
#   bits 56-58  mode-dependent, 0 for cells
#   bits 52-55  RESOLUTION, 0:15               <- `level(c)` reads this
#   bits 45-51  base cell, 0:121
#   bits 0-44   fifteen 3-bit digits, digit `k` at bit `45 - 3k`
#
# Digits below the cell's own resolution are padded with 7, which is what makes
# the raw integer order come out right: two cells at the same resolution share
# their padding, so comparing the raw `UInt64`s compares
# `(base cell, digit path)` lexicographically — exactly the canonical order
# their positions in a level grid are laid out in. Across resolutions the
# resolution field outranks the base cell, so the raw order is
# resolution-major, which is the order `has_sorted_subtrees` is declared
# against.
# ---------------------------------------------------------------------------

"""
    H3Cell(id::UInt64) <: AbstractCellIndex

The canonical cell id of [`H3System`](@ref): one `UInt64` in H3's own index
encoding, carrying its resolution in bits 52-55.

`level(c)` reads those bits, so it needs no system and no table, and
[`rawid`](@ref) hands back the `UInt64` to print in hex or pass to libh3.

# Order

`isless` is unsigned comparison of the raw index, which is:

  - **within one level**, `(base cell, digit path)` lexicographic order — the
    same order the level grid's positions run in, so a level grid's cells come
    out of [`cellindex`](@ref) sorted, as the base interface requires;
  - **across levels**, resolution-major, because the resolution field sits
    above the base-cell field.

Pentagon cells have a deleted digit, so some digit paths name no cell; those
gaps are absent from both orders rather than reordering anything around them.

```jldoctest
julia> c = DiscreteGlobalGrids.H3.H3Cell(0x8928308280fffff)
H3Cell(0x08928308280fffff, res 9)

julia> DiscreteGlobalGrids.level(c)
9
```
"""
struct H3Cell <: AbstractCellIndex
    id::UInt64
end

H3Cell(id::Integer) = H3Cell(UInt64(id))
H3Cell(c::H3Cell) = c

# The resolution field, bits 52-55. Total on every `H3Cell`, valid or not:
# `level` is a read of the encoding, not an assertion about it.
level(c::H3Cell) = Int((c.id >> 52) & 0x0f)

rawid(c::H3Cell) = c.id

Base.isless(a::H3Cell, b::H3Cell) = isless(a.id, b.id)

Base.show(io::IO, c::H3Cell) =
    print(io, "H3Cell(0x", string(c.id; base=16, pad=16), ", res ", level(c), ")")

# ---------------------------------------------------------------------------
# The bit vocabulary the hierarchy arithmetic and the border automaton share.
# Ported from the old `src/H3/H3Kernel.jl`.
# ---------------------------------------------------------------------------

const MAX_RESOLUTION = H3Native.MAX_RESOLUTION      # 15

const _H3_RESOLUTION_MASK = UInt64(0x0f) << 52

# Bit offset of digit `k`, `k in 1:15`.
_h3_digit_shift(k::Int) = 45 - 3k

_h3_resolution(id::UInt64) = Int((id >> 52) & 0x0f)

# The resolution field rewritten, digits left alone. Used when walking a
# subtree: the id keeps its prefix and gains the target resolution once.
_h3_with_resolution(id::UInt64, res::Int) =
    (id & ~_H3_RESOLUTION_MASK) | (UInt64(res) << 52)

"""
    isvalid(c::H3Cell) -> Bool

Whether libh3 recognises `c` as a real cell index.

Worth spelling out because libh3 validates almost nothing on its own:
`cellToChildren` of a malformed index cheerfully returns a subtree of
malformed indices.

The entry points that can be handed an arbitrary id and would otherwise
enumerate or place cells that do not exist — [`cellposition`](@ref) and
`subtree_border` — check this first. [`children`](@ref) and
[`descendants`](@ref) deliberately do **not**: they sit in tree-descent inner
loops, where their caller already holds a cell it got from this system, and a
validity ccall per node is a real cost for a case that cannot arise there.
Garbage in, garbage out is the contract for those two.

Malformed does not mean "random bits": clearing a padding digit, writing a 7
into an active digit slot, naming base cell 122, or taking the deleted K-axis
child of a pentagon all produce indices that look plausible and are not cells.
"""
Base.isvalid(c::H3Cell) = H3Native.is_valid_cell(c.id)

"""
    ispentagon(c::H3Cell) -> Bool

Whether `c` is one of the twelve pentagons at its resolution — the cells with
five neighbours and six children instead of six and seven.
"""
ispentagon(c::H3Cell) = H3Native.is_pentagon(c.id)

"""
    basecell(c::H3Cell) -> Int

The `0:121` base cell `c` descends from: bits 45-51 of the index, and the
major key of the canonical order.
"""
basecell(c::H3Cell) = H3Native.get_base_cell(c.id)
