# `H3Cell` uses H3's self-describing `UInt64` layout:
#
#   bit  63     reserved, 0
#   bits 59-62  mode (1 for a cell index)
#   bits 56-58  mode-dependent, 0 for cells
#   bits 52-55  RESOLUTION, 0:15               <- `level(c)` reads this
#   bits 45-51  base cell, 0:121
#   bits 0-44   fifteen 3-bit digits, digit `k` at bit `45 - 3k`
#
# Digits below the resolution are padded with 7. Raw order is `(base cell,
# digit path)` within a level and resolution-major across levels.

"""
    H3Cell(id::UInt64) <: AbstractCellIndex

H3's `UInt64` cell encoding, with resolution in bits 52–55. [`level`](@ref)
reads that field and [`rawid`](@ref) returns the libh3-compatible value.
Unsigned comparison is `(base cell, digit path)` order within a level and
resolution-major across levels; deleted pentagon paths are absent.
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

Whether libh3 recognises `c` as a valid cell. This checks mode, base cell,
active and padding digits, and deleted pentagon branches. [`cellposition`](@ref)
and `subtree_border` check it first; [`children`](@ref) and
[`descendants`](@ref) deliberately do not, and libh3 will happily descend a
malformed index into malformed ones.
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
