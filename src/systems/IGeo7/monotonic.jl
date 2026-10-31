# The base-7 structural line of Z7 ids. Pure integer code on top of `z7.jl`.
#
# A Z7 id read as a base-7 numeral, with the base cell as the leading digit,
# gives a position on `[0, 12·7^r)`: every digit string of resolution `r` has a
# slot, including the strings below a pentagon's deleted digit that name no
# cell. Every parent at resolution `r` owns exactly the block
# `[m·7^k, (m+1)·7^k)` at resolution `r + k`, and this block property is what
# the IGEO7 range stores written by z7py and xdggs-dggrid4py are built on.
#
# `cell_to_index` is the other line: the dense rank among the `10·7^r + 2`
# cells that exist. Both lines are strictly increasing in the packed `UInt64`,
# so they agree on order and differ on spacing, and only across a pentagon's
# deleted branch. A stored `[start, stop]` interval is a statement about
# spacing, which is why a reader has to know which line a store counts on.

"""
    z7_to_monotonic(z7) -> UInt64
    z7_to_monotonic(z7, res) -> UInt64

Position of `z7` on the base-7 structural line,
`base·7^res + Σₖ dₖ·7^(res−k)` over `[0, 12·7^res)`: the digit string read in
base 7 with the base cell as the leading digit. This is the value z7py's
`z7_to_monotonic_int` computes for the IGEO7 range stores.

It is not [`cell_to_index`](@ref), which ranks a cell among the cells that
exist and skips the pentagons' deleted digits. Both give the same order, but
this line keeps a slot for every digit string, so the two differ wherever an
interval crosses a deleted branch.

The two-argument form requires `res` to be the id's own resolution, so a
mixed-resolution id vector fails instead of being projected onto a line it does
not belong to. Throws [`InvalidZ7Error`](@ref) for a malformed padding tail
(`:bad_padding`), for `res ∉ 0:$(Z7_MAX_RESOLUTION)` (`:resolution_range`) and
for a resolution other than the id's (`:monotonic_res`).
"""
@inline function z7_to_monotonic(z7::UInt64, res::Integer)
    r = Int(res)
    0 <= r <= Z7_MAX_RESOLUTION || throw(InvalidZ7Error(
        :resolution_range, z7, _z7_int(res), Z7_MAX_RESOLUTION))
    own = z7_resolution(z7)
    own == r || throw(InvalidZ7Error(:monotonic_res, z7, own, r))
    value = UInt64(z7_base_cell(z7))
    for k in 1:r
        value = value * UInt64(7) + UInt64(_z7_digit(z7, k))
    end
    return value
end

@inline z7_to_monotonic(z7::UInt64) = z7_to_monotonic(z7, z7_resolution(z7))
@inline z7_to_monotonic(z7::Unsigned, res::Integer) = z7_to_monotonic(UInt64(z7), res)
@inline z7_to_monotonic(z7::Unsigned) = z7_to_monotonic(UInt64(z7))

"`MONOTONIC_LIMIT[res+1] = 12·7^res`, one past the last position at `res`."
const MONOTONIC_LIMIT = ntuple(i -> UInt64(Z7_NUM_BASES) * UInt64(7)^(i - 1),
    Z7_MAX_RESOLUTION + 1)

"""
    z7_from_monotonic(value, res) -> UInt64

Inverse of [`z7_to_monotonic`](@ref): the packed id at position `value` of the
base-7 structural line at resolution `res`, padded with the sentinel `7`.

The result is a well-formed digit string, but not necessarily a cell: `value`
can name a string below a pentagon's deleted digit, which
[`is_valid_z7`](@ref) rejects. The line has more positions than the grid has
cells, and a caller that needs cells validates. Throws
[`InvalidZ7Error`](@ref) for `res ∉ 0:$(Z7_MAX_RESOLUTION)`
(`:resolution_range`) and for `value ≥ 12·7^res` (`:monotonic_range`).
"""
function z7_from_monotonic(value::UInt64, res::Integer)
    r = Int(res)
    0 <= r <= Z7_MAX_RESOLUTION || throw(InvalidZ7Error(
        :resolution_range, zero(UInt64), _z7_int(res), Z7_MAX_RESOLUTION))
    limit = @inbounds MONOTONIC_LIMIT[r+1]
    value < limit || throw(InvalidZ7Error(
        :monotonic_range, zero(UInt64), _z7_int(value), _z7_int(limit - one(UInt64))))
    # Peel the digits from the lowest, writing each into its slot; what is left
    # after `r` divisions is the base cell.
    z = Z7_PAD_MASK
    remainder = value
    for k in r:-1:1
        z = _z7_set_digit(z, k, remainder % UInt64(7))
        remainder ÷= UInt64(7)
    end
    return (remainder << Z7_BASE_SHIFT) | z
end

z7_from_monotonic(value::Unsigned, res::Integer) = z7_from_monotonic(UInt64(value), res)
