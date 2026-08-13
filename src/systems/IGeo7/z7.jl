# Z7 index layer: the UInt64 cell id, its string/hex codecs, prefix
# (hierarchy) operations and validation. Pure integer code — no geometry, no
# floating point, no allocation on any scalar path.
#
# Provenance
#   [contract] spec/interface-contract.md ("Z7 index" section), distilled from
#              the untainted ../IGeo7/README.md and ../IGeo7/test/runtests.jl.
#   [a7]       spec/aperture7-indexing-spec.md §4 (bit format, prefix ops,
#              subtree counting).
#   [design]   spec/design.md §3 (the Z7 UInt64 *is* the cell id) and §4.3.
#   [fitted]   test/IGeo7/vectors/{pentagon_chains,hierarchy}.csv.
#
# `SmallList` comes from the package-wide `Helpers` module, bound by
# `IGeo7.jl`'s `import ..Helpers`.

"""
    MAX_RESOLUTION

Finest resolution with geometry (`0:19`) **[contract]**. The Z7 bit format
holds one more level: prefix/string/hex operations are valid through
`Z7_MAX_RESOLUTION == 20`, which geometry constructors must reject.
"""
const MAX_RESOLUTION = 19

"""
    Z7_MAX_RESOLUTION

Number of 3-bit digit slots in the Z7 `UInt64`, i.e. the coarsest-to-finest
resolution range the bit format can represent (`0:20`) **[a7 §4.1]**.
"""
const Z7_MAX_RESOLUTION = 20

# Bit layout [contract, a7 §4.1]:
#   bits [63:60]                base cell 0:11
#   digit k (k = 1..20)         bits [62-3k : 60-3k]  (digit 1 = bits 59:57)
#   digit values 0:6 active, 7 = padding for every slot past the resolution.
#   Verification: base 8, digits 0,0,4,3,3 -> 0x80237fffffffffff.
const Z7_BASE_SHIFT = 60
# All twenty digit slots set to 7 (the resolution-0 digit pattern).
const Z7_PAD_MASK = UInt64(0x0fffffffffffffff)
# One bit per digit slot, at the slot's low bit (positions 0, 3, ..., 57):
# (2^60 - 1) / 7. Used to fold each 3-bit slot down to a single flag bit.
const Z7_SLOT_LSB = Z7_PAD_MASK ÷ UInt64(7)
# = `ISEA.NBASE`; restated here because this is the pure integer layer and must
# not depend on the shared geometry module. Pinned equal by
# test/IGeo7/test_icosahedron.jl "cross-file constant consistency".
const Z7_NUM_BASES = 12

"""
    Z7_DELETED_DIGIT

Per-base deleted (missing) child digit of the pentagon chain: digit `2` for
bases `0:5`, digit `5` for bases `6:11`; indexed by `base + 1`.

The deletion applies at every level for which the digit prefix is still all
zero — i.e. while the cell is still a pentagon — and only there.

Provenance: `test/IGeo7/vectors/pentagon_chains.csv` (all 12 bases, chain depths
1..6, `missing_digits` column) **[fitted]**, matching the hemisphere rule
predicted in `spec/z7-paper-spec.md` §5.4 and the `"002"` / `"065"` rejections
recorded in `spec/interface-contract.md` **[contract]**.

This is the package's single definition of the table; every other consumer
(`grid.jl`'s pentagon collapse and cone wrap, the codecs below) reaches it
through [`z7_deleted_digit`](@ref). See `PROVENANCE.md`.
"""
const Z7_DELETED_DIGIT = (2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5)

# --- error type ------------------------------------------------------------
#
# Every validation rejection in this file and in the id-validation paths of
# `grid.jl` throws one `InvalidZ7Error`, written out at the site as a plain
# `throw(InvalidZ7Error(:reason, ...))`. The struct carries only the *facts*
# (a reason tag plus the offending values); the human-readable sentence is
# built in `Base.showerror`, i.e. when the error is printed, never when it is
# thrown. So a failure edge stores four immediates (plus, for the codecs, the
# caller's own string) and calls `throw` — no `print_to_string` / `String`
# formatting machinery lands in `z7_resolution`, `z7_from_string`,
# `_geometry_checked`, ... Laziness comes from the *struct*, not from where
# the throw is written.
#
# There are deliberately no `_throw_*` helper functions: a function whose only
# job is to throw is indirection for its own sake, and the direct form keeps
# each site's reason and payload readable next to the predicate that rejects.
#
# `BoundsError` (dense index out of range) and the plain `ArgumentError`s of
# the lookups/trees wiring layers are deliberately left alone — they are cold
# API-shape contracts, not hot validators.

"""
    InvalidZ7Error(reason, value, got, limit[, input]) <: Exception

A Z7 index, digit, level, resolution or codec string was rejected. The
message is constructed lazily by [`Base.showerror`](@ref) from:

| field    | meaning |
|:---------|:--------|
| `reason` | tag from the taxonomy below; selects the message |
| `value`  | the offending Z7 index (`0` when the reason is not about an id) |
| `got`    | the offending scalar — digit, level, resolution, base cell, string position (truncated to `Int`) |
| `limit`  | the bound that made it invalid — a maximum, the id's own resolution, the base cell, or an offending character's code point |
| `input`  | the offending text for the string/hex codecs, `""` otherwise |

Reasons: `:base_range`, `:bad_padding`, `:invalid_index`, `:deleted_digit`,
`:digit_level`, `:child_digit`, `:no_children`, `:no_parent`, `:parent_res`,
`:bad_length`, `:bad_base_char`, `:bad_digit_char`, `:bad_hex`,
`:res20_geometry`, `:resolution_range`, `:descendant_res`,
`:no_child_geometry`.

The type is not `isbits` because `input` is a `String`; that costs nothing on
the non-throwing path (the field is set from a string the caller already
holds, and only on the failure edge that is about to throw).
"""
struct InvalidZ7Error <: Exception
    reason::Symbol
    value::UInt64
    got::Int
    limit::Int
    input::String
end

InvalidZ7Error(reason::Symbol, value::UInt64, got::Int, limit::Int) =
    InvalidZ7Error(reason, value, got, limit, "")

# Offending scalars reach the struct as `Int` without ever throwing a
# conversion error of their own (which would mask the validation error).
@inline _z7_int(x::Integer) = x % Int
@inline _z7_int(x::Char) = Int(x)

_z7_hex(z::UInt64) = "0x" * string(z, base=16, pad=16)

function Base.showerror(io::IO, e::InvalidZ7Error)
    r = e.reason
    if r === :base_range
        isempty(e.input) || print(io, "invalid Z7 string \"", e.input, "\": ")
        print(io, "Z7 base cell must be in 0:", Z7_NUM_BASES - 1, ", got ", e.got)
    elseif r === :bad_padding
        print(io, "malformed Z7 index ", _z7_hex(e.value),
            ": active digit after the padding sentinel 7")
    elseif r === :invalid_index
        print(io, "invalid Z7 index ", _z7_hex(e.value))
        isempty(e.input) || print(io, " (parsed from \"", e.input, "\")")
        print(io, ": the base cell must be in 0:", Z7_NUM_BASES - 1,
            ", the padding tail must be all-7 to slot ", Z7_MAX_RESOLUTION,
            ", and the first nonzero digit must not be the base's deleted ",
            "pentagon-chain digit")
    elseif r === :deleted_digit
        isempty(e.input) || print(io, "invalid Z7 string \"", e.input, "\": ")
        print(io, "digit ", e.got,
            " is deleted below the pentagon chain of base cell ", e.limit)
        isempty(e.input) && print(io, " (index ", _z7_hex(e.value), ")")
    elseif r === :digit_level
        print(io, "Z7 digit level must be in 1:", Z7_MAX_RESOLUTION, ", got ", e.got)
    elseif r === :child_digit
        print(io, "Z7 child digit must be in 0:6 (7 is the padding sentinel), got ", e.got)
    elseif r === :no_children
        print(io, "Z7 index ", _z7_hex(e.value), " is at resolution ", Z7_MAX_RESOLUTION,
            ", the finest the Z7 bit format holds, so it has no children")
    elseif r === :no_parent
        print(io, "Z7 index ", _z7_hex(e.value),
            " is a resolution-0 base cell and has no parent")
    elseif r === :parent_res
        print(io, "Z7 ancestor resolution must be in 0:", e.limit, ", got ", e.got,
            " (index ", _z7_hex(e.value), ")")
    elseif r === :bad_length
        print(io, "invalid Z7 string \"", e.input, "\": expected 2 base characters and 0:",
            Z7_MAX_RESOLUTION, " digits, got ", e.got,
            e.got == 1 ? " character" : " characters")
    elseif r === :bad_base_char
        print(io, "invalid Z7 string \"", e.input, "\": base character ", e.got, " is '",
            Char(e.limit), "'; the base cell must be two decimal digits")
    elseif r === :bad_digit_char
        print(io, "invalid Z7 string \"", e.input, "\": digit ", e.got, " is '",
            Char(e.limit), "'; Z7 string digits must be in 0:6 (7 is the padding sentinel)")
    elseif r === :bad_hex
        print(io, "invalid Z7 hex \"", e.input,
            "\": expected 16 hexadecimal characters, optionally prefixed with \"0x\"")
    elseif r === :res20_geometry
        print(io, "Z7 index ", _z7_hex(e.value), " is at resolution ", Z7_MAX_RESOLUTION,
            ", which has no geometry (MAX_RESOLUTION = ", e.limit,
            "); it is valid for prefix arithmetic only")
    elseif r === :resolution_range
        print(io, "resolution must be in 0:", e.limit, ", got ", e.got)
    elseif r === :descendant_res
        print(io, "descendant resolution must be in ", e.limit, ":", MAX_RESOLUTION,
            ", got ", e.got, " (cell ", _z7_hex(e.value), ")")
    elseif r === :no_child_geometry
        print(io, "cell ", _z7_hex(e.value), " is at resolution ", e.limit,
            " (MAX_RESOLUTION), so its children would be resolution ",
            Z7_MAX_RESOLUTION, " and have no geometry")
    else
        print(io, "invalid Z7 value ", _z7_hex(e.value), " (", r, ", got ", e.got,
            ", limit ", e.limit, ")")
        isempty(e.input) || print(io, " from \"", e.input, "\"")
    end
    return nothing
end

"""
    z7_deleted_digit(base) -> Int

Child digit that does not exist below the pentagon chain of `base`
(see `Z7_DELETED_DIGIT`). Throws [`InvalidZ7Error`](@ref) (`:base_range`) for
`base ∉ 0:11`.
"""
@inline function z7_deleted_digit(base::Integer)
    0 <= base < Z7_NUM_BASES ||
        throw(InvalidZ7Error(:base_range, zero(UInt64), _z7_int(base), Z7_NUM_BASES - 1))
    return @inbounds Z7_DELETED_DIGIT[Int(base)+1]
end

# --- raw slot arithmetic ---------------------------------------------------

# Low bit position of digit slot `k` (k = 1..20); also the width of the tail
# below it.
@inline _z7_shift(k::Integer) = Z7_BASE_SHIFT - 3 * Int(k)

# Mask of every slot strictly finer than resolution `res` (all-7 there for a
# well-formed id). `res == 20` gives 0.
@inline _z7_tail_mask(res::Integer) = (UInt64(1) << (Z7_BASE_SHIFT - 3 * Int(res))) - UInt64(1)

@inline _z7_digit(z::UInt64, k::Integer) = Int((z >> _z7_shift(k)) & UInt64(7))

@inline _z7_set_digit(z::UInt64, k::Integer, d::Integer) =
    (z & ~(UInt64(7) << _z7_shift(k))) | (UInt64(d) << _z7_shift(k))

# Number of leading non-7 digit slots, without checking that the padding tail
# is well formed. O(1): XOR against the all-7 pattern turns a padded slot into
# a zero slot, `~(y | y>>1 | y>>2)` folds each zero slot to a 1 bit at the
# slot's low bit, and `leading_zeros` locates the first one [design §3].
@inline function _z7_leading_resolution(z::UInt64)
    y = z ⊻ Z7_PAD_MASK
    m = ~(y | (y >> 1) | (y >> 2)) & Z7_SLOT_LSB
    # Slot k folds to bit 60-3k, so leading_zeros == 6 + 3*(k-1) for the first
    # padded slot k; no padded slot at all means resolution 20.
    return m == zero(UInt64) ? Z7_MAX_RESOLUTION : (leading_zeros(m) - 6) ÷ 3
end

# --- accessors -------------------------------------------------------------

"""
    z7_base_cell(z7) -> Int

Base cell (`0:11`) held in bits `[63:60]`.
"""
@inline z7_base_cell(z7::UInt64) = Int(z7 >> Z7_BASE_SHIFT)
@inline z7_base_cell(z7::Unsigned) = z7_base_cell(UInt64(z7))

"""
    z7_resolution(z7) -> Int

Number of leading active (non-`7`) digit slots, `0:20`.

Throws [`InvalidZ7Error`](@ref) (`:bad_padding`) when the padding tail is
malformed (a non-`7` slot after the first `7` slot) **[contract]**; that is the
only validity check performed here, so the happy path stays branch-predictable
and allocation-free.
"""
@inline function z7_resolution(z7::UInt64)
    res = _z7_leading_resolution(z7)
    tail = _z7_tail_mask(res)
    (z7 & tail) == tail || throw(InvalidZ7Error(:bad_padding, z7, 0, 0))
    return res
end
@inline z7_resolution(z7::Unsigned) = z7_resolution(UInt64(z7))

"""
    z7_digit(z7, k) -> Int

Digit at 1-based level `k`, i.e. the child digit taken at resolution `k`.
Returns `7` for every level past the index's resolution **[contract]**.
Throws [`InvalidZ7Error`](@ref) (`:digit_level`) for `k ∉ 1:20`.
"""
@inline function z7_digit(z7::UInt64, k::Integer)
    1 <= k <= Z7_MAX_RESOLUTION ||
        throw(InvalidZ7Error(:digit_level, zero(UInt64), _z7_int(k), Z7_MAX_RESOLUTION))
    return _z7_digit(z7, k)
end
@inline z7_digit(z7::Unsigned, k::Integer) = z7_digit(UInt64(z7), k)

"""
    z7_is_pentagon(z7) -> Bool

`true` when every active digit is `0` — the twelve base cells and their center
descendants, one pentagon per icosahedron vertex per resolution **[contract]**.
O(1): mask off the padding tail and compare the active region with zero.
"""
@inline function z7_is_pentagon(z7::UInt64)
    res = z7_resolution(z7)
    active = Z7_PAD_MASK ⊻ _z7_tail_mask(res)
    return (z7 & active) == zero(UInt64)
end
@inline z7_is_pentagon(z7::Unsigned) = z7_is_pentagon(UInt64(z7))

# --- validation ------------------------------------------------------------

"""
    is_valid_z7(z7) -> Bool

Full structural validation **[contract]**:

1. base cell `< 12`;
2. no active digit after the first `7` (padding is solid to slot 20);
3. the first nonzero digit is not the base's deleted digit — pentagon-chain
   subsequences such as `"002"` and `"065"` do not exist.

Active digits are necessarily `0:6` because the resolution is defined as the
run of leading non-`7` slots.
"""
function is_valid_z7(z7::UInt64)
    base = z7_base_cell(z7)
    base < Z7_NUM_BASES || return false
    res = _z7_leading_resolution(z7)
    tail = _z7_tail_mask(res)
    (z7 & tail) == tail || return false
    deleted = @inbounds Z7_DELETED_DIGIT[base+1]
    for k in 1:res
        d = _z7_digit(z7, k)
        d == 0 && continue
        return d != deleted   # first nonzero digit decides the chain
    end
    return true
end
is_valid_z7(z7::Unsigned) = is_valid_z7(UInt64(z7))

# --- prefix operations -----------------------------------------------------

"""
    z7_parent(z7) -> UInt64
    z7_parent(z7, res) -> UInt64

Ancestor at resolution `res` (default: one coarser), formed by filling every
finer digit slot with the padding sentinel `7` **[a7 §4.2]**.

Throws [`InvalidZ7Error`](@ref) when `res` is outside `0:z7_resolution(z7)`
(`:parent_res`), in particular when a resolution-0 base cell is asked for its
parent (`:no_parent`).
"""
@inline function z7_parent(z7::UInt64, res::Integer)
    current = z7_resolution(z7)
    0 <= res <= current ||
        throw(InvalidZ7Error(:parent_res, z7, _z7_int(res), current))
    return z7 | _z7_tail_mask(res)
end

@inline function z7_parent(z7::UInt64)
    current = z7_resolution(z7)
    current > 0 || throw(InvalidZ7Error(:no_parent, z7, 0, 0))
    return z7 | _z7_tail_mask(current - 1)
end

@inline z7_parent(z7::Unsigned) = z7_parent(UInt64(z7))
@inline z7_parent(z7::Unsigned, res::Integer) = z7_parent(UInt64(z7), res)

"""
    z7_child(z7, digit) -> UInt64

Append `digit` (`0:6`) to the digit string.

Throws [`InvalidZ7Error`](@ref) when `digit` is outside `0:6` (`:child_digit`;
`7` is the padding sentinel, never a child), when `z7` is a pentagon and
`digit` is that base's deleted digit (`:deleted_digit`) **[contract]**, or when
`z7` is already at resolution $(Z7_MAX_RESOLUTION) and has no room for another
digit (`:no_children`).
"""
@inline function z7_child(z7::UInt64, digit::Integer)
    res = z7_resolution(z7)
    res < Z7_MAX_RESOLUTION ||
        throw(InvalidZ7Error(:no_children, z7, 0, Z7_MAX_RESOLUTION))
    0 <= digit <= 6 ||
        throw(InvalidZ7Error(:child_digit, zero(UInt64), _z7_int(digit), 6))
    base = z7_base_cell(z7)
    base < Z7_NUM_BASES ||
        throw(InvalidZ7Error(:base_range, zero(UInt64), base, Z7_NUM_BASES - 1))
    active = Z7_PAD_MASK ⊻ _z7_tail_mask(res)
    if (z7 & active) == zero(UInt64) && digit == @inbounds Z7_DELETED_DIGIT[base+1]
        throw(InvalidZ7Error(:deleted_digit, z7, _z7_int(digit), base))
    end
    return _z7_set_digit(z7, res + 1, digit)
end
@inline z7_child(z7::Unsigned, digit::Integer) = z7_child(UInt64(z7), digit)

"""
    z7_children(z7) -> SmallList{7,UInt64}

Immediate children in ascending id order: seven for a hexagon, six for a
pentagon (its deleted digit is skipped) **[contract, fitted:
test/IGeo7/vectors/pentagon_chains.csv]**. Ascending digit order is ascending
`UInt64` order because the digit slots are the high bits of the tail.

Throws [`InvalidZ7Error`](@ref) (`:no_children`) at resolution
$(Z7_MAX_RESOLUTION) **[contract]**.
"""
function z7_children(z7::UInt64)
    res = z7_resolution(z7)
    res < Z7_MAX_RESOLUTION ||
        throw(InvalidZ7Error(:no_children, z7, 0, Z7_MAX_RESOLUTION))
    base = z7_base_cell(z7)
    base < Z7_NUM_BASES ||
        throw(InvalidZ7Error(:base_range, zero(UInt64), base, Z7_NUM_BASES - 1))
    active = Z7_PAD_MASK ⊻ _z7_tail_mask(res)
    pentagon = (z7 & active) == zero(UInt64)
    deleted = @inbounds Z7_DELETED_DIGIT[base+1]
    shift = _z7_shift(res + 1)
    cleared = z7 & ~(UInt64(7) << shift)
    out = Helpers.empty_small_list(Val(7), zero(UInt64))
    for digit in 0:6
        pentagon && digit == deleted && continue
        out = Helpers.small_push(out, cleared | (UInt64(digit) << shift))
    end
    return out
end
z7_children(z7::Unsigned) = z7_children(UInt64(z7))

"""
    z7_is_descendant(z7, ancestor) -> Bool

`true` when `z7` lies in the subtree rooted at `ancestor`, i.e. when it is at
least as fine and shares the ancestor's digit prefix **[a7 §4.2]**. Reflexive:
an index is a descendant of itself.
"""
@inline function z7_is_descendant(z7::UInt64, ancestor::UInt64)
    res = z7_resolution(ancestor)
    z7_resolution(z7) >= res || return false
    tail = _z7_tail_mask(res)
    return (z7 | tail) == (ancestor | tail)
end
@inline z7_is_descendant(z7::Unsigned, ancestor::Unsigned) =
    z7_is_descendant(UInt64(z7), UInt64(ancestor))

# --- string codec ----------------------------------------------------------

"""
    z7_to_string(z7) -> String

The Z7 spec's `Z7_STRING` form (`spec/z7-paper-spec.md`): the two-character
zero-padded base cell followed by one character per active digit, e.g.
`"00"`, `"0800433"` **[contract]**.
"""
function z7_to_string(z7::UInt64)
    base = z7_base_cell(z7)
    base < Z7_NUM_BASES ||
        throw(InvalidZ7Error(:base_range, zero(UInt64), base, Z7_NUM_BASES - 1))
    res = z7_resolution(z7)
    buffer = Vector{UInt8}(undef, 2 + res)
    @inbounds buffer[1] = UInt8('0') + (base ÷ 10) % UInt8
    @inbounds buffer[2] = UInt8('0') + (base % 10) % UInt8
    for k in 1:res
        @inbounds buffer[2+k] = UInt8('0') + _z7_digit(z7, k) % UInt8
    end
    return String(buffer)
end
z7_to_string(z7::Unsigned) = z7_to_string(UInt64(z7))

"""
    z7_from_string(text) -> UInt64

Parse a `Z7_STRING`. Validates length (2 base characters plus `0:20` digits),
base cell `< 12`, digit range `0:6` (`7` is padding and never appears in the
string form), and the deleted pentagon subsequence — so `"0"`, `"12"`, `"002"`,
`"065"` and `"0107"` all throw [`InvalidZ7Error`](@ref) (`:bad_length`,
`:base_range`, `:deleted_digit`, `:bad_base_char`, `:bad_digit_char`)
**[contract]**. The loop below contains no string interpolation: each failure
edge hands the caller's own text to the error struct, which formats it only if
the error is printed.
"""
function z7_from_string(text::AbstractString)
    n = length(text)
    2 <= n <= 2 + Z7_MAX_RESOLUTION || throw(InvalidZ7Error(
        :bad_length, zero(UInt64), _z7_int(n), 2 + Z7_MAX_RESOLUTION, String(text)))
    z = Z7_PAD_MASK
    base = 0
    deleted = 0
    res = 0
    pentagon = true
    for (i, char) in enumerate(text)
        value = Int(char) - Int('0')
        if i <= 2
            0 <= value <= 9 || throw(InvalidZ7Error(
                :bad_base_char, zero(UInt64), i, _z7_int(char), String(text)))
            base = 10 * base + value
            if i == 2
                base < Z7_NUM_BASES || throw(InvalidZ7Error(
                    :base_range, zero(UInt64), base, Z7_NUM_BASES - 1, String(text)))
                deleted = @inbounds Z7_DELETED_DIGIT[base+1]
            end
        else
            0 <= value <= 6 || throw(InvalidZ7Error(
                :bad_digit_char, zero(UInt64), i - 2, _z7_int(char), String(text)))
            if pentagon
                value == deleted && throw(InvalidZ7Error(
                    :deleted_digit, zero(UInt64), value, base, String(text)))
                pentagon = value == 0
            end
            res += 1
            z = _z7_set_digit(z, res, value)
        end
    end
    return z | (UInt64(base) << Z7_BASE_SHIFT)
end

# --- hexadecimal codec -----------------------------------------------------

"""
    z7_to_hex(z7; prefix=false) -> String

Sixteen lowercase hexadecimal characters, optionally prefixed with `"0x"`
**[contract]**.
"""
function z7_to_hex(z7::UInt64; prefix::Bool=false)
    text = string(z7, base=16, pad=16)
    return prefix ? "0x" * text : text
end
z7_to_hex(z7::Unsigned; kwargs...) = z7_to_hex(UInt64(z7); kwargs...)

"""
    z7_from_hex(text) -> UInt64

Parse the 16-character hexadecimal form, with an optional `"0x"` / `"0X"`
prefix, and validate the result with `is_valid_z7` **[contract]**. Throws
[`InvalidZ7Error`](@ref) (`:bad_hex`, `:invalid_index`).
"""
function z7_from_hex(text::AbstractString)
    body = (startswith(text, "0x") || startswith(text, "0X")) ? SubString(text, 3) : text
    # `length` + `isxdigit` rather than trusting `tryparse`, which tolerates
    # surrounding whitespace and would silently accept a short index.
    (length(body) == 16 && all(isxdigit, body)) ||
        throw(InvalidZ7Error(:bad_hex, zero(UInt64), 0, 16, String(text)))
    z7 = tryparse(UInt64, body; base=16)
    z7 === nothing && throw(InvalidZ7Error(:bad_hex, zero(UInt64), 0, 16, String(text)))
    is_valid_z7(z7) ||
        throw(InvalidZ7Error(:invalid_index, z7, 0, 0, String(text)))
    return z7
end
