# A compact, prefix ordered aperture-3 identifier.  The external string uses
# DGGRID's two-decimal-digit root followed by ternary digits.  Internally the
# root occupies four bits and thirty base-3 digits occupy two bits each; 3 is
# padding.  Root 00 and root 11 are the two polar, all-zero chains.

const Z3_DIGITS = 30
const Z3_ROOT_SHIFT = 60
const Z3_PAD = UInt64(0x0fff_ffff_ffff_ffff)

"""
    Z3Cell(text)

Compact `Z3Cell` identifier for `ISEA3HSystem`. Text is the two-digit decimal
root followed by zero to thirty base-3 prefix digits, for example `"071201"`.
"""
struct Z3Cell <: DGG.AbstractCellIndex
    id::UInt64
end

Z3Cell(x::Unsigned) = Z3Cell(UInt64(x))

@inline _z3shift(k::Integer) = 60 - 2 * Int(k)
@inline _z3digit(z::UInt64, k::Integer) = Int((z >> _z3shift(k)) & 0x03)
@inline _z3root(z::UInt64) = Int(z >> Z3_ROOT_SHIFT)

function _z3level(z::UInt64)
    for k in 1:Z3_DIGITS
        _z3digit(z, k) == 3 && return k - 1
    end
    return Z3_DIGITS
end

DGG.level(c::Z3Cell) = _z3level(c.id)
DGG.rawid(c::Z3Cell) = c.id
Base.isless(a::Z3Cell, b::Z3Cell) = isless(a.id, b.id)

function _valid_z3(z::UInt64)
    root = _z3root(z)
    root <= 11 || return false
    level = _z3level(z)
    for k in 1:level
        _z3digit(z, k) <= 2 || return false
    end
    for k in (level + 1):Z3_DIGITS
        _z3digit(z, k) == 3 || return false
    end
    if root == 0 || root == 11
        for k in 1:level
            _z3digit(z, k) == 0 || return false
        end
    end
    return true
end

function Z3Cell(s::AbstractString)
    ncodeunits(s) >= 2 || throw(ArgumentError("Z3 text needs a two-digit root"))
    root = tryparse(Int, s[1:2])
    root === nothing && throw(ArgumentError("invalid Z3 root in $s"))
    0 <= root <= 11 || throw(ArgumentError("Z3 root must be 00:11"))
    level = ncodeunits(s) - 2
    level <= Z3_DIGITS || throw(ArgumentError("Z3 supports at most $Z3_DIGITS digits"))
    z = (UInt64(root) << Z3_ROOT_SHIFT) | Z3_PAD
    for k in 1:level
        ch = s[k + 2]
        '0' <= ch <= '2' || throw(ArgumentError("Z3 digits must be 0:2"))
        sh = _z3shift(k)
        z = (z & ~(UInt64(3) << sh)) | (UInt64(ch - '0') << sh)
    end
    _valid_z3(z) || throw(ArgumentError("$s does not name an ISEA3H cell"))
    return Z3Cell(z)
end

function z3_string(c::Z3Cell)
    _valid_z3(c.id) || throw(ArgumentError("invalid Z3 cell 0x$(string(c.id, base=16))"))
    io = IOBuffer()
    print(io, lpad(_z3root(c.id), 2, '0'))
    for k in 1:DGG.level(c)
        print(io, _z3digit(c.id, k))
    end
    return String(take!(io))
end

Base.show(io::IO, c::Z3Cell) = print(io, "Z3Cell(\"", z3_string(c), "\")")

@inline function _z3from(root::Int, path::Int64, level::Int)
    z = (UInt64(root) << Z3_ROOT_SHIFT) | Z3_PAD
    for k in level:-1:1
        d = path % 3
        path ÷= 3
        sh = _z3shift(k)
        z = (z & ~(UInt64(3) << sh)) | (UInt64(d) << sh)
    end
    return Z3Cell(z)
end

@inline function _z3path(c::Z3Cell)
    p = Int64(0)
    for k in 1:DGG.level(c)
        p = 3p + _z3digit(c.id, k)
    end
    return p
end
