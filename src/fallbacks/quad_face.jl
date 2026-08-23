# The radix-4 quad-face family: one aligned `2^l x 2^l` lattice per base face,
# ids `face * 4^l + curvecode`, positions `id + 1`. Everything that follows from
# that identity alone is written once here and dispatches on
# `AbstractQuadFaceGridSystem`; the per-system files keep their charts,
# adjacency, extents, and the three declarations below.

"""
    nbasefaces(sys) -> Int
    systemname(sys) -> String
    idname(sys) -> String

What an [`AbstractQuadFaceGridSystem`](@ref) declares about itself: how many base
faces it has, the name its errors call the system by, and the noun its errors
call a canonical id by (`"nested id"`, `"scaffold ordinal"`, …).
"""
function nbasefaces end

@doc (@doc nbasefaces)
function systemname end

@doc (@doc nbasefaces)
function idname end

"""
    subtree_curve(sys) -> curve
    subtree_orientation(sys, c) -> UInt8

The space-filling curve a subtree's square block is walked in, and the curve
state that block's own root is read under. They default to [`MortonCurve`](@ref)
and `0x0`; a system whose curve advances with depth overrides both.
"""
subtree_curve(::AbstractQuadFaceGridSystem) = MortonCurve()

@doc (@doc subtree_curve)
subtree_orientation(::AbstractQuadFaceGridSystem, ::LevelIndex) = 0x0

"""
    nside(level) -> Int64

The face side `2^level`, in cells.
"""
@inline nside(level::Integer) = Int64(1) << Int(level)

"""
    checked_id(sys, c) -> Int64
    checked_id(grid, c) -> Int64

The raw id of `c` after checking it against the system's id range at `c`'s own
level, and in the grid form against the grid's level too. Chart decoders read
this rather than [`rawid`](@ref): an out-of-range id would otherwise decode to
the geometry of a cell that does not exist.
"""
@inline function checked_id(sys::AbstractQuadFaceGridSystem, c::LevelIndex)
    l = level(c)
    n = ncells(sys, l)
    0 <= c.index < n || throw(ArgumentError(
        "$(idname(sys)) $(c.index) is out of range 0:$(n - 1) at level $l"))
    return c.index
end

@doc (@doc checked_id)
@inline function checked_id(g::HierarchicalLevelGrid{<:AbstractQuadFaceGridSystem},
        c::LevelIndex)
    level(c) == g.level || throw(ArgumentError(
        "cell $c is at level $(level(c)), not the grid's level $(g.level)"))
    return checked_id(g.system, c)
end

# ===========================================================================
# Identity and hierarchy
# ===========================================================================

cellindextype(::AbstractQuadFaceGridSystem) = LevelIndex
has_sorted_subtrees(::AbstractQuadFaceGridSystem) = true

"""
    rootcells(sys::AbstractQuadFaceGridSystem)

The base faces as level-0 cells `LevelIndex(0, 0:nbasefaces(sys)-1)`: at level 0
the curve code is empty, so the id *is* the face number.
"""
rootcells(sys::AbstractQuadFaceGridSystem) =
    [LevelIndex(0, i) for i in 0:(nbasefaces(sys) - 1)]

"""
    parent(sys::AbstractQuadFaceGridSystem, c) -> LevelIndex

The radix-4 parent `index ÷ 4`, one level up. Throws an `ArgumentError` on a
level-0 cell, which is a root.
"""
function Base.parent(sys::AbstractQuadFaceGridSystem, c::LevelIndex)
    l = level(c)
    l > 0 || throw(ArgumentError(
        "level-0 $(systemname(sys)) cell $c is a root and has no parent"))
    return LevelIndex(l - 1, c.index >> 2)
end

"""
    children(sys::AbstractQuadFaceGridSystem, c)

The four children `4*index .+ (0:3)`, ascending at the next level. Refinement is
a uniform quadtree, so every cell has exactly four. Throws an `ArgumentError` at
`maxlevel`.
"""
function children(sys::AbstractQuadFaceGridSystem, c::LevelIndex)
    l = level(c)
    l < maxlevel(sys) || throw(ArgumentError(
        "$(systemname(sys)) cell $c is at maxlevel $(maxlevel(sys)) and has no children"))
    base = c.index << 2
    return [LevelIndex(l + 1, base + k) for k in 0:3]
end

"""
    ancestor(sys::AbstractQuadFaceGridSystem, c, l) -> LevelIndex

The ancestor at level `l`, in one shift: `index >> 2Δ`.

Sound because the id is `face * 4^level + curvecode` and the curve code is
positional — dropping its low `2Δ` bits drops `Δ` bits from each lattice
coordinate, which is exactly `Δ` steps up the quadtree, and the face term divides
through untouched.
"""
function ancestor(sys::AbstractQuadFaceGridSystem, c::LevelIndex, l::Integer)
    target = Int(l)
    lc = level(c)
    target <= lc || throw(ArgumentError(
        "ancestor level $target is deeper than the cell's own level $lc"))
    target >= 0 || throw(ArgumentError(
        "ancestor level $target is above the root level 0"))
    return LevelIndex(target, c.index >> (2 * (lc - target)))
end

"""
    descendant_range(sys::AbstractQuadFaceGridSystem, c, l) -> UnitRange{Int}

The contiguous **positions** in `levelgrid(sys, l)` occupied by `c`'s level-`l`
descendants: ids `index * 4^Δ` through `(index + 1) * 4^Δ - 1`, shifted into
1-based positions.

Exact and hole-free in both directions, which is what
`has_sorted_subtrees(sys) == true` asserts: curve order is depth-first order by
construction, the id space is dense, and sibling ranges tile the parent's
exactly.
"""
function descendant_range(sys::AbstractQuadFaceGridSystem, c::LevelIndex, l::Integer)
    target = Int(l)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "descendant level $target is above the cell's own level $lc"))
    target <= maxlevel(sys) || throw(ArgumentError(
        "descendant level $target is past maxlevel $(maxlevel(sys))"))
    shift = 2 * (target - lc)
    lo = c.index << shift
    hi = ((c.index + 1) << shift) - 1
    return Int(lo + 1):Int(hi + 1)
end

"""
    descendants(sys::AbstractQuadFaceGridSystem, c, l)

Every level-`l` descendant of `c`, ascending: the dense, subtree-contiguous
[`descendant_range`](@ref) read off as consecutive ids, with no `children`
recursion and no sort.
"""
function descendants(sys::AbstractQuadFaceGridSystem, c::LevelIndex, l::Integer)
    r = descendant_range(sys, c, l)          # validates `l` both ways
    target = Int(l)
    return [LevelIndex(target, i - 1) for i in r]
end

# ===========================================================================
# The level grid: size, and positions <-> ids
# ===========================================================================

ncells(sys::AbstractQuadFaceGridSystem, l::Integer) =
    Int(nbasefaces(sys) * (Int64(1) << (2 * Int(l))))

# The grid bounds-checks `i`, so this is the bijection and nothing else.
cellindex(::AbstractQuadFaceGridSystem, l::Integer, i::Int) = LevelIndex(l, i - 1)

"""
    cellposition(sys::AbstractQuadFaceGridSystem, c) -> Union{Int,Nothing}

`index + 1` for an in-range id, `nothing` otherwise. The grid has already
rejected a cell from another level and converted any alternate scheme.
"""
function cellposition(sys::AbstractQuadFaceGridSystem, c::LevelIndex)
    0 <= c.index < ncells(sys, level(c)) || return nothing
    return Int(c.index + 1)
end

# ===========================================================================
# Subtree engines
#
# A subtree at depth Δ is the aligned `2^Δ x 2^Δ` lattice block on one face, and
# a descendant's offset within the subtree's contiguous id range IS its position
# along the curve inside that block. Ascending id order over the subtree is
# therefore curve order over the block, and the border is the block's perimeter,
# `4*2^Δ - 4` cells. Both connectivities agree: a perimeter cell has an axis
# neighbour outside the block already, and a block flush against a face edge is
# no special case, because its neighbours across the seam are on another face and
# so outside the subtree regardless.
# ===========================================================================

# `descendant_range` validates the target level and yields the block's first id.
function _quad_block(sys::AbstractQuadFaceGridSystem, c::LevelIndex, target::Int)
    r = descendant_range(sys, c, target)
    checked_id(sys, c)
    return Int64(first(r)) - 1, Int64(1) << (target - level(c))
end

function border_engine(sys::AbstractQuadFaceGridSystem, c::LevelIndex, target::Int,
        connectivity::Connectivity)
    lo, side = _quad_block(sys, c, target)
    return SquareBorderEngine(subtree_curve(sys), lo, target, side,
        subtree_orientation(sys, c))
end

function interior_engine(sys::AbstractQuadFaceGridSystem, c::LevelIndex, target::Int,
        connectivity::Connectivity)
    lo, side = _quad_block(sys, c, target)
    return SquareInteriorEngine(subtree_curve(sys), lo, target, side,
        subtree_orientation(sys, c))
end

# The family's `halo_engine` is `Engine`'s: it is the band walk's wiring, and
# the band walk lives beside the other halo engines.

# ===========================================================================
# Chart-independent geometry
# ===========================================================================

"""
    chart_perimeter(chart, ix, iy, face, nside, nseg) -> Vector{UnitSphericalPoint}

The cell's boundary walked in `nseg` equal chart steps per edge, starting at the
`(x+, y+)` corner and running `(x+,y+) → (x-,y+) → (x-,y-) → (x+,y-)`. `chart` is
the system's `(u, v, face) -> point` map on the unit face square.

Implicitly closed: each edge contributes its start vertex and its interior
points, never its end vertex, so the next edge's start is not duplicated.
`nseg == 1` reproduces the four corners.
"""
function chart_perimeter(chart, ix::Integer, iy::Integer, face::Integer,
        nside::Integer, nseg::Integer)
    n = nside
    x0 = Int64(ix)
    y0 = Int64(iy)
    pts = Vector{USPoint}(undef, 4 * nseg)
    k = 0
    for i in 0:(nseg - 1)          # (x+,y+) -> (x-,y+), along y = (iy+1)/n
        t = i / nseg
        pts[k += 1] = chart((x0 + 1 - t) / n, (y0 + 1) / n, face)
    end
    for i in 0:(nseg - 1)          # (x-,y+) -> (x-,y-), along x = ix/n
        t = i / nseg
        pts[k += 1] = chart(x0 / n, (y0 + 1 - t) / n, face)
    end
    for i in 0:(nseg - 1)          # (x-,y-) -> (x+,y-), along y = iy/n
        t = i / nseg
        pts[k += 1] = chart((x0 + t) / n, y0 / n, face)
    end
    for i in 0:(nseg - 1)          # (x+,y-) -> (x+,y+), along x = (ix+1)/n
        t = i / nseg
        pts[k += 1] = chart((x0 + 1) / n, (y0 + t) / n, face)
    end
    return pts
end

"""
    sampled_cap(center, pts) -> SphericalCap

A cap about `center` covering the region whose perimeter `pts` samples: the
sampled maximum radius, plus half the largest gap between consecutive samples,
plus one outward ULP.

For a chart-square cell this bounds the whole subtree, since children tile the
parent's square exactly and the distance from the centre is maximised on the
perimeter — in fact at a corner, and every corner is a sample. `gap/2` is
conservative measured slack rather than a formal Lipschitz bound, because `gap`
is a geodesic chord rather than chart-edge arc length.
"""
function sampled_cap(center, pts)
    rmax = 0.0
    gap = 0.0
    prev = pts[end]
    for p in pts
        rmax = max(rmax, US.spherical_distance(center, p))
        gap = max(gap, US.spherical_distance(prev, p))
        prev = p
    end
    radius = min(Float64(π), rmax + gap / 2)
    return SphericalCap(center, nextfloat(radius))
end

"""
    morton_encode(ix, iy) -> Int64

Interleave a lattice coordinate pair into a Z-order code: `ix` occupies the even
bits and `iy` the odd bits. That positional layout is what makes `code ÷ 4` the
parent and `4code .+ (0:3)` the children.

Defined for `0 <= ix, iy < 2^31`, which covers every representable refinement
level; the low 32 bits of each coordinate are the ones interleaved.
"""
@inline function morton_encode(ix::Integer, iy::Integer)
    return Int64(_spread_bits(UInt64(Int64(ix))) | (_spread_bits(UInt64(Int64(iy))) << 1))
end

# Spread the low 32 bits of `x` into the even bit positions of a 64-bit word by
# repeated halving of the gaps: the classic branchless bit-interleave. Constant
# time, unlike a per-bit loop whose cost grows with the refinement level.
@inline function _spread_bits(x::UInt64)
    x &= 0x00000000ffffffff
    x = (x | (x << 16)) & 0x0000ffff0000ffff
    x = (x | (x << 8))  & 0x00ff00ff00ff00ff
    x = (x | (x << 4))  & 0x0f0f0f0f0f0f0f0f
    x = (x | (x << 2))  & 0x3333333333333333
    x = (x | (x << 1))  & 0x5555555555555555
    return x
end

# Inverse of `_spread_bits`: gather the even bit positions back down into the
# low 32 bits, halving the gaps in the opposite order.
@inline function _gather_bits(x::UInt64)
    x &= 0x5555555555555555
    x = (x | (x >> 1))  & 0x3333333333333333
    x = (x | (x >> 2))  & 0x0f0f0f0f0f0f0f0f
    x = (x | (x >> 4))  & 0x00ff00ff00ff00ff
    x = (x | (x >> 8))  & 0x0000ffff0000ffff
    x = (x | (x >> 16)) & 0x00000000ffffffff
    return x
end

"""
    morton_decode(code) -> (ix, iy)

Inverse of [`morton_encode`](@ref): even bits rebuild `ix`, odd bits rebuild
`iy`.
"""
@inline function morton_decode(code::Integer)
    c = UInt64(Int64(code))
    return (Int64(_gather_bits(c)), Int64(_gather_bits(c >> 1)))
end
