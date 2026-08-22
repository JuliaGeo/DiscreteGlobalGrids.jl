# Helpers shared by more than one space, tree, or method.

# Chunk-local index maps

# Contiguous index sets use arithmetic; other index sets use a lookup table.

struct OffsetIndexMap
    offset::Int
    n::Int
end

struct LookupIndexMap
    lookup::Dict{Int,Int}
    n::Int
end

"""
    indexmap(inds) -> map

Return a map from global cell positions to their one-based positions within
`inds`. Query it with [`localindex`](@ref); `length` is `length(inds)`.
"""
indexmap(inds::AbstractUnitRange{<:Integer}) =
    OffsetIndexMap(Int(first(inds)) - 1, length(inds))
indexmap(inds) =
    LookupIndexMap(Dict{Int,Int}(Int(p) => k for (k, p) in enumerate(inds)), length(inds))

"""
    localindex(map, i::Integer) -> Int

Return `i`'s one-based position within the mapped index set, or `0` when `i` is
not in it.
"""
@inline localindex(m::OffsetIndexMap, i::Integer) =
    (k = Int(i) - m.offset; 1 <= k <= m.n ? k : 0)
@inline localindex(m::LookupIndexMap, i::Integer) = get(m.lookup, Int(i), 0)

Base.length(m::OffsetIndexMap) = m.n
Base.length(m::LookupIndexMap) = m.n

# Spherical caps

const _WHOLE_SPHERE = SphericalCap(USPoint(0.0, 0.0, 1.0), nextfloat(Float64(pi)))

# Nudged outward past dot-product rounding noise, so containment stays closed.
@inline _padcap(r::Real) = nextfloat(Float64(r) * 1.0001 + 1e-12)

# Area of a unit-sphere cap of angular radius `r`.
_caparea(r::Float64) = r >= pi ? 4pi : 2pi * (1 - cos(r))

# Area of the part of a cap of angular radius `r` lying beyond a great circle
# `h` radians from the cap's centre. `h <= -r` keeps the whole cap, `h >= r`
# none of it, and `h == 0` exactly half.
function _capsegment(r::Float64, h::Float64)
    h >= r && return 0.0
    h <= -r && return _caparea(r)
    s = clamp(sin(h) / sin(r), -1.0, 1.0)
    t = clamp(tan(h) / tan(r), -1.0, 1.0)
    return 2 * (acos(s) - cos(r) * acos(t))
end

# The lens two caps share, for radii of at most a quarter turn. The great
# circle through the two crossing points splits it: each cap contributes the
# segment on the far side of that circle, `h1` radians from its centre.
function _caplens(r1::Float64, r2::Float64, d::Float64)
    d >= r1 + r2 && return 0.0
    d <= abs(r1 - r2) && return _caparea(min(r1, r2))
    # `atan`'s two-argument form puts `h1` in the half-turn where `cos(h1)`
    # matches `cos(r1)`, which is what the right spherical triangle requires.
    h1 = atan(cos(r2) - cos(r1) * cos(d), cos(r1) * sin(d))
    return _capsegment(r1, h1) + _capsegment(r2, d - h1)
end

"""
    _capoverlap(r1, r2, d) -> Float64

Return the area two unit-sphere caps share, given their angular radii and the
angle `d` between their centres. A cap wider than a quarter turn goes through
its complement — itself a cap, of radius `pi - r` about the antipode — because
the lens construction only holds while both caps are at most that wide.
"""
function _capoverlap(r1::Float64, r2::Float64, d::Float64)
    r1 >= pi && return _caparea(r2)
    r2 >= pi && return _caparea(r1)
    half = Float64(pi) / 2
    if r1 > half && r2 > half
        # Inclusion-exclusion on the two complements, which share the angle `d`.
        return 4pi - _caparea(pi - r1) - _caparea(pi - r2) + _caplens(pi - r1, pi - r2, d)
    elseif r1 > half
        # The second cap, less the part of it outside the first.
        return _caparea(r2) - _caplens(pi - r1, r2, pi - d)
    elseif r2 > half
        return _caparea(r1) - _caplens(r1, pi - r2, pi - d)
    end
    return _caplens(r1, r2, d)
end

# The angle between two caps' centres.
_capdistance(a, b) = Float64(US.spherical_distance(a.point, b.point))

# Maximum cells per tree leaf.

# Threading policy

# True while an outer loop already runs one task per block, so nested weight
# builds must not spawn their own.
const OUTER_PARALLEL = ScopedValue(false)

# Outer parallelism wins: thread inner weight builds only at top level.
_innerthreaded() =
    Threads.nthreads() > 1 && !OUTER_PARALLEL[] ? GOCore.True() : GOCore.False()
