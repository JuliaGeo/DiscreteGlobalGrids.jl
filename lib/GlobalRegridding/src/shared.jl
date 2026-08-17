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

# Maximum cells per tree leaf.
const _CELL_TREE_LEAF = 16
