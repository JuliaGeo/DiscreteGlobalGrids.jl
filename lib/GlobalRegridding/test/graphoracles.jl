# Dependency-graph oracles, shared by every site that checks the relation.
#
# Three call sites need the same definitions and had drifted into three
# hand-copied ones:
#
#   - `lib/GlobalRegridding/test/test_chunkgraph.jl`  (this subpackage's suite)
#   - `test/systems/crosssystem/regrid.jl`            (the root suite)
#   - `benchmark/chunk_graph_gates.jl`                (the G1 harness)
#
# They live here because this is the innermost of the three: the root suite and
# `benchmark/` both already depend on GlobalRegridding, so including a file from
# its test directory adds no dependency edge in the wrong direction. Follows the
# `test/helpers.jl` convention — a module, `include`d by path from the suites
# that want it, defined inside the including module.
#
# Nothing here may consult a cap, a chunk extent or a chunk index except where
# the docstring says so: that independence is what makes `contributing_pairs` an
# oracle rather than a second opinion.

module ChunkGraphOracles

import GlobalRegridding as GR
import GeoInterface as GI
import GeometryOps as GO

const US = GO.UnitSpherical
const USPoint = GO.UnitSphericalPoint{Float64}

export contributing_pairs, graph_pairs, demanded_pairs, capjoin_pairs, eager_pairs

"""
    graph_pairs(g) -> Set{Tuple{Int,Int}}

The `(dstchunk, srcchunk)` relation a `ChunkDependencyGraph` holds, read out of
its destination-major rows.
"""
graph_pairs(g) = Set((d, Int(s)) for d in 1:GR.ndestinationchunks(g)
                     for s in GR.sourcesof(g, d))

"""
    demanded_pairs(dst, src; radius = 0.0) -> Set{Tuple{Int,Int}}

Every pair a lazy read can ask for: one `candidatechunks!` per destination cap
on the source space's own chunk index. A relation that does not contain this set
cannot back a refcount, because a source it retired is still going to be
demanded.

Post-#69 `chunk_dependency_graph` is built from exactly these queries, so with
no `refine` the graph and this set must be *equal*, not merely nested. That
equality is the contract a builder swap has to preserve, which is why the tests
assert it rather than a containment.
"""
function demanded_pairs(dst, src; radius = 0.0)
    index = GR.chunkindex(src)
    out = Set{Tuple{Int,Int}}()
    buf = Int[]
    for (d, cap) in pairs(GR.chunkextents(dst))
        GR.candidatechunks!(buf, index, cap; radius)
        for s in buf
            push!(out, (d, s))
        end
    end
    return out
end

"""
    capjoin_pairs(dst, src; radius = 0.0) -> Set{Tuple{Int,Int}}

Brute-force `DilatedIntersects` over `chunkextents(dst) × chunkextents(src)` —
the definition of the relation the deleted latitude-sorted join computed, and a
second superset of the geometric truth. It is *not* a bound on the graph in
either direction on the shipped native hierarchies; see the testsets that pin
that.
"""
capjoin_pairs(dst, src; radius = 0.0) =
    let dcaps = GR.chunkextents(dst), scaps = GR.chunkextents(src),
        hits = GR.DilatedIntersects(Float64(radius))

        Set((d, s) for d in eachindex(dcaps), s in eachindex(scaps)
            if hits(dcaps[d], scaps[s]))
    end

"""
    eager_pairs(dst, src; method = GR.Conservative()) -> Set{Tuple{Int,Int}}

Every `(dstchunk, srcchunk)` on which the **eager** path puts a nonzero weight:
the whole-domain block `method` builds for this pair, read entry by entry and
mapped to chunks.

This is the relation a lazy read must be free to reproduce. It differs from
[`contributing_pairs`](@ref) in what it is evidence about: that one is the
geometric truth, computed independently of the weight builder, while this one is
the weight builder's own answer. A dependency relation that does not contain it
would let a lazy read miss a source chunk the eager path weighted — the exact
way lazy and eager values could disagree.

Reads no cap, no chunk extent and no chunk index, so it is an oracle for the
relation rather than a second opinion. `O(ncells(dst) * ncells(src))` lookups;
for small spaces only.
"""
function eager_pairs(dst, src; method = GR.Conservative())
    W = GR.wholeblock(method, dst, src).weights
    out = Set{Tuple{Int,Int}}()
    for j in axes(W, 2), i in axes(W, 1)
        iszero(W[i, j]) && continue
        d, s = GR.chunkat(dst, i), GR.chunkat(src, j)
        # A partial space reports `nothing` for a position it does not cover.
        (d === nothing || s === nothing) && continue
        push!(out, (Int(d), Int(s)))
    end
    return out
end

"""
    contributing_pairs(dst, src; radius = 0.0) -> Set{Tuple{Int,Int}}

The geometric truth a dependency relation must dominate: every `(dstchunk,
srcchunk)` holding at least one cell pair with positive spherical intersection
area, plus — at nonzero support — every pair whose cells come within `radius`
radians of each other.

Built from real cell geometry through the same intersection kernel the
conservative weight builder uses, so a pair outside this set is a pair no weight
could be nonzero on. Nothing here consults a cap, a chunk extent or a chunk
index, which is what makes it an oracle rather than a second opinion.

`O(ncells(dst) * ncells(src))`; for small spaces only.
"""
function contributing_pairs(dst, src; radius = 0.0)
    op = GR._intersectionoperator(GR.manifold(dst))
    nd, ns = Int(GR.ncells(dst)), Int(GR.ncells(src))
    dcells = [GR.getcell(dst, i) for i in 1:nd]
    scells = [GR.getcell(src, j) for j in 1:ns]
    dchunk = [GR.chunkat(dst, i) for i in 1:nd]
    schunk = [GR.chunkat(src, j) for j in 1:ns]
    r = Float64(radius)
    dring = r > 0 ? [cellvertices(dst, i) for i in 1:nd] : nothing
    sring = r > 0 ? [cellvertices(src, j) for j in 1:ns] : nothing
    out = Set{Tuple{Int,Int}}()
    for i in 1:nd, j in 1:ns
        # A partial space reports `nothing` for a position it does not cover.
        (dchunk[i] === nothing || schunk[j] === nothing) && continue
        key = (Int(dchunk[i]), Int(schunk[j]))
        # Rows repeat heavily; skipping a known pair skips the clip too.
        key in out && continue
        hit = op(dcells[i], scells[j]) > 0
        if !hit && r > 0
            hit = ringdistance(dring[i], sring[j]) <= r
        end
        hit && push!(out, key)
    end
    return out
end

"Cell rings carry unit-sphere points, not geographic ones."
cellvertices(space, i) = [USPoint(GI.x(p), GI.y(p), GI.z(p))
                          for p in GI.getpoint(GI.getexterior(GR.getcell(space, i)))]

"""
    ringdistance(a, b) -> Float64

Minimum spherical distance, in radians, between two cell boundary rings given as
vectors of unit-sphere points.

This is a **cell-to-cell** distance, not a vertex-to-vertex one. Two cells can
lie within a support radius edge-to-edge with no vertex pair that close, and a
vertex-only test would then under-approximate the truth exactly in the cases the
nonzero-support oracle exists to cover. For two rings with no positive
intersection area — the only case `contributing_pairs` reaches this from — the
minimum is attained at a vertex of one and a point of an edge of the other, so
testing every vertex against every edge in both directions is exact.
"""
ringdistance(a::AbstractVector, b::AbstractVector) =
    min(_vertices_to_edges(a, b), _vertices_to_edges(b, a))

_vertices_to_edges(pts, ring) = minimum(_arcdistance(p, ring) for p in pts)

function _arcdistance(p, ring)
    n = length(ring)
    best = Inf
    for i in 1:n
        best = min(best, _segmentdistance(p, ring[i], ring[i == n ? 1 : i+1]))
    end
    return best
end

_dot(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
_cross(a, b) = USPoint(a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])

# Distance from `p` to the minor great-circle arc `a -> b`. The endpoints always
# bound it; the perpendicular foot improves on them only when it falls inside
# the arc.
function _segmentdistance(p, a, b)
    best = min(US.spherical_distance(p, a), US.spherical_distance(p, b))
    n = _cross(a, b)
    nn = sqrt(_dot(n, n))
    nn > 0 || return best                       # degenerate edge: a == ±b
    n = n / nn
    h = _dot(p, n)                              # sine of the signed angle to the plane
    q = p - h * n
    qq = sqrt(_dot(q, q))
    qq > 0 || return best                       # p is a pole of the great circle
    q = q / qq
    inside = _dot(_cross(a, q), n) >= 0 && _dot(_cross(q, b), n) >= 0
    return inside ? min(best, abs(asin(clamp(h, -1.0, 1.0)))) : best
end

end # module
