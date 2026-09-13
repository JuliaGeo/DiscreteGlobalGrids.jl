# Spherical query descent: prune by conservative caps, then apply exact per-cell
# predicates. `Intersects` also uses centroid acceptance and a boundary-distance
# sandwich. Input ring edges are great-circle arcs; densify intended parallels
# (`GO.segmentize`).

# Leaf bucket size balancing cursor overhead against cap selectivity. 16 measured
# well across all four systems of the previous design on box, antimeridian, polar
# and quarter-sphere queries.
const QUERY_BUCKET_SIZE = 16

# ===========================================================================
# Targets
# ===========================================================================

abstract type QueryTarget end

"""
    GeometryTarget(geom)

A prepared query target: the geometry lifted to the unit sphere once, a
`RelateNG` preparation of it, a bounding cap for tree pruning, its boundary
arcs for the border sandwich, and an [`ArcTree`](@ref) over those arcs for the
descent's boundary frontier.
"""
struct GeometryTarget{P,G,A,E} <: QueryTarget
    prepared::P
    geom::G
    cap::Cap
    arcs::A
    edges::E
end

GeometryTarget(prepared, geom, cap::Cap, arcs) =
    GeometryTarget(prepared, geom, cap, arcs, _arc_tree(arcs))

"""
    CapTarget(cap)

A spherical cap as a query target. Handled natively — exactly, and without any
polygonisation — for the predicates whose answer a cap can give directly.
"""
struct CapTarget <: QueryTarget
    cap::Cap
end

_query_target(x::US.SphericalCap) = CapTarget(SphericalCap(
    USPoint(x.point[1], x.point[2], x.point[3]), Float64(x.radius)))

_query_target(x::Extents.Extent) = _query_target(_extent_target(x))

function _query_target(geom)
    GI.isgeometry(geom) || GI.trait(geom) !== nothing || throw(ArgumentError(
        "a query target must be a GeoInterface geometry, an Extents.Extent or a " *
        "SphericalCap; got $(typeof(geom))"))
    spherical = _to_unit_sphere(geom)
    prepared = GO.prepare(GO.RelateNG(; manifold=GO.Spherical()), spherical)
    arcs = _boundary_arcs(spherical)
    return GeometryTarget(prepared, spherical, _geometry_cap(prepared, spherical),
        arcs, _arc_tree(arcs))
end

_query_target(::Nothing) = throw(ArgumentError(
    "the predicate carries no target; write e.g. `Intersects(geometry)`"))

# The target's coordinates on the unit sphere, once, at the boundary of the
# call: 2-D coordinates are lon/lat degrees, 3-D coordinates are xyz as-is.
# Rebuilt through GeoInterface so that any input geometry type lands in the
# concrete form the spherical kernels want.
_to_unit_sphere(geom) = _to_unit_sphere(GI.trait(geom), geom)
_to_unit_sphere(::GI.PointTrait, geom) = query_point(geom)
_to_unit_sphere(::GI.MultiPointTrait, geom) =
    GI.MultiPoint([query_point(p) for p in GI.getpoint(geom)])
_to_unit_sphere(::GI.LineStringTrait, geom) =
    GI.LineString([query_point(p) for p in GI.getpoint(geom)])
_to_unit_sphere(::GI.LinearRingTrait, geom) =
    GI.LinearRing([query_point(p) for p in GI.getpoint(geom)])
_to_unit_sphere(::GI.MultiLineStringTrait, geom) =
    GI.MultiLineString([_to_unit_sphere(GI.LineStringTrait(), l) for l in GI.getgeom(geom)])
_to_unit_sphere(::GI.PolygonTrait, geom) =
    GI.Polygon([GI.LinearRing([query_point(p) for p in GI.getpoint(r)])
                for r in GI.getring(geom)])
_to_unit_sphere(::GI.MultiPolygonTrait, geom) =
    GI.MultiPolygon([_to_unit_sphere(GI.PolygonTrait(), p) for p in GI.getpolygon(geom)])
_to_unit_sphere(t, geom) = throw(ArgumentError(
    "query targets of trait $t are not supported; pass a point, line, polygon, " *
    "multi-geometry, Extents.Extent or SphericalCap"))

# How finely a lon/lat box outline is sampled before becoming a spherical
# polygon. Its edges are great-circle arcs, so a parallel has to be densified
# to be traced; 2 degrees keeps the sag under ~150 m on Earth.
const EXTENT_STEP_DEGREES = 2.0

# Convert a lon/lat extent to a polygon or polar cap. A full-longitude,
# non-polar extent is an annulus and cannot be represented by either.
function _extent_target(ext::Extents.Extent)
    hasproperty(ext, :X) && hasproperty(ext, :Y) || throw(ArgumentError(
        "a query extent needs X and Y bounds in lon/lat degrees, got $(ext)"))
    x0, x1 = Float64.(ext.X)
    y0, y1 = Float64.(ext.Y)
    y0 <= y1 || throw(ArgumentError("extent Y bounds must be ascending, got $(ext.Y)"))
    if x1 - x0 >= 359.999
        y1 >= 89.999 && return SphericalCap(USPoint(0.0, 0.0, 1.0), deg2rad(90.0 - y0))
        y0 <= -89.999 && return SphericalCap(USPoint(0.0, 0.0, -1.0), deg2rad(y1 + 90.0))
        throw(ArgumentError(
            "an extent spanning all longitudes between two parallels is an annulus, " *
            "not a polygon; query it as the difference of two SphericalCaps"))
    end
    x0 <= x1 || throw(ArgumentError(
        "extent X bounds must be ascending; an antimeridian-crossing box is two boxes"))
    steps = max(1, ceil(Int, (x1 - x0) / EXTENT_STEP_DEGREES))
    points = USPoint[]
    # Counter-clockwise seen from outside: east along the south edge, north,
    # west along the north edge, south. A polar edge collapses to one vertex.
    if y0 <= -89.999
        push!(points, unit_point(0.0, -90.0))
    else
        for k in 0:steps
            push!(points, unit_point(x0 + (x1 - x0) * k / steps, y0))
        end
    end
    if y1 >= 89.999
        push!(points, unit_point(0.0, 90.0))
    else
        for k in steps:-1:0
            push!(points, unit_point(x0 + (x1 - x0) * k / steps, y1))
        end
    end
    push!(points, points[1])
    return GI.Polygon([GI.LinearRing(points)])
end

# A radius-`≤ π/2` vertex cap contains all great-circle boundary arcs. Its
# connected, boundary-free complement is wholly inside or outside the geometry,
# decided at the cap antipode. Wider caps or an interior antipode use the full
# sphere.

function _geometry_cap(prepared, geom)
    points = USPoint[]
    for p in GI.getpoint(geom)
        push!(points, query_point(p))
    end
    cap = points_cap(points)
    cap.radius >= Float64(pi) && return full_sphere_cap()
    antipode = USPoint(-cap.point[1], -cap.point[2], -cap.point[3])
    GO.relate_predicate(prepared, GO.pred_intersects(), antipode) && return full_sphere_cap()
    return cap
end

# A point target has no interior to test the antipode against.
_geometry_cap(_, geom::GO.UnitSphericalPoint) = SphericalCap(
    USPoint(geom[1], geom[2], geom[3]), 0.0)

# The border sandwich compares centroid-to-target-boundary distance `d` with cell
# radius bounds:
#
#   * `d > r_out` proves disjointness when the centroid is outside the target.
#   * `d < r_in` proves intersection by placing a target boundary point inside.
#
# The remaining annulus requires the exact predicate. `r_out` is the maximum
# centroid-to-vertex distance; `r_in` is the minimum distance to an edge's
# carrying great circle. `SANDWICH_SLACK` moves both bounds conservatively.

const SANDWICH_SLACK = 1e-6

"""
    BoundaryArc(a, b)

One great-circle edge of the target's boundary, with the derived quantities the
distance test needs: the unnormalised great-circle normal and its squared
length. Precomputed once per query, scanned once per candidate cell.
"""
struct BoundaryArc
    a::USPoint
    b::USPoint
    nx::Float64
    ny::Float64
    nz::Float64
    nn::Float64
end

function BoundaryArc(a::USPoint, b::USPoint)
    nx = a[2] * b[3] - a[3] * b[2]
    ny = a[3] * b[1] - a[1] * b[3]
    nz = a[1] * b[2] - a[2] * b[1]
    return BoundaryArc(a, b, nx, ny, nz, nx * nx + ny * ny + nz * nz)
end

# Enumerate every topological-boundary arc, including holes and all multipart
# rings. Return `nothing` when the geometry cannot provide a safe boundary set.
# Rings are walked one at a time (`GI.getring`, not `GI.getpoint`): an edge
# invented between consecutive points of two different rings is a wrong accept.
_boundary_arcs(geom) = _boundary_arcs!(BoundaryArc[], GI.trait(geom), geom)

_boundary_arcs!(arcs, ::Any, geom) = nothing

function _boundary_arcs!(arcs, ::GI.PolygonTrait, geom)
    _rings_arcs!(arcs, geom)
    return _finish_arcs!(arcs)
end

function _boundary_arcs!(arcs, ::GI.MultiPolygonTrait, geom)
    for poly in GI.getpolygon(geom)
        _rings_arcs!(arcs, poly)
    end
    return _finish_arcs!(arcs)
end

function _boundary_arcs!(arcs, ::GI.LineStringTrait, geom)
    _chain_arcs!(arcs, geom, false)
    return _finish_arcs!(arcs)
end

function _boundary_arcs!(arcs, ::GI.MultiLineStringTrait, geom)
    for line in GI.getgeom(geom)
        _chain_arcs!(arcs, line, false)
    end
    return _finish_arcs!(arcs)
end

function _boundary_arcs!(arcs, ::GI.PointTrait, geom)
    point = query_point(geom)
    push!(arcs, BoundaryArc(point, point))
    return _finish_arcs!(arcs)
end

function _boundary_arcs!(arcs, ::GI.MultiPointTrait, geom)
    for p in GI.getpoint(geom)
        point = query_point(p)
        push!(arcs, BoundaryArc(point, point))
    end
    return _finish_arcs!(arcs)
end

function _rings_arcs!(arcs, poly)
    for r in GI.getring(poly)
        _chain_arcs!(arcs, r, true)
    end
    return arcs
end

# Consecutive-vertex edges of one ring or linestring; `close` adds the edge back
# to the first vertex, which a polygon ring always has whether or not its
# coordinates repeat it (a duplicated closing vertex just makes that edge
# degenerate, and degenerate arcs are handled).
function _chain_arcs!(arcs, chain, close::Bool)
    first_point = nothing
    previous = nothing
    count = 0
    for p in GI.getpoint(chain)
        point = query_point(p)
        count += 1
        previous === nothing ? (first_point = point) :
        push!(arcs, BoundaryArc(previous, point))
        previous = point
    end
    # A one-point chain is that point; a degenerate arc keeps it in the set.
    (close || count == 1) && previous !== nothing &&
        push!(arcs, BoundaryArc(previous, first_point))
    return arcs
end

# Disable the sandwich for empty or near-antipodal edges. Near antipodality
# makes the great-circle normal ill-conditioned and endpoint distance unsafe for
# rejection. Repeated points remain safe and are distinguished by dot product.
# The guard is on `nn = |a × b|²` rather than on exact antipodality, so it covers
# every edge within 1e-6 rad of a half circle, where the normal's direction has
# already lost its accuracy.
const ANTIPODAL_NN = 1e-12

function _finish_arcs!(arcs)
    isempty(arcs) && return nothing
    for arc in arcs
        arc.nn > ANTIPODAL_NN && continue
        arc.a[1] * arc.b[1] + arc.a[2] * arc.b[2] + arc.a[3] * arc.b[3] < 0 && return nothing
    end
    return arcs
end

# `cos` of the distance from `c` to the great-circle arc — the nearer endpoint
# unless the foot of the perpendicular falls on the arc itself. Everything is
# dot and cross products on the raw coordinates: no normalisation, one `sqrt`,
# no inverse trigonometry.
@inline function _arc_cos_distance(arc::BoundaryArc, c)
    ca = c[1] * arc.a[1] + c[2] * arc.a[2] + c[3] * arc.a[3]
    cb = c[1] * arc.b[1] + c[2] * arc.b[2] + c[3] * arc.b[3]
    best = max(ca, cb)
    arc.nn > 0 || return best        # degenerate arc: its endpoints are all of it
    # The near foot is `f ~ c - ((c.n)/(n.n)) n`, and that second term drops out
    # of both triple products below, so `a -> f -> b` in `n`'s orientation is
    # decided by `c` directly. If it fails, the arc's minimum is at an endpoint.
    ((arc.a[2] * c[3] - arc.a[3] * c[2]) * arc.nx +
     (arc.a[3] * c[1] - arc.a[1] * c[3]) * arc.ny +
     (arc.a[1] * c[2] - arc.a[2] * c[1]) * arc.nz) > 0 || return best
    ((c[2] * arc.b[3] - c[3] * arc.b[2]) * arc.nx +
     (c[3] * arc.b[1] - c[1] * arc.b[3]) * arc.ny +
     (c[1] * arc.b[2] - c[2] * arc.b[1]) * arc.nz) > 0 || return best
    cn = c[1] * arc.nx + c[2] * arc.ny + c[3] * arc.nz
    # sin(distance to the great circle) = |c.n| / |n|
    return max(best, sqrt(max(0.0, 1.0 - cn * cn / arc.nn)))
end

"""
    ArcTree(arcs)

The boundary arcs sorted by Morton key under a binary tree of caps. The descent
carries a frontier of its nodes — the ones whose arcs may still reach the grid
node in hand — and narrows it cap by cap. An empty frontier under a convex cap
proves that cap free of the target's boundary.
"""
struct ArcTree
    arcs::Vector{BoundaryArc}
    caps::Vector{Cap}
    first::Vector{Int}
    last::Vector{Int}
    left::Vector{Int}
    right::Vector{Int}
end

const ARC_TREE_LEAF_SIZE = 8

function ArcTree(arcs::Vector{BoundaryArc})
    caps = [points_cap((a.a, a.b)) for a in arcs]
    order = sortperm([_morton_key(c.point) for c in caps])
    tree = ArcTree(arcs[order], Cap[], Int[], Int[], Int[], Int[])
    _arc_node!(tree, caps[order], 1, length(arcs))
    return tree
end

function _arc_node!(tree::ArcTree, caps::Vector{Cap}, lo::Int, hi::Int)
    n = length(tree.caps) + 1
    push!(tree.caps, full_sphere_cap())
    push!(tree.first, lo)
    push!(tree.last, hi)
    push!(tree.left, 0)
    push!(tree.right, 0)
    if hi - lo < ARC_TREE_LEAF_SIZE
        tree.caps[n] = _merge_range(caps, lo, hi)
    else
        mid = (lo + hi) ÷ 2
        l = _arc_node!(tree, caps, lo, mid)
        r = _arc_node!(tree, caps, mid + 1, hi)
        tree.left[n] = l
        tree.right[n] = r
        tree.caps[n] = Extents.union(tree.caps[l], tree.caps[r])
    end
    return n
end

_arc_tree(::Nothing) = nothing
_arc_tree(arcs::Vector{BoundaryArc}) = ArcTree(arcs)

_target_edges(::QueryTarget) = nothing
_target_edges(target::GeometryTarget) = target.edges

# Narrow frontier node `n` to `cap`. A tree node smaller than the cap is kept
# whole: splitting it costs more than the arcs it would exclude. A leaf is
# settled by exact arc distance, so a frontier never holds an arc that cannot
# reach the cap.
function _edge_frontier!(out::Vector{Int}, tree::ArcTree, n::Int, cap)
    Extents.intersects(tree.caps[n], cap) || return out
    if tree.left[n] == 0
        threshold = cos(min(Float64(pi),
            Float64(cap.radius) * (1 + SANDWICH_SLACK) + 1e-12))
        for k in tree.first[n]:tree.last[n]
            if _arc_cos_distance(tree.arcs[k], cap.point) >= threshold
                push!(out, n)
                break
            end
        end
    elseif tree.caps[n].radius > cap.radius
        _edge_frontier!(out, tree, tree.left[n], cap)
        _edge_frontier!(out, tree, tree.right[n], cap)
    else
        push!(out, n)
    end
    return out
end

function _narrow!(out::Vector{Int}, tree::ArcTree, frontier::Vector{Int}, cap)
    empty!(out)
    for n in frontier
        _edge_frontier!(out, tree, n, cap)
    end
    return out
end

# `1` the cell provably meets the target, `-1` it provably does not, `0`
# undecided — see the note above for why each arm is a proof.
function _sandwich(arcs::Vector{BoundaryArc}, centroid, ring)
    bounds = _sandwich_bounds(centroid, ring)
    bounds === nothing && return 0
    cos_out, cos_in = bounds
    best = -1.0
    for arc in arcs
        value = _arc_cos_distance(arc, centroid)
        value > cos_in && return 1          # d < r_in: a boundary point is inside the cell
        best = max(best, value)
    end
    return best < cos_out ? -1 : 0          # d > r_out: no boundary point is near
end

# The same verdict from the frontier's arcs alone: an arc excluded at an
# ancestor cap cannot reach a cell inside it.
function _sandwich(tree::ArcTree, frontier::Vector{Int}, centroid, ring)
    bounds = _sandwich_bounds(centroid, ring)
    bounds === nothing && return 0
    cos_out, cos_in = bounds
    best = -1.0
    for n in frontier, k in tree.first[n]:tree.last[n]
        value = _arc_cos_distance(tree.arcs[k], centroid)
        value > cos_in && return 1
        best = max(best, value)
    end
    return best < cos_out ? -1 : 0
end

# `(cos r_out, cos r_in)` for one cell, or `nothing` past a quarter sphere.
function _sandwich_bounds(centroid, ring)
    vertex = 0.0
    for p in ring
        vertex = max(vertex, US.spherical_distance(centroid, USPoint(p[1], p[2], p[3])))
    end
    r_out = vertex * (1 + SANDWICH_SLACK)
    # A quarter sphere is where both the concavity argument for r_out and the
    # sin/cos comparisons below stop holding. No cell of any system in scope
    # comes near it, and giving up is free.
    r_out < Float64(pi) / 2 || return nothing
    cos_out = cos(r_out)
    # Every edge has to be in the minimum or the bound could come out too
    # large, so the loop wraps rather than trusting the ring to repeat its
    # first vertex (where it does, that edge is degenerate and drops out).
    sin2_in = 1.0
    @inbounds for i in eachindex(ring)
        a = ring[i]
        b = ring[i == lastindex(ring) ? firstindex(ring) : i+1]
        nx = a[2] * b[3] - a[3] * b[2]
        ny = a[3] * b[1] - a[1] * b[3]
        nz = a[1] * b[2] - a[2] * b[1]
        nn = nx * nx + ny * ny + nz * nz
        nn > 0 || continue
        cn = centroid[1] * nx + centroid[2] * ny + centroid[3] * nz
        sin2_in = min(sin2_in, cn * cn / nn)
    end
    r_in = asin(sqrt(min(1.0, sin2_in)))
    # A ring that yields no usable inradius (all-degenerate edges, or one wider
    # than the cap it sits in) keeps only the reject arm: `cos_in > 1` can never
    # be exceeded.
    cos_in = 0 < r_in < r_out ? cos(r_in * (1 - SANDWICH_SLACK)) : 2.0
    return (cos_out, cos_in)
end

# ===========================================================================
# Proving a whole subtree outside the target
# ===========================================================================

# The sandwich's reject arm asked of a NODE EXTENT rather than one cell's ring,
# so the answer covers the node's whole subtree.
#
# A cap no boundary arc reaches is connected and free of the target's boundary,
# hence wholly interior or wholly exterior; its centre says which. Exterior
# prunes: `node_extent`'s covering law carries that from the cap to every
# descendant.
#
# The polygon analogue of the wide-cap complement move. A target's bounding cap
# is the cheap prune, and `_geometry_cap` answers the whole sphere whenever the
# target's antipode is interior — anything wider than a hemisphere — leaving no
# prune at all. It pays on ordinary targets too: California fills a small part
# of the 6.6-degree disc that bounds it, and the rest is pruned here instead of
# descended into.
#
# `nothing` arcs (empty or near-antipodal boundary, see `_finish_arcs!`) means
# no proof, and no proof means no prune.
_subtree_outside(::QueryTarget, extent) = false

function _subtree_outside(target::GeometryTarget, extent)
    arcs = target.arcs
    arcs === nothing && return false
    # `cos` decreases on `[0, pi]`, so "farther than the radius" is "cosine
    # below the radius's". The slack widens the radius, which is conservative: a
    # cap that only just clears the boundary is not trusted. The clamp keeps a
    # full-sphere extent — what `AuthalicSystem` returns where its warp bound
    # cannot contain a node — from wrapping past the antipode.
    threshold = cos(min(Float64(pi), Float64(extent.radius) * (1 + SANDWICH_SLACK)))
    for arc in arcs
        _arc_cos_distance(arc, extent.point) >= threshold && return false
    end
    return !GO.relate_predicate(target.prepared, GO.pred_contains(), extent.point)
end

# ===========================================================================
# The per-cell exact tests
# ===========================================================================

"""
    CentroidCovered()

Internal predicate: the cell's centroid lies on or inside the target — the
raster `boundary=:center` rule. A point rule rather than a DE9IM relation of
the cell polygon, so `query`'s public predicate set leaves it out.
"""
struct CentroidCovered end

const QueryPredicate = Union{DE9IM.DE9IMPredicate,CentroidCovered}

# On a point, `pred_intersects` is `covers`, so a centroid on the boundary
# counts; `pred_contains` would drop it.
_covers_centroid(target::GeometryTarget, grid, c) =
    GO.relate_predicate(target.prepared, GO.pred_intersects(), cell_centroid(grid, c))
_covers_centroid(target::CapTarget, grid, c) = cap_contains(target.cap, cell_centroid(grid, c))

_matches(::CentroidCovered, target::QueryTarget, grid, c) = _covers_centroid(target, grid, c)

# A boundary-free cell inside the target meets it, lies within it and is covered
# by it; it can contain, cover, touch, overlap or equal the target only through
# a boundary it does not hold.
_accepts_interior(::Union{DE9IM.Intersects,DE9IM.Within,DE9IM.CoveredBy,CentroidCovered}) = true
_accepts_interior(::QueryPredicate) = false

# Which `pred_*` functor answers the predicate with the TARGET prepared as the
# `A` side. A DE9IM predicate wraps its second argument (`Intersects(t)` means
# "cell intersects t"), so every asymmetric relation is spelled here as its
# converse. `nothing` means "not implemented".
_converse_predicate(::DE9IM.Intersects) = GO.pred_intersects()
_converse_predicate(::DE9IM.Within) = GO.pred_contains()      # cell within t <=> t contains cell
_converse_predicate(::DE9IM.Contains) = GO.pred_within()
_converse_predicate(::DE9IM.Covers) = GO.pred_coveredby()     # cell covers t <=> t coveredby cell
_converse_predicate(::DE9IM.CoveredBy) = GO.pred_covers()
_converse_predicate(::DE9IM.Touches) = GO.pred_touches()      # symmetric
_converse_predicate(::DE9IM.Overlaps) = GO.pred_overlaps()    # symmetric
_converse_predicate(::DE9IM.Equals) = GO.pred_equalstopo()    # symmetric
_converse_predicate(pred::DE9IM.DE9IMPredicate) = nothing

function _check_predicate(pred::DE9IM.DE9IMPredicate)
    pred isa DE9IM.Disjoint && return nothing
    _converse_predicate(pred) === nothing && throw(ArgumentError(
        "$(typeof(pred).name.name) is not implemented by this query engine; " *
        "implemented: Intersects, Disjoint, Contains, Within, Covers, CoveredBy, " *
        "Touches, Overlaps, Equals"))
    return nothing
end

# `Intersects` against a prepared geometry: the centroid fast accept, then the
# border sandwich, then the exact polygon predicate. The other predicates go
# straight to their exact test — a cell whose centroid is inside the target is
# not thereby within/covering/touching it.
function _matches(pred::DE9IM.Intersects, target::GeometryTarget, grid, c)
    centroid = cell_centroid(grid, c)
    GO.relate_predicate(target.prepared, GO.pred_contains(), centroid) && return true
    ring = closed_ring(cell_boundary(grid, c))
    if target.arcs !== nothing
        verdict = _sandwich(target.arcs, centroid, ring)
        verdict == 0 || return verdict > 0
    end
    return GO.relate_predicate(target.prepared, GO.pred_intersects(),
        GI.Polygon([GI.LinearRing(ring)]))
end

# A cell within the target holds its centroid inside it, so a centroid the
# target does not even cover settles `Within` without the polygon test.
function _matches(::DE9IM.Within, target::GeometryTarget, grid, c)
    _covers_centroid(target, grid, c) || return false
    return GO.relate_predicate(target.prepared, GO.pred_contains(), cell_polygon(grid, c))
end

_matches(pred::DE9IM.DE9IMPredicate, target::GeometryTarget, grid, c) =
    GO.relate_predicate(target.prepared, _converse_predicate(pred), cell_polygon(grid, c))

# --- cap targets, answered exactly and without polygonising the cap ---------

_matches(::DE9IM.Intersects, target::CapTarget, grid, c) =
    _cell_meets_cap(target.cap, grid, c, false)

# Exact cell/cap overlap for convex caps. `strict` tests the open cap used as a
# closed cap's complement. Wider caps reach this helper only via their convex
# complements.
function _cell_meets_cap(cap, grid, c, strict::Bool)
    ring, n = open_ring(cell_boundary(grid, c))
    radius = Float64(cap.radius)
    holds(p) = strict ? US.spherical_distance(cap.point, p) < radius :
               cap_contains(cap, p)
    holds(cell_centroid(grid, c)) && return true
    for i in 1:n
        holds(ring[i]) && return true
    end
    # Conservatively accept undecidable centre containment to avoid false misses:
    # by here the cap holds no vertex, so it is wholly inside the cell or wholly
    # outside it, and "outside" would drop a cap sitting deep inside one.
    radius > 0 && point_in_cell(ring, cap.point) !== false && return true
    threshold = cos(min(Float64(pi), radius))
    for i in 1:n
        a = ring[i]
        b = ring[i == n ? 1 : i+1]
        arc = BoundaryArc(USPoint(a[1], a[2], a[3]), USPoint(b[1], b[2], b[3]))
        d = _arc_cos_distance(arc, cap.point)
        (strict ? d > threshold : d >= threshold) && return true
    end
    return false
end

# For `r ≤ π/2`, convexity makes vertex containment sufficient. For wider caps,
# test the open convex complement about the antipode:
#
#     cell within C(p, r)  <=>  cell does not meet C(-p, pi - r)
#
# The complement test is strict so an internally tangent cell remains within.
function _matches(::DE9IM.Within, target::CapTarget, grid, c)
    cap = target.cap
    radius = Float64(cap.radius)
    radius >= Float64(pi) && return true          # the whole sphere holds everything
    if radius <= Float64(pi) / 2
        for p in cell_boundary(grid, c)
            cap_contains(cap, p) || return false
        end
        return true
    end
    point = cap.point
    complement = SphericalCap(USPoint(-point[1], -point[2], -point[3]),
        Float64(pi) - radius)
    return !_cell_meets_cap(complement, grid, c, true)
end

_matches(pred::DE9IM.DE9IMPredicate, ::CapTarget, grid, c) = throw(ArgumentError(
    "$(typeof(pred).name.name) is not implemented for a SphericalCap target; " *
    "implemented for caps: Intersects, Disjoint, Within. Pass a polygon instead."))

# ===========================================================================
# The descent
# ===========================================================================

"""
    query(grid, pred::DE9IM.DE9IMPredicate) -> Vector{<:AbstractCellIndex}
    query(sys, pred::DE9IM.DE9IMPredicate; level) -> Vector{<:AbstractCellIndex}

Return every matching cell as a sorted vector of typed ids.

Implemented predicates: `Intersects`, `Disjoint`, `Contains`, `Within`,
`Covers`, `CoveredBy`, `Touches`, `Overlaps`, `Equals`. `Crosses` throws, since
inverting it is not a matter of naming its converse. Targets: a GeoInterface
geometry, an `Extents.Extent` in lon/lat degrees, or a
`GO.UnitSpherical.SphericalCap`.

`Disjoint` is computed as the full-grid complement of `Intersects` and cannot
prune the output traversal.

The traversal runs over [`treeify`](@ref)`(grid)`, so a system's own tree
serves its queries; a grid that reports a system but no level (stored cells
of mixed levels) gets an [`IndexTree`](@ref).

Caps are handled without polygonization. `Within` uses the complement for caps
wider than a hemisphere. `Intersects` conservatively accepts undecidable cap-
centre containment. Cap targets support only `Intersects`, `Disjoint`, and
`Within`.
"""
function query(grid::AbstractGrid, pred::DE9IM.DE9IMPredicate)
    # Validate the grid, predicate, and target before the empty-grid shortcut.
    n = ncells(grid)
    _check_predicate(pred)
    target = _query_target(Base.parent(pred))
    n == 0 && return _empty_ids(grid)
    return _run_query(grid, pred, target)
end

# The id type of a grid's cells, without asking for a cell.
_id_type(grid::AbstractGrid) = (sys = system(grid);
sys === nothing ? typeof(cellindex(grid, 1)) : cellindextype(sys))

# Empty standalone grids have no discoverable concrete id type.
_empty_ids(grid::AbstractGrid) = (sys = system(grid);
sys === nothing ? AbstractCellIndex[] : cellindextype(sys)[])

function query(sys::AbstractHierarchicalGridSystem, pred::DE9IM.DE9IMPredicate;
        level::Integer)
    # `level` the keyword shadows `level` the function throughout this body.
    return query(levelgrid(sys, level), pred)
end

function _run_query(grid::AbstractGrid, pred::DE9IM.DE9IMPredicate, target::QueryTarget)
    indices = _query_indices(grid, pred, target)
    out = Vector{_id_type(grid)}(undef, length(indices))
    for (k, i) in enumerate(indices)
        out[k] = cellindex(grid, i)
    end
    # Indices ascend in canonical id order for every grid a system produces,
    # but a standalone grid's own order need not, and the contract is sorted
    # BY ID.
    return sort!(out)
end

# The complement, over every cell of the grid. Sound by construction and the
# only honest way to answer it: a cell that no tree node's extent reaches is
# still disjoint from the target, so there is nothing to prune away.
function _run_query(grid::AbstractGrid, ::DE9IM.Disjoint, target::QueryTarget)
    hits = _query_indices(grid, DE9IM.Intersects(nothing), target)
    out = Vector{_id_type(grid)}()
    j = 1
    for i in 1:ncells(grid)
        if j <= length(hits) && hits[j] == i
            j += 1
            continue
        end
        push!(out, cellindex(grid, i))
    end
    return sort!(out)
end

function _query_indices(grid::AbstractGrid, pred::QueryPredicate, target::QueryTarget)
    ncells(grid) == 0 && return Int[]
    out = Int[]
    _descend!(out, _query_tree(grid), grid, pred, target, Int[1], Int[])
    return sort!(out)
end

# `treeify`'s tree, with the cursor's leaf bucket widened for cap selectivity.
# A level-less grid names stored cells of mixed levels, which only the index
# tree can hold.
function _query_tree(grid::AbstractGrid)
    (system(grid) === nothing || level(grid) === nothing) &&
        return IndexTreeNode(IndexTree(grid), 1)
    return _rebucket(treeify(grid), QUERY_BUCKET_SIZE)
end

_rebucket(tree, ::Int) = tree
_rebucket(c::HierarchicalGridCursor, bucket::Int) = typeof(c)(c.grid, c.system, c.top_level,
    c.leaf_level, max(bucket, c.bucket_size), c.level, c.id, c.first_index, c.last_index,
    c.complete_subtree, c.selection)

# `frontier` holds the arc-tree nodes whose arcs may still reach this node.
# Only a convex node cap narrows it: its intersection with every convex
# ancestor cap is connected, and every descendant's geometry lies in that
# intersection by `node_extent`'s covering law. An empty frontier means no
# boundary crosses that region, so one real descendant centroid decides the
# whole node. The node's own cap centre need not lie in the region and is never
# the witness. `_narrow!` empties its output first, so a frontier is never
# narrowed into itself: each node gets a fresh vector, cells share `scratch`.
function _descend!(out, node, grid, pred, target, frontier::Vector{Int}, scratch::Vector{Int})
    extent = STI.node_extent(node)
    Extents.intersects(target.cap, extent) || return nothing
    edges = _target_edges(target)
    if edges !== nothing && Float64(extent.radius) < Float64(pi) / 2
        frontier = _narrow!(Int[], edges, frontier, extent)
        if isempty(frontier)
            index = _witness_index(node)
            index == 0 && return nothing
            _accepts_interior(pred) && _covers_centroid(target, grid, cellindex(grid, index)) &&
                _emit!(out, node)
            return nothing
        end
    end
    if STI.isleaf(node)
        # A node whose whole extent sits inside the target's cap has no per-cell
        # cap prune left to win — every cell cap would pass — so it skips
        # straight to the exact tests.
        interior = Extents.contains(target.cap, extent)
        _leaf_scan!(out, node, grid, pred, target, frontier, scratch, interior)
        return nothing
    end
    for child in STI.getchild(node)
        _descend!(out, child, grid, pred, target, frontier, scratch)
    end
    return nothing
end

# The first cell under a node, or `0` for an empty one.
_witness_index(node::HierarchicalGridCursor) =
    (indices = node_indices(node); isempty(indices) ? 0 : first(indices))

function _witness_index(node)
    if STI.isleaf(node)
        cells = STI.child_indices_extents(node)
        return isempty(cells) ? 0 : first(cells)[1]
    end
    for child in STI.getchild(node)
        index = _witness_index(child)
        index == 0 || return index
    end
    return 0
end

_emit!(out, node::HierarchicalGridCursor) = (append!(out, node_indices(node)); nothing)

function _emit!(out, node)
    if STI.isleaf(node)
        for (index, _) in STI.child_indices_extents(node)
            push!(out, index)
        end
    else
        for child in STI.getchild(node)
            _emit!(out, child)
        end
    end
    return nothing
end

# The leaf scan, once per tree kind: the hierarchical cursor derives a cell's
# cap from its boundary, the index tree has it stored, and any other tree
# lists its leaf caps.
function _leaf_scan!(out, node::HierarchicalGridCursor, grid, pred, target, frontier, scratch, interior::Bool)
    for index in node_indices(node)
        _cell_matches(pred, target, grid, cellindex(grid, index), frontier, scratch, interior) &&
            push!(out, index)
    end
    return nothing
end

function _leaf_scan!(out, node::IndexTreeNode, grid, pred, target, frontier, scratch, interior::Bool)
    tree = node.tree
    for k in tree.node_first[node.index]:tree.node_last[node.index]
        interior || Extents.intersects(target.cap, tree.caps[k]) || continue
        index = tree.order[k]
        _cell_matches(pred, target, grid, cellindex(grid, index), frontier, scratch, true) &&
            push!(out, index)
    end
    return nothing
end

function _leaf_scan!(out, node, grid, pred, target, frontier, scratch, interior::Bool)
    for (index, cap) in STI.child_indices_extents(node)
        interior || Extents.intersects(target.cap, cap) || continue
        _cell_matches(pred, target, grid, cellindex(grid, index), frontier, scratch, true) &&
            push!(out, index)
    end
    return nothing
end

# One cell's answer; `capped` says the cell already passed the target-cap
# prune. With an arc tree, the cell narrows the frontier to its own cap first:
# no arc left means the cell is wholly inside or outside and its centroid
# decides, and otherwise only the surviving arcs feed the sandwich.
function _cell_matches(pred::DE9IM.DE9IMPredicate, target::GeometryTarget, grid, c,
        frontier::Vector{Int}, scratch::Vector{Int}, capped::Bool)
    edges = target.edges
    if edges === nothing
        capped || Extents.intersects(target.cap, cell_cap(grid, c)) || return false
        return _matches(pred, target, grid, c)
    end
    # The cap prune goes first: systems with an analytical cap answer it
    # without the boundary the frontier step needs.
    capped || Extents.intersects(target.cap, cell_cap(grid, c)) || return false
    ring = closed_ring(cell_boundary(grid, c))
    cap = points_cap(ring)
    Float64(cap.radius) < Float64(pi) / 2 || return _matches(pred, target, grid, c)
    _narrow!(scratch, edges, frontier, cap)
    centroid = cell_centroid(grid, c)
    if isempty(scratch)
        return _accepts_interior(pred) &&
               GO.relate_predicate(target.prepared, GO.pred_intersects(), centroid)
    end
    return _matches_near(pred, target, grid, c, centroid, ring, scratch)
end

function _cell_matches(pred::CentroidCovered, target::QueryTarget, grid, c,
        ::Vector{Int}, ::Vector{Int}, capped::Bool)
    capped || Extents.intersects(target.cap, cell_cap(grid, c)) || return false
    return _matches(pred, target, grid, c)
end

function _cell_matches(pred::DE9IM.DE9IMPredicate, target::CapTarget, grid, c,
        ::Vector{Int}, ::Vector{Int}, capped::Bool)
    capped || Extents.intersects(target.cap, cell_cap(grid, c)) || return false
    return _matches(pred, target, grid, c)
end

# The exact tests for a cell the boundary does reach, with the ring in hand.
function _matches_near(::DE9IM.Intersects, target::GeometryTarget, grid, c, centroid, ring,
        frontier::Vector{Int})
    GO.relate_predicate(target.prepared, GO.pred_contains(), centroid) && return true
    verdict = _sandwich(target.edges, frontier, centroid, ring)
    verdict == 0 || return verdict > 0
    return GO.relate_predicate(target.prepared, GO.pred_intersects(),
        GI.Polygon([GI.LinearRing(ring)]))
end

function _matches_near(::DE9IM.Within, target::GeometryTarget, grid, c, centroid, ring,
        ::Vector{Int})
    GO.relate_predicate(target.prepared, GO.pred_intersects(), centroid) || return false
    return GO.relate_predicate(target.prepared, GO.pred_contains(),
        GI.Polygon([GI.LinearRing(ring)]))
end

_matches_near(pred::DE9IM.DE9IMPredicate, target::GeometryTarget, grid, c, centroid, ring,
        ::Vector{Int}) = GO.relate_predicate(target.prepared, _converse_predicate(pred),
    GI.Polygon([GI.LinearRing(ring)]))
