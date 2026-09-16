# Spherical query descent: prune by conservative caps, then apply exact per-cell
# predicates. `Intersects` also uses centroid acceptance and a boundary-distance
# sandwich. Input ring edges are great-circle arcs; densify intended parallels
# (`GO.segmentize`).

# Leaf bucket size balancing cursor overhead against cap selectivity.
const QUERY_BUCKET_SIZE = 16

# ===========================================================================
# Targets
# ===========================================================================

abstract type QueryTarget end

"""
    GeometryTarget(geom)

A geometry prepared once for querying.

- `prepared`: the `RelateNG` form of `geom`; its spherical kernel reads 2-D
  coordinates as lon/lat degrees and 3-D coordinates as unit-sphere xyz.
- `cap`: bounds the geometry for tree pruning.
- `arcs`: boundary arcs for the border sandwich, in the leaf order of `edges`,
  an [`ArcTree`](@ref) over them for the descent's frontier.
"""
struct GeometryTarget{P,A,E} <: QueryTarget
    prepared::P
    cap::Cap
    arcs::A
    edges::E
end

function GeometryTarget(prepared, cap::Cap, arcs)
    edges = _arc_tree(arcs)
    return GeometryTarget(prepared, cap, edges === nothing ? arcs : edges.data, edges)
end

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

const TargetTrait = Union{GI.PointTrait,GI.MultiPointTrait,GI.LineStringTrait,
    GI.LinearRingTrait,GI.MultiLineStringTrait,GI.PolygonTrait,GI.MultiPolygonTrait}

function _query_target(geom)
    GI.trait(geom) isa TargetTrait || throw(ArgumentError(
        "a query target must be a point, line, polygon, multi-geometry, " *
        "Extents.Extent or SphericalCap; got $(typeof(geom))"))
    prepared = GO.prepare(GO.RelateNG(; manifold=GO.Spherical()), geom)
    points = USPoint[]
    arcs = _boundary_arcs!(BoundaryArc[], points, GI.trait(geom), geom)
    return GeometryTarget(prepared, _geometry_cap(prepared, geom, points), arcs)
end

_query_target(::Nothing) = throw(ArgumentError(
    "the predicate carries no target; write e.g. `Intersects(geometry)`"))

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
    # Counter-clockwise seen from outside, sampled in lon/lat so parallels stay
    # parallels. A polar edge lifts to one repeated vertex, kept once.
    box = Extents.Extent(X=(x0, x1), Y=(y0 <= -89.999 ? -90.0 : y0, y1 >= 89.999 ? 90.0 : y1))
    outline = GO.segmentize(GO.Planar(), GO.extent_to_polygon(box); max_distance=EXTENT_STEP_DEGREES)
    points = USPoint[]
    for p in GI.getpoint(outline)
        q = unit_point(GI.x(p), GI.y(p))
        (isempty(points) || q != last(points)) && push!(points, q)
    end
    return GI.Polygon([GI.LinearRing(points)])
end

# A radius-`≤ π/2` vertex cap contains all great-circle boundary arcs. Its
# connected, boundary-free complement is wholly inside or outside the geometry,
# decided at the cap antipode. Wider caps or an interior antipode use the full
# sphere.

# `points` are the vertices the arc walk lifted, in `GI.getpoint` order; a
# geometry without arcs (a bare ring) is lifted here.
_geometry_cap(prepared, geom) = _geometry_cap(prepared, geom, USPoint[])

function _geometry_cap(prepared, geom, points::Vector{USPoint})
    GI.trait(geom) isa GI.PointTrait && return SphericalCap(query_point(geom), 0.0)
    isempty(points) && (points = USPoint[query_point(p) for p in GI.getpoint(geom)])
    cap = points_cap(points)
    cap.radius >= Float64(pi) && return full_sphere_cap()
    antipode = USPoint(-cap.point[1], -cap.point[2], -cap.point[3])
    GO.relate_predicate(prepared, GO.pred_intersects(), antipode) && return full_sphere_cap()
    return cap
end

# Border sandwich: `d` (centroid to target boundary) against `r_out` (farthest
# vertex) and `r_in` (nearest edge great circle), both slackened by
# `SANDWICH_SLACK`. `d > r_out` with the centroid outside proves disjoint;
# `d < r_in` puts a boundary point inside the cell; between them, the exact test.

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
# rings, lifting each vertex once into `points` on the way. Return `nothing`
# when the geometry cannot provide a safe boundary set. Rings are walked one
# at a time (`GI.getring`, not `GI.getpoint`): an edge invented between
# consecutive points of two different rings is a wrong accept.
_boundary_arcs(geom) = _boundary_arcs!(BoundaryArc[], USPoint[], GI.trait(geom), geom)

_boundary_arcs!(arcs, points, ::Any, geom) = nothing

function _boundary_arcs!(arcs, points, ::GI.PolygonTrait, geom)
    _rings_arcs!(arcs, points, geom)
    return _finish_arcs!(arcs)
end

function _boundary_arcs!(arcs, points, ::GI.MultiPolygonTrait, geom)
    for poly in GI.getpolygon(geom)
        _rings_arcs!(arcs, points, poly)
    end
    return _finish_arcs!(arcs)
end

function _boundary_arcs!(arcs, points, ::GI.LineStringTrait, geom)
    _chain_arcs!(arcs, points, geom, false)
    return _finish_arcs!(arcs)
end

function _boundary_arcs!(arcs, points, ::GI.MultiLineStringTrait, geom)
    for line in GI.getgeom(geom)
        _chain_arcs!(arcs, points, line, false)
    end
    return _finish_arcs!(arcs)
end

function _boundary_arcs!(arcs, points, ::GI.PointTrait, geom)
    point = query_point(geom)
    push!(points, point)
    push!(arcs, BoundaryArc(point, point))
    return _finish_arcs!(arcs)
end

function _boundary_arcs!(arcs, points, ::GI.MultiPointTrait, geom)
    for p in GI.getpoint(geom)
        point = query_point(p)
        push!(points, point)
        push!(arcs, BoundaryArc(point, point))
    end
    return _finish_arcs!(arcs)
end

function _rings_arcs!(arcs, points, poly)
    for r in GI.getring(poly)
        _chain_arcs!(arcs, points, r, true)
    end
    return arcs
end

# Consecutive-vertex edges of one ring or linestring; `close` adds the edge back
# to the first vertex, which a polygon ring always has whether or not its
# coordinates repeat it (a duplicated closing vertex just makes that edge
# degenerate, and degenerate arcs are handled).
function _chain_arcs!(arcs, points, chain, close::Bool)
    first_point = nothing
    previous = nothing
    count = 0
    for p in GI.getpoint(chain)
        point = query_point(p)
        push!(points, point)
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
    ArcTree

A packed R-tree over the target's boundary arcs, one exact 3-D box per arc.

- Built by `GeometryOps.FlexibleRTrees.RTree`; `data` holds the
  arcs in leaf order, so leaf slot `k` is `data[k]`.
- The descent narrows a frontier of [`ArcNode`](@ref)s box by box.
- An empty frontier under a convex cap proves that cap free of the boundary.
"""
const ArcTree = GO.FlexibleRTrees.RTree{GO.FlexibleRTrees.Unsorted,XYZExtent,
    Vector{BoundaryArc},Base.OneTo{Int}}
"""
    ArcNode

One node of an [`ArcTree`](@ref): the tree, its 0-based level, its index in that
level, and its box. Frontiers are vectors of these.
"""
const ArcNode = GO.FlexibleRTrees.RTreeNode{ArcTree,XYZExtent}

const ARC_TREE_CAPACITY = 8

# Sorting the arcs once up front lets the tree alias its leaf order, so the
# sandwich reads `data` sequentially and never goes through a permutation.
function _arc_tree(arcs::Vector{BoundaryArc})
    isempty(arcs) && return nothing
    extents = XYZExtent[GO.UnitSpherical.spherical_arc_extent(a.a, a.b) for a in arcs]
    order = GO.FlexibleRTrees.loadorder(GO.FlexibleRTrees.STR(), extents, ARC_TREE_CAPACITY)
    return GO.FlexibleRTrees.RTree(GO.FlexibleRTrees.Unsorted(), arcs[order];
        nodecapacity=ARC_TREE_CAPACITY, extents=extents[order])
end

_arc_tree(::Nothing) = nothing

_target_edges(::QueryTarget) = nothing
_target_edges(target::GeometryTarget) = target.edges

# Narrow frontier node `node` to `cap`. A node smaller than the cap is kept
# whole (splitting costs more than it excludes); a leaf survives only if an arc
# can reach the cap by exact distance.
function _edge_frontier!(out::Vector{ArcNode}, node::ArcNode, cap, threshold::Float64)
    Extents.intersects(cap, node.extent) || return out
    if STI.isleaf(node)
        arcs = node.tree.data
        for k in _leaf_range(node)
            if _arc_cos_distance(arcs[k], cap.point) >= threshold
                push!(out, node)
                break
            end
        end
    # The box overshoots the node's chord diameter by up to `sqrt(3)`, so this
    # splits some nodes already smaller than the cap: extra work, never a lost arc.
    elseif _box_half_diagonal(node.extent) > Float64(cap.radius)
        for i in 1:STI.nchild(node)
            _edge_frontier!(out, STI.getchild(node, i), cap, threshold)
        end
    else
        push!(out, node)
    end
    return out
end

function _narrow!(out::Vector{ArcNode}, frontier::Vector{ArcNode}, cap)
    empty!(out)
    # The widened cap serves both the box test and the exact leaf test. The
    # subtracted ulp keeps the widening conservative at radii so small that
    # `cos` of the widened radius rounds back to `cos(r)`.
    wide = SphericalCap(cap.point, min(Float64(pi), Float64(cap.radius) * (1 + SANDWICH_SLACK)))
    threshold = wide.radiuslike - 1e-15
    for node in frontier
        _edge_frontier!(out, node, wide, threshold)
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
function _sandwich(tree::ArcTree, frontier::Vector{ArcNode}, centroid, ring)
    bounds = _sandwich_bounds(centroid, ring)
    bounds === nothing && return 0
    cos_out, cos_in = bounds
    best = -1.0
    arcs = tree.data
    for node in frontier, k in _leaf_range(node)
        value = _arc_cos_distance(arcs[k], centroid)
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

# The sandwich's reject arm on a node extent. A cap no boundary arc reaches is
# wholly interior or exterior, and its centre says which; exterior prunes every
# descendant by `node_extent`'s covering law. `nothing` arcs (empty or
# near-antipodal boundary, see `_finish_arcs!`) prove nothing.
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
    CentroidCovered(region)

A [`query`](@ref) predicate selecting every cell whose centroid,
[`cell_centroid`](@ref), lies on or inside `region`. A centroid on the region's
boundary counts, so the rule is the raster `boundary = :center` rule. Read it as
**cell centroid COVERED BY region**: `Within(region)` selects a subset of
it, and `Intersects(region)` a superset.

`region` is any target `query` accepts: a GeoInterface geometry, an
`Extents.Extent` in lon/lat degrees, or a `GO.UnitSpherical.SphericalCap`.
`Base.parent(pred)` gives it back. The predicate works through `query` and as a
[`Cells`](@ref DiscreteGlobalGrids.CellLookups.Cells) selector:

```julia
query(grid, CentroidCovered(county))    # cells centred in a polygon
A[Cells(CentroidCovered(cap))]          # cells centred in a SphericalCap
```
"""
struct CentroidCovered{T}
    val::T
end
CentroidCovered() = CentroidCovered(nothing)
Base.parent(pred::CentroidCovered) = pred.val

const QueryPredicate = Union{DE9IM.DE9IMPredicate,CentroidCovered}

# On a point, `pred_intersects` is `covers`, so a centroid on the boundary
# counts; `pred_contains` would drop it.
_covers_centroid(target::GeometryTarget, grid, c) =
    GO.relate_predicate(target.prepared, GO.pred_intersects(), cell_centroid(grid, c))
_covers_centroid(target::CapTarget, grid, c) = cap_contains(target.cap, cell_centroid(grid, c))

_matches(::CentroidCovered, target::QueryTarget, grid, c) = _covers_centroid(target, grid, c)

# Predicates a cell wholly inside the target satisfies; the rest need the
# target's boundary, which such a cell does not hold.
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

_ring_polygon(ring) = GI.Polygon([GI.LinearRing(ring)])

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
    return GO.relate_predicate(target.prepared, GO.pred_intersects(), _ring_polygon(ring))
end

# `cell_centroid` lies strictly inside its cell, so a centroid the target does
# not cover rejects `Within` before the polygon test.
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
    query(grid, pred) -> Vector{<:AbstractCellIndex}
    query(sys, pred; level) -> Vector{<:AbstractCellIndex}

Return every cell matching the predicate as a sorted vector of typed ids.

Implemented predicates: `Intersects`, `Disjoint`, `Contains`, `Within`,
`Covers`, `CoveredBy`, `Touches`, `Overlaps`, `Equals`, and the centroid rule
[`CentroidCovered`](@ref). `Crosses` throws, since inverting it is not a matter
of naming its converse. Targets: a GeoInterface geometry, an `Extents.Extent`
in lon/lat degrees, or a `GO.UnitSpherical.SphericalCap`.

`Disjoint` is computed as the full-grid complement of `Intersects` and cannot
prune the output traversal.

Traversal and cap targets:

- The traversal runs over [`treeify`](@ref)`(grid)`; a hierarchical cursor is
  rebuilt at the wider of the query's leaf bucket and the grid's own, and a
  level-less grid (stored cells of mixed levels) gets the engine's index tree.
- Cap targets are answered without polygonisation and support only
  `Intersects`, `Disjoint`, `Within`, and `CentroidCovered`; `Within` uses
  the complement past a hemisphere, and `Intersects` accepts undecidable
  cap-centre containment.
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

function query(sys::AbstractHierarchicalGridSystem, pred::QueryPredicate;
        level::Integer)
    # `level` the keyword shadows `level` the function throughout this body.
    return query(levelgrid(sys, level), pred)
end

function query(grid::AbstractGrid, pred::CentroidCovered)
    n = ncells(grid)
    pred.val === nothing && throw(ArgumentError(
        "CentroidCovered needs a region; write e.g. `CentroidCovered(geometry)`"))
    target = _query_target(pred.val)
    n == 0 && return _empty_ids(grid)
    return _run_query(grid, pred, target)
end

function _run_query(grid::AbstractGrid, pred::QueryPredicate, target::QueryTarget)
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
    edges = _target_edges(target)
    frontier = edges === nothing ? ArcNode[] : ArcNode[_root_node(edges)]
    _descend!(out, _query_tree(grid), grid, pred, target, frontier, ArcNode[],
        Vector{ArcNode}[], 1)
    return sort!(out)
end

# A hierarchical cursor is rebuilt at the wider bucket for cap selectivity; a
# level-less grid (stored cells of mixed levels) fits only the index tree.
function _query_tree(grid::AbstractGrid)
    (system(grid) === nothing || level(grid) === nothing) &&
        return IndexTreeNode(IndexTree(grid))
    tree = treeify(grid)
    tree isa HierarchicalGridCursor || return tree
    return HierarchicalGridCursor(grid;
        bucket_size=max(QUERY_BUCKET_SIZE, _grid_bucket_size(grid)))
end

# `frontier`: arc-tree nodes whose arcs may still reach this node. Only a convex
# cap narrows it: its intersection with the convex ancestor caps is connected
# and holds every descendant (covering law), so an empty frontier lets one real
# descendant centroid decide the node; the cap centre may lie outside and is
# never the witness. `buffers[depth]` never aliases a live ancestor's frontier.
function _descend!(out, node, grid, pred, target, frontier::Vector{ArcNode},
        scratch::Vector{ArcNode}, buffers::Vector{Vector{ArcNode}}, depth::Int)
    extent = STI.node_extent(node)
    Extents.intersects(target.cap, extent) || return nothing
    edges = _target_edges(target)
    if edges !== nothing && Float64(extent.radius) < Float64(pi) / 2
        frontier = _narrow!(_depth_buffer(buffers, depth), frontier, extent)
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
        _descend!(out, child, grid, pred, target, frontier, scratch, buffers, depth + 1)
    end
    return nothing
end

# Grown on demand: a descent touches as many depths as the tree is deep.
function _depth_buffer(buffers::Vector{Vector{ArcNode}}, depth::Int)
    while length(buffers) < depth
        push!(buffers, ArcNode[])
    end
    return buffers[depth]
end

# The first cell under a node, or `0` for an empty one.
_witness_index(node::HierarchicalGridCursor) =
    (indices = node_indices(node); isempty(indices) ? 0 : first(indices))

# Any other tree: the first leaf index the search reaches.
function _witness_index(node)
    hit = STI.depth_first_search(_first_hit, Returns(true), node)
    return hit isa GO.LoopStateMachine.Action ? hit.x : 0
end

_first_hit(index::Int) = GO.LoopStateMachine.Action(:full_return, index)

# An index-tree node's leaf run spans its whole subtree, leaf or not.
_node_window(node::IndexTreeNode) = @view node.tree.rtree.indices[_leaf_range(node.node)]

_witness_index(node::IndexTreeNode) =
    (window = _node_window(node); isempty(window) ? 0 : first(window))

_emit!(out, node::HierarchicalGridCursor) = (append!(out, node_indices(node)); nothing)

_emit!(out, node::IndexTreeNode) = (append!(out, _node_window(node)); nothing)

_emit!(out, node) = (STI.depth_first_search(Base.Fix1(push!, out), Returns(true), node); nothing)

# One leaf scan per tree kind: the cursor leaves the cap to `_cell_matches`;
# any other tree lists its leaf extents.
function _leaf_scan!(out, node::HierarchicalGridCursor, grid, pred, target, frontier, scratch, interior::Bool)
    for index in node_indices(node)
        _cell_matches(pred, target, grid, cellindex(grid, index), frontier, scratch, interior) &&
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

# One cell's answer; `capped` says the target-cap prune already passed. The
# frontier narrows to `cell_cap` (sound for any convex cap, hence the `pi/2`
# guard); an empty frontier lets the centroid decide, else the survivors feed
# the sandwich.
function _cell_matches(pred::DE9IM.DE9IMPredicate, target::GeometryTarget, grid, c,
        frontier::Vector{ArcNode}, scratch::Vector{ArcNode}, capped::Bool)
    # One cap serves both the prune and the narrowing, so a system with an
    # analytical `cell_cap` decides most cells without a boundary at all.
    cap = cell_cap(grid, c)
    capped || Extents.intersects(target.cap, cap) || return false
    edges = target.edges
    edges === nothing && return _matches(pred, target, grid, c)
    Float64(cap.radius) < Float64(pi) / 2 || return _matches(pred, target, grid, c)
    _narrow!(scratch, frontier, cap)
    centroid = cell_centroid(grid, c)
    if isempty(scratch)
        return _accepts_interior(pred) &&
               GO.relate_predicate(target.prepared, GO.pred_intersects(), centroid)
    end
    # Most cells leave here with an empty frontier, so the ring is built only
    # for the few the boundary actually reaches.
    ring = closed_ring(cell_boundary(grid, c))
    return _matches_near(pred, target, grid, c, centroid, ring, scratch)
end

# Every other pair — `CentroidCovered`, and any predicate against a cap — reads
# a centroid or the cap itself, so the frontier has nothing to narrow.
function _cell_matches(pred, target::QueryTarget, grid, c,
        ::Vector{ArcNode}, ::Vector{ArcNode}, capped)
    capped || Extents.intersects(target.cap, cell_cap(grid, c)) || return false
    return _matches(pred, target, grid, c)
end

# The exact tests for a cell the boundary does reach, with the ring in hand.
function _matches_near(::DE9IM.Intersects, target::GeometryTarget, grid, c, centroid, ring,
        frontier::Vector{ArcNode})
    GO.relate_predicate(target.prepared, GO.pred_contains(), centroid) && return true
    verdict = _sandwich(target.edges, frontier, centroid, ring)
    verdict == 0 || return verdict > 0
    return GO.relate_predicate(target.prepared, GO.pred_intersects(), _ring_polygon(ring))
end

function _matches_near(::DE9IM.Within, target::GeometryTarget, grid, c, centroid, ring,
        ::Vector{ArcNode})
    GO.relate_predicate(target.prepared, GO.pred_intersects(), centroid) || return false
    return GO.relate_predicate(target.prepared, GO.pred_contains(), _ring_polygon(ring))
end

_matches_near(pred::DE9IM.DE9IMPredicate, target::GeometryTarget, grid, c, centroid, ring,
        ::Vector{ArcNode}) = GO.relate_predicate(target.prepared, _converse_predicate(pred),
    _ring_polygon(ring))
