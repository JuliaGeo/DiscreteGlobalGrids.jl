# ---------------------------------------------------------------------------
# The query engine
#
# One spherical tree descent, two stages per node:
#
#   1. **Cap prune** (always sound): the node's `SphericalCap` extent — which
#      bounds every cell under it by the covering law — against a cap bounding
#      the target (`_geometry_cap`). Dot products, no polygon touched. Pruning
#      is an optimisation and can only ever over-select, so it never appears in
#      the answer.
#   2. **Exact per-cell tests**: the target is `GO.prepare`d once per query, so
#      a cell's centroid test is a cached point location. For `Intersects` that
#      doubles as a fast accept — a centroid is interior to its cell by
#      contract, so "centroid inside the target" already proves intersection —
#      and only boundary-grazing cells pay the full polygon predicate.
#   3. **The rim sandwich** (`_sandwich`, `Intersects` only): before paying that
#      predicate, decide the cell from its distance to the target's *boundary*.
#      Both arms are proofs, not heuristics — see the note above `_sandwich`.
#
# All predicates run on the unit sphere (`GO.RelateNG(manifold=Spherical())`),
# consuming `UnitSphericalPoint` geometry directly — no lon/lat round trip,
# hence no antimeridian or pole special cases. The flip side is that ring edges
# are **great-circle arcs**: a long east-west edge does not follow its parallel.
# Densify (`GO.segmentize`) geometries whose edges are meant to trace parallels.
# ---------------------------------------------------------------------------

# Leaf-bucket granularity of the query tree. Single-cell leaves pay cursor
# construction and two binary searches per cell; a bucket amortises both across
# its cells and prunes them with one cap first — but too big a bucket makes
# that union cap its own hot spot. 16 measured well across all four systems of
# the previous design on box, antimeridian, polar and quarter-sphere queries.
const QUERY_BUCKET_SIZE = 16

# ===========================================================================
# Targets
# ===========================================================================

abstract type QueryTarget end

"""
    GeometryTarget(geom)

A prepared query target: the geometry lifted to the unit sphere once, a
`RelateNG` preparation of it, a bounding cap for tree pruning, and its boundary
arcs for the rim sandwich.
"""
struct GeometryTarget{P,G,A} <: QueryTarget
    prepared::P
    geom::G
    cap::Cap
    arcs::A
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

function _query_target(geom)
    GI.isgeometry(geom) || GI.trait(geom) !== nothing || throw(ArgumentError(
        "a query target must be a GeoInterface geometry, an Extents.Extent or a " *
        "SphericalCap; got $(typeof(geom))"))
    spherical = _to_unit_sphere(geom)
    prepared = GO.prepare(GO.RelateNG(; manifold=GO.Spherical()), spherical)
    return GeometryTarget(prepared, spherical, _geometry_cap(prepared, spherical),
        _boundary_arcs(spherical))
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

# A lon/lat extent as a query target. A box that spans every longitude is not a
# polygon at all — it is a cap (when it reaches a pole) or an annular band
# (when it does not), so the two are separated here rather than silently
# producing a degenerate ring.
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

# ===========================================================================
# Bounding cap of a target geometry
#
# Two facts make a rigorous cap out of nothing but the vertices and one
# prepared point query:
#
#  1. a cap of radius <= pi/2 is geodesically convex, so the cap centred on the
#     normalised vertex mean with radius = max vertex distance contains every
#     great-circle edge between consecutive vertices — the whole *boundary*.
#  2. the cap's complement is an open cap — connected — containing no boundary
#     point, so it lies entirely inside or entirely outside the geometry, and
#     one point decides which: the cap centre's antipode, the deepest point of
#     that complement. Outside the geometry => the whole complement is outside
#     => the geometry (interior, holes, multi-parts and all) is inside the cap.
#
# A vertex radius past pi/2 breaks the convexity argument and an antipode
# inside the geometry means the region really does reach into the complement;
# both give up and return the full sphere — no pruning, never a wrong prune.
# ===========================================================================

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

# ===========================================================================
# The rim sandwich
#
# Measured on the previous design's descent, the polygon predicate is where an
# intersection query spends ~90% of its time, and MOST of those calls answer
# "no": the cap prune can only bound cells against the whole geometry's cap, so
# every cell in that cap whose centroid falls outside the region pays a full
# polygon-polygon relate (~10 us cold, against ~0.06 us for a prepared centroid
# query). On a long thin geometry the ratio was 6629 wasted calls out of 6730.
#
# What separates those populations is one number: `d`, the distance from the
# cell centroid to the target's BOUNDARY. Sandwich it between two radii of the
# cell:
#
#   * `d > r_out`, the cell's circumradius => DISJOINT. No boundary point lies
#     in the cell, so the cell — connected — is wholly inside or wholly outside
#     the region; its centroid is outside (the centroid test just said so), so
#     all of it is.
#   * `d < r_in`, a lower bound on the cell's inradius => INTERSECTS. The
#     nearest boundary point is then strictly inside the cell, and it belongs
#     to the geometry (a closed region contains its boundary), so they share it.
#
# Only the annulus `r_in <= d <= r_out` — cells the boundary genuinely grazes —
# still needs the exact predicate. `r_out` is the ring's max centroid-to-vertex
# distance, which bounds the whole great-circle polygon (along an arc,
# `t -> cos(dist(centroid, p(t)))` solves `g'' = -k g`, so where positive it is
# concave and dips no lower than at the arc's ends). `r_in` is the distance to
# the nearest great circle *carrying* a ring edge, which is at most the
# distance to the edge itself and therefore a sound lower bound. Each is nudged
# by a relative `SANDWICH_SLACK` in the safe direction so floating-point noise
# cannot turn a grazing cell into a wrong verdict.
# ===========================================================================

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

# Every great-circle edge of the point set the sandwich argument needs, or
# `nothing` for a geometry kind this walk cannot enumerate — in which case the
# sandwich is off and the exact predicate runs as before.
#
# The set must satisfy two properties, and both arms rest on exactly these: its
# points all belong to the geometry (so a point of it inside a cell proves
# intersection), and a cell meeting none of it cannot straddle inside and
# outside (so a cell with an outside centroid that misses it is disjoint). For
# a polygon that set is the topological boundary — ALL rings of ALL parts,
# holes included, which is why the walk goes through `GI.getring` rather than
# `GI.getpoint`: consecutive points across two rings are not an edge, and
# inventing one could produce a wrong accept.
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
    for p in GI.getpoint(chain)
        point = query_point(p)
        previous === nothing ? (first_point = point) :
        push!(arcs, BoundaryArc(previous, point))
        previous = point
    end
    close && previous !== nothing && push!(arcs, BoundaryArc(previous, first_point))
    return arcs
end

# An empty edge set would let the reject arm fire on every cell; there is no
# geometry to be near, but there is also nothing to gain, so switch off.
#
# A NEAR-ANTIPODAL edge switches it off too. A vanishing `nn` has two causes: a
# repeated point, where the arc really is its endpoints and `_arc_cos_distance`
# is exact; and `a == -b`, where the cross product also vanishes but the edge is
# a whole half great circle whose direction is undefined. There the endpoint
# distance UNDER-estimates how close the boundary comes, which is the unsafe
# direction for the reject arm — it could call a cell disjoint that a boundary
# point sits inside. Such a geometry is pathological, and giving up the
# optimisation for it costs only time.
#
# The guard is on `nn` rather than on exact antipodality because the normal is
# already ill-conditioned before it vanishes: `nn` is |a x b|^2, so the
# threshold below covers every edge within 1e-6 rad of a half circle, well past
# the point where the normal's direction carries any accuracy. The dot-product
# test keeps repeated points — where `nn` is just as small but the endpoint
# distance is exact — on the fast path.
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

# `1` the cell provably meets the target, `-1` it provably does not, `0`
# undecided — see the note above for why each arm is a proof.
function _sandwich(arcs::Vector{BoundaryArc}, centroid, ring)
    vertex = 0.0
    for p in ring
        vertex = max(vertex, US.spherical_distance(centroid, USPoint(p[1], p[2], p[3])))
    end
    r_out = vertex * (1 + SANDWICH_SLACK)
    # A quarter sphere is where both the concavity argument for r_out and the
    # sin/cos comparisons below stop holding. No cell of any system in scope
    # comes near it, and giving up is free.
    r_out < Float64(pi) / 2 || return 0
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
    best = -1.0
    for arc in arcs
        value = _arc_cos_distance(arc, centroid)
        value > cos_in && return 1          # d < r_in: a boundary point is inside the cell
        best = max(best, value)
    end
    return best < cos_out ? -1 : 0          # d > r_out: no boundary point is near
end

# ===========================================================================
# Proving a whole subtree outside the target
# ===========================================================================

# The reject arm of the sandwich again, but asked of a NODE EXTENT rather than
# of one cell's ring, so that the answer covers the node's entire subtree.
#
# A cap that no boundary arc reaches is a connected set containing no point of
# the target's boundary, so it lies wholly in the target's interior or wholly in
# its exterior, and its centre says which. Exterior means no cell of the subtree
# can meet the target — the covering law of `node_extent` carries the conclusion
# from the cap to every descendant — so the subtree can be pruned.
#
# This is the polygon analogue of the wide-cap complement move. A target's own
# bounding cap is the cheap prune, and it is the WHOLE SPHERE whenever the
# target's antipode is interior to it (`_geometry_cap`); a target wider than a
# hemisphere therefore prunes nothing by cap alone and the traversal would visit
# every cell of its deepest level. It also pays on ordinary targets: California's
# bounding cap is a 6.6-degree disc and the state fills a small part of it, and
# everything in between is pruned here instead of being descended into.
#
# `nothing` arcs (an empty or near-antipodal boundary — see `_finish_arcs!`)
# means no proof is available, and no proof means no prune.
_subtree_outside(::QueryTarget, extent) = false

function _subtree_outside(target::GeometryTarget, extent)
    arcs = target.arcs
    arcs === nothing && return false
    # `cos` is decreasing on `[0, pi]`, so "farther than the radius" is "cosine
    # below the radius's cosine". The slack widens the radius, which is the
    # conservative direction: a cap that only just clears the boundary is not
    # trusted. Clamping at `pi` keeps a full-sphere extent — what
    # `AuthalicSystem` returns for a node its warp bound cannot contain — from
    # wrapping past the antipode and pruning everything.
    threshold = cos(min(Float64(pi), Float64(extent.radius) * (1 + SANDWICH_SLACK)))
    for arc in arcs
        _arc_cos_distance(arc, extent.point) >= threshold && return false
    end
    return !GO.relate_predicate(target.prepared, GO.pred_contains(), extent.point)
end

# ===========================================================================
# The per-cell exact tests
# ===========================================================================

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
# rim sandwich, then the exact polygon predicate. The other predicates go
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

_matches(pred::DE9IM.DE9IMPredicate, target::GeometryTarget, grid, c) =
    GO.relate_predicate(target.prepared, _converse_predicate(pred), cell_polygon(grid, c))

# --- cap targets, answered exactly and without polygonising the cap ---------

_matches(::DE9IM.Intersects, target::CapTarget, grid, c) =
    _cell_meets_cap(target.cap, grid, c, false)

# The exact cell/cap overlap test, at any radius. `strict` asks about the OPEN
# cap (distance strictly below the radius) instead of the closed one, which is
# what the complement of a closed cap actually is — see `Within` below.
#
# Exact because: if the cap's centre lies in the cell the two meet, and
# otherwise any point of the intersection can be joined to the centre by a path
# inside the cap, which must cross the cell's boundary on the way out. So a cell
# meeting the cap without any boundary point in the cap is impossible. That
# argument needs the cap to be connected along those paths, i.e. convex, which
# holds up to radius pi/2; a wider cap is only ever asked about here through its
# complement, which is narrower than pi/2 by construction.
function _cell_meets_cap(cap, grid, c, strict::Bool)
    ring, n = open_ring(cell_boundary(grid, c))
    radius = Float64(cap.radius)
    holds(p) = strict ? US.spherical_distance(cap.point, p) < radius :
               cap_contains(cap, p)
    holds(cell_centroid(grid, c)) && return true
    for i in 1:n
        holds(ring[i]) && return true
    end
    # The cap's centre inside the cell, or the cell's boundary reaching into
    # the cap: between them these cover every remaining way to overlap.
    #
    # An undecidable containment is taken as a hit rather than a miss. By this
    # point the cap holds no vertex and (below) may hold no boundary point at
    # all, so the cap is either wholly inside the cell or wholly outside it —
    # and reading "undecidable" as "outside" would drop a cap that sits deep
    # inside a cell, which is precisely the configuration that makes the
    # containment test give up. The `query` docstring states the trade.
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

# A cell is inside a cap of radius <= pi/2 exactly when its vertices are: such a
# cap is convex, so it carries the great-circle edge between any two points it
# holds. A wider cap is not convex and that argument collapses — but its
# complement is, being the OPEN cap of radius pi - r < pi/2 about the antipode,
# and
#
#     cell within C(p, r)  <=>  cell does not meet C(-p, pi - r)
#
# which `_cell_meets_cap` decides with no restriction on r at all. The
# complement is open, so the test there is strict: a cell tangent to the rim of
# `cap` from the inside is within it, and must not read as meeting the
# complement. ("Everything except a region" caps are ordinary enough — a
# coverage of the whole world bar one country — so this is not a corner.)
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

Every cell satisfying the spatial predicate, as a sorted `Vector` of typed ids —
see the interface docstring for the contract.

Implemented predicates: `Intersects`, `Disjoint`, `Contains`, `Within`,
`Covers`, `CoveredBy`, `Touches`, `Overlaps`, `Equals`. `Crosses` throws, since
inverting it is not a matter of naming its converse. Targets: a GeoInterface
geometry, an `Extents.Extent` in lon/lat degrees, or a
`GO.UnitSpherical.SphericalCap`.

`Disjoint` is the complement of `Intersects` and is therefore the one predicate
that cannot prune: it visits every cell of the grid by construction.

A `SphericalCap` target is answered from the cap itself rather than by
polygonising it, at any radius — `Within` decides a cap wider than a hemisphere
through its complement, which is narrower than one. Two consequences worth
knowing: where the containment of the cap's *centre* in a cell is numerically
undecidable, `Intersects` counts it as a hit, since at that point the cap is
either wholly inside the cell or wholly outside and reading "undecidable" as
"outside" would drop a cap sitting deep inside a cell; and only `Intersects`,
`Disjoint` and `Within` are implemented for caps, the rest throwing.
"""
function query(grid::AbstractGrid, pred::DE9IM.DE9IMPredicate)
    # The grid is asked its size first, deliberately: a grid that implements
    # nothing must bounce off the interface, not off the target parser. The
    # predicate and target are then validated BEFORE the empty short-circuit,
    # so that an unusable query is an error whether or not the grid has cells.
    n = ncells(grid)
    _check_predicate(pred)
    target = _query_target(Base.parent(pred))
    n == 0 && return _empty_ids(grid)
    return _run_query(grid, pred, target)
end

# The id type of a grid's cells, without asking for a cell.
_id_type(grid::AbstractGrid) = (sys = system(grid);
sys === nothing ? typeof(cellindex(grid, 1)) : cellindextype(sys))

# A standalone grid names its id type only through a cell, and an empty grid has
# none to ask — so the empty answer for one is abstractly typed. Unavoidable
# without adding an id-type primitive to the interface, and harmless: the vector
# is empty, so nothing is ever stored in it.
_empty_ids(grid::AbstractGrid) = (sys = system(grid);
sys === nothing ? AbstractCellIndex[] : cellindextype(sys)[])

function query(sys::AbstractHierarchicalGridSystem, pred::DE9IM.DE9IMPredicate;
        level::Integer)
    # `level` the keyword shadows `level` the function throughout this body.
    return query(levelgrid(sys, level), pred)
end

function _run_query(grid::AbstractGrid, pred::DE9IM.DE9IMPredicate, target::QueryTarget)
    positions = _query_positions(grid, pred, target)
    out = Vector{_id_type(grid)}(undef, length(positions))
    for (k, i) in enumerate(positions)
        out[k] = cellindex(grid, i)
    end
    # Positions ascend in canonical id order for every grid a system produces,
    # but a standalone grid's own order need not, and the contract is sorted
    # BY ID.
    return sort!(out)
end

# The complement, over every cell of the grid. Sound by construction and the
# only honest way to answer it: a cell that no tree node's extent reaches is
# still disjoint from the target, so there is nothing to prune away.
function _run_query(grid::AbstractGrid, ::DE9IM.Disjoint, target::QueryTarget)
    hits = _query_positions(grid, DE9IM.Intersects(nothing), target)
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

function _query_positions(grid::AbstractGrid, pred::DE9IM.DE9IMPredicate,
        target::QueryTarget)
    ncells(grid) == 0 && return Int[]
    out = Int[]
    _descend!(out, _query_tree(grid), grid, pred, target)
    return sort!(out)
end

_query_tree(grid::AbstractGrid) = system(grid) === nothing ?
                                  PositionTreeNode(PositionTree(grid), 1) :
                                  HierarchicalGridCursor(grid;
    bucket_size=max(QUERY_BUCKET_SIZE, _grid_bucket_size(grid)))

function _descend!(out, node, grid, pred, target)
    extent = STI.node_extent(node)
    intersects_cap(target.cap, extent) || return nothing
    if STI.isleaf(node)
        # A node whose whole extent sits inside the target's cap has no per-cell
        # cap prune left to win — every cell cap would pass — so it skips
        # straight to the exact tests.
        interior = US.spherical_distance(target.cap.point, extent.point) +
                   extent.radius <= target.cap.radius
        _leaf_scan!(out, node, grid, pred, target, interior)
        return nothing
    end
    for child in STI.getchild(node)
        _descend!(out, child, grid, pred, target)
    end
    return nothing
end

# The leaf scan, once per tree kind: the hierarchical cursor derives a cell's
# cap from its boundary, the position tree has it stored.
function _leaf_scan!(out, node::HierarchicalGridCursor, grid, pred, target, interior::Bool)
    for index in node_indices(node)
        c = cellindex(grid, index)
        interior || intersects_cap(target.cap, cell_cap(grid, c)) || continue
        _matches(pred, target, grid, c) && push!(out, index)
    end
    return nothing
end

function _leaf_scan!(out, node::PositionTreeNode, grid, pred, target, interior::Bool)
    tree = node.tree
    for k in tree.node_first[node.index]:tree.node_last[node.index]
        interior || intersects_cap(target.cap, tree.caps[k]) || continue
        index = tree.order[k]
        _matches(pred, target, grid, cellindex(grid, index)) && push!(out, index)
    end
    return nothing
end
