# ---------------------------------------------------------------------------
# Lookup-level operations: the neighbor halo table, `stencil`, and `zonal`
#
# These used to be HEALPix-only functions in `HealpixLookups`. Each one is a
# composition of things the kernel already answers generically — neighbors
# come from `cell_neighbors`, positions from `cell_position`, spatial queries
# from the `SpatialTreeInterface` tree over `DGGSPartialGrid` — so the
# functions live here, are
# exported once from `DiscreteGlobalGrids`, and every system whose kernel is
# wired gets them for free. `HealpixLookups` re-exports the same bindings (as
# it always exported these names), so a `using` of both namespaces can never
# make them ambiguous — the `treeify` re-export pattern.
#
# Two lookup accessors anchor the group: a lookup knows which system and level
# its ids live at, but each concrete type spells that differently (`level`,
# `resolution`), so `dggs_system` / `dggs_level` are the one place generic
# code asks. Wired next to `DGGSPartialGrid(l)` in the per-system kernel
# files; an unwired lookup type is a MethodError, not a silent default.
# ---------------------------------------------------------------------------

"""
    dggs_system(l::AbstractDGGSLookup) -> AbstractDGGS

The grid-system singleton whose cells the lookup holds, e.g. `HEALPixDGGS()`
for a `HealpixLookup`. Wired per lookup type next to `DGGSPartialGrid(l)` in
the system's kernel file.
"""
function dggs_system end

"""
    dggs_level(l::AbstractDGGSLookup) -> Int

The refinement level the lookup's ids live at — the field the concrete types
call `level` (HEALPix) or `resolution` (H3, IGEO7, A5). Wired per lookup type
next to [`dggs_system`](@ref).
"""
function dggs_level end

# The unique dimension of `A` holding a DGGS lookup — how `stencil` and
# `zonal` find their cell axis without naming any system's dim type (`Cells`,
# `H3Cells`, ...).
function _dggs_dim(A::DD.AbstractDimArray)
    dims = DD.dims(A)
    hits = [i for i in 1:length(dims) if DD.val(dims[i]) isa AbstractDGGSLookup]
    length(hits) == 1 || throw(ArgumentError(
        "expected exactly one dimension holding an AbstractDGGSLookup, found $(length(hits))"))
    return dims[hits[1]]
end

# --------------------------------------------------------------------------
# Neighbor halo table
# --------------------------------------------------------------------------

"""
    neighbor_indices(l::AbstractDGGSLookup) -> Vector{SmallVector{max_neighbors(system),Int}}

For each stored cell, the positions (into the lookup) of its edge neighbors —
[`cell_neighbors`](@ref) mapped through [`cell_position`](@ref) — ascending by
neighbor id, `0` where the neighbor cell is not stored (coverage boundary).
Computed once; this is the static "halo table" a [`stencil`](@ref) sweep
resolves values through (the DLWP-HPX pattern), so pass it back as `nbidx` to
amortize it across many stencils.
"""
function neighbor_indices(l::AbstractDGGSLookup)
    system = dggs_system(l)
    level = dggs_level(l)
    ids = DD.parent(l)
    return map(ids) do id
        _neighbor_positions(ids, cell_neighbors(system, level, id))
    end
end

function _neighbor_positions(ids, neighbors::SmallVector{N}) where {N}
    out = SmallVector{N,Int}()
    for neighbor in neighbors
        out = SmallCollections.push(out, something(cell_position(ids, neighbor), 0))
    end
    return out
end

# --------------------------------------------------------------------------
# Stencil sweep
# --------------------------------------------------------------------------

"""
    stencil(f, A::DD.AbstractDimArray; nbidx=nothing)

Apply `f(center_value, neighbor_values::SmallVector)` over every stored cell
of the 1-D DGGS-dimensioned array `A`. Neighbors outside the stored coverage
are simply absent from the values container (partial-coverage semantics:
reductions skip, like `nanmean`). Pass a precomputed
[`neighbor_indices`](@ref)`(lookup)` as `nbidx` to amortize the halo table
across many stencils; with it in hand the sweep itself does not allocate
beyond the output array.
"""
function stencil(f, A::DD.AbstractDimArray; nbidx=nothing)
    l = DD.val(_dggs_dim(A))
    ndims(A) == 1 || throw(ArgumentError(
        "stencil expects a 1-D array over the DGGS cell dimension, got $(ndims(A)) dimensions"))
    nbi = nbidx === nothing ? neighbor_indices(l) : nbidx
    data = parent(A)
    length(nbi) == length(data) ||
        throw(DimensionMismatch("neighbor index table and array must have the same length"))
    return DD.rebuild(A; data=_stencil_sweep(f, data, nbi))
end

# Function barrier: `N` enters as a type parameter here, so the loop below is
# fully inferred and the per-cell values container never touches the heap.
function _stencil_sweep(f, data::AbstractVector,
        nbi::AbstractVector{SmallVector{N,Int}}) where {N}
    T = eltype(data)
    out = similar(data, Base.promote_op(f, T, SmallVector{N,T}))
    @inbounds for i in eachindex(data, nbi)
        values = SmallVector{N,T}()
        for j in nbi[i]
            j > 0 && (values = SmallCollections.push(values, data[j]))
        end
        out[i] = f(data[i], values)
    end
    return out
end

# --------------------------------------------------------------------------
# Zonal statistics
# --------------------------------------------------------------------------

"""
    zonal(f, A::DD.AbstractDimArray; of, boundary=:center, skipmissing=true)

Zonal statistics over the DGGS cell dimension of `A`. `of` is a geometry,
feature(collection), or vector thereof; `boundary=:center` selects stored
cells whose center lies in the geometry's interior, any other value selects
cells whose spherical cell polygon intersects it. Returns one value per
geometry; `missing` where no stored cell matches.

On the equal-area systems (HEALPix, IGEO7) `zonal(mean, ...)` is the true
(unweighted) areal mean — no latitude weighting.

The cell query is the spherical tree descent below (`_query_positions` /
`_tree_query`): all predicates are evaluated on the unit sphere by the
spherical `GeometryOps.RelateNG` engine, so geometries crossing the
antimeridian or enclosing a pole need no special handling. The flip side is
that ring edges are **great-circle arcs**: a long east–west edge does not
follow its parallel. Densify (`GO.segmentize`) geometries whose edges are
meant to trace parallels. 2-D coordinates are lon/lat degrees; 3-D
coordinates are taken as unit-sphere points as-is.
"""
function zonal(f, A::DD.AbstractDimArray; of, boundary::Symbol=:center, skipmissing::Bool=true)
    dim = _dggs_dim(A)
    l = DD.val(dim)
    geoms = _geometries(of)
    mode = boundary === :center ? :center : :touches
    map(geoms) do g
        idx = _query_positions(l, g, mode)
        isempty(idx) && return missing
        sub = A[DD.rebuild(dim, idx)]
        vals = skipmissing ? Base.skipmissing(sub) : sub
        isempty(vals) ? missing : f(vals)
    end
end

_geometries(of) = GI.isgeometry(of) ? [of] :
    GI.trait(of) isa GI.AbstractFeatureCollectionTrait ? [GI.geometry(f) for f in GI.getfeature(of)] :
    GI.trait(of) isa GI.AbstractFeatureTrait ? [GI.geometry(of)] :
    of isa AbstractVector ? map(g -> GI.trait(g) isa GI.AbstractFeatureTrait ? GI.geometry(g) : g, of) :
    throw(ArgumentError("cannot extract geometries from $(typeof(of))"))

# Positions (into the lookup) of the stored cells matching `geom` under
# `mode` (`:center` / `:touches`), ascending. One spherical tree descent for
# every lookup type — the per-lookup dispatch point stays so a system with a
# genuinely better native query could still override it.
_query_positions(l::AbstractDGGSLookup, geom, mode::Symbol) = _tree_query(l, geom, mode)

#=
## Spherical tree query

The cell query behind the `Touching` selectors and `zonal`, over the
`SpatialTreeInterface` tree the lookup already exposes
(`treeify(DGGSPartialGrid(l))`, `core/generic_cursor.jl`). One descent, two
stages per node:

1. **Cap prune** (always sound): the node's `STI.node_extent` — a
   `SphericalCap` bounding every stored leaf under it, by the cap contracts
   the kernel test suites measure — against a cap bounding the query geometry
   (`_geometry_cap` below). Dot-product arithmetic, no polygon touched.
2. **Exact leaf tests**: the geometry is `GO.prepare`d once per query, so a
   leaf's center test is a cached point-location. That doubles as a fast
   accept for `:touches` — every wired system's `cell_center` is interior to
   its cell, so "center inside the geometry" already proves intersection and
   only boundary-grazing cells pay the full polygon predicate. Leaves come
   in `QUERY_BUCKET_SIZE` buckets (each with its own cap, pruned per cell
   unless the bucket's cap already sits inside the geometry's).
3. **The rim sandwich** (`:touches` only, `_sandwich` below): before paying
   the polygon predicate, decide the cell from its distance to the query
   geometry's *boundary*. Both arms are proofs, not heuristics.

All predicates run on the unit sphere (`GO.RelateNG(manifold=Spherical())`),
consuming the kernel's `UnitSphericalPoint` geometry directly — no lon/lat
round trip, hence no antimeridian or pole special cases. `GO.prepare`
validates the query geometry on the spherical manifold by default (a ring
whose edges cross as great-circle arcs would invert containment); repair such
geometry with `GO.CrossingEdgeSplit` as the error suggests.

### Design note: why there is no polygon prune / bulk-accept stage

Where a parent geographically contains its descendants
([`has_congruent_geometry`](@ref) — HEALPix), classifying internal nodes
against the exact subtree outline ([`subtree_polygon_unitsphere`](@ref)) —
`pred_disjoint` to prune, `pred_covers` to accept every stored leaf at once —
is sound, and is what the old planar HEALPix descent did. Measured against
this file's engine, it loses at every scale tried (levels 6 and 8, box and
quarter-sphere queries, both modes, 2–7× slower): the prepared point query
makes an interior leaf cost well under a microsecond, so bulk-accept saves
almost nothing, while every node the geometry's *boundary* crosses pays a
polygon predicate that scales with the densified outline — and a failing
`pred_covers` on a 4·2^Δ-vertex outline is the engine's worst case, not its
best. The outline API stays wired and tested (HEALPix), so a traversal with
different economics — a dual-tree sweep, a cheaper prepared-B engine — can
pick the classification back up without re-deriving the geometry.
=#

# Leaf-bucket granularity of the query tree: single-cell leaves pay cursor
# construction and two binary searches per stored cell, a bucket amortizes
# both across its cells and prunes them with one union cap first — but too
# big a bucket makes that union cap (O(bucket) boundary calls on the systems
# without exact subtree caps) its own hot spot. 16 measured well across all
# four systems on box, antimeridian, polar and quarter-sphere queries;
# 64 was already pathological for an H3 globe lookup (49-cell buckets of
# native-call boundaries).
const QUERY_BUCKET_SIZE = 16

function _tree_query(l::AbstractDGGSLookup, geom, mode::Symbol)
    system = dggs_system(l)
    leaf = dggs_level(l)
    ids = DD.parent(l)
    tree = treeify(DGGSPartialGrid(l; bucket_size=QUERY_BUCKET_SIZE))
    prep = GO.prepare(GO.RelateNG(; manifold=GO.Spherical()), geom)
    cap = _geometry_cap(prep, geom)
    # Only `:touches` ever reaches the polygon predicate, so only `:touches`
    # needs the boundary arcs the sandwich measures against.
    arcs = mode === :center ? nothing : _boundary_arcs(geom)
    out = Int[]
    _tree_query!(out, tree, ids, system, leaf, prep, cap, arcs, mode)
    # The range-backed cursors already yield ascending positions; only a
    # selection-cursor system without descendant ranges (A5) can interleave
    # across children, and a near-sorted sort! is cheap.
    return sort!(out)
end

function _tree_query!(out, node, ids, system, leaf, prep, cap, arcs, mode)
    extent = STI.node_extent(node)
    intersects_cap(cap, extent) || return nothing
    if STI.isleaf(node)
        indices = node_indices(node)
        if length(indices) == 1
            # A single-cell leaf's node extent IS its cell cap — already
            # tested on entry, so go straight to the exact test.
            index = first(indices)
            _leaf_hit(prep, arcs, system, leaf, ids[index], mode) && push!(out, index)
        else
            # A bucket whose cap sits entirely inside the geometry's cap has
            # no per-cell cap prune left to win — every cell cap would pass —
            # so skip straight to the exact tests.
            interior = GO.UnitSpherical.spherical_distance(cap.point, extent.point) +
                       extent.radius <= cap.radius
            for index in indices
                interior || intersects_cap(cap, cell_cap(system, leaf, ids[index])) ||
                    continue
                _leaf_hit(prep, arcs, system, leaf, ids[index], mode) && push!(out, index)
            end
        end
        return nothing
    end
    for child in STI.getchild(node)
        _tree_query!(out, child, ids, system, leaf, prep, cap, arcs, mode)
    end
    return nothing
end

# The exact per-cell test. `:center` is "cell center in the geometry's
# interior" (`pred_contains`, the DE-9IM sense the old planar
# `GO.contains(geom, center)` had). For `:touches` the same cached point
# query is a fast accept — the center is interior to its cell in every wired
# system, so center ∈ geom ⟹ cell ∩ geom ≠ ∅ — and cells whose center lands
# outside go to the rim sandwich first, paying the polygon–polygon
# `pred_intersects` only when neither of its arms proves the answer.
function _leaf_hit(prep, arcs, system, leaf, id, mode::Symbol)
    center = cell_center(system, leaf, id)
    GO.relate_predicate(prep, GO.pred_contains(), center) && return true
    mode === :center && return false
    if arcs !== nothing
        # One boundary fetch feeds both the sandwich and the polygon it falls
        # through to: `cell_polygon_unitsphere` is exactly this ring wrapped,
        # since every system wires `cell_boundary` and none overrides the
        # polygon itself — pinned per system in `test_tree_queries.jl` so that
        # a future override cannot silently make this a different test.
        ring = cell_boundary(system, leaf, id; closed=true)
        verdict = _sandwich(arcs, center, ring)
        verdict == 0 || return verdict > 0
        return GO.relate_predicate(prep, GO.pred_intersects(),
            GI.Polygon([GI.LinearRing(ring)]))
    end
    return GO.relate_predicate(prep, GO.pred_intersects(),
        cell_polygon_unitsphere(system, leaf, id))
end

#=
## The rim sandwich

Measured on the shipped descent, the polygon predicate is where a `:touches`
query spends ~90% of its time, and **most of those calls answer "no"**: on a
HEALPix level-6 globe and a 10°×10° box, 155 cells reach `pred_intersects`
and 131 of them are disjoint — the cap prune can only bound cells against the
*whole geometry's* cap, so every cell in that cap whose center falls outside
the region pays a full polygon–polygon relate (~10 µs cold, against ~0.06 µs
for the prepared center query). On a long thin geometry the ratio is worse
still: 6629 of 6730 calls answer "no".

What separates those two populations is one number: `d`, the spherical
distance from the cell center to the geometry's **boundary** (`_boundary_arcs`
below — every great-circle edge of every ring of every part). Sandwich it
between two radii of the cell:

* `d > r_out`, the cell polygon's circumradius ⟹ **disjoint**. No boundary
  point lies in the cell, so the cell — connected — is wholly inside or
  wholly outside the region; its center is outside (the `pred_contains` above
  just said so), therefore all of it is.
* `d < r_in`, a lower bound on the cell's inradius ⟹ **intersects**. The
  nearest boundary point is then strictly inside the cell, and it belongs to
  the geometry (a closed region contains its boundary), so the two share it.

Only the annulus `r_in ≤ d ≤ r_out` — cells the boundary genuinely grazes —
still needs the exact predicate. Both arms are proofs: results are identical
to running `pred_intersects` on every cell, which
`test/core/test_tree_queries.jl` asserts against its brute-force oracle.

Both radii are measured against the polygon `pred_intersects` is being spared,
i.e. [`cell_polygon_unitsphere`](@ref)'s **great-circle** ring — deliberately
not against [`cell_cap`](@ref), whose radius is a different quantity per
system (HEALPix overrides it with an exact cap over the *curved* chart
boundary; the aperture-7 systems inflate by [`cell_cap_inflation`](@ref) so
that a leaf cap also covers its descendants). Neither is the bound this
argument needs, and both are looser than it.

`r_out` is the ring's max center-to-vertex distance, which bounds the whole
polygon: along a great-circle arc `t ↦ cos(dist(center, p(t)))` solves
`g'' = -Ω²g`, so where it is positive it is concave and dips no lower than at
the arc's ends — the boundary's farthest point from the center is therefore a
vertex — and the cap of that radius, being under a quarter sphere, has a
connected complement of area over `2π` that no cell region can contain, so
the interior is inside the cap too. `r_in` is the distance from the center to
the nearest great circle *carrying* a ring edge, which is at most the distance
to the edge itself and therefore a sound lower bound on the inradius (for a
spherically convex cell the two coincide). Each is nudged by a relative
`SANDWICH_SLACK` in the safe direction so that floating-point noise cannot
turn a grazing cell into a wrong verdict.

The scan is `O(edges of the query geometry)` per cell with no transcendentals
in the loop, which is why it pays: ~0.08 µs for a 5-edge box and ~2.6 µs for a
415-vertex US-state ring, against the ~10 µs it stands in for. Measured
end to end on the `:touches` query it never loses — 1.0–2.2× on the systems
whose native boundary calls dominate the descent anyway (H3, IGEO7, A5) and
1.2–31× on HEALPix, best on the shapes with the most wasted predicates. It
stays a plain scan on purpose: the 415-vertex ring is the worst case measured
and still wins, and indexing the query boundary to shorten the scan would be a
different (and much larger) change.
=#

# Relative margin by which each sandwich radius is nudged the safe way — the
# circumradius out, the inradius in — before a verdict may rest on it. 1 ppm of
# a level-6 cell's radius is ~5 cm on Earth: six orders of magnitude above the
# ~1e-14 rad noise of the dot products below, and small enough that it costs no
# verdicts in practice.
const SANDWICH_SLACK = 1e-6

# One great-circle edge of the query geometry's boundary, with the two derived
# quantities the distance test needs: the (unnormalized) great-circle normal
# `n = a × b` and its squared length. Precomputed once per query, scanned once
# per candidate cell.
struct BoundaryArc
    a::GO.UnitSphericalPoint{Float64}
    b::GO.UnitSphericalPoint{Float64}
    nx::Float64
    ny::Float64
    nz::Float64
    nn::Float64
end

function BoundaryArc(a::GO.UnitSphericalPoint{Float64}, b::GO.UnitSphericalPoint{Float64})
    nx = a[2] * b[3] - a[3] * b[2]
    ny = a[3] * b[1] - a[1] * b[3]
    nz = a[1] * b[2] - a[2] * b[1]
    return BoundaryArc(a, b, nx, ny, nz, nx * nx + ny * ny + nz * nz)
end

#=
Every great-circle edge of the point set the sandwich argument needs, or
`nothing` for a geometry kind this walk cannot enumerate — in which case the
sandwich is simply off and the exact predicate runs as before.

The set must satisfy two properties, and both arms rest on exactly these:
its points all belong to the geometry (so a point of it inside a cell proves
intersection), and a cell meeting none of it cannot straddle inside and
outside (so a cell with an outside center that misses it is disjoint). For a
polygonal geometry that set is the topological boundary — **all** rings of
**all** parts, holes included, which is why the walk goes through
`GI.getring` rather than `GI.getpoint`: consecutive points across two rings
are not an edge, and inventing one could produce a wrong accept. For a
curve or a point set the geometry *is* its own such set.

`LinearRing` and `GeometryCollection` are deliberately unhandled: a bare ring
is ambiguous about closure and a collection can mix dimensions, and neither is
worth a soundness argument when the fallback is just the existing exact test.
=#
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
    point = _query_point(geom)
    push!(arcs, BoundaryArc(point, point))
    return _finish_arcs!(arcs)
end

function _boundary_arcs!(arcs, ::GI.MultiPointTrait, geom)
    for p in GI.getpoint(geom)
        point = _query_point(p)
        push!(arcs, BoundaryArc(point, point))
    end
    return _finish_arcs!(arcs)
end

function _rings_arcs!(arcs, poly)
    for ring in GI.getring(poly)
        _chain_arcs!(arcs, ring, true)
    end
    return arcs
end

# Consecutive-vertex edges of one ring or linestring; `close` adds the edge
# back to the first vertex, which a polygon ring always has whether or not its
# coordinates repeat it (a duplicated closing vertex just makes that edge
# degenerate, and degenerate arcs are handled).
function _chain_arcs!(arcs, chain, close::Bool)
    first_point = nothing
    previous = nothing
    for p in GI.getpoint(chain)
        point = _query_point(p)
        previous === nothing ? (first_point = point) :
            push!(arcs, BoundaryArc(previous, point))
        previous = point
    end
    close && previous !== nothing && push!(arcs, BoundaryArc(previous, first_point))
    return arcs
end

# An empty edge set would let the reject arm fire on every cell; there is no
# geometry to be near, but there is also nothing to gain, so switch off.
_finish_arcs!(arcs) = isempty(arcs) ? nothing : arcs

# `cos` of the distance from `c` to the great-circle arc — the nearer endpoint
# unless the foot of the perpendicular falls on the arc itself, in which case
# the foot is nearer. Everything is dot and cross products on the raw
# coordinates: no normalization, one `sqrt`, no inverse trigonometry.
@inline function _arc_cos_distance(arc::BoundaryArc, c)
    ca = c[1] * arc.a[1] + c[2] * arc.a[2] + c[3] * arc.a[3]
    cb = c[1] * arc.b[1] + c[2] * arc.b[2] + c[3] * arc.b[3]
    best = max(ca, cb)
    arc.nn > 0 || return best        # degenerate arc: its endpoints are all of it
    # The near foot is `f ∝ c - ((c⋅n)/n⋅n) n`, and that second term drops out
    # of both triple products below, so `a → f → b` in `n`'s orientation is
    # decided by `c` directly. If it fails the arc's minimum is at an endpoint.
    ((arc.a[2] * c[3] - arc.a[3] * c[2]) * arc.nx +
     (arc.a[3] * c[1] - arc.a[1] * c[3]) * arc.ny +
     (arc.a[1] * c[2] - arc.a[2] * c[1]) * arc.nz) > 0 || return best
    ((c[2] * arc.b[3] - c[3] * arc.b[2]) * arc.nx +
     (c[3] * arc.b[1] - c[1] * arc.b[3]) * arc.ny +
     (c[1] * arc.b[2] - c[2] * arc.b[1]) * arc.nz) > 0 || return best
    cn = c[1] * arc.nx + c[2] * arc.ny + c[3] * arc.nz
    # sin(distance to the great circle) = |c⋅n| / |n|
    return max(best, sqrt(max(0.0, 1.0 - cn * cn / arc.nn)))
end

# `1` the cell provably meets the geometry, `-1` it provably does not, `0`
# undecided — see the note above for why each arm is a proof.
function _sandwich(arcs::Vector{BoundaryArc}, center, ring)
    # r_out: the ring's max center-to-vertex distance, which bounds the whole
    # great-circle cell polygon (see the note above).
    vertex = 0.0
    for point in ring
        vertex = max(vertex, GO.UnitSpherical.spherical_distance(center, point))
    end
    r_out = vertex * (1 + SANDWICH_SLACK)
    # A quarter sphere is where both the concavity argument for r_out and the
    # sin/cos comparisons below stop holding. No cell of any wired system comes
    # near it, and giving up is free.
    r_out < Float64(pi) / 2 || return 0
    cos_out = cos(r_out)
    # r_in: distance from the center to the nearest great circle carrying a
    # ring edge — at most the distance to the edge, hence a sound lower bound
    # on the inradius. Every edge has to be in the minimum or the bound could
    # come out too large, so the loop wraps rather than trusting `closed=true`
    # to have repeated the first vertex (where it did, that edge is degenerate
    # and drops out below).
    sin2_in = 1.0
    @inbounds for i in eachindex(ring)
        a = ring[i]
        b = ring[i == lastindex(ring) ? firstindex(ring) : i + 1]
        nx = a[2] * b[3] - a[3] * b[2]
        ny = a[3] * b[1] - a[1] * b[3]
        nz = a[1] * b[2] - a[2] * b[1]
        nn = nx * nx + ny * ny + nz * nz
        nn > 0 || continue
        cn = center[1] * nx + center[2] * ny + center[3] * nz
        sin2_in = min(sin2_in, cn * cn / nn)
    end
    r_in = asin(sqrt(min(1.0, sin2_in)))
    # A ring that yields no usable inradius (all-degenerate edges, or one
    # wider than the cap it sits in) keeps only the reject arm: `cos_in > 1`
    # can never be exceeded.
    cos_in = 0 < r_in < r_out ? cos(r_in * (1 - SANDWICH_SLACK)) : 2.0
    best = -1.0
    for arc in arcs
        value = _arc_cos_distance(arc, center)
        value > cos_in && return 1          # d < r_in: boundary point inside the cell
        best = max(best, value)
    end
    return best < cos_out ? -1 : 0          # d > r_out: nothing of the boundary is near
end

#=
Bounding cap of a spherical geometry, for the conservative node prune. Two
facts make a rigorous cap out of nothing but the geometry's vertices and one
prepared point query:

1. A spherical cap of radius ≤ π/2 is geodesically convex, so a cap centered
   on the normalized vertex mean with radius = max distance to any vertex
   contains every great-circle edge between consecutive vertices — i.e. the
   geometry's entire *boundary*.
2. The cap's complement is an open cap — connected — containing no boundary
   point, so it lies entirely inside or entirely outside the geometry, and
   one point decides which: the cap center's antipode, the deepest point of
   the complement. Outside the geometry ⟹ the whole complement is outside ⟹
   the geometry (interior, boundary, holes and all, multi-parts included) is
   inside the cap.

A vertex radius past π/2 breaks the convexity argument, and an antipode
inside the geometry means the region really does extend into the complement;
both give up and return the full sphere — no pruning, never a wrong prune.
The radius is inflated as in `cells_cap`, which keeps the containments strict
against rounding while staying compatible with `intersects_cap`'s closed-cap
arithmetic.
=#
function _geometry_cap(prep, geom)
    points = GO.UnitSphericalPoint{Float64}[]
    for p in GI.getpoint(geom)
        push!(points, _query_point(p))
    end
    isempty(points) && return full_sphere_extent()
    mean = reduce(+, points) / length(points)
    norm = sqrt(sum(abs2, mean))
    norm <= eps(Float64) && return full_sphere_extent()
    center = GO.UnitSphericalPoint(mean[1] / norm, mean[2] / norm, mean[3] / norm)
    radius = 0.0
    for point in points
        radius = max(radius, GO.UnitSpherical.spherical_distance(center, point))
    end
    radius > pi / 2 && return full_sphere_extent()
    antipode = GO.UnitSphericalPoint(-center[1], -center[2], -center[3])
    GO.relate_predicate(prep, GO.pred_intersects(), antipode) && return full_sphere_extent()
    return GO.UnitSpherical.SphericalCap(center,
        nextfloat(min(Float64(pi), radius * 1.0001 + 1e-12)))
end

# Query-geometry vertex on the unit sphere: 2-D coordinates are lon/lat
# degrees, 3-D coordinates are unit-sphere xyz as-is — the same convention
# the spherical RelateNG kernel ingests vertices with.
_query_point(p) = GI.is3d(p) ?
    GO.UnitSphericalPoint(Float64(GI.x(p)), Float64(GI.y(p)), Float64(GI.z(p))) :
    _unit_sphere_point(GI.x(p), GI.y(p))

function _unit_sphere_point(lon::Real, lat::Real)
    # Float64 like the 3-D branch above: a Float32-coordinate geometry (GeoJSON
    # reads Natural Earth that way) would otherwise carry its element type into
    # `UnitSphericalPoint` and miss the concrete `BoundaryArc` fields.
    lon64, lat64 = Float64(lon), Float64(lat)
    coslat = cosd(lat64)
    return GO.UnitSphericalPoint(coslat * cosd(lon64), coslat * sind(lon64), sind(lat64))
end
