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
    out = Int[]
    _tree_query!(out, tree, ids, system, leaf, prep, cap, mode)
    # The range-backed cursors already yield ascending positions; only a
    # selection-cursor system without descendant ranges (A5) can interleave
    # across children, and a near-sorted sort! is cheap.
    return sort!(out)
end

function _tree_query!(out, node, ids, system, leaf, prep, cap, mode)
    extent = STI.node_extent(node)
    intersects_cap(cap, extent) || return nothing
    if STI.isleaf(node)
        indices = node_indices(node)
        if length(indices) == 1
            # A single-cell leaf's node extent IS its cell cap — already
            # tested on entry, so go straight to the exact test.
            index = first(indices)
            _leaf_hit(prep, system, leaf, ids[index], mode) && push!(out, index)
        else
            # A bucket whose cap sits entirely inside the geometry's cap has
            # no per-cell cap prune left to win — every cell cap would pass —
            # so skip straight to the exact tests.
            interior = GO.UnitSpherical.spherical_distance(cap.point, extent.point) +
                       extent.radius <= cap.radius
            for index in indices
                interior || intersects_cap(cap, cell_cap(system, leaf, ids[index])) ||
                    continue
                _leaf_hit(prep, system, leaf, ids[index], mode) && push!(out, index)
            end
        end
        return nothing
    end
    for child in STI.getchild(node)
        _tree_query!(out, child, ids, system, leaf, prep, cap, mode)
    end
    return nothing
end

# The exact per-cell test. `:center` is "cell center in the geometry's
# interior" (`pred_contains`, the DE-9IM sense the old planar
# `GO.contains(geom, center)` had). For `:touches` the same cached point
# query is a fast accept — the center is interior to its cell in every wired
# system, so center ∈ geom ⟹ cell ∩ geom ≠ ∅ — and only cells whose center
# lands outside pay the polygon–polygon `pred_intersects`.
function _leaf_hit(prep, system, leaf, id, mode::Symbol)
    GO.relate_predicate(prep, GO.pred_contains(), cell_center(system, leaf, id)) &&
        return true
    mode === :center && return false
    return GO.relate_predicate(prep, GO.pred_intersects(),
        cell_polygon_unitsphere(system, leaf, id))
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
    coslat = cosd(lat)
    return GO.UnitSphericalPoint(coslat * cosd(lon), coslat * sind(lon), sind(lat))
end
