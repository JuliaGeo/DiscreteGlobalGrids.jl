# Raster membership is separate from allocation/reduction. Edge-tree frontiers
# are private to one traversal; GeometryOps prepares point location only once.
struct _RasterEdgeTree
    arcs::Vector{Engine.BoundaryArc}
    caps::Vector{Engine.Cap}
    first::Vector{Int}
    last::Vector{Int}
    left::Vector{Int}
    right::Vector{Int}
end

function _RasterEdgeTree(arcs::Vector{Engine.BoundaryArc})
    caps = [Fallbacks.points_cap((a.a, a.b)) for a in arcs]
    order = sortperm([Engine._morton_key(c.point) for c in caps])
    tree = _RasterEdgeTree(arcs[order], Engine.Cap[], Int[], Int[], Int[], Int[])
    _raster_edge_node!(tree, caps[order], 1, length(arcs))
    return tree
end

function _raster_edge_node!(tree, caps, lo, hi)
    n = length(tree.caps) + 1
    push!(tree.caps, Fallbacks.full_sphere_cap())
    push!(tree.first, lo); push!(tree.last, hi)
    push!(tree.left, 0); push!(tree.right, 0)
    if hi - lo < 8
        tree.caps[n] = Engine._merge_range(caps, lo, hi)
    else
        mid = (lo + hi) ÷ 2
        l = _raster_edge_node!(tree, caps, lo, mid)
        r = _raster_edge_node!(tree, caps, mid + 1, hi)
        tree.left[n] = l; tree.right[n] = r
        tree.caps[n] = Extents.union(tree.caps[l], tree.caps[r])
    end
    return n
end

struct _PreparedRasterGeometry{T,E}
    target::T
    edges::E
    polygon::Bool
end

function _raster_boundary(boundary)
    boundary === :touches && return :intersects
    boundary in (:center, :inside, :intersects) || throw(ArgumentError(
        "boundary must be :center, :inside, :intersects, or :touches"))
    return boundary
end

function _raster_shape(shape)
    shape === nothing && return nothing
    shape in (:point, :points) && return :point
    shape in (:line, :lines, :linestring, :linestrings) && return :line
    shape in (:polygon, :polygons) && return :polygon
    throw(ArgumentError("shape must be :point, :line, :polygon, or nothing"))
end

function _raster_ring_line(ring)
    points = collect(GI.getpoint(ring))
    isempty(points) || first(points) == last(points) || push!(points, first(points))
    return GI.LineString(points)
end

function _raster_shaped_geometry(geom, shape)
    shape = _raster_shape(shape)
    shape === nothing && return geom
    trait = GI.trait(geom)
    shape === :point && return trait isa GI.PointTrait ? geom : GI.MultiPoint(collect(GI.getpoint(geom)))
    if shape === :line
        trait isa GI.PolygonTrait && return GI.MultiLineString([
            _raster_ring_line(r) for r in GI.getring(geom)])
        trait isa GI.MultiPolygonTrait && return GI.MultiLineString([
            _raster_ring_line(r) for p in GI.getpolygon(geom)
            for r in GI.getring(p)])
        trait isa GI.MultiLineStringTrait && return geom
        return GI.LineString(collect(GI.getpoint(geom)))
    end
    trait isa Union{GI.PolygonTrait,GI.MultiPolygonTrait} && return geom
    return GI.Polygon([GI.LinearRing(collect(GI.getpoint(geom)))])
end

function _prepare_raster_geometry(geom; shape=nothing)
    geom = _raster_shaped_geometry(geom, shape)
    trait = GI.trait(geom)
    polygon = trait isa Union{GI.PolygonTrait,GI.MultiPolygonTrait}
    target = Engine._query_target(geom)
    singleton = trait isa GI.MultiLineStringTrait && any(l -> GI.npoint(l) == 1, GI.getgeom(geom))
    edges = target.arcs === nothing || singleton ? nothing : _RasterEdgeTree(target.arcs)
    return _PreparedRasterGeometry(target, edges, polygon)
end

# Retain whole edge nodes when splitting the grid is cheaper; expand only
# boundary nodes larger than the grid extent. Terminal edge buckets use exact
# arc distance, never a full feature-boundary scan.
function _raster_edge_frontier!(out, tree, n, cap)
    Extents.intersects(tree.caps[n], cap) || return out
    if tree.left[n] == 0
        threshold = cos(min(Float64(pi), cap.radius * (1 + Engine.SANDWICH_SLACK) + 1e-12))
        for k in tree.first[n]:tree.last[n]
            if Engine._arc_cos_distance(tree.arcs[k], cap.point) >= threshold
                push!(out, n)
                break
            end
        end
    elseif tree.caps[n].radius > cap.radius
        _raster_edge_frontier!(out, tree, tree.left[n], cap)
        _raster_edge_frontier!(out, tree, tree.right[n], cap)
    else
        push!(out, n)
    end
    return out
end

_raster_node_indices(node::Engine.HierarchicalGridCursor) = Engine.node_indices(node)
_raster_node_indices(node::Engine.IndexTreeNode) =
    @view node.tree.order[node.tree.node_first[node.index]:node.tree.node_last[node.index]]
_raster_leaf_indices(node::Union{Engine.HierarchicalGridCursor,Engine.IndexTreeNode}) =
    _raster_node_indices(node)
_raster_leaf_indices(node) = (i for (i, _) in Engine.STI.child_indices_extents(node))

function _raster_witness_index(node)
    if node isa Union{Engine.HierarchicalGridCursor,Engine.IndexTreeNode}
        return first(_raster_node_indices(node))
    end
    Engine.STI.isleaf(node) && return first(_raster_leaf_indices(node))
    return _raster_witness_index(first(Engine.STI.getchild(node)))
end

function _raster_emit!(out, node)
    if node isa Union{Engine.HierarchicalGridCursor,Engine.IndexTreeNode}
        append!(out, _raster_node_indices(node))
    elseif Engine.STI.isleaf(node)
        append!(out, _raster_leaf_indices(node))
    else
        for child in Engine.STI.getchild(node)
            _raster_emit!(out, child)
        end
    end
    return out
end

function _raster_leaf_matches(prepared, grid, i, boundary)
    c = cellindex(grid, i)
    if prepared.polygon
        covered = _raster_relate(prepared.target, GO.pred_intersects,
            cell_centroid(grid, c))
        boundary === :center && return covered
        boundary === :intersects && covered && return true
        boundary === :inside && !covered && return false
    end
    predicate = prepared.polygon && boundary === :inside ?
        GO.pred_contains : GO.pred_intersects
    return _raster_relate(prepared.target, predicate, cell_polygon(grid, c))
end

# The intersection sandwich scans only edge buckets that survived both trees.
# Excluded ancestral arcs cannot cross this cell, by the extent covering law.
function _raster_candidate_matches(prepared, grid, i, boundary, frontier)
    (boundary === :intersects || !prepared.polygon) ||
        return _raster_leaf_matches(prepared, grid, i, boundary)
    prepared.edges === nothing && return _raster_leaf_matches(prepared, grid, i, boundary)
    c = cellindex(grid, i)
    centroid = cell_centroid(grid, c)
    prepared.polygon && _raster_relate(prepared.target,
        GO.pred_intersects, centroid) && return true
    ring = Fallbacks.closed_ring(cell_boundary(grid, c))
    cap = Fallbacks.points_cap(ring)
    localfrontier = Int[]
    for n in frontier
        _raster_edge_frontier!(localfrontier, prepared.edges, n, cap)
    end
    arcs = Engine.BoundaryArc[]
    for n in localfrontier, k in prepared.edges.first[n]:prepared.edges.last[n]
        push!(arcs, prepared.edges.arcs[k])
    end
    verdict = Engine._sandwich(arcs, centroid, ring)
    verdict == 0 || return verdict > 0
    return _raster_relate(prepared.target, GO.pred_intersects,
        GI.Polygon([GI.LinearRing(ring)]))
end

function _raster_descend!(out, node, grid, prepared, frontier, boundary)
    cap = Engine.STI.node_extent(node)
    Extents.intersects(prepared.target.cap, cap) || return out
    nextfrontier = frontier
    # Only convex caps contribute edge exclusions. Their intersection with all
    # earlier convex caps is connected, even when child caps protrude. Every
    # descendant geometry lies in that intersection by node_extent's contract.
    # Consequently a real descendant centroid is a certified interior witness;
    # the current cap centre is NOT, and would cause false bulk acceptance.
    if prepared.edges !== nothing && cap.radius < Float64(pi) / 2
        nextfrontier = Int[]
        for n in frontier
            _raster_edge_frontier!(nextfrontier, prepared.edges, n, cap)
        end
        if isempty(nextfrontier)
            prepared.polygon || return out
            i = _raster_witness_index(node)
            covered = _raster_relate(prepared.target, GO.pred_intersects,
                cell_centroid(grid, cellindex(grid, i)))
            covered && _raster_emit!(out, node)
            return out
        end
    end
    if Engine.STI.isleaf(node)
        for i in _raster_leaf_indices(node)
            _raster_candidate_matches(prepared, grid, i, boundary, nextfrontier) && push!(out, i)
        end
    else
        for child in Engine.STI.getchild(node)
            _raster_descend!(out, child, grid, prepared, nextfrontier, boundary)
        end
    end
    return out
end

function _raster_indices(grid::AbstractGrid, prepared::_PreparedRasterGeometry;
        boundary=:center, shape=nothing)
    boundary = _raster_boundary(boundary)
    shape === nothing || throw(ArgumentError("apply shape when preparing geometry"))
    ncells(grid) == 0 && return Int[]
    # A level-less adapter names stored multi-order cells, not a complete level.
    node = level(grid) === nothing ? Engine.IndexTreeNode(Engine.IndexTree(grid), 1) : treeify(grid)
    return sort!(unique!(_raster_descend!(Int[], node, grid, prepared, [1], boundary)))
end

function _raster_empty_geometry(geom)
    (ismissing(geom) || geom === nothing) && return false
    trait = GI.trait(geom)
    trait === nothing && return false
    trait isa GI.FeatureTrait && return _raster_empty_geometry(GI.geometry(geom))
    if trait isa Union{GI.GeometryCollectionTrait,GI.FeatureCollectionTrait}
        parts = trait isa GI.FeatureCollectionTrait ? GI.getfeature(geom) : GI.getgeom(geom)
        return all(_raster_empty_geometry, parts)
    end
    GI.isempty(geom) && return true
    trait isa GI.PointTrait && return false
    # Flattening getpoint can fail for a polygon with no rings. Count its
    # points through the ring hierarchy before asking for coordinates.
    trait isa Union{GI.PolygonTrait,GI.MultiPolygonTrait} && return GI.npoint(geom) == 0
    return isempty(GI.getpoint(geom))
end
_raster_empty_geometry(::Union{Extents.Extent,GO.UnitSpherical.SphericalCap}) = false

function _raster_indices(grid::AbstractGrid, geom; boundary=:center, shape=nothing)
    boundary = _raster_boundary(boundary)
    shape = _raster_shape(shape)
    (ismissing(geom) || geom === nothing) && return Int[]
    trait = GI.trait(geom)
    trait isa GI.FeatureTrait && return _raster_indices(grid, GI.geometry(geom); boundary, shape)
    if trait isa Union{GI.GeometryCollectionTrait,GI.FeatureCollectionTrait}
        out = Int[]
        parts = trait isa GI.FeatureCollectionTrait ? GI.getfeature(geom) : GI.getgeom(geom)
        for part in parts
            append!(out, _raster_indices(grid, part; boundary, shape))
        end
        return sort!(unique!(out))
    end
    _raster_empty_geometry(geom) && return Int[]
    geom = _raster_shaped_geometry(geom, shape)
    trait = GI.trait(geom)
    if trait isa Union{GI.PointTrait,GI.MultiPointTrait}
        out = Int[]
        points = trait isa GI.PointTrait ? (geom,) : GI.getpoint(geom)
        for point in points
            c = cellat(grid, Fallbacks.query_point(point))
            c === nothing && continue
            i = localindex(grid, c)
            i === nothing || push!(out, i)
        end
        return sort!(unique!(out))
    end
    _raster_empty_geometry(geom) && return Int[]
    trait isa GI.LinearRingTrait && (geom = _raster_ring_line(geom))
    return _raster_indices(grid, _prepare_raster_geometry(geom); boundary)
end

# Extents use the query API's densified longitude/latitude outline. Polar
# full-longitude extents have a native spherical-cap representation instead.
function _raster_indices(grid::AbstractGrid, extent::Extents.Extent;
        boundary=:center,shape=nothing)
    shape === nothing || throw(ArgumentError("shape overrides require a GeoInterface geometry"))
    return _raster_indices(grid,Engine._extent_target(extent);boundary)
end
function _raster_indices(grid::AbstractGrid, cap::GO.UnitSpherical.SphericalCap;
        boundary=:center,shape=nothing)
    boundary=_raster_boundary(boundary)
    shape === nothing || throw(ArgumentError("shape overrides require a GeoInterface geometry"))
    ncells(grid)==0 && return Int[]
    target=Engine._query_target(cap)
    predicate=boundary===:inside ? DE9IM.Within(nothing) : DE9IM.Intersects(nothing)
    indices=Engine._query_indices(grid,predicate,target)
    if boundary===:center
        filter!(indices) do i
            GO.UnitSpherical.spherical_distance(cell_centroid(grid,cellindex(grid,i)),target.cap.point) <= target.cap.radius
        end
    end
    return indices
end
