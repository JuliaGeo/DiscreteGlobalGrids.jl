# test/core/test_node_extent_containment.jl — the one contract a node extent
# has, checked directly on the wired systems rather than inferred from the rule
# `STI.node_extent` implements (`src/core/generic_cursor.jl`):
#
#   the cap a node reports must contain every point of every
#   `cell_polygon_unitsphere` of every leaf `node_indices(node)` names.
#
# A superset is fine and expected — a node extent is a prune filter, and the
# exact geometry is tested downstream. A cap that misses a stored leaf is not:
# the dual depth-first search behind a `Regridder` prunes that pair and the
# intersection quietly disappears from the matrix, which is a wrong answer, not
# an error. Nothing else in the suite would catch it either, so this file walks
# whole trees and measures the worst overshoot in radians.
#
# Since the internal-node extent is the node's *own* cell cap (O(1) from the
# hierarchy — `cell_cap_inflation` is what makes that sound where children
# overhang their parent), the containment being asserted here is exactly the
# per-system CAP-VALIDATION measurement, read through the trees that depend on
# it and including the paths those sweeps do not reach: the sparse partial
# nodes that keep a union cap, bucketed leaves, the pentagon subtrees, and the
# `cell_parent` descent mode (A5) whose nodes carry materialized selections.
#
# Sampled points per cell: the ring vertices, plus three slerp points inside
# each edge. The polygon a traversal consumes is its vertices joined by
# geodesics, so those samples are on the geometry itself, not near it; the caps
# in play are far under a quarter turn and therefore convex, so vertices alone
# would already imply the edges, and the interior samples are insurance
# against that argument rather than the argument itself.
module NodeExtentContainmentTests

using Test
using Printf

using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
import GeometryOps as GO
import GeometryOps: SpatialTreeInterface as STI
import ConservativeRegridding: Trees

const SD = GO.UnitSpherical.spherical_distance

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

"Depth-first walk over every node of `tree`, root included."
function walk_nodes(f::Function, node)
    f(node)
    STI.isleaf(node) && return nothing
    for child in STI.getchild(node)
        walk_nodes(f, child)
    end
    return nothing
end

"Canonical id behind tree-level leaf index `index`."
leaf_id_at(tree::DGG.DGGSCursor{<:DGGSPartialGrid}, index::Int) = tree.grid.ids[index]
leaf_id_at(tree::DGG.DGGSCursor{<:DGGSGrid}, index::Int) =
    ordinal_to_cell(tree.grid.system, tree.grid.level, index)

"Points of the cell's polygon: every vertex, plus three inside every edge."
function polygon_samples(system, level, id)
    ring = collect(GO.UnitSpherical.UnitSphericalPoint{Float64},
        DGG.cell_boundary(system, level, id; closed=true))
    points = similar(ring, 0)
    sizehint!(points, 4 * length(ring))
    for i in 1:(length(ring) - 1)
        push!(points, ring[i])
        for t in (0.25, 0.5, 0.75)
            push!(points, GO.UnitSpherical.slerp(ring[i], ring[i + 1], t))
        end
    end
    return points
end

"""
    containment(tree) -> (worst, cell_worst, nodes, checks)

Worst `distance(extent.point, p) - extent.radius` over every node of `tree` and
every sampled polygon point `p` of every leaf that node owns — negative means
contained, with the magnitude as the margin. `cell_worst` is the same over the
nodes that stand for a cell: the synthetic whole-sphere root reports a cap of
radius `nextfloat(pi)`, which any antipodal sample sits one ulp inside, so it
would otherwise be the reported margin of every whole-globe tree. Leaf geometry
is sampled once per leaf index and reused, so the walk costs one
`cell_boundary` per leaf.
"""
function containment(tree)
    system = tree.grid.system
    leaf_level = tree.grid.level
    cache = Dict{Int,Vector{GO.UnitSpherical.UnitSphericalPoint{Float64}}}()
    samples(index) = get!(cache, index) do
        polygon_samples(system, leaf_level, leaf_id_at(tree, index))
    end
    worst = -Inf
    cell_worst = -Inf
    nodes = 0
    checks = 0
    walk_nodes(tree) do node
        extent = STI.node_extent(node)
        nodes += 1
        for index in node_indices(node), point in samples(index)
            overshoot = SD(extent.point, point) - extent.radius
            worst = max(worst, overshoot)
            node_level(node) >= 0 && (cell_worst = max(cell_worst, overshoot))
            checks += 1
        end
        return nothing
    end
    return (worst, cell_worst, nodes, checks)
end

"Every node's extent contains every polygon point of every leaf it owns."
function test_contains(label, tree; min_nodes=2)
    worst, cell_worst, nodes, checks = containment(tree)
    @test nodes >= min_nodes
    @test checks > 0
    @test worst <= 0
    @test cell_worst <= 0
    @printf("  %-42s %6d nodes %9d point checks, margin %.3e rad\n",
        label, nodes, checks, -cell_worst)
    return nothing
end

"The `level` cell reached by taking the first (center) child of `id` all the way down."
function first_descendant(system, id, level)
    for l in 0:(Int(level) - 1)
        id = first(cell_children(system, l, id))
    end
    return id
end

"""
    root_pair(system, level) -> (id, id)

Two cells at `level`, one under a pentagon root and one under a hexagon root
where the system has both (H3, A5); two different roots where every root is the
same shape (IGEO7's twelve pentagons, HEALPix's twelve pixels). Descending by
the center child keeps a pentagon a pentagon, which is the case whose child
count, subtree size and geometry all differ from the rest of the grid.
"""
function root_pair(system, level)
    roots = collect(DGG.root_ids(system))
    counts = [length(cell_children(system, 0, id)) for id in roots]
    small, large = roots[argmin(counts)], roots[argmax(counts)]
    small == large && (large = last(roots))
    return (first_descendant(system, small, level), first_descendant(system, large, level))
end

"Sparse partial grids over the ids of a full chunk: a stride and a scatter."
function sparse_grids(system, chunk)
    ids = chunk.ids
    stride = max(1, length(ids) ÷ 40)
    scattered = ids[1:stride:end]                       # a few per subtree
    clustered = ids[1:min(length(ids), 25)]             # one corner of it
    return (("stride", DGGSPartialGrid(system, chunk.level, scattered)),
        ("cluster", DGGSPartialGrid(system, chunk.level, clustered)),
        ("stride/bucket", DGGSPartialGrid(system, chunk.level, scattered; bucket_size=8)))
end

# --------------------------------------------------------------------------
# The systems, and the trees each of them is checked on.
#
# `root_level` / `leaf_level` are chosen so the deepest chunk is a few thousand
# leaves: enough depth that a node's extent is asked to bound several levels of
# descendants, small enough that the walk stays a second or two. Each system
# contributes a dense whole-globe grid (the `subtree_cap` path), a full subtree
# chunk (every internal node stores its whole subtree — the O(1) path), and the
# sparse grids above (the union-cap path and the bucketed leaves).
# --------------------------------------------------------------------------

const CASES = (
    # (system, dense level, chunk root_level, chunk leaf_level)
    (IGEO7DGGS(), 2, 0, 4),
    (IGEO7DGGS(), 3, 2, 6),
    (H3DGGS(), 1, 0, 3),
    (H3DGGS(), 2, 1, 5),
    (HEALPixDGGS(), 3, 0, 5),
    (A5DGGS(), 2, 1, 5),
)

@testset "node extents contain the leaves they own" begin
    println("\n  node-extent containment (margin = how far inside the cap the ",
        "worst leaf point sits)")
    for (system, dense_level, root_level, leaf_level) in CASES
        name = system_name(system)
        @testset "$name (dense $dense_level, chunk $root_level->$leaf_level)" begin
            test_contains("$name dense grid, level $dense_level",
                treeify(DGGSGrid(system, dense_level)))

            for root in root_pair(system, root_level)
                chunk = subtree_grid(system, root; root_level=root_level, leaf_level=leaf_level)
                test_contains("$name chunk $root_level->$leaf_level, root $root",
                    treeify(chunk))
                for (kind, grid) in sparse_grids(system, chunk)
                    test_contains("$name chunk $kind, root $root", treeify(grid))
                end
            end
        end
    end
    println()
end

# --------------------------------------------------------------------------
# ...and the rule that containment is bought with: which cap a node reports.
# The containment above is what must never break; this is what makes it O(1).
# --------------------------------------------------------------------------

@testset "internal nodes report their own cell's cap" begin
    for (system, _, root_level, leaf_level) in CASES
        chunk = subtree_grid(system, first(root_pair(system, root_level));
            root_level=root_level, leaf_level=leaf_level)
        tree = treeify(chunk)
        internal = 0
        walk_nodes(tree) do node
            node_level(node) < leaf_level || return nothing
            internal += 1
            # A node that stores its whole subtree — which is every internal
            # node of a chunk — takes the hierarchy's O(1) answer. `cells_cap`
            # over the same leaves would differ (it is a different center and
            # a 1.0001 rather than a `cell_cap_inflation` radius), so this
            # equality is what pins the enumeration out of the traversal.
            expected = has_exact_subtree_cap(system) ?
                       DGG.subtree_cap(system, node_level(node), node_id(node), leaf_level) :
                       cell_cap(system, node_level(node), node_id(node))
            @test STI.node_extent(node) == expected
            return nothing
        end
        @test internal > 1
    end
end

@testset "sparse nodes still get the tighter union cap" begin
    # The other half of the rule: a node storing a handful of a big cell's
    # leaves bounds those leaves, not its cell — that tightness is the pruning
    # power a sparse partial grid is built for, and it is bounded work per node
    # (`STORED_UNION_CAP_LIMIT`) rather than work proportional to the subtree.
    @test DGG.STORED_UNION_CAP_LIMIT >= 1
    for (system, _, root_level, leaf_level) in CASES
        has_exact_subtree_cap(system) && continue   # HEALPix: the cell cap is exact
        chunk = subtree_grid(system, first(root_pair(system, root_level));
            root_level=root_level, leaf_level=leaf_level)
        grid = DGGSPartialGrid(system, chunk.level, chunk.ids[1:3])
        tree = treeify(grid)
        tighter = 0
        walk_nodes(tree) do node
            0 <= node_level(node) < leaf_level || return nothing   # not the synthetic root
            stored = [grid.ids[i] for i in node_indices(node)]
            length(stored) <= DGG.STORED_UNION_CAP_LIMIT || return nothing
            length(stored) < subtree_leaf_count(system, node_level(node),
                node_id(node), leaf_level) || return nothing
            extent = STI.node_extent(node)
            @test extent == DGG.cells_cap(system, leaf_level, stored)
            extent.radius < cell_cap(system, node_level(node), node_id(node)).radius &&
                (tighter += 1)
            return nothing
        end
        # ...and it really is tighter, on at least the nodes near the root of a
        # three-cell grid, where the cell they sit in is enormous beside them.
        @test tighter > 0
    end
end

end # module NodeExtentContainmentTests
