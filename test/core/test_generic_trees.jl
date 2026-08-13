# test/core/test_generic_trees.jl — the generic DGGS cursor
# (src/core/generic_cursor.jl) against a self-contained mock system.
#
# The mock is an equirectangular quadtree: level 0 is 4 longitude sectors x 2
# latitude bands, every cell splits into 2x2, and `cell_boundary` returns the
# real unit-sphere corners of the lon/lat box. It comes in three trait
# configurations so every cursor path is exercised on identical geometry:
#
#   * `QuadMock{true,false}` — `has_ordinal_ids`, hence
#     `has_descendant_ranges`: the `descendant_range` + `searchsorted`
#     index-bound path (HEALPix, and H3 / IGeo7 once their kernels are wired).
#   * `QuadMock{false,false}` — hierarchy wired natively,
#     `has_descendant_ranges` false: the `cell_parent` membership-filter path
#     with materialized index vectors (A5-style).
#   * `QuadMock{true,true}` — adds `has_exact_subtree_cap`: partial internal
#     nodes take the O(1) `subtree_cap` instead of the O(stored) union cap
#     (HEALPix, where a parent pixel contains its children).
#
# Both must produce the *same* leaf indices for the same query, which is the
# strongest available cross-check of the two implementations.
#
# The mock overrides `cells_cap`/`subtree_cap` with a cap-of-caps fold plus a
# hair of slack, so a node's extent provably contains its descendants' leaf
# caps. That is what makes the brute-force query comparisons *equalities*: with
# containing extents the traversal cannot over-prune, so the returned index set
# is exactly the set of leaves whose `cell_cap` meets the query. (The kernel's
# own vertex-union caps contain descendant *cells*, not descendant *caps*, so
# against them only a subset relation would be assertable.) The nodes that
# report their own `cell_cap` — a node storing its whole subtree takes the
# hierarchy's O(1) answer — carry the same property here for a reason of the
# mock's geometry rather than of the fold: a child box halves its parent's
# radius, so a depth-`d` leaf cap reaches `1 + 0.2 * 2^-d` of the parent's own
# radius against the `1.2` the cap is inflated by. The "mock system
# self-consistency" testset measures exactly that, so the premise is checked
# rather than argued.
module GenericTreeTests

using Test
using Random

using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
import GeometryOps: SpatialTreeInterface as STI
import ConservativeRegridding: Trees

# --------------------------------------------------------------------------
# Mock system
# --------------------------------------------------------------------------

struct QuadMock{ORDINAL,EXACT_CAP} <: AbstractDGGS end

const ROOT_COUNT = 8
const RADIX = 4
const ORDINAL_MOCK = QuadMock{true,false}()
const STRUCTURAL_MOCK = QuadMock{false,false}()
# Third configuration: same grid, but `has_exact_subtree_cap` — parents contain
# their descendants, so partial internal nodes take the O(1) `subtree_cap`
# instead of the O(stored) union cap. HEALPix is the real instance.
const EXACT_CAP_MOCK = QuadMock{true,true}()
const MOCKS = (ORDINAL_MOCK, STRUCTURAL_MOCK)

DGG.system_name(::QuadMock) = :QuadMock
DGG.root_count(::QuadMock) = ROOT_COUNT
DGG.radix(::QuadMock) = RADIX
DGG.max_level(::QuadMock) = 12
DGG.has_ordinal_ids(::QuadMock{ORDINAL}) where {ORDINAL} = ORDINAL
DGG.has_exact_subtree_cap(::QuadMock{ORDINAL,EXACT_CAP}) where {ORDINAL,EXACT_CAP} = EXACT_CAP

# The non-ordinal variant wires the same arithmetic natively, so the two
# variants describe one grid and differ only in which kernel path the cursor
# takes. `has_descendant_ranges` stays false for it (the trait's default).
const Structural = QuadMock{false}

DGG.root_ids(::Structural) = collect(Int64, 0:(ROOT_COUNT - 1))
DGG.cell_children(::Structural, level::Integer, id) =
    collect(Int64, (Int64(id) * RADIX):(Int64(id) * RADIX + RADIX - 1))
DGG.cell_parent(::Structural, level::Integer, id, parent_level::Integer) =
    Int64(id) ÷ Int64(RADIX)^(Int(level) - Int(parent_level))
DGG.num_cells(::Structural, level::Integer) = Int64(ROOT_COUNT) * Int64(RADIX)^Int(level)
DGG.subtree_leaf_count(::Structural, level::Integer, id, leaf_level::Integer) =
    Int64(RADIX)^(Int(leaf_level) - Int(level))
DGG.cell_to_ordinal(::Structural, level::Integer, id) = Int(id) + 1
DGG.ordinal_to_cell(::Structural, level::Integer, ordinal::Integer) = Int64(ordinal - 1)

"Longitude/latitude box (degrees) of cell `(level, id)`."
function cell_box(level::Integer, id::Integer)
    span = RADIX^Int(level)
    root, within = divrem(Int(id), span)
    x = 0
    y = 0
    for bit in 0:(Int(level) - 1)
        x |= ((within >> (2bit)) & 1) << bit
        y |= ((within >> (2bit + 1)) & 1) << bit
    end
    width = 90.0 / (1 << Int(level))
    lon = -180.0 + 90.0 * (root % 4) + x * width
    lat = -90.0 + 90.0 * (root ÷ 4) + y * width
    return (lon, lat, lon + width, lat + width)
end

function unit_point(lon::Real, lat::Real)
    lambda = deg2rad(lon)
    phi = deg2rad(lat)
    cosphi = cos(phi)
    return GO.UnitSphericalPoint(cosphi * cos(lambda), cosphi * sin(lambda), sin(phi))
end

function DGG.cell_boundary(::QuadMock, level::Integer, id; closed::Bool=false)
    lon0, lat0, lon1, lat1 = cell_box(level, id)
    points = [unit_point(lon0, lat0), unit_point(lon1, lat0),
        unit_point(lon1, lat1), unit_point(lon0, lat1)]
    closed && push!(points, points[1])
    return points
end

# Native center (box midpoint): real systems wire one, and it keeps `cell_cap`
# down to a single boundary evaluation.
function DGG.cell_center(::QuadMock, level::Integer, id)
    lon0, lat0, lon1, lat1 = cell_box(level, id)
    return unit_point((lon0 + lon1) / 2, (lat0 + lat1) / 2)
end

"Fold `cell_cap`s into one cap that provably contains all of them."
function DGG.cells_cap(system::QuadMock, level::Integer, ids)
    isempty(ids) && return DGG.full_sphere_extent()
    length(ids) > DGG.SUBTREE_CAP_EXACT_LIMIT && return DGG.full_sphere_extent()
    cap = cell_cap(system, level, first(ids))
    for id in Iterators.drop(ids, 1)
        cap = GO.UnitSpherical._merge(cap, cell_cap(system, level, id))
    end
    # A hair of slack so the containment is robust to the merge's rounding.
    return GO.UnitSpherical.SphericalCap(cap.point,
        nextfloat(min(Float64(pi), cap.radius * (1 + 1e-12) + 1e-12)))
end

function DGG.subtree_cap(system::QuadMock, level::Integer, id, leaf_level::Integer)
    level == leaf_level && return cell_cap(system, level, id)
    subtree_leaf_count(system, level, id, leaf_level) > DGG.SUBTREE_CAP_EXACT_LIMIT &&
        return DGG.full_sphere_extent()
    return DGG.cells_cap(system, leaf_level, cell_descendants(system, level, id, leaf_level))
end

# --------------------------------------------------------------------------
# Traversal helpers
#
# `treeify` below is the package's re-export of `Trees.treeify`, used in its
# one-argument form throughout: the manifold is `best_manifold(grid)`, which is
# `Spherical()` for every DGGS grid. The manifold-explicit form is exercised in
# the "treeify and manifold" testset.
# --------------------------------------------------------------------------

"Every leaf index the tree yields, in traversal order."
all_leaf_indices(tree) = STI.depth_first_search(Returns(true), tree)

"Leaf indices whose extent meets `cap`, via the public depth-first search."
query(tree, cap) =
    sort!(STI.depth_first_search(extent -> GO.UnitSpherical._intersects(cap, extent), tree))

"Ground truth: the leaves whose own `cell_cap` meets `cap`."
brute_force(system, level, ids, cap) =
    [i for i in eachindex(ids)
     if GO.UnitSpherical._intersects(cap, cell_cap(system, level, ids[i]))]

"Depth-first walk over every node (including the root)."
function walk_nodes(f::Function, node)
    f(node)
    STI.isleaf(node) && return nothing
    for child in STI.getchild(node)
        walk_nodes(f, child)
    end
    return nothing
end

descend_to_leaf(node) = STI.isleaf(node) ? node : descend_to_leaf(first(STI.getchild(node)))

ring_points(polygon) = collect(GI.getpoint(GI.getexterior(polygon)))

"Cursors are compared by what they mean, not by identity of their index vectors."
node_key(node) = (node.level, node.id, Trees.ncells(node))

"Canonical id behind tree-level leaf index `index`."
leaf_id_at(tree::DGG.DGGSCursor{<:DGGSPartialGrid}, index::Int) = tree.grid.ids[index]
leaf_id_at(tree::DGG.DGGSCursor{<:DGGSGrid}, index::Int) =
    ordinal_to_cell(tree.grid.system, tree.grid.level, index)

"Sorted, strictly ascending sample of `count` cell ids at `level`."
function sample_ids(rng, system, level, count)
    n = Int(DGG.num_cells(system, level))
    return sort!(Int64.(randperm(rng, n)[1:count]) .- 1)
end

const QUERY_CAPS = (
    GO.UnitSpherical.SphericalCap(unit_point(20.0, 35.0), 0.25),
    GO.UnitSpherical.SphericalCap(unit_point(-140.0, -70.0), 0.4),
    GO.UnitSpherical.SphericalCap(unit_point(0.0, 0.0), 0.05),
    GO.UnitSpherical.SphericalCap(unit_point(179.0, 12.5), 0.6),
    GO.UnitSpherical.SphericalCap(unit_point(-90.0, 89.0), 1.1),
)

const RNG = MersenneTwister(20260805)
const LEVEL = 4
const IDS = sample_ids(RNG, ORDINAL_MOCK, LEVEL, 300)
const SPARSE_IDS = sample_ids(RNG, ORDINAL_MOCK, LEVEL, 17)

# --------------------------------------------------------------------------

@testset "mock system self-consistency" begin
    for system in MOCKS
        @test DGG.cell_id_type(system) === Int64
        @test DGG.num_cells(system, 3) == ROOT_COUNT * RADIX^3
        @test root_ids(system) == collect(Int64, 0:(ROOT_COUNT - 1))
        for level in 0:2, id in 0:(DGG.num_cells(system, level) - 1)
            children = cell_children(system, level, id)
            @test issorted(children)
            @test length(children) == RADIX
            @test all(child -> cell_parent(system, level + 1, child, level) == id, children)
        end
        @test subtree_leaf_count(system, 1, 5, 4) == RADIX^3
        @test length(cell_descendants(system, 1, 5, 4)) == RADIX^3
    end
    @test has_descendant_ranges(ORDINAL_MOCK)
    @test !has_descendant_ranges(STRUCTURAL_MOCK)

    # Child boxes tile the parent box exactly.
    lon0, lat0, lon1, lat1 = cell_box(1, 5)
    child_boxes = [cell_box(2, child) for child in cell_children(ORDINAL_MOCK, 1, 5)]
    @test minimum(b -> b[1], child_boxes) == lon0
    @test minimum(b -> b[2], child_boxes) == lat0
    @test maximum(b -> b[3], child_boxes) == lon1
    @test maximum(b -> b[4], child_boxes) == lat1

    # The overridden caps really do contain their descendants' leaf caps —
    # this is the premise of every brute-force equality below.
    batch = collect(Int64, 40:75)
    batch_cap = DGG.cells_cap(ORDINAL_MOCK, LEVEL, batch)
    for id in batch
        @test GO.UnitSpherical._contains(batch_cap, cell_cap(ORDINAL_MOCK, LEVEL, id))
    end
    node_cap = DGG.subtree_cap(ORDINAL_MOCK, 1, 5, 4)
    for id in cell_descendants(ORDINAL_MOCK, 1, 5, 4)
        @test GO.UnitSpherical._contains(node_cap, cell_cap(ORDINAL_MOCK, LEVEL, id))
    end

    # ...and so does the plain `cell_cap`, which is what a node storing its
    # whole subtree reports. The 1.2 inflation covers a descendant's *cap*
    # here, not merely its vertices, at every depth these trees use — the other
    # half of the premise behind the brute-force equalities (see the header).
    for level in 0:3, id in (0, 5, 21, 87, DGG.num_cells(ORDINAL_MOCK, level) - 1)
        id < DGG.num_cells(ORDINAL_MOCK, level) || continue
        parent_cap = cell_cap(ORDINAL_MOCK, level, id)
        for leaf in cell_descendants(ORDINAL_MOCK, level, id, LEVEL)
            @test GO.UnitSpherical._contains(parent_cap, cell_cap(ORDINAL_MOCK, LEVEL, leaf))
        end
    end
end

# `max_level` is a fact about the canonical index, so the grid constructors
# reject levels the id encoding cannot hold — no cursor ever addresses one.
# `subtree_grid` and the keyword form both funnel through the inner
# constructors, which is where the check lives, so every path is covered.
# Systems that impose no bound (`max_level === nothing`) are untouched.
@testset "constructor level bounds" begin
    @test max_level(ORDINAL_MOCK) == 12
    @test DGGSGrid(ORDINAL_MOCK, 12).level == 12
    @test_throws ArgumentError DGGSGrid(ORDINAL_MOCK, 13)
    @test_throws ArgumentError DGGSGrid(ORDINAL_MOCK, -1)
    @test DGGSPartialGrid(ORDINAL_MOCK, 12, Int64[0, 1]).level == 12
    @test_throws ArgumentError DGGSPartialGrid(ORDINAL_MOCK, 13, Int64[0, 1])
    @test_throws ArgumentError subtree_grid(ORDINAL_MOCK, 0; root_level=12, leaf_level=13)

    # The two bounds this milestone verified, on the real systems.
    @test DGGSGrid(HEALPixDGGS(), 29).level == 29
    @test_throws ArgumentError DGGSGrid(HEALPixDGGS(), 30)
    @test DGGSPartialGrid(HEALPixDGGS(), 29, Int64[0, 5]).level == 29
    @test_throws ArgumentError DGGSPartialGrid(HEALPixDGGS(), 30, Int64[0, 5])
    @test DGGSGrid(IGEO7DGGS(), 19).level == 19
    @test_throws ArgumentError DGGSGrid(IGEO7DGGS(), 20)

    # Systems without a bound construct as they always have.
    @test max_level(RHEALPixDGGS()) === nothing
    @test DGGSGrid(RHEALPixDGGS(), 40).level == 40
end

# Ordinal ids at `level` are exactly `0:num_cells - 1`, and `ids` is already
# known sorted, so its two endpoints bound every entry: the range check is O(1)
# *and* complete. Without it an out-of-range id built a structurally valid grid
# and threw only later, from inside the chart math of a traversal.
@testset "partial grid id range" begin
    @test_throws ArgumentError DGGSPartialGrid(HEALPixDGGS(), 0, Int64[100, 200])
    @test_throws ArgumentError DGGSPartialGrid(HEALPixDGGS(), 0, Int64[-1, 3])
    @test_throws ArgumentError DGGSPartialGrid(ORDINAL_MOCK, 1,
        Int64[0, Int64(DGG.num_cells(ORDINAL_MOCK, 1))])

    # The valid path is untouched, empty grids included (no endpoints to check).
    @test DGGSPartialGrid(HEALPixDGGS(), 0, Int64[0, 11]).ids == Int64[0, 11]
    @test isempty(DGGSPartialGrid(HEALPixDGGS(), 0, Int64[]).ids)
    @test DGGSPartialGrid(ORDINAL_MOCK, LEVEL, IDS).ids === IDS

    # Gated on `has_ordinal_ids`: a structural-id system has no cheap complete
    # test, and the grid docstring says so rather than pretending otherwise.
    @test !has_ordinal_ids(STRUCTURAL_MOCK)
    @test DGGSPartialGrid(STRUCTURAL_MOCK, LEVEL, IDS).ids === IDS
end

# `descendant_range`'s ordinal default is pure radix arithmetic, so a
# nonexistent root id used to produce a perfectly well-formed interval of
# nonexistent cells — and `subtree_grid` turned that into a grid whose
# membership check passes trivially. Two integer comparisons close it, the
# ordinal counterpart of H3's and IGEO7's encoded-resolution guards.
@testset "subtree roots must exist" begin
    @test_throws ArgumentError subtree_grid(HEALPixDGGS(), 50; root_level=0, leaf_level=2)
    @test_throws ArgumentError subtree_grid(ORDINAL_MOCK, ROOT_COUNT; root_level=0, leaf_level=2)
    @test_throws ArgumentError subtree_grid(ORDINAL_MOCK, -1; root_level=0, leaf_level=2)

    # ...and the last valid root still builds the same 16-leaf chunk.
    chunk = subtree_grid(HEALPixDGGS(), 11; root_level=0, leaf_level=2)
    @test chunk.ids == collect(Int64, 176:191)
    @test Trees.ncells(treeify(chunk)) == 16
    @test Trees.ncells(treeify(subtree_grid(ORDINAL_MOCK, ROOT_COUNT - 1;
        root_level=0, leaf_level=2))) == RADIX^2
end

# The cursor is exported and its raw six-argument form used to be open beside
# the constructors that build only valid states: `node_level` reported a forged
# level verbatim, and `Trees.getcell`'s bounds check is computed FROM the
# window, so a forged window is not caught downstream — it *is* the downstream
# authority. Every testset in this file exercises the valid side of the same
# constructor: `treeify` and all three descent modes go through it unchanged.
@testset "cursor construction invariants" begin
    dense = DGGSGrid(ORDINAL_MOCK, 3)
    partial = DGGSPartialGrid(ORDINAL_MOCK, LEVEL, IDS)
    total = Int(DGG.num_cells(ORDINAL_MOCK, 3))

    # The audit's forged cursor, then each of its lies on its own.
    @test_throws ArgumentError DGG.DGGSCursor(dense, 99, Int64(-5), -7, 10^9, nothing)
    @test_throws ArgumentError DGG.DGGSCursor(dense, 4, Int64(0), 1, 1, nothing)
    @test_throws ArgumentError DGG.DGGSCursor(dense, -2, Int64(0), 1, total, nothing)
    @test_throws ArgumentError DGG.DGGSCursor(dense, 2, Int64(0), 1, total + 1, nothing)
    @test_throws ArgumentError DGG.DGGSCursor(dense, 2, Int64(0), 0, 5, nothing)
    @test_throws ArgumentError DGG.DGGSCursor(dense, 2, Int64(0), 5, 3, nothing)
    # Level -1 is the synthetic root: it owns the whole grid or it is no root.
    @test_throws ArgumentError DGG.DGGSCursor(dense, -1, Int64(0), 1, total - 1, nothing)
    # Exactly one of the window and the selection carries the leaf set.
    @test_throws ArgumentError DGG.DGGSCursor(dense, 1, Int64(0), 1, 1, [1])
    @test_throws ArgumentError DGG.DGGSCursor(partial, 1, Int64(0), 1, 2, [1])
    @test_throws ArgumentError DGG.DGGSCursor(partial, 1, Int64(0), 1, 1, [length(IDS) + 1])
    @test_throws ArgumentError DGG.DGGSCursor(partial, 1, Int64(0), 1, length(IDS) + 1, nothing)

    # Every form the descent modes actually build still constructs — including
    # the dense "no interval" marker and a partial node's empty window.
    @test DGG.DGGSCursor(dense, -1, Int64(0), 1, total, nothing).level == -1
    @test DGG.DGGSCursor(dense, 3, Int64(7), 8, 8, nothing).first_index == 8
    @test !DGG._has_interval(DGG.DGGSCursor(dense, 1, Int64(2), 0, -1, nothing))
    @test Trees.ncells(DGG.DGGSCursor(partial, LEVEL, IDS[1], 1, 0, nothing)) == 0
    @test Trees.ncells(DGG.DGGSCursor(partial, 2, Int64(0), 1, 3, [4, 5, 6])) == 3
    # A forged state is now unreachable, so the descent's own children are
    # bitwise what they were: same node keys as the traversal above builds.
    @test node_key(first(STI.getchild(treeify(dense)))) ==
        node_key(DGG.DGGSCursor(dense, 0, Int64(0), 1, RADIX^3, nothing))
end

@testset "treeify and manifold" begin
    for system in MOCKS
        dense = DGGSGrid(system, 3)
        partial = DGGSPartialGrid(system, LEVEL, IDS)
        @test GOCore.best_manifold(dense) === GO.Spherical()
        @test GOCore.best_manifold(partial) === GO.Spherical()

        dense_tree = treeify(dense)
        partial_tree = treeify(partial)
        @test dense_tree isa DGG.DGGSCursor
        @test partial_tree isa DGG.DGGSCursor
        @test GOCore.best_manifold(dense_tree) === GO.Spherical()
        # Pass-through: treeifying a cursor is the identity.
        @test treeify(dense_tree) === dense_tree
        @test treeify(partial_tree) === partial_tree
        # `treeify` must not collect or reorder the ids — leaf index i has to
        # stay position i of the lookup vector the grid was built from.
        @test partial_tree.grid.ids === IDS

        # Roots: synthetic level -1 whole-sphere node.
        @test dense_tree.level == -1
        @test partial_tree.level == -1
        @test STI.node_extent(dense_tree) == DGG.full_sphere_extent()
        @test STI.node_extent(partial_tree) == DGG.full_sphere_extent()
        @test STI.nchild(dense_tree) == ROOT_COUNT

        # The one-argument form is `Trees.treeify` resolving through
        # `best_manifold`, not a second function — so the manifold-explicit
        # call a `Regridder` makes lands on the same methods.
        @test treeify === Trees.treeify
        explicit = Trees.treeify(GO.Spherical(), dense)
        @test node_key(dense_tree) == node_key(explicit)
        @test Trees.treeify(GO.Spherical(), partial).grid.ids === partial_tree.grid.ids
    end

    # `ncells` / `getcell` are re-exported too, so a tree consumer needs
    # neither `Trees` nor `ConservativeRegridding` in scope.
    @test ncells === Trees.ncells
    @test getcell === Trees.getcell
    tree = treeify(DGGSGrid(ORDINAL_MOCK, 2))
    @test ncells(tree) == Int(DGG.num_cells(ORDINAL_MOCK, 2))
    @test ring_points(getcell(tree, 3)) ==
        DGG.cell_boundary(ORDINAL_MOCK, 2, ordinal_to_cell(ORDINAL_MOCK, 2, 3); closed=true)
end

# The `(level, id)` a node stands for is the pair every kernel operation takes,
# so it is public; the index bookkeeping beside it is not.
@testset "node accessors" begin
    tree = treeify(DGGSPartialGrid(ORDINAL_MOCK, LEVEL, IDS))
    @test node_level(tree) == -1
    @test node_id(tree) == 0            # synthetic root: level is what tells you
    walk_nodes(tree) do node
        @test node_level(node) == node.level
        @test node_id(node) == node.id
        node.level < 0 && return nothing
        # The pair is directly usable as kernel arguments.
        @test cell_cap(ORDINAL_MOCK, node_level(node), node_id(node)) isa
            GO.UnitSpherical.SphericalCap{Float64}
        return nothing
    end
    chunk = treeify(subtree_grid(ORDINAL_MOCK, 5; root_level=1, leaf_level=4))
    @test (node_level(chunk), node_id(chunk)) == (1, 5)
    leaf = descend_to_leaf(chunk)
    @test node_level(leaf) == 4
    @test cell_parent(ORDINAL_MOCK, 4, node_id(leaf), 1) == 5
end

# The public replacement for `Base.Fix1(GO.UnitSpherical._intersects, cap)`.
@testset "intersects_cap" begin
    cap = QUERY_CAPS[1]
    other = cell_cap(ORDINAL_MOCK, LEVEL, IDS[1])
    @test intersects_cap(cap, other) == GO.UnitSpherical._intersects(cap, other)
    predicate = intersects_cap(cap)
    @test predicate(other) == intersects_cap(cap, other)
    # Symmetric, and true of a cap against itself.
    @test intersects_cap(cap, cap)
    @test intersects_cap(other, cap) == intersects_cap(cap, other)

    # ...and it is the predicate the traversals take, answering exactly what
    # the hand-rolled closure did.
    for system in MOCKS, query_cap in QUERY_CAPS
        tree = treeify(DGGSPartialGrid(system, LEVEL, IDS))
        @test sort!(STI.query(tree, intersects_cap(query_cap))) ==
            brute_force(system, LEVEL, IDS, query_cap)
    end
end

@testset "SpatialTreeInterface contract" begin
    @test STI.isspatialtree(DGG.DGGSCursor)
    for system in MOCKS, grid in (DGGSGrid(system, 3),
            DGGSPartialGrid(system, LEVEL, IDS),
            DGGSPartialGrid(system, LEVEL, IDS; bucket_size=8))
        tree = treeify(grid)
        @test STI.isspatialtree(tree)
        @test !STI.isleaf(tree)

        children = collect(STI.getchild(tree))
        @test length(children) == STI.nchild(tree)
        @test all(child -> child isa DGG.DGGSCursor, children)
        @test [node_key(STI.getchild(tree, i)) for i in 1:STI.nchild(tree)] ==
            node_key.(children)
        @test_throws BoundsError STI.getchild(tree, 0)
        @test_throws BoundsError STI.getchild(tree, STI.nchild(tree) + 1)

        # Node extents are spherical caps everywhere; leaf lists are
        # materialized vectors of (index, cap) pairs, never lazy or tuples.
        leaf = descend_to_leaf(tree)
        @test STI.node_extent(leaf) isa GO.UnitSpherical.SphericalCap{Float64}
        entries = STI.child_indices_extents(leaf)
        @test entries isa Vector{Tuple{Int,GO.UnitSpherical.SphericalCap{Float64}}}
        @test length(entries) == Trees.ncells(leaf)
        @test_throws ArgumentError STI.child_indices_extents(tree)

        # Structural invariants over the whole tree.
        node_count = 0
        leaf_total = 0
        walk_nodes(tree) do node
            node_count += 1
            @test STI.node_extent(node) isa GO.UnitSpherical.SphericalCap{Float64}
            if STI.isleaf(node)
                leaf_total += Trees.ncells(node)
            else
                @test STI.nchild(node) > 0
                # A node's leaf count is exactly its children's.
                @test sum(Trees.ncells(child) for child in STI.getchild(node)) ==
                    Trees.ncells(node)
            end
        end
        @test node_count > 1
        @test leaf_total == Trees.ncells(tree)
    end
end

@testset "dense tree indexing" begin
    for system in MOCKS
        tree = treeify(DGGSGrid(system, 3))
        total = Int(DGG.num_cells(system, 3))
        @test Trees.ncells(tree) == total

        indices = all_leaf_indices(tree)
        # Dense leaves come out in ordinal order, exactly once each.
        @test indices == collect(1:total)

        for i in (1, 2, 57, total)
            id = ordinal_to_cell(system, 3, i)
            @test ring_points(Trees.getcell(tree, i)) ==
                DGG.cell_boundary(system, 3, id; closed=true)
        end
        @test_throws BoundsError Trees.getcell(tree, 0)
        @test_throws BoundsError Trees.getcell(tree, total + 1)
        @test length(collect(Trees.getcell(tree))) == total
    end
end

@testset "dense ordinal contiguity" begin
    # The identity the dense O(1) leaf interval rests on: with
    # `has_descendant_ranges`, a node's interval width is its subtree size, so
    # `descendant_range`, `cell_to_ordinal` and `subtree_leaf_count` all agree.
    tree = treeify(DGGSGrid(ORDINAL_MOCK, 4))
    checked = 0
    walk_nodes(tree) do node
        node.level < 0 && return nothing
        width = node.last_index - node.first_index + 1
        @test width == subtree_leaf_count(ORDINAL_MOCK, node.level, node.id, 4)
        @test width == Trees.ncells(node)
        lo, hi = descendant_range(ORDINAL_MOCK, node.level, node.id, 4)
        @test node.first_index == cell_to_ordinal(ORDINAL_MOCK, 4, lo)
        @test node.last_index == cell_to_ordinal(ORDINAL_MOCK, 4, hi)
        checked += 1
        return nothing
    end
    @test checked == sum(ROOT_COUNT * RADIX^level for level in 0:4)

    # Without the trait, internal dense nodes carry no interval and fall back
    # to `subtree_leaf_count`; the leaf indices must still be the ordinals.
    fallback = treeify(DGGSGrid(STRUCTURAL_MOCK, 3))
    internal = first(STI.getchild(fallback))
    @test internal.first_index == 0
    @test Trees.ncells(internal) == subtree_leaf_count(STRUCTURAL_MOCK, 0, internal.id, 3)
    @test all_leaf_indices(fallback) == collect(1:Int(DGG.num_cells(STRUCTURAL_MOCK, 3)))
end

@testset "partial tree leaf ordering" begin
    for system in MOCKS
        tree = treeify(DGGSPartialGrid(system, LEVEL, IDS))
        @test Trees.ncells(tree) == length(IDS)
        # Leaf indices are positions in `ids`, emitted in `ids` order — this is
        # what lets a Regridder matrix line up with a DimensionalData lookup.
        @test all_leaf_indices(tree) == collect(1:length(IDS))
        for i in (1, 2, 150, length(IDS))
            @test ring_points(Trees.getcell(tree, i)) ==
                DGG.cell_boundary(system, LEVEL, IDS[i]; closed=true)
        end
        @test_throws BoundsError Trees.getcell(tree, length(IDS) + 1)

        # A one-cell grid degenerates gracefully.
        single = treeify(DGGSPartialGrid(system, LEVEL, IDS[3:3]))
        @test Trees.ncells(single) == 1
        @test all_leaf_indices(single) == [1]

        # So does an empty one: the root is a leaf with no entries.
        empty_tree = treeify(DGGSPartialGrid(system, LEVEL, Int64[]))
        @test Trees.ncells(empty_tree) == 0
        @test STI.isleaf(empty_tree)
        @test isempty(STI.child_indices_extents(empty_tree))
        @test isempty(all_leaf_indices(empty_tree))
    end
end

@testset "bucket_size" begin
    for system in MOCKS, bucket in (1, 4, 32)
        tree = treeify(DGGSPartialGrid(system, LEVEL, IDS; bucket_size=bucket))
        @test all_leaf_indices(tree) == collect(1:length(IDS))
        walk_nodes(tree) do node
            STI.isleaf(node) || return nothing
            # A leaf is either a single cell or a bucket under the limit.
            @test node.level == LEVEL || Trees.ncells(node) <= bucket
            return nothing
        end
    end

    # Buckets stop descent: fewer nodes than the unbucketed tree.
    count_nodes(tree) = (n = 0; walk_nodes(_ -> (n += 1), tree); n)
    plain = count_nodes(treeify(DGGSPartialGrid(ORDINAL_MOCK, LEVEL, IDS)))
    bucketed = count_nodes(treeify(DGGSPartialGrid(ORDINAL_MOCK, LEVEL, IDS; bucket_size=16)))
    @test bucketed < plain
end

@testset "subtree rooting" begin
    # Keyword-only levels: there is no positional form to mis-order, and the
    # two levels cannot be supplied by accident.
    @test_throws UndefKeywordError subtree_grid(ORDINAL_MOCK, 5; root_level=1)
    @test_throws UndefKeywordError subtree_grid(ORDINAL_MOCK, 5; leaf_level=4)
    @test_throws MethodError subtree_grid(ORDINAL_MOCK, 1, 5, 4)
    @test subtree_grid(ORDINAL_MOCK, 5; root_level=1, leaf_level=4,
        bucket_size=8).bucket_size == 8

    for system in MOCKS
        grid = subtree_grid(system, 5; root_level=1, leaf_level=4)
        tree = treeify(grid)
        expected = cell_descendants(system, 1, 5, 4)
        @test grid.root_level == 1
        @test grid.root_id == 5
        @test Trees.ncells(tree) == subtree_leaf_count(system, 1, 5, 4)
        @test Trees.ncells(tree) == length(expected)

        # The cursor roots at the cell, not at the sphere, and its extent is
        # that cell's own cap rather than the full sphere: a chunk stores its
        # whole subtree, and the hierarchy bounds a whole subtree in O(1).
        @test tree.level == 1
        @test tree.id == 5
        @test STI.node_extent(tree) == cell_cap(system, 1, 5)
        @test STI.node_extent(tree) != DGG.full_sphere_extent()
        # ...and it really does bound the chunk: every vertex of every leaf.
        @test all(expected) do id
            all(point -> GO.UnitSpherical._contains(STI.node_extent(tree), point),
                DGG.cell_boundary(system, 4, id))
        end

        # Leaves are numbered 1:subtree_leaf_count in ascending id order.
        @test all_leaf_indices(tree) == collect(1:length(expected))
        @test ring_points(Trees.getcell(tree, 7)) ==
            DGG.cell_boundary(system, 4, expected[7]; closed=true)

        # A subtree rooted at the leaf level is a single-cell leaf.
        leafroot = treeify(subtree_grid(system, 300; root_level=4, leaf_level=4))
        @test Trees.ncells(leafroot) == 1
        @test STI.isleaf(leafroot)
        @test STI.child_indices_extents(leafroot) ==
            [(1, cell_cap(system, 4, Int64(300)))]
    end
end

@testset "query correctness vs brute force" begin
    for system in MOCKS
        # Dense: leaf index space is the ordinals.
        dense_tree = treeify(DGGSGrid(system, 3))
        dense_ids = collect(Int64, 0:(DGG.num_cells(system, 3) - 1))
        for cap in QUERY_CAPS
            @test query(dense_tree, cap) == brute_force(system, 3, dense_ids, cap)
        end

        # Partial, with and without buckets, dense-ish and sparse.
        for ids in (IDS, SPARSE_IDS), bucket in (0, 5)
            tree = treeify(DGGSPartialGrid(system, LEVEL, ids; bucket_size=bucket))
            for cap in QUERY_CAPS
                @test query(tree, cap) == brute_force(system, LEVEL, ids, cap)
            end
        end

        # Subtree-rooted chunks answer in chunk-local indices.
        chunk = subtree_grid(system, 5; root_level=1, leaf_level=4)
        chunk_tree = treeify(chunk)
        for cap in QUERY_CAPS
            @test query(chunk_tree, cap) == brute_force(system, 4, chunk.ids, cap)
        end
    end

    # The two descent modes agree exactly, on identical geometry.
    for ids in (IDS, SPARSE_IDS), cap in QUERY_CAPS
        ordinal_tree = treeify(DGGSPartialGrid(ORDINAL_MOCK, LEVEL, ids))
        structural_tree = treeify(DGGSPartialGrid(STRUCTURAL_MOCK, LEVEL, ids))
        @test query(ordinal_tree, cap) == query(structural_tree, cap)
    end
end

@testset "dual traversal" begin
    # What `ConservativeRegridding` actually runs: candidate pairs from a dual
    # depth-first search over two DGGS cursors.
    dense_tree = treeify(DGGSGrid(ORDINAL_MOCK, 2))
    partial_tree = treeify(DGGSPartialGrid(STRUCTURAL_MOCK, LEVEL, SPARSE_IDS))
    dense_ids = collect(Int64, 0:(DGG.num_cells(ORDINAL_MOCK, 2) - 1))

    pairs = Tuple{Int,Int}[]
    STI.dual_depth_first_search(GO.UnitSpherical._intersects, dense_tree, partial_tree) do i, j
        push!(pairs, (i, j))
    end

    expected = Tuple{Int,Int}[]
    for i in eachindex(dense_ids), j in eachindex(SPARSE_IDS)
        left = cell_cap(ORDINAL_MOCK, 2, dense_ids[i])
        right = cell_cap(STRUCTURAL_MOCK, LEVEL, SPARSE_IDS[j])
        GO.UnitSpherical._intersects(left, right) && push!(expected, (i, j))
    end

    @test !isempty(expected)
    @test sort(pairs) == sort(expected)
end

@testset "should_parallelize policy" begin
    chunks = 32 # DGG.PARALLELIZE_CHUNKS_PER_THREAD
    @test DGG.PARALLELIZE_CHUNKS_PER_THREAD == chunks

    for grid in (DGGSGrid(ORDINAL_MOCK, 4), DGGSPartialGrid(ORDINAL_MOCK, LEVEL, IDS),
            subtree_grid(ORDINAL_MOCK, 5; root_level=1, leaf_level=4))
        tree = treeify(grid)
        total = grid isa DGGSGrid ? Int(DGG.num_cells(grid.system, grid.level)) :
                length(grid.ids)
        threshold = max(1, total ÷ (Threads.nthreads() * chunks))

        # The method that answers must be ours, not ConservativeRegridding's
        # `::Any` cap-area fallback — dispatching on the node type is the whole
        # point of carrying the grid in the cursor.
        method = which(Trees.should_parallelize,
            Tuple{typeof(tree),GO.UnitSpherical.SphericalCap{Float64}})
        @test method.module === DiscreteGlobalGrids
        @test method.sig.parameters[2] !== Any

        # Leaf-count formula, on the root and on every node of a walk.
        @test Trees.should_parallelize(tree, STI.node_extent(tree)) ==
            (Trees.ncells(tree) <= threshold)
        walk_nodes(tree) do node
            extent = STI.node_extent(node)
            @test Trees.should_parallelize(node, extent) == (Trees.ncells(node) <= threshold)
            return nothing
        end

        # The root of a grid this size is above the threshold; a single cell is
        # always below it.
        @test total > threshold
        @test !Trees.should_parallelize(tree, STI.node_extent(tree))
        leaf = descend_to_leaf(tree)
        @test Trees.should_parallelize(leaf, STI.node_extent(leaf))
    end
end

@testset "ncells/getcell consistency" begin
    for system in MOCKS, grid in (DGGSGrid(system, 3),
            DGGSPartialGrid(system, LEVEL, IDS),
            DGGSPartialGrid(system, LEVEL, IDS; bucket_size=6),
            subtree_grid(system, 5; root_level=1, leaf_level=4))
        tree = treeify(grid)
        system = grid.system
        level = grid.level
        # Node-level indexing is a window on the tree-level index space: the
        # k-th cell of a leaf is the cell at the index that leaf reports.
        walk_nodes(tree) do node
            STI.isleaf(node) || return nothing
            entries = STI.child_indices_extents(node)
            @test length(entries) == Trees.ncells(node)
            for k in eachindex(entries)
                index, extent = entries[k]
                @test 1 <= index <= Trees.ncells(tree)
                @test ring_points(Trees.getcell(node, k)) ==
                    ring_points(Trees.getcell(tree, index))
                @test extent == cell_cap(system, level, leaf_id_at(tree, index))
            end
            return nothing
        end
        # An internal node's window is likewise a slice of the tree's.
        internal = first(STI.getchild(tree))
        under = all_leaf_indices(internal)
        @test length(under) == Trees.ncells(internal)
        @test ring_points(Trees.getcell(internal, 1)) ==
            ring_points(Trees.getcell(tree, first(under)))
    end
end

@testset "node extents come from the kernel" begin
    for system in MOCKS
        # A partial internal node reports its own cell's cap — O(1), straight
        # from the hierarchy — unless it stores a *proper* subset of its
        # subtree and few enough leaves that bounding them directly is bounded
        # work (`STORED_UNION_CAP_LIMIT`); that union is the tighter cap a
        # sparse chunk exists for. Dense internal nodes bound the whole
        # subtree, and leaves use the cell cap.
        @test !has_exact_subtree_cap(system)
        partial = DGGSPartialGrid(system, LEVEL, IDS)
        tree = treeify(partial)
        union_nodes = 0
        cell_nodes = 0
        walk_nodes(tree) do node
            node.level < 0 && return nothing
            stored = [IDS[i] for i in all_leaf_indices(node)]
            extent = STI.node_extent(node)
            if node.level == LEVEL
                @test extent == cell_cap(system, LEVEL, node.id)
            elseif length(stored) <= DGG.STORED_UNION_CAP_LIMIT &&
                   length(stored) < subtree_leaf_count(system, node.level, node.id, LEVEL)
                @test extent == DGG.cells_cap(system, LEVEL, stored)
                union_nodes += 1
            else
                @test extent == cell_cap(system, node.level, node.id)
                cell_nodes += 1
            end
            # Whichever cap it is, it bounds what the node owns — the one
            # thing a node extent may never get wrong, since a leaf outside it
            # is silently dropped from every traversal, not reported.
            @test all(stored) do id
                all(point -> GO.UnitSpherical._contains(extent, point),
                    DGG.cell_boundary(system, LEVEL, id))
            end
            return nothing
        end
        # Both halves of the rule are live on this grid (300 of 2,048 cells).
        @test union_nodes > 0
        @test cell_nodes > 0

        dense = treeify(DGGSGrid(system, 3))
        walk_nodes(dense) do node
            node.level < 0 && return nothing
            expected = node.level == 3 ? cell_cap(system, 3, node.id) :
                       DGG.subtree_cap(system, node.level, node.id, 3)
            @test STI.node_extent(node) == expected
            return nothing
        end
    end
end

@testset "has_exact_subtree_cap" begin
    # Third mock configuration: identical grid and geometry, `subtree_cap`
    # declared O(1) and geographically tight. Partial internal nodes must then
    # take that cap instead of the O(stored) union over their stored ids —
    # which is what restores the O(1) node extents the old per-system HEALPix
    # tree had.
    @test has_exact_subtree_cap(EXACT_CAP_MOCK)
    @test !has_exact_subtree_cap(ORDINAL_MOCK)
    @test !has_exact_subtree_cap(STRUCTURAL_MOCK)
    # The trait must not change any other trait's answer.
    @test has_ordinal_ids(EXACT_CAP_MOCK) == has_ordinal_ids(ORDINAL_MOCK)
    @test has_descendant_ranges(EXACT_CAP_MOCK) == has_descendant_ranges(ORDINAL_MOCK)

    for ids in (IDS, SPARSE_IDS), bucket in (0, 5)
        exact_tree = treeify(DGGSPartialGrid(EXACT_CAP_MOCK, LEVEL, ids; bucket_size=bucket))
        plain_tree = treeify(DGGSPartialGrid(ORDINAL_MOCK, LEVEL, ids; bucket_size=bucket))

        differs = false
        walk_nodes(exact_tree) do node
            node.level < 0 && return nothing
            extent = STI.node_extent(node)
            if node.level == LEVEL
                # Leaf extents are untouched by the trait.
                @test extent == cell_cap(EXACT_CAP_MOCK, LEVEL, node.id)
            else
                @test extent == DGG.subtree_cap(EXACT_CAP_MOCK, node.level, node.id, LEVEL)
                stored = [ids[i] for i in all_leaf_indices(node)]
                differs |= extent != DGG.cells_cap(EXACT_CAP_MOCK, LEVEL, stored)
            end
            return nothing
        end
        # The branch is really taken: at least one internal node's extent is
        # not the union cap it would have been with the trait off.
        @test differs

        # Same tree shape and same leaf indexing — only the extents move.
        @test Trees.ncells(exact_tree) == Trees.ncells(plain_tree)
        @test all_leaf_indices(exact_tree) == all_leaf_indices(plain_tree)
        exact_keys = Tuple{Int,Int64,Int}[]
        plain_keys = Tuple{Int,Int64,Int}[]
        walk_nodes(node -> push!(exact_keys, node_key(node)), exact_tree)
        walk_nodes(node -> push!(plain_keys, node_key(node)), plain_tree)
        @test exact_keys == plain_keys

        # And the traversal still answers exactly: `subtree_cap` covers the
        # whole subtree, hence every stored leaf's cap.
        for cap in QUERY_CAPS
            @test query(exact_tree, cap) == brute_force(EXACT_CAP_MOCK, LEVEL, ids, cap)
            @test query(exact_tree, cap) == query(plain_tree, cap)
        end
    end

    # Subtree-rooted grids root on the cell's own cap too.
    chunk = subtree_grid(EXACT_CAP_MOCK, 5; root_level=1, leaf_level=4)
    chunk_tree = treeify(chunk)
    @test STI.node_extent(chunk_tree) == DGG.subtree_cap(EXACT_CAP_MOCK, 1, 5, 4)
    for cap in QUERY_CAPS
        @test query(chunk_tree, cap) == brute_force(EXACT_CAP_MOCK, 4, chunk.ids, cap)
    end

    # A non-rooted partial root stays the whole sphere, trait or not.
    @test STI.node_extent(treeify(DGGSPartialGrid(EXACT_CAP_MOCK, LEVEL, IDS))) ==
        DGG.full_sphere_extent()

    # The dense path is untouched: it already used `subtree_cap`.
    dense_extents(system) = (caps = GO.UnitSpherical.SphericalCap{Float64}[];
        walk_nodes(node -> push!(caps, STI.node_extent(node)),
            treeify(DGGSGrid(system, 3))); caps)
    @test dense_extents(EXACT_CAP_MOCK) == dense_extents(ORDINAL_MOCK)
end

@testset "type stability" begin
    for grid in (DGGSGrid(ORDINAL_MOCK, 3), DGGSGrid(STRUCTURAL_MOCK, 3),
            DGGSPartialGrid(ORDINAL_MOCK, LEVEL, IDS),
            DGGSPartialGrid(STRUCTURAL_MOCK, LEVEL, IDS),
            DGGSPartialGrid(EXACT_CAP_MOCK, LEVEL, IDS))
        tree = treeify(grid)
        child = first(STI.getchild(tree))
        # The node type is a fixed point of descent: one concrete cursor type
        # for the whole traversal, no per-level widening.
        @test typeof(child) === typeof(first(STI.getchild(child)))
        @inferred STI.isleaf(child)
        @inferred STI.nchild(child)
        @inferred STI.getchild(child)
        @inferred STI.node_extent(child)
        @inferred Trees.ncells(child)
        @inferred Trees.should_parallelize(child, STI.node_extent(child))
        leaf = descend_to_leaf(tree)
        @inferred STI.child_indices_extents(leaf)
        @inferred Trees.getcell(leaf, 1)
    end
end

end # module GenericTreeTests
