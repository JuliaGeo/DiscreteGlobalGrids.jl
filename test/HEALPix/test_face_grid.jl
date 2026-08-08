module HealpixFaceGridTestSuite

# Tests for `src/HEALPix/face_grid.jl`: dense HEALPix grids (all `12 nside²`
# pixels of one resolution) as spatial trees, under a swappable data ordering.
#
# Three things are being pinned down here, in rising order of load-bearingness:
#
# 1. *Orderings are bijections.* `data_index` / `lattice_index` must be exact
#    mutual inverses over `1:12nside²`, because that pair is the whole
#    alignment contract: column `j` of a `Regridder` is data position `j`.
#
# 2. *Node extents contain their subtree's geometry.* This is what makes tree
#    pruning sound. It is asserted exhaustively — every node of every face
#    tree, against both the block's lattice vertices and a dense sampling of
#    the block's chart rectangle — and it is the same invariant
#    `test_healpix_kernel.jl` asserts for the id-hierarchy path ("HEALPix
#    subtree caps contain descendants": descendant *vertices* inside the
#    parent cap).
#
#    Note what is deliberately NOT claimed: a node's extent does not contain
#    its descendants' *caps*. It cannot, for any extent that tightly bounds
#    only the node's own geometry — a leaf's circumscribed cap sticks out past
#    the block boundary by roughly its own radius, whatever the parent does.
#    So the query comparison below is a two-sided sandwich rather than an
#    equality: `query ⊇ {leaves that really meet the cap}` (soundness, the
#    load-bearing half) and `query ⊆ {leaves whose leaf extent meets the cap}`
#    (no spurious extras). The generic-tree suite gets an equality only
#    because its mock deliberately folds a cap-of-caps at every internal node,
#    which is O(subtree) and throws away the O(1) block cap this layer exists
#    for.
#
# 3. *The face-grid path and the id-hierarchy path agree.* At `nside = 2^level`
#    with `NestedOrder`, column identity is the same on both sides (nested data
#    position `p + 1` == DGGS ordinal `p + 1`), so their `Regridder` matrices
#    are directly comparable entry by entry.
#
# Regridders here are built on the *unit* sphere (`GO.Spherical(radius = 1)`)
# rather than the default Earth-radius manifold, so that conservation reads as
# `4π` and per-entry tolerances are absolute numbers rather than fractions of
# 5.1e14 m².

using Test
using Printf
using Random
import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
import GeometryOps: SpatialTreeInterface as STI
import ConservativeRegridding as CR
import ConservativeRegridding: Trees
import SparseArrays

using DiscreteGlobalGrids
using DiscreteGlobalGrids.HEALPix
using DiscreteGlobalGrids.HEALPix: FaceChartGrid, HealpixFaceRoot,
    num_pixels, pixel_polygon, data_index, lattice_index, validate_ordering,
    xyf_to_point, pixel_corners, xyf_to_ring, xyf_to_nested, nested_to_xyf,
    nested_to_ring

const US = GO.UnitSpherical
# Unit-radius sphere: areas come out in steradians, so conservation is `4π`.
const UNIT = GO.Spherical(radius=1.0)

# Recorded so the numbers land in the test log (and in the milestone report).
# `-Inf` rather than `0.0` as the neutral element: several of the quantities
# below are *negative* by design (a cap overhang inside the cap), and a `0.0`
# floor would quietly report them as zero.
const MEASURED = Dict{String,Float64}()
record!(key, value) = (MEASURED[key] = max(get(MEASURED, key, -Inf), value))

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

cross3(a, b) = (a[2] * b[3] - a[3] * b[2],
                a[3] * b[1] - a[1] * b[3],
                a[1] * b[2] - a[2] * b[1])

# Signed area of the 4-gon seen from *outside* the sphere; positive ⇔ CCW.
# Same measure as `test_chart.jl` uses on `pixel_corners`.
function ccw_measure(corners)
    acc = (0.0, 0.0, 0.0)
    n = length(corners)
    for i in 1:n
        acc = acc .+ cross3(corners[i], corners[i % n + 1])
    end
    outward = reduce((a, b) -> a .+ Tuple(b), corners; init=(0.0, 0.0, 0.0))
    return sum(acc .* outward)
end

ring_points(poly) = collect(GI.getpoint(GI.getexterior(poly)))
# `getcell` closes the ring, so the four distinct corners are points 1:4.
open_ring(poly) = ring_points(poly)[1:4]

"Depth-first walk over every node of a face cursor (the node itself included)."
function walk_nodes(f::Function, node)
    f(node)
    STI.isleaf(node) && return nothing
    for child in STI.getchild(node)
        walk_nodes(f, child)
    end
    return nothing
end

"""
Every `(data index, leaf extent)` pair the traversal can ever test, harvested
from `STI.child_indices_extents` — i.e. the literal leaf-level predicate input,
not a reconstruction of it.
"""
function leaf_extents(root)
    exts = Dict{Int,US.SphericalCap{Float64}}()
    for f in 1:STI.nchild(root)
        walk_nodes(STI.getchild(root, f)) do node
            STI.isleaf(node) || return nothing
            for (index, extent) in STI.child_indices_extents(node)
                exts[index] = extent
            end
            return nothing
        end
    end
    return exts
end

"Random unit vector / random cap, from a seeded RNG."
function random_point(rng)
    v = randn(rng, 3)
    v ./= sqrt(sum(abs2, v))
    return GO.UnitSphericalPoint(v[1], v[2], v[3])
end
random_cap(rng) = US.SphericalCap(random_point(rng), 0.01 + 1.2 * rand(rng))

"Is any point of pixel `(ix, iy, face)` inside `cap`? (dense chart sampling)"
function pixel_meets_cap(cap, ix, iy, face, nside, samples=4)
    for a in 0:samples, b in 0:samples
        p = xyf_to_point((ix + a / samples) / nside, (iy + b / samples) / nside, face)
        US.spherical_distance(cap.point, p) <= cap.radius && return true
    end
    return false
end

# --------------------------------------------------------------------------
# 1. Construction and validation
# --------------------------------------------------------------------------

@testset "HealpixFaceSpace / HealpixFaceGrid construction" begin
    @test_throws ArgumentError HealpixFaceSpace(0)
    @test_throws ArgumentError HealpixFaceSpace(-3)
    @test HealpixFaceSpace(1).nside == 1
    @test num_pixels(HealpixFaceSpace(3)) == 108
    @test num_pixels(HealpixFaceSpace(4)) == 192

    # Upper bound: `12 * nside^2` and `ring_first` overflow Int64 past 2^29
    # (the same bound as `HEALPixDGGS`'s `max_level == 29`).
    @test HealpixFaceSpace(2^29).nside == 2^29
    @test_throws ArgumentError HealpixFaceSpace(2^30)

    # The nested index is a Morton code, so it exists only on a 2^k x 2^k face.
    # It must be refused at *construction*, not from inside a traversal.
    @test_throws ArgumentError HealpixFaceGrid(3; ordering=NestedOrder())
    @test_throws ArgumentError HealpixFaceGrid(5; ordering=NestedOrder())
    @test_throws ArgumentError HealpixFaceGrid(HealpixFaceSpace(6), NestedOrder())
    @test_throws ArgumentError HealpixFaceGrid(0)
    @test_throws ArgumentError HealpixFaceGrid(0; ordering=NestedOrder())

    # RING carries no such restriction — that is the point of this layer.
    for nside in (1, 2, 3, 4, 5, 7)
        grid = HealpixFaceGrid(nside)
        @test grid.ordering isa RingOrder            # the default
        @test grid.space.nside == nside
        @test num_pixels(grid) == 12nside^2
    end
    @test HealpixFaceGrid(4; ordering=NestedOrder()).ordering isa NestedOrder
    @test RingOrder() isa AbstractHealpixOrdering
    @test NestedOrder() isa AbstractHealpixOrdering
    # The optional third contract method defaults to "accept everything". It
    # takes the *space*, not a loose `nside`: the resolution reaches this layer
    # only through the type that has already checked it.
    @test validate_ordering(RingOrder(), HealpixFaceSpace(3)) === nothing
    @test validate_ordering(NestedOrder(), HealpixFaceSpace(4)) === nothing
    @test_throws ArgumentError validate_ordering(NestedOrder(), HealpixFaceSpace(3))

    @test sprint(show, HealpixFaceSpace(4)) == "HealpixFaceSpace(4)"
    @test occursin("nside=4", sprint(show, HealpixFaceGrid(4)))
end

# `HealpixFaceRoot` is independently constructible — it is the REPL/test
# convenience form — so it must not be a way *around* the checks
# `HealpixFaceGrid` runs. Its resolution is a `HealpixFaceSpace` (so `nside` is
# checked by construction) and it re-runs both ordering checks itself.
@testset "HealpixFaceRoot construction validates like the grid" begin
    # (a) `nside`, through the space the convenience form builds.
    @test_throws ArgumentError HealpixFaceRoot(0)
    @test_throws ArgumentError HealpixFaceRoot(-5, RingOrder())
    @test_throws ArgumentError HealpixFaceRoot(2^30, RingOrder())
    @test HealpixFaceRoot(2^29, RingOrder()).space.nside == 2^29

    # (b) another system's ordering: an `ArgumentError` here rather than a
    # `MethodError` from inside a traversal.
    @test_throws ArgumentError HealpixFaceRoot(4, DiscreteGlobalGrids.S2.RowMajorOrder())
    @test_throws ArgumentError HealpixFaceRoot(
        GO.Spherical(), HealpixFaceSpace(4), DiscreteGlobalGrids.ISEA4R.MortonOrder())

    # (c) `validate_ordering`, which now runs at root construction too.
    @test_throws ArgumentError HealpixFaceRoot(3, NestedOrder())
    @test_throws ArgumentError HealpixFaceRoot(GO.Spherical(), HealpixFaceSpace(5), NestedOrder())
    @test HealpixFaceRoot(4, NestedOrder()).ordering isa NestedOrder
end

@testset "treeify and tree toplevel" begin
    for (nside, ordering) in ((4, NestedOrder()), (5, RingOrder()), (3, RingOrder()))
        grid = HealpixFaceGrid(nside; ordering)
        # One-argument `treeify` resolves through `best_manifold`.
        @test GOCore.best_manifold(grid) == GO.Spherical()
        root = treeify(grid)
        @test root isa HealpixFaceRoot
        @test root.space == HealpixFaceSpace(nside)
        @test root.space.nside == nside
        @test root.ordering === ordering
        @test GOCore.best_manifold(root) == GO.Spherical()
        # Explicit-manifold form, and idempotent passthrough on the tree.
        @test treeify(GO.Spherical(), grid) isa HealpixFaceRoot
        @test treeify(root) === root
        @test treeify(GO.Spherical(), root) === root
        @test Trees.treeify(UNIT, root) === root
        # The direct constructor defaults the manifold for REPL use.
        @test HealpixFaceRoot(nside, ordering) == root

        @test STI.isspatialtree(typeof(root))
        @test !STI.isleaf(root)
        @test STI.nchild(root) == 12
        @test Trees.ncells(root) == 12nside^2
        @test length(collect(STI.getchild(root))) == 12
        @test occursin("nside=$nside", sprint(show, root))

        # The root bounds nothing tighter than the whole sphere: the twelve
        # faces tile it.
        @test STI.node_extent(root).radius >= Float64(pi)

        # Children are stock quadtree cursors over per-face chart grids.
        for f in 1:12
            child = STI.getchild(root, f)
            @test child isa Trees.TopDownQuadtreeCursor{<:FaceChartGrid}
            @test child.grid.face == f - 1
            @test child.grid.space == HealpixFaceSpace(nside)
            @test Trees.ncells(child.grid, 1) == nside
            @test Trees.ncells(child.grid, 2) == nside
            @test GOCore.manifold(child.grid) == GO.Spherical()
        end

        @test_throws BoundsError Trees.getcell(root, 0)
        @test_throws BoundsError Trees.getcell(root, 12nside^2 + 1)
        @test length(collect(Trees.getcell(root))) == 12nside^2
    end
end

# --------------------------------------------------------------------------
# 2. The ordering contract
# --------------------------------------------------------------------------

@testset "ordering bijections (nside = $nside)" for nside in (1, 2, 3, 4, 5, 7, 8)
    npix = 12nside^2
    space = HealpixFaceSpace(nside)
    orderings = ispow2(nside) ? (RingOrder(), NestedOrder()) : (RingOrder(),)
    for ordering in orderings
        # `lattice_index` is onto the whole lattice and `data_index` inverts it.
        seen = Set{NTuple{3,Int}}()
        for j in 1:npix
            ix, iy, face = lattice_index(ordering, space, j)
            @test 0 <= ix < nside && 0 <= iy < nside && 0 <= face <= 11
            push!(seen, (ix, iy, face))
            @test data_index(ordering, space, ix, iy, face) == j
        end
        @test length(seen) == npix

        # ... and the other direction, over the lattice.
        for face in 0:11, ix in 0:(nside - 1), iy in 0:(nside - 1)
            j = data_index(ordering, space, ix, iy, face)
            @test 1 <= j <= npix
            @test lattice_index(ordering, space, j) == (ix, iy, face)
        end
    end

    # The two shipped orderings are exactly the chart's two index maps, with
    # the 0-based/1-based conventions of `chart.jl` reconciled to data
    # positions (RING is already 1-based; NESTED is 0-based, hence the `+ 1`).
    ix, iy = min(1, nside - 1), min(2, nside - 1)
    @test data_index(RingOrder(), space, ix, iy, 5) == xyf_to_ring(ix, iy, 5, nside)
    if ispow2(nside)
        @test data_index(NestedOrder(), space, ix, iy, 5) ==
            xyf_to_nested(ix, iy, 5, nside) + 1
        @test lattice_index(NestedOrder(), space, 7) == nested_to_xyf(6, nside)
    end
end

@testset "getcell round-trips through the ordering (nside = $nside, $ordering)" for
        (nside, ordering) in ((4, NestedOrder()), (4, RingOrder()), (5, RingOrder()), (3, RingOrder()))

    space = HealpixFaceSpace(nside)
    root = treeify(HealpixFaceGrid(nside; ordering))
    npix = 12nside^2
    for j in 1:npix
        ix, iy, face = lattice_index(ordering, space, j)
        # `Trees.getcell(root, j)` is the geometry side of the alignment rule.
        @test ring_points(Trees.getcell(root, j)) ==
            ring_points(pixel_polygon(ix, iy, face, nside))
        # Closed CCW 4-gon, with the corners `chart.jl` emits, in its order.
        pts = ring_points(Trees.getcell(root, j))
        @test length(pts) == 5 && pts[1] == pts[5]
        @test Tuple(pts[1:4]) == pixel_corners(ix, iy, face, nside)
    end
    # The iterator form must agree with the indexed one, position by position.
    @test collect(Trees.getcell(root)) == [Trees.getcell(root, j) for j in 1:npix]

    # The per-face grid's own index maps target the *global* data layout.
    for face in 0:11
        g = FaceChartGrid(GO.Spherical(), space, face, ordering)
        for ix in 0:(nside - 1), iy in 0:(nside - 1)
            j = Trees.cartesian_to_linear_idx(g, CartesianIndex(ix + 1, iy + 1))
            @test j == data_index(ordering, space, ix, iy, face)
            @test Trees.linear_to_cartesian_idx(g, j) == CartesianIndex(ix + 1, iy + 1)
            @test ring_points(Trees.getcell(g, ix + 1, iy + 1)) ==
                ring_points(Trees.getcell(root, j))
            @test Trees.getvertex(g, ix + 1, iy + 1) ==
                xyf_to_point(ix / nside, iy / nside, face)
        end
    end
end

# NestedOrder is the id space `HealpixKernel.jl` wires, so a nested-ordered
# grid at `nside = 2^level` must be column-for-column the DGGS grid at `level`.
@testset "NestedOrder position j is nested id j - 1" begin
    for level in 0:3
        nside = 2^level
        root = treeify(HealpixFaceGrid(nside; ordering=NestedOrder()))
        system = HEALPixDGGS()
        for j in 1:(12nside^2)
            id = Int64(j - 1)
            @test DiscreteGlobalGrids.cell_to_ordinal(system, level, id) == j
            @test lattice_index(NestedOrder(), HealpixFaceSpace(nside), j) ==
                nested_to_xyf(id, nside)
        end
    end
end

# --------------------------------------------------------------------------
# 3. Node extents contain their subtree's geometry
#
# The soundness proof for pruning, by exhaustion. Two probes per node: the
# block's lattice vertices (the polygon corners of every pixel it owns) and a
# dense sampling of the block's chart rectangle (which also covers the
# constant-z / constant-φ bulge of the pixel edges between lattice points).
# --------------------------------------------------------------------------

@testset "node extents contain their block geometry (nside = $nside)" for nside in (3, 4, 5, 8)
    root = treeify(HealpixFaceGrid(nside; ordering=RingOrder()))
    worst_vertex = -Inf
    worst_sample = -Inf
    nodes = 0
    for f in 1:12
        cursor = STI.getchild(root, f)
        grid = cursor.grid
        walk_nodes(cursor) do node
            nodes += 1
            extent = STI.node_extent(node)
            imin, imax = extrema(node.leafranges[1]); imax += 1
            jmin, jmax = extrema(node.leafranges[2]); jmax += 1
            for i in imin:imax, j in jmin:jmax
                p = Trees.getvertex(grid, i, j)
                worst_vertex = max(worst_vertex,
                    US.spherical_distance(extent.point, p) - extent.radius)
            end
            x0, x1 = (imin - 1) / nside, (imax - 1) / nside
            y0, y1 = (jmin - 1) / nside, (jmax - 1) / nside
            for a in 0:8, b in 0:8
                p = xyf_to_point(x0 + (x1 - x0) * a / 8, y0 + (y1 - y0) * b / 8, grid.face)
                worst_sample = max(worst_sample,
                    US.spherical_distance(extent.point, p) - extent.radius)
            end
            return nothing
        end
    end
    @test nodes > 12                                   # the walk actually descended
    # Non-positive overhang: every probe sits inside the cap.
    @test worst_vertex <= 0
    @test worst_sample <= 0
    record!("node-extent overhang, lattice vertices (rad)", worst_vertex)
    record!("node-extent overhang, dense chart samples (rad)", worst_sample)

    # Leaf extents likewise contain their own pixel, and the leaf index space
    # is exactly `1:12nside²` with no index emitted twice.
    exts = leaf_extents(root)
    @test sort!(collect(keys(exts))) == collect(1:(12nside^2))
    worst_leaf = -Inf
    for (j, extent) in exts
        for p in open_ring(Trees.getcell(root, j))
            worst_leaf = max(worst_leaf,
                US.spherical_distance(extent.point, p) - extent.radius)
        end
    end
    @test worst_leaf <= 0
    record!("leaf-extent overhang, pixel corners (rad)", worst_leaf)
end

# --------------------------------------------------------------------------
# 4. Cap queries against brute force
# --------------------------------------------------------------------------

@testset "STI queries vs brute force (nside = $nside, $ordering)" for
        (nside, ordering) in ((4, NestedOrder()), (4, RingOrder()), (5, RingOrder()))

    root = treeify(HealpixFaceGrid(nside; ordering))
    npix = 12nside^2
    exts = leaf_extents(root)
    lattice = [lattice_index(ordering, HealpixFaceSpace(nside), j) for j in 1:npix]

    rng = MersenneTwister(20260805 + nside)
    pruned_total = 0
    leaf_positive_total = 0
    tested_total = 0
    for _ in 1:50
        cap = random_cap(rng)
        # Instrumented predicate: `tested_total` counts every extent the
        # traversal actually evaluates, internal blocks included. See the
        # assertion after the loop.
        answer = sort!(STI.query(root, function (extent)
            tested_total += 1
            return intersects_cap(cap, extent)
        end))

        @test allunique(answer)
        @test all(j -> 1 <= j <= npix, answer)

        # Upper bound: the traversal only ever *returns* leaves whose own leaf
        # extent passes the predicate, so no result can be outside this set.
        leaf_positives = [j for j in 1:npix if intersects_cap(cap, exts[j])]
        @test issubset(answer, leaf_positives)

        # Lower bound — the load-bearing half. Any pixel that genuinely meets
        # the cap must survive every ancestor's pruning test. (Ground truth by
        # dense chart sampling of the pixel: a sampled point inside the cap is
        # a witness that the pixel really intersects it.)
        truth = [j for j in 1:npix if pixel_meets_cap(cap, lattice[j]..., nside)]
        @test issubset(truth, answer)

        leaf_positive_total += length(leaf_positives)
        pruned_total += length(leaf_positives) - length(answer)
    end
    # Internal pruning has to actually happen: a regression that inflated the
    # block caps to trivial full-sphere extents would still satisfy both subset
    # assertions above (it only ever *adds* leaves to the answer) and would
    # silently pass the rest of this testset. Under that regression every
    # internal node passes, so the traversal degenerates to testing every node
    # of every face tree — strictly MORE extents than there are leaves
    # (`50 * npix`). Real pruning comes in far under that: measured 4028 vs
    # 9600 at nside = 4 and 6907 vs 15000 at nside = 5, against 12000 / 24600
    # for the fully unpruned walk.
    @test tested_total < 50 * npix

    # Informational: how much of the leaf-extent-positive set internal pruning
    # removes. Every removed leaf is a false positive of the inflated leaf cap,
    # never a real intersection — which is exactly what the `truth ⊆ answer`
    # assertion above pins down. NOT asserted positive: a leaf whose own cap
    # meets the query cap while an ancestor's does not is a rare accident of
    # cap slack, not an invariant (0 of 50 caps here at nside = 4, 5 of 500;
    # 0 of 500 at nside = 2 and nside = 16), so `pruned_total > 0` would be a
    # flaky assertion. `tested_total` above is the load-bearing one.
    record!("leaf-extent positives pruned by ancestors (fraction)",
        pruned_total / max(leaf_positive_total, 1))
end

# --------------------------------------------------------------------------
# 5. Regridder vs the id-hierarchy path
# --------------------------------------------------------------------------

@testset "face grid vs id-hierarchy Regridder (nside = 4, level = 2)" begin
    face_grid(nside) = treeify(HealpixFaceGrid(nside; ordering=NestedOrder()))
    dggs_grid(level) = treeify(DGGSGrid(HEALPixDGGS(), level))

    # Column identity is the same on both sides — nested data position `p + 1`
    # on the face-grid side, DGGS ordinal `p + 1` on the id-hierarchy side — so
    # the matrices are comparable entry by entry with no permutation.
    #
    # Self-regridding is diagonal (neighbouring pixels share edges, so their
    # intersection area is zero and is dropped), so a cross-resolution pair is
    # run alongside it: that one has genuine off-diagonal structure and is the
    # test that actually exercises polygon clipping on both paths.
    cases = (("self, level 2", face_grid(4), face_grid(4), dggs_grid(2), dggs_grid(2)),
             ("level 2 <- level 1", face_grid(4), face_grid(2), dggs_grid(2), dggs_grid(1)))
    for (label, fdst, fsrc, ddst, dsrc) in cases
        A = CR.Regridder(UNIT, fdst, fsrc; threaded=false, normalize=false)
        B = CR.Regridder(UNIT, ddst, dsrc; threaded=false, normalize=false)

        @test size(A.intersections) == size(B.intersections)
        # Same sparsity pattern: the same cell pairs are found to overlap.
        @test A.intersections.colptr == B.intersections.colptr
        @test A.intersections.rowval == B.intersections.rowval

        # Exact, not approximate: `HealpixKernel.jl`'s `cell_boundary` evaluates
        # the same `chart.jl` closed forms this grid does (via `nested_to_xyf`),
        # so the two paths clip *literally the same* polygons and every entry
        # agrees bitwise. An `isapprox` here would hide a fork between them.
        deviation = maximum(abs, A.intersections - B.intersections)
        @test A.intersections == B.intersections
        record!("face grid vs id hierarchy, max |ΔA| ($label)", deviation)

        # Conservation on both paths, on the unit sphere.
        @test isapprox(sum(A.intersections), 4π; rtol=1e-10)
        @test isapprox(sum(B.intersections), 4π; rtol=1e-10)
        @test isapprox(A.dst_areas, B.dst_areas; rtol=1e-12)
        @test isapprox(A.src_areas, B.src_areas; rtol=1e-12)
    end
end

# --------------------------------------------------------------------------
# 6. Ring vs nested: one grid, two orderings, one permutation
# --------------------------------------------------------------------------

@testset "ring/nested Regridders differ by exactly the index permutation" begin
    nside = 4
    npix = 12nside^2
    nested_tree = treeify(HealpixFaceGrid(nside; ordering=NestedOrder()))
    ring_tree = treeify(HealpixFaceGrid(nside; ordering=RingOrder()))
    N = CR.Regridder(UNIT, nested_tree, nested_tree; threaded=false, normalize=false)
    R = CR.Regridder(UNIT, ring_tree, ring_tree; threaded=false, normalize=false)

    # σ maps a nested data position to the ring data position of the same
    # pixel: position i holds nested id i - 1, whose ring index is
    # `nested_to_ring(i - 1, nside)` (already 1-based).
    sigma = [nested_to_ring(i - 1, nside) for i in 1:npix]
    @test sort(sigma) == collect(1:npix)               # it is a permutation

    # Both orderings hand the clipper *bit-identical* polygons — same chart
    # kernel, same corner order, only the data slot differs — and assembly with
    # `threaded = false` is deterministic. So this is exact equality, not
    # `isapprox`: any drift here would mean the ordering leaked into the
    # geometry, which is precisely what this layer is designed to prevent.
    @test N.intersections == R.intersections[sigma, sigma]
    @test N.dst_areas == R.dst_areas[sigma]
    @test N.src_areas == R.src_areas[sigma]
    # And the geometry itself, cell by cell.
    for i in 1:npix
        @test ring_points(Trees.getcell(nested_tree, i)) ==
            ring_points(Trees.getcell(ring_tree, sigma[i]))
    end
end

# --------------------------------------------------------------------------
# 7. Non-power-of-two conservation
#
# The whole reason this layer is separate from the nested kernel: `nside = 3`
# and `nside = 5` have no nested id space at all, but they are perfectly good
# HEALPix grids and must regrid conservatively.
# --------------------------------------------------------------------------

@testset "non-power-of-two conservation (nside = $nside)" for nside in (3, 5)
    @test !ispow2(nside)
    npix = 12nside^2
    grid = HealpixFaceGrid(nside; ordering=RingOrder())
    R = CR.Regridder(UNIT, treeify(grid), treeify(grid); threaded=false, normalize=false)

    @test size(R.intersections) == (npix, npix)
    @test isapprox(sum(R.intersections), 4π; rtol=1e-10)
    @test isapprox(sum(R.dst_areas), 4π; rtol=1e-10)
    @test isapprox(sum(R.src_areas), 4π; rtol=1e-10)
    record!("|sum(intersections) - 4π|, nside=$nside", abs(sum(R.intersections) - 4π))

    # Every pixel overlaps itself: no cell may be missed by the traversal.
    @test all(i -> R.intersections[i, i] > 0, 1:npix)
    # The *chart* is equal-area, but the cell handed to the clipper is the
    # 4-corner spherical polygon and HEALPix edges are not great circles (they
    # follow constant-z / constant-φ chart lines). So an individual 4-gon is
    # only within ~10% of `4π / npix` — while the *sum* is exact to the last
    # bits above, because neighbours share edges and the bulges cancel pairwise.
    # (Densifying the edges through `xyf_to_point` would shrink the per-cell
    # gap; nothing downstream needs it, since conservation is the property that
    # matters.)
    per_cell = maximum(abs.(R.dst_areas .- 4π / npix) ./ (4π / npix))
    @test per_cell < 0.15
    @test isapprox(sum(R.dst_areas) / npix, 4π / npix; rtol=1e-10)
    record!("max |cell 4-gon area - 4π/npix| / (4π/npix), nside=$nside", per_cell)

    # A constant field regrids to itself — the defining property of a
    # conservative (mean-preserving) operator.
    src = fill(3.5, npix)
    dst = CR.regrid!(zeros(npix), R, src)
    @test maximum(abs, dst .- 3.5) <= 1e-12
    record!("constant-field regrid deviation, nside=$nside", maximum(abs, dst .- 3.5))
end

# --------------------------------------------------------------------------
# 8. CCW discipline
#
# The convex-clip kernel clips a clockwise ring to EMPTY, so a reversed ring
# yields silent zero intersection areas rather than an error. Every polygon
# this layer emits must wind CCW as seen from outside the sphere.
# --------------------------------------------------------------------------

@testset "getcell polygons are CCW from outside (nside = $nside)" for nside in (3, 4, 5)
    ordering = RingOrder()
    root = treeify(HealpixFaceGrid(nside; ordering))
    worst = Inf
    for j in 1:(12nside^2)
        worst = min(worst, ccw_measure(open_ring(Trees.getcell(root, j))))
    end
    @test worst > 0
    # Same for the per-face grid's `(i, j)` accessor, which is the interface
    # method the cursor machinery calls.
    for face in 0:11
        g = FaceChartGrid(GO.Spherical(), HealpixFaceSpace(nside), face, ordering)
        for i in 1:nside, j in 1:nside
            @test ccw_measure(open_ring(Trees.getcell(g, i, j))) > 0
        end
    end
    record!("min CCW measure (nside=$nside)", worst)
end

@printf("[HEALPix face grid] measured:\n")
for key in sort!(collect(keys(MEASURED)))
    @printf("[HEALPix face grid]   %-52s %+.3e\n", key, MEASURED[key])
end

end # module HealpixFaceGridTestSuite
