module S2FaceGridTestSuite

# Tests for `src/S2/face_grid.jl`: dense S2 face grids (all `6 nside²` cells of
# one resolution) as spatial trees, under a swappable data ordering.
#
# Three things are being pinned down here, in rising order of load-bearingness:
#
# 1. *Orderings are bijections.* `data_index` / `lattice_index` must be exact
#    mutual inverses over `1:6nside²`, because that pair is the whole alignment
#    contract: column `j` of a `Regridder` is data position `j`.
#
# 2. *Node extents contain their subtree's geometry.* This is what makes tree
#    pruning sound. It is asserted exhaustively — every node of every face tree,
#    against both the block's lattice vertices and a dense sampling of the
#    block's chart rectangle. For S2 the containment is provable rather than
#    measured (see the argument on `STI.node_extent` in `face_grid.jl`: block
#    edges are geodesics, and no boundary point is within a quarter turn of the
#    cap centre, so the farthest boundary point is a corner); these assertions
#    are the check on the proof, not a substitute for it.
#
#    Note what is deliberately NOT claimed: a node's extent does not contain its
#    descendants' *caps*. It cannot, for any extent that tightly bounds only the
#    node's own geometry — a leaf's circumscribed cap sticks out past the block
#    boundary by roughly its own radius, whatever the parent does. So the query
#    comparison below is a two-sided sandwich rather than an equality:
#    `query ⊇ {leaves that really meet the cap}` (soundness, the load-bearing
#    half) and `query ⊆ {leaves whose leaf extent meets the cap}` (no spurious
#    extras).
#
# 3. *The grid regrids conservatively, and the Hilbert ordering nests.* There is
#    no id-hierarchy *tree* to compare against — `DGGSGrid(S2DGGS(), level)`
#    needs `has_ordinal_ids`, which stays false while the native `s2_cellid`
#    encoding is unported — so its place is taken by cross-resolution and
#    cross-system `Regridder` checks: a Hilbert `nside = 4` grid refines a
#    Hilbert `nside = 2` grid in contiguous blocks of four, and an S2 grid
#    regrids against a HEALPix grid conservatively. The DGGS *geometry* path
#    (`cell_polygon(S2DGGS(), level, id)`, wired in `src/S2/S2Kernel.jl`) does
#    exist now; a test below smoke-checks that it agrees with this grid, and
#    `test/S2/test_s2_kernel.jl` sweeps it bitwise.
#
# Naming: `S2.cell_polygon` / `S2.cell_center` are functions in the `S2`
# namespace and are NOT methods of the top-level `cell_polygon(::AbstractDGGS,
# level, id)` that `using DiscreteGlobalGrids` brings into scope here. Both are
# used below, so the S2 ones are always written qualified (`S2.cell_polygon`)
# and the DGGS one bare.
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
using DiscreteGlobalGrids.S2
using DiscreteGlobalGrids.HEALPix                 # for the cross-system case
using DiscreteGlobalGrids.S2: FaceChartGrid, S2FaceRoot,
    num_cells, data_index, lattice_index, validate_ordering,
    stf_to_point, cell_corners, xyf_to_rowmajor, xyf_to_hilbert, hilbert_to_xyf

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
# Same measure as `test_chart.jl` uses on `cell_corners`.
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

"Is any point of cell `(ix, iy, face)` inside `cap`? (dense chart sampling)"
function cell_meets_cap(cap, ix, iy, face, nside, samples=4)
    for a in 0:samples, b in 0:samples
        p = stf_to_point((ix + a / samples) / nside, (iy + b / samples) / nside, face)
        US.spherical_distance(cap.point, p) <= cap.radius && return true
    end
    return false
end

# --------------------------------------------------------------------------
# 1. Construction and validation
# --------------------------------------------------------------------------

@testset "S2FaceSpace / S2FaceGrid construction" begin
    @test_throws ArgumentError S2FaceSpace(0)
    @test_throws ArgumentError S2FaceSpace(-3)
    @test S2FaceSpace(1).nside == 1
    @test num_cells(S2FaceSpace(3)) == 54
    @test num_cells(S2FaceSpace(4)) == 96

    # Upper bound: `6 * nside^2` is `3 * 2^61 < typemax(Int64)` at `nside = 2^30`
    # and overflows silently at `2^31`. `2^30` is also `nside` at the deepest
    # native S2 level (`max_level(S2DGGS()) == 30`).
    @test S2FaceSpace(2^30).nside == 2^30
    @test_throws ArgumentError S2FaceSpace(2^30 + 1)
    @test_throws ArgumentError S2FaceSpace(2^31)
    @test 6 * Int64(2^30)^2 < typemax(Int64)

    # The Hilbert position is built two bits per level, so it exists only on a
    # 2^k x 2^k face. It must be refused at *construction*, not from inside a
    # traversal.
    @test_throws ArgumentError S2FaceGrid(3; ordering=HilbertOrder())
    @test_throws ArgumentError S2FaceGrid(5; ordering=HilbertOrder())
    @test_throws ArgumentError S2FaceGrid(S2FaceSpace(6), HilbertOrder())
    @test_throws ArgumentError S2FaceGrid(0)
    @test_throws ArgumentError S2FaceGrid(0; ordering=HilbertOrder())

    # Row-major carries no such restriction — that is the point of this layer.
    for nside in (1, 2, 3, 4, 5, 7)
        grid = S2FaceGrid(nside)
        @test grid.ordering isa RowMajorOrder          # the default
        @test grid.space.nside == nside
        @test num_cells(grid) == 6nside^2
    end
    @test S2FaceGrid(4; ordering=HilbertOrder()).ordering isa HilbertOrder
    @test RowMajorOrder() isa AbstractS2Ordering
    @test HilbertOrder() isa AbstractS2Ordering
    # The optional third contract method defaults to "accept everything". It
    # takes the *space*, not a loose `nside`: the resolution reaches this layer
    # only through the type that has already checked it.
    @test validate_ordering(RowMajorOrder(), S2FaceSpace(3)) === nothing
    @test validate_ordering(HilbertOrder(), S2FaceSpace(4)) === nothing
    @test_throws ArgumentError validate_ordering(HilbertOrder(), S2FaceSpace(3))

    @test sprint(show, S2FaceSpace(4)) == "S2FaceSpace(4)"
    @test occursin("nside=4", sprint(show, S2FaceGrid(4)))
    @test occursin("96 cells", sprint(show, S2FaceGrid(4)))
end

# `S2FaceRoot` is independently constructible — it is the REPL/test convenience
# form — so it must not be a way *around* the checks `S2FaceGrid` runs. Its
# resolution is an `S2FaceSpace` (so `nside` is checked by construction) and it
# re-runs both ordering checks itself.
@testset "S2FaceRoot construction validates like the grid" begin
    # (a) `nside`, through the space the convenience form builds.
    @test_throws ArgumentError S2FaceRoot(0)
    @test_throws ArgumentError S2FaceRoot(-5, RowMajorOrder())
    @test_throws ArgumentError S2FaceRoot(2^31, RowMajorOrder())
    @test S2FaceRoot(2^30, RowMajorOrder()).space.nside == 2^30

    # (b) another system's ordering: an `ArgumentError` here rather than a
    # `MethodError` from inside a traversal.
    @test_throws ArgumentError S2FaceRoot(4, DiscreteGlobalGrids.HEALPix.RingOrder())
    @test_throws ArgumentError S2FaceRoot(
        GO.Spherical(), S2FaceSpace(4), DiscreteGlobalGrids.ISEA4R.MortonOrder())

    # (c) `validate_ordering`, which now runs at root construction too.
    @test_throws ArgumentError S2FaceRoot(3, HilbertOrder())
    @test_throws ArgumentError S2FaceRoot(GO.Spherical(), S2FaceSpace(5), HilbertOrder())
    @test S2FaceRoot(4, HilbertOrder()).ordering isa HilbertOrder
end

@testset "treeify and tree toplevel" begin
    for (nside, ordering) in ((4, HilbertOrder()), (5, RowMajorOrder()), (3, RowMajorOrder()))
        grid = S2FaceGrid(nside; ordering)
        # One-argument `treeify` resolves through `best_manifold`.
        @test GOCore.best_manifold(grid) == GO.Spherical()
        root = treeify(grid)
        @test root isa S2FaceRoot
        @test root.space == S2FaceSpace(nside)
        @test root.space.nside == nside
        @test root.ordering === ordering
        @test GOCore.best_manifold(root) == GO.Spherical()
        # Explicit-manifold form, and idempotent passthrough on the tree.
        @test treeify(GO.Spherical(), grid) isa S2FaceRoot
        @test treeify(root) === root
        @test treeify(GO.Spherical(), root) === root
        @test Trees.treeify(UNIT, root) === root
        # The direct constructor defaults the manifold for REPL use.
        @test S2FaceRoot(nside, ordering) == root

        @test STI.isspatialtree(typeof(root))
        @test !STI.isleaf(root)
        @test STI.nchild(root) == 6
        @test Trees.ncells(root) == 6nside^2
        @test length(collect(STI.getchild(root))) == 6
        @test occursin("nside=$nside", sprint(show, root))

        # The root bounds nothing tighter than the whole sphere: the six faces
        # tile it.
        @test STI.node_extent(root).radius >= Float64(pi)

        # Children are stock quadtree cursors over per-face chart grids.
        for f in 1:6
            child = STI.getchild(root, f)
            @test child isa Trees.TopDownQuadtreeCursor{<:FaceChartGrid}
            @test child.grid.face == f - 1
            @test child.grid.space == S2FaceSpace(nside)
            @test Trees.ncells(child.grid, 1) == nside
            @test Trees.ncells(child.grid, 2) == nside
            @test GOCore.manifold(child.grid) == GO.Spherical()
        end

        @test_throws BoundsError Trees.getcell(root, 0)
        @test_throws BoundsError Trees.getcell(root, 6nside^2 + 1)
        @test length(collect(Trees.getcell(root))) == 6nside^2
    end
end

# --------------------------------------------------------------------------
# 2. The ordering contract
# --------------------------------------------------------------------------

@testset "ordering bijections (nside = $nside)" for nside in (1, 2, 3, 4, 5, 7, 8)
    ncell = 6nside^2
    space = S2FaceSpace(nside)
    orderings = ispow2(nside) ? (RowMajorOrder(), HilbertOrder()) : (RowMajorOrder(),)
    for ordering in orderings
        # `lattice_index` is onto the whole lattice and `data_index` inverts it.
        seen = Set{NTuple{3,Int}}()
        for j in 1:ncell
            ix, iy, face = lattice_index(ordering, space, j)
            @test 0 <= ix < nside && 0 <= iy < nside && 0 <= face <= 5
            push!(seen, (ix, iy, face))
            @test data_index(ordering, space, ix, iy, face) == j
        end
        @test length(seen) == ncell

        # ... and the other direction, over the lattice.
        for face in 0:5, ix in 0:(nside - 1), iy in 0:(nside - 1)
            j = data_index(ordering, space, ix, iy, face)
            @test 1 <= j <= ncell
            @test lattice_index(ordering, space, j) == (ix, iy, face)
        end
    end

    # The two shipped orderings are exactly the chart's two index maps, with the
    # single 0-based/1-based reconciliation both need (unlike HEALPix, where
    # RING is already 1-based and only NESTED takes the `+ 1`).
    ix, iy = min(1, nside - 1), min(2, nside - 1)
    @test data_index(RowMajorOrder(), space, ix, iy, 3) ==
        xyf_to_rowmajor(ix, iy, 3, nside) + 1
    if ispow2(nside)
        @test data_index(HilbertOrder(), space, ix, iy, 3) ==
            xyf_to_hilbert(ix, iy, 3, nside) + 1
        # `min` because `nside = 1` has only 6 cells in total (HEALPix's 12
        # would have absorbed a bare `7` here).
        j = min(7, ncell)
        @test lattice_index(HilbertOrder(), space, j) == hilbert_to_xyf(j - 1, nside)
    end
end

@testset "getcell round-trips through the ordering (nside = $nside, $ordering)" for
        (nside, ordering) in ((4, HilbertOrder()), (4, RowMajorOrder()),
                              (5, RowMajorOrder()), (3, RowMajorOrder()))

    space = S2FaceSpace(nside)
    root = treeify(S2FaceGrid(nside; ordering))
    ncell = 6nside^2
    for j in 1:ncell
        ix, iy, face = lattice_index(ordering, space, j)
        # `Trees.getcell(root, j)` is the geometry side of the alignment rule.
        @test ring_points(Trees.getcell(root, j)) ==
            ring_points(S2.cell_polygon(ix, iy, face, nside))
        # Closed CCW 4-gon, with the corners `chart.jl` emits, in its order.
        pts = ring_points(Trees.getcell(root, j))
        @test length(pts) == 5 && pts[1] == pts[5]
        @test Tuple(pts[1:4]) == cell_corners(ix, iy, face, nside)
    end
    # The iterator form must agree with the indexed one, position by position.
    @test collect(Trees.getcell(root)) == [Trees.getcell(root, j) for j in 1:ncell]

    # The per-face grid's own index maps target the *global* data layout.
    for face in 0:5
        g = FaceChartGrid(GO.Spherical(), space, face, ordering)
        for ix in 0:(nside - 1), iy in 0:(nside - 1)
            j = Trees.cartesian_to_linear_idx(g, CartesianIndex(ix + 1, iy + 1))
            @test j == data_index(ordering, space, ix, iy, face)
            @test Trees.linear_to_cartesian_idx(g, j) == CartesianIndex(ix + 1, iy + 1)
            @test ring_points(Trees.getcell(g, ix + 1, iy + 1)) ==
                ring_points(Trees.getcell(root, j))
            @test Trees.getvertex(g, ix + 1, iy + 1) ==
                stf_to_point(ix / nside, iy / nside, face)
        end
    end
end

# `HilbertOrder` is written against the scaffold ordinal the `S2DGGS` registry
# entry records: `face * 4^level + hilbert_position` at `nside = 2^level`. What
# is pinned here is the ordinal arithmetic alone; the geometry comparison
# against the DGGS kernel over the same ordinals is `test/S2/test_s2_kernel.jl`.
@testset "HilbertOrder position j is scaffold ordinal j - 1" begin
    for level in 0:3
        nside = 2^level
        for face in 0:5, ix in 0:(nside - 1), iy in 0:(nside - 1)
            j = data_index(HilbertOrder(), S2FaceSpace(nside), ix, iy, face)
            ordinal = j - 1
            position = ordinal - face * 4^level
            @test 0 <= position < 4^level
            @test ordinal == face * 4^level + position
            @test hilbert_to_xyf(ordinal, nside) == (ix, iy, face)
        end
        @test root_count(S2DGGS()) == 6                # the 6 in `6 * 4^level`
        @test radix(S2DGGS()) == 4
    end
end

# --------------------------------------------------------------------------
# 3. Node extents contain their subtree's geometry
#
# The soundness check for pruning, by exhaustion. Two probes per node: the
# block's lattice vertices (the polygon corners of every cell it owns) and a
# dense sampling of the block's chart rectangle.
# --------------------------------------------------------------------------

@testset "node extents contain their block geometry (nside = $nside)" for nside in (3, 4, 5, 8)
    root = treeify(S2FaceGrid(nside; ordering=RowMajorOrder()))
    worst_vertex = -Inf
    worst_sample = -Inf
    nodes = 0
    for f in 1:6
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
            s0, s1 = (imin - 1) / nside, (imax - 1) / nside
            t0, t1 = (jmin - 1) / nside, (jmax - 1) / nside
            for a in 0:8, b in 0:8
                p = stf_to_point(s0 + (s1 - s0) * a / 8, t0 + (t1 - t0) * b / 8, grid.face)
                worst_sample = max(worst_sample,
                    US.spherical_distance(extent.point, p) - extent.radius)
            end
            return nothing
        end
    end
    @test nodes > 6                                    # the walk actually descended
    # Non-positive overhang: every probe sits inside the cap.
    @test worst_vertex <= 0
    @test worst_sample <= 0
    record!("node-extent overhang, lattice vertices (rad)", worst_vertex)
    record!("node-extent overhang, dense chart samples (rad)", worst_sample)

    # Leaf extents likewise contain their own cell, and the leaf index space is
    # exactly `1:6nside²` with no index emitted twice.
    exts = leaf_extents(root)
    @test sort!(collect(keys(exts))) == collect(1:(6nside^2))
    worst_leaf = -Inf
    for (j, extent) in exts
        for p in open_ring(Trees.getcell(root, j))
            worst_leaf = max(worst_leaf,
                US.spherical_distance(extent.point, p) - extent.radius)
        end
    end
    @test worst_leaf <= 0
    record!("leaf-extent overhang, cell corners (rad)", worst_leaf)
end

# --------------------------------------------------------------------------
# 4. Cap queries against brute force
# --------------------------------------------------------------------------

@testset "STI queries vs brute force (nside = $nside, $ordering)" for
        (nside, ordering) in ((4, HilbertOrder()), (4, RowMajorOrder()), (5, RowMajorOrder()))

    root = treeify(S2FaceGrid(nside; ordering))
    ncell = 6nside^2
    exts = leaf_extents(root)
    lattice = [lattice_index(ordering, S2FaceSpace(nside), j) for j in 1:ncell]

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
        @test all(j -> 1 <= j <= ncell, answer)

        # Upper bound: the traversal only ever *returns* leaves whose own leaf
        # extent passes the predicate, so no result can be outside this set.
        leaf_positives = [j for j in 1:ncell if intersects_cap(cap, exts[j])]
        @test issubset(answer, leaf_positives)

        # Lower bound — the load-bearing half. Any cell that genuinely meets the
        # cap must survive every ancestor's pruning test. (Ground truth by dense
        # chart sampling of the cell: a sampled point inside the cap is a
        # witness that the cell really intersects it.)
        truth = [j for j in 1:ncell if cell_meets_cap(cap, lattice[j]..., nside)]
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
    # (`50 * ncell`). Real pruning comes in comfortably under that.
    @test tested_total < 50 * ncell

    # Informational: how much of the leaf-extent-positive set internal pruning
    # removes. Every removed leaf is a false positive of the inflated leaf cap,
    # never a real intersection — which is exactly what the `truth ⊆ answer`
    # assertion above pins down. NOT asserted positive: a leaf whose own cap
    # meets the query cap while an ancestor's does not is a rare accident of cap
    # slack, not an invariant, so `pruned_total > 0` would be flaky.
    # `tested_total` above is the load-bearing one.
    record!("leaf-extent positives pruned by ancestors (fraction)",
        pruned_total / max(leaf_positive_total, 1))
end

# --------------------------------------------------------------------------
# 5. Regridder structure
#
# The place `test/HEALPix/test_face_grid.jl` compares the face-grid path against
# the id-hierarchy path. S2 has no id-hierarchy *grid* (that needs
# `has_ordinal_ids`), so (a) checks what it does have — the DGGS geometry path
# over the same ordinals — and records what is still missing, while (b)-(d)
# stand in for the matrix comparison.
# --------------------------------------------------------------------------

@testset "the DGGS geometry path agrees; the id hierarchy is not ported" begin
    # `cell_polygon` here is the top-level `cell_polygon(::AbstractDGGS, level,
    # id)` — NOT `S2.cell_polygon(ix, iy, face, nside)`, which is what the rest
    # of this file exercises. It used to throw `NotPortedError`; since
    # `src/S2/S2Kernel.jl` it answers over the scaffold ordinal, and the ordinal
    # it answers for is exactly this grid's `HilbertOrder` data position minus
    # one. A smoke check of that correspondence lives here; the full bitwise
    # sweep, the caps and the still-unported hierarchy group are
    # `test/S2/test_s2_kernel.jl`'s business.
    root = treeify(S2FaceGrid(8; ordering=HilbertOrder()))
    for id in (0, 17, 6 * 64 - 1)
        @test ring_points(cell_polygon(S2DGGS(), 3, id)) ==
              ring_points(Trees.getcell(root, id + 1))
    end
    # What is still NOT ported: the native `s2_cellid` id hierarchy, hence the
    # whole hierarchy/ordinal/pruning group.
    @test_throws NotPortedError cell_children(S2DGGS(), 0, 0)
    @test_throws NotPortedError descendant_range(S2DGGS(), 0, 0, 1)
    @test max_level(S2DGGS()) == 30                    # nside = 2^30 at the leaf level
end

@testset "self-regridding is diagonal (nside = 4)" begin
    nside = 4
    ncell = 6nside^2
    tree = treeify(S2FaceGrid(nside; ordering=RowMajorOrder()))
    R = CR.Regridder(UNIT, tree, tree; threaded=false, normalize=false)

    @test size(R.intersections) == (ncell, ncell)
    # Neighbouring cells share *geodesic* edges, so their intersection is a
    # zero-area sliver and is dropped: the matrix is exactly diagonal. (The same
    # holds for HEALPix, whose shared edges are shared chart lines.)
    @test SparseArrays.nnz(R.intersections) == ncell
    @test SparseArrays.diag(R.intersections) ≈ R.dst_areas
    @test isapprox(sum(R.intersections), 4π; rtol=1e-10)
end

@testset "cross-resolution Hilbert nesting (nside 4 <- 2)" begin
    fine = treeify(S2FaceGrid(4; ordering=HilbertOrder()))
    coarse = treeify(S2FaceGrid(2; ordering=HilbertOrder()))
    R = CR.Regridder(UNIT, fine, coarse; threaded=false, normalize=false)
    M = R.intersections
    @test size(M) == (6 * 16, 6 * 4)

    # Every fine cell lies in exactly one coarse cell, and every coarse cell is
    # exactly four fine ones: the Hilbert prefix property (`pos >> 2` is the
    # parent's position), made observable through the clipper.
    @test SparseArrays.nnz(M) == 6 * 16
    for j in 1:size(M, 2)
        rows = sort(M.rowval[M.colptr[j]:(M.colptr[j + 1] - 1)])
        @test rows == collect((4j - 3):(4j))
    end
    for i in 1:size(M, 1)
        @test count(==(i), M.rowval) == 1
    end

    # Areas: each nonzero is the whole fine cell, each column sums to the whole
    # coarse cell, and the total is the sphere.
    @test isapprox(vec(sum(M, dims=2)), R.dst_areas; rtol=1e-12)
    @test isapprox(vec(sum(M, dims=1)), R.src_areas; rtol=1e-12)
    @test isapprox(sum(M), 4π; rtol=1e-10)
    record!("|sum(intersections) - 4π|, Hilbert 4<-2", abs(sum(M) - 4π))
end

@testset "cross-system: S2 against HEALPix (nside = 4)" begin
    s2_tree = treeify(S2FaceGrid(4; ordering=RowMajorOrder()))
    hp_tree = treeify(HealpixFaceGrid(4; ordering=RingOrder()))
    R = CR.Regridder(UNIT, s2_tree, hp_tree; threaded=false, normalize=false)

    @test size(R.intersections) == (6 * 16, 12 * 16)
    # Two completely different tessellations of the same sphere: the shared
    # total is the real check that both sets of polygons are wound the same way
    # and cover everything exactly once.
    @test isapprox(sum(R.intersections), 4π; rtol=1e-10)
    @test isapprox(sum(R.src_areas), 4π; rtol=1e-10)
    @test isapprox(vec(sum(R.intersections, dims=2)), R.dst_areas; rtol=1e-12)
    record!("|sum(intersections) - 4π|, S2 <- HEALPix", abs(sum(R.intersections) - 4π))

    # A constant field regrids to itself — the defining property of a
    # conservative (mean-preserving) operator, here across systems.
    dst = CR.regrid!(zeros(6 * 16), R, fill(3.5, 12 * 16))
    @test maximum(abs, dst .- 3.5) <= 1e-12
    record!("constant-field regrid deviation, S2 <- HEALPix", maximum(abs, dst .- 3.5))
end

# --------------------------------------------------------------------------
# 6. Row-major vs Hilbert: one grid, two orderings, one permutation
# --------------------------------------------------------------------------

@testset "row-major/Hilbert Regridders differ by exactly the index permutation" begin
    nside = 4
    ncell = 6nside^2
    hilbert_tree = treeify(S2FaceGrid(nside; ordering=HilbertOrder()))
    rowmajor_tree = treeify(S2FaceGrid(nside; ordering=RowMajorOrder()))
    H = CR.Regridder(UNIT, hilbert_tree, hilbert_tree; threaded=false, normalize=false)
    R = CR.Regridder(UNIT, rowmajor_tree, rowmajor_tree; threaded=false, normalize=false)

    # σ maps a Hilbert data position to the row-major data position of the same
    # cell: position i holds Hilbert id i - 1, whose lattice cell has row-major
    # id `xyf_to_rowmajor(...)`, hence data position `+ 1`.
    sigma = [Int(xyf_to_rowmajor(hilbert_to_xyf(i - 1, nside)..., nside)) + 1 for i in 1:ncell]
    @test sort(sigma) == collect(1:ncell)              # it is a permutation

    # Both orderings hand the clipper *bit-identical* polygons — same chart
    # kernel, same corner order, only the data slot differs — and assembly with
    # `threaded = false` is deterministic. So this is exact equality, not
    # `isapprox`: any drift here would mean the ordering leaked into the
    # geometry, which is precisely what this layer is designed to prevent.
    @test H.intersections == R.intersections[sigma, sigma]
    @test H.dst_areas == R.dst_areas[sigma]
    @test H.src_areas == R.src_areas[sigma]
    # And the geometry itself, cell by cell.
    for i in 1:ncell
        @test ring_points(Trees.getcell(hilbert_tree, i)) ==
            ring_points(Trees.getcell(rowmajor_tree, sigma[i]))
    end
end

# --------------------------------------------------------------------------
# 7. Non-power-of-two conservation
#
# The whole reason this layer is separate from the (future) id hierarchy:
# `nside = 3` and `nside = 5` have no S2 id space at all, but they are perfectly
# good cube-face grids and must regrid conservatively.
# --------------------------------------------------------------------------

@testset "each face carries exactly 4π/6" begin
    # True for ANY monotone `st_to_uv` — it follows from cube symmetry (the six
    # face frames are in SO(3) and the six gnomonic patches are congruent), NOT
    # from the quadratic transform. The quadratic transform only narrows the
    # spread *within* a face; S2 remains `is_equal_area == false`.
    @test !is_equal_area(S2DGGS())
    grid = S2FaceGrid(1)
    R = CR.Regridder(UNIT, treeify(grid), treeify(grid); threaded=false, normalize=false)
    @test length(R.dst_areas) == 6
    for a in R.dst_areas
        @test isapprox(a, 4π / 6; rtol=1e-12)
    end
    record!("max |face area - 4π/6| / (4π/6)",
        maximum(abs.(R.dst_areas .- 4π / 6) ./ (4π / 6)))
end

@testset "non-power-of-two conservation (nside = $nside)" for nside in (3, 5)
    @test !ispow2(nside)
    ncell = 6nside^2
    grid = S2FaceGrid(nside; ordering=RowMajorOrder())
    R = CR.Regridder(UNIT, treeify(grid), treeify(grid); threaded=false, normalize=false)

    @test size(R.intersections) == (ncell, ncell)
    @test isapprox(sum(R.intersections), 4π; rtol=1e-10)
    @test isapprox(sum(R.dst_areas), 4π; rtol=1e-10)
    @test isapprox(sum(R.src_areas), 4π; rtol=1e-10)
    record!("|sum(intersections) - 4π|, nside=$nside", abs(sum(R.intersections) - 4π))

    # Every cell overlaps itself: no cell may be missed by the traversal.
    @test all(i -> R.intersections[i, i] > 0, 1:ncell)

    # Unlike HEALPix, the 4-gon handed to the clipper IS the cell (geodesic
    # edges), so there is no per-cell 4-gon deficit to bound. What varies is the
    # cell area itself: the quadratic transform leaves a within-level spread
    # that grows toward ~2.08x asymptotically. Measured here: ~1.10 at nside=3,
    # ~1.42 at nside=5.
    ratio = maximum(R.dst_areas) / minimum(R.dst_areas)
    @test ratio < 2.2
    record!("max/min cell area, nside=$nside", ratio)

    # A constant field regrids to itself — the defining property of a
    # conservative (mean-preserving) operator.
    src = fill(3.5, ncell)
    dst = CR.regrid!(zeros(ncell), R, src)
    @test maximum(abs, dst .- 3.5) <= 1e-12
    record!("constant-field regrid deviation, nside=$nside", maximum(abs, dst .- 3.5))
end

# --------------------------------------------------------------------------
# 8. CCW discipline
#
# The convex-clip kernel clips a clockwise ring to EMPTY, so a reversed ring
# yields silent zero intersection areas rather than an error. Every polygon this
# layer emits must wind CCW as seen from outside the sphere.
# --------------------------------------------------------------------------

@testset "getcell polygons are CCW from outside (nside = $nside)" for nside in (3, 4, 5)
    ordering = RowMajorOrder()
    root = treeify(S2FaceGrid(nside; ordering))
    worst = Inf
    for j in 1:(6nside^2)
        worst = min(worst, ccw_measure(open_ring(Trees.getcell(root, j))))
    end
    @test worst > 0
    # Same for the per-face grid's `(i, j)` accessor, which is the interface
    # method the cursor machinery calls.
    for face in 0:5
        g = FaceChartGrid(GO.Spherical(), S2FaceSpace(nside), face, ordering)
        for i in 1:nside, j in 1:nside
            @test ccw_measure(open_ring(Trees.getcell(g, i, j))) > 0
        end
    end
    record!("min CCW measure (nside=$nside)", worst)
end

@printf("[S2 face grid] measured:\n")
for key in sort!(collect(keys(MEASURED)))
    @printf("[S2 face grid]   %-52s %+.3e\n", key, MEASURED[key])
end

end # module S2FaceGridTestSuite
