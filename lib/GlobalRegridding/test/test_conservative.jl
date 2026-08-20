# Conservative weight construction.

import ConservativeRegridding as CR
import DimensionalData as DD
import GeometryOps: SpatialTreeInterface as STI

"""
    DensifiedCellSpace(lon, lat, n)

A one-cell graticule box whose parallel edges have `n` segments, producing a
non-convex spherical ring for clipping tests.
"""
struct DensifiedCellSpace <: RegridSpace
    lon::Tuple{Float64,Float64}
    lat::Tuple{Float64,Float64}
    n::Int
end

ncells(::DensifiedCellSpace) = 1
nchunks(::DensifiedCellSpace) = 1
cellindices(::DensifiedCellSpace, ::Int) = 1:1
manifold(::DensifiedCellSpace) = GOCore.Spherical(; radius = 1.0)
celltree(space::DensifiedCellSpace) = GR.CellCapTree(space, 1:1)

function getcell(space::DensifiedCellSpace, i::Int)
    i == 1 || throw(BoundsError(space, i))
    lon0, lon1 = space.lon
    lat0, lat1 = space.lat
    ring = USPoint[]
    for k in 0:space.n          # south edge, west to east
        push!(ring, toy_point(lon0 + (lon1 - lon0) * k / space.n, lat0))
    end
    for k in 0:space.n          # north edge, east to west
        push!(ring, toy_point(lon1 - (lon1 - lon0) * k / space.n, lat1))
    end
    push!(ring, ring[1])
    return GI.Polygon([GI.LinearRing(ring)])
end

# Helpers

conservative_block(dst, dst_inds, src, src_inds) =
    WeightBlock(
        build_weights!(WeightCOO(length(dst_inds)), Conservative(),
            dst, dst_inds, src, src_inds),
        length(dst_inds), length(src_inds))

cellareas(space, inds) = [GO.area(manifold(space), getcell(space, i)) for i in inds]

@testset "Conservative weights" begin
    @test Conservative() isa AbstractRegriddingMethod

    @testset "identity" begin
        space = ToyLonLatSpace(4, 2)
        inds = cellindices(space, 1)
        block = conservative_block(space, inds, space, inds)
        areas = cellareas(space, inds)
        W = Matrix(block.weights)

        # Self-regridding produces an area-valued diagonal.
        @test size(W) == (8, 8)
        @test LinearAlgebra.diag(W) ≈ areas rtol = 1e-8
        @test vec(sum(W; dims = 2)) ≈ areas rtol = 1e-8

        # Full-coverage denominators equal destination areas.
        @test block.denom ≈ areas rtol = 1e-8
        @test sum(W) ≈ 4pi rtol = 1e-8
    end

    @testset "refinement" begin
        coarse = ToyLonLatSpace(4, 2)
        fine = ToyLonLatSpace(8, 4)
        dst_inds = cellindices(fine, 1)
        src_inds = cellindices(coarse, 1)
        block = conservative_block(fine, dst_inds, coarse, src_inds)
        W = Matrix(block.weights)
        areas = cellareas(fine, dst_inds)

        # Each fine cell lies in one coarse cell.
        @test size(W) == (32, 8)
        @test all(count(>(1e-12), view(W, j, :)) == 1 for j in 1:32)
        @test vec(sum(W; dims = 2)) ≈ areas rtol = 1e-8

        # Refinement preserves total area.
        @test sum(W) ≈ 4pi rtol = 1e-8
        @test block.denom ≈ areas rtol = 1e-8
    end

    @testset "chunked source partitions the block" begin
        fine = ToyLonLatSpace(8, 4)
        whole = ToyLonLatSpace(4, 2)
        chunked = ToyLonLatSpace(4, 2; chunks = (2, 2))
        dst_inds = cellindices(fine, 1)

        reference = conservative_block(fine, dst_inds, whole, cellindices(whole, 1))

        # Non-contiguous chunks still use local block indices.
        @test nchunks(chunked) == 2
        parts = [cellindices(chunked, c) for c in 1:nchunks(chunked)]
        @test !(parts[1] isa AbstractUnitRange)
        @test sort(vcat(parts...)) == collect(1:8)

        assembled = zeros(32, 8)
        denom = zeros(32)
        for inds in parts
            block = conservative_block(fine, dst_inds, chunked, inds)
            @test size(block) == (32, length(inds))
            assembled[:, inds] .+= Matrix(block.weights)
            denom .+= block.denom
        end

        @test assembled ≈ Matrix(reference.weights) rtol = 1e-8
        @test denom ≈ reference.denom rtol = 1e-8
    end

    @testset "cells across the antimeridian" begin
        fine = ToyLonLatSpace(8, 4)
        # Unit-sphere geometry handles cells crossing 180°.
        band = ToyLonLatSpace(2, 2; lon = (135.0, 225.0), lat = (-45.0, 45.0))
        dst_inds = cellindices(fine, 1)
        src_inds = cellindices(band, 1)
        block = conservative_block(fine, dst_inds, band, src_inds)
        W = Matrix(block.weights)
        areas = cellareas(band, src_inds)

        # Band cells coincide with four destination cells.
        @test count(>(1e-12), W) == 4
        for k in 1:4
            j = cellat(fine, cellcentroid(band, k))
            @test W[j, k] ≈ areas[k] rtol = 1e-8
        end
        @test sum(W) ≈ sum(areas) rtol = 1e-8
    end

    @testset "non-convex cells" begin
        dense = DensifiedCellSpace((0.0, 90.0), (0.0, 45.0), 12)
        tiling = ToyLonLatSpace(16, 8)
        exact = GO.area(manifold(dense), getcell(dense, 1))
        tiling_inds = cellindices(tiling, 1)

        # A non-convex source intersects a full tiling to its own area.
        as_source = conservative_block(tiling, tiling_inds, dense, 1:1)
        @test sum(as_source.weights) ≈ exact rtol = 1e-6

        # GeometryOps currently requires the destination clip ring to be convex.
        as_destination = conservative_block(dense, 1:1, tiling, tiling_inds)
        @test_broken isapprox(sum(as_destination.weights), exact; rtol = 1e-6)
    end

    @testset "the block is the assembled matrix, entry for entry" begin
        # Direct sparse-block copying preserves entries exactly.
        coarse = ToyLonLatSpace(4, 2)
        fine = ToyLonLatSpace(8, 4)
        dst_inds, src_inds = cellindices(coarse, 1), cellindices(fine, 1)
        block = conservative_block(coarse, dst_inds, fine, src_inds)
        areas = CR.intersection_areas(manifold(coarse), GOCore.False(),
            GR.subtree(coarse, dst_inds), GR.subtree(fine, src_inds); progress = false)
        @test block.weights == areas
        # Denominators use the same stored entries.
        @test block.denom == vec(sum(areas; dims = 2))
    end

    @testset "banded chunk extents" begin
        # Full-longitude row chunks retain bounded caps.
        banded = ToyLonLatSpace(16, 8; chunks = (16, 2))
        caps = chunktree(banded).caps
        @test all(cap.radius < Float64(pi) for cap in caps)

        # Caps cover geodesic edges, not only corners.
        function samples(space, i, nseg = 6)
            ring = collect(GI.getpoint(GI.getexterior(getcell(space, i))))
            out = USPoint[]
            for k in 1:(length(ring)-1)
                p, q = ring[k], ring[k+1]
                for t in range(0, 1; length = nseg + 1)
                    v = (1 - t) .* Tuple(p) .+ t .* Tuple(q)
                    n = sqrt(sum(v .^ 2))
                    n > 0 && push!(out, USPoint(v[1] / n, v[2] / n, v[3] / n))
                end
            end
            return out
        end
        @test all(US._contains(caps[c], p)
                  for c in 1:nchunks(banded)
                  for i in cellindices(banded, c)
                  for p in samples(banded, i))

        # Opposite poles lie outside the band caps.
        @test !US._contains(caps[1], toy_point(0.0, 90.0))
        @test !US._contains(caps[4], toy_point(0.0, -90.0))

        # Regional chunks retain convex caps.
        region = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (4, 2))
        @test all(c -> c.radius < Float64(pi) / 2, chunktree(region).caps)
    end

    @testset "threaded and serial builds agree bit for bit" begin
        # Threaded and serial intersection assembly must be bit-identical.
        dst = RasterGrid(DD.DimArray(zeros(30, 15), (DD.X(range(-174.0, 174.0; length = 30)),
            DD.Y(range(-84.0, 84.0; length = 15)))))
        src = RasterGrid(DD.DimArray(zeros(48, 24), (DD.X(range(-176.25, 176.25; length = 48)),
            DD.Y(range(-86.25, 86.25; length = 24)))))
        ndst, nsrc = ncells(dst), ncells(src)
        block = conservative_block(dst, 1:ndst, src, 1:nsrc)

        m = manifold(dst)
        # One operator per call: the inner one carries a mutable clipping cache.
        blockop() = GR.BlockAreaOperator(CR.DefaultIntersectionOperator(m),
            GR.indexmap(1:ndst), GR.indexmap(1:nsrc))
        assemble(threaded) = WeightBlock(
            GR._fillcoo!(WeightCOO(ndst),
                CR.intersection_areas(m, threaded,
                    GR.subtree(dst, 1:ndst), GR.subtree(src, 1:nsrc);
                    intersection_operator = blockop(), progress = false)),
            ndst, nsrc)
        reference = assemble(GOCore.False())

        @test GR.SparseArrays.nnz(block.weights) > 2000

        # The explicitly threaded build and the production build — whichever way
        # the thread count sends the latter — both equal the serial reference.
        for other in (assemble(GOCore.True()), block)
            @test other.weights.colptr == reference.weights.colptr
            @test other.weights.rowval == reference.weights.rowval
            @test all(other.weights.nzval .=== reference.weights.nzval)
            @test all(other.denom .=== reference.denom)
        end
    end

    @testset "a pair whose descent finds nothing" begin
        # Deep disjoint trees assemble an empty block.
        north = RasterGrid(DD.DimArray(zeros(8, 8),
            (DD.X(range(-3.5, 3.5; length = 8)), DD.Y(range(66.5, 73.5; length = 8)))))
        south = RasterGrid(DD.DimArray(zeros(8, 8),
            (DD.X(range(-3.5, 3.5; length = 8)), DD.Y(range(-73.5, -66.5; length = 8)))))
        @test !STI.isleaf(celltree(north))
        @test !STI.isleaf(celltree(south))

        block = conservative_block(north, 1:64, south, 1:64)
        @test GR.SparseArrays.nnz(block.weights) == 0
        @test block.denom == zeros(64)
    end

    @testset "one tile's cell geometry, synthesized once" begin
        # Cached tile geometry matches the wrapped space exactly.
        space = RasterGrid(DD.DimArray(zeros(12, 6),
            (DD.X(range(-165.0, 165.0; length = 12)),
                DD.Y(range(-75.0, 75.0; length = 6)))); chunks = ([1:6, 7:12], [1:3, 4:6]))
        inds = cellindices(space, 3)
        tc = GR.TileCells(space, inds)
        tree = GR.subtree(tc, inds)
        @test tree isa GR.CachedCellTree
        @test all(getcell(tree, i) == getcell(space, i) for i in inds)

        # Indices outside the tile fall through to the wrapped tree.
        outside = setdiff(1:ncells(space), inds)
        @test all(getcell(tree, i) == getcell(space, i) for i in outside)

        # Only the exact tile subtree uses cached geometry.
        @test all(getcell(tc, i) == getcell(space, i) for i in 1:ncells(space))
        @test !(GR.subtree(tc, cellindices(space, 1)) isa GR.CachedCellTree)

        # Cached geometry preserves weights bit for bit.
        other = cellindices(space, 2)
        plain = conservative_block(space, inds, space, other)
        stored = conservative_block(tc, inds, space, other)
        @test plain.weights.colptr == stored.weights.colptr
        @test plain.weights.rowval == stored.weights.rowval
        @test all(plain.weights.nzval .=== stored.weights.nzval)
        @test all(plain.denom .=== stored.denom)
    end

    @testset "cap tree split weights are the node windows" begin
        space = ToyLonLatSpace(8, 4)
        tree = GR.CellCapTree(space, 1:ncells(space))

        # Split weights are the node's own window: children partition the parent.
        @test CR.Trees.split_weight(tree) == ncells(space)
        @test sum(CR.Trees.split_weight, STI.getchild(tree)) == CR.Trees.split_weight(tree)
        @test all(c -> CR.Trees.split_weight(c) < CR.Trees.split_weight(tree), STI.getchild(tree))
    end

    @testset "one tile's restricted tree, built once" begin
        # A non-chunk range on a fallback space forces the CellCapTree path.
        xd = DD.X(-168.75:22.5:168.75)
        yd = DD.Y(-78.75:22.5:78.75)
        src = RasterGrid(DD.DimArray(zeros(16, 8), (xd, yd));
            chunks = ([1:8, 9:16], [1:4, 5:8]))
        dst = CountingSpace(ToyLonLatSpace(8, 4; chunks = (8, 2)))
        inds = 5:20
        s1, s2 = cellindices(src, 1), cellindices(src, 4)

        # One shared wrapper builds the destination tree once for both blocks.
        tc = GR.TileCells(dst, inds)
        b0 = dst.builds[]
        sh1 = conservative_block(tc, inds, src, s1)
        sh2 = conservative_block(tc, inds, src, s2)
        @test dst.builds[] - b0 == 1

        # Fresh wrappers rebuild per block: the uncached reference.
        b1 = dst.builds[]
        fr1 = conservative_block(GR.TileCells(dst, inds), inds, src, s1)
        fr2 = conservative_block(GR.TileCells(dst, inds), inds, src, s2)
        @test dst.builds[] - b1 == 2

        # The memoized tree changes nothing, bit for bit.
        for (a, b) in ((sh1, fr1), (sh2, fr2))
            @test a.weights.colptr == b.weights.colptr
            @test a.weights.rowval == b.weights.rowval
            @test all(a.weights.nzval .=== b.weights.nzval)
            @test all(a.denom .=== b.denom)
        end
    end

    @testset "cell memos change nothing, and hits return the built value" begin
        space = ToyLonLatSpace(16, 8)
        memo = GR._cellmemo(space, 1:ncells(space))
        @test memo isa GR.CellMemo
        tree = GR.subtree(space, 1:ncells(space))
        a = GR._memocell(memo, tree, 5)
        @test a == CR.Trees.getcell(tree, 5)
        # A hit returns the stored object; a colliding key rebuilds.
        @test GR._memocell(memo, tree, 5) === a
        b = GR._memocell(memo, tree, 5 + length(memo.keys))
        @test b == CR.Trees.getcell(tree, 5 + length(memo.keys))
        @test GR._memocell(memo, tree, 5) == a
        # The memo-free operator spelling still assembles (covered above by the
        # threaded-vs-serial testset, whose reference uses it).
        @test GR.BlockAreaOperator(1, 2, 3).srcmemo === nothing
    end

    @testset "the area-only clip matches the default operator, bit for bit" begin
        m = GOCore.Spherical(; radius = 1.0)
        default = CR.task_local_operator(CR.DefaultIntersectionOperator(m))
        streaming = CR.task_local_operator(GR.SphericalClipAreaOperator(m))
        @test GR._intersectionoperator(m) isa GR.SphericalClipAreaOperator
        @test GR._intersectionoperator(GOCore.Planar()) isa CR.DefaultIntersectionOperator

        coarse = ToyLonLatSpace(4, 2)
        fine = ToyLonLatSpace(8, 4)
        dense = DensifiedCellSpace((10.0, 30.0), (20.0, 40.0), 5)
        cells = vcat([getcell(coarse, i) for i in 1:8],
            [getcell(fine, i) for i in 1:32], [getcell(dense, 1)])
        # Overlapping, contained, identical, and disjoint pairs, one value each:
        # every branch of the streaming kernel against the materializing one.
        checked = 0
        for a in cells, b in cells
            va = default(a, b)
            vb = streaming(a, b)
            @test va === vb
            checked += 1
        end
        @test checked == length(cells)^2
    end

    @testset "the eager whole block adopts the assembled matrix unchanged" begin
        dst = ToyLonLatSpace(4, 2)
        src = ToyLonLatSpace(8, 4)
        fast = GR.wholeblock(Conservative(), dst, src)
        slow = invoke(GR.wholeblock,
            Tuple{AbstractRegriddingMethod,RegridSpace,RegridSpace},
            Conservative(), dst, src)
        @test fast.weights.colptr == slow.weights.colptr
        @test fast.weights.rowval == slow.weights.rowval
        @test all(fast.weights.nzval .=== slow.weights.nzval)
        @test all(fast.denom .=== slow.denom)
        @test GR.hasdenom(fast) == GR.hasdenom(slow) == true
    end

    @testset "disjoint chunks keep zero denominators" begin
        north = ToyLonLatSpace(2, 1; lat = (60.0, 90.0))
        south = ToyLonLatSpace(2, 1; lat = (-90.0, -60.0))
        block = conservative_block(north, 1:2, south, 1:2)

        # Disjoint conservative blocks retain zero denominators.
        @test GR.SparseArrays.nnz(block.weights) == 0
        @test block.denom == zeros(2)
    end
end
