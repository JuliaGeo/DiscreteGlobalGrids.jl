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
ownedindices(::DensifiedCellSpace, ::Int) = 1:1
manifold(::DensifiedCellSpace) = GOCore.Spherical(; radius = 1.0)
celltree(space::DensifiedCellSpace) = GR.CellSpaceRTree(space, 1:1)

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

"A geometry-free space used to pin the destination tile cache threshold."
struct TileCacheBoundSpace <: RegridSpace
    n::Int
    reads::Base.RefValue{Int}
end

ncells(space::TileCacheBoundSpace) = space.n
nchunks(::TileCacheBoundSpace) = 1
ownedindices(space::TileCacheBoundSpace, ::Int) = 1:space.n
function getcell(space::TileCacheBoundSpace, i::Int)
    1 <= i <= space.n || throw(BoundsError(space, i))
    space.reads[] += 1
    return i
end
GR.subtree(space::TileCacheBoundSpace, inds) = space

# Helpers

conservative_block(dst, dst_inds, src, src_inds) =
    WeightBlock(
        buildweights!(WeightCOO(length(dst_inds)), Conservative(),
            dst, dst_inds, src, src_inds),
        length(dst_inds), length(src_inds))

# The generic coordinate-list assembly, reached past a method's own `pairblock`.
generic_pairblock(method, dst, dst_inds, src, src_inds) =
    invoke(GR.pairblock,
        Tuple{AbstractRegriddingMethod,RegridSpace,Any,RegridSpace,Any},
        method, dst, dst_inds, src, src_inds)

cellareas(space, inds) = [GO.area(manifold(space), getcell(space, i)) for i in inds]

"""
    T3AdoptMethod(inner)

Forward `pairblock` to `inner`, counting the forwards, and count the generic
assemblies that reach `buildweights!` anyway.
"""
struct T3AdoptMethod{M<:AbstractRegriddingMethod} <: AbstractRegriddingMethod
    inner::M
    builds::Base.RefValue{Int}
    fallbacks::Base.RefValue{Int}
end

T3AdoptMethod(inner::AbstractRegriddingMethod) =
    T3AdoptMethod(inner, Ref(0), Ref(0))

function GR.pairblock(method::T3AdoptMethod, dst_space::RegridSpace, dst_inds,
    src_space::RegridSpace, src_inds)
    method.builds[] += 1
    return GR.pairblock(method.inner, dst_space, dst_inds, src_space, src_inds)
end

function buildweights!(coo::WeightCOO, method::T3AdoptMethod,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    method.fallbacks[] += 1
    return buildweights!(coo, method.inner, dst_space, dst_inds, src_space, src_inds)
end

"""
    T3CooMethod(inner)

Forward only `buildweights!` to `inner`, counting the calls, as a third-party
emitter that knows nothing but that hook.
"""
struct T3CooMethod{M<:AbstractRegriddingMethod} <: AbstractRegriddingMethod
    inner::M
    builds::Base.RefValue{Int}
end

T3CooMethod(inner::AbstractRegriddingMethod) = T3CooMethod(inner, Ref(0))

function buildweights!(coo::WeightCOO, method::T3CooMethod,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    method.builds[] += 1
    return buildweights!(coo, method.inner, dst_space, dst_inds, src_space, src_inds)
end

@testset "Conservative weights" begin
    @test Conservative() isa AbstractRegriddingMethod

    @testset "identity" begin
        space = ToyLonLatSpace(4, 2)
        inds = ownedindices(space, 1)
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
        dst_inds = ownedindices(fine, 1)
        src_inds = ownedindices(coarse, 1)
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
        dst_inds = ownedindices(fine, 1)

        reference = conservative_block(fine, dst_inds, whole, ownedindices(whole, 1))

        # Non-contiguous chunks still use chunk-local block indices.
        @test nchunks(chunked) == 2
        parts = [ownedindices(chunked, c) for c in 1:nchunks(chunked)]
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

        # The seam every builder reaches partitions the same way. Each entry is
        # computed in exactly one part, so the columns reassemble exactly.
        seamwhole = GR.weightblock(Conservative(), fine, dst_inds, whole,
            ownedindices(whole, 1))
        seam = zeros(32, 8)
        seamdenom = zeros(32)
        for inds in parts
            block = GR.weightblock(Conservative(), fine, dst_inds, chunked, inds)
            @test size(block) == (32, length(inds))
            seam[:, inds] .+= Matrix(block.weights)
            seamdenom .+= block.denom
        end
        @test seam == Matrix(seamwhole.weights)
        @test seamdenom ≈ seamwhole.denom rtol = 1e-12
        @test Matrix(seamwhole.weights) == Matrix(reference.weights)
        @test seamwhole.denom == reference.denom
    end

    @testset "cells across the antimeridian" begin
        fine = ToyLonLatSpace(8, 4)
        # Unit-sphere geometry handles cells crossing 180°.
        band = ToyLonLatSpace(2, 2; lon = (135.0, 225.0), lat = (-45.0, 45.0))
        dst_inds = ownedindices(fine, 1)
        src_inds = ownedindices(band, 1)
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
        tiling_inds = ownedindices(tiling, 1)

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
        dst_inds, src_inds = ownedindices(coarse, 1), ownedindices(fine, 1)
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
        caps = GR.chunkextents(banded)
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
                  for i in ownedindices(banded, c)
                  for p in samples(banded, i))

        # Opposite poles lie outside the band caps.
        @test !US._contains(caps[1], toy_point(0.0, 90.0))
        @test !US._contains(caps[4], toy_point(0.0, -90.0))

        # Regional chunks retain convex caps.
        region = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (4, 2))
        @test all(c -> c.radius < Float64(pi) / 2, GR.chunkextents(region))
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

        # Sequential borrowers reuse retained, correctly typed buffers.
        cache = GR._with_sparse_assembly_cache(Float64) do cache
            push!(cache.candidate_pairs, (0, 0))
            push!(cache.rows, 0); push!(cache.cols, 0); push!(cache.vals, 0.0)
            cache
        end
        buffers = (cache.candidate_pairs, cache.rows, cache.cols, cache.vals)
        production = GR._intersectionareas(m,
            GR.subtree(dst, 1:ndst), GR.subtree(src, 1:nsrc), blockop())
        @test all(isempty, buffers)
        @test buffers === (cache.candidate_pairs, cache.rows, cache.cols, cache.vals)
        @test cache === GR._with_sparse_assembly_cache(identity, Float64)

        # Two live borrowers can never own the same cache.  Holding both until
        # the parent observes them also proves an empty pool grows instead of
        # making the second borrower wait.
        ready = Channel{Any}(2)
        release = Channel{Nothing}(2)
        jobs = map(1:2) do _
            Threads.@spawn GR._with_sparse_assembly_cache(Float64) do borrowed
                put!(ready, borrowed)
                take!(release)
                borrowed
            end
        end
        concurrent = (take!(ready), take!(ready))
        @test concurrent[1] !== concurrent[2]
        put!(release, nothing); put!(release, nothing)
        fetch.(jobs)

        # Pools are separated by the exact output value type.
        cache32 = GR._with_sparse_assembly_cache(identity, Float32)
        @test eltype(cache32.vals) === Float32
        @test cache32 !== cache

        # A failed borrower returns its cache through the production helper's
        # `finally`, ready for the next acquire.
        failed = Ref{Any}()
        @test_throws ErrorException GR._with_sparse_assembly_cache(Float64) do borrowed
            failed[] = borrowed
            error("deliberate assembly failure")
        end
        @test failed[] === GR._with_sparse_assembly_cache(identity, Float64)

        @test GR.SparseArrays.nnz(block.weights) > 2000

        # The explicitly threaded build and the production build — whichever way
        # the thread count sends the latter — both equal the serial reference.
        for other in (assemble(GOCore.True()), WeightBlock(production,
                GR._blockdenom(production, ndst)), block)
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
        inds = ownedindices(space, 3)
        tc = GR.TileCells(space, inds)
        tree = GR.subtree(tc, inds)
        @test tree isa GR.CachedCellTree
        @test all(getcell(tree, i) == getcell(space, i) for i in inds)

        # Indices outside the tile fall through to the wrapped tree.
        outside = setdiff(1:ncells(space), inds)
        @test all(getcell(tree, i) == getcell(space, i) for i in outside)

        # Only the exact tile subtree uses cached geometry.
        @test all(getcell(tc, i) == getcell(space, i) for i in 1:ncells(space))
        @test !(GR.subtree(tc, ownedindices(space, 1)) isa GR.CachedCellTree)

        # Cached geometry preserves weights bit for bit.
        other = ownedindices(space, 2)
        plain = conservative_block(space, inds, space, other)
        stored = conservative_block(tc, inds, space, other)
        @test plain.weights.colptr == stored.weights.colptr
        @test plain.weights.rowval == stored.weights.rowval
        @test all(plain.weights.nzval .=== stored.weights.nzval)
        @test all(plain.denom .=== stored.denom)
    end

    @testset "the destination tile cell cache is bounded" begin
        limit = GR._TILE_CELL_CACHE_MAX

        # The documented boundary is inclusive: a tile at the limit is
        # synthesized once, then shared by later subtree requests.
        cached_space = TileCacheBoundSpace(limit, Ref(0))
        cached = GR.TileCells(cached_space, 1:limit)
        @test GR.subtree(cached, 1:limit) isa GR.CachedCellTree
        @test GR.subtree(cached, 1:limit) isa GR.CachedCellTree
        @test cached_space.reads[] == limit
        @test length(cached.cells) == limit

        # One cell beyond the limit keeps geometry on demand and allocates no
        # tile-wide cell vector. Repeated access must not revisit initialization.
        uncached_space = TileCacheBoundSpace(limit + 1, Ref(0))
        uncached = GR.TileCells(uncached_space, 1:(limit + 1))
        @test GR.subtree(uncached, 1:(limit + 1)) === uncached_space
        @test GR.subtree(uncached, 1:(limit + 1)) === uncached_space
        @test uncached.cells === nothing
        @test uncached_space.reads[] == 0
    end

    @testset "packed cell fallback preserves indices and packing results" begin
        space = ToyLonLatSpace(8, 4)
        inds = [1, 4, 7, 10, 13, 18, 23, 28, 31]

        function leafindices!(out, node)
            if STI.isleaf(node)
                append!(out, first(e) for e in STI.child_indices_extents(node))
            else
                foreach(child -> leafindices!(out, child), STI.getchild(node))
            end
            return out
        end

        function candidates(tree, other)
            out = Tuple{Int,Int}[]
            STI.dual_depth_first_search(GO.Extents.intersects, tree, other) do i, j
                push!(out, (i, j))
            end
            return sort!(out)
        end

        other = celltree(ToyLonLatSpace(16, 8))
        reference_pairs = nothing
        reference_weights = nothing
        for capacity in (2, 3, 16)
            tree = GR.CellSpaceRTree(space, inds; nodecapacity = capacity)
            @test tree.space === space
            @test tree.indices == inds
            @test sort!(leafindices!(Int[], tree)) == sort(inds)
            @test CR.Trees.ncells(tree) == length(inds)
            @test CR.Trees.cell_index_count(tree) == ncells(space)
            @test CR.Trees.split_weight(tree) == length(inds)
            STI.isleaf(tree) || @test sum(CR.Trees.split_weight, STI.getchild(tree)) ==
                                      CR.Trees.split_weight(tree)

            pairs = candidates(tree, other)
            weights = CR.intersection_areas(manifold(space), GOCore.False(), other, tree;
                progress = false)
            if reference_pairs === nothing
                reference_pairs = pairs
                reference_weights = weights
            else
                @test pairs == reference_pairs
                @test weights == reference_weights
            end
        end
    end

    @testset "packed fallback and native cursor blocks are identical" begin
        dst = ToyLonLatSpace(8, 4)
        src = ToyLonLatSpace(16, 8)
        native_dst, native_src = celltree(dst), celltree(src)
        fallback_dst = GR.CellSpaceRTree(dst, 1:ncells(dst); nodecapacity = 3)
        fallback_src = GR.CellSpaceRTree(src, 1:ncells(src); nodecapacity = 2)

        m = manifold(dst)
        reference = CR.intersection_areas(m, GOCore.False(), native_dst, native_src;
            progress = false)
        serial = CR.intersection_areas(m, GOCore.False(), fallback_dst, fallback_src;
            progress = false)
        threaded = CR.intersection_areas(m, GOCore.True(), fallback_dst, fallback_src;
            progress = false)
        @test serial == reference
        @test threaded == serial
    end

    @testset "one tile's restricted tree, built once" begin
        # A non-chunk range on a fallback space forces the packed R-tree path.
        xd = DD.X(-168.75:22.5:168.75)
        yd = DD.Y(-78.75:22.5:78.75)
        src = RasterGrid(DD.DimArray(zeros(16, 8), (xd, yd));
            chunks = ([1:8, 9:16], [1:4, 5:8]))
        dst = CountingSpace(ToyLonLatSpace(8, 4; chunks = (8, 2)))
        inds = 5:20
        s1, s2 = ownedindices(src, 1), ownedindices(src, 4)

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
        areaonly = CR.task_local_operator(GR.IntersectionAreaOperator(m))
        @test GR._intersectionoperator(m) isa GR.IntersectionAreaOperator
        @test GR._intersectionoperator(GOCore.Planar()) isa CR.DefaultIntersectionOperator

        coarse = ToyLonLatSpace(4, 2)
        fine = ToyLonLatSpace(8, 4)
        dense = DensifiedCellSpace((10.0, 30.0), (20.0, 40.0), 5)
        cells = vcat([getcell(coarse, i) for i in 1:8],
            [getcell(fine, i) for i in 1:32], [getcell(dense, 1)])
        # Overlapping, contained, identical, and disjoint pairs, one value each:
        # every branch of the area-only kernel against the materializing one.
        checked = 0
        for a in cells, b in cells
            va = default(a, b)
            vb = areaonly(a, b)
            @test va === vb
            checked += 1
        end
        @test checked == length(cells)^2
    end

    @testset "every conservative block adopts the assembled matrix unchanged" begin
        dst = ToyLonLatSpace(4, 2)
        src = ToyLonLatSpace(8, 4)

        # The seam, against the generic coordinate-list assembly it replaces.
        fast = GR.pairblock(Conservative(), dst, 1:8, src, 1:32)
        slow = generic_pairblock(Conservative(), dst, 1:8, src, 1:32)
        @test fast.weights.colptr == slow.weights.colptr
        @test fast.weights.rowval == slow.weights.rowval
        @test all(fast.weights.nzval .=== slow.weights.nzval)
        @test all(fast.denom .=== slow.denom)
        # Both routes report denominators, and both reference them rather than
        # the row sums a method reporting none would leave.
        @test fast.reference === fast.denom
        @test slow.reference === slow.denom

        # The eager whole domain is that seam and nothing else.
        eager = GR.wholeblock(Conservative(), dst, src)
        @test eager.weights.colptr == fast.weights.colptr
        @test eager.weights.rowval == fast.weights.rowval
        @test all(eager.weights.nzval .=== fast.weights.nzval)
        @test all(eager.denom .=== fast.denom)

        # So is a chunked plan's pair, over a source chunk that is neither the
        # whole domain nor contiguous.
        chunked = ToyLonLatSpace(8, 4; chunks = (4, 2))
        plan = ChunkedPlan(Conservative(), Weighted(0.5), dst, chunked;
            storage = PerChunk())
        for s in 1:nchunks(chunked)
            sinds = ownedindices(chunked, s)
            pair = GR.buildblock(plan, 1:8, sinds)
            reference = generic_pairblock(Conservative(), dst, 1:8, chunked, sinds)
            @test pair.weights.colptr == reference.weights.colptr
            @test pair.weights.rowval == reference.weights.rowval
            @test all(pair.weights.nzval .=== reference.weights.nzval)
            @test all(pair.denom .=== reference.denom)
            @test pair.reference === pair.denom
        end
    end

    @testset "a conservative pair assembles no coordinate list" begin
        dst = ToyLonLatSpace(24, 12)
        src = ToyLonLatSpace(96, 48)
        di, si = 1:Int(ncells(dst)), 1:Int(ncells(src))
        # Serial assembly, so the two routes differ by what they allocate and
        # by nothing else.
        Base.ScopedValues.@with GR.OUTER_PARALLEL => true begin
            adopted = GR.pairblock(Conservative(), dst, di, src, si)
            generic_pairblock(Conservative(), dst, di, src, si)
            nz = GR.SparseArrays.nnz(adopted.weights)
            adoptbytes = @allocated GR.pairblock(Conservative(), dst, di, src, si)
            coobytes = @allocated generic_pairblock(Conservative(), dst, di, src, si)
            # A coordinate list holds a row, a column and a value for every
            # entry, before the second matrix assembled from it.
            @test coobytes - adoptbytes >= 24 * nz
        end
    end

    @testset "a wrapper takes the build it forwards" begin
        dst = ToyLonLatSpace(4, 2)
        src = ToyLonLatSpace(8, 4)
        field = collect(reshape(1.0:32.0, 8, 4))
        reference = GR.pairblock(Conservative(), dst, 1:8, src, 1:32)

        # Forwarding `pairblock` takes the inner method's own assembly and
        # never reaches the generic route.
        adopt = T3AdoptMethod(Conservative())
        block = GR.wholeblock(adopt, dst, src)
        @test block.weights.colptr == reference.weights.colptr
        @test block.weights.rowval == reference.weights.rowval
        @test all(block.weights.nzval .=== reference.weights.nzval)
        @test all(block.denom .=== reference.denom)
        @test adopt.builds[] == 1
        @test adopt.fallbacks[] == 0

        # Forwarding `buildweights!` alone keeps the generic route, for the
        # same values.
        emit = T3CooMethod(Conservative())
        generic = GR.wholeblock(emit, dst, src)
        @test emit.builds[] == 1
        @test generic.weights == reference.weights
        @test generic.denom == reference.denom
        @test generic.reference === generic.denom

        # One block carries one denominator, computed where it was built: a
        # reused plan applies the stored one and builds nothing further.
        counted = T3AdoptMethod(Conservative())
        plan = plan_regrid(field; to = dst, from = src, method = counted, lazy = false)
        @test counted.builds[] == 1
        @test plan.block.reference === plan.block.denom
        first = regrid(field, plan)
        @test regrid(field, plan) == first
        @test counted.builds[] == 1
        @test counted.fallbacks[] == 0
        @test first == regrid(field; to = dst, from = src, lazy = false)

        # A chunked plan makes the same choice for both wrappers.
        chunkedadopt = T3AdoptMethod(Conservative())
        chunkedemit = T3CooMethod(Conservative())
        adoptpair = GR.buildblock(ChunkedPlan(chunkedadopt, Weighted(0.5), dst, src;
            storage = PerChunk()), 1, 1)
        emitpair = GR.buildblock(ChunkedPlan(chunkedemit, Weighted(0.5), dst, src;
            storage = PerChunk()), 1, 1)
        @test adoptpair.weights == emitpair.weights
        @test adoptpair.denom == emitpair.denom
        @test chunkedadopt.builds[] == 1
        @test chunkedadopt.fallbacks[] == 0
        @test chunkedemit.builds[] == 1
    end

    @testset "empty sides keep the denominator asymmetry" begin
        dst = ToyLonLatSpace(4, 2)
        src = ToyLonLatSpace(8, 4)

        # No destination returns before denominators are declared, no source
        # after, so the two sides report differently. Both are the generic
        # route's answer, entry for entry.
        nodst = GR.weightblock(Conservative(), dst, 1:0, src, 1:32)
        @test size(nodst) == (0, 32)
        @test nodst.denom === nothing
        @test isempty(nodst.reference)

        nosrc = GR.weightblock(Conservative(), dst, 1:8, src, 1:0)
        @test size(nosrc) == (8, 0)
        @test nosrc.denom == zeros(8)
        @test nosrc.reference === nosrc.denom

        for (di, si) in ((1:0, 1:32), (1:8, 1:0), (1:0, 1:0))
            seam = GR.weightblock(Conservative(), dst, di, src, si)
            plain = generic_pairblock(Conservative(), dst, di, src, si)
            @test size(seam) == size(plain)
            @test typeof(seam.denom) === typeof(plain.denom)
            @test seam.denom == plain.denom
            @test seam.weights == plain.weights
        end
    end

    @testset "disjoint chunks keep zero denominators" begin
        north = ToyLonLatSpace(2, 1; lat = (60.0, 90.0))
        south = ToyLonLatSpace(2, 1; lat = (-90.0, -60.0))
        block = conservative_block(north, 1:2, south, 1:2)

        # Disjoint conservative blocks retain zero denominators.
        @test GR.SparseArrays.nnz(block.weights) == 0
        @test block.denom == zeros(2)

        # An adopted assembly with nothing in it reports them too, rather than
        # falling back on the row sums of no weights.
        seam = GR.weightblock(Conservative(), north, 1:2, south, 1:2)
        @test GR.SparseArrays.nnz(seam.weights) == 0
        @test seam.denom == zeros(2)
        @test seam.reference === seam.denom
    end
end
