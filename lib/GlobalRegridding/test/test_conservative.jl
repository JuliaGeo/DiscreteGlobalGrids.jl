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

"""
    CellCountSpace(space)

Wrap a lon/lat space and count both restricted-tree builds and cell
syntheses. The tile's tree is rooted at the wrapper and its node extents come
from cell corners, so the count is exactly what a weight build synthesizes.
Block builders may run on `Threads.@spawn` tasks, hence the atomic counters.
"""
struct CellCountSpace{S<:ToyLonLatSpace} <: RegridSpace
    space::S
    cells::Threads.Atomic{Int}
    builds::Threads.Atomic{Int}
end

CellCountSpace(space::ToyLonLatSpace) =
    CellCountSpace(space, Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))

ncells(cs::CellCountSpace) = ncells(cs.space)
nchunks(cs::CellCountSpace) = nchunks(cs.space)
ownedindices(cs::CellCountSpace, chunk::Int) = ownedindices(cs.space, chunk)
manifold(cs::CellCountSpace) = manifold(cs.space)
hascellchart(cs::CellCountSpace) = hascellchart(cs.space)
cellcentroid(cs::CellCountSpace, i::Int) = cellcentroid(cs.space, i)
cellat(cs::CellCountSpace, p) = cellat(cs.space, p)
chunkextents(cs::CellCountSpace) = GR.chunkextents(cs.space)

getcell(cs::CellCountSpace, i::Int) =
    (Threads.atomic_add!(cs.cells, 1); getcell(cs.space, i))

celltree(cs::CellCountSpace) = _countedtree(cs, 1:ncells(cs.space))

GR.subtree(cs::CellCountSpace, inds) =
    (Threads.atomic_add!(cs.builds, 1); _countedtree(cs, inds))

_countedtree(cs::CellCountSpace, inds) =
    ToyCapTree(cs, collect(Int, inds),
        [toy_cap(cellcorners(cs.space, Int(i))) for i in inds])

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

    @testset "a prepared destination changes no weight" begin
        # Both destination kinds: a raster lattice and a cell space.
        raster = RasterGrid(DD.DimArray(zeros(12, 6),
            (DD.X(range(-165.0, 165.0; length = 12)),
                DD.Y(range(-75.0, 75.0; length = 6)))); chunks = ([1:6, 7:12], [1:3, 4:6]))
        cells = ToyLonLatSpace(12, 6; chunks = (6, 3))
        for (dst, src) in ((raster, raster),
            (cells, ToyLonLatSpace(16, 8; chunks = (8, 4))))
            inds = ownedindices(dst, 3)
            sinds = ownedindices(src, 2)
            cache = GR.DestinationCache(dst, inds)
            prepared = GR.pairblock(Conservative(), dst, cache, src, sinds)
            plain = GR.pairblock(Conservative(), dst, inds, src, sinds)
            @test prepared.weights == plain.weights
            @test prepared.weights.colptr == plain.weights.colptr
            @test GR.SparseArrays.rowvals(prepared.weights) ==
                  GR.SparseArrays.rowvals(plain.weights)
            @test all(GR.SparseArrays.nzrange(prepared.weights, c) ==
                      GR.SparseArrays.nzrange(plain.weights, c)
                      for c in axes(prepared.weights, 2))
            @test all(prepared.weights.nzval .=== plain.weights.nzval)
            @test all(prepared.denom .=== plain.denom)
            @test prepared.reference === prepared.denom

            # And the generic coordinate-list route, which prepares nothing.
            reference = generic_pairblock(Conservative(), dst, inds, src, sinds)
            @test prepared.weights.rowval == reference.weights.rowval
            @test all(prepared.weights.nzval .=== reference.weights.nzval)

            # An empty source keeps the degenerate contract of the index form.
            empty_prepared = GR.pairblock(Conservative(), dst, cache, src, 1:0)
            empty_plain = GR.pairblock(Conservative(), dst, inds, src, 1:0)
            @test empty_prepared.denom == empty_plain.denom
            @test size(empty_prepared) == size(empty_plain)
        end
    end

    @testset "a prepared destination is built once and synthesized once" begin
        # A non-chunk range on a fallback space forces the restricted tree.
        xd = DD.X(-168.75:22.5:168.75)
        yd = DD.Y(-78.75:22.5:78.75)
        src = RasterGrid(DD.DimArray(zeros(16, 8), (xd, yd));
            chunks = ([1:8, 9:16], [1:4, 5:8]))
        dst = CellCountSpace(ToyLonLatSpace(8, 4; chunks = (8, 2)))
        inds = 5:20
        # The source's northern half, so the tile holds southern destination
        # cells neither chunk reaches.
        s1, s2 = ownedindices(src, 3), ownedindices(src, 4)

        cache = GR.DestinationCache(dst, inds)
        @test cache.inds === inds
        @test length(cache.polygons) == length(inds)
        @test !any(cache.filled)
        @test dst.builds[] == 1
        built = dst.cells[]

        sh1 = GR.pairblock(Conservative(), dst, cache, src, s1)
        sh2 = GR.pairblock(Conservative(), dst, cache, src, s2)

        # One restricted tree for the tile, whatever its block count.
        @test dst.builds[] == 1
        # One synthesis per prepared row, over both blocks together.
        @test dst.cells[] - built == count(cache.filled)
        # Every row a block weighs was prepared.
        @test all(cache.filled[r] for r in GR.SparseArrays.rowvals(sh1.weights))
        @test all(cache.filled[r] for r in GR.SparseArrays.rowvals(sh2.weights))
        # Rows the tile holds but neither source chunk reaches are never
        # synthesized.
        @test count(cache.filled) == 12
        @test !any(cache.filled[1:4])

        # Fresh preparations per block: the unprepared reference.
        b1 = dst.builds[]
        fr1 = GR.pairblock(Conservative(), dst, inds, src, s1)
        fr2 = GR.pairblock(Conservative(), dst, inds, src, s2)
        @test dst.builds[] - b1 == 2

        for (a, b) in ((sh1, fr1), (sh2, fr2))
            @test a.weights.colptr == b.weights.colptr
            @test a.weights.rowval == b.weights.rowval
            @test all(a.weights.nzval .=== b.weights.nzval)
            @test all(a.denom .=== b.denom)
        end
    end

    @testset "concurrent blocks share one prepared destination" begin
        xd = DD.X(-168.75:22.5:168.75)
        yd = DD.Y(-78.75:22.5:78.75)
        src = RasterGrid(DD.DimArray(zeros(16, 8), (xd, yd));
            chunks = ([1:8, 9:16], [1:4, 5:8]))
        dst = CellCountSpace(ToyLonLatSpace(8, 4; chunks = (8, 2)))
        inds = 1:32
        chunkinds = [ownedindices(src, s) for s in 1:nchunks(src)]

        cache = GR.DestinationCache(dst, inds)
        built = dst.cells[]
        tasks = map(chunkinds) do sinds
            Base.ScopedValues.@with GR.OUTER_PARALLEL => true begin
                Threads.@spawn GR.pairblock(Conservative(), dst, cache, src, sinds)
            end
        end
        blocks = [fetch(t)::WeightBlock for t in tasks]

        @test dst.builds[] == 1
        @test dst.cells[] - built == count(cache.filled)
        # The whole tile is covered here, so every row was prepared exactly once.
        @test count(cache.filled) == length(inds)

        for (block, sinds) in zip(blocks, chunkinds)
            reference = GR.pairblock(Conservative(), dst, inds, src, sinds)
            @test block.weights.colptr == reference.weights.colptr
            @test block.weights.rowval == reference.weights.rowval
            @test all(block.weights.nzval .=== reference.weights.nzval)
            @test all(block.denom .=== reference.denom)
        end
    end

    @testset "prepared geometry is charged to the budget, never to a tile's size" begin
        cells = ToyLonLatSpace(32, 16)
        raster = RasterGrid(DD.DimArray(zeros(12, 6),
            (DD.X(range(-165.0, 165.0; length = 12)),
                DD.Y(range(-75.0, 75.0; length = 6)))))
        ncell = Int(ncells(cells))

        # A raster cell is cheap to synthesize; an unknown space's is not.
        @test !GR.expensivecellgeometry(raster)
        @test GR.expensivecellgeometry(cells)
        @test GR.preparesdestination(Conservative(), cells)
        @test !GR.preparesdestination(Conservative(), raster)
        @test !GR.preparesdestination(BarycentricPoint(), cells)

        # Only a destination whose geometry is kept charges for it.
        percell = GR.destcellbytes(Conservative(), cells, ncell)
        @test percell > GR.DEST_BYTES_PER_CELL
        @test GR.destcellbytes(Conservative(), raster, Int(ncells(raster))) ==
              GR.DEST_BYTES_PER_CELL
        @test GR.destcellbytes(BarycentricPoint(), cells, ncell) ==
              GR.DEST_BYTES_PER_CELL

        # A tile is sized by the executor's own per-cell cost alone, whatever
        # the destination's geometry costs to keep.
        budget = 1 << 20
        @test GR._defaulttilesizes(1 << 20, 1, budget) ==
              [fld(GR.destcellbudget(budget), GR.DEST_BYTES_PER_CELL)]

        # Preparing takes what is left of the tile's share, and is refused
        # where it does not fit. A tile at the size cap cannot hold it.
        capped = only(GR._defaulttilesizes(1 << 20, 1, budget))
        @test !GR.destcellsfit(Conservative(), cells, capped, budget)
        @test GR.destcellsfit(Conservative(), cells,
            fld(GR.destcellbudget(budget), percell), budget)

        # The same rule decides the index set a build is handed.
        @test GR.preparedestination(Conservative(), cells, 1:ncell, 1 << 10) === 1:ncell
        @test GR.preparedestination(Conservative(), cells, 1:ncell,
            GR.DEFAULT_BUDGET) isa GR.DestinationCache
        # A cheap destination and a point method prepare nothing at any budget.
        @test GR.preparedestination(Conservative(), raster, 1:24,
            GR.DEFAULT_BUDGET) === 1:24
        @test GR.preparedestination(BarycentricPoint(), cells, 1:ncell,
            GR.DEFAULT_BUDGET) === 1:ncell
        # An empty tile has no cell to probe, and prepares nothing.
        @test GR.DestinationCache(cells, 1:0) === nothing
        @test GR.preparedestination(Conservative(), cells, 1:0,
            GR.DEFAULT_BUDGET) === 1:0
    end

    @testset "a method that reads no prepared geometry takes the tile's cells" begin
        dst = ToyLonLatSpace(4, 2)
        src = ToyLonLatSpace(8, 4)
        inds, sinds = 1:8, 1:32
        cache = GR.DestinationCache(dst, inds)

        method = T3CooMethod(Conservative())
        block = GR.pairblock(method, dst, cache, src, sinds)
        reference = GR.pairblock(method, dst, inds, src, sinds)
        @test block.weights == reference.weights
        @test all(block.denom .=== reference.denom)
        # It names the tile's cells through the space, and fills no slot.
        @test !any(cache.filled)
    end

    @testset "a prepared destination adds nothing to the clipping loop" begin
        dst = ToyLonLatSpace(8, 4)
        src = ToyLonLatSpace(16, 8)
        inds, sinds = 1:Int(ncells(dst)), 1:Int(ncells(src))
        m = manifold(dst)
        src_tree = GR.subtree(src, sinds)

        cache = GR.DestinationCache(dst, inds)
        pairs = [(i1, i2) for i1 in sinds, i2 in inds]
        CR.work_items(GR.BlockAreaOperator(GR._intersectionoperator(m), cache.map,
            GR.indexmap(sinds), nothing, cache), pairs)
        @test count(cache.filled) == length(inds)

        cachedop = CR.task_local_operator(
            GR.BlockAreaOperator(GR._intersectionoperator(m), cache.map,
                GR.indexmap(sinds), GR._cellmemo(src, sinds), cache))
        memoop = CR.task_local_operator(
            GR.BlockAreaOperator(GR._intersectionoperator(m), GR.indexmap(inds),
                GR.indexmap(sinds), GR._cellmemo(src, sinds),
                GR._cellmemo(dst, inds)))
        rows, cols, vals = Int[], Int[], Float64[]
        sizehint!(rows, 8)
        sizehint!(cols, 8)
        sizehint!(vals, 8)
        item = (1, 1)
        # Warm both operators' memos and both vectors' capacity.
        cachedop(rows, cols, vals, item, src_tree, cache.tree)
        memoop(rows, cols, vals, item, src_tree, cache.tree)
        cachedbytes = @allocated cachedop(rows, cols, vals, item, src_tree, cache.tree)
        memobytes = @allocated memoop(rows, cols, vals, item, src_tree, cache.tree)
        @test cachedbytes == 0
        @test cachedbytes <= memobytes
        @test length(rows) == 4
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
