# Nearest-cell and bilinear weight construction, and the build path the output
# sampling selects.

import DimensionalData as DD

# Cell centres in native degrees.
GR.chartaxes(space::ToyLonLatSpace) =
    ([space.lon0 + (ix - 0.5) * dlon(space) for ix in 1:space.nlon],
        [space.lat0 + (iy - 0.5) * dlat(space) for iy in 1:space.nlat])

# Return longitude on the chart axes' branch.
function GR.chartcoords(space::ToyLonLatSpace, p)
    lon, lat = toy_lonlat(p)
    return (space.lon0 + mod(lon - space.lon0, 360.0), lat)
end

GR.chartlocalindex(space::ToyLonLatSpace, ix::Int, iy::Int) = localindex(space, ix, iy)

GR.chartperiod(space::ToyLonLatSpace) =
    (space.lon1 - space.lon0 >= 360 ? 360.0 : nothing, nothing)

# Δλ bounds adjacent distance along a parallel.
GR.chartspacing(space::ToyLonLatSpace) = (deg2rad(dlon(space)), deg2rad(dlat(space)))

# Helpers

function t4_build(method, dst, dst_inds, src, src_inds)
    coo = WeightCOO(length(dst_inds))
    buildweights!(coo, method, dst, dst_inds, src, src_inds)
    return coo
end

"""
    t4_entries(coo, dst_inds, src_inds) -> Dict{(dst index, src index), weight}

Return block entries keyed by the spaces' local destination and source
indices.
"""
function t4_entries(coo, dst_inds, src_inds)
    entries = Dict{Tuple{Int,Int},Float64}()
    for k in eachindex(coo.vals)
        key = (Int(dst_inds[coo.rows[k]]), Int(src_inds[coo.cols[k]]))
        entries[key] = get(entries, key, 0.0) + coo.vals[k]
    end
    return entries
end

t4_rowsum(entries, dst_index) =
    sum(w for ((j, _), w) in entries if j == dst_index; init = 0.0)

# Parent coarse cell of a twice-refined fine cell.
function t4_parent(coarse, fine, i)
    ix, iy = cellsubscript(fine, i)
    return localindex(coarse, cld(ix, 2), cld(iy, 2))
end

struct T4NoChartSpace <: RegridSpace end

struct T4BareChartSpace <: RegridSpace end
GR.hascellchart(::T4BareChartSpace) = true

# Two methods whose sampling is a field, so one method type can be put on either
# build path. Nothing in the split may notice which type they are.

"""
    T4SplitMethod(sampling)

Report `sampling`, count `buildweights!` calls, and give each destination weight
one to the source cell of the same local index. It supplies no weights outside
that hook.
"""
struct T4SplitMethod{S<:DD.Lookups.Sampling} <: AbstractRegriddingMethod
    sampling::S
    builds::Threads.Atomic{Int}
end

T4SplitMethod(sampling::DD.Lookups.Sampling) =
    T4SplitMethod(sampling, Threads.Atomic{Int}(0))

GR.outputsampling(method::T4SplitMethod) = method.sampling

function buildweights!(coo::WeightCOO, method::T4SplitMethod,
    ::RegridSpace, dst_inds, ::RegridSpace, src_inds)
    Threads.atomic_add!(method.builds, 1)
    chunklocal_of = Dict{Int,Int}(p => k for (k, p) in enumerate(src_inds))
    for (j, p) in enumerate(dst_inds)
        k = get(chunklocal_of, p, 0)
        k == 0 && continue
        addweight!(coo, j, k, 1.0)
        adddenom!(coo, j, 1.0)
    end
    return coo
end

"""
    T4TileMethod(sampling)

Report `sampling` and supply a constant block on the point path. It implements
no `buildweights!`, so it answers only where the point path runs.
"""
struct T4TileMethod{S<:DD.Lookups.Sampling} <: AbstractRegriddingMethod
    sampling::S
end

GR.outputsampling(method::T4TileMethod) = method.sampling

GR.weightblock(::DD.Lookups.Points, ::T4TileMethod, ::RegridSpace, dst_inds,
    ::RegridSpace, src_inds) =
    WeightBlock(fill(0.5, length(dst_inds), length(src_inds)), nothing)

"""
    T4PlaceCount()

Report `Points()`, record how many destination cells each `buildweights!` call
is asked to place, and emit [`BilinearPoint`](@ref) stencils. `placed` is
therefore the number of destination point locations the builder performs, and
`calls` the number of blocks it was asked for. Its support radius is
`BilinearPoint`'s, so a plan around it discovers the same source chunks.
"""
struct T4PlaceCount <: AbstractRegriddingMethod
    placed::Threads.Atomic{Int}
    calls::Threads.Atomic{Int}
end

T4PlaceCount() = T4PlaceCount(Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))

GR.outputsampling(::T4PlaceCount) = DD.Lookups.Points()

GR.supportradius(::T4PlaceCount, src_space::RegridSpace) =
    supportradius(BilinearPoint(), src_space)

function buildweights!(coo::WeightCOO, method::T4PlaceCount,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    Threads.atomic_add!(method.calls, 1)
    Threads.atomic_add!(method.placed, length(dst_inds))
    return buildweights!(coo, BilinearPoint(), dst_space, dst_inds, src_space, src_inds)
end

@testset "Interpolation weights" begin

    @testset "NearestCell" begin
        # Self-regridding is the identity.
        space = ToyLonLatSpace(8, 4)
        inds = ownedindices(space, 1)
        block = WeightBlock(t4_build(NearestCell(), space, inds, space, inds),
            length(inds), length(inds))
        @test Matrix(block.weights) == Matrix(LinearAlgebra.I, 32, 32)

        # Point samples have no denominator.
        @test block.denom === nothing

        # Each fine centroid selects its coarse parent.
        coarse = ToyLonLatSpace(4, 2)
        fine = ToyLonLatSpace(8, 4)
        dst_inds, src_inds = ownedindices(fine, 1), ownedindices(coarse, 1)
        entries = t4_entries(t4_build(NearestCell(), fine, dst_inds, coarse, src_inds),
            dst_inds, src_inds)
        @test length(entries) == ncells(fine)
        @test all(entries[(i, t4_parent(coarse, fine, i))] == 1.0
                  for i in 1:ncells(fine))

        # Centroids outside coverage emit no weight.
        north = ToyLonLatSpace(4, 1; lat = (0.0, 90.0))
        global_dst = ToyLonLatSpace(4, 2)
        dst_inds = ownedindices(global_dst, 1)
        entries = t4_entries(
            t4_build(NearestCell(), global_dst, dst_inds, north, ownedindices(north, 1)),
            dst_inds, ownedindices(north, 1))
        @test sort(unique(first.(keys(entries)))) == collect(5:8)
    end

    @testset "BilinearPoint stencils" begin
        # Source centres are at ±135°/±45° longitude and ±45° latitude.
        src = ToyLonLatSpace(4, 2)
        src_inds = ownedindices(src, 1)

        # Asymmetric fractional position verifies axis order and weights.
        dst = ToyLonLatSpace(1, 1; lon = (-27.5, -17.5), lat = (17.5, 27.5))
        entries = t4_entries(t4_build(BilinearPoint(), dst, [1], src, src_inds),
            [1], src_inds)
        @test length(entries) == 4
        @test entries[(1, localindex(src, 2, 1))] ≈ 0.1875
        @test entries[(1, localindex(src, 3, 1))] ≈ 0.0625
        @test entries[(1, localindex(src, 2, 2))] ≈ 0.5625
        @test entries[(1, localindex(src, 3, 2))] ≈ 0.1875
        @test t4_rowsum(entries, 1) ≈ 1.0

        # A centroid at 180° interpolates across the periodic seam.
        seam = ToyLonLatSpace(1, 1; lon = (175.0, 185.0), lat = (-5.0, 5.0))
        entries = t4_entries(t4_build(BilinearPoint(), seam, [1], src, src_inds),
            [1], src_inds)
        @test length(entries) == 4
        @test all(entries[(1, localindex(src, ix, iy))] ≈ 0.25
                  for ix in (1, 4), iy in (1, 2))

        # Outside latitude centres, clamp latitude and interpolate longitude.
        polar = ToyLonLatSpace(1, 1; lon = (-5.0, 5.0), lat = (75.0, 85.0))
        entries = t4_entries(t4_build(BilinearPoint(), polar, [1], src, src_inds),
            [1], src_inds)
        @test length(entries) == 2
        @test entries[(1, localindex(src, 2, 2))] ≈ 0.5
        @test entries[(1, localindex(src, 3, 2))] ≈ 0.5

        # Non-periodic corners clamp to one point.
        patch = ToyLonLatSpace(4, 2; lon = (-40.0, 40.0), lat = (-20.0, 20.0))
        corner = ToyLonLatSpace(1, 1; lon = (38.0, 40.0), lat = (18.0, 20.0))
        entries = t4_entries(
            t4_build(BilinearPoint(), corner, [1], patch, ownedindices(patch, 1)),
            [1], ownedindices(patch, 1))
        @test entries == Dict((1, localindex(patch, 4, 2)) => 1.0)

        # Bilinear interpolation reproduces a linear chart field.
        fsrc = ToyLonLatSpace(36, 18)
        fdst = ToyLonLatSpace(17, 8; lon = (-170.0, 170.0), lat = (-80.0, 80.0))
        fdst_inds, fsrc_inds = ownedindices(fdst, 1), ownedindices(fsrc, 1)
        function linear(p)
            lon, lat = toy_lonlat(p)
            return 2.0 + 0.01 * lon + 0.03 * lat
        end
        field = [linear(cellcentroid(fsrc, i)) for i in fsrc_inds]
        weights = WeightBlock(
            t4_build(BilinearPoint(), fdst, fdst_inds, fsrc, fsrc_inds),
            length(fdst_inds), length(fsrc_inds)).weights
        @test weights * field ≈ [linear(cellcentroid(fdst, i)) for i in fdst_inds]
        @test all(≈(1.0), sum(weights; dims = 2))

        # Bilinear interpolation requires a complete chart interface.
        @test_throws "hascellchart" buildweights!(
            WeightCOO(1), BilinearPoint(), src, [1], T4NoChartSpace(), [1])
        @test_throws "chartaxes" buildweights!(
            WeightCOO(1), BilinearPoint(), src, [1], T4BareChartSpace(), [1])
    end

    @testset "BilinearPoint sampling policy" begin
        # Interior: the stencil is a tensor product, so it reproduces the whole
        # Q1 space `a + bx + cy + dxy` on the chart, not merely affine fields.
        # Source centres lie 20° apart; every destination centroid falls
        # strictly inside the lattice hull, so nothing here is clamped.
        qsrc = ToyLonLatSpace(6, 5; lon = (-60.0, 60.0), lat = (-50.0, 50.0))
        qdst = ToyLonLatSpace(7, 4; lon = (-49.0, 49.0), lat = (-39.0, 39.0))
        qsrc_inds, qdst_inds = ownedindices(qsrc, 1), ownedindices(qdst, 1)
        q1(x, y) = 2.0 + 0.01x + 0.03y + 0.0007x * y
        chart(space, i) = GR.chartcoords(qsrc, cellcentroid(space, i))
        qfield = [q1(chart(qsrc, i)...) for i in qsrc_inds]
        qweights = WeightBlock(
            t4_build(BilinearPoint(), qdst, qdst_inds, qsrc, qsrc_inds),
            length(qdst_inds), length(qsrc_inds)).weights
        @test qweights * qfield ≈ [q1(chart(qdst, i)...) for i in qdst_inds]
        @test all(≈(1.0), sum(qweights; dims = 2))

        # Four entries per interior destination, never more.
        @test all(count(!iszero, view(qweights, r, :)) == 4
                  for r in axes(qweights, 1))

        # Periodic seam: a point 35° east of the last centre interpolates
        # between that centre and the first across the 90° seam, so the two
        # longitude weights are 11/18 and 7/18 rather than a clamp to one.
        src = ToyLonLatSpace(4, 2)
        src_inds = ownedindices(src, 1)
        seam = ToyLonLatSpace(1, 1; lon = (169.0, 171.0), lat = (-1.0, 1.0))
        entries = t4_entries(t4_build(BilinearPoint(), seam, [1], src, src_inds),
            [1], src_inds)
        @test length(entries) == 4
        @test all(entries[(1, localindex(src, 4, iy))] ≈ 11 / 36 for iy in (1, 2))
        @test all(entries[(1, localindex(src, 1, iy))] ≈ 7 / 36 for iy in (1, 2))
        @test t4_rowsum(entries, 1) ≈ 1.0

        # Non-periodic rim: past the last centre of one axis the stencil drops
        # to that axis' rim line and the other axis keeps interpolating, so two
        # entries carry the other axis' weights and the row still sums to one.
        patch = ToyLonLatSpace(4, 2; lon = (-40.0, 40.0), lat = (-20.0, 20.0))
        patch_inds = ownedindices(patch, 1)
        latrim = ToyLonLatSpace(1, 1; lon = (-6.0, -4.0), lat = (14.0, 16.0))
        entries = t4_entries(t4_build(BilinearPoint(), latrim, [1], patch, patch_inds),
            [1], patch_inds)
        @test length(entries) == 2
        @test entries[(1, localindex(patch, 2, 2))] ≈ 0.75
        @test entries[(1, localindex(patch, 3, 2))] ≈ 0.25
        @test t4_rowsum(entries, 1) ≈ 1.0

        lonrim = ToyLonLatSpace(1, 1; lon = (34.0, 36.0), lat = (2.0, 4.0))
        entries = t4_entries(t4_build(BilinearPoint(), lonrim, [1], patch, patch_inds),
            [1], patch_inds)
        @test length(entries) == 2
        @test entries[(1, localindex(patch, 4, 1))] ≈ 0.35
        @test entries[(1, localindex(patch, 4, 2))] ≈ 0.65
        @test t4_rowsum(entries, 1) ≈ 1.0

        # Outside coverage: the chart is located but never tested for coverage,
        # so a point 40° north of the source's northmost cell is served by the
        # rim line with a row summing to one. `NearestCell`, which asks
        # `cellat`, emits nothing for the same point.
        outside = ToyLonLatSpace(1, 1; lon = (-1.0, 1.0), lat = (59.0, 61.0))
        @test cellat(patch, cellcentroid(outside, 1)) === nothing
        entries = t4_entries(t4_build(BilinearPoint(), outside, [1], patch, patch_inds),
            [1], patch_inds)
        @test length(entries) == 2
        @test entries[(1, localindex(patch, 2, 2))] ≈ 0.5
        @test entries[(1, localindex(patch, 3, 2))] ≈ 0.5
        @test t4_rowsum(entries, 1) ≈ 1.0
        @test isempty(t4_entries(
            t4_build(NearestCell(), outside, [1], patch, patch_inds),
            [1], patch_inds))

        # The same holds along a non-periodic longitude axis: whatever branch
        # the chart reports, a coordinate past the axis span is clamped rather
        # than refused.
        west = ToyLonLatSpace(1, 1; lon = (-171.0, -169.0), lat = (-1.0, 1.0))
        @test GR.chartcoords(patch, cellcentroid(west, 1))[1] ≈ 190.0
        @test cellat(patch, cellcentroid(west, 1)) === nothing
        entries = t4_entries(t4_build(BilinearPoint(), west, [1], patch, patch_inds),
            [1], patch_inds)
        @test length(entries) == 2
        @test all(entries[(1, localindex(patch, 4, iy))] ≈ 0.5 for iy in (1, 2))
        @test t4_rowsum(entries, 1) ≈ 1.0
    end

    @testset "stencils partition across source chunks" begin
        # Every stencil crosses the east/west source-chunk boundary.
        src = ToyLonLatSpace(4, 2; chunks = (2, 2))
        dst = ToyLonLatSpace(3, 2; lon = (-30.0, 30.0), lat = (-30.0, 30.0))
        dst_inds = ownedindices(dst, 1)
        @test nchunks(src) == 2

        whole_inds = ownedindices(ToyLonLatSpace(4, 2), 1)
        whole = t4_entries(t4_build(BilinearPoint(), dst, dst_inds, src, whole_inds),
            dst_inds, whole_inds)

        blocks = [t4_entries(
                      t4_build(BilinearPoint(), dst, dst_inds, src, ownedindices(src, c)),
                      dst_inds, ownedindices(src, c)) for c in 1:nchunks(src)]

        # Each non-empty block emits only its own source cells.
        @test all(!isempty(b) for b in blocks)
        @test all(key[2] in ownedindices(src, c)
                  for (c, b) in enumerate(blocks) for key in keys(b))

        # Chunked stencils merge to the unchunked stencil.
        merged = merge(+, blocks...)
        @test keys(merged) == keys(whole)
        @test all(merged[key] ≈ whole[key] for key in keys(whole))

        # Merged stencils still sum to one.
        @test all(t4_rowsum(merged, i) ≈ 1.0 for i in dst_inds)
    end

    @testset "support radius" begin
        # Bilinear support reaches beyond geometric overlap.
        space = ToyLonLatSpace(8, 4)
        @test supportradius(NearestCell(), space) == 0.0
        @test supportradius(BilinearPoint(), space) > 0
        @test supportradius(BilinearPoint(), space) ≈ deg2rad(45.0)

        # The larger axis spacing bounds both directions.
        oblong = ToyLonLatSpace(36, 6)
        @test supportradius(BilinearPoint(), oblong) >= deg2rad(dlon(oblong))
        @test supportradius(BilinearPoint(), oblong) >= deg2rad(dlat(oblong))
    end

    @testset "output sampling selects the build path" begin
        src = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (4, 2))
        dst = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (8, 2))
        field = collect(reshape(1.0:32.0, 8, 4))
        ndst, nsrc = Int(ncells(dst)), Int(ncells(src))

        # A point method that supplies no weights of its own still reaches
        # `buildweights!`, once per block, and keeps every entry it emits.
        point = T4SplitMethod(DD.Lookups.Points())
        area = T4SplitMethod(DD.Lookups.Intervals(DD.Lookups.Center()))
        pointweights = GR.wholeblock(point, dst, src).weights
        @test Matrix(pointweights) == Matrix(LinearAlgebra.I, ndst, nsrc)
        @test point.builds[] == 1
        @test pointweights == GR.wholeblock(area, dst, src).weights

        # Eager and chunked agree for it, so the split moved no value.
        eager = regrid(field; to = dst, from = src, method = point, lazy = false)
        dest = Vector{Float64}(undef, ndst)
        regrid!(dest, field,
            ChunkedPlan(point, Weighted(0.5), dst, src; storage = PerChunk()))
        @test dest == eager

        # One method type on two traits. Only the point one takes the point
        # path, in the eager builder and in the chunk-pair builder alike; the
        # area one falls to `buildweights!`, which it does not implement.
        tilepoint = T4TileMethod(DD.Lookups.Points())
        tilearea = T4TileMethod(DD.Lookups.Intervals(DD.Lookups.Center()))
        @test Matrix(GR.wholeblock(tilepoint, dst, src).weights) ==
              fill(0.5, ndst, nsrc)
        @test_throws "buildweights! is not implemented" GR.wholeblock(
            tilearea, dst, src)

        tileplan = ChunkedPlan(tilepoint, Weighted(0.5), dst, src;
            storage = PerChunk())
        @test Matrix(GR.buildblock(tileplan, 1, 1).weights) ==
              fill(0.5, length(ownedindices(dst, 1)), length(ownedindices(src, 1)))
        @test_throws "buildweights! is not implemented" GR.buildblock(
            ChunkedPlan(tilearea, Weighted(0.5), dst, src; storage = PerChunk()),
            1, 1)

        # The shipped point methods answer the same eagerly and by chunk.
        pdst = ToyLonLatSpace(5, 3; lon = (-30.0, 30.0), lat = (-15.0, 15.0),
            chunks = (5, 2))
        for method in (NearestCell(), BilinearPoint())
            reference = regrid(field; to = pdst, from = src, method, lazy = false)
            out = Vector{Float64}(undef, Int(ncells(pdst)))
            regrid!(out, field, ChunkedPlan(method, Weighted(0.5), pdst, src;
                storage = PerChunk()))
            @test out == reference
        end
    end

    @testset "point lookup repeats once per candidate source chunk" begin
        # One destination tile against a chunked source. The tile is the whole
        # destination, so it is destination chunk 1 and the plan's relation
        # names its candidate source chunks in one row. No destination centroid
        # shares a source centre, so every stencil is a full four entries.
        src = ToyLonLatSpace(36, 18; chunks = (6, 6))
        dst = ToyLonLatSpace(20, 10; lon = (-19.0, 21.0), lat = (-20.0, 20.0))
        @test nchunks(dst) == 1
        tile = ownedindices(dst, 1)

        method = T4PlaceCount()
        plan = ChunkedPlan(method, Weighted(0.5), dst, src; storage = PerChunk())
        candidates = GR.sourcesof(GR.dependencies(plan), 1)
        k = length(candidates)

        # The tile reaches several source chunks but not the whole source, so
        # `k` measures a candidate set rather than the chunk count.
        @test 1 < k < nchunks(src)

        # Driving every pair the relation names builds one block per candidate.
        contributing = 0
        nonzeros = 0
        for s in candidates
            block = GR.blockfor(plan, (1, Int(s)), tile).block
            n = count(!iszero, block.weights)
            nonzeros += n
            n > 0 && (contributing += 1)
        end
        @test method.calls[] == k

        # Every pair places every destination cell of the tile: the point
        # location a stencil needs is repeated once per candidate chunk.
        @test method.placed[] == length(tile) * k

        # What the repeated work produces is one four-entry stencil per
        # destination, partitioned across the chunks that own its sample sites
        # — fewer chunks than were searched.
        @test nonzeros == 4 * length(tile)
        @test 0 < contributing < k
    end
end
