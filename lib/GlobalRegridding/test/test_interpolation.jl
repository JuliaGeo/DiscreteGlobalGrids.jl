# Point-sampling construction and output-sampling dispatch.

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

# Storing sampling in the method tests value-based `weightblock` dispatch.

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

Test sampling-specific [`weightblock`](@ref) dispatch with a constant `Points()`
block. The type intentionally omits `buildweights!` and [`sampler`](@ref): its
specialized sampling succeeds, other samplings reach the generic error, and
plans retain chunk-pair construction.
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

Emit [`BarycentricPoint`](@ref) stencils while counting construction work.
`placed` counts destination locations and `calls` counts requested blocks. The
matching support radius preserves source-chunk discovery.
"""
struct T4PlaceCount <: AbstractRegriddingMethod
    placed::Threads.Atomic{Int}
    calls::Threads.Atomic{Int}
end

T4PlaceCount() = T4PlaceCount(Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))

GR.outputsampling(::T4PlaceCount) = DD.Lookups.Points()

GR.supportradius(::T4PlaceCount, src_space::RegridSpace) =
    supportradius(BarycentricPoint(), src_space)

function buildweights!(coo::WeightCOO, method::T4PlaceCount,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    Threads.atomic_add!(method.calls, 1)
    Threads.atomic_add!(method.placed, length(dst_inds))
    return buildweights!(coo, BarycentricPoint(), dst_space, dst_inds, src_space, src_inds)
end

"""
    t5_bracket(space::ToyLonLatSpace, p) -> ((i0, w0), (i1, w1))

Return the two longitude sample sites bracketing `p` in its latitude row, with
linear weights. Longitude wrapping joins the first and last columns; a
one-column source therefore names the same cell twice.
"""
function t5_bracket(space::ToyLonLatSpace, p)
    lon, lat = toy_lonlat(p)
    iy = clamp(floor(Int, (lat - space.lat0) / dlat(space)) + 1, 1, space.nlat)
    # Centres sit half a cell in from the axis start, so `t` is where the point
    # lies between the two of them that bracket it.
    t = mod((lon - space.lon0) / dlon(space) - 0.5, Float64(space.nlon))
    ix = floor(Int, t)
    frac = t - ix
    return ((localindex(space, mod(ix, space.nlon) + 1, iy), 1.0 - frac),
        (localindex(space, mod(ix + 1, space.nlon) + 1, iy), frac))
end

"""
    T5Bracket(space, placed, yielding)

The sampler state [`T5PlaceCount`](@ref) prepares once per source space: the
space its stencils are written on, the method's location counter, and whether
the first location of a pass yields.
"""
struct T5Bracket
    space::ToyLonLatSpace
    placed::Threads.Atomic{Int}
    yielding::Bool
end

"""
    T5PlaceCount(; yielding = false)

Provide a point sampler that emits the two source sites bracketing each
destination. `placed` counts stencil queries across all produced blocks.
`yielding = true` pauses the first query so another task can observe an active
build.
"""
struct T5PlaceCount <: AbstractRegriddingMethod
    placed::Threads.Atomic{Int}
    yielding::Bool
end

T5PlaceCount(; yielding::Bool = false) = T5PlaceCount(Threads.Atomic{Int}(0), yielding)

GR.outputsampling(::T5PlaceCount) = DD.Lookups.Points()

GR.supportradius(::T5PlaceCount, src_space::RegridSpace) =
    GR.chartradius(src_space)

GR.sampler(method::T5PlaceCount, space::ToyLonLatSpace) =
    GR.Sampler(space, GR.samplesites(space), T5Bracket(space, method.placed,
        method.yielding), method)

function GR.weightsat!(row::GR.WeightRow,
    s::GR.Sampler{<:RegridSpace,<:AbstractVector,T5Bracket}, p)
    empty!(row)
    state = s.state
    Threads.atomic_add!(state.placed, 1) == 0 && state.yielding && yield()
    for (i, w) in t5_bracket(state.space, p)
        w > 0 || continue
        push!(row.indices, i)
        push!(row.weights, w)
    end
    return GR.WeightsMapped
end

# Eager and uncached pair builds need the same stencil oracle.
function buildweights!(coo::WeightCOO, method::T5PlaceCount, dst_space::RegridSpace,
    dst_inds, src_space::RegridSpace, src_inds)
    smp = GR.sampler(method, src_space)
    sites = GR.samplesites(dst_space)
    indexer = GR.indexmap(src_inds)
    row = GR.WeightRow()
    for (j, i) in enumerate(dst_inds)
        GR.ismapped(GR.weightsat!(row, smp, sites[Int(i)])) || continue
        for k in eachindex(row.indices)
            c = GR.localindex(indexer, row.indices[k])
            c == 0 && continue
            addweight!(coo, j, c, row.weights[k])
        end
    end
    return coo
end

"""
    T6LocateCount(space)

Wrap a source space and count [`cellat`](@ref) queries while delegating the
remaining interface. The atomic `located` counter covers builds running on
spawned tasks.
"""
struct T6LocateCount{S<:RegridSpace} <: RegridSpace
    space::S
    located::Threads.Atomic{Int}
end

T6LocateCount(space::RegridSpace) = T6LocateCount(space, Threads.Atomic{Int}(0))

function cellat(cs::T6LocateCount, p)
    Threads.atomic_add!(cs.located, 1)
    return cellat(cs.space, p)
end

ncells(cs::T6LocateCount) = ncells(cs.space)
getcell(cs::T6LocateCount, i::Int) = getcell(cs.space, i)
manifold(cs::T6LocateCount) = manifold(cs.space)
hascellchart(cs::T6LocateCount) = hascellchart(cs.space)
cellcentroid(cs::T6LocateCount, i::Int) = cellcentroid(cs.space, i)
nchunks(cs::T6LocateCount) = nchunks(cs.space)
ownedindices(cs::T6LocateCount, chunk::Int) = ownedindices(cs.space, chunk)
celltree(cs::T6LocateCount) = celltree(cs.space)
chunkextents(cs::T6LocateCount) = chunkextents(cs.space)

"""
    t6_owners(dst, tile, src) -> Vector{Int}

The source chunks holding the cells `tile`'s destinations sit in, read from
`cellat` directly so no build takes part in the answer.
"""
t6_owners(dst, tile, src) =
    sort(unique(Int(GR.chunkat(src, cellat(src, cellcentroid(dst, Int(i)))))
                for i in tile))

# Direct construction keeps plan storage outside this oracle.
t5_weights(method, dst, tile, src) =
    GR.tileweights(method, dst, tile, src,
        GR.sampler(method, src))

"""
    t5_entries(dst, tile, src) -> Dict{(tile row, source local index), weight}

The stencils `tile`'s destinations take from `src`, computed from the bracket
directly so a build has nothing to do with them.
"""
function t5_entries(dst, tile, src)
    out = Dict{Tuple{Int,Int},Float64}()
    for (j, i) in enumerate(tile), (s, w) in t5_bracket(src, cellcentroid(dst, i))
        w > 0 || continue
        out[(j, s)] = get(out, (j, s), 0.0) + w
    end
    return out
end

# The entries a tile's blocks hold, read back in the source space's own indices.
function t5_blockentries(weights, src)
    out = Dict{Tuple{Int,Int},Float64}()
    for (k, s) in enumerate(weights.sourcechunks)
        inds = ownedindices(src, s)
        M = Matrix(weights.blocks[k].weights)
        for row in axes(M, 1), col in axes(M, 2)
            iszero(M[row, col]) && continue
            out[(row, Int(inds[col]))] = M[row, col]
        end
    end
    return out
end

# The source chunks a tile's stencils name, under whatever chunking `src` has.
t5_owners(dst, tile, src) =
    sort(unique(GR.chunkat(src, s) for (_, s) in keys(t5_entries(dst, tile, src))))

@testset "Interpolation weights" begin

    @testset "source sampling is independent of output sampling" begin
        @test GR.sourcesampling(Conservative()) isa DD.Lookups.Intervals
        @test GR.sourcesampling(UnimplementedMethod()) isa DD.Lookups.Intervals
        for method in (NearestCell(), DirectNearest(), BarycentricPoint())
            @test GR.sourcesampling(method) === DD.Lookups.Points()
        end
        # Methods may override either trait independently.
        @test GR.sourcesampling !== GR.outputsampling

        # A source with one presentation resolves identically for every method.
        space = ToyLonLatSpace(8, 4)
        for method in (Conservative(), NearestCell(), DirectNearest(),
                       BarycentricPoint(), UnimplementedMethod())
            @test GR.sourcespacefor(space, method) === space
        end
        @test GR.sourcespacefor(space, NearestCell()) ===
              GR._asspace(space, "from")

        plan = plan_regrid(collect(reshape(1.0:32.0, 8, 4)); to = ToyLonLatSpace(4, 2),
            from = space, method = NearestCell(), lazy = false)
        @test plan.src_space === space
    end

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

    @testset "chart Q1 stencils" begin
        # Source centres are at ±135°/±45° longitude and ±45° latitude.
        src = ToyLonLatSpace(4, 2)
        src_inds = ownedindices(src, 1)

        # Asymmetric fractional position verifies axis order and weights.
        dst = ToyLonLatSpace(1, 1; lon = (-27.5, -17.5), lat = (17.5, 27.5))
        entries = t4_entries(t4_build(BarycentricPoint(), dst, [1], src, src_inds),
            [1], src_inds)
        @test length(entries) == 4
        @test entries[(1, localindex(src, 2, 1))] ≈ 0.1875
        @test entries[(1, localindex(src, 3, 1))] ≈ 0.0625
        @test entries[(1, localindex(src, 2, 2))] ≈ 0.5625
        @test entries[(1, localindex(src, 3, 2))] ≈ 0.1875
        @test t4_rowsum(entries, 1) ≈ 1.0

        # A centroid at 180° interpolates across the periodic seam.
        seam = ToyLonLatSpace(1, 1; lon = (175.0, 185.0), lat = (-5.0, 5.0))
        entries = t4_entries(t4_build(BarycentricPoint(), seam, [1], src, src_inds),
            [1], src_inds)
        @test length(entries) == 4
        @test all(entries[(1, localindex(src, ix, iy))] ≈ 0.25
                  for ix in (1, 4), iy in (1, 2))

        # Past the outermost latitude centres there is no cell to interpolate
        # in, on a periodic longitude axis as anywhere else: the stencil is
        # empty and the missing policy decides the destination.
        polar = ToyLonLatSpace(1, 1; lon = (-5.0, 5.0), lat = (75.0, 85.0))
        @test isempty(t4_entries(t4_build(BarycentricPoint(), polar, [1], src, src_inds),
            [1], src_inds))

        # And likewise past a non-periodic corner, where a clamp would have had
        # the corner cell answer for ground it does not reach.
        patch = ToyLonLatSpace(4, 2; lon = (-40.0, 40.0), lat = (-20.0, 20.0))
        patch_inds = ownedindices(patch, 1)
        corner = ToyLonLatSpace(1, 1; lon = (38.0, 40.0), lat = (18.0, 20.0))
        @test isempty(t4_entries(
            t4_build(BarycentricPoint(), corner, [1], patch, patch_inds),
            [1], patch_inds))

        # A point outside the source's coverage takes no weights either way:
        # `NearestCell` asks `cellat`, the chart point brackets its axes, and
        # neither invents a source for it.
        outside = ToyLonLatSpace(1, 1; lon = (-1.0, 1.0), lat = (59.0, 61.0))
        @test cellat(patch, cellcentroid(outside, 1)) === nothing
        @test isempty(t4_entries(
            t4_build(BarycentricPoint(), outside, [1], patch, patch_inds),
            [1], patch_inds))
        @test isempty(t4_entries(
            t4_build(NearestCell(), outside, [1], patch, patch_inds),
            [1], patch_inds))

        # Chart interpolation reproduces a linear chart field.
        fsrc = ToyLonLatSpace(36, 18)
        fdst = ToyLonLatSpace(17, 8; lon = (-170.0, 170.0), lat = (-80.0, 80.0))
        fdst_inds, fsrc_inds = ownedindices(fdst, 1), ownedindices(fsrc, 1)
        function linear(p)
            lon, lat = toy_lonlat(p)
            return 2.0 + 0.01 * lon + 0.03 * lat
        end
        field = [linear(cellcentroid(fsrc, i)) for i in fsrc_inds]
        weights = WeightBlock(
            t4_build(BarycentricPoint(), fdst, fdst_inds, fsrc, fsrc_inds),
            length(fdst_inds), length(fsrc_inds)).weights
        @test weights * field ≈ [linear(cellcentroid(fdst, i)) for i in fdst_inds]
        @test all(≈(1.0), sum(weights; dims = 2))

        # Chart interpolation requires a complete chart interface.
        @test_throws "hascellchart" buildweights!(
            WeightCOO(1), BarycentricPoint(), src, [1], T4NoChartSpace(), [1])
        @test_throws "chartaxes" buildweights!(
            WeightCOO(1), BarycentricPoint(), src, [1], T4BareChartSpace(), [1])
    end

    @testset "stencils partition across source chunks" begin
        # Every stencil crosses the east/west source-chunk boundary.
        src = ToyLonLatSpace(4, 2; chunks = (2, 2))
        dst = ToyLonLatSpace(3, 2; lon = (-30.0, 30.0), lat = (-30.0, 30.0))
        dst_inds = ownedindices(dst, 1)
        @test nchunks(src) == 2

        whole_inds = ownedindices(ToyLonLatSpace(4, 2), 1)
        whole = t4_entries(t4_build(BarycentricPoint(), dst, dst_inds, src, whole_inds),
            dst_inds, whole_inds)

        blocks = [t4_entries(
                      t4_build(BarycentricPoint(), dst, dst_inds, src,
                          ownedindices(src, c)),
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
        # Chart support reaches beyond geometric overlap.
        space = ToyLonLatSpace(8, 4)
        @test supportradius(NearestCell(), space) == 0.0
        @test GR.chartradius(space) > 0
        @test GR.chartradius(space) ≈ deg2rad(45.0)
        @test supportradius(BarycentricPoint(), space) == GR.chartradius(space)

        # The larger axis spacing bounds both directions.
        oblong = ToyLonLatSpace(36, 6)
        @test GR.chartradius(oblong) >= deg2rad(dlon(oblong))
        @test GR.chartradius(oblong) >= deg2rad(dlat(oblong))

        # A source with no chart has no chart spacing to bound.
        @test_throws "hascellchart" GR.chartradius(T4NoChartSpace())

        # Nearest support is zero on a chart source too: the stencil is the
        # cell the destination point is already inside, so it does not follow
        # the source's spacing anywhere.
        @test supportradius(NearestCell(), oblong) == 0.0
    end

    @testset "a sampling may specialise the weight block seam" begin
        src = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (4, 2))
        dst = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (8, 2))
        field = collect(reshape(1.0:32.0, 8, 4))
        ndst, nsrc = Int(ncells(dst)), Int(ncells(src))

        # A method that specialises no seam of its own reaches `buildweights!`
        # under either sampling, once per block, and keeps every entry it
        # emits.
        point = T4SplitMethod(DD.Lookups.Points())
        area = T4SplitMethod(DD.Lookups.Intervals(DD.Lookups.Center()))
        pointweights = GR.wholeblock(point, dst, src).weights
        @test Matrix(pointweights) == Matrix(LinearAlgebra.I, ndst, nsrc)
        @test point.builds[] == 1
        @test pointweights == GR.wholeblock(area, dst, src).weights

        # Eager and chunked agree for it, so the seam moved no value.
        eager = regrid(field; to = dst, from = src, method = point, lazy = false)
        dest = Vector{Float64}(undef, ndst)
        regrid!(dest, field,
            ChunkedPlan(point, Weighted(0.5), dst, src; storage = PerChunk()))
        @test dest == eager

        # One method type on two samplings. Only the `Points()` one is
        # answered by the fixture's own `weightblock` method, in the eager
        # builder and in the chunk-pair builder alike; the other reaches the
        # generic `weightblock` and so `buildweights!`, which the fixture does
        # not implement.
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

        # The shipped point methods answer the same eagerly, by chunk, lazily
        # and against a differently chunked source, entry for entry rather than
        # approximately. `odst` reaches past the source on both axes, so the
        # cells outside coverage take no weight and the missing policy decides
        # them on every route.
        rows = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (8, 1))
        pdst = ToyLonLatSpace(5, 3; lon = (-30.0, 30.0), lat = (-15.0, 15.0),
            chunks = (5, 2))
        odst = ToyLonLatSpace(6, 3; lon = (-60.0, 60.0), lat = (-33.0, 33.0),
            chunks = (3, 2))
        @test count(i -> cellat(src, cellcentroid(odst, i)) === nothing,
            1:Int(ncells(odst))) > 0
        # A stencil of several sites is summed in as many partial sums as the
        # chunking gives it blocks, so a second chunking may reassociate it by
        # an ulp; `NearestCell` takes one source value whole and cannot. A
        # destination the source cannot map is blanked either way, and a blank
        # compares by where it is rather than by value.
        for (method, same) in ((NearestCell(), isequal), (BarycentricPoint(), isapprox)),
            tdst in (pdst, odst)

            nd = Int(ncells(tdst))
            reference = regrid(field; to = tdst, from = src, method, lazy = false)
            out = Vector{Float64}(undef, nd)
            regrid!(out, field, ChunkedPlan(method, Weighted(0.5), tdst, src;
                storage = PerChunk()))
            @test all(isequal(out[i], reference[i]) for i in 1:nd)

            lazy = LazyRegridArray(field, ChunkedPlan(method, Weighted(0.5), tdst,
                src; storage = PerChunk()))[1:nd]
            @test all(isequal(lazy[i], reference[i]) for i in 1:nd)

            rechunked = Vector{Float64}(undef, nd)
            regrid!(rechunked, field, ChunkedPlan(method, Weighted(0.5), tdst, rows;
                storage = PerChunk()))
            @test all(isnan(rechunked[i]) == isnan(reference[i]) for i in 1:nd)
            mapped = findall(i -> !isnan(reference[i]), 1:nd)
            @test same([rechunked[i] for i in mapped], [reference[i] for i in mapped])
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
    @testset "a point tile is one build with an exact manifest" begin
        # The counting fixture above, for a method that supplies a sampler: one
        # destination tile against an 18-chunk source, with stencils that
        # bracket each destination across source-chunk seams.
        src = ToyLonLatSpace(36, 18; chunks = (6, 6))
        dst = ToyLonLatSpace(20, 10; lon = (-19.0, 21.0), lat = (-20.0, 20.0))
        tile = ownedindices(dst, 1)
        method = T5PlaceCount()
        plan = ChunkedPlan(method, Weighted(0.5), dst, src; storage = PerChunk())
        candidates = GR.sourcesof(GR.dependencies(plan), 1)
        k = length(candidates)
        @test 1 < k < nchunks(src)

        # Driving every pair the relation names locates each destination once,
        # where the pair route above locates it once per candidate.
        blocks = [GR.blockfor(plan, (1, Int(s)), tile).block for s in candidates]
        @test method.placed[] == length(tile)
        @test method.placed[] < length(tile) * k

        # What that one pass produces is one two-entry stencil per destination,
        # partitioned across fewer chunks than were searched.
        nonzeros = sum(count(!iszero, b.weights) for b in blocks)
        contributing = count(b -> any(!iszero, b.weights), blocks)
        @test nonzeros == 2 * length(tile)
        @test 0 < contributing < k

        # The tile is the cache entry, and no chunk pair is one.
        @test GR.nblocks(plan.storage) == 1
        @test isempty(plan.storage.blocks)

        # The manifest is exactly the chunks owning a nonzero entry: ascending,
        # strictly increasing, inside the relation's row, and shorter than it.
        weights = t5_weights(method, dst, tile, src)
        @test weights.sourcechunks == t5_owners(dst, tile, src)
        @test issorted(weights.sourcechunks) && allunique(weights.sourcechunks)
        @test length(weights.blocks) == length(weights.sourcechunks)
        @test issubset(weights.sourcechunks, candidates)
        @test length(weights.sourcechunks) < k
        @test all(any(!iszero, b.weights) for b in weights.blocks)

        # Every entry sits at its chunk-local column of its owner's block, and
        # some stencil straddles a seam, so the split is a real one.
        expected = t5_entries(dst, tile, src)
        got = t5_blockentries(weights, src)
        @test length(got) == length(expected)
        @test all(got[key] ≈ expected[key] for key in keys(expected))
        @test count(1:length(tile)) do j
            (a, _), (b, _) = t5_bracket(src, cellcentroid(dst, tile[j]))
            GR.chunkat(src, a) != GR.chunkat(src, b)
        end > 0
    end

    @testset "a nearest tile is one build with an exact manifest" begin
        # The counting fixture the two routes above share, for `NearestCell`:
        # one destination tile against a 72-chunk source, straddling source
        # chunk seams on both axes. The relation is cap overlap at radius zero,
        # so it names more chunks than the tile's cells sit in.
        src = ToyLonLatSpace(36, 18; chunks = (3, 3))
        dst = ToyLonLatSpace(20, 10; lon = (-19.0, 21.0), lat = (-20.0, 20.0))
        @test nchunks(dst) == 1
        tile = ownedindices(dst, 1)

        # The chunk-pair route, which `buildblock` still takes, locates every
        # destination of the tile once per candidate source chunk.
        paired = T6LocateCount(src)
        pairplan = ChunkedPlan(NearestCell(), Weighted(0.5), dst, paired;
            storage = PerChunk())
        candidates = GR.sourcesof(GR.dependencies(pairplan), 1)
        k = length(candidates)
        @test 1 < k < nchunks(src)
        for s in candidates
            GR.buildblock(pairplan, tile, ownedindices(src, Int(s)))
        end
        @test paired.located[] == length(tile) * k

        # A chunked plan takes the tile route, and it locates each destination
        # once however many pairs ask for it.
        counting = T6LocateCount(src)
        plan = ChunkedPlan(NearestCell(), Weighted(0.5), dst, counting;
            storage = PerChunk())
        @test GR.tilesampler(plan) !== nothing
        blocks = [GR.blockfor(plan, (1, Int(s)), tile).block for s in candidates]
        @test counting.located[] == length(tile)
        @test counting.located[] < length(tile) * k

        # One entry of weight exactly one per destination, on fewer chunks than
        # were searched.
        nonzeros = sum(count(!iszero, b.weights) for b in blocks)
        contributing = count(b -> any(!iszero, b.weights), blocks)
        @test nonzeros == length(tile)
        @test all(sum(b.weights) == count(!iszero, b.weights) for b in blocks)
        @test 0 < contributing < k

        # The tile is the cache entry, and no chunk pair is one.
        @test GR.nblocks(plan.storage) == 1
        @test isempty(plan.storage.blocks)

        # The manifest is exactly the chunks owning those cells, and the tile
        # straddles a source chunk seam, so the split is a real one.
        weights = t5_weights(NearestCell(), dst, tile, src)
        @test weights.sourcechunks == t6_owners(dst, tile, src)
        @test length(weights.sourcechunks) > 1
        @test issubset(weights.sourcechunks, candidates)
        @test length(weights.sourcechunks) < k
    end

    @testset "point tile values match eager, lazy and another chunking" begin
        # One lattice under two chunkings: patchwork chunks, whose chunk-local
        # index is a lookup, and full-width ones, whose chunk-local index is
        # arithmetic.
        src = ToyLonLatSpace(36, 18; chunks = (6, 6))
        rows = ToyLonLatSpace(36, 18; chunks = (36, 6))
        dst = ToyLonLatSpace(20, 10; lon = (-19.0, 21.0), lat = (-20.0, 20.0))
        tile = ownedindices(dst, 1)
        field = collect(reshape(1.0:648.0, 36, 18))
        ndst = Int(ncells(dst))
        @test !(ownedindices(src, 1) isa AbstractUnitRange) &&
              ownedindices(rows, 1) isa AbstractUnitRange

        # The eager whole domain locates each destination once too, and gives
        # the bracket's own interpolation.
        method = T5PlaceCount()
        eager = regrid(field; to = dst, from = src, method, lazy = false)
        @test method.placed[] == ndst
        @test eager ≈ [sum(w * field[i] for (i, w) in t5_bracket(src, cellcentroid(dst, j)))
                       for j in 1:ndst]

        chunked = Vector{Float64}(undef, ndst)
        regrid!(chunked, field, ChunkedPlan(T5PlaceCount(), Weighted(0.5), dst, src;
            storage = PerChunk()))
        @test chunked ≈ eager

        lazy = LazyRegridArray(field, ChunkedPlan(T5PlaceCount(), Weighted(0.5), dst,
            src; storage = PerChunk()))[1:ndst]
        @test lazy == chunked

        rechunked = Vector{Float64}(undef, ndst)
        regrid!(rechunked, field, ChunkedPlan(T5PlaceCount(), Weighted(0.5), dst, rows;
            storage = PerChunk()))
        @test rechunked ≈ eager

        # The stencils do not move with the chunking; the manifest does.
        @test t5_weights(method, dst, tile, rows).sourcechunks == t5_owners(dst, tile, rows)
        @test t5_owners(dst, tile, rows) != t5_owners(dst, tile, src)
    end

    @testset "a stencil naming one source cell twice keeps one entry" begin
        # A one-column source: the bracket wraps onto the same cell twice.
        src = ToyLonLatSpace(1, 4)
        dst = ToyLonLatSpace(1, 1; lon = (80.0, 100.0), lat = (10.0, 20.0))
        (a, wa), (b, wb) = t5_bracket(src, cellcentroid(dst, 1))
        @test a == b && 0 < wa < 1 && wa + wb ≈ 1.0

        plan = ChunkedPlan(T5PlaceCount(), Weighted(0.5), dst, src; storage = PerChunk())
        block = GR.blockfor(plan, (1, GR.chunkat(src, a)), ownedindices(dst, 1)).block
        @test count(!iszero, block.weights) == 1
        @test sum(block.weights) ≈ 1.0
    end

    @testset "a point tile builds once, evicts and spills whole" begin
        src = ToyLonLatSpace(36, 18; chunks = (6, 6))
        dst = ToyLonLatSpace(20, 10; lon = (-19.0, 21.0), lat = (-20.0, 20.0),
            chunks = (20, 5))
        @test nchunks(dst) == 2
        first_tile, second_tile = ownedindices(dst, 1), ownedindices(dst, 2)

        # Two tasks asking for two source chunks of one tile: the second meets
        # the build in flight and waits for it instead of starting a second.
        method = T5PlaceCount(; yielding = true)
        plan = ChunkedPlan(method, Weighted(0.5), dst, src; storage = PerChunk())
        rows = (GR.sourcesof(GR.dependencies(plan), 1),
            GR.sourcesof(GR.dependencies(plan), 2))
        tasks = [Threads.@spawn GR.blockfor(plan, (1, Int(s)), first_tile)
                 for s in rows[1][1:2]]
        foreach(fetch, tasks)
        @test method.placed[] == length(first_tile)
        @test GR.nblocks(plan.storage) == 1
        manifest = copy(plan.storage.tiles[1].sourcechunks)

        # A budget too small for two tiles evicts the first, and the next
        # request rebuilds it once, the same tile it was.
        s1, s2 = Int(first(rows[1])), Int(first(rows[2]))
        tight = T5PlaceCount()
        tightplan = ChunkedPlan(tight, Weighted(0.5), dst, src;
            storage = PerChunk(; maxbytes = 1))
        once = Matrix(GR.blockfor(tightplan, (1, s1), first_tile).block.weights)
        @test tight.placed[] == length(first_tile)
        GR.blockfor(tightplan, (2, s2), second_tile)
        @test GR.nblocks(tightplan.storage) == 1
        again = Matrix(GR.blockfor(tightplan, (1, s1), first_tile).block.weights)
        @test tight.placed[] == 2 * length(first_tile) + length(second_tile)
        @test again == once

        # A spilled tile comes back off disk with its manifest, and no
        # destination pass runs to recover either.
        spill = T5PlaceCount()
        storage = Spilled(mktempdir(); maxbytes = 1)
        spillplan = ChunkedPlan(spill, Weighted(0.5), dst, src; storage)
        before = Matrix(GR.blockfor(spillplan, (1, s1), first_tile).block.weights)
        @test spill.placed[] == length(first_tile)
        @test length(GR.spilledfiles(storage)) == 1
        GR.blockfor(spillplan, (2, s2), second_tile)
        recovered = Matrix(GR.blockfor(spillplan, (1, s1), first_tile).block.weights)
        @test spill.placed[] == length(first_tile) + length(second_tile)
        @test recovered == before
        @test storage.memory.tiles[1].sourcechunks == manifest
    end

    @testset "the conservative build unit is still one chunk pair" begin
        src = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (4, 2))
        cdst = ToyLonLatSpace(4, 2; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (4, 1))
        field = collect(reshape(1.0:32.0, 8, 4))
        plan = ChunkedPlan(Conservative(), Weighted(0.5), cdst, src; storage = PerChunk())
        out = Vector{Float64}(undef, Int(ncells(cdst)))
        regrid!(out, field, plan)
        @test out ≈ regrid(field; to = cdst, from = src, method = Conservative(),
            lazy = false)

        # An area method has no sampler, so it keeps the pair key and holds no
        # tile.
        @test GR.tilesampler(plan) === nothing
        @test isempty(plan.storage.tiles)
        @test GR.nblocks(plan.storage) == length(plan.storage.blocks) > 1
        @test all(key isa Tuple{Int,Int} for key in keys(plan.storage.blocks))
    end
end
