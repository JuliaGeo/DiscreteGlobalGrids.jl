# Reproducible planning baseline for the production Copernicus DEM -> IGeo7
# workload. It downloads only the official tile-name list; no raster data or
# provider metadata is read.
#
#     julia -t 4 --project=benchmark benchmark/regridding_plan_baseline.jl
#
# Baseline on 2026-08-22, Julia 1.12.6, 4 threads: 26,475 source chunks,
# 66,175 destination chunks, 326,386 edges, 3,352,520-byte graph, 0.0589 s
# median graph construction over five samples.

include(joinpath(@__DIR__, "..", "scripts", "copdem_production.jl"))

const TILE_LIST_URL =
    "https://copernicus-dem-90m.s3.amazonaws.com/tileList.txt"
const EXPECTED_TILES = 26_475
const NSAMPLES = parse(Int, get(ENV, "DGG_REGRID_PLAN_SAMPLES", "5"))

function production_spaces()
    sys = DGG.CopernicusDEMSystem(90)
    sys7 = DGG.IGeo7System()
    tilelist = Downloads.download(TILE_LIST_URL)
    tiles = listedtiles(sys, tilelist, nothing)
    length(tiles) == EXPECTED_TILES || error(
        "official GLO-90 workload changed: expected $EXPECTED_TILES tiles, got $(length(tiles))")

    srcids = TileIds(sys, tiles)
    srcgrid = DGG.PartialGrid(sys, 1, srcids)
    srcspace = DGG.DGGSpace(srcgrid; chunklevel = 0)

    chunks = covering_chunks(sys7, sys, tiles, 5;
        nthreads = max(1, Threads.nthreads()))
    chunkgrid = DGG.levelgrid(sys7, 5)
    dstids = SubtreeIds(sys7,
        [DGG.cellindex(chunkgrid, chunk) for chunk in chunks], 12)
    timedspace = @timed DGG.DGGSpace(DGG.PartialGrid(sys7, 12, dstids);
        chunklevel = 5)
    return (; srcspace, dstspace = timedspace.value, tiles, chunks, timedspace)
end

function main()
    spaces = production_spaces()
    radius = Float64(GR.support_radius(DGG.Conservative(), spaces.srcspace))

    # Compile the graph builder against the production concrete space types
    # without using the production destination size as the warm-up.
    GR.chunk_dependency_graph(spaces.dstspace, spaces.srcspace; radius)

    times = Float64[]
    allocations = Int[]
    graph = nothing
    for _ in 1:NSAMPLES
        GC.gc()
        timed = @timed GR.chunk_dependency_graph(
            spaces.dstspace, spaces.srcspace; radius)
        push!(times, timed.time)
        push!(allocations, timed.bytes)
        graph = timed.value
    end

    edges = sum(d -> GR.sourcedegree(graph, d), 1:length(spaces.chunks))
    println((
        source_chunks = GR.nchunks(spaces.srcspace),
        destination_chunks = GR.nchunks(spaces.dstspace),
        edges,
        destination_space_seconds = spaces.timedspace.time,
        destination_space_allocated = spaces.timedspace.bytes,
        graph_seconds_min = minimum(times),
        graph_seconds_median = Statistics.median(times),
        graph_allocated_min = minimum(allocations),
        graph_summarysize = Base.summarysize(graph),
        samples = NSAMPLES,
    ))
end

main()
