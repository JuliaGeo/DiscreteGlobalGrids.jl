# Regrid the Copernicus DEM onto IGeo7 and write it to a Zarr store.
#
#     one machine   julia --project=scripts/copdem_regrid -p 4 -t 2 scripts/copdem_regrid/regrid.jl
#     Slurm         sbatch scripts/copdem_regrid/regrid.sbatch
#
# The globe is cut into IGeo7 cells at `chunklevel`. Each one is a unit of work
# and one Zarr chunk, so workers never write to the same file. Dagger hands
# batches of chunks to the workers. Running the script again resumes the store.

using Distributed

# Inside a Slurm allocation, start one worker per Slurm task.
if haskey(ENV, "SLURM_NTASKS") && nprocs() == 1
    using SlurmClusterManager
    addprocs(SlurmManager(); exeflags = "--threads=$(get(ENV, "SLURM_CPUS_PER_TASK", "1"))")
end

@everywhere begin
    import DiscreteGlobalGrids as DGG
    import GlobalRegridding as GR
    import Dagger
    import Zarr                                     # loads the DGG Zarr store
    using CopernicusUtils
    using Base.ScopedValues: with
end

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

res = parse(Int, get(ENV, "COPDEM_RES", "30"))      # 30 for GLO-30, 90 for GLO-90
level = res == 30 ? 13 : 12                         # first IGeo7 level finer than the DEM's posts
data = get(ENV, "COPDEM_DATA", joinpath(homedir(), "copdem"))
synthetic = get(ENV, "COPDEM_SYNTHETIC", "0") != "0"

settings = (;
    res, level,
    synthetic,                                      # a trial run: every tile is fabricated and no DEM is read
    chunklevel = level - 7,                         # 7^7 = 823 543 cells per chunk
    method     = DGG.Conservative(),                # or BarycentricPoint(), NearestCell()
    policy     = DGG.Weighted(0.5),                 # a cell under half covered by data is NaN; Weighted(1) for BarycentricPoint()
    tiledir    = get(ENV, "COPDEM_TILES", joinpath(data, "glo$res")), # see "Where the tiles are"
    # true: fetch missing tiles from AWS to the path `locator` gives. false: the
    # tiles on disk are complete, and a tile with no file is ocean.
    download   = get(ENV, "COPDEM_DOWNLOAD", "1") != "0",
    store      = get(ENV, "COPDEM_STORE",
                     joinpath(data, "glo$res$(synthetic ? "-synthetic" : "")-igeo7-l$level.zarr")),
    # `nothing` for the globe, or boxes of whole degrees, [(west, east, south, north)].
    # COPDEM_REGION=5,8,45,48 gives one box.
    region     = haskey(ENV, "COPDEM_REGION") ?
                 [Tuple(parse.(Float64, split(ENV["COPDEM_REGION"], ',')))] : nothing,
    batch      = parse(Int, get(ENV, "COPDEM_BATCH", "8")), # chunks per Dagger task
    cachetiles = res == 30 ? 32 : 256,              # decoded tiles each worker keeps: 52 MB each at 30 m, 6 MB at 90 m
)

# ---------------------------------------------------------------------------
# Where the tiles are
# ---------------------------------------------------------------------------

# `locator(s)` returns a function from a tile's south-west corner, in whole
# degrees, to the path of its GeoTIFF, or to `nothing` for a tile that does not
# exist. Replace it to read any layout.
#
# The default reads AWS file names under `tiledir`, flat
# (`<tiledir>/<stem>.tif`) or laid out like the bucket
# (`<tiledir>/<stem>/<stem>.tif`), where `<stem>` is
# `Copernicus_DSM_COG_10_N45_00_E006_00_DEM`.
@everywhere locator(s) = TileDirectory(s.tiledir, s.res)

# A custom layout: one directory per hemisphere and latitude band, such as
# `/mnt/dem/north/45/<stem>.tif`.
#
# @everywhere locator(s) = (lat, lon) -> joinpath("/mnt/dem", lat < 0 ? "south" : "north",
#     string(abs(lat)), tilestem(s.res, lat, lon) * ".tif")

# ---------------------------------------------------------------------------
# What each worker does
# ---------------------------------------------------------------------------

# The destination geometry. The store is indexed by plain IGeo7 cells.
@everywhere igeo7() = DGG.AuthalicSystem(DGG.IGeo7System())

# A worker's view of the run: the DEM as one lazy array, its grid, and the store.
@everywhere function openworker(s, tiles)
    # glibc otherwise keeps freed tile buffers: resident memory settles near three times the live heap.
    Sys.islinux() && ccall(:mallopt, Cint, (Cint, Cint), -1, 32 * 2^20)
    sys = DGG.CopernicusDEMSystem(s.res)
    source = s.synthetic ? SyntheticTiles(sys) :
        CopernicusTiles(sys, tiles; locate = locator(s), download = s.download)
    dem = TiledDEM(source, tiles; slots = s.cachetiles)
    from = DGG.DGGSpace(DGG.PartialGrid(sys, 1, dem.ids); chunklevel = 0)
    return (; dem, from, store = DGG.subzonestore(s.store))
end

# Regrid and write `chunks`; returns them with the id of the worker that ran them.
@everywhere function regridbatch(worker, chunks, s)
    # Dagger runs one batch per thread, so each regrid stays on its own thread.
    with(GR.OUTER_PARALLEL => true) do
        for chunk in chunks
            cell = DGG.columncell(worker.store.layout, chunk)
            to = DGG.DGGSpace(DGG.subtree(igeo7(), cell, s.level); chunklevel = s.chunklevel)
            elevation = GR.regrid(worker.dem; to, from = worker.from, method = s.method,
                missingpolicy = s.policy, lazy = true, budget = 2^30)
            DGG.dggwrite!(worker.store, chunk, Float32.(vec(collect(elevation))))
        end
    end
    return (; chunks, worker = myid())
end

# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------

sys = DGG.CopernicusDEMSystem(res)
tiles = settings.download || synthetic ?
    landtiles(sys, tilelist(data, res), settings.region) :
    landtiles(sys, locator(settings), settings.region)
isempty(tiles) && error("no tiles found for region $(settings.region)")

chunks = covering_chunks(igeo7(), sys, tiles, settings.chunklevel; nthreads = Threads.nthreads())

isdir(settings.store) || DGG.subzonestore(settings.store, DGG.IGeo7System(), level;
    ancestor_level = settings.chunklevel, capacity = 7^(level - settings.chunklevel),
    layers = ("elevation" => Float32,), fill_value = NaN, ancestor_coordinate = true,
    attrs = Dict{String,Any}("title" => "Copernicus DEM GLO-$res$(synthetic ? " (synthetic)" : "") on IGeo7 level $level",
        "regridding_method" => string(settings.method)))

# The ledger of finished chunks, which is what a rerun skips.
ledger = settings.store * ".done.txt"
done = isfile(ledger) ? Set(parse.(Int, readlines(ledger))) : Set{Int}()
todo = filter(!in(done), chunks)
println("$(length(tiles)) tiles, $(length(chunks)) chunks, $(length(done)) already done, ",
    "$(nworkers()) workers")

# Ascending chunk indices are neighbours on the sphere, so a batch reads the same few tiles.
batches = collect(Iterators.partition(todo, settings.batch))
state = Dagger.shard(() -> openworker(settings, tiles); workers = workers())
tasks = [Dagger.@spawn scope = Dagger.scope(; workers = workers()) regridbatch(state, collect(b), settings)
         for b in batches]

started = time()
failed = 0
open(ledger, "a") do io
    for (i, task) in enumerate(tasks)
        try
            result = fetch(task)
            foreach(c -> println(io, c), result.chunks)
            flush(io)
            println("batch $i/$(length(tasks)) on worker $(result.worker), ",
                round(Int, time() - started), " s")
        catch err
            global failed += 1
            @error "batch $i failed; rerun to retry chunks $(batches[i])" exception = err
        end
    end
end

failed == 0 || error("$failed of $(length(tasks)) batches failed; run again to retry them")
println("done: $(settings.store)")
