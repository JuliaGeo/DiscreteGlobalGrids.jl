# Focused network/integration smoke test for the lazy Copernicus GLO-90 source.
# Downloads at most three real tiles, then regrids one IGeo7 level-5 column.
#
#   RASTERDATASOURCES_PATH=/path/to/bench/data \
#     nice -n 10 julia --project=benchmark --threads=8 scripts/copdem_lazy_smoke.jl

include("copdem_production.jl")

check_smoke(name, condition) =
    condition ? println("PASS  ", name) : error("FAILED: $name")

function smoke()
    dataroot = get(ENV, "RASTERDATASOURCES_PATH",
        joinpath(@__DIR__, "..", "bench", "data"))
    listpath = joinpath(dataroot, "CopernicusDEM", "tileList-glo90.txt")
    cache = get(ENV, "COPDEM_TILE_CACHE",
        joinpath(@__DIR__, "..", "bench", "data", "CopernicusDEM", "tiles"))
    sys = DGG.CopernicusDEMSystem(90)
    listed = listedtiles(sys, listpath, nothing)
    check_smoke("the production list has 26,475 unique tiles", length(listed) == 26_475)

    stems = [
        "Copernicus_DSM_COG_30_N00_00_E006_00_DEM",
        "Copernicus_DSM_COG_30_N60_00_E010_00_DEM",
        "Copernicus_DSM_COG_30_S90_00_E000_00_DEM",
    ]
    ordinals = [Int(stemtile(sys, stem).index) for stem in stems]
    listedset = Set(listed)
    check_smoke("all three smoke tiles are listed", all(in(listedset), ordinals))

    provider = LazyCopernicusTiles(sys, listed; cachedir = cache)

    # The ocean contract is checked before any network operation. `loadtile`
    # has to return nodata while the provider's successful-GET count stays put.
    ocean = first(t for t in 0:(DGG.ncells(sys, 0) - 1) if !(t in listedset))
    before_ocean = provider.ndownloads[]
    oceanvals = loadtile(provider, ocean)
    check_smoke("unlisted ocean tile is all NaN", all(isnan, oceanvals))
    check_smoke("unlisted ocean tile makes no network request",
        provider.ndownloads[] == before_ocean)

    # Four workers ask for the same uncached tile. The per-tile lock must turn
    # that into exactly one GET (or zero when this is a cache-resume run).
    equatorial = ordinals[1]
    wascached = isfile(tilecachepath(provider, equatorial))
    before = provider.ndownloads[]
    before_cold = provider.ncold[]
    tasks = [Threads.@spawn tilepath!(provider, equatorial) for _ in 1:4]
    paths = fetch.(tasks)
    check_smoke("concurrent requests return one cache path", all(==(paths[1]), paths))
    check_smoke("concurrent first access downloads exactly once",
        provider.ndownloads[] - before == (wascached ? 0 : 1))
    check_smoke("concurrent demand records one cold download at most",
        provider.ncold[] - before_cold == (wascached ? 0 : 1))

    println("\nDecoded real tiles:")
    for (stem, ordinal) in zip(stems, ordinals)
        vals = loadtile(provider, ordinal)
        tile = DGG.LevelIndex(0, ordinal)
        lat, _ = CD.tilecorner(sys, tile)
        dims = (Int(CD.ncols_at(sys, lat)), Int(CD.lat_intervals(sys)))
        finite = count(isfinite, vals)
        check_smoke("$stem has expected $(dims[1])x$(dims[2]) dimensions",
            length(vals) == prod(dims))
        check_smoke("$stem contains finite elevations", finite > 0)
        lo, hi = extrema(Iterators.filter(isfinite, vals))
        println("  ", stem, "  dims=", dims, " finite=", finite,
            "/", length(vals), " range=[", lo, ", ", hi, "] path=",
            tilecachepath(provider, ordinal))
        check_smoke("$stem left no trusted partial object",
            !isfile(tilecachepath(provider, ordinal) * ".part"))
    end

    # Normal regridding path: the real tile is one source chunk and the output
    # is one complete level-5 -> level-12 rooted destination column.
    tile = equatorial
    ids = TileIds(sys, [tile])
    builder = TileBuilder(sys, [tile], Dict{Int,String}(), provider, NOMASK)
    tilecache = StripedLRUCache{Vector{Float32}}(k -> buildtile(builder, k);
        slots = 4, stripes = 1)
    dem = TiledDEM(ids, builder, tilecache)
    srcgrid = DGG.PartialGrid(sys, 1, ids)
    srcspace = DGG.DGGSpace(srcgrid; chunklevel = 0)

    sys7 = DGG.IGeo7System()
    g5 = DGG.levelgrid(sys7, 5)
    columns = covering_chunks(sys7, sys, [tile], 5)
    target = (6.5, 0.5)
    distance(cpos) = let c = DGG.cellindex(g5, cpos),
                         ll = US.GeographicFromUnitSphere()(DGG.cell_centroid(g5, c))
        (ll[1] - target[1])^2 + (ll[2] - target[2])^2
    end
    column = columns[argmin(distance.(columns))]
    ancestor = DGG.cellindex(g5, column)
    dstgrid = DGG.subtree(sys7, ancestor, 12)
    dstspace = DGG.DGGSpace(dstgrid; chunklevel = 5)
    vals = Float32[]
    elapsed = @elapsed begin
        output = GR.regrid(dem; to = dstspace, from = srcspace,
            method = DGG.Conservative(), missingpolicy = DGG.Weighted(0.5),
            lazy = true, budget = 2^29)
        vals = Float32.(vec(collect(output)))
    end
    nfinite = count(isfinite, vals)
    check_smoke("one level-5 column produces 7^7 cells", length(vals) == 7^7)
    check_smoke("real-tile regrid produces non-NaN cells", nfinite > 0)
    lo, hi = extrema(Iterators.filter(isfinite, vals))
    println("\nRegrid: column_index=", column, " values=", length(vals),
        " finite=", nfinite, " NaN=", count(isnan, vals),
        " range=[", lo, ", ", hi, "] elapsed=", round(elapsed; digits = 2), " s")
    println("Successful GETs this process: ", provider.ndownloads[])
    println("ALL LAZY COPDEM SMOKE CHECKS PASSED")
end

smoke()
