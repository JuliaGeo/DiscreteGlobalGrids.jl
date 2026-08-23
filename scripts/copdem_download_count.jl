# Count the distinct real GLO-90 source tiles demanded by the production chunk
# graph, then estimate transfer size from 18 stratified HTTP HEAD requests.
# Extents and headers only: this script downloads no GeoTIFF bodies.
#
#   RASTERDATASOURCES_PATH=/path/to/bench/data \
#   COPDEM_COLUMNS=/path/to/copdem90-igeo7-l12.columns.txt \
#     nice -n 10 julia --project=benchmark --threads=4 scripts/copdem_download_count.jl

import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import Graphs
import Downloads
using Printf

const CDCOUNT = DGG.CopernicusDEM
const BASEURLCOUNT = "https://copernicus-dem-90m.s3.amazonaws.com"

# Lazy concatenation of complete-grid position ranges, matching the production
# TileIds representation without materializing tens of billions of cell ids.
struct CountConcatIds{G<:DGG.AbstractGrid,ID} <: AbstractVector{ID}
    grid::G
    starts::Vector{Int}
    offsets::Vector{Int}
end

function CountConcatIds(grid::DGG.AbstractGrid, ranges::Vector{UnitRange{Int}})
    starts = [first(r) for r in ranges]
    offsets = Vector{Int}(undef, length(ranges) + 1)
    offsets[1] = 0
    for (k, r) in enumerate(ranges)
        offsets[k + 1] = offsets[k] + length(r)
    end
    ID = typeof(DGG.cellindex(grid, first(starts)))
    return CountConcatIds{typeof(grid),ID}(grid, starts, offsets)
end

Base.size(v::CountConcatIds) = (v.offsets[end],)
Base.IndexStyle(::Type{<:CountConcatIds}) = Base.IndexLinear()
Base.@propagate_inbounds function Base.getindex(v::CountConcatIds, i::Int)
    @boundscheck checkbounds(v, i)
    b = searchsortedlast(v.offsets, i - 1)
    return DGG.cellindex(v.grid, v.starts[b] + (i - 1 - v.offsets[b]))
end
DGG.Helpers.strictly_increasing(::CountConcatIds) = true

function parsestem(sys, stem)
    m = match(r"_([NS])(\d{2})_00_([EW])(\d{3})_00_DEM$", stem)
    m === nothing && error("unparseable Copernicus stem $stem")
    lat = parse(Int, m[2]) * (m[1] == "S" ? -1 : 1)
    lon = parse(Int, m[4]) * (m[3] == "W" ? -1 : 1)
    return CDCOUNT.tilecell(sys, lat, lon)
end

function tileurl(stem)
    return string(BASEURLCOUNT, "/", stem, "/", stem, ".tif")
end

# Same conservative lon/lat-box narrow phase used in the merged DAG proof.
const COUNTPAD = 0.01
function count_lonhalfwidth(latdeg::Float64, rdeg::Float64)
    abs(latdeg) + rdeg >= 90 - COUNTPAD && return 180.0
    s = sind(rdeg) / cosd(latdeg)
    s >= 1 && return 180.0
    return rad2deg(asin(s))
end
count_circulardlon(a::Float64, b::Float64) =
    abs(mod(a - b + 180.0, 360.0) - 180.0)
function count_boxesoverlap(dlat::Float64, dlon::Float64, dlathalf::Float64,
        dlonhalf::Float64, tlat::Float64, tlon::Float64)
    abs(dlat - tlat) <= dlathalf + 0.5 + COUNTPAD || return false
    dlonhalf >= 180 - COUNTPAD && return true
    return count_circulardlon(dlon, tlon) <= dlonhalf + 0.5 + COUNTPAD
end

bandindex(lat::Real) = abs(lat) < 30 ? 1 : abs(lat) < 60 ? 2 : 3
const BANDLABELS = ("|lat| 0-30", "|lat| 30-60", "|lat| 60-90")

function selectsample(indices, latitudes, longitudes)
    targets = [
        (-25.0, -150.0), (-15.0, -90.0), (-5.0, -30.0),
        (5.0, 30.0), (15.0, 90.0), (25.0, 150.0),
        (-55.0, -150.0), (-45.0, -90.0), (-35.0, -30.0),
        (35.0, 30.0), (45.0, 90.0), (55.0, 150.0),
        (-85.0, -150.0), (-75.0, -90.0), (-65.0, -30.0),
        (65.0, 30.0), (75.0, 90.0), (85.0, 150.0),
    ]
    remaining = Set(indices)
    picked = Int[]
    for (targetlat, targetlon) in targets
        b = bandindex(targetlat)
        candidates = (i for i in remaining if bandindex(latitudes[i]) == b)
        best = argmin(i -> 20 * abs(latitudes[i] - targetlat) +
                           count_circulardlon(longitudes[i], targetlon), candidates)
        push!(picked, best)
        delete!(remaining, best)
    end
    return picked
end

function headlength(stem)
    response = Downloads.request(tileurl(stem); method = "HEAD", timeout = 60.0)
    response.status == 200 || error("HEAD for listed tile $stem returned HTTP $(response.status)")
    entry = findfirst(p -> lowercase(first(p)) == "content-length", response.headers)
    entry === nothing && error("HEAD for $stem returned no Content-Length")
    return parse(Int, last(response.headers[entry]))
end

function countreport(label, graph, latitudes)
    used = [s for s in 1:GR.nsourcechunks(graph) if GR.consumerdegree(graph, s) > 0]
    bands = [count(s -> bandindex(latitudes[s]) == b, used) for b in 1:3]
    @printf("%-10s distinct source chunks with consumers: %d\n", label, length(used))
    @printf("%-10s net listed tiles to download:          %d\n", label, length(used))
    for b in 1:3
        @printf("  %-15s %6d (%5.1f%%)\n", BANDLABELS[b], bands[b],
            100 * bands[b] / length(used))
    end
    return used, bands
end

function main_count()
    dataroot = get(ENV, "RASTERDATASOURCES_PATH",
        joinpath(@__DIR__, "..", "bench", "data"))
    listpath = joinpath(dataroot, "CopernicusDEM", "tileList-glo90.txt")
    columnpath = get(ENV, "COPDEM_COLUMNS",
        "/home/asinghvi17/geo/dggstores/copdem90-igeo7-l12-synthetic.zarr.columns.txt")

    sys = DGG.CopernicusDEMSystem(90)
    stems = String[]
    for line in eachline(listpath)
        stem = strip(line)
        isempty(stem) || push!(stems, String(stem))
    end
    tilecells = parsestem.(Ref(sys), stems)
    ordinals = [Int(c.index) for c in tilecells]
    perm = sortperm(ordinals)
    stems = stems[perm]
    ordinals = ordinals[perm]
    length(stems) == 26_475 || error("expected 26,475 listed tiles, got $(length(stems))")
    allunique(ordinals) || error("tile list contains duplicate ordinals")

    g1 = DGG.levelgrid(sys, 1)
    ranges = [let r = DGG.descendant_range(sys, DGG.LevelIndex(0, o), 1)
                  Int(first(r)):Int(last(r))
              end for o in ordinals]
    srcids = CountConcatIds(g1, ranges)
    srcspace = DGG.DGGSpace(DGG.PartialGrid(sys, 1, srcids); chunklevel = 0)

    columns = Int[]
    for line in eachline(columnpath)
        text = strip(line)
        (isempty(text) || startswith(text, "#")) && continue
        push!(columns, parse(Int, text))
    end
    issorted(columns) && allunique(columns) || error("column list is not sorted/unique")
    sys7 = DGG.IGeo7System()
    g5, g12 = DGG.levelgrid(sys7, 5), DGG.levelgrid(sys7, 12)
    dstranges = [let r = DGG.descendant_range(sys7, DGG.cellindex(g5, c), 12)
                     Int(first(r)):Int(last(r))
                 end for c in columns]
    dstids = CountConcatIds(g12, dstranges)
    dstspace = DGG.DGGSpace(DGG.PartialGrid(sys7, 12, dstids); chunklevel = 5)

    @printf("source chunks=%d destination chunks=%d threads=%d\n",
        GR.nchunks(srcspace), GR.nchunks(dstspace), Threads.nthreads())
    radius = Float64(GR.support_radius(DGG.Conservative(), srcspace))
    GR.chunk_dependency_graph(dstspace, srcspace; radius) # compile warm-up
    capseconds = @elapsed caps = GR.chunk_dependency_graph(dstspace, srcspace; radius)

    dstcaps, srccaps = GR.chunkextents(dstspace), GR.chunkextents(srcspace)
    latof(c) = rad2deg(atan(c.point[3], hypot(c.point[1], c.point[2])))
    lonof(c) = rad2deg(atan(c.point[2], c.point[1]))
    dstlat, dstlon = map(latof, dstcaps), map(lonof, dstcaps)
    dstrdeg = [rad2deg(c.radius) for c in dstcaps]
    dstlonhalf = [count_lonhalfwidth(dstlat[d], dstrdeg[d]) for d in eachindex(dstcaps)]
    tilelat, tilelon = Float64[], Float64[]
    for ordinal in ordinals
        lat, lon = CDCOUNT.tilecorner(sys, DGG.LevelIndex(0, ordinal))
        push!(tilelat, lat + 0.5)
        push!(tilelon, lon + 0.5)
    end
    refine(d, s) = count_boxesoverlap(dstlat[d], dstlon[d], dstrdeg[d],
        dstlonhalf[d], tilelat[s], tilelon[s])
    # A narrow phase is an argument to plan construction and to nothing else
    # (Task G4), so the refined relation comes from a plan that owns it. The
    # plan's radius is `support_radius(Conservative(), srcspace)`, which is the
    # same `radius` the cap graph above was built at.
    refinedgraph() = GR.dependencies(GR.ChunkedPlan(DGG.Conservative(),
        GR.Weighted(0.5), dstspace, srcspace;
        dependencies = true, refine, narrow = :copdem_tile_lonlat_box))
    refinedgraph() # compile warm-up
    refinedseconds = @elapsed refined = refinedgraph()

    @printf("caps graph:    %.3f s, %d edges\n", capseconds, Graphs.ne(caps))
    @printf("refined graph: %.3f s, %d edges\n\n", refinedseconds, Graphs.ne(refined))
    capused, capbands = countreport("caps", caps, tilelat)
    println()
    refinedused, refinedbands = countreport("refined", refined, tilelat)

    # Source chunks in this production space are exactly the listed tile set;
    # retain the explicit intersection assertion because the executor's net-I/O
    # accounting depends on that identity.
    listedordinals = Set(ordinals)
    all(ordinals[s] in listedordinals for s in capused) || error("caps used an unlisted tile")
    all(ordinals[s] in listedordinals for s in refinedused) || error("refined used an unlisted tile")

    unionused = sort!(collect(union(Set(capused), Set(refinedused))))
    samples = selectsample(unionused, tilelat, tilelon)
    lengths = Dict{Int,Int}()
    println("\nHTTP HEAD sample (decimal MB):")
    println("| band | tile | center lat | pixel columns | Content-Length (bytes) | MB |")
    println("|---|---|---:|---:|---:|---:|")
    for s in samples
        bytes = headlength(stems[s])
        lengths[s] = bytes
        lat_s, _ = CDCOUNT.tilecorner(sys, DGG.LevelIndex(0, ordinals[s]))
        ncols = CDCOUNT.ncols_at(sys, lat_s)
        @printf("| %s | `%s` | %.1f | %d | %d | %.3f |\n",
            BANDLABELS[bandindex(tilelat[s])], stems[s], tilelat[s], ncols,
            bytes, bytes / 1e6)
    end

    means = [sum(lengths[s] for s in samples if bandindex(tilelat[s]) == b) /
             count(s -> bandindex(tilelat[s]) == b, samples) for b in 1:3]
    println("\nBand model and extrapolation (decimal GB):")
    println("| band | sample n | mean bytes/tile | caps tiles | caps GB | refined tiles | refined GB |")
    println("|---|---:|---:|---:|---:|---:|---:|")
    capbytes = refinedbytes = 0.0
    for b in 1:3
        n = count(s -> bandindex(tilelat[s]) == b, samples)
        cb = capbands[b] * means[b]
        rb = refinedbands[b] * means[b]
        capbytes += cb
        refinedbytes += rb
        @printf("| %s | %d | %.0f | %d | %.3f | %d | %.3f |\n",
            BANDLABELS[b], n, means[b], capbands[b], cb / 1e9,
            refinedbands[b], rb / 1e9)
    end
    @printf("| **total** | **%d** | — | **%d** | **%.3f** | **%d** | **%.3f** |\n",
        length(samples), length(capused), capbytes / 1e9,
        length(refinedused), refinedbytes / 1e9)

    println("\nLand-mask-restricted count: skipped. The production space is built before the " *
        "15-arcsecond Natural Earth mask is applied per decoded pixel, and it carries no " *
        "column-level non-triviality metadata; deriving that subset is a separate spatial pass.")
end

main_count()
