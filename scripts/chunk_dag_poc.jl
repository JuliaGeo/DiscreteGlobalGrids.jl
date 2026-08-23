# Build the real CopDEM GLO-90 x IGeo7-L12 chunk dependency graph and report on
# it. Extents only: no pixel is read, no weight is built, no regrid is run.
#
#   RASTERDATASOURCES_PATH=.../bench/data \
#     julia --project=benchmark -t 4 --gcthreads=2 scripts/chunk_dag_poc.jl
#
# The source space is one chunk per listed 1-degree tile; the destination is one
# chunk per level-5 ancestor column of the level-12 grid. Both are built through
# the ordinary `DGGSpace` constructor, so the graph sees exactly the spaces the
# production run regrids between.

import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import Graphs
using Printf
const CD = DGG.CopernicusDEM

const COLUMNS = "/home/asinghvi17/geo/dggstores/copdem90-igeo7-l12-synthetic.zarr.columns.txt"
const OUTDIR = "/home/asinghvi17/geo/DiscreteGlobalGrids.jl/regrid-notes"
const ADJACENCY = joinpath(OUTDIR, "2026-08-21-chunk-dag-adjacency.ndjson")
const OTHER = joinpath(OUTDIR, "2026-08-21-chunk-dag-adjacency-mine.ndjson.gz")
const REFINED = joinpath(OUTDIR, "2026-08-21-chunk-dag-adjacency-refined.ndjson")
const MEASUREMENTS = joinpath(OUTDIR, "2026-08-21-chunk-dag-api.ndjson")

# A lazy concatenation of consecutive runs of a complete level grid: the ids of
# every descendant of a sorted set of ancestors, without materializing 10^10
# cell indices. This is the same trick the production run's `TileIds` uses.
struct ConcatIds{G<:DGG.AbstractGrid,ID} <: AbstractVector{ID}
    grid::G
    starts::Vector{Int}    # first position of each run in `grid`
    offsets::Vector{Int}   # cumulative run lengths, `offsets[1] == 0`
end

function ConcatIds(grid::DGG.AbstractGrid, ranges::Vector{UnitRange{Int}})
    starts = [first(r) for r in ranges]
    offsets = Vector{Int}(undef, length(ranges) + 1)
    offsets[1] = 0
    for (k, r) in enumerate(ranges)
        offsets[k+1] = offsets[k] + length(r)
    end
    ID = typeof(DGG.cellindex(grid, first(starts)))
    return ConcatIds{typeof(grid),ID}(grid, starts, offsets)
end

Base.size(v::ConcatIds) = (v.offsets[end],)
Base.IndexStyle(::Type{<:ConcatIds}) = Base.IndexLinear()
Base.@propagate_inbounds function Base.getindex(v::ConcatIds, i::Int)
    @boundscheck checkbounds(v, i)
    b = searchsortedlast(v.offsets, i - 1)
    return DGG.cellindex(v.grid, v.starts[b] + (i - 1 - v.offsets[b]))
end
# The runs are ascending and disjoint by construction, so the O(n) ascent check
# `PartialGrid` would otherwise run over 10^10 ids has nothing to find.
DGG.Helpers.strictly_increasing(::ConcatIds) = true

check(name, cond) = (cond || error("FAILED: $name"); println("  ok: $name"))

# Write one line per destination chunk. Rows come out in source-chunk order,
# which is ascending tile *ordinal*; that is not ascending tile *name* (the name
# sorts by hemisphere letter before latitude), and the artifact asks for sorted
# names, so sort each row here.
function writeadjacency(path, g, cols, tilestems)
    open(path, "w") do io
        buf = IOBuffer()
        names = String[]
        for d in 1:GR.ndestinationchunks(g)
            empty!(names)
            for s in GR.sourcesof(g, d)
                push!(names, tilestems[s])
            end
            sort!(names)
            print(buf, "{\"col\":", cols[d], ",\"tiles\":[")
            sep = false
            for nm in names
                sep && print(buf, ",")
                sep = true
                print(buf, "\"", nm, "\"")
            end
            print(buf, "]}\n")
            buf.size > (1 << 20) && write(io, take!(buf))
        end
        write(io, take!(buf))
    end
    return path
end

# A conservative narrow phase for this specific pair, to exercise the `refine`
# hook. A Copernicus tile is exactly a 1-degree lon/lat box (always 1x1 in
# extent — it is the pixel count that shrinks poleward, not the footprint), and
# a spherical cap has an exact lon/lat bounding box, so two boxes that do not
# overlap cannot intersect. This is far tighter than cap-versus-cap, where a
# circle circumscribing a square already inflates by its half-diagonal.
#
# `PAD` absorbs the half-pixel outer frame `node_extent` adds to a tile box and
# any rounding, so the test only ever errs toward keeping an edge.
const PAD = 0.01

# Half-width in longitude degrees of a cap's bounding box; 180 when it reaches
# a pole and therefore spans every meridian.
function lonhalfwidth(latdeg::Float64, rdeg::Float64)
    abs(latdeg) + rdeg >= 90 - PAD && return 180.0
    s = sind(rdeg) / cosd(latdeg)
    s >= 1 && return 180.0
    return rad2deg(asin(s))
end

# Circular longitude separation, in degrees.
circulardlon(a::Float64, b::Float64) = abs(mod(a - b + 180.0, 360.0) - 180.0)

function boxesoverlap(dlat::Float64, dlon::Float64, dlathalf::Float64,
        dlonhalf::Float64, tlat::Float64, tlon::Float64)
    abs(dlat - tlat) <= dlathalf + 0.5 + PAD || return false
    (dlonhalf >= 180 - PAD) && return true
    return circulardlon(dlon, tlon) <= dlonhalf + 0.5 + PAD
end

function main()
    nt = Threads.nthreads()
    @printf("threads = %d\n", nt)
    records = Dict{String,Any}[]

    # Source: one chunk per listed GLO-90 tile.
    println("\n== source space ==")
    sys = DGG.CopernicusDEMSystem(90)
    listpath = joinpath(ENV["RASTERDATASOURCES_PATH"], "CopernicusDEM",
        "tileList-glo90.txt")
    rx = r"_([NS])(\d{2})_00_([EW])(\d{3})_00_DEM$"
    stems = String[]
    for line in eachline(listpath)
        s = strip(line)
        isempty(s) && continue
        push!(stems, String(s))
    end
    tilecells = map(stems) do s
        m = match(rx, s)
        m === nothing && error("unparseable tile stem $s")
        lat = parse(Int, m[2]) * (m[1] == "S" ? -1 : 1)
        lon = parse(Int, m[4]) * (m[3] == "W" ? -1 : 1)
        CD.tilecell(sys, lat, lon)
    end
    ordinals = [Int(c.index) for c in tilecells]
    perm = sortperm(ordinals)
    tilestems = stems[perm]           # source chunk k <-> tilestems[k]
    tileordinals = ordinals[perm]
    check("tile ordinals are unique", allunique(tileordinals))

    g1 = DGG.levelgrid(sys, 1)
    tsrc = @elapsed begin
        srcranges = [let r = DGG.descendant_range(sys, DGG.LevelIndex(0, o), 1)
                         Int(first(r)):Int(last(r))
                     end for o in tileordinals]
        srcids = ConcatIds(g1, srcranges)
        srcgrid = DGG.PartialGrid(sys, 1, srcids)
        srcspace = DGG.DGGSpace(srcgrid; chunklevel = 0)
    end
    @printf("built in %.2f s: %s\n", tsrc, srcspace)
    check("one source chunk per listed tile",
        GR.nchunks(srcspace) == length(tilestems) == 26475)

    # Destination: one chunk per level-5 ancestor column of the level-12 grid.
    println("\n== destination space ==")
    sys7 = DGG.IGeo7System()
    cols = Int[]
    for line in eachline(COLUMNS)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        push!(cols, parse(Int, s))
    end
    check("columns file is ascending and unique",
        issorted(cols) && allunique(cols))
    g5, g12 = DGG.levelgrid(sys7, 5), DGG.levelgrid(sys7, 12)
    tdst = @elapsed begin
        dstranges = [let r = DGG.descendant_range(sys7, DGG.cellindex(g5, c), 12)
                         Int(first(r)):Int(last(r))
                     end for c in cols]
        dstids = ConcatIds(g12, dstranges)
        dstgrid = DGG.PartialGrid(sys7, 12, dstids)
        dstspace = DGG.DGGSpace(dstgrid; chunklevel = 5)
    end
    @printf("built in %.2f s: %s\n", tdst, dstspace)
    check("one destination chunk per listed column",
        GR.nchunks(dstspace) == length(cols) == 66178)
    colof = [DGG.cellposition(g5, id) for id in dstspace.chunkids]
    check("destination chunk k is column k of the columns file", colof == cols)

    # The graph. Conservative regridding has zero support radius; the relation
    # is still a superset because the caps cover the cells.
    println("\n== chunk dependency graph ==")
    radius = Float64(GR.support_radius(DGG.Conservative(), srcspace))
    @printf("support radius = %g rad\n", radius)
    GR.chunk_dependency_graph(dstspace, srcspace; radius)   # warm up
    tbuild = @elapsed graph = GR.chunk_dependency_graph(dstspace, srcspace; radius)
    @printf("built in %.2f s (%d threads): %s\n", tbuild, nt, graph)

    # Structural invariants.
    println("\n== invariants ==")
    check("bipartition sizes", GR.nsourcechunks(graph) == 26475 &&
        GR.ndestinationchunks(graph) == 66178)
    check("Graphs.nv is the vertex total",
        Graphs.nv(graph) == 26475 + 66178)
    check("every edge joins a source to a destination",
        all(e -> Graphs.src(e) <= 26475 < Graphs.dst(e), Graphs.edges(graph)))
    check("Graphs agrees the graph is bipartite", Graphs.is_bipartite(graph))
    check("the two CSR directions have equal edge counts",
        sum(GR.consumerdegree(graph, s) for s in 1:GR.nsourcechunks(graph)) ==
        Graphs.ne(graph))
    check("every destination column needs at least one tile",
        minimum(GR.sourcedegree(graph, d) for d in 1:GR.ndestinationchunks(graph)) >= 1)
    check("every tile has at least one consumer",
        minimum(GR.consumerdegree(graph, s) for s in 1:GR.nsourcechunks(graph)) >= 1)

    # Predictive completeness, on the real pair rather than on toy spaces: every
    # source chunk a lazy read can ask for has to be an edge, or `consumersof`
    # is not a refcount. This is the check the 2026-08-23 full run failed.
    tdemand = @elapsed unpredicted = let index = GR.chunkindex(srcspace),
        caps = GR.chunkextents(dstspace), buf = Int[], n = 0

        for d in 1:GR.nchunks(dstspace)
            GR.candidatechunks!(buf, index, caps[d]; radius)
            row = GR.sourcesof(graph, d)
            n += count(s -> !insorted(Int32(s), row), buf)
        end
        n
    end
    check("the graph holds every pair the source index answers ($unpredicted missing)",
        unpredicted == 0)
    @printf("demand sweep %.2f s\n", tdemand)

    dstcaps = GR.chunkextents(dstspace)
    srccaps = GR.chunkextents(srcspace)
    latof(c) = rad2deg(atan(c.point[3], hypot(c.point[1], c.point[2])))
    lonof(c) = rad2deg(atan(c.point[2], c.point[1]))
    dstlat = map(latof, dstcaps)
    srclat = map(latof, srccaps)

    # The `refine` hook, exercised on the real pair. Precompute each side's
    # bounding box once, then reject pairs whose boxes are disjoint.
    println("\n== refined graph (lon/lat box narrow phase) ==")
    dstlon = map(lonof, dstcaps)
    dstrdeg = [rad2deg(c.radius) for c in dstcaps]
    dstlonhalf = [lonhalfwidth(dstlat[d], dstrdeg[d]) for d in eachindex(dstcaps)]
    # A tile's box comes from its name, not its cap: exact by construction.
    tilelat = Float64[]
    tilelon = Float64[]
    for o in tileordinals
        lat_s, lon_w = CD.tilecorner(sys, DGG.LevelIndex(0, o))
        push!(tilelat, lat_s + 0.5)
        push!(tilelon, lon_w + 0.5)
    end
    boxrefine(d, s) = boxesoverlap(dstlat[d], dstlon[d], dstrdeg[d], dstlonhalf[d],
        tilelat[s], tilelon[s])
    # Task G4: a narrow phase is an argument to plan construction and to
    # nothing else, so the refined relation comes from a plan that owns it. The
    # plan's radius is `support_radius(Conservative(), srcspace)`, the same
    # `radius` the unrefined graph above was built at.
    refinedgraph() = GR.dependencies(GR.ChunkedPlan(DGG.Conservative(),
        GR.Weighted(0.5), dstspace, srcspace;
        dependencies = true, refine = boxrefine,
        narrow = :copdem_tile_lonlat_box))
    refinedgraph()   # warm up
    trefined = @elapsed refined = refinedgraph()
    @printf("built in %.2f s: %s\n", trefined, refined)
    check("refinement only removes edges", all(
        issubset(GR.sourcesof(refined, d), GR.sourcesof(graph, d))
        for d in 1:GR.ndestinationchunks(graph)))

    # Degree distributions, split by destination latitude.
    println("\n== degree distributions ==")
    coldeg = [GR.sourcedegree(graph, d) for d in 1:GR.ndestinationchunks(graph)]
    tiledeg = [GR.consumerdegree(graph, s) for s in 1:GR.nsourcechunks(graph)]

    function stats(name, v)
        s = sort(v)
        n = length(s)
        q(p) = s[clamp(ceil(Int, p * n), 1, n)]
        @printf("%-28s n=%-7d min=%-4d med=%-5d p90=%-5d p99=%-5d max=%d mean=%.2f\n",
            name, n, s[1], q(0.5), q(0.9), q(0.99), s[end], sum(s) / n)
        return Dict{String,Any}("group" => name, "n" => n, "min" => s[1],
            "median" => q(0.5), "p90" => q(0.9), "p99" => q(0.99), "max" => s[end],
            "mean" => sum(s) / n, "total" => sum(s))
    end

    degrecords = Dict{String,Any}[]
    push!(degrecords, stats("tiles per column (all)", coldeg))
    push!(degrecords, stats("tiles per column |lat|<=60", coldeg[abs.(dstlat).<=60]))
    push!(degrecords, stats("tiles per column |lat|>60", coldeg[abs.(dstlat).>60]))
    push!(degrecords, stats("columns per tile (all)", tiledeg))
    push!(degrecords, stats("columns per tile |lat|<=60", tiledeg[abs.(srclat).<=60]))
    push!(degrecords, stats("columns per tile |lat|>60", tiledeg[abs.(srclat).>60]))
    for r in degrecords
        r["kind"] = "degree"
    end
    append!(records, degrecords)

    ncomp = length(Graphs.connected_components(graph))
    @printf("connected components (independent work groups): %d\n", ncomp)

    # The adjacency artifact.
    println("\n== writing adjacency ==")
    # `ADJACENCY` is what the API returns with no refinement — the library's own
    # answer. `REFINED` is the same relation narrowed by the box test above, and
    # is the better input for a scheduler; both are valid supersets.
    twrite = @elapsed writeadjacency(ADJACENCY, graph, cols, tilestems)
    writeadjacency(REFINED, refined, cols, tilestems)
    @printf("wrote %s in %.2f s (%d lines, %.1f MiB)\n", ADJACENCY, twrite,
        GR.ndestinationchunks(graph), filesize(ADJACENCY) / 2^20)
    @printf("wrote %s (%.1f MiB)\n", REFINED, filesize(REFINED) / 2^20)

    # Cross-check against the independently built adjacency. That one uses the
    # run's own exact `MultiOrderCoverage` rule; this one uses covering caps, so
    # the correct relation between them is CONTAINMENT, not equality: every
    # exact edge must appear here, and the excess is this graph's conservatism.
    println("\n== cross-check ==")
    other = Dict{Int,Vector{String}}()
    rowrx = r"^\{\"col\":(\d+),\"tiles\":\[(.*)\]\}$"
    if isfile(OTHER)
        open(`zcat $OTHER`) do io
            for line in eachline(io)
                m = match(rowrx, strip(line))
                m === nothing && error("unparseable cross-check row: $line")
                names = isempty(m[2]) ? String[] :
                        [String(strip(t, '"')) for t in split(m[2], ",")]
                other[parse(Int, m[1])] = names
            end
        end
        check("both files cover the same columns",
            Set(keys(other)) == Set(cols))
        otheredges = sum(length, values(other))

        # Containment is the acceptance test: the exact relation must be a
        # subset of the conservative one. A missing edge is a correctness bug;
        # an extra edge is a wasted tile load.
        function compare(label, g)
            missingedges = extraedges = exactrows = 0
            for d in 1:GR.ndestinationchunks(g)
                theirs = Set(other[cols[d]])
                ours = Set(tilestems[s] for s in GR.sourcesof(g, d))
                missingedges += length(setdiff(theirs, ours))
                extraedges += length(setdiff(ours, theirs))
                theirs == ours && (exactrows += 1)
            end
            check("$label: no exact edge is missing", missingedges == 0)
            @printf("%-12s edges %6d  false positives %6d (%.2fx exact)  identical rows %5d/%d (%.1f%%)\n",
                label, Graphs.ne(g), extraedges, Graphs.ne(g) / otheredges,
                exactrows, length(other), 100 * exactrows / length(other))
            return Dict{String,Any}("kind" => "crosscheck", "graph" => label,
                "exact_edges" => otheredges, "edges" => Graphs.ne(g),
                "missing_edges" => missingedges, "false_positive_edges" => extraedges,
                "inflation" => Graphs.ne(g) / otheredges,
                "identical_rows" => exactrows, "rows" => length(other))
        end
        push!(records, compare("caps", graph))
        push!(records, compare("refined", refined))
    else
        println("  skipped: $OTHER not found")
    end

    push!(records, Dict{String,Any}("kind" => "build", "threads" => nt,
        "nsrcchunks" => GR.nsourcechunks(graph),
        "ndstchunks" => GR.ndestinationchunks(graph),
        "edges" => Graphs.ne(graph), "radius" => radius,
        "graph_seconds" => tbuild, "demand_sweep_seconds" => tdemand,
        "unpredicted_pairs" => unpredicted,
        "srcspace_seconds" => tsrc, "dstspace_seconds" => tdst,
        "write_seconds" => twrite, "refined_seconds" => trefined, "refined_edges" => Graphs.ne(refined), "components" => ncomp,
        "adjacency_bytes" => filesize(ADJACENCY)))

    open(MEASUREMENTS, "w") do io
        for r in records
            print(io, "{")
            print(io, join(("\"$k\":" * (v isa AbstractString ? "\"$v\"" : string(v))
                            for (k, v) in sort(collect(r); by = first)), ","))
            print(io, "}\n")
        end
    end
    println("wrote ", MEASUREMENTS)
    println("\nALL CHECKS PASSED")
end

main()
