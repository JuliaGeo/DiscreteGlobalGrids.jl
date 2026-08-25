# Correctness and performance gates for the chunk dependency graph (Task G1).
#
#     julia -t 8 --project=benchmark benchmark/chunk_graph_gates.jl
#
# Three questions, one harness, over one case matrix:
#
#   1. Does the relation the builder produces hold every pair that a weight
#      builder would find a nonzero weight between? The oracle builds REAL cell
#      geometry for both spaces and finds the pairs with genuine positive
#      spherical intersection area — no cap, no chunk index and no
#      `chunkextents` anywhere in its construction. The same oracle runs in the
#      test suites (`lib/GlobalRegridding/test/test_chunkgraph.jl` and
#      `test/systems/crosssystem/regrid.jl`); here it also covers the cases too
#      large to assert on every CI run.
#
#   2. What does each row builder cost? The `:indexed` arm is production's
#      `chunk_dependency_graph`. The `:latjoin` arm reimplements the
#      latitude-sorted cap join that PR #69 deleted, verbatim from
#      `ba2bbfa^:lib/GlobalRegridding/src/chunkgraph.jl`, so the comparison that
#      #69 bypassed can still be run. It is reimplemented HERE, not imported, so
#      that this harness keeps working as later tasks change the production
#      builder. `:latjoin_raw` is the same with the latitude prefilter off; it
#      exists only to attribute the prefilter's share and is opt-in behind
#      `DGG_GRAPH_GATE_RAW=1`.
#
#   3. What does a per-column plan pay for its rows? Task G3 added `restrict`,
#      which hands back a row view sharing the parent's destination-major CSR.
#      The alternative is a rebuild, and since #69 a rebuild costs a fresh
#      `chunkindex(src_space)` per column. `restrict_measurement` times both and
#      asserts they agree row for row.
#
# The two arms answer different relations on purpose: `:indexed` asks the source
# space's own chunk index the question a lazy read asks, `:latjoin` joins the
# caps `chunkextents` reports. Neither dominates the other (measured; see the
# `only` columns), and only the first can back a refcount. Compare their COST
# here, and read their correctness from the oracle columns.
#
# SUNSET CONDITION for the `:latjoin` arm. It reconstructs a builder that no
# longer exists, so it is not maintenance the repo owes indefinitely. Delete it
# — and this file's whole second question with it — at whichever comes first:
#
#   - G2's waived performance gate is retired (the waiver in
#     `regrid-notes/2026-08-23-g1-graph-oracles.md` §5 is what the arm exists to
#     keep auditable, and a retired waiver has nothing left to audit); or
#   - the production builder stops being comparable to a flat cap join at all,
#     so that "same relation, different cost" is no longer the question being
#     asked and the timing is measuring two unrelated things.
#
# Question 1 — the oracles — has no sunset: G3 and G4 use them to prove they did
# not change the relation.
#
# Environment
#
#   DGG_GRAPH_GATE_SAMPLES   timed samples per arm (default 5)
#   DGG_GRAPH_GATE_CASES     comma-separated case-name filter (default: all)
#   DGG_GRAPH_GATE_NDJSON    path to append one JSON line per (case, arm)
#   DGG_GRAPH_GATE_RAW       set to 1 to add the unprefiltered `:latjoin_raw`
#                            arm (default off)
#   DGG_COPDEM_TILELIST      LOCAL Copernicus tile list. The production case is
#                            SKIPPED without it; this harness never downloads.
#
# Every row is stamped with the Julia version, the thread count, and the git
# revisions of GeometryOps, GeometryOpsCore and ConservativeRegridding, all read
# from the active manifest. Clipping cost is the largest single term in this
# workload and those three are pinned to branches rather than releases, so a
# number recorded here is only comparable against another with the same stamp.

import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import GeometryOps as GO
import DimensionalData as DD
import Extents
import Graphs
import Statistics
import Pkg
using Printf: @printf

const SphCap = GO.UnitSpherical.SphericalCap

const NSAMPLES = parse(Int, get(ENV, "DGG_GRAPH_GATE_SAMPLES", "5"))
const NDJSON = get(ENV, "DGG_GRAPH_GATE_NDJSON", "")
const CASEFILTER = let s = get(ENV, "DGG_GRAPH_GATE_CASES", "")
    isempty(s) ? nothing : Set(strip.(split(s, ',')))
end
const TILELIST = get(ENV, "DGG_COPDEM_TILELIST", "")
const WITHRAW = get(ENV, "DGG_GRAPH_GATE_RAW", "") in ("1", "true", "yes")

# The oracles are shared with both test suites; see `graphoracles.jl` for why
# they live there. Do not re-spell any of them in this file.
include(joinpath(@__DIR__, "..", "lib", "GlobalRegridding", "test",
    "graphoracles.jl"))
using .ChunkGraphOracles: contributing_pairs, graph_pairs, demanded_pairs

# ===========================================================================
# Provenance
# ===========================================================================

"""
    provenance() -> NamedTuple

The identity of the measurement: interpreter, parallelism, and the exact
revisions of the geometry stack the timings are dominated by. Read from the
active manifest rather than from `Project.toml`, so a re-run after a pin bump
relabels itself instead of inheriting a stale SHA.
"""
function provenance()
    revs = Dict{String,String}()
    trees = Dict{String,String}()
    for (_, p) in Pkg.dependencies()
        p.name in ("GeometryOps", "GeometryOpsCore", "ConservativeRegridding",
            "DiscreteGlobalGrids", "GlobalRegridding") || continue
        revs[p.name] = something(p.git_revision, "")
        trees[p.name] = p.tree_hash === nothing ? "" : string(p.tree_hash)
    end
    # A branch name is not an identity; the tree hash beside it is.
    stamp(name) = (get(revs, name, ""), get(trees, name, ""))
    go, gotree = stamp("GeometryOps")
    goc, goctree = stamp("GeometryOpsCore")
    cr, crtree = stamp("ConservativeRegridding")
    return (; julia = string(VERSION), threads = Threads.nthreads(),
        gcthreads = Threads.ngcthreads(),
        geometryops_rev = go, geometryops_tree = gotree,
        geometryopscore_rev = goc, geometryopscore_tree = goctree,
        conservativeregridding_rev = cr, conservativeregridding_tree = crtree,
        repo_head = repohead())
end

function repohead()
    try
        return strip(read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String))
    catch
        return ""
    end
end

# ===========================================================================
# Arm 1 — the deleted latitude-sorted cap join, kept runnable
# ===========================================================================
#
# Verbatim in behaviour from `ba2bbfa^`, with the private helpers copied in so
# this arm cannot drift when `chunkgraph.jl` changes. It joins the caps
# `chunkextents` reports on BOTH sides, after a latitude-band prefilter over the
# source caps sorted by centre latitude.

# `atan(z, hypot(x, y))` is the robust spelling: `asin(z)` loses precision at
# the poles, which is exactly where the widest chunk caps are.
_caplatitude(c::SphCap) =
    atan(c.point[3], hypot(c.point[1], c.point[2]))

function latjoin_graph(dst_space, src_space; radius::Real = 0.0, refine = nothing,
        prefilter::Bool = true)
    r = Float64(radius)
    (isfinite(r) && r >= 0) || throw(ArgumentError("radius must be finite and >= 0"))
    # G3 gave the graph an identity; this arm stamps the same one the production
    # builder would, so the two arms stay comparable as objects and not just as
    # relations. The stamping cost is inside the timed region for both.
    id = GR.dependency_identity(dst_space, src_space; radius = r)
    return _latjoin(id, GR.chunkextents(dst_space), GR.chunkextents(src_space), r,
        refine, prefilter)
end

function _latjoin(id, dstcaps::Vector{<:SphCap}, srccaps::Vector{<:SphCap},
        radius::Float64, refine, prefilter::Bool)
    T = Int32
    ndst, nsrc = length(dstcaps), length(srccaps)
    lats = map(_caplatitude, srccaps)
    order = sortperm(lats)
    sortedlats = lats[order]
    maxsrcradius = isempty(srccaps) ? 0.0 : maximum(c -> c.radius, srccaps)

    rows = Vector{Vector{T}}(undef, ndst)
    nblocks = max(1, min(ndst, 8 * Threads.nthreads()))
    Threads.@sync for b in _blockranges(ndst, nblocks)
        let b = b
            Threads.@spawn begin
                buf = T[]
                for d in b
                    _latrow!(buf, rows, d, dstcaps[d], srccaps, order, sortedlats,
                        maxsrcradius, radius, refine, prefilter)
                end
            end
        end
    end
    return _assemble(id, rows, nsrc, ndst, T)
end

function _latrow!(buf::Vector{T}, rows::Vector{Vector{T}}, d::Int, dcap,
        srccaps::Vector, order::Vector{Int}, sortedlats::Vector{Float64},
        maxsrcradius::Float64, radius::Float64, refine, prefilter::Bool) where {T}
    empty!(buf)
    intersects = GR.DilatedIntersects(radius)
    band = maxsrcradius + dcap.radius + radius
    n = length(order)
    # A band at or above pi reaches every cap and the latitude bound degenerates.
    lo, hi = if !prefilter || !(band < Float64(pi))
        1, n
    else
        dlat = _caplatitude(dcap)
        (searchsortedfirst(sortedlats, prevfloat(dlat - band)),
         searchsortedlast(sortedlats, nextfloat(dlat + band)))
    end
    @inbounds for j in lo:hi
        s = order[j]
        intersects(dcap, srccaps[s]) || continue
        (refine === nothing || refine(d, s)::Bool) || continue
        push!(buf, T(s))
    end
    sort!(buf)
    rows[d] = copy(buf)
    return nothing
end

function _blockranges(n::Int, k::Int)
    out = UnitRange{Int}[]
    n == 0 && return out
    lo = 1
    for i in 1:k
        hi = lo + cld(n - lo + 1, k - i + 1) - 1
        push!(out, lo:hi)
        lo = hi + 1
        lo > n && break
    end
    return out
end

function _assemble(id, rows::Vector{Vector{T}}, nsrc::Int, ndst::Int,
        ::Type{T}) where {T}
    dstoff = Vector{Int}(undef, ndst + 1)
    dstoff[1] = 1
    for d in 1:ndst
        dstoff[d+1] = dstoff[d] + length(rows[d])
    end
    srcof = Vector{T}(undef, dstoff[end] - 1)
    for d in 1:ndst
        copyto!(srcof, dstoff[d], rows[d], 1, length(rows[d]))
    end
    srcoff = zeros(Int, nsrc + 1)
    @inbounds for k in eachindex(srcof)
        srcoff[srcof[k]+1] += 1
    end
    srcoff[1] = 1
    @inbounds for s in 2:(nsrc+1)
        srcoff[s] += srcoff[s-1]
    end
    dstof = Vector{T}(undef, length(srcof))
    cursor = copy(srcoff)
    @inbounds for d in 1:ndst, k in dstoff[d]:(dstoff[d+1]-1)
        s = srcof[k]
        dstof[cursor[s]] = T(d)
        cursor[s] += 1
    end
    return GR.ChunkDependencyGraph(id, dstoff, srcof, srcoff, dstof)
end

# ===========================================================================
# Measurement
# ===========================================================================

struct Timed
    seconds_min::Float64
    seconds_median::Float64
    bytes_min::Int
    gc_seconds_median::Float64
    peak_rss_growth::Int
end

function timeit(f, samples::Int; warm::Bool = true)
    value = warm ? f() : nothing       # compile against these concrete types
    times, bytes, gcs = Float64[], Int[], Float64[]
    GC.gc(true)
    rss0 = Sys.maxrss()
    for _ in 1:samples
        GC.gc(true)
        t = @timed f()
        push!(times, t.time)
        push!(bytes, t.bytes)
        push!(gcs, t.gctime)
        value = t.value
    end
    # `maxrss` is a high-water mark and never falls, so its growth across the
    # samples is the memory this arm added at its own peak.
    growth = Int(max(Sys.maxrss(), rss0) - rss0)
    return value, Timed(minimum(times), Statistics.median(times), minimum(bytes),
        Statistics.median(gcs), growth)
end

# ===========================================================================
# Case matrix
# ===========================================================================

_axis(D, centres, step) = D(DD.Sampled(collect(centres); span = DD.Regular(step),
    sampling = DD.Intervals(DD.Center()), order = DD.ForwardOrdered()))

"A global lon/lat raster of `nx x ny` cells, chunked `cx x cy` cells per chunk."
function gridraster(nx, ny, cx, cy)
    dlon, dlat = 360 / nx, 180 / ny
    lon = (-180 + dlon / 2):dlon:180
    lat = (-90 + dlat / 2):dlat:90
    array = DD.DimArray(zeros(Float32, nx, ny),
        (_axis(DD.X, lon, dlon), _axis(DD.Y, lat, dlat)))
    return GR.RasterGrid(array;
        chunks = ([lo:min(lo + cx - 1, nx) for lo in 1:cx:nx],
            [lo:min(lo + cy - 1, ny) for lo in 1:cy:ny]))
end

"A regional lon/lat raster, chunked in one block per `cx x cy` cells."
function boxraster(nx, ny, cx, cy, lon::Tuple, lat::Tuple)
    dlon = (lon[2] - lon[1]) / nx
    dlat = (lat[2] - lat[1]) / ny
    xs = (lon[1] + dlon / 2):dlon:lon[2]
    ys = (lat[1] + dlat / 2):dlat:lat[2]
    array = DD.DimArray(zeros(Float32, nx, ny),
        (_axis(DD.X, xs, dlon), _axis(DD.Y, ys, dlat)))
    return GR.RasterGrid(array;
        chunks = ([lo:min(lo + cx - 1, nx) for lo in 1:cx:nx],
            [lo:min(lo + cy - 1, ny) for lo in 1:cy:ny]))
end

const SYS7 = DGG.IGeo7System()
const S2 = DGG.S2System()

dgg(sys, level, chunklevel) = DGG.DGGSpace(DGG.levelgrid(sys, level); chunklevel)

function rooted_dgg(level, chunklevel, ordinal = 20)
    root = DGG.cellindex(DGG.levelgrid(SYS7, 1), ordinal)
    return DGG.DGGSpace(DGG.subtree(SYS7, root, level); chunklevel)
end

function sparse_dgg(level, chunklevel, extent)
    cells = DGG.covering(DGG.CellVector(DGG.levelgrid(SYS7, level)), extent)
    return DGG.DGGSpace(DGG.PartialGrid(cells); chunklevel)
end

# `oracle = true` runs the O(ncells^2) actual-cell sweep. Keep it off for cases
# whose cell counts make that a minutes-long job; those cases are covered by the
# demand-domination column instead.
const CASES = [
    (name = "raster-small", oracle = true, radius = 0.0,
     make = () -> (dgg(SYS7, 2, 1), gridraster(36, 18, 7, 5))),
    (name = "raster-support", oracle = true, radius = 0.05,
     make = () -> (dgg(SYS7, 2, 1), gridraster(36, 18, 7, 5))),
    (name = "raster-nonuniform", oracle = true, radius = 0.0,
     make = () -> (dgg(SYS7, 2, 1), gridraster(74, 38, 11, 6))),
    (name = "dgg-complete", oracle = true, radius = 0.0,
     make = () -> (dgg(SYS7, 2, 1), dgg(SYS7, 3, 2))),
    (name = "dgg-crosssystem", oracle = true, radius = 0.0,
     make = () -> (dgg(SYS7, 2, 1), dgg(S2, 3, 1))),
    (name = "dgg-rooted", oracle = true, radius = 0.0,
     make = () -> (rooted_dgg(4, 2), dgg(SYS7, 2, 1))),
    (name = "dgg-sparse", oracle = true, radius = 0.0,
     make = () -> (sparse_dgg(3, 2, Extents.Extent(X = (-140.0, 40.0), Y = (-20.0, 60.0))),
         dgg(S2, 3, 1))),
    (name = "polar-source", oracle = true, radius = 0.0,
     make = () -> (dgg(SYS7, 2, 1), boxraster(24, 6, 5, 2, (-180.0, 180.0), (60.0, 90.0)))),
    (name = "antimeridian-source", oracle = true, radius = 0.0,
     make = () -> (dgg(SYS7, 2, 1), boxraster(12, 16, 5, 5, (150.0, 180.0), (-40.0, 40.0)))),
    # The residual risk PR #69 left open: a raster source deep enough that the
    # quadtree descent per destination cap costs far more than a cap join over
    # the same chunks. 4320x2160 with 162 chunks is the shape that measurement
    # was taken on; the 1800-chunk row varies chunk count at fixed raster size.
    (name = "raster-4320-162chunks", oracle = false, radius = 0.0,
     make = () -> (dgg(SYS7, 4, 3), gridraster(4320, 2160, 240, 240))),
    (name = "raster-4320-1800chunks", oracle = false, radius = 0.0,
     make = () -> (dgg(SYS7, 4, 3), gridraster(4320, 2160, 72, 72))),
    (name = "dgg-large", oracle = false, radius = 0.0,
     make = () -> (dgg(SYS7, 5, 3), dgg(SYS7, 4, 2))),
]

# The production pair needs `listedtiles`, `TileIds`, `SubtreeIds` and
# `covering_chunks`. Loading that file pulls in ArchGDAL and Zarr, so only do it
# when the case is actually going to run.
if !isempty(TILELIST)
    isfile(TILELIST) || error("DGG_COPDEM_TILELIST=$TILELIST is not a file")
    include(joinpath(@__DIR__, "..", "scripts", "copdem_production.jl"))
end

"""
    production_case() -> case or nothing

The real Copernicus GLO-90 x IGeo7-L12 pair, built from a LOCAL tile list. The
list names tiles; no raster is opened and nothing is downloaded. Without
`DGG_COPDEM_TILELIST` the case is skipped rather than fetched.
"""
function production_case()
    isempty(TILELIST) && return nothing
    return (name = "copdem90-igeo7-l12", oracle = false, radius = 0.0,
        make = function ()
            sys = DGG.CopernicusDEMSystem(90)
            tiles = listedtiles(sys, TILELIST, nothing)
            src = DGG.DGGSpace(DGG.PartialGrid(sys, 1, TileIds(sys, tiles));
                chunklevel = 0)
            chunks = covering_chunks(SYS7, sys, tiles, 5;
                nthreads = max(1, Threads.nthreads()))
            chunkgrid = DGG.levelgrid(SYS7, 5)
            ids = SubtreeIds(SYS7, [DGG.cellindex(chunkgrid, c) for c in chunks], 12)
            dst = DGG.DGGSpace(DGG.PartialGrid(SYS7, 12, ids); chunklevel = 5)
            return (dst, src)
        end)
end

# ===========================================================================
# Runner
# ===========================================================================

const ARMS = let base = (
        (:indexed, (dst, src, r) -> GR.chunk_dependency_graph(dst, src; radius = r)),
        (:latjoin, (dst, src, r) -> latjoin_graph(dst, src; radius = r, prefilter = true)),
    )
    # The unprefiltered join is an attribution aid, not a candidate builder, and
    # it is the slowest arm in the matrix by an order of magnitude. Opt in.
    WITHRAW ? (base...,
        (:latjoin_raw,
            (dst, src, r) -> latjoin_graph(dst, src; radius = r, prefilter = false))) :
        base
end

# ===========================================================================
# Question 3 — what a row view saves against rebuilding (Task G3)
# ===========================================================================
#
# A per-column plan needs the dependency rows of one column of destinations.
# Two ways to get them: `restrict` the graph that already exists, or build a
# graph over just those destinations.
#
# `rebuild_graph` is that second path. It does every piece of work
# `chunk_dependency_graph` does for a column-sized destination — stamp both
# spaces, take the destination caps, build the source chunk index, query it once
# per destination, assemble the destination-major CSR and transpose it — inside
# the timed region, and it hands back a real `ChunkDependencyGraph`. Both sides
# of the comparison therefore produce the same kind of object and pay the same
# `O(nsourcechunks)` transpose; the difference is exactly the geometry work
# `restrict` does not repeat.
#
# What it leaves OUT, and therefore understates, is constructing the per-column
# destination space itself. A real per-column rebuild pays that too.

"Build a one-column dependency graph the way a caller without `restrict` must."
function rebuild_graph(dst, src, ds::Vector{Int}, radius::Float64)
    # Mirrors `_builddependencies`: each side's caps are taken once and serve
    # both the stamp and the relation. Task E1 made the graph keep them, so a
    # rebuild that took them twice would no longer be what a rebuild costs.
    dcaps = GR.chunkextents(dst)
    scaps = GR.chunkextents(src)
    id = GR.DependencyIdentity(GR._spacestamp(dst, dcaps), GR._spacestamp(src, scaps),
        radius, :none)
    return GR._chunkgraph(id, [dcaps[d] for d in ds], GR.chunkindex(src),
        scaps, radius, nothing)
end

rows_of(g) = [Int.(collect(GR.sourcesof(g, d))) for d in 1:GR.ndestinationchunks(g)]

"""
    restrict_measurement(graph, dst, src, radius, ndst) -> NamedTuple

Time `restrict` against a rebuild, for one destination and for a column-sized
block, and check that the two agree row for row and in both CSR directions. A
disagreement would mean a row view is not the relation it claims to be, so it is
an assertion and not a printed column.
"""
function restrict_measurement(graph, dst, src, radius, ndst)
    one = ndst == 0 ? Int[] : [max(1, ndst ÷ 2)]
    # A "column" of a 1/16th of the destination space, contiguous, which is the
    # shape a per-column plan actually asks for.
    lo = ndst == 0 ? 1 : max(1, ndst ÷ 2)
    block = collect(lo:min(ndst, lo + max(0, cld(ndst, 16) - 1)))

    vone, tone = timeit(() -> GR.restrict(graph, one), NSAMPLES)
    vblk, tblk = timeit(() -> GR.restrict(graph, block), NSAMPLES)
    bone, rone = timeit(() -> rebuild_graph(dst, src, one, radius), NSAMPLES)
    bblk, rblk = timeit(() -> rebuild_graph(dst, src, block, radius), NSAMPLES)

    # The relation a view reports must be the parent's rows for those chunks,
    # must be what a rebuild would have produced, and must transpose the same.
    consumers(g) = [Int.(collect(GR.consumersof(g, s))) for s in 1:GR.nsourcechunks(g)]
    ok = rows_of(vone) == rows_of(bone) && rows_of(vblk) == rows_of(bblk) &&
         consumers(vone) == consumers(bone) && consumers(vblk) == consumers(bblk) &&
         GR.globaldestinations(vone) == one && GR.globaldestinations(vblk) == block
    ok || error("restrict disagreed with a rebuild on $(GR.nsourcechunks(graph)) sources")

    @printf("   restrict:    1 dst %9.6f s vs rebuild %9.6f s (%6.1f×)   \
%d dst %9.6f s vs rebuild %9.6f s (%6.1f×)\n",
        tone.seconds_median, rone.seconds_median,
        rone.seconds_median / max(tone.seconds_median, eps()),
        length(block), tblk.seconds_median, rblk.seconds_median,
        rblk.seconds_median / max(tblk.seconds_median, eps()))

    return (; restrict_one_seconds = tone.seconds_median,
        restrict_one_bytes = tone.bytes_min,
        rebuild_one_seconds = rone.seconds_median,
        rebuild_one_bytes = rone.bytes_min,
        restrict_block_destinations = length(block),
        restrict_block_seconds = tblk.seconds_median,
        restrict_block_bytes = tblk.bytes_min,
        rebuild_block_seconds = rblk.seconds_median,
        rebuild_block_bytes = rblk.bytes_min)
end

function runcase(case, prov)
    (dst, src), spacetime = timeit(case.make, 1; warm = false)
    ndst, nsrc = Int(GR.nchunks(dst)), Int(GR.nchunks(src))
    @printf("\n== %s  (dst %d cells / %d chunks  <-  src %d cells / %d chunks, radius %g)\n",
        case.name, GR.ncells(dst), ndst, GR.ncells(src), nsrc, case.radius)
    @printf("   spaces built in %.3f s\n", spacetime.seconds_min)

    _, plantime = timeit(1) do
        GR.ChunkedPlan(GR.Conservative(), GR.Weighted(0.5), dst, src)
    end

    truth = case.oracle ? contributing_pairs(dst, src; radius = case.radius) : nothing
    demanded = demanded_pairs(dst, src; radius = case.radius)

    # The identity a graph now carries, and what stamping it costs. Both spaces
    # are stamped inside `chunk_dependency_graph`, so this is a component of the
    # `:indexed` timing below, not something added beside it.
    _, idtime = timeit(NSAMPLES) do
        GR.dependency_identity(dst, src; radius = case.radius)
    end
    _, indextime = timeit(() -> GR.chunkindex(src), NSAMPLES)

    rows = NamedTuple[]
    relations = Dict{Symbol,Set{Tuple{Int,Int}}}()
    measured = Dict{Symbol,Any}()
    for (arm, build) in ARMS
        graph, t = timeit(() -> build(dst, src, case.radius), NSAMPLES)
        pairs = graph_pairs(graph)
        relations[arm] = pairs
        measured[arm] = (graph, t, pairs)
    end
    restriction = restrict_measurement(measured[:indexed][1], dst, src,
        Float64(case.radius), ndst)
    for (arm, _) in ARMS
        graph, t, pairs = measured[arm]
        push!(rows, (; case = case.name, arm = String(arm),
            destination_chunks = ndst, source_chunks = nsrc,
            destination_cells = Int(GR.ncells(dst)), source_cells = Int(GR.ncells(src)),
            radius = case.radius,
            edges = Graphs.ne(graph),
            space_seconds = spacetime.seconds_min,
            plan_seconds = plantime.seconds_min,
            graph_seconds_min = t.seconds_min,
            graph_seconds_median = t.seconds_median,
            complete_plan_seconds = plantime.seconds_min + t.seconds_median,
            graph_allocated_bytes = t.bytes_min,
            graph_gc_seconds_median = t.gc_seconds_median,
            graph_summarysize_bytes = Base.summarysize(graph),
            peak_rss_growth_bytes = t.peak_rss_growth,
            per_destination_microseconds = 1e6 * t.seconds_median / max(ndst, 1),
            demanded_pairs = length(demanded),
            # A relation that misses a demanded pair cannot back a refcount.
            demand_missing = length(setdiff(demanded, pairs)),
            oracle_pairs = truth === nothing ? -1 : length(truth),
            oracle_missing = truth === nothing ? -1 : length(setdiff(truth, pairs)),
            # How this arm differs from the production relation, both ways.
            only_here = length(setdiff(pairs, relations[:indexed])),
            missing_here = length(setdiff(relations[:indexed], pairs)),
            # G3: what the identity costs, and what a row view saves.
            identity_seconds = idtime.seconds_median,
            identity_bytes = idtime.bytes_min,
            source_index_seconds = indextime.seconds_median,
            restriction...,
            samples = NSAMPLES, prov...))
    end

    # How the two relations differ, in both directions. Neither dominates the
    # other; that is the fact that made the cap join unusable as a refcount.
    ix, lj = relations[:indexed], relations[:latjoin]
    @printf("   relation: indexed %d edges, latjoin %d edges; indexed-only %d, latjoin-only %d\n",
        length(ix), length(lj), length(setdiff(ix, lj)), length(setdiff(lj, ix)))
    if haskey(relations, :latjoin_raw) && relations[:latjoin] != relations[:latjoin_raw]
        @printf("   WARNING: the latitude prefilter changed the relation by %d pairs\n",
            length(symdiff(relations[:latjoin], relations[:latjoin_raw])))
    end
    if truth !== nothing
        for (arm, _) in ARMS
            m = length(setdiff(truth, relations[arm]))
            @printf("   oracle: %-12s holds %d of %d contributing pairs%s\n",
                String(arm), length(truth) - m, length(truth), m == 0 ? "" : "  <-- MISS")
        end
    end
    for r in rows
        @printf("   %-12s %9.4f s  %10d B alloc  %10d B graph  %8.2f us/dst  %8d edges\n",
            r.arm, r.graph_seconds_median, r.graph_allocated_bytes,
            r.graph_summarysize_bytes, r.per_destination_microseconds, r.edges)
    end
    return rows, (; indexed_only = length(setdiff(ix, lj)),
        latjoin_only = length(setdiff(lj, ix)))
end

# Minimal JSON: every value here is a number, a string or a bool.
_json(v::AbstractString) = '"' * replace(String(v), '\\' => "\\\\", '"' => "\\\"") * '"'
_json(v::Symbol) = _json(String(v))
_json(v::Bool) = v ? "true" : "false"
_json(v::Real) = isfinite(v) ? string(v) : "null"
jsonline(nt::NamedTuple) = '{' * join((_json(String(k)) * ':' * _json(v)
                                       for (k, v) in pairs(nt)), ',') * '}'

function main()
    prov = provenance()
    println("chunk dependency graph gates")
    println("  julia ", prov.julia, "  threads ", prov.threads,
        "  gcthreads ", prov.gcthreads)
    println("  GeometryOps ", prov.geometryops_rev, " (tree ", prov.geometryops_tree, ")")
    println("  GeometryOpsCore ", prov.geometryopscore_rev,
        " (tree ", prov.geometryopscore_tree, ")")
    println("  ConservativeRegridding ", prov.conservativeregridding_rev,
        " (tree ", prov.conservativeregridding_tree, ")")
    println("  repo ", prov.repo_head, "  samples ", NSAMPLES)

    cases = copy(CASES)
    pc = production_case()
    pc === nothing ? println("  production case SKIPPED (set DGG_COPDEM_TILELIST)") :
        push!(cases, pc)

    allrows = NamedTuple[]
    for case in cases
        CASEFILTER === nothing || case.name in CASEFILTER || continue
        rows, _ = runcase(case, prov)
        append!(allrows, rows)
    end

    if !isempty(NDJSON)
        open(NDJSON, "a") do io
            for r in allrows
                println(io, jsonline(r))
            end
        end
        println("\nwrote ", length(allrows), " rows to ", NDJSON)
    end

    return report(allrows)
end

"""
    report(allrows) -> Int

The indexed arm's correctness verdict, and — the point of separating this out —
what it is a verdict *about*.

`oracle_missing` is `-1`, not 0, when the `O(ncells²)` sweep was skipped as too
large, and `demand_missing` is 0 for `:indexed` by construction: post-#69 the
builder issues exactly the `candidatechunks!` queries `demanded_pairs` replays,
so that column can only ever read 0 on this arm and is not evidence of anything.
Counting either as a pass is how a run over nothing but skipped cases printed an
unqualified PASS. So: count the cases that got a real geometric verdict, name
the ones that did not, and refuse to call a run with no checked case a pass.
"""
function report(allrows)
    indexed = filter(r -> r.arm == "indexed", allrows)
    checked = filter(r -> r.oracle_pairs >= 0, indexed)
    skipped = filter(r -> r.oracle_pairs < 0, indexed)
    misses = count(r -> r.oracle_missing > 0, checked)

    println("\nindexed-arm gate")
    @printf("  cases run: %d.  Oracle-checked: %d.  Oracle skipped (space too large \
             for the O(ncells^2) sweep): %d.\n",
        length(indexed), length(checked), length(skipped))
    isempty(skipped) || println("    skipped: ",
        join((r.case for r in skipped), ", "))
    for r in checked
        @printf("    %-26s %6d contributing pairs, %d missing\n",
            r.case, r.oracle_pairs, r.oracle_missing)
    end
    println("  demand_missing is 0 on this arm by construction (the builder IS the \
             `candidatechunks!` query `demanded_pairs` replays), so it is not part \
             of this verdict; read it on the `:latjoin` rows.")

    verdict = if isempty(checked)
        "NOT CHECKED — no case in this run had a geometric oracle. This run \
         proves nothing about the relation; add a case with `oracle = true`."
    elseif misses == 0
        "PASS on $(length(checked)) oracle-checked case(s)" *
        (isempty(skipped) ? "" : "; $(length(skipped)) case(s) unchecked")
    else
        "FAIL ($misses of $(length(checked)) oracle-checked cases miss pairs)"
    end
    println("  verdict: ", verdict)
    return misses
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main() == 0 ? 0 : 1)
end
