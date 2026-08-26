# What a point stencil costs against nearest and conservative on real
# Copernicus DEM tiles regridded to IGeo7.
#
#     julia -t auto --project=benchmark benchmark/copdem_point.jl \
#         res=90 box=10,11,46,47
#
# `benchmark/copdem_nearest.jl` prices `NearestCell` on this workload and owns
# the source, the destination, the tile cache and the counters; this file
# includes it and reuses all of them, so the two agree on what a "tile", a
# "chunk read" and a "sweep" are by construction rather than by convention.
#
# ARMS. One destination and one source shape, and three methods in ONE process,
# so the ratios between them belong to one machine state:
#
#   1. `destination`   the IGeo7 cells the box covers at `level`, as a
#                      `PartialGrid` tiled at `chunk`. Paid once, shared.
#   for each of `NearestCell`, `BarycentricPoint`, `Conservative`:
#   2. `plan`          `plan_regrid`, which builds the dependency relation and
#                      reads no data.
#   3. `cold sweep`    every destination tile read once, in order, on a plan
#                      never read and a source whose tile cache is empty.
#   4. `warm sweep`    the same sweep against the plan arm 3 left warm.
#
# Each method gets a source object of its own, so no tile one arm decoded is
# resident for the next; the GeoTIFFs stay in the OS page cache either way, so
# what a cold sweep pays is decode and weights, not disk.
#
# WHAT IS COUNTED: source chunk reads and pixels copied, from the counting
# `DiskArrays` wrapper; GeoTIFF decodes from `TileBuilder`'s own atomics; peak
# RSS from `Sys.maxrss()`, a high-water mark for the whole process and therefore
# including the warm-up.
#
# CORRECTNESS. `check=N` verifies N random destination cells of the barycentric
# arm against the stencil recomputed here and the source pixels read directly:
# the value must be the weighted sum of the posts `weightsat!` names, and a cell
# whose stencil the source cannot complete must be missing. This checks the
# executor against the stencil, not the stencil against geometry — the CopDEM
# section of `test/systems/CopernicusDEM/runtests.jl` does that.
#
# WARM-UP AND STATISTIC. Every arm is a one-shot measurement of a phase running
# for seconds; there is no min-of-n. Compilation is removed instead: each method
# first runs one complete miniature regrid over a small box at the same level.
#
# CONFIGURATION is `copdem_nearest.jl`'s, without `method` and `profile`.
# `sweeps=0` runs the plan arm alone.
#
# No data is committed and no path here is absolute.

include(joinpath(@__DIR__, "copdem_nearest.jl"))

const ARMS = ("nearest", "barycentric", "conservative")

armmethod(name) = name == "barycentric" ? GR.BarycentricPoint() : buildmethod(name)

"`copdem_nearest.jl`'s options, with the two this file does not take removed."
pointoptions(args = ARGS) =
    options(filter(a -> !startswith(a, "method=") && !startswith(a, "profile="), args))

"""
    stencilpass(smp, sites, n) -> (mapped, entries)

One `weightsat!` per destination sample site, into a row reused across all of
them and nothing else: the stencil's own cost, with no source read and no weight
stored. `mapped` counts the sites the source could place and `entries` the
nonzeros they emitted.
"""
function stencilpass(smp, sites, n::Int)
    row = GR.WeightRow()
    mapped = 0
    entries = 0
    for i in 1:n
        GR.ismapped(GR.weightsat!(row, smp, sites[i])) || continue
        mapped += 1
        entries += length(row)
    end
    return (mapped, entries)
end

"One method's plan and its two sweeps, on a source whose tile cache is empty."
function arm(name, opt, dir, dst, ndst)
    method = armmethod(name)
    src = sourceside(opt, dir)
    println("\n$name")
    p = phase("plan", () -> GR.plan_regrid(src.dem; to = dst.space, from = src.space,
        method = method, missingpolicy = DGG.Weighted(0.5), lazy = true,
        budget = opt.budget))
    plan = p.value
    graph = GR.dependencies(plan)
    cand = [length(GR.sourcesof(graph, c)) for c in 1:Int(GR.nchunks(dst.space))]
    @printf("  %-26s min %d, median %.1f, max %d\n", "candidates per tile",
        minimum(cand), Statistics.median(cand), maximum(cand))
    opt.sweeps || return (; name, cold = p, warm = p, cand, check = (0, 0, 0, 0))

    if name == "barycentric"
        smp = GR.sampler(GR.BarycentricPoint(), src.space)
        sites = GR.samplesites(dst.space)
        stencilpass(smp, sites, min(ndst, 1024))          # warm the row's buffers
        st = phase("stencils", () -> stencilpass(smp, sites, ndst))
        @printf("  %-26s %.1f ns and %.1f B each, %s entries over %s placed\n",
            "one stencil", 1e9 * st.time / ndst, st.bytes / ndst,
            commas(st.value[2]), commas(st.value[1]))
        # The same stencils against the COMPLETE level, where a post's index is
        # arithmetic rather than a search of a holding's ids. The difference is
        # what `localindex` on the holding costs, four times a query.
        whole = GR.sampler(GR.BarycentricPoint(),
            DGG.DGGSpace(DGG.levelgrid(src.sys, 1); chunklevel = 0))
        stencilpass(whole, sites, min(ndst, 1024))
        sw = phase("stencils, complete level", () -> stencilpass(whole, sites, ndst))
        @printf("  %-26s %.1f ns each, %s entries over %s placed\n",
            "one stencil", 1e9 * sw.time / ndst, commas(sw.value[2]),
            commas(sw.value[1]))
        whole = nothing
    end

    out = Vector{Float32}(undef, ndst)
    A = GR.regrid(src.dem, plan)
    cold = phase("cold sweep", () -> sweep!(out, A, dst.space))
    coldreads = sum(values(src.dem.reads); init = 0)
    coldpixels = src.dem.pixels[]
    colddecodes = src.builder.nreal[]
    resetreads!(src.dem)
    warm = phase("warm sweep", () -> sweep!(out, A, dst.space))
    nfinite = count(isfinite, out)
    @printf("  %-26s %s reads, %s pixels, %d GeoTIFF decodes\n", "cold source",
        commas(coldreads), commas(coldpixels), colddecodes)
    @printf("  %-26s %s of %s (%.2f%%), %s resident weight B of %s\n", "finite values",
        commas(nfinite), commas(ndst), 100nfinite / ndst,
        commas(weightstate(plan).accounted), commas(GR.weightbudget(opt.budget)))
    # Verified here, while this arm's source and values are still alive, so no
    # arm holds another's plan or weights and the peak below is one arm's.
    check = name == "barycentric" ?
            verifypoint(out, dst, src, src.dem, opt.check) : (0, 0, 0, 0)
    name == "barycentric" && @printf("  %-26s %s cells, %s mapped, %s outside, %s wrong\n",
        "check", commas(check[1]), commas(check[2]), commas(check[3]), commas(check[4]))
    @printf("  %-26s %.2f GiB\n", "peak RSS after this arm", Sys.maxrss() / 2^30)
    flush(stdout)
    return (; name, cold, warm, cand, coldreads, coldpixels, colddecodes, nfinite,
        check)
end

"""
    verifypoint(out, dst, src, dem, n) -> (checked, mapped, missing, wrong)

`n` random destination cells of the barycentric arm, against the stencil
recomputed here and the source pixels read directly.
"""
function verifypoint(out, dst, src, dem, n::Int; seed = 20260826)
    n <= 0 && return (0, 0, 0, 0)
    rng = Random.MersenneTwister(seed)
    smp = GR.sampler(GR.BarycentricPoint(), src.space)
    row = GR.WeightRow()
    sites = GR.samplesites(dst.space)
    checked = mapped = absent = wrong = 0
    for i in Random.rand(rng, 1:DGG.ncells(dst.grid), n)
        checked += 1
        if !GR.ismapped(GR.weightsat!(row, smp, sites[i]))
            absent += 1
            isfinite(out[i]) && (wrong += 1)
            continue
        end
        mapped += 1
        want = 0.0
        good = true
        for k in 1:length(row)
            v = dem[row.indices[k]]
            isfinite(v) || (good = false; break)
            want += row.weights[k] * v
        end
        good || continue
        isapprox(out[i], Float32(want); rtol = 1e-6) || (wrong += 1)
    end
    return (checked, mapped, absent, wrong)
end

function run(opt = pointoptions())
    dir = tiledir(opt)
    sys7 = DGG.IGeo7System()
    println(provenance())
    @printf("\nres=GLO-%d  box=(%d,%d,%d,%d)  budget=%s B  tiles from %s\n",
        opt.res, opt.box..., commas(opt.budget), dir)

    probe = sourceside(opt, dir)
    level = opt.level === nothing ? DGG.levelfor(sys7, probe.space) : opt.level
    chunklevel = opt.chunk === nothing ? max(first(DGG.levels(sys7)), level - 7) : opt.chunk
    @printf("source   %d tiles, %s pixels, %s B of GeoTIFF, %.2f m pixels\n",
        length(probe.tiles), commas(DGG.ncells(probe.grid)), commas(probe.bytes),
        DGG.cellsize(probe.space))
    @printf("dest     IGeo7 level %d, chunk level %d, %.2f m cells\n",
        level, chunklevel, DGG.cellsize(sys7, level))
    flush(stdout)

    # Compilation, on sources of their own, so no arm below starts warm.
    warm = @timed for name in ARMS
        w = sourceside(opt, dir)
        ex = Extents.Extent(X = (Float64(opt.box[1]) + 0.40, Float64(opt.box[1]) + 0.42),
            Y = (Float64(opt.box[3]) + 0.40, Float64(opt.box[3]) + 0.42))
        set = DGG.query(sys7, DGG.MultiOrderCoverage(ex); level = level)
        g = DGG.PartialGrid(DGG.CellVector(set; level = level))
        sp = DGG.DGGSpace(g; chunklevel = chunklevel)
        p = GR.plan_regrid(w.dem; to = sp, from = w.space, method = armmethod(name),
            missingpolicy = DGG.Weighted(0.5), lazy = true, budget = opt.budget)
        collect(GR.regrid(w.dem, p))
    end
    @printf("\nwarm-up  %.3f s for all three methods (compilation; in no arm below)\n",
        warm.time)

    println("\narms")
    a1 = phase("destination", () -> destinationside(sys7, opt.box, level, chunklevel))
    dst = a1.value
    ndst = DGG.ncells(dst.grid)
    @printf("  %-26s %s cells in %s tiles\n", "destination",
        commas(ndst), commas(Int(GR.nchunks(dst.space))))
    flush(stdout)

    probe = nothing                       # the sizing source is not an arm's
    results = [arm(name, opt, dir, dst, ndst) for name in ARMS]

    if opt.sweeps
        println("\n$(Threads.nthreads()) threads, $(commas(ndst)) destination cells")
        @printf("%-14s %10s %12s %12s %10s %12s %12s\n",
            "method", "cold s", "cold ns/cell", "cold cells/s",
            "warm s", "warm ns/cell", "warm cells/s")
        for r in results
            @printf("%-14s %10.3f %12.1f %12s %10.3f %12.1f %12s\n", r.name,
                r.cold.time, 1e9 * r.cold.time / ndst,
                commas(round(Int, ndst / r.cold.time)),
                r.warm.time, 1e9 * r.warm.time / ndst,
                commas(round(Int, ndst / r.warm.time)))
        end
        base = results[1].cold.time
        println()
        for r in results
            @printf("%-14s cold sweep %.2fx nearest's\n", r.name, r.cold.time / base)
        end
    end
    @printf("\npeak RSS   %.2f GiB\n", Sys.maxrss() / 2^30)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
