# South-pole throughput shakedown for the production CopDEM DAG executor.
#
# This is a measurement harness, not another driver.  It includes the validated
# production driver, names a deterministic Antarctic canary, samples this
# process from /proc every five seconds, and writes only to the campaign's
# scratch-store directory and ndjson record.
#
#     nice -n 10 julia --project=benchmark -t 1  --gcthreads=1 scripts/southpole_shakedown.jl select
#     nice -n 10 julia --project=benchmark -t 26 --gcthreads=4 scripts/southpole_shakedown.jl run A
#     nice -n 10 julia --project=benchmark -t 26 --gcthreads=8 scripts/southpole_shakedown.jl run B
#     # C is launched with whichever first-field mark count wins A/B:
#     nice -n 10 julia --project=benchmark -t 26 --gcthreads=N scripts/southpole_shakedown.jl run C
#     nice -n 10 julia --project=benchmark -t 26 --gcthreads=16 scripts/southpole_shakedown.jl run D
#     nice -n 10 julia --project=benchmark -t 1  --gcthreads=1 scripts/southpole_shakedown.jl crosscheck
#     nice -n 10 julia --project=benchmark -t 1  --gcthreads=1 scripts/southpole_shakedown.jl diagnose
#
# Never use a nonzero second --gcthreads field.  copdem_production.jl's gcguard
# independently refuses that known-bad concurrent-sweeper configuration.

include(joinpath(@__DIR__, "copdem_production.jl"))

import JSON3
import SHA

const SCRATCH = "/home/asinghvi17/geo/scratch-stores"
const PRODUCTION =
    "/home/asinghvi17/geo/dggstores/copdem90-igeo7-l12-synthetic.zarr"
const DATA = "/home/asinghvi17/geo/DiscreteGlobalGrids.jl/bench/data"
const RECORD = joinpath(
    "/home/asinghvi17/geo/DiscreteGlobalGrids.jl/regrid-notes",
    "2026-08-21-southpole-shakedown.ndjson")
const TARGET_COLUMNS = 2_000
const SAMPLE_SECONDS = 5.0
const RUN_TAGS = Dict(
    "A" => (mark = 4, shape = :outer),
    "B" => (mark = 8, shape = :outer),
    "D" => (mark = 16, shape = :outer),
)

utcstamp() = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS.sss") * "Z"

function appendrecord(record::Dict{String,Any})
    mkpath(dirname(RECORD))
    open(RECORD, "a") do io
        JSON3.write(io, record)
        write(io, '\n')
        flush(io)
    end
    return record
end

function emit(kind::AbstractString; kw...)
    record = Dict{String,Any}("kind" => String(kind), "ts" => utcstamp())
    for (key, value) in kw
        record[String(key)] = value
    end
    appendrecord(record)
end

function productioncolumns()
    columns = load_chunklist(PRODUCTION * ".columns.txt")
    columns === nothing && error("production column cache is missing")
    return columns
end

function productionledger()
    rows = Dict{Int,@NamedTuple{cells::Int,nan::Int}}()
    for line in eachline(PRODUCTION * ".done.ndjson")
        row = try
            JSON3.read(line)
        catch
            continue
        end
        rows[Int(row.col)] = (cells = Int(row.cells), nan = Int(row.nan))
    end
    return rows
end

"Select exactly n evenly spaced positions, including both ends when n > 1."
function evensample(values::AbstractVector, n::Integer)
    n == 0 && return eltype(values)[]
    n <= length(values) || error("cannot take $n values from $(length(values))")
    n == 1 && return [values[cld(length(values), 2)]]
    indices = [1 + fld((k - 1) * (length(values) - 1), n - 1) for k in 1:n]
    allunique(indices) || error("even-sample indices are not unique")
    return values[indices]
end

"The production-covering Antarctic canary and the facts that define it."
function columnset()
    sys = DGG.IGeo7System()
    layout = DGG.SubzoneLayout(sys, CONFIG.level, CONFIG.ancestor)
    grid = DGG.levelgrid(sys, CONFIG.ancestor)
    latitude(column) = asind(DGG.cell_centroid(
        grid, DGG.columncell(layout, column))[3])
    ispentagon(column) = DGG.IGeo7.z7_is_pentagon(
        DGG.columncell(layout, column).id)

    covering = productioncolumns()
    ledger = productionledger()
    eligible = sort!([c for c in covering if latitude(c) < -60.0];
        by = c -> (latitude(c), c))

    # This is the complete full-land block in the dress rehearsal's failure
    # band, not a sample of it.  `nan == 0` is the production ledger's direct
    # observation that every real cell in the column was land-fed.
    deep_full_land = sort!([c for c in eligible
        if latitude(c) < -80.0 && get(ledger, c, (cells = -1, nan = -1)).nan == 0])

    # IGeo7's two southernmost topological pentagons have centroids at
    # -58.2825 degrees: just outside the strict band, but explicitly retained
    # as the polar-pentagon exception required by this campaign.  They are not
    # in the production land covering, so both stores should read as fill.
    pentagons = Int[]
    for c in 1:layout.ncolumns
        ispentagon(c) || continue
        latitude(c) < 0 || continue
        push!(pentagons, c)
    end
    southmost = minimum(latitude, pentagons)
    polar_pentagons = sort!([c for c in pentagons
        if isapprox(latitude(c), southmost; atol = 1e-12)])
    length(polar_pentagons) == 2 || error(
        "expected two southernmost pentagons, found $(polar_pentagons)")

    mandatory = Set([deep_full_land; polar_pentagons])
    length(mandatory) <= TARGET_COLUMNS || error(
        "mandatory Antarctic set has $(length(mandatory)) columns")
    remaining = [c for c in eligible if !(c in mandatory)]
    sampled = evensample(remaining, TARGET_COLUMNS - length(mandatory))
    selected = sort!([collect(mandatory); sampled])
    length(selected) == TARGET_COLUMNS || error(
        "selector produced $(length(selected)) columns")
    allunique(selected) || error("selector produced duplicate columns")
    all(c -> c in selected, deep_full_land) || error("deep full-land block was sampled")
    all(c -> c in selected, polar_pentagons) || error("polar pentagon omitted")

    digest = bytes2hex(SHA.sha256(join(selected, ',')))
    return (columns = selected, sha256 = digest, covering = length(covering),
        eligible = length(eligible), deep_full_land = deep_full_land,
        polar_pentagons = polar_pentagons,
        latitude = Dict(c => latitude(c) for c in selected),
        eligible_latitude = extrema(latitude.(eligible)))
end

function recordselection()
    set = columnset()
    emit("column_set"; definition =
        "production-covering level-5 columns with centroid latitude < -60 deg; " *
        "retain every production-ledger nan==0 column below -80 deg and the two " *
        "southernmost IGeo7 pentagons, then evenly sample the remainder by " *
        "(latitude,column) rank to exactly 2000",
        target = TARGET_COLUMNS, production_covering = set.covering,
        eligible_south_of_60 = set.eligible,
        eligible_latitude_min = set.eligible_latitude[1],
        eligible_latitude_max = set.eligible_latitude[2],
        deep_full_land_count = length(set.deep_full_land),
        deep_full_land_columns = set.deep_full_land,
        polar_pentagon_columns = set.polar_pentagons,
        polar_pentagon_latitudes = [set.latitude[c] for c in set.polar_pentagons],
        columns_sha256 = set.sha256, columns = set.columns)
    println("selected $(length(set.columns)) columns; sha256=$(set.sha256)")
    println("eligible=$(set.eligible), deep-full-land=$(length(set.deep_full_land)), " *
            "polar-pentagons=$(set.polar_pentagons)")
    return set
end

"Process utime+stime in kernel clock ticks, parsed after the parenthesized comm."
function proccputicks()
    stat = read("/proc/self/stat", String)
    closeparen = findlast(==(')'), stat)
    closeparen === nothing && error("malformed /proc/self/stat")
    fields = split(SubString(stat, closeparen + 2)) # field 3 is fields[1]
    return parse(Int64, fields[12]) + parse(Int64, fields[13])
end

const CLOCK_TICKS = let ticks = ccall(:sysconf, Clong, (Cint,), 2) # _SC_CLK_TCK
    ticks > 0 || error("sysconf(_SC_CLK_TCK) failed")
    Float64(ticks)
end

function procmemory()
    rss_kib = -1
    hwm_kib = -1
    for line in eachline("/proc/self/status")
        startswith(line, "VmRSS:") && (rss_kib = parse(Int, split(line)[2]))
        startswith(line, "VmHWM:") && (hwm_kib = parse(Int, split(line)[2]))
    end
    rss_kib >= 0 || error("VmRSS missing from /proc/self/status")
    hwm_kib >= 0 || error("VmHWM missing from /proc/self/status")
    return (rss_gib = rss_kib / 2^20, hwm_gib = hwm_kib / 2^20)
end

mutable struct Telemetry
    tag::String
    started::Float64
    startticks::Int64
    lasttime::Float64
    lastticks::Int64
    stop::Threads.Atomic{Bool}
    cores::Vector{Float64}
    rss::Vector{Float64}
    hwm::Vector{Float64}
end

function Telemetry(tag)
    now = time()
    ticks = proccputicks()
    Telemetry(tag, now, ticks, now, ticks, Threads.Atomic{Bool}(false),
        Float64[], Float64[], Float64[])
end

function takesample!(telemetry::Telemetry; initial = false, final = false)
    now = time()
    ticks = proccputicks()
    elapsed = now - telemetry.started
    interval = now - telemetry.lasttime
    cores = interval > 0 ? (ticks - telemetry.lastticks) / CLOCK_TICKS / interval : 0.0
    memory = procmemory()
    initial || push!(telemetry.cores, cores)
    push!(telemetry.rss, memory.rss_gib)
    push!(telemetry.hwm, memory.hwm_gib)
    emit("sample"; config = telemetry.tag, elapsed_s = elapsed,
        interval_s = interval, cpu_ticks = ticks,
        instantaneous_cores = initial ? nothing : cores,
        rss_gib = memory.rss_gib, peak_rss_gib = memory.hwm_gib,
        gc_live_gib = Base.gc_live_bytes() / 2^30,
        initial = initial, final = final)
    telemetry.lasttime = now
    telemetry.lastticks = ticks
    return nothing
end

function startsampler(tag)
    telemetry = Telemetry(tag)
    takesample!(telemetry; initial = true)
    task = Threads.@spawn begin
        while !telemetry.stop[]
            Base.timedwait(() -> telemetry.stop[], SAMPLE_SECONDS; pollint = 0.05) ==
                :ok && break
            try
                takesample!(telemetry)
            catch err
                emit("sampler_error"; config = telemetry.tag,
                    error = sprint(showerror, err, catch_backtrace()))
            end
        end
    end
    return telemetry, task
end

function stopsampler!(telemetry, task)
    telemetry.stop[] = true
    wait(task)
    takesample!(telemetry; final = true)
    return telemetry
end

function ledgerrows(path)
    rows = NamedTuple[]
    isfile(path) || return rows
    for line in eachline(path)
        row = try
            JSON3.read(line)
        catch
            continue
        end
        push!(rows, (column = Int(row.col), cells = Int(row.cells),
            nan = Int(row.nan), seconds = Float64(row.secs), worker = Int(row.w),
            completed = Dates.datetime2unix(Dates.DateTime(String(row.t)))))
    end
    return rows
end

function runtimes(rows)
    isempty(rows) && return (compute_wall = NaN, last10_span = NaN)
    compute_start = minimum(row.completed - row.seconds for row in rows)
    compute_stop = maximum(row.completed for row in rows)
    ordered = sort(rows; by = row -> row.completed)
    tail = ordered[max(1, end - 9):end]
    return (compute_wall = compute_stop - compute_start,
        last10_span = maximum(row.completed for row in tail) -
                      minimum(row.completed for row in tail))
end

function runconfig(tag::String)
    nsweep = Base.JLOptions().nsweepthreads
    nsweep == 0 || error("refusing unsafe concurrent sweeper: nsweepthreads=$nsweep")
    Threads.nthreads() == 26 || error("campaign runs require -t 26")
    spec = if tag == "C"
        Base.JLOptions().nmarkthreads in (4, 8, 16) || error(
            "C must reuse a measured first-field mark count")
        (mark = Int(Base.JLOptions().nmarkthreads), shape = :inner)
    else
        get(RUN_TAGS, tag, nothing)
    end
    spec === nothing && error("run tag must be A, B, C, or D")
    Base.JLOptions().nmarkthreads == spec.mark || error(
        "$tag requires --gcthreads=$(spec.mark), got $(Base.JLOptions().nmarkthreads)")

    set = columnset()
    config_name = "$tag-gc$(spec.mark)-$(spec.shape)"
    store = joinpath(SCRATCH, "southpole-$config_name.zarr")
    mkpath(SCRATCH)
    islink(SCRATCH) && error("scratch root must not be a symlink")
    abspath(dirname(store)) == SCRATCH || error("store escaped scratch root")
    for path in (store, donelogpath(store), chunklistpath(store))
        ispath(path) && error("fresh-store requirement: $path already exists")
    end
    save_chunklist(chunklistpath(store), CONFIG.ancestor, set.columns)

    config = merge(CONFIG, (source = :synthetic, store = store, region = nothing,
        real = :none, data = DATA, workers = 0, cores = 24, shape = spec.shape,
        schedule = :affinity, cachepolicy = :refcount, taper = true,
        prefetch = 0, fetchdelay = 0.0, resume = false, checks = false,
        heartbeat = 300, maxchunks = 0, chunks = set.columns, dryrun = false))
    workers, resolved_shape = workercount(config)
    expected_workers = spec.shape === :outer ? 23 : 8
    workers == expected_workers || error(
        "expected $expected_workers workers for $(spec.shape), resolved $workers")

    FAILURES[] = 0
    LASTCACHE[] = nothing
    emit("run_start"; config = config_name, campaign_tag = tag,
        pid = getpid(), julia = string(VERSION), threads = Threads.nthreads(),
        gc_mark_threads = spec.mark, gc_sweep_threads = nsweep,
        shape = String(resolved_shape), workers = workers, cores_budget = config.cores,
        source = String(config.source), prefetch = config.prefetch,
        store = store, columns = length(set.columns), columns_sha256 = set.sha256,
        nice = try Base.Libc.getpriority(0, 0) catch; 10 end)

    telemetry, task = startsampler(config_name)
    failures = -1
    caught = nothing
    try
        failures = main(config)
    catch err
        caught = (err, catch_backtrace())
    finally
        stopsampler!(telemetry, task)
    end
    if caught !== nothing
        emit("run_error"; config = config_name,
            wall_s = time() - telemetry.started,
            error = sprint(showerror, caught[1], caught[2]))
        Base.display_error(Base.stderr, caught[1], caught[2])
        return 1
    end

    rows = ledgerrows(donelogpath(store))
    timings = runtimes(rows)
    completed = length(rows)
    cells = sum(row.cells for row in rows; init = 0)
    cpu_seconds = (proccputicks() - telemetry.startticks) / CLOCK_TICKS
    wall = time() - telemetry.started
    cache = LASTCACHE[]
    memory = procmemory()
    core_s_per_column = cpu_seconds / max(completed, 1)
    samplecores = filter(isfinite, telemetry.cores)
    regime = abs(core_s_per_column - 8.3) <= abs(core_s_per_column - 16.5) ?
        "epoch-B-like" : "epoch-C-like"
    emit("run_final"; config = config_name, campaign_tag = tag,
        exit_failures = failures, columns = completed, cells = cells,
        wall_s = wall, compute_wall_s = timings.compute_wall,
        aggregate_cells_s = cells / timings.compute_wall,
        cpu_s = cpu_seconds, core_s_per_column = core_s_per_column,
        cores_mean = cpu_seconds / wall,
        cores_sample_median = isempty(samplecores) ? nothing : Statistics.median(samplecores),
        cores_sample_p10 = isempty(samplecores) ? nothing : Statistics.quantile(samplecores, 0.1),
        cores_sample_p90 = isempty(samplecores) ? nothing : Statistics.quantile(samplecores, 0.9),
        peak_rss_gib = memory.hwm_gib,
        sampled_peak_rss_gib = maximum(telemetry.rss),
        end_rss_gib = memory.rss_gib,
        cache_peak_gib = cache === nothing ? nothing : cache.peakbytes / 2^30,
        cache_peak_tiles = cache === nothing ? nothing : cache.peaktiles,
        cache_loads = cache === nothing ? nothing : cache.loads,
        cache_uncredited = cache === nothing ? nothing : cache.uncredited,
        cache_live_end = cache === nothing ? nothing : cache.live,
        cold_downloads = 0, last10_columns_span_s = timings.last10_span,
        workers = workers, shape = String(resolved_shape), gc_mark_threads = spec.mark,
        regime = regime, columns_sha256 = set.sha256, store = store)
    println("FINAL $config_name: $(round(core_s_per_column; digits=3)) core-s/column, " *
            "$(round(cpu_seconds / wall; digits=2)) cores, " *
            "$(round(memory.hwm_gib; digits=2)) GiB peak, $regime")
    return failures == 0 && completed == TARGET_COLUMNS ? 0 : 1
end

function crosscheck()
    set = columnset()
    astore = joinpath(SCRATCH, "southpole-A-gc4-outer.zarr")
    isdir(astore) || error("config A store is missing: $astore")
    isdir(PRODUCTION) || error("production store is missing")

    pentagons = set.polar_pentagons
    deep = sort!([c for c in set.columns if set.latitude[c] < -80.0];
        by = c -> (set.latitude[c], c))
    band = sort!([c for c in set.columns
        if set.latitude[c] >= -80.0 && !(c in pentagons)];
        by = c -> (set.latitude[c], c))
    sampled = sort!(unique([pentagons; evensample(deep, 6); evensample(band, 16)]))
    length(sampled) >= 20 || error("cross-check selected only $(length(sampled)) columns")
    count(c -> c in pentagons, sampled) >= 2 || error("cross-check lacks pentagons")
    count(c -> set.latitude[c] < -80.0, sampled) >= 5 || error("cross-check lacks deep columns")

    production = Zarr.zopen(PRODUCTION, "r")["elevation"]
    candidate = Zarr.zopen(astore, "r")["elevation"]
    total_bits = 0
    total_masks = 0
    total_cells = 0
    maximum_residual = 0.0
    for column in sampled
        expected = Vector{Float32}(production[:, column])
        observed = Vector{Float32}(candidate[:, column])
        length(expected) == length(observed) || error("column $column length mismatch")
        masks = count(i -> isnan(expected[i]) != isnan(observed[i]), eachindex(expected))
        bits = count(i -> reinterpret(UInt32, expected[i]) !=
                          reinterpret(UInt32, observed[i]), eachindex(expected))
        residual = 0.0
        for i in eachindex(expected)
            if isfinite(expected[i]) && isfinite(observed[i])
                residual = max(residual,
                    abs(Float64(expected[i]) - Float64(observed[i])))
            end
        end
        category = column in pentagons ? "pentagon" :
            set.latitude[column] < -80.0 ? "south-of-80" : "band"
        emit("crosscheck_column"; column = column,
            latitude = set.latitude[column], category = category,
            cells = length(expected), nan_mask_mismatches = masks,
            non_bit_equal_float32 = bits, max_abs_residual_m = residual)
        total_cells += length(expected)
        total_masks += masks
        total_bits += bits
        maximum_residual = max(maximum_residual, residual)
        println("crosscheck $column ($category): residual=$residual, bits=$bits, masks=$masks")
    end
    emit("crosscheck_final"; config = "A-gc4-outer", columns = sampled,
        column_count = length(sampled), cells = total_cells,
        pentagon_count = count(c -> c in pentagons, sampled),
        south_of_80_count = count(c -> set.latitude[c] < -80.0, sampled),
        nan_mask_mismatches = total_masks, non_bit_equal_float32 = total_bits,
        max_abs_residual_m = maximum_residual,
        threshold_m = 1e-3, passed = total_masks == 0 && maximum_residual <= 1e-3,
        production_store = PRODUCTION, candidate_store = astore)
    println("CROSSCHECK: $(length(sampled)) columns, $total_cells cells, " *
            "max residual $maximum_residual m, $total_bits non-bit-equal, " *
            "$total_masks NaN-mask mismatches")
    return total_masks == 0 && maximum_residual <= 1e-3 ? 0 : 1
end

"Explain the two pole-column residuals without mutating either store."
function polediagnostic()
    stores = Dict(
        "production" => PRODUCTION,
        "A" => joinpath(SCRATCH, "southpole-A-gc4-outer.zarr"),
        "B" => joinpath(SCRATCH, "southpole-B-gc8-outer.zarr"),
        "C" => joinpath(SCRATCH, "southpole-C-gc4-inner.zarr"),
    )
    all(isdir, values(stores)) || error("pole diagnostic needs production and A/B/C stores")
    arrays = Dict(tag => Zarr.zopen(path, "r")["elevation"]
        for (tag, path) in stores)
    sys = DGG.IGeo7System()
    layout = DGG.SubzoneLayout(sys, CONFIG.level, CONFIG.ancestor)
    grid = DGG.levelgrid(sys, CONFIG.level)

    failed = false
    for column in (123203, 123204)
        values = Dict(tag => Vector{Float32}(array[:, column])
            for (tag, array) in arrays)
        a = values["A"]
        production = values["production"]
        ab = count(i -> reinterpret(UInt32, a[i]) !=
                        reinterpret(UInt32, values["B"][i]), eachindex(a))
        ac = count(i -> reinterpret(UInt32, a[i]) !=
                        reinterpret(UInt32, values["C"][i]), eachindex(a))
        differing = [i for i in eachindex(a)
            if reinterpret(UInt32, a[i]) != reinterpret(UInt32, production[i])]
        finite = [i for i in differing if isfinite(a[i]) && isfinite(production[i])]
        above = count(i -> abs(Float64(a[i]) - Float64(production[i])) > 1e-3, finite)
        worst = isempty(finite) ? nothing : finite[argmax(
            [abs(Float64(a[i]) - Float64(production[i])) for i in finite])]
        maximum_residual = worst === nothing ? 0.0 :
            abs(Float64(a[worst]) - Float64(production[worst]))
        production_out_of_range = count(v -> isfinite(v) && !( -1400.0 <= v <= 1600.0),
            production)
        candidate_out_of_range = count(v -> isfinite(v) && !( -1400.0 <= v <= 1600.0), a)

        row = worst === nothing ? nothing : worst
        cell = row === nothing ? nothing : DGG.cellindex(grid,
            first(DGG.columnpositions(layout, column)) + row - 1)
        point = cell === nothing ? nothing : DGG.cell_centroid(grid, cell)
        lon = point === nothing ? nothing : atand(point[2], point[1])
        lat = point === nothing ? nothing : asind(point[3])
        oracle = point === nothing ? nothing : SYNTHETIC(point)
        emit("pole_diagnostic"; column = column,
            a_vs_b_non_bit_equal = ab, a_vs_c_non_bit_equal = ac,
            production_non_bit_equal = length(differing),
            production_residual_gt_1e3 = above,
            max_abs_residual_m = maximum_residual,
            production_outside_source_range = production_out_of_range,
            candidate_outside_source_range = candidate_out_of_range,
            source_range_m = [-1400.0, 1600.0], worst_row = row,
            worst_longitude = lon, worst_latitude = lat,
            candidate_value_m = row === nothing ? nothing : a[row],
            production_value_m = row === nothing ? nothing : production[row],
            analytic_centroid_m = oracle,
            candidate_centroid_residual_m = row === nothing ? nothing :
                abs(Float64(a[row]) - oracle),
            production_centroid_residual_m = row === nothing ? nothing :
                abs(Float64(production[row]) - oracle))
        println("diagnose $column: A/B bits=$ab, A/C bits=$ac, " *
            "production bits=$(length(differing)), >1e-3=$above, " *
            "max=$maximum_residual m, production out-of-range=$production_out_of_range")
        failed |= ab != 0 || ac != 0 || candidate_out_of_range != 0
    end
    emit("pole_diagnostic_final"; passed_fresh_determinism = !failed,
        conclusion = failed ? "fresh stores disagree or leave the analytic range" :
        "A/B/C are bit-identical at both pole columns; production alone contains " *
        "values outside the synthetic field's exact global range")
    return failed ? 1 : 0
end

function main_shakedown()
    isempty(ARGS) && error(
        "usage: southpole_shakedown.jl select | run A|B|C|D | crosscheck | diagnose")
    if ARGS[1] == "select"
        length(ARGS) == 1 || error("select takes no arguments")
        recordselection()
        return 0
    elseif ARGS[1] == "run"
        length(ARGS) == 2 || error("run requires one tag: A, B, C, or D")
        return runconfig(uppercase(ARGS[2]))
    elseif ARGS[1] == "crosscheck"
        length(ARGS) == 1 || error("crosscheck takes no arguments")
        return crosscheck()
    elseif ARGS[1] == "diagnose"
        length(ARGS) == 1 || error("diagnose takes no arguments")
        return polediagnostic()
    end
    error("unknown mode $(repr(ARGS[1]))")
end

exit(main_shakedown())
