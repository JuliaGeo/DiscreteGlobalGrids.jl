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
# Follow-up thread-scaling campaign (fixed 24-core sizing budget):
#
#     nice -n 10 julia --project=benchmark -t 1  --gcthreads=1 scripts/southpole_shakedown.jl scale-select
#     nice -n 10 julia --project=benchmark -t 32 --gcthreads=4 scripts/southpole_shakedown.jl scale-run outer-t32
#     nice -n 10 julia --project=benchmark -t 32 --gcthreads=4 scripts/southpole_shakedown.jl scale-run inner-t32
#     nice -n 10 julia --project=benchmark -t 32 --gcthreads=4 scripts/southpole_shakedown.jl scale-run inner-t32-Wboost
#     nice -n 10 julia --project=benchmark -t 21 --gcthreads=4 scripts/southpole_shakedown.jl scale-run outer-t21
#     nice -n 10 julia --project=benchmark -t 21 --gcthreads=4 scripts/southpole_shakedown.jl scale-run inner-t21
#     nice -n 10 julia --project=benchmark -t 1  --gcthreads=1 scripts/southpole_shakedown.jl scale-crosscheck
#
# Follow-up t21 utilization-attribution campaign (the command supervises a
# fresh t21 child so /proc sampling is independent of Julia-pool load):
#
#     nice -n 10 julia --project=benchmark -t 1 --gcthreads=1 scripts/southpole_shakedown.jl util-select
#     nice -n 10 julia --project=benchmark -t 1 --gcthreads=1 scripts/southpole_shakedown.jl util-run baseline
#     nice -n 10 julia --project=benchmark -t 1 --gcthreads=1 scripts/southpole_shakedown.jl util-run W32-gc4
#     nice -n 10 julia --project=benchmark -t 1 --gcthreads=1 scripts/southpole_shakedown.jl util-crosscheck
#
# Never use a nonzero second --gcthreads field.  copdem_production.jl's gcguard
# independently refuses that known-bad concurrent-sweeper configuration.

include(joinpath(@__DIR__, "copdem_production.jl"))

import JSON3
import Profile
import SHA

const SCRATCH = "/home/asinghvi17/geo/scratch-stores"
const PRODUCTION =
    "/home/asinghvi17/geo/dggstores/copdem90-igeo7-l12-synthetic.zarr"
const DATA = "/home/asinghvi17/geo/DiscreteGlobalGrids.jl/bench/data"
const RECORD = joinpath(
    "/home/asinghvi17/geo/DiscreteGlobalGrids.jl/regrid-notes",
    "2026-08-21-southpole-shakedown.ndjson")
const SCALING_RECORD = joinpath(
    "/home/asinghvi17/geo/DiscreteGlobalGrids.jl/regrid-notes",
    "2026-08-21-inner-scaling.ndjson")
const UTIL_RECORD = joinpath(
    "/home/asinghvi17/geo/DiscreteGlobalGrids.jl/regrid-notes",
    "2026-08-22-t21-utilization.ndjson")
const ACTIVE_RECORD = Ref(RECORD)
const TARGET_COLUMNS = 2_000
const UTIL_TARGET_COLUMNS = 800
const TARGET_COLUMNS_SHA256 =
    "62cef2b46fcac83f308f4a6314e8acd80ebae71f9b96b72fe8d3380184dd1635"
const SAMPLE_SECONDS = 5.0
const UTIL_PROC_SECONDS = 3.0
const UTIL_MAX_WALL_SECONDS = 40 * 60.0
const UTIL_MAX_RSS_GIB = 40.0
const RUN_TAGS = Dict(
    "A" => (mark = 4, shape = :outer),
    "B" => (mark = 8, shape = :outer),
    "D" => (mark = 16, shape = :outer),
)
const SCALING_RUNS = Dict(
    "outer-t32" => (threads = 32, shape = :outer, workers = 0),
    "inner-t32" => (threads = 32, shape = :inner, workers = 0),
    # The driver-supported override is round(2/3 * outer W) = round(2/3 * 23).
    "inner-t32-Wboost" => (threads = 32, shape = :inner, workers = 15),
    "outer-t21" => (threads = 21, shape = :outer, workers = 0),
    "inner-t21" => (threads = 21, shape = :inner, workers = 0),
)
const UTIL_RUNS = Dict(
    "baseline" => (threads = 21, mark = 4, shape = :outer, workers = 0,
        batch = 8, taper = true, schedule = :affinity),
    "W32-gc4" => (threads = 21, mark = 4, shape = :outer, workers = 32,
        batch = 8, taper = true, schedule = :affinity),
    "W40-gc4" => (threads = 21, mark = 4, shape = :outer, workers = 40,
        batch = 8, taper = true, schedule = :affinity),
    "W32-gc8" => (threads = 21, mark = 8, shape = :outer, workers = 32,
        batch = 8, taper = true, schedule = :affinity),
    "W40-gc8" => (threads = 21, mark = 8, shape = :outer, workers = 40,
        batch = 8, taper = true, schedule = :affinity),
    # Reserved targeted existing-knob variant if attribution implicates taper.
    "W40-batch1-gc4" => (threads = 21, mark = 4, shape = :outer, workers = 40,
        batch = 1, taper = true, schedule = :affinity),
)

utcstamp() = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS.sss") * "Z"

function appendrecord(record::Dict{String,Any})
    path = ACTIVE_RECORD[]
    mkpath(dirname(path))
    buffer = IOBuffer()
    JSON3.write(buffer, record)
    write(buffer, '\n')
    payload = take!(buffer)
    open(path, "a") do io
        # One append write keeps parent /proc samples and child run records from
        # interleaving when both processes share the utilization record.
        write(io, payload)
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

"Recover the already-recorded canonical selection if its source sidecars are gone."
function recordedcolumnset()
    canonical = nothing
    for line in eachline(RECORD)
        row = try
            JSON3.read(line)
        catch
            continue
        end
        if hasproperty(row, :kind) && String(row.kind) == "column_set"
            canonical = row
        end
    end
    canonical === nothing && error("canonical column_set record is missing from $RECORD")
    columns = Int.(canonical.columns)
    deep_full_land = Int.(canonical.deep_full_land_columns)
    polar_pentagons = Int.(canonical.polar_pentagon_columns)
    digest = bytes2hex(SHA.sha256(join(columns, ',')))
    recorded_digest = String(canonical.columns_sha256)
    digest == recorded_digest == TARGET_COLUMNS_SHA256 || error(
        "canonical column set digest mismatch: computed=$digest, " *
        "recorded=$recorded_digest, required=$TARGET_COLUMNS_SHA256")
    length(columns) == TARGET_COLUMNS || error(
        "canonical record has $(length(columns)) columns, expected $TARGET_COLUMNS")
    allunique(columns) || error("canonical record has duplicate columns")

    sys = DGG.IGeo7System()
    layout = DGG.SubzoneLayout(sys, CONFIG.level, CONFIG.ancestor)
    grid = DGG.levelgrid(sys, CONFIG.ancestor)
    latitude(column) = asind(DGG.cell_centroid(
        grid, DGG.columncell(layout, column))[3])
    ispentagon(column) = DGG.IGeo7.z7_is_pentagon(
        DGG.columncell(layout, column).id)
    all(ispentagon, polar_pentagons) || error(
        "canonical polar exceptions are not topological pentagons")
    all(c -> c in columns, deep_full_land) || error(
        "canonical deep full-land block is not contained in the selection")
    all(c -> c in columns, polar_pentagons) || error(
        "canonical polar pentagons are not contained in the selection")
    return (columns = columns, sha256 = digest,
        covering = Int(canonical.production_covering),
        eligible = Int(canonical.eligible_south_of_60),
        deep_full_land = deep_full_land,
        polar_pentagons = polar_pentagons,
        latitude = Dict(c => latitude(c) for c in columns),
        eligible_latitude = (Float64(canonical.eligible_latitude_min),
            Float64(canonical.eligible_latitude_max)),
        selection_source = "canonical column_set in $RECORD; original source sidecars absent")
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
    if !isfile(PRODUCTION * ".columns.txt") || !isfile(PRODUCTION * ".done.ndjson")
        return recordedcolumnset()
    end
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
        eligible_latitude = extrema(latitude.(eligible)),
        selection_source = "production column and completion sidecars")
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
        columns_sha256 = set.sha256, columns = set.columns,
        selection_source = set.selection_source)
    println("selected $(length(set.columns)) columns; sha256=$(set.sha256)")
    println("eligible=$(set.eligible), deep-full-land=$(length(set.deep_full_land)), " *
            "polar-pentagons=$(set.polar_pentagons)")
    return set
end

"The expensive deterministic 800-column subset for t21 utilization tuning."
function utilcolumnset()
    canonical = columnset()
    length(canonical.deep_full_land) == 1_087 || error(
        "expected 1,087 canonical deep full-land columns, got " *
        "$(length(canonical.deep_full_land))")
    deep = sort(copy(canonical.deep_full_land);
        by = c -> (canonical.latitude[c], c))
    selected_deep = evensample(deep, UTIL_TARGET_COLUMNS -
        length(canonical.polar_pentagons))
    selected = sort!([selected_deep; canonical.polar_pentagons])
    length(selected) == UTIL_TARGET_COLUMNS || error(
        "utilization selector produced $(length(selected)) columns")
    allunique(selected) || error("utilization selector produced duplicate columns")
    all(c -> c in selected, canonical.polar_pentagons) || error(
        "utilization selector omitted a polar pentagon")
    all(c -> c in canonical.deep_full_land, selected_deep) || error(
        "utilization selector left the deep full-land block")
    digest = bytes2hex(SHA.sha256(join(selected, ',')))
    return (columns = selected, sha256 = digest,
        deep_full_land = selected_deep,
        canonical_deep_full_land_count = length(deep),
        polar_pentagons = canonical.polar_pentagons,
        latitude = Dict(c => canonical.latitude[c] for c in selected),
        selection_source = canonical.selection_source)
end

function recordutilselection()
    set = utilcolumnset()
    emit("util_column_set"; definition =
        "the two canonical southernmost IGeo7 pentagons plus 798 evenly " *
        "spaced (including endpoint) ranks after sorting the canonical 1,087 " *
        "production-ledger nan==0 columns below -80 degrees by " *
        "(latitude,column)",
        target = UTIL_TARGET_COLUMNS,
        canonical_deep_full_land_count = set.canonical_deep_full_land_count,
        selected_deep_full_land_count = length(set.deep_full_land),
        selected_deep_full_land_columns = set.deep_full_land,
        polar_pentagon_columns = set.polar_pentagons,
        polar_pentagon_latitudes = [set.latitude[c] for c in set.polar_pentagons],
        columns_sha256 = set.sha256, columns = set.columns,
        selection_source = set.selection_source)
    println("selected $(length(set.columns)) utilization columns; sha256=$(set.sha256)")
    println("deep-full-land=$(length(set.deep_full_land)), " *
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

function procnice(pid::Integer = getpid())
    stat = read("/proc/$pid/stat", String)
    closeparen = findlast(==(')'), stat)
    closeparen === nothing && error("malformed /proc/$pid/stat")
    fields = split(SubString(stat, closeparen + 2)) # field 3 is fields[1]
    return parse(Int, fields[17]) # /proc stat field 19
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

function gcdeltas(new::Base.GC_Num, old::Base.GC_Num)
    diff = Base.GC_Diff(new, old)
    return (time_s = diff.total_time / 1e9,
        pauses = diff.pause, full_sweeps = diff.full_sweep,
        mark_time_s = (new.total_mark_time - old.total_mark_time) / 1e9,
        sweep_time_s = (new.total_sweep_time - old.total_sweep_time) / 1e9,
        safepoint_time_s = (new.total_time_to_safepoint -
            old.total_time_to_safepoint) / 1e9,
        allocated_gib = diff.allocd / 2^30)
end

mutable struct Telemetry
    tag::String
    started::Float64
    startticks::Int64
    lasttime::Float64
    lastticks::Int64
    startgc::Base.GC_Num
    lastgc::Base.GC_Num
    stop::Threads.Atomic{Bool}
    cores::Vector{Float64}
    rss::Vector{Float64}
    hwm::Vector{Float64}
    gctime::Vector{Float64}
end

function Telemetry(tag)
    now = time()
    ticks = proccputicks()
    gc = Base.gc_num()
    Telemetry(tag, now, ticks, now, ticks, gc, gc, Threads.Atomic{Bool}(false),
        Float64[], Float64[], Float64[], Float64[])
end

function takesample!(telemetry::Telemetry; initial = false, final = false)
    now = time()
    ticks = proccputicks()
    elapsed = now - telemetry.started
    interval = now - telemetry.lasttime
    cores = interval > 0 ? (ticks - telemetry.lastticks) / CLOCK_TICKS / interval : 0.0
    gc = Base.gc_num()
    gcdiff = gcdeltas(gc, telemetry.lastgc)
    memory = procmemory()
    initial || push!(telemetry.cores, cores)
    initial || push!(telemetry.gctime, gcdiff.time_s)
    push!(telemetry.rss, memory.rss_gib)
    push!(telemetry.hwm, memory.hwm_gib)
    emit("sample"; config = telemetry.tag, elapsed_s = elapsed,
        interval_s = interval, cpu_ticks = ticks,
        instantaneous_cores = initial ? nothing : cores,
        rss_gib = memory.rss_gib, peak_rss_gib = memory.hwm_gib,
        gc_live_gib = Base.gc_live_bytes() / 2^30,
        gc_time_interval_s = initial ? nothing : gcdiff.time_s,
        gc_wall_fraction = initial || interval <= 0 ? nothing : gcdiff.time_s / interval,
        gc_pause_count = initial ? nothing : gcdiff.pauses,
        gc_full_sweeps = initial ? nothing : gcdiff.full_sweeps,
        gc_mark_time_interval_s = initial ? nothing : gcdiff.mark_time_s,
        gc_sweep_time_interval_s = initial ? nothing : gcdiff.sweep_time_s,
        gc_safepoint_time_interval_s = initial ? nothing : gcdiff.safepoint_time_s,
        gc_allocated_interval_gib = initial ? nothing : gcdiff.allocated_gib,
        initial = initial, final = final)
    telemetry.lasttime = now
    telemetry.lastticks = ticks
    telemetry.lastgc = gc
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

"Map Julia's default-pool thread ids to Linux tids before the measured run."
function juliathreadtids()
    Sys.islinux() || return Int[]
    tids = zeros(Int, Threads.nthreads())
    # The thread that enters @threads is not assigned one of its static loop
    # bodies on Julia 1.12, so record it explicitly first.
    tids[Threads.threadid()] = Int(ccall(:gettid, Cint, ()))
    Threads.@threads :static for _ in 1:Threads.nthreads()
        id = Threads.threadid()
        1 <= id <= length(tids) || continue
        tids[id] = Int(ccall(:gettid, Cint, ()))
    end
    all(>(0), tids) || error("failed to map every Julia thread to a Linux tid: $tids")
    allunique(tids) || error("Julia thread/tid map is not unique: $tids")
    return tids
end

function configureutilprofile(tag)
    tag == "baseline" || return false
    Profile.init(n = 5_000_000, delay = 0.01)
    Profile.set_peek_duration(60.0)
    Profile.peek_report[] = () -> begin
        println(stderr, "UTIL_PROFILE_BEGIN config=$tag duration_s=60 delay_s=0.01")
        Profile.print(stderr; format = :flat, sortedby = :count,
            combine = true, mincount = 20)
        println(stderr, "UTIL_PROFILE_END config=$tag")
        flush(stderr)
    end
    return true
end

function runconfig(tag::String; scaling = false, utilization = false)
    scaling && utilization && error("a run cannot be both scaling and utilization")
    nsweep = Base.JLOptions().nsweepthreads
    nsweep == 0 || error("refusing unsafe concurrent sweeper: nsweepthreads=$nsweep")
    spec = if utilization
        util_spec = get(UTIL_RUNS, tag, nothing)
        util_spec === nothing && error(
            "utilization run must be one of $(join(sort!(collect(keys(UTIL_RUNS))), ", "))")
        Threads.nthreads() == util_spec.threads || error(
            "$tag requires -t $(util_spec.threads), got $(Threads.nthreads())")
        Base.JLOptions().nmarkthreads == util_spec.mark || error(
            "$tag requires --gcthreads=$(util_spec.mark), got " *
            "$(Base.JLOptions().nmarkthreads)")
        util_spec
    elseif scaling
        scaling_spec = get(SCALING_RUNS, tag, nothing)
        scaling_spec === nothing && error(
            "scaling run must be one of $(join(sort!(collect(keys(SCALING_RUNS))), ", "))")
        Threads.nthreads() == scaling_spec.threads || error(
            "$tag requires -t $(scaling_spec.threads), got $(Threads.nthreads())")
        Base.JLOptions().nmarkthreads == 4 || error(
            "$tag requires --gcthreads=4, got $(Base.JLOptions().nmarkthreads)")
        (mark = 4, shape = scaling_spec.shape, workers = scaling_spec.workers,
            batch = CONFIG.batch, taper = true, schedule = :affinity)
    else
        Threads.nthreads() == 26 || error("shakedown runs require -t 26")
        shakedown_spec = if tag == "C"
            Base.JLOptions().nmarkthreads in (4, 8, 16) || error(
                "C must reuse a measured first-field mark count")
            (mark = Int(Base.JLOptions().nmarkthreads), shape = :inner)
        else
            get(RUN_TAGS, tag, nothing)
        end
        shakedown_spec === nothing && error("run tag must be A, B, C, or D")
        (mark = shakedown_spec.mark, shape = shakedown_spec.shape, workers = 0,
            batch = CONFIG.batch, taper = true, schedule = :affinity)
    end
    Base.JLOptions().nmarkthreads == spec.mark || error(
        "$tag requires --gcthreads=$(spec.mark), got $(Base.JLOptions().nmarkthreads)")

    set = utilization ? utilcolumnset() : columnset()
    config_name = scaling || utilization ? tag : "$tag-gc$(spec.mark)-$(spec.shape)"
    store = joinpath(SCRATCH,
        utilization ? "t21util-$config_name.zarr" :
        scaling ? "innerscale-$config_name.zarr" : "southpole-$config_name.zarr")
    mkpath(SCRATCH)
    islink(SCRATCH) && error("scratch root must not be a symlink")
    abspath(dirname(store)) == SCRATCH || error("store escaped scratch root")
    for path in (store, donelogpath(store), chunklistpath(store))
        ispath(path) && error("fresh-store requirement: $path already exists")
    end
    save_chunklist(chunklistpath(store), CONFIG.ancestor, set.columns)

    config = merge(CONFIG, (source = :synthetic, store = store, region = nothing,
        real = :none, data = DATA, workers = spec.workers, cores = 24,
        shape = spec.shape, batch = spec.batch,
        schedule = spec.schedule, cachepolicy = :refcount, taper = spec.taper,
        prefetch = 0, fetchdelay = 0.0, resume = false, checks = false,
        heartbeat = 300, maxchunks = 0, chunks = set.columns, dryrun = false))
    workers, resolved_shape = workercount(config)
    expected_workers = spec.workers > 0 ? spec.workers :
        (spec.shape === :outer ? 23 : 8)
    workers == expected_workers || error(
        "expected $expected_workers workers for $(spec.shape), resolved $workers")

    FAILURES[] = 0
    LASTCACHE[] = nothing
    profile_enabled = utilization && configureutilprofile(tag)
    thread_tids = utilization ? juliathreadtids() : Int[]
    emit("run_start"; config = config_name, campaign_tag = tag,
        pid = getpid(), julia = string(VERSION), threads = Threads.nthreads(),
        gc_mark_threads = spec.mark, gc_sweep_threads = nsweep,
        shape = String(resolved_shape), workers = workers, cores_budget = config.cores,
        workers_override = spec.workers,
        batch = config.batch, taper = config.taper, schedule = String(config.schedule),
        source = String(config.source), prefetch = config.prefetch,
        store = store, columns = length(set.columns), columns_sha256 = set.sha256,
        julia_thread_tids = thread_tids,
        profile_enabled = profile_enabled,
        profile_duration_s = profile_enabled ? 60.0 : 0.0,
        nice = Sys.islinux() ? procnice() : 10)

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
    gctotal = gcdeltas(Base.gc_num(), telemetry.startgc)
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
        gc_time_s = gctotal.time_s, gc_wall_fraction = gctotal.time_s / wall,
        gc_pause_count = gctotal.pauses, gc_full_sweeps = gctotal.full_sweeps,
        gc_mark_time_s = gctotal.mark_time_s, gc_sweep_time_s = gctotal.sweep_time_s,
        gc_safepoint_time_s = gctotal.safepoint_time_s,
        gc_allocated_gib = gctotal.allocated_gib,
        peak_rss_gib = memory.hwm_gib,
        sampled_peak_rss_gib = maximum(telemetry.rss),
        end_rss_gib = memory.rss_gib,
        cache_peak_gib = cache === nothing ? nothing : cache.peakbytes / 2^30,
        cache_peak_tiles = cache === nothing ? nothing : cache.peaktiles,
        cache_loads = cache === nothing ? nothing : cache.loads,
        cache_demanded = cache === nothing ? nothing : cache.demanded,
        cache_reloads = cache === nothing ? nothing : cache.loads - cache.demanded,
        cache_joined_loads = cache === nothing ? nothing : cache.waits,
        cache_uncredited = cache === nothing ? nothing : cache.uncredited,
        cache_live_end = cache === nothing ? nothing : cache.live,
        cache_pinned_end = cache === nothing ? nothing : cache.pinned,
        cold_downloads = 0, last10_columns_span_s = timings.last10_span,
        workers = workers, workers_override = spec.workers,
        batch = config.batch, taper = config.taper, schedule = String(config.schedule),
        shape = String(resolved_shape), gc_mark_threads = spec.mark,
        regime = regime, columns_sha256 = set.sha256, store = store)
    println("FINAL $config_name: $(round(core_s_per_column; digits=3)) core-s/column, " *
            "$(round(cpu_seconds / wall; digits=2)) cores, " *
            "$(round(memory.hwm_gib; digits=2)) GiB peak, $regime")
    expected_columns = utilization ? UTIL_TARGET_COLUMNS : TARGET_COLUMNS
    return failures == 0 && completed == expected_columns ? 0 : 1
end

"Bit-compare the scaling stores on a deterministic 12-column polar sample."
function scalingcrosscheck()
    set = columnset()
    baseline_tag = "outer-t32"
    baseline_store = joinpath(SCRATCH, "innerscale-$baseline_tag.zarr")
    isdir(baseline_store) || error("scaling baseline store is missing: $baseline_store")
    candidate_tags = [tag for tag in keys(SCALING_RUNS) if tag != baseline_tag]
    sort!(candidate_tags)
    candidate_stores = Dict(tag => joinpath(SCRATCH, "innerscale-$tag.zarr")
        for tag in candidate_tags)
    all(isdir, values(candidate_stores)) || error(
        "scaling cross-check needs every configured store: $(candidate_stores)")

    pentagons = set.polar_pentagons
    deep = sort!([c for c in set.columns if set.latitude[c] < -80.0];
        by = c -> (set.latitude[c], c))
    band = sort!([c for c in set.columns
        if set.latitude[c] >= -80.0 && !(c in pentagons)];
        by = c -> (set.latitude[c], c))
    sampled = sort!(unique([pentagons; evensample(deep, 5); evensample(band, 5)]))
    length(sampled) >= 10 || error("scaling cross-check selected only $(length(sampled)) columns")
    all(c -> c in sampled, pentagons) || error("scaling cross-check lacks a pentagon")

    baseline = Zarr.zopen(baseline_store, "r")["elevation"]
    total_differences = 0
    comparisons = 0
    all_passed = true
    for tag in candidate_tags
        candidate = Zarr.zopen(candidate_stores[tag], "r")["elevation"]
        candidate_differences = 0
        for column in sampled
            expected = Vector{Float32}(baseline[:, column])
            observed = Vector{Float32}(candidate[:, column])
            length(expected) == length(observed) || error(
                "$tag column $column length mismatch")
            differences = count(i -> reinterpret(UInt32, expected[i]) !=
                reinterpret(UInt32, observed[i]), eachindex(expected))
            emit("scaling_bitcheck_column"; baseline = baseline_tag,
                candidate = tag, column = column,
                pentagon = column in pentagons, cells = length(expected),
                non_bit_equal_float32 = differences)
            candidate_differences += differences
            comparisons += length(expected)
        end
        passed = candidate_differences == 0
        emit("scaling_bitcheck_candidate"; baseline = baseline_tag,
            candidate = tag, columns = sampled, column_count = length(sampled),
            non_bit_equal_float32 = candidate_differences, passed = passed)
        println("bitcheck $baseline_tag vs $tag: $(length(sampled)) columns, " *
                "$candidate_differences differing Float32 values")
        total_differences += candidate_differences
        all_passed &= passed
    end
    emit("scaling_bitcheck_final"; baseline = baseline_tag,
        candidates = candidate_tags, columns = sampled,
        column_count = length(sampled), pentagons = pentagons,
        compared_float32 = comparisons,
        non_bit_equal_float32 = total_differences, passed = all_passed)
    return all_passed ? 0 : 1
end

"Read one Linux stat file, robust to spaces and parentheses in comm."
function linuxstat(path)
    stat = read(path, String)
    openparen = findfirst(==('('), stat)
    closeparen = findlast(==(')'), stat)
    (openparen === nothing || closeparen === nothing) && error("malformed $path")
    fields = split(SubString(stat, closeparen + 2)) # field 3 is fields[1]
    return (comm = String(SubString(stat, openparen + 1, closeparen - 1)),
        state = String(fields[1]),
        ticks = parse(Int64, fields[12]) + parse(Int64, fields[13]),
        nice = parse(Int, fields[17]),
        processor = parse(Int, fields[37])) # /proc stat field 39
end

function procmemory(pid::Integer)
    rss_kib = -1
    hwm_kib = -1
    for line in eachline("/proc/$pid/status")
        startswith(line, "VmRSS:") && (rss_kib = parse(Int, split(line)[2]))
        startswith(line, "VmHWM:") && (hwm_kib = parse(Int, split(line)[2]))
    end
    rss_kib >= 0 || error("VmRSS missing from /proc/$pid/status")
    hwm_kib >= 0 || error("VmHWM missing from /proc/$pid/status")
    return (rss_gib = rss_kib / 2^20, hwm_gib = hwm_kib / 2^20)
end

function procthreads(pid::Integer)
    out = Dict{Int,Any}()
    taskdir = "/proc/$pid/task"
    isdir(taskdir) || return out
    for name in readdir(taskdir)
        tid = tryparse(Int, name)
        tid === nothing && continue
        stat = try
            linuxstat(joinpath(taskdir, name, "stat"))
        catch err
            (isa(err, SystemError) || isa(err, IOError)) && continue
            rethrow()
        end
        out[tid] = stat
    end
    return out
end

function recordedjuliatids(pid::Integer, tag::String)
    isfile(UTIL_RECORD) || return Dict{Int,Int}()
    found = Dict{Int,Int}()
    for line in eachline(UTIL_RECORD)
        row = try
            JSON3.read(line)
        catch
            continue
        end
        hasproperty(row, :kind) && String(row.kind) == "run_start" || continue
        hasproperty(row, :pid) && Int(row.pid) == pid || continue
        hasproperty(row, :config) && String(row.config) == tag || continue
        empty!(found)
        for (julia_thread, tid) in enumerate(Int.(row.julia_thread_tids))
            found[tid] = julia_thread
        end
    end
    return found
end

function emitprocsample(tag, pid, started, lasttime, lastproc, lastthreads, juliamap;
        initial = false, final = false)
    now = time()
    proc = linuxstat("/proc/$pid/stat")
    threads = procthreads(pid)
    memory = procmemory(pid)
    interval = now - lasttime
    process_cores = initial || interval <= 0 ? nothing :
        (proc.ticks - lastproc.ticks) / CLOCK_TICKS / interval
    rows = Any[]
    julia_busy = Float64[]
    for tid in sort!(collect(keys(threads)))
        stat = threads[tid]
        old = get(lastthreads, tid, nothing)
        busy = initial || old === nothing || interval <= 0 ? nothing :
            (stat.ticks - old.ticks) / CLOCK_TICKS / interval
        julia_thread = get(juliamap, tid, nothing)
        julia_thread === nothing || busy === nothing || push!(julia_busy, busy)
        push!(rows, (tid = tid, julia_thread = julia_thread,
            comm = stat.comm, state = stat.state, processor = stat.processor,
            cpu_ticks = stat.ticks, busy_fraction = busy))
    end
    emit("proc_thread_sample"; config = tag, child_pid = pid,
        elapsed_s = now - started, interval_s = interval,
        process_cpu_ticks = proc.ticks, instantaneous_cores = process_cores,
        rss_gib = memory.rss_gib, peak_rss_gib = memory.hwm_gib,
        task_count = length(threads),
        runnable_task_count = count(s -> s.state == "R", values(threads)),
        mapped_julia_threads = length(julia_busy),
        julia_busy_p10 = isempty(julia_busy) ? nothing : Statistics.quantile(julia_busy, 0.1),
        julia_busy_median = isempty(julia_busy) ? nothing : Statistics.median(julia_busy),
        julia_busy_p90 = isempty(julia_busy) ? nothing : Statistics.quantile(julia_busy, 0.9),
        julia_busy_lt10_count = count(<(0.1), julia_busy),
        julia_busy_gt80_count = count(>=(0.8), julia_busy),
        threads = rows, initial = initial, final = final)
    return (time = now, proc = proc, threads = threads, memory = memory,
        process_cores = process_cores, julia_busy = julia_busy)
end

function stopchild!(process; grace_seconds = 300.0)
    process_running(process) || return
    Base.kill(process, Base.SIGINT)
    deadline = time() + grace_seconds
    while process_running(process) && time() < deadline
        sleep(1.0)
    end
    process_running(process) || return
    Base.kill(process, Base.SIGTERM)
    deadline = time() + 10.0
    while process_running(process) && time() < deadline
        sleep(0.25)
    end
    process_running(process) && Base.kill(process, Base.SIGKILL)
    return
end

"Supervise one t21 child and sample it externally from /proc every three seconds."
function superviseutilrun(tag::String)
    Threads.nthreads() == 1 || error("util-run supervisor requires -t1")
    Base.JLOptions().nmarkthreads == 1 || error(
        "util-run supervisor requires --gcthreads=1")
    Base.JLOptions().nsweepthreads == 0 || error("refusing unsafe concurrent sweeper")
    spec = get(UTIL_RUNS, tag, nothing)
    spec === nothing && error(
        "util-run must be one of $(join(sort!(collect(keys(UTIL_RUNS))), ", "))")
    priority = procnice()
    priority >= 10 || error("util-run must be launched under nice -n 10 (got $priority)")
    mkpath(SCRATCH)
    logpath = joinpath(SCRATCH, "t21util-$tag.log")
    ispath(logpath) && error("fresh-log requirement: $logpath already exists")
    project = Base.active_project()
    project === nothing && error("util-run requires --project=benchmark")
    cmd = `$(Base.julia_cmd()) --project=$(dirname(project)) -t 21 --gcthreads=$(spec.mark) $(@__FILE__) util-child $tag`
    started = time()
    process = nothing
    logio = open(logpath, "w")
    try
        process = run(pipeline(cmd, stdout = logio, stderr = logio); wait = false)
        pid = getpid(process)
        emit("proc_supervisor_start"; config = tag, child_pid = pid,
            command = string(cmd), sample_seconds = UTIL_PROC_SECONDS,
            max_wall_s = UTIL_MAX_WALL_SECONDS, max_rss_gib = UTIL_MAX_RSS_GIB,
            log = logpath, nice = priority)
        lasttime = time()
        lastproc = linuxstat("/proc/$pid/stat")
        lastthreads = procthreads(pid)
        juliamap = Dict{Int,Int}()
        firstsample = emitprocsample(tag, pid, started, lasttime, lastproc,
            lastthreads, juliamap; initial = true)
        lasttime, lastproc, lastthreads =
            firstsample.time, firstsample.proc, firstsample.threads
        profile_sent = false
        abort_reason = nothing
        max_hwm = firstsample.memory.hwm_gib
        last_progress = started
        while process_running(process)
            sleep(UTIL_PROC_SECONDS)
            process_running(process) || break
            isempty(juliamap) && (juliamap = recordedjuliatids(pid, tag))
            sample = try
                emitprocsample(tag, pid, started, lasttime, lastproc,
                    lastthreads, juliamap)
            catch err
                (isa(err, SystemError) || isa(err, IOError)) && !process_running(process) && break
                rethrow()
            end
            lasttime, lastproc, lastthreads = sample.time, sample.proc, sample.threads
            max_hwm = max(max_hwm, sample.memory.hwm_gib)
            elapsed = time() - started
            if tag == "baseline" && !profile_sent && !isempty(juliamap) && elapsed >= 180.0
                # Julia 1.12 exposes no Base.SIGUSR1 constant. This supervisor
                # is Linux-only (/proc is its data source), where SIGUSR1 is 10.
                ccall(:kill, Cint, (Cint, Cint), Cint(pid), Cint(10)) == 0 ||
                    error("kill(SIGUSR1) failed for child $pid")
                profile_sent = true
                emit("profile_signal"; config = tag, child_pid = pid,
                    elapsed_s = elapsed, duration_s = 60.0)
                println("PROFILE $tag: sent 60-second SIGUSR1 capture at $(round(elapsed; digits=1)) s")
                flush(stdout)
            end
            if sample.memory.hwm_gib > UTIL_MAX_RSS_GIB
                abort_reason = "peak RSS $(sample.memory.hwm_gib) GiB exceeded $(UTIL_MAX_RSS_GIB) GiB"
                break
            elseif elapsed > UTIL_MAX_WALL_SECONDS
                abort_reason = "wall $(elapsed) s exceeded $(UTIL_MAX_WALL_SECONDS) s"
                break
            end
            if time() - last_progress >= 30.0
                busy = sample.julia_busy
                med = isempty(busy) ? NaN : Statistics.median(busy)
                println("PROC $tag: $(round(elapsed; digits=0)) s, " *
                    "$(round(something(sample.process_cores, NaN); digits=2)) cores, " *
                    "Julia-thread median=$(round(med; digits=2)), " *
                    "RSS=$(round(sample.memory.rss_gib; digits=2)) GiB " *
                    "(peak $(round(max_hwm; digits=2)))")
                flush(stdout)
                last_progress = time()
            end
        end
        if abort_reason !== nothing
            emit("proc_guard_abort"; config = tag, child_pid = pid,
                reason = abort_reason, elapsed_s = time() - started, peak_rss_gib = max_hwm)
            println("ABORT $tag: $abort_reason")
            flush(stdout)
            stopchild!(process)
        end
        wait(process)
        emit("proc_supervisor_final"; config = tag, child_pid = pid,
            elapsed_s = time() - started, exitcode = process.exitcode,
            profile_sent = profile_sent, aborted = abort_reason !== nothing,
            abort_reason = abort_reason, peak_rss_gib = max_hwm)
        return success(process) && abort_reason === nothing ? 0 : 1
    finally
        process !== nothing && process_running(process) && stopchild!(process)
        close(logio)
    end
end

"Bit-compare completed utilization stores on six deterministic polar columns."
function utilcrosscheck()
    set = utilcolumnset()
    baseline_tag = "baseline"
    baseline_store = joinpath(SCRATCH, "t21util-$baseline_tag.zarr")
    isdir(baseline_store) || error("utilization baseline store is missing: $baseline_store")
    candidate_tags = sort!([tag for tag in keys(UTIL_RUNS)
        if tag != baseline_tag && isdir(joinpath(SCRATCH, "t21util-$tag.zarr")) &&
        length(ledgerrows(donelogpath(joinpath(SCRATCH, "t21util-$tag.zarr")))) ==
            UTIL_TARGET_COLUMNS])
    isempty(candidate_tags) && error("no completed utilization variant store exists")
    sampled = sort!(unique([set.polar_pentagons;
        evensample(sort(set.deep_full_land; by = c -> (set.latitude[c], c)), 4)]))
    length(sampled) >= 5 || error("utilization bitcheck selected only $(length(sampled)) columns")
    any(c -> c in set.polar_pentagons, sampled) || error("utilization bitcheck lacks a pentagon")
    baseline = Zarr.zopen(baseline_store, "r")["elevation"]
    total_differences = 0
    comparisons = 0
    all_passed = true
    for tag in candidate_tags
        candidate = Zarr.zopen(joinpath(SCRATCH, "t21util-$tag.zarr"), "r")["elevation"]
        differences = 0
        for column in sampled
            expected = Vector{Float32}(baseline[:, column])
            observed = Vector{Float32}(candidate[:, column])
            length(expected) == length(observed) || error("$tag column $column length mismatch")
            column_differences = count(i -> reinterpret(UInt32, expected[i]) !=
                reinterpret(UInt32, observed[i]), eachindex(expected))
            emit("util_bitcheck_column"; baseline = baseline_tag, candidate = tag,
                column = column, pentagon = column in set.polar_pentagons,
                cells = length(expected), non_bit_equal_float32 = column_differences)
            differences += column_differences
            comparisons += length(expected)
        end
        passed = differences == 0
        emit("util_bitcheck_candidate"; baseline = baseline_tag, candidate = tag,
            columns = sampled, column_count = length(sampled),
            non_bit_equal_float32 = differences, passed = passed)
        println("bitcheck $baseline_tag vs $tag: $(length(sampled)) columns, " *
                "$differences differing Float32 values")
        total_differences += differences
        all_passed &= passed
    end
    emit("util_bitcheck_final"; baseline = baseline_tag, candidates = candidate_tags,
        columns = sampled, column_count = length(sampled),
        pentagons = set.polar_pentagons, compared_float32 = comparisons,
        non_bit_equal_float32 = total_differences, passed = all_passed)
    return all_passed ? 0 : 1
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
        "usage: southpole_shakedown.jl select | run A|B|C|D | crosscheck | " *
        "diagnose | scale-select | scale-run CONFIG | scale-crosscheck | " *
        "util-select | util-run CONFIG | util-child CONFIG | util-crosscheck")
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
    elseif ARGS[1] == "scale-select"
        length(ARGS) == 1 || error("scale-select takes no arguments")
        ACTIVE_RECORD[] = SCALING_RECORD
        recordselection()
        return 0
    elseif ARGS[1] == "scale-run"
        length(ARGS) == 2 || error("scale-run requires one config name")
        ACTIVE_RECORD[] = SCALING_RECORD
        return runconfig(ARGS[2]; scaling = true)
    elseif ARGS[1] == "scale-crosscheck"
        length(ARGS) == 1 || error("scale-crosscheck takes no arguments")
        ACTIVE_RECORD[] = SCALING_RECORD
        return scalingcrosscheck()
    elseif ARGS[1] == "util-select"
        length(ARGS) == 1 || error("util-select takes no arguments")
        ACTIVE_RECORD[] = UTIL_RECORD
        recordutilselection()
        return 0
    elseif ARGS[1] == "util-run"
        length(ARGS) == 2 || error("util-run requires one config name")
        ACTIVE_RECORD[] = UTIL_RECORD
        return superviseutilrun(ARGS[2])
    elseif ARGS[1] == "util-child"
        length(ARGS) == 2 || error("util-child requires one config name")
        ACTIVE_RECORD[] = UTIL_RECORD
        return runconfig(ARGS[2]; utilization = true)
    elseif ARGS[1] == "util-crosscheck"
        length(ARGS) == 1 || error("util-crosscheck takes no arguments")
        ACTIVE_RECORD[] = UTIL_RECORD
        return utilcrosscheck()
    end
    error("unknown mode $(repr(ARGS[1]))")
end

exit(main_shakedown())
