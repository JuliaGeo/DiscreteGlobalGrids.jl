# The store side of the run: the Zarr store, its two sidecar files, and what
# "already done" means. Included by `copdem_production.jl`.
#
# A CHUNK is one level-`ancestor` IGeo7 cell together with all its
# level-`level` descendants: one work unit, one Zarr chunk, one file on disk.
# The DGG store API calls the same thing a "column", so `DGG.columncell`,
# `DGG.columnindex` and `DGG.columnlength` keep that name.

# ---------------------------------------------------------------------------
# Sidecar files
# ---------------------------------------------------------------------------

"The ledger of finished chunks, appended as the run goes."
donelogpath(store) = store * ".done.ndjson"

"""
The cached covering, beside the store.

The `.columns.txt` name predates the column -> chunk rename and is kept as it
is, so a store written by an earlier run still finds its cache and resumes.
"""
chunklistpath(store) = store * ".columns.txt"

# ---------------------------------------------------------------------------
# The chunk list
# ---------------------------------------------------------------------------

"""
    load_chunklist(path) -> Vector{Int} or nothing

The cached covering, or `nothing` if there is none. Cached because computing it
is ~26 000 coverage queries and it comes out the same on every resume.
"""
function load_chunklist(path)
    isfile(path) || return nothing
    out = Int[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        push!(out, parse(Int, s))
    end
    return out
end

"Write the covering to `path`, one chunk index per line."
function save_chunklist(path, ancestor, chunks)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# level-$(ancestor) column indices covering the listed tiles")
        for c in chunks
            println(io, c)
        end
    end
    return path
end

# ---------------------------------------------------------------------------
# The done ledger
# ---------------------------------------------------------------------------

"""
    DoneLog(path)

Append-only ndjson, one line per finished chunk, flushed as written, so a hard
kill loses at most the line in flight.

The chunk index goes under the key `"col"`: that name predates the
column -> chunk rename, and [`donechunks`](@ref) still reads it, so existing
ledgers stay valid.
"""
struct DoneLog
    io::IOStream
    lock::ReentrantLock
end

DoneLog(path::AbstractString) = DoneLog(open(path, "a"), ReentrantLock())
Base.close(l::DoneLog) = close(l.io)

function record!(l::DoneLog, chunk, ncells, nnan, secs, worker)
    lock(l.lock) do
        println(l.io, @sprintf("{\"col\":%d,\"cells\":%d,\"nan\":%d,\"secs\":%.2f,\"w\":%d,\"t\":\"%s\"}",
            chunk, ncells, nnan, secs, worker, stamp()))
        flush(l.io)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Resume
# ---------------------------------------------------------------------------

"""
    doneledger(logpath) -> Dict{Int,@NamedTuple{cells::Int, nan::Int}}

Every chunk the ledger has a line for, mapped to the counts that line recorded.

One entry per CHUNK, not per line: a chunk recomputed after a resume — because a
crash lost the file but not the line, or because the ledger was truncated — has
two lines, and the last one is the one whose values are in the store, so it
wins. Counting lines instead would double a chunk's cells.

Both [`donechunks`](@ref) and the store-wide totals of the final banner read
this, so "which chunks are done" and "how many cells that is" cannot disagree.
A field the line does not carry comes back as `-1`, which
[`storetotals`](@ref) reports as unaccounted rather than as zero.
"""
function doneledger(logpath)
    out = Dict{Int,@NamedTuple{cells::Int, nan::Int}}()
    isfile(logpath) || return out
    field(line, key) = (m = match(Regex("\"$key\":(\\d+)"), line);
                        m === nothing ? -1 : parse(Int, m[1]))
    for line in eachline(logpath)
        col = field(line, "col")
        col < 0 && continue
        out[col] = (cells = field(line, "cells"), nan = field(line, "nan"))
    end
    return out
end

"""
    storetotals(logpath, written) -> (; chunks, cells, nan, unaccounted)

What the STORE holds, across every session that ever wrote to it — as against
the [`Progress`](@ref) counters, which only ever describe the session printing
them.

`written` is the resume union from [`donechunks`](@ref), so `chunks` includes
chunks this run never touched, and an all-`NaN` chunk — which Zarr stores as no
file at all — still counts if the ledger saw it. Cells and NaNs can only be
counted for chunks the ledger has a line for; `unaccounted` is how many of
`written` it does not, so totals read off a truncated ledger are never mistaken
for the whole store.
"""
function storetotals(logpath, written)
    led = doneledger(logpath)
    cells = 0
    nan = 0
    unaccounted = 0
    for ch in written
        r = get(led, ch, nothing)
        if r === nothing || r.cells < 0 || r.nan < 0
            unaccounted += 1
        else
            cells += r.cells
            nan += r.nan
        end
    end
    return (chunks = length(written), cells = cells, nan = nan,
            unaccounted = unaccounted)
end

"""
    donechunks(logpath, storepath, layer; label = "resume") -> Set{Int}

Which chunks are already written: the done ledger UNION the chunk files the
store itself holds.

Neither source is complete alone. The ledger is what this script did and can be
deleted or truncated; the file listing is what Zarr HAS, and Zarr writes no file
for a chunk whose every value is the fill value — so an all-`NaN` chunk, a
normal outcome on the ocean side of the covering, leaves nothing behind and looks
exactly like one nobody computed. Skipping the union is therefore the safe
choice: recomputing a written chunk only wastes time, skipping an unwritten one
leaves a hole.

`label` names the caller in the notes this prints; the final banner reads the
same union the resume path did, and calling both of them "resume" would read as
if the run were starting over.
"""
function donechunks(logpath, storepath, layer; label = "resume")
    fromlog = Set(keys(doneledger(logpath)))
    fromdisk = storechunks(storepath, layer)
    if fromdisk !== nothing
        only_log = length(setdiff(fromlog, fromdisk))
        only_disk = length(setdiff(fromdisk, fromlog))
        (only_log == 0 && only_disk == 0) ||
            say("$label: $only_log logged chunks have no file (an all-NaN chunk is " *
                "stored as nothing at all), $only_disk files have no ledger line; " *
                "taking the union")
        return union(fromlog, fromdisk)
    end
    say("$label: no chunk listing available at $storepath, using the ledger alone")
    return fromlog
end

"""
    storechunks(path, layer) -> Set{Int} or nothing

The chunks a directory store already has a file for, or `nothing` if the listing
does not parse — in which case the caller falls back to the ledger.

Zarr v2 names a chunk by its indices in ZARR order, the reverse of Julia's, and
this layout is `(capacity, nchunks)` in Julia: a chunk's file is
`"<chunk-1>.0"`, leading field first.
"""
function storechunks(path, layer)
    dir = joinpath(path, layer)
    isdir(dir) || return nothing
    out = Set{Int}()
    for f in readdir(dir)
        startswith(f, ".") && continue
        m = match(r"^(\d+)\.(\d+)$", f)
        m === nothing && return nothing
        parse(Int, m[2]) == 0 || return nothing
        push!(out, parse(Int, m[1]) + 1)
    end
    return out
end

# ---------------------------------------------------------------------------
# The store itself
# ---------------------------------------------------------------------------

"""
    openstore(config, sys7, capacity; geometry_tag=nothing) -> store

Reopen the subzone store named by `config.store`, or create it if it is not
there. One `Float32` layer, `elevation`, filled with `NaN`; one Zarr chunk per
level-`config.ancestor` chunk, `capacity` cells wide.

A new store records `config.method` as `regridding_method`, so a file says which
semantics produced its elevations — area means, post samples, or the nearest
post. Reopening one written under a different method is refused: one store holds
one reading of the data. A store from before the attribute existed carries none
and is accepted as it is.
"""
function openstore(config, sys7, capacity;
        geometry_tag::Union{Nothing,String} = nothing)
    path = config.store
    method_tag = String(config.method)
    if isdir(path)
        store = DGG.subzonestore(path)
        observed = get(store.group.attrs, "destination_geometry", nothing)
        observed == geometry_tag || error(
            "store destination_geometry is $(repr(observed)), expected " *
            "$(repr(geometry_tag)); refusing to mix spherical and authalic columns")
        seen = get(store.group.attrs, "regridding_method", nothing)
        seen === nothing || seen == method_tag || error(
            "store regridding_method is $(repr(seen)), expected " *
            "$(repr(method_tag)); refusing to mix two readings of the data in " *
            "one store")
        say("store: reopened $path")
        return store
    end
    mkpath(dirname(path))
    t0 = time()
    sourcelabel = config.source === :real ? "REAL" : "SYNTHETIC"
    sourcedescription = config.source === :real ?
        "Copernicus DEM GLO-$(config.res) COGs from the public AWS Open Data bucket" :
        "synthetic analytic field over the real Copernicus GLO-$(config.res) tile list"
    producerattrs = Dict{String,Any}(
        "title" => "Copernicus DEM GLO-$(config.res) ($sourcelabel) on IGEO7 level $(config.level)",
        "source" => sourcedescription,
        "regridding_method" => method_tag,
        "created" => stamp())
    geometry_tag === nothing ||
        (producerattrs["destination_geometry"] = geometry_tag)
    store = DGG.subzonestore(path, sys7, config.level;
        ancestor_level = config.ancestor,
        layers = ("elevation" => Float32,), capacity = capacity,
        fill_value = NaN, ancestor_coordinate = true,
        attrs = producerattrs)
    say("store: created $path, $(DGG.ncells(sys7, config.ancestor)) chunks of " *
        "$capacity, $(secs(time() - t0))")
    return store
end
