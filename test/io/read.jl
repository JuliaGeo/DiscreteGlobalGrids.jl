# `dggread`: opening a store, describing it, and assembling the lazy DimStack.
#
# Fixtures are real Zarr v2 stores written here with Zarr.jl and hand-stamped
# with each convention's attrs, so what is under test is the whole path from
# bytes to cube and not a snapshot someone typed. Every store is a few hundred
# cells of IGEO7 level 4 or the whole of HEALPix level 2.
#
# The suite self-skips when Zarr.jl is absent: `dggread` is then the stub.

module DGGSIOReadTests

using Test
import DiscreteGlobalGrids as DGG

const HAS_ZARR = try
    @eval using Zarr
    true
catch err
    @warn "Zarr.jl is not loadable: the dggread suite is skipped." exception = err
    false
end

if HAS_ZARR

import DimensionalData as DD
using DimensionalData: Lookups
using DiscreteGlobalGrids: IGeo7System, HEALPixSystem, Z7Cell, LevelIndex,
    levelgrid, ncells, cellindex, Cells, CellVector, CellLookup, DGGSFormatError,
    dggread, dggwrite, StoreDescription
using DiscreteGlobalGrids.Encodings: idselect, idranges, DenseEncoding,
    RangesEncoding, ImplicitEncoding

const ZarrExt = Base.get_extension(DGG, :DiscreteGlobalGridsZarrExt)

d(pairs...) = Dict{String,Any}(pairs...)

# ---------------------------------------------------------------------------
# A store that counts the chunks read from it
# ---------------------------------------------------------------------------

# What it proves is laziness: which arrays `dggread` touches and which it does
# not. Metadata keys are not chunks and are not counted; every other key is a
# chunk of the array whose name prefixes it.
struct CountingStore{S<:Zarr.AbstractStore} <: Zarr.AbstractStore
    parent::S
    reads::Dict{String,Int}
end
CountingStore(p::Zarr.AbstractStore) = CountingStore(p, Dict{String,Int}())

const METADATA_KEYS = (".zarray", ".zattrs", ".zgroup", ".zmetadata")

function Base.getindex(s::CountingStore, k::String)
    if !any(m -> endswith(k, m), METADATA_KEYS)
        name = first(split(k, '/'))
        s.reads[name] = get(s.reads, name, 0) + 1
    end
    return s.parent[k]
end
Base.setindex!(s::CountingStore, v, k::String) = (s.parent[k] = v)
Base.delete!(s::CountingStore, k::String) = delete!(s.parent, k)
Zarr.subdirs(s::CountingStore, p) = Zarr.subdirs(s.parent, p)
Zarr.subkeys(s::CountingStore, p) = Zarr.subkeys(s.parent, p)
Zarr.isinitialized(s::CountingStore, k::AbstractString) = Zarr.isinitialized(s.parent, k)
Zarr.storagesize(s::CountingStore, p) = Zarr.storagesize(s.parent, p)

counting(path) = Zarr.zopen(CountingStore(Zarr.DirectoryStore(path)), "r")
reads(g::Zarr.ZGroup, name) = get(g.storage.reads, name, 0)
resetreads!(g::Zarr.ZGroup) = empty!(g.storage.reads)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

const GRID = levelgrid(IGeo7System(), 4)
# Three index runs: the shape of a regional store, and enough runs that the
# ranges twin has more than one row.
const RANKS = [50:149; 399:399; 699:848]
const IDS = UInt64[idselect(GRID, r) for r in RANKS]
const CELLS = [Z7Cell(x) for x in IDS]
const N = length(IDS)
const ELEVATION = Float32.(10 .* (1:N))
const SLOPE = Float64.(inv.(1:N))
const COORD_CHUNK = 64

# Base cell 0's pentagon chain deletes digit 2, so `00002` is well formed and
# names no cell. It sits between the level-4 cells of ranks 1 and 2.
const PHANTOM = (UInt64(2) << 48) | ((UInt64(1) << 48) - 1)

# The same axis shape in base cell 8, whose ids set the top bit: written into an
# `Int64` sidecar they are NEGATIVE integers that wrap back to exactly these
# ids, which is the one case where reinterpreting a foreign width would silently
# work and is why the width guard refuses to.
const HIGH_IDS = UInt64[idselect(GRID, r) for r in 16008:(16008+N-1)]

zarr_conventions() = Any[deepcopy(DGG.ZARR_DGGS_DECLARATION)]

function dggs_attrs(; level=4, compression="none", coordinate="cell_ids")
    dggs = d("name" => "igeo7", "refinement_level" => level,
        "spatial_dimension" => "cell_ids", "dggs_vert0_lon" => 11.2)
    coordinate === nothing || (dggs["coordinate"] = coordinate;
    dggs["compression"] = compression)
    return d("zarr_conventions" => zarr_conventions(), "dggs" => dggs)
end

xdggs_attrs(; level=4) = d("grid_name" => "igeo7", "level" => level,
    "igeo7_dggs_vert0_lon" => 11.2)

adims(names...) = d("_ARRAY_DIMENSIONS" => Any[names...])

# The data variables every IGEO7 fixture carries, chunked unlike the coordinate
# so that a store whose arrays disagree about chunking is the ordinary case.
function write_data!(g, ids=IDS)
    n = length(ids)
    e = Zarr.zcreate(Float32, g, "elevation", n; chunks=(50,), attrs=adims("cell_ids"))
    e[:] = Float32.(10 .* (1:n))
    s = Zarr.zcreate(Float64, g, "slope", n; chunks=(n,),
        attrs=merge(adims("cell_ids"), d("units" => "m/m")))
    s[:] = Float64.(inv.(1:n))
    return g
end

"""
A dense store: one stored id per cell, `zarr-conventions/dggs` only. `T` is the
integer width the coordinate itself is written at, which a store chooses and a
reader does not: `Int64` holds the same ids as `UInt64` with the top bit read as
a sign.
"""
function dense_store(dir; name="dense.zarr", group=dggs_attrs(), coord=adims("cell_ids"),
    ids=IDS, T=UInt64)
    path = joinpath(dir, name)
    g = Zarr.zgroup(path; attrs=group)
    c = Zarr.zcreate(T, g, "cell_ids", length(ids); chunks=(COORD_CHUNK,), attrs=coord)
    c[:] = map(x -> x % T, ids)
    write_data!(g, ids)
    return path
end

"The same cells as `dense_store`, as an `(n, 2)` inclusive-range coordinate."
function ranges_store(dir; name="ranges.zarr", ids=IDS, T=UInt64)
    path = joinpath(dir, name)
    g = Zarr.zgroup(path; attrs=dggs_attrs(compression="ranges",
        coordinate="cell_id_ranges"))
    rows = idranges(GRID, ids)
    # Zarr declares shapes outermost-first, Julia innermost-first: an `(m, 2)`
    # array is `(2, m)` here.
    c = Zarr.zcreate(T, g, "cell_id_ranges", 2, size(rows, 1);
        chunks=(2, size(rows, 1)), attrs=adims("ranges", "bounds"))
    c[:, :] = map(x -> x % T, permutedims(rows))
    write_data!(g, ids)
    return path
end

"The cube the dense fixtures hold, over this package's own compressed lookup."
cube() = DD.DimStack((elevation=copy(ELEVATION), slope=copy(SLOPE)),
    (Cells(CellLookup(CellVector(IGeo7System(), 4, CELLS))),))

"""
A dense store carrying a chunk-manifest sidecar whose MARKER the caller shapes:
`marker` overrides fields and `drop` removes them. The rows are otherwise the
truth about the ids, so what a fixture varies is trust alone; `rows` is the one
exception and rewrites the table itself — as Julia holds it, `(2, n_chunks)` —
and `T` the integer width it is written at.
"""
function manifest_store(dir; name, ids=IDS, chunk=COORD_CHUNK,
    marker=d(), drop=String[], rows=identity, T=UInt64)
    path = dense_store(dir; name=name, ids=ids)
    n = length(ids)
    nc = cld(n, chunk)
    table = rows(permutedims(hcat(UInt64[ids[(c-1)*chunk+1] for c in 1:nc],
        UInt64[ids[min(c * chunk, n)] for c in 1:nc])))
    data = map(x -> x % T, table)
    m = merge(d("writer" => "DiscreteGlobalGrids.jl", "format" => 1,
            "validated" => "strict", "spatial_dimension" => "cell_ids",
            "chunk_length" => chunk, "length" => n,
            "level" => 4, "grid" => "igeo7"), marker)
    for k in drop
        delete!(m, k)
    end
    g = Zarr.zopen(path, "w")
    z = Zarr.zcreate(T, g, "cell_chunk_manifest", size(data)...; chunks=size(data),
        attrs=merge(adims("chunks", "bounds"), d("dggs_chunk_manifest" => m)))
    z[:, :] = data
    return path
end

"Chunk `c`'s first and last id, the wrong way round."
function swaprow(R, c)
    S = copy(R)
    S[1, c], S[2, c] = S[2, c], S[1, c]
    return S
end

"A nextGEMS-style store: a `crs` variable, no cell array, index is the id."
function dkrz_store(dir; name="dkrz.zarr")
    path = joinpath(dir, name)
    grid = levelgrid(HEALPixSystem(), 2)
    n = ncells(grid)
    g = Zarr.zgroup(path; attrs=d())
    crs = Zarr.zcreate(Int8, g, "crs", 1; chunks=(1,),
        attrs=merge(adims("crs"), d("grid_mapping_name" => "healpix",
            "healpix_nside" => 4, "healpix_order" => "nest")))
    crs[:] = Int8[0]
    tas = Zarr.zcreate(Float32, g, "tas", n, 3; chunks=(64, 1),
        attrs=merge(adims("time", "cell"), d("grid_mapping" => "crs")))
    tas[:, :] = Float32.(reshape(1:3n, n, 3))
    t = Zarr.zcreate(Int64, g, "time", 3; chunks=(3,), attrs=adims("time"))
    t[:] = Int64[2020, 2021, 2022]
    return path
end

# ---------------------------------------------------------------------------
# The dense path, end to end
# ---------------------------------------------------------------------------

@testset "a dense store reads back as the cube it was written from" begin
    mktempdir() do dir
        st = dggread(dense_store(dir))
        @test st isa DD.DimStack
        @test collect(keys(st)) == [:elevation, :slope]

        lk = DD.lookup(st[:elevation], Cells)
        @test lk isa DGG.ChunkedCellLookup
        @test length(lk) == N
        @test collect(lk) == CELLS

        # Values, addressed by cell rather than by index: the axis and the
        # data have to agree for this to hold anywhere but index 1.
        @test st[:elevation][Cells(DD.At(CELLS[137]))] == ELEVATION[137]
        @test st[:elevation][Cells(DD.At(CELLS[end]))] == ELEVATION[end]
        @test st[:slope][Cells(DD.At(CELLS[137]))] == SLOPE[137]
        @test collect(parent(st[:elevation])) == ELEVATION

        md = DD.metadata(st)
        @test md["encoding"] == "none"
        @test md["conventions"] == ["zarr-conventions/dggs"]
        @test occursin("dense.zarr", md["source"])
        # The group attrs verbatim, which is what regenerating the store needs.
        @test md["attrs"]["dggs"]["refinement_level"] == 4
        @test md["description"].level == 4
        # Array attrs ride on the layer they came from.
        @test DD.metadata(st[:slope])["units"] == "m/m"
    end
end

@testset "the data is not read until it is indexed" begin
    mktempdir() do dir
        path = dense_store(dir)
        g = counting(path)
        st = dggread(g)
        # The scan reads the coordinate, once per chunk, and nothing else.
        @test reads(g, "cell_ids") == cld(N, COORD_CHUNK)
        @test reads(g, "elevation") == 0
        @test reads(g, "slope") == 0
        @test parent(st[:elevation]) isa Zarr.ZArray

        resetreads!(g)
        @test st[:elevation][Cells(DD.At(CELLS[137]))] == ELEVATION[137]
        @test reads(g, "elevation") == 1
        @test reads(g, "slope") == 0
        # The axis inherited the store's own chunk grid, so resolving the cell
        # reads the one chunk the manifest names it in and no other.
        @test reads(g, "cell_ids") == 1

        # `lazy = false` is the opposite promise.
        g2 = counting(path)
        st2 = dggread(g2; lazy=false)
        @test parent(st2[:elevation]) isa Array
        @test reads(g2, "elevation") > 0
        @test collect(parent(st2[:slope])) == SLOPE
    end
end

# ---------------------------------------------------------------------------
# The persisted manifest
# ---------------------------------------------------------------------------

@testset "a store we wrote opens from its manifest instead of scanning" begin
    # The design's trust model: our own marker means the chunk grid on disk IS
    # the chunk grid, so opening a store of tens of millions of cells reads no
    # ids at all. Kills a reader that scans anyway, and a manifest whose rows
    # are read but whose chunk boundaries are not the store's.
    mktempdir() do dir
        path = joinpath(dir, "written.zarr")
        dggwrite(path, cube(); encoding=:dense, chunks=COORD_CHUNK)

        g = counting(path)
        st = dggread(g)
        @test reads(g, "cell_ids") == 0
        @test length(DD.lookup(st[:elevation], Cells)) == N

        resetreads!(g)
        @test st[:elevation][Cells(DD.At(CELLS[137]))] == ELEVATION[137]
        # Exactly the one chunk the manifest named — and, being decoded, the one
        # chunk the spot check was able to verify.
        @test reads(g, "cell_ids") == 1
        @test collect(DD.lookup(st[:elevation], Cells)) == CELLS
    end
end

@testset "a stale manifest is caught on the chunk it lies about" begin
    # A sidecar is only as true as the ids under it. Another writer rewrote them
    # and left the manifest behind: the reader trusts it at open, because that
    # is the whole point, and refuses on the first chunk it actually decodes.
    # Kills trust without the spot check, and a check that fires at open.
    mktempdir() do dir
        path = joinpath(dir, "stale.zarr")
        dggwrite(path, cube(); encoding=:dense, chunks=COORD_CHUNK)
        z = Zarr.zopen(path, "w")["cell_ids"]
        z[:] = UInt64[idselect(GRID, r) for r in 2000:(2000+N-1)]

        g = counting(path)
        st = dggread(g)
        @test reads(g, "cell_ids") == 0

        err = try
            st[:elevation][Cells(DD.At(CELLS[137]))]
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :stale_manifest
        @test occursin("stale.zarr", err.store)
    end
end

@testset "a manifest is trusted only as far as its marker says" begin
    # Everything the marker has to agree about before a scan is skipped. Each
    # row is a store whose manifest ROWS are correct and whose marker is not
    # ours to trust: the manifest is then an extra array and nothing else, so
    # the reader scans and gets the same axis. Kills a trust check that reads
    # the rows and forgets to read the marker.
    mktempdir() do dir
        nc = cld(N, COORD_CHUNK)
        # Trusted, as the reference the fallbacks are compared against.
        good = counting(manifest_store(dir; name="good.zarr"))
        trusted = dggread(good)
        @test reads(good, "cell_ids") == 0

        cases = (("foreign", d("writer" => "xarray-dggs"), String[]),
            ("format", d("format" => 2), String[]),
            ("sampled", d("validated" => "lazy"), String[]),
            ("unvalidated", d(), ["validated"]),
            ("otherdim", d("spatial_dimension" => "cells"), String[]),
            ("regridded", d("chunk_length" => 32), String[]),
            ("shortened", d("length" => N - 1), String[]),
            ("otherlevel", d("level" => 3), String[]),
            ("othergrid", d("grid" => "healpix"), String[]))
        for (name, marker, drop) in cases
            g = counting(manifest_store(dir; name="$name.zarr", marker, drop))
            st = dggread(g)
            @test reads(g, "cell_ids") == nc
            @test DD.lookup(st[:elevation], Cells) ==
                  DD.lookup(trusted[:elevation], Cells)
        end
        # The rows are read with the same suspicion. The two-level search binary
        # searches them, so chunk intervals that do not ascend and stay disjoint
        # would resolve an id into the wrong chunk — or into none — rather than
        # fail; this reader scans instead.
        g = counting(manifest_store(dir;
            name="descending.zarr", rows=R -> R[:, end:-1:1]))
        st = dggread(g)
        @test reads(g, "cell_ids") == nc
        @test DD.lookup(st[:elevation], Cells) == DD.lookup(trusted[:elevation], Cells)

        # A foreign manifest that sorts BEFORE ours does not shadow it. The
        # sidecar is found by the marker it carries, so the search is for a
        # marker this reader can trust and not for the first array wearing one;
        # kills a first-hit lookup, which costs the store its fast open on the
        # strength of somebody else's array name.
        shadowed = manifest_store(dir; name="shadowed.zarr")
        h = Zarr.zopen(shadowed, "w")
        foreign = Zarr.zcreate(UInt64, h, "aaa_manifest", 2, 1; chunks=(2, 1),
            attrs=merge(adims("chunks", "bounds"),
                d("dggs_chunk_manifest" => d("writer" => "xarray-dggs"))))
        foreign[:, :] = zeros(UInt64, 2, 1)
        g = counting(shadowed)
        st = dggread(g)
        @test reads(g, "cell_ids") == 0
        @test DD.lookup(st[:elevation], Cells) == DD.lookup(trusted[:elevation], Cells)

        # And the trusted axis is the axis a scan would have built.
        @test collect(DD.lookup(trusted[:elevation], Cells)) == CELLS
    end
end

@testset "a sidecar the rows themselves disqualify is ignored, not an error" begin
    # Three ways a table can fail before it is believed, one fixture each. All
    # three are the same verdict — a manifest this reader cannot use is an extra
    # array, so the scan runs and the store opens — and each is the ONLY thing
    # wrong with its store, so deleting any one of the three guards is caught
    # here and nowhere else.
    mktempdir() do dir
        nc = cld(N, COORD_CHUNK)
        trusted = DD.lookup(dggread(manifest_store(dir; name="ref.zarr"))[:elevation],
            Cells)

        # One row more than the chunk grid has chunks, and a row that ascends
        # away from the last one so that nothing but the SHAPE is wrong with it:
        # the table describes some other store's axis, however well formed.
        tall = counting(manifest_store(dir; name="tallrows.zarr",
            rows=R -> hcat(R, [R[2, end] + 0x1, R[2, end] + 0x2])))

        # One row whose two entries are the wrong way round. The rows still
        # ascend and stay disjoint ACROSS rows — the existing descending fixture
        # is what kills that condition — so only the per-row `first <= last`
        # comparison rejects this one, and a chunk whose interval runs backwards
        # answers `nothing` for every id genuinely inside it.
        swapped = counting(manifest_store(dir; name="swapped.zarr",
            rows=R -> swaprow(R, 2)))

        # A sidecar written at another width. These ids sit in base cell 8, so
        # as `Int64` they are negative — and wrap back to exactly the right
        # `UInt64` ids, which is what makes reinterpreting them look harmless.
        @test all(>=(UInt64(1) << 63), HIGH_IDS)
        narrow = counting(manifest_store(dir; name="narrow.zarr", ids=HIGH_IDS,
            T=Int64))

        for g in (tall, swapped)
            st = dggread(g)
            @test reads(g, "cell_ids") == nc
            @test DD.lookup(st[:elevation], Cells) == trusted
        end
        st = dggread(narrow)
        @test reads(narrow, "cell_ids") == nc
        @test collect(DD.lookup(st[:elevation], Cells)) == [Z7Cell(x) for x in HIGH_IDS]
    end
end

@testset "a coordinate written as Int64 is the same axis as one written UInt64" begin
    # The one normative rule about integer width: a reader REINTERPRETS a stored
    # id, it does not convert it. Base-cell-8 Z7 ids set the top bit, so a store
    # that wrote its coordinate as `Int64` — which xarray does whenever the ids
    # went through a signed dtype — holds negative integers that are exactly
    # these cells' bits. Kills a reader that converts and raises `InexactError`
    # on a store half the earth's cells cannot avoid.
    #
    # The sidecar is the deliberate asymmetry: `narrow.zarr` above declines a
    # MANIFEST of another width, because declining one costs a scan and nothing
    # else, while declining the coordinate would make the store unreadable.
    mktempdir() do dir
        @test all(>=(UInt64(1) << 63), HIGH_IDS)
        for (wide, narrow) in ((dense_store(dir; name="wide_dense.zarr", ids=HIGH_IDS),
                dense_store(dir; name="narrow_dense.zarr", ids=HIGH_IDS, T=Int64)),
            (ranges_store(dir; name="wide_ranges.zarr", ids=HIGH_IDS),
                ranges_store(dir; name="narrow_ranges.zarr", ids=HIGH_IDS, T=Int64)))
            w, n = dggread(wide), dggread(narrow)
            @test collect(DD.lookup(n[:elevation], Cells)) == [Z7Cell(x) for x in HIGH_IDS]
            @test DD.lookup(n[:elevation], Cells) == DD.lookup(w[:elevation], Cells)
            # And a cell still resolves to its own value, which is the half of
            # the bijection a sign-flipped id would break silently.
            @test n[:elevation][Cells(DD.At(Z7Cell(HIGH_IDS[137])))] == ELEVATION[137]
        end
    end
end

@testset "a coordinate value that is no id of this width is named, not InexactError" begin
    # The other side of the same rule. Reinterpretation is for ids that fit the
    # grid's width; a stored value that fits neither that nor a plain conversion
    # is a malformed coordinate, and saying so beats an `InexactError` from four
    # frames down with no store in it.
    mktempdir() do dir
        grid = levelgrid(HEALPixSystem(), 2)
        n = ncells(grid)
        path = joinpath(dir, "fractional.zarr")
        g = Zarr.zgroup(path; attrs=d())
        c = Zarr.zcreate(Float64, g, "cell_ids", n; chunks=(64,),
            attrs=merge(adims("cell_ids"), d("grid_name" => "healpix",
                "level" => 2, "indexing_scheme" => "nested")))
        values = Float64.(0:n-1)
        values[5] = 3.5
        c[:] = values
        v = Zarr.zcreate(Float32, g, "t2m", n; chunks=(64,), attrs=adims("cell_ids"))
        v[:] = Float32.(1:n)

        err = try
            dggread(path)
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :coordinate_width
        msg = sprint(showerror, err)
        @test occursin("3.5", msg) && occursin("Int64", msg)
        @test occursin("fractional.zarr", err.store)
    end
end

@testset "the stub says so when the extension is loaded and the call is wrong" begin
    # `dggread`/`dggwrite` are stubs in the main package and methods in the
    # extension, so everything that does not match a method falls back to the
    # stub. Telling a caller to run `using Zarr` when Zarr is already loaded
    # sends them to fix the one thing that is not wrong; kills a stub message
    # that never looks.
    @test Base.get_extension(DGG, :DiscreteGlobalGridsZarrExt) !== nothing
    for call in (() -> dggread(42), () -> dggwrite(42, cube()))
        err = try
            call()
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("extension is loaded", err.msg)
        @test !occursin("using Zarr", err.msg)
    end
end

@testset "a marker that does not name this store's level is not its manifest" begin
    # The store `dggwrite` produced, with both halves of the dual stamp edited
    # from level 4 to level 3 — the Ifremer failure under a sidecar. Every
    # GEOMETRY field of the marker still agrees: same dimension, same chunk
    # length, same length. A marker that carries no level and no grid is
    # therefore still trusted, and the store opens claiming a level none of its
    # ids names. Kills exactly that: the marker has to say what it validated.
    mktempdir() do dir
        path = joinpath(dir, "relevelled.zarr")
        dggwrite(path, cube(); encoding=:dense, chunks=COORD_CHUNK)
        group = joinpath(path, ".zattrs")
        coord = joinpath(path, "cell_ids", ".zattrs")
        write(group, replace(read(group, String),
            "\"refinement_level\":4" => "\"refinement_level\":3"))
        write(coord, replace(read(coord, String), "\"level\":4" => "\"level\":3"))
        # The edits landed: a fixture that quietly stopped lying would pass for
        # the wrong reason.
        @test occursin("\"refinement_level\":3", read(group, String))
        @test occursin("\"level\":3", read(coord, String))

        g = counting(path)
        err = try
            dggread(g)
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :id_names_no_cell
        @test err.declared == 3
        # It fell back rather than raising on the marker: the ids were read.
        @test reads(g, "cell_ids") > 0
        # And the error names the level the ids themselves are at.
        @test occursin("level 4", sprint(showerror, err))
    end
end

@testset "validate = :scan declines a sidecar this reader would trust" begin
    # The switch that puts a store this package wrote back on the ordinary
    # scanning path: the sidecar is ignored, every chunk of ids is read, and
    # every id is checked. Kills a `validate` that only ever chooses between
    # sample counts, leaving a trusted store with no way to be checked at all.
    mktempdir() do dir
        path = joinpath(dir, "scanned.zarr")
        dggwrite(path, cube(); encoding=:dense, chunks=COORD_CHUNK)

        g = counting(path)
        st = dggread(g; validate=:scan)
        @test reads(g, "cell_ids") == cld(N, COORD_CHUNK)
        @test collect(DD.lookup(st[:elevation], Cells)) == CELLS

        # The same store on the default setting reads no id at all, which is
        # what makes the two settings a choice rather than a spelling.
        h = counting(path)
        @test length(DD.lookup(dggread(h)[:elevation], Cells)) == N
        @test reads(h, "cell_ids") == 0
    end
end

# ---------------------------------------------------------------------------
# The ranges path
# ---------------------------------------------------------------------------

@testset "a ranges store is the same cube, read with no data IO" begin
    mktempdir() do dir
        dense = dggread(dense_store(dir))
        g = counting(ranges_store(dir))
        ranged = dggread(g)

        # The range array is small and read eagerly; no data chunk is touched,
        # and the axis it produced is complete before any of them is.
        @test reads(g, "elevation") == 0
        @test reads(g, "slope") == 0
        @test length(DD.lookup(ranged[:elevation], Cells)) == N

        # Twin equivalence: the axis is arithmetic on one side and a scan on
        # the other, and they name the same cells in the same order.
        @test DD.lookup(ranged[:elevation], Cells) == DD.lookup(dense[:elevation], Cells)
        @test collect(DD.lookup(ranged[:elevation], Cells)) == CELLS
        for v in (:elevation, :slope)
            @test collect(parent(ranged[v])) == collect(parent(dense[v]))
        end
        @test DD.metadata(ranged)["encoding"] == "ranges"
    end
end

# ---------------------------------------------------------------------------
# The implicit path
# ---------------------------------------------------------------------------

@testset "an implicit store has an axis and a time dimension and no cell array" begin
    mktempdir() do dir
        g = counting(dkrz_store(dir))
        st = dggread(g)
        grid = levelgrid(HEALPixSystem(), 2)

        @test collect(keys(st)) == [:tas]
        @test DD.metadata(st)["encoding"] == "implicit"
        # An implicit axis is arithmetic: the store is opened without reading
        # one byte of the array it describes.
        @test reads(g, "tas") == 0

        lk = DD.lookup(st[:tas], Cells)
        @test length(lk) == ncells(grid)
        @test collect(lk) == [cellindex(grid, i) for i in 1:ncells(grid)]

        # A non-cell dimension is an ordinary Dim, and Julia's dimension order
        # is the reverse of the store's.
        @test DD.name.(DD.dims(st[:tas])) == (:Cells, :time)
        @test collect(DD.lookup(st[:tas], DD.Dim{:time})) == [2020, 2021, 2022]
        @test st[:tas][Cells(DD.At(cellindex(grid, 7))), DD.Dim{:time}(DD.At(2021))] ==
              Float32(ncells(grid) + 7)
    end
end

@testset "a healpix store is read nested, and refused in ring" begin
    mktempdir() do dir
        grid = levelgrid(HEALPixSystem(), 2)
        n = ncells(grid)
        for scheme in ("nested", "ring")
            g = Zarr.zgroup(joinpath(dir, scheme * ".zarr"); attrs=d())
            c = Zarr.zcreate(Int64, g, "cell_ids", n; chunks=(64,),
                attrs=merge(adims("cell_ids"), d("grid_name" => "healpix",
                    "level" => 2, "indexing_scheme" => scheme)))
            c[:] = Int64.(0:n-1)
            v = Zarr.zcreate(Float32, g, "t2m", n; chunks=(64,), attrs=adims("cell_ids"))
            v[:] = Float32.(1:n)
        end

        st = dggread(joinpath(dir, "nested.zarr"))
        @test collect(DD.lookup(st[:t2m], Cells)) == [cellindex(grid, i) for i in 1:n]
        @test st[:t2m][Cells(DD.At(cellindex(grid, 42)))] == 42.0f0

        # The same integers in another order are another axis, and reading them
        # as nested would misplace every cell without saying so.
        err = try
            dggread(joinpath(dir, "ring.zarr"))
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :unsupported_indexing_scheme
        @test occursin("nested", sprint(showerror, err))
    end
end

# ---------------------------------------------------------------------------
# Several conventions on one store
# ---------------------------------------------------------------------------

@testset "a dual-stamped store merges, and a contradiction does not" begin
    mktempdir() do dir
        both = dense_store(dir; name="dual.zarr",
            coord=merge(adims("cell_ids"), xdggs_attrs()))
        st = dggread(both)
        @test DD.metadata(st)["conventions"] == ["zarr-conventions/dggs", "xdggs"]
        @test collect(DD.lookup(st[:elevation], Cells)) == CELLS

        # The same store with the two halves disagreeing about the level: the
        # "attrs lie" failure, refused rather than guessed, and the error names
        # the store it came from.
        lying = dense_store(dir; name="lying.zarr",
            coord=merge(adims("cell_ids"), xdggs_attrs(level=5)))
        err = try
            dggread(lying)
        catch e
            e
        end
        @test err isa DGGSFormatError
        @test err.check === :level_disagreement
        @test occursin("lying.zarr", err.store)
        @test err.conventions == ["zarr-conventions/dggs", "xdggs"]
    end
end

# ---------------------------------------------------------------------------
# Asserting the description instead of detecting it
# ---------------------------------------------------------------------------

@testset "an attribute-less store is read from a supplied description" begin
    mktempdir() do dir
        path = dense_store(dir; name="bare.zarr", group=d())
        err = try
            dggread(path)
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :no_convention_detected

        desc = StoreDescription(gridname="igeo7", system=IGeo7System(), level=4,
            encoding=DenseEncoding(), coordinate="cell_ids",
            spatial_dimension="cell_ids")
        st = dggread(path; description=desc)
        @test collect(DD.lookup(st[:elevation], Cells)) == CELLS
        @test collect(keys(st)) == [:elevation, :slope]
        @test DD.metadata(st)["conventions"] == String[]

        # Mechanical checks still run: the ids are scanned against the level the
        # caller asserted, and a wrong one is caught.
        wrong = StoreDescription(gridname="igeo7", system=IGeo7System(), level=5,
            encoding=DenseEncoding(), coordinate="cell_ids",
            spatial_dimension="cell_ids")
        err = try
            dggread(path; description=wrong)
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :id_names_no_cell
        @test occursin("bare.zarr", err.store)
    end
end

# ---------------------------------------------------------------------------
# The validation pair
# ---------------------------------------------------------------------------

@testset "strict verifies every id where lazy samples, on a scanned store" begin
    mktempdir() do dir
        # One chunk of 100 ids with a phantom at slot 3 — where an even sample
        # spread does not look. The array stays sorted, so only validity can
        # reject it. This store carries no manifest, so the scan is what runs;
        # a store whose sidecar this reader trusts checks no id at all until
        # `:scan` asks it to.
        ids = UInt64[idselect(GRID, r) for r in 0:99]
        ids[3] = PHANTOM
        @test issorted(ids) && allunique(ids)
        path = dense_store(dir; name="phantom.zarr", ids=ids)

        # `:scan` is `:strict` on a store with nothing to decline, which is what
        # makes it safe to reach for: it never checks less.
        for validate in (:strict, :scan)
            err = try
                dggread(path; validate)
            catch e
                e
            end
            @test err isa DGGSFormatError && err.check === :id_names_no_cell
            @test occursin("phantom.zarr", err.store)
            @test err.conventions == ["zarr-conventions/dggs"]
        end

        st = dggread(path; validate=:lazy)
        @test length(DD.lookup(st[:elevation], Cells)) == 100

        @test_throws ArgumentError dggread(path; validate=:sometimes)
    end
end

@testset "a duplicated id is refused whatever the validation setting, when scanned" begin
    mktempdir() do dir
        ids = copy(IDS)
        ids[130] = ids[129]
        path = dense_store(dir; name="dup.zarr", ids=ids)
        for validate in (:strict, :lazy, :scan)
            err = try
                dggread(path; validate)
            catch e
                e
            end
            @test err isa DGGSFormatError && err.check === :duplicate_ids
            @test occursin("dup.zarr", err.store)
        end
    end
end

# ---------------------------------------------------------------------------
# What the store is not
# ---------------------------------------------------------------------------

@testset "a path that names an array is not a store" begin
    # A DGGS store is a GROUP of arrays over one cell axis, and pointing at a
    # single array is the everyday typo. Kills a reader that lets a ZArray
    # through to fail somewhere less legible.
    mktempdir() do dir
        path = joinpath(dir, "lonely.zarr")
        z = Zarr.zcreate(Float32, Zarr.DirectoryStore(path), 8; chunks=(8,))
        z[:] = Float32.(1:8)

        err = try
            dggread(path)
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :not_a_group
        @test occursin("lonely.zarr", err.store)
    end
end

@testset "a coordinate the store does not hold is named as missing" begin
    # Both ways a dense layout can end up without its ids: a description that
    # names an array which is not there, and one that names none.
    mktempdir() do dir
        path = dense_store(dir; name="nocoord.zarr")

        wrongname = StoreDescription(gridname="igeo7", system=IGeo7System(), level=4,
            encoding=DenseEncoding(), coordinate="cell_index",
            spatial_dimension="cell_ids")
        err = try
            dggread(path; description=wrongname)
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :missing_coordinate
        msg = sprint(showerror, err)
        @test occursin("cell_index", msg) && occursin("cell_ids", msg)

        anonymous = StoreDescription(gridname="igeo7", system=IGeo7System(), level=4,
            encoding=DenseEncoding(), spatial_dimension="cell_ids")
        err = try
            dggread(path; description=anonymous)
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :missing_coordinate
    end
end

@testset "two layers disagreeing about a dimension's length is a format error" begin
    # The cell axis already gets this check; a `time` two arrays give different
    # lengths is the same lie about the same store, and has to read as one.
    # Kills a dimension cache keyed on the name alone.
    mktempdir() do dir
        path = joinpath(dir, "ragged.zarr")
        g = Zarr.zgroup(path; attrs=dggs_attrs())
        c = Zarr.zcreate(UInt64, g, "cell_ids", N; chunks=(COORD_CHUNK,),
            attrs=adims("cell_ids"))
        c[:] = IDS
        for (name, steps) in (("pr", 4), ("tas", 3))
            v = Zarr.zcreate(Float32, g, name, N, steps; chunks=(64, steps),
                attrs=adims("time", "cell_ids"))
            v[:, :] = zeros(Float32, N, steps)
        end

        err = try
            dggread(path)
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :dimension_length_disagreement
        msg = sprint(showerror, err)
        @test occursin("time", msg) && occursin("3", msg) && occursin("4", msg)
    end
end

# ---------------------------------------------------------------------------
# Variable selection
# ---------------------------------------------------------------------------

@testset "vars selects layers and an unknown one lists what there is" begin
    mktempdir() do dir
        path = dense_store(dir)
        @test collect(keys(dggread(path; vars=(:elevation,)))) == [:elevation]
        @test collect(keys(dggread(path; vars=[:slope, :elevation]))) == [:slope, :elevation]

        err = try
            dggread(path; vars=(:nope,))
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :unknown_variable
        msg = sprint(showerror, err)
        @test occursin("nope", msg) && occursin("elevation", msg) && occursin("slope", msg)

        # The single-variable method is a DimArray, not a one-layer stack.
        A = dggread(path, :elevation)
        @test A isa DD.DimArray
        @test collect(parent(A)) == ELEVATION
        @test collect(DD.lookup(A, Cells)) == CELLS
    end
end

# ---------------------------------------------------------------------------
# URLs
# ---------------------------------------------------------------------------

@testset "a gs:// URL is the https one, and s3:// says so" begin
    # String level only: no request is made by either assertion.
    @test ZarrExt.normalize_store_url("gs://geo-assets/igeo7-zarr/pori.zarr") ==
          "https://storage.googleapis.com/geo-assets/igeo7-zarr/pori.zarr"
    @test ZarrExt.normalize_store_url("gs://bucket") ==
          "https://storage.googleapis.com/bucket"
    for u in ("https://storage.googleapis.com/geo-assets/x.zarr",
        "http://localhost:8000/x.zarr", "/local/path/x.zarr")
        @test ZarrExt.normalize_store_url(u) == u
    end

    err = try
        dggread("s3://some-bucket/some.zarr")
    catch e
        e
    end
    @test err isa DGGSFormatError && err.check === :unsupported_store_scheme
    msg = sprint(showerror, err)
    @test occursin("AWSS3", msg) && occursin("https", msg)
end

# ---------------------------------------------------------------------------
# Compacted stores
# ---------------------------------------------------------------------------

@testset "a compacted store round-trips a coarsened field" begin
    # The whole mixed-level path on both radices: coarsen leaf data, write,
    # read back, query, expand. The field is exactly flat where it merges, so
    # the expansion of the read-back cube reproduces the leaf values
    # bit-for-bit — which kills any reordering of cells against values, a lost
    # or misread levels column, and a signed/unsigned id-column mixup.
    for sys in (HEALPixSystem(), IGeo7System())
        L = 3
        grid = levelgrid(sys, L)
        cv = CellVector(grid)
        n = ncells(grid)
        leafvals = Float64[k <= n ÷ 2 ? 1.0 : k for k in 1:n]
        mov, vals = DGG.coarsen(cv, leafvals; atol=0.5)
        @test length(mov) < n

        M = DD.DimArray(vals, Cells(DGG.MultiOrderLookup(mov)); name=:field)
        path = joinpath(mktempdir(), "compacted.zarr")
        dggwrite(path, M)

        A = dggread(path, :field; lazy=false)
        lk = DD.lookup(A, Cells)
        @test lk isa DGG.MultiOrderLookup
        @test collect(parent(lk)) == collect(mov)
        @test collect(A) == vals

        # Contains resolves a leaf cell to its stored ancestor's value.
        leaf = cv[1]
        k = DGG.covering_index(parent(lk), leaf)
        @test A[Cells(DD.Contains(leaf))] == vals[k]

        # The read-back cube expands to the leaf level it was coarsened from.
        @test collect(DGG.expand(A, L)) == leafvals
    end
end

end # if HAS_ZARR

end # module DGGSIOReadTests
