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
    levelgrid, ncells, cellindex, Cells, DGGSFormatError, dggread,
    StoreDescription
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
# Three position runs: the shape of a regional store, and enough runs that the
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

"A dense store: one stored id per cell, `zarr-conventions/dggs` only."
function dense_store(dir; name="dense.zarr", group=dggs_attrs(), coord=adims("cell_ids"),
    ids=IDS)
    path = joinpath(dir, name)
    g = Zarr.zgroup(path; attrs=group)
    c = Zarr.zcreate(UInt64, g, "cell_ids", length(ids); chunks=(COORD_CHUNK,), attrs=coord)
    c[:] = ids
    write_data!(g, ids)
    return path
end

"The same cells as `dense_store`, as an `(n, 2)` inclusive-range coordinate."
function ranges_store(dir; name="ranges.zarr")
    path = joinpath(dir, name)
    g = Zarr.zgroup(path; attrs=dggs_attrs(compression="ranges",
        coordinate="cell_id_ranges"))
    rows = idranges(GRID, IDS)
    # Zarr declares shapes outermost-first, Julia innermost-first: an `(m, 2)`
    # array is `(2, m)` here.
    c = Zarr.zcreate(UInt64, g, "cell_id_ranges", 2, size(rows, 1);
        chunks=(2, size(rows, 1)), attrs=adims("ranges", "bounds"))
    c[:, :] = permutedims(rows)
    write_data!(g)
    return path
end

"A nextGEMS-style store: a `crs` variable, no cell array, position is the id."
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

        # Values, addressed by cell rather than by position: the axis and the
        # data have to agree for this to hold anywhere but position 1.
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

@testset "strict verifies every id where lazy samples" begin
    mktempdir() do dir
        # One chunk of 100 ids with a phantom at slot 3 — where an even sample
        # spread does not look. The array stays sorted, so only validity can
        # reject it.
        ids = UInt64[idselect(GRID, r) for r in 0:99]
        ids[3] = PHANTOM
        @test issorted(ids) && allunique(ids)
        path = dense_store(dir; name="phantom.zarr", ids=ids)

        err = try
            dggread(path)
        catch e
            e
        end
        @test err isa DGGSFormatError && err.check === :id_names_no_cell
        @test occursin("phantom.zarr", err.store)
        @test err.conventions == ["zarr-conventions/dggs"]

        st = dggread(path; validate=:lazy)
        @test length(DD.lookup(st[:elevation], Cells)) == 100

        @test_throws ArgumentError dggread(path; validate=:sometimes)
    end
end

@testset "a duplicated id is refused whatever the validation setting" begin
    mktempdir() do dir
        ids = copy(IDS)
        ids[130] = ids[129]
        path = dense_store(dir; name="dup.zarr", ids=ids)
        for validate in (:strict, :lazy)
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

end # if HAS_ZARR

end # module DGGSIOReadTests
