# The store corpus: one cell set written every way this package can write it,
# and one trap fixture per §5 check that only real bytes can spring.
#
# What this file adds to `read.jl` and `write.jl` is the CROSS-STORE claim. Those
# two each exercise one direction against one store; here the same cells are
# written dense, as `:rank` ranges, as `:step` ranges, and — where the ids permit
# — as an implicit axis with no coordinate at all, and every spelling has to read
# back as the same cube with the same chunk grid. An expansion-semantics mutant,
# a rank/select mutant and a scan mutant cannot all agree, so one comparison
# kills whichever of them is wrong.
#
# The corpus is BUILT, not committed: `dggwrite` is under test in `write.jl` and
# trustworthy by the time this file runs, and blosc is not byte-deterministic
# across versions, so a committed `.zarr` would pin bytes nobody promised. The
# generators below are the fixture.
#
# The suite self-skips when Zarr.jl is absent, and the public-store suite is
# gated on `DGG_IO_NETWORK_TESTS=1`.

module DGGSIOStoreTests

using Test
import DiscreteGlobalGrids as DGG

const HAS_ZARR = try
    @eval using Zarr
    true
catch err
    @warn "Zarr.jl is not loadable: the store-corpus suite is skipped." exception = err
    false
end

if HAS_ZARR

import DimensionalData as DD
using DiscreteGlobalGrids: IGeo7System, HEALPixSystem, Z7Cell, levelgrid, ncells,
    cellindex, CellVector, CellLookup, Cells, DGGSFormatError, dggread, dggwrite,
    chunkmanifest
using DiscreteGlobalGrids.Encodings: idselect, idranges

d(pairs...) = Dict{String,Any}(pairs...)
adims(names...) = d("_ARRAY_DIMENSIONS" => Any[names...])

"`f()`'s exception, or `nothing` when it returns."
caught(f) = try
    f()
    nothing
catch e
    e
end

# ---------------------------------------------------------------------------
# The corpus
# ---------------------------------------------------------------------------

const SYS = IGeo7System()
const LEVEL = 4
const GRID = levelgrid(SYS, LEVEL)
# Three position runs: the shape of a regional store, and each long enough to
# span several digit rollovers, which is what makes `:rank` and `:step` disagree
# about how many runs there are.
const RANKS = [100:199; 500:519; 900:929]
const IDS = UInt64[idselect(GRID, r) for r in RANKS]
const CELLS = [Z7Cell(x) for x in IDS]
const N = length(IDS)
const ELEVATION = Float32.(10 .* (1:N))
const SLOPE = Float64.(inv.(1:N))
# Three chunks over the cell axis, so a manifest has interior rows to disagree
# about and a selector has a wrong chunk to land in.
const CHUNK = 50

# The whole of HEALPix level 2. An implicit axis is a COMPLETE level and nothing
# less, and the DKRZ dialect that writes them is healpix-only; the same claim on
# the IGEO7 corpus would need all 24 012 cells of level 4 rather than 150. So the
# implicit twins get their own small stack.
const HSYS = HEALPixSystem()
const HGRID = levelgrid(HSYS, 2)
const HN = ncells(HGRID)
const HCELLS = [cellindex(HGRID, i) for i in 1:HN]
const TAS = Float32.(1:HN)

const CORPUS = Ref{Any}(nothing)

"""
    corpus() -> NamedTuple

The twin corpus, written once per test run into a temporary directory: the same
150 IGEO7 cells as a dense store, a `:rank`-merged ranges store, a `:step`-merged
one and one written with no `merge` keyword at all, and the same complete HEALPix
level four ways — ranges, dense, this package's own implicit store, and a
hand-stamped DKRZ one.

Every path names a store this package wrote except `hpx_dkrz`, which is another
producer's dialect and the reason the implicit read path exists.
"""
function corpus()
    CORPUS[] === nothing || return CORPUS[]
    dir = mktempdir()
    stack = DD.DimStack((elevation=copy(ELEVATION), slope=copy(SLOPE)),
        (Cells(CellLookup(CellVector(SYS, LEVEL, CELLS))),))
    hpx = DD.DimArray(copy(TAS), (Cells(CellLookup(CellVector(HSYS, 2, HCELLS))),);
        name=:tas)

    paths = (dense=dggwrite(joinpath(dir, "dense.zarr"), stack;
            encoding=:dense, chunks=CHUNK),
        ranges=dggwrite(joinpath(dir, "ranges.zarr"), stack;
            encoding=:ranges, merge=:rank, chunks=CHUNK),
        step=dggwrite(joinpath(dir, "step.zarr"), stack;
            encoding=:ranges, merge=:step, chunks=CHUNK),
        default=dggwrite(joinpath(dir, "default.zarr"), stack;
            encoding=:ranges, chunks=CHUNK),
        hpx_ranges=dggwrite(joinpath(dir, "healpix.zarr"), hpx; chunks=64),
        hpx_dense=dggwrite(joinpath(dir, "healpix_dense.zarr"), hpx;
            encoding=:dense, chunks=64),
        hpx_implicit=dggwrite(joinpath(dir, "healpix_implicit.zarr"), hpx;
            encoding=:implicit, chunks=64),
        hpx_dkrz=dkrz_store(joinpath(dir, "dkrz.zarr")))
    return CORPUS[] = paths
end

"A nextGEMS-style store of the same cells: a `crs` variable, no cell array."
function dkrz_store(path)
    g = Zarr.zgroup(path; attrs=d())
    crs = Zarr.zcreate(Int8, g, "crs", 1; chunks=(1,),
        attrs=merge(adims("crs"), d("grid_mapping_name" => "healpix",
            "healpix_nside" => 2^2, "healpix_order" => "nest")))
    crs[:] = Int8[0]
    tas = Zarr.zcreate(Float32, g, "tas", HN; chunks=(64,),
        attrs=merge(adims("cell"), d("grid_mapping" => "crs")))
    tas[:] = TAS
    return path
end

"""
Every layer of `a` and `b`, cell for cell and value for value.

Values are compared with `isequal`, because a nodata `NaN` is a value a store
holds and a twin has to hold the same one; `==` calls two of them different.
"""
function sametwin(a, b)
    @test collect(keys(a)) == collect(keys(b))
    for k in keys(a)
        la, lb = DD.lookup(a[k], Cells), DD.lookup(b[k], Cells)
        @test la == lb
        @test collect(la) == collect(lb)
        @test isequal(collect(parent(a[k])), collect(parent(b[k])))
    end
end

# ---------------------------------------------------------------------------
# Twin equivalence
# ---------------------------------------------------------------------------

@testset "the dense and the ranges store are one cube" begin
    # The load-bearing corpus assertion. One side scanned 150 stored ids, the
    # other computed them from three intervals by rank/select, and they have to
    # name the same cells in the same order against the same values. Kills every
    # expansion-semantics mutant at once — an exclusive `stop`, a step that walks
    # integers instead of cells, a validity filter applied to interval interiors
    # — because any of them shifts the ranges axis against the dense one.
    c = corpus()
    dense, ranged = dggread(c.dense), dggread(c.ranges)
    @test DD.metadata(dense)["encoding"] == "none"
    @test DD.metadata(ranged)["encoding"] == "ranges"

    sametwin(dense, ranged)
    @test collect(DD.lookup(dense[:elevation], Cells)) == CELLS
    # Addressed by cell rather than by position, which is the only way the axis
    # and the data can be caught disagreeing anywhere but position 1.
    for k in (1, 101, N)
        @test ranged[:elevation][Cells(DD.At(CELLS[k]))] == ELEVATION[k]
        @test ranged[:slope][Cells(DD.At(CELLS[k]))] == SLOPE[k]
    end
end

@testset "the chunk grid is the same from arithmetic and from a scan" begin
    # The dense manifest is what the open-time scan recorded; the ranges manifest
    # is closed-form rank/select over three intervals, computed without reading a
    # byte. Kills rank/select and scan mutants against each other, and does it at
    # a chunk length the store was NOT written in as well as the one it was, so a
    # manifest that only ever reproduces its own store's boundaries fails.
    c = corpus()
    dense = DD.lookup(dggread(c.dense)[:elevation], Cells)
    ranged = DD.lookup(dggread(c.ranges)[:elevation], Cells)
    for cl in (CHUNK, 37, N)
        @test chunkmanifest(dense, cl) == chunkmanifest(ranged, cl)
    end
    m = chunkmanifest(dense, CHUNK)
    @test m.lengths == [CHUNK, CHUNK, CHUNK]
    @test m.firstids == IDS[1:CHUNK:N]

    # And what the writer persisted is what the reader rebuilds: the sidecar is
    # the manifest a store of tens of millions of cells will be opened on instead
    # of scanning, so the two must be the same table.
    sidecar = permutedims(Zarr.zopen(c.dense)["cell_chunk_manifest"][:, :])
    @test sidecar[:, 1] == m.firstids
    @test sidecar[:, 2] == m.lastids
end

@testset "the merge rule is a spelling, not an axis" begin
    # `:step` cannot enclose an id naming no cell and `:rank` can, so the two
    # stores hold different intervals — 24 rows against 3 for these cells. The
    # reader counts cells inside an interval either way, so both are the same
    # axis. Kills a reader that recovers the axis by stepping integers, which
    # works on the `:step` store and misplaces every cell after the first
    # rollover on the `:rank` one.
    c = corpus()
    # `(n, 2)` in the store is `(2, n)` here, so the row count is the second
    # extent: 24 intervals against 3 for the same 150 cells. Pinned rather than
    # compared, because `>` also holds for a `:step` rule that merges nothing at
    # all — 150 rows of one cell each — which is not what it does.
    @test size(Zarr.zopen(c.step)["cell_id_ranges"], 2) == 24
    @test size(Zarr.zopen(c.ranges)["cell_id_ranges"], 2) == 3

    stepped = dggread(c.step)
    sametwin(dggread(c.ranges), stepped)
    @test collect(DD.lookup(stepped[:elevation], Cells)) == CELLS
end

@testset "a write that names no merge rule is step-merged" begin
    # The default is interop, not compactness: a caller who names no rule gets
    # the intervals a structural-count reader also counts correctly, which is
    # what the published stores hold. Kills a default of `:rank`, which for these
    # cells writes 3 rows enclosing ids that name no cell — same cube, but only
    # for a reader that expands an interval by rank.
    c = corpus()
    @test Zarr.zopen(c.default)["cell_id_ranges"][:, :] ==
          Zarr.zopen(c.step)["cell_id_ranges"][:, :]
    @test size(Zarr.zopen(c.default)["cell_id_ranges"], 2) == 24
    sametwin(dggread(c.dense), dggread(c.default))
end

@testset "an implicit axis is the same complete level, ours and theirs" begin
    # The whole of HEALPix level 2 written four ways: as intervals, as stored
    # ids, as this package's own implicit store, and as a DKRZ store from another
    # producer's dialect. The last two write no axis at all and let position be
    # the id, which kills an off-by-one there — rank is zero-based and position
    # is one-based — that nothing else catches, since every other encoding reads
    # its ids back from the store.
    c = corpus()
    ranged = dggread(c.hpx_ranges)
    dense = dggread(c.hpx_dense)
    ours = dggread(c.hpx_implicit)
    theirs = dggread(c.hpx_dkrz)
    @test DD.metadata(ours)["encoding"] == "implicit"
    @test DD.metadata(theirs)["encoding"] == "implicit"
    @test DD.metadata(theirs)["conventions"] == ["dkrz-healpix"]

    sametwin(ranged, ours)
    sametwin(dense, ours)
    sametwin(ours, theirs)
    @test collect(DD.lookup(theirs[:tas], Cells)) == HCELLS
    for cl in (64, 50)
        @test chunkmanifest(DD.lookup(dense[:tas], Cells), cl) ==
              chunkmanifest(DD.lookup(theirs[:tas], Cells), cl)
    end

    # An implicit store is one whose coordinate is ABSENT: the convention says so
    # by leaving the key out, and there is no array to say it with. Kills a write
    # path that stamps a coordinate name the store does not hold, which reads
    # back as a dense store missing its ids.
    g = Zarr.zopen(c.hpx_implicit)
    @test !haskey(g.attrs["dggs"], "coordinate")
    @test !any(startswith("cell_id"), keys(g.arrays))
end

# ---------------------------------------------------------------------------
# Trap fixtures: hand-built, because `dggwrite` refuses to produce any of them
# ---------------------------------------------------------------------------

zarr_conventions() = Any[deepcopy(DGG.ZARR_DGGS_DECLARATION)]

function dggs_attrs(; level=LEVEL, compression="none", coordinate="cell_ids")
    dggs = d("name" => "igeo7", "refinement_level" => level,
        "spatial_dimension" => "cell_ids")
    coordinate === nothing || (dggs["coordinate"] = coordinate;
    dggs["compression"] = compression)
    return d("zarr_conventions" => zarr_conventions(), "dggs" => dggs)
end

"A dense store, `zarr-conventions/dggs` on the group and whatever is asked on the coordinate."
function dense_store(dir, name; ids=IDS, group=dggs_attrs(), coord=adims("cell_ids"),
    declared=length(ids), chunks=CHUNK)
    path = joinpath(dir, name)
    g = Zarr.zgroup(path; attrs=group)
    c = Zarr.zcreate(UInt64, g, "cell_ids", length(ids); chunks=(chunks,), attrs=coord)
    c[:] = ids
    e = Zarr.zcreate(Float32, g, "elevation", declared; chunks=(chunks,),
        attrs=adims("cell_ids"))
    e[:] = Float32.(1:declared)
    return path, c
end

"A ranges store holding `rows` verbatim, however malformed they are."
function ranges_store(dir, name, rows; declared=N)
    path = joinpath(dir, name)
    g = Zarr.zgroup(path; attrs=dggs_attrs(compression="ranges",
        coordinate="cell_id_ranges"))
    # Zarr declares shapes outermost-first, Julia innermost-first: an `(m, 2)`
    # array is `(2, m)` here.
    c = Zarr.zcreate(UInt64, g, "cell_id_ranges", 2, size(rows, 1);
        chunks=(2, size(rows, 1)), attrs=adims("ranges", "bounds"))
    c[:, :] = permutedims(rows)
    e = Zarr.zcreate(Float32, g, "elevation", declared; chunks=(CHUNK,),
        attrs=adims("cell_ids"))
    e[:] = Float32.(1:declared)
    return path
end

@testset "a duplicated cell axis is counted and located" begin
    # Two duplicates, one of them STRADDLING a chunk boundary: the pair is the
    # last id of chunk 1 and the first of chunk 2. Kills a scan that compares
    # within a chunk and starts each block fresh, which sees neither duplicate as
    # adjacent and opens a store that names one cell twice.
    mktempdir() do dir
        ids = copy(IDS)
        ids[CHUNK+1] = ids[CHUNK]
        ids[100] = ids[99]
        @test issorted(ids)

        err = caught(() -> dggread(dense_store(dir, "dup.zarr"; ids)[1]))
        @test err isa DGGSFormatError && err.check === :duplicate_ids
        @test err.observed == 2                 # how many
        @test err.declared == CHUNK + 1         # where the first one is
        @test occursin("dup.zarr", err.store)   # and which store holds them
        @test err.conventions == ["zarr-conventions/dggs"]
    end
end

@testset "attributes that lie about the level are refused, not reinterpreted" begin
    # The Ifremer pattern: a store whose attributes declare one level and whose
    # ids are written at another. Every id is then structurally well formed and
    # names nothing at the declared level, so only the validity check can catch
    # it — and it must, because reading level-4 ids as level 3 would answer every
    # selector with the wrong cell rather than with an error.
    mktempdir() do dir
        path, _ = dense_store(dir, "ifremer.zarr"; group=dggs_attrs(level=3))
        for validate in (:strict, :lazy)
            err = caught(() -> dggread(path; validate))
            @test err isa DGGSFormatError && err.check === :id_names_no_cell
            @test err.declared == 3             # the level the attributes claim
            @test err.observed == IDS[1]        # the id that contradicts it
            @test occursin("ifremer.zarr", err.store)
        end
        # Sampled or not, the FIRST id already contradicts the attributes, which
        # is why `:lazy` catches this one and misses a lone phantom.
        msg = sprint(showerror, caught(() -> dggread(path)))
        @test occursin("level 3", msg) && occursin("position 1", msg)
        # Design §5 asks for an error naming BOTH levels and pointing at the
        # reference-level pattern. The id's own level is read back out of the id
        # where the scheme carries one, which for Z7 it does.
        @test occursin("level 4", msg)
        @test occursin("reference level", msg)
    end
end

@testset "padding past the declared shape is invisible" begin
    # 150 cells in chunks of 64: the last chunk holds 22 of them and 42 slots of
    # pad, three of which this store filled with zeros before shrinking its
    # declared shape — Zarr's own semantics, and what a producer leaves behind
    # when it rewrites a store's length without rewriting its chunks. Kills a
    # reader that derives the axis length from the CHUNK GRID rather than from
    # the declared shape: three chunks of 64 is 192, and the zeros are neither
    # cells nor above the last id, so such a reader rejects the store as
    # unsorted instead of reading it.
    mktempdir() do dir
        path, coord = dense_store(dir, "padded.zarr";
            ids=[IDS; zeros(UInt64, 3)], declared=N, chunks=64)
        @test coord[N+1:N+3] == zeros(UInt64, 3)    # the garbage is really there

        meta = joinpath(path, "cell_ids", ".zarray")
        write(meta, replace(read(meta, String),
            "\"shape\":[$(N + 3)]" => "\"shape\":[$N]"))

        st = dggread(path)
        lk = DD.lookup(st[:elevation], Cells)
        @test length(lk) == N
        @test collect(lk) == CELLS
        @test st[:elevation][Cells(DD.At(CELLS[end]))] == Float32(N)
        # The manifest ends where the axis does: the final chunk is short, not
        # padded out to the chunk length the store still has room for.
        m = chunkmanifest(lk, 64)
        @test m.lengths == [64, 64, 22]
        @test m.lastids[end] == IDS[end]
    end
end

@testset "malformed ranges name the store that holds them" begin
    # Both structural checks, sprung by bytes rather than by a matrix argument,
    # which is what makes them a test of the extension: the range array is read
    # and transposed out of the store, and the encoding layer's rejection has to
    # cross the API boundary with the store filled in. Kills a read path that
    # validates ranges outside `with_store_context`, leaving a lazy cube's error
    # with no way to say which store it came from.
    mktempdir() do dir
        rows = idranges(GRID, IDS)
        @test size(rows, 1) == 3

        # Rows in descending order: they do not ascend, and the same comparison
        # is what a genuinely overlapping pair fails.
        err = caught(() -> dggread(ranges_store(dir, "descending.zarr", rows[3:-1:1, :])))
        @test err isa DGGSFormatError && err.check === :overlapping_ranges
        @test occursin("descending.zarr", err.store)
        @test occursin("ascend", sprint(showerror, err))

        # A row whose stop is below its start holds no cells at all.
        backwards = copy(rows)
        backwards[2, 1], backwards[2, 2] = backwards[2, 2], backwards[2, 1]
        err = caught(() -> dggread(ranges_store(dir, "backwards.zarr", backwards)))
        @test err isa DGGSFormatError && err.check === :empty_range_row
        @test occursin("backwards.zarr", err.store)
    end
end

@testset "the trust boundary a persisted manifest buys, pinned as intended" begin
    # NOT a bug report: this is the documented ceiling of design §4. The sidecar
    # is ours, the marker is ours, and the spot check the reader makes on every
    # chunk it decodes — first id, last id, length — passes, because the
    # duplicate is INSIDE a chunk and moves none of the three. `validated =
    # "strict"` is the writer's attestation that no such duplicate exists, and a
    # trusted store takes it: the axis opens holding the cell twice.
    #
    # What the trust does NOT buy is a wrong answer. `_axisposition` searches the
    # decoded chunk and compares the id it lands on, so every cell still resolves
    # to a position that really holds it. And `validate = :scan` declines the
    # sidecar and finds the duplicate, which is the way back.
    mktempdir() do dir
        path = joinpath(dir, "interior_dup.zarr")
        stack = DD.DimStack((elevation=copy(ELEVATION),),
            (Cells(CellLookup(CellVector(SYS, LEVEL, CELLS))),))
        dggwrite(path, stack; encoding=:dense, chunks=CHUNK)
        ids = copy(IDS)
        ids[3] = ids[2]                 # interior to chunk 1; 1 and 50 untouched
        Zarr.zopen(path, "w")["cell_ids"][:] = ids

        st = dggread(path)
        lk = DD.lookup(st[:elevation], Cells)
        cells = collect(lk)
        @test length(cells) == N
        @test cells[2] == cells[3]      # the duplicate is there to be seen

        # No cell is misplaced: whatever position each one resolves to, that
        # position holds it. `_axisposition` searches the decoded chunk and
        # compares the id it lands on, so the duplicate costs position 3 its
        # cell and costs no cell its position.
        @test all(eachindex(cells)) do k
            p = DGG.cellposition(lk, cells[k])
            p !== nothing && cells[p] == cells[k]
        end

        err = caught(() -> dggread(path; validate=:scan))
        @test err isa DGGSFormatError && err.check === :duplicate_ids
        @test err.declared == 3         # where the first one is
        @test occursin("interior_dup.zarr", err.store)
    end
end

@testset "two conventions that name different grids are refused" begin
    # The other half of the "attrs lie" failure: the group says igeo7 and the
    # coordinate says healpix, and the two are not the same cells under different
    # names — they are different tessellations whose ids happen to be integers.
    # `read.jl` pins the same refusal on the level; this pins that reconciliation
    # covers the field that decides which id arithmetic runs.
    mktempdir() do dir
        path, _ = dense_store(dir, "namefight.zarr";
            coord=merge(adims("cell_ids"), d("grid_name" => "healpix", "level" => LEVEL)))
        err = caught(() -> dggread(path))
        @test err isa DGGSFormatError && err.check === :gridname_disagreement
        @test err.declared == "igeo7" && err.observed == "healpix"
        @test err.conventions == ["zarr-conventions/dggs", "xdggs"]
        @test occursin("namefight.zarr", err.store)
    end
end

# ---------------------------------------------------------------------------
# The published stores
# ---------------------------------------------------------------------------

const NETWORK = get(ENV, "DGG_IO_NETWORK_TESTS", "0") in ("1", "true", "yes")
const PUBLIC = "https://storage.googleapis.com/geo-assets/igeo7-zarr/"

if !NETWORK
    @info "DGG_IO_NETWORK_TESTS is unset: the public-store suite is skipped. " *
          "Set DGG_IO_NETWORK_TESTS=1 to read the Pori stores over HTTPS."
else
    @testset "the published Pori stores read as themselves" begin
        # The same assertions as the synthetic corpus, against the three stores
        # this format was reverse-engineered from. They were written by another
        # implementation entirely (xdggs + dggrid4py), so what passes here is
        # cross-implementation agreement and not this package agreeing with
        # itself: the counts, the level and the layer names are the stores', and
        # a reader that only ever satisfies its own writer fails them.
        dense = dggread(PUBLIC * "pori_z7_r10.zarr"; validate=:strict)
        @test collect(keys(dense)) == [:elevation, :slope_geodesic, :slope_lookup]
        lk = DD.lookup(dense[:elevation], Cells)
        @test length(lk) == 3101
        @test DGG.level(lk) == 10
        @test DD.metadata(dense)["encoding"] == "none"
        @test DD.metadata(dense)["conventions"] == ["zarr-conventions/dggs", "xdggs"]

        # The real-world twin: 3 101 cells recovered from 136 intervals by
        # arithmetic, against 3 101 stored ids.
        ranged = dggread(PUBLIC * "pori_z7_r10_ranges.zarr")
        @test DD.metadata(ranged)["encoding"] == "ranges"
        @test DD.lookup(ranged[:elevation], Cells) == lk
        # `isequal`: the elevation model's nodata is NaN, and the twins hold the
        # same NaNs in the same places.
        for k in keys(ranged)
            @test isequal(collect(parent(ranged[k])), collect(parent(dense[k])))
        end
        @test chunkmanifest(DD.lookup(ranged[:elevation], Cells), 1000) ==
              chunkmanifest(lk, 1000)

        # 158 430 cells from 1 073 intervals. Opening it is already the normative
        # length check — the closed-form count of the intervals has to equal the
        # length the data arrays declare — and the data stays where it is.
        big = dggread(PUBLIC * "pori_z7_r12_ranges.zarr")
        blk = DD.lookup(big[:elevation], Cells)
        @test length(blk) == 158430
        @test DGG.level(blk) == 12
        @test parent(big[:elevation]) isa Zarr.ZArray
        @test !(parent(big[:elevation]) isa Array)
        # Forward and inverse agree in the middle of a level-12 axis, which is
        # further than any synthetic fixture here reaches. Position 100 000 holds
        # a real elevation rather than the model's NaN nodata, so `==` is the
        # comparison that means something.
        k = 100_000
        @test big[:elevation][Cells(DD.At(blk[k]))] == parent(big[:elevation])[k]
    end
end

end # if HAS_ZARR

end # module DGGSIOStoreTests
