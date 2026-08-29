# Enumerated canonical cells and well-formed phantom ids provide the Z7 oracle.

using Test
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: IGeo7System, HEALPixSystem, Z7Cell, LevelIndex,
    levelgrid, ncells, cellindex, localindex, rawid, level

function structural_z7_ids(L::Int)
    p = 7^L
    tail = (UInt64(1) << (3 * (20 - L))) - UInt64(1)
    out = Vector{UInt64}(undef, 12 * p)
    for k in 0:(12*p-1)
        b, r = divrem(k, p)
        x = UInt64(b) << 60
        w = p
        for j in 1:L
            w ÷= 7
            d, r = divrem(r, w)
            x |= UInt64(d) << (60 - 3j)
        end
        out[k+1] = x | tail
    end
    return out
end

@testset "Z7 existence rank/select" begin
    for L in 0:4
        grid = levelgrid(IGeo7System(), L)
        n = ncells(grid)
        @test n == 10 * 7^L + 2

        # `idselect` is the level's canonical order and nothing else.
        ids = [Encodings.idselect(grid, r) for r in 0:(n-1)]
        @test ids == [rawid(cellindex(grid, i)) for i in 1:n]

        # rank inverts select over every cell of the level.
        @test all(r -> Encodings.idrank(grid, ids[r+1]) == r, 0:(n-1))
        @test all(x -> Encodings.idvalid(grid, x), ids)

        phantoms = setdiff(structural_z7_ids(L), ids)
        @test length(phantoms) == 2 * 7^L - 2
        @test !any(x -> Encodings.idvalid(grid, x), phantoms)
        # A phantom ranks where it sits, and is counted by nobody.
        @test all(p -> Encodings.idrank(grid, p) == searchsortedfirst(ids, p) - 1,
            phantoms)
        @test all(p -> Encodings.idcount_between(grid, p, p) == 0, phantoms)
        @test all(x -> Encodings.idcount_between(grid, x, x) == 1, ids)
    end
end

@testset "Z7 count_between against enumerated truth" begin
    L = 3
    grid = levelgrid(IGeo7System(), L)
    n = ncells(grid)
    ids = [Encodings.idselect(grid, r) for r in 0:(n-1)]
    phantoms = setdiff(structural_z7_ids(L), ids)

    # Deterministic probes: cells, phantoms, the ids either side of a cell, the
    # ends of the integer range, and a base cell past the twelve that exist.
    probes = sort!(unique!(UInt64[
        ids[1:37:end]...,
        phantoms[1:13:end]...,
        (ids[1:53:end] .- 1)...,
        (ids[1:53:end] .+ 1)...,
        UInt64(0), typemax(UInt64), UInt64(12) << 60, UInt64(15) << 60,
    ]))
    @test length(probes) > 100

    bad = 0
    for lo in probes, hi in probes
        truth = hi < lo ? 0 :
                searchsortedlast(ids, hi) - searchsortedfirst(ids, lo) + 1
        Encodings.idcount_between(grid, lo, hi) == max(truth, 0) || (bad += 1)
    end
    @test bad == 0
end

@testset "HEALPix nested arithmetic" begin
    grid = levelgrid(HEALPixSystem(), 2)
    n = ncells(grid)
    @test n == 12 * 4^2

    @test [Encodings.idselect(grid, r) for r in 0:(n-1)] == collect(0:(n-1))
    @test all(r -> Encodings.idrank(grid, Encodings.idselect(grid, r)) == r, 0:(n-1))
    @test Encodings.idcount_between(grid, 5, 9) == 5
    @test Encodings.idcount_between(grid, 9, 5) == 0
    # Clamping at both ends, so an interval reaching outside the level counts
    # only what is inside it.
    @test Encodings.idcount_between(grid, -3, 4) == 5
    @test Encodings.idcount_between(grid, n - 2, 10n) == 2
    @test Encodings.idvalid(grid, n - 1)
    @test !Encodings.idvalid(grid, n)
    @test !Encodings.idvalid(grid, -1)
end

@testset "encoding registry" begin
    @test Encodings.ENCODING_REGISTRY["none"] === Encodings.DenseEncoding()
    @test Encodings.ENCODING_REGISTRY["ranges"] === Encodings.RangesEncoding()
    @test Encodings.ENCODING_REGISTRY["implicit"] === Encodings.ImplicitEncoding()
    @test all(e -> e isa Encodings.CellEncoding, values(Encodings.ENCODING_REGISTRY))
    # The registry key round-trips, which is what a writer stamps into attrs.
    @test all(k -> Encodings.encodingname(Encodings.ENCODING_REGISTRY[k]) == k,
        keys(Encodings.ENCODING_REGISTRY))
end

struct GriddedRunsEncoding <: Encodings.CellEncoding end
Encodings.encodingname(::GriddedRunsEncoding) = "gridded_runs"

@testset "a downstream encoding joins the registry" begin
    # The extension point the design promises: one pair added to the table and a
    # store's vocabulary resolves to it. Kills a registry that only ever holds
    # the three shipped layouts — and one typed loosely enough to hold anything.
    @test Encodings.ENCODING_REGISTRY isa Dict{String,Encodings.CellEncoding}
    n = length(Encodings.ENCODING_REGISTRY)
    try
        @test DGG.register_encoding!("gridded_runs", GriddedRunsEncoding()) ===
              Encodings.ENCODING_REGISTRY
        @test Encodings.ENCODING_REGISTRY["gridded_runs"] === GriddedRunsEncoding()
        @test length(Encodings.ENCODING_REGISTRY) == n + 1
        # And a store that names it decodes to it, which is the whole point.
        @test DGG.encoding_for("gridded_runs") === GriddedRunsEncoding()
    finally
        delete!(Encodings.ENCODING_REGISTRY, "gridded_runs")
    end
    @test length(Encodings.ENCODING_REGISTRY) == n
    # Only encodings: the table is what every `encodingname` and `cellaxis` call
    # is dispatched from, so a bare string in it is a MethodError somewhere else.
    @test_throws MethodError DGG.register_encoding!("bogus", "ranges")
end

@testset "the merge rule chooses between rank runs and step runs" begin
    L = 3
    grid = levelgrid(IGeo7System(), L)
    n = ncells(grid)
    ids = [Encodings.idselect(grid, r) for r in 0:(n-1)]
    structural = structural_z7_ids(L)

    ranked = Encodings.idranges(grid, ids)
    stepped = Encodings.idranges(grid, ids; merge=:step)

    # Rank adjacency holds across every digit rollover and every deleted
    # branch, so the whole level is one interval. A reader that counted
    # well-formed digit strings in it would find 12*7^L, not the 10*7^L + 2
    # cells there are: the divergence is exactly the phantoms.
    @test size(ranked, 1) == 1
    @test Encodings.idcount_between(grid, ranked[1, 1], ranked[1, 2]) == n
    @test count(x -> ranked[1, 1] <= x <= ranked[1, 2], structural) == 12 * 7^L
    @test 12 * 7^L - n == 2 * 7^L - 2

    # Step adjacency breaks at both, so no interval can hold a phantom and no
    # run outlives a sibling set: the committed 504-run structure of the whole
    # res-3 earth, longest run 7.
    runlength(i) = Encodings.idcount_between(grid, stepped[i, 1], stepped[i, 2])
    structuralcount(i) = searchsortedlast(structural, stepped[i, 2]) -
                         searchsortedfirst(structural, stepped[i, 1]) + 1
    @test size(stepped, 1) == 504
    @test maximum(runlength, axes(stepped, 1)) == 7
    @test sum(runlength, axes(stepped, 1)) == n
    # The property the row count is a consequence of: on step runs a structural
    # reader and an existence reader agree, interval for interval.
    @test all(i -> structuralcount(i) == runlength(i), axes(stepped, 1))

    # Two spellings of one axis.
    fromranked = Encodings.cellaxis(Encodings.RangesEncoding(), grid, ranked;
        declared_length=n)
    fromstepped = Encodings.cellaxis(Encodings.RangesEncoding(), grid, stepped;
        declared_length=n)
    @test collect(fromranked) == collect(fromstepped)
    @test_throws ArgumentError Encodings.idranges(grid, ids; merge=:maximal)
end

@testset "write eligibility and range construction" begin
    grid = levelgrid(IGeo7System(), 4)
    runs = [201:400; 900:900; 1500:2000]
    ids = [Encodings.idselect(grid, r) for r in runs]

    @test Encodings.write_eligible(Encodings.RangesEncoding(), grid, ids)
    @test Encodings.write_eligible(Encodings.DenseEncoding(), grid, ids)
    @test !Encodings.write_eligible(Encodings.RangesEncoding(), grid, reverse(ids))
    @test !Encodings.write_eligible(Encodings.RangesEncoding(), grid,
        sort!([ids; ids[7]]))
    # A cell from another level is not a cell of this one.
    other = Encodings.idselect(levelgrid(IGeo7System(), 3), 5)
    @test !Encodings.write_eligible(Encodings.RangesEncoding(), grid,
        sort!([ids; other]))
    @test Encodings.write_eligible(Encodings.DenseEncoding(), grid, reverse(ids))

    ranges = Encodings.idranges(grid, ids)
    @test size(ranges, 2) == 2
    # Runs break where the ids stop being adjacent, so the run structure of the
    # INDICES is what the ranges array records.
    @test size(ranges, 1) == 3
    @test sum(Encodings.idcount_between(grid, ranges[i, 1], ranges[i, 2])
              for i in axes(ranges, 1)) == length(ids)
    @test Encodings.validate_ranges(grid, ranges, length(ids)) === nothing
    # The normative check: a length the ranges do not add up to.
    err = try
        Encodings.validate_ranges(grid, ranges, length(ids) + 1)
    catch e
        e
    end
    @test err isa DGGSFormatError && err.check === :count_mismatch
    # A rejection from this layer carries no store: it never saw one. The layer
    # that did opened it adds that on the way out.
    @test err.store === nothing
    @test occursin("cell axis holds", sprint(showerror, err))
    enriched = try
        with_store_context("gs://bucket/store.zarr"; conventions=["xdggs"]) do
            Encodings.validate_ranges(grid, ranges, length(ids) + 1)
        end
    catch e
        e
    end
    @test enriched isa DGGSFormatError && enriched.check === :count_mismatch
    @test enriched.store == "gs://bucket/store.zarr"
    @test occursin("conventions fired: xdggs", sprint(showerror, enriched))
    # Descending or overlapping ranges are a malformed store, not a subset.
    @test_throws DGGSFormatError Encodings.validate_ranges(grid, ranges[[2, 1, 3], :],
        length(ids))
    @test_throws DGGSFormatError Encodings.validate_ranges(grid,
        [ranges[1:1, 1] ranges[1:1, 1] .- 1], 0)
end

@testset "the compacted axis: two columns to a validated MultiOrderVector" begin
    for sys in (HEALPixSystem(), IGeo7System())
        l1 = levelgrid(sys, 1)
        a = cellindex(l1, 2)
        b = cellindex(l1, 4)
        kids = collect(DGG.children(sys, b))
        mov = DGG.MultiOrderVector(sys, [a; kids])
        lv = Int8[level(c) for c in mov]
        ids = [rawid(c) for c in mov]

        m2 = Encodings.cellaxis(CompactedEncoding(), sys, lv, ids;
            declared_length=length(mov))
        @test m2 isa DGG.MultiOrderVector
        @test collect(m2) == collect(mov) && m2 == mov

        grab(f) = try
            f()
            nothing
        catch e
            e
        end
        checkof(f) = (e = grab(f); e isa DGGSFormatError ? e.check : e)

        @test checkof(() -> Encodings.cellaxis(CompactedEncoding(), sys,
            lv[1:end-1], ids)) === :compacted_column_mismatch
        @test checkof(() -> Encodings.cellaxis(CompactedEncoding(), sys, lv, ids;
            declared_length=length(mov) + 1)) === :count_mismatch
        @test checkof(() -> Encodings.cellaxis(CompactedEncoding(), sys,
            reverse(lv), reverse(ids))) === :compacted_axis_order
        @test checkof(() -> Encodings.cellaxis(CompactedEncoding(), sys,
            [lv; lv[end]], [ids; ids[end]])) === :invalid_compacted_axis
        @test checkof(() -> Encodings.cellaxis(CompactedEncoding(), sys,
            [lv; Int8(level(b))], [ids; rawid(b)])) === :invalid_compacted_axis
        @test checkof(() -> Encodings.cellaxis(CompactedEncoding(), sys,
            Int8[99; lv[2:end]], ids)) === :invalid_stored_level
        @test checkof(() -> Encodings.cellaxis(CompactedEncoding(), sys, lv,
            [typemax(Int64); ids[2:end]])) === :id_names_no_cell
        # Validate the level before narrowing an out-of-range UInt to Int.
        @test checkof(() -> Encodings.cellaxis(CompactedEncoding(), sys,
            [typemax(UInt64); UInt64.(lv[2:end])], ids)) === :invalid_stored_level
    end
end
