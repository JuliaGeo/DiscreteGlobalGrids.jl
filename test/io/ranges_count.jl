# How a ranges store counts its intervals: on the cells that exist, or on the
# base-7 structural line of Z7 ids that keeps a slot for every digit string. The
# two agree except across a pentagon's deleted branch, and this suite pins both
# the agreement and the difference.

module DGGSIORangesCountTests

using Test
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: IGeo7System, HEALPixSystem, AuthalicSystem, Z7Cell,
    CellVector, CellLookup, Cells, DGGSFormatError, levelgrid, descendants, rawid,
    dggread
using DiscreteGlobalGrids.Encodings: structural_count_between, idcount_between,
    checkrangescount, idranges, cellaxis, RangesEncoding, RANGES_COUNTS
import DimensionalData as DD

const I = DGG.IGeo7
const SYS = IGeo7System()

function thrown(f)
    try
        f()
        return nothing
    catch e
        return e
    end
end

@testset "ranges count" begin

    @testset "the base-7 structural line" begin
        for r in (0, 1, 3, 7, 12, 19)
            g = levelgrid(SYS, r)
            n = DGG.ncells(g)
            for k in unique(round.(Int, range(1, n; length=min(n, 400))))
                z = rawid(DGG.cellindex(g, k))
                m = I.z7_to_monotonic(z)
                @test m < I.MONOTONIC_LIMIT[r+1]
                @test I.z7_from_monotonic(m, r) == z
                @test I.z7_to_monotonic(z, r) == m
            end
        end
        # a parent owns one block of 7^k positions at every finer level
        p = Z7Cell("0941054")
        mp = I.z7_to_monotonic(rawid(p))
        for c in descendants(SYS, p, 7)
            @test I.z7_to_monotonic(rawid(c)) ÷ UInt64(7)^2 == mp
        end
        # order agrees with the packed ids
        ids = sort!(rawid.(collect(descendants(SYS, Z7Cell("0812"), 6))))
        @test issorted(I.z7_to_monotonic.(ids))
        # phantoms have positions too, and are not cells
        phantom = I.z7_from_monotonic(UInt64(2), 1)          # base 0, deleted digit 2
        @test !I.is_valid_z7(phantom)
        @test I.z7_to_monotonic(phantom, 1) == 2
        e = thrown(() -> I.z7_to_monotonic(rawid(p), 4))
        @test e isa I.InvalidZ7Error && e.reason === :monotonic_res
        e = thrown(() -> I.z7_from_monotonic(I.MONOTONIC_LIMIT[3], 2))
        @test e isa I.InvalidZ7Error && e.reason === :monotonic_range
        @test sprint(showerror, e) isa String
    end

    @testset "structural vs existence count" begin
        for L in (1, 3, 5)
            g = levelgrid(SYS, L)
            lo = DGG.idselect(g, 0)
            hi = DGG.idselect(g, DGG.ncells(g) - 1)
            # the whole level: every digit string vs the cells that exist
            @test structural_count_between(g, lo, hi) == 12 * 7^L
            @test idcount_between(g, lo, hi) == 10 * 7^L + 2
            # one row per interval, against the place-value difference
            for (a, b) in ((0, 41), (17, 300), (5, DGG.ncells(g) - 1))
                b < DGG.ncells(g) || continue
                x, y = DGG.idselect(g, a), DGG.idselect(g, b)
                @test structural_count_between(g, x, y) ==
                      Int(I.z7_to_monotonic(y) - I.z7_to_monotonic(x)) + 1
            end
            @test structural_count_between(g, hi, lo) == 0
            @test structural_count_between(g, lo, typemax(UInt64)) == 12 * 7^L
        end
        # off the pentagon chains both count the same
        g = levelgrid(SYS, 4)
        cells = sort!(collect(descendants(SYS, Z7Cell("0812"), 4)))
        x, y = rawid(first(cells)), rawid(last(cells))
        @test structural_count_between(g, x, y) == idcount_between(g, x, y) == 7^2
        # the authalic wrapper and HEALPix answer the same arithmetic
        @test structural_count_between(levelgrid(AuthalicSystem(SYS), 4), x, y) == 49
        hg = levelgrid(HEALPixSystem(), 2)
        @test structural_count_between(hg, 3, 40) == idcount_between(hg, 3, 40) == 38
    end

    @testset "checkrangescount" begin
        g = levelgrid(SYS, 3)
        whole = idranges(g, [DGG.cellindex(g, k) for k in 1:DGG.ncells(g)]; merge=:rank)
        @test size(whole, 1) == 1
        starts, stops = whole[:, 1], whole[:, 2]
        total = idcount_between(g, starts[1], stops[1])
        @test total == 3432
        # existence: the plain length check
        @test checkrangescount(g, starts, stops, total, 3432, :existence) === nothing
        @test checkrangescount(g, starts, stops, total, 3432, nothing) === nothing
        e = thrown(() -> checkrangescount(g, starts, stops, total, 3000, :existence))
        @test e isa DGGSFormatError && e.check === :count_mismatch
        # structural with phantom slots inside: not readable yet, said so
        e = thrown(() -> checkrangescount(g, starts, stops, total, 4116, :structural))
        @test e isa DGGSFormatError && e.check === :structural_phantom_slots
        @test e.observed == 3432
        # undeclared, and the length is the structural count: the same cause
        e = thrown(() -> checkrangescount(g, starts, stops, total, 4116, nothing))
        @test e isa DGGSFormatError && e.check === :structural_phantom_slots
        # undeclared, and neither count: a plain mismatch
        e = thrown(() -> checkrangescount(g, starts, stops, total, 4000, nothing))
        @test e isa DGGSFormatError && e.check === :count_mismatch
        # step rows never cross a phantom, so structural reads them
        step = idranges(g, [DGG.cellindex(g, k) for k in 1:DGG.ncells(g)]; merge=:step)
        @test checkrangescount(g, step[:, 1], step[:, 2], 3432, 3432, :structural) === nothing
        @test_throws ArgumentError checkrangescount(g, starts, stops, total, 3432, :slots)
        # the axis constructor runs the same check
        e = thrown(() -> cellaxis(RangesEncoding(), g, whole; declared_length=4116,
            ranges_count=:structural))
        @test e isa DGGSFormatError && e.check === :structural_phantom_slots
        @test length(cellaxis(RangesEncoding(), g, step; declared_length=3432,
            ranges_count=:structural)) == 3432
        @test RANGES_COUNTS == (:existence, :structural)
    end

    ZARR_LOADED = try
        @eval using Zarr
        true
    catch
        false
    end

    if !ZARR_LOADED
        @info "Zarr.jl is unavailable: the ranges_count round trips are skipped."
    else
        @testset "ranges_count is declared on write and read back" begin
            L = 3
            cells = sort!(reduce(vcat, [collect(descendants(SYS, Z7Cell(c), L))
                                        for c in ("0812", "0911", "1004")]))
            stack = DD.DimStack((elevation=Float32.(1:length(cells)),),
                (Cells(CellLookup(CellVector(SYS, L, cells))),))
            for (merge, count) in ((:step, "structural"), (:rank, "existence"))
                path = joinpath(mktempdir(), "rc_$merge.zarr")
                DGG.dggwrite(path, stack; merge)
                z = Zarr.zopen(path)
                @test z.attrs["dggs"]["ranges_count"] == count
                back = dggread(path)
                @test DD.metadata(back)["description"].ranges_count === Symbol(count)
                @test collect(DD.lookup(back, DD.dims(back)[1])) == cells
            end
            # a dense store says nothing about ranges
            path = joinpath(mktempdir(), "rc_dense.zarr")
            DGG.dggwrite(path, stack; encoding=:dense)
            @test !haskey(Zarr.zopen(path).attrs["dggs"], "ranges_count")
        end

        @testset "an unknown ranges_count is refused" begin
            cells = sort!(collect(descendants(SYS, Z7Cell("0812"), 3)))
            stack = DD.DimStack((elevation=Float32.(1:length(cells)),),
                (Cells(CellLookup(CellVector(SYS, 3, cells))),))
            path = joinpath(mktempdir(), "rc_bad.zarr")
            DGG.dggwrite(path, stack)
            f = joinpath(path, ".zattrs")
            s = read(f, String)
            t = replace(s, r"(\"ranges_count\":\s*)\"structural\"" => s"\g<1>\"slots\"")
            @test t != s
            write(f, t)
            e = thrown(() -> dggread(path))
            @test e isa DGGSFormatError && e.check === :unknown_ranges_count
        end
    end
end

end # module DGGSIORangesCountTests
