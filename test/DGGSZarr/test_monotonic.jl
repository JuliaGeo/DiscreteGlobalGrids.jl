# test_monotonic.jl — the base-7 number line and the lazy id vectors on it
# (src/IGeo7/z7_ranges.jl). No data on disk is touched: every id here is built
# from the Z7 codecs, so this is the tier that runs everywhere.
#
# The single most important assertion in this file is that `z7_to_monotonic`
# and `cell_to_index` are *not* the same function. They agree on order, which
# is why a spot check on a contiguous sample passes under either, and they
# disagree on spacing, which is what makes a range table expand to the wrong
# cell count when the wrong one is used. Everything else here is round trips.

using Test

module TestMonotonic

using Test

using DimensionalData

using DiscreteGlobalGrids: Helpers, IGeo7
using DiscreteGlobalGrids.IGeo7: Z7CachedIds, Z7RangeIds, cell_to_index, index_to_cell,
    is_valid_z7, num_cells, z7_from_monotonic, z7_from_string, z7_level, z7_materialize!,
    z7_nranges, z7_range_position, z7_resolution, z7_to_monotonic, z7_cached_position,
    z7_is_materialized
import DiscreteGlobalGrids as DGG

"Every cell of resolution `res`, in canonical ascending order."
all_cells(res) = [index_to_cell(i, res) for i in 1:num_cells(res)]

@testset "z7_to_monotonic / z7_from_monotonic" begin
    @testset "the published example" begin
        # base 8, digits 0,0,4,3,3 -> the place value read in base 7.
        z = z7_from_string("0800433")
        expected = ((((8 * 7 + 0) * 7 + 0) * 7 + 4) * 7 + 3) * 7 + 3
        @test z7_to_monotonic(z, 5) == UInt64(expected)
        @test z7_to_monotonic(z) == UInt64(expected)
        @test z7_from_monotonic(UInt64(expected), 5) == z
    end

    @testset "round trip over whole resolutions" begin
        for res in 0:4
            cells = all_cells(res)
            values = z7_to_monotonic.(cells, res)
            @test all(z7_from_monotonic.(values, res) .== cells)
            # Strictly increasing: the line preserves canonical id order, which
            # is the property the ranges form is built on.
            @test Helpers.strictly_increasing(values)
            # ... but it is *not* dense over the cells: it indexes all 12·7^res
            # digit strings, the deleted-digit ones included.
            @test last(values) - first(values) + 1 >= UInt64(length(cells))
        end
    end

    @testset "distinct from cell_to_index" begin
        # Order agrees at every resolution — this is why the two are confusable.
        for res in 0:4
            cells = all_cells(res)
            @test sortperm(z7_to_monotonic.(cells, res)) == sortperm(cell_to_index.(cells))
        end
        # Spacing does not. At res 2 the grid has 10·49+2 = 492 cells while the
        # line has 12·49 = 588 positions, so the two differ on all but a prefix.
        cells = all_cells(2)
        @test num_cells(2) == 492
        monotonic = z7_to_monotonic.(cells, 2)
        ranks = cell_to_index.(cells)
        @test monotonic .+ 1 != ranks
        # Concretely: expanding [first, last] on the number line yields more
        # positions than the grid has cells, and that gap is exactly what would
        # be lost by expanding a range table with `cell_to_index`.
        @test Int(last(monotonic) - first(monotonic)) + 1 > num_cells(2)
    end

    @testset "rejections" begin
        z = z7_from_string("0800433")
        @test_throws IGeo7.InvalidZ7Error z7_to_monotonic(z, 4)     # wrong resolution
        @test_throws IGeo7.InvalidZ7Error z7_to_monotonic(z, 21)    # out of range
        @test_throws IGeo7.InvalidZ7Error z7_from_monotonic(UInt64(12 * 7^2), 2)
        @test z7_from_monotonic(UInt64(12 * 7^2 - 1), 2) isa UInt64
    end

    @testset "positions with no cell" begin
        # The line has more positions than the grid has cells; those extra
        # positions decode to well-formed digit strings that are not cells.
        # Callers that need cells validate — this documents that they must.
        invalid = filter(v -> !is_valid_z7(z7_from_monotonic(UInt64(v), 2)), 0:(12*49-1))
        @test !isempty(invalid)
        @test length(invalid) == 12 * 49 - num_cells(2)
    end
end

@testset "Z7RangeIds" begin
    res = 3
    cells = all_cells(res)

    "Build the (starts, ends) range table covering `subset` (sorted cells)."
    function table_for(subset, res)
        values = z7_to_monotonic.(subset, res)
        breaks = findall(i -> values[i] != values[i-1] + 1, 2:length(values)) .+ 1
        starts_idx = [1; breaks]
        ends_idx = [breaks .- 1; length(values)]
        return subset[starts_idx], subset[ends_idx]
    end

    @testset "a contiguous whole resolution collapses to few ranges" begin
        starts, ends = table_for(cells, res)
        ids = Z7RangeIds(starts, ends, res)
        @test length(ids) == length(cells)
        @test collect(ids) == cells
        # The grid's deleted digits are the only breaks, so R is tiny next to N.
        @test z7_nranges(ids) < length(cells)
        @test z7_level(ids) == res
    end

    @testset "a scattered subset round trips" begin
        subset = cells[1:7:end]
        starts, ends = table_for(subset, res)
        ids = Z7RangeIds(starts, ends, res)
        @test collect(ids) == subset
        @test length(ids) == length(subset)
        # getindex and cell_position are mutual inverses over the whole vector.
        @test all(i -> z7_range_position(ids, ids[i]) == i, eachindex(ids))
        # Every cell of the resolution that is *not* in the subset is absent —
        # the gaps between ranges are real, not rounded into a neighbour.
        absent = setdiff(cells, subset)
        @test all(id -> z7_range_position(ids, id) === nothing, absent)
    end

    @testset "matrix constructor and its transposition" begin
        subset = cells[1:11:end]
        starts, ends = table_for(subset, res)
        matrix = hcat(starts, ends)                      # R × 2, the on-disk layout
        @test collect(Z7RangeIds(matrix, res)) == subset
        # A 2 × R matrix is the transposed (Zarr-order) form and must be rejected
        # rather than silently read as two ranges.
        @test_throws ArgumentError Z7RangeIds(permutedims(matrix), res)
    end

    @testset "rejections" begin
        a, b = cells[1], cells[5]
        @test_throws ArgumentError Z7RangeIds([a, b], [b], res)          # length mismatch
        @test_throws ArgumentError Z7RangeIds([b], [a], res)             # ends before starts
        @test_throws ArgumentError Z7RangeIds([a, a], [b, b], res)       # overlapping
        @test_throws IGeo7.InvalidZ7Error Z7RangeIds([a], [b], 25)       # bad level
    end

    @testset "state is O(R), not O(N)" begin
        starts, ends = table_for(cells, res)
        ids = Z7RangeIds(starts, ends, res)
        # Two vectors of R (+1) words and an Int — nothing proportional to N.
        bytes = Base.summarysize(ids)
        @test bytes < 64 * (z7_nranges(ids) + 2)
        @test bytes < 8 * length(ids)
    end

    @testset "show does not materialize" begin
        starts, ends = table_for(cells, res)
        ids = Z7RangeIds(starts, ends, res)
        text = sprint(show, MIME"text/plain"(), ids)
        @test occursin("Z7RangeIds", text)
        @test occursin("ranges=$(z7_nranges(ids))", text)
        @test !occursin("0x", text)          # no ids were printed
    end
end

@testset "Z7CachedIds" begin
    res = 3
    cells = all_cells(res)

    @testset "reads only when asked" begin
        ids = Z7CachedIds(cells, res)
        @test !z7_is_materialized(ids)
        @test length(ids) == length(cells)
        # The two endpoints are the documented exception: they are what the
        # lookup constructor and DimensionalData's compact display probe.
        @test ids[1] == cells[1]
        @test ids[end] == cells[end]
        @test !z7_is_materialized(ids)
        @test z7_cached_position(ids, cells[42]) == 42
        @test z7_is_materialized(ids)        # the position query pulled them
        @test z7_materialize!(ids) == cells
    end

    @testset "an interior element pulls the whole array" begin
        # Not an optimization detail but the contract: serving interior indices
        # one at a time from a chunked store is what cost 49 GiB on the res-12
        # archive, so the second element is already a full read.
        ids = Z7CachedIds(cells, res)
        @test ids[2] == cells[2]
        @test z7_is_materialized(ids)
        # ... and iterating is therefore one read, not `length(ids)` of them.
        @test collect(Z7CachedIds(cells, res)) == cells
    end

    @testset "validates on materialization, not before" begin
        # Descending ids: construction is fine (nothing was read), the pull fails.
        bad = Z7CachedIds(reverse(cells), res)
        @test bad isa Z7CachedIds
        @test_throws ArgumentError z7_materialize!(bad)

        # Mixed resolutions, same deal.
        mixed = Z7CachedIds([cells[1], IGeo7.cell_to_parent(cells[2])], res, 2)
        @test_throws IGeo7.InvalidZ7Error z7_materialize!(mixed)
    end

    @testset "absent ids answer nothing" begin
        subset = cells[1:2:end]
        ids = Z7CachedIds(subset, res)
        @test all(i -> z7_cached_position(ids, subset[i]) == i, eachindex(subset))
        @test all(id -> z7_cached_position(ids, id) === nothing, setdiff(cells, subset))
    end
end

@testset "lookup integration" begin
    res = 3
    cells = all_cells(res)
    subset = cells[1:5:end]
    values = z7_to_monotonic.(subset, res)
    breaks = findall(i -> values[i] != values[i-1] + 1, 2:length(values)) .+ 1
    ids = Z7RangeIds(subset[[1; breaks]], subset[[breaks .- 1; length(subset)]], res)

    lookup = IGeo7.IGeo7Lookups.IGeo7Lookup(ids)
    @test length(lookup) == length(subset)
    @test lookup.resolution == res

    # `cell_position` dispatches to the range search rather than the generic
    # binary search over stored ids — the wiring this whole layer rests on.
    @test DGG.cell_position(ids, subset[3]) == 3
    @test DGG.cell_position(ids, setdiff(cells, subset)[1]) === nothing

    # `rebuild` must not materialize: the returned lookup still wraps the lazy
    # ids, not a Vector. (Without the constructor fast path this silently
    # expands the whole dimension.)
    rebuilt = DimensionalData.rebuild(lookup)
    @test DimensionalData.parent(rebuilt) isa Z7RangeIds
end

end # module TestMonotonic
