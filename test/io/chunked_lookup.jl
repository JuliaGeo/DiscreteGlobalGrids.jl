# The chunk manifest, the computed and chunk-cached axis vectors, and the
# two-level selector resolution over them.
#
# Two independent implementations are checked against each other throughout: the
# arithmetic axis built from ranges with no IO, and the scanned axis built from
# stored ids. The oracle for selection is the package's own `CellLookup` over the
# same cells, which resolves selectors by a route this file shares nothing with.

using Test
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: IGeo7System, HEALPixSystem, Z7Cell, LevelIndex,
    levelgrid, ncells, cellindex, rawid, level, ancestor,
    cell_extent, CellVector, CellLookup, Cells, Covering
import DimensionalData as DD
using DimensionalData: Lookups

# An AbstractVector that counts block reads, standing in for a lazy DiskArray:
# what it proves is that the manifest prunes, i.e. that a selector touches one
# chunk and not the array.
struct CountingIds <: AbstractVector{UInt64}
    values::Vector{UInt64}
    reads::Base.RefValue{Int}
end
CountingIds(v::Vector{UInt64}) = CountingIds(v, Ref(0))
Base.size(c::CountingIds) = size(c.values)
Base.getindex(c::CountingIds, i::Int) = (c.reads[] += 1; c.values[i])
Base.getindex(c::CountingIds, r::AbstractUnitRange) = (c.reads[] += 1; c.values[r])

# A few hundred cells at IGEO7 level 4 in three index runs: the shape of a
# regional store, small enough to hold the whole truth in memory.
axis_grid = levelgrid(IGeo7System(), 4)
axis_runs = [201:400; 900:900; 1500:2000]
axis_ids = [Encodings.idselect(axis_grid, r) for r in axis_runs]
axis_cells = [Z7Cell(x) for x in axis_ids]
axis_chunk = 64
axis_ranges = Encodings.idranges(axis_grid, axis_ids)
# Base cell 0's pentagon chain deletes digit 2, so `00002` is well formed and
# names nothing. It sits between the level-4 cells of ranks 1 and 2.
axis_phantom = (UInt64(2) << 48) | ((UInt64(1) << 48) - 1)

@testset "a range may span an id that names no cell" begin
    ids = [Encodings.idselect(axis_grid, r) for r in 0:5]
    @test ids[2] < axis_phantom < ids[3]
    @test !Encodings.idvalid(axis_grid, axis_phantom)

    ranges = Encodings.idranges(axis_grid, ids)
    @test size(ranges, 1) == 1
    # Six cells in the interval, not the seven integers a structural count of
    # well-formed digit strings would find.
    @test Encodings.idcount_between(axis_grid, ranges[1, 1], ranges[1, 2]) == 6

    axis = Encodings.cellaxis(Encodings.RangesEncoding(), axis_grid, ranges;
        declared_length=6)
    @test collect(axis) == [Z7Cell(x) for x in ids]
    @test ChunkedLookups.axisindex(axis, axis_phantom) === nothing
end

@testset "a ranges row that holds no cell is skipped on read" begin
    # A well-formed interval containing only a deleted branch: structurally a
    # row, arithmetically empty. Indices must step straight over it.
    ids = [Encodings.idselect(axis_grid, r) for r in 0:5]
    ranges = UInt64[ids[1] ids[2]; axis_phantom axis_phantom; ids[3] ids[6]]
    axis = Encodings.cellaxis(Encodings.RangesEncoding(), axis_grid, ranges;
        declared_length=6)
    @test collect(axis) == [Z7Cell(x) for x in ids]
    @test ChunkedLookups.axisindex(axis, ids[3]) == 3
    @test ChunkedLookups.axisindex(axis, axis_phantom) === nothing
end

@testset "a manifest whose chunk length does not describe it is refused" begin
    # `chunkof` divides by `chunklength` rather than searching `offsets`, so a
    # manifest loaded from a sidecar with a foreign chunk length would resolve
    # every index into the wrong chunk. Only the last chunk may be short.
    ok = ChunkedLookups.ChunkManifest(UInt64[1, 5], UInt64[4, 8], [4, 3], [0, 4], 4)
    @test ChunkedLookups.nchunks(ok) == 2
    @test length(ok) == 7
    @test_throws DGGSFormatError ChunkedLookups.ChunkManifest(
        UInt64[1, 5], UInt64[4, 8], [4, 4], [0, 4], 5)
    @test_throws DGGSFormatError ChunkedLookups.ChunkManifest(
        UInt64[1, 5, 9], UInt64[4, 8], [4, 4], [0, 4], 4)
end

@testset "the axis is the id vector, however it is stored" begin
    ranged = Encodings.cellaxis(Encodings.RangesEncoding(), axis_grid, axis_ranges;
        declared_length=length(axis_ids))
    dense = Encodings.cellaxis(Encodings.DenseEncoding(), axis_grid, axis_ids;
        chunklength=axis_chunk)

    @test length(ranged) == length(axis_ids)
    @test collect(ranged) == axis_cells
    @test collect(dense) == axis_cells
    # The inverse, on both, including ids the axis does not hold.
    @test all(k -> ChunkedLookups.axisindex(ranged, axis_ids[k]) == k,
        eachindex(axis_ids))
    @test all(k -> ChunkedLookups.axisindex(dense, axis_ids[k]) == k,
        eachindex(axis_ids))
    absent = Encodings.idselect(axis_grid, 700)          # between two runs
    @test ChunkedLookups.axisindex(ranged, absent) === nothing
    @test ChunkedLookups.axisindex(dense, absent) === nothing
end

@testset "manifest from ranges == manifest from a dense scan" begin
    ranged = Encodings.cellaxis(Encodings.RangesEncoding(), axis_grid, axis_ranges;
        declared_length=length(axis_ids))
    dense = Encodings.cellaxis(Encodings.DenseEncoding(), axis_grid, axis_ids;
        chunklength=axis_chunk)

    fromranges = ChunkedLookups.ChunkManifest(ranged, axis_chunk)
    fromscan = ChunkedLookups.ChunkManifest(dense, axis_chunk)
    @test fromranges == fromscan

    n = length(axis_ids)
    nc = cld(n, axis_chunk)
    @test ChunkedLookups.nchunks(fromranges) == nc
    @test fromranges.firstids == [axis_ids[(c-1)*axis_chunk+1] for c in 1:nc]
    @test fromranges.lastids == [axis_ids[min(c * axis_chunk, n)] for c in 1:nc]
    @test sum(fromranges.lengths) == n
    @test fromranges.lengths[end] == n - (nc - 1) * axis_chunk
    @test fromranges.offsets == [(c - 1) * axis_chunk for c in 1:nc]

    # A data array chunked differently from the coordinate gets its own
    # manifest, from the same arithmetic and still without touching the store.
    other = ChunkedLookups.ChunkManifest(ranged, 97)
    @test ChunkedLookups.nchunks(other) == cld(n, 97)
    @test sum(other.lengths) == n
    @test other.firstids[3] == axis_ids[2*97+1]
end

@testset "implicit axis needs no store at all" begin
    grid = levelgrid(HEALPixSystem(), 2)
    axis = Encodings.cellaxis(Encodings.ImplicitEncoding(), grid, ncells(grid))
    @test length(axis) == ncells(grid)
    @test collect(axis) == [cellindex(grid, i) for i in 1:ncells(grid)]
    manifest = ChunkedLookups.ChunkManifest(axis, 50)
    @test manifest.lengths == [50, 50, 50, 42]
    @test manifest.firstids == [0, 50, 100, 150]
    @test manifest.lastids == [49, 99, 149, 191]
end

@testset "the dense scan is where a broken axis is caught" begin
    duplicated = copy(axis_ids)
    duplicated[130] = duplicated[129]
    duplicated[300] = duplicated[299]
    err = try
        Encodings.cellaxis(Encodings.DenseEncoding(), axis_grid, duplicated;
            chunklength=axis_chunk)
        nothing
    catch e
        e
    end
    @test err isa DGGSFormatError && err.check === :duplicate_ids
    # The count and the FIRST offending index, both, so the message locates
    # the damage rather than merely reporting it.
    msg = sprint(showerror, err)
    @test occursin("2", msg)
    @test occursin("130", msg)
    @test occursin("duplicate", msg)

    unsorted = copy(axis_ids)
    unsorted[10], unsorted[11] = unsorted[11], unsorted[10]
    @test_throws DGGSFormatError Encodings.cellaxis(Encodings.DenseEncoding(),
        axis_grid, unsorted; chunklength=axis_chunk)

    # A well-formed id that names no cell (a pentagon's deleted branch) is not a
    # cell of the level and must not pass ingest. It is placed where it keeps
    # the array sorted, so validity and not order is what rejects it.
    phantoms = copy(axis_ids)
    phantoms[1] = axis_phantom
    err = try
        Encodings.cellaxis(Encodings.DenseEncoding(), axis_grid, phantoms;
            chunklength=axis_chunk)
        nothing
    catch e
        e
    end
    @test err isa DGGSFormatError && err.check === :id_names_no_cell
    @test occursin("names no cell", sprint(showerror, err))

    # The manifest publishes each chunk's first and last id, so those two in
    # particular must be cells. Here the last slot of the first chunk holds a
    # level-5 id — the "attrs lie about the level" failure — placed so that
    # order and uniqueness still hold and only validity can reject it.
    mislevelled = copy(axis_ids)
    mislevelled[axis_chunk] = axis_ids[axis_chunk] & ~(UInt64(1) << 47)
    @test issorted(mislevelled) && allunique(mislevelled)
    @test !Encodings.idvalid(axis_grid, mislevelled[axis_chunk])
    err = try
        Encodings.cellaxis(Encodings.DenseEncoding(), axis_grid, mislevelled;
            chunklength=axis_chunk)
        nothing
    catch e
        e
    end
    @test err isa DGGSFormatError && err.check === :id_names_no_cell
    @test occursin("names no cell", sprint(showerror, err))

    # Zarr pads the final chunk with the fill value; the declared length, not
    # the stored one, says where the axis ends.
    n = length(axis_ids)
    padded = [axis_ids; zeros(UInt64, axis_chunk * cld(n, axis_chunk) - n)]
    axis = Encodings.cellaxis(Encodings.DenseEncoding(), axis_grid, padded;
        chunklength=axis_chunk, declared_length=n)
    @test length(axis) == n
    @test last(axis) == Z7Cell(last(axis_ids))
    @test ChunkedLookups.ChunkManifest(axis, axis_chunk).lastids[end] == last(axis_ids)
end

@testset "the dense scan verifies every id, not a sample" begin
    # One chunk of 100 cells. Sampling three ids visits slots 1, 50 and 100, so
    # slot 3 is exactly where a phantom hides from it. It is the deleted-branch
    # id between the cells of ranks 1 and 2, so the array stays sorted.
    ids = [Encodings.idselect(axis_grid, r) for r in 0:99]
    ids[3] = axis_phantom
    @test issorted(ids) && allunique(ids)

    err = try
        Encodings.cellaxis(Encodings.DenseEncoding(), axis_grid, ids;
            chunklength=100)
        nothing
    catch e
        e
    end
    @test err isa DGGSFormatError && err.check === :id_names_no_cell
    @test occursin("index 3", sprint(showerror, err))

    # Sampling is the explicit opt-out, and is what misses it.
    sampled = Encodings.cellaxis(Encodings.DenseEncoding(), axis_grid, ids;
        chunklength=100, samples=3)
    @test length(sampled) == 100
end

@testset "two-level resolution answers what a full search answers" begin
    naive = CellLookup(CellVector(IGeo7System(), 4, axis_cells))
    ranged = ChunkedLookups.ChunkedCellLookup(
        Encodings.cellaxis(Encodings.RangesEncoding(), axis_grid, axis_ranges;
            declared_length=length(axis_ids)))
    counting = CountingIds(copy(axis_ids))
    dense = ChunkedLookups.ChunkedCellLookup(
        Encodings.cellaxis(Encodings.DenseEncoding(), axis_grid, counting;
            chunklength=axis_chunk))

    hits = (DD.At(axis_cells[1]), DD.At(axis_cells[257]), DD.At(axis_cells[end]))
    misses = (DD.At(Z7Cell(Encodings.idselect(axis_grid, 700))),
        DD.At(Z7Cell(Encodings.idselect(axis_grid, 0))))

    for lookup in (ranged, dense)
        @test length(lookup) == length(axis_ids)
        @test Lookups.order(lookup) isa Lookups.ForwardOrdered
        for sel in hits
            @test Lookups.selectindices(lookup, sel) ==
                  Lookups.selectindices(naive, sel)
            @test Lookups.hasselection(lookup, sel)
        end
        for sel in misses
            @test !Lookups.hasselection(lookup, sel)
            @test_throws Lookups.SelectorError Lookups.selectindices(lookup, sel)
        end
        # A stored axis reports itself in prose, not as its parameterised type,
        # and names the level mismatch when there is one.
        @test_throws "not the axis's level 4" Lookups.selectindices(lookup,
            DD.At(ancestor(IGeo7System(), axis_cells[1], 3)))
        message = try
            Lookups.selectindices(lookup, first(misses))
        catch err
            sprint(showerror, err)
        end
        @test !occursin("ChunkedCellLookup{", message)
        @test occursin("ChunkedCellLookup(IGeo7System, level=4", message)
        # `Near` is refused here for the same reason it is on `CellLookup`.
        @test_throws "nearest id is not the nearest cell" Lookups.selectindices(
            lookup, DD.Near(axis_cells[1]))
    end

    # One chunk read per resolved selector: the manifest prunes, so the axis is
    # never scanned.
    counting.reads[] = 0
    Lookups.selectindices(dense, hits[2])
    @test counting.reads[] == 1

    # A region: the coverage of a coarse ancestor of one of the stored cells.
    target = cell_extent(levelgrid(IGeo7System(), 2),
        ancestor(IGeo7System(), axis_cells[300], 2))
    expected = Lookups.selectindices(naive, Covering(target))
    @test !isempty(expected)
    for lookup in (ranged, dense)
        @test Lookups.selectindices(lookup, Covering(target)) == expected
    end
end

@testset "the stored axis is a cube dimension" begin
    lookup = ChunkedLookups.ChunkedCellLookup(
        Encodings.cellaxis(Encodings.RangesEncoding(), axis_grid, axis_ranges;
            declared_length=length(axis_ids)))
    data = collect(1.0:length(axis_ids))
    cube = DD.DimArray(data, Cells(lookup))
    @test DD.lookup(cube, Cells) === lookup
    @test cube[Cells(DD.At(axis_cells[257]))] == 257.0

    target = cell_extent(levelgrid(IGeo7System(), 2),
        ancestor(IGeo7System(), axis_cells[300], 2))
    expected = Lookups.selectindices(lookup, Covering(target))
    sub = cube[Cells(Covering(target))]
    @test parent(sub) == data[expected]
    # A subset is no longer a stored axis: it becomes the package's own
    # compressed lookup, naming the same cells.
    @test DD.lookup(sub, Cells) isa CellLookup
    @test collect(DD.lookup(sub, Cells)) == axis_cells[expected]
end
