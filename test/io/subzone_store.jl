# The ancestor-subzone store on disk: what the two-dimensional bytes are, what
# the incremental writer touches, and that the lazy cell-axis view hands back
# exactly what went in.
#
# The layout's arithmetic is `subzones.jl`'s and is not repeated. What is under
# test here is the store: the array shapes and chunk grid read back with Zarr
# directly rather than through `dggread`, the fill semantics of a column nobody
# wrote, the claim that a column write touches one file, and the round trip
# through `dggread`'s DiskArrays view — including a regrid straight into it.
#
# The suite self-skips when Zarr.jl is absent, as the other extension suites do.

module DGGSubzoneStoreTests

using Test
import DiscreteGlobalGrids as DGG

const HAS_ZARR = try
    @eval using Zarr
    true
catch err
    @warn "Zarr.jl is not loadable: the subzone-store suite is skipped." exception = err
    false
end

if HAS_ZARR

import DimensionalData as DD
import DiskArrays
import Extents
using DiscreteGlobalGrids: IGeo7System, Cells, CellLookup, CellVector,
    DGGSFormatError, SubzoneLayout, cellindex, children, columncell, columnindex,
    columnlength, columnpositions, dggread, dggwrite, dggwrite!, levelgrid,
    ncells, rootcells, subzone_cellvector, subzonestore

const SYS = IGeo7System()
const LEVEL = 4
const ANCESTOR = 2
const LAYOUT = SubzoneLayout(SYS, LEVEL, ANCESTOR)
const GRID = levelgrid(SYS, LEVEL)
const CAPACITY = 49                     # 7^2, spelled out rather than imported
const PENT = 41                         # (5*7^2 + 1) / 6

# Three hexagon-rooted columns away from any pentagon chain, and the pentagon
# under root cell 0 — the two cases the padding rule distinguishes.
const HEXCOLS = [5, 6, 7]
const PENTCOL = columnindex(LAYOUT,
    foldl((c, _) -> first(children(SYS, c)), 1:ANCESTOR; init=first(rootcells(SYS))))

@assert columnlength(LAYOUT, PENTCOL) == PENT
@assert all(i -> columnlength(LAYOUT, i) == CAPACITY, HEXCOLS)

dest(name) = joinpath(mktempdir(), name)

"A cube over the cells of `cols`, values ascending from 1."
function democube(cols; name=:elevation)
    cv = subzone_cellvector(LAYOUT, cols)
    return DD.DimArray(Float32.(1:length(cv)), (Cells(CellLookup(cv)),); name=name)
end

"`f()`'s exception, or `nothing` when it returns."
caught(f) = try
    f()
    nothing
catch e
    e
end

# ---------------------------------------------------------------------------

@testset "the bytes a one-shot write leaves" begin
    cube = democube(HEXCOLS)
    path = dest("hex.zarr")
    @test dggwrite(path, cube; layout=:subzones, ancestor_level=ANCESTOR) == path

    g = Zarr.zopen(path)
    z = g["elevation"]
    # Julia sees (subzone, ancestor) and the store declares (ancestor, subzone):
    # one chunk is one whole column, which is the point of the layout.
    @test size(z) == (CAPACITY, LAYOUT.ncolumns)
    @test z.metadata.chunks == (CAPACITY, 1)
    @test isnan(z.metadata.fill_value)
    @test z.attrs["_ARRAY_DIMENSIONS"] == ["ancestor", "subzone"]

    block = g.attrs["dggs"]["subzone_layout"]
    @test g.attrs["dggs"]["name"] == "igeo7"
    @test g.attrs["dggs"]["refinement_level"] == LEVEL
    @test block["layout"] == "ancestor_subzone"
    @test block["ancestor_level"] == ANCESTOR
    @test block["subzone_count"] == CAPACITY
    @test block["ancestor_count"] == LAYOUT.ncolumns
    @test block["variables"] == ["elevation"]
    @test block["ancestor_coordinate"] == "ancestor_cell_ids"

    # The ancestor coordinate is the level-2 ids, in order, and carries the
    # xdggs spelling of what it is.
    ids = g["ancestor_cell_ids"][:]
    @test length(ids) == LAYOUT.ncolumns
    @test ids == [DGG.rawid(columncell(LAYOUT, i)) for i in 1:LAYOUT.ncolumns]
    @test g["ancestor_cell_ids"].attrs["level"] == ANCESTOR

    # The written columns hold the cube, column by column, and nothing else was
    # written at all.
    values = DD.data(cube)
    for (k, col) in pairs(HEXCOLS)
        @test z[:, col] == values[((k-1)*CAPACITY+1):(k*CAPACITY)]
    end
    @test all(isnan, z[:, 1])
    @test all(isnan, z[:, 200])
end

@testset "a pentagon column is padded, not truncated" begin
    cube = democube([PENTCOL, HEXCOLS[1]])
    path = dest("pent.zarr")
    dggwrite(path, cube; layout=:subzones, ancestor_level=ANCESTOR)

    z = Zarr.zopen(path)["elevation"]
    column = z[:, PENTCOL]
    @test length(column) == CAPACITY
    @test column[1:PENT] == DD.data(cube)[1:PENT]
    @test all(isnan, column[(PENT+1):CAPACITY])

    # And the axis that comes back has the padding dropped: p(d) cells, not 7^d.
    A = dggread(path; ancestors=[PENTCOL, HEXCOLS[1]])[:elevation]
    @test length(A) == PENT + CAPACITY
    @test Array(A) == DD.data(cube)
    @test collect(DD.lookup(A, Cells)) == collect(DD.lookup(cube, Cells))
end

@testset "the round trip" begin
    cube = democube(HEXCOLS)
    path = dest("trip.zarr")
    dggwrite(path, cube; layout=:subzones, ancestor_level=ANCESTOR)

    # The whole store: a cell axis over the COMPLETE level, because a column
    # nobody wrote is not absent, it is fill.
    stack = dggread(path)
    @test keys(stack) == (:elevation,)
    A = stack[:elevation]
    @test length(A) == ncells(GRID)
    @test DGG.level(DD.lookup(A, Cells)) == LEVEL
    @test parent(A) isa DiskArrays.AbstractDiskArray
    @test DD.metadata(stack)["layout"] == LAYOUT

    values = DD.data(cube)
    for (k, col) in pairs(HEXCOLS)
        window = columnpositions(LAYOUT, col)
        @test A[window] == values[((k-1)*CAPACITY+1):(k*CAPACITY)]
    end
    # Unwritten columns read as NaN, at zero storage cost.
    @test all(isnan, A[columnpositions(LAYOUT, 200)])
    @test count(!isnan, Array(A)) == length(values)

    # The restricted view is the cube that was written, cell for cell.
    B = dggread(path; ancestors=[columncell(LAYOUT, i) for i in HEXCOLS])[:elevation]
    @test Array(B) == values
    @test collect(DD.lookup(B, Cells)) == collect(DD.lookup(cube, Cells))
    # Columns may be named by index as readily as by cell.
    @test Array(dggread(path; ancestors=HEXCOLS)[:elevation]) == values
    # `lazy = false` materializes the same numbers.
    @test dggread(path; ancestors=HEXCOLS, lazy=false)[:elevation] == B

    @test caught(() -> dggread(path; vars=(:slope,))) isa DGGSFormatError

    # A cube read back out of a store is written again as readily: its axis is a
    # `ChunkedCellLookup` and its data are lazy, which is the path a production
    # rewrite takes and not the one the cube above took.
    plain = dest("plain-src.zarr")
    dggwrite(plain, cube)
    again = dest("again.zarr")
    dggwrite(again, dggread(plain)[:elevation]; layout=:subzones,
        ancestor_level=ANCESTOR)
    @test Array(dggread(again; ancestors=HEXCOLS)[:elevation]) == values
end

@testset "the lazy view is chunked by subtree" begin
    cube = democube([PENTCOL; HEXCOLS])
    path = dest("chunks.zarr")
    dggwrite(path, cube; layout=:subzones, ancestor_level=ANCESTOR)
    A = parent(dggread(path; ancestors=[PENTCOL; HEXCOLS])[:elevation])

    @test DiskArrays.haschunks(A) isa DiskArrays.Chunked
    chunks = DiskArrays.eachchunk(A)
    # Irregular by construction: Zarr's own chunk grid could not hold a 41 next
    # to three 49s, and the view publishes exactly the subtree boundaries.
    @test [length(only(c)) for c in chunks] == [PENT, CAPACITY, CAPACITY, CAPACITY]
    @test DiskArrays.eachchunk(A) === chunks       # built once, memoized

    values = DD.data(cube)
    # A read inside one column, a read spanning two, and the whole thing.
    @test A[3:20] == values[3:20]
    @test A[(PENT-2):(PENT+3)] == values[(PENT-2):(PENT+3)]
    @test A[:] == values
    @test A[end] == values[end]

    # The whole-level view chunks by every column of the level.
    whole = parent(dggread(path)[:elevation])
    @test length(DiskArrays.eachchunk(whole)) == LAYOUT.ncolumns
    @test sum(length(only(c)) for c in DiskArrays.eachchunk(whole)) == ncells(GRID)

    # Writing through the view is refused by name rather than by MethodError:
    # a column is written whole, padding rule included.
    @test caught(() -> (A[1] = 0.0f0)) isa ArgumentError
end

@testset "the incremental path is the one-shot path" begin
    cols = [PENTCOL; HEXCOLS]
    cube = democube(cols)
    values = DD.data(cube)

    oneshot = dest("one.zarr")
    dggwrite(oneshot, cube; layout=:subzones, ancestor_level=ANCESTOR)

    # Column by column, from a store created once and never re-stamped.
    piecemeal = dest("many.zarr")
    store = subzonestore(piecemeal, SYS, LEVEL; ancestor_level=ANCESTOR,
        layers=("elevation" => Float32,))
    @test keys(store) == ["elevation"]
    offset = 0
    for col in reverse(cols)            # out of order, as workers finish
        h = columnlength(LAYOUT, col)
        lo = sum(columnlength(LAYOUT, c) for c in cols if c < col; init=0)
        dggwrite!(store, columncell(LAYOUT, col), values[(lo+1):(lo+h)])
    end

    # Byte for byte in the arrays, and attribute for attribute in the group.
    a = Zarr.zopen(oneshot)["elevation"][:, :]
    b = Zarr.zopen(piecemeal)["elevation"][:, :]
    @test isequal(a, b)
    @test Zarr.zopen(oneshot).attrs["dggs"] == Zarr.zopen(piecemeal).attrs["dggs"]
    @test isequal(Array(dggread(oneshot; ancestors=cols)[:elevation]),
        Array(dggread(piecemeal; ancestors=cols)[:elevation]))

    # A reopened store writes more columns without disturbing the ones there.
    again = subzonestore(piecemeal)
    @test again.layout == LAYOUT
    extra = 300
    dggwrite!(again, columncell(LAYOUT, extra), fill(7.0f0, columnlength(LAYOUT, extra)))
    @test all(==(7.0f0), Zarr.zopen(piecemeal)["elevation"][:, extra])
    @test isequal(Zarr.zopen(piecemeal)["elevation"][:, cols[1]], a[:, cols[1]])
end

@testset "a column write touches one file" begin
    path = dest("files.zarr")
    store = subzonestore(path, SYS, LEVEL; ancestor_level=ANCESTOR,
        layers=("elevation" => Float32,))
    before = Dict(f => mtime(joinpath(path, "elevation", f))
                  for f in readdir(joinpath(path, "elevation")))
    sleep(0.01)
    dggwrite!(store, columncell(LAYOUT, 9), Float32.(1:CAPACITY))

    after = readdir(joinpath(path, "elevation"))
    fresh = setdiff(after, keys(before))
    # Exactly one new chunk file, and nothing shared was rewritten: that is what
    # makes concurrent writes of disjoint columns safe.
    @test length(fresh) == 1
    @test all(f -> mtime(joinpath(path, "elevation", f)) == before[f], keys(before))
    @test !isempty(readdir(path))       # .zgroup/.zattrs/.zmetadata still there
end

@testset "several layers, and the errors" begin
    path = dest("two.zarr")
    store = subzonestore(path, SYS, LEVEL; ancestor_level=ANCESTOR,
        layers=("elevation" => Float32, "slope" => Float64))
    @test keys(store) == ["elevation", "slope"]

    # Two layers means the layer has to be named.
    @test caught(() -> dggwrite!(store, 9, Float32.(1:CAPACITY))) isa ArgumentError
    dggwrite!(store, 9, Float32.(1:CAPACITY); var=:elevation)
    dggwrite!(store, 9, Float64.(1:CAPACITY); var="slope")
    @test caught(() -> dggwrite!(store, 9, Float32.(1:3); var=:elevation)) isa ArgumentError
    @test caught(() -> dggwrite!(store, 9, Float32.(1:CAPACITY); var=:aspect)) isa ArgumentError
    # A pentagon column takes p(d) values, not 7^d.
    @test caught(() -> dggwrite!(store, PENTCOL, Float32.(1:CAPACITY);
        var=:elevation)) isa ArgumentError
    dggwrite!(store, PENTCOL, Float32.(1:PENT); var=:elevation)

    stack = dggread(path; ancestors=[9])
    @test keys(stack) == (:elevation, :slope)
    @test Array(stack[:slope]) == Float64.(1:CAPACITY)

    # A NamedTuple writes both layers of one column.
    dggwrite!(store, 10, (elevation=fill(1.0f0, CAPACITY), slope=fill(2.0, CAPACITY)))
    @test Array(dggread(path; ancestors=[10])[:slope]) == fill(2.0, CAPACITY)

    # The one-shot path takes a stack, and each layer keeps its element type.
    cube = democube(HEXCOLS[1:2])
    values = DD.data(cube)
    stacked = DD.DimStack((elevation=copy(values), slope=Float64.(values)),
        (DD.dims(cube, Cells),))
    both = dest("stack.zarr")
    dggwrite(both, stacked; layout=:subzones, ancestor_level=ANCESTOR)
    read_back = dggread(both; ancestors=HEXCOLS[1:2])
    @test keys(read_back) == (:elevation, :slope)
    @test eltype(read_back[:elevation]) === Float32
    @test Array(read_back[:slope]) == Float64.(values)
end

@testset "what the writer refuses" begin
    cube = democube(HEXCOLS)

    # A cube whose coverage is not whole columns.
    partial = let cv = subzone_cellvector(LAYOUT, [5])
        cells = collect(cv)[2:end]
        DD.DimArray(Float32.(1:length(cells)),
            (Cells(CellLookup(CellVector(SYS, LEVEL, cells))),); name=:elevation)
    end
    err = caught(() -> dggwrite(dest("partial.zarr"), partial;
        layout=:subzones, ancestor_level=ANCESTOR))
    @test err isa DGGSFormatError && err.check === :incomplete_subtree

    # The keywords of the other layout, and no ancestor level at all.
    @test caught(() -> dggwrite(dest("a.zarr"), cube; layout=:subzones,
        ancestor_level=ANCESTOR, chunks=64)) isa ArgumentError
    @test caught(() -> dggwrite(dest("b.zarr"), cube; layout=:subzones)) isa ArgumentError
    @test caught(() -> dggwrite(dest("c.zarr"), cube; layout=:tiles)) isa ArgumentError
    @test caught(() -> dggwrite("gs://bucket/x.zarr", cube; layout=:subzones,
        ancestor_level=ANCESTOR)) isa ArgumentError

    # A layer with a second dimension: the layout spends both of its own on the
    # cell axis.
    wide = DD.DimArray(zeros(Float32, length(DD.lookup(cube, Cells)), 2),
        (DD.dims(cube, Cells), DD.Dim{:month}(1:2)); name=:elevation)
    @test caught(() -> dggwrite(dest("d.zarr"), wide; layout=:subzones,
        ancestor_level=ANCESTOR)) isa ArgumentError

    # `ancestors` names columns of a store that has them.
    plain = dest("plain.zarr")
    dggwrite(plain, cube)
    @test caught(() -> dggread(plain; ancestors=[1])) isa ArgumentError
    # And a subzone store is not read through a `description`.
    sub = dest("sub.zarr")
    dggwrite(sub, cube; layout=:subzones, ancestor_level=ANCESTOR)
    @test caught(() -> dggread(sub;
        description=DGG.StoreDescription(gridname="igeo7"))) isa ArgumentError
end

@testset "a regrid, written and read back" begin
    # A coarse global raster with abutting cell edges, so a conservative regrid
    # off it is defined everywhere the target reaches.
    step = 15.0
    lon = (-180+step/2):step:180
    lat = (-90+step/2):step:90
    axis(D, centres) = D(DD.Sampled(collect(centres); span=DD.Regular(step),
        sampling=DD.Intervals(DD.Center()), order=DD.ForwardOrdered()))
    raster = DD.DimArray([10 + sind(x) + cosd(2y) for x in lon, y in lat],
        (axis(DD.X, lon), axis(DD.Y, lat)); name=:elevation)

    # The target is ancestor-snapped by construction: complete level-2 subtrees.
    cols = [PENTCOL, 40, 41]
    target = subzone_cellvector(LAYOUT, cols)
    direct = DGG.regrid(raster; to=target)
    @test DD.lookup(direct, Cells) == CellLookup(target)

    path = dest("regrid.zarr")
    dggwrite(path, DD.rebuild(direct; name=:elevation); layout=:subzones,
        ancestor_level=ANCESTOR)
    back = dggread(path; ancestors=cols)[:elevation]

    @test isequal(Array(back), Array(DD.data(direct)))
    @test collect(DD.lookup(back, Cells)) == collect(DD.lookup(direct, Cells))
end

end # if HAS_ZARR

end # module DGGSubzoneStoreTests
