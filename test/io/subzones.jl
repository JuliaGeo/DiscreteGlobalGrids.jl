# The ancestor-subzone layout's arithmetic and vocabulary, with no store in
# sight: `src/io/subzones.jl` sees cells, indices and attribute dictionaries
# and nothing else, so this suite needs neither Zarr nor a temporary directory.
# What a real store does with all of it is `subzone_store.jl`.
#
# The claims here are the ones a wrong offset would break silently: that a
# column is exactly one subtree, that the twelve pentagon columns are short by
# the amount the aperture-7 count says, that the cell axis a set of columns
# spells is the inverse of the runs a cell axis maps onto, and that a store's
# attributes read back as the layout that wrote them.

module DGGSubzoneTests

using Test
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: IGeo7System, HEALPixSystem, A5System, Z7Cell,
    CellLookup, CellVector, DGGSFormatError, SubzoneLayout, cellindex, children,
    columncell, columnindex, columnlength, columnindices, descendants,
    descendant_range, issubzonestore, levelgrid, ncells, columnrow,
    rootcells, subzone_attrs, subzone_capacity, subzone_cellvector,
    subzone_columns, subzone_depth, subzone_layout, subzone_runs, subzoneindex

const SYS = IGeo7System()
const LEVEL = 4
const ANCESTOR = 2
const LAYOUT = SubzoneLayout(SYS, LEVEL, ANCESTOR)
const GRID = levelgrid(SYS, LEVEL)

# The closed forms the layout is measured against, spelled out here rather than
# imported: a test that read `POW7` and `PENT_COUNT` out of the implementation
# would not notice either of them changing.
pow7(d) = 7^d
pent(d) = (5 * 7^d + 1) ÷ 6

"The level-`l` pentagon of each root cell: the digit-0 chain, all the way down."
function pentagons(sys, l)
    out = eltype(rootcells(sys))[]
    for r in rootcells(sys)
        c = r
        for _ in 1:l
            c = first(children(sys, c))
        end
        push!(out, c)
    end
    return out
end

@testset "the layout" begin
    @test subzone_depth(LAYOUT) == 2
    @test LAYOUT.capacity == pow7(2)
    @test LAYOUT.ncolumns == ncells(levelgrid(SYS, ANCESTOR))
    @test LAYOUT.level == LEVEL
    @test LAYOUT.ancestor_level == ANCESTOR
    @test LAYOUT.gridname == "igeo7"

    # The measured capacity is the closed form, and a capacity handed in is
    # taken rather than measured again.
    @test subzone_capacity(SYS, ANCESTOR, LEVEL) == pow7(2)
    @test SubzoneLayout(SYS, LEVEL, ANCESTOR; capacity=pow7(2)) == LAYOUT

    @test_throws ArgumentError SubzoneLayout(SYS, 4, 5)
    @test_throws ArgumentError SubzoneLayout(SYS, 4, -1)
    # A5's descendants are not one run of their level, so no cell of it is a
    # column: the layout refuses the system rather than the store.
    @test_throws ArgumentError SubzoneLayout(A5System(), 4, 2)
end

@testset "columns are subtrees" begin
    for i in (1, 2, 42, 100, LAYOUT.ncolumns)
        c = columncell(LAYOUT, i)
        @test DGG.level(c) == ANCESTOR
        @test columnindex(LAYOUT, c) == i
        @test columnindices(LAYOUT, i) == descendant_range(SYS, c, LEVEL)
        @test columnlength(LAYOUT, i) == length(collect(descendants(SYS, c, LEVEL)))
    end
    # Every column together is the whole level, exactly once.
    @test sum(columnlength(LAYOUT, i) for i in 1:LAYOUT.ncolumns) == ncells(GRID)
    @test_throws BoundsError columncell(LAYOUT, 0)
    @test_throws BoundsError columncell(LAYOUT, LAYOUT.ncolumns + 1)
end

@testset "the twelve short columns are the pentagons" begin
    short = [i for i in 1:LAYOUT.ncolumns if columnlength(LAYOUT, i) != LAYOUT.capacity]
    @test length(short) == 12
    @test all(i -> columnlength(LAYOUT, i) == pent(2), short)
    # Which twelve is derivable from the grid, and this is the derivation: the
    # digit-0 chain from each root cell.
    @test short == sort!([columnindex(LAYOUT, p) for p in pentagons(SYS, ANCESTOR)])
    # p(d) < 7^d by exactly the deleted branch, which is what the padding covers.
    @test pent(2) < pow7(2)
end

@testset "cell <-> (column, row)" begin
    for p in (1, 2, 41, 42, 1000, ncells(GRID))
        c = cellindex(GRID, p)
        col, row = subzoneindex(LAYOUT, c)
        @test (col, row) == columnrow(LAYOUT, p)
        r = columnindices(LAYOUT, col)
        @test 1 <= row <= length(r)
        # The inverse: row `row` of column `col` is this cell and no other.
        @test cellindex(GRID, first(r) + row - 1) == c
    end
    # A cell of another level has no place in the store.
    @test_throws ArgumentError subzoneindex(LAYOUT, cellindex(levelgrid(SYS, 3), 1))
end

@testset "runs and the axis are inverse" begin
    cols = [1, 5, 6, 200]
    cv = subzone_cellvector(LAYOUT, cols)
    @test length(cv) == sum(columnlength(LAYOUT, i) for i in cols)
    runs = subzone_runs(LAYOUT, CellLookup(cv))
    @test [r.column for r in runs] == cols
    @test all(r -> r.rows == 1:columnlength(LAYOUT, r.column), runs)
    @test all(r -> length(r.rows) == length(r.axis), runs)
    # The axis pieces tile 1:length(cv) in order.
    @test reduce(vcat, [collect(r.axis) for r in runs]) == collect(1:length(cv))
    # And the cells really are the ones the runs claim.
    for r in runs
        window = columnindices(LAYOUT, r.column)
        @test collect(cv[r.axis]) == [cellindex(GRID, p) for p in window]
    end

    # The whole store: every column, which is the complete level in one window.
    whole = subzone_cellvector(LAYOUT, nothing)
    @test length(whole) == ncells(GRID)
    @test whole == CellVector(GRID)

    # An explicit id vector is the same answer by a slower road.
    ids = collect(cv)
    @test subzone_runs(LAYOUT, ids) == runs

    @test subzone_columns(LAYOUT, nothing) == collect(1:LAYOUT.ncolumns)
    @test subzone_columns(LAYOUT, [columncell(LAYOUT, 6), columncell(LAYOUT, 1)]) == [1, 6]
    @test_throws ArgumentError subzone_columns(LAYOUT, [0])
end

@testset "a partly covered column is refused" begin
    # One cell short of a complete subtree, at both ends and in the middle.
    full = collect(columnindices(LAYOUT, 5))
    for indices in (full[2:end], full[1:end-1], [full[1:10]; full[12:end]])
        cells = [cellindex(GRID, p) for p in indices]
        err = try
            subzone_runs(LAYOUT, cells)
            nothing
        catch e
            e
        end
        @test err isa DGGSFormatError
        @test err.check === :incomplete_subtree
        # The same cells are describable when nobody asks for whole columns.
        @test !isempty(subzone_runs(LAYOUT, cells; complete=false))
    end
    # A pentagon column is complete at p(d) cells and NOT at 7^d: the padding is
    # not part of the axis.
    pentcol = columnindex(LAYOUT, first(pentagons(SYS, ANCESTOR)))
    cells = [cellindex(GRID, p) for p in columnindices(LAYOUT, pentcol)]
    @test length(cells) == pent(2)
    @test length(subzone_runs(LAYOUT, cells)) == 1
end

@testset "attributes" begin
    attrs = subzone_attrs(LAYOUT; variables=["elevation"],
        coordinate="ancestor_cell_ids")
    @test issubzonestore(attrs)
    @test subzone_layout(attrs) == LAYOUT
    @test DGG.subzone_coordinate(attrs) == "ancestor_cell_ids"

    dggs = attrs["dggs"]
    @test dggs["name"] == "igeo7"
    @test dggs["refinement_level"] == LEVEL
    block = dggs["subzone_layout"]
    @test block["layout"] == DGG.SUBZONE_LAYOUT
    @test block["ancestor_level"] == ANCESTOR
    @test block["subzone_count"] == LAYOUT.capacity
    @test block["ancestor_count"] == LAYOUT.ncolumns
    @test block["subzone_order"] == "ascending_id"
    @test block["padding"] == "trailing_fill"
    @test block["padding_fill_value"] == "NaN"
    @test block["chunk_shape"] == [1, LAYOUT.capacity]
    @test block["variables"] == ["elevation"]
    # No `zarr_conventions` declaration: this is not the one-dimensional layout
    # that convention describes, and saying so would send its reader down a path
    # that cannot open the store.
    @test !haskey(attrs, "zarr_conventions")

    # Silence, and other people's stores, are not subzone stores.
    @test !issubzonestore(Dict{String,Any}())
    @test !issubzonestore(Dict{String,Any}("dggs" => Dict{String,Any}("name" => "igeo7")))
    @test !issubzonestore("dggs")
    @test_throws DGGSFormatError subzone_layout(Dict{String,Any}())

    # The declared counts are checked against the grid, not believed.
    bad = deepcopy(attrs)
    bad["dggs"]["subzone_layout"]["ancestor_count"] = 7
    @test_throws DGGSFormatError subzone_layout(bad)
    bad = deepcopy(attrs)
    bad["dggs"]["subzone_layout"]["version"] = 99
    @test_throws DGGSFormatError subzone_layout(bad)
    bad = deepcopy(attrs)
    bad["dggs"]["subzone_layout"]["subzone_order"] = "morton"
    @test_throws DGGSFormatError subzone_layout(bad)
    bad = deepcopy(attrs)
    bad["dggs"]["name"] = "not-a-grid"
    @test_throws DGGSFormatError subzone_layout(bad)
end

@testset "another system's columns" begin
    # HEALPix nested: every subtree is 4^d and no column is short, so the
    # padding rule is real but unused. The arithmetic is the same arithmetic.
    hpx = SubzoneLayout(HEALPixSystem(), 4, 2)
    @test hpx.capacity == 4^2
    @test hpx.ncolumns == ncells(levelgrid(HEALPixSystem(), 2))
    @test all(i -> columnlength(hpx, i) == 16, 1:hpx.ncolumns)
    @test subzone_layout(subzone_attrs(hpx)) == hpx
end

end # module DGGSubzoneTests
