# The pyramid store on disk: what the bytes are, which of them a half-covered
# write touches, and that the lazy cell-axis view hands back exactly what went
# in at every level.
#
# The layout's arithmetic is `pyramid.jl`'s and is not repeated. What is under
# test here is the store: the array shapes and chunk grid read back with Zarr
# directly rather than through `dggread`, the claim that a chunk nothing landed
# in is never stored, the claim that the coarse level says which of the fine
# level's chunks exist, and the round trip through the DiskArrays view.
#
# The suite self-skips when Zarr.jl is absent, as the other extension suites do.

module DGGPyramidStoreTests

using Test
import DiscreteGlobalGrids as DGG

const HAS_ZARR = try
    @eval using Zarr
    true
catch err
    @warn "Zarr.jl is not loadable: the pyramid-store suite is skipped." exception = err
    false
end

if HAS_ZARR

import DimensionalData as DD
import DiskArrays
using Statistics: mean
using DiscreteGlobalGrids: IGeo7System, Cells, CellLookup, CellVector,
    DGGSFormatError, PyramidLayout, StorePyramid, ancestor, cellindex, cellvalues,
    children, chunkcount, chunkindices, chunklength, chunkrootlevel, dggread,
    dggwrite, descendant_range, globalindex, holdsdata, levelgrid, levels,
    levelslots, ncells, rootcells, slotcell, slotcount, slotindex, system

const SYS = IGeo7System()
const LEVEL = 4
const EXPONENT = 2
const LAYOUT = PyramidLayout(SYS, LEVEL; chunkexponent=EXPONENT)
const GRID = levelgrid(SYS, LEVEL)
const N = ncells(SYS, LEVEL)

dest(name) = joinpath(mktempdir(), name)

"A cube over `cells`, values ascending from 1."
democube(cv; name=:elevation, T=Float64) =
    DD.DimArray(T.(1:length(cv)), (Cells(CellLookup(cv)),); name=name)

wholeearth(; kw...) = democube(CellVector(GRID); kw...)

"Every level-`l` cell under `roots`, as one cell vector."
subtree(roots, l) = CellVector(SYS, l,
    sort!(reduce(vcat, [collect(DGG.descendants(SYS, r, l)) for r in roots])))

@testset "what a whole-earth store is" begin
    d = dest("whole")
    @test dggwrite(d, wholeearth(); layout=:pyramid, chunks=49) == d

    g = Zarr.zopen(d)
    @test sort!(collect(keys(g.groups))) == ["elevation"]
    v = g.groups["elevation"]
    @test sort!(collect(keys(v.arrays))) ==
          ["level0", "level1", "level2", "level3", "level4"]

    # Each level is its SLOT space, not its cell count, and is chunked at a
    # power of seven that divides it.
    for J in 0:LEVEL
        z = v.arrays["level$J"]
        @test size(z) == (12 * 7^J,)
        @test only(z.metadata.chunks) == 7^min(EXPONENT, J)
        @test isnan(z.metadata.fill_value)
    end

    # The absent addresses really are fill on disk, and the cells really are
    # where the slot arithmetic says.
    z4 = v.arrays["level4"]
    for i in 1:N
        @test z4[slotindex(SYS, cellindex(GRID, i))] == Float64(i)
    end
    absent = [i for i in 1:slotcount(SYS, LEVEL) if slotcell(SYS, LEVEL, i) === nothing]
    @test length(absent) == 2 * (7^LEVEL - 1)
    @test all(isnan, z4[absent])

    # No `zarr_conventions` and no cell coordinate: this layout carries no array
    # of ids for a generic reader to fingerprint.
    @test !haskey(g.attrs, "zarr_conventions")
    @test !any(n -> occursin("cell", n), keys(v.arrays))
    @test g.attrs["dggs"]["pyramid_layout"]["variables"] == ["elevation"]
    @test g.attrs["dggs"]["pyramid_layout"]["aggregate"] == "mean"
end

@testset "a chunk nothing landed in is never stored" begin
    d = dest("whole")
    dggwrite(d, wholeearth(); layout=:pyramid, chunks=49)
    for J in 0:LEVEL
        files = filter(!startswith('.'), readdir(joinpath(d, "elevation", "level$J")))
        # One file per chunk that holds a cell, and none for the chunks under an
        # absent address: those are entirely on a deleted branch.
        @test length(files) == ncells(SYS, chunkrootlevel(LAYOUT, J))
        @test length(files) < chunkcount(LAYOUT, J) || J <= EXPONENT
    end

    # Half an earth costs half an earth. Six of the twelve base subtrees, so
    # roughly half the chunks, and none of the other half.
    half = dest("half")
    roots = collect(rootcells(SYS))[1:6]
    dggwrite(half, democube(subtree(roots, LEVEL)); layout=:pyramid, chunks=49)
    whole = length(filter(!startswith('.'),
        readdir(joinpath(d, "elevation", "level4"))))
    part = length(filter(!startswith('.'),
        readdir(joinpath(half, "elevation", "level4"))))
    @test part == whole ÷ 2
    # And the levels above it are half-covered in the same proportion, which is
    # what makes them usable as the index.
    @test length(filter(!startswith('.'),
        readdir(joinpath(half, "elevation", "level0")))) == 6
end

@testset "the round trip" begin
    d = dest("whole")
    dggwrite(d, wholeearth(); layout=:pyramid, chunks=49)
    pyr = dggread(d)

    @test pyr isa StorePyramid
    @test levels(pyr) == 0:LEVEL
    @test keys(pyr) == (:elevation,)
    @test system(pyr) == SYS
    @test pyr.layout == LAYOUT
    @test occursin("elevation", sprint(show, pyr))

    # A level is a cube over the COMPLETE level in canonical order: the slot
    # space stopped at the store.
    for J in 0:LEVEL
        A = pyr[J][:elevation]
        @test length(A) == ncells(SYS, J)
        @test DD.lookup(A, Cells) == CellLookup(CellVector(levelgrid(SYS, J)))
    end

    # The finest level is exactly what went in, read whole, in pieces, and one
    # cell at a time.
    A = pyr[LEVEL][:elevation]
    @test collect(A) == Float64.(1:N)
    @test A[1000:2000] == Float64.(1000:2000)
    @test only(A[7]) == 7.0
    # Including across a pentagon's subtree, where the slots and the cells part
    # company.
    pentagon = foldl((c, _) -> first(children(SYS, c)), 1:EXPONENT; init=first(rootcells(SYS)))
    r = descendant_range(SYS, pentagon, LEVEL)
    @test A[r] == Float64.(r)

    # Every coarser level is the mean of the children that carry a value —
    # checked against the level below as the store holds it, not against a
    # formula.
    for J in LEVEL:-1:1
        fine = collect(pyr[J][:elevation])
        coarse = collect(pyr[J-1][:elevation])
        pgrid = levelgrid(SYS, J - 1)
        for i in eachindex(coarse)
            kids = fine[descendant_range(SYS, cellindex(pgrid, i), J)]
            @test coarse[i] ≈ mean(kids)
        end
    end

    @test DiskArrays.haschunks(parent(A)) isa DiskArrays.Chunked
    chunks = DiskArrays.eachchunk(parent(A))
    # One chunk per cell of the chunk root level, each the subtree it stands for.
    @test length(chunks) == ncells(SYS, chunkrootlevel(LAYOUT, LEVEL))
    @test sum(c -> length(only(c)), chunks) == N
    @test reduce(vcat, [collect(only(c)) for c in chunks]) == collect(1:N)

    @test_throws ArgumentError parent(A)[1:2] = [0.0, 0.0]
end

@testset "reading a level that was only half written" begin
    d = dest("half")
    roots = collect(rootcells(SYS))[1:6]
    cv = subtree(roots, LEVEL)
    dggwrite(d, democube(cv); layout=:pyramid, chunks=49)
    pyr = dggread(d)[:elevation]

    A = pyr[LEVEL][:elevation]
    @test length(A) == N
    written = Set(globalindex(GRID, c) for c in cv)
    values = collect(A)
    @test all(i -> isnan(values[i]), setdiff(1:N, written))
    @test count(!isnan, values) == length(cv)

    # `holdsdata` is exact here, not an estimate: the stored coarse value says
    # whether anything under the cell was written.
    for r in collect(rootcells(SYS))[1:6]
        @test holdsdata(pyr, r)
    end
    for r in collect(rootcells(SYS))[7:12]
        @test !holdsdata(pyr, r)
    end
    # And the same all the way down.
    deep = cellindex(GRID, first(descendant_range(SYS, roots[1], LEVEL)))
    @test holdsdata(pyr, deep)
    outside = cellindex(GRID, first(descendant_range(SYS, collect(rootcells(SYS))[8], LEVEL)))
    @test !holdsdata(pyr, outside)

    # A parent is non-fill whenever ANY cell under it is, which is the invariant
    # the chunk index rests on.
    for J in 0:LEVEL
        coarse = collect(pyr[J][:elevation])
        grid = levelgrid(SYS, J)
        for i in eachindex(coarse)
            under = any(c -> globalindex(GRID, c) in written,
                DGG.descendants(SYS, cellindex(grid, i), LEVEL))
            @test isnan(coarse[i]) != under
        end
    end
end

@testset "cellvalues" begin
    d = dest("half")
    roots = collect(rootcells(SYS))[1:6]
    cv = subtree(roots, LEVEL)
    dggwrite(d, democube(cv); layout=:pyramid, chunks=49)
    pyr = dggread(d)[:elevation]

    # Scattered cells, in no order, some of them over nothing.
    picked = [cellindex(GRID, i) for i in [17, 3, 20000, 9001, 2, 24012]]
    got = cellvalues(pyr, LEVEL, picked)
    @test length(got) == length(picked)
    expected = collect(pyr[LEVEL][:elevation])
    for (k, c) in pairs(picked)
        i = globalindex(GRID, c)
        @test isnan(expected[i]) ? ismissing(got[k]) : got[k] == expected[i]
    end
    @test any(ismissing, got)
    @test cellvalues(pyr, LEVEL, typeof(first(picked))[]) == Union{Float64,Missing}[]

    # A coarse level answers for its own cells, which is what a zoomed-out frame
    # asks for.
    coarse = [cellindex(levelgrid(SYS, 2), i) for i in 1:10]
    @test cellvalues(pyr, 2, coarse) ==
          [isnan(v) ? missing : v for v in collect(pyr[2][:elevation])[1:10]]
    @test_throws ArgumentError cellvalues(pyr, 2, [cellindex(GRID, 1)])
    @test_throws ArgumentError cellvalues(pyr, LEVEL + 1, [cellindex(GRID, 1)])
end

@testset "several variables" begin
    d = dest("two")
    cv = CellVector(GRID)
    stack = DD.DimStack((elevation=Float64.(1:N), slope=Float64.(N:-1:1)),
        (Cells(CellLookup(cv)),))
    dggwrite(d, stack; layout=:pyramid, chunks=49)

    pyr = dggread(d)
    @test keys(pyr) == (:elevation, :slope)
    @test collect(pyr[LEVEL][:slope]) == Float64.(N:-1:1)
    # `holdsdata` and `cellvalues` answer for one variable, and say so.
    @test_throws ArgumentError holdsdata(pyr, first(rootcells(SYS)))
    @test_throws ArgumentError cellvalues(pyr, LEVEL, [cellindex(GRID, 1)])
    one = pyr[:slope]
    @test keys(one) == (:slope,)
    @test holdsdata(one, first(rootcells(SYS)))
    @test cellvalues(one, LEVEL, [cellindex(GRID, 1)]) == [Float64(N)]
    @test_throws ArgumentError pyr[:missing_variable]

    # A single variable can be read by name, and selected on the way in.
    @test keys(dggread(d; vars=(:slope,))) == (:slope,)
end

@testset "the chunk cache" begin
    d = dest("whole")
    dggwrite(d, wholeearth(); layout=:pyramid, chunks=49)
    plain = dggread(d)[:elevation]
    cached = dggread(d; cache=true)[:elevation]
    budget = dggread(d; cache=8)[:elevation]

    # The cache changes what a read costs, never what it answers.
    for p in (cached, budget)
        @test collect(p[LEVEL][:elevation]) == collect(plain[LEVEL][:elevation])
        @test cellvalues(p, 2, [cellindex(levelgrid(SYS, 2), 3)]) ==
              cellvalues(plain, 2, [cellindex(levelgrid(SYS, 2), 3)])
    end
    @test parent(cached[LEVEL][:elevation]) isa DiskArrays.CachedDiskArray
    @test !(parent(plain[LEVEL][:elevation]) isa DiskArrays.CachedDiskArray)
    # Reading the same run twice is the access pattern the cache is for.
    @test cached[LEVEL][:elevation][1:100] == plain[LEVEL][:elevation][1:100]
    @test cached[LEVEL][:elevation][1:100] == plain[LEVEL][:elevation][1:100]

    @test_throws ArgumentError dggread(d; cache=0)
    @test_throws ArgumentError dggread(d; cache=:yes)
    # And it is a pyramid keyword: a flat store says so rather than ignoring it.
    flat = dest("flat")
    dggwrite(flat, wholeearth())
    @test_throws ArgumentError dggread(flat; cache=true)
end

@testset "materialized" begin
    d = dest("whole")
    dggwrite(d, wholeearth(); layout=:pyramid, chunks=49)
    lazy = dggread(d)
    eager = dggread(d; lazy=false)
    @test eager isa StorePyramid
    @test parent(eager[LEVEL][:elevation]) isa Array
    @test collect(eager[LEVEL][:elevation]) == collect(lazy[LEVEL][:elevation])
    @test holdsdata(eager[:elevation], first(rootcells(SYS)))
end

@testset "what the writer refuses" begin
    cv = CellVector(GRID)

    # An element type with no NaN needs a fill value named, and a `mean` that
    # does not land back in it needs an aggregate named.
    ints = democube(cv; T=Int32)
    @test_throws ArgumentError dggwrite(dest("i"), ints; layout=:pyramid, chunks=49)
    @test_throws ArgumentError dggwrite(dest("i"), ints; layout=:pyramid, chunks=49,
        fill_value=typemin(Int32))
    d = dest("i")
    dggwrite(d, ints; layout=:pyramid, chunks=49, fill_value=typemin(Int32),
        aggregate=vs -> round(Int32, mean(vs)))
    @test eltype(dggread(d)[0][:elevation]) == Int32
    # An anonymous reduction has no name worth recording.
    @test !haskey(Zarr.zopen(d).attrs["dggs"]["pyramid_layout"], "aggregate")

    # A chunk that is not a power of the aperture is not a subtree.
    @test_throws ArgumentError dggwrite(dest("c"), wholeearth(); layout=:pyramid,
        chunks=1000)
    @test_throws ArgumentError dggwrite(dest("c"), wholeearth(); layout=:pyramid,
        chunks=7^5)

    # The one-dimensional layout's keywords mean nothing here.
    @test_throws ArgumentError dggwrite(dest("k"), wholeearth(); layout=:pyramid,
        encoding=:dense)
    @test_throws ArgumentError dggwrite(dest("k"), wholeearth(); layout=:pyramid,
        chunk_target=1000)

    # An unnamed array has nothing to call its subgroup.
    unnamed = DD.DimArray(Float64.(1:N), (Cells(CellLookup(cv)),))
    @test_throws ArgumentError dggwrite(dest("u"), unnamed; layout=:pyramid)
    # And a layer named like a level array would collide with one.
    @test_throws ArgumentError dggwrite(dest("u"),
        democube(cv; name=:level3); layout=:pyramid)

    # Writing over a store that already holds the variable, and the way past it.
    # A path names the whole store, so `overwrite` clears the whole store; a
    # group names one variable of it.
    d = dest("again")
    dggwrite(d, wholeearth(); layout=:pyramid, chunks=49)
    @test_throws Exception dggwrite(d, wholeearth(); layout=:pyramid, chunks=49)
    dggwrite(d, wholeearth(); layout=:pyramid, chunks=49, overwrite=true)
    @test collect(dggread(d)[LEVEL][:elevation]) == Float64.(1:N)

    g = Zarr.zopen(d, "w")
    @test_throws DGGSFormatError dggwrite(g, wholeearth(); layout=:pyramid, chunks=49)
    dggwrite(g, wholeearth(); layout=:pyramid, chunks=49, overwrite=true)
    @test collect(dggread(d)[LEVEL][:elevation]) == Float64.(1:N)

    @test_throws ArgumentError dggwrite("s3://bucket/store", wholeearth();
        layout=:pyramid)
    @test_throws ArgumentError dggwrite(dest("x"), wholeearth(); layout=:nonsense)
end

@testset "what the reader refuses" begin
    d = dest("whole")
    dggwrite(d, wholeearth(); layout=:pyramid, chunks=49)
    @test_throws ArgumentError dggread(d; ancestors=[1])

    # A level removed from a variable is a store that cannot index itself.
    broken = dest("broken")
    dggwrite(broken, wholeearth(); layout=:pyramid, chunks=49)
    rm(joinpath(broken, "elevation", "level2"); recursive=true)
    rm(joinpath(broken, ".zmetadata"))
    @test_throws DGGSFormatError dggread(broken)
end

end # if HAS_ZARR

end # module DGGPyramidStoreTests
