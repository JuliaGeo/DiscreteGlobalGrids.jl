# The pyramid layout's arithmetic and vocabulary, with no store in sight:
# `src/io/pyramid.jl` sees cells, slots, values and attribute dictionaries and
# nothing else. What a real store does with all of it is `pyramid_store.jl`.
#
# The claims here are the ones a wrong offset would break silently: that the
# slot space is the level's ids with the absent addresses left in, that a chunk
# is exactly one subtree at a coarser level, that the coarse level therefore
# indexes the fine level's chunks position for position, and that a reduction up
# the pyramid leaves a parent absent exactly when everything under it is.

module DGGPyramidTests

using Test
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: IGeo7System, A5System, Z7Cell, CellLookup, CellVector,
    DGGSFormatError, PyramidLayout, PyramidRun, ancestor, cellindex, cellslot,
    children, chunkcount, chunkexponent, chunkindices, chunklength, chunkroot,
    chunkrootlevel, chunkslots, descendant_range, globalindex, ispyramidstore,
    levelfromname, levelgrid, levelname, levelslots, levels, ncells,
    pyramid_attrs, pyramid_layout, pyramid_reduce, pyramid_runs, rootcells,
    slotcell, slotchunk, slotcount, slotindex

using Statistics: mean

const SYS = IGeo7System()

# The closed forms the layout is measured against, spelled out here rather than
# imported: a test that read the implementation's own constants would not notice
# either of them changing.
slots(l) = 12 * 7^l
cells(l) = 10 * 7^l + 2
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

@testset "the slot space" begin
    for L in 0:4
        @test slotcount(SYS, L) == slots(L)
        @test ncells(SYS, L) == cells(L)
        # Every absent address is on a pentagon's deleted branch, and there are
        # `2·(7^L - 1)` of them however the deletion is spelled.
        @test slotcount(SYS, L) - ncells(SYS, L) == 2 * (7^L - 1)
    end

    # `slotindex` is the dense order with the gaps left in: strictly ascending
    # over the level's own order, and `slotcell` inverts it.
    for L in 0:3
        grid = levelgrid(SYS, L)
        taken = Int[]
        previous = 0
        for i in 1:ncells(SYS, L)
            c = cellindex(grid, i)
            s = slotindex(SYS, c)
            @test s > previous
            @test slotcell(SYS, L, s) == c
            previous = s
            push!(taken, s)
        end
        # And the slots nothing took are exactly the ones naming no cell.
        absent = setdiff(1:slotcount(SYS, L), taken)
        @test all(i -> slotcell(SYS, L, i) === nothing, absent)
        @test length(absent) == slotcount(SYS, L) - ncells(SYS, L)
    end

    @test_throws Exception slotcell(SYS, 2, 0)
    @test_throws Exception slotcell(SYS, 2, slots(2) + 1)
end

@testset "the layout" begin
    l = PyramidLayout(SYS, 4; chunkexponent=2)
    @test levels(l) == 0:4
    @test l.aperture == 7
    @test levelname(l, 0) == "level0"
    @test levelname(l, 13) == "level13"
    @test levelfromname("level13") == 13
    @test levelfromname("level0") == 0
    # A second spelling of a level is no level at all: two arrays claiming one
    # level is a store no reader can resolve.
    @test levelfromname("level07") === nothing
    @test levelfromname("levels") === nothing
    @test levelfromname("cell_ids") === nothing

    for J in levels(l)
        @test levelslots(l, J) == slots(J)
        # A level too shallow for the nominal exponent is chunked at its own
        # depth rather than at one that would not divide it.
        @test chunkexponent(l, J) == min(2, J)
        @test chunklength(l, J) == 7^min(2, J)
        @test chunkrootlevel(l, J) == J - min(2, J)
        @test chunkcount(l, J) == slots(J - min(2, J))
        # The chunk grid divides the slot space exactly. This is the property
        # the whole layout is for.
        @test chunkcount(l, J) * chunklength(l, J) == levelslots(l, J)
    end

    @test_throws ArgumentError PyramidLayout(SYS, -1)
    @test_throws ArgumentError PyramidLayout(SYS, 3; chunkexponent=-1)
    # A system whose subtrees are scattered has no chunk that is a subtree.
    @test_throws ArgumentError PyramidLayout(A5System(), 3)
    @test_throws ArgumentError levelslots(PyramidLayout(SYS, 3), 4)

    @test occursin("igeo7", sprint(show, l))
    @test occursin("7^2", sprint(show, l))
    @test PyramidLayout(SYS, 4; chunkexponent=2) == l
    @test PyramidLayout(SYS, 4; chunkexponent=3) != l
end

@testset "a chunk is a subtree" begin
    l = PyramidLayout(SYS, 4; chunkexponent=2)
    for J in levels(l)
        w = chunklength(l, J)
        grid = levelgrid(SYS, J)
        covered = 0
        for c in 1:chunkcount(l, J)
            root = chunkroot(l, J, c)
            r = chunkindices(l, J, c)
            off = chunkslots(l, J, c)
            @test length(off) == length(r)
            if root === nothing
                # A chunk under an absent address holds no cell at all, and is
                # therefore a chunk that is never written.
                @test isempty(r)
                continue
            end
            @test r == descendant_range(SYS, root, J)
            covered += length(r)
            # The offsets land inside the chunk, ascend, and are the slots of
            # the cells they stand for.
            @test all(1 .<= off .<= w)
            @test issorted(off) && allunique(off)
            for (t, p) in enumerate(r)
                cell = cellindex(grid, p)
                @test slotindex(SYS, cell) == (c - 1) * w + off[t]
                @test slotchunk(l, J, slotindex(SYS, cell)) == (c, off[t])
                @test cellslot(l, cell) == slotindex(SYS, cell)
            end
            # A complete subtree is the identity, which is why only the twelve
            # pentagon-rooted chunks cost anything to lay out.
            if length(r) == w
                @test off === 1:w
            else
                @test length(r) == pent(chunkexponent(l, J))
            end
        end
        # Every cell of the level is in exactly one chunk.
        @test covered == ncells(SYS, J)
    end
end

@testset "the coarse level indexes the fine level's chunks" begin
    # The claim the read side rests on: chunk `c` of level `J` is slot `c` of
    # level `chunkrootlevel`, position for position and absence for absence.
    l = PyramidLayout(SYS, 5; chunkexponent=2)
    for J in levels(l)
        rl = chunkrootlevel(l, J)
        @test chunkcount(l, J) == levelslots(l, rl)
        for c in 1:chunkcount(l, J)
            @test chunkroot(l, J, c) === slotcell(SYS, rl, c)
        end
    end
end

@testset "a cell axis maps onto chunks" begin
    l = PyramidLayout(SYS, 4; chunkexponent=2)
    J = 4
    grid = levelgrid(SYS, J)

    # The whole level: every chunk that holds a cell, each covered whole, and
    # the runs concatenate to the axis in order.
    runs = pyramid_runs(l, J, CellLookup(CellVector(grid)))
    @test sum(length(r.axis) for r in runs) == ncells(SYS, J)
    @test reduce(vcat, [collect(r.axis) for r in runs]) == collect(1:ncells(SYS, J))
    @test issorted([r.chunk for r in runs])
    @test allunique([r.chunk for r in runs])
    for r in runs
        @test length(chunkindices(l, J, r.chunk)) == length(r.rows)
        @test first(r.rows) == 1
    end

    # A subset, and one that deliberately covers a chunk only in part: a chunk
    # is written from a prefilled buffer, so a partial chunk is ordinary here
    # where the ancestor-subzone layout refuses it.
    coarse = cellindex(levelgrid(SYS, 2), 40)
    r = descendant_range(SYS, coarse, J)
    picked = [cellindex(grid, p) for p in (first(r)+3):(first(r)+9)]
    runs = pyramid_runs(l, J, picked)
    @test length(runs) == 1
    @test only(runs).rows == 4:10
    @test only(runs).axis == 1:7
    @test only(runs).chunk == slotindex(SYS, coarse)

    # Two disjoint pieces of one chunk come back as two runs of that chunk.
    picked = [cellindex(grid, p) for p in [first(r), first(r) + 5, first(r) + 6]]
    runs = pyramid_runs(l, J, picked)
    @test [x.chunk for x in runs] == [slotindex(SYS, coarse), slotindex(SYS, coarse)]
    @test [x.rows for x in runs] == [1:1, 6:7]
    @test [x.axis for x in runs] == [1:1, 2:3]

    # A pentagon-rooted chunk is short, and its runs still cover it whole.
    pentagon = pentagons(SYS, 2)[1]
    runs = pyramid_runs(l, J, CellLookup(CellVector(SYS, J, DGG.descendants(SYS, pentagon, J))))
    @test sum(length(x.rows) for x in runs) == pent(2)

    @test_throws ArgumentError pyramid_runs(l, J, reverse(picked))
    @test_throws ArgumentError pyramid_runs(l, 3, CellLookup(CellVector(grid)))
end

@testset "reducing up the pyramid" begin
    l = PyramidLayout(SYS, 3; chunkexponent=1)
    J = 3
    grid = levelgrid(SYS, J)
    n = ncells(SYS, J)

    # A complete level of ones reduces to a complete level of ones, all the way
    # to the twelve root cells.
    idx = collect(1:n)
    val = ones(Float64, n)
    for j in J:-1:1
        idx, val = pyramid_reduce(l, j, idx, val, mean, NaN)
        @test idx == collect(1:ncells(SYS, j - 1))
        @test all(==(1.0), val)
    end
    @test length(val) == 12

    # A parent is the mean of the children that carry a value, and a child
    # holding the fill value is absent rather than NaN-poisoning its parent.
    parent = cellindex(levelgrid(SYS, J - 1), 5)
    r = descendant_range(SYS, parent, J)
    val = fill(NaN, n)
    val[first(r)] = 2.0
    val[first(r)+1] = 4.0
    pidx, pval = pyramid_reduce(l, J, collect(1:n), val, mean, NaN)
    @test pidx == [5]
    @test pval == [3.0]

    # A parent with nothing under it is left out entirely: that is the
    # invariant the chunk index rests on.
    pidx, pval = pyramid_reduce(l, J, collect(1:n), fill(NaN, n), mean, NaN)
    @test isempty(pidx) && isempty(pval)

    # Pentagons reduce over six children, hexagons over seven; the absent
    # addresses simply are not there to reduce.
    counts = Int[]
    pidx, pval = pyramid_reduce(l, J, collect(1:n), ones(Float64, n),
        vs -> (push!(counts, length(vs)); mean(vs)), NaN)
    @test count(==(6), counts) == 12
    @test count(==(7), counts) == ncells(SYS, J - 1) - 12

    # An aggregate that lands on the fill value would mark a populated subtree
    # as empty, and is refused rather than written.
    @test_throws ArgumentError pyramid_reduce(l, J, collect(1:n), ones(Float64, n),
        _ -> NaN, NaN)

    # An aggregate that cannot land back in the element type is refused with the
    # element type named, not with an `InexactError`. A mean over a mix of ones
    # and twos is never a whole number over six or seven children.
    mixed = Int32[i % 2 + 1 for i in 1:n]
    @test_throws ArgumentError pyramid_reduce(l, J, collect(1:n), mixed, mean,
        typemin(Int32))
    # One that can is fine, and every level keeps the element type it started in.
    pidx, pval = pyramid_reduce(l, J, collect(1:n), mixed,
        vs -> round(Int32, mean(vs)), typemin(Int32))
    @test eltype(pval) == Int32
    @test all(v -> v in Int32[1, 2], pval)

    @test_throws ArgumentError pyramid_reduce(l, 0, [1], [1.0], mean, NaN)
    @test_throws ArgumentError pyramid_reduce(l, J, [1, 2], [1.0], mean, NaN)
end

@testset "the attributes" begin
    l = PyramidLayout(SYS, 4; chunkexponent=3)
    attrs = pyramid_attrs(l; variables=["elevation"], aggregate="mean")
    @test ispyramidstore(attrs)
    block = attrs["dggs"]["pyramid_layout"]
    @test block["layout"] == "implicit_pyramid"
    @test block["levels"] == collect(0:4)
    @test block["aperture"] == 7
    @test block["chunk_exponent"] == 3
    @test block["variables"] == ["elevation"]
    @test block["aggregate"] == "mean"
    @test attrs["dggs"]["refinement_level"] == 4
    @test attrs["dggs"]["indexing_scheme"] == "z7int"
    # No convention declaration and no coordinate: this layout carries no array
    # of ids for a generic reader to fingerprint.
    @test !haskey(attrs, "zarr_conventions")

    @test pyramid_layout(attrs) == l

    @test !ispyramidstore(Dict{String,Any}())
    @test !ispyramidstore(Dict{String,Any}("dggs" => Dict{String,Any}()))
    @test_throws DGGSFormatError pyramid_layout(Dict{String,Any}())

    bad = deepcopy(attrs)
    bad["dggs"]["pyramid_layout"]["version"] = 99
    @test_throws DGGSFormatError pyramid_layout(bad)

    bad = deepcopy(attrs)
    bad["dggs"]["pyramid_layout"]["aperture"] = 4
    @test_throws DGGSFormatError pyramid_layout(bad)

    bad = deepcopy(attrs)
    bad["dggs"]["pyramid_layout"]["levels"] = [2, 3, 4]
    @test_throws DGGSFormatError pyramid_layout(bad)

    bad = deepcopy(attrs)
    bad["dggs"]["pyramid_layout"]["slot_order"] = "hilbert"
    @test_throws DGGSFormatError pyramid_layout(bad)

    bad = deepcopy(attrs)
    bad["dggs"]["name"] = "not_a_grid"
    @test_throws DGGSFormatError pyramid_layout(bad)
end

end # module DGGPyramidTests
