# Cross-system laws for the DimensionalData cell lookup. The tests compare its
# vector interface and selectors with independent cell-set expansions.

module DimensionalDataTests

using Test
using Statistics
import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GeoInterface as GI
import GeometryOps as GO
import SmallCollections

include(joinpath(@__DIR__, "..", "..", "helpers.jl"))
using .DGGTestHelpers: syslabel, sweepcovers

const CL = DGG.CellLookups


const REGION = GI.Polygon([GI.LinearRing([(6.0, 45.8), (10.5, 45.8), (10.5, 47.8),
    (6.0, 47.8), (6.0, 45.8)])])
const ZURICH = GI.Polygon([GI.LinearRing([(8.3, 47.2), (8.8, 47.2), (8.8, 47.6),
    (8.3, 47.6), (8.3, 47.2)])])
const FARAWAY = (-25.0, -40.0)

const LONLAT = GO.UnitSpherical.GeographicFromUnitSphere()

const SWEEP = [
    (DGG.IGeo7System(), 6, 3),
    (DGG.H3System(), 5, 3),
    (DGG.HEALPixSystem(), 9, 3),
    (DGG.A5System(), 9, 2),
    (DGG.S2System(), 9, 3),
    (DGG.ISEA4RSystem(), 9, 3),
    (DGG.AuthalicSystem(DGG.IGeo7System()), 6, 3),
]

@testset "the sweep covers every registered system" begin
    sweepcovers(SWEEP)
end

function expand(sys, set, l::Int)
    grid = DGG.levelgrid(sys, l)
    DGG.has_sorted_subtrees(sys) &&
        return [DGG.cellindex(grid, p) for r in DGG.level_ranges(set, l) for p in r]
    return sort!(reduce(vcat, [DGG.descendants(sys, c, l) for c in set]))
end

# ---------------------------------------------------------------------------
# The laws, once per system
# ---------------------------------------------------------------------------

@testset "a multi-order cell axis: $(syslabel(sys))" for (sys, leaf, deeper) in SWEEP
    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf)
    grid = DGG.levelgrid(sys, leaf)
    lk = DGG.CellLookup(set)
    ids = expand(sys, set, leaf)

    @testset "what it is" begin
        @test lk isa DD.Lookups.Lookup
        @test DGG.level(lk) == leaf
        @test DGG.system(lk) == sys
        # `parent` is the VALUES, per DimensionalData's contract, and the
        # backing has its own name. The whole of the DD machinery reads the
        # former, so the two must not be confused.
        @test parent(lk) isa AbstractVector{DGG.cellindextype(sys)}
        @test parent(lk) == ids
        @test DGG.cellset(lk) === set
        @test length(lk) == length(ids)
        @test eltype(lk) === DGG.cellindextype(sys)
        # A mixed-level backing standing for a single-level axis is the whole
        # premise; a set that happened to be flat would make the rest vacuous.
        @test length(unique(DGG.level.(collect(set)))) > 1
        @test all(==(leaf), DGG.level.(collect(lk)))
    end

    @testset "the lazy form equals the materialised leaf vector" begin
        @test collect(lk) == ids
        @test [c for c in lk] == ids
        @test first(lk) == first(ids)
        @test last(lk) == last(ids)
        @test lk[end] == last(ids)
        @test lk[begin] == first(ids)
        @test all(lk[k] == ids[k] for k in eachindex(ids))
        @test collect(lk[2:5]) == ids[2:5]
        @test lk[:] === lk
    end

    @testset "index <-> id round trips" begin
        @test all(DGG.localindex(lk, lk[k]) == k for k in eachindex(ids))
        outside = DGG.cellat(grid, FARAWAY...)
        @test outside !== nothing
        @test DGG.localindex(lk, outside) === nothing
        # An id from another level names no index here, rather than the
        # index of a cell with the same raw bits.
        coarser = DGG.ancestor(sys, first(ids), leaf - 1)
        @test DGG.localindex(lk, coarser) === nothing
    end

    @testset "the degenerate cases answer identically" begin
        # A whole level: one window, and the lookup's indices ARE the grid's.
        complete = DGG.CellLookup(grid)
        @test length(complete) == DGG.ncells(grid)
        @test CL.nwindows(CL.windows(complete)) == 1
        @test all(complete[k] == DGG.cellindex(grid, k) for k in (1, 7, DGG.ncells(grid)))

        # PARTIAL ⊂ COMPLETE: the subset lookup and the whole-level lookup
        # agree about every cell of the subset, and the whole-level lookup's
        # indices are the grid's own.
        @test all(DGG.localindex(complete, c) == DGG.globalindex(grid, c) for c in lk)
        @test all(DGG.localindex(complete, c) !== nothing for c in lk)

        # An arbitrary ascending subset, the "partial lookup": same content,
        # same answers, and equal to the multi-order form it was built from.
        partial = DGG.CellLookup(DGG.PartialGrid(sys, leaf, ids))
        @test collect(partial) == ids
        @test partial == lk
        @test all(DGG.localindex(partial, ids[k]) == k for k in eachindex(ids))

        # A rooted subtree, which is the one shape both backings can hold.
        root = DGG.ancestor(sys, first(ids), leaf - 1)
        rooted = DGG.CellLookup(DGG.subtree(sys, root, leaf))
        @test collect(rooted) == DGG.descendants(sys, root, leaf)
        @test DGG.localindex(rooted, first(ids)) !== nothing
    end

    @testset "the lookup read as a grid" begin
        pg = DGG.PartialGrid(lk)
        @test DGG.ncells(pg) == length(lk)
        @test all(DGG.cellindex(pg, k) == lk[k] for k in eachindex(ids))
        @test all(DGG.localindex(pg, lk[k]) == k for k in eachindex(ids))
        @test DGG.level(pg) == leaf
    end
end

# ---------------------------------------------------------------------------
# Selectors, and the cube they run against
# ---------------------------------------------------------------------------

@testset "selectors on a cell axis: $(syslabel(sys))" for (sys, leaf, _) in SWEEP
    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf)
    grid = DGG.levelgrid(sys, leaf)
    lk = DGG.CellLookup(set)
    ids = expand(sys, set, leaf)
    A = DD.DimArray(Float64.(eachindex(ids)), DGG.Cells(lk); name=:demo)

    @testset "the cube is a cube" begin
        @test DD.lookup(A, DGG.Cells) === lk
        @test length(A) == length(ids)
        @test DD.dims(A, DGG.Cells) isa DGG.Cells
    end

    @testset "Where sees the axis, not the backing" begin
        @test length(ids) != length(set)
        @test length(A[DGG.Cells(DD.Where(c -> true))]) == length(ids)
        @test length(A[DGG.Cells(DD.Where(c -> DGG.level(c) == leaf))]) == length(ids)
        # A predicate false at RUNTIME, not one Julia can fold: a literally
        # constant `false` makes the comprehension in DimensionalData's `Where`
        # infer `Vector{Union{}}`, which hits an ambiguity in its own
        # `getindex` — reproducible with a plain `Sampled` lookup and nothing
        # to do with this one.
        @test isempty(A[DGG.Cells(DD.Where(c -> DGG.level(c) == leaf + 1))])
        @test parent(A[DGG.Cells(DD.Where(c -> c == ids[2]))]) == [2.0]
    end

    @testset "At(id) selects that id's index" begin
        for k in (1, length(ids) ÷ 2, length(ids))
            @test A[DGG.Cells(DD.At(ids[k]))] == Float64(k)
        end
        outside = DGG.cellat(grid, FARAWAY...)
        @test_throws DD.Lookups.SelectorError A[DGG.Cells(DD.At(outside))]
    end

    @testset "Contains(lon, lat) selects the cell the point is in" begin
        for k in (1, length(ids) ÷ 3, length(ids) ÷ 2, length(ids))
            lon, lat = LONLAT(DGG.cell_centroid(grid, ids[k]))
            @test A[DGG.Cells(DD.Contains(lon, lat))] == Float64(k)
        end
        @test_throws DD.Lookups.SelectorError A[DGG.Cells(DD.Contains(FARAWAY...))]
    end

    # The polygon selector IS the coverage expansion intersected with the
    # lookup; the hand-rolled right-hand side is the same sentence spelled out.
    @testset "Covering(region) is the coverage expansion, intersected" begin
        for target in (ZURICH, REGION)
            byhand = Int[]
            for c in expand(sys, DGG.query(sys, DGG.MultiOrderCoverage(target); level=leaf), leaf)
                k = DGG.localindex(lk, c)
                k === nothing || push!(byhand, k)
            end
            sort!(byhand)
            sub = A[DGG.Cells(DGG.Covering(target))]
            @test parent(sub) == Float64.(byhand)
            @test length(sub) == length(byhand)
            # Subsetting stays in the compact form, and the sub-lookup is a
            # cell axis in its own right.
            sublk = DD.lookup(sub, DGG.Cells)
            @test sublk isa DGG.CellLookup
            @test DGG.level(sublk) == leaf
            @test collect(sublk) == [lk[k] for k in byhand]
            @test all(DGG.localindex(sublk, sublk[j]) == j for j in eachindex(byhand))
        end
        # The whole region selects the whole axis: the coverage that built the
        # lookup, run against the lookup, is the identity.
        @test length(A[DGG.Cells(DGG.Covering(REGION))]) == length(ids)
        # And a region the axis does not reach selects nothing at all.
        empty_target = GI.Polygon([GI.LinearRing([(-26.0, -41.0), (-24.0, -41.0),
            (-24.0, -39.0), (-26.0, -41.0)])])
        @test isempty(A[DGG.Cells(DGG.Covering(empty_target))])
    end
end

# ---------------------------------------------------------------------------
# Memory: the reason the type exists
#
# One set, two leaf levels. The deeper one names `aperture^deeper` times as
# many cells and stores the same windows, so `Base.summarysize` must not move.
# The comparison against the materialised vector is the same statement with a
# number on it.
# ---------------------------------------------------------------------------

@testset "memory is O(#entries), not O(#leaf cells): $(syslabel(sys))" for
    (sys, leaf, deeper) in SWEEP

    if !DGG.has_sorted_subtrees(sys)
        @test !DGG.has_sorted_subtrees(sys)
        continue
    end

    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf)
    shallow = DGG.CellLookup(set)
    deep = DGG.CellLookup(set; level=leaf + deeper)

    @test length(deep) > length(shallow)
    @test CL.nwindows(CL.windows(deep)) == CL.nwindows(CL.windows(shallow))
    @test Base.summarysize(deep) == Base.summarysize(shallow)

    # Against the thing it replaces. One word per leaf id is the floor for a
    # materialised vector; the lookup is under it, backing included.
    @test Base.summarysize(deep) < 8 * length(deep)

    # Construction does not walk the leaves either: the two builds allocate the
    # same, three levels apart.
    DGG.CellLookup(set)
    DGG.CellLookup(set; level=leaf + deeper)
    @test abs(@allocated(DGG.CellLookup(set; level=leaf + deeper)) -
              @allocated(DGG.CellLookup(set))) < 1024

    # And reading is O(1) in the number of leaves, not merely O(1) amortised.
    deep[1]
    @test @allocated(deep[length(deep)]) <= 64
end

# ---------------------------------------------------------------------------
# A5: the selection-mode decision, exercised
# ---------------------------------------------------------------------------

@testset "A5 stores indices because it has no descendant ranges" begin
    sys = DGG.A5System()
    leaf = 9
    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf)

    # The premise, stated rather than assumed.
    @test !DGG.has_sorted_subtrees(sys)
    @test_throws MethodError DGG.descendant_range(sys, first(set), leaf)
    @test_throws ArgumentError DGG.level_ranges(set, leaf)

    # The decision: the lookup exists anyway, and it is exactly the
    # `descendants` expansion — the pattern `subtree(sys, cell, level)`
    # already uses.
    lk = DGG.CellLookup(set)
    ids = sort!(reduce(vcat, [DGG.descendants(sys, c, leaf) for c in set]))
    @test collect(lk) == ids
    @test all(DGG.localindex(lk, lk[k]) == k for k in eachindex(ids))

    DGG.CellLookup(set)
    DGG.CellLookup(set; level=leaf + 2)
    @test @allocated(DGG.CellLookup(set; level=leaf + 2)) >
          4 * @allocated(DGG.CellLookup(set))
    deep = DGG.CellLookup(set; level=leaf + 2)
    @test length(deep) > length(lk)
    @test Base.summarysize(deep) <= 8 * length(deep)

    A = DD.DimArray(Float64.(eachindex(ids)), DGG.Cells(lk))
    inner = DGG.query(sys, DGG.MultiOrderCoverage(ZURICH); level=leaf)
    byhand = sort!(filter!(!isnothing,
        [DGG.localindex(lk, c) for c in sort!(reduce(vcat,
            [DGG.descendants(sys, c, leaf) for c in inner]))]))
    @test parent(A[DGG.Cells(DGG.Covering(ZURICH))]) == Float64.(byhand)
end

# ---------------------------------------------------------------------------
# The parts of the DimensionalData contract that are not per-system
# ---------------------------------------------------------------------------

@testset "lookup mechanics" begin
    sys = DGG.IGeo7System()
    leaf = 5
    grid = DGG.levelgrid(sys, leaf)
    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf)
    lk = DGG.CellLookup(set)
    ids = collect(lk)

    # Ascending indices keep the windowed form; anything else is not a set of
    # windows and falls back to an ordinary lookup rather than lying about it.
    @test lk[[1, 3, 5]] isa DGG.CellLookup
    @test collect(lk[[1, 3, 5]]) == ids[[1, 3, 5]]
    @test !(lk[[3, 1]] isa DGG.CellLookup)
    @test collect(lk[[3, 1]]) == ids[[3, 1]]
    mask = falses(length(lk));
    mask[2] = mask[4] = true
    @test collect(lk[mask]) == ids[[2, 4]]
    # A neighbour list is a `SmallVector`, and indexing an axis by one is
    # ambiguous unless the tie against SmallCollections' own method is broken.
    @test collect(lk[SmallCollections.SmallVector{8,Int}([3, 4, 5])]) == ids[3:5]
    @test !any(Test.detect_ambiguities(DGG; recursive=true)) do pair
        any(m -> occursin("dimensionaldata.jl", string(m.file)) ||
                     occursin("stencil.jl", string(m.file)), pair)
    end

    # `parent` is the values, so the whole collection surface DimensionalData
    # derives from it is correct without a per-method override.
    @test parent(lk) == ids
    @test DGG.cellset(lk) === set
    @test size(lk) == (length(ids),)
    @test axes(lk) == (Base.OneTo(length(ids)),)
    @test keys(lk) == LinearIndices((Base.OneTo(length(ids)),))
    @test first(lk) == ids[1] && last(lk) == ids[end]
    @test DD.Lookups.order(lk) === DD.Lookups.ForwardOrdered()
    @test DD.Lookups.val(lk) === parent(lk)
    @test DD.Lookups.metadata(lk) === DD.Lookups.NoMetadata()

    # `bounds` is the one caller that asks about an empty axis rather than
    # indexing it, and `first`/`last` of nothing is a BoundsError.
    @test DD.Lookups.bounds(lk) == (ids[1], ids[end])
    @test DD.Lookups.bounds(lk[Int[]]) == (nothing, nothing)

    # Membership without selecting.
    @test DD.Lookups.hasselection(lk, DD.At(ids[3]))
    @test !DD.Lookups.hasselection(lk, DD.At(DGG.cellat(grid, FARAWAY...)))
    @test DD.Lookups.hasselection(lk, DD.Contains(ids[3]))

    @test issorted(ids)
    @test searchsortedfirst(lk, ids[7]) == 7
    @test searchsortedlast(lk, ids[7]) == 7

    @test DD.Lookups.rebuild(lk) === lk
    @test DD.Lookups.rebuild(lk; data=parent(lk)) === lk
    @test DD.Lookups.rebuild(lk; data=ids) == lk
    @test DD.Lookups.rebuild(lk; data=ids) isa DGG.CellLookup
    @test DD.Lookups.rebuild(lk; data=reverse(ids)) isa DD.Lookups.Categorical
    @test_throws ArgumentError DD.Lookups.rebuild(lk; data=[1, 2, 3])
    @test occursin("cat", sprint(showerror,
        try
            DD.Lookups.rebuild(lk; data=[1, 2, 3])
        catch err
            err
        end))

    # Equality is about content, not about which backing produced it.
    @test lk == DGG.CellLookup(DGG.PartialGrid(sys, leaf, ids))
    @test lk != DGG.CellLookup(grid)
    @test DGG.CellLookup(lk) === lk

    # Two levels of the same region are different axes even though the set is
    # the same object.
    @test DGG.CellLookup(set; level=leaf + 1) != lk

    # Show is a summary, not the sixty thousand ids.
    s = sprint(show, lk)
    @test occursin("CellLookup", s)
    @test occursin("IGeo7System", s)
    @test occursin(string(length(ids)), s)
    @test occursin("Covering", sprint(show, DGG.Covering(REGION)))

    # A cube shows its axis the same way.
    A = DD.DimArray(Float64.(eachindex(ids)), DGG.Cells(lk))
    @test occursin("CellLookup", sprint(show, MIME"text/plain"(), DD.dims(A, DGG.Cells)))

    # Constructing a cube whose data does not match the axis is an error, not a
    # silently ragged array.
    @test_throws DimensionMismatch DD.DimArray(zeros(3), DGG.Cells(lk))

    # A level deeper than the set's coarsest cell is refused by the expansion,
    # and a level above it by `level_ranges`' own bound.
    @test_throws ArgumentError DGG.CellLookup(set; level=minimum(DGG.level, set) - 1)
end


@testset "cube operations over a cell axis" begin
    sys = DGG.IGeo7System()
    leaf = 5
    grid = DGG.levelgrid(sys, leaf)
    lk = DGG.CellLookup(DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf))
    ids = collect(lk)
    A = DD.DimArray(Float64.(eachindex(ids)), DGG.Cells(lk); name=:demo)

    # DimensionalData reverses the lookups it knows and falls through to
    # `Base.reverse` on the rest, which answers with a bare vector where a
    # lookup belongs; every selector on the result then recursed on it.
    @testset "reverse" begin
        revA = reverse(A; dims=DGG.Cells)
        revlk = DD.lookup(revA, DGG.Cells)
        @test revlk isa DD.Lookups.Lookup
        @test collect(revlk) == reverse(ids)
        @test parent(revA) == reverse(parent(A))
        @test revA[DGG.Cells(DD.At(ids[3]))] == 3.0
        # A descending axis is not a window set, so it comes back as the same
        # `Categorical` fallback `lk[[3, 1]]` takes. Reversing twice therefore
        # restores the VALUES, not the type: the data and the ids round trip,
        # `A == revrevA` does not, because DimensionalData's `==` on lookups
        # asks for `basetypeof` first.
        revrevA = reverse(revA; dims=DGG.Cells)
        @test parent(revrevA) == parent(A)
        @test collect(DD.lookup(revrevA, DGG.Cells)) == ids
        @test revlk isa DD.Lookups.Categorical
    end

    # A reduction collapses the axis to one element that no cell id names.
    @testset "reductions over the cell axis" begin
        for (f, want) in ((sum, sum(parent(A))), (mean, mean(parent(A))),
            (maximum, maximum(parent(A))), (minimum, minimum(parent(A))))
            r = f(A; dims=DGG.Cells)
            @test size(r) == (1,)
            @test r[1] == want
            @test DD.lookup(r, DGG.Cells) isa DD.Lookups.NoLookup
        end
        @test DD.Lookups.reducelookup(lk) isa DD.Lookups.NoLookup
    end

    # Concatenation hands the joined id vector back through `rebuild`; two
    # ascending disjoint axes join into one, and the result is a cell axis.
    @testset "concatenation" begin
        n = length(ids)
        lo, hi = A[1:(n÷2)], A[(n÷2+1):n]
        for joined in (vcat(lo, hi), cat(lo, hi; dims=DGG.Cells))
            @test length(joined) == n
            joinedlk = DD.lookup(joined, DGG.Cells)
            @test joinedlk isa DGG.CellLookup
            @test collect(joinedlk) == ids
            @test joinedlk == lk
            @test parent(joined) == parent(A)
            @test joined[DGG.Cells(DD.At(ids[end]))] == Float64(n)
        end
    end

    # Replacing the axis wholesale is `set`, and the replacement that works is
    # `NoLookup` — the documented collateral of a lookup with no free fields.
    @testset "replacing the axis" begin
        plain = DD.set(A, DGG.Cells => DD.NoLookup())
        @test DD.lookup(plain, DGG.Cells) isa DD.Lookups.NoLookup
        @test parent(plain) == parent(A)
        @test plain[DGG.Cells(3)] == 3.0
    end
end

# ---------------------------------------------------------------------------
# What a failed selection says. One system is enough: the messages are written
# once, in the lookup, and read the axis through the same accessors everywhere.
# ---------------------------------------------------------------------------

@testset "a failed cell selection is a sentence, not a type" begin
    sys = DGG.HEALPixSystem()
    grid = DGG.levelgrid(sys, 4)
    lk = DGG.CellLookup(DGG.CellVector(grid))
    A = DD.DimArray(collect(1.0:DGG.ncells(grid)), DGG.Cells(lk))
    coarse = DGG.cellindex(DGG.levelgrid(sys, 3), 1)
    foreign = DGG.cellindex(DGG.levelgrid(DGG.IGeo7System(), 4), 5)

    @test_throws "not the axis's level 4" A[DGG.Cells(DD.At(coarse))]
    @test_throws "names cells as LevelIndex" A[DGG.Cells(DD.At(foreign))]
    # The regression: `SelectorError`'s own message prints the lookup's full
    # parameterised type instead of what it holds.
    message = try
        A[DGG.Cells(DD.At(coarse))]
    catch err
        sprint(showerror, err)
    end
    @test !occursin("CellLookup{", message)

    # A cell of another system reaching a geometry verb names both systems
    # rather than surfacing as a `MethodError` about a verb nobody called.
    @test_throws "is a cell of IGeo7System, not of HEALPixSystem" DGG.cell_centroid(grid, foreign)
    @test_throws "is a cell of IGeo7System, not of HEALPixSystem" DGG.children(sys, foreign)

    # `Near` is refused rather than answered in id order.
    @test_throws "nearest id is not the nearest cell" A[DGG.Cells(DD.Near(lk[3]))]
end

end # module DimensionalDataTests
