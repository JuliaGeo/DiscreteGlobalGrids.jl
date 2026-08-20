# Cross-system laws for `CellVector`, including expansion, indexing, selection,
# iteration, and compact storage of multi-order cell sets.

module CellVectorTests

using Test
import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import GeometryOps as GO

include(joinpath(@__DIR__, "..", "..", "helpers.jl"))
using .DGGTestHelpers: syslabel, sweepcovers

const FB = DGG.Fallbacks
const EN = DGG.Engine

# The same Switzerland/Zurich fixture the DimensionalData suite uses, so the
# two files can be read against each other.
const REGION = GI.Polygon([GI.LinearRing([(6.0, 45.8), (10.5, 45.8), (10.5, 47.8),
    (6.0, 47.8), (6.0, 45.8)])])
const ZURICH = GI.Polygon([GI.LinearRing([(8.3, 47.2), (8.8, 47.2), (8.8, 47.6),
    (8.3, 47.6), (8.3, 47.2)])])
const NOWHERE = GI.Polygon([GI.LinearRing([(-26.0, -41.0), (-24.0, -41.0),
    (-24.0, -39.0), (-26.0, -41.0)])])
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

# The oracle, computed without the type under test and through the same trait
# branch it takes, but by a different route.
function expand(sys, set, l::Int)
    grid = DGG.levelgrid(sys, l)
    DGG.has_sorted_subtrees(sys) &&
        return [DGG.cellindex(grid, p) for r in DGG.level_ranges(set, l) for p in r]
    return sort!(reduce(vcat, [DGG.descendants(sys, c, l) for c in set]))
end

nwin(cv) = EN.nwindows(EN.windows(cv))

@testset "the sweep covers every registered system" begin
    sweepcovers(SWEEP)
end

# ---------------------------------------------------------------------------
# The core is reachable, and reachable without a cube
# ---------------------------------------------------------------------------

@testset "the compression is DimensionalData-free" begin
    sys, leaf = DGG.IGeo7System(), 5
    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf)
    cv = DGG.CellVector(set)

    # This module never imported it, and never had to.
    @test !isdefined(@__MODULE__, :DimensionalData)
    @test !isdefined(@__MODULE__, :DD)

    # The type and its verbs live in the engine substrate, not in the DD layer.
    @test parentmodule(DGG.CellVector) === DGG.Engine
    @test DGG.CellVector <: AbstractVector
    @test !(DGG.CellVector <: DGG.CellLookup)
    for f in (DGG.covering, DGG.covering_positions, DGG.cellset)
        @test any(m -> occursin("cell_vector.jl", string(m.file)), methods(f))
    end
    # `covering`/`covering_positions` are answered ONLY in the core: the cube
    # layer selects by calling them, never by reimplementing them.
    for f in (DGG.covering, DGG.covering_positions)
        @test !any(m -> occursin("dimensionaldata.jl", string(m.file)), methods(f))
    end

    # And nothing the core hands back is a DimensionalData type.
    for T in (typeof(cv), eltype(cv), typeof(DGG.PartialGrid(cv)),
        typeof(DGG.covering(cv, ZURICH)), typeof(DGG.cellset(cv)),
        typeof(intersect(cv, cv)))
        @test !startswith(string(parentmodule(T)), "DimensionalData")
    end

    # The cube face is the wrapper, in that direction and not the other: a
    # `CellLookup`'s VALUES are exactly one of these.
    @test parent(DGG.CellLookup(set)) isa DGG.CellVector
    @test parent(DGG.CellLookup(set)) == cv
end

# ---------------------------------------------------------------------------
# The laws, once per system
# ---------------------------------------------------------------------------

@testset "a compressed cell vector: $(syslabel(sys))" for (sys, leaf, _) in SWEEP
    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf)
    grid = DGG.levelgrid(sys, leaf)
    cv = DGG.CellVector(set)
    ids = expand(sys, set, leaf)

    @testset "what it is" begin
        @test cv isa AbstractVector{DGG.cellindextype(sys)}
        @test DGG.level(cv) == leaf
        @test DGG.system(cv) == sys
        @test DGG.cellset(cv) === set
        @test length(cv) == length(ids)
        @test eltype(cv) === DGG.cellindextype(sys)
        # A mixed-level backing standing for a single-level vector is the whole
        # premise; a set that happened to be flat would make the rest vacuous.
        @test length(unique(DGG.level.(collect(set)))) > 1
        @test all(==(leaf), DGG.level.(collect(cv)))
    end

    @testset "the lazy form equals the materialised leaf vector" begin
        @test collect(cv) == ids
        @test [c for c in cv] == ids
        @test first(cv) == first(ids)
        @test last(cv) == last(ids)
        @test cv[end] == last(ids)
        @test cv[begin] == first(ids)
        @test all(cv[k] == ids[k] for k in eachindex(ids))
        @test collect(cv[2:5]) == ids[2:5]
        @test cv[:] === cv
        # An ascending index set stays compressed; anything else is not a set of
        # windows and comes back as the plain vector that can hold it.
        @test cv[[1, 3, 5]] isa DGG.CellVector
        @test collect(cv[[1, 3, 5]]) == ids[[1, 3, 5]]
        @test !(cv[[3, 1]] isa DGG.CellVector)
        @test cv[[3, 1]] == ids[[3, 1]]
        mask = falses(length(cv))
        mask[2] = mask[4] = true
        @test collect(cv[mask]) == ids[[2, 4]]
    end

    @testset "position <-> id round trips" begin
        @test all(DGG.cellposition(cv, cv[k]) == k for k in eachindex(ids))
        outside = DGG.cellat(grid, FARAWAY...)
        @test outside !== nothing
        @test DGG.cellposition(cv, outside) === nothing
        @test !(outside in cv)
        @test all(cv[k] in cv for k in (1, length(ids) ÷ 2, length(ids)))
        # An id from another level names no position here, rather than the
        # position of a cell with the same raw bits.
        coarser = DGG.ancestor(sys, first(ids), leaf - 1)
        @test DGG.cellposition(cv, coarser) === nothing
        @test !(coarser in cv)
    end

    @testset "the five ways in agree" begin
        # A whole level: one window, and the vector's positions ARE the grid's.
        complete = DGG.CellVector(grid)
        @test length(complete) == DGG.ncells(grid)
        @test nwin(complete) == 1
        @test all(complete[k] == DGG.cellindex(grid, k) for k in (1, 7, DGG.ncells(grid)))
        @test all(DGG.cellposition(complete, c) == DGG.cellposition(grid, c) for c in cv)

        # An arbitrary ascending subset, by both spellings.
        @test DGG.CellVector(DGG.PartialGrid(sys, leaf, ids)) == cv
        @test DGG.CellVector(sys, leaf, ids) == cv
        @test collect(DGG.CellVector(sys, leaf, ids)) == ids

        # A rooted subtree, which is the one shape both backings can hold.
        root = DGG.ancestor(sys, first(ids), leaf - 1)
        rooted = DGG.CellVector(DGG.subtree(sys, root, leaf))
        @test collect(rooted) == DGG.descendants(sys, root, leaf)
        @test DGG.cellposition(rooted, first(ids)) !== nothing

        # And the identity.
        @test DGG.CellVector(cv) === cv
    end

    @testset "read as a grid, and back" begin
        pg = DGG.PartialGrid(cv)
        @test DGG.ncells(pg) == length(cv)
        @test DGG.level(pg) == leaf
        @test all(DGG.cellindex(pg, k) == cv[k] for k in eachindex(ids))
        @test all(DGG.cellposition(pg, cv[k]) == k for k in eachindex(ids))
        # The round trip is exact, and it is the same windows rather than a
        # re-derivation of them.
        back = DGG.CellVector(pg)
        @test back == cv
        @test collect(back) == ids
        @test EN.windows(back) === EN.windows(cv)
    end

    @testset "a point, and a region" begin
        # `cell_centroid` is interior to its cell by contract, so asking about a
        # cell's own centroid must come back with that cell — the `cellat` law
        # read through the compression.
        for k in (1, length(ids) ÷ 3, length(ids) ÷ 2, length(ids))
            lon, lat = LONLAT(DGG.cell_centroid(grid, ids[k]))
            @test DGG.cellposition(cv, lon, lat) == k
            @test DGG.cellat(cv, lon, lat) == ids[k]
        end
        @test DGG.cellat(grid, FARAWAY...) !== nothing
        @test DGG.cellat(cv, FARAWAY...) === nothing
        @test DGG.cellposition(cv, FARAWAY...) === nothing

        # `covering` IS the coverage expansion intersected with the vector; the
        # hand-rolled right-hand side is the same sentence spelled out.
        for target in (ZURICH, REGION)
            byhand = Int[]
            for c in expand(sys, DGG.query(sys, DGG.MultiOrderCoverage(target); level=leaf), leaf)
                k = DGG.cellposition(cv, c)
                k === nothing || push!(byhand, k)
            end
            sort!(byhand)
            sub = DGG.covering(cv, target)
            @test sub isa DGG.CellVector
            @test DGG.level(sub) == leaf
            @test DGG.system(sub) == sys
            @test collect(sub) == [cv[k] for k in byhand]
            @test DGG.covering_positions(cv, target) == byhand
            @test all(DGG.cellposition(sub, sub[j]) == j for j in eachindex(byhand))
            # A derived vector has no origin of its own and says so.
            @test DGG.cellset(sub) isa DGG.PartialGrid
        end
        # The coverage that built the vector, run against it, is the identity;
        # and a region it does not reach selects nothing at all.
        @test DGG.covering(cv, REGION) == cv
        @test isempty(DGG.covering(cv, NOWHERE))
        @test isempty(DGG.covering_positions(cv, NOWHERE))
    end
end

# ---------------------------------------------------------------------------
# Set arithmetic over the windows
# ---------------------------------------------------------------------------

@testset "intersect and issubset are O(#windows)" begin
    sys, leaf = DGG.IGeo7System(), 6
    big = DGG.CellVector(DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf))
    small = DGG.CellVector(DGG.query(sys, DGG.MultiOrderCoverage(ZURICH); level=leaf))

    both = intersect(big, small)
    @test both isa DGG.CellVector
    @test DGG.level(both) == leaf
    # Base's `intersect` keeps the left operand's order, which for an ascending
    # vector is ascending — so the compressed answer and the generic one are the
    # same vector, not merely the same set.
    @test collect(both) == intersect(collect(big), collect(small))
    @test !isempty(both)

    @test both ⊆ big
    @test both ⊆ small
    @test big ⊆ big
    @test !(big ⊆ small)
    @test intersect(big, big) == big
    @test isempty(intersect(big, DGG.CellVector(
        DGG.query(sys, DGG.MultiOrderCoverage(NOWHERE); level=leaf))))

    deepbig = DGG.CellVector(DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf); level=leaf + 3)
    deepsmall = DGG.CellVector(DGG.query(sys, DGG.MultiOrderCoverage(ZURICH); level=leaf); level=leaf + 3)
    @test length(deepbig) > 300 * length(big)
    intersect(big, small)
    intersect(deepbig, deepsmall)
    @test abs(@allocated(intersect(deepbig, deepsmall)) -
              @allocated(intersect(big, small))) < 1024

    # Two vectors that do not live in the same space are an error, not a
    # silently empty answer.
    @test_throws ArgumentError intersect(big, DGG.CellVector(
        DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf); level=leaf + 1))
    @test_throws ArgumentError intersect(big, DGG.CellVector(
        DGG.query(DGG.HEALPixSystem(), DGG.MultiOrderCoverage(REGION); level=9)))
end

# ---------------------------------------------------------------------------
# Memory: the reason the type exists, restated without a cube
# ---------------------------------------------------------------------------

@testset "memory is O(#windows): $(syslabel(sys))" for (sys, leaf, deeper) in SWEEP
    if !DGG.has_sorted_subtrees(sys)
        @test !DGG.has_sorted_subtrees(sys)
        continue
    end

    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf)
    shallow = DGG.CellVector(set)
    deep = DGG.CellVector(set; level=leaf + deeper)

    @test length(deep) > length(shallow)
    @test nwin(deep) == nwin(shallow)
    @test Base.summarysize(deep) == Base.summarysize(shallow)

    # Against the thing it replaces. One word per leaf id is the floor for a
    # materialised vector; the compressed form is under it, backing included.
    @test Base.summarysize(deep) < 8 * length(deep)

    # Construction does not walk the leaves either: the two builds allocate the
    # same, `deeper` levels apart.
    DGG.CellVector(set)
    DGG.CellVector(set; level=leaf + deeper)
    @test abs(@allocated(DGG.CellVector(set; level=leaf + deeper)) -
              @allocated(DGG.CellVector(set))) < 1024

    # Reading is O(1) in the number of leaves, not merely O(1) amortised.
    deep[1]
    @test @allocated(deep[length(deep)]) <= 64

    DGG.PartialGrid(shallow)
    DGG.PartialGrid(deep)
    @test abs(@allocated(DGG.PartialGrid(deep)) -
              @allocated(DGG.PartialGrid(shallow))) < 512
    @test @allocated(DGG.PartialGrid(deep)) < 4096
    @test DGG.Helpers.strictly_increasing(deep)
    @test Base.summarysize(DGG.PartialGrid(deep)) < 8 * length(deep)
end

@testset "A5 stores positions because it has no descendant ranges" begin
    sys, leaf = DGG.A5System(), 9
    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf)

    @test !DGG.has_sorted_subtrees(sys)
    @test_throws ArgumentError DGG.level_ranges(set, leaf)

    # The decision: the vector exists anyway, and it is exactly the
    # `descendants` expansion — the pattern `subtree(sys, cell, level)`
    # already uses.
    cv = DGG.CellVector(set)
    ids = sort!(reduce(vcat, [DGG.descendants(sys, c, leaf) for c in set]))
    @test collect(cv) == ids
    @test all(DGG.cellposition(cv, cv[k]) == k for k in eachindex(ids))

    # The cost is the CONSTRUCTION, not the stored form.
    DGG.CellVector(set)
    DGG.CellVector(set; level=leaf + 2)
    @test @allocated(DGG.CellVector(set; level=leaf + 2)) >
          4 * @allocated(DGG.CellVector(set))
    deep = DGG.CellVector(set; level=leaf + 2)
    @test length(deep) > length(cv)
    @test Base.summarysize(deep) <= 8 * length(deep)
end


@testset "the compressed grid route costs O(#windows), not O(#cells)" begin
    sys, leaf = DGG.IGeo7System(), 9
    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf)

    compressed = DGG.PartialGrid(DGG.CellVector(set))
    materialised = DGG.PartialGrid(sys, leaf, DGG.cellindices(set, leaf))

    # Same grid, in the position order a regridder lines up against.
    n = DGG.ncells(materialised)
    @test DGG.ncells(compressed) == n
    @test DGG.level(compressed) == DGG.level(materialised) == leaf
    @test all(DGG.cellindex(compressed, i) == DGG.cellindex(materialised, i)
              for i in (1, 2, n ÷ 3, n ÷ 2, n))
    @test collect(DGG.CellVector(compressed)) == DGG.cellindices(set, leaf)

    windows = nwin(DGG.CellVector(set))
    @test n > 50_000
    @test windows < n ÷ 50

    big = Base.summarysize(materialised)
    small = Base.summarysize(compressed)
    @test small < 200 * windows          # O(#windows), with room for the header
    @test big > 8 * n                    # O(#cells), one word per id at least
    @test big / small > 20               # measured 27.0x on this fixture

    deep = DGG.PartialGrid(DGG.CellVector(set; level=leaf + 3))
    @test DGG.ncells(deep) == 343 * n
    @test Base.summarysize(deep) == small
end


@testset "indexing keeps Base's contract" begin
    sys, leaf = DGG.IGeo7System(), 5
    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf)
    cv = DGG.CellVector(set)
    lk = DGG.CellLookup(set)
    ids = collect(cv)
    n = length(cv)

    @testset "a logical mask must match the axis" begin
        mask = falses(n)
        mask[2] = mask[4] = true
        @test collect(cv[mask]) == ids[[2, 4]]
        @test collect(lk[mask]) == ids[[2, 4]]
        # A mask of every position is the whole vector, not an error.
        @test collect(cv[trues(n)]) == ids
        @test isempty(cv[falses(n)])

        # Short and long both throw, on both faces. The short one is the
        # dangerous shape: `findall` would have answered with a prefix.
        for bad in (trues(n - 1), falses(n - 1), trues(n + 1), falses(n + 1))
            @test_throws BoundsError cv[bad]
            @test_throws BoundsError lk[bad]
        end
        @test (try
            cv[trues(n + 1)]
        catch e
            e
        end).a isa DGG.CellVector
        @test (try
            lk[trues(n + 1)]
        catch e
            e
        end).a isa DGG.CellLookup
        # A mask of the wrong RANK is the same mistake and lands the same way.
        @test_throws BoundsError cv[trues(2, 2)]
    end

    @testset "a shaped index answers the same way whatever its values" begin
        up = [1 3; 2 4]           # ascending down each column
        down = [4 2; 3 1]         # not ascending anywhere

        # The core falls through to Base's generic, which preserves the index's
        # shape — so both matrices come back as matrices, and the answer no
        # longer depends on whether the values happened to ascend.
        for m in (up, down)
            got = cv[m]
            @test got isa AbstractMatrix
            @test size(got) == (2, 2)
            @test got == ids[m]
        end

        # The cube face falls through to DimensionalData's own `getindex`, which
        # rebuilds the lookup around the result — and a `Lookup` is
        # one-dimensional, so there is no lookup a matrix could be. It refuses,
        # with the message `rebuild` gives. The law is the CONSISTENCY: the
        # ascending matrix and the descending one are answered identically,
        # where before one flattened to a `CellLookup` and the other did not.
        for m in (up, down)
            @test_throws ArgumentError lk[m]
        end

        @test cv[[1, 3, 5]] isa DGG.CellVector
        @test lk[[1, 3, 5]] isa DGG.CellLookup
        @test cv[2:5] isa DGG.CellVector
        @test lk[2:5] isa DGG.CellLookup
    end

    @testset "an out-of-range index reports the collection, as Base does" begin
        err = try
            cv[[1, n + 1]]
        catch e
            e
        end
        @test err isa BoundsError
        @test err.i == ([1, n + 1],)
    end
end

end # module CellVectorTests
