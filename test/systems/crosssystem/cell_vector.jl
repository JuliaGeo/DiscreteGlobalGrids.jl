# ---------------------------------------------------------------------------
# T19 — the compressed cell collection, WITHOUT DimensionalData.
#
# `CellVector` is what `CellLookup` is made of: a strictly ascending run of
# cells at one level, stored as leaf position windows and answering as the id
# vector. T16 shipped that arithmetic inside the DimensionalData layer, which
# meant a regridder, a chunker or any plain-`Array` caller had to materialise
# the ids to get at it — 892,568 bytes where the cube got 33,024. This file is
# the law that says the compression is available with no cube in sight.
#
# The module below deliberately does NOT import DimensionalData, and asserts
# that it has not: everything is reached through `DiscreteGlobalGrids`' own
# exports. `DimensionalData` is still a hard dependency of the package, so it IS
# loaded in the session — what is being pinned is that no *caller* needs it, and
# that no verb here is answered by a method in `src/dimensionaldata.jl`.
#
# The laws, in order:
#
#   * DD-FREEDOM — the type, its verbs and everything they return live under
#     `src/fallbacks/`, and none of their methods come from the DD layer.
#   * CONSTRUCTION — the five ways in (a multi-order set, a whole level, a
#     rooted subtree, a `PartialGrid`, an explicit id vector) all agree, on
#     every system, including A5's selection mode.
#   * EQUIVALENCE — the lazy form answers exactly what the materialised leaf
#     vector would, and `cellposition` is its inverse.
#   * POINT AND REGION — `cellat`/`cellposition` on a lon/lat pair, `in` for
#     membership, `covering` for a region: the three selector questions, asked
#     without a selector.
#   * MEMORY — re-expanding one set three levels deeper multiplies the cells by
#     the aperture cubed and must not move `Base.summarysize`. Excluded on A5
#     for the reason T16 states: selection mode walks the leaves to build.
#   * THE PROBE — the motivating measurement, as an assertion: the compressed
#     grid route costs O(#windows) where the `cellindices` route costs
#     O(#cells), and the ratio on the fixture is pinned.
# ---------------------------------------------------------------------------

module CellVectorTests

using Test
import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import GeometryOps as GO

const FB = DGG.Fallbacks

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

sysname(sys) = sys isa DGG.AuthalicSystem ?
               "Authalic($(nameof(typeof(parent(sys)))))" : string(nameof(typeof(sys)))

# The oracle, computed without the type under test and through the same trait
# branch it takes, but by a different route.
function expand(sys, set, l::Int)
    grid = DGG.levelgrid(sys, l)
    DGG.has_sorted_subtrees(sys) &&
        return [DGG.cellindex(grid, p) for r in DGG.level_ranges(set, l) for p in r]
    return sort!(reduce(vcat, [DGG.descendants(sys, c, l) for c in set]))
end

nwin(cv) = FB.nwindows(FB.windows(cv))

@testset "the sweep covers every registered system" begin
    swept = Set(typeof(s) for (s, _, _) in SWEEP)
    for s in DGG.systems()
        @test typeof(s) in swept
    end
    @test any(s -> s isa DGG.AuthalicSystem, first.(SWEEP))
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

    # The type and its verbs live in the fallback substrate, not in the DD layer.
    @test parentmodule(DGG.CellVector) === DGG.Fallbacks
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

@testset "a compressed cell vector: $(sysname(sys))" for (sys, leaf, _) in SWEEP
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
        rooted = DGG.CellVector(DGG.PartialGrid(sys, root, leaf))
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
        @test FB.windows(back) === FB.windows(cv)
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
        # A point the vector does not reach answers `nothing`, not the cell of
        # the level grid that happens to hold it.
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

    # The point of doing it over intervals: the same operands three levels
    # deeper name 343 times as many cells and cost the same to intersect.
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

@testset "memory is O(#windows): $(sysname(sys))" for (sys, leaf, deeper) in SWEEP
    if !DGG.has_sorted_subtrees(sys)
        # EXCLUDED for the reason T16 pins on the lookup: selection mode
        # materialises one position per leaf to BUILD, so nothing bounds its
        # construction by the entry count even where the compression afterwards
        # happens to. The A5 testset below states the decision from the other
        # side.
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

    # And so is the handshake, which is the whole point of it being O(1): the
    # ascent check `PartialGrid` runs on an arbitrary id vector short-circuits
    # here, so building a grid over twenty million cells costs what building one
    # over sixty thousand costs. The invariance is the law; the absolute figure
    # is a few hundred bytes of constructor and varies with the id codec, so it
    # is bounded rather than pinned.
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
    # `descendants` expansion — the pattern `PartialGrid(sys, cell, level)`
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

# ---------------------------------------------------------------------------
# The motivating probe, as a law
#
# Before T19 the only way to hand a coverage to something that is not a cube was
# `cellindices(set, l)` — one typed id per leaf cell, and a `PartialGrid` over
# it. The compressed route now gets there through `CellVector`, which is the
# same grid, cell for cell, in a fraction of the memory. Both halves are
# asserted: same answers, and the ratio.
# ---------------------------------------------------------------------------

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

    # The fixture, pinned so that a regression in the compression is visible as
    # a number rather than as a slowdown: ~666 windows for ~60,000 leaf cells.
    windows = nwin(DGG.CellVector(set))
    @test n > 50_000
    @test windows < n ÷ 50

    big = Base.summarysize(materialised)
    small = Base.summarysize(compressed)
    @test small < 200 * windows          # O(#windows), with room for the header
    @test big > 8 * n                    # O(#cells), one word per id at least
    @test big / small > 20               # measured 27.0x on this fixture

    # The same claim in the direction that matters most: re-expanding to a level
    # deep enough that materialising is out of the question costs the same.
    deep = DGG.PartialGrid(DGG.CellVector(set; level=leaf + 3))
    @test DGG.ncells(deep) == 343 * n
    @test Base.summarysize(deep) == small
end

end # module CellVectorTests
