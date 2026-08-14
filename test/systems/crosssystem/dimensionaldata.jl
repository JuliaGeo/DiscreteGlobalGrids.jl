# ---------------------------------------------------------------------------
# T16 — the DimensionalData cell axis, on every system.
#
# A `CellLookup` claims to BE the leaf id vector while storing the multi-order
# set it came from. Everything below is one of the two halves of that claim:
#
#   * EQUIVALENCE — the lazy form answers exactly what the materialised leaf
#     vector would. `collect(lk)`, `lk[k]`, iteration, `first`/`last`, and the
#     three degenerate constructions (a whole level, a rooted subset, an
#     arbitrary ascending subset) are all checked against an expansion computed
#     independently of the lookup, through the trait branch the lookup itself
#     is not allowed to skip.
#   * ROUND TRIP — `cellposition(lk, lk[k]) == k` for every k, and `nothing`
#     for a cell the lookup does not hold, including one at another level. That
#     is the bijection every selector lands on.
#   * PARTIAL ⊂ COMPLETE — a lookup over a subset and a lookup over the whole
#     level agree about every cell they share: same ids, and positions related
#     by the complete grid's own `cellposition`.
#   * SELECTORS — `At(c)` is the position of `c`; `Contains(lon, lat)` on a
#     cell's own centroid is that cell's position, which is the `cellat`
#     contract read through the lookup; `Covering(region)` is exactly the
#     coverage expansion intersected with the lookup, computed here by hand.
#   * MEMORY — the point of the type. Re-expanding one set to three levels
#     deeper multiplies the cells by the aperture cubed and must not move
#     `Base.summarysize` at all, because the stored windows are the same
#     windows. That law is EXCLUDED on A5 with a reason; see below.
#
# A5 has `has_sorted_subtrees == false`, so `level_ranges` throws and the
# windowed backing is unavailable. The decision T16 took is SELECTION MODE:
# the lookup is built by naming the leaves through `descendants`, resolving
# them to positions and compressing what comes out. Every law above still holds
# there and is run there; the memory law is the one that cannot, because the
# construction walks the leaves whatever the compression then finds. Its
# testset states that as an exclusion and pins the decision from the other
# side — that the A5 lookup IS the `descendants` expansion, and that building
# it at sixteen times the depth costs an order of magnitude more.
# ---------------------------------------------------------------------------

module DimensionalDataTests

using Test
import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GeoInterface as GI
import GeometryOps as GO
const CL = DGG.CellLookups

# ---------------------------------------------------------------------------
# Fixtures
#
# A Switzerland-sized box and a Zurich-sized one strictly inside it, so that
# `Covering` has both a non-empty and a proper-subset answer to give, plus a
# point in the South Atlantic that no lookup here can hold.
# ---------------------------------------------------------------------------

const REGION = GI.Polygon([GI.LinearRing([(6.0, 45.8), (10.5, 45.8), (10.5, 47.8),
    (6.0, 47.8), (6.0, 45.8)])])
const ZURICH = GI.Polygon([GI.LinearRing([(8.3, 47.2), (8.8, 47.2), (8.8, 47.6),
    (8.3, 47.6), (8.3, 47.2)])])
const FARAWAY = (-25.0, -40.0)

const LONLAT = GO.UnitSpherical.GeographicFromUnitSphere()

# Levels chosen so that a leaf cell is a few kilometres across on all seven —
# the apertures are 7, 7, 4, 4, 4, 4, so a fixed level is not comparable. The
# third column is how much deeper the memory law re-expands the same set.
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

@testset "the sweep covers every registered system" begin
    swept = Set(typeof(s) for (s, _, _) in SWEEP)
    for s in DGG.systems()
        @test typeof(s) in swept
    end
    @test any(s -> s isa DGG.AuthalicSystem, first.(SWEEP))
end

# The oracle: a multi-order set expanded to leaf ids WITHOUT the lookup, taking
# the same trait branch the lookup takes but by a different route — positions
# from `level_ranges` where subtrees are sorted, ids from `descendants` where
# they are not. `sort` is what makes the two routes comparable: only the first
# is in curve order already.
function expand(sys, set, l::Int)
    grid = DGG.levelgrid(sys, l)
    DGG.has_sorted_subtrees(sys) &&
        return [DGG.cellindex(grid, p) for r in DGG.level_ranges(set, l) for p in r]
    return sort!(reduce(vcat, [DGG.descendants(sys, c, l) for c in set]))
end

# ---------------------------------------------------------------------------
# The laws, once per system
# ---------------------------------------------------------------------------

@testset "a multi-order cell axis: $(sysname(sys))" for (sys, leaf, deeper) in SWEEP
    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf)
    grid = DGG.levelgrid(sys, leaf)
    lk = DGG.CellLookup(set)
    ids = expand(sys, set, leaf)

    @testset "what it is" begin
        @test lk isa DD.Lookups.Lookup
        @test DGG.level(lk) == leaf
        @test DGG.system(lk) == sys
        @test parent(lk) === set
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

    @testset "position <-> id round trips" begin
        @test all(DGG.cellposition(lk, lk[k]) == k for k in eachindex(ids))
        outside = DGG.cellat(grid, FARAWAY...)
        @test outside !== nothing
        @test DGG.cellposition(lk, outside) === nothing
        # An id from another level names no position here, rather than the
        # position of a cell with the same raw bits.
        coarser = DGG.ancestor(sys, first(ids), leaf - 1)
        @test DGG.cellposition(lk, coarser) === nothing
    end

    @testset "the degenerate cases answer identically" begin
        # A whole level: one window, and the lookup's positions ARE the grid's.
        complete = DGG.CellLookup(grid)
        @test length(complete) == DGG.ncells(grid)
        @test CL.nwindows(CL.windows(complete)) == 1
        @test all(complete[k] == DGG.cellindex(grid, k) for k in (1, 7, DGG.ncells(grid)))

        # PARTIAL ⊂ COMPLETE: the subset lookup and the whole-level lookup
        # agree about every cell of the subset, and the whole-level lookup's
        # positions are the grid's own.
        @test all(DGG.cellposition(complete, c) == DGG.cellposition(grid, c) for c in lk)
        @test all(DGG.cellposition(complete, c) !== nothing for c in lk)

        # An arbitrary ascending subset, the "partial lookup": same content,
        # same answers, and equal to the multi-order form it was built from.
        partial = DGG.CellLookup(DGG.PartialGrid(sys, leaf, ids))
        @test collect(partial) == ids
        @test partial == lk
        @test all(DGG.cellposition(partial, ids[k]) == k for k in eachindex(ids))

        # A rooted subtree, which is the one shape both backings can hold.
        root = DGG.ancestor(sys, first(ids), leaf - 1)
        rooted = DGG.CellLookup(DGG.PartialGrid(sys, root, leaf))
        @test collect(rooted) == DGG.descendants(sys, root, leaf)
        @test DGG.cellposition(rooted, first(ids)) !== nothing
    end

    @testset "the lookup read as a grid" begin
        pg = DGG.PartialGrid(lk)
        @test DGG.ncells(pg) == length(lk)
        @test all(DGG.cellindex(pg, k) == lk[k] for k in eachindex(ids))
        @test all(DGG.cellposition(pg, lk[k]) == k for k in eachindex(ids))
        @test DGG.level(pg) == leaf
    end
end

# ---------------------------------------------------------------------------
# Selectors, and the cube they run against
# ---------------------------------------------------------------------------

@testset "selectors on a cell axis: $(sysname(sys))" for (sys, leaf, _) in SWEEP
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

    @testset "At(id) selects that id's position" begin
        for k in (1, length(ids) ÷ 2, length(ids))
            @test A[DGG.Cells(DD.At(ids[k]))] == Float64(k)
        end
        outside = DGG.cellat(grid, FARAWAY...)
        @test_throws DD.Lookups.SelectorError A[DGG.Cells(DD.At(outside))]
    end

    # `cell_centroid` is interior to its cell by contract, so a point selector
    # asked about a cell's own centroid must come back with that cell — this is
    # the `cellat` law read through the lookup, and the composition is where it
    # could be lost.
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
                k = DGG.cellposition(lk, c)
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
            @test all(DGG.cellposition(sublk, sublk[j]) == j for j in eachindex(byhand))
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

@testset "memory is O(#entries), not O(#leaf cells): $(sysname(sys))" for
    (sys, leaf, deeper) in SWEEP

    if !DGG.has_sorted_subtrees(sys)
        # EXCLUDED, with the reason and the decision both pinned rather than
        # skipped: selection mode materialises one position per leaf to build
        # the lookup, so nothing bounds its cost by the entry count even where
        # the compression afterwards happens to. What CAN be asserted is that
        # the decision was taken — see the A5 testset below.
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

@testset "A5 stores positions because it has no descendant ranges" begin
    sys = DGG.A5System()
    leaf = 9
    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf)

    # The premise, stated rather than assumed.
    @test !DGG.has_sorted_subtrees(sys)
    @test_throws MethodError DGG.descendant_range(sys, first(set), leaf)
    @test_throws ArgumentError DGG.level_ranges(set, leaf)

    # The decision: the lookup exists anyway, and it is exactly the
    # `descendants` expansion — the pattern `PartialGrid(sys, cell, level)`
    # already uses.
    lk = DGG.CellLookup(set)
    ids = sort!(reduce(vcat, [DGG.descendants(sys, c, leaf) for c in set]))
    @test collect(lk) == ids
    @test all(DGG.cellposition(lk, lk[k]) == k for k in eachindex(ids))

    # The cost of the decision, pinned so that a later change to it is visible.
    # It is the CONSTRUCTION that walks the leaves, not the stored form: sixteen
    # times the cells cost an order of magnitude more to build, and then
    # compress to the same windows, because a coverage of a connected region
    # happens to have contiguous runs even where the subtrees are unsorted.
    # That "happens to" is the whole reason the memory law above is excluded
    # here rather than merely restated: nothing guarantees it.
    DGG.CellLookup(set)
    DGG.CellLookup(set; level=leaf + 2)
    @test @allocated(DGG.CellLookup(set; level=leaf + 2)) >
          4 * @allocated(DGG.CellLookup(set))
    deep = DGG.CellLookup(set; level=leaf + 2)
    @test length(deep) > length(lk)
    @test Base.summarysize(deep) <= 8 * length(deep)

    # The consequence the decision inherits rather than causes: an A5 cell's
    # descendants need not lie inside its own footprint, so a coverage's leaf
    # expansion — and therefore a `Covering` selection — over-covers. The
    # selector is still exactly that expansion, which is what is asserted here;
    # the over-covering itself belongs to `MultiOrderCoverage` and is measured
    # in test/systems/crosssystem/multiorder_polygons.jl.
    A = DD.DimArray(Float64.(eachindex(ids)), DGG.Cells(lk))
    inner = DGG.query(sys, DGG.MultiOrderCoverage(ZURICH); level=leaf)
    byhand = sort!(filter!(!isnothing,
        [DGG.cellposition(lk, c) for c in sort!(reduce(vcat,
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

    # `parent` is the backing, deliberately, so the collection surface that
    # DimensionalData derives from `parent` is answered by the lookup instead.
    @test parent(lk) === set
    @test size(lk) == (length(ids),)
    @test axes(lk) == (Base.OneTo(length(ids)),)
    @test DD.Lookups.order(lk) === DD.Lookups.ForwardOrdered()
    @test DD.Lookups.val(lk) === lk
    @test DD.Lookups.metadata(lk) === DD.Lookups.NoMetadata()

    # Ids ascend with position on a complete level grid, which is what makes
    # the binary searches sound; `searchsortedfirst` is that fact, exposed.
    @test issorted(ids)
    @test searchsortedfirst(lk, ids[7]) == 7
    @test searchsortedlast(lk, ids[7]) == 7

    # A lookup has no properties to vary, so a rebuild that would change its
    # values is refused rather than silently ignored.
    @test DD.Lookups.rebuild(lk) === lk
    @test_throws ArgumentError DD.Lookups.rebuild(lk; data=ids)

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

end # module DimensionalDataTests
