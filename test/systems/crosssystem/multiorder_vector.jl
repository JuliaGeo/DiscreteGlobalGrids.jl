# ---------------------------------------------------------------------------
# `MultiOrderVector`, the mixed-level cell container. Laws checked here:
#
#   * round trip: a coverage read as storage keeps its cells, order and
#     reference level, and expands to the set's own `CellVector` at every level.
#   * membership: `cellposition` is exact; `covering_position` resolves every
#     covered leaf (or deeper cell) to its stored ancestor.
#   * point location: a stored cell's centroid resolves back to that cell.
#   * set algebra: against a brute force over leaf-position `Set`s, which
#     shares no code with the interval index.
#   * normalization: results are the coarsest cells that tile them — no
#     complete sibling family survives in any result.
#   * validation: overlap, duplicates and a too-shallow reference level throw;
#     unsorted input is sorted.
#
# A5 has no `has_sorted_subtrees`, hence no container; the last testset asserts
# the refusal.
# ---------------------------------------------------------------------------

module MultiOrderVectorTests

using Test
import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import GeometryOps as GO

const EN = DGG.Engine

# Same fixture regions as the `CellVector` suite.
const REGION = GI.Polygon([GI.LinearRing([(6.0, 45.8), (10.5, 45.8), (10.5, 47.8),
    (6.0, 47.8), (6.0, 45.8)])])
const ZURICH = GI.Polygon([GI.LinearRing([(8.3, 47.2), (8.8, 47.2), (8.8, 47.6),
    (8.3, 47.6), (8.3, 47.2)])])
const NOWHERE = GI.Polygon([GI.LinearRing([(-26.0, -41.0), (-24.0, -41.0),
    (-24.0, -39.0), (-26.0, -41.0)])])
const FARAWAY = (-25.0, -40.0)

# Overlapping continent-sized boxes for the set algebra — big enough to be
# non-trivial at the shallow oracle level.
const BOXA = GI.Polygon([GI.LinearRing([(-20.0, 20.0), (40.0, 20.0), (40.0, 60.0),
    (-20.0, 60.0), (-20.0, 20.0)])])
const BOXB = GI.Polygon([GI.LinearRing([(10.0, 0.0), (70.0, 0.0), (70.0, 40.0),
    (10.0, 40.0), (10.0, 0.0)])])

const LONLAT = GO.UnitSpherical.GeographicFromUnitSphere()

# Per-system leaf levels (a leaf is a few km across on each), deep enough that
# the REGION coverage comes back mixed-level. A5 has no container; see the last
# testset.
const SWEEP = [
    (DGG.IGeo7System(), 5),
    (DGG.H3System(), 4),
    (DGG.HEALPixSystem(), 9),
    (DGG.S2System(), 9),
    (DGG.ISEA4RSystem(), 9),
    (DGG.AuthalicSystem(DGG.IGeo7System()), 5),
]

sysname(sys) = sys isa DGG.AuthalicSystem ?
               "Authalic($(nameof(typeof(parent(sys)))))" : string(nameof(typeof(sys)))

# Shallowest level with at least `atleast` cells, so the oracle size is
# comparable across apertures.
function shallow_level(sys; atleast=400)
    for l in DGG.levels(sys)
        DGG.ncells(DGG.levelgrid(sys, l)) >= atleast && return l
    end
    return last(DGG.levels(sys))
end

# Leaf positions a container names, via `descendant_range` — independent of
# the container's own index.
leafset(mov, l) = Set(p for c in mov for p in DGG.descendant_range(DGG.system(mov), c, l))

# Intervals sorted and pairwise disjoint.
function disjoint_and_sorted(mov)
    ivs = [DGG.descendant_range(DGG.system(mov), c, EN.reference_level(mov)) for c in mov]
    return issorted(ivs; by=first) &&
           all(k -> first(ivs[k]) > last(ivs[k-1]), 2:length(ivs))
end

# True if a complete sibling family is stored, i.e. the container is not
# coarsest.
function has_complete_family(sys, mov)
    stored = Set(collect(mov))
    top = first(DGG.levels(sys))
    for c in mov
        DGG.level(c) > top || continue
        all(k -> k in stored, DGG.children(sys, parent(sys, c))) && return true
    end
    return false
end

# Deterministic probes: ends, middles, and the shallowest and deepest stored
# cells.
function probes(mov)
    n = length(mov)
    return unique([1, cld(n, 3), cld(n, 2), n,
        argmin(k -> DGG.level(mov[k]), eachindex(mov)),
        argmax(k -> DGG.level(mov[k]), eachindex(mov))])
end

@testset "the sweep covers every registered system" begin
    swept = Set(typeof(s) for (s, _) in SWEEP)
    for s in DGG.systems()
        # The sweep must track `has_sorted_subtrees` exactly.
        @test (typeof(s) in swept) == DGG.has_sorted_subtrees(s)
    end
    @test any(s -> s isa DGG.AuthalicSystem, first.(SWEEP))
end

# ---------------------------------------------------------------------------
# The container's own laws, once per system
# ---------------------------------------------------------------------------

@testset "a mixed-level container: $(sysname(sys))" for (sys, leaf) in SWEEP
    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=leaf)
    cells = collect(set)
    mov = DGG.MultiOrderVector(set)
    grid = DGG.levelgrid(sys, leaf)
    top = first(DGG.levels(sys))
    ks = probes(mov)

    @testset "what it is" begin
        @test mov isa AbstractVector{DGG.cellindextype(sys)}
        @test DGG.system(mov) === sys
        @test EN.reference_level(mov) == leaf
        @test length(mov) == length(set)
        @test collect(mov) == cells
        @test mov[:] === mov
        @test mov == mov
        # The fixture must be genuinely mixed-level.
        @test length(unique(DGG.level.(cells))) > 1
        @test maximum(DGG.level, cells) == leaf
        @test disjoint_and_sorted(mov)
        # `offsets` is the cumulative leaf count at the reference level and
        # must agree with the expansion.
        @test mov.offsets[end] == length(DGG.CellVector(mov))
        @test mov.offsets == cumsum(mov.stops .- mov.starts .+ 1)
        @test occursin("MultiOrderVector", sprint(show, mov))
        @test occursin("ref $leaf", sprint(show, mov))
    end

    @testset "the CellVector bridge" begin
        for l in (leaf, leaf + 1)
            bridged = DGG.CellVector(mov; level=l)
            direct = DGG.CellVector(set; level=l)
            @test bridged == direct
            # Same windows, not just same cells: unmerged adjacent subtrees
            # would still compare equal.
            @test EN.windows(bridged).starts == EN.windows(direct).starts
            @test EN.windows(bridged).stops == EN.windows(direct).stops
        end
        @test DGG.cellset(DGG.CellVector(mov)) === mov
        # Expanding above the deepest stored cell would need ancestors, which
        # cover more.
        @test_throws ArgumentError DGG.CellVector(mov; level=leaf - 1)
    end

    @testset "membership is exact" begin
        @test all(DGG.cellposition(mov, mov[k]) == k for k in eachindex(mov))
        @test all(c in mov for c in mov)

        # A descendant of a stored cell is covered but NOT stored — kills
        # `cellposition` answering with the covering position.
        coarse = argmin(k -> DGG.level(mov[k]), eachindex(mov))
        kid = first(DGG.children(sys, mov[coarse]))
        @test DGG.level(kid) <= leaf
        @test DGG.cellposition(mov, kid) === nothing
        @test !(kid in mov)

        # An ancestor of a stored cell is neither stored nor covered.
        deepest = argmax(k -> DGG.level(mov[k]), eachindex(mov))
        anc = DGG.ancestor(sys, mov[deepest], DGG.level(mov[deepest]) - 1)
        @test DGG.cellposition(mov, anc) === nothing
        @test EN.covering_position(mov, anc) === nothing

        # A leaf-grid cell outside the container entirely.
        outside = DGG.cellat(grid, FARAWAY...)
        @test outside !== nothing
        @test DGG.cellposition(mov, outside) === nothing
        @test EN.covering_position(mov, outside) === nothing
    end

    @testset "the covering ancestor" begin
        # Every leaf under a stored cell resolves to it; sampled at each
        # probe's range ends and middle.
        for i in ks
            r = DGG.descendant_range(sys, mov[i], leaf)
            for p in (first(r), (first(r) + last(r)) ÷ 2, last(r))
                @test EN.covering_position(mov, DGG.cellindex(grid, p)) == i
            end
        end
        # A cell deeper than the reference level resolves through its
        # reference-level ancestor.
        deep = first(DGG.children(sys, DGG.cellindex(grid, mov.stops[end])))
        @test DGG.level(deep) == leaf + 1
        @test EN.covering_position(mov, deep) == length(mov)
        @test DGG.cellposition(mov, deep) === nothing
    end

    @testset "a point lands in the cell holding it" begin
        # `cell_centroid` is strictly interior by contract, so a stored cell's
        # centroid must come back as that cell.
        for i in ks
            c = mov[i]
            lon, lat = LONLAT(DGG.cell_centroid(DGG.levelgrid(sys, DGG.level(c)), c))
            @test DGG.cellat(mov, lon, lat) == c
            @test DGG.cellposition(mov, lon, lat) == i
        end
        @test DGG.cellat(grid, FARAWAY...) !== nothing
        @test DGG.cellat(mov, FARAWAY...) === nothing
        @test DGG.cellposition(mov, FARAWAY...) === nothing
    end

    @testset "a region selects whole cells" begin
        for target in (ZURICH, REGION)
            covered = Set(p for r in DGG.level_ranges(
                DGG.query(sys, DGG.MultiOrderCoverage(target); level=leaf), leaf) for p in r)
            byhand = [i for i in eachindex(mov)
                      if any(in(covered), DGG.descendant_range(sys, mov[i], leaf))]
            @test DGG.covering_positions(mov, target) == byhand
            sub = DGG.covering(mov, target)
            @test sub isa DGG.MultiOrderVector
            @test EN.reference_level(sub) == leaf
            # The selection is a sub-vector of stored cells, never a re-cut of
            # the region.
            @test collect(sub) == [mov[i] for i in byhand]
            @test all(DGG.cellposition(sub, sub[j]) == j for j in eachindex(sub))
        end
        @test DGG.covering(mov, REGION) == mov
        @test isempty(DGG.covering(mov, NOWHERE))
        @test isempty(DGG.covering_positions(mov, NOWHERE))
    end

    @testset "indexing takes CellVector's fork" begin
        @test mov[[1, 3, 5]] isa DGG.MultiOrderVector
        @test collect(mov[[1, 3, 5]]) == cells[[1, 3, 5]]
        @test !(mov[[3, 1]] isa DGG.MultiOrderVector)
        @test mov[[3, 1]] == cells[[3, 1]]
        mask = falses(length(mov))
        mask[2] = mask[4] = true
        @test collect(mov[mask]) == cells[[2, 4]]
        # A wrong-length mask is a `BoundsError`, not a shorter answer.
        @test_throws BoundsError mov[falses(length(mov) - 1)]
        @test_throws BoundsError mov[[1, length(mov) + 1]]
        # A subset is a container in its own right, not a view with stale keys.
        sub = mov[2:4]
        @test sub isa DGG.MultiOrderVector
        @test all(DGG.cellposition(sub, sub[j]) == j for j in eachindex(sub))
        @test disjoint_and_sorted(sub)
    end

    @testset "the reference level is a unit, not content" begin
        deeper = DGG.MultiOrderVector(sys, cells; reference_level=leaf + 1)
        @test EN.reference_level(deeper) == leaf + 1
        @test deeper.starts != mov.starts        # ... different integers ...
        @test deeper == mov && mov == deeper     # ... same container
        @test collect(deeper) == cells
        @test DGG.CellVector(deeper) == DGG.CellVector(set; level=leaf + 1)
    end

    @testset "construction from a loose vector" begin
        @test DGG.MultiOrderVector(sys, cells) == mov
        @test DGG.MultiOrderVector(mov) === mov
        # Unsorted input is sorted, not refused.
        @test collect(DGG.MultiOrderVector(sys, reverse(cells))) == cells
        @test collect(DGG.MultiOrderVector(sys, vcat(cells[2:2:end], cells[1:2:end]))) == cells

        # Overlap is refused: ancestor/descendant pairs and duplicates.
        coarse = argmin(k -> DGG.level(mov[k]), eachindex(mov))
        @test_throws ArgumentError DGG.MultiOrderVector(sys,
            vcat(cells, [first(DGG.children(sys, mov[coarse]))]))
        @test_throws ArgumentError DGG.MultiOrderVector(sys, vcat(cells, [cells[1]]))
        # A reference level shallower than a stored cell has no interval for it.
        @test_throws ArgumentError DGG.MultiOrderVector(sys, cells; reference_level=leaf - 1)
        @test_throws ArgumentError DGG.MultiOrderVector(sys, cells;
            reference_level=last(DGG.levels(sys)) + 1)
    end

    @testset "geometry of a mixed-level container" begin
        polys = DGG.cell_polygons(mov)
        @test length(polys) == length(mov)
        for i in ks
            c = mov[i]
            @test GI.coordinates(polys[i]) ==
                  GI.coordinates(DGG.cell_polygon(DGG.levelgrid(sys, DGG.level(c)), c))
        end
    end

    @testset "an empty container is legal" begin
        empty = DGG.MultiOrderVector(sys, DGG.cellindextype(sys)[])
        @test isempty(empty)
        @test EN.reference_level(empty) == top
        @test DGG.cellposition(empty, cells[1]) === nothing
        @test EN.covering_position(empty, cells[1]) === nothing
        @test isempty(DGG.CellVector(empty))
        @test isempty(DGG.covering(empty, REGION))
        @test occursin("0 cells", sprint(show, empty))
    end
end

# ---------------------------------------------------------------------------
# Set algebra, against a brute force over leaf position sets
# ---------------------------------------------------------------------------

@testset "set algebra: $(sysname(sys))" for (sys, _) in SWEEP
    l = shallow_level(sys)
    n = DGG.ncells(DGG.levelgrid(sys, l))
    a = DGG.MultiOrderVector(DGG.query(sys, DGG.MultiOrderCoverage(BOXA); level=l))
    b = DGG.MultiOrderVector(DGG.query(sys, DGG.MultiOrderCoverage(BOXB); level=l))
    A, B = leafset(a, l), leafset(b, l)

    # The boxes must genuinely overlap or three of the four laws are vacuous.
    @test !isempty(intersect(A, B))
    @test !isempty(setdiff(A, B))
    @test !isempty(setdiff(B, A))

    @testset "$name" for (name, got, want) in (
        ("union", union(a, b), union(A, B)),
        ("intersect", intersect(a, b), intersect(A, B)),
        ("setdiff", setdiff(a, b), setdiff(A, B)),
        ("complement", EN.complement(a), setdiff(Set(1:n), A)))

        @test got isa DGG.MultiOrderVector
        @test DGG.system(got) === sys
        @test EN.reference_level(got) == l
        @test leafset(got, l) == want
        @test disjoint_and_sorted(got)
        # Coarsest decomposition: no complete sibling family survives.
        @test !has_complete_family(sys, got)
    end

    c = EN.complement(a)
    @test isempty(intersect(a, c))
    @test leafset(union(a, c), l) == Set(1:n)
    # The whole sphere normalizes to the root cells.
    @test collect(union(a, c)) == collect(DGG.rootcells(sys))
    # Double complement normalizes rather than round-trips: same leaves, as the
    # coarsest cells that name them.
    @test leafset(EN.complement(c), l) == A
    @test EN.complement(c) == union(a, a)
    @test !has_complete_family(sys, EN.complement(c))
    # Identity only when the operand was already minimal; aperture-7 coverages
    # can carry complete sibling families, which normalize away.
    if has_complete_family(sys, a)
        @test length(EN.complement(c)) < length(a)
    else
        @test EN.complement(c) == a
    end

    # The n-ary forms fold the binary one (Base's would build a `Set`).
    @test union(a, b, a) == union(a, b)
    @test intersect(a, b, a) == intersect(a, b)
    @test setdiff(a, b, b) == setdiff(a, b)

    @testset "operands are re-keyed to the deeper reference level" begin
        deep = DGG.MultiOrderVector(DGG.query(sys, DGG.MultiOrderCoverage(BOXB); level=l + 1))
        @test EN.reference_level(deep) == l + 1
        u = union(a, deep)
        @test EN.reference_level(u) == l + 1
        @test leafset(u, l + 1) == union(leafset(a, l + 1), leafset(deep, l + 1))
        @test leafset(intersect(a, deep), l + 1) ==
              intersect(leafset(a, l + 1), leafset(deep, l + 1))
    end

    @testset "a complete sibling family is stored as its parent" begin
        p = DGG.cellindex(DGG.levelgrid(sys, first(DGG.levels(sys)) + 1), 1)
        kids = DGG.MultiOrderVector(sys, collect(DGG.children(sys, p)))
        @test length(kids) > 1
        @test has_complete_family(sys, kids)
        @test collect(union(kids, kids)) == [p]
        @test union(kids, kids) == DGG.MultiOrderVector(sys, [p])
    end

    @testset "the empty container and its complement" begin
        rootgrid = DGG.levelgrid(sys, first(DGG.levels(sys)))
        empty = DGG.MultiOrderVector(sys, DGG.cellindextype(sys)[])
        @test collect(EN.complement(empty)) == collect(DGG.rootcells(sys))
        @test isempty(EN.complement(DGG.MultiOrderVector(sys, collect(DGG.rootcells(sys)))))
        @test DGG.ncells(rootgrid) == length(DGG.rootcells(sys))
    end

    # Mixed-system operands throw, never a silently empty answer.
    other = DGG.MultiOrderVector(DGG.query(sys isa DGG.HEALPixSystem ? DGG.S2System() :
                                           DGG.HEALPixSystem(),
        DGG.MultiOrderCoverage(BOXA); level=3))
    @test_throws ArgumentError union(a, other)
    @test_throws ArgumentError intersect(a, other)
    @test_throws ArgumentError setdiff(a, other)
end

# ---------------------------------------------------------------------------
# A5: no descendant ranges, so no container
# ---------------------------------------------------------------------------

@testset "A5 has no intervals to index" begin
    sys = DGG.A5System()
    @test !DGG.has_sorted_subtrees(sys)

    # Both constructors throw `ArgumentError`, not `descendant_range`'s
    # `MethodError`.
    @test_throws ArgumentError DGG.MultiOrderVector(sys, collect(DGG.rootcells(sys)))
    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=5)
    @test !isempty(set)
    @test_throws ArgumentError DGG.MultiOrderVector(set)
    @test_throws ArgumentError DGG.level_ranges(set, 5)

    # Single-level `CellVector` still works on A5.
    @test DGG.CellVector(set) isa DGG.CellVector
end

end # module MultiOrderVectorTests
