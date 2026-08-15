# ---------------------------------------------------------------------------
# MOC storage, wave 1: `MultiOrderVector`, the mixed-level cell container.
#
# `MultiOrderCellSet` is a coverage; this is storage — the axis an adaptively
# refined mesh's values hang off. What has to be true of it is not what has to
# be true of a query result, so this file states the container's own laws:
#
#   * ROUND TRIP — a coverage read as storage keeps its cells, in its order, at
#     its reference level, and expands to exactly the `CellVector` the set
#     itself expands to, at every level.
#   * TWO KINDS OF MEMBERSHIP — `cellposition` is EXACT (a descendant of a
#     stored cell is not stored) and `covering_position` is the compression
#     verb (every leaf under a stored cell resolves to it). Confusing the two
#     is the mistake the type exists to make impossible, so both directions are
#     asserted against each other.
#   * POINT LOCATION — a cell's own centroid comes back as that cell, whatever
#     level it sits at. That is `cellat`'s contract read through a mixture of
#     levels, and it is the only verb here that touches geometry.
#   * SET ALGEBRA — union, intersect, setdiff and complement against a brute
#     force over leaf POSITION SETS at a shallow level. The oracle is
#     `Base.union` on a `Set{Int}`, which shares no code with the interval
#     arithmetic under test.
#   * NORMALIZATION — the results are the COARSEST cells that tile them. Stated
#     twice, because it is the one claim a correct-looking implementation
#     silently drops: no complete sibling family survives in any result, a
#     family unioned with itself comes back as its parent, and the whole sphere
#     comes back as the root cells.
#   * VALIDATION — an ancestor/descendant pair, a duplicate and a reference
#     level shallower than a stored cell are rejected; an unsorted vector is
#     sorted rather than refused.
#
# `has_sorted_subtrees(A5System())` is `false`, so A5 has no descendant ranges,
# no intervals, and no container at all. It is excluded from the sweep with
# that reason and asserted to throw, the way `level_ranges` is.
# ---------------------------------------------------------------------------

module MultiOrderVectorTests

using Test
import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import GeometryOps as GO

const FB = DGG.Fallbacks

# The Switzerland/Zurich fixture the `CellVector` suite uses, so the two files
# can be read against each other.
const REGION = GI.Polygon([GI.LinearRing([(6.0, 45.8), (10.5, 45.8), (10.5, 47.8),
    (6.0, 47.8), (6.0, 45.8)])])
const ZURICH = GI.Polygon([GI.LinearRing([(8.3, 47.2), (8.8, 47.2), (8.8, 47.6),
    (8.3, 47.6), (8.3, 47.2)])])
const NOWHERE = GI.Polygon([GI.LinearRing([(-26.0, -41.0), (-24.0, -41.0),
    (-24.0, -39.0), (-26.0, -41.0)])])
const FARAWAY = (-25.0, -40.0)

# Two overlapping continent-sized boxes, for the set algebra: each has cells the
# other does not, and they share a quadrant. Big, because the algebra is checked
# at a level shallow enough to expand the whole sphere.
const BOXA = GI.Polygon([GI.LinearRing([(-20.0, 20.0), (40.0, 20.0), (40.0, 60.0),
    (-20.0, 60.0), (-20.0, 20.0)])])
const BOXB = GI.Polygon([GI.LinearRing([(10.0, 0.0), (70.0, 0.0), (70.0, 40.0),
    (10.0, 40.0), (10.0, 0.0)])])

const LONLAT = GO.UnitSpherical.GeographicFromUnitSphere()

# The leaf level per system — the `CellVector` suite's levels, chosen so that a
# leaf is a few kilometres across on all of them. The apertures are 7 and 4, so
# a fixed level would be two very different depths, and a coverage that never
# fits a whole cell inside REGION comes back flat, which would make every law
# below a single-level law wearing another name. A5 is absent by construction
# and the testset at the bottom says why.
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

# The shallowest level with enough cells to make a set-algebra oracle
# interesting. Computed rather than tabulated: the apertures are 4 and 7, so a
# fixed level would mean six different problem sizes.
function shallow_level(sys; atleast=400)
    for l in DGG.levels(sys)
        DGG.ncells(DGG.levelgrid(sys, l)) >= atleast && return l
    end
    return last(DGG.levels(sys))
end

# The oracle's currency: which leaf positions a container names. Built from
# `descendant_range` directly, never from the container's own index.
leafset(mov, l) = Set(p for c in mov for p in DGG.descendant_range(DGG.system(mov), c, l))

# The two halves of "sorted and pairwise disjoint", read off the intervals.
function disjoint_and_sorted(mov)
    ivs = [DGG.descendant_range(DGG.system(mov), c, FB.reference_level(mov)) for c in mov]
    return issorted(ivs; by=first) &&
           all(k -> first(ivs[k]) > last(ivs[k-1]), 2:length(ivs))
end

# The normalization law, from the other side: a complete sibling family means a
# coarser cell would have said the same thing in fewer words.
function has_complete_family(sys, mov)
    stored = Set(collect(mov))
    top = first(DGG.levels(sys))
    for c in mov
        DGG.level(c) > top || continue
        all(k -> k in stored, DGG.children(sys, parent(sys, c))) && return true
    end
    return false
end

# Deterministic probes: the ends, the middles, and one cell at each extreme of
# the level mixture — the levels are what makes this container different from a
# `CellVector`, so both extremes are always sampled.
function probes(mov)
    n = length(mov)
    return unique([1, cld(n, 3), cld(n, 2), n,
        argmin(k -> DGG.level(mov[k]), eachindex(mov)),
        argmax(k -> DGG.level(mov[k]), eachindex(mov))])
end

@testset "the sweep covers every registered system" begin
    swept = Set(typeof(s) for (s, _) in SWEEP)
    for s in DGG.systems()
        # The container IS interval arithmetic, so the sweep is exactly the
        # systems that have intervals. A system gaining or losing the trait
        # without this file noticing fails here.
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
        @test FB.reference_level(mov) == leaf
        @test length(mov) == length(set)
        @test collect(mov) == cells
        @test mov[:] === mov
        @test mov == mov
        # Mixed levels are the whole premise; a coverage that came back flat
        # would make every law below a `CellVector` law wearing another name.
        @test length(unique(DGG.level.(cells))) > 1
        @test maximum(DGG.level, cells) == leaf
        @test disjoint_and_sorted(mov)
        # `offsets` is the cumulative leaf count at the reference level — the
        # thing wave 2's lazy expansion binary-searches — so it has to agree
        # with the expansion itself.
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
            # Not merely the same cells: the same WINDOWS, so a container that
            # forgot to merge adjacent subtrees is caught here rather than
            # surviving as a slower equal answer.
            @test FB.windows(bridged).starts == FB.windows(direct).starts
            @test FB.windows(bridged).stops == FB.windows(direct).stops
        end
        @test DGG.cellset(DGG.CellVector(mov)) === mov
        # A level above the deepest stored cell would need an ancestor, which
        # covers more than the container does.
        @test_throws ArgumentError DGG.CellVector(mov; level=leaf - 1)
    end

    @testset "membership is exact" begin
        @test all(DGG.cellposition(mov, mov[k]) == k for k in eachindex(mov))
        @test all(c in mov for c in mov)

        # A descendant of a stored cell is covered but NOT stored. This is the
        # line between the two verbs, and the mutant it kills — `cellposition`
        # answering with the covering position — is the plausible one.
        coarse = argmin(k -> DGG.level(mov[k]), eachindex(mov))
        kid = first(DGG.children(sys, mov[coarse]))
        @test DGG.level(kid) <= leaf
        @test DGG.cellposition(mov, kid) === nothing
        @test !(kid in mov)

        # So is an ancestor of one, from the other direction.
        deepest = argmax(k -> DGG.level(mov[k]), eachindex(mov))
        anc = DGG.ancestor(sys, mov[deepest], DGG.level(mov[deepest]) - 1)
        @test DGG.cellposition(mov, anc) === nothing
        @test FB.covering_position(mov, anc) === nothing

        # And a cell of the leaf grid the container does not reach at all.
        outside = DGG.cellat(grid, FARAWAY...)
        @test outside !== nothing
        @test DGG.cellposition(mov, outside) === nothing
        @test FB.covering_position(mov, outside) === nothing
    end

    @testset "the covering ancestor" begin
        # Every leaf position under a stored cell resolves to that cell — the
        # compression law, sampled at both ends and the middle of each probe's
        # range rather than swept over millions.
        for i in ks
            r = DGG.descendant_range(sys, mov[i], leaf)
            for p in (first(r), (first(r) + last(r)) ÷ 2, last(r))
                @test FB.covering_position(mov, DGG.cellindex(grid, p)) == i
            end
        end
        # A cell DEEPER than the reference level is keyed through its ancestor
        # there, so the container answers about leaves it was never keyed for.
        deep = first(DGG.children(sys, DGG.cellindex(grid, mov.stops[end])))
        @test DGG.level(deep) == leaf + 1
        @test FB.covering_position(mov, deep) == length(mov)
        @test DGG.cellposition(mov, deep) === nothing
    end

    @testset "a point lands in the cell holding it" begin
        # `cell_centroid` is strictly interior to its cell by contract, so a
        # stored cell's own centroid must come back as that cell — at whatever
        # level it sits, which is the mixture this type exists for.
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
            @test FB.reference_level(sub) == leaf
            # Kept whole, never clipped: the selection is a sub-vector of the
            # container's own cells, not a re-cut of the region.
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
        # A mask names positions by index, so a short one is a bounds error
        # rather than a silently shorter answer.
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
        @test FB.reference_level(deeper) == leaf + 1
        @test deeper.starts != mov.starts        # ... different integers ...
        @test deeper == mov && mov == deeper     # ... same container
        @test collect(deeper) == cells
        @test DGG.CellVector(deeper) == DGG.CellVector(set; level=leaf + 1)
    end

    @testset "construction from a loose vector" begin
        @test DGG.MultiOrderVector(sys, cells) == mov
        @test DGG.MultiOrderVector(mov) === mov
        # Unsorted input is SORTED, not refused: order is derived from the
        # cells, so there is nothing for a caller to get wrong.
        @test collect(DGG.MultiOrderVector(sys, reverse(cells))) == cells
        @test collect(DGG.MultiOrderVector(sys, vcat(cells[2:2:end], cells[1:2:end]))) == cells

        # Overlapping subtrees are refused, in both spellings of overlap.
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
        @test FB.reference_level(empty) == top
        @test DGG.cellposition(empty, cells[1]) === nothing
        @test FB.covering_position(empty, cells[1]) === nothing
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

    # The fixture has to be a real overlap or three of the four laws are vacuous.
    @test !isempty(intersect(A, B))
    @test !isempty(setdiff(A, B))
    @test !isempty(setdiff(B, A))

    @testset "$name" for (name, got, want) in (
        ("union", union(a, b), union(A, B)),
        ("intersect", intersect(a, b), intersect(A, B)),
        ("setdiff", setdiff(a, b), setdiff(A, B)),
        ("complement", FB.complement(a), setdiff(Set(1:n), A)))

        @test got isa DGG.MultiOrderVector
        @test DGG.system(got) === sys
        @test FB.reference_level(got) == l
        @test leafset(got, l) == want
        @test disjoint_and_sorted(got)
        # Minimal: a complete sibling family would mean the parent could have
        # stood for it, which is the whole content of "coarsest decomposition".
        @test !has_complete_family(sys, got)
    end

    c = FB.complement(a)
    @test isempty(intersect(a, c))
    @test leafset(union(a, c), l) == Set(1:n)
    # The whole sphere IS the root cells — the decomposition's strongest single
    # statement, and the one a climb that stops early cannot make.
    @test collect(union(a, c)) == collect(DGG.rootcells(sys))
    # Applying it twice NORMALIZES rather than round-trips: the same leaves, as
    # the coarsest cells that name them.
    @test leafset(FB.complement(c), l) == A
    @test FB.complement(c) == union(a, a)
    @test !has_complete_family(sys, FB.complement(c))
    # Whether that is the identity is a fact about the OPERAND, and both
    # regimes occur here: a coverage under congruent refinement is already
    # minimal, while an aperture-7 one emits a complete sibling family whenever
    # the parent itself was not inside the target — so normalizing it is
    # strictly coarser. The leaves are the same either way, which is what the
    # docstring promises and the line above already asserted.
    if has_complete_family(sys, a)
        @test length(FB.complement(c)) < length(a)
    else
        @test FB.complement(c) == a
    end

    # The n-ary forms fold the binary one. Base's own would build a `Set` and
    # lose both the order and the normalization.
    @test union(a, b, a) == union(a, b)
    @test intersect(a, b, a) == intersect(a, b)
    @test setdiff(a, b, b) == setdiff(a, b)

    @testset "operands are re-keyed to the deeper reference level" begin
        deep = DGG.MultiOrderVector(DGG.query(sys, DGG.MultiOrderCoverage(BOXB); level=l + 1))
        @test FB.reference_level(deep) == l + 1
        u = union(a, deep)
        @test FB.reference_level(u) == l + 1
        # Without the re-key the two operands' intervals would be compared in
        # different units and the answer would be arbitrary.
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
        @test collect(FB.complement(empty)) == collect(DGG.rootcells(sys))
        @test isempty(FB.complement(DGG.MultiOrderVector(sys, collect(DGG.rootcells(sys)))))
        @test DGG.ncells(rootgrid) == length(DGG.rootcells(sys))
    end

    # Two containers that do not live in the same space are an error, not a
    # silently empty answer.
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

    # Both ways in refuse, with the reason `level_ranges` gives rather than the
    # `MethodError` that `descendant_range` would raise on its own.
    @test_throws ArgumentError DGG.MultiOrderVector(sys, collect(DGG.rootcells(sys)))
    set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=5)
    @test !isempty(set)
    @test_throws ArgumentError DGG.MultiOrderVector(set)
    @test_throws ArgumentError DGG.level_ranges(set, 5)

    # The single-level compression is still available there, and says so: what
    # A5 lacks is the INTERVAL, not the coverage.
    @test DGG.CellVector(set) isa DGG.CellVector
end

end # module MultiOrderVectorTests
