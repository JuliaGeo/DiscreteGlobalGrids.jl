# ---------------------------------------------------------------------------
# T18 — the BUDGET mode of multi-order coverage, on every system.
#
# `query(sys, MultiOrderCoverage(target); level = l)` refines everything the
# boundary crosses to a fixed depth and returns however many cells that takes.
# `query(sys, MultiOrderCoverage(target); maxcells = n)` answers the other
# question — "ten cells that cover California, or a hundred" — by refining the
# crossing cells coarsest first and stopping when the next replacement would not
# fit. Same predicates, same prunes, same result type; a different schedule and
# a different thing held fixed.
#
# This file is the T15 treatment of that second mode: the same committed
# `test/fixtures/california.txt` outline, the same seven systems (six registered
# plus an `AuthalicSystem` wrap), the same dense-plus-adversarial point sampler.
#
# THE LAWS, and what each is worth:
#
#   * THE BUDGET IS A BOUND — `length(set) <= maxcells`, at every budget, with
#     one documented exception: a SEED bigger than the budget is returned whole,
#     because the seed is the coarsest covering the traversal can express and
#     there is nothing to refine away. Both arms are pinned, the exception with
#     a target that produces it.
#   * DETERMINISM — same inputs, same output, cells and flags and keys. The
#     schedule is a priority order over (level, position within the level) and
#     nothing else, so this is a statement about the implementation containing
#     no `Dict` iteration and no ties.
#   * SIZE IS MONOTONE IN THE BUDGET — a measurement of the schedule, not a
#     theorem. Raising the budget cannot make the answer smaller in any case
#     swept here, and the testset says which it swept.
#   * SHAPE — the result is a `MultiOrderCellSet` in every sense the accuracy
#     mode's is: no emitted cell has an emitted ancestor, the reference-level
#     intervals ascend and are disjoint, `cell_polygons` and
#     `coarsest_contained` work, and A5 takes the documented `(level, id)`
#     fallback.
#   * PROVENANCE, AND WHY IT IS SHARPER HERE — the budget traversal asks
#     `Within` of every cell it admits, because contained and crossing is the
#     classification the schedule runs on. So `is_contained` is EXACT in both
#     directions on a budget set, unless refinement actually reached the
#     `maxlevel` cap, where the accuracy mode's blind spot returns verbatim.
#     Both arms are asserted.
#   * COVERING — and this is the one that needs care, because the two readings
#     of "cover" do not agree here the way they do for the accuracy mode.
#
#       - THE UNION READING: every point of the target lies inside some emitted
#         cell. This is what "ten cells that cover California" means, and it is
#         EXACT — zero misses at every budget — on the three systems whose four
#         children tile their parent. Elsewhere it degrades exactly as
#         `MultiOrderCoverage`'s own warning describes, because it is the same
#         phenomenon: replacing a cell by its children swaps one footprint for
#         another. Bounded, not asserted zero.
#       - THE LEAF READING: every cell of the reference level that meets the
#         target is a member or the descendant of one. The accuracy mode has
#         this on all seven — it earns it by descending into cells that MISS the
#         target, because a child can overhang its parent, and carrying that all
#         the way down to `level`. A budget has no fixed depth to carry it to;
#         stopping early is the whole point of the mode, and a branch stopped at
#         a cell that misses the target is a branch whose overhang was not
#         followed. So this too is a law on the congruent three and a bounded
#         measurement on the other four.
#
#     Neither is a defect of the schedule. A bounded look past a child that
#     misses the target was measured while this was written: two extra levels
#     take IGEO7, H3 and the authalic wrap to zero misses on these fixtures and
#     A5 to a handful, three extra levels leave A5 alone still failing. There is
#     no depth at which it becomes a proof — `_subtree_outside` never fires on a
#     cell hugging the target's boundary from outside, at any level — so the
#     constant would be a number that happened to pass, and the mode ships
#     without one. The numbers below are what it does instead.
#
#   * EQUIVALENCE WITH THE ACCURACY MODE — give the budget exactly the number of
#     cells `level = L` produces, and cap it at `L`, and the two modes return
#     the SAME SET, cell for cell and flag for flag. That is the sharpest
#     statement available that the budget mode is the same traversal: it holds
#     on the three congruent systems at every level swept. On the other four the
#     two diverge for the covering reason above, and the testset states that as
#     an exclusion with its reason rather than weakening the law.
#   * COMPOSITION — a budget set backs a `CellLookup` at its own reference level
#     and at any deeper one, the positions round-trip, and the expansion equals
#     the brute-force `descendants` walk. Ten cells chosen for the picture still
#     name every leaf under them for the data.
#
# A5 runs all of it. The budget traversal needs `children` and `cellposition`
# and nothing else — no `descendant_range` — and the testset at the bottom
# asserts that from the other side.
# ---------------------------------------------------------------------------

module MultiOrderBudgetTests

using Test
import DiscreteGlobalGrids as DGG
const FB = DGG.Fallbacks
import GeoInterface as GI
import GeometryOps as GO
import DimensionalData as DD

# ---------------------------------------------------------------------------
# The fixture — the same committed outline the accuracy mode is tested against
# ---------------------------------------------------------------------------

const FIXTURE = joinpath(@__DIR__, "..", "..", "fixtures", "california.txt")

# `P` opens a polygon, `R` opens a ring, everything else is `lon lat`.
function load_parts(path)
    parts = Vector{Vector{Vector{Tuple{Float64,Float64}}}}()
    ring = Tuple{Float64,Float64}[]
    for line in eachline(path)
        (isempty(line) || startswith(line, '#')) && continue
        if line == "P"
            push!(parts, Vector{Vector{Tuple{Float64,Float64}}}())
        elseif line == "R"
            ring = Tuple{Float64,Float64}[]
            push!(parts[end], ring)
        else
            lon, lat = split(line)
            push!(ring, (parse(Float64, lon), parse(Float64, lat)))
        end
    end
    return parts
end

topolygon(rings) = GI.Polygon([GI.LinearRing(r) for r in rings])

const PARTS = load_parts(FIXTURE)
const MAINLAND = topolygon(PARTS[1])
const CALIFORNIA = GI.MultiPolygon([topolygon(p) for p in PARTS])

# A San Francisco city block: smaller than one cell of any system's coarsest
# levels, which is the edge case where a budget buys depth rather than breadth.
const BLOCK = GI.Polygon([GI.LinearRing([(-122.42, 37.77), (-122.40, 37.77),
    (-122.40, 37.79), (-122.42, 37.79), (-122.42, 37.77)])])

# Two polar caps whose union is 66% of the sphere. Every root cell of every
# system meets it, so its SEED is the whole top level and a small budget cannot
# refine any of it — the documented over-budget return.
parallel_ring(lat) = GI.Polygon([GI.LinearRing([(lon, lat) for lon in 0.0:5.0:360.0])])
const WIDE = GI.MultiPolygon([parallel_ring(20.0), parallel_ring(-20.0)])

# Somewhere no cell of anything meets: a degenerate target is not swept here,
# but an empty answer has to stay a well-formed set, so one is built below.

# ---------------------------------------------------------------------------
# Systems
#
# The last column is CONGRUENCE — whether a cell's children exactly tile it —
# for the same reason `multiorder_polygons.jl` carries it: it is not a trait,
# nothing in the interface asks for it, and it decides which arm two of the
# laws take. HEALPix, S2 and ISEA4R are aperture-4 quadtrees on a chart and
# four children tile their parent; IGEO7 and H3 are aperture 7 and their seven
# children are a rotated rosette with the parent's area and not its footprint;
# A5's four Hilbert children do not even stay inside it.
#
# The middle column is the level the equivalence law sweeps up to. It is
# per-system because the apertures differ and the set sizes have to stay small
# enough to run the accuracy mode repeatedly.
# ---------------------------------------------------------------------------

const SWEEP = [
    (DGG.IGeo7System(), 5, false),
    (DGG.H3System(), 4, false),
    (DGG.HEALPixSystem(), 6, true),
    (DGG.A5System(), 6, false),
    (DGG.S2System(), 7, true),
    (DGG.ISEA4RSystem(), 6, true),
    (DGG.AuthalicSystem(DGG.IGeo7System()), 5, false),
]

# The budgets the sampled laws run at, and the wider ladder the counting laws
# use. 10 and 100 are the owner's own two numbers.
const BUDGETS = (10, 40, 100)
const LADDER = (1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377)

sysname(sys) = sys isa DGG.AuthalicSystem ?
               "Authalic($(nameof(typeof(parent(sys)))))" : string(nameof(typeof(sys)))

# ---------------------------------------------------------------------------
# How much of the target a budget covering is allowed to miss, per system and
# per reading.
#
# These are NOT slack. `query`'s docstring ships these very numbers as the
# measured cost of non-congruence, so a uniform bound would let IGEO7 drift to
# A5's figure with the suite still green and the docstring quietly false. Each
# bound is that system's own claim, and the sweep below asserts it.
#
# Calibration: the worst fraction seen over budgets 10/40/100 on this file's
# 0.1-degree sampler, and separately on an independent 7505-point sampler
# during review. Each bound is at least 1.3x the worse of the two, rounded to a
# round number — wide enough that no reshuffle of the sample grid trips it,
# narrow enough that a change in KIND does. The congruent three are exact and
# get no entry: zero is their bound, and it is a theorem rather than a
# measurement.
#
#   system            union: here / review / bound     leaf: here / review / bound
#   IGEO7               1.01% / 1.25% / 2%                 0% / 0.19% / 1%
#   Authalic(IGEO7)     1.36% /    -  / 3%              0.31% /    -  / 2%
#   H3                 11.09% /  7.4% / 15%                 0% /    0% / 2%
#   A5                 22.77% / 16.3% / 30%             11.72% /  8.7% / 18%
#
# H3's union figure is worst at the SMALLEST budget, not the largest: ten cells
# over a state is a handful of aperture-7 rosettes, and one rosette's footprint
# mismatch is a large share of the answer. That is why the bound is not simply
# the large-budget number.
# ---------------------------------------------------------------------------

const UNION_BOUND = Dict("IGeo7System" => 0.02, "Authalic(IGeo7System)" => 0.03,
    "H3System" => 0.15, "A5System" => 0.30)
const LEAF_BOUND = Dict("IGeo7System" => 0.01, "Authalic(IGeo7System)" => 0.02,
    "H3System" => 0.02, "A5System" => 0.18)

@testset "the budget sweep covers every registered system" begin
    swept = Set(typeof(s) for (s, _, _) in SWEEP)
    for s in DGG.systems()
        @test typeof(s) in swept
    end
    @test any(s -> s isa DGG.AuthalicSystem, first.(SWEEP))
end

# ---------------------------------------------------------------------------
# Oracles — the engine's own preparation, asked about POINTS
# ---------------------------------------------------------------------------

prepare(geom) = FB._query_target(geom)
inside(t, lon, lat) = GO.relate_predicate(t.prepared, GO.pred_contains(),
    FB.unit_point(lon, lat))

function interior_samples(geom, t, n::Int)
    ext = GI.extent(geom)
    x0, x1 = ext.X
    y0, y1 = ext.Y
    out = Tuple{Float64,Float64}[]
    for i in 0:n, j in 0:n
        lon = x0 + (x1 - x0) * i / n
        lat = y0 + (y1 - y0) * j / n
        inside(t, lon, lat) && push!(out, (lon, lat))
    end
    return out
end

# A hair inside the outline, where a covering that rounds the wrong way loses
# points. Offsets straddle the boundary and the oracle keeps what lands inside.
function boundary_samples(geom, t; stride::Int=1)
    corners = Tuple{Float64,Float64}[]
    for r in GI.getring(geom), p in GI.getpoint(r)
        push!(corners, (GI.x(p), GI.y(p)))
    end
    out = Tuple{Float64,Float64}[]
    for k in 1:stride:length(corners)
        a = corners[k]
        b = corners[k == length(corners) ? 1 : k+1]
        for base in (a, ((a[1] + b[1]) / 2, (a[2] + b[2]) / 2)),
            e in (1e-3, 1e-4, 1e-5, -1e-4),
            d in ((e, 0.0), (0.0, e), (e, e), (-e, e))

            q = (base[1] + d[1], base[2] + d[2])
            inside(t, q...) && push!(out, q)
        end
    end
    return out
end

const TARGET = prepare(MAINLAND)
const SAMPLES = vcat(interior_samples(MAINLAND, TARGET, 40),
    boundary_samples(MAINLAND, TARGET; stride=5))

# THE LEAF READING: a point is covered when its own reference-level cell is
# emitted, or an ancestor of it is.
function leaf_misses(sys, set, samples)
    leaf = set.reference_level
    grid = DGG.levelgrid(sys, leaf)
    emitted = Set(collect(set))
    n = 0
    for (lon, lat) in samples
        c = DGG.cellat(grid, lon, lat)
        if c === nothing
            n += 1
            continue
        end
        hit = c in emitted
        for l in DGG.level(c)-1:-1:first(DGG.levels(sys))
            hit && break
            hit = DGG.ancestor(sys, c, l) in emitted
        end
        hit || (n += 1)
    end
    return n
end

# THE UNION READING: a point is covered when it lies in some emitted cell's
# polygon. One cap per cell so the scan is a distance comparison per pair and a
# point-in-polygon only where the cap admits it.
function union_misses(set, samples)
    rings = [FB.open_ring(DGG.cell_boundary(set, c)) for c in set]
    caps = [FB.points_cap(view(r[1], 1:r[2])) for r in rings]
    n = 0
    for (lon, lat) in samples
        p = FB.unit_point(lon, lat)
        hit = false
        for k in eachindex(rings)
            FB.cap_contains(caps[k], p) || continue
            if FB.point_in_cell(view(rings[k][1], 1:rings[k][2]), p) !== false
                hit = true
                break
            end
        end
        hit || (n += 1)
    end
    return n
end

function emitted_ancestors(sys, set)
    emitted = Set(collect(set))
    return [(c, l) for c in set for l in first(DGG.levels(sys)):DGG.level(c)-1
            if DGG.ancestor(sys, c, l) in emitted]
end

expand(sys, set, l) = DGG.has_sorted_subtrees(sys) ? DGG.cellindices(set, l) :
                      sort!(reduce(vcat, [DGG.descendants(sys, c, l) for c in set]))

# ---------------------------------------------------------------------------
# The keyword rules — one system is enough, they are not per-system
# ---------------------------------------------------------------------------

@testset "the two modes are modes, not a bound and a hint" begin
    sys = DGG.S2System()
    cov = DGG.MultiOrderCoverage(MAINLAND)

    @test_throws ArgumentError DGG.query(sys, cov)
    @test_throws ArgumentError DGG.query(sys, cov; level=5, maxcells=10)
    @test_throws ArgumentError DGG.query(sys, cov; level=5, maxlevel=7)
    @test_throws ArgumentError DGG.query(sys, cov; maxcells=0)
    @test_throws ArgumentError DGG.query(sys, cov; maxcells=-3)
    @test_throws ArgumentError DGG.query(sys, cov; maxcells=10, maxlevel=99)
    @test_throws ArgumentError DGG.query(sys, cov; maxcells=10, maxlevel=-1)
    # `level` mode is untouched by any of it.
    @test length(DGG.query(sys, cov; level=5)) > 0
    # The type constructor is the same two modes and the same errors.
    @test FB.MultiOrderCellSet(sys, cov; maxcells=12) isa DGG.MultiOrderCellSet
    @test_throws ArgumentError FB.MultiOrderCellSet(sys, cov)
    @test_throws ArgumentError FB.MultiOrderCellSet(sys, cov; level=5, maxcells=10)
end

# ---------------------------------------------------------------------------
# The laws, per system
# ---------------------------------------------------------------------------

@testset "a budget covering of a real outline: $(sysname(sys))" for
    (sys, toplevel, congruent) in SWEEP

    sets = Dict(b => DGG.query(sys, DGG.MultiOrderCoverage(MAINLAND); maxcells=b)
                for b in BUDGETS)

    @testset "shape" begin
        for b in BUDGETS
            set = sets[b]
            @test set isa DGG.MultiOrderCellSet
            @test DGG.system(set) === sys
            @test eltype(set) === DGG.cellindextype(sys)
            @test !isempty(set)
            cells = collect(set)
            @test allunique(cells)
            @test isempty(emitted_ancestors(sys, set))
            # The reference level is the deepest level the budget reached, not a
            # number the caller supplied — that is what makes it a *budget*.
            @test set.reference_level == maximum(DGG.level, cells)
            @test length(DGG.cell_polygons(set)) == length(set)
            if DGG.has_sorted_subtrees(sys)
                intervals = [DGG.descendant_range(sys, c, set.reference_level) for c in cells]
                @test issorted(intervals; by=first)
                @test all(k -> first(intervals[k]) > last(intervals[k-1]), 2:length(intervals))
                @test FB.curve_keys(set) == first.(intervals)
            else
                # EXCLUDED with its reason: no descendant ranges on A5, so the
                # documented `(level, id)` fallback is what is asserted.
                @test !DGG.has_sorted_subtrees(sys)
                @test issorted(cells; by=c -> (DGG.level(c), c))
            end
        end
    end

    @testset "the budget is a bound, and size climbs with it" begin
        sizes = Int[]
        for b in LADDER
            set = DGG.query(sys, DGG.MultiOrderCoverage(MAINLAND); maxcells=b)
            push!(sizes, length(set))
            # The seed on this outline is one or two root cells, so the
            # exception below cannot fire here; the over-budget case has its own
            # testset with a target that produces it.
            @test length(set) <= max(b, 2)
        end
        @test issorted(sizes)                 # measurement of the schedule
        @test sizes[end] == LADDER[end]       # and the budget is actually spent
    end

    @testset "determinism" begin
        for b in (10, 137)
            a = DGG.query(sys, DGG.MultiOrderCoverage(MAINLAND); maxcells=b)
            c = DGG.query(sys, DGG.MultiOrderCoverage(MAINLAND); maxcells=b)
            @test collect(a) == collect(c)
            @test a.contained == c.contained
            @test FB.curve_keys(a) == FB.curve_keys(c)
            @test a.reference_level == c.reference_level
        end
    end

    @testset "provenance is exact, because the schedule runs on it" begin
        set = sets[100]
        # Refinement never reached the system's deepest level, so `Within` was
        # asked of every admitted cell and the flag is exact BOTH ways — which
        # the accuracy mode's is not, at its own maximum depth.
        @test set.reference_level < last(DGG.levels(sys))
        @test all(eachindex(set)) do i
            c = set[i]
            DGG.is_contained(set, i) == FB._matches(DGG.Within(nothing), TARGET,
                DGG.levelgrid(sys, DGG.level(c)), c)
        end
        best = DGG.coarsest_contained(set)
        @test best !== nothing
        @test FB._matches(DGG.Within(nothing), TARGET,
            DGG.levelgrid(sys, DGG.level(best)), best)
        @test DGG.level(best) ==
              minimum(DGG.level(set[i]) for i in eachindex(set) if DGG.is_contained(set, i))
    end

    @testset "a maxlevel the refinement reaches restores the blind spot" begin
        # Cap the descent shallow enough that the budget cannot spend itself,
        # and the cells stranded at the cap are flagged unproven without being
        # asked — the accuracy mode's contract, verbatim.
        #
        # The cap is `first + 6` and not something shallower because the second
        # assertion needs cells that ACTUALLY fit inside California. At a cap of
        # `first + 2` no cell of any of the seven does — level-2 cells are
        # continent-sized — so the flag being `false` there would be a true
        # statement rather than a blind one, and would prove nothing.
        cap = first(DGG.levels(sys)) + 6
        set = DGG.query(sys, DGG.MultiOrderCoverage(MAINLAND); maxcells=10_000, maxlevel=cap)
        @test set.reference_level == cap
        stranded = [i for i in eachindex(set) if DGG.level(set[i]) == cap]
        @test !isempty(stranded)
        @test all(i -> !DGG.is_contained(set, i), stranded)
        # And it IS a blind spot rather than a negative fact: ask `Within` of
        # the same cells and some of them say yes. The flag does not, because
        # nobody asked it.
        grid = DGG.levelgrid(sys, cap)
        fits = count(i -> FB._matches(DGG.Within(nothing), TARGET, grid, set[i]), stranded)
        @test fits > 0
    end

    @testset "covering: the union reading" begin
        @test length(SAMPLES) > 2000
        # Zero for the congruent three, and this system's OWN docstring figure
        # for the rest — see `UNION_BOUND` for why it is not one number.
        bound = congruent ? 0.0 : UNION_BOUND[sysname(sys)]
        for b in BUDGETS
            missed = union_misses(sets[b], SAMPLES)
            if congruent
                # Four children tile their parent, so replacing a cell by the
                # children that meet the target loses nothing.
                @test missed == 0
            else
                @test missed / length(SAMPLES) < bound
            end
        end
    end

    @testset "covering: the leaf reading" begin
        bound = congruent ? 0.0 : LEAF_BOUND[sysname(sys)]
        for b in BUDGETS
            missed = leaf_misses(sys, sets[b], SAMPLES)
            if congruent
                @test missed == 0
            else
                # Same cause, measured through ancestry rather than geometry.
                @test missed / length(SAMPLES) < bound
            end
        end
    end

    @testset "a budget the size of a level reproduces that level" begin
        for l in (first(DGG.levels(sys))+1):toplevel
            accurate = DGG.query(sys, DGG.MultiOrderCoverage(MAINLAND); level=l)
            budgeted = DGG.query(sys, DGG.MultiOrderCoverage(MAINLAND);
                maxcells=length(accurate), maxlevel=l)
            if congruent
                # The sharpest available statement that this is one traversal
                # with two schedules: same cells, same order, same flags.
                @test collect(budgeted) == collect(accurate)
                @test budgeted.contained == accurate.contained
                @test budgeted.reference_level == accurate.reference_level
            else
                # EXCLUDED, with its reason: the accuracy mode descends into
                # cells that miss the target because a child can overhang its
                # parent, and emits what it finds beneath them. A budget cannot,
                # so it is a subset of the same size or smaller. Asserting
                # equality here would be asserting congruence.
                @test length(budgeted) <= length(accurate)
                @test !isempty(budgeted)
            end
        end
    end

    @testset "composition: a budget set backs a cell axis" begin
        set = sets[40]
        leaf = set.reference_level + 2
        lk = DGG.CellLookup(set; level=leaf)
        @test DGG.level(lk) == leaf
        @test DGG.system(lk) === sys
        brute = expand(sys, set, leaf)
        @test collect(lk) == brute
        @test length(lk) == length(brute)
        for k in (1, max(1, length(lk) ÷ 3), length(lk))
            @test DGG.cellposition(lk, lk[k]) == k
        end
        # Its own reference level is the default, and it is a lookup too.
        @test length(DGG.CellLookup(set)) == length(expand(sys, set, set.reference_level))
        A = DD.DimArray(Float64.(eachindex(lk)), DGG.Cells(lk))
        @test A[DGG.Cells(DD.At(lk[3]))] == 3.0
        @test DGG.cellset(lk) === set
    end
end

# ---------------------------------------------------------------------------
# The documented edges
# ---------------------------------------------------------------------------

@testset "a seed bigger than the budget comes back whole: $(sysname(sys))" for
    (sys, _, _) in SWEEP

    # Every root cell of every system meets a target covering 66% of the sphere,
    # so the seed IS the top level and no refinement of it can fit in three
    # cells. Returning the seed over budget is the documented answer: it is the
    # coarsest covering this traversal can express, and there is nothing to
    # refine away.
    set = DGG.query(sys, DGG.MultiOrderCoverage(WIDE); maxcells=3)
    roots = collect(DGG.rootcells(sys))
    @test length(set) > 3
    @test length(set) <= length(roots)
    @test all(c -> DGG.level(c) == first(DGG.levels(sys)), set)
    @test set.reference_level == first(DGG.levels(sys))
    @test isempty(emitted_ancestors(sys, set))
    # And it does cover: at the top level the two readings coincide.
    samples = [(lon, lat) for lon in -170.0:20.0:170.0, lat in -85.0:5.0:85.0]
    t = prepare(WIDE)
    samples = [q for q in vec(samples) if inside(t, q...)]
    @test length(samples) > 200
    @test leaf_misses(sys, set, samples) == 0
end

@testset "a target smaller than one cell: $(sysname(sys))" for (sys, _, _) in SWEEP
    # A budget of one cannot split anything, so refinement can only follow the
    # single crossing child down — the chain stops at the first cell the target
    # crosses into two of.
    one = DGG.query(sys, DGG.MultiOrderCoverage(BLOCK); maxcells=1)
    @test length(one) == 1
    @test !DGG.is_contained(one, 1)
    @test DGG.coarsest_contained(one) === nothing
    @test one.reference_level == DGG.level(one[1])
    t = prepare(BLOCK)
    @test FB._matches(DGG.Intersects(nothing), t,
        DGG.levelgrid(sys, DGG.level(one[1])), one[1])
    # More budget can only go deeper or wider, never fewer.
    more = DGG.query(sys, DGG.MultiOrderCoverage(BLOCK); maxcells=8)
    @test length(more) >= length(one)
    @test length(more) <= 8
end

@testset "an empty answer is still a set" begin
    # A target the traversal admits nothing for: a cap of zero radius over a
    # point is met by the cell containing it, so instead the emptiness is built
    # from a level cap the seed cannot pass — a `maxlevel` at the top with a
    # target no root meets is not reachable, so this asserts the shape of the
    # zero-cell path directly through the constructor instead.
    sys = DGG.S2System()
    empty = FB._sorted_cell_set(sys, DGG.cellindextype(sys)[], falses(0),
        first(DGG.levels(sys)))
    @test isempty(empty)
    @test length(empty) == 0
    @test DGG.coarsest_contained(empty) === nothing
    @test isempty(DGG.cell_polygons(empty))
end

# ---------------------------------------------------------------------------
# A5, from the other side
# ---------------------------------------------------------------------------

@testset "the budget traversal needs children, not descendant ranges" begin
    sys = DGG.A5System()
    @test !DGG.has_sorted_subtrees(sys)
    set = DGG.query(sys, DGG.MultiOrderCoverage(MAINLAND); maxcells=100)
    @test length(set) == 100
    # The schedule's tie-break is `cellposition` within a level, which every
    # system answers; `descendant_range` has no method here at all, and the
    # traversal never reaches for one.
    @test_throws MethodError DGG.descendant_range(sys, first(collect(set)),
        set.reference_level)
    @test_throws ArgumentError DGG.level_ranges(set, set.reference_level)
    # Everything downstream of the set still works, through the selection path.
    @test issorted(collect(set); by=c -> (DGG.level(c), c))
    lk = DGG.CellLookup(set)
    @test length(lk) == length(expand(sys, set, set.reference_level))
end

@testset "multipolygon and multipart budgets" begin
    # The offshore islands are separate components, and a budget spent on the
    # mainland must not lose them: every part has to meet some emitted cell.
    for (sys, _, _) in SWEEP
        set = DGG.query(sys, DGG.MultiOrderCoverage(CALIFORNIA); maxcells=200)
        @test length(set) <= 200
        @test isempty(emitted_ancestors(sys, set))
        emitted = Set(collect(set))
        for (k, part) in enumerate(PARTS)
            poly = topolygon(part)
            t = prepare(poly)
            # The seed cell of each part is coarse enough that SOME emitted cell
            # meets it, whatever depth the budget reached elsewhere.
            hit = any(set) do c
                FB._matches(DGG.Intersects(nothing), t,
                    DGG.levelgrid(sys, DGG.level(c)), c)
            end
            @test hit || error("part $k of $(sysname(sys)) met no emitted cell")
        end
    end
end

end # module MultiOrderBudgetTests
