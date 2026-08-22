# Cross-system laws for depth-limited multi-order polygon coverage. Oracles use
# expanded target-level cells rather than the compressed representation.

module MultiOrderPolygonTests

using Test
import DiscreteGlobalGrids as DGG
const FB = DGG.Fallbacks
const EN = DGG.Engine
import GeoInterface as GI
import GeometryOps as GO

include(joinpath(@__DIR__, "..", "..", "helpers.jl"))
using .DGGTestHelpers: syslabel, isquadface, sweepcovers

# ---------------------------------------------------------------------------
# The fixture
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

# One box in the Central Valley with a square hole in it, and one box straddling
# the antimeridian south of Fiji. Both synthetic, both the shape of a real
# defect rather than of a real place.
const DONUT = GI.Polygon([
    GI.LinearRing([(-122.0, 35.0), (-118.0, 35.0), (-118.0, 39.0), (-122.0, 39.0),
        (-122.0, 35.0)]),
    GI.LinearRing([(-121.0, 36.0), (-121.0, 38.0), (-119.0, 38.0), (-119.0, 36.0),
        (-121.0, 36.0)])])
const HOLE = GI.Polygon([GI.LinearRing([(-121.0, 36.0), (-119.0, 36.0), (-119.0, 38.0),
    (-121.0, 38.0), (-121.0, 36.0)])])
const SEAM = GI.Polygon([GI.LinearRing([(176.0, -19.0), (-178.0, -19.0), (-178.0, -16.0),
    (176.0, -16.0), (176.0, -19.0)])])

parallel_ring(lat) = GI.Polygon([GI.LinearRing([(lon, lat) for lon in 0.0:5.0:360.0])])
const WIDE = GI.MultiPolygon([parallel_ring(20.0), parallel_ring(-20.0)])


# (system, leaf level, coarse level). Congruence — whether a cell's children
# tile it exactly, which two of the laws below branch on — is `isquadface`.
const SWEEP = [
    (DGG.IGeo7System(), 7, 5),
    (DGG.H3System(), 6, 4),
    (DGG.HEALPixSystem(), 10, 7),
    (DGG.A5System(), 10, 7),
    (DGG.S2System(), 10, 7),
    (DGG.ISEA4RSystem(), 10, 7),
    (DGG.AuthalicSystem(DGG.IGeo7System()), 7, 5),
]

# Every registered system is swept, so a system added to `systems()` without
# being added here fails this rather than being silently untested.
@testset "the sweep covers every registered system" begin
    sweepcovers(SWEEP)
end


prepare(geom) = EN._query_target(geom)
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

# The adversarial half: points a hair inside the outline, where a coverage that
# rounds the boundary the wrong way loses them. Offsets straddle the boundary in
# both directions and the oracle keeps whichever land inside.
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

function uncovered(sys, set, leaf, samples)
    grid = DGG.levelgrid(sys, leaf)
    emitted = Set(collect(set))
    misses = Tuple{Float64,Float64}[]
    for (lon, lat) in samples
        c = DGG.cellat(grid, lon, lat)
        if c === nothing
            push!(misses, (lon, lat))
            continue
        end
        hit = c in emitted
        for l in DGG.level(c)-1:-1:first(DGG.levels(sys))
            hit && break
            hit = DGG.ancestor(sys, c, l) in emitted
        end
        hit || push!(misses, (lon, lat))
    end
    return misses
end

function emitted_ancestors(sys, set)
    emitted = Set(collect(set))
    return [(c, l) for c in set for l in first(DGG.levels(sys)):DGG.level(c)-1
            if DGG.ancestor(sys, c, l) in emitted]
end

# The parents of every complete sibling family emitted below the maximum depth.
function complete_families(sys, set, leaf)
    emitted = Set(collect(set))
    seen = Set{eltype(set)}()
    out = eltype(set)[]
    for c in set
        (DGG.level(c) <= first(DGG.levels(sys)) || DGG.level(c) >= leaf) && continue
        p = parent(sys, c)
        p in seen && continue
        push!(seen, p)
        all(ch -> ch in emitted, DGG.children(sys, p)) && push!(out, p)
    end
    return out
end

# The set expanded to `l`, whichever way the system can be asked.
expand(sys, set, l) = DGG.has_sorted_subtrees(sys) ? DGG.cellindices(set, l) :
                      sort!(reduce(vcat, [DGG.descendants(sys, c, l) for c in set]))

# ---------------------------------------------------------------------------
# The laws, per system
# ---------------------------------------------------------------------------

@testset "coverage of a real outline: $(syslabel(sys))" for (sys, leaf, coarse) in SWEEP
    congruent = isquadface(sys)
    target = prepare(MAINLAND)
    set = DGG.query(sys, DGG.MultiOrderCoverage(MAINLAND); level=leaf)

    @test set isa DGG.MultiOrderCellSet
    @test DGG.system(set) === sys
    @test eltype(set) === DGG.cellindextype(sys)
    @test !isempty(set)
    cells = collect(set)
    @test allunique(cells)
    @test minimum(DGG.level, cells) < leaf
    @test maximum(DGG.level, cells) == leaf

    @testset "covering" begin
        samples = vcat(interior_samples(MAINLAND, target, 48),
            boundary_samples(MAINLAND, target; stride=3))
        @test length(samples) > 2000          # ... and the sampler is not trivial
        @test isempty(uncovered(sys, set, leaf, samples))
    end

    @testset "no emitted cell has an emitted ancestor" begin
        @test isempty(emitted_ancestors(sys, set))
    end

    @testset "curve order" begin
        if DGG.has_sorted_subtrees(sys)
            intervals = [DGG.descendant_range(sys, c, leaf) for c in cells]
            @test issorted(intervals; by=first)
            @test all(k -> first(intervals[k]) > last(intervals[k-1]), 2:length(intervals))
            @test EN.curve_keys(set) == first.(intervals)
        else
            # EXCLUDED, with its reason: no descendant ranges, so no curve
            # intervals to order by. The documented fallback is `(level, id)`,
            # and that much is still asserted.
            @test !DGG.has_sorted_subtrees(sys)
            @test issorted(cells; by=c -> (DGG.level(c), c))
            @test_throws MethodError DGG.descendant_range(sys, first(cells), leaf)
        end
    end

    @testset "sibling compaction" begin
        # A complete family below the maximum depth means the traversal split a
        # cell into children it then emitted whole — which is right only if the
        # parent itself was not inside the target. Under a congruent refinement
        # that never happens; under aperture 7, and under A5's Hilbert
        # refinement, it does, because the children cover the parent's area
        # without covering its footprint.
        grid_of(c) = DGG.levelgrid(sys, DGG.level(c))
        families = complete_families(sys, set, leaf)
        if congruent
            @test isempty(families)
        end
        for p in families
            @test !EN._matches(DGG.Within(nothing), target, grid_of(p), p)
        end
    end

    @testset "contained vs crossed" begin
        @test length(EN.curve_keys(set)) == length(set)
        contained = [i for i in eachindex(set) if DGG.iscontained(set, i)]
        @test !isempty(contained)
        @test all(i -> DGG.iscontained(set, i) == (DGG.level(set[i]) < leaf), eachindex(set))
        # And where the flag is exact, it agrees with the predicate it stands for.
        for i in Iterators.take(contained, 12)
            c = set[i]
            @test EN._matches(DGG.Within(nothing), target, DGG.levelgrid(sys, DGG.level(c)), c)
        end
        best = DGG.coarsest_contained(set)
        @test best !== nothing
        @test EN._matches(DGG.Within(nothing), target,
            DGG.levelgrid(sys, DGG.level(best)), best)
        @test DGG.level(best) == minimum(DGG.level(set[i]) for i in contained)
    end

    @testset "geometry of a mixed-level set" begin
        polys = DGG.cell_polygons(set)
        @test length(polys) == length(set)
        for i in (1, length(set) ÷ 2, length(set))
            c = set[i]
            g = DGG.levelgrid(sys, DGG.level(c))
            @test DGG.cell_boundary(set, c) == DGG.cell_boundary(g, c)
            @test DGG.cell_centroid(set, c) == DGG.cell_centroid(g, c)
            @test GI.coordinates(DGG.cell_polygon(set, c)) == GI.coordinates(polys[i])
        end
    end

    @testset "expansion against the single-level query" begin
        coarse_set = DGG.query(sys, DGG.MultiOrderCoverage(MAINLAND); level=coarse)
        brute = DGG.query(DGG.levelgrid(sys, coarse), DGG.Intersects(MAINLAND))
        @test !isempty(brute)
        expanded = expand(sys, coarse_set, coarse)
        # The law: nothing that meets the target is lost. Equality is the
        # stronger statement, and on this outline it holds everywhere — a
        # member's descendants can leave the target only where the refinement is
        # non-congruent AND the target is locally concave, which California's
        # coastline at these depths is not, but a hole in the target would be.
        @test brute ⊆ expanded
        @test expanded == brute
        if DGG.has_sorted_subtrees(sys)
            ranges = DGG.level_ranges(coarse_set, coarse)
            @test issorted(first.(ranges))
            @test all(k -> first(ranges[k]) > last(ranges[k-1]) + 1, 2:length(ranges))
            @test sum(length, ranges) == length(brute)
            # Expanding to a level above the set's own cells would need an
            # ancestor, which covers more than the set does.
            @test_throws ArgumentError DGG.level_ranges(coarse_set, first(DGG.levels(sys)))
        else
            @test_throws ArgumentError DGG.level_ranges(coarse_set, coarse)
        end
    end
end


# Fraction of the target's interior that lies in NO emitted cell's polygon.
function sliver_fraction(sys, set, geom, n::Int)
    t = prepare(geom)
    rings = [FB.open_ring(DGG.cell_boundary(set, c)) for c in set]
    caps = [FB.points_cap(view(r[1], 1:r[2])) for r in rings]
    ext = GI.extent(geom)
    x0, x1 = ext.X
    y0, y1 = ext.Y
    total = 0
    missed = 0
    for i in 0:n, j in 0:n
        lon = x0 + (x1 - x0) * i / n
        lat = y0 + (y1 - y0) * j / n
        inside(t, lon, lat) || continue
        total += 1
        p = FB.unit_point(lon, lat)
        hit = false
        for k in eachindex(rings)
            FB.cap_contains(caps[k], p) || continue
            if FB.point_in_cell(view(rings[k][1], 1:rings[k][2]), p) !== false
                hit = true
                break
            end
        end
        hit || (missed += 1)
    end
    return missed, total
end

@testset "the drawn cells tile the target only where the refinement does: $(syslabel(sys))" for
    (sys, leaf, _) in SWEEP

    congruent = isquadface(sys)

    set = DGG.query(sys, DGG.MultiOrderCoverage(MAINLAND); level=leaf)
    missed, total = sliver_fraction(sys, set, MAINLAND, 90)
    @test total > 2000
    if congruent
        # Four children tile their parent, so a subtree and its root have the
        # same footprint and the mixed-level set is a genuine cover.
        @test missed == 0
    else
        @test missed > 0
        @test missed / total < 0.25
    end
end


@testset "multipolygon: the offshore islands: $(syslabel(sys))" for (sys, leaf, _) in SWEEP
    set = DGG.query(sys, DGG.MultiOrderCoverage(CALIFORNIA); level=leaf)
    @test isempty(emitted_ancestors(sys, set))
    # Every one of the eight parts has to be covered, not just the big one. The
    # smallest is Santa Barbara Island, 2.6 km across — smaller than a leaf cell
    # on some of these systems, which is exactly the case that gets dropped.
    for (k, part) in enumerate(PARTS)
        poly = topolygon(part)
        t = prepare(poly)
        samples = vcat(interior_samples(poly, t, 24), boundary_samples(poly, t))
        @test !isempty(samples)
        @test isempty(uncovered(sys, set, leaf, samples)) || error("part $k uncovered")
    end
    # The islands add cells rather than replacing them.
    @test length(set) > length(DGG.query(sys, DGG.MultiOrderCoverage(MAINLAND); level=leaf))
end

@testset "hole: $(syslabel(sys))" for (sys, leaf, _) in SWEEP
    t = prepare(DONUT)
    set = DGG.query(sys, DGG.MultiOrderCoverage(DONUT); level=leaf)
    @test isempty(emitted_ancestors(sys, set))

    # The ring between the two boundaries is covered ...
    samples = vcat(interior_samples(DONUT, t, 64), boundary_samples(DONUT, t))
    @test length(samples) > 500
    @test isempty(uncovered(sys, set, leaf, samples))

    # ... and no emitted cell lies inside the hole. Cells that STRADDLE the
    # hole's edge are emitted, and must be: they meet the target.
    for c in set
        g = DGG.levelgrid(sys, DGG.level(c))
        @test !EN._matches(DGG.Within(nothing), prepare(HOLE), g, c)
    end
end

@testset "antimeridian: $(syslabel(sys))" for (sys, leaf, _) in SWEEP
    # Nothing in the engine works in longitude — a ring becomes unit-sphere
    # points at the boundary of the call and its edges are great-circle arcs —
    # so a seam-crossing target is not a special case. This asserts that it
    # stays that way.
    t = prepare(SEAM)
    @test inside(t, 179.0, -17.5)
    @test inside(t, -179.0, -17.5)
    @test !inside(t, 170.0, -17.5)

    set = DGG.query(sys, DGG.MultiOrderCoverage(SEAM); level=leaf)
    @test !isempty(set)
    @test isempty(emitted_ancestors(sys, set))
    samples = [(lon > 180 ? lon - 360 : lon, lat)
               for lon in 176.0:0.2:182.0, lat in -19.0:0.2:-16.0]
    samples = [q for q in samples if inside(t, q...)]
    @test length(samples) > 100
    @test isempty(uncovered(sys, set, leaf, samples))
end

@testset "a target larger than a hemisphere: $(syslabel(sys))" for (sys, _, coarse) in SWEEP
    # 66% of the sphere. Its bounding cap is the whole sphere — the antipode of
    # any cap enclosing the boundary is interior to the target — so the cheap
    # cap prune prunes nothing, and the traversal stays output-sensitive only
    # because of the boundary-arc prune. Without that this testset would visit
    # every cell of `coarse`.
    t = prepare(WIDE)
    @test t.cap.radius >= Float64(pi)
    @test inside(t, 0.0, 50.0) && inside(t, 0.0, -50.0)
    @test !inside(t, 0.0, 0.0)

    set = DGG.query(sys, DGG.MultiOrderCoverage(WIDE); level=coarse)
    @test !isempty(set)
    @test isempty(emitted_ancestors(sys, set))
    @test length(set) < DGG.ncells(DGG.levelgrid(sys, coarse))

    samples = [(lon, lat) for lon in -170.0:20.0:170.0, lat in -85.0:5.0:85.0]
    samples = [q for q in samples if inside(t, q...)]
    @test length(samples) > 200
    @test isempty(uncovered(sys, set, coarse, samples))

    # And the equatorial band, which is outside the target, is not covered.
    equator = [DGG.cellat(DGG.levelgrid(sys, coarse), lon, 0.0) for lon in -180.0:10.0:170.0]
    emitted = Set(collect(set))
    for c in equator
        c === nothing && continue
        @test !(c in emitted)
    end
end

# ---------------------------------------------------------------------------
# The prune that makes the two previous testsets affordable
# ---------------------------------------------------------------------------

@testset "the boundary-arc prune changes speed, not answers" begin
    # `_subtree_outside` is a proof, so switching it off must reproduce the
    # coverage exactly. This runs the traversal both ways on the mainland and
    # compares, on the two systems whose refinements are furthest apart.
    function unpruned(sys, geom, maxlevel)
        target = EN._query_target(geom)
        cells = DGG.cellindextype(sys)[]
        top = first(DGG.levels(sys))
        grids = [DGG.levelgrid(sys, l) for l in top:maxlevel]
        function visit(c)
            GO.Extents.intersects(target.cap, DGG.node_extent(sys, c)) || return nothing
            lc = DGG.level(c)
            grid = grids[lc-top+1]
            meets = EN._matches(DGG.Intersects(nothing), target, grid, c)
            if lc >= maxlevel
                meets && push!(cells, c)
                return nothing
            end
            if meets && EN._matches(DGG.Within(nothing), target, grid, c)
                push!(cells, c)
                return nothing
            end
            for ch in DGG.children(sys, c)
                visit(ch)
            end
            return nothing
        end
        for c in DGG.rootcells(sys)
            visit(c)
        end
        return sort!(cells; by=c -> (DGG.level(c), c))
    end

    for (sys, level) in ((DGG.IGeo7System(), 5), (DGG.S2System(), 7))
        set = DGG.query(sys, DGG.MultiOrderCoverage(MAINLAND); level)
        @test sort(collect(set); by=c -> (DGG.level(c), c)) == unpruned(sys, MAINLAND, level)
    end

    # And it really does prune: the node extent of a cell on the far side of the
    # planet is provably outside California.
    target = EN._query_target(MAINLAND)
    far = DGG.cellat(DGG.levelgrid(DGG.S2System(), 3), 60.0, 20.0)
    @test EN._subtree_outside(target, DGG.node_extent(DGG.S2System(), far))
    near = DGG.cellat(DGG.levelgrid(DGG.S2System(), 6), -120.0, 37.0)
    @test !EN._subtree_outside(target, DGG.node_extent(DGG.S2System(), near))
    # A cap target carries no boundary arcs, so no proof is available and the
    # traversal keeps its cap prune alone.
    cap = GO.UnitSpherical.SphericalCap(FB.unit_point(0.0, 0.0), 0.1)
    @test !EN._subtree_outside(EN._query_target(cap),
        DGG.node_extent(DGG.S2System(), far))
end

end # module MultiOrderPolygonTests
