# test/core/test_tree_queries.jl — the spherical tree query behind the
# `Touching` selectors and `zonal` (src/core/lookup_ops.jl): equivalence
# against an independent brute-force spherical oracle on all four lookup
# types — antimeridian-crossing and pole-enclosing geometries included, the
# cases the old planar predicates documented away — plus the query's
# supporting API: the `node_indices` cursor accessor, the spherical geometry
# cap, the rim sandwich that keeps `:touches` off the polygon predicate, the
# congruence trait, and a descent-vs-brute-force perf sanity bound.

module TreeQueryTests

using Test
using Random
import DimensionalData as DD
using DimensionalData: Lookups
import GeometryOps as GO
import GeometryOps.SpatialTreeInterface as STI
import GeoInterface as GI

using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
using DiscreteGlobalGrids.HEALPix.HealpixLookups: HealpixLookup, HealpixLookups
using DiscreteGlobalGrids.H3.H3Lookups: H3Lookup
using DiscreteGlobalGrids.A5.A5Lookups: A5Lookup
using DiscreteGlobalGrids.IGeo7.IGeo7Lookups: IGeo7Lookup

const SPHERICAL = GO.RelateNG(; manifold=GO.Spherical())

# The independent slow oracle: one unprepared-per-geometry prepared engine,
# every stored cell tested exactly — no caps, no tree, no shared code with
# the descent beyond the kernel geometry itself.
function brute_positions(l, geom, mode::Symbol)
    system = dggs_system(l)
    level = dggs_level(l)
    prep = GO.prepare(SPHERICAL, geom)
    ids = DD.parent(l)
    if mode === :center
        return findall(i -> GO.relate_predicate(prep, GO.pred_contains(),
            DGG.cell_center(system, level, ids[i])), eachindex(ids))
    end
    return findall(i -> GO.relate_predicate(prep, GO.pred_intersects(),
        cell_polygon_unitsphere(system, level, ids[i])), eachindex(ids))
end

lonlat_ring(points) = GI.Polygon([GI.LinearRing(points)])

box_ring(lon0, lat0, w, h) = lonlat_ring([
    (lon0, lat0), (lon0 + w, lat0), (lon0 + w, lat0 + h), (lon0, lat0 + h), (lon0, lat0)])

# A polygon straddling ±180 (written continuously past 180, as spherical
# geometry allows) and one enclosing the north pole — the two geometry
# families the pre-spherical planar predicates could not answer.
const ANTIMERIDIAN = box_ring(170.0, -12.0, 22.0, 25.0)
const POLAR_CAP = lonlat_ring([(lon, 76.0) for lon in 0.0:15.0:360.0])
# A triangle with deliberately long great-circle edges: planar and spherical
# predicates disagree on much of its rim, so agreement with the oracle here
# is agreement about *spherical* semantics, not just about small shapes.
const BIG_TRIANGLE = lonlat_ring([(-40.0, 10.0), (35.0, 25.0), (-5.0, 65.0), (-40.0, 10.0)])
# A 70°-long, 0.3°-tall sliver: nearly every cell its bounding cap admits is
# far from the region, the worst case for the cap prune and the best case for
# the rim sandwich (`_sandwich`) that backs it up.
const THIN_BAND = box_ring(-20.0, 5.0, 70.0, 0.3)
# A 200-vertex ring — the many-edge query geometry, where the sandwich's
# per-cell boundary scan costs the most relative to the predicate it saves.
const WOBBLY = lonlat_ring([
    let angle = 2 * pi * i / 200, radius = 9.0 + 2.5 * sin(5 * angle)
        (-25.0 + radius * cos(angle), 15.0 + 0.7 * radius * sin(angle))
    end for i in 0:200])

# One moderately sized lookup per system — small enough that the O(n · cost)
# oracle stays cheap, global or near-global so the fixed geometries hit.
function test_lookups()
    hp = HealpixLookup(collect(Int64, 0:(12 * 4^4 - 1)); level=4)
    h3 = H3Lookup(DGGSGlobeIds(H3DGGS(), 1))
    ig = IGeo7Lookup(sort!(reduce(vcat,
        [collect(cell_descendants(IGEO7DGGS(), 0, r, 2)) for r in root_ids(IGEO7DGGS())]));
        resolution=2)
    a5 = A5Lookup(sort!(reduce(vcat,
        [collect(cell_descendants(A5DGGS(), 0, r, 2)) for r in root_ids(A5DGGS())]));
        resolution=2)
    return (hp, h3, ig, a5)
end

# The corpus the equivalence testsets share: the fixed geometries whose
# families the planar predicates could not answer, plus random boxes. `rng` is
# threaded rather than seeded here so that one stream spans the whole loop over
# lookups and every system draws *different* boxes; the testsets below walk the
# same lookups in the same order from the same seed, so they still see the same
# corpus as each other.
function fuzz_geometries(rng)
    geoms = Any[ANTIMERIDIAN, POLAR_CAP, BIG_TRIANGLE, THIN_BAND, WOBBLY]
    for _ in 1:12
        lon0 = rand(rng) * 340 - 170
        lat0 = rand(rng) * 130 - 75
        push!(geoms, box_ring(lon0, lat0, rand(rng) * 40 + 2, rand(rng) * 25 + 2))
    end
    return geoms
end

# The descent with the rim sandwich switched off, i.e. every cell whose center
# falls outside the geometry paying the exact `pred_intersects` — the code path
# that shipped before `_sandwich` existed, and the reference the sandwich must
# reproduce position for position.
function sandwich_off(l, geom, mode::Symbol)
    system = dggs_system(l)
    leaf = dggs_level(l)
    ids = DD.parent(l)
    tree = treeify(DGGSPartialGrid(l; bucket_size=DGG.QUERY_BUCKET_SIZE))
    prep = GO.prepare(SPHERICAL, geom)
    cap = DGG._geometry_cap(prep, geom)
    out = Int[]
    DGG._tree_query!(out, tree, ids, system, leaf, prep, cap, nothing, mode)
    return sort!(out)
end

@testset "spherical tree queries" begin
    @testset "descent == spherical brute force, all four systems" begin
        rng = MersenneTwister(1234)
        for l in test_lookups()
            for geom in fuzz_geometries(rng), mode in (:center, :touches)
                got = DGG._query_positions(l, geom, mode)
                want = brute_positions(l, geom, mode)
                @test got == want
            end
            # ...and the fixed geometries really do select something, on every
            # system — an empty-vs-empty agreement would prove nothing.
            for geom in (ANTIMERIDIAN, POLAR_CAP, BIG_TRIANGLE, THIN_BAND, WOBBLY)
                @test !isempty(DGG._query_positions(l, geom, :touches))
            end
        end
    end

    @testset "the query IS the selector and zonal path" begin
        # Same seed geometry through the three public doors: Touching,
        # Contains(geometry), zonal — one query, three surfaces.
        hp = HealpixLookup(collect(Int64, 0:(12 * 4^4 - 1)); level=4)
        touching = Lookups.selectindices(hp, HealpixLookups.Touching(BIG_TRIANGLE))
        @test touching == DGG._query_positions(hp, BIG_TRIANGLE, :touches)
        centerhits = Lookups.selectindices(hp, Lookups.Contains(BIG_TRIANGLE))
        @test centerhits == DGG._query_positions(hp, BIG_TRIANGLE, :center)
        @test issubset(centerhits, touching)
        A = DD.DimArray(ones(length(hp)), HealpixLookups.Cells(hp))
        @test zonal(sum, A; of=[BIG_TRIANGLE])[1] == length(centerhits)
        @test zonal(sum, A; of=[BIG_TRIANGLE], boundary=:touches)[1] == length(touching)
    end

    @testset "A5 zonal works through the selection cursor" begin
        # A5 has no descendant ranges; the old kernel descent threw
        # NotPortedError here. The tree query runs on the selection-cursor
        # fallback instead, so A5 gets zonal like everyone else.
        a5 = A5Lookup(sort!(reduce(vcat,
            [collect(cell_descendants(A5DGGS(), 0, r, 1)) for r in root_ids(A5DGGS())]));
            resolution=1)
        A = DD.DimArray(ones(length(a5)), DD.Dim{:cells}(a5))
        z = zonal(sum, A; of=[POLAR_CAP], boundary=:touches)
        @test z[1] == length(brute_positions(a5, POLAR_CAP, :touches)) > 0
    end

    @testset "node_indices" begin
        # partial grid, range cursors: the root owns everything, children
        # partition it, and the indices are the same leaf index space
        # Trees.getcell addresses.
        ids = Int64[3, 5, 9, 17, 33, 100, 101]
        tree = treeify(DGGSPartialGrid(HEALPixDGGS(), 3, ids))
        @test node_indices(tree) == 1:7
        seen = Int[]
        for child in STI.getchild(tree)
            append!(seen, node_indices(child))
        end
        @test sort(seen) == 1:7
        # selection cursors (A5 path) expose their materialized vector
        a5ids = sort(root_ids(A5DGGS()))
        a5tree = treeify(DGGSPartialGrid(A5DGGS(), 0, a5ids))
        @test sort(reduce(vcat, [collect(node_indices(c)) for c in STI.getchild(a5tree)])) ==
              collect(1:length(a5ids))
        # dense grids: ordinal windows
        dense = treeify(DGGSGrid(HEALPixDGGS(), 2))
        @test node_indices(dense) == 1:Int(DGG.num_cells(HEALPixDGGS(), 2))
        firstchild = first(STI.getchild(dense))
        @test node_indices(firstchild) == 1:16
    end

    @testset "spherical geometry cap is conservative" begin
        rng = MersenneTwister(23)
        slerp(a, b, t) = begin
            v = a * (1 - t) + b * t
            n = sqrt(sum(abs2, v))
            GO.UnitSphericalPoint(v[1] / n, v[2] / n, v[3] / n)
        end
        for _ in 1:200
            lon0 = rand(rng) * 300 - 150
            lat0 = rand(rng) * 140 - 80
            w = rand(rng) * 25 + 0.1
            h = min(rand(rng) * 15 + 0.1, 89 - lat0)
            b = box_ring(lon0, lat0, w, h)
            prep = GO.prepare(SPHERICAL, b)
            cap = DGG._geometry_cap(prep, b)
            in_cap(p) = GO.UnitSpherical.spherical_distance(cap.point, p) <= cap.radius
            # the geometry's boundary — vertices and points along the
            # great-circle edges — lies inside the cap
            vertices = [DGG._unit_sphere_point(GI.x(p), GI.y(p))
                        for p in GI.getpoint(GI.getexterior(b))]
            for i in 1:(length(vertices) - 1), t in (0.0, 0.25, 0.5, 0.75)
                @test in_cap(slerp(vertices[i], vertices[i + 1], t))
            end
            # ...and so does the interior, which is the half of the claim the
            # prune actually rests on: the cap may not clip any part of the
            # region. Sampled points are only *asked* about when the engine
            # agrees they are in the geometry, so this is exactly the
            # implication "in the geometry ⟹ in the cap" and never a guess
            # about where a great-circle edge runs relative to a parallel.
            for _ in 1:8
                p = DGG._unit_sphere_point(lon0 + rand(rng) * w, lat0 + rand(rng) * h)
                @test !GO.relate_predicate(prep, GO.pred_intersects(), p) || in_cap(p)
            end
        end
        # A ring whose lon-lat coordinates *look* like a 260°-wide box is not
        # wide on the sphere: great-circle edges take the short arcs across
        # the antimeridian and the interior is the smaller region, so this is
        # really the ~100°-wide patch around lon 180 — and the cap is a
        # genuine one (radius < π/2) that covers that patch, not a give-up.
        wide = lonlat_ring([(-130.0, -40.0), (130.0, -40.0), (130.0, 40.0),
                            (-130.0, 40.0), (-130.0, -40.0)])
        wprep = GO.prepare(SPHERICAL, wide)
        wcap = DGG._geometry_cap(wprep, wide)
        @test wcap.radius < Float64(pi) / 2
        @test GO.UnitSpherical.spherical_distance(wcap.point,
            DGG._unit_sphere_point(180.0, 0.0)) <= wcap.radius
        # A geometry whose vertex spread genuinely passes a quarter sphere
        # gives up — the convexity argument needs radius <= π/2. A thin
        # equatorial band densified across 240° of longitude has its vertex
        # mean near (0°, 0°) and vertices 120° away from it.
        band = lonlat_ring(vcat(
            [(lon, 2.0) for lon in -120.0:10.0:120.0],
            [(lon, -2.0) for lon in 120.0:-10.0:-120.0],
            [(-120.0, 2.0)]))
        bprep = GO.prepare(SPHERICAL, band)
        @test DGG._geometry_cap(bprep, band).radius >= Float64(pi)
        # a polar ring's cap covers the pole its interior encloses, so cap
        # pruning cannot drop polar cells (the old planar cap could not
        # express this region at all)
        pprep = GO.prepare(SPHERICAL, POLAR_CAP)
        pcap = DGG._geometry_cap(pprep, POLAR_CAP)
        @test pcap.radius < Float64(pi) / 2   # a real cap, not the give-up answer
        north = GO.UnitSphericalPoint(0.0, 0.0, 1.0)
        @test GO.UnitSpherical.spherical_distance(pcap.point, north) <= pcap.radius
    end

    @testset "rim sandwich" begin
        # 1. The whole point: turning the prefilter on changes no answer, on
        #    the same corpus and every system, in both modes. The equality is
        #    exact — `_sandwich` may only skip the polygon predicate when it
        #    has proved what that predicate would have said.
        rng = MersenneTwister(1234)
        for l in test_lookups()
            for geom in fuzz_geometries(rng), mode in (:center, :touches)
                @test DGG._query_positions(l, geom, mode) == sandwich_off(l, geom, mode)
            end
        end

        # 2. Verdict by verdict, against the predicate it replaces — over
        #    *every* stored cell, not just the ones the cap prune admits, so
        #    the reject arm is exercised on distant cells too. Counted rather
        #    than asserted per cell to keep the summary readable.
        for l in test_lookups()
            system = dggs_system(l)
            level = dggs_level(l)
            ids = DD.parent(l)
            total_accepted = 0
            for geom in (ANTIMERIDIAN, POLAR_CAP, BIG_TRIANGLE, THIN_BAND, WOBBLY,
                         box_ring(5.0, 40.0, 10.0, 10.0), box_ring(-3.0, -60.0, 30.0, 12.0))
                prep = GO.prepare(SPHERICAL, geom)
                arcs = DGG._boundary_arcs(geom)
                @test arcs !== nothing
                wrong = accepted = rejected = undecided = 0
                for id in ids
                    center = DGG.cell_center(system, level, id)
                    # the sandwich only ever runs behind a failed center test
                    GO.relate_predicate(prep, GO.pred_contains(), center) && continue
                    ring = DGG.cell_boundary(system, level, id; closed=true)
                    verdict = DGG._sandwich(arcs, center, ring)
                    if verdict == 0
                        undecided += 1
                        continue
                    end
                    verdict > 0 ? (accepted += 1) : (rejected += 1)
                    truth = GO.relate_predicate(prep, GO.pred_intersects(),
                        cell_polygon_unitsphere(system, level, id))
                    (verdict > 0) == truth || (wrong += 1)
                end
                total_accepted += accepted
                @test wrong == 0
                # the reject arm may never be vacuous — most of the globe is
                # far from any of these geometries...
                @test rejected > 0
                # ...and the annulus it leaves behind must be a small
                # remainder, or the prefilter would not be worth its cost
                @test undecided < rejected
            end
            # The accept arm needs the query boundary to pass within a cell's
            # inradius, which on these deliberately coarse lookups some single
            # geometry can miss entirely; over the corpus it must fire.
            @test total_accepted > 0
        end

        # 3. The sandwich holds the cell's boundary ring anyway, so the exact
        #    test it falls through to builds the cell polygon from that ring
        #    instead of asking twice. That is only the same test while no
        #    system overrides `cell_polygon_unitsphere` away from its
        #    `cell_boundary`-derived default — pin it.
        for l in test_lookups()
            system = dggs_system(l)
            level = dggs_level(l)
            for id in DD.parent(l)[1:37:end]
                ring = DGG.cell_boundary(system, level, id; closed=true)
                @test collect(GI.getpoint(GI.Polygon([GI.LinearRing(ring)]))) ==
                      collect(GI.getpoint(cell_polygon_unitsphere(system, level, id)))
            end
        end

        # 4. The edge set is the *boundary*, all rings of all parts. A hole's
        #    ring counts (a cell straddling it must still reach the exact
        #    test), and points across two rings are never joined into an edge.
        square(lon0, lat0, w) = GI.LinearRing([
            (lon0, lat0), (lon0 + w, lat0), (lon0 + w, lat0 + w), (lon0, lat0 + w), (lon0, lat0)])
        holed = GI.Polygon([square(0.0, 0.0, 10.0), square(3.0, 3.0, 4.0)])
        # a closed n-point ring contributes n arcs: n-1 edges plus a
        # degenerate closing one
        @test length(DGG._boundary_arcs(holed)) == 10
        multi = GI.MultiPolygon([GI.Polygon([square(0.0, 0.0, 10.0)]),
                                 GI.Polygon([square(40.0, 0.0, 10.0)])])
        @test length(DGG._boundary_arcs(multi)) == 10
        @test length(DGG._boundary_arcs(GI.Point(1.0, 2.0))) == 1
        # A geometry whose coordinates are not Float64 — GeoJSON reads Natural
        # Earth as Float32 — must still land in the concrete `BoundaryArc`
        # fields rather than throwing a MethodError.
        square32(lon0, lat0, w) = GI.LinearRing([
            (Float32(lon0), Float32(lat0)), (Float32(lon0 + w), Float32(lat0)),
            (Float32(lon0 + w), Float32(lat0 + w)), (Float32(lon0), Float32(lat0 + w)),
            (Float32(lon0), Float32(lat0))])
        holed32 = GI.Polygon([square32(0.0, 0.0, 10.0), square32(3.0, 3.0, 4.0)])
        @test length(DGG._boundary_arcs(holed32)) == 10
        @test eltype(DGG._boundary_arcs(holed32)) === DGG.BoundaryArc
        # kinds the walk does not vouch for switch the prefilter off rather
        # than guess, and the query still answers exactly
        @test DGG._boundary_arcs(square(0.0, 0.0, 10.0)) === nothing
        hp = HealpixLookup(collect(Int64, 0:(12 * 4^4 - 1)); level=4)
        @test DGG._query_positions(hp, holed, :touches) ==
              brute_positions(hp, holed, :touches)
        @test DGG._query_positions(hp, multi, :touches) ==
              brute_positions(hp, multi, :touches)
        # the hole really is respected: cells strictly inside it are selected
        # by `:touches` only where the hole's rim crosses them
        @test length(DGG._query_positions(hp, holed, :touches)) >
              length(DGG._query_positions(hp, holed, :center))
    end

    @testset "congruence trait" begin
        @test has_congruent_geometry(HEALPixDGGS())
        for system in (H3DGGS(), IGEO7DGGS(), A5DGGS(), S2DGGS())
            @test !has_congruent_geometry(system)
        end
        # ...and the default answer for the outline that goes with it
        @test subtree_polygon_unitsphere(H3DGGS(), 1, first(root_ids(H3DGGS())), 3) === nothing
    end

    @testset "perf sanity: descent beats brute force" begin
        # Full-globe HEALPix at level 6 (49,152 cells), a country-sized box.
        # The point is asymptotics, not a benchmark: the descent touches the
        # cells near the geometry, brute force touches every stored cell.
        # Measured ~100x on this workload; 3x keeps CI noise out of it.
        l = HealpixLookup(collect(Int64, 0:(12 * 4^6 - 1)); level=6)
        geom = box_ring(5.0, 40.0, 15.0, 12.0)
        DGG._query_positions(l, geom, :touches)             # compile
        t_descent = @elapsed hits = DGG._query_positions(l, geom, :touches)
        brute_positions(l[1:100], geom, :touches)           # compile oracle path
        t_brute = @elapsed truth = brute_positions(l, geom, :touches)
        @test hits == truth
        @test !isempty(hits)
        @test t_descent < t_brute / 3
    end
end

end # module TreeQueryTests
