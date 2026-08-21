module RHEALPixSystemTests

using Test
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
import GeometryOps as GO
using DiscreteGlobalGridsConformanceTesting

const RH = DiscreteGlobalGrids.RHEALPix
const SYS = RH.RHEALPixSystem()
const ORACLES = normpath(joinpath(@__DIR__, "..", "..", "oracles", "rhealpix"))

csvrows(name) = (split(line, ','; keepempty=true) for
    line in readlines(joinpath(ORACLES, name))[2:end])

function json_string(line, name)
    match(Regex("\\\"" * name * "\\\":\\\"([^\\\"]+)\\\""), line).captures[1]
end

function unitpoint(lon, lat)
    c = cos(lat)
    return GO.UnitSphericalPoint(c * cos(lon), c * sin(lon), sin(lat))
end

cross3(a, b) = (a[2] * b[3] - a[3] * b[2],
                a[3] * b[1] - a[1] * b[3],
                a[1] * b[2] - a[2] * b[1])

function ccw_measure(points)
    acc = (0.0, 0.0, 0.0)
    outward = (0.0, 0.0, 0.0)
    for i in eachindex(points)
        acc = acc .+ cross3(points[i], points[mod1(i + 1, length(points))])
        outward = outward .+ Tuple(points[i])
    end
    return sum(acc .* outward)
end

@testset "SUID codec and hierarchy" begin
    for text in ("N", "O0", "P28", "Q444", "R7751215231", "S08")
        c = RH.parse_suid(text)
        @test RH.suid(c) == text
        @test RH.RHEALPixCell(text) == c
        @test isbitstype(typeof(c))
    end
    @test_throws ArgumentError RH.parse_suid("")
    @test_throws ArgumentError RH.parse_suid("A0")
    @test_throws ArgumentError RH.parse_suid("N9")

    worst = 0.0
    for a in csvrows("hierarchy.csv")
        l, ordinal, text = parse(Int, a[1]), parse(Int64, a[2]), a[3]
        c = RH.RHEALPixCell(l, ordinal)
        @test RH.suid(c) == text
        @test cellposition(levelgrid(SYS, l), c) == ordinal + 1
        if l == 0
            @test isempty(a[4])
        else
            @test RH.suid(parent(SYS, c)) == a[4]
        end
        @test RH.suid.(children(SYS, c)) == split(a[5])
        got = RH.cell_rectangle(SYS, c)
        expected = parse.(Float64, a[6:8])
        worst = max(worst, maximum(abs.(collect(got) .- expected)))
    end
    @test worst < 1e-14
    @test ncells(levelgrid(SYS, 19)) == 6 * Int64(9)^19
    @test last(levels(SYS)) == 19
    deepest = RH.RHEALPixCell("S" * repeat("8", 19))
    @test cellat(levelgrid(SYS, 19), cell_centroid(SYS, deepest)) == deepest
end

@testset "unit-sphere projection oracle: all polar placements" begin
    forward_error = inverse_error = 0.0
    for a in csvrows("projection_unit_all_polar_placements.csv")
        north, south = parse(Int, a[1]), parse(Int, a[2])
        lon, lat, ex, ey, elon, elat = parse.(Float64, a[4:9])
        x, y = RH.rhealpix_forward(lon, lat;
            north_square=north, south_square=south)
        ilon, ilat = RH.rhealpix_inverse(x, y;
            north_square=north, south_square=south)
        forward_error = max(forward_error, abs(x - ex), abs(y - ey))
        inverse_error = max(inverse_error, abs(ilon - elon), abs(ilat - elat))
    end
    @test forward_error < 2e-15
    @test inverse_error < 2e-15
    @test RH.in_rhealpix_image(0, 0)
    @test !RH.in_rhealpix_image(0, pi / 2)
    @test_throws DomainError RH.rhealpix_inverse(0, pi / 2)
    @test_throws DomainError RH.healpix_inverse(4pi, 0)
    @test_throws DomainError RH.healpix_inverse(0, 1.4)
end

@testset "AusPIX WGS84 projection profile" begin
    forward_error = inverse_error = 0.0
    for a in csvrows("projection_auspix_wgs84.csv")
        lon, lat, ex, ey, elon, elat = parse.(Float64, a[2:7])
        x, y = RH.auspix_forward(lon, lat)
        ilon, ilat = RH.auspix_inverse(x, y)
        forward_error = max(forward_error, abs(x - ex), abs(y - ey))
        inverse_error = max(inverse_error, abs(ilon - elon), abs(ilat - elat))
    end
    @test forward_error < 5e-8
    @test inverse_error < 5e-13

    auspix = RH.AusPIXSystem()
    @test cellindextype(auspix) == RH.RHEALPixCell
    @test RH.suid.(rootcells(auspix)) == collect(string.(RH.ROOT_NAMES))
    @test parent(auspix) == SYS
end

@testset "point lookup oracle and boundary ownership" begin
    bad = 0
    for a in csvrows("lookup.csv")
        profile, l = a[1], parse(Int, a[2])
        lon, lat = parse(Float64, a[3]), parse(Float64, a[4])
        if profile == "rhealpix_unit_003_00"
            got = cellat(levelgrid(SYS, l), unitpoint(lon, lat))
        else
            got = cellat(levelgrid(RH.AusPIXSystem(), l), rad2deg(lon), rad2deg(lat))
        end
        bad += RH.suid(got) != a[5]
    end
    @test bad == 0

    # Non-exact normal/quadrant probes seal both sides of every root and level-1
    # edge.  Exact points are deliberately excluded: the implementation uses
    # the paper's depth-independent half-open convention, whereas the oracle's
    # exact tie nudge depends on its configured artificial maximum depth.
    bad = total = 0
    for a in csvrows("boundary_ties.csv")
        a[14] == "true" || continue
        isempty(a[19]) && continue
        a[9] == "exact" && continue
        l = parse(Int, a[2])
        lon, lat = parse(Float64, a[15]), parse(Float64, a[16])
        if a[1] == "rhealpix_unit_003_00"
            got = cellat(levelgrid(SYS, l), unitpoint(lon, lat))
        else
            got = cellat(levelgrid(RH.AusPIXSystem(), l), rad2deg(lon), rad2deg(lat))
        end
        total += 1
        bad += RH.suid(got) != a[19]
    end
    @test total == 9408
    @test bad == 0
end

@testset "cell geometry, edge topology, and vertex topology" begin
    edge_bad = 0
    geometry_error = 0.0
    polygon_area_error = 0.0
    for line in eachline(joinpath(ORACLES, "cells.jsonl"))
        json_string(line, "profile") == "rhealpix_unit_003_00" || continue
        text = json_string(line, "cell_id")
        m = match(r"\"edge_neighbors\":\{\"down\":\"([^\"]+)\",\"left\":\"([^\"]+)\",\"right\":\"([^\"]+)\",\"up\":\"([^\"]+)\"\}", line)
        expected = (m.captures[4], m.captures[2], m.captures[1], m.captures[3])
        c = RH.parse_suid(text)
        g = levelgrid(SYS, level(c))
        edge_bad += RH.suid.(RH.edge_neighbors(g, c)) != expected
        @test cellat(g, cell_centroid(g, c)) == c
        boundary = cell_boundary(g, c)
        @test length(boundary) == 32
        @test ccw_measure(boundary) > 0
        @test cell_area(g, c) ≈ 2pi / (3 * 9^level(c)) rtol=2e-15
        polygon_area = GO.area(GO.Spherical(radius=1.0), cell_polygon(g, c))
        polygon_area_error = max(polygon_area_error,
            abs(polygon_area - cell_area(g, c)) / cell_area(g, c))

        bm = match(r"\"boundary_unit_directions_authalic\":\[(.*?)\],\"boundary_winding\"", line)
        expected_boundary = [parse.(Float64, point.captures) for point in
            eachmatch(r"\[([-+0-9.eE]+),([-+0-9.eE]+),([-+0-9.eE]+)\]", bm.captures[1])]
        @test length(expected_boundary) == 16
        # The sealed corpus has four samples per edge; production geometry has
        # eight. Every second production vertex is therefore an oracle point.
        for (got, expected) in zip(boundary[1:2:end], expected_boundary)
            geometry_error = max(geometry_error,
                maximum(abs.(collect(Tuple(got)) .- expected)))
        end
        nm = match(r"\"nucleus_unit_direction_authalic\":\[([-+0-9.eE]+),([-+0-9.eE]+),([-+0-9.eE]+)\]", line)
        expected_nucleus = parse.(Float64, nm.captures)
        geometry_error = max(geometry_error,
            maximum(abs.(collect(Tuple(cell_centroid(g, c))) .- expected_nucleus)))
    end
    @test edge_bad == 0
    @test geometry_error < 2e-14
    @test polygon_area_error < 7e-3

    vertex_bad = 0
    for line in eachline(joinpath(ORACLES, "vertex_neighbors.jsonl"))
        text = json_string(line, "cell_id")
        m = match(r"\"vertex_neighbor_ids_sorted\":\[([^\]]*)\]", line)
        expected = sort([x.captures[1] for
            x in eachmatch(r"\"([^\"]+)\"", m.captures[1])])
        c = RH.parse_suid(text)
        got = sort(RH.suid.(collect(neighbors(levelgrid(SYS, level(c)), c))))
        vertex_bad += got != expected
    end
    @test vertex_bad == 0
end

@testset "all 16 placement topologies" begin
    for north in 0:3, south in 0:3, l in 0:2
        sys = RH.RHEALPixSystem(north, south)
        g = levelgrid(sys, l)
        for i in 1:ncells(g)
            c = cellindex(g, i)
            adjacent = neighbors(g, c; connectivity=Edge())
            @test length(adjacent) == 4
            @test length(unique(adjacent)) == 4
            @test c ∉ adjacent
            @test all(c ∈ neighbors(g, x; connectivity=Edge()) for x in adjacent)
        end
    end
end

@testset "package interface conformance" begin
    for l in 0:2
        test_grid_interface(levelgrid(SYS, l); label="rHEALPix level $l")
    end
    test_hierarchical_system(SYS)
end

end # module
