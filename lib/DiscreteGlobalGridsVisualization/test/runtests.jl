using Test
import DiscreteGlobalGrids as DGG
import DiscreteGlobalGridsVisualization as DGGV
using DiscreteGlobalGridsVisualization
using GeometryBasics: Point2d, Point3d
import GeometryBasics
using Makie
using CairoMakie
using GeoMakie
import Proj
import DimensionalData as DD

CairoMakie.activate!(type = "png")

const SYS = DGG.IGeo7System()
const ALPS = DGG.Extents.Extent(X = (10.0, 11.0), Y = (46.0, 47.0))

patch(level) = DGG.CellVector(DGG.query(SYS, DGG.MultiOrderCoverage(ALPS); level))
planar(cut = 180.0) = DGGV.PlanarTarget(identity, cut)

lons(mesh) = [p[1] for p in mesh.positions]
lats(mesh) = [p[2] for p in mesh.positions]

"A one-cell grid whose cell is a hexagon centred on the north pole."
struct PolarSource end

DGG.cell_boundary(::PolarSource, ::Int) =
    [DGG.UnitSphericalPoint(cosd(80) * cosd(a), cosd(80) * sind(a), sind(80))
        for a in 0.0:60.0:300.0]

"A unit-sphere point from longitude and latitude in degrees."
sphere(lon, lat) = DGG.UnitSphericalPoint(cosd(lat) * cosd(lon), cosd(lat) * sind(lon),
    sind(lat))

"""
A one-cell grid whose cell has a corner *at* the north pole, listed first: the
shape HEALPix puts around each pole, where four cells meet at a point.
"""
struct PoleCornerSource end

DGG.cell_boundary(::PoleCornerSource, ::Int) =
    [sphere(0.0, 90.0), sphere(0.0, 85.0), sphere(45.0, 85.0), sphere(90.0, 85.0)]

"The signed area of a spherical triangle, from its three unit-sphere corners."
function spherical_area(a, b, c)
    triple = a[1] * (b[2] * c[3] - b[3] * c[2]) -
        a[2] * (b[1] * c[3] - b[3] * c[1]) +
        a[3] * (b[1] * c[2] - b[2] * c[1])
    return 2 * atan(triple, 1 + a' * b + b' * c + c' * a)
end

"The signed area of a plane triangle."
plane_area(p, q, r) = 0.5 * ((q[1] - p[1]) * (r[2] - p[2]) - (r[1] - p[1]) * (q[2] - p[2]))

"Every triangle of `mesh`, as its three positions."
corners(mesh) = ((mesh.positions[f[1]], mesh.positions[f[2]], mesh.positions[f[3]])
    for f in mesh.faces)

"Each triangle as its three vertices, sorted, so that a repeat is visible."
triangle_set(mesh) = Set(Tuple(sort([Int(GeometryBasics.value(f[j])) for j in 1:3]))
    for f in mesh.faces)

globe() = DGGV.GlobeTarget(identity, 1.0, 0.0, 0.0)

"Save `figure` to a temporary file and say whether it landed."
function saves(figure)
    path = tempname() * ".png"
    Makie.save(path, figure)
    return isfile(path)
end

@testset "DiscreteGlobalGridsVisualization" begin

    @testset "cellset" begin
        cells = patch(7)
        # Every container that names cells reads as the same set.
        @test length(DGGV.cellset(cells)) == length(cells)
        @test length(DGGV.cellset(DGG.CellLookup(cells))) == length(cells)
        @test length(DGGV.cellset(DGG.PartialGrid(cells))) == length(cells)

        grid = DGG.levelgrid(SYS, 2)
        @test length(DGGV.cellset(grid)) == DGG.ncells(grid)

        # A system with a bare vector of ids needs no level agreement.
        mixed = [DGG.cellindex(DGG.levelgrid(SYS, 1), 1), DGG.cellindex(DGG.levelgrid(SYS, 3), 1)]
        @test length(DGGV.cellset(SYS, mixed)) == 2
    end

    @testset "a regional patch" begin
        cells = patch(8)
        mesh = DGGV.tessellate(planar(), cells)

        @test mesh.ncells == length(cells)
        @test mesh.ndropped == 0
        # Away from the poles and the cut, one hexagon is six corners and four
        # triangles, and nothing is shared between cells.
        @test DGGV.nrings(mesh) == length(cells)
        @test length(mesh.positions) == 6 * length(cells)
        @test length(mesh.faces) == 4 * length(cells)
        @test length(mesh.vertex_cell) == length(mesh.positions)
        @test mesh.ring_start[end] == length(mesh.positions) + 1

        # The mesh lands where the cells are.
        @test all(10.0 - 1 .<= lons(mesh) .<= 11.0 + 1)
        @test all(46.0 - 1 .<= lats(mesh) .<= 47.0 + 1)

        # Every vertex of a ring names the same cell, and cells are in order.
        for i in 1:DGGV.nrings(mesh)
            span = mesh.ring_start[i]:(mesh.ring_start[i + 1] - 1)
            @test allequal(mesh.vertex_cell[span])
        end
        @test mesh.vertex_cell == repeat(Int32.(1:length(cells)), inner = 6)

        # Faces index inside the buffer, and each names three distinct vertices.
        for f in mesh.faces
            indices = convert.(Int, f)
            @test all(1 .<= indices .<= length(mesh.positions))
            @test allunique(indices)
        end
    end

    @testset "threading is not observable" begin
        cells = patch(9)
        serial = DGGV.tessellate(planar(), cells; ntasks = 1)
        parallel = DGGV.tessellate(planar(), cells; ntasks = 8)
        @test serial.positions == parallel.positions
        @test serial.faces == parallel.faces
        @test serial.vertex_cell == parallel.vertex_cell
        @test serial.ring_start == parallel.ring_start
    end

    @testset "empty input" begin
        cells = DGG.CellVector(SYS, 5, DGG.Z7Cell[])
        mesh = DGGV.tessellate(planar(), cells)
        @test isempty(mesh)
        @test mesh.ncells == 0
        @test DGGV.nrings(mesh) == 0
    end

    @testset "a whole grid stays inside the map" begin
        # A global grid exercises every hard case at once: cells across the cut,
        # cells containing a pole, and cells with an edge running over one.
        for level in 0:4
            grid = DGG.levelgrid(SYS, level)
            mesh = DGGV.tessellate(planar(), grid)

            @test mesh.ndropped == 0
            @test mesh.ncells == DGG.ncells(grid)
            @test all(-180.0 - 1.0e-9 .<= lons(mesh) .<= 180.0 + 1.0e-9)
            @test all(-90.0 .<= lats(mesh) .<= 90.0)
            # The poles are reached, so the map has no hole at the top or bottom.
            @test maximum(lats(mesh)) == 90.0
            @test minimum(lats(mesh)) == -90.0
            # Every cell is drawn, and splitting only ever adds rings.
            @test DGGV.nrings(mesh) >= DGG.ncells(grid)
            @test length(unique(mesh.vertex_cell)) == DGG.ncells(grid)
        end
    end

    @testset "cells are split at the cut, not smeared" begin
        grid = DGG.levelgrid(SYS, 3)
        mesh = DGGV.tessellate(planar(), grid)
        # A ring smeared across the map instead of split is the artefact this
        # whole file exists to prevent, so no drawn ring may be wider than the
        # cell it came from — and no IGEO7 cell reaches a quarter of a turn.  The
        # exception is a cell reaching a pole, which is genuinely that wide,
        # because it closes along the top of the map.
        for i in 1:DGGV.nrings(mesh)
            span = mesh.ring_start[i]:(mesh.ring_start[i + 1] - 1)
            width = maximum(lons(mesh)[span]) - minimum(lons(mesh)[span])
            @test width < 90.0 || any(abs.(lats(mesh)[span]) .== 90.0)
        end
    end

    @testset "cells with an edge over a pole" begin
        # IGEO7 puts each pole on an edge shared by two cells rather than inside
        # one.  Such an edge has no unambiguous longitude, and the four cells it
        # affects are the reason the ring is traced rather than merely unwrapped.
        for level in 0:4
            grid = DGG.levelgrid(SYS, level)
            buffer = Point2d[]
            traced = 0
            for i in 1:DGG.ncells(grid)
                ring = DGG.cell_boundary(grid, DGG.cellindex(grid, i))
                DGGV.ring_lonlat!(buffer, ring, length(ring))
                length(buffer) > length(ring) && (traced += 1)
            end
            @test traced == 4  # two cells meeting at each pole
        end

        # The traced ring reaches the pole and closes along the top of the map,
        # rather than guessing a side and spanning half the world.
        mesh = DGGV.tessellate(planar(), DGG.levelgrid(SYS, 0))
        at_pole = findall(p -> abs(p[2]) == 90.0, mesh.positions)
        @test length(unique(mesh.vertex_cell[at_pole])) == 4
    end

    @testset "a cell with a corner at a pole" begin
        # HEALPix meets four cells at a point on each pole.  A corner there has
        # no longitude of its own, and reading the `0` that `atand(0, 0)`
        # returns as one tilted every step of the ring after it: a quarter-turn
        # cell came out winding three quarters of the way around the pole, was
        # taken for a cell *containing* one, and was drawn as a band across the
        # whole map — the streak from the centre of an azimuthal map to its
        # pole.
        source = PoleCornerSource()
        ring = DGG.cell_boundary(source, 1)
        buffer = Point2d[]
        @test abs(DGGV.ring_lonlat!(buffer, ring, length(ring))) < 1.0e-9
        # The pole corner is drawn as two, one on each meridian the cell's edges
        # arrive and leave on, so that the ring closes along the top of the map.
        @test length(buffer) == length(ring) + 1
        @test count(p -> p[2] == 90.0, buffer) == 2

        mesh = DGGV.tessellate(planar(), DGGV.CellSet(source, [1]))
        @test DGGV.nrings(mesh) == 1
        @test maximum(lats(mesh)) == 90.0
        @test extrema(lons(mesh)) == (0.0, 90.0)

        # And on the grid the case comes from: no HEALPix cell contains a pole,
        # so none may be drawn as a band, at any level.
        for level in 0:3
            grid = DGG.levelgrid(DGG.HEALPixSystem(), level)
            mesh = DGGV.tessellate(planar(), grid)
            @test length(unique(mesh.vertex_cell)) == DGG.ncells(grid)
            @test maximum(lats(mesh)) == 90.0
            @test minimum(lats(mesh)) == -90.0
            for i in 1:DGGV.nrings(mesh)
                span = mesh.ring_start[i]:(mesh.ring_start[i + 1] - 1)
                width = maximum(lons(mesh)[span]) - minimum(lons(mesh)[span])
                @test width <= 90.0 + 1.0e-9   # a base cell's own width
            end
        end
    end

    @testset "a cell containing a pole" begin
        # No IGEO7 cell contains a pole, so the case gets its own grid:
        # `PolarSource` is a single hexagon centred on the north pole, whose ring
        # winds a full turn.
        mesh = DGGV.tessellate(planar(), DGGV.CellSet(PolarSource(), [1]))
        @test mesh.ndropped == 0
        @test DGGV.nrings(mesh) == 1
        # The cell is drawn as the whole band between its ring and the pole.
        @test extrema(lons(mesh)) == (-180.0, 180.0)
        @test maximum(lats(mesh)) == 90.0
        @test minimum(lats(mesh)) ≈ 80.0
        @test all(f -> all(1 .<= convert.(Int, f) .<= length(mesh.positions)), mesh.faces)
    end

    @testset "the cut follows the projection" begin
        grid = DGG.levelgrid(SYS, 3)
        # Shifting the central meridian shifts the seam with it: with lon_0 = 90
        # the map runs from -90 to 270.
        mesh = DGGV.tessellate(DGGV.PlanarTarget(identity, 270.0), grid)
        @test all(-90.0 - 1.0e-9 .<= lons(mesh) .<= 270.0 + 1.0e-9)

        # And a real projection publishes it.
        transform = Proj.Transformation("+proj=longlat +datum=WGS84", "+proj=moll +lon_0=150"; always_xy = true)
        target = DGGV.plot_target(transform)
        @test target isa DGGV.PlanarTarget
        @test target.cut ≈ -30.0 || target.cut ≈ 330.0
    end

    @testset "wrap = false leaves the ring alone" begin
        grid = DGG.levelgrid(SYS, 2)
        target = DGGV.uncut(planar())
        @test !DGGV.needs_cutting(target)
        mesh = DGGV.tessellate(target, grid)
        # Nothing is split, so ring and cell counts agree.
        @test DGGV.nrings(mesh) == DGG.ncells(grid)
    end

    @testset "projection is applied to the whole buffer" begin
        cells = patch(7)
        transform = Proj.Transformation("+proj=longlat +datum=WGS84", "+proj=merc"; always_xy = true)
        mesh = DGGV.tessellate(DGGV.plot_target(transform), cells)
        plain = DGGV.tessellate(planar(), cells)
        @test length(mesh.positions) == length(plain.positions)
        for i in eachindex(mesh.positions)
            expected = transform(Tuple(plain.positions[i]))
            @test mesh.positions[i][1] ≈ expected[1] atol = 1.0e-6
            @test mesh.positions[i][2] ≈ expected[2] atol = 1.0e-6
        end
    end

    @testset "globe target" begin
        transform = GeoMakie.create_globe_transform(GeoMakie.Geodesy.wgs84_ellipsoid, "+proj=longlat +datum=WGS84", 0.0)
        target = DGGV.plot_target(transform)
        @test target isa DGGV.GlobeTarget
        @test target.a ≈ 6378137.0 rtol = 1.0e-9
        @test target.e2 ≈ 0.00669437999014 rtol = 1.0e-6

        cells = patch(7)
        mesh = DGGV.tessellate(target, cells)
        @test eltype(mesh.positions) === Point3d
        # A globe has no seam, so no cell is ever split.
        @test DGGV.nrings(mesh) == length(cells)
        @test length(mesh.positions) == 6 * length(cells)

        # And the shortcut reproduces what the axis's own transform would do.
        grid = DGG.levelgrid(SYS, DGG.level(cells))
        for (k, cell) in enumerate(cells)
            k > 20 && break
            corner = first(DGG.cell_boundary(grid, cell))
            reference = transform(Makie.Point2d(atand(corner[2], corner[1]), asind(corner[3])))
            ours = mesh.positions[mesh.ring_start[k]]
            @test all(isapprox.(ours, reference; rtol = 1.0e-9))
        end
    end

    @testset "colour follows the cells" begin
        cells = patch(7)
        mesh = DGGV.tessellate(planar(), cells)
        values = Float64.(1:length(cells))

        vertex = DGGV.vertex_colors(mesh, values)
        @test length(vertex) == length(mesh.positions)
        @test vertex == values[mesh.vertex_cell]

        ring = DGGV.ring_colors(mesh, values)
        @test length(ring) == DGGV.nrings(mesh)
        @test ring == values

        # A single colour passes straight through.
        @test DGGV.vertex_colors(mesh, :red) === :red
        @test DGGV.ring_colors(mesh, :red) === :red

        # A wrong-length vector is a mistake worth naming.
        @test_throws ArgumentError DGGV.vertex_colors(mesh, values[1:3])
        @test_throws ArgumentError DGGV.ring_colors(mesh, values[1:3])
    end

    @testset "outlines close every ring" begin
        cells = patch(7)
        mesh = DGGV.tessellate(planar(), cells)
        outline = DGGV.outline_points(mesh)
        @test length(outline) == length(mesh.positions) + 2 * DGGV.nrings(mesh)
        @test count(p -> isnan(p[1]), outline) == DGGV.nrings(mesh)
        # Each loop returns to where it started.
        @test outline[1] == outline[7]
        @test isnan(outline[8][1])
    end

    @testset "plots" begin
        cells = patch(8)
        values = Float64.(1:length(cells))

        # CairoMakie is a vector renderer, so it gets paths, not triangles.
        figure, axis, plot = dggpoly(cells; color = values, strokewidth = 1)
        @test any(p -> p isa Makie.Poly, plot.plots)
        @test saves(figure)

        # The mesh path can be asked for explicitly.
        figure, axis, plot = dggpoly(cells; color = values, primitive = :mesh)
        @test any(p -> p isa Makie.Mesh, plot.plots)
        @test saves(figure)

        # A GeoAxis projects, a GlobeAxis does not.
        figure = Figure()
        geo = GeoAxis(figure[1, 1]; dest = "+proj=moll")
        plot = dggpoly!(geo, DGG.levelgrid(SYS, 2); color = 1:DGG.ncells(DGG.levelgrid(SYS, 2)))
        @test eltype(plot.cellmesh[].positions) === Point2d
        @test saves(figure)

        figure = Figure()
        globe = GlobeAxis(figure[1, 1])
        plot = dggpoly!(globe, DGG.levelgrid(SYS, 2); color = 1:DGG.ncells(DGG.levelgrid(SYS, 2)))
        @test eltype(plot.cellmesh[].positions) === Point3d
        @test saves(figure)
    end


    # ------------------------------------------------------------------
    # dggsurface
    # ------------------------------------------------------------------

    @testset "cellregion" begin
        cells = patch(7)
        # Every container that names cells *and* their adjacency reads the same.
        for x in (cells, DGG.PartialGrid(cells), DGG.CellLookup(cells))
            @test length(DGGV.cellregion(x)) == length(cells)
        end
        grid = DGG.levelgrid(SYS, 2)
        @test length(DGGV.cellregion(grid)) == DGG.ncells(grid)

        # A surface needs topology, and these two do not carry it.
        multiorder = DGG.query(SYS, DGG.MultiOrderCoverage(ALPS); level = 7)
        @test_throws ArgumentError DGGV.cellregion(multiorder)
        @test_throws ArgumentError DGGV.cellregion([DGG.cellindex(grid, 1)])
    end

    @testset "a set is its cells, not the object holding them" begin
        cells = patch(7)
        # Two sets built separately over the same cells: different objects,
        # equal contents.  `===` says no, which is Julia's default for a struct
        # over an array, and is why the guards ask with `==`.
        a, b = DGGV.cellset(cells), DGGV.cellset(DGG.CellVector(DGG.query(SYS,
            DGG.MultiOrderCoverage(ALPS); level = 7)))
        @test !(a === b)
        @test a == b
        @test hash(a) == hash(b)
        @test DGGV.cellregion(cells) == DGGV.cellregion(cells)

        # A grid answers in one comparison rather than in `ncells` of them.
        grid = DGG.levelgrid(SYS, 3)
        @test DGGV.cellset(grid) == DGGV.cellset(grid)
        @test DGGV.GridCells(grid) == DGGV.GridCells(grid)

        # And a different set is still a different set.
        @test DGGV.cellset(patch(6)) != a
        @test DGGV.cellset(DGG.levelgrid(SYS, 4)) != DGGV.cellset(grid)
    end

    @testset "the dual of a whole level" begin
        # Three cells meet at every corner of an IGeo7 grid, so its dual is a
        # triangulation with `V - E + F == 2`, `3V == 2E` and one face per
        # corner: exactly `2 * ncells - 4` triangles, whatever the level.
        for level in 1:4
            grid = DGG.levelgrid(SYS, level)
            mesh = DGGV.triangulate(globe(), grid)
            @test mesh.ncells == DGG.ncells(grid)
            @test mesh.nsplit == 0
            # One vertex per cell and nothing else: a globe has no seam.
            @test length(mesh.positions) == DGG.ncells(grid)
            # Nothing is cut, so the cells are the vertices and a value vector
            # passes through untouched.
            @test DGGV.spread(1:DGG.ncells(grid), mesh) === 1:DGG.ncells(grid)
            @test DGGV.ntriangles(mesh) == 2 * DGG.ncells(grid) - 4
            # No corner is drawn twice.
            @test length(triangle_set(mesh)) == DGGV.ntriangles(mesh)
        end
    end

    @testset "the dual tiles the sphere" begin
        # The real test of the anti-duplication rule, and the one that would
        # catch a gap as readily as an overlap: the triangles' signed areas add
        # up to the whole sphere, so none is missing, none is doubled, and all
        # of them wind the same way.
        for system in DGG.systems(), level in 2:3
            grid = DGG.levelgrid(system, level)
            DGG.ncells(grid) > 60_000 && continue
            mesh = DGGV.triangulate(globe(), grid)
            areas = [spherical_area(t...) for t in corners(mesh)]
            @test sum(areas) ≈ 4π rtol = 1.0e-9
            @test all(>(0), areas)
            @test length(triangle_set(mesh)) == DGGV.ntriangles(mesh)
        end
    end

    @testset "the dual tiles the map" begin
        # The same statement in longitude/latitude, where it also says that the
        # cut was split cleanly and that both polar caps were filled: the
        # triangles cover the whole 360 by 180 rectangle exactly once.
        for system in DGG.systems(), level in 2:3, cut in (180.0, 0.0, -75.0)
            grid = DGG.levelgrid(system, level)
            DGG.ncells(grid) > 60_000 && continue
            mesh = DGGV.triangulate(DGGV.PlanarTarget(identity, cut), grid)
            areas = [plane_area(t...) for t in corners(mesh)]
            @test sum(areas) ≈ 360 * 180 rtol = 1.0e-9
            # A clip can leave a triangle with no area at all; none may have
            # negative area, which would mean it was drawn inside out.
            @test all(>=(-1.0e-9), areas)
            @test all(p -> -90 <= p[2] <= 90, mesh.positions)
            @test all(p -> cut - 360 - 1.0e-9 <= p[1] <= cut + 1.0e-9, mesh.positions)
        end
    end

    @testset "four cells to a corner" begin
        # A quadrilateral grid puts four cells around each corner, and its dual
        # face is therefore a quadrilateral to be split in two — not the four
        # overlapping triangles that taking every adjacent pair would give.
        # Both readings are duplicate-free, so only the count tells them apart:
        # about two triangles per cell rather than about three.
        for system in (DGG.HEALPixSystem(), DGG.S2System(), DGG.ISEA4RSystem())
            grid = DGG.levelgrid(system, 3)
            mesh = DGGV.triangulate(globe(), grid)
            @test DGGV.ntriangles(mesh) / DGG.ncells(grid) < 2.05
        end
    end

    @testset "a partial grid stops where the data does" begin
        cells = patch(6)
        mesh = DGGV.triangulate(globe(), cells)
        @test mesh.ncells == length(cells)
        @test length(triangle_set(mesh)) == DGGV.ntriangles(mesh)

        # Every triangle of the patch is a triangle of the whole level whose
        # three cells all survived, and every such triangle is there.  So the
        # surface neither invents a triangle across a hole nor drops one at the
        # edge — it is the complete dual, restricted.
        grid = DGG.levelgrid(SYS, DGG.level(cells))
        whole = DGGV.triangulate(globe(), grid)
        index = Dict(DGG.cellindex(grid, p) => p for p in 1:DGG.ncells(grid))
        member = Dict(cells[p] => p for p in 1:length(cells))
        restricted = Set(
            Tuple(sort([member[DGG.cellindex(grid, Int(GeometryBasics.value(f[j])))] for j in 1:3]))
                for f in whole.faces
                if all(DGG.cellindex(grid, Int(GeometryBasics.value(f[j]))) in keys(member) for j in 1:3)
        )
        @test triangle_set(mesh) == restricted
        @test index isa Dict  # the level grid is indexed by array index, as assumed above
    end

    @testset "threading is not observable" begin
        cells = patch(7)
        one = DGGV.triangulate(planar(), cells; ntasks = 1)
        many = DGGV.triangulate(planar(), cells; ntasks = 8)
        @test one.positions == many.positions
        @test triangle_set(one) == triangle_set(many)
        @test one.nsplit == many.nsplit
    end

    @testset "empty input" begin
        cells = DGG.CellVector(SYS, 5, DGG.Z7Cell[])
        mesh = DGGV.triangulate(planar(), cells)
        @test isempty(mesh)
        @test mesh.ncells == 0
        @test DGGV.ntriangles(mesh) == 0
    end

    @testset "wrap = false leaves the seam alone" begin
        grid = DGG.levelgrid(SYS, 2)
        mesh = DGGV.triangulate(DGGV.uncut(planar()), grid)
        # Nothing is cut and nothing is capped, so every triangle is three
        # centroids and the mesh carries no vertex of its own.
        @test mesh.nsplit == 0
        @test length(mesh.positions) == DGG.ncells(grid)
    end

    @testset "surface colour is the cell vector itself" begin
        grid = DGG.levelgrid(SYS, 2)
        values = Float64.(1:DGG.ncells(grid))

        # On a globe there is one vertex per cell, in cell order, so a per-cell
        # vector needs no gather at all.
        mesh = DGGV.triangulate(globe(), grid)
        @test DGGV.vertex_colors(mesh, values) === values

        # On a map the seam and the poles add vertices, which do need one.  A
        # vertex the seam created is a point inside a triangle, and takes that
        # triangle's three values mixed — so it lands between them, never
        # outside.
        cut = DGGV.triangulate(planar(), grid)
        @test cut.nsplit > 0
        vertex = DGGV.vertex_colors(cut, values)
        @test length(vertex) == length(cut.positions)
        @test vertex[1:DGG.ncells(grid)] == values
        for k in 1:length(cut.extra_tri)
            v = vertex[DGG.ncells(grid) + k]
            corners = values[collect(cut.extra_tri[k])]
            @test minimum(corners) - 1.0e-9 <= v <= maximum(corners) + 1.0e-9
        end

        @test DGGV.vertex_colors(cut, :red) === :red
        @test_throws ArgumentError DGGV.vertex_colors(cut, values[1:3])

        # A cube axis is unwrapped: a `DimArray` would still be one after
        # Makie's colour conversion, and the backends will not draw that.
        A = DD.DimArray(values, DGG.Cells(DGG.CellLookup(grid)))
        @test DGGV.vertex_colors(mesh, A) == values
        @test DGGV.vertex_colors(mesh, A) isa Vector{Float64}
    end

    @testset "ntasks follows the session, not the precompile worker" begin
        # A recipe's attribute defaults are evaluated once, where the `@recipe`
        # block is read — during precompilation, in a worker with one thread.
        # A default of `Threads.nthreads()` therefore bakes `1` into the package
        # image and no session ever gets more, which is what this pins.
        cells = patch(6)
        for make in (dggpoly, dggsurface, dggresample)
            _, _, plot = make(cells)
            @test plot.ntasks[] === Makie.automatic
        end
        @test DGGV.task_count(Makie.automatic) == Threads.nthreads()
        @test DGGV.task_count(4) == 4
    end

    @testset "heights are part of the geometry" begin
        grid = DGG.levelgrid(SYS, 2)
        n = DGG.ncells(grid)
        zs = Float64.(1:n)

        # Left out, there is no third coordinate to carry, and the fill that
        # stands in for the heights stores nothing but its length.
        flat = DGGV.triangulate(planar(), grid)
        @test flat.nsplit > 0
        @test eltype(flat.positions) === Point2d
        @test DGGV.triangulate(planar(), grid, DGGV.ZeroHeights(n)).positions ==
            flat.positions
        @test all(iszero, DGGV.ZeroHeights(n))
        @test (@allocated DGGV.ZeroHeights(n)) == 0

        # Given, the map is the same map with a height on top of it.  The seam
        # and the poles add vertices of their own, and each takes the height of
        # the cell it takes its colour from.
        raised = DGGV.triangulate(planar(), grid, zs)
        @test eltype(raised.positions) === Point3d
        @test [Point2d(p[1], p[2]) for p in raised.positions] == flat.positions
        @test [p[3] for p in raised.positions] == DGGV.spread(zs, raised)

        @test_throws ArgumentError DGGV.triangulate(planar(), grid, zs[1:3])
    end

    @testset "a height on a globe is a height above the ellipsoid" begin
        # Not a third coordinate but a lift straight out along the cell's own
        # centroid, which on an oblate ellipsoid is not the direction the
        # projected vertex points in: lifting along that instead lands 3.4 m
        # away in the 1 km below, against the 0.1 m `≈` allows at this size.
        grid = DGG.levelgrid(SYS, 2)
        n = DGG.ncells(grid)
        target = DGGV.GlobeTarget(identity, 6378137.0, 0.00669437999014, 0.0)
        flat = DGGV.triangulate(target, grid)
        raised = DGGV.triangulate(target, grid, fill(1000.0, n))
        @test all(1:n) do p
            u = DGG.cell_centroid(grid, DGG.cellindex(grid, p))
            raised.positions[p] ≈ flat.positions[p] + 1000.0 * Point3d(u...)
        end
    end

    @testset "surface plots" begin
        cells = patch(8)
        values = Float64.(1:length(cells))

        figure, axis, plot = dggsurface(cells; color = values)
        @test any(p -> p isa Makie.Mesh, plot.plots)
        @test plot.surfacemesh[].ncells == length(cells)
        @test saves(figure)

        grid = DGG.levelgrid(SYS, 2)
        figure = Figure()
        geo = GeoAxis(figure[1, 1]; dest = "+proj=moll")
        plot = dggsurface!(geo, grid; color = 1:DGG.ncells(grid))
        @test eltype(plot.surfacemesh[].positions) === Point2d
        @test saves(figure)

        figure = Figure()
        globeaxis = GlobeAxis(figure[1, 1])
        plot = dggsurface!(globeaxis, grid; color = 1:DGG.ncells(grid))
        @test eltype(plot.surfacemesh[].positions) === Point3d
        # A globe has no seam, so the mesh is exactly one vertex per cell.
        @test length(plot.surfacemesh[].positions) == DGG.ncells(grid)
        @test saves(figure)
    end

    @testset "surface heights through the recipe" begin
        cells = patch(7)
        n = length(cells)
        values = Float64.(1:n)

        # No heights given is the flat surface, and a flat surface keeps the
        # two coordinates it had.
        figure, axis, plot = dggsurface(cells; color = values)
        @test plot.zs[] isa DGGV.ZeroHeights
        @test eltype(plot.surfacemesh[].positions) === Point2d

        figure, axis, plot = dggsurface(cells, values; color = values)
        @test [p[3] for p in plot.surfacemesh[].positions] ==
            DGGV.spread(values, plot.surfacemesh[])
        @test saves(figure)

        # A one-dimensional cube axis is both at once: its lookup names the
        # cells and its values are their heights.
        A = DD.DimArray(values, DGG.Cells(DGG.CellLookup(cells)))
        region, zs = Makie.convert_arguments(DGGV.DGGSurface, A)
        @test length(region) == n
        @test zs === values

        figure, axis, plot = dggsurface(A; color = A)
        @test plot.surfacemesh[].ncells == n
        @test [p[3] for p in plot.surfacemesh[].positions] ==
            DGGV.spread(values, plot.surfacemesh[])
        @test saves(figure)

        # A dimension that is not cells names no surface.
        @test_throws ArgumentError dggsurface(DD.DimArray(values, DD.X(1:n)))
    end

    @testset "recolouring a surface keeps the geometry" begin
        cells = patch(7)
        figure, axis, plot = dggsurface(cells; color = Float64.(1:length(cells)))
        before = plot.surfacemesh[]
        plot.color = Float64.(length(cells):-1:1)
        @test plot.surfacemesh[] === before
        @test plot.mesh_color[] == DGGV.spread(Float64.(length(cells):-1:1), before)
    end

    @testset "a recolour and a re-raise write over their own buffers" begin
        # A whole level crosses the map's cut and both poles, so the surface has
        # vertices the cells did not give it and the colours need a buffer.
        cells = DGG.CellVector(DGG.levelgrid(SYS, 2))
        n = length(cells)
        values = Float64.(1:n)

        figure, axis, plot = dggsurface(cells, values; color = values)
        positions, colors = plot.mesh_positions[], plot.mesh_color[]
        @test length(colors) > n

        plot.color = reverse(values)
        @test plot.mesh_color[] === colors
        @test plot.mesh_color[] == DGGV.spread(reverse(values), plot.surfacemesh[])

        Makie.update!(plot, cells, values .* 2)
        @test plot.mesh_positions[] === positions
        @test [p[3] for p in plot.mesh_positions[]] ==
            DGGV.spread(2 .* values, plot.surfacemesh[])

        # What a plot does not own it does not write over.  A flat map at no
        # height draws the topology's own vertices, and a surface with nothing
        # cut draws the caller's own colours.
        figure, axis, flat = dggsurface(cells; color = values)
        @test flat.mesh_positions[] === flat.surfacetopology[].positions

        figure = Figure(size = (400, 300))
        onglobe = dggsurface!(GlobeAxis(figure[1, 1]), cells; color = values)
        @test onglobe.mesh_color[] === values
    end

    @testset "a projection set after the plot re-projects rather than rebuilds" begin
        # What a tutorial does: plot first, then hand the axis its projection.
        # The topology has to survive that — the transform function is no part
        # of its type — and only the vertices may move.
        cells = patch(6)
        values = Float64.(1:length(cells))
        ortho(lat) = GeoMakie.create_transform(
            "+proj=ortho +lon_0=10.5 +lat_0=$lat +datum=WGS84",
            "+proj=longlat +datum=WGS84")

        figure, axis, plot = dggsurface(cells; color = values)
        flat = copy(plot.mesh_positions[])
        plot.transformation.transform_func[] = ortho(46.5)
        Makie.update_state_before_display!(figure)

        built = plot.surfacetopology[]
        @test plot.mesh_positions[] != flat
        @test plot.mesh_positions[] == plot.projectedtopology[].positions
        # The build keeps its longitudes and latitudes, or there would be
        # nothing left to project the second time.
        @test built.positions == flat

        # Another projection about the same central meridian has the same seam,
        # so the expensive half is reused and the vertices are projected again
        # into the buffer that node already owns.
        projected = plot.mesh_positions[]
        plot.transformation.transform_func[] = ortho(0.0)
        Makie.update_state_before_display!(figure)
        @test plot.surfacetopology[] === built
        @test plot.mesh_positions[] === projected
        @test built.positions == flat
        @test saves(figure)

        # Heights are read after the projection, not before it.
        figure, axis, raised = dggsurface(cells, values)
        raised.transformation.transform_func[] = ortho(46.5)
        Makie.update_state_before_display!(figure)
        @test [p[3] for p in raised.mesh_positions[]] ==
            DGGV.spread(values, raised.surfacemesh[])
    end

    @testset "a surface can be coloured by its own heights" begin
        cells = patch(7)
        values = Float64.(1:length(cells))

        figure, axis, plot = dggsurface(cells, values; color = nothing)
        @test plot.mesh_color[] == DGGV.spread(values, plot.surfacemesh[])
        @test saves(figure)

        # Explicit colour still wins, and so does the default, which is a
        # colour rather than a field.
        figure, axis, plot = dggsurface(cells, values)
        @test plot.mesh_color[] != DGGV.spread(values, plot.surfacemesh[])
    end

    @testset "a cube axis names cells for every recipe" begin
        cells = patch(7)
        n = length(cells)
        values = Float64.(1:n)
        A = DD.DimArray(values, DGG.Cells(DGG.CellLookup(cells)))

        # `dggpoly` and `dggresample` read only the cells: a cube axis is a
        # cell set for them, and its values are a colour like any other.
        @test length(DGGV.cellset(A)) == n
        @test length(only(Makie.convert_arguments(DGGV.DGGPoly, A))) == n
        @test length(first(Makie.convert_arguments(DGGV.DGGResample, A))) == n

        figure, axis, plot = dggpoly(A; color = A)
        @test plot.cellmesh[].ncells == n
        @test saves(figure)
    end

    @testset "a mesh to hand somewhere else" begin
        cells = patch(7)
        n = length(cells)
        values = Float64.(1:n)
        region = DGGV.cellregion(cells)

        # The default space is the unit sphere, and a height is a height above
        # it — the same reading the recipe gives a globe.
        flat = GeometryBasics.mesh(region)
        @test length(GeometryBasics.coordinates(flat)) == n
        @test all(p -> sqrt(sum(abs2, p)) ≈ 1, GeometryBasics.coordinates(flat))

        raised = GeometryBasics.mesh(region, values ./ 100)
        @test [sqrt(sum(abs2, p)) for p in GeometryBasics.coordinates(raised)] ≈
            1 .+ values ./ 100
        @test GeometryBasics.faces(raised) == DGGV.triangulate(globe(), cells).faces

        # A value per cell becomes a value per vertex, including where the map's
        # cut has made more vertices than there are cells.
        whole = DGGV.cellregion(DGG.levelgrid(SYS, 2))
        heights = Float64.(1:length(whole))
        cut = GeometryBasics.mesh(whole, heights; target = planar(), color = heights)
        @test length(GeometryBasics.coordinates(cut)) > length(whole)
        @test cut.color ==
            DGGV.spread(heights, DGGV.triangulate(planar(), whole, heights))

        # A cube axis is cells and heights at once, but it takes both verbs to
        # say so — `mesh` is `GeometryBasics`', and only our own types reach it.
        A = DD.DimArray(values ./ 100, DGG.Cells(DGG.CellLookup(cells)))
        @test GeometryBasics.coordinates(GeometryBasics.mesh(DGGV.cellregion(A), A)) ==
            GeometryBasics.coordinates(raised)

        # The patch mesh is the other one the package draws, and carries no
        # heights.
        patches = GeometryBasics.mesh(DGGV.cellset(cells); color = values)
        @test length(patches.color) == length(GeometryBasics.coordinates(patches))
        @test patches.color[1] == values[1]
        @test_throws ArgumentError GeometryBasics.mesh(DGGV.cellset(cells), values)

        # A set spanning levels has no surface, but it does have patches.
        spanning = DGG.MultiOrderCellSet(SYS, DGG.MultiOrderCoverage(ALPS); level = 7)
        @test GeometryBasics.mesh(DGGV.cellset(spanning)) isa GeometryBasics.Mesh
    end

    # ------------------------------------------------------------------
    # dggresample
    # ------------------------------------------------------------------

    "Bring `figure` up to date and hand back the frame `plot` settled on."
    function frame(figure, plot)
        Makie.update_state_before_display!(figure)
        return plot.resampled[]
    end

    "Set an axis's limits to a square `width` degrees across, centred on the Alps."
    function look!(axis, width)
        xlims!(axis, 10.5 - width / 2, 10.5 + width / 2)
        ylims!(axis, 46.5 - width / 2, 46.5 + width / 2)
        return axis
    end

    "The arithmetic mean.  A reduction is handed to the plot, so it need not be one the package knows."
    mean(v) = sum(v) / length(v)

    "The middle value of `v`, or the lower of the two middle ones."
    median(v) = partialsort(collect(v), (length(v) + 1) ÷ 2)

    @testset "a pyramid over a cell set" begin
        cells = patch(9)
        pyramid = DGGV.CellPyramid(DGGV.cellset(cells))

        @test pyramid.leaflevel == 9
        @test pyramid.rootlevel == first(DGG.levels(SYS))
        @test pyramid.ncells == length(cells)

        # The cap covers the data it was sampled from.
        for c in cells[1:97:end]
            @test DGGV._angle(pyramid.capcentre, DGG.cell_centroid(SYS, c)) <= pyramid.capradius
        end

        # A cell of the set finds itself; a cell on the far side of the world
        # finds nothing, which is what makes a coarse view of a partial grid
        # show the grid rather than its bounding box.
        @test DGGV.nearest(pyramid, cells[1]) == 1
        @test DGGV.nearest(pyramid, cells[end]) == length(cells)
        antipode = DGG.cellat(DGG.levelgrid(SYS, 9), -169.5, -46.5)
        @test DGGV.nearest(pyramid, antipode) == 0

        # A coarse cell over the patch resolves to one of the leaves under it.
        coarse = DGG.ancestor(SYS, cells[end ÷ 2], 6)
        @test DGGV.nearest(pyramid, coarse) in 1:length(cells)

        @test_throws ArgumentError DGGV.CellPyramid(DGGV.CellSet(SYS, empty(cells)))
    end

    @testset "locating a cell in its set" begin
        cells = patch(7)
        # Every way of naming a set can say where one of its cells sits.
        @test DGGV.celllocator(cells)(cells[3]) == 3
        @test DGGV.celllocator(collect(cells))(cells[3]) == 3
        @test DGGV.celllocator(collect(cells))(DGG.cellat(DGG.levelgrid(SYS, 7), -169.5, -46.5)) == 0

        grid = DGG.levelgrid(SYS, 2)
        gridcells = DGGV.GridCells(grid)
        @test DGGV.celllocator(gridcells)(DGG.cellindex(grid, 5)) == 5
    end

    @testset "the rim of a cap" begin
        p = DGG.UnitSphericalPoint(cosd(20) * cosd(35), cosd(20) * sind(35), sind(20))
        a, b = DGGV.rimpoints(p, 0.3)
        @test DGGV._angle(p, a) ≈ 0.3
        @test DGGV._angle(p, b) ≈ 0.3
        # A quarter turn apart around the cap, so between them they see both
        # of the ways a projection can stretch a circle.
        @test DGGV._angle(a, b) ≈ acos(cos(0.3)^2) atol = 1e-12
    end

    @testset "the level follows the zoom" begin
        cells = patch(11)
        pyramid = DGGV.CellPyramid(DGGV.cellset(cells))

        figure = Figure(size = (600, 400))
        axis = Axis(figure[1, 1])
        dggpoly!(axis, patch(4); color = :red)   # something to fix the limits
        Makie.update_state_before_display!(figure)

        levels = Int[]
        counts = Int[]
        for width in (2.0, 0.5, 0.125)
            look!(axis, width)
            Makie.update_state_before_display!(figure)
            view = DGGV.ScreenView(planar(), axis.scene, 1.6)
            drawn, level = DGGV.resample(pyramid, view)
            push!(levels, level)
            push!(counts, length(drawn))
        end

        # Zooming in shows a finer level, and never past the level stored.
        @test issorted(levels)
        @test levels[1] < levels[end] <= pyramid.leaflevel
        # And the number of cells on screen stays of a size with the screen,
        # rather than growing with the zoom.
        @test all(n -> n < 200_000, counts)

        # How many threads did the descent is not observable in its answer.
        look!(axis, 0.5)
        Makie.update_state_before_display!(figure)
        view = DGGV.ScreenView(planar(), axis.scene, 1.6)
        @test DGGV.resample(pyramid, view; ntasks = 1) ==
            DGGV.resample(pyramid, view; ntasks = 4)
    end

    @testset "a resampled plot" begin
        cells = patch(11)
        values = Float64.(1:length(cells))

        figure = Figure(size = (600, 400))
        axis = Axis(figure[1, 1])
        plot = dggresample!(axis, cells; color = values)
        wide = frame(figure, plot)

        # Far less is drawn than was handed in, and every drawn cell carries a
        # value from the set.
        @test 0 < length(wide) < length(cells) ÷ 5
        @test wide.level <= 11
        @test length(wide.cells.cells) == length(wide.index)
        @test all(i -> 1 <= i <= length(cells), wide.index)
        @test plot.cellcolor[] == values[wide.index]
        @test saves(figure)

        # Zooming in refines without drawing more than a screen's worth.
        look!(axis, 0.05)
        close = frame(figure, plot)
        @test close.level > wide.level
        @test length(close) < 200_000
        @test saves(figure)
    end

    @testset "the limits do not follow the resampling" begin
        # The camera drives the resampling, so if the resampling drove the
        # limits the plot would chase itself.  What the axis is told is the
        # whole data set, whatever is on screen.
        figure = Figure(size = (600, 400))
        axis = Axis(figure[1, 1])
        plot = dggresample!(axis, patch(11); color = :grey)
        wide = frame(figure, plot)
        before = plot.datalimits[]

        look!(axis, 0.05)
        close = frame(figure, plot)
        @test close !== wide
        @test plot.datalimits[] == before
    end

    @testset "a build covers more than the viewport" begin
        figure = Figure(size = (600, 400))
        axis = Axis(figure[1, 1])
        plot = dggresample!(axis, patch(10); color = :grey, buffer = 3.0)
        look!(axis, 0.4)
        built = frame(figure, plot)

        # A pan well inside the buffer reuses what is already there.
        xlims!(axis, 10.5 - 0.2 + 0.02, 10.5 + 0.2 + 0.02)
        @test frame(figure, plot) === built

        # Leaving it does not.
        look!(axis, 0.02)
        @test frame(figure, plot) !== built
    end

    @testset "dynamic = false builds once" begin
        figure = Figure(size = (600, 400))
        axis = Axis(figure[1, 1])
        plot = dggresample!(axis, patch(10); color = :grey, dynamic = false)
        built = frame(figure, plot)
        look!(axis, 0.02)
        @test frame(figure, plot) === built
    end

    @testset "resampling in every axis" begin
        cells = DGG.CellVector(DGG.levelgrid(SYS, 5))
        values = Float64.(1:length(cells))

        figure = Figure(size = (400, 300))
        geo = GeoAxis(figure[1, 1]; dest = "+proj=moll")
        plot = dggresample!(geo, cells; color = values)
        @test 0 < length(frame(figure, plot)) <= length(cells)
        @test saves(figure)

        figure = Figure(size = (400, 300))
        globe = GlobeAxis(figure[1, 1])
        plot = dggresample!(globe, cells; color = values)
        drawn = frame(figure, plot)
        # A globe only ever shows one side of itself.
        @test 0 < length(drawn) < length(cells)
        @test saves(figure)
    end

    @testset "recolouring a resampled plot keeps the frame" begin
        cells = patch(10)
        figure = Figure(size = (600, 400))
        axis = Axis(figure[1, 1])
        plot = dggresample!(axis, cells; color = Float64.(1:length(cells)))
        built = frame(figure, plot)

        plot.color = Float64.(length(cells):-1:1)
        @test plot.resampled[] === built
        @test plot.cellcolor[] == Float64.(length(cells):-1:1)[built.index]
    end

    @testset "a resampled plot can carry heights" begin
        cells = patch(9)
        n = length(cells)
        values = Float64.(1:n)

        # Without heights the frame is drawn as patches, as it always was.
        figure, axis, plot = dggresample(cells; color = values)
        @test plot.zs[] === nothing
        @test any(p -> p isa DGGV.DGGPoly, plot.plots)

        figure, axis, plot = dggresample(cells, values; color = values)
        @test any(p -> p isa DGGV.DGGSurface, plot.plots)
        drawn = frame(figure, plot)
        region, index = DGGV.surfaceframe(drawn)

        # Heights and colours are both indexed by the cells handed in, so both
        # are the same gather over the frame.
        @test length(region) == length(drawn)
        @test plot.cellheights[] == values[index]
        @test plot.cellcolor[] == values[index]

        # The frame is one level, which is what lets it have a surface at all.
        @test DGG.level(region.region) == drawn.level
        @test saves(figure)


        # Re-raising is the same gather again: the cells are the cells the
        # frame was built from, so the frame stands and only the heights move.
        Makie.update!(plot, cells, values .* 2)
        @test frame(figure, plot) === drawn
        @test plot.cellheights[] == 2 .* values[index]

        # A source and a vector of ids is still a source and its ids.
        ids = collect(cells)
        @test only(Makie.convert_arguments(DGGV.DGGResample, SYS, ids)[2:2]) === nothing

        # A globe raises above the ellipsoid here as it does anywhere else.
        globefigure = Figure(size = (400, 300))
        globe = GlobeAxis(globefigure[1, 1])
        onglobe = dggresample!(globe, cells, values ./ 50; color = values)
        @test length(frame(globefigure, onglobe)) > 0
        @test saves(globefigure)

        # Looking away from the data leaves a frame with no cells in it, and a
        # surface over no cells is a surface with no triangles, not an error.
        away, awayaxis, awayplot = dggresample(cells, values; color = values)
        xlims!(awayaxis, -170, -160)
        ylims!(awayaxis, -80, -70)
        @test length(frame(away, awayplot)) == 0
        @test length(awayplot.displayed[]) == 0
        @test saves(away)
    end

    @testset "the leaves under a drawn cell" begin
        cells = patch(9)
        pyramid = DGGV.CellPyramid(DGGV.cellset(cells))
        leaves = DGGV.subtreeranges(pyramid)

        # Every level-7 ancestor of the patch claims a run of it; the runs
        # ascend, do not overlap, and between them are the whole patch.  That
        # is the property that lets a range stand in for a list of leaves.
        coarse = unique(DGG.ancestor(SYS, c, 7) for c in cells)
        runs = [leaves(c) for c in coarse]
        @test issorted(runs; by = first)
        @test all(i -> last(runs[i]) < first(runs[i + 1]), 1:length(runs) - 1)
        @test sum(length, runs) == length(cells)
        # And a run holds the cells under its own ancestor, not its neighbour's.
        @test all(k -> DGG.ancestor(SYS, cells[k], 7) == coarse[1], runs[1])

        # A cell of the leaf level is its own one leaf, and a cell over nothing
        # gets an empty run rather than a wrong one.
        @test leaves(cells[5]) == 5:5
        away = DGG.ancestor(SYS, DGG.cellat(DGG.levelgrid(SYS, 9), -169.5, -46.5), 7)
        @test isempty(leaves(away))
    end

    @testset "a summarised frame reads every leaf" begin
        cells = patch(11)
        values = Float64.(1:length(cells))

        figure = Figure(size = (400, 300))
        axis = Axis(figure[1, 1])
        plot = dggresample!(axis, cells; color = values, aggregate = mean)
        drawn = frame(figure, plot)

        # The frame is far coarser than the data, and the runs partition it:
        # every value handed in is read by exactly one drawn cell.
        @test 0 < length(drawn) < length(cells) ÷ 10
        @test sum(length, drawn.index) == length(cells)
        @test plot.cellcolor[] == [mean(values[g]) for g in drawn.index]
        @test saves(figure)

        # The reduction is the API, so anything of the shape works — including
        # two the package has never heard of.
        for f in (sum, maximum, minimum, mean, median)
            plot.aggregate = f
            @test plot.cellcolor[] == [f(values[g]) for g in drawn.index]
        end

        # A height is indexed by the cells handed in exactly as a colour is, so
        # it is summarised by the same call over the same runs.
        raised, _, surface = dggresample(cells, values; color = values, aggregate = mean)
        lifted = frame(raised, surface)
        _, index = DGGV.surfaceframe(lifted)
        @test surface.cellheights[] == [mean(values[g]) for g in index]
        @test surface.cellcolor[] == surface.cellheights[]
        @test saves(raised)
    end

    @testset "one reduction in place of another re-reads the frame" begin
        cells = patch(10)
        values = Float64.(1:length(cells))
        figure = Figure(size = (400, 300))
        axis = Axis(figure[1, 1])
        plot = dggresample!(axis, cells; color = values, aggregate = mean)
        built = frame(figure, plot)

        # Which values lie under a drawn cell does not depend on what they are,
        # so another reduction re-reads the frame rather than descending again.
        plot.aggregate = maximum
        @test plot.resampled[] === built
        @test plot.cellcolor[] == [maximum(values[g]) for g in built.index]

        # Nearest neighbour is a different grouping, though, and that one is
        # built.
        plot.aggregate = nothing
        nearest = frame(figure, plot)
        @test nearest !== built
        @test eltype(nearest.index) == Int32
        @test plot.cellcolor[] == values[nearest.index]
    end

    @testset "aggregation needs the leaves to sit together" begin
        cells = patch(8)
        values = Float64.(1:length(cells))
        ids = collect(cells)

        # A bare list of ids promises no order, so there is no run to reduce
        # over — and saying so is better than a slower answer or a wrong one.
        # Nearest neighbour still has a leaf under every centre.
        figure = Figure(size = (400, 300))
        @test_throws ArgumentError dggresample!(Axis(figure[1, 1]), SYS, ids;
            color = values, aggregate = mean)
        sampled = dggresample!(Axis(figure[2, 1]), SYS, ids; color = values)
        @test length(frame(figure, sampled)) > 0

        # A5 scatters a subtree through its level, which is the same refusal
        # for the other of the two reasons.
        a5 = DGG.A5System()
        a5cells = DGG.CellVector(DGG.query(a5, DGG.MultiOrderCoverage(ALPS); level = 5))
        a5values = Float64.(1:length(a5cells))
        other = Figure(size = (400, 300))
        @test_throws ArgumentError dggresample!(Axis(other[1, 1]), a5cells;
            color = a5values, aggregate = mean)
        onA5 = dggresample!(Axis(other[2, 1]), a5cells; color = a5values)
        @test length(frame(other, onA5)) > 0
    end

    @testset "draw chooses the plot, not the arguments" begin
        cells = patch(9)
        values = Float64.(1:length(cells))

        # A flat field drawn as a surface: no heights, so `ZeroHeights` stands
        # in for them, which is what keeps the vertex buffer two-dimensional.
        figure, axis, plot = dggresample(cells; color = values, draw = :surface)
        drawn = frame(figure, plot)
        @test any(p -> p isa DGGV.DGGSurface, plot.plots)
        @test plot.cellheights[] isa DGGV.ZeroHeights
        @test length(plot.cellheights[]) == length(drawn) > 0
        @test saves(figure)

        # Heights as patches is the one pairing with no picture behind it, and
        # a name that is neither is not quietly ignored.
        @test_throws ArgumentError dggresample(cells, values; color = values, draw = :patches)
        @test_throws ArgumentError dggresample(cells; color = values, draw = :blobs)

        # Naming what would have been chosen anyway changes nothing.
        figure, axis, plot = dggresample(cells; color = values, draw = :patches)
        @test any(p -> p isa DGGV.DGGPoly, plot.plots)
    end

    @testset "recolouring keeps the geometry" begin
        cells = patch(8)
        figure, axis, plot = dggpoly(cells; color = Float64.(1:length(cells)), primitive = :mesh)
        before = plot.cellmesh[]
        plot.color = Float64.(length(cells):-1:1)
        @test plot.cellmesh[] === before
        @test plot.mesh_color[] == Float64.(length(cells):-1:1)[before.vertex_cell]
    end
end
