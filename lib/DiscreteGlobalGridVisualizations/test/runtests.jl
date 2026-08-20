using Test
import DiscreteGlobalGrids as DGG
import DiscreteGlobalGridVisualizations as DGGV
using DiscreteGlobalGridVisualizations
using GeometryBasics: Point2d, Point3d
using Makie
using CairoMakie
using GeoMakie
import Proj

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

"Save `figure` to a temporary file and say whether it landed."
function saves(figure)
    path = tempname() * ".png"
    Makie.save(path, figure)
    return isfile(path)
end

@testset "DiscreteGlobalGridVisualizations" begin

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

    @testset "recolouring keeps the geometry" begin
        cells = patch(8)
        figure, axis, plot = dggpoly(cells; color = Float64.(1:length(cells)), primitive = :mesh)
        before = plot.cellmesh[]
        plot.color = Float64.(length(cells):-1:1)
        @test plot.cellmesh[] === before
        @test plot.mesh_color[] == Float64.(length(cells):-1:1)[before.vertex_cell]
    end
end
