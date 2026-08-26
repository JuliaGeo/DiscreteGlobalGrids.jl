isdefined(@__MODULE__, :DGGSTalkFigures) ||
    include(joinpath(@__DIR__, "00-dggs-theme.jl"))

using .DGGSTalkFigures
using DiscreteGlobalGrids, DiscreteGlobalGridsVisualization, Makie
using GeometryBasics
using LinearAlgebra

const CHUNK_LEVEL, LOCAL_CURVE_LEVEL = 5, 9
const GLOBAL_CURVE_LEVEL = 3

function z7_curve_fixture()
    system = IGeo7System()
    chunks = CellVector(levelgrid(system, CHUNK_LEVEL))
    chunk = chunks[286]
    chunk_grid = levelgrid(system, CHUNK_LEVEL)
    chunk_center = Point3f(cell_centroid(chunk_grid, chunk)...)
    local_grid = levelgrid(system, LOCAL_CURVE_LEVEL)
    local_interval = descendant_range(system, chunk, LOCAL_CURVE_LEVEL)
    local_cells = CellVector(local_grid)[collect(local_interval)]
    footprint_cells = collect(border(
        subtree(system, chunk, LOCAL_CURVE_LEVEL); cells = true))
    global_grid = levelgrid(system, GLOBAL_CURVE_LEVEL)
    global_cells = collect(CellVector(global_grid))
    return (; system, chunk, chunk_center, local_grid, local_cells,
        footprint_cells,
        global_grid, global_cells)
end

centroid_points(grid, cells) = Point3f.(cell_centroid.(Ref(grid), cells))

function dense_spherical_path(points; angular_step = deg2rad(0.20))
    dense = Point3f[]
    for i in 1:(length(points) - 1)
        a, b = points[i], points[i + 1]
        angle = acos(clamp(dot(a, b), -1, 1))
        samples = max(2, ceil(Int, angle / angular_step))
        for k in 0:(samples - 1)
            t = k / samples
            # Normalized interpolation traces the short spherical arc while
            # avoiding the below-surface chords produced by sparse 3D lines.
            push!(dense, Point3f(normalize((1 - t) * a + t * b)))
        end
    end
    push!(dense, last(points))
    dense
end

function tangent_points(center, points)
    east = normalize(Point3f(-center[2], center[1], 0))
    north = normalize(cross(center, east))
    Point2f[Point2f(dot(p, east), dot(p, north)) for p in points]
end

function z7_curve_figure()
    data = z7_curve_fixture()
    local_centroids = centroid_points(data.local_grid, data.local_cells)
    local_path = local_centroids
    local_polygons = [GeometryBasics.Polygon(tangent_points(data.chunk_center,
        Point3f.(cell_boundary(data.system, cell)))) for cell in data.local_cells]
    footprint_path = tangent_points(data.chunk_center,
        centroid_points(data.local_grid, data.footprint_cells))
    sort!(footprint_path; by = p -> atan(p[2], p[1]))
    push!(footprint_path, first(footprint_path))
    global_centroids = centroid_points(data.global_grid, data.global_cells)
    global_path = dense_spherical_path(global_centroids)

    with_theme(JG_THEME) do
        fig, body = slide_figure()

        left = GridLayout(body[1, 1])
        Label(left[1, 1], "ONE LEVEL-5 CHUNK · LEVEL-9 Z7 WALK";
            fontsize = 11, font = FONT_BODY, color = JG.caption,
            tellwidth = false)
        local_axis = Axis(left[2, 1]; aspect = DataAspect(),
            backgroundcolor = JG.paper)
        hidedecorations!(local_axis); hidespines!(local_axis)
        local_xy = tangent_points(data.chunk_center, local_path)
        poly!(local_axis, footprint_path; color = JG.green_100,
            strokecolor = JG.green_dark, strokewidth = 1.6)
        poly!(local_axis, local_polygons; color = (JG.paper, 0.0),
            strokecolor = (JG.green_dark, 0.22), strokewidth = 0.24)
        lines!(local_axis, local_xy; color = (JG.red, 0.76), linewidth = 0.62)
        scatter!(local_axis, local_xy; color = (JG.red, 0.72), markersize = 1.7)
        scatter!(local_axis, local_xy[[1, end]];
            color = [JG.green_dark, JG.red], markersize = 6)
        rowgap!(left, 3)

        right = GridLayout(body[1, 2])
        Label(right[1, 1], "Z7 INDEX CURVE · GLOBAL VIEW";
            fontsize = 11, font = FONT_BODY, color = JG.caption,
            tellwidth = false)
        globe = globe_axis(right[2, 1]; camera_longlat = (12, 24),
            camera_altitude = 1.28)
        dggpoly!(globe, data.global_grid; color = (JG.green_100, 0.30),
            strokecolor = (JG.green_dark, 0.26), strokewidth = 0.35)
        lines!(globe, global_path; source = DGGSTalkFigures.CARTESIAN_SPHERE,
            color = (JG.red, 0.84), linewidth = 0.78, zlevel = 0.004)
        scatter!(globe, global_centroids;
            source = DGGSTalkFigures.CARTESIAN_SPHERE,
            color = (JG.red, 0.70), markersize = 1.8, zlevel = 0.004)
        coastlines!(globe)
        rowgap!(right, 3)

        colsize!(body, 1, Relative(0.46))
        colsize!(body, 2, Relative(0.54))
        colgap!(body, 16)
        fig
    end
end

z7_curve_app() = static_app(z7_curve_figure())

export_z7_curve(dir = joinpath(@__DIR__, "html")) =
    export_html(joinpath(dir, "11c-z7-curve.html"), z7_curve_app())

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_z7_curve()
