isdefined(@__MODULE__, :DaggerTalkFigures) ||
    include(joinpath(@__DIR__, "11-dagger-theme.jl"))

using .DGGSTalkFigures, .DaggerTalkFigures
using DiscreteGlobalGrids, GeometryBasics, Makie
using LinearAlgebra

function tangent_frame(center)
    n = normalize(collect(Tuple(center)))
    reference = abs(n[3]) < 0.9 ? [0.0, 0.0, 1.0] : [0.0, 1.0, 0.0]
    east = normalize(cross(reference, n))
    north = normalize(cross(n, east))
    return n, east, north
end

function tangent_point(point, frame; scale = 1.0)
    _, east, north = frame
    p = collect(Tuple(point))
    Point2f(scale * dot(p, east), scale * dot(p, north))
end

function tangent_polygon(system, cell, frame; scale = 1.0)
    GeometryBasics.Polygon(tangent_point.(cell_boundary(system, cell), Ref(frame);
        scale))
end

function ordered_border_points(system, root, leaf, frame)
    cells = collect(border(subtree(system, root, leaf); cells = true))
    points = tangent_point.(cell_centroid.(Ref(levelgrid(system, leaf)), cells),
        Ref(frame))
    sort!(points; by = p -> atan(p[2], p[1]))
    push!(points, first(points))
    return points
end

function storage_subtree_app()
    system = IGeo7System()
    parent_grid = levelgrid(system, 5)
    root = cellat(parent_grid, -105.0, 40.0)
    center = cell_centroid(parent_grid, root)
    frame = tangent_frame(center)
    interval = descendant_range(system, root, 12)
    borderline = ordered_border_points(system, root, 12, frame)

    with_theme(JG_THEME) do
        fig, body = slide_figure()
        diagram = clean_axis(body[1, 1]; limits = (-1.12, 1.12, -0.72, 0.72))

        cell_scale = 48.0
        parent_poly = tangent_polygon(system, root, frame)
        parent_draw = GeometryBasics.Polygon([
            Point2f(-0.58, 0.08) + cell_scale * p for p in parent_poly.exterior])
        poly!(diagram, parent_draw; color = JG.green_100,
            strokecolor = JG.green_dark, strokewidth = 2.2)
        text!(diagram, -0.58, 0.54; text = "LEVEL 5",
            font = FONT_BOLD, fontsize = 13, color = JG.green_dark,
            align = (:center, :center))
        text!(diagram, -0.58, 0.44; text = "one chunk key",
            font = FONT_BODY, fontsize = 11, color = JG.caption,
            align = (:center, :center))

        lines!(diagram, [Point2f(p[1] * cell_scale + 0.58, p[2] * cell_scale + 0.08)
            for p in borderline]; color = JG.purple, linewidth = 2.3)
        parent_overlay = GeometryBasics.Polygon([
            Point2f(0.58, 0.08) + cell_scale * p for p in parent_poly.exterior])
        poly!(diagram, parent_overlay; color = (JG.paper, 0.0),
            strokecolor = JG.green_dark, strokewidth = 1.0)
        text!(diagram, 0.58, 0.54; text = "LEVEL 12",
            font = FONT_BOLD, fontsize = 13, color = JG.purple,
            align = (:center, :center))
        text!(diagram, 0.58, 0.44; text = "823,543 descendant cells",
            font = FONT_BODY, fontsize = 11, color = JG.caption,
            align = (:center, :center))

        draw_arrow!(diagram, (-0.19, 0.08), (0.19, 0.08);
            color = JG.ink, linewidth = 1.8)
        text!(diagram, 0.0, 0.16; text = "descendant_range",
            font = FONT_BODY, fontsize = 10, color = JG.caption,
            align = (:center, :bottom))

        strip = Rect2f(-0.87, -0.66, 1.74, 0.075)
        poly!(diagram, strip; color = JG.paper_off,
            strokecolor = JG.hairline, strokewidth = 1)
        poly!(diagram, Rect2f(-0.79, -0.645, 1.58, 0.045);
            color = JG.green_100, strokecolor = JG.green_dark, strokewidth = 0.8)
        text!(diagram, -0.79, -0.565; text = string(first(interval)),
            font = FONT_BODY, fontsize = 9, color = JG.caption,
            align = (:left, :bottom))
        text!(diagram, 0.79, -0.565; text = string(last(interval)),
            font = FONT_BODY, fontsize = 9, color = JG.caption,
            align = (:right, :bottom))
        text!(diagram, 0.0, -0.622; text = "one contiguous address interval",
            font = FONT_BOLD, fontsize = 10, color = JG.green_dark,
            align = (:center, :center))

        static_app(fig)
    end
end

export_storage_subtree(dir = joinpath(@__DIR__, "html")) =
    export_html(joinpath(dir, "11a-storage-subtree.html"), storage_subtree_app())

abspath(PROGRAM_FILE) == (@__FILE__) && export_storage_subtree()
