isdefined(@__MODULE__, :DaggerTalkFigures) ||
    include(joinpath(@__DIR__, "11-dagger-theme.jl"))

using .DGGSTalkFigures, .DaggerTalkFigures
using DiscreteGlobalGrids
using DiscreteGlobalGridsVisualization, Makie
using LinearAlgebra
using SparseArrays
import GlobalRegridding as GR
import DimensionalData as DD

"True when one conservative chunk-pair block contains at least one cell intersection."
function _nonempty_conservative_pair(dstspace, srcspace, d::Int, s::Int)
    block = GR.weightblock(Conservative(), dstspace,
        GR.ownedindices(dstspace, d), srcspace, GR.ownedindices(srcspace, s))
    return nnz(block.weights) > 0
end

"A small, exact chunk graph over the same analytic raster as figure 12."
function chunkgraph_fixture()
    source = synthetic_cube(; nt = 4)
    xchunks = [i:(i + 2) for i in 1:3:24]
    ychunks = [j:(j + 1) for j in 1:2:12]
    srcspace = GR.RasterGrid(DD.dims(source.data);
        chunks = (xchunks, ychunks))
    grid = levelgrid(IGeo7System(), 2)
    dstspace = DGGSpace(grid; chunkcells = 8)
    method = Conservative()
    # The raster chunk index intentionally answers a conservative broad phase:
    # a row can include a nearby chunk whose exact weight block is empty.  This
    # presentation is about data actually read for regridding, so give the plan
    # a narrow phase that keeps only nonempty conservative intersection blocks.
    refine = (d, s) -> _nonempty_conservative_pair(dstspace, srcspace, d, s)
    plan = GR.ChunkedPlan(method, Weighted(0.01), dstspace, srcspace;
        refine, narrow = :nonempty_conservative_intersections)
    graph = GR.dependencies(plan)

    # Prefer an exact three-edge row near the centre of the map.  This is D23 in
    # the current IGeo7 level-2 grid and reaches the three geographically
    # adjacent source chunks S12, S19 and S20.
    degrees = [GR.sourcedegree(graph, d) for d in 1:GR.ndestinationchunks(graph)]
    degree_three = findall(==(3), degrees)
    isempty(degree_three) && error("the fixture needs a destination with three exact sources")
    destination = argmin(degree_three) do d
        cap = GR.chunkextents(dstspace)[d]
        lon = rad2deg(atan(cap.point[2], cap.point[1]))
        lat = rad2deg(asin(cap.point[3]))
        abs(lat + 20) + 0.15abs(lon + 30)
    end
    return (; source, srcspace, grid, dstspace, graph, destination)
end

function _chunk_rect(source, srcspace, chunk)
    xr, yr = GR.chunkranges(srcspace, chunk, (length(source.lon), length(source.lat)))
    dx = source.lon[2] - source.lon[1]
    dy = source.lat[2] - source.lat[1]
    x0, y0 = source.lon[first(xr)] - dx / 2, source.lat[first(yr)] - dy / 2
    return Rect2f(x0, y0, length(xr) * dx, length(yr) * dy)
end

"A closed longitude–latitude outline for a source-raster chunk."
function _chunk_outline(source, srcspace, chunk)
    rect = _chunk_rect(source, srcspace, chunk)
    x, y = rect.origin
    w, h = rect.widths
    bottom = Point2f.(range(x, x + w; length = 25), y)
    right = Point2f.(x + w, range(y, y + h; length = 13))
    top = Point2f.(range(x + w, x; length = 25), y + h)
    left = Point2f.(x, range(y + h, y; length = 13))
    vcat(bottom, right[2:end], top[2:end], left[2:end])
end

function _edge!(axis, a, b; color, width = 2.0)
    lines!(axis, Point2f[a, b]; color, linewidth = width)
    direction = Point2f(b[1] - a[1], b[2] - a[2])
    scatter!(axis, [Point2f(b)]; marker = :rtriangle,
        rotation = atan(direction[2], direction[1]), markersize = 8, color)
end

function regrid_chunkgraph_figure()
    data = chunkgraph_fixture()
    g, d = data.graph, data.destination
    needed = Int.(collect(GR.sourcesof(g, d)))
    selected = Set(needed)
    output_cells = GR.ownedindices(data.dstspace, d)
    output_cell_objects = CellVector(data.grid)[collect(output_cells)]
    output_direction = normalize(reduce(+,
        collect.(Tuple.(cell_centroid.(Ref(data.grid), output_cell_objects)))))
    output_camera = (rad2deg(atan(output_direction[2], output_direction[1])),
        rad2deg(asin(output_direction[3])))

    with_theme(JG_THEME) do
        fig, body = slide_figure()

        # Source chunks stay rectilinear: this is deliberately a memory/storage
        # view, not a second geographic map.
        source = Axis(body[1, 1];
            title = "SOURCE MAP · $(GR.nsourcechunks(g)) CHUNKS",
            titlesize = 12, titlefont = FONT_BODY, xlabel = "longitude",
            ylabel = "latitude", aspect = 2)
        heatmap!(source, data.source.lon, data.source.lat,
            parent(data.source.data)[:, :, 1]; colormap = FIELD_COLORMAP,
            colorrange = FIELD_RANGE, alpha = 0.34)
        for s in 1:GR.nsourcechunks(g)
            active = s in selected
            poly!(source, _chunk_rect(data.source, data.srcspace, s);
                color = active ? (JG.green, 0.28) : (:white, 0.12),
                strokecolor = active ? JG.green_dark : (JG.ink, 0.30),
                strokewidth = active ? 2.2 : 0.75)
            r = _chunk_rect(data.source, data.srcspace, s)
            active && text!(source,
                r.origin[1] + r.widths[1] / 2,
                r.origin[2] + r.widths[2] / 2;
                text = "S$s", align = (:center, :center), fontsize = 9,
                font = FONT_BOLD, color = JG.green_dark)
        end
        xlims!(source, -180, 180); ylims!(source, -90, 90)

        # The central panel is a selected real row of the bipartite relation.
        # It provides the readable graph-level abstraction of the two maps.
        relation = clean_axis(body[1, 2]; limits = (0, 1, 0, 1))
        text!(relation, 0.5, 0.96; text = "EXACT OVERLAP GRAPH",
            align = (:center, :center), fontsize = 10, font = FONT_BODY, color = JG.ink)
        text!(relation, 0.08, 0.86; text = "SOURCES", align = (:left, :center),
            fontsize = 10, font = FONT_BOLD, color = JG.ink)
        text!(relation, 0.92, 0.86; text = "OUTPUT", align = (:right, :center),
            fontsize = 10, font = FONT_BOLD, color = JG.ink)

        node_y = collect(range(0.72, 0.32; length = max(length(needed), 1)))
        for (i, s) in enumerate(needed)
            a, b = Point2f(0.22, node_y[i]), Point2f(0.74, 0.53)
            _edge!(relation, a, b; color = (JG.green_dark, 0.76))
            scatter!(relation, [a]; marker = :circle, markersize = 20,
                color = JG.green, strokecolor = JG.green_dark, strokewidth = 1)
            text!(relation, 0.22, node_y[i]; text = "S$s", align = (:center, :center),
                fontsize = 9, font = FONT_BOLD, color = :white)
        end
        scatter!(relation, [Point2f(0.74, 0.53)]; marker = :rect, markersize = 30,
            color = JG.purple, strokecolor = JG.ink, strokewidth = 1.2)
        text!(relation, 0.74, 0.53; text = "D$d", align = (:center, :center),
            fontsize = 10, font = FONT_BOLD, color = :white)
        text!(relation, 0.5, 0.16;
            text = "D$d overlaps\nS$(join(needed, " · S"))",
            align = (:center, :center), fontsize = 10, font = FONT_BOLD,
            color = JG.green_dark)

        # Destination is geographic and hexagonal. Only Dd is highlighted;
        # all other cells remain visible as context for the exact source
        # chunks whose conservative weight blocks are nonempty.
        target_panel = GridLayout(body[1, 3])
        Label(target_panel[1, 1], "DESTINATION · IGEO7 HEXAGON CHUNKS";
            fontsize = 12, font = FONT_BODY, color = JG.ink, tellwidth = false)
        target = globe_axis(target_panel[2, 1]; camera_longlat = output_camera,
            camera_altitude = 1.55)
        cell_colors = fill(JG.green_50, ncells(data.grid))
        cell_colors[output_cells] .= JG.purple_100
        dggpoly!(target, data.grid; color = cell_colors,
            strokecolor = (JG.ink, 0.35), strokewidth = 0.7)
        # Keep the source chunks geographically legible on the destination:
        # they are faint red outlines, not a competing raster layer.
        for s in needed
            lines!(target, _chunk_outline(data.source, data.srcspace, s);
                source = "+proj=longlat +R=1", color = (JG.red, 0.48),
                linewidth = 1.35, zlevel = 0.025)
        end
        coastlines!(target)
        Label(target_panel[3, 1],
            "D$d owns $(length(output_cells)) hex cells · $(length(needed)) source overlaps";
            fontsize = 11, font = FONT_BODY, color = JG.ink, tellwidth = false)
        rowgap!(target_panel, 2)

        colsize!(body, 1, Relative(0.24)); colsize!(body, 2, Relative(0.18)); colsize!(body, 3, Relative(0.58))
        colgap!(body, 12)
        fig
    end
end

regrid_chunkgraph_app() = static_app(regrid_chunkgraph_figure())

export_regrid_chunkgraph(dir = joinpath(@__DIR__, "html")) =
    export_html(joinpath(dir, "14-regrid-chunkgraph.html"), regrid_chunkgraph_app())

abspath(PROGRAM_FILE) == (@__FILE__) && export_regrid_chunkgraph()
