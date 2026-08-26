isdefined(@__MODULE__, :DaggerTalkFigures) ||
    include(joinpath(@__DIR__, "11-dagger-theme.jl"))

using .DGGSTalkFigures, .DaggerTalkFigures
using DiscreteGlobalGrids, DiscreteGlobalGridsVisualization
using FlyThroughPaths, GLMakie, LinearAlgebra, Makie, Statistics
import DimensionalData as DD

const STENCIL_ORDER_OUTPUT = joinpath(@__DIR__, "video", "15-stencil-order.mp4")
const STENCIL_CHUNK_SIZE = 128

function stencil_fixture()
    grid = levelgrid(H3System(), 3)
    cells = CellVector(grid)
    adjacency_table = adjacency(grid)
    centroids = cell_centroid.(Ref(grid), collect(cells))
    input = [begin
        x, y, z = Tuple(p)
        0.60sin(8atan(y, x)) + 0.36cos(7asin(z)) + 0.24(isodd(i) ? -1 : 1)
    end for (i, p) in enumerate(centroids)]
    output = mapneighbors(cells, input; threaded = false) do _, center, neighbors
        mean((center, neighbors...))
    end

    axis = DD.DimArray(Float32.(input),
        (Cells(CellLookup(cells)),); name = :synthetic)
    plan = chunkplan(axis; chunks = STENCIL_CHUNK_SIZE, halo = 1)
    score(chunk) = length(unique(cld.(chunkhalo(chunk), STENCIL_CHUNK_SIZE)))
    function angular_spread(chunk)
        indices = collect(ownedindices(chunk))
        directions = collect.(Tuple.(centroids[indices]))
        mean_direction = normalize(reduce(+, directions))
        maximum(acos(clamp(dot(mean_direction, direction), -1, 1))
            for direction in directions)
    end
    scores = score.(collect(plan))
    # Keep a dependency-rich example, but among the top two supplier counts
    # choose the most spatially compact interval.  Its full storage-order walk
    # then remains visible in one fixed camera view at the finer H3 level.
    candidates = findall(>=(maximum(scores) - 1), scores)
    spreads = angular_spread.([plan[k] for k in candidates])
    chunk_number = candidates[argmin(spreads)]
    chunk = plan[chunk_number]
    owned = collect(ownedindices(chunk))
    halo = collect(chunkhalo(chunk))
    supplying = sort!(unique(cld.(halo, STENCIL_CHUNK_SIZE)))
    direction = normalize(reduce(+, collect.(Tuple.(centroids[owned]))))
    camera_longlat = (rad2deg(atan(direction[2], direction[1])),
        rad2deg(asin(direction[3])))
    return (; grid, cells, adjacency_table, input, output, plan, chunk_number,
        chunk, owned, halo, supplying, camera_longlat)
end

function field_colors(values, lo, hi)
    gradient = cgrad(FIELD_COLORMAP)
    [gradient[clamp((value - lo) / (hi - lo), 0, 1)] for value in values]
end

function render_stencil_order(path = STENCIL_ORDER_OUTPUT)
    GLMakie.activate!()
    data = stencil_fixture()
    lo, hi = extrema(vcat(data.input, data.output))

    hold_in, hold_out = 18, 30
    frames = vcat(fill(0, hold_in), collect(eachindex(data.owned)),
        fill(length(data.owned), hold_out))
    shown = copy(data.input)

    with_theme(JG_THEME) do
        fig, body = slide_figure()
        memory_panel = GridLayout(body[1, 2])
        axis = globe_axis(body[1, 1]; camera_longlat = data.camera_longlat,
            camera_altitude = 0.82)

        field = Observable(field_colors(shown, lo, hi))
        active_pos = Observable(data.owned[1])
        active_cell = @lift([data.cells[$active_pos]])
        active_neighbors = @lift(data.cells[collect(data.adjacency_table[$active_pos])])
        active_alpha = Observable(0.0)

        dggpoly!(axis, data.grid; color = field,
            strokecolor = (JG.ink, 0.20), strokewidth = 0.55)
        coastlines!(axis)
        dggpoly!(axis, system(data.grid), active_neighbors;
            color = @lift((JG.green_100, 0.30 * $active_alpha)),
            strokecolor = @lift((JG.green_dark, 0.92 * $active_alpha)),
            strokewidth = 2.0, zlevel = 0.020)
        dggpoly!(axis, system(data.grid), active_cell;
            color = @lift((JG.purple_100, 0.92 * $active_alpha)),
            strokecolor = @lift((JG.purple, $active_alpha)),
            strokewidth = 3.2, zlevel = 0.026)
        memory = clean_axis(memory_panel[1, 1]; limits = (0, 1, 0, 1))
        nchunks_total = nchunks(data.plan)
        chunk_rects = [Rect2f(0.04 + (i - 1) * 0.92 / nchunks_total,
            0.68, 0.92 / nchunks_total - 0.002, 0.12) for i in 1:nchunks_total]
        chunk_colors = Observable(fill(JG.paper_off, nchunks_total))
        poly!(memory, chunk_rects; color = chunk_colors,
            strokecolor = JG.hairline, strokewidth = 0.4)
        text!(memory, 0.04, 0.84; text = "STORAGE CHUNKS",
            align = (:left, :bottom), font = FONT_BOLD,
            fontsize = 11, color = JG.ink)

        owned_rects = [Rect2f(0.04 + (i - 1) * 0.92 / length(data.owned),
            0.34, 0.92 / length(data.owned), 0.12) for i in eachindex(data.owned)]
        owned_colors = Observable(fill(JG.paper_off, length(data.owned)))
        poly!(memory, owned_rects; color = owned_colors,
            strokecolor = (JG.hairline, 0.55), strokewidth = 0.25)
        cursor_x = Observable(0.04)
        vlines!(memory, cursor_x; color = JG.purple, linewidth = 3.0)
        text!(memory, 0.04, 0.50; text = "OWNED INTERVAL",
            align = (:left, :bottom), font = FONT_BOLD,
            fontsize = 11, color = JG.ink)
        text!(memory, 0.04, 0.29; text = string(first(data.owned)),
            align = (:left, :top), font = FONT_BODY,
            fontsize = 10, color = JG.caption)
        text!(memory, 0.96, 0.29; text = string(last(data.owned)),
            align = (:right, :top), font = FONT_BODY,
            fontsize = 10, color = JG.caption)

        progress = Observable("read halo")
        active_label = Observable("chunk $(data.chunk_number)")
        text!(memory, 0.04, 0.15; text = progress,
            align = (:left, :center), font = FONT_BOLD,
            fontsize = 16, color = JG.purple)
        text!(memory, 0.04, 0.07; text = active_label,
            align = (:left, :center), font = FONT_BODY,
            fontsize = 11, color = JG.caption)

        colsize!(body, 1, Relative(0.64))
        colsize!(body, 2, Relative(0.36))
        colgap!(body, 12)

        mkpath(dirname(path))
        record(fig, path, frames; framerate = 24, backend = GLMakie,
            compression = 18, profile = "high", pixel_format = "yuv420p",
            px_per_unit = 1, visible = false, loglevel = "warning") do step
            k = step == 0 ? 1 : step
            pos = data.owned[k]
            active_pos[] = pos
            active_alpha[] = step == 0 ? 0.0 : 1.0

            if step > 0
                shown[pos] = data.output[pos]
                field[] = field_colors(shown, lo, hi)
            end

            chunks = fill(JG.paper_off, nchunks_total)
            chunks[data.supplying] .= JG.green_100
            chunks[data.chunk_number] = JG.purple_100
            for neighbor in data.adjacency_table[pos]
                chunks[cld(neighbor, STENCIL_CHUNK_SIZE)] = JG.green
            end
            chunks[data.chunk_number] = JG.purple
            chunk_colors[] = chunks

            cells = fill(JG.paper_off, length(data.owned))
            step > 1 && (cells[1:(step - 1)] .= JG.green_100)
            step > 0 && (cells[step] = JG.purple)
            owned_colors[] = cells
            cursor_x[] = 0.04 + (k - 0.5) * 0.92 / length(data.owned)
            progress[] = step == 0 ? "read halo" :
                step == length(data.owned) ? "write owned results" : "compute → advance"
            active_label[] = "global cell $(pos) · $(k) / $(length(data.owned))"
        end
    end
    println("wrote ", path)
    return path
end

abspath(PROGRAM_FILE) == (@__FILE__) && render_stencil_order()
