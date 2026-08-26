isdefined(@__MODULE__, :DaggerTalkFigures) ||
    include(joinpath(@__DIR__, "11-dagger-theme.jl"))

using .DGGSTalkFigures, .DaggerTalkFigures
using DiscreteGlobalGrids, DiscreteGlobalGridsVisualization, Makie
using LinearAlgebra
import DimensionalData as DD

# Match the storage geometry used by the stencil animation immediately before
# this slide. The scheduler-facing lookup is:
# one owned interval -> its exact halo cells -> the chunks that own those cells.
const HALO_CHUNK_CELLS = 128
const HALO_CHUNK_NUMBER = 130

function halo_mapping_fixture()
    grid = levelgrid(H3System(), 3)
    cells = CellVector(grid)
    field = DD.DimArray(Float32.(1:length(cells)),
        (Cells(CellLookup(cells)),); name = :synthetic)
    plan = chunkplan(field; chunks = HALO_CHUNK_CELLS, halo = 1)
    chunk = plan[HALO_CHUNK_NUMBER]
    owned = collect(ownedindices(chunk))
    halo = collect(chunkhalo(chunk))

    # Ask the plan which interval owns each halo index rather than assuming a
    # formula for chunk numbering. This is the mapping a scheduler consumes.
    ownerof(index) = something(findfirst(k -> index in ownedindices(plan[k]),
        1:nchunks(plan)))
    suppliers = sort!(unique(ownerof.(halo)))
    supplied = [filter(index -> ownerof(index) == supplier, halo)
        for supplier in suppliers]

    directions = collect.(Tuple.(cell_centroid.(Ref(grid), cells[owned])))
    center = normalize(reduce(+, directions))
    camera_longlat = (rad2deg(atan(center[2], center[1])),
        rad2deg(asin(center[3])))
    return (; grid, cells, plan, chunk, owned, halo, suppliers, supplied,
        camera_longlat)
end

function chunk_strip!(axis, fixture)
    total = nchunks(fixture.plan)
    supplier_set = Set(fixture.suppliers)
    rows = (1:60, 61:120, 121:total)
    ypositions = (0.46, 0.29, 0.12)
    x0, width = 0.10, 0.82

    text!(axis, 0.06, 0.58; text = "HALO CELL → STORAGE OWNER",
        align = (:left, :center), font = FONT_BOLD, fontsize = 11,
        color = JG.ink)

    for (range, y) in zip(rows, ypositions)
        n = length(range)
        boxwidth = width / n
        rects = [Rect2f(x0 + (j - 1) * boxwidth, y, boxwidth - 0.0015, 0.075)
            for j in 1:n]
        colors = [i == fixture.chunk.index ? JG.purple :
            i in supplier_set ? JG.green : JG.paper_off for i in range]
        strokes = [i == fixture.chunk.index ? JG.purple :
            i in supplier_set ? JG.green_dark : JG.hairline for i in range]
        poly!(axis, rects; color = colors, strokecolor = strokes,
            strokewidth = 0.55)
        text!(axis, x0 - 0.018, y + 0.0375; text = string(first(range)),
            align = (:right, :center), font = FONT_BODY, fontsize = 9,
            color = JG.caption)
        text!(axis, x0 + width + 0.018, y + 0.0375; text = string(last(range)),
            align = (:left, :center), font = FONT_BODY, fontsize = 9,
            color = JG.caption)
    end

    text!(axis, 0.06, 0.035;
        text = "$(length(fixture.suppliers)) supplying chunks · only chunk $(fixture.chunk.index) is written",
        align = (:left, :center), font = FONT_BODY, fontsize = 11,
        color = JG.muted)
    return axis
end

function halo_mapping_app()
    fixture = halo_mapping_fixture()
    system_ = system(fixture.grid)
    supplier_cells = reduce(vcat,
        [collect(ownedindices(fixture.plan[s])) for s in fixture.suppliers])

    with_theme(JG_THEME) do
        fig, body = slide_figure()

        spatial_panel = GridLayout(body[1, 1])
        Label(spatial_panel[1, 1], "THE CHUNK AND THE CELLS JUST OUTSIDE IT";
            fontsize = 11, font = FONT_BOLD, color = JG.ink, tellwidth = false)
        spatial = globe_axis(spatial_panel[2, 1];
            camera_longlat = fixture.camera_longlat, camera_altitude = 0.82)
        spatial.tellwidth[] = false
        # Full supplying chunks are faint context; saturated green cells are
        # the exact subset that this one-ring stencil reads.
        dggpoly!(spatial, system_, fixture.cells[supplier_cells];
            color = (JG.green_100, 0.18), strokecolor = (JG.green_dark, 0.18),
            strokewidth = 0.45, zlevel = 0.010)
        dggpoly!(spatial, system_, fixture.cells[fixture.owned];
            color = (JG.purple_100, 0.94), strokecolor = JG.purple,
            strokewidth = 1.15, zlevel = 0.018)
        dggpoly!(spatial, system_, fixture.cells[fixture.halo];
            color = (JG.green_100, 0.98), strokecolor = JG.green_dark,
            strokewidth = 1.55, zlevel = 0.026)
        coastlines!(spatial)
        Label(spatial_panel[3, 1], "PURPLE: 128 OWNED · GREEN: 106 EXACT HALO CELLS";
            fontsize = 10, font = FONT_BODY, color = JG.muted, tellwidth = false)
        rowgap!(spatial_panel, 2)

        mapping = clean_axis(body[1, 2]; limits = (0, 1, 0, 1))
        mapping.tellwidth[] = false
        draw_box!(mapping, 0.05, 0.78, 0.30, 0.13,
            "chunk $(fixture.chunk.index)\n128 owned";
            fill = JG.purple_100, stroke = JG.purple, linewidth = 1.6,
            fontsize = 12)
        draw_arrow!(mapping, Point2f(0.37, 0.845), Point2f(0.52, 0.845);
            color = JG.ink, linewidth = 1.4)
        text!(mapping, 0.445, 0.885; text = "one-ring lookup",
            align = (:center, :bottom), font = FONT_BODY, fontsize = 10,
            color = JG.muted)
        draw_box!(mapping, 0.54, 0.78, 0.38, 0.13,
            "106 exact\nhalo cells";
            fill = JG.green_100, stroke = JG.green_dark, linewidth = 1.6,
            fontsize = 12)
        chunk_strip!(mapping, fixture)

        colsize!(body, 1, Relative(0.54))
        colsize!(body, 2, Relative(0.46))
        colgap!(body, 10)
        static_app(fig)
    end
end

export_halo_mapping(dir = joinpath(@__DIR__, "html")) =
    export_html(joinpath(dir, "16-halo-mapping.html"), halo_mapping_app())

abspath(PROGRAM_FILE) == (@__FILE__) && export_halo_mapping()
