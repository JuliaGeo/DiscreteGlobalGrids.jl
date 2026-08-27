isdefined(@__MODULE__, :DGGSTalkFigures) || include(joinpath(@__DIR__, "00-dggs-theme.jl"))

using .DGGSTalkFigures, DiscreteGlobalGrids, DiscreteGlobalGridsVisualization, GLMakie
using ConservativeRegridding, Oceananigans, RasterDataSources, Rasters
import ArchGDAL
import GeometryOps as GO

const FPS, LEVEL = 30, 2
const OUTPUT = joinpath(@__DIR__, "video", "05-h3-stencil.mp4")

smoothstep(t) = t^2 * (3 - 2t)

function hamiltonian_patch(grid, seedpos, nbidx)
    patch = Set(neighbors(grid, seedpos, 2))
    push!(patch, seedpos)
    path, used = Int[seedpos], Set([seedpos])
    function visit!()
        length(path) == length(patch) && return true
        choices = [j for j in nbidx[last(path)] if j in patch && j ∉ used]
        sort!(choices; by = j -> count(k -> k in patch && k ∉ used, nbidx[j]))
        for j in choices
            push!(path, j); push!(used, j)
            visit!() && return true
            pop!(path); delete!(used, j)
        end
        return false
    end
    visit!() || error("could not find an adjacent path through the stencil patch")
    return path
end

function worldclim_on_h3(grid)
    ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH",
        joinpath(DEPOT_PATH[1], "rasterdatasources")))
    raster = Raster(getraster(WorldClim{Climate}, :tavg; month = 1, res = "10m"))
    raw = reverse(Array(raster); dims = 2)
    values = Float64.(coalesce.(raw, 0)); coverage = Float64.((.!ismissing.(raw)))
    source = Oceananigans.LatitudeLongitudeGrid(; size = (2160, 1080, 1),
        longitude = (-180, 180), latitude = (-90, 90), z = (0, 1))
    regridder = ConservativeRegridding.Regridder(GO.Spherical(), grid,
        ConservativeRegridding.Trees.treeify(GO.Spherical(), source))
    numerator, denominator = zeros(ncells(grid)), zeros(ncells(grid))
    ConservativeRegridding.regrid!(numerator, regridder, vec(values); normalize = false)
    ConservativeRegridding.regrid!(denominator, regridder, vec(coverage); normalize = false)
    return numerator, denominator
end

function stencil_data()
    grid = levelgrid(H3System(), LEVEL)
    nbidx = adjacency(grid)
    heat, land = worldclim_on_h3(grid)
    input = [w > 1e-12 ? h / w : NaN for (h, w) in zip(heat, land)]
    packed = [(; heat = h, land = w) for (h, w) in zip(heat, land)]
    function land_mean(center, neighbors)
        center.land > 1e-12 || return NaN
        h, w = center.heat, center.land
        for neighbor in neighbors
            h += neighbor.heat; w += neighbor.land
        end
        return h / w
    end
    output = mapneighbors(CellVector(grid), packed) do _, center, neighbors
        land_mean(center, neighbors)
    end
    seed = cellat(grid, 85.0, 28.0)
    seedpos = something(localindex(grid, seed))
    path = hamiltonian_patch(grid, seedpos, nbidx)
    @assert length(path) == 19 && all(path[i + 1] in nbidx[path[i]] for i in 1:18)
    return grid, nbidx, path, input, output
end

function animation_frames(path, input, output)
    frames = NamedTuple[]
    values = copy(input)
    append!(frames, [(values = copy(values), active = 1, alpha = 0.0) for _ in 1:21])
    for alpha in smoothstep.(range(0, 1; length = 13)[2:end])
        push!(frames, (; values = copy(values), active = 1, alpha))
    end
    for (active, pos) in enumerate(path)
        for k in 1:8
            t = k <= 3 ? 0.0 : smoothstep((k - 3) / 5)
            shown = copy(values); shown[pos] = (1 - t) * input[pos] + t * output[pos]
            push!(frames, (; values = shown, active, alpha = 1.0))
        end
        values[pos] = output[pos]
    end
    for alpha in 1 .- smoothstep.(range(0, 1; length = 13)[2:end])
        push!(frames, (; values = copy(values), active = length(path), alpha))
    end
    append!(frames, [(values = copy(values), active = length(path), alpha = 0.0) for _ in 1:36])
    return frames
end

function render_stencil(path = OUTPUT)
    grid, nbidx, sweep, input, output = stencil_data()
    cells = CellVector(grid)
    frames = animation_frames(sweep, input, output)
    cmap = cgrad([JG.purple, JG.purple_100, JG.paper,
        JG.green_100, colorant"#7dda71"])
    paint(values) = [isfinite(v) ? cmap[clamp((v + 20) / 40, 0, 1)] : JG.paper_off
        for v in values]

    GLMakie.activate!()
    with_theme(DGGSTalkFigures.JG_THEME) do
        fig, body = DGGSTalkFigures.slide_figure()
        axis = DGGSTalkFigures.globe_axis(body[1, 1];
            camera_longlat = (86, 28), camera_altitude = 0.65)
        colors = Observable(paint(input))
        active = Observable(1)
        opacity = Observable(0.0)
        center = @lift([cells[sweep[$active]]])
        neighbors = @lift(cells[collect(nbidx[sweep[$active]])])

        dggpoly!(axis, grid;
            color = colors, strokecolor = (JG.ink, 0.22), strokewidth = 0.7)
        DGGSTalkFigures.coastlines!(axis)
        dggpoly!(axis, system(grid), neighbors;
            color = (JG.paper, 0), strokecolor = @lift((JG.paper, 0.95 * $opacity)),
            strokewidth = 3.2, zlevel = 0.017)
        dggpoly!(axis, system(grid), neighbors;
            color = @lift((JG.paper, 0.03 * $opacity)),
            strokecolor = @lift((JG.green_dark, 0.88 * $opacity)),
            strokewidth = 1.5, zlevel = 0.018)
        dggpoly!(axis, system(grid), center;
            color = (JG.paper, 0), strokecolor = @lift((JG.paper, 0.98 * $opacity)),
            strokewidth = 4.6, zlevel = 0.021)
        dggpoly!(axis, system(grid), center;
            color = @lift((JG.paper, 0.05 * $opacity)),
            strokecolor = @lift((JG.green_dark, $opacity)),
            strokewidth = 2.5, zlevel = 0.022)

        mkpath(dirname(path))
        record(fig, path, frames; framerate = FPS, backend = GLMakie,
            compression = 18, profile = "high", pixel_format = "yuv420p",
            px_per_unit = 1, visible = false, loglevel = "warning") do frame
            colors[] = paint(frame.values)
            active[] = frame.active
            opacity[] = frame.alpha
        end
    end
    println("wrote ", path)
    return path
end

abspath(PROGRAM_FILE) == (@__FILE__) && render_stencil()
