# # The resampling recipe
#
# `dggresample` puts the descent in `resample.jl` behind a plot that follows the
# camera.  What it draws is a `dggpoly` of the resampled cells, so everything
# already true of that plot — the targets, the cut meridian, the poles, the
# CairoMakie path — stays true here; the only new thing is *which* cells.

"""
    Resampled(cells, index, level)

One frame's worth of resampling: the cells to draw, which entry of the user's
value vector each of them takes, and the level they came from.

The three travel together because a plot that had the cells of one frame and the
colours of another would be wrong for as long as it took the second to arrive.
"""
struct Resampled{C}
    cells::C
    index::Vector{Int32}
    level::Int
end

Base.length(r::Resampled) = length(r.index)

"""
    resample_frame(pyramid, view; cellpixels, maxcells, ntasks) -> Resampled

Resample the pyramid for one view: descend to the level worth drawing, then keep
the cells that have a value underneath them.

Cells over a hole in a partial grid are dropped rather than drawn in the NaN
colour, so a coarse view of a ragged set shows the set rather than its bounding
box.
"""
function resample_frame(pyr::CellPyramid, view::ScreenView;
        cellpixels::Real = 3.0, maxcells::Integer = 400_000, ntasks::Integer = 1)
    cells, level = resample(pyr, view; cellpixels, maxcells, ntasks)
    kept, index = withdata(pyr, cells, ntasks)
    return Resampled(CellSet(pyr.system, kept), index, level)
end

"""
    withdata(pyramid, cells, ntasks) -> (kept, index)

Look every cell up in the data and keep the ones that hit something.

This is the nearest-neighbour resample, and it is the one part of a frame that
is proportional to the number of cells drawn, so it is worth spreading over the
threads the plot was given.
"""
function withdata(pyr::CellPyramid, cells::AbstractVector, ntasks::Integer)
    n = length(cells)
    n == 0 && return similar(cells, 0), Int32[]

    nt = clamp(Int(ntasks), 1, max(1, cld(n, 4096)))
    bounds = round.(Int, range(0, n; length = nt + 1))
    parts = Vector{Tuple{typeof(cells), Vector{Int32}}}(undef, nt)
    Threads.@sync for t in 1:nt
        Threads.@spawn parts[t] = _withdata(pyr, cells, bounds[t] + 1, bounds[t + 1])
    end

    kept = similar(cells, 0)
    index = Int32[]
    total = sum(length(p[2]) for p in parts)
    sizehint!(kept, total)
    sizehint!(index, total)
    for (k, ix) in parts
        append!(kept, k)
        append!(index, ix)
    end
    return kept, index
end

# One task's share, in its own function so that its buffers are its own: two
# `similar` calls in one function body are one variable as far as a closure is
# concerned, and two tasks growing one vector is a data race.
function _withdata(pyr::CellPyramid, cells::AbstractVector, lo::Int, hi::Int)
    kept = similar(cells, 0)
    index = Int32[]
    for i in lo:hi
        c = cells[i]
        j = nearest(pyr, c)
        j == 0 && continue
        push!(kept, c)
        push!(index, j % Int32)
    end
    return kept, index
end

"""
    BuildState()

What the last build was looking at, so that the next camera event can tell
whether it still applies.

A build covers a padded rectangle of pixels; while the viewport stays inside it
and the zoom stays within `hysteresis` of what it was, the cells on screen are
already the right cells and nothing is done.
"""
mutable struct BuildState
    pv::Union{Nothing, Makie.Mat4d}
    resolution::Makie.Vec2f
    xmin::Float64
    ymin::Float64
    xmax::Float64
    ymax::Float64
end

BuildState() = BuildState(nothing, Makie.Vec2f(0, 0), 0.0, 0.0, 0.0, 0.0)

function record!(state::BuildState, view::ScreenView)
    state.pv = view.pv
    state.resolution = view.resolution
    state.xmin, state.ymin = view.xmin, view.ymin
    state.xmax, state.ymax = view.xmax, view.ymax
    return state
end

"""
    stale(state, scene, hysteresis) -> Bool

Whether the view has moved far enough to need building again.

The current viewport's corners are carried back into the space the last build
measured in, which works the same way for a map and for a globe: a pan shifts
the rectangle, a zoom scales it, and turning a globe throws the corners across
it.  A rebuild is due when the rectangle leaves what was built, or when the zoom
has changed by more than `hysteresis` — the second test being what keeps the
plot from switching levels back and forth on the boundary between two of them.
"""
function stale(state::BuildState, scene, hysteresis::Real)
    state.pv === nothing && return true
    res = scene.camera.resolution[]
    res == state.resolution || return true

    now = Makie.Mat4d(scene.camera.projectionview[])
    back = try
        inv(now)
    catch
        return true
    end

    xlo = ylo = Inf
    xhi = yhi = -Inf
    for x in (-1.0, 1.0), y in (-1.0, 1.0)
        w = back * Makie.Vec4d(x, y, 0.0, 1.0)
        abs(w[4]) > 1e-30 || return true
        p, _ = to_pixels(state.pv, res, Point3d(w[1] / w[4], w[2] / w[4], w[3] / w[4]))
        isfinite(p[1]) && isfinite(p[2]) || return true
        xlo, xhi = min(xlo, p[1]), max(xhi, p[1])
        ylo, yhi = min(ylo, p[2]), max(yhi, p[2])
    end

    inside = xlo >= state.xmin && ylo >= state.ymin &&
        xhi <= state.xmax && yhi <= state.ymax
    inside || return true

    zoom = (xhi - xlo) / max(Float64(res[1]), 1.0)
    return !(inv(hysteresis) <= zoom <= hysteresis)
end

"""
    dggresample(cells; color = ..., kwargs...)
    dggresample(source, ids; color = ..., kwargs...)

Draw a set of DGGS cells at whatever level of their own hierarchy the current
zoom can actually show, resampled nearest neighbour.

Takes the same arguments as [`dggpoly`](@ref) and draws through it.  The
difference is what reaches the mesh: rather than every cell handed in,
`dggresample` descends the system's hierarchy from its root cells, keeping only
branches that are on screen and hold data, and stops at the level whose cells
come out about `cellpixels` across.  Each of those cells is coloured by the
value of the leaf cell under its centre.

```julia
sys = DGG.IGeo7System()
cells = DGG.CellVector(DGG.query(sys, DGG.MultiOrderCoverage(extent); level = 13))
dggresample(cells; color = elevation)   # draws a few thousand, not sixteen million
```

Neither half of a frame is proportional to the number of cells given: the
descent costs what is on screen, and the resampling costs one point location per
cell drawn.  Zooming in refines, zooming out coarsens, and both keep the picture
about the same size.

Because rebuilding on every camera event would be its own kind of slow, a build
covers a `buffer` factor more than the viewport and stands until the view leaves
it — so panning within the buffer, and zooming within `hysteresis`, cost nothing
at all.

!!! note "Nearest neighbour"
    A drawn cell shows one leaf value, not a summary of the leaves under it, so
    a coarse view of noisy data shows a sample of the noise rather than its
    mean.  This is what makes a frame cost what it does; averaging would have to
    read every leaf cell.
"""
@recipe DGGResample (cells,) begin
    """
    Sets the colour of the cells: a single colour, or one value per cell **of
    the set handed in** — the resampling picks from it, so it is indexed by the
    original cells and not by the cells drawn.
    """
    color = @inherit patchcolor
    """
    How wide a drawn cell should come out, in pixels.  Smaller draws more cells
    and resolves more detail; three is about where a hexagon stops being
    distinguishable from the pixels under it, and each step finer costs seven
    times as many cells.
    """
    cellpixels = 3.0
    """
    How much more than the viewport to build, as a factor.  `1` builds exactly
    what is visible and rebuilds on the smallest pan; larger values trade memory
    and build time for fewer rebuilds.
    """
    buffer = 1.6
    """
    How far the zoom may drift from the zoom a build was made at before another
    is due.  Keep it below the ratio between two levels of the system, or the
    plot will show a level coarser than it needs to.
    """
    hysteresis = 1.5
    "The most cells to draw.  The descent stops rather than refine past this."
    maxcells = 400_000
    """
    Whether to follow the camera.  `false` builds once for the view at the time
    and leaves it, which is what a figure being saved to a file wants.
    """
    dynamic = true
    "Controls whether lights affect the cells."
    shading = Makie.NoShading
    "How the cells are handed to the backend; see [`dggpoly`](@ref)."
    primitive = Makie.automatic
    "Whether to split cells that straddle the map's cut meridian."
    wrap = true
    "How many tasks build the mesh and do the resampling."
    ntasks = Threads.nthreads()
    "The resampled frame currently drawn.  Set by the plot; read it to see which level is on screen."
    resampled = nothing
    cycle = [:color => :patchcolor]
    Makie.mixin_generic_plot_attributes()...
    Makie.mixin_colormap_attributes()...
end

Makie.convert_arguments(::Type{<:DGGResample}, cs::CellSet) = (cs,)
Makie.convert_arguments(::Type{<:DGGResample}, x) = (cellset(x),)
Makie.convert_arguments(::Type{<:DGGResample}, source, ids::AbstractVector) =
    (cellset(source, ids),)

# The axis's limits must not depend on what is currently drawn: they drive the
# camera, the camera drives the resampling, and the resampling would then drive
# the limits.  Reporting the whole data set, always, breaks that loop.
# The box is in the space the plot draws in — the target has already projected
# it — so only the model matrix may be applied to it, never the axis's transform
# a second time.  That is the same contract the children carry through
# `transformation = :inherit_model`.
Makie.data_limits(plot::DGGResample) = plot.datalimits[]
Makie.boundingbox(plot::DGGResample, space::Symbol = :data) =
    Makie.apply_transform_and_model(plot.model[], identity, plot.datalimits[], Point3d)

function Makie.plot!(plot::DGGResample{<:Tuple{<:CellSet}})
    Makie.map!(cs -> (CellPyramid(cs),), plot, [:cells], [:pyramid])
    Makie.map!(plot, [:pyramid, :transform_func], [:datalimits]) do pyr, transform_func
        return (extentbox(pyr, plot_target(transform_func)),)
    end
    Makie.map!(r -> (r.cells,), plot, [:resampled], [:displayed])
    Makie.map!((r, c) -> (pick(c, r.index),), plot, [:resampled, :color], [:cellcolor])

    scene = Makie.parent_scene(plot)
    state = BuildState()
    inflight = Ref(false)

    function rebuild!(force::Bool)
        inflight[] && return nothing
        (force || plot.dynamic[]) || return nothing
        force || stale(state, scene, Float64(plot.hysteresis[])) || return nothing

        target = plot_target(Makie.to_value(plot.transform_func))
        plot.wrap[] || (target = uncut(target))
        view = ScreenView(target, scene, Float64(plot.buffer[]))
        frame = resample_frame(plot.pyramid[], view;
            cellpixels = plot.cellpixels[], maxcells = plot.maxcells[],
            ntasks = plot.ntasks[])

        inflight[] = true
        try
            plot.resampled[] = frame
        finally
            inflight[] = false
        end
        record!(state, view)
        return nothing
    end

    rebuild!(true)
    for attribute in (:pyramid, :cellpixels, :maxcells, :buffer, :wrap, :ntasks)
        on(_ -> rebuild!(true), getproperty(plot, attribute))
    end
    on(_ -> rebuild!(false), scene.camera.projectionview)

    dggpoly!(
        plot, plot.displayed;
        color = plot.cellcolor,
        colormap = plot.colormap,
        colorscale = plot.colorscale,
        colorrange = plot.colorrange,
        lowclip = plot.lowclip,
        highclip = plot.highclip,
        nan_color = plot.nan_color,
        alpha = plot.alpha,
        shading = plot.shading,
        primitive = plot.primitive,
        wrap = plot.wrap,
        ntasks = plot.ntasks,
        strokewidth = 0,
        visible = plot.visible,
        transparency = plot.transparency,
        inspectable = plot.inspectable,
        space = plot.space,
        depth_shift = plot.depth_shift,
        clip_planes = plot.clip_planes,
        fxaa = plot.fxaa,
    )

    return plot
end

"""
    pick(color, index) -> color

The colour of each drawn cell: `color[index]` when there is one value per cell of
the original set, and `color` itself otherwise.

Recolouring therefore costs a gather, exactly as it does in [`dggpoly`](@ref) —
the resampling does not run again.
"""
function pick(color, index::Vector{Int32})
    color isa AbstractVector || return color
    length(color) >= maximum(index; init = Int32(0)) || throw(ArgumentError(
        "color has $(length(color)) entries, which does not cover the cells given"))
    return color[index]
end

"""
    extentbox(pyramid, target) -> Rect3d

A box around the whole data set in the space the target draws in, taken from the
sample the pyramid measured its cap with.

It is what the axis is shown instead of the cells actually drawn, which is what
keeps a resampled plot from moving the limits that decide how it is resampled.
"""
function extentbox(pyr::CellPyramid, target)
    lo = Point3d(Inf, Inf, Inf)
    hi = Point3d(-Inf, -Inf, -Inf)
    for p in pyr.samplepoints
        q = project_probe(target, p)
        z = length(q) == 3 ? q[3] : 0.0
        lo = Point3d(min(lo[1], q[1]), min(lo[2], q[2]), min(lo[3], z))
        hi = Point3d(max(hi[1], q[1]), max(hi[2], q[2]), max(hi[3], z))
    end
    all(isfinite, lo) && all(isfinite, hi) ||
        return Makie.Rect3d(Makie.Vec3d(0), Makie.Vec3d(0))
    return Makie.Rect3d(Makie.Vec3d(lo...), Makie.Vec3d((hi .- lo)...))
end
