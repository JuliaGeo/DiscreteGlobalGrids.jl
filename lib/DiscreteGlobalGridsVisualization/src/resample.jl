# # Resampling to the screen
#
# `dggpoly` makes drawing `n` cells as cheap as it can be made.  `dggresample`
# asks a different question: how few cells can be drawn without the picture
# changing?  A figure is on the order of a million pixels, so past a certain
# level every extra cell lands under a pixel that another cell already owns, and
# the work is spent on a difference nobody can see.
#
# The answer is the hierarchy the data already lives in.  Rather than draw a
# level-13 cell set, draw the level of the *same* system whose cells come out a
# couple of pixels across, and colour each of those cells by nearest neighbour —
# the value of the leaf cell under its centre.  Which level that is depends on
# the zoom, so it is chosen per view, and the cells to draw are found by
# descending the hierarchy from its roots, keeping only branches that are both
# on screen and hold data.  Neither half of the cost depends on how many cells
# the user handed in.
#
# The view moves constantly and rebuilding is not free, so what is built covers
# more than the viewport — a `buffer` factor around it — and stands until the
# view leaves that window or the zoom changes enough to want another level.

# WGS84 authalic radius, in metres: the sphere `DiscreteGlobalGrids.cellsize`
# lays its areas on, so the two agree about what a metre is.
const EARTH_RADIUS = 6371007.180918475

"""
    to_pixels(pv, resolution, p) -> (Point2d, depth)

Where a point of the plot's own space lands on screen, and how deep it is.

`pv` is the whole way from the plot's space to clip space — the camera, and
before it the Float32 rescaling an axis applies so that a deep zoom into large
coordinates still has precision left.  Leaving that factor out is not a small
error: it is the difference between reading the camera an axis is using and
reading the one it would use if it never zoomed.

The matrix is applied here rather than through `Makie.project` because these
points are already in the space the recipe draws in — the target has done the
projecting — and going through the scene would apply the axis's transform to
them a second time.
"""
@inline function to_pixels(pv, resolution, p)
    z = length(p) == 3 ? p[3] : 0.0
    v = pv * Makie.Vec4d(p[1], p[2], z, 1.0)
    w = v[4]
    x = (v[1] / w * 0.5 + 0.5) * resolution[1]
    y = (v[2] / w * 0.5 + 0.5) * resolution[2]
    return Point2d(x, y), v[3] / w
end

"""
    ScreenView(target, scene, buffer)

Everything a descent needs to know about where the plot is being looked at: the
camera, the size of the picture, and the padded rectangle of pixels a cell has
to touch to be worth keeping.

`buffer` widens that rectangle beyond the viewport, so panning a little reuses
what is already built instead of asking for it again.
"""
struct ScreenView{T}
    target::T
    pv::Makie.Mat4d
    resolution::Makie.Vec2f
    xmin::Float64
    ymin::Float64
    xmax::Float64
    ymax::Float64
    centredepth::Float64
end

function ScreenView(target, scene, buffer::Real)
    pv = Makie.Mat4d(scene.camera.projectionview[]) *
        Makie.f32_convert_matrix(scene, :data)
    res = scene.camera.resolution[]
    pad = (max(buffer, 1.0) - 1) / 2
    w, h = Float64(res[1]), Float64(res[2])
    # The depth of the origin is where the far side of a globe begins: a point
    # deeper than the centre of the earth is behind it.  On a planar target the
    # comparison never fires, because every vertex shares the origin's depth.
    _, centredepth = to_pixels(pv, res, Point3d(0, 0, 0))
    return ScreenView(target, pv, res, -pad * w, -pad * h, (1 + pad) * w,
        (1 + pad) * h, centredepth)
end

"""
    capdisc(view, cap) -> (radius_px, scale, visible)

The screen footprint of a spherical cap: how many pixels its rim reaches, how
many pixels a radian comes out as there, and whether any of it is inside the
padded viewport.

A cap is the right thing to test because `node_extent` gives one that covers a
cell *and every cell beneath it*, so a cap that misses the screen rules out the
whole subtree in one go.

The rim is measured in two perpendicular directions rather than one, because on
a map a circle on the sphere is not a circle on the screen — in a lon/lat axis
at fifty degrees north it is nearly twice as wide as it is tall.  Visibility
takes the larger of the two, so nothing is culled that should not be; the scale
takes their geometric mean, which is the one that makes a *level* come out at
the right number of cells for the area of the screen.
"""
function capdisc(view::ScreenView, cap)
    centre, depth = to_pixels(view.pv, view.resolution, project_probe(view.target, cap.point))
    a, b = rimpoints(cap.point, cap.radius)
    ra = _pixeldistance(view, centre, a)
    rb = _pixeldistance(view, centre, b)

    # Two samples of a linear map bound its largest stretch to within √2; the
    # margin buys back the directions they missed, and being generous here costs
    # candidates rather than correctness.
    radius = 1.5 * max(ra, rb)
    scale = cap.radius > 0 ? sqrt(ra * rb) / cap.radius : 0.0

    # On a globe the far hemisphere projects onto the near one; depth separates
    # them, and a cap large enough to reach around the edge is kept either way.
    behind = depth > view.centredepth && cap.radius < 1.0
    visible = !behind &&
        centre[1] + radius >= view.xmin && centre[1] - radius <= view.xmax &&
        centre[2] + radius >= view.ymin && centre[2] - radius <= view.ymax
    return radius, scale, visible
end

function _pixeldistance(view::ScreenView, centre, p)
    q, _ = to_pixels(view.pv, view.resolution, project_probe(view.target, p))
    return hypot(q[1] - centre[1], q[2] - centre[2])
end

"""
    rimpoints(p, r) -> (q1, q2)

Two points on the rim of the cap centred at `p` with angular radius `r`, a
quarter turn apart around it.
"""
function rimpoints(p, r::Real)
    ax = abs(p[1]) < 0.9 ? (1.0, 0.0, 0.0) : (0.0, 1.0, 0.0)
    ux = ax[2] * p[3] - ax[3] * p[2]
    uy = ax[3] * p[1] - ax[1] * p[3]
    uz = ax[1] * p[2] - ax[2] * p[1]
    n = sqrt(ux^2 + uy^2 + uz^2)
    n > 0 || return p, p
    ux, uy, uz = ux / n, uy / n, uz / n
    # The second tangent completes a right-handed frame with `p` and the first.
    vx = p[2] * uz - p[3] * uy
    vy = p[3] * ux - p[1] * uz
    vz = p[1] * uy - p[2] * ux
    s, c = sincos(min(r, Float64(pi) / 2))
    return (DGG.UnitSphericalPoint(c * p[1] + s * ux, c * p[2] + s * uy, c * p[3] + s * uz),
        DGG.UnitSphericalPoint(c * p[1] + s * vx, c * p[2] + s * vy, c * p[3] + s * vz))
end

"""
    resample(pyramid, view; cellpixels, maxcells) -> (cells, level)

Choose a level and the cells of it worth drawing.

The descent starts at the system's root cells and refines one level at a time,
dropping every branch whose covering cap is off screen or holds none of the
data.  It stops when the surviving cells are about `cellpixels` across, when the
data's own level is reached, or when refining once more would cost more than
`maxcells` — whichever comes first.

The work at each level is proportional to what survived the level above, so a
zoomed-in view costs the same on a sixteen-million-cell set as on a small one.
"""
function resample(pyr::CellPyramid, view::ScreenView;
        cellpixels::Real = 3.0, maxcells::Integer = 400_000, ntasks::Integer = 1)
    cells, scales = keepvisible(pyr, view, collect(DGG.rootcells(pyr.system)), ntasks)
    level = pyr.rootlevel

    while level < pyr.leaflevel
        isempty(cells) && break
        # Stop at the level nearest the target rather than the first level
        # under it: cells shrink by a fixed ratio, so the level above the target
        # is often the closer of the two, and the one below costs seven times as
        # much to draw.
        cellpixels_at(pyr, level, scales) <= cellpixels * LEVEL_MIDPOINT && break

        next = similar(cells, 0)
        for c in cells
            append!(next, DGG.children(pyr.system, c))
        end
        finer, finerscales = keepvisible(pyr, view, next, ntasks)
        (isempty(finer) || length(finer) > maxcells) && break
        cells, scales, level = finer, finerscales, level + 1
    end

    return cells, level
end

"""
    keepvisible(pyramid, view, candidates, ntasks) -> (cells, scales)

The candidates whose covering caps are on screen and can hold data, each with
the scale its cap came out at — pixels per radian, measured where the user is
looking.

This is the whole cost of a descent — every candidate is tested here, and a
level has seven times as many as the one above it — and the test on one cell
says nothing about the test on another, so it is done in parallel.  Chunks are
merged in order, so the answer does not depend on how many threads ran it.
"""
function keepvisible(pyr::CellPyramid, view::ScreenView, candidates::AbstractVector,
        ntasks::Integer = 1)
    n = length(candidates)
    n == 0 && return similar(candidates, 0), Float64[]

    nt = clamp(Int(ntasks), 1, max(1, cld(n, 2048)))
    if nt == 1
        return _keepvisible(pyr, view, candidates, 1, n)
    end

    bounds = round.(Int, range(0, n; length = nt + 1))
    parts = Vector{Tuple{typeof(candidates), Vector{Float64}}}(undef, nt)
    Threads.@sync for t in 1:nt
        Threads.@spawn parts[t] = _keepvisible(pyr, view, candidates,
            bounds[t] + 1, bounds[t + 1])
    end

    cells = similar(candidates, 0)
    scales = Float64[]
    for (c, s) in parts
        append!(cells, c)
        append!(scales, s)
    end
    return cells, scales
end

function _keepvisible(pyr::CellPyramid, view::ScreenView, candidates::AbstractVector,
        lo::Int, hi::Int)
    cells = similar(candidates, 0)
    scales = Float64[]
    for i in lo:hi
        c = candidates[i]
        cap = DGG.node_extent(pyr.system, c)
        _angle(pyr.capcentre, cap.point) <= cap.radius + pyr.capradius || continue
        _, scale, visible = capdisc(view, cap)
        visible || continue
        push!(cells, c)
        push!(scales, scale)
    end
    return cells, scales
end

"""
    cellpixels_at(pyramid, level, scales) -> Float64

How wide a cell of `level` comes out on screen, in pixels.

The caps the descent measures are inflated covers rather than cells, so what is
taken from them is only the scale — pixels per radian — and the width itself
comes from the level's own cell size.  The median is what makes that robust: a
handful of caps land across a projection's cut meridian and project to nonsense.
"""
function cellpixels_at(pyr::CellPyramid, level::Integer, scales::Vector{Float64})
    isempty(scales) && return 0.0
    scale = length(scales) == 1 ? scales[1] :
        partialsort(copy(scales), (length(scales) + 1) ÷ 2)
    return scale * levelsize(pyr, level) / EARTH_RADIUS
end

"""
    levelsize(pyramid, level) -> Float64

The typical width of a cell of `level`, in metres, remembered between frames
because a level's size is a property of the system and a view asks for it every
time it changes.
"""
levelsize(pyr::CellPyramid, level::Integer) =
    get!(() -> DGG.cellsize(pyr.system, level; samples = 64),
        LEVEL_SIZES, (nameof(typeof(pyr.system)), Int(level)))

const LEVEL_SIZES = Dict{Tuple{Symbol, Int}, Float64}()

# Halfway between two levels, in the ratio sense.  IGEO7 refines sevenfold in
# area, so a level's cells are √7 times narrower than its parent's, and the
# level nearest a target width is the one whose width is within √(√7) of it.
const LEVEL_MIDPOINT = 7^0.25
