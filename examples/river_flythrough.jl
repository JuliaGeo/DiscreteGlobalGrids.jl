# Fly a camera down a river valley over an IGEO7 DEM.
#
# The terrain is the Copernicus 30 m DEM of the Alps regridded onto IGEO7 and
# drawn by `dggsurface` — one vertex per hexagon, lifted to its own elevation —
# so the surface in the video is the cell geometry itself, not a resampled
# raster. It is coloured by a geomorphometric indicator computed on the DGGS
# neighbourhood, and the flight follows the valley's own main stem: the camera
# path is traced from the flow accumulation field, upstream from the outlet,
# so the route is derived from the data rather than hand-placed.
#
# The scene is built in real coordinates. `dggsurface` reads the plot's
# `transform_func`, so an orthographic PROJ transform centred on the tile puts
# x/y in kilometres on the ground while z stays elevation in metres — the same
# construction as `docs/src/tutorials/hydrology.jl`, in kilometres so a camera
# a few hundred metres off the deck keeps its depth precision.
#
# Run in the docs environment, which has GLMakie, Geomorphometry and
# FlyThroughPaths:
#
#     julia -t auto --project=docs examples/river_flythrough.jl
#
# Everything is configurable, by `--flag value` or by the matching environment
# variable. The defaults render the Vinschgau reach of the Adige at native
# resolution, which is ~1.6M cells:
#
#     --indicator NAME   what colours the terrain (see INDICATORS below)
#     --level N          IGEO7 level; `auto` matches the raster's cell size
#     --extent W,S,E,N   lon/lat box to fly, in degrees
#     --exaggeration X   vertical exaggeration of the relief
#     --seconds S        flight duration
#     --fps N            frame rate
#     --width / --height video size in pixels
#     --output PATH      output file
#     --cache PATH       file holding the regridded DEM values; written on the
#                        first run and reused after, so changing only the
#                        indicator or the camera skips the regrid entirely
#     --presmooth N      neighbour-averaging passes over the DEM before use
#     --simplify-tol KM  how much of a corner the flight path keeps
#     --palette NAME     the whole look: colormap, sky and water (see PALETTES)
#     --colormap NAME    the terrain colormap on its own
#     --sky B / --haze X  the gradient behind the scene, and how far into it
#                        the distant ridges dissolve
#     --streamlines B    run water down the channel network
#     --flow-style NAME  `dots`, `comets` or `pulses`
#     --stream-threshold N  upstream cells a channel is taken to start at
#     --flow-speed KM    how fast the water runs on the trunk; `auto` sets it
#                        against the speed of the flight. Everywhere else the
#                        speed is Manning's, off the reach slope
#     --stills           write three PNGs instead of a video, to frame a shot
#
# The video has no text, no axes and no colourbar — it is only the flight.

import DiscreteGlobalGrids as DGG
import Geomorphometry as GM
import GeoInterface as GI
import GeometryOps as GO
import ArchGDAL
import DimensionalData
import Extents
using Rasters, RasterDataSources
using GLMakie, GeoMakie
using DiscreteGlobalGridsVisualization: dggsurface!
using FlyThroughPaths
using LinearAlgebra, Statistics

const US = GO.UnitSpherical

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# `--flag value` on the command line, else `$FLAG` in the environment, else the
# default. One lookup so a flag and its variable can never drift apart.
const ARGS_DICT = let d = Dict{String,String}()
    i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        if startswith(a, "--")
            key = replace(a[3:end], '-' => '_')
            if occursin('=', key)
                k, v = split(key, '=', limit = 2)
                d[k] = v
            elseif i < length(ARGS) && !startswith(ARGS[i + 1], "--")
                d[key] = ARGS[i + 1]; i += 1
            else
                d[key] = "true"
            end
        end
        i += 1
    end
    d
end

function setting(name::String, default)
    raw = get(ARGS_DICT, name, get(ENV, uppercase(name), nothing))
    raw === nothing && return default
    return _parsesetting(default, raw)
end
_parsesetting(::AbstractString, raw) = raw
_parsesetting(::Symbol, raw) = Symbol(raw)
_parsesetting(::Bool, raw) = raw in ("1", "true", "yes")
_parsesetting(d::Integer, raw) = parse(typeof(d), raw)
_parsesetting(d::AbstractFloat, raw) = parse(typeof(d), raw)
_parsesetting(::Extents.Extent, raw) = begin
    w, s, e, n = parse.(Float64, split(raw, ','))
    Extents.Extent(X = (w, e), Y = (s, n))
end

# The Vinschgau: the Adige runs the length of it, west to east, with the Ortler
# group to the south. One Copernicus tile covers it.
const EXTENT = setting("extent", Extents.Extent(X = (10.45, 10.95), Y = (46.58, 46.78)))
const LEVEL = setting("level", "auto")
const EXAGGERATION = setting("exaggeration", 2.5)
const INDICATOR = setting("indicator", :flowaccumulation)
const SECONDS = setting("seconds", 24.0)
const FPS = setting("fps", 30)
const WIDTH = setting("width", 1920)
const HEIGHT = setting("height", 1080)
const BACKGROUND = setting("background", "black")
const STILLS = setting("stills", false)
const CACHE = setting("cache", "")
const OUTPUT = setting("output", "")

# Camera, in real units. Heights are metres above the river bed and are
# exaggerated with the relief, so the framing holds when `--exaggeration`
# changes; the horizontal distances are ground kilometres and do not.
const HEIGHT_START = setting("height_start", 1600.0)  # m above the bed, at the source
const HEIGHT_END = setting("height_end", 640.0)       # m above the bed, at the mouth
const TRAIL = setting("trail", 2.0)                   # km the eye sits behind its subject
const LOOKAHEAD = setting("lookahead", 5.5)           # km down-valley the eye looks
const FOV = setting("fov", 50.0)                      # degrees
const TRIM_HEAD = setting("trim_head", 0.10)          # skip the headwater cirque
const TRIM_TAIL = setting("trim_tail", 0.02)          # skip the stub at the outlet
const SMOOTHING = setting("smoothing", 12)            # thalweg moving-average half-width
const SIMPLIFY_TOL = setting("simplify_tol", 0.25)    # km a corner must be worth keeping
const SAMPLES = setting("samples", 160)               # camera waypoints along the track
const HOLD_START = setting("hold_start", 0.0)         # seconds held before the flight
const HOLD_END = setting("hold_end", 0.0)             # seconds held after it
# Hold the camera at one second of the flight while the water keeps running:
# the only way to see what the water is doing without the terrain sliding past.
const FREEZE = setting("freeze", -1.0)

# Neighbour-averaging passes over the elevation before anything reads it. What
# they take off is the surface rather than the terrain: Copernicus is a DSM, so
# across cells 24.7 m wide every one of them carries canopy, buildings and the
# resampling of a posting no finer than itself. Two passes halve the cell-to-cell
# roughness — median |z - ring mean| 1.02 m to 0.47 m, 4.15 m to 1.88 m at the
# 95th. Zero leaves the DEM exactly as it was regridded.
const PRESMOOTH = setting("presmooth", 2)

# `height_above_nearest_drainage` needs to be told what counts as a channel.
const HAND_THRESHOLD = setting("hand_threshold", 100)

# ---------------------------------------------------------------------------
# The look
# ---------------------------------------------------------------------------

"""
    Palette

One coherent look: what the terrain is coloured by, what is behind it, and what
the water is.

Sky and water cannot be chosen independently of the terrain — water the colour
of the bright end of the colormap disappears into the channel it is running
down, and a sky that does not share the terrain's cast reads as a hole cut in
the frame. So they travel together, and `--palette` swaps all four at once.
Each part is still a flag of its own for anyone who wants to break the set.
"""
struct Palette
    terrain::Symbol         # colormap; `:auto` keeps whatever the indicator uses
    zenith::String
    horizon::String
    dot::String
    halo::String
    accent::String          # the pulse: a hue the terrain colormap does not own
end

const PALETTES = Dict{Symbol,Palette}(
    :ice => Palette(:devon, "#050a18", "#3d5f8c", "#eaf8ff", "#03101c", "#ffd27f"),
    :ember => Palette(:lipari, "#180c16", "#c66a45", "#dcf6ff", "#170a10", "#7fe4ff"),
    :glacier => Palette(:davos, "#071320", "#7fa6bd", "#ffffff", "#04131f", "#ff9ecb"),
    :dusk => Palette(:acton, "#120a1e", "#a67fae", "#e6fbff", "#0a0614", "#ffd98a"),
    :auto => Palette(:auto, "#05070d", "#33405e", "#eaf8ff", "#03101c", "#ffd27f"),
)

const PALETTE = let name = Symbol(setting("palette", "ice"))
    haskey(PALETTES, name) || error("unknown palette $name; one of \
        $(join(sort(collect(keys(PALETTES))), ", "))")
    PALETTES[name]
end

# The terrain colormap: the palette's, unless the palette defers to the
# indicator's own, and `--colormap` beats both.
const COLORMAP = setting("colormap",
    PALETTE.terrain === :auto ? "" : string(PALETTE.terrain))

# A hard edge between a ridge and a flat background is the one thing that says
# "render" loudest, so the background is a gradient and the terrain dissolves
# into it with distance. `HAZE` is how much of the horizon colour the furthest
# ridge takes on; between `HAZE_NEAR` and `HAZE_FAR` kilometres it ramps.
const SKY = setting("sky", true)
const SKY_TOP = setting("sky_top", PALETTE.zenith)
const SKY_HORIZON = setting("sky_horizon", PALETTE.horizon)
const HAZE = setting("haze", 0.82)
const HAZE_NEAR = setting("haze_near", 2.0)
const HAZE_FAR = setting("haze_far", 24.0)
const HAZE_POWER = setting("haze_power", 1.3)

# The water. Points run down the channel network at the river's own speed, and
# what sets that speed is the gradient they are falling down, not how much of
# them there is: a torrent off a headwall runs faster than the Adige does on
# its own floodplain, however much more water the Adige is carrying. Manning
# has velocity going as the square root of the slope, which is `SLOPE_EXPONENT`
# — his other term, the one that grows with depth, is `FLOW_EXPONENT` and is
# off by default, because at its real value of about 0.27 it very nearly
# cancels the slope and rivers come out running at one speed everywhere, which
# is true and shows nothing.
#
# Since the dots are spaced by travel time rather than by distance, slowing the
# trunk crowds them along it: the big river ends up with more of them, larger
# (they are sized by flow) and moving slower, and the headwaters with a few
# small fast ones. That is the potential turning into velocity, drawn.
#
# They are objects in the scene rather than an overlay: occluded by the
# terrain, and shrinking with distance like everything else in it.
#
# `--flow-style` picks how they are drawn: `dots` are beads, `comets` give each
# a tail so the line reads as motion even in a still, and `pulses` drops the
# beads for the channels themselves with a bright wave running down them.
const FLOW_STYLE = Symbol(setting("flow_style", "dots"))
const STREAMLINES = setting("streamlines", true)
const STREAM_THRESHOLD = setting("stream_threshold", 600.0) # upstream cells a channel needs
const STREAM_BRANCHES = setting("stream_branches", 500)     # 0 keeps every branch
const STREAM_MIN_CELLS = setting("stream_min_cells", 8)     # shorter branches are dropped
const STREAM_LINES = setting("stream_lines", false)         # draw the channels themselves
const FLOW_SPEED = setting("flow_speed", "auto")      # km per second, at the trunk
const SLOPE_EXPONENT = setting("slope_exponent", 0.5) # Manning: speed goes as sqrt(slope)
const FLOW_EXPONENT = setting("flow_exponent", 0.0)   # Manning's depth term; 0.27 is real
const SLOPE_WINDOW = setting("slope_window", 20)      # cells either side of a reach slope
const SPEED_RANGE = setting("speed_range", 1.8)       # the fastest and slowest, over the trunk
const DOT_GAP = setting("dot_gap", 0.25)              # km between dots, at the trunk
const DOT_GAP_MAX = setting("dot_gap_max", 0.35)      # and the widest they ever get
const DOT_RADIUS = setting("dot_radius", 25.0)        # m across the ground, at the trunk
const DOT_LIFT = setting("dot_lift", 8.0)             # m above the bed, to clear the mesh
const DOT_FADE = setting("dot_fade", 0.5)             # seconds a dot takes to appear
const DOT_COLOR = setting("dot_color", PALETTE.dot)
const DOT_HALO = setting("dot_halo", 2.0)             # dark rim, as a multiple of the dot
const HALO_COLOR = setting("halo_color", PALETTE.halo)
const MAX_DOTS = setting("max_dots", 60000)
const COMET_TAIL = setting("comet_tail", 5)           # ghosts trailing each dot
const COMET_GAP = setting("comet_gap", 0.2)           # their spacing, as a share of the dots'
const LINE_WIDTH = setting("line_width", 1.8)         # px, where the river is biggest
const PULSE_PERIOD = setting("pulse_period", 1.8)     # seconds between one wave and the next
const PULSE_WIDTH = setting("pulse_width", 0.16)      # how much of the period a wave fills
const PULSE_SWELL = setting("pulse_swell", 3.0)       # how much wider a wave makes the line
const PULSE_COLOR = setting("pulse_color", PALETTE.accent)

# ---------------------------------------------------------------------------
# The indicators
# ---------------------------------------------------------------------------

# Everything Geomorphometry computes from a DGGS neighbourhood rather than from
# a square raster window, which is what an IGEO7 cell axis offers it. `scale`
# is applied before the colour range is taken, and `symmetric` centres a
# diverging map on zero instead of on the data.
struct Indicator
    compute::Function     # (elevation, accumulation) -> Vector
    colormap::Symbol
    scale::Function
    symmetric::Bool
end

const INDICATORS = Dict{Symbol,Indicator}(
    # Upstream area, in cell equivalents: the drainage network picks itself out.
    :flowaccumulation => Indicator(
        (elev, acc) -> vec(parent(acc)) ./ GM.cellarea(elev, first(eachindex(elev))),
        :devon, log10, false),
    :elevation => Indicator((elev, acc) -> vec(parent(elev)), :bone, identity, false),
    # Maximum downward gradient — the only slope method that reads a DGGS ring.
    :slope => Indicator((elev, acc) -> vec(parent(GM.slope(elev; method = GM.MDG()))),
        :magma, identity, false),
    :tpi => Indicator((elev, acc) -> vec(parent(GM.topographic_position_index(elev))),
        :balance, identity, true),
    :roughness => Indicator((elev, acc) -> vec(parent(GM.roughness(elev))),
        :lajolla, identity, false),
    :tri => Indicator((elev, acc) -> vec(parent(GM.terrain_ruggedness_index(elev))),
        :acton, identity, false),
    :twi => Indicator(
        (elev, acc) -> vec(parent(GM.topographic_wetness_index(elev; method = GM.D8()))),
        :haline, identity, false),
    :spi => Indicator(
        (elev, acc) -> vec(parent(GM.stream_power_index(elev; method = GM.D8()))),
        :dense, identity, false),
    :drainage_potential => Indicator(
        (elev, acc) -> vec(parent(GM.drainage_potential(elev; method = GM.D8()))),
        :tempo, identity, false),
    # Height above the nearest channel: the valley floor reads flat and dark,
    # the walls bright. IGEO7 only — it needs relative-cell arithmetic.
    :hand => Indicator(
        (elev, acc) -> vec(parent(GM.height_above_nearest_drainage(elev;
            method = GM.D8(), threshold = HAND_THRESHOLD))),
        :oslo, identity, false),
    :prominence => Indicator((elev, acc) -> Float32.(vec(parent(GM.prominence(elev)))),
        :berlin, identity, false),
    :depression_depth => Indicator(
        (elev, acc) -> vec(parent(GM.depression_depth(elev))), :turbid, identity, false),
)

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

"""
    elevation_cube() -> Raster

The DEM over `EXTENT`, regridded onto IGEO7 and returned over a `Cells` axis.

Regridding a level-13 tile is by far the slowest step here, and nothing about
it depends on the indicator or the camera, so `CACHE` names a file to keep the
result in. Only the values go in it: the axis is rebuilt by re-running the
`MultiOrderCoverage` query, which is a second's work and is a pure function of
the extent and the level, and a rebuilt axis is a plain `CellLookup` — which is
what `mapneighbors`, Geomorphometry and `dggsurface` all take.
"""
function elevation_cube()
    sys = DGG.IGeo7System()
    cached = cached_values()
    level = cached === nothing ?
        (LEVEL == "auto" ? nothing : parse(Int, LEVEL)) : cached[1]

    # With cached values and a known level the raster is never touched.
    dem = nothing
    if level === nothing || cached === nothing
        # AWS's `copernicus-dem-30m` bucket: the DGED packaging of GLO-30 as
        # COGs, so elevations are `Float32` and continuous — nothing downstream
        # is working around a quantised source.
        ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH",
            joinpath(tempdir(), "rasterdatasources")))
        centre = GI.extent(((EXTENT.X[1] + EXTENT.X[2]) / 2, (EXTENT.Y[1] + EXTENT.Y[2]) / 2))
        path = only(skipmissing(RasterDataSources.getraster(CopernicusDEM; extent = centre)))
        Sys.isapple() && Rasters.checkmem!(false)
        dem = view(Raster(path; lazy = false), EXTENT)
        level === nothing && (level = DGG.levelfor(sys, dem))
    end

    region = DGG.query(sys, DGG.MultiOrderCoverage(EXTENT); level)
    axis = DGG.Cells(DGG.CellLookup(region))
    ncells = length(axis)

    if cached !== nothing
        values = cached[2]
        length(values) == ncells ||
            error("cache $CACHE holds $(length(values)) cells but the level-$level \
                coverage of this extent has $ncells; delete it or pass another --extent")
        @info "reusing the regridded DEM" file = CACHE level cells = ncells
        return Raster(DimensionalData.DimArray(values, axis))
    end

    @info "regridding" level cells = ncells
    cube = DGG.regrid(dem; to = region)
    isempty(CACHE) || write_cache(level, Float32.(vec(parent(cube))))
    return Raster(DimensionalData.DimArray(Float32.(vec(parent(cube))), axis))
end

# The cache is a level and a run of `Float32`s — the axis is not in it, because
# the query rebuilds it and a stored axis reads back chunked.
const CACHE_MAGIC = 0x44474746   # "DGGF"

function cached_values()
    (isempty(CACHE) || !isfile(CACHE)) && return nothing
    return open(CACHE, "r") do io
        read(io, UInt32) == CACHE_MAGIC || return nothing
        level = Int(read(io, Int64))
        n = Int(read(io, Int64))
        values = Vector{Float32}(undef, n)
        read!(io, values)
        (level, values)
    end
end

function write_cache(level, values)
    open(CACHE, "w") do io
        write(io, CACHE_MAGIC, Int64(level), Int64(length(values)), values)
    end
    @info "cached the regridded DEM" file = CACHE MB = round(
        filesize(CACHE) / 2^20, digits = 1)
end

"""
    fill_gaps!(elev) -> Raster

Replace the elevations the coverage did not reach with the mean of their
finite neighbours, and return a cube with no `NaN` left.

The IGEO7 coverage of a lon/lat box overhangs it slightly, so the outermost
ring of cells regrids to `NaN`. A hole in a value map is honest; a hole in a
surface being flown through is a tear in the geometry.
"""
function fill_gaps!(elev)
    z = vec(parent(elev))
    for _ in 1:8
        count(isnan, z) == 0 && break
        patched = DGG.mapneighbors(elev) do cell, nbrs
            here = z[DGG.localindex(cell)]
            isnan(here) || return here
            total, seen = 0.0, 0
            for nb in nbrs
                zn = z[DGG.localindex(nb)]
                isnan(zn) || (total += zn; seen += 1)
            end
            (seen == 0 ? NaN32 : Float32(total / seen))::Float32
        end
        z .= vec(parent(patched))
    end
    return rebuild(elev; data = z)
end

"""
    presmooth(elev) -> Raster

`elev` averaged with its neighbours, `PRESMOOTH` times over.

Each pass replaces a cell with the mean of itself and its one-ring, which on a
hexagonal grid is an isotropic low-pass filter — there is no axis for it to
smear along, which is the reason a square-window mean would need a shape.
"""
function presmooth(elev)
    PRESMOOTH <= 0 && return elev
    z = vec(parent(elev))
    for _ in 1:PRESMOOTH
        averaged = DGG.mapneighbors(elev) do cell, nbrs
            total, seen = Float64(z[DGG.localindex(cell)]), 1
            for nb in nbrs
                total += z[DGG.localindex(nb)]; seen += 1
            end
            Float32(total / seen)::Float32
        end
        z = vec(parent(averaged))
    end
    @info "pre-smoothed the relief" passes = PRESMOOTH
    return rebuild(elev; data = z)
end

# ---------------------------------------------------------------------------
# The route
# ---------------------------------------------------------------------------

"""
    thalweg(elevation, accumulation) -> Vector{Int}

Cell indices along the valley's main stem, from the headwater to the outlet.

Walks upstream from the cell of greatest flow accumulation, taking at each step
the neighbour carrying the most water of those carrying less than the cell
itself. Accumulation is the field to climb rather than elevation: it strictly
decreases upstream and, unlike elevation, has no flats, so the walk crosses a
braided valley floor — where no neighbour is any lower — without stalling.
"""
function thalweg(elev, acc)
    a = Float64.(vec(parent(acc)))
    z = vec(parent(elev))
    up = vec(parent(DGG.mapneighbors(elev) do cell, nbrs
        here = a[DGG.localindex(cell)]
        best, bestflow = 0, -Inf
        for nb in nbrs
            j = DGG.localindex(nb)
            if a[j] < here && a[j] > bestflow
                best, bestflow = j, a[j]
            end
        end
        best::Int
    end))

    outlet = argmax(a)
    floor_ = 0.005 * a[outlet]     # stop before the trace frays into hillslopes
    chain = [outlet]
    visited = falses(length(a))
    while true
        here = chain[end]
        visited[here] = true
        next = up[here]
        (next == 0 || visited[next] || a[next] < floor_ || isnan(z[next])) && break
        push!(chain, next)
    end
    return reverse!(chain)         # headwater first, so the flight runs downstream
end

"""
    projector(elev, grid, project) -> (index, lift = 0) -> Point3d

A function from a cell's index in the cube to its centre in the scene.

Everything that is not the terrain — the flight path, the water — is built as
plain geometry rather than as a plot over cells, so each of them has to land in
the same coordinates the surface does: the tile's orthographic projection in
kilometres, carrying elevation in metres through the same exaggeration. `lift`
is metres above the ground, for the things that must not sit inside it.
"""
function projector(elev, grid, project)
    z = vec(parent(elev))
    cells = parent(DimensionalData.lookup(elev, DGG.Cells))
    geographic = US.GeographicFromUnitSphere()
    return function (p, lift = 0.0)
        lonlat = geographic(DGG.cell_centroid(grid, cells[p]))
        x, y = project((lonlat[1], lonlat[2]))
        return Point3d(x, y, (z[p] + lift) * EXAGGERATION / 1000)
    end
end

"""
    track(indices, place) -> Vector{Point3d}

The thalweg as a smooth, equally spaced polyline in scene coordinates.

Cell centres zigzag between hexagons, which a camera following them reads as
a shudder, so the line is averaged along its length before it is resampled to
constant spacing — constant spacing is what makes constant-speed interpolation
between waypoints a constant ground speed.
"""
function track(indices, place)
    raw = map(place, indices)

    n = length(raw)
    w = min(SMOOTHING, max(0, (n - 1) ÷ 2))
    smoothed = [mean(@view raw[max(1, i - w):min(n, i + w)]) for i in 1:n]

    along = cumsum([0.0; [norm(smoothed[i + 1] - smoothed[i]) for i in 1:n - 1]])
    return map(range(0, along[end], length = SAMPLES)) do d
        i = clamp(searchsortedlast(along, d), 1, n - 1)
        u = (d - along[i]) / max(along[i + 1] - along[i], eps())
        smoothed[i] + u * (smoothed[i + 1] - smoothed[i])
    end
end

"""
    flightpath(track) -> FlyThroughPaths.Path

A camera path down `track`, curved through its corners.

The eye trails its subject and looks down-valley at a point ahead of it, which
is what makes the terrain open out in front of the camera rather than sweep
past it. It starts high enough to see over the headwall and descends into the
valley as the flight goes on.

A move per sample would be a polyline, and a polyline is only C⁰: the heading
changes discontinuously at every waypoint and the flight reads as a stutter,
one tick per cell the thalweg stepped through. So the track is first reduced to
the corners that actually matter — `simplify` keeps the vertices whose removal
would move the line, and drops the rest — and consecutive corners are joined by
a `BezierMove` whose control points are the Catmull-Rom tangents. Sharing a
tangent is what makes the joins C¹, and C¹ is what the eye reads as smooth.

Segments are timed by their length rather than given an equal share, so
straightening the line does not also speed the camera up along it.
"""
function flightpath(tr)
    k = length(tr)
    rise(s) = HEIGHT_END + (HEIGHT_START - HEIGHT_END) * (1 - smoothstep(s))^2
    spacing = norm(tr[end] - tr[1]) / max(k - 1, 1)
    ahead = max(1, round(Int, LOOKAHEAD / max(spacing, eps())))

    # The eye and its target at sample `i`, as a pair of points. They are curved
    # separately and reassembled, so the camera turns as smoothly as it travels.
    function station(i)
        s = (i - 1) / max(k - 1, 1)
        here = tr[i]
        target = tr[clamp(i + ahead, 1, k)]
        forward = target - here
        forward = norm(forward) < 1e-9 ? Vec3d(1, 0, 0) : normalize(forward)
        lift = rise(s) * EXAGGERATION / 1000
        return (eye = Vec3d(here - TRAIL * forward + Vec3d(0, 0, lift)),
                at = Vec3d(target + Vec3d(0, 0, 0.3 * lift)))
    end

    knots = corners(tr)
    stations = [station(i) for i in knots]
    state(st) = ViewState(eyeposition = st.eye, lookat = st.at,
        upvector = Vec3d(0, 0, 1), fov = FOV)

    n = length(stations)
    eyes = [st.eye for st in stations]
    ats = [st.at for st in stations]
    spans = [max(norm(eyes[i + 1] - eyes[i]), eps()) for i in 1:n - 1]
    moving = SECONDS - HOLD_START - HOLD_END
    times = [moving * s / sum(spans) for s in spans]

    # Catmull-Rom, in seconds rather than in index. The velocity at a knot is
    # the chord through its neighbours over the time that chord takes, and a
    # cubic Bézier whose control arm is a third of `velocity * duration` leaves
    # each segment at exactly that velocity. Matching them is what makes the
    # joins C¹ *in time*; arms scaled in index space are C¹ only where adjacent
    # segments last equally long, which after `simplify` is nowhere.
    function velocity(v, i)
        i == 1 && return (v[2] - v[1]) / times[1]
        i == n && return (v[n] - v[n - 1]) / times[n - 1]
        return (v[i + 1] - v[i - 1]) / (times[i - 1] + times[i])
    end

    control(eye, at) = ViewState(eyeposition = eye, lookat = at,
        upvector = Vec3d(0, 0, 1), fov = FOV)

    path = Path(state(stations[1]))
    HOLD_START > 0 && (path = path * Pause(HOLD_START))
    for i in 1:n - 1
        t = times[i]
        arm = t / 3
        path = path * BezierMove(t, state(stations[i + 1]), [
            control(eyes[i] + velocity(eyes, i) * arm,
                ats[i] + velocity(ats, i) * arm),
            control(eyes[i + 1] - velocity(eyes, i + 1) * arm,
                ats[i + 1] - velocity(ats, i + 1) * arm),
        ])
    end
    HOLD_END > 0 && (path = path * Pause(HOLD_END))
    @info "path" corners = n moves = n - 1
    return path
end

"""
    corners(track) -> Vector{Int}

The indices of `track` that survive a Douglas-Peucker simplification of it.

`simplify` returns geometry rather than indices and drops the third coordinate,
so the run is done on the ground track and the survivors are matched back by
index — they are verbatim copies of the vertices that went in.
"""
function corners(tr)
    line = GI.LineString([(p[1], p[2]) for p in tr])
    kept = GO.simplify(GO.DouglasPeucker(; tol = SIMPLIFY_TOL), line)
    index = Dict((p[1], p[2]) => i for (i, p) in enumerate(tr))
    out = unique(sort([index[(GI.x(p), GI.y(p))] for p in GI.getpoint(kept)]))
    # Two corners cannot be curved, and a Bézier needs a neighbour on each side
    # to take a tangent from; fall back to the ends plus the middle.
    length(out) >= 3 && return out
    return unique([1, (1 + length(tr)) ÷ 2, length(tr)])
end

smoothstep(x) = (y = clamp(x, 0.0, 1.0); y * y * (3 - 2y))

"""
    warp(t) -> Float64

`t` bent so the flight eases in and out once, over the whole of it.

Easing per move would put a beat at every corner, so the moves all run at
constant speed and the clock is bent instead. The holds are outside the ease
and pass through untouched — an eased clock over the whole duration would map
the flight onto the moving part and never sample the pauses at all.
"""
function warp(t)
    moving = SECONDS - HOLD_START - HOLD_END
    moving <= 0 && return t
    t <= HOLD_START && return t
    t >= SECONDS - HOLD_END && return t
    return HOLD_START + smoothstep((t - HOLD_START) / moving) * moving
end

# ---------------------------------------------------------------------------
# The water
# ---------------------------------------------------------------------------

"""
    downstream(flow, elev) -> Vector{Int}

For every cell, the index of the neighbour its water leaves by, or `0` where
none does.

A cell drains to the neighbour carrying the most water: accumulation at a cell
counts everything upstream of it, so the cell below is larger than the cell
itself and larger than any of its own tributaries. Reading the routing back out
of the accumulation field instead of out of the D8 direction codes keeps this
free of the convention a system encodes directions in, and makes the map
acyclic by construction — flow strictly increases at every step, so no trace
can return to a cell it has left.
"""
function downstream(flow, elev)
    return vec(parent(DGG.mapneighbors(elev) do cell, nbrs
        best, most = 0, flow[DGG.localindex(cell)]
        for nb in nbrs
            j = DGG.localindex(nb)
            f = flow[j]
            if isfinite(f) && f > most
                best, most = j, f
            end
        end
        best::Int
    end))
end

"""
    upstream(down) -> (starts, kids)

`down` inverted: the cells draining into each cell, as one flat list and the
index its run begins at.

A cell has one downstream neighbour and any number of upstream ones, so the
inverse is a ragged list — held flat, because a vector of vectors over a couple
of million cells is a couple of million allocations.
"""
function upstream(down)
    n = length(down)
    starts = ones(Int, n + 1)
    for d in down
        d == 0 || (starts[d + 1] += 1)
    end
    for i in 2:n + 1
        starts[i] += starts[i - 1] - 1
    end
    at = copy(starts)
    kids = Vector{Int}(undef, starts[n + 1] - 1)
    for (i, d) in enumerate(down)
        d == 0 && continue
        kids[at[d]] = i
        at[d] += 1
    end
    return starts, kids
end

"""
    channels(flow, down, starts, kids) -> Vector{Vector{Int}}

The channel network as a set of branches, each running downstream.

Cells carrying less than `STREAM_THRESHOLD` cells' worth of water are hillslope
and are left out; what remains is walked from the bottom up, the same way the
flight's own route is. Starting at the cell with the most water and climbing
into its largest tributary at every fork traces the main stem from the outlet
to its headwater; the next branch starts at the largest cell that walk did not
reach, which is the mouth of the biggest tributary, and so on until the network
is used up. Every channel cell therefore lands on exactly one branch, and each
branch is prefixed with the cell it drains into so it meets the trunk rather
than stopping one cell short of it.
"""
function channels(flow, down, starts, kids)
    network = findall(>=(STREAM_THRESHOLD), flow)
    isempty(network) && return Vector{Int}[]
    sort!(network; by = i -> flow[i], rev = true)

    taken = falses(length(down))
    branches = Vector{Int}[]
    for seed in network
        taken[seed] && continue
        branch = Int[]
        mouth = down[seed]
        (mouth != 0 && flow[mouth] >= STREAM_THRESHOLD) && push!(branch, mouth)
        node = seed
        while true
            push!(branch, node)
            taken[node] = true
            best, most = 0, -Inf
            for k in starts[node]:starts[node + 1] - 1
                up = kids[k]
                (taken[up] || flow[up] < STREAM_THRESHOLD) && continue
                if flow[up] > most
                    best, most = up, flow[up]
                end
            end
            best == 0 && break
            node = best
        end
        length(branch) >= STREAM_MIN_CELLS && push!(branches, reverse!(branch))
    end

    sort!(branches; by = length, rev = true)
    if STREAM_BRANCHES > 0 && length(branches) > STREAM_BRANCHES
        @info "kept the largest branches" found = length(branches) kept = STREAM_BRANCHES
        branches = branches[1:STREAM_BRANCHES]
    end
    return branches
end

"""
    Water

The channel network, parameterised by how long water takes to run down it.

`points`, `times` and `weights` are the branches laid end to end, `starts`
holding where each begins. `times` restarts at zero on every branch and counts
seconds of travel from its head, which is what makes the animation a lookup:
a dot is a branch and a departure time, and where it is at `t` is where that
branch is `offset + t` seconds along.
"""
struct Water
    starts::Vector{Int}
    points::Vector{Point3d}
    times::Vector{Float64}
    weights::Vector{Float32}
    branch::Vector{Int32}
    offset::Vector{Float64}
    tail::Vector{Float32}
end

"""
    water(branches, flow, place, speed) -> Water

`branches` measured out into dots, spaced by travel time rather than by
distance.

The velocity at a point is Manning's: `speed` times the reach slope over the
trunk's, to `SLOPE_EXPONENT`, optionally scaled by the flow to
`FLOW_EXPONENT`, and held inside `SPEED_RANGE` either way so that no reach
stalls and none of them bolts. `speed` is therefore what the *trunk* runs at,
not the fastest anything goes — which is the useful end to pin, because the
trunk is the reach the flight follows and the one whose pace the eye judges
everything else against.

Slope is taken over a reach rather than between neighbours: on the valley floor
the DSM's roughness across one 24.7 m cell runs to 0.3 m against the channel's
0.09 m of drop, so a cell-to-cell gradient there is surface cover and not river.
Smoothing does not rescue it — `PRESMOOTH` 2 moves the floor's median neighbour
slope from 12.8 % to 11.3 %, against a true gradient of 0.36 %. `SLOPE_WINDOW`
cells either side is what makes it a number rather than noise.

Spacing the dots evenly in duration rather than along the ground is what a
steady release of them would look like: they crowd where the water is slow and
string out where it runs, which is the difference between a river and a row of
beads. `DOT_GAP_MAX` bounds how far that stringing out is allowed to go, since
past a certain sparseness a mountain torrent stops reading as a torrent and
starts reading as nothing at all.
"""
function water(branches, flow, place, speed)
    isempty(branches) && return Water(Int[1], Point3d[], Float64[], Float32[],
        Int32[], Float64[], Float32[])
    biggest = maximum(flow)
    starts = Int[1]
    points = Point3d[]
    weights = Float32[]
    for branch in branches
        for p in branch
            push!(points, place(p, DOT_LIFT))
            push!(weights, Float32(clamp(flow[p] / biggest, 1e-6, 1.0)))
        end
        push!(starts, length(points) + 1)
    end

    # True elevation and true along-channel distance: `place` exaggerates the
    # vertical, and a slope read off the exaggerated geometry is the
    # exaggeration, not the ground.
    height = [p[3] / EXAGGERATION for p in points]
    along = zeros(Float64, length(points))
    step = zeros(Float64, length(points))
    for c in 1:length(starts) - 1
        lo, hi = starts[c], starts[c + 1] - 1
        for i in lo:hi - 1
            d = points[i + 1] - points[i]
            step[i] = sqrt(d[1]^2 + d[2]^2 + (d[3] / EXAGGERATION)^2)
            along[i + 1] = along[i] + step[i]
        end
    end

    slope = zeros(Float64, length(points))
    for c in 1:length(starts) - 1
        lo, hi = starts[c], starts[c + 1] - 1
        for j in lo:hi
            a = max(lo, j - SLOPE_WINDOW)
            b = min(hi, j + SLOPE_WINDOW)
            slope[j] = max(height[a] - height[b], 0.0) / max(along[b] - along[a], eps())
        end
    end

    # What `speed` pins is the main stem, so the reference has to be the main
    # stem's own slope and not the network's. Taking the median weighted by
    # flow does that without a threshold to pick: the trunk carries most of the
    # water, so it decides the median, and a network with no trunk degrades to
    # its own middle instead of to an empty selection.
    refslope = max(flowmedian(slope, weights), 1e-6)
    refflow = max(flowmedian(weights, weights), 1e-6)

    times = zeros(Float64, length(points))
    pace = similar(times)
    for i in eachindex(pace)
        pace[i] = clamp((max(slope[i], 1e-9) / refslope)^SLOPE_EXPONENT *
            (weights[i] / refflow)^FLOW_EXPONENT, 1 / SPEED_RANGE, SPEED_RANGE)
    end
    # A reference slope is only ever approximately the trunk's, so pinning the
    # trunk to `speed` is done by measuring where the field actually came out
    # and scaling it there. Without this `--flow-speed` is a number the water
    # is near, which is no use when what it has to beat is the camera. The
    # cells it is measured over are the main stem proper — the reach the flight
    # actually follows — and not the whole network by weight, which the smaller
    # and much faster tributaries drag upward.
    stem = findall(>=(0.3f0), weights)
    pace .*= speed / max(isempty(stem) ? flowmedian(pace, weights) :
        median(@view pace[stem]), eps())
    for c in 1:length(starts) - 1
        lo, hi = starts[c], starts[c + 1] - 1
        for i in lo:hi - 1
            times[i + 1] = times[i] + step[i] / max((pace[i] + pace[i + 1]) / 2, eps())
        end
    end
    @info "the river's pace" trunk_km_per_s = round(speed, digits = 2) fastest = round(
        maximum(pace), digits = 2) slowest = round(minimum(pace), digits = 2)

    # `pulses` lights the channels themselves, so it needs the geometry and
    # none of the dots.
    FLOW_STYLE === :pulses &&
        return Water(starts, points, times, weights, Int32[], Float64[], Float32[])

    ghosts = FLOW_STYLE === :comets ? COMET_TAIL : 0
    total = sum(times[starts[c + 1] - 1] for c in eachindex(branches))
    spacing = DOT_GAP / speed
    if total / spacing * (ghosts + 1) > MAX_DOTS
        spacing = total * (ghosts + 1) / MAX_DOTS
        @info "thinned the water to stay under --max-dots" dot_gap = round(
            spacing * speed, digits = 2)
    end

    branch_of = Int32[]
    offset = Float64[]
    tail = Float32[]
    for c in eachindex(branches)
        lo, hi = starts[c], starts[c + 1] - 1
        span = times[hi]
        span <= 0 && continue
        # One release interval for the whole network puts the dots `DOT_GAP`
        # apart on the trunk and four times that on a torrent, which is a
        # torrent nobody can follow. So each branch releases at its own rate,
        # fast enough that its dots are never further apart on the ground than
        # `DOT_GAP_MAX` — the spacing still opens out where the water runs, but
        # only so far. Within a branch the interval is constant, which is what
        # keeps the dots behaving like a steady release rather than a pattern
        # of their own travelling down it.
        gap = min(spacing, DOT_GAP_MAX / max(median(@view pace[lo:hi]), eps()))
        # Branches are dephased against each other, or every one of them would
        # release its first dot on the same frame.
        t = gap * mod(c * 0.6180339887498949, 1.0)
        while t < span
            # A comet is its head plus a few ghosts released a moment earlier,
            # which is why they trail it: they are the same dot, behind in time.
            for g in 0:ghosts
                push!(branch_of, c)
                push!(offset, t - g * COMET_GAP * gap)
                push!(tail, Float32((1 - g / (ghosts + 1))^1.5))
            end
            t += gap
        end
    end
    return Water(starts, points, times, weights, branch_of, offset, tail)
end

"""
    flowmedian(values, weights) -> eltype(values)

The median of `values` with each one counting for the water it carries.

A plain median over a channel network is a median over its headwaters, because
that is where nearly all of the cells are; weighting by flow asks instead what
the water sees, which is the trunk.
"""
function flowmedian(values, weights)
    isempty(values) && return zero(eltype(values))
    order = sortperm(values)
    half = sum(weights) / 2
    running = zero(eltype(weights))
    for i in order
        running += weights[i]
        running >= half && return values[i]
    end
    return values[order[end]]
end

# Where in `times[lo:hi]` — sorted, and a window of a longer vector that is not
# — the value `t` falls.
@inline function _seek(times, lo, hi, t)
    a, b = lo, hi
    while b - a > 1
        m = (a + b) >>> 1
        times[m] <= t ? (a = m) : (b = m)
    end
    return a
end

"""
    advance!(w, t, positions, sizes, colors, tint)

Put every dot where it is `t` seconds into the flight, and size and fade it.

Dots wrap round their branch, so water leaving the bottom of one arrives at its
head; a dot fades in and out over `DOT_FADE` at the ends, because a branch that
began or stopped mid-air would read as a glitch and not as water. The clock
here is the real one, not the eased clock the camera flies on — the river runs
at its own rate whatever the camera is doing.
"""
function advance!(w::Water, t, positions, sizes, colors, tint)
    @inbounds for k in eachindex(w.branch)
        c = w.branch[k]
        lo, hi = w.starts[c], w.starts[c + 1] - 1
        span = w.times[hi]
        along = mod(w.offset[k] + t, span)
        i = _seek(w.times, lo, hi, along)
        i = min(i, hi - 1)
        u = (along - w.times[i]) / max(w.times[i + 1] - w.times[i], eps())
        positions[k] = Point3f(w.points[i] + u * (w.points[i + 1] - w.points[i]))
        carried = w.weights[i] + u * (w.weights[i + 1] - w.weights[i])
        behind = w.tail[k]
        sizes[k] = Float32(DOT_RADIUS / 1000 * (0.45 + 0.55 * carried^0.25) * behind)
        edge = min(DOT_FADE, span / 3)
        fade = clamp(min(along, span - along) / edge, 0.0, 1.0)
        colors[k] = Makie.RGBAf(tint.r, tint.g, tint.b, tint.alpha * fade * behind)
    end
    return
end

# ---------------------------------------------------------------------------
# The scene
# ---------------------------------------------------------------------------

"""
    colorrange(values) -> Tuple

A robust range for `values`, clipped to its 2nd and 98th percentiles.

A geomorphometric field is usually long-tailed — one gorge holds the whole top
of a slope histogram — and a range taken from the extrema spends the entire
colormap on cells that are a few pixels wide.
"""
function colorrange(values, symmetric)
    finite = filter(isfinite, values)
    isempty(finite) && return (0.0, 1.0)
    lo, hi = quantile(finite, 0.02), quantile(finite, 0.98)
    lo == hi && ((lo, hi) = (lo - 1, hi + 1))
    if symmetric
        m = max(abs(lo), abs(hi))
        return (-m, m)
    end
    return (lo, hi)
end

"""
    pulses!(scene, river) -> (t -> nothing)

Draw the channels themselves, with a wave of light running down each of them.

The alternative to beads: the network is one `lines` whose width is the flow it
carries, lit steadily along its length and brightened where a wave is passing.
Because the wave's position is `times - t` — the same travel-time coordinate
the dots move on — it runs downstream at the water's speed and slows in the
headwaters exactly as they do. Only the colour buffer changes per frame.
"""
function pulses!(scene, river::Water)
    tint = Makie.RGBAf(Makie.to_color(PULSE_COLOR))
    rim = Makie.RGBAf(Makie.to_color(HALO_COLOR))
    line = Point3f[]
    width = Float32[]
    carried = Float32[]
    when = Float64[]
    for c in 1:length(river.starts) - 1
        for i in river.starts[c]:river.starts[c + 1] - 1
            push!(line, Point3f(river.points[i]))
            push!(width, Float32(LINE_WIDTH * (0.3 + 0.7 * river.weights[i]^0.3)))
            push!(carried, river.weights[i])
            push!(when, river.times[i])
        end
        # A break in the line, so one branch does not join the next.
        push!(line, Point3f(NaN32, NaN32, NaN32))
        push!(width, 1.0f0)
        push!(carried, 0.0f0)
        push!(when, 0.0)
    end

    n = length(line)
    colors = Observable(Vector{Makie.RGBAf}(undef, n))
    widths = Observable(copy(width))
    shadows = Observable(Vector{Makie.RGBAf}(undef, n))
    shadowwidths = Observable(copy(width))
    glow!(colors[], widths[], shadows[], shadowwidths[], width, carried, when,
        0.0, tint, rim)
    # The dark line goes down first and the bright one rides in it, for the same
    # reason the dots are rimmed: what the wave has to stand out against is a
    # channel the colormap has already painted its brightest.
    lines!(scene, line; color = shadows, linewidth = shadowwidths,
        transparency = true)
    lines!(scene, line; color = colors, linewidth = widths, transparency = true)
    @info "water" style = FLOW_STYLE vertices = n branches = length(
        river.starts) - 1

    return function (t)
        glow!(colors[], widths[], shadows[], shadowwidths[], width, carried,
            when, t, tint, rim)
        notify(colors); notify(widths); notify(shadows); notify(shadowwidths)
        return
    end
end

"""
    glow!(colors, widths, shadows, shadowwidths, rest, carried, when, t, tint, rim)

Light every vertex of the network for time `t`.

A wave has to differ from the channel it is on in something other than opacity.
The bright end of a flow colormap is already white, and the channel is exactly
where it lives, so a white line brightening from most-of-the-way to fully opaque
is a change nothing can see — which is what the first version of this did. So a
wave does three things at once: it swells the line to `PULSE_SWELL` times its
width, takes it from a dim thread to full `PULSE_COLOR`, and darkens the rim
under it. The swelling and the rim are what carry it; a silhouette in a hue the
colormap does not own reads against any terrain, where a brightness does not.

The wave itself is a gaussian in phase, written twice — once on each side of
the wrap — so it crosses from one period into the next without a flicker.
"""
function glow!(colors, widths, shadows, shadowwidths, rest, carried, when, t,
        tint, rim)
    @inbounds for i in eachindex(colors)
        here = carried[i]
        if here <= 0
            colors[i] = Makie.RGBAf(0, 0, 0, 0)
            shadows[i] = Makie.RGBAf(0, 0, 0, 0)
            widths[i] = rest[i]
            shadowwidths[i] = rest[i]
            continue
        end
        phase = mod(when[i] - t, PULSE_PERIOD) / PULSE_PERIOD
        wave = clamp(exp(-(phase / PULSE_WIDTH)^2) +
            exp(-((1 - phase) / PULSE_WIDTH)^2), 0.0, 1.0)
        size = here^0.3
        w = rest[i] * (1 + PULSE_SWELL * wave)
        widths[i] = Float32(w)
        shadowwidths[i] = Float32(w * 2.0 + 1.0)
        lit = 0.4 + 0.6 * wave
        colors[i] = Makie.RGBAf(tint.r * lit, tint.g * lit, tint.b * lit,
            tint.alpha * clamp(0.10 + 0.18 * size + 0.78 * wave, 0.0, 1.0))
        shadows[i] = Makie.RGBAf(rim.r, rim.g, rim.b,
            rim.alpha * clamp(0.18 + 0.55 * wave, 0.0, 1.0) * size)
    end
    return
end

"""
    flow!(scene, river) -> (t -> nothing)

Put the water in `scene`, and return the function that runs it.

The dots are `meshscatter` rather than `scatter` so they are the size of real
things — shrinking with distance and hidden by the ridge in front of them,
where fixed-pixel markers would read as fireflies laid over the terrain.
Nothing is lit: the water is a colour, not a surface. Every frame rewrites the
same three buffers, because the geometry never changes — only which part of it
each dot is on.
"""
function flow!(scene, river::Water)
    isempty(river.points) && return t -> nothing
    FLOW_STYLE === :pulses && return pulses!(scene, river)
    FLOW_STYLE in (:dots, :comets) ||
        error("unknown --flow-style $FLOW_STYLE; one of comets, dots, pulses")
    isempty(river.branch) && return t -> nothing
    tint = Makie.RGBAf(Makie.to_color(DOT_COLOR))

    if STREAM_LINES
        line = Point3f[]
        colour = Makie.RGBAf[]
        for c in 1:length(river.starts) - 1
            for i in river.starts[c]:river.starts[c + 1] - 1
                push!(line, Point3f(river.points[i]))
                push!(colour, Makie.RGBAf(tint.r, tint.g, tint.b,
                    0.10 + 0.35 * river.weights[i]^0.25))
            end
            push!(line, Point3f(NaN32, NaN32, NaN32))
            push!(colour, Makie.RGBAf(0, 0, 0, 0))
        end
        lines!(scene, line; color = colour, linewidth = 1.2, transparency = true)
    end

    n = length(river.branch)
    positions = Observable(Vector{Point3f}(undef, n))
    sizes = Observable(Vector{Float32}(undef, n))
    colors = Observable(Vector{Makie.RGBAf}(undef, n))
    advance!(river, 0.0, positions[], sizes[], colors[], tint)

    # A pale dot on a pale channel is not a dot, and the channel is exactly
    # where the bright end of a flow colormap lives, so each one is rimmed:
    # a larger, darker sphere at the same place, drawn transparently so the
    # core sits inside it rather than behind it. That reads on any colormap,
    # which a single colour chosen against one of them does not.
    if DOT_HALO > 1
        rim = Makie.RGBAf(Makie.to_color(HALO_COLOR))
        meshscatter!(scene, positions;
            markersize = map(v -> v .* Float32(DOT_HALO), sizes),
            color = map(v -> [Makie.RGBAf(rim.r, rim.g, rim.b, 0.55 * c.alpha) for c in v],
                colors),
            shading = Makie.NoShading, transparency = true)
    end
    meshscatter!(scene, positions; markersize = sizes, color = colors,
        shading = Makie.NoShading, transparency = true)
    @info "water" dots = n branches = length(river.starts) - 1

    return function (t)
        advance!(river, t, positions[], sizes[], colors[], tint)
        notify(positions); notify(sizes); notify(colors)
        return
    end
end

"""
    sky!(figure)

A vertical gradient behind the whole frame, in place of a flat background.

Four vertices in the figure's own pixel space, coloured zenith at the top and
horizon at the bottom: the scene is drawn over it without clearing, so what the
terrain does not cover is sky. A flat colour behind a ridgeline reads as a hole
cut in the frame however dark it is; a gradient reads as air.
"""
function sky!(figure)
    Makie.campixel!(figure.scene)
    zenith = Makie.RGBAf(Makie.to_color(SKY_TOP))
    horizon = Makie.RGBAf(Makie.to_color(SKY_HORIZON))
    quad = Point2f[(0, 0), (WIDTH, 0), (WIDTH, HEIGHT), (0, HEIGHT)]
    # The scene over this one does not clear the frame, so it shares this
    # plot's depth buffer: without pushing the quad to the far plane it is
    # written at the middle of the range and everything behind that — which is
    # most of the terrain — fails the depth test and never appears.
    mesh!(figure.scene, quad, [1 2 3; 1 3 4];
        color = [horizon, horizon, zenith, zenith], shading = Makie.NoShading,
        depth_shift = 1.0f0, fxaa = false)
    return
end

"""
    shade(values, colormap, colorrange) -> Vector{RGBAf}

`values` put through the colormap here rather than on the GPU.

Aerial perspective mixes each cell's colour towards the sky by how far away it
is, and a distance is not something a colorrange can express — so the mapping
is done once in Julia and the mixing happens on top of it, per frame.
"""
function shade(values, colormap, colorrange)
    grad = Makie.to_colormap(colormap)
    n = length(grad)
    lo, hi = colorrange
    span = max(hi - lo, eps())
    return map(values) do v
        u = isfinite(v) ? clamp((v - lo) / span, 0.0, 1.0) : 0.0
        @inbounds grad[clamp(round(Int, 1 + u * (n - 1)), 1, n)]
    end
end

"""
    mist!(out, base, vertices, eye, air)

Mix the terrain's colours towards `air` by how far each vertex is from `eye`.

This is the whole of the depth in the picture: without it a ridge twenty
kilometres off is exactly as saturated as the bank in front of the camera, and
the eye reads the two as the same distance. `HAZE_NEAR` is where it starts,
`HAZE_FAR` where it is full, and `HAZE` how much of the sky the furthest ridge
ends up wearing.
"""
function mist!(out, base, vertices, eye, air)
    span = max(HAZE_FAR - HAZE_NEAR, eps())
    ex, ey, ez = eye[1], eye[2], eye[3]
    Threads.@threads for i in eachindex(out)
        @inbounds begin
            p = vertices[i]
            d = sqrt((p[1] - ex)^2 + (p[2] - ey)^2 + (p[3] - ez)^2)
            f = HAZE * clamp((d - HAZE_NEAR) / span, 0.0, 1.0)^HAZE_POWER
            c = base[i]
            out[i] = Makie.RGBAf(c.r + (air.r - c.r) * f, c.g + (air.g - c.g) * f,
                c.b + (air.b - c.b) * f, c.alpha)
        end
    end
    return
end

function main()
    haskey(INDICATORS, INDICATOR) || error("unknown indicator $INDICATOR; \
        one of $(join(sort(collect(keys(INDICATORS))), ", "))")
    spec = INDICATORS[INDICATOR]

    elev = presmooth(fill_gaps!(elevation_cube()))
    cells = DimensionalData.lookup(elev, DGG.Cells)
    vector = parent(cells)
    grid = DGG.levelgrid(DGG.system(vector), DGG.level(vector))
    @info "cube ready" cells = length(elev) level = DGG.level(vector)

    # D8 rather than DInf: it is the one method every system implements, so the
    # script keeps working if `IGeo7System` above is swapped for another.
    accumulation, _ = GM.flowaccumulation(elev; method = GM.D8())
    # Upstream area is metres squared; upstream cells is what a channel
    # threshold is quoted in and what the routing below compares.
    upstream_cells = Float64.(vec(parent(accumulation))) ./
        GM.cellarea(elev, first(eachindex(elev)))

    values = Float64.(spec.compute(elev, accumulation))
    scaled = spec.scale.(values)
    range_ = colorrange(scaled, spec.symmetric)
    palette = isempty(COLORMAP) ? spec.colormap : Symbol(COLORMAP)
    @info "coloured by $INDICATOR" colormap = palette colorrange = range_

    stem = thalweg(elev, accumulation)
    centre = ((EXTENT.X[1] + EXTENT.X[2]) / 2, (EXTENT.Y[1] + EXTENT.Y[2]) / 2)
    # Kilometres on the ground, elevation in metres: the same orthographic
    # construction the hydrology tutorial uses, in units a close camera can
    # keep its depth precision in.
    transform = GeoMakie.create_transform(
        "+proj=ortho +lon_0=$(centre[1]) +lat_0=$(centre[2]) +datum=WGS84 +units=km",
        "+proj=longlat +datum=WGS84")

    place = projector(elev, grid, transform)
    full = track(stem, place)
    lo = max(1, round(Int, TRIM_HEAD * length(full)))
    hi = min(length(full), length(full) - round(Int, TRIM_TAIL * length(full)))
    tr = full[lo:hi]
    flown = sum(norm(tr[i + 1] - tr[i]) for i in 1:length(tr) - 1)
    @info "flight" waypoints = length(tr) km = round(flown, digits = 1)

    branches = Vector{Int}[]
    if STREAMLINES
        down = downstream(upstream_cells, elev)
        branches = channels(upstream_cells, down, upstream(down)...)
        @info "channel network" branches = length(branches) cells = sum(
            length, branches; init = 0)
    end
    # Referenced to the flight rather than given in absolute terms, so that it
    # survives a change of `--seconds` or of extent: a third of the speed the
    # camera covers the same ground at. That is slower than the camera, so the
    # water does drift backwards across the screen — which was worth avoiding
    # right up until the alternative was tried, because at speeds fast enough
    # to outrun the camera nothing coming off the mountains can be followed.
    river = water(branches, upstream_cells, place,
        FLOW_SPEED == "auto" ? 0.35 * flown / max(SECONDS - HOLD_START - HOLD_END, eps()) :
            parse(Float64, FLOW_SPEED))

    figure = Figure(size = (WIDTH, HEIGHT),
        backgroundcolor = SKY ? SKY_HORIZON : BACKGROUND, figure_padding = 0)
    SKY && sky!(figure)
    scene = LScene(figure[1, 1]; show_axis = false,
        scenekw = (backgroundcolor = SKY ? :transparent : BACKGROUND, clear = !SKY))

    # With haze on, the colours are mapped here rather than on the GPU, because
    # what goes to the GPU changes every frame: the mix towards the sky depends
    # on where the camera is.
    air = Makie.RGBAf(Makie.to_color(SKY ? SKY_HORIZON : BACKGROUND))
    ground = HAZE > 0 ? shade(scaled, palette, range_) : nothing
    tinted = ground === nothing ? nothing : Observable(copy(ground))
    surface = dggsurface!(scene, cells,
        Float64.(vec(parent(elev))) .* (EXAGGERATION / 1000);
        color = tinted === nothing ? scaled : tinted,
        colormap = palette, colorrange = range_,
        # `shading` because here the surface is a shape, not a value; `wrap`
        # off because a regional patch never reaches the cut meridian, and
        # without it the mesh is exactly one vertex per cell.
        shading = true, wrap = false,
        transformation = Makie.Transformation(transform))

    # The mesh's own vertex buffer, already projected: the distances the haze
    # needs are to these, and computing them again would be the same work twice.
    fog! = (eye -> nothing)
    if tinted !== nothing
        vertices = surface.mesh_positions[]
        if length(vertices) == length(ground)
            fog! = function (eye)
                mist!(tinted[], ground, vertices, eye, air)
                notify(tinted)
                return
            end
        else
            @warn "no haze: the surface has one vertex per something other than \
                a cell" vertices = length(vertices) cells = length(ground)
        end
    end
    # Without this, `save` and `record` refit the camera on the scene's bounding
    # box and discard every view the path sets.
    cameracontrols(scene.scene).settings.center[] = false

    run! = flow!(scene, river)

    path = flightpath(tr)
    if STILLS
        for (name, t) in (("start", 0.0), ("middle", SECONDS / 2), ("end", SECONDS))
            view = path(warp(FREEZE >= 0 ? FREEZE : t))
            set_view!(scene.scene, view)
            fog!(view.eyeposition)
            run!(t)
            save("river_flythrough_$(INDICATOR)_$name.png", figure)
        end
        @info "wrote stills"
        return
    end

    out = isempty(OUTPUT) ? "river_flythrough_$(INDICATOR).mp4" : OUTPUT
    frames = round(Int, SECONDS * FPS)
    record(figure, out, range(0, SECONDS, length = frames);
        framerate = FPS, px_per_unit = 1, compression = 12,
        profile = "high", pixel_format = "yuv420p") do t
        view = path(warp(FREEZE >= 0 ? FREEZE : t))
        set_view!(scene.scene, view)
        fog!(view.eyeposition)
        run!(t)
    end
    @info "wrote $out" seconds = SECONDS frames MB = round(filesize(out) / 2^20, digits = 1)
    return
end

GLMakie.activate!(visible = false, fxaa = true, ssao = true, render_on_demand = true)
# GLFW can hand back a null primary monitor on an unattended macOS session and
# take the process down before the first frame. The video's geometry is fixed by
# `size` and `px_per_unit`, so the window's DPI is not needed.
@static if Sys.isapple()
    function Makie.window_area(scene::Scene, screen::GLMakie.Screen)
        Makie.disconnect!(screen, Makie.window_area)
        scene.events.window_dpi[] = 96.0
        return
    end
end

(abspath(PROGRAM_FILE) == @__FILE__) && main()
