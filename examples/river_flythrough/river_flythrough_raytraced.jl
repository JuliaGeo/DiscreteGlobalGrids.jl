# Fly a camera down a river valley over an IGEO7 DEM, path-traced.
#
# This is `river_flythrough.jl` under a physically based renderer. It `include`s
# that script rather than restating it, so the two fly the same route: the same
# Copernicus tile regridded onto the same IGEO7 coverage, the same thalweg
# traced out of the same flow accumulation field, the same channel network with
# the same Manning timings running down it. Everything above the drawing is
# shared, and `scene_inputs` over there is the seam.
#
# Three devices in the raster version stand in for light, and none survives
# here. Its gradient-quad background and its per-vertex haze both become one
# environment map, which is at once the background and what lights the ridge.
# Its dark rim around each water dot goes: the channel surface is a dielectric
# at index 1.33 and is visible the way a river is, by refracting what is behind
# it and reflecting the sky.
#
# The sky's colour is the one deliberate departure from physics. `skymap` runs
# Hosek-Wilkie for real and substitutes only its chroma, so the palettes survive
# on a sky whose angular radiance distribution is still a real one.
#
# Run in `docs/raytracing`, which is the docs environment resolved against
# Makie's `sd/lava` branch and the Vulkan stack under it. It carries its own
# manifest and is deliberately not a workspace member — see the comment at the
# top of `docs/raytracing/Project.toml` for why the two stacks cannot share one:
#
#     julia --project=docs/raytracing -e 'using Pkg; Pkg.instantiate()'
#     julia -t auto --project=docs/raytracing examples/river_flythrough/river_flythrough_raytraced.jl
#
# Every flag of the raster script still applies — `--extent`, `--level`,
# `--indicator`, `--palette`, `--cache`, `--exaggeration`, and the whole of the
# river's parameterisation. What this file adds is the renderer and the light:
#
#     --spp N            samples per pixel; the whole quality/time dial
#     --max-depth N      path depth
#     --denoise B        à-trous denoise; on, and worth leaving on
#     --denoise-iterations N   filter passes, each doubling the radius.
#                        1 is what the edge-stopping guides were tuned
#                        against; 4 carries the filter across a ridge
#     --sun-altitude D   degrees above the horizon
#     --sun-azimuth D    degrees clockwise from north
#     --turbidity X      2 is clear air, 10 is hazy
#     --sun-intensity X / --sky-intensity X   the two lights, against each other
#     --sky-contrast X   how much of the palette's own darkness the sky takes on
#     --albedo X         the brightest the terrain's colormap is allowed to be
#     --terrain-coat B   dielectric coat over the terrain — see `terrain_material`
#     --terrain-roughness X   how broad the coat's highlight is; 0 is a mirror
#     --range-high Q     quantile the colormap's top end clips at; turn it
#                        down to spend the map on the hillslopes rather than
#                        on the two per cent of cells that are channels
#     --water-ior X      1.33 is water; 1.31 is ice
#     --water-width M    the trunk's width across the ground, in metres
#     --bead-radius M    the droplets running down it
#     --bead-material    pearl or water — see `bead_material`
#     --bead-color C     the beads' colour; the palette's accent by default
#     --at T             seconds into the flight, for a single still
#     --stills           three stills instead, to frame a shot
#     --video            render a clip of the flight
#     --clip-start S / --clip-seconds S   which part of it
#
# Measured on one RTX A5000 with hardware ray tracing, over the full level-13
# cube of 1.6M cells, warm — a run's first frame also pays scene build and
# compilation:
#
#     640×360    64 spp    0.6 s/frame
#     1920×1080  128 spp   5.2 s/frame    (24 s at 24 fps = 576 frames, ~55 min)
#
# Cost is sub-linear in pixels: the second line is 18× the samples of the first
# for 9× the time, because 1080p keeps the GPU better occupied.
#
# And on one RX 9060 XT (RDNA4, RADV) at 1920×1080 over the same cube, warm,
# which is where the defaults below come from:
#
#     128 spp  --denoise false   7.2 s/frame   (24 s at 30 fps = 720 frames, ~87 min)
#      32 spp  --denoise true    1.9 s/frame   (the same 720 frames, ~23 min)
#
# The cheaper line is also the better picture. Its grain is gone where the
# 128-spp frame still speckles the shaded slopes, and a quarter of the samples
# is what pays for the four-fold speedup. Measured against a 256-spp reference
# the raw 128-spp frame is nearer in RMS — 0.015 against 0.023 — so the filter
# is trading a little bias for the noise it removes, and on this terrain that
# trade reads as an improvement.
#
# Nearly all of a frame is the trace itself: `rt_indirect` is 79% of GPU time,
# and the GPU is busy 89% of the wall clock, so the sample count is the dial
# that matters and host-side work is not worth tuning.
#
# `--video` goes through `record_longrunning`, which writes every frame to
# `<output>_frames/` and skips the ones already there when re-run — so an
# interrupted clip costs the frame it was on rather than all of them, and a
# clip can be re-encoded without re-rendering.

# The raster script reads settings from the command line first and the
# environment second, so a default this backend wants differently is set here,
# before the include. A `--flag` still wins.
#
# A dot there is a marker; here it is a sphere with its own material and its own
# place in the acceleration structure.
get!(ENV, "MAX_DOTS", "4000")

# Flow accumulation is skewed enough that at the 98th percentile the top two per
# cent of cells — the channels — take most of the colormap, leaving every
# hillslope crowded into the bottom of it. Clipping at the 80th spreads the
# hillslopes across the map; the channels saturate at the top.
get!(ENV, "RANGE_HIGH", "0.80")

include(joinpath(@__DIR__, "river_flythrough.jl"))

using RayMakie, Hikari
import Lava
using Colors: RGB
using FileIO: save

# ---------------------------------------------------------------------------
# The renderer
# ---------------------------------------------------------------------------

const HW_ACCEL = setting("hw_accel", true)

"""
    device() -> Lava.LavaBackend

The Vulkan backend this renders on, having checked it can trace rays in
hardware.

Named directly rather than obtained through `GPUSelect`, which RayDemo's scenes
use: `GPUSelect.Backend(:Lava)` looks for `Manifest.toml` beside the active
project, which is not where a Julia workspace keeps it, and returns
`KernelAbstractions.CPU()` when it finds none.

`rt_pipeline_properties` is non-`nothing` exactly when the device carries
`VK_KHR_ray_tracing_pipeline` and `VK_KHR_acceleration_structure`, the pair
`hw_accel` needs. Without them Lava traverses the BVH in a compute shader —
correct, and much slower on geometry this dense — so a missing extension is an
error rather than a silent downgrade. `--hw-accel false` takes that path on
purpose.
"""
function device()
    backend = Lava.LavaBackend()
    context = Lava.vk_context()
    rt = context.rt_pipeline_properties
    if rt === nothing
        HW_ACCEL && error("this GPU has no VK_KHR_ray_tracing_pipeline, so \
            `--hw-accel true` cannot be honoured: $(context.device_name). Pass \
            `--hw-accel false` to trace on Lava's software BVH instead.")
        @warn "tracing on Lava's software BVH" gpu = context.device_name
    else
        @info "Lava" gpu = context.device_name hardware_ray_tracing = HW_ACCEL max_ray_recursion =
            Int(rt.max_ray_recursion_depth)
    end
    return backend
end

const DEVICE = device()

const SPP = setting("spp", 24)
const MAX_DEPTH = setting("max_depth", 8)
# On by default: at one iteration the filter removes more noise than the bias
# it adds, so it buys back most of the samples it lets you drop.
const DENOISE = setting("denoise", true)
# Each pass doubles the filter radius, so 4 — the `DenoiseConfig` default —
# reaches about sixteen pixels and takes the ridges with it. The aux normal and
# depth guides that stop the filter at an edge were tuned at one pass.
const DENOISE_ITERATIONS = setting("denoise_iterations", 1)
# 0.86 rather than 0.5 because the terrain is coated: a clearcoat over a bright
# diffuse base returns about 58% of what the bare base returns, so the default
# exposure carries the 1/0.58 that puts the valley back where it was. With
# `--terrain-coat false`, 0.5 is the matching exposure.
const EXPOSURE = setting("exposure", 0.86)
const TONEMAP = setting("tonemap", :aces)
const GAMMA = setting("gamma", 2.2)
const ISO = setting("iso", 100.0)
const WHITEBALANCE = setting("whitebalance", 6500.0)

# What to render. A still by default: at 1080p and a real sample count a frame
# is minutes, and the flight is `--seconds` × `--fps` of them.
const AT = setting("at", -1.0)                   # seconds in; the middle if unset
const VIDEO = setting("video", false)
const CLIP_START = setting("clip_start", -1.0)   # seconds in; centred if unset
const CLIP_SECONDS = setting("clip_seconds", 4.0)

# ---------------------------------------------------------------------------
# The light
# ---------------------------------------------------------------------------

# Altitude and azimuth rather than a vector: a valley's appearance is mostly a
# question of the sun's angle to the way it runs. Azimuth is clockwise from
# north, so 135° rakes an east-west valley from the side and puts a shadow in
# every side gully.
const SUN_ALTITUDE = setting("sun_altitude", 34.0)
const SUN_AZIMUTH = setting("sun_azimuth", 138.0)
# Hikari's clear-air sun-to-sky ratio is 5:1. This scene narrows it: the
# terrain's colour is the indicator, and a hard sun against a dim sky drops half
# the valley into shadow too deep to read the field in.
const SUN_INTENSITY = setting("sun_intensity", 2.8)
const SKY_INTENSITY = setting("sky_intensity", 2.2)
const TURBIDITY = setting("turbidity", 2.6)
const SKY_RESOLUTION = setting("sky_resolution", 512)

# How much of the palette's own darkness the sky takes on top of the model's.
# The palettes are near-black at the zenith, which obeyed literally leaves the
# valley unlit. 0 keeps the model's brightness and only the palette's hue; 1 is
# as dark as the palette says.
const SKY_CONTRAST = setting("sky_contrast", 0.45)

# A near-white sun; the palette lives in the sky and its accent goes to the
# beads. A saturated sun multiplies every lit surface by its own hue, which
# repaints the indicator the terrain is carrying. A sky tints by adding instead,
# so the look survives without the data being recoloured.
const SUN_COLOR = setting("sun_color", "#fff3e0")

# What fills the environment map below the horizon. The terrain is one tile, so
# between ridges the camera sees past its edge. The horizon's colour, dimmed,
# reads as distance; anything darker reads as a hole in the frame.
const GROUND_DIM = setting("ground_dim", 0.55)

"""
    sundirection() -> Vec3f

The unit vector pointing at the sun, in the scene's own frame.

The scene is x east, y north, z up — the orthographic projection puts ground
kilometres in x and y and elevation in z — so an altitude and an azimuth
clockwise from north land here without a change of basis.
"""
function sundirection()
    alt, azi = deg2rad(SUN_ALTITUDE), deg2rad(SUN_AZIMUTH)
    return Makie.Vec3f(cos(alt) * sin(azi), cos(alt) * cos(azi), sin(alt))
end

"""
    skymap(sun) -> Matrix{RGBf}

The Hosek-Wilkie sky over this valley, wearing the palette's colours.

The background of a path-traced scene is also its light, so this is at once
what is behind the ridge and what lights it.

The model is run for real and only its chroma replaced: each pixel keeps the
radiance Hosek-Wilkie gives it — the aureole around the sun, the bright horizon
band, the falloff to the zenith — and takes its hue from the palette's
horizon-to-zenith ramp at its own elevation. The angular distribution of
radiance survives, which is what reads as sky and what the environment light's
importance sampling needs to find the bright parts.

`--sky-contrast` is the one dial here that is art rather than physics; see
`SKY_CONTRAST`.

The map is in Hikari's equal-area octahedral square, the parameterisation
`EnvironmentMap` reads. `sunsky_to_envlight` returns one already in it, so
directions are recovered through the same mapping that wrote it.

The sun disc is not in the map — `sunsky_to_envlight` returns it separately as a
delta light, pbrt's arrangement, which keeps a hundred-thousand-to-one pixel
from throwing fireflies.
"""
function skymap(sun)
    env, _ = Hikari.sunsky_to_envlight(;
        direction = sun,
        intensity = 1f0,
        turbidity = Float32(TURBIDITY),
        # No ground hemisphere: with it, everything below the horizon is one
        # flat albedo and the seam between that and the sky is a visible
        # scalloped edge where the equal-area square folds. Without it the model
        # clamps the zenith angle at the horizon, so a downward direction gets
        # the radiance of the horizon it is under — which is both smooth and the
        # right order of magnitude for haze.
        ground_enabled = false,
        resolution = SKY_RESOLUTION,
    )
    sky = env.env_map.data
    rows, cols = size(sky)

    horizon = Makie.RGBf(Makie.to_color(SKY_HORIZON))
    zenith = Makie.RGBf(Makie.to_color(SKY_TOP))

    # A palette colour carries a hue and a brightness, and they are wanted for
    # different things: the hue outright, the brightness only as far as
    # `SKY_CONTRAST` allows. Splitting it at its own maximum channel is what
    # separates them — `chroma` is the colour at full brightness, `value` is how
    # bright the palette wanted it.
    split(c) = (v = max(c.r, c.g, c.b, 1f-4); (Makie.RGBf(c.r / v, c.g / v, c.b / v), v))
    horizon_chroma, horizon_value = split(horizon)

    out = Matrix{Makie.RGBf}(undef, rows, cols)
    for v in 1:rows, u in 1:cols
        uv = Makie.Point2f((u - 0.5f0) / cols, (v - 0.5f0) / rows)
        direction = Hikari.equal_area_square_to_sphere(uv)
        Y = Hikari.to_Y(sky[v, u])
        if direction[3] <= 0
            chroma, value = horizon_chroma, horizon_value * Float32(GROUND_DIM)
        else
            # `smoothstep` is the raster script's, and the ramp is over the
            # sine of the elevation rather than the angle: that is the
            # coordinate the sky's own gradient runs in.
            t = smoothstep(direction[3])
            mixed = Makie.RGBf(horizon.r + (zenith.r - horizon.r) * t,
                horizon.g + (zenith.g - horizon.g) * t,
                horizon.b + (zenith.b - horizon.b) * t)
            chroma, value = split(mixed)
        end
        scale = Float32(Y * value^SKY_CONTRAST)
        out[v, u] = Makie.RGBf(chroma.r * scale, chroma.g * scale, chroma.b * scale)
    end
    return out
end

"""
    lights(sun) -> Vector

The sky and the sun, as the two lights this scene has.

Nothing else is lit by anything: there is no fill, no ambient term and no light
following the camera, because a valley at ten in the morning has a sun and a sky
and that is the whole of it. What stands in for a fill is the ground albedo in
the sky map — the far wall of the valley is lit by the near one, which is what
bounced light is, and a path tracer does that on its own.

The sun is near-white and the palette lives in the sky, which is the opposite of
where it first went. See `SUN_COLOR`: a saturated sun multiplies every lit
surface by its own hue and the terrain stops showing the indicator and starts
showing the sun. A sky tints by adding rather than by multiplying, so the look
survives there without the data being repainted on its way to the eye.

The sun's magnitude is a ratio against the sky's, not a number in a vacuum —
Hikari's own figure for clear air is 5:1 and `SUN_INTENSITY` explains why this
scene does not use it.
"""
function lights(sun)
    tint = Makie.RGBf(Makie.to_color(SUN_COLOR))
    peak = max(tint.r, tint.g, tint.b, 1f-4)
    scale = Float32(SUN_INTENSITY) / peak
    return [
        Makie.EnvironmentLight(Float32(SKY_INTENSITY), skymap(sun)),
        # Makie's `DirectionalLight` takes the direction light *travels*, which
        # is away from the sun.
        Makie.DirectionalLight(
            Makie.RGBf(tint.r * scale, tint.g * scale, tint.b * scale),
            Makie.Vec3f(-sun)),
    ]
end

# ---------------------------------------------------------------------------
# The materials
# ---------------------------------------------------------------------------

const ALBEDO = setting("albedo", 0.70)      # the brightest the terrain gets

# The terrain carries both a colour and a material, and the two have to combine:
# the colormap is the indicator and the material is what makes it look like a
# surface rather than a chart. RayMakie merges them by putting the vertex colours
# into the material's diffuse slot — `CoatedDiffuse.reflectance` here.
#
# That merge needs MakieOrg/Makie.jl#5756. `extract_material` decided whether to
# merge by testing whether `:material` is in `plot.attributes.inputs`, reading a
# miss as a deliberate post-construction override; but a recipe forwards
# `Computed` nodes, which are graph edges and never land in `inputs`, so every
# recipe-drawn mesh given a material lost its colormap. (The override it was
# protecting is not reachable that way in any case: `plot.material = mat` on such
# a plot throws `Cannot attach input with name material - already exists!`.)
# Without that fix this surface renders in flat coat-grey.

const TERRAIN_COAT = setting("terrain_coat", true)
const TERRAIN_ROUGHNESS = setting("terrain_roughness", 0.0)

"""
    terrain_material() -> Union{Hikari.Material, Nothing}

What the ground is made of: a dielectric coat over the colormap.

`Diffuse` alone — which is what a plot with a colour and no material gets — is
Lambertian, and a Lambertian surface has no highlight at all. It is the right
model for chalk and the wrong one for a mountainside, which is mineral and damp
and throws a broad sheen back at the sky. `CoatedDiffuse` layers a smooth
dielectric over the diffuse base, so the colormap still reads as the albedo
underneath while the coat adds the specular lobe on top.

`TERRAIN_ROUGHNESS` sets how wide that lobe is: 0 is a mirror-flat varnish that
returns a sharp image of the sky, and larger values spread the highlight out over
the slope. 0 is the default for two measured reasons. It holds the colormap: a
rough coat scatters a broad *white* highlight over the whole slope, and measured
over this frame that costs a fifth of the image's mean saturation by roughness
0.05 and a third by 0.3 — the greying this script refuses to do to the colormap
elsewhere. And it is not the noisier choice it sounds like: per-pixel variance on
the shadowed slope is the same at roughness 0 as with no coat at all (0.0282 vs
0.0281), because what is grainy there is the dim slope, not the coat.

`--terrain-coat false` drops back to the plain Lambertian surface.

The coat darkens what it covers, and by a lot: about 42% of the terrain's
luminance, converged — 8x the layered estimator's walk budget moves it by 0.2 of
a percentage point, so this is the answer and not a truncated one. It is also the
right answer. At eta 1.5 the critical angle is 41.8 degrees, so a bit over half of
the light leaving the diffuse base is totally internally reflected back down into
it and partly re-absorbed; a clearcoat over a bright base really does darken it.
`EXPOSURE` carries the compensation.
"""
function terrain_material()
    TERRAIN_COAT || return nothing
    return Hikari.CoatedDiffuse(roughness = Float32(TERRAIN_ROUGHNESS))
end

"""
    groundcolors(scaled, colormap, colorrange) -> Vector{RGBAf}

The indicator through the colormap, rescaled to albedos a surface could have.

`shade` — the raster script's, reused — maps the values to colours. What it
cannot know is that those colours are about to be *reflectances*: a colormap
runs to white at one end, and a white mountainside in direct sun blows out to
the top of the sensor and takes its neighbours with it. Snow is 0.8, dry rock is
0.2, fresh basalt is 0.04. `ALBEDO` is the brightest end of that.

The scaling must be multiplicative — one factor on all three channels. Mapping
each channel into `[floor, ALBEDO]` separately adds a constant to red, green and
blue alike, and a constant added equally to three channels is grey: it
desaturates the whole map towards neutral in proportion to darkness. Dimming the
colormap is allowed, greying it is not.
"""
function groundcolors(scaled, colormap, colorrange)
    scale = Float32(ALBEDO)
    return map(shade(scaled, colormap, colorrange)) do c
        Makie.RGBAf(c.r * scale, c.g * scale, c.b * scale, 1f0)
    end
end

const WATER_IOR = setting("water_ior", 1.33)
const WATER_COLOR = setting("water_color", PALETTE.dot)

"""
    water_material() -> Hikari.Material

What the river is made of.

A smooth dielectric at the index of water. It neither glows nor overlays the
terrain: it is visible because it refracts the valley behind it and reflects the
sky above it.

`Kt` is a transmission tint, not a surface colour, so it acts on what comes
through and grows with path length through the water.
"""
water_material() = Hikari.Dielectric(index = Float32(WATER_IOR))

# What a bead is made of. The beads carry the water's speed, and at this range
# each is a few pixels across — too small to read as a lens, which is why the
# default is opaque rather than the channel's dielectric.
#
#   pearl  an opaque diffuse bead. Most legible against any colormap.
#   water  the channel's dielectric. Physically of a piece with the river and
#          more subtle; `Kt` takes the bead colour, so a warm `--bead-color`
#          gives amber glass rather than water.
#
# Emissive beads — the raster version's glowing dots — are not available:
# Hikari has no per-instance emissive material, and one `mesh!` per bead would
# be thousands of plots rebuilt every frame.
#
# Their colour is the palette's `accent`, the hue it reserves as one the terrain
# colormap does not own. A bead in the water's own pale colour disappears into
# the channel, which is what the raster version's dark rim was there to prevent.
#
# Both per-instance material paths need JuliaGraphics/Hikari.jl#40 and
# SimonDanisch/Lava.jl#9: without the first a dielectric bead renders black, and
# without the second every bead loses its material after the first frame of an
# animation.
const BEAD_MATERIAL = Symbol(setting("bead_material", "pearl"))
const BEAD_COLOR = setting("bead_color", PALETTE.accent)

function bead_material()
    BEAD_MATERIAL === :pearl && return Hikari.Diffuse(σ = 0f0)
    BEAD_MATERIAL === :water && return water_material()
    error("unknown --bead-material $BEAD_MATERIAL; one of pearl, water \
        (`glow` is not available: Hikari has no per-instance emissive)")
end

# ---------------------------------------------------------------------------
# The river
# ---------------------------------------------------------------------------

# The water is drawn wider than the Adige actually is, and deliberately so. At
# the height this camera flies, a channel at its true thirty to fifty metres is
# one or two pixels of a 1280-wide frame — which is honest and invisible, and an
# invisible river is not what this is for. The trunk is drawn at roughly twice
# life size and the network narrows from there by the same flow law, so what the
# width still tells you truthfully is the *ordering*: which reach carries more.
const WATER_WIDTH = setting("water_width", 90.0)   # m across, at the trunk
const WATER_LIFT = setting("water_lift", 6.0)      # m above the bed
const BEADS = setting("beads", true)
const BEAD_RADIUS = setting("bead_radius", 34.0)   # m at the trunk

"""
    channels!(scene, river) -> Plot

The channel network, as the surface of the water running down it.

A flat strip, not a tube, and the shape is what makes a river visible from the
air. `Makie.tubes!` would build this for free, but a cylinder's visible normals
fan through every direction, so most of what a camera sees of one is at
near-normal incidence — where a dielectric reflects two per cent and transmits
the rest into an unlit riverbed. A level surface viewed from above is seen at
sixty to eighty degrees off normal instead, and Fresnel there returns a tenth to
a half of the sky.

The strip is offset in the horizontal plane, so it stays level across its width
however steeply the reach falls. Where the valley is a narrow V it cuts into
both walls, which is what a river in a gorge does.

Width goes as the flow, so the network narrows upstream — the quantity the
raster version put into `linewidth`, in ground metres rather than screen
pixels.
"""
function channels!(scene, river::Water)
    isempty(river.points) && return nothing
    GB = Makie.GeometryBasics
    # `river.points` sit `DOT_LIFT` above the bed, which is where a bead of
    # spray belongs and is too high for a water surface; the difference is taken
    # back off here, through the same exaggeration the heights carry.
    drop = Float32((DOT_LIFT - WATER_LIFT) * EXAGGERATION / 1000)

    verts = Makie.Point3f[]
    faces = GB.GLTriangleFace[]
    for c in 1:length(river.starts) - 1
        lo, hi = river.starts[c], river.starts[c + 1] - 1
        hi - lo < 1 && continue
        base = length(verts)
        for j in lo:hi
            p = river.points[j]
            # The centred difference along the branch, flattened to the ground
            # plane: the strip is laid across the direction the water runs, and
            # how fast it is falling is not part of which way that is.
            ahead = river.points[min(j + 1, hi)]
            behind = river.points[max(j - 1, lo)]
            tx, ty = ahead[1] - behind[1], ahead[2] - behind[2]
            len = sqrt(tx * tx + ty * ty)
            sx, sy = len < 1e-12 ? (1.0, 0.0) : (-ty / len, tx / len)
            half = WATER_WIDTH / 2000 * (0.22 + 0.78 * river.weights[j]^0.35)
            z = Float32(p[3]) - drop
            push!(verts, Makie.Point3f(p[1] + sx * half, p[2] + sy * half, z))
            push!(verts, Makie.Point3f(p[1] - sx * half, p[2] - sy * half, z))
        end
        for k in 0:(hi - lo - 1)
            a, b = base + 2k + 1, base + 2k + 2
            push!(faces, GB.GLTriangleFace(a, b, a + 2))
            push!(faces, GB.GLTriangleFace(b, b + 2, a + 2))
        end
    end
    isempty(faces) && return nothing
    @info "water" surface_vertices = length(verts) triangles = length(faces) branches =
        length(river.starts) - 1
    return Makie.mesh!(scene, GB.normal_mesh(GB.Mesh(verts, faces));
        color = WATER_COLOR, material = water_material())
end

"""
    beads!(w, t, positions, sizes)

Put every droplet where it is `t` seconds into the flight, and size it.

The same lookup `advance!` does in the raster script, with the fade moved from
the colour into the radius: a dielectric has no opacity to fade, so a droplet
reaching the end of its branch shrinks away instead.

Each bead is lifted by its own radius so it rests on the water rather than in
it. `DOT_LIFT` is eight metres and a bead is tens across, so a bead at the
channel's own height sits mostly below ground, where rays entering it meet
opaque terrain immediately.
"""
function beads!(w::Water, t, positions, sizes)
    surface = Float32((DOT_LIFT - WATER_LIFT) * EXAGGERATION / 1000)
    @inbounds for k in eachindex(w.branch)
        c = w.branch[k]
        lo, hi = w.starts[c], w.starts[c + 1] - 1
        span = w.times[hi]
        along = mod(w.offset[k] + t, span)
        i = min(_seek(w.times, lo, hi, along), hi - 1)
        u = (along - w.times[i]) / max(w.times[i + 1] - w.times[i], eps())
        p = w.points[i] + u * (w.points[i + 1] - w.points[i])
        carried = w.weights[i] + u * (w.weights[i + 1] - w.weights[i])
        edge = min(DOT_FADE, span / 3)
        fade = clamp(min(along, span - along) / edge, 0.0, 1.0)
        r = Float32(BEAD_RADIUS / 1000 * (0.45 + 0.55 * carried^0.25) *
            w.tail[k] * fade)
        sizes[k] = r
        positions[k] = Makie.Point3f(p[1], p[2], Float32(p[3]) - surface + r)
    end
    return
end

"""
    droplets!(scene, river) -> (t -> nothing)

Put the beads of water in `scene`, and return the function that runs them.

They are the raster version's dots, made of the same thing the channels are.
What they are for is unchanged and is the one thing a still cannot show: the
water's *speed*, which is Manning's off the reach slope, so the torrent off a
headwall outruns the Adige on its floodplain by the ratio of the square roots of
their gradients. Spacing them by travel time rather than by distance is what
draws that — they crowd where the water is slow and string out where it runs.

Only the positions and the radii change per frame. The geometry is one sphere
and the material is one dielectric; every droplet is an instance of them, so a
frame costs a transform update and not a rebuild.
"""
function droplets!(scene, river::Water)
    (!BEADS || isempty(river.branch)) && return t -> nothing
    n = length(river.branch)
    positions = Makie.Observable(Vector{Makie.Point3f}(undef, n))
    sizes = Makie.Observable(Vector{Float32}(undef, n))
    beads!(river, 0.0, positions[], sizes[])
    Makie.meshscatter!(scene, positions; markersize = sizes,
        color = BEAD_COLOR, material = bead_material())
    @info "water" droplets = n
    return function (t)
        beads!(river, t, positions[], sizes[])
        notify(positions); notify(sizes)
        return
    end
end

# ---------------------------------------------------------------------------
# The scene
# ---------------------------------------------------------------------------

"""
    build() -> (scene, path, run!)

The valley, its river and the camera's route through it.

A bare `Scene` rather than a `Figure` with an `LScene` in it. The raster version
needs the figure because its sky is a plot in the figure's pixel space; there is
no such plot here, and a figure would only bring a layout and a second camera
for a picture that is one 3D scene edge to edge.
"""
function build()
    (; elev, cells, transform, scaled, palette, range_, tr, river) = scene_inputs()

    sun = sundirection()
    scene = Makie.Scene(size = (WIDTH, HEIGHT), lights = lights(sun))
    Makie.cam3d!(scene)

    dggsurface!(scene, cells,
        Float64.(vec(parent(elev))) .* (EXAGGERATION / 1000);
        color = groundcolors(scaled, palette, range_),
        material = terrain_material(),
        # `wrap` off because a regional patch never reaches the cut meridian,
        # and without it the mesh is exactly one vertex per cell. `shading` is
        # nominally on, but the ray tracer takes the material and not the
        # attribute: what shades this surface is the light in the scene.
        shading = true, wrap = false,
        transformation = Makie.Transformation(transform))

    channels!(scene, river)
    run! = droplets!(scene, river)

    # Without this, the first `colorbuffer` refits the camera on the scene's
    # bounding box and discards every view the path sets.
    Makie.cameracontrols(scene).settings.center[] = false

    return scene, flightpath(tr), run!
end

"""
    integrator() -> Hikari.VolPath

The path tracer, on hardware ray tracing.

`hw_accel` asks Lava for a `VK_KHR_ray_tracing_pipeline` acceleration structure
instead of a BVH walked in a compute shader. On this scene it is the difference
worth having: a couple of million triangles of terrain plus the water surface is exactly
the geometry-heavy case hardware traversal is for.
"""
function integrator()
    vp = Hikari.VolPath(samples = SPP, max_depth = MAX_DEPTH, hw_accel = HW_ACCEL)
    # ISO and white balance are the sensor's, not the tone mapper's: they act on
    # the spectral radiance before it is ever RGB, which is where a camera acts.
    vp.sensor = Hikari.PixelSensor(iso = Float32(ISO), whitebalance = Float32(WHITEBALANCE))
    return vp
end

"""
    renderkw() -> NamedTuple

The screen's whole configuration, built once per run.

Once, because Makie keys its cached screen on this: a fresh `integrator()` per
call is a different config, which is a different screen, which is the
acceleration structure over three million triangles built again. Three stills
would pay for it three times.
"""
renderkw() = (device = DEVICE, integrator = integrator(), exposure = Float32(EXPOSURE),
    tonemap = TONEMAP === :none ? nothing : TONEMAP, gamma = Float32(GAMMA),
    denoise = DENOISE,
    denoise_config = DENOISE ? Hikari.DenoiseConfig(iterations = DENOISE_ITERATIONS) : nothing,
    update = false)

"""
    frame!(scene, path, run!, t)

Put the camera and the water where they are `t` seconds into the flight.

The camera flies on the eased clock and the water runs on the real one, which is
the raster version's arrangement and is right: the ease is a property of the
shot, and the river does not know about the shot.
"""
function frame!(scene, path, run!, t)
    set_view!(scene, path(warp(FREEZE >= 0 ? FREEZE : t)))
    run!(t)
    return
end

function render()
    RayMakie.activate!()

    scene, path, run! = build()
    kw = renderkw()
    @info "renderer" hw_accel = HW_ACCEL spp = SPP max_depth = MAX_DEPTH size = (WIDTH, HEIGHT)

    if VIDEO
        start = CLIP_START >= 0 ? CLIP_START : max(0.0, (SECONDS - CLIP_SECONDS) / 2)
        stop = min(SECONDS, start + CLIP_SECONDS)
        frames = max(2, round(Int, (stop - start) * FPS))
        out = isempty(OUTPUT) ? "river_flythrough_raytraced_$(INDICATOR).mp4" : OUTPUT
        @info "clip" from = start to = stop frames
        # `record_longrunning` writes every frame to disk beside the video and
        # skips the ones already there on a re-run. At minutes a frame that is
        # not a convenience, it is the difference between an interrupted render
        # costing the frame it was on and costing all of them.
        Makie.record_longrunning(scene, out, range(start, stop, length = frames);
            framerate = FPS, compression = 12, profile = "high",
            pixel_format = "yuv420p", kw...) do t
            frame!(scene, path, run!, t)
        end
        @info "wrote $out" MB = round(filesize(out) / 2^20, digits = 1)
        return out
    end

    stills = STILLS ?
        [("start", 0.0), ("middle", SECONDS / 2), ("end", SECONDS)] :
        [("", AT >= 0 ? AT : SECONDS / 2)]
    written = String[]
    for (name, t) in stills
        frame!(scene, path, run!, t)
        image = @time("render", Makie.colorbuffer(scene; kw...))
        out = if !isempty(OUTPUT) && length(stills) == 1
            OUTPUT
        else
            "river_flythrough_raytraced_$(INDICATOR)$(isempty(name) ? "" : "_" * name).png"
        end
        save(out, image)
        @info "wrote $out" at_seconds = round(t, digits = 2)
        push!(written, out)
    end
    return written
end

if abspath(PROGRAM_FILE) == @__FILE__
    render()
end
