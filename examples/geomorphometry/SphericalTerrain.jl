"""
    SphericalTerrain

Terrain analysis on a DGGS, written directly against `DiscreteGlobalGrids`'
`neighbors` / `subtree_halo` / `halo` API.

The functions here are ports of the generic (non-`AbstractMatrix`) methods
proposed in <https://github.com/Deltares/Geomorphometry.jl/pull/83>: maximum
downward gradient, slope, aspect, curvature, TPI, TRI, roughness, prominence,
D8 flow direction and flow accumulation.

Deliberately NOT ported: that PR's "relative cell" abstraction
(`RelativeZ7Cell`, `neighbor - center`, `cell + dir`). Everything here keeps
absolute cell ids and resolves them to array positions through `cellposition`.

Everything runs against one small interface, so that the *same* function body
runs over a whole level grid and over a chunk-plus-halo read, and the two can
be compared cell for cell:

    positions(f)     # grid positions we produce output for, in output order
    value(f, p)      # elevation at grid position `p`, or NaN if not held
    ctx(f)           # shared per-grid stencil cache

`value` counts every unavailable position in `f.misses[]`. That counter is the
sharp halo test: over a chunk field it must stay at zero, because the halo is
by definition every outside cell a one-ring stencil can reach.
"""
module SphericalTerrain

import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: Vertex, Edge, Connectivity

export GridCtx, WholeField, ChunkField, value, positions, ctx, reset_misses!,
    chunk_range
export max_downward_gradient, slope_mdg, roughness, tpi, tri, prominence,
    aspect, planefit_slope, curvature, flow_direction, flow_accumulation
export gcdistance, bearing_deg, METRICS, same, first_difference
export FIELD_KINDS, make_field

# ---------------------------------------------------------------------------
# geometry helpers
# ---------------------------------------------------------------------------

"Great-circle distance, in radians, between two unit-sphere points."
@inline function gcdistance(a, b)
    dx = a[1] - b[1]; dy = a[2] - b[2]; dz = a[3] - b[3]
    chord = sqrt(dx * dx + dy * dy + dz * dz)
    return 2 * asin(min(1.0, chord / 2))
end

@inline dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]

"""
    local_frame(a) -> (east, north)

A deterministic right-handed ENU tangent frame at unit-sphere point `a`.
Falls back to a fixed frame at the poles rather than producing NaNs, so polar
cells need no special case in the callers.
"""
function local_frame(a)
    ex, ey, ez = -a[2], a[1], 0.0          # ẑ × a
    n = sqrt(ex * ex + ey * ey + ez * ez)
    if n < 1e-12                            # a is (numerically) a pole
        ex, ey, ez = 1.0, 0.0, 0.0
        d = dot3((ex, ey, ez), a)
        ex -= d * a[1]; ey -= d * a[2]; ez -= d * a[3]
        n = sqrt(ex * ex + ey * ey + ez * ez)
        if n < 1e-12
            ex, ey, ez, n = 0.0, 1.0, 0.0, 1.0
        end
    end
    east = (ex / n, ey / n, ez / n)
    north = (a[2] * east[3] - a[3] * east[2],
             a[3] * east[1] - a[1] * east[3],
             a[1] * east[2] - a[2] * east[1])
    return east, north
end

"Bearing in degrees clockwise from local north, from unit point `a` to `b`."
function bearing_deg(a, b)
    east, north = local_frame(a)
    d = dot3(a, b)
    v = (b[1] - d * a[1], b[2] - d * a[2], b[3] - d * a[3])
    return mod(atand(dot3(v, east), dot3(v, north)), 360.0)
end

# ---------------------------------------------------------------------------
# the per-grid stencil cache
#
# This is the first thing a real consumer has to write, and it is the first
# API friction: `neighbors` speaks cell ids, arrays speak positions, and there
# is no positional one-ring table that spans a chunk PLUS its halo
# (`halo_table` is in-set only). So we build one.
# ---------------------------------------------------------------------------

struct GridCtx{G,K,C}
    grid::G
    connectivity::K
    cells::Vector{C}                     # position -> cell id
    centroids::Vector{NTuple{3,Float64}} # position -> unit-sphere centroid
    nbr::Vector{Int}                     # CSR neighbour positions
    ptr::Vector{Int}                     # CSR row pointers, length n+1
end

function GridCtx(sys, level::Integer; connectivity::Connectivity = Vertex())
    g = DGG.levelgrid(sys, level)
    n = DGG.ncells(g)
    cs = [DGG.cellindex(g, p) for p in 1:n]
    cent = Vector{NTuple{3,Float64}}(undef, n)
    for p in 1:n
        q = DGG.cell_centroid(g, cs[p])
        cent[p] = (q[1], q[2], q[3])
    end
    nbr = Int[]; ptr = Vector{Int}(undef, n + 1)
    for p in 1:n
        ptr[p] = length(nbr) + 1
        for m in DGG.neighbors(g, cs[p], 1; connectivity)
            push!(nbr, DGG.cellposition(g, m))
        end
    end
    ptr[n + 1] = length(nbr) + 1
    return GridCtx(g, connectivity, cs, cent, nbr, ptr)
end

@inline nrange(k::GridCtx, p::Int) = k.ptr[p]:(k.ptr[p + 1] - 1)
Base.length(k::GridCtx) = length(k.cells)

"""
    chunk_range(sys, root, level) -> UnitRange{Int}

The chunk's contiguous position block.

`descendant_range` is the documented way to get this, but it only exists where
`has_sorted_subtrees(sys)` is true — and it is NOT true for `A5System`, which
still has a perfectly good `subtree_halo`. So a consumer who wants the
chunk-plus-halo read to work on all six systems has to write this fallback:
materialise `descendants`, take the position hull, and check it really is a
hull. That is API friction worth naming, not a bug.
"""
function chunk_range(sys, root, level::Integer)
    if DGG.has_sorted_subtrees(sys)
        return DGG.descendant_range(sys, root, level)
    end
    g = DGG.levelgrid(sys, level)
    ps = sort!([DGG.cellposition(g, c) for c in DGG.descendants(sys, root, level)])
    r = first(ps):last(ps)
    length(r) == length(ps) || throw(ArgumentError(
        "subtree of $root at level $level is not a contiguous position block"))
    return r
end

const CTX_CACHE = Dict{Any,Any}()

"Memoised `GridCtx`; building one is O(ncells) and the fuzzer reuses them."
function gridctx(sys, level::Integer, connectivity::Connectivity)
    get!(CTX_CACHE, (sys, Int(level), connectivity)) do
        GridCtx(sys, level; connectivity)
    end
end

# ---------------------------------------------------------------------------
# fields
# ---------------------------------------------------------------------------

"""
    WholeField(ctx, z)

Elevation over an entire level grid: `z[p]` is the value at position `p`.
"""
struct WholeField{K}
    ctx::K
    z::Vector{Float64}
    misses::Base.RefValue{Int}
end
WholeField(k::GridCtx, z::Vector{Float64}) = WholeField(k, z, Ref(0))

"""
    ChunkField(ctx, sys, root, z_whole; halo_connectivity = ctx.connectivity)

A chunk-plus-halo read, built exactly the way the halo design doc advertises:

* `descendant_range(sys, root, level)` is the contiguous position block a
  chunked store hands back for the data;
* `subtree_halo(sys, root, level; connectivity)` is the extra fetch list, whose
  cell ids we must map to positions ourselves.

Only those two blocks of `z_whole` are copied in. Everything else is
unavailable and shows up as a miss.

`halo_connectivity` exists so that a deliberately-too-small halo can be used as
a negative control: an `Edge()` halo under a `Vertex()` stencil must produce
misses, or the miss detector proves nothing.
"""
struct ChunkField{K,C}
    ctx::K
    root::C
    range::UnitRange{Int}
    z::Vector{Float64}     # chunk values, z[p - offset]
    offset::Int
    halopos::Vector{Int}   # ascending grid positions of the halo cells
    haloz::Vector{Float64}
    halocells::Vector{C}
    misses::Base.RefValue{Int}
end

function ChunkField(k::GridCtx, sys, root, level::Integer, z_whole::Vector{Float64};
        halo_connectivity::Connectivity = k.connectivity, halo = nothing)
    r = chunk_range(sys, root, level)
    hc = halo === nothing ?
        DGG.subtree_halo(sys, root, level; connectivity = halo_connectivity) :
        collect(halo)
    hp = [DGG.cellposition(k.grid, c) for c in hc]
    perm = sortperm(hp)
    hp = hp[perm]; hc = hc[perm]
    return ChunkField(k, root, r, z_whole[r], first(r) - 1, hp, z_whole[hp], hc, Ref(0))
end

const Field = Union{WholeField,ChunkField}

ctx(f::Field) = f.ctx
positions(f::WholeField) = 1:length(f.ctx)
positions(f::ChunkField) = f.range
Base.length(f::Field) = length(positions(f))
reset_misses!(f::Field) = (f.misses[] = 0; f)

@inline value(f::WholeField, p::Int) = @inbounds f.z[p]

@inline function value(f::ChunkField, p::Int)
    if first(f.range) <= p <= last(f.range)
        return @inbounds f.z[p - f.offset]
    end
    i = searchsortedfirst(f.halopos, p)
    if i <= length(f.halopos) && @inbounds(f.halopos[i]) == p
        return @inbounds f.haloz[i]
    end
    f.misses[] += 1
    return NaN
end

# ---------------------------------------------------------------------------
# 1. Maximum downward gradient  (the headline function)
# ---------------------------------------------------------------------------

"""
    max_downward_gradient(f; degrees = true)

For each cell, the steepest DOWNHILL gradient to any one-ring neighbour:

    max over neighbours n of  (z(c) - z(n)) / d(c, n)

with `d` the great-circle distance between cell centroids on the unit sphere,
in radians. Pits (no lower neighbour) give exactly `0`. With `degrees = true`
the result is `atand` of that, i.e. a slope angle.

Neighbour degree varies — 5 at pentagons, 6 for hexagons, 3..11 on A5, 7 or 8
at HEALPix and S2 corners — and this needs no knowledge of that: it iterates
whatever `neighbors` returns. Non-uniform spacing is handled by dividing by the
actual centroid distance rather than by a nominal cell size.
"""
function max_downward_gradient(f::Field; degrees::Bool = true)
    k = f.ctx
    out = Vector{Float64}(undef, length(f))
    @inbounds for (i, p) in enumerate(positions(f))
        zc = value(f, p); pc = k.centroids[p]
        g = 0.0
        for j in nrange(k, p)
            q = k.nbr[j]
            zn = value(f, q)
            d = gcdistance(pc, k.centroids[q])
            d == 0 && continue
            cand = (zc - zn) / d
            cand > g && (g = cand)
            isnan(zn) && (g = NaN)
        end
        out[i] = degrees ? atand(g) : g
    end
    return out
end

# ---------------------------------------------------------------------------
# 2. The rest of the PR's local metrics
# ---------------------------------------------------------------------------

"""
    slope_mdg(f)

PR #83's `slope!(::MaximumDownwardGradient, ...)` literally: the maximum
*absolute* inter-cell gradient, in degrees. (Uphill neighbours count too, so
this is not a downhill slope; kept for fidelity to the PR.)
"""
function slope_mdg(f::Field)
    k = f.ctx
    out = Vector{Float64}(undef, length(f))
    @inbounds for (i, p) in enumerate(positions(f))
        zc = value(f, p); pc = k.centroids[p]; g = 0.0
        for j in nrange(k, p)
            q = k.nbr[j]; zn = value(f, q)
            d = gcdistance(pc, k.centroids[q])
            d == 0 && continue
            cand = abs(zn - zc) / d
            cand > g && (g = cand)
            isnan(zn) && (g = NaN)
        end
        out[i] = atand(g)
    end
    return out
end

"Largest absolute inter-cell difference (Wilson et al. roughness)."
function roughness(f::Field)
    k = f.ctx
    out = Vector{Float64}(undef, length(f))
    @inbounds for (i, p) in enumerate(positions(f))
        zc = value(f, p); v = 0.0
        for j in nrange(k, p)
            zn = value(f, k.nbr[j]); d = abs(zn - zc)
            d > v && (v = d)
            isnan(zn) && (v = NaN)
        end
        out[i] = v
    end
    return out
end

"Topographic Position Index: centre minus the mean of its one-ring."
function tpi(f::Field)
    k = f.ctx
    out = Vector{Float64}(undef, length(f))
    @inbounds for (i, p) in enumerate(positions(f))
        zc = value(f, p); total = 0.0; n = 0
        for j in nrange(k, p)
            total += value(f, k.nbr[j]); n += 1
        end
        out[i] = n == 0 ? NaN : zc - total / n
    end
    return out
end

"Terrain Ruggedness Index."
function tri(f::Field; normalize::Bool = false, squared::Bool = true)
    k = f.ctx
    out = Vector{Float64}(undef, length(f))
    @inbounds for (i, p) in enumerate(positions(f))
        zc = value(f, p); v = 0.0; n = 0
        for j in nrange(k, p)
            d = abs(value(f, k.nbr[j]) - zc)
            v += squared ? d * d : d
            n += 1
        end
        normalize && n > 0 && (v /= n)
        out[i] = squared ? sqrt(v) : v
    end
    return out
end

"Number of one-ring neighbours no higher than the centre (0 = pit)."
function prominence(f::Field)
    k = f.ctx
    out = Vector{Int}(undef, length(f))
    @inbounds for (i, p) in enumerate(positions(f))
        zc = value(f, p); n = 0
        for j in nrange(k, p)
            zn = value(f, k.nbr[j])
            isnan(zn) && (n -= 1000)         # poison, so a miss is visible
            zn <= zc && (n += 1)
        end
        out[i] = n
    end
    return out
end

"""
    _plane_fit(f, p) -> (east_gradient, north_gradient, ok)

Least-squares fit of a tilted plane to the one-ring, in the local ENU tangent
frame — PR #83's `_localaspect`, with great-circle distance and bearing in
place of a raster's cell size. `ok = false` marks a degenerate normal matrix.
"""
@inline function _plane_fit(f::Field, p::Int)
    k = f.ctx
    xx = xy = yy = xz = yz = 0.0
    zc = value(f, p); pc = k.centroids[p]; bad = false
    for j in nrange(k, p)
        q = k.nbr[j]; pn = k.centroids[q]
        d = gcdistance(pc, pn)
        d == 0 && continue
        b = bearing_deg(pc, pn)
        e = d * sind(b); nn = d * cosd(b)
        zn = value(f, q)
        isnan(zn) && (bad = true)
        dz = zn - zc
        xx += e * e; xy += e * nn; yy += nn * nn
        xz += e * dz; yz += nn * dz
    end
    det = xx * yy - xy * xy
    (det == 0 || !isfinite(det)) && return (NaN, NaN, false)
    bad && return (NaN, NaN, true)
    return ((xz * yy - yz * xy) / det, (yz * xx - xz * xy) / det, true)
end

"Aspect in degrees clockwise from north, pointing DOWNSLOPE. `NaN` when flat."
function aspect(f::Field)
    out = Vector{Float64}(undef, length(f))
    @inbounds for (i, p) in enumerate(positions(f))
        eg, ng, ok = _plane_fit(f, p)
        out[i] = (!ok || (eg == 0 && ng == 0)) ? NaN : mod(atand(-eg, -ng), 360.0)
    end
    return out
end

"Slope in degrees from the same least-squares plane fit `aspect` uses."
function planefit_slope(f::Field)
    out = Vector{Float64}(undef, length(f))
    @inbounds for (i, p) in enumerate(positions(f))
        eg, ng, ok = _plane_fit(f, p)
        out[i] = ok ? atand(hypot(eg, ng)) : NaN
    end
    return out
end

"""
    curvature(f)

A discrete spherical Laplacian: `4 * mean over neighbours of (z(n) - z(c)) / d²`.
Positive is concave up (a valley), negative convex (a ridge). Uses actual
centroid distances, so it degrades gracefully where cells are unequal.
"""
function curvature(f::Field)
    k = f.ctx
    out = Vector{Float64}(undef, length(f))
    @inbounds for (i, p) in enumerate(positions(f))
        zc = value(f, p); pc = k.centroids[p]; acc = 0.0; n = 0
        for j in nrange(k, p)
            q = k.nbr[j]
            d = gcdistance(pc, k.centroids[q])
            d == 0 && continue
            acc += (value(f, q) - zc) / (d * d)
            n += 1
        end
        out[i] = n == 0 ? NaN : 4 * acc / n
    end
    return out
end

# ---------------------------------------------------------------------------
# 3. Flow
# ---------------------------------------------------------------------------

"""
    flow_direction(f) -> Vector{Int}

D8-style single flow direction generalised to a DGGS: the GRID POSITION of the
neighbour with the steepest downhill gradient. A pit points at itself. Ties
break on ascending position, so the answer does not depend on the rotational
phase of `neighbors`.

Absolute positions rather than a packed direction code — the deliberate
alternative to PR #83's relative-index type.
"""
function flow_direction(f::Field)
    k = f.ctx
    out = Vector{Int}(undef, length(f))
    @inbounds for (i, p) in enumerate(positions(f))
        zc = value(f, p); pc = k.centroids[p]
        best = p; bestg = 0.0
        for j in nrange(k, p)
            q = k.nbr[j]; zn = value(f, q)
            d = gcdistance(pc, k.centroids[q])
            d == 0 && continue
            g = (zc - zn) / d
            if g > bestg || (g == bestg && g > 0 && q < best)
                best = q; bestg = g
            end
            isnan(zn) && (best = -1)
        end
        out[i] = best
    end
    return out
end

"""
    flow_accumulation(f::WholeField) -> Vector{Float64}

Single-direction (D8) flow accumulation: every cell starts with one unit and
pushes its total downhill in descending-elevation order.

Only defined for a whole field, because it is a GLOBAL algorithm — a cell's
value depends on its whole upslope basin, not on a one-ring. That is itself a
finding: a one-ring halo does not make the interesting half of PR #83
chunkable.
"""
function flow_accumulation(f::WholeField)
    dir = flow_direction(f)
    n = length(f)
    acc = ones(Float64, n)
    for i in sortperm(f.z; rev = true)
        d = dir[i]
        d == i && continue
        acc[d] += acc[i]
    end
    return acc
end

# ---------------------------------------------------------------------------
# 4. Synthetic elevation fields
# ---------------------------------------------------------------------------

"Degree-3 spherical-harmonic-ish smooth field; nothing flat, nothing jumping."
function smooth_harmonic(k::GridCtx)
    n = length(k)
    z = Vector{Float64}(undef, n)
    for p in 1:n
        x, y, zz = k.centroids[p]
        z[p] = 3 * (x * x - y * y) * zz + 2 * x * y * (3 * zz * zz - 1) + 0.5 * zz
    end
    return z
end

"A single 1000-unit spike at position `at`, zero elsewhere."
function spike_field(k::GridCtx, at::Int)
    z = zeros(Float64, length(k)); z[at] = 1000.0; return z
end

"A 0/1000 step across the plane `n · x = 0`."
function step_field(k::GridCtx, normal = (0.3, 0.7, 0.6))
    n = length(k)
    z = Vector{Float64}(undef, n)
    for p in 1:n
        z[p] = dot3(k.centroids[p], normal) > 0 ? 1000.0 : 0.0
    end
    return z
end

"Linear ramp in latitude, degrees."
function lat_ramp(k::GridCtx)
    n = length(k)
    z = Vector{Float64}(undef, n)
    for p in 1:n
        z[p] = asind(clamp(k.centroids[p][3], -1.0, 1.0))
    end
    return z
end

"Reproducible white noise — maximally sensitive to a wrong neighbour value."
white_noise(k::GridCtx, rng) = rand(rng, length(k)) .* 1000

const FIELD_KINDS = (:harmonic, :spike, :step, :ramp, :noise)

"""
    make_field(kind, ctx, rng; spike_at = 1) -> Vector{Float64}
"""
function make_field(kind::Symbol, k::GridCtx, rng; spike_at::Int = 1)
    kind === :harmonic && return smooth_harmonic(k)
    kind === :spike    && return spike_field(k, spike_at)
    kind === :step     && return step_field(k)
    kind === :ramp     && return lat_ramp(k)
    kind === :noise    && return white_noise(k, rng)
    throw(ArgumentError("unknown field kind $kind"))
end

# ---------------------------------------------------------------------------
# 5. Registry + comparison used by the checkers
# ---------------------------------------------------------------------------

const METRICS = (
    (:max_downward_gradient, max_downward_gradient),
    (:slope_mdg,             slope_mdg),
    (:roughness,             roughness),
    (:tpi,                   tpi),
    (:tri,                   tri),
    (:prominence,            prominence),
    (:aspect,                aspect),
    (:planefit_slope,        planefit_slope),
    (:curvature,             curvature),
    (:flow_direction,        flow_direction),
)

isnan_(x) = x isa AbstractFloat && isnan(x)

"NaN-aware exact equality, elementwise."
function same(a::AbstractVector, b::AbstractVector)
    length(a) == length(b) || return false
    for i in eachindex(a, b)
        (isnan_(a[i]) && isnan_(b[i])) && continue
        a[i] == b[i] || return false
    end
    return true
end

"Index and values of the first elementwise difference, or `(0, nothing, nothing)`."
function first_difference(a::AbstractVector, b::AbstractVector)
    length(a) == length(b) || return (-1, length(a), length(b))
    for i in eachindex(a, b)
        (isnan_(a[i]) && isnan_(b[i])) && continue
        a[i] == b[i] || return (i, a[i], b[i])
    end
    return (0, nothing, nothing)
end

end # module
