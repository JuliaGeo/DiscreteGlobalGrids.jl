# # Turning cells into one surface
#
# `tessellate.jl` draws the cells.  This file draws the field they sample: a
# vertex per cell **centroid**, carrying that cell's value, joined by the
# triangles of the grid's dual so that the GPU interpolates between cell
# centres.  Where `n` cells meet at a grid corner their centroids ring it and
# fan into `n - 2` triangles; over every corner, those triangles tile the sphere
# exactly once.
#
# ## Finding the corners
#
# There is no verb that lists corners, and cell boundaries are floating-point
# rings that neighbours do not agree on bit for bit, so coordinates cannot be
# matched either.  Adjacency is stated exactly, though, and the cells at a
# corner all touch one another, which puts them in each other's rings as a run
# of consecutive neighbours.  For a cell `a` and a consecutive ring pair
# `(b, c)`:
#
#  1. `b` and `c` must touch each other, or the three share no corner.  On a
#     partial grid this also stops a triangle bridging a hole, where two members
#     either side of a missing one come out consecutive.
#  2. Widen `(b, c)` to the largest run of `a`'s neighbours that all touch one
#     another.  That run plus `a` is the corner's cell set.
#  3. Emit only if `a` is the smallest member of that set.
#
# Step 3 gives each corner one owner, so no triangle is emitted twice without a
# shared hash set or a later `unique`.  Emitting per pair rather than per run
# keeps that true where two maximal runs overlap.
#
# Step 2 is what makes one rule right for square cells as well as hexagons.  A
# corner of `n` cells is a run of `n - 1`, fanned into `n - 2` triangles: one at
# the three-cell corners of IGeo7 and H3, two at the four-cell corners of
# HEALPix, S2 and ISEA4R, three where ISEA4R's diamonds meet an icosahedral
# vertex.  Without the widening, a four-cell corner would be emitted as four
# overlapping triangles instead of two.

"""
    ZeroHeights(n)

`n` heights, all zero, stored as the number `n` and nothing else.

The stand-in for the `zs` a caller did not give.  It is an ordinary
`AbstractVector{Float64}`, so anything that reads heights reads it, and
[`triangulate`](@ref) recognises the type: a surface at height zero everywhere
keeps the two coordinates it would have had rather than carrying a third that is
always `0.0`.
"""
struct ZeroHeights <: AbstractVector{Float64}
    n::Int
end

Base.size(z::ZeroHeights) = (z.n,)

Base.@propagate_inbounds function Base.getindex(z::ZeroHeights, i::Int)
    @boundscheck checkbounds(z, i)
    return 0.0
end


# `adjacency` gives counter-clockwise rings addressed by in-region position,
# which is also the index the vertex buffer and the colour vector use.

@inline function inring(r, b::Int)
    @inbounds for x in r
        x == b && return true
    end
    return false
end

"""
    touches_all(adj, r, k, lo, hi, m) -> Bool

Does ring member `m` touch every member of the cyclic run `lo:hi`?

The run is already pairwise touching, so only the new member has to be checked.
The scan starts at the far end, which is the one likely to answer `false`: in a
hexagonal ring the near end always touches.
"""
@inline function touches_all(adj, r, k::Int, lo::Int, hi::Int, m::Int)
    @inbounds y = r[mod1(m, k)]
    step = m > hi ? (lo:hi) : (hi:-1:lo)
    @inbounds for i in step
        inring(adj[r[mod1(i, k)]], y) || return false
    end
    return true
end

"""
    corner_run(adj, r, k, i) -> (lo, hi)

The corner that the consecutive ring pair `(i, i + 1)` belongs to: the widest
cyclic run of the ring `r` whose members all touch one another.

The run grows outwards from the pair, forwards first, so where two maximal runs
overlap this picks one of them.  The caller emits per pair, which keeps that
choice from mattering.
"""
@inline function corner_run(adj, r, k::Int, i::Int)
    lo, hi = i, i + 1
    while hi - lo + 1 < k && touches_all(adj, r, k, lo, hi, hi + 1)
        hi += 1
    end
    while hi - lo + 1 < k && touches_all(adj, r, k, lo, hi, lo - 1)
        lo -= 1
    end
    return lo, hi
end

"""
    owns_corner(r, k, lo, hi, a) -> Bool

Is `a` the smallest cell of the corner whose other members are `r[lo:hi]`?
"""
@inline function owns_corner(r, k::Int, lo::Int, hi::Int, a::Int)
    @inbounds for j in lo:hi
        a < r[mod1(j, k)] || return false
    end
    return true
end

# A task owns a contiguous block of cells and writes triangles for the corners
# those cells own.  Face indices address `[all centroids; this task's extra
# vertices]`, so an index of `n` or less names a centroid and needs no fixing up
# at merge time, and a larger one is task-local and gets shifted.

"""
    CornerWeights

How much of each of a triangle's three cells a vertex is made of.

A vertex at a corner is that corner outright — `(1, 0, 0)` and its rotations.  A
vertex the seam created lies partway along an edge, and is the two ends mixed in
exactly the proportion the cut fell at, which is the value the GPU would have
interpolated there had the triangle not been cut.  Carrying the mix rather than
the nearer end is what keeps a cut triangle continuous with the one across the
seam from it, in height as well as in colour.
"""
const CornerWeights = NTuple{3, Float64}

const CORNER1 = (1.0, 0.0, 0.0)
const CORNER2 = (0.0, 1.0, 0.0)
const CORNER3 = (0.0, 0.0, 1.0)

@inline mix(a::CornerWeights, b::CornerWeights, t::Float64) =
    (a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2]), a[3] + t * (b[3] - a[3]))

struct SurfaceChunk{P}
    faces::Vector{NTuple{3, Int32}}
    extra::Vector{P}                    # vertices private to a cut or polar triangle
    extra_tri::Vector{NTuple{3, Int32}} # the three cells that vertex's triangle spans
    extra_weight::Vector{NTuple{3, Float32}}  # how much of each it is
    poly::Vector{Point2d}     # scratch: a triangle being cut
    polyw::Vector{CornerWeights}
    clip::Vector{Point2d}     # scratch: the surviving side of a cut
    clipw::Vector{CornerWeights}
    ring::Vector{Point2d}     # scratch: a triangle's outline, traced
    ringw::Vector{CornerWeights}
end

function SurfaceChunk{P}(nhint::Int) where {P}
    faces = NTuple{3, Int32}[]
    sizehint!(faces, 2 * nhint)
    return SurfaceChunk{P}(faces, P[], NTuple{3, Int32}[], NTuple{3, Float32}[],
        Point2d[], CornerWeights[], Point2d[], CornerWeights[],
        Point2d[], CornerWeights[])
end

@inline function push_face!(chunk::SurfaceChunk, a::Int32, b::Int32, c::Int32)
    push!(chunk.faces, (a, b, c))
    return
end

"""
    emit_private_fan!(chunk, pts, ws, m, n, tri) -> ntriangles

Add `m` vertices of a cut or polar triangle as this task's own, and fan them.

`n` is the number of shared centroid vertices, where task-local numbering starts.
`tri` is the three cells the triangle spans, which the vertices' weights are
weights of.
"""
function emit_private_fan!(chunk::SurfaceChunk{P}, pts::Vector{P},
        ws::Vector{CornerWeights}, m::Int, n::Int, tri::NTuple{3, Int32}) where {P}
    m < 3 && return 0
    base = Int32(n + length(chunk.extra))
    @inbounds for k in 1:m
        w = ws[k]
        push!(chunk.extra, pts[k])
        push!(chunk.extra_tri, tri)
        push!(chunk.extra_weight, (Float32(w[1]), Float32(w[2]), Float32(w[3])))
    end
    @inbounds for k in 2:(m - 1)
        push_face!(chunk, base + Int32(1), base + Int32(k), base + Int32(k + 1))
    end
    return m - 2
end

# Nothing to cut and nothing to cap, so every triangle is three shared
# centroids and the task writes faces only.

function fill_surface!(chunk::SurfaceChunk, ::GlobeTarget, adj, lo::Int, hi::Int, n::Int)
    @inbounds for a in lo:hi
        r = adj[a]
        k = length(r)
        k < 2 && continue
        # A ring of two has one pair: `(r₁, r₂)` and `(r₂, r₁)` are the same
        # corner, and walking both would draw its triangle twice.
        npairs = k == 2 ? 1 : k
        for i in 1:npairs
            b = r[i]
            c = r[mod1(i + 1, k)]
            inring(adj[b], c) || continue
            rlo, rhi = corner_run(adj, r, k, i)
            owns_corner(r, k, rlo, rhi, a) || continue
            push_face!(chunk, Int32(a), Int32(b), Int32(c))
        end
    end
    return 0
end

# The centroids already sit in the map's window, so a triangle is shared-vertex
# work unless longitude says otherwise.  How far longitude turns in going once
# round tells the two exceptions apart: `0` for one straddling the map's cut,
# `±360` for the one at each pole that contains it.

"""
    outline_triangle!(pts, ws, u1, u2, u3) -> winding

Fill `pts` with the triangle's outline in longitude/latitude degrees and return
how far longitude turned in going once around.

The outline is the three corners and nothing else: every edge is shared, the
triangle across it draws that edge as the straight lon/lat segment between the
same two corners, and bending one copy and not the other opens a gap.  Only the
longitude sweep along each edge is worked out rather than read off.

For a sweep of a quarter turn or more, the direction is the sign of `(a × b)ᶻ`,
longitude being monotonic along any great circle that is not a meridian.  The
sign is anti-symmetric — the triangle across the edge gets it negated bit for
bit, however nearly zero it is — so the two never disagree about the edge.
Below a quarter turn the sign is noise and the short way is right anyway.

`(a × b)ᶻ` is exactly zero when the edge's great circle passes through both
poles; if the corners are then half a turn apart, the arc runs over one, and the
two pole corners are inserted so that the triangle closes along the top of the
map.  Grids whose cells meet four to a corner put a pole on a dual edge exactly.
The jump follows the triangle's reference longitude, so the triangle across the
edge takes the other half turn and between them they cover it once.
"""
function outline_triangle!(pts::Vector{Point2d}, ws::Vector{CornerWeights}, u1, u2, u3)
    empty!(pts)
    empty!(ws)
    us = (u1, u2, u3)
    cs = (CORNER1, CORNER2, CORNER3)
    reference = reference_longitude(us, 3)
    previous = atand(u1[2], u1[1])
    first_lon = previous
    @inbounds for k in 1:3
        j = k == 3 ? 1 : k + 1
        a = us[k]
        b = us[j]
        push!(pts, Point2d(previous, asind(clamp(a[3], -1.0, 1.0))))
        push!(ws, cs[k])

        lon_b = atand(b[2], b[1])
        step = lon_delta(lon_b, previous)
        if abs(step) >= 90.0
            nz = a[1] * b[2] - a[2] * b[1]   # (a × b)ᶻ
            if nz > 0
                step < 0 && (step += 360.0)  # eastward
            elseif nz < 0
                step > 0 && (step -= 360.0)  # westward
            else
                # A meridian great circle, and the corners half a turn apart:
                # the arc runs over a pole.
                pole = a[3] + b[3] >= 0 ? 90.0 : -90.0
                forward = lon_distance(previous + 90.0, reference) <=
                    lon_distance(previous - 90.0, reference)
                step = forward ? 180.0 : -180.0
                push!(pts, Point2d(previous, pole))
                push!(ws, cs[k])
                push!(pts, Point2d(previous + step, pole))
                push!(ws, cs[j])
            end
        end
        previous += step
    end
    return previous - first_lon
end

"""
    emit_polar_band!(chunk, m, n, north, left, cut, tri) -> ntriangles

Draw a triangle that *contains* a pole as the polar cap it is.

Such a triangle winds a full turn in longitude, so cutting it at the seam gives
no cap.  It is instead the band between its outline — the `m` points of
`chunk.ring`, monotonic in longitude, increasing at the north pole and
decreasing at the south — and the pole, spanning the map and closed along the
top or bottom edge.

One trapezoid per outline segment rather than a fan from one corner: a trapezoid
stays convex whatever the latitudes do, where a fan can turn a triangle inside
out.  Each is emitted at two turns of longitude and clipped to the window, which
leaves the band's lower edge on the longitudes its neighbours draw.
"""
function emit_polar_band!(chunk::SurfaceChunk{Point2d}, m::Int, n::Int, north::Bool,
        left::Float64, cut::Float64, tri::NTuple{3, Int32})
    ring = chunk.ring
    ringw = chunk.ringw
    # Longitude increases around the north pole and decreases around the south.
    # Turning the southern outlines round lets one strip serve both.
    if !north
        reverse!(ring)
        reverse!(ringw)
    end
    # There is no data at the pole itself, so its corners take the value of
    # whichever cell reaches nearest to it.
    polew = ringw[1]
    highest = abs(ring[1][2])
    @inbounds for k in 2:m
        if abs(ring[k][2]) > highest
            highest = abs(ring[k][2])
            polew = ringw[k]
        end
    end
    pole = north ? 90.0 : -90.0

    # Put the first corner in the window, so that the two turns emitted below
    # are certain to cover the whole of it between them.
    @inbounds base_shift = (left + mod(ring[1][1] - left, 360.0)) - ring[1][1]

    pts = chunk.poly
    ws = chunk.polyw
    ntri = 0
    for turn in (-360.0, 0.0)
        @inbounds for k in 1:m
            j = k == m ? 1 : k + 1
            p = ring[k]
            q = ring[j]
            x0 = p[1] + base_shift + turn
            # The segment that closes the outline lands one full turn along.
            x1 = (k == m ? q[1] + 360.0 : q[1]) + base_shift + turn
            x0 == x1 && continue
            empty!(pts)
            empty!(ws)
            # Wound counter-clockwise either way: at the north pole the cap lies
            # above its lower edge, at the south pole below it.
            if north
                push!(pts, Point2d(x0, p[2]));  push!(ws, ringw[k])
                push!(pts, Point2d(x1, q[2]));  push!(ws, ringw[j])
                push!(pts, Point2d(x1, pole));  push!(ws, polew)
                push!(pts, Point2d(x0, pole));  push!(ws, polew)
            else
                push!(pts, Point2d(x1, q[2]));  push!(ws, ringw[j])
                push!(pts, Point2d(x0, p[2]));  push!(ws, ringw[k])
                push!(pts, Point2d(x0, pole));  push!(ws, polew)
                push!(pts, Point2d(x1, pole));  push!(ws, polew)
            end
            ntri += emit_window_polygon!(chunk, 4, n, left, cut, tri)
        end
    end
    return ntri
end

"""
    emit_traced_polygon!(chunk, m, n, left, cut, tri) -> ntriangles

Draw an ordinary triangle whose outline, in `chunk.ring`, does not fit the map
in one piece.

The outline is emitted at three turns of longitude and each copy clipped to the
window, so a triangle straddling the cut comes out as its two pieces without
either having to be identified as the far side.  Three, because the outline can
run either way in longitude from the corner the turns are measured off.
"""
function emit_traced_polygon!(chunk::SurfaceChunk{Point2d}, m::Int, n::Int,
        left::Float64, cut::Float64, tri::NTuple{3, Int32})
    ring = chunk.ring
    ringw = chunk.ringw
    @inbounds base_shift = (left + mod(ring[1][1] - left, 360.0)) - ring[1][1]
    pts = chunk.poly
    ws = chunk.polyw
    ntri = 0
    for turn in (-360.0, 0.0, 360.0)
        empty!(pts)
        empty!(ws)
        @inbounds for k in 1:m
            push!(pts, Point2d(ring[k][1] + base_shift + turn, ring[k][2]))
            push!(ws, ringw[k])
        end
        ntri += emit_window_polygon!(chunk, m, n, left, cut, tri)
    end
    return ntri
end

"""
    emit_window_polygon!(chunk, m, n, left, cut, tri) -> ntriangles

Clip the `m`-gon in `chunk.poly` to the map's window `[left, cut]` and fan what
survives.  The two clips ping-pong between the chunk's scratch buffers, leaving
the result in `chunk.poly`.
"""
function emit_window_polygon!(chunk::SurfaceChunk{Point2d}, m::Int, n::Int,
        left::Float64, cut::Float64, tri::NTuple{3, Int32})
    kept = clip_weighted!(chunk.clip, chunk.clipw, chunk.poly, chunk.polyw,
        m, left, false, 0.0)
    kept < 3 && return 0
    kept = clip_weighted!(chunk.poly, chunk.polyw, chunk.clip, chunk.clipw,
        kept, cut, true, 0.0)
    return emit_private_fan!(chunk, chunk.poly, chunk.polyw, kept, n, tri)
end

"""
    clip_weighted!(out, outw, pts, ws, m, cut, keep_below, shift) -> nkept

Sutherland–Hodgman against one vertical half-plane in longitude, carrying each
vertex's [`CornerWeights`](@ref) along with it, with `shift` added to the
surviving longitudes so that the far side of a cut triangle lands on the far edge
of the map.

A crossing landing exactly on a corner is skipped: the corner is already there,
and repeating it would put a zero-area triangle in the fan.

A vertex created *on* the cut is its edge's two ends mixed at the crossing
parameter, which is the value the triangle carried at that point before it was
cut.  The mix costs one lerp per created vertex, and there are only as many of
those as the seam and the poles make.
"""
function clip_weighted!(out::Vector{Point2d}, outw::Vector{CornerWeights},
        pts::Vector{Point2d}, ws::Vector{CornerWeights}, m::Int, cut::Float64,
        keep_below::Bool, shift::Float64)
    empty!(out)
    empty!(outw)
    @inbounds for i in 1:m
        j = i == m ? 1 : i + 1
        a = pts[i]; b = pts[j]
        wa = ws[i]; wb = ws[j]
        a_in = keep_below ? (a[1] <= cut) : (a[1] >= cut)
        b_in = keep_below ? (b[1] <= cut) : (b[1] >= cut)
        if a_in
            push!(out, Point2d(a[1] + shift, a[2]))
            push!(outw, wa)
        end
        if a_in != b_in
            t = (cut - a[1]) / (b[1] - a[1])
            # `t` at either end means the crossing IS one of the corners, which
            # is pushed as a corner in its own right; adding it again would put
            # a zero-area triangle in the fan.
            if 0.0 < t < 1.0
                push!(out, Point2d(cut + shift, a[2] + t * (b[2] - a[2])))
                push!(outw, mix(wa, wb, t))
            end
        end
    end
    return length(out)
end

function fill_surface!(chunk::SurfaceChunk{Point2d}, target::PlanarTarget, adj,
        lo::Int, hi::Int, n::Int, positions::Vector{Point2d}, cr::CellRegion)
    cut = target.cut
    left = cut - 360.0
    wrapping = needs_cutting(target)
    source = cr.source
    cells = cr.cells
    nsplit = 0
    @inbounds for a in lo:hi
        r = adj[a]
        k = length(r)
        k < 2 && continue
        # A ring of two has one pair: `(r₁, r₂)` and `(r₂, r₁)` are the same
        # corner, and walking both would draw its triangle twice.
        npairs = k == 2 ? 1 : k
        for i in 1:npairs
            b = r[i]
            c = r[mod1(i + 1, k)]
            inring(adj[b], c) || continue
            rlo, rhi = corner_run(adj, r, k, i)
            owns_corner(r, k, rlo, rhi, a) || continue

            if !wrapping
                push_face!(chunk, Int32(a), Int32(b), Int32(c))
                continue
            end

            l2, l3, span = unwrap_triangle(positions, a, b, c)
            # Short in longitude, with all three corners already together in
            # the window: three shared vertices and nothing to clip.
            if span < 90.0 && l2 == positions[b][1] && l3 == positions[c][1]
                push_face!(chunk, Int32(a), Int32(b), Int32(c))
                continue
            end

            # Everything else is drawn from its outline, and the winding says
            # whether it straddles the cut or lies over a pole.
            winding = outline_triangle!(chunk.ring, chunk.ringw,
                DGG.cell_centroid(source, cells[a]),
                DGG.cell_centroid(source, cells[b]),
                DGG.cell_centroid(source, cells[c]))
            m = length(chunk.ring)
            tri = (Int32(a), Int32(b), Int32(c))
            if abs(winding) > 180.0
                # Longitude increases around the north pole and decreases around
                # the south.
                emit_polar_band!(chunk, m, n, winding > 0, left, cut, tri)
            else
                emit_traced_polygon!(chunk, m, n, left, cut, tri)
            end
            nsplit += 1
        end
    end
    return nsplit
end

"""
    unwrap_triangle(positions, a, b, c) -> (l2, l3, span)

The triangle's second and third longitudes made continuous with its first, and
the longest step taken to get around.

`positions` holds centroids in the map's window, so a corner whose continuous
longitude is not the one it is stored at belongs to a triangle straddling the
cut.  A step of a quarter turn or more means an edge beside a pole, which
[`outline_triangle!`](@ref) has to sweep from the edge itself.
"""
@inline function unwrap_triangle(positions::Vector{Point2d}, a::Int, b::Int, c::Int)
    @inbounds l1 = positions[a][1]
    @inbounds d2 = lon_delta(positions[b][1], l1)
    l2 = l1 + d2
    @inbounds d3 = lon_delta(positions[c][1], l2)
    l3 = l2 + d3
    d1 = lon_delta(l1, l3)
    span = max(abs(d1), abs(d2), abs(d3))
    return l2, l3, span
end

"""
    surface_vertex(target, p, z) -> Point2d or Point3d

A cell centroid `p`, a unit-sphere point, raised to height `z`, in the target's
build space.

On a [`GlobeTarget`](@ref) the build space is the unit sphere itself, and `z` is
carried no further: a globe vertex is finished by [`vertex_positions`](@ref),
which is what lets the heights change without the topology being rebuilt.

On a [`PlanarTarget`](@ref) it is longitude and latitude in degrees inside the
map's window `[cut - 360, cut]`, which [`surface_topology`](@ref) projects in
bulk at the end; `z` is not part of that space either, because every seam and
polar cut is worked out in longitude and latitude alone.
"""
@inline surface_vertex(::GlobeTarget, p) = Point3d(p[1], p[2], p[3])

@inline function surface_vertex(target::PlanarTarget, p)
    lon = atand(p[2], p[1])
    lat = asind(clamp(p[3], -1.0, 1.0))
    needs_cutting(target) || return Point2d(lon, lat)
    left = target.cut - 360.0
    return Point2d(left + mod(lon - left, 360.0), lat)
end

function surface_vertices(target::PlotTarget, cr::CellRegion, ntasks::Int)
    P = pointtype(target)
    n = length(cr)
    positions = Vector{P}(undef, n)
    cells = cr.cells
    source = cr.source
    inparallel(n, ntasks) do lo, hi
        @inbounds for p in lo:hi
            positions[p] = surface_vertex(target, DGG.cell_centroid(source, cells[p]))
        end
    end
    return positions
end

"""
    SurfaceTopology(target, cells, positions, faces, extra_tri, extra_weight,
                    ncells, nsplit)

Everything about a DGGS surface that its heights do not touch.

  * `target` — the space it was built in.
  * `cells` — the [`CellRegion`](@ref) it was built from.  Held so that a plot
    can tell a new cell set from the one it already has: the pipeline reports
    every update as a change, and rebuilding this is the expensive half.
  * `positions` — the vertices in **build space**: longitude and latitude
    projected into the axis's plane on a flat map, and the unit-sphere direction
    on a globe, which is as far as a vertex gets before its height is known.
  * `faces` — the triangles of the grid's dual.
  * `extra_tri`, `extra_weight` — for each vertex past the first `ncells`, the
    three cells of the triangle that created it and its [`CornerWeights`](@ref)
    within it.  The first `ncells` vertices are the cells themselves and need
    neither.
  * `ncells` — the length a value vector must have.
  * `nsplit` — triangles cut at the seam or capped at a pole.  Zero on a globe.

Built by [`surface_topology`](@ref); finished into a [`SurfaceMesh`](@ref) by
[`vertex_positions`](@ref).
"""
struct SurfaceTopology{P, T <: PlotTarget, R <: CellRegion}
    target::T
    cells::R
    positions::Vector{P}
    faces::Vector{GLTriangleFace}
    extra_tri::Vector{NTuple{3, Int32}}
    extra_weight::Vector{NTuple{3, Float32}}
    ncells::Int
    nsplit::Int
end

Base.isempty(t::SurfaceTopology) = isempty(t.faces)

ntriangles(t::SurfaceTopology) = length(t.faces)

nvertices(t::SurfaceTopology) = length(t.positions)

"""
    samebuild(top, target, cells) -> Bool

Would rebuilding `top` from `target` and `cells` produce the same topology?

The compute graph cannot answer it: for a value that is not `isbits` it reports
`a === b` as *changed*, on the grounds that the same object may have been mutated
behind its back, so the surface asks for itself.  A target is either immutable to
its leaves or an object whose identity is the thing that matters, and `===` is
the test for it.  The cells are compared with `==`, which for a
[`CellRegion`](@ref) is a value comparison: a plot handed its arguments a second
time may get a set rebuilt from the same cells rather than the same object, and
`===` would call that a different surface and rebuild the expensive half.
"""
samebuild(top::SurfaceTopology, target::PlotTarget, cells::CellRegion) =
    top.target === target && top.cells == cells

samebuild(::SurfaceTopology, ::PlotTarget, ::Any) = false

"""
    blend(v1, v2, v3, w) -> value

One vertex's value, from the three cells of its triangle and its weights.

Numbers are mixed.  Anything else — a colour, a category, a string — takes the
cell it is most of, because there is no meaning to two thirds of a colour that
the backends and the colormap would agree on.
"""
@inline blend(v1::Number, v2::Number, v3::Number, w::NTuple{3, Float32}) =
    Float64(w[1]) * v1 + Float64(w[2]) * v2 + Float64(w[3]) * v3

@inline function blend(v1, v2, v3, w::NTuple{3, Float32})
    w[1] >= w[2] && w[1] >= w[3] && return v1
    return w[2] >= w[3] ? v2 : v3
end

"""
    blendtype(T) -> Type

What a mix of `T`s has to be stored as.

A float stays itself: mixing loses no more than the values already carried.  Any
other number becomes a `Float64`, because two thirds of the way between two
integers is not an integer — a colour vector of `1:ncells` would otherwise fail
to hold its own seam vertices.  Anything that is not mixed at all keeps its type.
"""
blendtype(::Type{T}) where {T <: AbstractFloat} = T
blendtype(::Type{T}) where {T <: Number} = Float64
blendtype(::Type{T}) where {T} = T

"""
    spread(values, top) -> AbstractVector

One value per cell spread over one value per vertex.

The first `ncells` vertices are the cells in order, so a set with nothing cut —
a globe, or `wrap = false` — hands `values` straight back.  Each vertex past
them is a point inside a triangle, and takes that triangle's three values mixed
by its [`CornerWeights`](@ref): exactly the value the GPU would have interpolated
there had the triangle not been cut, so the two halves of a split triangle meet
without a step.
"""
function spread(values::AbstractVector, top::SurfaceTopology)
    nextra = length(top.extra_tri)
    nextra == 0 && return values
    length(values) == top.ncells ||
        throw(ArgumentError("got $(length(values)) values, but there are \
            $(top.ncells) cells"))
    out = Vector{blendtype(eltype(values))}(undef, top.ncells + nextra)
    copyto!(out, 1, values, firstindex(values), top.ncells)
    @inbounds for k in 1:nextra
        c = top.extra_tri[k]
        out[top.ncells + k] = blend(values[c[1]], values[c[2]], values[c[3]],
            top.extra_weight[k])
    end
    return out
end

"""
    vertex_positions(top, zs) -> Vector{Point2d} or Vector{Point3d}

The drawable vertex buffer: the topology's vertices at the heights `zs`.

On a [`PlanarTarget`](@ref) the height is the vertex's third coordinate, and
[`ZeroHeights`](@ref) leaves the two-coordinate buffer alone.  On a
[`GlobeTarget`](@ref) the direction and the height turn into one earth-centred
point together, which is the step a globe's build stops short of.

This is the only part of a surface that the heights touch, and it is linear in
the number of vertices — no adjacency, no centroids, no projection.
"""
vertex_positions(top::SurfaceTopology{Point2d, <:PlanarTarget}, ::ZeroHeights) =
    top.positions

vertex_positions(top::SurfaceTopology, zs::AbstractVector, ntasks::Int = 1) =
    vertex_positions!(Vector{Point3d}(undef, length(top.positions)), top, zs, ntasks)

"""
    vertex_positions!(out, top, zs, ntasks = 1) -> out

[`vertex_positions`](@ref) into a buffer that already exists.

A plot at new heights writes over its old vertex buffer rather than allocating
another, which at a few million cells is the difference between a redraw that
allocates a hundred megabytes and one that allocates none.  The height of each
vertex is read where it is needed instead of through a [`spread`](@ref) of its
own, so the pass allocates nothing at all.
"""
function vertex_positions!(out::Vector{Point3d},
        top::SurfaceTopology{Point2d, <:PlanarTarget}, zs::AbstractVector,
        ntasks::Int = 1)
    checkheights(zs, top)
    n = top.ncells
    inparallel(n, ntasks) do lo, hi
        @inbounds for i in lo:hi
            p = top.positions[i]
            out[i] = Point3d(p[1], p[2], zs[i])
        end
    end
    @inbounds for k in eachindex(top.extra_tri)
        c = top.extra_tri[k]
        p = top.positions[n + k]
        out[n + k] = Point3d(p[1], p[2],
            blend(zs[c[1]], zs[c[2]], zs[c[3]], top.extra_weight[k]))
    end
    return out
end

function vertex_positions!(out::Vector{Point3d},
        top::SurfaceTopology{Point3d, <:GlobeTarget}, zs::AbstractVector,
        ntasks::Int = 1)
    checkheights(zs, top)
    target = top.target
    n = top.ncells
    inparallel(n, ntasks) do lo, hi
        @inbounds for i in lo:hi
            out[i] = globe_vertex(target, top.positions[i], zs[i])
        end
    end
    @inbounds for k in eachindex(top.extra_tri)
        c = top.extra_tri[k]
        z = blend(zs[c[1]], zs[c[2]], zs[c[3]], top.extra_weight[k])
        out[n + k] = globe_vertex(target, top.positions[n + k], z)
    end
    return out
end

function checkheights(zs::AbstractVector, top::SurfaceTopology)
    length(zs) == top.ncells || throw(ArgumentError("zs has $(length(zs)) entries, \
        but there are $(top.ncells) cells"))
    return nothing
end

"""
    surface_topology(target::PlotTarget, cells;
                     ntasks = Threads.nthreads(),
                     connectivity = DiscreteGlobalGrids.Vertex()) -> SurfaceTopology

Build everything about the surface over `cells` that its heights do not touch.

`cells` is anything [`cellregion`](@ref) accepts: a set of cells at one level
whose adjacency the package can answer.  That is more than [`tessellate`](@ref)
needs, because a surface is built out of which cells touch which.

The result has one vertex per cell and the triangles of the grid's dual between
them.  A cell whose neighbours are absent takes part in fewer triangles, so a
partial grid comes out with a ragged edge rather than a wrong one.

`connectivity` is what counts as a neighbour.  `Vertex()`, the default, is the
one that gives a complete surface; `Edge()` would miss every corner of a
square-celled grid, whose four cells touch only diagonally.

On a [`PlanarTarget`](@ref) two cases need more than three centroids, both
counted in `nsplit`: a triangle straddling the cut is split against it, and the
one at each pole is drawn as the cap it covers.

This is the expensive half of a surface — adjacency, a centroid per cell, the
corner scan, and one bulk projection.  [`vertex_positions`](@ref) is the other
half, and the only one the heights reach.
"""
function surface_topology(target::PlotTarget, cr::CellRegion;
        ntasks::Int = Threads.nthreads(), connectivity = DGG.Vertex())
    P = pointtype(target)
    n = length(cr)
    n == 0 && return SurfaceTopology(target, cr, P[], GLTriangleFace[],
        NTuple{3, Int32}[], NTuple{3, Float32}[], 0, 0)

    positions = surface_vertices(target, cr, ntasks)
    adj = DGG.adjacency(cr.region; connectivity)

    nt = clamp(ntasks, 1, max(1, cld(n, 512)))
    chunks = Vector{SurfaceChunk{P}}(undef, nt)
    splits = zeros(Int, nt)
    Threads.@sync for t in 1:nt
        Threads.@spawn begin
            lo = 1 + div((t - 1) * n, nt)
            hi = div(t * n, nt)
            chunk = SurfaceChunk{P}(hi - lo + 1)
            chunks[t] = chunk
            splits[t] = fill_surface!(chunk, target, adj, lo, hi, n, positions, cr)
        end
    end

    return merge_surface(chunks, positions, n, sum(splits), target, cr)
end

surface_topology(target::PlotTarget, x; kwargs...) =
    surface_topology(target, cellregion(x); kwargs...)

"""
    SurfaceMesh(topology, positions)

A set of DGGS cells as one interpolated surface in an axis's data space: a
[`SurfaceTopology`](@ref) and the vertex buffer it comes to at some set of
heights.

`positions` are the drawable vertices — `Point2d` on a flat map, `Point3d` on a
globe or on a flat map carrying heights.  `faces`, `ncells` and `nsplit` read
through to the topology, and the arrays are shared with it rather than copied,
so a mesh at new heights costs one vertex buffer and nothing else.

Built by [`triangulate`](@ref).  Compare [`CellMesh`](@ref), which draws the
same cells as flat patches.
"""
struct SurfaceMesh{V, Top <: SurfaceTopology}
    topology::Top
    positions::Vector{V}
end

function Base.getproperty(m::SurfaceMesh, name::Symbol)
    name === :topology && return getfield(m, :topology)
    name === :positions && return getfield(m, :positions)
    return getproperty(getfield(m, :topology), name)
end

Base.propertynames(::SurfaceMesh) =
    (:topology, :positions, :faces, :ncells, :nsplit, :target, :cells,
        :extra_tri, :extra_weight)

Base.isempty(m::SurfaceMesh) = isempty(m.topology)

ntriangles(m::SurfaceMesh) = ntriangles(m.topology)

spread(values::AbstractVector, m::SurfaceMesh) = spread(values, m.topology)

"""
    triangulate(target::PlotTarget, cells, zs = ZeroHeights(length(cells));
                ntasks = Threads.nthreads(),
                connectivity = DiscreteGlobalGrids.Vertex()) -> SurfaceMesh

Build the interpolated surface over `cells` in `target`'s coordinate space, at
the heights `zs`.

The two halves of the work, [`surface_topology`](@ref) and
[`vertex_positions`](@ref), run back to back.  A plot keeps them apart, so that
changing the heights costs only the second; this is the one-shot form, for
everything that is not a live plot.

`zs` is one height per cell, and it lands in the geometry: on a
[`PlanarTarget`](@ref) as each vertex's third coordinate, in the axis's data
space, and on a [`GlobeTarget`](@ref) as a height above the ellipsoid, in the
axis's units, which lifts the vertex straight out from the centre.  The default,
[`ZeroHeights`](@ref), is the flat surface, and keeps a planar mesh's vertices
two-dimensional.
"""
function triangulate(target::PlotTarget, cr::CellRegion,
        zs::AbstractVector = ZeroHeights(length(cr)); kwargs...)
    length(zs) == length(cr) || throw(ArgumentError("zs has $(length(zs)) entries, \
        but there are $(length(cr)) cells"))
    top = surface_topology(target, cr; kwargs...)
    return SurfaceMesh(top, vertex_positions(top, zs))
end

triangulate(target::PlotTarget, x, zs::AbstractVector; kwargs...) =
    triangulate(target, cellregion(x), zs; kwargs...)

triangulate(target::PlotTarget, x; kwargs...) = triangulate(target, cellregion(x); kwargs...)

# The globe path ignores the planar-only arguments, so that the driver has one
# call site.
fill_surface!(chunk::SurfaceChunk, target::GlobeTarget, adj, lo, hi, n, positions, cr) =
    fill_surface!(chunk, target, adj, lo, hi, n)

function merge_surface(chunks::Vector{SurfaceChunk{P}}, positions::Vector{P},
        n::Int, nsplit::Int, target::PlotTarget, cr::CellRegion) where {P}
    nt = length(chunks)
    face_offset = zeros(Int, nt + 1)
    extra_offset = zeros(Int, nt + 1)
    for t in 1:nt
        chunk = chunks[t]
        face_offset[t + 1] = face_offset[t] + length(chunk.faces)
        extra_offset[t + 1] = extra_offset[t] + length(chunk.extra)
    end

    nextra = extra_offset[end]
    faces = Vector{GLTriangleFace}(undef, face_offset[end])
    extra_tri = Vector{NTuple{3, Int32}}(undef, nextra)
    extra_weight = Vector{NTuple{3, Float32}}(undef, nextra)
    nextra > 0 && resize!(positions, n + nextra)

    Threads.@sync for t in 1:nt
        Threads.@spawn begin
            chunk = chunks[t]
            foff = face_offset[t]
            # An index at or below `n` names a centroid and is already right;
            # above it is this task's own vertex.
            shift = Int32(extra_offset[t])
            cutoff = Int32(n)
            @inbounds for j in eachindex(chunk.faces)
                x, y, z = chunk.faces[j]
                faces[foff + j] = GLTriangleFace(
                    x > cutoff ? x + shift : x,
                    y > cutoff ? y + shift : y,
                    z > cutoff ? z + shift : z,
                )
            end
            eoff = extra_offset[t]
            copyto!(positions, n + eoff + 1, chunk.extra, 1, length(chunk.extra))
            copyto!(extra_tri, eoff + 1, chunk.extra_tri, 1, length(chunk.extra_tri))
            copyto!(extra_weight, eoff + 1, chunk.extra_weight, 1, length(chunk.extra_weight))
        end
    end

    project!(target, positions)
    return SurfaceTopology(target, cr, positions, faces, extra_tri, extra_weight,
        n, nsplit)
end
