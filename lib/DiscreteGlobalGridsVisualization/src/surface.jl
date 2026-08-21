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
    SurfaceMesh(positions, faces, vertex_cell, ncells, nsplit)

A set of DGGS cells as one interpolated surface in an axis's data space.

  * `positions` — the vertices.  The first `ncells` are the cells' centroids in
    cell order; the rest belong to triangles cut at the map's edge.
  * `faces` — the triangles of the grid's dual, over `positions`.
  * `vertex_cell` — the cell each vertex takes its value from; the identity over
    the first `ncells`.
  * `ncells` — the length a colour vector must have.
  * `nsplit` — triangles cut at the seam or capped at a pole.  Zero on a globe.

Built by [`triangulate`](@ref).  Compare [`CellMesh`](@ref), which draws the
same cells as flat patches.
"""
struct SurfaceMesh{P}
    positions::Vector{P}
    faces::Vector{GLTriangleFace}
    vertex_cell::Vector{Int32}
    ncells::Int
    nsplit::Int
end

Base.isempty(m::SurfaceMesh) = isempty(m.faces)

ntriangles(m::SurfaceMesh) = length(m.faces)

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

struct SurfaceChunk{P}
    faces::Vector{NTuple{3, Int32}}
    extra::Vector{P}          # vertices private to a cut or polar triangle
    extra_cell::Vector{Int32} # the cell each private vertex takes its value from
    poly::Vector{Point2d}     # scratch: a triangle being cut
    polytag::Vector{Int32}
    clip::Vector{Point2d}     # scratch: the surviving side of a cut
    cliptag::Vector{Int32}
    ring::Vector{Point2d}     # scratch: a triangle's outline, traced
    ringtag::Vector{Int32}
end

function SurfaceChunk{P}(nhint::Int) where {P}
    faces = NTuple{3, Int32}[]
    sizehint!(faces, 2 * nhint)
    return SurfaceChunk{P}(faces, P[], Int32[], Point2d[], Int32[], Point2d[], Int32[],
        Point2d[], Int32[])
end

@inline function push_face!(chunk::SurfaceChunk, a::Int32, b::Int32, c::Int32)
    push!(chunk.faces, (a, b, c))
    return
end

"""
    emit_private_fan!(chunk, pts, tags, m, n) -> ntriangles

Add `m` vertices of a cut or polar triangle as this task's own, and fan them.

`n` is the number of shared centroid vertices, where task-local numbering starts.
"""
function emit_private_fan!(chunk::SurfaceChunk{P}, pts::Vector{P}, tags::Vector{Int32},
        m::Int, n::Int) where {P}
    m < 3 && return 0
    base = Int32(n + length(chunk.extra))
    @inbounds for k in 1:m
        push!(chunk.extra, pts[k])
        push!(chunk.extra_cell, tags[k])
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
    outline_triangle!(pts, tags, u1, u2, u3, c1, c2, c3) -> winding

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
function outline_triangle!(pts::Vector{Point2d}, tags::Vector{Int32}, u1, u2, u3,
        c1::Int32, c2::Int32, c3::Int32)
    empty!(pts)
    empty!(tags)
    us = (u1, u2, u3)
    cs = (c1, c2, c3)
    reference = reference_longitude(us, 3)
    previous = atand(u1[2], u1[1])
    first_lon = previous
    @inbounds for k in 1:3
        j = k == 3 ? 1 : k + 1
        a = us[k]
        b = us[j]
        push!(pts, Point2d(previous, asind(clamp(a[3], -1.0, 1.0))))
        push!(tags, cs[k])

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
                push!(tags, cs[k])
                push!(pts, Point2d(previous + step, pole))
                push!(tags, cs[j])
            end
        end
        previous += step
    end
    return previous - first_lon
end

"""
    emit_polar_band!(chunk, m, n, north, left, cut) -> ntriangles

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
        left::Float64, cut::Float64)
    ring = chunk.ring
    ringtag = chunk.ringtag
    # Longitude increases around the north pole and decreases around the south.
    # Turning the southern outlines round lets one strip serve both.
    if !north
        reverse!(ring)
        reverse!(ringtag)
    end
    # The pole corners belong to whichever cell reaches nearest the pole.
    polecell = ringtag[1]
    highest = abs(ring[1][2])
    @inbounds for k in 2:m
        if abs(ring[k][2]) > highest
            highest = abs(ring[k][2])
            polecell = ringtag[k]
        end
    end
    pole = north ? 90.0 : -90.0

    # Put the first corner in the window, so that the two turns emitted below
    # are certain to cover the whole of it between them.
    @inbounds base_shift = (left + mod(ring[1][1] - left, 360.0)) - ring[1][1]

    pts = chunk.poly
    tags = chunk.polytag
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
            empty!(tags)
            # Wound counter-clockwise either way: at the north pole the cap lies
            # above its lower edge, at the south pole below it.
            if north
                push!(pts, Point2d(x0, p[2]));  push!(tags, ringtag[k])
                push!(pts, Point2d(x1, q[2]));  push!(tags, ringtag[j])
                push!(pts, Point2d(x1, pole));  push!(tags, polecell)
                push!(pts, Point2d(x0, pole));  push!(tags, polecell)
            else
                push!(pts, Point2d(x1, q[2]));  push!(tags, ringtag[j])
                push!(pts, Point2d(x0, p[2]));  push!(tags, ringtag[k])
                push!(pts, Point2d(x0, pole));  push!(tags, polecell)
                push!(pts, Point2d(x1, pole));  push!(tags, polecell)
            end
            ntri += emit_window_polygon!(chunk, 4, n, left, cut)
        end
    end
    return ntri
end

"""
    emit_traced_polygon!(chunk, m, n, left, cut) -> ntriangles

Draw an ordinary triangle whose outline, in `chunk.ring`, does not fit the map
in one piece.

The outline is emitted at three turns of longitude and each copy clipped to the
window, so a triangle straddling the cut comes out as its two pieces without
either having to be identified as the far side.  Three, because the outline can
run either way in longitude from the corner the turns are measured off.
"""
function emit_traced_polygon!(chunk::SurfaceChunk{Point2d}, m::Int, n::Int,
        left::Float64, cut::Float64)
    ring = chunk.ring
    ringtag = chunk.ringtag
    @inbounds base_shift = (left + mod(ring[1][1] - left, 360.0)) - ring[1][1]
    pts = chunk.poly
    tags = chunk.polytag
    ntri = 0
    for turn in (-360.0, 0.0, 360.0)
        empty!(pts)
        empty!(tags)
        @inbounds for k in 1:m
            push!(pts, Point2d(ring[k][1] + base_shift + turn, ring[k][2]))
            push!(tags, ringtag[k])
        end
        ntri += emit_window_polygon!(chunk, m, n, left, cut)
    end
    return ntri
end

"""
    emit_window_polygon!(chunk, m, n, left, cut) -> ntriangles

Clip the `m`-gon in `chunk.poly` to the map's window `[left, cut]` and fan what
survives.  The two clips ping-pong between the chunk's scratch buffers, leaving
the result in `chunk.poly`.
"""
function emit_window_polygon!(chunk::SurfaceChunk{Point2d}, m::Int, n::Int,
        left::Float64, cut::Float64)
    kept = clip_tagged!(chunk.clip, chunk.cliptag, chunk.poly, chunk.polytag,
        m, left, false, 0.0)
    kept < 3 && return 0
    kept = clip_tagged!(chunk.poly, chunk.polytag, chunk.clip, chunk.cliptag,
        kept, cut, true, 0.0)
    return emit_private_fan!(chunk, chunk.poly, chunk.polytag, kept, n)
end

"""
    clip_tagged!(out, outtag, pts, tags, m, cut, keep_below, shift) -> nkept

Sutherland–Hodgman against one vertical half-plane in longitude, carrying each
vertex's cell along with it, with `shift` added to the surviving longitudes so
that the far side of a cut triangle lands on the far edge of the map.

A crossing landing exactly on a corner is skipped: the corner is already there,
and repeating it would put a zero-area triangle in the fan.

A vertex created *on* the cut takes the cell of the nearer end of its edge.  This
is the one approximate value in the mesh, confined to a strip one cell wide along
the seam; carrying a blend of two cells instead would cost every vertex a lerp.
"""
function clip_tagged!(out::Vector{Point2d}, outtag::Vector{Int32}, pts::Vector{Point2d},
        tags::Vector{Int32}, m::Int, cut::Float64, keep_below::Bool, shift::Float64)
    empty!(out)
    empty!(outtag)
    @inbounds for i in 1:m
        j = i == m ? 1 : i + 1
        a = pts[i]; b = pts[j]
        ta = tags[i]; tb = tags[j]
        a_in = keep_below ? (a[1] <= cut) : (a[1] >= cut)
        b_in = keep_below ? (b[1] <= cut) : (b[1] >= cut)
        if a_in
            push!(out, Point2d(a[1] + shift, a[2]))
            push!(outtag, ta)
        end
        if a_in != b_in
            t = (cut - a[1]) / (b[1] - a[1])
            # `t` at either end means the crossing IS one of the corners, which
            # is pushed as a corner in its own right; adding it again would put
            # a zero-area triangle in the fan.
            if 0.0 < t < 1.0
                push!(out, Point2d(cut + shift, a[2] + t * (b[2] - a[2])))
                push!(outtag, t <= 0.5 ? ta : tb)
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
            winding = outline_triangle!(chunk.ring, chunk.ringtag,
                DGG.cell_centroid(source, cells[a]),
                DGG.cell_centroid(source, cells[b]),
                DGG.cell_centroid(source, cells[c]),
                Int32(a), Int32(b), Int32(c))
            m = length(chunk.ring)
            if abs(winding) > 180.0
                # Longitude increases around the north pole and decreases around
                # the south.
                emit_polar_band!(chunk, m, n, winding > 0, left, cut)
            else
                emit_traced_polygon!(chunk, m, n, left, cut)
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
    surface_vertex(target, p) -> Point2d or Point3d

A cell centroid `p`, a unit-sphere point, in the target's build space.

On a [`GlobeTarget`](@ref) that is the finished vertex.  On a
[`PlanarTarget`](@ref) it is longitude and latitude in degrees inside the map's
window `[cut - 360, cut]`, which [`triangulate`](@ref) projects in bulk at the
end.
"""
@inline surface_vertex(target::GlobeTarget, p) = globe_vertex(target, p)

@inline function surface_vertex(target::PlanarTarget, p)
    lon = atand(p[2], p[1])
    lat = asind(clamp(p[3], -1.0, 1.0))
    needs_cutting(target) || return Point2d(lon, lat)
    left = target.cut - 360.0
    return Point2d(left + mod(lon - left, 360.0), lat)
end

function surface_vertices(target::PlotTarget, cr::CellRegion)
    P = pointtype(target)
    n = length(cr)
    positions = Vector{P}(undef, n)
    cells = cr.cells
    source = cr.source
    Threads.@threads for p in 1:n
        @inbounds positions[p] = surface_vertex(target, DGG.cell_centroid(source, cells[p]))
    end
    return positions
end

"""
    triangulate(target::PlotTarget, cells; ntasks = Threads.nthreads(),
                connectivity = DiscreteGlobalGrids.Vertex()) -> SurfaceMesh

Build the interpolated surface over `cells` in `target`'s coordinate space.

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
"""
function triangulate(target::PlotTarget, cr::CellRegion;
        ntasks::Int = Threads.nthreads(), connectivity = DGG.Vertex())
    P = pointtype(target)
    n = length(cr)
    n == 0 && return SurfaceMesh(P[], GLTriangleFace[], Int32[], 0, 0)

    positions = surface_vertices(target, cr)
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

    return merge_surface(chunks, positions, n, sum(splits), target)
end

triangulate(target::PlotTarget, x; kwargs...) = triangulate(target, cellregion(x); kwargs...)

# The globe path ignores the planar-only arguments, so that the driver has one
# call site.
fill_surface!(chunk::SurfaceChunk, target::GlobeTarget, adj, lo, hi, n, positions, cr) =
    fill_surface!(chunk, target, adj, lo, hi, n)

function merge_surface(chunks::Vector{SurfaceChunk{P}}, positions::Vector{P},
        n::Int, nsplit::Int, target::PlotTarget) where {P}
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
    vertex_cell = Vector{Int32}(undef, n + nextra)
    @inbounds for p in 1:n
        vertex_cell[p] = Int32(p)
    end
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
            eoff = n + extra_offset[t]
            copyto!(positions, eoff + 1, chunk.extra, 1, length(chunk.extra))
            copyto!(vertex_cell, eoff + 1, chunk.extra_cell, 1, length(chunk.extra_cell))
        end
    end

    project!(target, positions)
    return SurfaceMesh(positions, faces, vertex_cell, n, nsplit)
end
