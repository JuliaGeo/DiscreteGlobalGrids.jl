# # Turning cells into one mesh
#
# This is the whole point of the package.  Makie's `poly` builds one polygon per
# cell, triangulates each one with a general-purpose algorithm, and concatenates
# the results on every update.  A DGGS cell is a convex ring of five to ten
# corners, so it can be fanned from its first corner without any triangulator at
# all, and every cell is independent of every other one — which makes the whole
# job a parallel `for` loop over disjoint output buffers.
#
# Vertices are *not* shared between neighbouring cells even though the grids
# share corners.  Sharing them would mean a vertex belonging to three cells, and
# a mesh vertex can carry only one colour; per-cell colour is the entire point of
# these plots.

"""
    CellMesh(positions, faces, vertex_cell, ring_start, ncells, ndropped)

A set of cells as a single triangle mesh in an axis's data space.

  * `positions` — one vertex per cell corner, in the axis's data coordinates.
  * `faces` — the triangle fan of each cell, over `positions`.
  * `vertex_cell` — for each vertex, the position of the cell it came from.
    This is what makes recolouring cheap: new per-cell values become per-vertex
    values with one gather, and no geometry is rebuilt.
  * `ring_start` — where each drawn ring begins in `positions`, with a trailing
    sentinel, so that `positions[ring_start[i]:ring_start[i+1]-1]` is one closed
    outline.  Strokes are drawn from this and nothing else.
  * `ncells` — how many cells went in, which is the length a colour vector must
    have.
  * `ndropped` — how many cells were skipped as undrawable (see
    [`tessellate`](@ref)).

A cell can contribute more than one ring: on a planar target a cell straddling
the map's cut is split into two pieces, both of which point back at the same
cell through `vertex_cell`.
"""
struct CellMesh{P}
    positions::Vector{P}
    faces::Vector{GLTriangleFace}
    vertex_cell::Vector{Int32}
    ring_start::Vector{Int32}
    ncells::Int
    ndropped::Int
end

Base.isempty(m::CellMesh) = isempty(m.faces)

nrings(m::CellMesh) = max(0, length(m.ring_start) - 1)

# ## Per-task output
#
# Each task writes into its own buffers with chunk-local vertex indices; the
# merge below shifts them into place.  Faces are held as plain integer triples
# because that shift is easier to express on integers than on
# `GeometryBasics`' offset-encoded face indices.

struct MeshChunk{P}
    positions::Vector{P}
    faces::Vector{NTuple{3, Int32}}
    vertex_cell::Vector{Int32}
    ring_start::Vector{Int32}
    ring::Vector{Point2d}     # scratch: the cell ring in lon/lat degrees
    piece::Vector{Point2d}    # scratch: one side of a cut cell
end

function MeshChunk{P}(nhint::Int) where {P}
    positions = P[]
    faces = NTuple{3, Int32}[]
    vertex_cell = Int32[]
    ring_start = Int32[]
    sizehint!(positions, 6 * nhint)
    sizehint!(faces, 4 * nhint)
    sizehint!(vertex_cell, 6 * nhint)
    sizehint!(ring_start, nhint)
    return MeshChunk{P}(positions, faces, vertex_cell, ring_start, Point2d[], Point2d[])
end

# Push a ring of `n` points, already in the target's build space, as one fan.
@inline function emit_fan!(chunk::MeshChunk, ring, n::Int, cell::Int32)
    n < 3 && return
    base = Int32(length(chunk.positions))
    push!(chunk.ring_start, base + Int32(1))
    @inbounds for k in 1:n
        push!(chunk.positions, ring[k])
        push!(chunk.vertex_cell, cell)
    end
    @inbounds for k in 2:(n - 1)
        push!(chunk.faces, (base + Int32(1), base + Int32(k), base + Int32(k + 1)))
    end
    return
end

# ## The globe path
#
# Nothing to cut, nothing to project: a unit-sphere corner goes straight to a
# vertex.

function fill_chunk!(chunk::MeshChunk, target::GlobeTarget, source, cells, lo::Int, hi::Int)
    ndropped = 0
    @inbounds for i in lo:hi
        ring = DGG.cell_boundary(source, cells[i])
        n = length(ring)
        if n < 3
            ndropped += 1
            continue
        end
        cell = Int32(i)
        base = Int32(length(chunk.positions))
        push!(chunk.ring_start, base + Int32(1))
        for k in 1:n
            push!(chunk.positions, globe_vertex(target, ring[k]))
            push!(chunk.vertex_cell, cell)
        end
        for k in 2:(n - 1)
            push!(chunk.faces, (base + Int32(1), base + Int32(k), base + Int32(k + 1)))
        end
    end
    return ndropped
end

# ## The planar path
#
# Corners become longitude/latitude, the ring is unwrapped so that it is
# contiguous in longitude rather than jumping a full turn, and a ring that ends
# up straddling the map's cut is split against it.  Projection happens later,
# over the finished buffer, so that it can be one bulk call.

"The signed difference `a - b` folded into `(-180, 180]`."
@inline lon_delta(a::Float64, b::Float64) = a - b - 360.0 * round((a - b) / 360.0)

"The unsigned circular distance between two longitudes."
@inline lon_distance(a::Float64, b::Float64) = abs(lon_delta(a, b))

"""
    ring_lonlat!(buf, ring, n) -> winding

Fill `buf` with the ring as longitude/latitude in degrees, unwrapped so that it
runs continuously, and return how far longitude turned in going once around.

Each corner is placed within half a turn of the one before it.  Summing those
steps around the closed ring gives `0` for an ordinary cell and `±360` for one
that *encircles* a pole — the sign says which, because the corners run
counter-clockwise seen from outside the sphere.  That total is the only reliable
way to tell a cell that crosses the map's cut from a cell that contains a pole,
and the two need completely different treatment.

The step from one corner to the next is a great-circle arc, and its longitude
sweep is the short way round — the arc's equatorial shadow stays inside the
wedge its endpoints span — so a step of less than a quarter turn is never
ambiguous.  A longer one can be, and rather than guess, [`trace_edge!`](@ref)
walks it.  Because that only happens near a pole, it is kept off the common
path: the loop below notes the longest step it took, and only if one was long
does the ring get traced again with care.
"""
function ring_lonlat!(buf::Vector{Point2d}, ring, n::Int)
    resize!(buf, n)
    previous = 0.0
    longest = 0.0
    @inbounds for k in 1:n
        p = ring[k]
        lon = atand(p[2], p[1])
        lat = asind(clamp(p[3], -1.0, 1.0))
        if k > 1
            step = lon_delta(lon, previous)
            longest = max(longest, abs(step))
            lon = previous + step
        end
        previous = lon
        buf[k] = Point2d(lon, lat)
    end
    # The step that closes the ring, back onto the first corner.
    @inbounds closing = lon_delta(buf[1][1], previous)
    longest = max(longest, abs(closing))
    longest >= 90.0 && return trace_ring!(buf, ring, n)
    @inbounds return (previous + closing) - buf[1][1]
end

"""
    trace_ring!(buf, ring, n) -> winding

[`ring_lonlat!`](@ref) for a ring with an edge long enough in longitude to be
ambiguous — which in practice means an edge running over or beside a pole.

Every such edge is walked by [`trace_edge!`](@ref), so `buf` can come out longer
than the ring: an edge that runs exactly across a pole gets the two pole corners
it passes through inserted, which is what turns a cell merely *touching* the
pole into a polygon that closes along the top of the map instead of one that
guesses a side and smears.
"""
function trace_ring!(buf::Vector{Point2d}, ring, n::Int)
    reference = reference_longitude(ring, n)
    empty!(buf)
    previous = atand(ring[1][2], ring[1][1])
    first_lon = previous
    @inbounds for k in 1:n
        a = ring[k]
        b = ring[k == n ? 1 : k + 1]
        push!(buf, Point2d(previous, asind(clamp(a[3], -1.0, 1.0))))
        previous = trace_edge!(buf, a, b, previous, reference, 24)
    end
    return previous - first_lon
end

"""
    reference_longitude(ring, n) -> Float64

The direction the ring sits in, as a longitude: the mean of its corners'
directions in the equatorial plane.

This is the tie-breaker for an edge that runs exactly over a pole, where both
ways round are the same arc and only the cell's own position says which sweep of
longitude belongs to it.
"""
function reference_longitude(ring, n::Int)
    x = 0.0
    y = 0.0
    @inbounds for k in 1:n
        p = ring[k]
        r = hypot(p[1], p[2])
        r > 1.0e-9 || continue
        x += p[1] / r
        y += p[2] / r
    end
    return atand(y, x)
end

"""
    trace_edge!(buf, a, b, previous, reference, depth) -> longitude of `b`

Follow one great-circle edge in longitude, pushing whatever points `buf` needs to
follow it faithfully, and return where `b` lands.

A step under a quarter turn is taken directly.  A longer one is halved at the
arc's midpoint and each half is followed in turn, which converges quickly because
only an edge passing close to a pole sweeps that far.  The exception is an edge
whose midpoint *is* a pole: there the two ways round are geometrically the same
arc, so the sweep is chosen to be the one that stays on the cell's side of the
sphere, and both pole corners are inserted so the drawn polygon closes along the
top or bottom of the map.
"""
function trace_edge!(buf::Vector{Point2d}, a, b, previous::Float64,
        reference::Float64, depth::Int)
    step = lon_delta(atand(b[2], b[1]), previous)
    (abs(step) < 90.0 || depth <= 0) && return previous + step

    mx, my, mz = a[1] + b[1], a[2] + b[2], a[3] + b[3]
    if hypot(mx, my) < 1.0e-9
        pole = mz >= 0 ? 90.0 : -90.0
        # Both sweeps end at `b`; take the one passing the cell's own side.
        forward = lon_distance(previous + 90.0, reference) <=
            lon_distance(previous - 90.0, reference)
        target = previous + (forward ? 180.0 : -180.0)
        push!(buf, Point2d(previous, pole))
        push!(buf, Point2d(target, pole))
        return target
    end

    scale = sqrt(mx * mx + my * my + mz * mz)
    mid = (mx / scale, my / scale, mz / scale)
    previous = trace_edge!(buf, a, mid, previous, reference, depth - 1)
    push!(buf, Point2d(previous, asind(clamp(mid[3], -1.0, 1.0))))
    return trace_edge!(buf, mid, b, previous, reference, depth - 1)
end

# Sutherland–Hodgman against one vertical half-plane in longitude, with `shift`
# added to the surviving longitudes so that the far side of a cut cell lands on
# the far edge of the map.  The clip is done on straight lon/lat segments, which
# is the same approximation the rest of Makie makes about a cell edge.
function clip_half!(out::Vector{Point2d}, ring::Vector{Point2d}, n::Int,
        cut::Float64, keep_below::Bool, shift::Float64)
    empty!(out)
    @inbounds for i in 1:n
        a = ring[i]
        b = ring[i == n ? 1 : i + 1]
        a_in = keep_below ? (a[1] <= cut) : (a[1] >= cut)
        b_in = keep_below ? (b[1] <= cut) : (b[1] >= cut)
        a_in && push!(out, Point2d(a[1] + shift, a[2]))
        if a_in != b_in
            t = (cut - a[1]) / (b[1] - a[1])
            push!(out, Point2d(cut + shift, a[2] + t * (b[2] - a[2])))
        end
    end
    return length(out)
end

@inline function shift_ring!(ring::Vector{Point2d}, n::Int, shift::Float64)
    shift == 0.0 && return ring
    @inbounds for k in 1:n
        ring[k] = Point2d(ring[k][1] + shift, ring[k][2])
    end
    return ring
end

"""
    emit_polar!(chunk, buf, n, cell, north, left)

Draw a cell that contains a pole.

Longitude is not a coordinate such a cell can be cut along: its ring winds a
full turn, so there is no meridian that separates an inside from an outside.
What is true instead is that the cell is everything between its ring and the
pole, and *that* is a shape longitude describes perfectly — a band spanning the
whole map, closed along the top or bottom edge.

The ring's corners already run monotonically in longitude (they went once around
the pole), so the band is drawn by sliding the ring's start to the map's left
edge, walking the corners across to the right edge, and fanning each step up to
the pole corner above the start.  One extra triangle closes the far top corner.
"""
function emit_polar!(chunk::MeshChunk, buf::Vector{Point2d}, n::Int, cell::Int32,
        north::Bool, left::Float64)
    piece = chunk.piece
    empty!(piece)
    # Corners increase in longitude around the north pole and decrease around the
    # south; reversing the southern ones lets one strip serve both.
    order = north ? (1:n) : (n:-1:1)
    @inbounds start = buf[first(order)][1]
    shift = left - start
    @inbounds for k in order
        push!(piece, Point2d(buf[k][1] + shift, buf[k][2]))
    end
    # The corner the ring closes onto, one full turn along.
    @inbounds push!(piece, Point2d(left + 360.0, piece[1][2]))
    pole = north ? 90.0 : -90.0
    push!(piece, Point2d(left, pole))          # index n + 2: the fan's apex
    push!(piece, Point2d(left + 360.0, pole))  # index n + 3

    base = Int32(length(chunk.positions))
    push!(chunk.ring_start, base + Int32(1))
    for p in piece
        push!(chunk.positions, p)
        push!(chunk.vertex_cell, cell)
    end
    apex = base + Int32(n + 2)
    @inbounds for k in 1:n
        push!(chunk.faces, (apex, base + Int32(k), base + Int32(k + 1)))
    end
    push!(chunk.faces, (apex, base + Int32(n + 1), base + Int32(n + 3)))
    return
end

function fill_chunk!(chunk::MeshChunk, target::PlanarTarget, source, cells, lo::Int, hi::Int)
    cut = target.cut
    left = cut - 360.0
    wrapping = needs_cutting(target)
    ndropped = 0
    @inbounds for i in lo:hi
        ring = DGG.cell_boundary(source, cells[i])
        n = length(ring)
        if n < 3
            ndropped += 1
            continue
        end
        cell = Int32(i)
        buf = chunk.ring
        winding = ring_lonlat!(buf, ring, n)
        # Tracing a pole-side edge can lengthen the ring, so the corner count
        # from here on is the buffer's, not the cell's.
        m = length(buf)

        # `wrap = false` means take the unwrapped ring as it is — no window to
        # move it into and no seam to split it on.
        if !wrapping
            emit_fan!(chunk, buf, m, cell)
            continue
        end

        if abs(winding) > 180.0
            emit_polar!(chunk, buf, m, cell, winding > 0, left)
            continue
        end

        lon_min = buf[1][1]
        lon_max = lon_min
        for k in 2:m
            l = buf[k][1]
            lon_min = min(lon_min, l)
            lon_max = max(lon_max, l)
        end

        # Move the ring bodily into the drawable window `[cut - 360, cut]`,
        # judged by its centre so that the ring stays whole.
        mid = 0.5 * (lon_min + lon_max)
        shift = (left + mod(mid - left, 360.0)) - mid
        shift_ring!(buf, m, shift)
        lon_min += shift
        lon_max += shift

        if lon_max > cut
            kept = clip_half!(chunk.piece, buf, m, cut, true, 0.0)
            emit_fan!(chunk, chunk.piece, kept, cell)
            kept = clip_half!(chunk.piece, buf, m, cut, false, -360.0)
            emit_fan!(chunk, chunk.piece, kept, cell)
        elseif lon_min < left
            kept = clip_half!(chunk.piece, buf, m, left, false, 0.0)
            emit_fan!(chunk, chunk.piece, kept, cell)
            kept = clip_half!(chunk.piece, buf, m, left, true, 360.0)
            emit_fan!(chunk, chunk.piece, kept, cell)
        else
            emit_fan!(chunk, buf, m, cell)
        end
    end
    return ndropped
end

# ## Driver

"""
    tessellate(target::PlotTarget, cells; ntasks = Threads.nthreads()) -> CellMesh

Build the mesh for `cells` in `target`'s coordinate space.

`cells` is anything [`cellset`](@ref) accepts.  Work is split into `ntasks`
contiguous chunks of cells, each of which fills its own buffers; the chunks are
then concatenated with their vertex indices shifted into place.

Cells that cannot be drawn are skipped and counted in the result's `ndropped`;
the only such cell is one whose boundary has fewer than three corners.

On a [`PlanarTarget`](@ref) a cell may end up drawn as more than one ring: two,
when it straddles the map's cut, and one that spans the full width of the map,
when it contains a pole.
"""
function tessellate(target::PlotTarget, cs::CellSet; ntasks::Int = Threads.nthreads())
    return _tessellate(target, cs.source, cs.cells, ntasks)
end

tessellate(target::PlotTarget, x; kwargs...) = tessellate(target, cellset(x); kwargs...)

function _tessellate(target::PlotTarget, source, cells, ntasks::Int)
    P = pointtype(target)
    n = length(cells)
    n == 0 && return CellMesh(P[], GLTriangleFace[], Int32[], Int32[Int32(1)], 0, 0)

    # Chunks below a few hundred cells cost more to spawn than to run.
    nt = clamp(ntasks, 1, max(1, cld(n, 512)))
    chunks = [MeshChunk{P}(cld(n, nt)) for _ in 1:nt]
    dropped = zeros(Int, nt)

    Threads.@sync for t in 1:nt
        lo = firstindex(cells) + div((t - 1) * n, nt)
        hi = firstindex(cells) + div(t * n, nt) - 1
        Threads.@spawn dropped[t] = fill_chunk!(chunks[t], target, source, cells, lo, hi)
    end

    return merge_chunks(chunks, n, sum(dropped), target)
end

function merge_chunks(chunks::Vector{MeshChunk{P}}, ncells::Int, ndropped::Int,
        target::PlotTarget) where {P}
    nt = length(chunks)
    vertex_offset = zeros(Int, nt + 1)
    face_offset = zeros(Int, nt + 1)
    ring_offset = zeros(Int, nt + 1)
    for t in 1:nt
        vertex_offset[t + 1] = vertex_offset[t] + length(chunks[t].positions)
        face_offset[t + 1] = face_offset[t] + length(chunks[t].faces)
        ring_offset[t + 1] = ring_offset[t] + length(chunks[t].ring_start)
    end

    positions = Vector{P}(undef, vertex_offset[end])
    vertex_cell = Vector{Int32}(undef, vertex_offset[end])
    faces = Vector{GLTriangleFace}(undef, face_offset[end])
    ring_start = Vector{Int32}(undef, ring_offset[end] + 1)
    ring_start[end] = Int32(vertex_offset[end] + 1)

    Threads.@sync for t in 1:nt
        Threads.@spawn begin
            chunk = chunks[t]
            voff = vertex_offset[t]
            foff = face_offset[t]
            roff = ring_offset[t]
            copyto!(positions, voff + 1, chunk.positions, 1, length(chunk.positions))
            copyto!(vertex_cell, voff + 1, chunk.vertex_cell, 1, length(chunk.vertex_cell))
            shift = Int32(voff)
            @inbounds for j in eachindex(chunk.faces)
                a, b, c = chunk.faces[j]
                faces[foff + j] = GLTriangleFace(a + shift, b + shift, c + shift)
            end
            @inbounds for j in eachindex(chunk.ring_start)
                ring_start[roff + j] = chunk.ring_start[j] + shift
            end
        end
    end

    project!(target, positions)
    return CellMesh(positions, faces, vertex_cell, ring_start, ncells, ndropped)
end
