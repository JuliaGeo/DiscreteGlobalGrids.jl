# # Turning cells into one surface
#
# `tessellate.jl` draws the cells: every cell becomes its own little patch of a
# single colour, and neighbouring cells share nothing.  This file draws the
# *field* the cells sample.  A vertex is a cell **centroid**, it carries that
# cell's value, and the triangles between centroids let the GPU interpolate —
# which is what turns a discrete global grid into a continuous surface.
#
# The mesh is the grid's topological dual.  Where `n` cells meet at a grid
# corner, their `n` centroids form a small polygon around that corner, and
# fanning it gives `n - 2` triangles.  Do that at every corner and the triangles
# tile the sphere exactly once — no gap, no overlap.
#
# ## Finding the corners without a corner index
#
# `DiscreteGlobalGrids` has no "list the corners" verb, and cell boundaries are
# floating-point rings that neighbouring cells do not agree on bit for bit, so
# corners cannot be recovered by matching coordinates.  What the package does
# state exactly is adjacency: `neighbors` returns a cell's neighbours in one
# counter-clockwise turn.  That is enough, because the cells meeting at a corner
# are precisely a set of cells that all touch one another — a clique in the
# adjacency graph — and each such clique appears in each of its members' rings
# as a **run of consecutive neighbours**.
#
# So for a cell `a` and a consecutive pair `(b, c)` in its ring:
#
#  1. `b` and `c` must touch each other, or `a`, `b` and `c` do not share a
#     corner at all.  On a partial grid this is also what stops a triangle from
#     bridging a hole: two ring members either side of a missing one come out
#     consecutive, and this test throws the pair away.
#  2. Widen `(b, c)` to the largest run of neighbours of `a` that all touch one
#     another.  That run plus `a` is the corner's cell set.
#  3. Emit the triangle `(a, b, c)` only if `a` is the *smallest* member of that
#     set.
#
# Step 3 is the anti-duplication rule, and it costs nothing: no shared hash set,
# no post-hoc `unique`, nothing that would serialise the loop.  Each corner is
# owned by one of its cells, that cell fans it, and every other cell around it
# stays quiet.  Emitting per *pair* rather than per run matters too — in a
# degenerate neighbourhood two different maximal runs can overlap, and a pair
# belongs to at most one triangle whatever the runs do.
#
# Step 2 is what makes the same code right for square cells as for hexagons.
# Three cells to a corner (IGeo7, H3) gives runs of two and one triangle each;
# four (HEALPix, S2, ISEA4R) gives runs of three and a fan of two; five, where
# ISEA4R's diamonds meet an icosahedral vertex, gives three.  Without step 2 a
# four-cell corner would emit all four of its triangles instead of two, and the
# surface would be drawn twice over.
#
# The result is checked, not asserted: the test suite sums the signed spherical
# area of every triangle of a whole level and requires `4π`.

"""
    SurfaceMesh(positions, faces, vertex_cell, ncells, nsplit)

A set of DGGS cells as one interpolated surface in an axis's data space.

  * `positions` — the vertices.  The first `ncells` are the cells' own
    centroids, in cell order; anything after them belongs to a triangle that had
    to be cut at the map's edge (see below).
  * `faces` — the triangles of the grid's dual, over `positions`.
  * `vertex_cell` — the cell each vertex takes its value from.  For the first
    `ncells` this is the identity, which is what lets a per-cell colour vector be
    handed to the GPU untouched.
  * `ncells` — how many cells went in, which is the length a colour vector must
    have.
  * `nsplit` — how many triangles needed cutting at the map's seam or filling in
    at a pole.  Zero on a globe, and zero on a plot that does not reach either.

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

"How many triangles the surface is drawn from."
ntriangles(m::SurfaceMesh) = length(m.faces)

# ## Reading the adjacency graph
#
# `adjacency` hands back a CSR table of counter-clockwise rings addressed by
# in-region position, which is the same index the vertex buffer and the colour
# vector use.  Everything below is integers into that table.

"Is `b` in the ring `r`?"
@inline function inring(r, b::Int)
    @inbounds for x in r
        x == b && return true
    end
    return false
end

"""
    touches_all(adj, r, k, lo, hi, m) -> Bool

Does ring member `m` touch every member of the cyclic run `lo:hi`?

The run is already known to be pairwise touching, so one member against the run
is all that is left to check.  The scan runs from the far end inwards because
the near end is the one most likely to touch — in a hexagonal ring the next
neighbour round always does — and the answer wanted here is usually `false`.
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

The corner that the consecutive ring pair `(i, i + 1)` belongs to, as the widest
cyclic run of `r` whose members all touch one another.

`r` is one cell's counter-clockwise ring and `k` its length.  The run is grown
outwards from the pair, forwards first; a neighbourhood dense enough for two
maximal runs to overlap therefore gets one of them rather than an argument, and
the caller's per-pair emission keeps that from mattering.
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

The one rule that keeps a corner from being drawn once per cell around it.
"""
@inline function owns_corner(r, k::Int, lo::Int, hi::Int, a::Int)
    @inbounds for j in lo:hi
        a < r[mod1(j, k)] || return false
    end
    return true
end

# ## Per-task output
#
# A task owns a contiguous block of cells and writes triangles for the corners
# those cells own.  Face indices are into a buffer that is `[all centroids;
# this task's extra vertices]`: an index of `n` or less names a centroid and
# needs no fixing up at merge time, and a larger one is task-local and gets
# shifted.  Since the overwhelming majority of triangles are three centroids,
# the overwhelming majority of the merge is a `copyto!` of faces that are
# already correct.

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

`n` is the number of shared centroid vertices, which is where task-local vertex
numbering starts.
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

# ## The globe path
#
# Nothing to cut and nothing to fill: every triangle is three centroids, so the
# whole task writes faces and no vertices at all.

function fill_surface!(chunk::SurfaceChunk, ::GlobeTarget, adj, lo::Int, hi::Int, n::Int)
    @inbounds for a in lo:hi
        r = adj[a]
        k = length(r)
        k < 2 && continue
        # A ring of two has one pair, not two: `(r₁, r₂)` and `(r₂, r₁)` are the
        # same corner, and walking both would draw its triangle twice.  A cell
        # left with two neighbours is a partial grid's edge, so this is the
        # ragged-boundary case, not a curiosity.
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

# ## The planar path
#
# The centroids are already in longitude and latitude inside the map's window,
# so a triangle is shared-vertex work unless longitude says otherwise.  Two
# things say otherwise, and they are told apart by how far longitude turns in
# going once round the triangle: `0` for one that merely straddles the map's cut,
# and `±360` for the one triangle at each pole that *contains* it.

"""
    outline_triangle!(pts, tags, u1, u2, u3, c1, c2, c3) -> winding

Fill `pts` with the triangle's outline in longitude/latitude degrees and return
how far longitude turned in going once around.

The outline is the triangle's three corners and **nothing else** — no point is
interpolated along an edge — because an edge here is shared with the triangle on
its other side, which draws it as the straight lon/lat segment between the same
two corners.  Bending one copy of a shared edge and not the other is what opens
a gap between them, so the only thing this does that reading the corners would
not is work out how far longitude sweeps along each edge, and insert the two
pole corners when an edge runs exactly over a pole.

## Which way an edge sweeps

An edge is a great-circle arc between two neighbouring cell centroids, so it is
short, and its longitude sweep is the short way round — unless it runs beside a
pole, where a short arc can sweep almost half a turn either way and where the
short way is not necessarily the true one.

Longitude runs monotonically along any great circle that is not a meridian, and
which way it runs is the sign of `(a × b)ᶻ`.  That is the exact answer, and more
to the point it is an **anti-symmetric** one: the triangle on the other side of
this edge computes `(b × a)ᶻ`, which is this number negated bit for bit however
nearly zero it is.  So the two of them always sweep the shared edge in opposite
directions and always agree on the ground it covers, which is the property that
makes the drawn triangles meet rather than overlap or part.

The sign is only consulted for a sweep of a quarter turn or more.  A short
north-south edge has a nearly-zero cross product whose sign is noise, and it
does not need it: the short way is right.

## Edges over a pole

`(a × b)ᶻ` is exactly zero when the edge's great circle passes through both
poles, and if the corners are then half a turn apart in longitude the arc runs
over one of them.  Longitude jumps half a turn at the pole and the two pole
corners are inserted, which is what lets such a triangle close along the top of
the map instead of leaving a wedge of it undrawn.  Grids whose cells meet four
to a corner put a pole on a dual edge exactly, so this is the ordinary case for
them, not an exotic one.

Which half turn is decided by the triangle's own reference direction, and the
triangle on the other side of the edge, lying the opposite way, chooses the
other — so between them they cover the whole turn once.
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

Such a triangle winds a full turn in longitude, so no meridian separates its
inside from its outside and there is nothing to cut it along.  What is true
instead is that it is everything between its outline — the `m` points in
`chunk.ring`, which went once around the pole and so run monotonically in
longitude, increasing at the north pole and decreasing at the south — and the
pole itself: a band spanning the whole width of the map, closed along the top or
bottom edge.

The band is built as one trapezoid per outline segment, each running from that
segment up to the pole.  A trapezoid is convex whatever the latitudes do, which
a fan from a single corner is not: with corners at 88°, 88° and 87° the fan
turns one of its triangles inside out.

Each trapezoid is emitted twice, one turn of longitude apart, and both copies
are clipped to the map's window.  That is what keeps the band's lower edge lying
exactly on top of the edges its neighbouring triangles draw: the corners stay at
the longitudes they actually have, and the window decides which turn of each is
the one on the map, rather than the band being slid bodily to the left margin.
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

The outline is emitted at three turns of longitude and each copy is clipped to
the window, which is how a triangle straddling the cut comes out as the two
pieces it is — one against each edge of the map — without either of them having
to be identified as "the far side" first.  Three rather than two because the
outline runs both ways from the corner the turns are measured off, so a copy
either side of it can reach into the window.
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
survives.

The two clips ping-pong between the chunk's two scratch buffers, so the finished
polygon is back in `chunk.poly` and nothing was allocated.
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

A crossing that lands exactly on a corner is not added, because that corner is
already there and the repeat would only put a zero-area triangle in the fan.

A vertex created *on* the cut is given the cell of whichever end of the edge it
came out nearer to.  It is the one place in this file where a vertex's value is
approximate rather than exact, it is confined to a strip one cell wide along the
map's seam, and the alternative — carrying a blend of two cells' values through
to the colour gather — would make every vertex in the mesh pay for it.
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
        # A ring of two has one pair, not two: `(r₁, r₂)` and `(r₂, r₁)` are the
        # same corner, and walking both would draw its triangle twice.  A cell
        # left with two neighbours is a partial grid's edge, so this is the
        # ragged-boundary case, not a curiosity.
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
            # The ordinary triangle: short in longitude, and its three corners
            # already sit together in the window.  Three shared vertices, no
            # geometry of its own, and nothing to clip.
            if span < 90.0 && l2 == positions[b][1] && l3 == positions[c][1]
                push_face!(chunk, Int32(a), Int32(b), Int32(c))
                continue
            end

            # Everything else is drawn from its outline in longitude: it
            # either straddles the map's cut, or it lies over a pole, and the
            # winding is what says which.
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

The longitudes of the triangle's second and third corners made continuous with
its first, and the longest longitude step it took to get around.

`positions` holds each cell's centroid in the map's window, so a corner whose
continuous longitude is not the one it is stored at belongs to a triangle
straddling the cut.  A step of a quarter turn or more means an edge running
beside a pole, where longitude is not to be trusted and
[`outline_triangle!`](@ref) works the sweep out from the edge itself instead.
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

Where a cell centroid `p` — a unit-sphere point — starts life in the mesh.

On a [`GlobeTarget`](@ref) this is the finished vertex.  On a
[`PlanarTarget`](@ref) it is longitude and latitude in degrees, moved into the
map's drawable window `[cut - 360, cut]`, and the bulk projection at the end of
[`triangulate`](@ref) is what turns it into the axis's coordinates.
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

# ## Driver

"""
    triangulate(target::PlotTarget, cells; ntasks = Threads.nthreads(),
                connectivity = DiscreteGlobalGrids.Vertex()) -> SurfaceMesh

Build the interpolated surface over `cells` in `target`'s coordinate space.

`cells` is anything [`cellregion`](@ref) accepts — a grid, a `CellVector`, a
`PartialGrid`, a `CellLookup` — which is to say, a set of cells that all sit at
one level and whose adjacency the package can answer.  A cell's *neighbours* are
what a surface is built out of, which is why this asks for more than
[`tessellate`](@ref) does: a bare list of ids has no topology.

The result has one vertex per cell, carrying that cell's value, and the
triangles of the grid's dual between them.  Cells whose neighbours are absent
from the set simply take part in fewer triangles, so a partial grid comes out
with a ragged edge rather than a wrong one.

`connectivity` is what counts as a neighbour.  The default,
`DiscreteGlobalGrids.Vertex()`, is the one that produces a complete surface:
`Edge()` would miss the corners of a square-celled grid entirely, because the
four cells at such a corner touch only diagonally.

On a [`PlanarTarget`](@ref) two things need more than a triangle between three
centroids, and both are counted in the result's `nsplit`: a triangle straddling
the map's cut is split against it, and the single triangle at each pole is drawn
as the polar cap it covers.
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

# The globe path has no `positions` to consult, but taking it keeps one call
# site in the driver.
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
            # A face index at or below `n` names a centroid and is already
            # right; anything above it is this task's own vertex.
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
