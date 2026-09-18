# Streaming a pyramid store: given a view, which level is worth drawing and
# which cells of it.
#
#   julia --project=docs scripts/pyramid_stream_poc.jl
#
# A prototype, in the terminal, of what an interactive plot would do on every
# camera change. The descent is the layout's own index used the way it was
# meant to be: start at the twelve root cells, keep the ones that are in view
# and that the store says hold something, refine, and stop when a cell is about
# a pixel across.
#
# That last rule has to be about PIXELS and not about a cell budget. A budget is
# wrong in both directions: over a sparsely covered earth it keeps descending
# until the count catches up -- drawing fifty thousand sub-pixel cells for a tile
# one pixel wide -- and over a densely covered one it stops early and blurs the
# frame.
#
# One `cellvalues` call per level, which is one read per CHUNK touched -- not
# per cell, and never a probe for a chunk that does not exist. The trace prints
# that count so a zoom that loads more than it needs is visible.
#
# Nothing here is part of the package yet. It is the shape the Makie recipe
# should have, written where it can be run and read first.

module PyramidStream

import DiscreteGlobalGrids as DGG
import GeometryOps as GO
import Extents
using Printf

export streamframe, FrameTrace

# ===========================================================================
# A view, as the descent sees it
# ===========================================================================

"""
    viewcap(extent) -> (centre, radius)

A spherical cap covering a longitude/latitude box, as a unit-sphere centre and
an angular radius in radians.

Conservative on purpose. The descent prunes with a cell's
[`node_extent`](@ref DiscreteGlobalGrids.node_extent), which covers its whole
subtree, and two caps that miss each other rule that subtree out; a cap over the
box rather than the box itself keeps the test to one dot product and errs
towards keeping.
"""
function viewcap(e::Extents.Extent)
    xs, ys = e.X, e.Y
    tosphere = GO.UnitSpherical.UnitSphereFromGeographic()
    corners = [DGG.UnitSphericalPoint(tosphere((x, y)))
               for x in (xs[1], (xs[1] + xs[2]) / 2, xs[2]),
                   y in (ys[1], (ys[1] + ys[2]) / 2, ys[2])]
    acc = reduce(+, corners)
    n = sqrt(sum(x -> x^2, acc))
    n > 1e-9 || return DGG.UnitSphericalPoint(0.0, 0.0, 1.0), Float64(pi)
    centre = DGG.UnitSphericalPoint(acc[1] / n, acc[2] / n, acc[3] / n)
    return centre, maximum(p -> angle(centre, p), corners)
end

@inline angle(a, b) = acos(clamp(a[1] * b[1] + a[2] * b[2] + a[3] * b[3], -1.0, 1.0))

"""
    metresperpixel(extent, pixels) -> Float64

How wide one pixel is on the ground, for a viewport `pixels` across showing
`extent`. The wider of the two spans decides it, so a window that is tall and
narrow is not drawn seven times finer than it can show.
"""
function metresperpixel(e::Extents.Extent, pixels::Integer)
    lat = (e.Y[1] + e.Y[2]) / 2
    wide = (e.X[2] - e.X[1]) * 111_320 * cosd(clamp(lat, -89.9, 89.9))
    tall = (e.Y[2] - e.Y[1]) * 110_540
    return max(wide, tall) / max(pixels, 1)
end

"Whether `cell`'s whole subtree could reach into the view cap."
function inview(sys, cell, centre, radius)
    cap = DGG.node_extent(sys, cell)
    return angle(centre, cap.point) <= cap.radius + radius
end

# ===========================================================================
# The descent
# ===========================================================================

"""
    FrameTrace

What one level of a descent cost: the candidates it considered, how many
survived the view and the store, and how many CHUNKS were read to find out.

The last number is the one to watch. It is the store reads the frame paid for,
and a descent that is working keeps it in the low tens however far in the camera
is.
"""
struct FrameTrace
    level::Int
    candidates::Int
    inview::Int
    present::Int
    chunks::Int
end

Base.show(io::IO, t::FrameTrace) = @printf(io,
    "  level %2d: %8d children -> %7d in view -> %7d with data   [%3d chunk%s read]",
    t.level, t.candidates, t.inview, t.present, t.chunks, t.chunks == 1 ? "" : "s")

"""
    streamframe(pyr, extent; pixels = 900, cellpixels = 3, maxcells = 200_000,
                verbose = true) -> (cells, level, values, trace)

The level of `pyr` worth drawing over `extent`, and the cells of it that are in
view and hold data.

`pyr` is a single-variable [`StorePyramid`](@ref DiscreteGlobalGrids.StorePyramid).
The descent stops when a cell is about `cellpixels` across on a viewport
`pixels` wide — that, and not the size of the data, is what decides the level, so
zooming in by a factor of seven moves it down one. `maxcells` is a safety net
underneath, not the rule.

Stopping on the cell BUDGET instead is the obvious mistake and it is wrong in
both directions: over a sparsely covered earth it descends until the count
catches up, drawing a hundred thousand sub-pixel cells for a tile a pixel wide,
and over a dense one it stops early and draws a blurred frame.

Returns the cells, the level they are at, their values, and one
[`FrameTrace`](@ref) per level visited.
"""
function streamframe(pyr, extent::Extents.Extent; pixels::Integer=900,
    cellpixels::Real=3, maxcells::Integer=200_000, verbose::Bool=true)

    sys = DGG.system(pyr)
    target = cellpixels * metresperpixel(extent, pixels)
    layout = pyr.layout
    centre, radius = viewcap(extent)
    trace = FrameTrace[]

    level = first(DGG.levels(pyr))
    roots = collect(DGG.rootcells(sys))
    seen = filter(c -> inview(sys, c, centre, radius), roots)
    cells, values = _keep(pyr, layout, level, seen, length(roots), trace, verbose)

    while level < last(DGG.levels(pyr)) && !isempty(cells)
        # Cells already at or below the target size: refining would draw seven
        # times as many of them into the same pixels.
        DGG.cellsize(sys, level) <= target && break
        born = 0
        next = eltype(cells)[]
        for c in cells, kid in DGG.children(sys, c)
            born += 1
            inview(sys, kid, centre, radius) && push!(next, kid)
        end
        # Stop BEFORE reading a level the screen cannot use. This is the whole
        # budget: nothing below here is fetched, so a zoomed-out view never
        # touches a deep chunk.
        (isempty(next) || length(next) > maxcells) && break
        newcells, newvalues = _keep(pyr, layout, level + 1, next, born, trace, verbose)
        isempty(newcells) && break
        level, cells, values = level + 1, newcells, newvalues
    end
    return cells, level, values, trace
end

# One level: which candidates the store has a value for, in one `cellvalues`
# call. The store answers `missing` where it holds nothing, so presence and
# value come out of the same read rather than out of a probe and then a read.
function _keep(pyr, layout, level, candidates, born, trace, verbose)
    chunks = _chunkstouched(layout, level, candidates)
    got = DGG.cellvalues(pyr, level, candidates)
    keep = findall(!ismissing, got)
    t = FrameTrace(level, born, length(candidates), length(keep), chunks)
    push!(trace, t)
    verbose && println(t)
    return candidates[keep], [got[i] for i in keep]
end

# How many store reads `cellvalues` is about to do: it groups its cells by
# chunk, and a chunk is the subtree of one cell of the chunk root level.
function _chunkstouched(layout, level, cells)
    sys = DGG.system(layout)
    rl = DGG.chunkrootlevel(layout, level)
    seen = Set{eltype(cells)}()
    for c in cells
        push!(seen, DGG.ancestor(sys, c, rl))
    end
    return length(seen)
end

# ===========================================================================
# The demo
# ===========================================================================

"""
    zoomtrace(pyr, centre, spans; maxcells)

Run [`streamframe`](@ref) over a sequence of ever-smaller boxes around `centre`,
which is what a camera zoom looks like to the store.
"""
function zoomtrace(pyr, centre::Tuple{Real,Real}, spans; pixels::Integer=900,
    maxcells::Integer=200_000)
    for span in spans
        e = Extents.Extent(X=(centre[1] - span / 2, centre[1] + span / 2),
            Y=(centre[2] - span / 2, centre[2] + span / 2))
        @printf("\nview %.4f deg across at (%.3f, %.3f), %.0f m per pixel\n",
            span, centre..., metresperpixel(e, pixels))
        t0 = time()
        cells, level, values, trace = streamframe(pyr, e; pixels=pixels, maxcells=maxcells)
        @printf("  => level %d, %d cells, %d chunk reads total, %.3f s\n",
            level, length(cells), sum(t.chunks for t in trace), time() - t0)
    end
    return nothing
end

end # module PyramidStream
