# # A pyramid over a cell set
#
# Drawing every cell of a large set is the wrong shape of work: a level-13 IGEO7
# tile is sixteen million hexagons, and a figure is a million pixels.  What a
# pyramid gives is a way of naming *a coarser level* of the same hierarchy and
# reading the data through it, so the cost of a frame follows the size of the
# screen rather than the size of the data.
#
# Two things are needed for that.  The first is a way down: from a cell at any
# level, its children, which every hierarchical system provides analytically.
# The second is a way back to the data: given a coarse cell, which entry of the
# user's value vector it should show.  Nearest neighbour makes that the leaf
# cell under the coarse cell's centre, so it is one `cellat` and one lookup, and
# — the point — its cost is per *drawn* cell, not per stored cell.

"""
    celllocator(cells) -> f

Build `f(cell) -> Int`, the position of `cell` in `cells`, or `0` when it is not
one of them.

This is what ties a resampled cell back to the value the user passed: `color[i]`
for `i = f(cell)`.  `CellVector`s and whole grids answer it themselves through
`DiscreteGlobalGrids.cellposition`; a bare vector of ids is indexed here.
"""
function celllocator end

celllocator(cv::DGG.CellVector) = c -> something(DGG.cellposition(cv, c), 0)
celllocator(gc::GridCells) = c -> something(DGG.cellposition(gc.grid, c), 0)

function celllocator(cells::AbstractVector)
    index = Dict{eltype(cells), Int}()
    sizehint!(index, length(cells))
    for (i, c) in enumerate(cells)
        get!(index, c, i)
    end
    return c -> get(index, convert(eltype(cells), c), 0)
end

"""
    CellPyramid(cellset; samples = 4096)

A cell set read as one level of a hierarchy, so that coarser levels of the same
hierarchy can stand in for it.

Holds the system, the level the data lives at, a spherical cap covering the data
(estimated from `samples` cells, so that the descent can prune whole branches
that hold nothing), and the locator that turns a cell back into a position in
the user's value vector.

Building one is `O(samples)`, not `O(ncells)`: nothing here walks the data.
"""
struct CellPyramid{S, G, L}
    system::S
    leafgrid::G
    leaflevel::Int
    rootlevel::Int
    capcentre::DGG.UnitSphericalPoint{Float64}
    capradius::Float64
    samplepoints::Vector{DGG.UnitSphericalPoint{Float64}}
    locate::L
    ncells::Int
end

# `source` may be a system already, or a grid drawn from one.
_system(sys::DGG.AbstractHierarchicalGridSystem) = sys
_system(x) = DGG.system(x)

function CellPyramid(cs::CellSet; samples::Integer = 4096)
    sys = _system(cs.source)
    sys isa DGG.AbstractHierarchicalGridSystem || throw(ArgumentError(
        "dggresample needs a hierarchical grid system to resample through, \
         and $(typeof(sys)) is not one"))
    n = length(cs.cells)
    n > 0 || throw(ArgumentError("cannot build a pyramid over an empty cell set"))

    # The leaf level is read off the ids: every cell index carries its own.
    step = max(1, cld(n, Int(samples)))
    picks = 1:step:n
    leaf = maximum(DGG.level(cs.cells[i]) for i in picks)

    points = [DGG.cell_centroid(cs.source, cs.cells[i]) for i in picks]
    centre, radius = _samplecap(sys, cs.cells[first(picks)], points)

    return CellPyramid(sys, DGG.levelgrid(sys, leaf), leaf, first(DGG.levels(sys)),
        centre, radius, points, celllocator(cs.cells), n)
end

"""
    _samplecap(system, cell, points) -> (centre, radius)

A spherical cap over the sampled cells, as a unit-sphere centre and an angular
radius in radians.

The cap is a *sample*, widened by one cell's own covering radius, so it is close
but not guaranteed to be a cover.  It is only used to skip branches of the
hierarchy that plainly hold nothing, and a slightly small cap costs a sliver at
the edge of a very ragged set rather than anything structural.
"""
function _samplecap(sys, cell, points)
    acc = DGG.UnitSphericalPoint(0.0, 0.0, 0.0)
    for p in points
        acc = acc + p
    end
    n = sqrt(sum(x -> x^2, acc))
    # Cells spread over a whole hemisphere or more sum to nothing in particular;
    # a cap over everything is the honest answer there.
    n > 1e-9 || return DGG.UnitSphericalPoint(0.0, 0.0, 1.0), Float64(pi)

    centre = DGG.UnitSphericalPoint(acc[1] / n, acc[2] / n, acc[3] / n)
    radius = 0.0
    for p in points
        radius = max(radius, _angle(centre, p))
    end
    # One cell's worth of slack, so a sampled centroid stands in for its cell.
    radius += DGG.node_extent(sys, cell).radius
    return centre, min(radius, Float64(pi))
end

@inline _angle(a, b) = acos(clamp(a[1] * b[1] + a[2] * b[2] + a[3] * b[3], -1.0, 1.0))

"""
    holdsdata(pyramid, cell) -> Bool

Whether the subtree under `cell` can hold any of the pyramid's cells.

`node_extent` covers every descendant of `cell` at every depth, so two caps that
miss each other rule the whole subtree out — which is what keeps the descent's
cost proportional to what is on screen rather than to the size of the world.
"""
function holdsdata(pyr::CellPyramid, cell)
    cap = DGG.node_extent(pyr.system, cell)
    return _angle(pyr.capcentre, cap.point) <= cap.radius + pyr.capradius
end

"""
    nearest(pyramid, cell) -> Int

The position in the value vector of the leaf cell under `cell`'s centre, or `0`
when there is none — which is how a resampled cell over a hole, or over the sea,
is recognised and dropped.

This is the nearest-neighbour resampling rule, and the reason a frame costs what
it does: one point location per *drawn* cell, whatever the level below holds.
"""
function nearest(pyr::CellPyramid, cell)
    centre = DGG.cell_centroid(pyr.system, cell)
    return pyr.locate(DGG.cellat(pyr.leafgrid, centre))
end
