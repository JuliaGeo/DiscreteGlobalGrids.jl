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
`DiscreteGlobalGrids.cellposition`, in no space at all; a bare vector of ids has
to be indexed here, which is the one place a large set costs something up front.
Hand a `CellVector` rather than a `collect` of one when the set is large.
"""
function celllocator end

celllocator(cv::DGG.AbstractCellVector) = c -> something(DGG.cellposition(cv, c), 0)
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
that hold nothing), the cells themselves, and the locator that turns one of them
back into a position in the user's value vector.

Building one is `O(samples)`, not `O(ncells)`: nothing here walks the data.
"""
struct CellPyramid{S, G, C, L}
    system::S
    leafgrid::G
    leaflevel::Int
    rootlevel::Int
    capcentre::DGG.UnitSphericalPoint{Float64}
    capradius::Float64
    samplepoints::Vector{DGG.UnitSphericalPoint{Float64}}
    cells::C
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
        centre, radius, points, cs.cells, celllocator(cs.cells), n)
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

# # Every leaf under a drawn cell
#
# Nearest neighbour asks the level below one question — which leaf lies under
# this cell's centre.  A summary asks after all of them, and listing them would
# cost what drawing them costs, which is the thing a resampling exists not to
# do.  Naming them costs nothing, though, wherever a system keeps a subtree
# together in its canonical order: `DiscreteGlobalGrids.descendant_range` is
# then the interval of leaf-level *positions* the subtree occupies, and a set
# stored in that same order meets the interval in one run, found by binary
# search.  So a grouping is two searches per cell drawn, whatever the subtree
# under it holds, and only the reduction that follows reads the leaves.

"""
    LeafPositions(cells, grid)

`cells` as the positions they occupy in `grid`, read one at a time rather than
built.

Both containers a frame can be summarised over — a cell vector and a grid — hold
their cells in the level's own order, so this is an ascending vector of integers
and [`subtreeranges`](@ref) can search it without materialising it.
"""
struct LeafPositions{C, G} <: AbstractVector{Int}
    cells::C
    grid::G
end

Base.size(lp::LeafPositions) = (length(lp.cells),)

Base.@propagate_inbounds function Base.getindex(lp::LeafPositions, k::Int)
    @boundscheck checkbounds(lp, k)
    return something(DGG.cellposition(lp.grid, lp.cells[k]))
end

"""
    leafpositions(cells, grid) -> LeafPositions or nothing

`cells` as leaf-level positions, or `nothing` where they are stored in no order
worth searching.

A `CellVector` is strictly ascending by construction, and a grid's positions
ascend in canonical id order, so either can be met by an interval.  A bare
vector of ids promises nothing of the kind — it is a list, and the honest answer
for it is that there is no answer.
"""
leafpositions(cells::DGG.AbstractCellVector, grid) = LeafPositions(cells, grid)
leafpositions(cells::GridCells, grid) = LeafPositions(cells, grid)
leafpositions(::AbstractVector, grid) = nothing

"""
    subtreeranges(pyramid) -> f

Build `f(cell) -> UnitRange{Int32}`, the positions in the user's value vector of
every leaf under `cell`, empty where the cell has none.

Throws where the pyramid cannot answer it, rather than falling back on something
slower: a system whose subtrees are scattered through their level has no such
range at all, and a bare list of ids does not say where in itself one would
fall.  Both are conditions on a summary alone — nearest neighbour asks for one
leaf at a time and works over any backing.
"""
function subtreeranges(pyr::CellPyramid)
    sys = pyr.system
    DGG.has_sorted_subtrees(sys) || throw(ArgumentError(
        "`aggregate` needs the leaves under a drawn cell to sit together in the \
         level below, and $(nameof(typeof(sys))) does not order its cells that \
         way — `has_sorted_subtrees` is false, so a subtree is scattered \
         through its level.  Leave `aggregate` unset to resample nearest \
         neighbour."))
    positions = leafpositions(pyr.cells, pyr.leafgrid)
    positions === nothing && throw(ArgumentError(
        "`aggregate` needs the cells in the level's own order, so that the \
         leaves under a drawn cell are one run of the values, and a \
         $(nameof(typeof(pyr.cells))) is a list that promises no order.  Hand a \
         `CellVector`, a `CellLookup` or a grid, or leave `aggregate` unset to \
         resample nearest neighbour."))
    level = pyr.leaflevel
    return function (cell)
        r = DGG.descendant_range(sys, cell, level)
        lo = searchsortedfirst(positions, first(r))
        hi = searchsortedlast(positions, last(r))
        return (lo % Int32):(hi % Int32)
    end
end
