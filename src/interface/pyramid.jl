# What a pyramid is, to anything that reads one.
#
# Drawing every cell of a large set is the wrong shape of work, and the way out
# is always the same: name a COARSER level of the same hierarchy and read the
# data through it, so the cost of a frame follows the size of the screen rather
# than the size of the data. Two implementations do that and share nothing else.
#
#   * a pyramid over data already in memory, whose coarse levels are SYNTHESIZED
#     on demand — `DiscreteGlobalGridsVisualization.CellPyramid`
#   * a pyramid over a store, whose coarse levels were WRITTEN and are read like
#     any other level — the pyramid layout's `StorePyramid`
#
# The descent is the same either way, which is why it is worth one vocabulary:
# ask which levels there are, prune the branches that hold nothing, and read a
# value per drawn cell.
#
# Include order: after `interface/system.jl`, whose `levels` this extends.

"""
    AbstractPyramid

A set of values on a hierarchical grid, readable at several levels.

Four verbs, and no fields:

    system(pyr)                      -> AbstractHierarchicalGridSystem
    levels(pyr)                      -> AbstractUnitRange{Int}
    holdsdata(pyr, cell)             -> Bool
    cellvalues(pyr, level, cells)    -> AbstractVector

A descent starts at [`rootcells`](@ref), keeps the cells that are on screen and
whose subtrees [`holdsdata`](@ref DiscreteGlobalGrids.holdsdata) admits, refines
what survives, and stops when the cells are about a pixel across — at which
point [`cellvalues`](@ref DiscreteGlobalGrids.cellvalues) answers the frame.
"""
abstract type AbstractPyramid end

"""
    levels(pyr::AbstractPyramid) -> AbstractUnitRange{Int}

The levels `pyr` can be read at, coarsest first. `first` is where a descent
starts and `last` is the level the data itself lives at.

A pyramid's range need not be its system's: a store written from level-9 data
holds `0:9` however deep the system goes.
"""
function levels(::AbstractPyramid) end

"""
    holdsdata(pyr::AbstractPyramid, cell) -> Bool

Whether the subtree under `cell` can hold any of `pyr`'s data.

This is what keeps a descent proportional to what is on screen: two cells whose
answer is `false` are never refined, and the branches below them are never
visited. `cell` may be at any level of [`levels`](@ref DiscreteGlobalGrids.levels)`(pyr)`.

**`false` is a promise and `true` is not.** An implementation that answers
`true` where nothing is stored costs a wasted refinement; one that answers
`false` over data drops it from the picture. A pyramid whose coarse levels are
stored answers exactly; one that estimates from a bounding cap errs towards
`true`.
"""
function holdsdata end

"""
    cellvalues(pyr::AbstractPyramid, level::Integer, cells) -> AbstractVector

The value `pyr` shows at each of `cells`, which are at `level`.

As long as `cells`, in the same order, with `missing` wherever the pyramid holds
no value for that cell — a cell over the sea, over a hole in the coverage, or
inside a subtree nothing was ever written to.

This is the verb a frame is drawn from, so it is asked for a whole level's worth
of cells at once rather than one at a time: a stored pyramid answers it with one
read per chunk touched, which a per-cell verb could not.
"""
function cellvalues end
