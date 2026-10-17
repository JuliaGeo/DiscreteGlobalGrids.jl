# The regridding face: what a `to` or `from` of this package resolves to, and
# how a result over one of its grids is labelled. The grids themselves are the
# regridding sources and destinations — `GlobalRegridding` reads them through
# GridInterface's verbs — so this file holds only the spellings that are not
# grids yet and the two output hooks.

import GlobalRegridding as GR
# The target-resolution hook, extended here for this package's spellings.
import GlobalRegridding: _asgrid
import DimensionalData as DD

# Resolving `to` and `from`

"""
    GlobalRegridding._asgrid(target, name)
    GlobalRegridding._asgrid(target, name, src_grid)

Return the grid over the cells a regridding target names. A grid stands for
itself; a [`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup), a
[`CellVector`](@ref) and a [`MultiOrderCellSet`](@ref) name the
[`PartialGrid`](@ref) of their cells.

A bare system names no cells until a level is chosen. As a destination it takes
the level whose cells are closest in size to the source's, which is the only
spelling that reads `src_grid`; as a source there is nothing to match against
and it is an error.
"""
_asgrid(lk::AbstractCellLookup, name::AbstractString) = PartialGrid(lk)

_asgrid(cv::AbstractCellVector, name::AbstractString) = PartialGrid(cv)

_asgrid(set::MultiOrderCellSet, name::AbstractString) = PartialGrid(CellVector(set))

_asgrid(sys::AbstractHierarchicalGridSystem, name::AbstractString) =
    throw(ArgumentError(
        "`$name = $(typeof(sys).name.name)()` names no cells until a level is " *
        "chosen. As a destination the level is matched to the source's cell " *
        "areas, but as a source you must name it with `levelgrid(sys, l)`."))

_asgrid(sys::AbstractHierarchicalGridSystem, name::AbstractString, src_grid) =
    levelgrid(sys, levelfor(sys, src_grid))

# A `Cells` axis already names the cells a regrid would otherwise look for a
# raster lattice in, so a source given no `from` can point at the grid itself.
GR.dimsource(lk::AbstractCellLookup) = cellset(lk)

# Labelling the output

"""
    GlobalRegridding.destinationdims(grid::AbstractGrid, sampling)

Return the single [`Cells`](@ref DiscreteGlobalGrids.CellLookups.Cells) dimension a result over this grid carries: a
[`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup) over the destination's own cells, in its local index order.

A cell holds one value however that value was measured, so the lookup is the
same whichever `sampling` the method asks for. Being the grid's only axis, it
is the whole of the destination's shape, and a regrid needs no reshape to put a
result on it.
"""
GR.destinationdims(grid::AbstractGrid, ::DD.Lookups.Sampling) = (Cells(CellLookup(grid)),)

# `CR.Regridder(dst, src)` on a grid of this package reads it as the
# GridInterface grid it is (`ncells`, `cell_geometry`, `cell_area`) and sweeps
# `treeify(grid; leaf_capacity = ForSweeps())`, this package's own cursors,
# whose task weight is SpatialTreeInterface's `split_weight`: no forwarding.
