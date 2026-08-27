# # Plotting a cube axis
#
# A one-dimensional `DimArray` over a cell dimension is a set of cells and a
# value for each of them in one object — which is exactly the pair `dggsurface`
# draws.  Reading it apart is therefore a `convert_arguments` and nothing more:
# the dimension's lookup is the cells, and the array's values are the heights.

"""
    cellregion(A::DimensionalData.AbstractDimArray)

The cells of a one-dimensional cube axis.

The array's single dimension carries a `CellLookup`, whose cells these are.  A
dimension of any other kind names no cells, and falls to the same error a bare
list of ids does.
"""
cellregion(A::DD.AbstractDimArray{<:Any, 1}) = cellregion(DD.lookup(A, 1))

"""
    cellset(A::DimensionalData.AbstractDimArray)

The cells of a one-dimensional cube axis, for the recipes that need no adjacency.

`dggpoly(A)` and `dggresample(A)` draw the cells `A` is indexed by.  Neither
reads `A`'s values — only [`dggsurface`](@ref) does, as heights — so pass
`color = A` to see them.
"""
cellset(A::DD.AbstractDimArray{<:Any, 1}) = cellset(DD.lookup(A, 1))

# See [`cellvalues`](@ref): the dimensions have been read by the time a cube
# axis reaches a height or a colour, and what goes on to the mesh is its values.
cellvalues(A::DD.AbstractDimArray{<:Any, 1}) = parent(A)

"""
    convert_arguments(::Type{<:DGGSurface}, A::DimensionalData.AbstractDimArray)

Expand a one-dimensional cube axis into the `(cells, zs)` a surface is built
from: the cells its dimension names, and its own values as their heights.

Colour is not implied — `dggsurface(A)` draws `A` as relief in one flat colour,
and `dggsurface(A; color = A)` colours it by the same field.
"""
Makie.convert_arguments(P::Type{<:DGGSurface}, A::DD.AbstractDimArray{<:Any, 1}) =
    Makie.convert_arguments(P, cellregion(A), parent(A))

"""
    convert_arguments(::Type{<:DGGResample}, A::DimensionalData.AbstractDimArray, zs)

A cube axis in the cells slot of a resampled plot names its cells, and the
vector after it is heights.

Without it a cube axis followed by a vector would fall to the `(source, ids)`
form, which is for a system paired with the ids of its cells and not for an
array that already knows both.
"""
Makie.convert_arguments(::Type{<:DGGResample}, A::DD.AbstractDimArray{<:Any, 1},
    zs::AbstractVector) = (cellset(A), cellvalues(zs))
