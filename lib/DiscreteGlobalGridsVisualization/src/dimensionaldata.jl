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
