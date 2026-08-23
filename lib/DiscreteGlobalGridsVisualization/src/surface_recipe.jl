# # The surface recipe
#
# `dggsurface` takes the same cells and the same colour vector `dggpoly` does,
# and draws the field they sample.  The mesh comes from `surface.jl`.

"""
    dggsurface(cells; color = ..., kwargs...)
    dggsurface(cells, zs; kwargs...)
    dggsurface(cell_dimarray; kwargs...)

Draw a set of DGGS cells as one interpolated surface.

Where [`dggpoly`](@ref) gives each cell its own patch of one flat colour, this
puts a vertex at each cell's **centroid**, gives it that cell's value, and joins
the centroids with the triangles of the grid's dual, so the value varies
continuously between cell centres instead of jumping at cell edges.  That is the
picture for a sampled field — elevation, temperature, a model output — and not
for a categorical one, where the cell boundaries are the point.

`cells` is anything [`cellregion`](@ref) accepts: an `AbstractGrid`, a
`PartialGrid`, a `CellVector` or a `CellLookup`.  Unlike `dggpoly` it takes
neither a bare list of ids nor a `MultiOrderCellSet`, neither of which says
which cells touch which.

`color` is a single colour, or one value or colour **per cell**.

```julia
sys = DGG.IGeo7System()
cells = DGG.CellVector(DGG.query(sys, DGG.MultiOrderCoverage(extent); level = 9))
dggsurface(cells; color = elevation)
```

## Raising the surface

`zs`, one height per cell, goes into the geometry the way it does in Makie's
`surface`: the vertex over a cell is that cell's centroid lifted to its height.
In a flat axis the height is the vertex's third coordinate, in the axis's data
space; on a `GeoMakie.GlobeAxis` it is a height above the ellipsoid, in the
axis's units, and lifts the vertex straight out from the centre.

```julia
dggsurface(cells, elevation; color = elevation)
```

Left out, the heights are [`ZeroHeights`](@ref) — the flat surface, at no cost
in memory and none in the vertex buffer, which stays two-dimensional.

A one-dimensional `DimArray` over a cell dimension gives both at once: its
lookup is the cells and its values are the heights, so `dggsurface(A)` is
`dggsurface(cells, values)`.  Colour is a separate matter either way; pass
`color = A` to colour the surface by the same field it is raised by.

A partial grid comes out with a ragged edge: a cell whose neighbours are missing
takes part in fewer triangles, and the surface stops where the data does — half
a cell short of `dggpoly`, since it can reach no further than the outermost
centroid.

The plot adapts to the axis it is placed in, as `dggpoly` does.  A
`GeoMakie.GlobeAxis` has nothing to cut.  An `Axis` or `GeoAxis` is built in
longitude/latitude, with triangles straddling the map's cut split against it and
the one over each pole drawn as the cap it covers.
"""
@recipe DGGSurface (cells, zs) begin
    """
    Sets the colour of the surface.  Either a single colour, or a vector with one
    entry per cell — numbers to be mapped through the colormap, or colours —
    interpolated across the triangles between cell centres.
    """
    color = @inherit patchcolor
    "Controls whether lights affect the surface.  Off by default: a DGGS surface carries a value, not a shape."
    shading = Makie.NoShading
    """
    Whether to split triangles that straddle the map's cut meridian and to fill
    the polar caps.  Only planar targets have a cut; on a globe this does
    nothing.
    """
    wrap = true
    "How many tasks build the mesh."
    ntasks = Threads.nthreads()
    cycle = [:color => :patchcolor]
    Makie.mixin_generic_plot_attributes()...
    Makie.mixin_colormap_attributes()...
end

"""
    cellvalues(x)

The plain vector inside whatever names one value per cell.

This is where a cube axis is unwrapped.  By the time one reaches a height or a
colour its dimensions have already been read — in the `zs` slot they are the
cells the surface is being built over — and only its values travel on to the
mesh.  They have to: a `DimArray` is a `DimArray` still after Makie's colour
conversion, and the backends will not draw that.
"""
cellvalues(x) = x

Makie.convert_arguments(::Type{<:DGGSurface}, cr::CellRegion) =
    (cr, ZeroHeights(length(cr)))
Makie.convert_arguments(P::Type{<:DGGSurface}, x) =
    Makie.convert_arguments(P, cellregion(x))
Makie.convert_arguments(::Type{<:DGGSurface}, cr::CellRegion, zs::AbstractVector) =
    (cr, cellvalues(zs))
Makie.convert_arguments(P::Type{<:DGGSurface}, x, zs::AbstractVector) =
    Makie.convert_arguments(P, cellregion(x), zs)

"""
    vertex_colors(mesh::SurfaceMesh, color)

Spread a per-cell colour over the mesh's vertices.

The first `ncells` vertices are the cells in order, so a vector already as long
as the vertex buffer is handed on untouched.  Only a mesh carrying extra
vertices from a seam or a pole needs the gather.  A cube axis is read for its
values first; see [`cellvalues`](@ref).
"""
function vertex_colors(mesh::SurfaceMesh, color)
    values = cellvalues(color)
    values isa AbstractVector || return values
    nvertices = length(mesh.positions)
    length(values) == nvertices && return values
    length(values) == mesh.ncells ||
        throw(ArgumentError("color has $(length(values)) entries, but there are \
            $(mesh.ncells) cells"))
    return @inbounds [values[Int(c)] for c in mesh.vertex_cell]
end

function Makie.plot!(plot::DGGSurface{<:Tuple{<:CellRegion, <:AbstractVector}})
    Makie.map!(
        plot, [:cells, :zs, :transform_func, :wrap, :ntasks], [:surfacemesh]
    ) do cells, zs, transform_func, wrap, ntasks
        target = plot_target(transform_func)
        wrap || (target = uncut(target))
        return (triangulate(target, cells, zs; ntasks),)
    end

    Makie.map!(m -> (m.positions,), plot, [:surfacemesh], [:mesh_positions])
    Makie.map!(m -> (m.faces,), plot, [:surfacemesh], [:mesh_faces])
    Makie.map!((m, c) -> (vertex_colors(m, c),), plot, [:surfacemesh, :color], [:mesh_color])

    Makie.mesh!(
        plot, plot.mesh_positions, plot.mesh_faces;
        color = plot.mesh_color,
        colormap = plot.colormap,
        colorscale = plot.colorscale,
        colorrange = plot.colorrange,
        lowclip = plot.lowclip,
        highclip = plot.highclip,
        nan_color = plot.nan_color,
        alpha = plot.alpha,
        shading = plot.shading,
        visible = plot.visible,
        transparency = plot.transparency,
        inspectable = plot.inspectable,
        space = plot.space,
        depth_shift = plot.depth_shift,
        clip_planes = plot.clip_planes,
        fxaa = plot.fxaa,
        transformation = :inherit_model,
    )
    return plot
end
