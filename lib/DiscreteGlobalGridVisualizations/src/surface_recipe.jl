# # The surface recipe
#
# `dggsurface` is `dggpoly`'s counterpart: same cells, same colour vector, but
# drawn as the field the cells sample rather than as the cells themselves.  The
# mesh comes from `surface.jl`; everything here is Makie plumbing.

"""
    dggsurface(cells; color = ..., kwargs...)

Draw a set of DGGS cells as one interpolated surface.

Where [`dggpoly`](@ref) gives each cell its own patch of one flat colour, this
puts a vertex at each cell's **centroid**, gives it that cell's value, and joins
the centroids with the triangles of the grid's dual — so the value varies
continuously between cell centres instead of jumping at cell edges.  It is the
right picture for a sampled field (elevation, temperature, a model output) and
the wrong one for a categorical one, where the cell boundaries are the point.

`cells` is anything [`cellregion`](@ref) accepts: an `AbstractGrid`, a
`PartialGrid`, a `CellVector` or a `CellLookup`.  Unlike `dggpoly` it will not
take a bare list of ids or a `MultiOrderCellSet`, because a surface is built out
of which cells touch which and neither of those says.

`color` is a single colour, or one value or colour **per cell**.

```julia
sys = DGG.IGeo7System()
cells = DGG.CellVector(DGG.query(sys, DGG.MultiOrderCoverage(extent); level = 9))
dggsurface(cells; color = elevation)
```

The mesh has one vertex per cell and about two triangles per cell, against
`dggpoly`'s six vertices and four triangles, so it is also the cheaper of the
two to draw — and a recolour costs nothing at all, because a per-cell colour
vector *is* the vertex buffer.

A partial grid comes out with a ragged edge: a cell whose neighbours are missing
takes part in fewer triangles, and the surface simply stops where the data does.
It stops half a cell short of where `dggpoly` would draw, because a surface can
only reach as far as the outermost centroid.

The plot adapts to the axis it is placed in, as `dggpoly` does.  In a
`GeoMakie.GlobeAxis` the surface closes over the whole sphere with nothing to
cut.  In an `Axis` or a `GeoAxis` it is built in longitude/latitude: triangles
straddling the map's cut are split against it, and the one triangle over each
pole is drawn as the polar cap it covers, so a global surface has no wedge
missing at the top.
"""
@recipe DGGSurface (cells,) begin
    """
    Sets the colour of the surface.  Either a single colour, or a vector with
    one entry per cell — numbers to be mapped through the colormap, or colours —
    which is interpolated across the triangles between cell centres.
    """
    color = @inherit patchcolor
    "Controls whether lights affect the surface.  Off by default: the height of a DGGS surface is not a shape, it is a colour."
    shading = Makie.NoShading
    """
    Whether to split triangles that straddle the map's cut meridian, and to fill
    the polar caps.  Only planar targets have a cut; on a globe this attribute
    does nothing.
    """
    wrap = true
    "How many tasks build the mesh."
    ntasks = Threads.nthreads()
    cycle = [:color => :patchcolor]
    Makie.mixin_generic_plot_attributes()...
    Makie.mixin_colormap_attributes()...
end

Makie.convert_arguments(::Type{<:DGGSurface}, cr::CellRegion) = (cr,)
Makie.convert_arguments(::Type{<:DGGSurface}, x) = (cellregion(x),)

"""
    vertex_colors(mesh::SurfaceMesh, color)

Spread a per-cell colour over the mesh's vertices.

The first `ncells` vertices *are* the cells, in order, so a colour vector as long
as the cell set is already the vertex buffer and is handed on untouched — which
is the usual case, and costs nothing.  It is only a mesh carrying extra vertices
from a seam or a pole that needs the gather, and then only for those.
"""
function vertex_colors(mesh::SurfaceMesh, color)
    color isa AbstractVector || return color
    nvertices = length(mesh.positions)
    length(color) == nvertices && return color
    length(color) == mesh.ncells ||
        throw(ArgumentError("color has $(length(color)) entries, but there are \
            $(mesh.ncells) cells"))
    return @inbounds [color[Int(c)] for c in mesh.vertex_cell]
end

function Makie.plot!(plot::DGGSurface{<:Tuple{<:CellRegion}})
    Makie.map!(
        plot, [:cells, :transform_func, :wrap, :ntasks], [:surfacemesh]
    ) do cells, transform_func, wrap, ntasks
        target = plot_target(transform_func)
        wrap || (target = uncut(target))
        return (triangulate(target, cells; ntasks),)
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
