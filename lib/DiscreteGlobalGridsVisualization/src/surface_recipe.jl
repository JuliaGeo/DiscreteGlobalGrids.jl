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
`color = A` to colour the surface by the same field it is raised by, or
`color = nothing`, which colours a surface by its own heights the way `surface`
does.

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
    interpolated across the triangles between cell centres.  `nothing` colours
    the surface by its own heights, the way Makie's `surface` does.
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
    "How many tasks build the mesh.  `automatic` is one per thread."
    ntasks = Makie.automatic
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
    reusable(cached, T, n, notthis) -> Vector{T} or nothing

The buffer a node computed last time, if it can be written over again.

It can when it is a `Vector{T}` of the length wanted and is not `notthis` — an
array the node does not own, which it would otherwise corrupt for whoever does.
"""
function reusable(cached, ::Type{T}, n::Int, notthis) where {T}
    cached === nothing && return nothing
    buffer = cached[1]
    buffer isa Vector{T} || return nothing
    (length(buffer) == n && buffer !== notthis) || return nothing
    return buffer
end

"""
    vertex_colors(top::SurfaceTopology, color)

Spread a per-cell colour over the surface's vertices.

A single colour passes through.  A vector one entry per cell is
[`spread`](@ref) over the vertices — which is a no-op unless the seam or a pole
added vertices of their own, and blends the triangle's three cells for each one
that did.  A cube axis is read for its values first; see [`cellvalues`](@ref).
"""
function vertex_colors(top::SurfaceTopology, color)
    values = cellvalues(color)
    values isa AbstractVector || return values
    length(values) == nvertices(top) && length(top.extra_tri) == 0 && return values
    return spread(values, top)
end

vertex_colors(mesh::SurfaceMesh, color) = vertex_colors(mesh.topology, color)

function Makie.plot!(plot::DGGSurface{<:Tuple{<:CellRegion, <:AbstractVector}})
    # The topology is the expensive half of a surface, and the heights do not
    # touch it — but the compute graph reports `cells` as changed on every
    # update, however identical it is (see `samebuild`), so the rebuild is
    # skipped here rather than by the edge.
    Makie.register_computation!(
        plot.attributes, [:cells, :transform_func, :wrap, :ntasks], [:surfacetopology]
    ) do inputs, changed, cached
        target = plot_target(inputs.transform_func)
        inputs.wrap || (target = uncut(target))
        if cached !== nothing
            previous = cached[1]
            previous isa SurfaceTopology &&
                samebuild(previous, target, inputs.cells) && return nothing
        end
        return (surface_topology(target, inputs.cells;
            ntasks = task_count(inputs.ntasks)),)
    end

    # The vertex buffer is the one thing a height change costs, so it is written
    # over rather than allocated again.  `reusable` is what keeps that from
    # scribbling on a buffer the surface does not own — the topology's own
    # positions, which a flat surface hands out as its vertices.
    Makie.register_computation!(
        plot.attributes, [:surfacetopology, :zs, :ntasks], [:mesh_positions]
    ) do inputs, changed, cached
        top, zs = inputs.surfacetopology, inputs.zs
        zs isa ZeroHeights && return (vertex_positions(top, zs),)
        buffer = reusable(cached, Point3d, nvertices(top), top.positions)
        buffer === nothing &&
            return (vertex_positions(top, zs, task_count(inputs.ntasks)),)
        return (vertex_positions!(buffer, top, zs, task_count(inputs.ntasks)),)
    end

    Makie.map!(top -> (top.faces,), plot, [:surfacetopology], [:mesh_faces])
    # `color = nothing` is `surface`'s way of saying "colour it by its heights".
    Makie.map!((top, c, zs) -> (vertex_colors(top, c === nothing ? zs : c),), plot,
        [:surfacetopology, :color, :zs], [:mesh_color])
    # Kept for introspection: an O(1) view over the two halves, sharing their
    # arrays rather than copying them.
    Makie.map!((top, p) -> (SurfaceMesh(top, p),), plot,
        [:surfacetopology, :mesh_positions], [:surfacemesh])

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
