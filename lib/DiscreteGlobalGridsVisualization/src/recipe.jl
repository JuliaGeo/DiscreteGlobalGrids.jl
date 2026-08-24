# # The recipe
#
# `dggpoly` reads like `poly`: cells in, one filled patch per cell out, coloured
# by a vector as long as the cell set.  What it does underneath depends on where
# the plot lives — see `targets.jl` for the space and `tessellate.jl` for the
# mesh — and on which backend is drawing it.

"""
    dggpoly(cells; color = ..., kwargs...)
    dggpoly(source, ids; color = ..., kwargs...)

Draw a set of DGGS cells as filled patches.

`cells` is anything [`cellset`](@ref) accepts: an `AbstractGrid`, a
`CellVector`, a `CellLookup`, a `MultiOrderCellSet`, a one-dimensional
`DimArray` over cells, or a system paired with a vector of cell ids.  Partial
grids are the normal case, not a special one — the cells drawn are exactly the
cells given.

`color` is a single colour, or one value or colour **per cell**.  A cube axis
names only its cells here — `dggpoly(A)` draws them in one flat colour, and
`dggpoly(A; color = A)` colours them by its values; it is [`dggsurface`](@ref)
that reads a cube axis as heights.

```julia
sys = DGG.IGeo7System()
cells = DGG.CellVector(DGG.query(sys, DGG.MultiOrderCoverage(extent); level = 9))
dggpoly(cells; color = elevation)
```

The plot adapts to the axis it is placed in.  In a `GeoMakie.GlobeAxis` the
cells are built directly on the globe from their unit-sphere corners; in an
`Axis` or a `GeoAxis` they are built in longitude/latitude, split where they
straddle the map's cut meridian, and projected in bulk.  In both cases the
projection is done here rather than by the axis, which is why the plot can be
one mesh with one draw call.
"""
@recipe DGGPoly (cells,) begin
    """
    Sets the colour of the cells.  Either a single colour, or a vector with one
    entry per cell — numbers to be mapped through the colormap, or colours.
    """
    color = @inherit patchcolor
    "Sets the colour of the outline around each cell."
    strokecolor = @inherit patchstrokecolor
    "Sets the width of the cell outlines.  `0` draws no outlines, which is the default because outlines cost as much as the fill."
    strokewidth = 0
    "Sets the dash pattern of the cell outlines."
    linestyle = nothing
    "Controls the rendering of outline corners."
    joinstyle = @inherit joinstyle
    "Sets the minimum inner join angle below which miter line joins truncate."
    miter_limit = @inherit miter_limit
    "Controls whether lights affect the cells.  Off by default: a DGGS patch is a flat value, not a surface."
    shading = Makie.NoShading
    """
    How the cells are handed to the backend.  `:mesh` builds one triangle mesh
    for the whole set, which is what GPU backends want.  `:polygons` emits one
    filled path per cell, which is what CairoMakie wants — a vector renderer
    draws a path faster than a fan of triangles, and without the hairline seams
    that show up between triangles of the same cell.  `automatic` picks
    `:polygons` under CairoMakie and `:mesh` otherwise.
    """
    primitive = Makie.automatic
    """
    Whether to split cells that straddle the map's cut meridian.  Only planar
    targets have a cut; on a globe this attribute does nothing.
    """
    wrap = true
    "How many tasks build the mesh.  `automatic` is one per thread."
    ntasks = Makie.automatic
    "Depth shift of the outlines, to keep them from z-fighting with the fill."
    stroke_depth_shift = -1.0f-5
    cycle = [:color => :patchcolor]
    Makie.mixin_generic_plot_attributes()...
    Makie.mixin_colormap_attributes()...
end

Makie.convert_arguments(::Type{<:DGGPoly}, cs::CellSet) = (cs,)
Makie.convert_arguments(::Type{<:DGGPoly}, x) = (cellset(x),)
Makie.convert_arguments(::Type{<:DGGPoly}, source, ids::AbstractVector) = (cellset(source, ids),)

# `wrap = false` is expressed by moving the cut out of reach rather than by a
# branch in the inner loop.
uncut(target::PlanarTarget) = PlanarTarget(target.projection, NaN)
uncut(target::PlotTarget) = target

"""
    task_count(ntasks) -> Int

How many tasks to build a mesh with: one per thread unless the plot names a
number.

The session's thread count cannot be written as the attribute's default.  A
recipe's defaults are evaluated once, where the `@recipe` block is read — which
is during precompilation, in a worker that has one thread — so
`ntasks = Threads.nthreads()` bakes `1` into the package image and every plot
gets it, however many threads the session was started with.  `automatic` is
resolved here instead, when the mesh is actually built.
"""
task_count(::Makie.Automatic) = Threads.nthreads()
task_count(ntasks::Integer) = Int(ntasks)

"""
    resolve_primitive(primitive, target) -> :mesh or :polygons

Decide how to hand the cells to the backend.

`automatic` asks CairoMakie for paths and every other backend for a mesh.  A
globe plot is three-dimensional, and a vector renderer has no way to draw a
3D path, so it gets the mesh whatever the backend.
"""
function resolve_primitive(primitive, target)
    primitive === Makie.automatic || return primitive
    (pointtype(target) === Point2d && is_cairo_backend()) || return :mesh
    return :polygons
end

function is_cairo_backend()
    backend = Makie.current_backend()
    backend === missing && return false
    return nameof(backend) === :CairoMakie
end

"""
    vertex_colors(mesh::CellMesh, color)

Spread a per-cell colour over the mesh's vertices.

A vector as long as the cell set is gathered through `mesh.vertex_cell`; anything
else — a single colour, or something already as long as the vertex buffer — is
passed through untouched.  The gather is the only work a recolour costs: the
geometry never moves.
"""
function vertex_colors(mesh::CellMesh, color)
    color isa AbstractVector || return color
    length(color) == length(mesh.positions) && return color
    length(color) == mesh.ncells ||
        throw(ArgumentError("color has $(length(color)) entries, but there are \
            $(mesh.ncells) cells"))
    return @inbounds [color[Int(c)] for c in mesh.vertex_cell]
end

"""
    cell_polygons(mesh::CellMesh, ntasks = Threads.nthreads()) -> Vector{Makie.Polygon}

The mesh's rings as one polygon each, for backends that draw paths.
"""
function cell_polygons(mesh::CellMesh{P}, ntasks::Int = Threads.nthreads()) where {P}
    polys = Vector{Makie.Polygon{2, Float64}}(undef, nrings(mesh))
    starts = mesh.ring_start
    inparallel(nrings(mesh), ntasks) do lo, hi
        @inbounds for i in lo:hi
            polys[i] = Makie.Polygon(mesh.positions[starts[i]:(starts[i + 1] - 1)])
        end
    end
    return polys
end

"""
    ring_colors(mesh::CellMesh, color)

`color` reduced to one entry per drawn ring, for the polygon path.
"""
function ring_colors(mesh::CellMesh, color)
    color isa AbstractVector || return color
    length(color) == mesh.ncells ||
        throw(ArgumentError("color has $(length(color)) entries, but there are \
            $(mesh.ncells) cells"))
    starts = mesh.ring_start
    return @inbounds [color[Int(mesh.vertex_cell[starts[i]])] for i in 1:nrings(mesh)]
end

"""
    outline_points(mesh::CellMesh) -> Vector

Every drawn ring as a closed loop, the loops separated by `NaN` points, ready
for a single `lines!`.
"""
function outline_points(mesh::CellMesh{GeometryBasics.Point{N, T}}) where {N, T}
    P = GeometryBasics.Point{N, T}
    n = nrings(mesh)
    out = Vector{P}(undef, length(mesh.positions) + 2n)
    nan = P(ntuple(_ -> T(NaN), N))
    k = 0
    starts = mesh.ring_start
    @inbounds for i in 1:n
        lo, hi = starts[i], starts[i + 1] - 1
        for j in lo:hi
            out[k += 1] = mesh.positions[j]
        end
        out[k += 1] = mesh.positions[lo]
        out[k += 1] = nan
    end
    return out
end

function Makie.plot!(plot::DGGPoly{<:Tuple{<:CellSet}})
    Makie.map!(
        plot, [:cells, :transform_func, :wrap, :ntasks], [:cellmesh]
    ) do cells, transform_func, wrap, ntasks
        target = plot_target(transform_func)
        wrap || (target = uncut(target))
        return (tessellate(target, cells; ntasks = task_count(ntasks)),)
    end

    # Resolved once, when the plot is created: the backend does not change under
    # a live plot, and keeping both drawing paths alive would double what a large
    # cell set costs in memory.
    target = plot_target(Makie.to_value(plot.transform_func))
    if resolve_primitive(Makie.to_value(plot.primitive), target) === :polygons
        draw_polygons!(plot)
    else
        draw_mesh!(plot)
    end

    draw_strokes!(plot)
    return plot
end

function draw_mesh!(plot::DGGPoly)
    Makie.map!(m -> (m.positions,), plot, [:cellmesh], [:mesh_positions])
    Makie.map!(m -> (m.faces,), plot, [:cellmesh], [:mesh_faces])
    Makie.map!((m, c) -> (vertex_colors(m, c),), plot, [:cellmesh, :color], [:mesh_color])

    return Makie.mesh!(
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
end

function draw_polygons!(plot::DGGPoly)
    Makie.map!((m, n) -> (cell_polygons(m, task_count(n)),), plot,
        [:cellmesh, :ntasks], [:polygons])
    Makie.map!((m, c) -> (ring_colors(m, c),), plot, [:cellmesh, :color], [:polygon_color])

    return Makie.poly!(
        plot, plot.polygons;
        color = plot.polygon_color,
        colormap = plot.colormap,
        colorscale = plot.colorscale,
        colorrange = plot.colorrange,
        lowclip = plot.lowclip,
        highclip = plot.highclip,
        nan_color = plot.nan_color,
        alpha = plot.alpha,
        strokewidth = 0,
        strokecolor = :transparent,
        shading = plot.shading,
        visible = plot.visible,
        transparency = plot.transparency,
        inspectable = plot.inspectable,
        space = plot.space,
        depth_shift = plot.depth_shift,
        clip_planes = plot.clip_planes,
        transformation = :inherit_model,
    )
end

function draw_strokes!(plot::DGGPoly)
    Makie.map!(plot, [:cellmesh, :strokewidth], [:outline]) do mesh, strokewidth
        # Outlines cost as much as the fill, so they are only built when they
        # will actually be drawn.
        strokewidth > 0 || return (similar(mesh.positions, 0),)
        return (outline_points(mesh),)
    end
    Makie.map!(
        (visible, strokewidth) -> (visible && strokewidth > 0,),
        plot, [:visible, :strokewidth], [:stroke_visible]
    )

    return Makie.lines!(
        plot, plot.outline;
        color = plot.strokecolor,
        linewidth = plot.strokewidth,
        linestyle = plot.linestyle,
        joinstyle = plot.joinstyle,
        miter_limit = plot.miter_limit,
        visible = plot.stroke_visible,
        transparency = plot.transparency,
        inspectable = plot.inspectable,
        space = plot.space,
        depth_shift = plot.stroke_depth_shift,
        clip_planes = plot.clip_planes,
        transformation = :inherit_model,
    )
end
