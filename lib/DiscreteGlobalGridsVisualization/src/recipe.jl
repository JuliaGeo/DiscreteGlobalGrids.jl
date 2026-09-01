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

Cells whose colour is `NaN` or `missing` are not drawn.

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
    surface_target(transform_func, wrap) -> PlotTarget

The target a plot draws into: its axis's transform function read once, with the
seam removed if the plot does not want one.

Every recipe here starts from these two attributes, and this is the one place
that turns them into a [`PlotTarget`](@ref).
"""
surface_target(transform_func, wrap::Bool) =
    wrap ? plot_target(transform_func) : uncut(plot_target(transform_func))

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
    isdrawn(x) -> Bool

Whether a cell coloured `x` should be drawn at all.

`false` for `missing` and for a `NaN` float; `true` for everything else,
including an ordinary number and a named colour.  `dggpoly` drops a cell's
faces — or, on the polygon path, its whole ring — rather than hand either kind
of blank value to a backend: GLMakie interpolates a `NaN` vertex colour to the
top of the colormap instead of leaving it unpainted (`nan_color` is never
consulted on that path), and Makie's colour conversion rejects `missing`
outright.
"""
isdrawn(x) = true
isdrawn(::Missing) = false
isdrawn(x::AbstractFloat) = !isnan(x)

"""
    _fill_value(color::AbstractVector)

A finite, in-range stand-in for an undrawn cell's colour.

The mesh keeps one colour slot per cell for every vertex (and, before
filtering, every ring) it built, whether or not that cell ends up drawn:
dropping a cell's faces stops it from being rasterised, but a backend that
scans the whole colour buffer to pick a default `colorrange` — GLMakie among
them — still sees the undrawn slots, and a `NaN` anywhere in that scan can
break the range for every cell, drawn or not.  Filling undrawn slots with the
first drawn colour in `color` (or `0` if none is drawn) keeps the scan finite;
which value is used does not matter, because no undrawn slot is ever
rasterised.
"""
function _fill_value(color::AbstractVector)
    for c in color
        isdrawn(c) && return c
    end
    S = nonmissingtype(eltype(color))
    return S <: Real ? zero(float(S)) : first(color)
end

"""
    drawn_faces(mesh::CellMesh, color) -> Vector{GLTriangleFace}

The mesh's faces with every triangle belonging to an undrawn cell removed.

Every vertex of a face belongs to the same cell — faces never span a cell
boundary, see [`CellMesh`](@ref) — so checking the first vertex is enough.
`color` that is not a per-cell vector (a single colour, or something already
spread to one entry per vertex) keeps every face, as before.
"""
function drawn_faces(mesh::CellMesh, color)
    (color isa AbstractVector && length(color) == mesh.ncells) || return mesh.faces
    return filter(f -> isdrawn(@inbounds color[Int(mesh.vertex_cell[f[1]])]), mesh.faces)
end

"""
    drawn_ring_indices(mesh::CellMesh, color) -> AbstractVector{Int}

Which of the mesh's rings belong to a drawn cell.

Every ring, unless `color` is a per-cell vector, in which case a ring is kept
only when its cell's colour is [`isdrawn`](@ref).  Used by the polygon path
and by the outline path, so that both agree with the mesh path about which
cells get drawn.
"""
function drawn_ring_indices(mesh::CellMesh, color)
    (color isa AbstractVector && length(color) == mesh.ncells) || return 1:nrings(mesh)
    starts = mesh.ring_start
    return findall(i -> isdrawn(@inbounds color[Int(mesh.vertex_cell[starts[i]])]), 1:nrings(mesh))
end

"""
    vertex_colors(mesh::CellMesh, color)

Spread a per-cell colour over the mesh's vertices.

A vector as long as the cell set is gathered through `mesh.vertex_cell`; anything
else — a single colour, or something already as long as the vertex buffer — is
passed through untouched.  The gather is the only work a recolour costs: the
geometry never moves.

A `missing` or `NaN` entry never reaches the output: the vertices of a cell
that fails [`isdrawn`](@ref) get [`_fill_value`](@ref) instead, because their
faces are dropped by [`drawn_faces`](@ref) but their colour slot is not — the
colour buffer this returns holds no `NaN` and no `missing` at all.
"""
function vertex_colors(mesh::CellMesh, color)
    color isa AbstractVector || return color
    length(color) == length(mesh.positions) && return color
    length(color) == mesh.ncells ||
        throw(ArgumentError("color has $(length(color)) entries, but there are \
            $(mesh.ncells) cells"))
    fill_value = _fill_value(color)
    T = typeof(fill_value)
    out = Vector{T}(undef, length(mesh.vertex_cell))
    @inbounds for (i, c) in enumerate(mesh.vertex_cell)
        v = color[Int(c)]
        out[i] = isdrawn(v) ? v : fill_value
    end
    return out
end

"""
    cell_polygons(mesh::CellMesh, color, ntasks = Threads.nthreads()) -> Vector{Makie.Polygon}

The mesh's drawn rings as one polygon each, for backends that draw paths.

A ring whose cell fails [`isdrawn`](@ref) — see [`drawn_ring_indices`](@ref) —
is left out entirely, matching [`ring_colors`](@ref).
"""
function cell_polygons(mesh::CellMesh{P}, color, ntasks::Int = Threads.nthreads()) where {P}
    starts = mesh.ring_start
    keep = drawn_ring_indices(mesh, color)
    m = length(keep)
    polys = Vector{Makie.Polygon{2, Float64}}(undef, m)
    inparallel(m, ntasks) do lo, hi
        @inbounds for j in lo:hi
            i = keep[j]
            polys[j] = Makie.Polygon(mesh.positions[starts[i]:(starts[i + 1] - 1)])
        end
    end
    return polys
end

"""
    ring_colors(mesh::CellMesh, color)

`color` reduced to one entry per drawn ring, for the polygon path.

Unlike [`vertex_colors`](@ref), an undrawn ring is not filled — it is left out
of the result altogether, in the same order [`cell_polygons`](@ref) leaves it
out of the polygon list, so the two stay aligned.
"""
function ring_colors(mesh::CellMesh, color)
    color isa AbstractVector || return color
    length(color) == mesh.ncells ||
        throw(ArgumentError("color has $(length(color)) entries, but there are \
            $(mesh.ncells) cells"))
    starts = mesh.ring_start
    keep = drawn_ring_indices(mesh, color)
    T = nonmissingtype(eltype(color))
    out = Vector{T}(undef, length(keep))
    @inbounds for (j, i) in enumerate(keep)
        out[j] = color[Int(mesh.vertex_cell[starts[i]])]
    end
    return out
end

"""
    outline_points(mesh::CellMesh) -> Vector
    outline_points(mesh::CellMesh, color) -> Vector

Every drawn ring as a closed loop, the loops separated by `NaN` points, ready
for a single `lines!`.

With `color`, a ring whose cell fails [`isdrawn`](@ref) is left out (see
[`drawn_ring_indices`](@ref)), so a positive `strokewidth` never outlines a
blank cell.  Without it, every ring is kept.
"""
function _outline_points(mesh::CellMesh{GeometryBasics.Point{N, T}}, keep) where {N, T}
    P = GeometryBasics.Point{N, T}
    starts = mesh.ring_start
    total = 0
    @inbounds for i in keep
        total += (starts[i + 1] - starts[i]) + 2
    end
    out = Vector{P}(undef, total)
    nan = P(ntuple(_ -> T(NaN), N))
    k = 0
    @inbounds for i in keep
        lo, hi = starts[i], starts[i + 1] - 1
        for j in lo:hi
            out[k += 1] = mesh.positions[j]
        end
        out[k += 1] = mesh.positions[lo]
        out[k += 1] = nan
    end
    return out
end

outline_points(mesh::CellMesh) = _outline_points(mesh, 1:nrings(mesh))
outline_points(mesh::CellMesh, color) = _outline_points(mesh, drawn_ring_indices(mesh, color))

function Makie.plot!(plot::DGGPoly{<:Tuple{<:CellSet}})
    Makie.map!(
        plot, [:cells, :transform_func, :wrap, :ntasks], [:cellmesh]
    ) do cells, transform_func, wrap, ntasks
        target = surface_target(transform_func, wrap)
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
    # `mesh_faces` depends on `color` too: an undrawn cell's faces are dropped
    # from the mesh, not just recoloured, because GLMakie interpolates a NaN
    # vertex colour into the colormap rather than leaving it unpainted.
    Makie.map!((m, c) -> (drawn_faces(m, c),), plot, [:cellmesh, :color], [:mesh_faces])
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
    Makie.map!((m, c, n) -> (cell_polygons(m, c, task_count(n)),), plot,
        [:cellmesh, :color, :ntasks], [:polygons])
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
    Makie.map!(plot, [:cellmesh, :color, :strokewidth], [:outline]) do mesh, color, strokewidth
        # Outlines cost as much as the fill, so they are only built when they
        # will actually be drawn.  An undrawn cell gets no outline either.
        strokewidth > 0 || return (similar(mesh.positions, 0),)
        return (outline_points(mesh, color),)
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
