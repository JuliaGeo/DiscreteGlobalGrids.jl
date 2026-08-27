# # Handing the mesh to somebody else
#
# Everything above builds a mesh in order to draw it.  The same mesh is worth
# having on its own — to write out as OBJ or glTF, to hand to a mesh library, or
# to plot with Makie's plain `mesh` — and `GeometryBasics.Mesh` is the type
# every one of those speaks.
#
# The translation is nearly free: a `GeometryBasics.Mesh` is a vertex-attribute
# `NamedTuple` plus faces, and both halves already exist.  What the entry point
# adds is the one thing a plot does that a bare mesh cannot: values arrive per
# *cell*, and a mesh wants them per *vertex*.
#
# `mesh` is `GeometryBasics`', so every method here dispatches on a type this
# package owns — a `CellRegion`, a `CellSet`, or one of the two meshes.  That is
# what keeps them from being piracy, and it costs the caller one call: whatever
# names the cells goes through `cellregion` or `cellset` first, the same two
# verbs the recipes use.

"""
    GeometryBasics.mesh(cells::CellRegion[, elevation]; target, ntasks,
                        connectivity, facetype, attributes...) -> GeometryBasics.Mesh

The cells as one interpolated surface, outside any plot.

`cells` is a [`CellRegion`](@ref) — `cellregion(x)` reads one out of anything
[`dggsurface`](@ref) takes — and the result is the surface that plot would draw:
one vertex per cell centroid, joined by the triangles of the grid's dual.
`elevation` is one height per cell and lands in the geometry exactly as it does
in the plot; without it the surface is flat.

Extra keyword arguments become **vertex attributes** of the mesh, under the
names given.  A vector of one value per cell is spread to the vertices the way a
colour is — the split vertices of a cut triangle take the mix their position
earns, so nothing steps at the seam — and a vector already one per vertex is
taken as it is.  `normal`, `uv` and `color` are the names most readers look for.

A one-dimensional `DimArray` over cells is cells and heights at once, as it is
for [`dggsurface`](@ref), but it takes both verbs to say so:

```julia
using GeometryBasics, FileIO
save("terrain.obj", GeometryBasics.mesh(cellregion(A), A))
GeometryBasics.mesh(cellregion(grid), heights; color = temperature)
```

`target` is the space the mesh lives in, a [`PlotTarget`](@ref).  The default is
the unit sphere, `GlobeTarget()`, on which a height is in units of that radius;
`PlanarTarget(identity)` gives longitude and latitude instead, with cells cut at
the antimeridian.  `ntasks` and `connectivity` are [`triangulate`](@ref)'s, and
`facetype` is `GeometryBasics`' own.

Pass a [`CellSet`](@ref) instead — `cellset(x)` — to get the flat patches
[`dggpoly`](@ref) draws.
"""
function GeometryBasics.mesh(cells::CellRegion, elevation = nothing;
        target::PlotTarget = GlobeTarget(), ntasks::Int = Threads.nthreads(),
        connectivity = DGG.Vertex(), facetype = GLTriangleFace, attributes...)
    zs = elevation === nothing ? ZeroHeights(length(cells)) : cellvalues(elevation)
    m = triangulate(target, cells, zs; ntasks = ntasks, connectivity = connectivity)
    return GeometryBasics.mesh(m; facetype = facetype, attributes...)
end

"""
    GeometryBasics.mesh(cells::CellSet; target, ntasks, facetype,
                        attributes...) -> GeometryBasics.Mesh

The cells as one mesh of flat patches: what [`dggpoly`](@ref) draws.

Each cell contributes its own ring of vertices, so the cells meet at their real
boundaries and a value belongs to the whole cell rather than fading into its
neighbours.  Per-cell attributes are handed to every vertex of the cell.

This is the form for a set with no adjacency to read — a `MultiOrderCellSet`
spanning several levels, or a bare list of ids — and the one that keeps cell
edges.  It carries no heights: every vertex of a cell would sit at the same
height, so neighbouring cells would meet at a step rather than a slope.
"""
function GeometryBasics.mesh(cells::CellSet; target::PlotTarget = GlobeTarget(),
        ntasks::Int = Threads.nthreads(), facetype = GLTriangleFace, attributes...)
    m = tessellate(target, cells; ntasks = ntasks)
    return GeometryBasics.mesh(m; facetype = facetype, attributes...)
end

GeometryBasics.mesh(::CellSet, ::AbstractVector; kwargs...) =
    throw(ArgumentError("a patch mesh carries no heights, because every vertex \
        of a cell would take the cell's own height and neighbouring cells would \
        meet at a step. `cellregion` instead of `cellset` gives the interpolated \
        surface, which is what heights are for."))

"""
    GeometryBasics.mesh(m::SurfaceMesh; facetype, attributes...)
    GeometryBasics.mesh(m::CellMesh; facetype, attributes...)

A mesh this package has already built, as a `GeometryBasics.Mesh`.

The vertex buffer and the faces are shared, not copied — the conversion costs
only what the attributes cost.
"""
function GeometryBasics.mesh(m::SurfaceMesh; facetype = GLTriangleFace, attributes...)
    return GeometryBasics.mesh(m.positions, m.faces;
        facetype = facetype, pervertex(m, values(attributes))...)
end

function GeometryBasics.mesh(m::CellMesh; facetype = GLTriangleFace, attributes...)
    return GeometryBasics.mesh(m.positions, m.faces;
        facetype = facetype, pervertex(m, values(attributes))...)
end

"""
    pervertex(mesh, values) -> value per vertex
    pervertex(mesh, attributes::NamedTuple) -> NamedTuple

One value per cell as one value per vertex, for a mesh of either kind.

A vector as long as the mesh has vertices is already what it needs to be — on a
globe, where nothing is cut, a surface's per-cell vector is both.  Anything else
is spread: [`spread`](@ref) blends a surface's split vertices, and a patch mesh
hands each cell's value to every vertex of its ring.  A value that is not a
vector at all — a single colour, a `FaceView` — passes through untouched.
"""
pervertex(m, attributes::NamedTuple) =
    NamedTuple{keys(attributes)}(map(v -> pervertex(m, v), values(attributes)))

pervertex(::Any, value) = value

pervertex(m::SurfaceMesh, values::AbstractVector) =
    length(values) == length(m.positions) ? cellvalues(values) :
        spread(cellvalues(values), m.topology)

function pervertex(m::CellMesh, values::AbstractVector)
    length(values) == length(m.positions) && return cellvalues(values)
    vs = cellvalues(values)
    length(vs) == m.ncells || throw(ArgumentError("got $(length(vs)) values, but \
        there are $(m.ncells) cells and $(length(m.positions)) vertices"))
    return vs[m.vertex_cell]
end
