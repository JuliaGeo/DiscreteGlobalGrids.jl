# # Plot targets
#
# A cell boundary is a ring on the unit sphere.  Where that ring belongs in the
# axis's data space is the axis's business, and the axis states it in exactly
# one place: its transform function, which every plot carries as the
# `:transform_func` attribute.  A `PlotTarget` is that transform function read
# once, up front, into the few numbers the mesh builder actually needs — so the
# per-vertex inner loop is arithmetic and never a dynamic dispatch or a `ccall`.

"""
    abstract type PlotTarget

Where a cell mesh is built: the axis's transform function, interrogated.

Subtypes answer three questions for [`tessellate`](@ref):

  * `pointtype(target)` — the coordinate type of a mesh vertex.
  * `needs_cutting(target)` — whether the target space has a seam that a cell
    can straddle, so that such cells must be split.
  * how a single cell corner becomes a vertex, and how the finished vertex
    buffer is projected.

Add support for a new axis by adding a [`plot_target`](@ref) method for its
transform function.
"""
abstract type PlotTarget end

"""
    plot_target(transform_func) -> PlotTarget

The [`PlotTarget`](@ref) for a Makie transform function.

The fallback treats an unrecognised transform function as a planar projection
from longitude/latitude with its cut at the antimeridian, which is right for a
plain `Axis` (`identity`) and for any point-wise callable.  `GeoMakie`'s axes
are handled by this package's GeoMakie extension.
"""
function plot_target end

# ## Planar targets

"""
    PlanarTarget(projection, cut)

A 2D target: cells are built in longitude/latitude degrees and then pushed
through `projection` in bulk.

`projection` maps `(lon, lat)` in degrees to the axis's data space; `identity`
means the axis plots longitude/latitude directly.

`cut` is the longitude of the map's seam — the meridian opposite the
projection's central meridian, `lon_0 + 180`.  The drawable domain is
`[cut - 360, cut]`, and a cell that straddles `cut` is split there and the
outside piece is moved a full turn so that it reappears on the far edge.  This
is the "extract the antimeridian from the projection" step, and it is why the
target is derived from the transform function rather than assumed.

`cut = NaN` says the target has no seam, and turns both the split and the
longitude wrapping off.
"""
struct PlanarTarget{TF} <: PlotTarget
    projection::TF
    cut::Float64
end

PlanarTarget(projection) = PlanarTarget(projection, 180.0)

pointtype(::PlanarTarget) = Point2d
needs_cutting(target::PlanarTarget) = !isnan(target.cut)

# The fallback: anything callable is assumed to be a lon/lat projection.  The
# antimeridian is the only cut we can name without knowing more.
plot_target(tf) = PlanarTarget(tf, 180.0)

# A plain `Axis` states its scales as one function per dimension, and the usual
# case is `(identity, identity)` — worth recognising, because it turns the
# projection pass into a no-op.
function plot_target(tf::Tuple)
    all(f -> f === identity, tf) && return PlanarTarget(identity, 180.0)
    return PlanarTarget(tf, 180.0)
end

"""
    project!(target::PlotTarget, positions) -> positions

Convert a vertex buffer from the target's build space into the axis's data
space, in place.

For a [`PlanarTarget`](@ref) the buffer arrives holding `(lon, lat)` in degrees.
The fallback applies the projection point by point; the Proj extension replaces
it with a single strided `proj_trans_generic` call over the whole buffer, which
is where the bulk of a `GeoAxis` plot's time would otherwise go.
"""
function project! end

project!(::PlanarTarget{typeof(identity)}, positions::Vector{Point2d}) = positions

function project!(target::PlanarTarget, positions::Vector{Point2d})
    f = target.projection
    # `apply_transform` is Makie's own dispatch over what a transform function
    # may be — a callable, a per-dimension tuple of scales, a `PointTrans` — so
    # going through it keeps this method as general as the axis is.
    @inbounds for i in eachindex(positions)
        p = Makie.apply_transform(f, positions[i])
        positions[i] = Point2d(p[1], p[2])
    end
    return positions
end

# ## Globe targets

"""
    GlobeTarget(transform, a, e2, height)

A 3D target: cells are built straight from their unit-sphere corners into the
axis's earth-centred coordinates, and nothing is ever split, because a globe has
no seam.

`a` is the ellipsoid's semi-major axis and `e2` its first eccentricity squared,
both in the axis's units, and `height` is the constant offset the axis's
transform adds (`GeoMakie`'s `zlevel`).  Together they turn a unit vector into
the same point the axis's own transform would produce, using three multiplies
and a square root instead of a `ccall`:

    N = a / sqrt(1 - e2 * z^2)
    (X, Y, Z) = ((N + h) * x, (N + h) * y, (N * (1 - e2) + h) * z)

which is exact because a unit vector's components *are* `cos(lat)cos(lon)`,
`cos(lat)sin(lon)` and `sin(lat)`.

See [`probe_ellipsoid`](@ref) for how `a` and `e2` are recovered from a
transform function that does not publish them.
"""
struct GlobeTarget{TF} <: PlotTarget
    transform::TF
    a::Float64
    e2::Float64
    height::Float64
end

"""
    GlobeTarget()

The unit sphere: a globe of radius one, with no flattening and no offset.

The default space for a mesh built outside a plot, where there is no axis to
take an ellipsoid from.  A height is then in units of that radius, so real
elevations want scaling — `zs ./ 6.371e6` puts metres on it — and any other
ellipsoid is the four-argument form.
"""
GlobeTarget() = GlobeTarget(identity, 1.0, 0.0, 0.0)

pointtype(::GlobeTarget) = Point3d
needs_cutting(::GlobeTarget) = false

project!(::GlobeTarget, positions) = positions

"""
    globe_vertex(target::GlobeTarget, p, above = 0.0) -> Point3d

The unit-sphere point `p`, raised `above` the ellipsoid, in the globe axis's
coordinates.

Height enters the formula only as `h`, and only additively, so raising a point
moves it straight out along `p`: `globe_vertex(target, p, above)` is
`globe_vertex(target, p) + above * p`.  That is what lets a surface carry a
height per cell without carrying a normal per vertex.
"""
@inline function globe_vertex(target::GlobeTarget, p, above = 0.0)
    x, y, z = p[1], p[2], p[3]
    N = target.a / sqrt(1 - target.e2 * z * z)
    h = target.height + above
    return Point3d((N + h) * x, (N + h) * y, (N * (1 - target.e2) + h) * z)
end

"""
    probe_ellipsoid(tf, height) -> (a, e2)

Recover the ellipsoid a globe transform projects onto by asking it where two
known points go.

A transform to earth-centred coordinates sends `(0°, 0°)` to `(a + h, 0, 0)` and
`(0°, 90°)` to `(0, 0, b + h)`, so two calls give both axes whatever library is
behind the transform and however it was configured.  This keeps
[`GlobeTarget`](@ref) free of any assumption about WGS84 — a globe on a sphere,
on Mars, or in kilometres all come out right.
"""
function probe_ellipsoid(tf, height::Real)
    equator = tf(Point2d(0.0, 0.0))
    pole = tf(Point2d(0.0, 90.0))
    a = abs(equator[1]) - height
    b = abs(pole[3]) - height
    (isfinite(a) && isfinite(b) && a > 0 && b > 0) ||
        throw(ArgumentError("could not probe the ellipsoid of $tf: got a = $a, b = $b"))
    return a, 1 - (b / a)^2
end

# ## Build space

"""
    buildspace(target::PlotTarget) -> PlotTarget

The target with its projection dropped: the part of it a mesh *build* reads.

A build never calls the projection.  It reads numbers out of the target — the
seam's longitude on a flat map, the ellipsoid and the offset on a globe — and
the projection is applied once at the end, over the finished vertex buffer.
The space is those numbers, and nothing else.

Keeping the two apart is what lets an axis re-project a live plot.  A
[`SurfaceTopology`](@ref) holds the space rather than the target, so its type
does not name the transform function: an axis that swaps `identity` for a
`Proj.Transformation` leaves the topology's type — and its contents — alone, and
costs one bulk [`project!`](@ref) instead of a rebuild.

The space is `isbits`, which is what makes it a value [`samebuild`](@ref) can
compare.
"""
buildspace(target::PlanarTarget) = PlanarTarget(identity, target.cut)
buildspace(target::GlobeTarget) =
    GlobeTarget(identity, target.a, target.e2, target.height)

"""
    needs_projection(target::PlotTarget) -> Bool

Does [`project!`](@ref) do anything on this target?

False for the two targets whose build space already *is* their data space: a
flat map plotting longitude and latitude, and any globe, which is finished by
[`globe_vertex`](@ref) instead.  Asked before a buffer to project into is
reached for, since those are the cases that need none.
"""
needs_projection(::PlotTarget) = true
needs_projection(::PlanarTarget{typeof(identity)}) = false
needs_projection(::GlobeTarget) = false

"""
    project_probe(target, p::UnitSphericalPoint) -> Point2d or Point3d

Where one unit-sphere point lands in the space the target draws in.

The mesh builder projects whole buffers at once, which is what makes a large
cell set affordable; this is the other case — a handful of probe points, asked
about one at a time, when something needs to know where a cell would appear
before deciding whether to build it.
"""
function project_probe end

project_probe(target::GlobeTarget, p) = globe_vertex(target, p)

function project_probe(target::PlanarTarget, p)
    buffer = [Point2d(atand(p[2], p[1]), asind(clamp(p[3], -1.0, 1.0)))]
    project!(target, buffer)
    return buffer[1]
end
