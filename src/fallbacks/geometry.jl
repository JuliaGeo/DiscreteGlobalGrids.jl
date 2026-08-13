# ---------------------------------------------------------------------------
# Geometry generics
#
# Everything a grid can say about a cell's shape, derived from the one required
# primitive `cell_boundary`. A system overrides any of these for speed; none of
# them may be overridden for a different answer.
#
# The boundary contract these are all written against: implicitly closed (the
# first vertex is NOT repeated), counter-clockwise seen from outside the
# sphere, unit-norm points, great-circle arcs between consecutive vertices.
# ---------------------------------------------------------------------------

"""
    closed_ring(points) -> Vector{UnitSphericalPoint{Float64}}

The implicitly closed [`cell_boundary`](@ref) ring made explicitly closed: the
first vertex appended at the end, which is what a GeoInterface `LinearRing`
wants. A ring that already repeats its first vertex is passed through unchanged
— boundaries are implicitly closed by contract, but a system that closes its
own is not thereby made to describe a doubled edge.
"""
function closed_ring(points)
    n = length(points)
    n == 0 && return USPoint[]
    out = Vector{USPoint}(undef, n + 1)
    # `enumerate` rather than `points[i]`: the boundary contract allows any
    # `AbstractVector`, including one that is not 1-based.
    for (i, p) in enumerate(points)
        @inbounds out[i] = USPoint(p[1], p[2], p[3])
    end
    @inbounds out[n+1] = out[1]
    # Already closed: drop the duplicate we just introduced.
    n > 1 && out[n] == out[1] && return out[1:n]
    return out
end

"""
    open_ring(points) -> (ring, n)

A boundary as a **1-based** vector together with its **open** length: `ring[1:n]`
are the distinct vertices, with the closing edge `ring[n] -> ring[1]` implied.

Two conversions in one, both needed by the spherical ring predicates, which
index `pts[1:n]` positionally and count the closing edge themselves:

  - a boundary that is not 1-based is collected (the contract allows any
    `AbstractVector`);
  - a boundary that repeats its first vertex at the end — which the contract
    says it should not, but which [`closed_ring`](@ref) also tolerates — has
    that vertex dropped from the count. Left in, it becomes a retraced
    zero-length edge, and a retraced edge silently breaks the crossing parity
    those predicates decide containment by.
"""
function open_ring(points)
    ring = firstindex(points) == 1 ? points : collect(points)
    n = length(ring)
    n > 1 && ring[n] == ring[1] && (n -= 1)
    return ring, n
end

"""
    point_in_cell(boundary, p) -> Union{Bool,Nothing}

Whether the unit-sphere point `p` lies in the cell bounded by `boundary`
(boundary points count as inside), or `nothing` when the question is
genuinely undecidable.

Two independent algorithms, because each has a degenerate case the other does
not: `spherical_ring_encloses` bootstraps from a definitionally exterior
anchor — the antipode of the ring's vertex mass — and gives up when `p` sits
near that mass, i.e. near the middle of the cell, which is exactly where a
caller most wants an answer. `spherical_ring_contains` bootstraps instead from
a wedge at one ring edge, valid here because the boundary contract fixes the
winding (counter-clockwise seen from outside). Only if both decline does this
return `nothing`.
"""
function point_in_cell(boundary, p)
    ring, n = open_ring(boundary)
    n >= 3 || return nothing
    verdict = US.spherical_ring_encloses(ring, n, p)
    verdict === nothing || return verdict
    return US.spherical_ring_contains(ring, n, p)
end

"""
    cell_polygon(grid, c) -> GI.Polygon

See the interface docstring. One `LinearRing`, explicitly closed, unit-sphere
`(x, y, z)` coordinates.
"""
cell_polygon(grid::AbstractGrid, c::AbstractCellIndex) =
    GI.Polygon([GI.LinearRing(closed_ring(cell_boundary(grid, c)))])

"""
    getcell(grid, i) -> GI.Polygon

`ConservativeRegridding.Trees.getcell`: position -> unit-sphere polygon.
Implemented once, here, as `cell_polygon(grid, cellindex(grid, i))`; grid
authors never write it.
"""
getcell(grid::AbstractGrid, i::Int) = cell_polygon(grid, cellindex(grid, i))

# `Trees` addresses a grid as a tree source too, and its lazy all-cells form is
# part of the surface every grid gets.
getcell(grid::AbstractGrid) = (getcell(grid, i) for i in 1:ncells(grid))

"""
    cell_area(grid, c) -> Float64

The spherical area of the cell in **steradians**, from the exact ring.

Delegates to `GeometryOps.area` on a unit-radius `Spherical` manifold, which
fan-triangulates the ring and sums signed spherical excess — right for cells
that span a large solid angle, where a planar formula is wrong in sign rather
than merely inaccurate.
"""
cell_area(grid::AbstractGrid, c::AbstractCellIndex) =
    Float64(GO.area(GO.Spherical(; radius=1.0), cell_polygon(grid, c)))

# ===========================================================================
# Longitude/latitude extent
# ===========================================================================

"""
    cell_extent(grid, c) -> Extents.Extent{(:X, :Y)}

The lon/lat bounding box of a cell, in degrees, conservative where a rectangle
cannot be exact (see the interface docstring). This is an interoperability
convenience — tree pruning uses [`node_extent`](@ref), a `SphericalCap`.

Three things a naive vertex bounding box gets wrong, and what is done instead:

  - **Arc bulge.** A great-circle edge between two vertices reaches a higher
    latitude than either endpoint. Each edge's extreme-latitude point is
    computed and included when it actually falls on the arc.
  - **The antimeridian.** An edge whose endpoints straddle ±180° makes the
    vertex longitude range meaningless, so `X` becomes `(-180, 180)`.
  - **Poles.** A cell enclosing a pole has no longitude bound at all, and its
    latitude range runs to ±90. Enclosure is decided by the ring's longitude
    winding, and which pole by a spherical containment test.
"""
function cell_extent(grid::AbstractGrid, c::AbstractCellIndex)
    points, n = open_ring(cell_boundary(grid, c))
    n == 0 && throw(ArgumentError("cell $c has an empty boundary"))

    latmin = latmax = 0.0
    lonmin = 180.0
    lonmax = -180.0
    winding = 0.0
    crosses = false
    first_lon = 0.0
    previous_lon = 0.0
    for i in 1:n
        p = points[i]
        lon, lat = lonlat(p)
        if i == 1
            latmin = latmax = lat
            first_lon = lon
        else
            latmin = min(latmin, lat)
            latmax = max(latmax, lat)
            delta = lon - previous_lon
            if delta > 180.0
                delta -= 360.0
                crosses = true
            elseif delta < -180.0
                delta += 360.0
                crosses = true
            end
            winding += delta
        end
        previous_lon = lon
        lonmin = min(lonmin, lon)
        lonmax = max(lonmax, lon)
        # Arc bulge: the closing edge is included by the wrap-around index.
        q = points[i == n ? 1 : i + 1]
        elat = _arc_lat_extremes(p, q)
        if elat !== nothing
            latmin = min(latmin, elat[1])
            latmax = max(latmax, elat[2])
        end
    end
    # The closing edge's longitude step, which the loop above never took.
    delta = first_lon - previous_lon
    if delta > 180.0
        delta -= 360.0
        crosses = true
    elseif delta < -180.0
        delta += 360.0
        crosses = true
    end
    winding += delta

    if abs(winding) > 180.0
        # The ring winds around a pole: no longitude bound exists, and the
        # enclosed pole is at the latitude limit. An undecidable ring gives up
        # on both poles rather than guessing one.
        north = point_in_cell(points, USPoint(0.0, 0.0, 1.0))
        south = point_in_cell(points, USPoint(0.0, 0.0, -1.0))
        ymin = south === false ? latmin : -90.0
        ymax = north === false ? latmax : 90.0
        return Extents.Extent(X=(-180.0, 180.0), Y=(ymin, ymax))
    end
    crosses && return Extents.Extent(X=(-180.0, 180.0), Y=(latmin, latmax))
    return Extents.Extent(X=(lonmin, lonmax), Y=(latmin, latmax))
end

# The latitude extremes of the great-circle arc `a -> b`, or `nothing` when
# neither extreme falls on the arc (the common case, where the endpoints
# already bound it).
#
# The great circle with unit normal `n` attains its extreme latitudes at
# `+-normalize(z - (z . n) n)`; that point is on the minor arc `a -> b` exactly
# when it lies strictly between them, which is the same pair of triple products
# the boundary-arc distance uses.
function _arc_lat_extremes(a, b)
    nx = a[2] * b[3] - a[3] * b[2]
    ny = a[3] * b[1] - a[1] * b[3]
    nz = a[1] * b[2] - a[2] * b[1]
    nn = nx * nx + ny * ny + nz * nz
    nn <= 1e-24 && return nothing            # degenerate or antipodal edge
    inv = 1.0 / sqrt(nn)
    ux, uy, uz = nx * inv, ny * inv, nz * inv
    # z - (z . u) u, the tangent direction of steepest latitude gain.
    vx, vy, vz = -uz * ux, -uz * uy, 1.0 - uz * uz
    vn = sqrt(vx * vx + vy * vy + vz * vz)
    # `u` parallel to the pole means the arc's great circle IS the equator,
    # where latitude is constant and there is no extreme to find.
    vn <= 1e-12 && return nothing
    px, py, pz = vx / vn, vy / vn, vz / vn
    lo = nothing
    hi = nothing
    for s in (1.0, -1.0)
        qx, qy, qz = s * px, s * py, s * pz
        # `a -> q -> b` in the orientation of `n`.
        ((a[2] * qz - a[3] * qy) * nx + (a[3] * qx - a[1] * qz) * ny +
         (a[1] * qy - a[2] * qx) * nz) > 0 || continue
        ((qy * b[3] - qz * b[2]) * nx + (qz * b[1] - qx * b[3]) * ny +
         (qx * b[2] - qy * b[1]) * nz) > 0 || continue
        lat = asind(clamp(qz, -1.0, 1.0))
        lo = lo === nothing ? lat : min(lo, lat)
        hi = hi === nothing ? lat : max(hi, lat)
    end
    lo === nothing && return nothing
    return (lo, hi)
end

# ===========================================================================
# The default node extent
# ===========================================================================

"""
    node_extent(sys, c) -> SphericalCap

The generic default: the cell's own bounding cap inflated by
[`cap_inflation(sys)`](@ref cap_inflation), which is O(1) at every depth and
sound for every system whose measured child overhang stays under that factor.

See the covering law in the interface docstring — this default is the reason
`cap_inflation` exists, and a system that overrides `node_extent` outright
(HEALPix's exact subtree cap) ignores it.
"""
node_extent(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex) =
    inflate_cap(cell_cap(levelgrid(sys, level(c)), c), cap_inflation(sys))

# ===========================================================================
# Manifold
# ===========================================================================

"""
    GeometryOpsCore.best_manifold(grid::AbstractGrid) -> GO.Spherical

The compute manifold of every grid in this package: the **unit** sphere, which
is the frame all its geometry is expressed in.

`radius = 1` is not a placeholder. Cell coordinates are unit-sphere `(x, y, z)`,
so an area computed against this manifold comes out in steradians — see
[`cell_area`](@ref) — and scaling to a physical sphere is the caller's one
multiplication by `R^2`. Handing `ConservativeRegridding` a manifold whose
radius disagreed with the coordinates would put that factor in twice.
"""
GOCore.best_manifold(grid::AbstractGrid) = GO.Spherical(; radius=1.0)
