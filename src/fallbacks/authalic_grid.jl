# Authalic wrappers convert geometry between authalic and geodetic latitude.
# Cell ids, indices, hierarchy, and ordering are unchanged. The two latitudes
# differ by up to 0.1283° for WGS84 — ~14.3 km along a meridian at ±45°, wider
# than a level-8 cell — so overlaying the frames without this is a silent
# misregistration.

# ===========================================================================
# The two bounds the wrapper is built on
# ===========================================================================

# Explicitly fixes the wrapped grid's potentially ambiguous manifold semantics.
"""
    authalic_shift(t::Helpers.AuthalicTransform) -> Float64

Upper bound in radians on the authalic-to-geodetic displacement,
`max_ξ |φ(ξ) - ξ|`, computed as `Σ|C'_j|` for
`φ(ξ) = ξ + Σ C'_j sin(2jξ)`. For WGS84 the bound is `2.2421e-3` radians.
"""
authalic_shift(t::Helpers.AuthalicTransform) = sum(abs ∘ Float64, t.inv)

"""
    authalic_stretch(t::Helpers.AuthalicTransform) -> Float64

Lipschitz constant of the authalic-to-geodetic warp `Φ` on the sphere:
`d(Φp, Φq) ≤ authalic_stretch(t) · d(p, q)` for every pair of points, with `d`
the great-circle distance. It equals `1 + Σ 2j|C'_j|` (`1 + 4.4886e-3` for
WGS84), obtained by bounding the north and east components of the warp
differential. It is exactly `1` when `e² = 0`.

Multiplicative rather than additive is what keeps it usable: inflating by
[`authalic_shift`](@ref) instead puts a 0.13° floor under every node extent,
which swamps the cell itself from about level 8 down.
"""
function authalic_stretch(t::Helpers.AuthalicTransform)
    total = 0.0
    for (j, c) in enumerate(t.inv)
        total += 2 * j * abs(Float64(c))
    end
    return 1 + total
end

"""
    geodetic_point(t::Helpers.AuthalicTransform, p) -> USPoint

Redraw authalic point `p` at its geodetic latitude, preserving longitude.
Poles are fixed, and an `e² = 0` transform returns bit-identical coordinates.
"""
@inline geodetic_point(t::Helpers.AuthalicTransform, p) =
    _warp_latitude(p, ξ -> Helpers.authalic_to_geodetic(t, ξ), t)

"""
    authalic_point(t::Helpers.AuthalicTransform, p) -> USPoint

The inverse of [`geodetic_point`](@ref): a point given in geodetic latitude,
re-drawn on the authalic sphere the grid is tessellated on. This is what
[`cellat`](@ref) runs on its argument before handing it to the wrapped grid.
"""
@inline authalic_point(t::Helpers.AuthalicTransform, p) =
    _warp_latitude(p, ξ -> Helpers.geodetic_to_authalic(t, ξ), t)

@inline function _warp_latitude(p, f, t::Helpers.AuthalicTransform)
    x, y, z = Float64(p[1]), Float64(p[2]), Float64(p[3])
    iszero(t.eccentricity_squared) && return USPoint(x, y, z)
    ρ = hypot(x, y)
    # A pole has no longitude to preserve and is a fixed point of both series.
    iszero(ρ) && return USPoint(x, y, z)
    s, c = sincos(Float64(f(atan(z, ρ))))
    scale = c / ρ
    return USPoint(x * scale, y * scale, s)
end

# ===========================================================================
# The grid wrapper
# ===========================================================================

"""
    AuthalicGrid(grid, ellipsoid = Helpers.WGS84_AUTHALIC) <: AbstractGrid

Read `grid` geometry at geodetic rather than authalic latitude. Cell ids,
indices, hierarchy, adjacency, and winding are unchanged. [`cellat`](@ref)
takes its query point in that same geodetic frame.

`ellipsoid` may be a [`Helpers.AuthalicTransform`](@ref) or a
`GeometryOpsCore.Manifold` (`Geodesic(; semimajor_axis, inv_flattening)` for an
ellipsoid, or `Spherical(; radius)` for the identity). `Planar` and
`AutoManifold` throw. The default is WGS84.

# Areas

[`cell_area`](@ref) returns the warped ring's unit-sphere area, not true
ellipsoidal area. For true ellipsoidal area, multiply the base grid's area by
`Helpers.authalic_radius(ellipsoid)^2`.

Wrapping an `AuthalicGrid` throws; use `parent(grid)` before changing ellipsoid.

A [`PartialGrid`](@ref) is also rejected. Wrap its system instead:

    PartialGrid(AuthalicSystem(sys), level, ids)

See also [`AuthalicSystem`](@ref).
"""
struct AuthalicGrid{G<:AbstractGrid,T<:AbstractFloat} <: AbstractGrid
    grid::G
    transform::Helpers.AuthalicTransform{T}

    function AuthalicGrid(grid::G, t::Helpers.AuthalicTransform{T}) where {G<:AbstractGrid,T}
        _check_wrappable(grid)
        return new{G,T}(grid, t)
    end
end

AuthalicGrid(grid::AbstractGrid) = AuthalicGrid(grid, Helpers.WGS84_AUTHALIC)

AuthalicGrid(grid::AbstractGrid, m::GOCore.Manifold) =
    AuthalicGrid(grid, Helpers.AuthalicTransform(m))

_check_wrappable(::AbstractGrid) = nothing

_check_wrappable(::AuthalicGrid) = throw(ArgumentError(
    "this grid is already an AuthalicGrid: its geometry is in geodetic latitude \
already, and warping it again would read those latitudes as authalic ones. Use \
`parent(grid)` to get the unwarped grid back if you meant to re-wrap it on a \
different ellipsoid."))

# The subset refusal is `Engine`'s: it needs the subset type, which lives there.

"""
    Base.parent(grid::AuthalicGrid) -> AbstractGrid
    Base.parent(sys::AuthalicSystem) -> AbstractHierarchicalGridSystem

Return the wrapped grid or system in authalic coordinates.
"""
Base.parent(grid::AuthalicGrid) = grid.grid

# --- the base grid interface ----------------------------------------------
#
# Identity forwards; geometry warps. Nothing here reorders or filters, so the
# `cellindex`/`localindex` bijection is the base grid's, unchanged.

ncells(grid::AuthalicGrid) = ncells(grid.grid)
cellindex(grid::AuthalicGrid, i::Int) = cellindex(grid.grid, i)
localindex(grid::AuthalicGrid, c::AbstractCellIndex) = localindex(grid.grid, c)
# Delegated rather than left to the `AbstractGrid` bridge: the bridge reads a
# grid's global index off its local one, which is only the same number when the
# wrapped grid is complete. This wrapper does not require that.
globalindex(grid::AuthalicGrid, c::AbstractCellIndex) = globalindex(grid.grid, c)
level(grid::AuthalicGrid) = level(grid.grid)

function system(grid::AuthalicGrid)
    sys = system(grid.grid)
    sys === nothing && return nothing
    return AuthalicSystem(sys, grid.transform)
end

# `map` rather than a comprehension: a system that returns a static vector of
# vertices keeps returning one, so a warped boundary allocates exactly where the
# unwarped one did.
cell_boundary(grid::AuthalicGrid, c::AbstractCellIndex) =
    map(p -> geodetic_point(grid.transform, p), cell_boundary(grid.grid, c))

cell_centroid(grid::AuthalicGrid, c::AbstractCellIndex) =
    geodetic_point(grid.transform, cell_centroid(grid.grid, c))

# The one input-side warp: the caller's point is in the geodetic frame this grid
# publishes, and the grid underneath speaks authalic.
cellat(grid::AuthalicGrid, p::GO.UnitSphericalPoint) =
    cellat(grid.grid, authalic_point(grid.transform, p))

# A monotone latitude warp preserves adjacency and winding.
neighbors(grid::AuthalicGrid, c::AbstractCellIndex, k::Integer=1;
    connectivity::Connectivity=Vertex()) =
    neighbors(grid.grid, c, k; connectivity)

# The latitude warp also preserves degree.
neighborcount(grid::AuthalicGrid, c::AbstractCellIndex;
    connectivity::Connectivity=Vertex()) =
    neighborcount(grid.grid, c; connectivity)

ring(grid::AuthalicGrid, c::AbstractCellIndex, k::Integer;
    connectivity::Connectivity=Vertex()) =
    ring(grid.grid, c, k; connectivity)

"""
    GeometryOpsCore.best_manifold(grid::AuthalicGrid) -> GO.Spherical

Return the unit sphere: the wrapper changes latitude, not coordinate scale.
Using the authalic radius here would apply `R_A²` twice to computed areas.
"""
# Dispatch-redundant with the `AbstractGrid` fallback, and kept anyway to pin
# these semantics against a future change to it. Not dead code.
GOCore.best_manifold(::AuthalicGrid) = GO.Spherical(; radius=1.0)

function Base.show(io::IO, grid::AuthalicGrid)
    print(io, "AuthalicGrid(", grid.grid, ", e² = ",
        grid.transform.eccentricity_squared, ")")
end

Base.show(io::IO, ::MIME"text/plain", grid::AuthalicGrid) = show(io, grid)

# ===========================================================================
# The system wrapper
# ===========================================================================

"""
    AuthalicSystem(sys, ellipsoid = Helpers.WGS84_AUTHALIC) <: AbstractHierarchicalGridSystem

Read every level grid's geometry at geodetic latitude. Hierarchical identities,
levels, ordering, and descendant ranges are forwarded unchanged.

`ellipsoid` is read exactly as [`AuthalicGrid`](@ref)'s is, and wrapping an
`AuthalicSystem` throws for the same reason.

[`node_extent`](@ref) is recomputed with the analytic
[`authalic_stretch`](@ref) bound because the warp is not an isometry.
"""
struct AuthalicSystem{S<:AbstractHierarchicalGridSystem,T<:AbstractFloat} <:
       AbstractHierarchicalGridSystem
    system::S
    transform::Helpers.AuthalicTransform{T}

    function AuthalicSystem(sys::S, t::Helpers.AuthalicTransform{T}) where {
            S<:AbstractHierarchicalGridSystem,T}
        sys isa AuthalicSystem && throw(ArgumentError(
            "this system is already an AuthalicSystem; its grids are in geodetic \
latitude already. Use `parent(sys)` to get the unwarped system back if you meant \
to re-wrap it on a different ellipsoid."))
        return new{S,T}(sys, t)
    end
end

AuthalicSystem(sys::AbstractHierarchicalGridSystem) =
    AuthalicSystem(sys, Helpers.WGS84_AUTHALIC)

AuthalicSystem(sys::AbstractHierarchicalGridSystem, m::GOCore.Manifold) =
    AuthalicSystem(sys, Helpers.AuthalicTransform(m))

Base.parent(sys::AuthalicSystem) = sys.system

# --- forwarded, because ids are not geometry -------------------------------

cellindextype(sys::AuthalicSystem) = cellindextype(sys.system)
cellindextypes(sys::AuthalicSystem) = cellindextypes(sys.system)
levels(sys::AuthalicSystem) = levels(sys.system)
maxlevel(sys::AuthalicSystem) = maxlevel(sys.system)
rootcells(sys::AuthalicSystem) = rootcells(sys.system)
children(sys::AuthalicSystem, c::AbstractCellIndex) = children(sys.system, c)
Base.parent(sys::AuthalicSystem, c::AbstractCellIndex) = Base.parent(sys.system, c)
has_sorted_subtrees(sys::AuthalicSystem) = has_sorted_subtrees(sys.system)
maxneighbors(sys::AuthalicSystem, connectivity::Connectivity) =
    maxneighbors(sys.system, connectivity)
ancestor(sys::AuthalicSystem, c::AbstractCellIndex, l::Integer) =
    ancestor(sys.system, c, l)
descendants(sys::AuthalicSystem, c::AbstractCellIndex, l::Integer) =
    descendants(sys.system, c, l)
descendant_range(sys::AuthalicSystem, c::AbstractCellIndex, l::Integer) =
    descendant_range(sys.system, c, l)
reindex(::Type{T}, sys::AuthalicSystem, c::AbstractCellIndex) where {T<:AbstractCellIndex} =
    reindex(T, sys.system, c)

# The warp preserves refinement geometry; `node_extent` is overridden below.
cap_inflation(sys::AuthalicSystem) = cap_inflation(sys.system)

levelgrid(sys::AuthalicSystem, l::Integer) =
    AuthalicGrid(levelgrid(sys.system, l), sys.transform)

"""
    node_extent(sys::AuthalicSystem, c) -> SphericalCap

Return the base node extent re-centred at its warped centre and inflated by
[`authalic_stretch`](@ref). For base cap `(c, r)` and warp `Φ`, the Lipschitz
bound gives

    d(Φc, Φv) ≤ L · d(c, v) ≤ L · r

so `(Φc, Lr)` covers every warped descendant vertex. If `Lr > π/2`, return the
full sphere because vertex containment no longer guarantees containment of the
great-circle boundary arcs.
"""
function node_extent(sys::AuthalicSystem, c::AbstractCellIndex)
    cap = node_extent(sys.system, c)
    radius = authalic_stretch(sys.transform) * Float64(cap.radius)
    # A cap wider than a quadrant is not convex, so containing the warped
    # vertices would no longer contain the arcs between them.
    radius > Float64(pi) / 2 && return full_sphere_cap()
    return SphericalCap(geodetic_point(sys.transform, cap.point), nextfloat(radius))
end

function Base.show(io::IO, sys::AuthalicSystem)
    print(io, "AuthalicSystem(", sys.system, ", e² = ",
        sys.transform.eccentricity_squared, ")")
end

Base.show(io::IO, ::MIME"text/plain", sys::AuthalicSystem) = show(io, sys)
