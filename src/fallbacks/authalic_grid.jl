# ---------------------------------------------------------------------------
# `AuthalicGrid` / `AuthalicSystem` — the ellipsoid wrapper
#
# Every DGGS in this package tiles a SPHERE, and an ellipsoidal DGGS is that
# same tessellation read on the *authalic* sphere: cell latitudes are authalic
# latitudes, and areas scaled by `R_A^2` are the true ellipsoidal areas (see
# `src/core/manifolds.jl`). Data that comes from anywhere else on Earth is in
# GEODETIC latitude. The two differ by up to 0.1283° for WGS84 — ~14.3 km along
# a meridian at ±45°, larger than a level-8 cell — so laying one over the other
# without a conversion is a silent misregistration.
#
# This pair is that conversion, applied at the geometry boundary and nowhere
# else:
#
#   * `AuthalicGrid(grid)` warps `cell_boundary`/`cell_centroid` OUTPUTS from
#     authalic to geodetic latitude and inverse-warps `cellat` INPUTS, so the
#     wrapped grid's geometry is in the frame a lon/lat dataset is in.
#   * `AuthalicSystem(sys)` does the same for the hierarchical surface, and owes
#     the covering law in the warped frame — which is the one thing here that is
#     not a forward, and the reason this file has a derivation in it.
#
# IDS ARE NOT TOUCHED. A cell is the same cell in both frames; only where its
# corners are drawn changes. So `ncells`, `cellindex`, `cellposition`,
# `parent`/`children`, `levels`, `rootcells` and every ordering fact forward
# verbatim, and a position stays a position in `1:ncells(grid)`.
# ---------------------------------------------------------------------------

# ===========================================================================
# The two bounds the wrapper is built on
# ===========================================================================

"""
    authalic_shift(t::Helpers.AuthalicTransform) -> Float64

An upper bound, in **radians**, on how far the warp moves a point:
`max_ξ |φ(ξ) − ξ|` for the authalic→geodetic series
`φ(ξ) = ξ + Σ_{j=1}^{6} C'_j sin(2jξ)` (Karney 2024 eq. A20).

Bounded term by term — `|sin(2jξ)| ≤ 1`, so the sum of the coefficient
magnitudes bounds the sum — which is sharp to three digits because `C'_1`
dominates the rest by three orders of magnitude and `sin 2ξ` does reach 1 (at
±45°, where the deviation is famously largest).

For WGS84 this is `2.2416e-3` rad = `0.12843°`, the ~14.3 km figure quoted in
`src/Helpers/authalic.jl`.

Not used by [`node_extent`](@ref) — [`authalic_stretch`](@ref) gives a strictly
better cap there — but it is the quantity "how wrong is it to skip the warp?"
is measured in, and the tests pin it against a swept maximum.
"""
authalic_shift(t::Helpers.AuthalicTransform) = sum(abs ∘ Float64, t.inv)

"""
    authalic_stretch(t::Helpers.AuthalicTransform) -> Float64

The **Lipschitz constant** of the authalic→geodetic warp `Φ` on the sphere:
`d(Φp, Φq) ≤ authalic_stretch(t) · d(p, q)` for every pair of points, with `d`
the great-circle distance.

`1 + Σ_{j=1}^{6} 2j |C'_j|`, which for WGS84 is `1 + 4.4886e-3`.

# Where that comes from

`Φ: (λ, ξ) ↦ (λ, φ(ξ))` with `φ(ξ) = ξ + δ(ξ)`, `δ(ξ) = Σ_j C'_j sin(2jξ)`.
Longitude is untouched, so in the orthonormal (east, north) frames at a point
and its image the differential is diagonal:

    dΦ = diag( cos φ(ξ) / cos ξ ,  φ'(ξ) )

and both entries are bounded by `L := 1 + Σ_j 2j |C'_j|`:

  - **north.** `φ'(ξ) = 1 + Σ_j 2j C'_j cos(2jξ) ≤ 1 + Σ_j 2j |C'_j| = L`.
  - **east.** `cos φ / cos ξ = cos δ − tan ξ · sin δ ≤ 1 + |tan ξ| · |δ(ξ)|`.
    Writing `θ = π/2 − ξ` gives `sin(2jξ) = ±sin(2jθ)`, and `|sin(2jθ)| ≤
    2j |sin θ| = 2j cos ξ` (the standard `|sin nθ| ≤ n|sin θ|`), so
    `|δ(ξ)| ≤ cos ξ · Σ_j 2j |C'_j|` and therefore
    `|tan ξ| · |δ(ξ)| ≤ |sin ξ| · Σ_j 2j |C'_j| ≤ L − 1`.

A map whose differential has operator norm ≤ `L` everywhere sends the geodesic
joining `p` and `q` to a path of length ≤ `L · d(p, q)` joining `Φp` and `Φq`,
which is the stated bound.

Both bounds are **analytic**, from the published series alone: no sampling, and
no dependence on the ellipsoid beyond its coefficients. The degenerate `e² = 0`
transform has every coefficient zero and gives exactly `1`, so a spherical grid
pays nothing.

Multiplicative rather than additive is what keeps the wrapper usable: the
alternative bound — every point moves by at most [`authalic_shift`](@ref), so a
cap of radius `r` about `c` sends its contents into a cap of radius `r + 2Δ`
about `Φc` — puts a 0.13° FLOOR under every node extent, which swamps the cell
itself from about level 8 down and destroys pruning there. `L·r` shrinks with
`r` and is tighter than `r + 2Δ` for every radius under ~1 rad, which is every
node extent any system in this package produces.
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

`p` read on the authalic sphere, re-drawn at its **geodetic** latitude:
longitude is kept, latitude `ξ = atan(z, √(x²+y²))` is replaced by
`φ = authalic_to_geodetic(t, ξ)`.

The poles are fixed points and are returned unchanged rather than rebuilt, both
because `atan` there is `atan(±1, 0)` and because the rebuilt point would
otherwise carry a meaningless longitude. The degenerate `e² = 0` transform is
the exact identity — every series coefficient is zero — and is short-circuited
so that a spherical grid's coordinates come back bit-identical rather than
merely equal.
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

`grid` with its geometry read on an ellipsoid: the same cells, with the same
ids and the same positions, drawn at **geodetic** latitude instead of the
authalic latitude the tessellation is computed in.

`ellipsoid` may be a [`Helpers.AuthalicTransform`](@ref) or a
`GeometryOpsCore.Manifold` (`Geodesic(; semimajor_axis, inv_flattening)` for an
ellipsoid, `Spherical(; radius)` for the degenerate identity); `Planar` and
`AutoManifold` throw, because a DGGS's ellipsoid is part of its reference system
and guessing one is the misregistration this type exists to prevent. The default
is WGS84 — the datum an unqualified terrestrial "lon/lat" means — and it is a
default rather than a guess: it is named, and it is the *shape*, not the kind of
frame.

# What changes and what does not

| | |
|---|---|
| [`cell_boundary`](@ref), [`cell_centroid`](@ref) | warped **outputs** |
| [`cellat`](@ref) (both methods) | warped **input**, then forwarded |
| [`ncells`](@ref), [`cellindex`](@ref), [`cellposition`](@ref), [`level`](@ref) | forwarded verbatim |
| [`neighbors`](@ref), [`ring`](@ref) | forwarded: the warp is a homeomorphism that maps shared vertices to shared vertices, so adjacency and winding are unchanged |
| [`system`](@ref) | the base system, wrapped in [`AuthalicSystem`](@ref) |

Positions still run `1:ncells(grid)` and mean exactly what they meant in the
base grid, so a data array laid out against one is laid out against the other.

# Areas

[`cell_area`](@ref) on a wrapped grid is the area of the **warped** ring, in
steradians of the coordinate sphere — *not* the cell's true ellipsoidal area,
and not equal across cells even for an equal-area system. That is not a defect:
it is the area consistent with the geometry this grid publishes, which is what
makes `ConservativeRegridding`'s overlap fractions add up on it. For the true
ellipsoidal area take the **base** grid's `cell_area` and multiply by
`Helpers.authalic_radius(ellipsoid)^2`.

# Wrapping a wrapped grid

Throws an `ArgumentError`. Composing two latitude warps is not a coordinate
transformation anyone means — the second one would be fed geodetic latitudes as
if they were authalic — and silently collapsing to one of the two transforms
would quietly discard the other, which is how a WGS84 grid ends up being read on
GRS80. Unwrap with `parent(grid)` and re-wrap if that is what was meant.

A [`PartialGrid`](@ref) is rejected for a different reason: a subset is a
property of the id set and the warp is a property of the system, so the two
compose the other way round —

    PartialGrid(AuthalicSystem(sys), level, ids)

— which is not merely tidier but the only spelling the tree cursor reads
correctly, since it dispatches its position windows on `PartialGrid` itself.

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

_check_wrappable(::PartialGrid) = throw(ArgumentError(
    "wrap the SYSTEM, not the subset: `PartialGrid(AuthalicSystem(sys), level, ids)`. \
A subset is a property of the id set and the warp is a property of the system, and \
only that order keeps the tree cursor's position windows correct."))

"""
    Base.parent(grid::AuthalicGrid) -> AbstractGrid
    Base.parent(sys::AuthalicSystem) -> AbstractHierarchicalGridSystem

The wrapped grid or system, with its geometry back on the authalic sphere. The
inverse of the constructor, and the way to re-wrap on a different ellipsoid
without composing two warps.
"""
Base.parent(grid::AuthalicGrid) = grid.grid

# --- the base grid interface ----------------------------------------------
#
# Identity forwards; geometry warps. Nothing here reorders or filters, so the
# `cellindex`/`cellposition` bijection is the base grid's, unchanged.

ncells(grid::AuthalicGrid) = ncells(grid.grid)
cellindex(grid::AuthalicGrid, i::Int) = cellindex(grid.grid, i)
cellposition(grid::AuthalicGrid, c::AbstractCellIndex) = cellposition(grid.grid, c)
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

# Adjacency is combinatorial. The warp is a homeomorphism of the sphere that
# fixes longitudes and is strictly monotone in latitude, so cells that share a
# vertex still share it and a counter-clockwise ring is still counter-clockwise:
# forwarding is not an approximation, it is the same answer without the
# geometric search.
neighbors(grid::AuthalicGrid, c::AbstractCellIndex, k::Integer=1;
    connectivity::Connectivity=Vertex()) =
    neighbors(grid.grid, c, k; connectivity)

ring(grid::AuthalicGrid, c::AbstractCellIndex, k::Integer;
    connectivity::Connectivity=Vertex()) =
    ring(grid.grid, c, k; connectivity)

"""
    GeometryOpsCore.best_manifold(grid::AuthalicGrid) -> GO.Spherical

The **unit** sphere, exactly as for every other grid here — the wrapper changes
coordinates, not the manifold, and its coordinates are still unit-sphere
`(x, y, z)`.

Emphatically *not* [`authalic_sphere`](@ref): that manifold's radius is the one
on which the *unwarped* grid's areas are true ellipsoidal areas, and a wrapped
grid's geometry is precisely the geometry that is no longer equal-area. Handing
`ConservativeRegridding` an `R_A`-radius manifold for coordinates that are unit
vectors would put a factor of `R_A²` in every area twice over.
"""
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

`sys` with its geometry read on an ellipsoid: every level grid is the
corresponding [`AuthalicGrid`](@ref), and every hierarchical fact — ids, levels,
[`rootcells`](@ref), [`parent`](@ref)/[`children`](@ref), the canonical order
and therefore [`descendant_range`](@ref) — is the base system's, forwarded
unchanged.

`ellipsoid` is read exactly as [`AuthalicGrid`](@ref)'s is, and wrapping an
`AuthalicSystem` throws for the same reason.

The one method that is not a forward is [`node_extent`](@ref), because the warp
is not an isometry and a cap is not a cap after it. See its docstring for the
bound, which is derived from the transform's own series
([`authalic_stretch`](@ref)) rather than measured.

    julia> sys = AuthalicSystem(HEALPixSystem());

    julia> grid = levelgrid(sys, 4);        # an AuthalicGrid

    julia> system(grid) === sys
    true
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
max_level(sys::AuthalicSystem) = max_level(sys.system)
rootcells(sys::AuthalicSystem) = rootcells(sys.system)
children(sys::AuthalicSystem, c::AbstractCellIndex) = children(sys.system, c)
Base.parent(sys::AuthalicSystem, c::AbstractCellIndex) = Base.parent(sys.system, c)
has_sorted_subtrees(sys::AuthalicSystem) = has_sorted_subtrees(sys.system)
max_neighbors(sys::AuthalicSystem, connectivity::Connectivity) =
    max_neighbors(sys.system, connectivity)
ancestor(sys::AuthalicSystem, c::AbstractCellIndex, l::Integer) =
    ancestor(sys.system, c, l)
descendants(sys::AuthalicSystem, c::AbstractCellIndex, l::Integer) =
    descendants(sys.system, c, l)
descendant_range(sys::AuthalicSystem, c::AbstractCellIndex, l::Integer) =
    descendant_range(sys.system, c, l)
subtree_border(sys::AuthalicSystem, c::AbstractCellIndex, l::Integer;
    connectivity::Connectivity=Vertex()) =
    subtree_border(sys.system, c, l; connectivity)
subtree_interior(sys::AuthalicSystem, c::AbstractCellIndex, l::Integer;
    connectivity::Connectivity=Vertex()) =
    subtree_interior(sys.system, c, l; connectivity)
reindex(::Type{T}, sys::AuthalicSystem, c::AbstractCellIndex) where {T<:AbstractCellIndex} =
    reindex(T, sys.system, c)

# The base system's headroom factor, forwarded so that the trait keeps saying
# something true about the refinement geometry — which the warp does not change,
# since it maps a cell's descendants to that cell's descendants. Nothing here
# reads it: `node_extent` below is an override.
cap_inflation(sys::AuthalicSystem) = cap_inflation(sys.system)

levelgrid(sys::AuthalicSystem, l::Integer) =
    AuthalicGrid(levelgrid(sys.system, l), sys.transform)

"""
    node_extent(sys::AuthalicSystem, c) -> SphericalCap

The base system's node extent, **inflated by [`authalic_stretch`](@ref)** and
re-centred on the warp of its centre.

# Why an inflation is owed at all

A spherical cap is not a spherical cap after the authalic latitude shift. The
warp `Φ` moves points along their meridians by up to 0.1283° (WGS84), and it
moves them by *different* amounts at different latitudes — that is the entire
content of an auxiliary latitude — so the image of a cap is an egg, not a cap,
and the base system's own extent does not contain it. Under-covering a node
extent silently drops cells from every pruned traversal in the package (see the
covering law in the interface docstring), so this is not a place to be
approximately right.

# The bound

Let `C = (c, r)` be `node_extent(sys.system, x)`. By the base system's covering
law every boundary vertex `v` of every descendant of `x`, at every depth,
satisfies `d(c, v) ≤ r`. `Φ` is Lipschitz with constant
`L = authalic_stretch(t)` ([derived there](@ref authalic_stretch) from the
series coefficients, no sampling), so

    d(Φc, Φv) ≤ L · d(c, v) ≤ L · r

for every one of them: the cap `(Φc, L·r)` contains every warped vertex.

The wrapped grid's boundary is the great-arc polygon through those warped
vertices — *not* the image of the base boundary, which is not a great-arc
polygon at all — so vertices alone are the whole story only while the cap is
geodesically convex. It is, at radius ≤ 90°, and a convex region containing two
points contains the arc between them. Past 90° that argument is gone and the
honest answer is the full sphere, which is what is returned there. (No system in
this package comes close: the widest node extent any of them produces is
HEALPix's 0.907 rad at level 0, which `L` moves to 0.911.)

For WGS84 `L − 1 = 4.4886e-3`, so this costs a node extent 0.45% of its radius —
against the 0.1283° *additive* floor the naive "inflate by the maximum
deviation" bound would put under every level, which at level 10 is thirty times
the cell.
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
