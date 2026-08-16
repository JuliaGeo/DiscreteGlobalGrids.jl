# The two-level hierarchy: 64 800 tiles at level 0, one tile's raster at level 1.
# Every relation here is arithmetic on the dense 0-based ordinals of `bands.jl`.

# `levelgrid(CopernicusDEMSystem(...), l)` is the package's `HierarchicalLevelGrid`:
# all 64 800 tiles, or every pixel on Earth, in ordinal order. The methods this system
# takes over from the generic grid hang off the alias, and the primitives they forward
# to are the `(sys, ...)` methods below.
const LevelGrid{N} = DGG.HierarchicalLevelGrid{CopernicusDEMSystem{N}}

# ===========================================================================
# System interface
# ===========================================================================

DGG.cellindextype(::CopernicusDEMSystem) = DGG.LevelIndex
DGG.levels(::CopernicusDEMSystem) = 0:1
DGG.has_sorted_subtrees(::CopernicusDEMSystem) = true

"""
    max_neighbors(CopernicusDEMSystem{N}(), connectivity) -> Int

`360 * ncols(pole row) + 2`, which is `36N + 2`, under `Vertex()`; `6` under `Edge()`.
Both bound every cell of both levels, and both are attained.

A pole row sets the `Vertex()` bound, not the lattice interior. Raster row 0 of a
`lat_s = 89` tile and row `N - 1` of a `lat_s = -90` tile are the spherical triangles
[`cell_boundary`](@ref) emits, and their apex is the one exact ±90 vertex every cell of
the row shares — so each such pixel touches the whole pole ring of `360 * ncols` pixels,
plus the three pixels below it. Everything else is far under that: eight for a pixel
interior to its tile, for a tile edge within a band and for the wide side of a band
boundary, and at most `360 + 2` for a level-0 pole tile.

A band boundary sets the `Edge()` bound, which is why it exceeds the von Neumann four.
The two sides of a boundary parallel carry different column counts, so a cell on the
wide side faces the narrow side across a segment that up to three of its cells divide:
with the reduced ratio `p:q` written wide side first, every one the band table produces
(3:2, 4:3, 3:2, 5:3, 2:1) leaves `ceil(p/q) <= 2` narrow-side breakpoints strictly
inside that segment, hence at most three overlaps of positive length. Two laterals and
the one cell facing the other way bring the total to six. A level-0 tile reaches five
the same way. Under `Edge()` a shared apex is a point rather than a segment, so a pole
cell has three and never approaches the `Vertex()` figure.

Both bounds are attained, not merely respected: `"max_neighbors is attained"` in
`test/systems/CopernicusDEM/runtests.jl` asserts equality at the pole rows and at the
±85 boundary of all three lattices, and sweeps the boundaries, the antimeridian, the
tile edges and both levels for anything larger.
"""
DGG.max_neighbors(sys::CopernicusDEMSystem, ::DGG.Vertex) =
    360 * Int(max(ncols(sys, 0), ncols(sys, NROWS - 1))) + 2
DGG.max_neighbors(::CopernicusDEMSystem, ::DGG.Edge) = 6

"Lazy 0-based ids at one level; `rootcells` and `children` are windows over it."
struct IdRange <: AbstractVector{DGG.LevelIndex}
    level::Int32
    first::Int64      # 0-based index of element 1
    n::Int
end

Base.size(v::IdRange) = (v.n,)
Base.IndexStyle(::Type{IdRange}) = Base.IndexLinear()
Base.@propagate_inbounds function Base.getindex(v::IdRange, i::Int)
    @boundscheck checkbounds(v, i)
    return DGG.LevelIndex(v.level, v.first + i - 1)
end

"""
    rootcells(CopernicusDEMSystem(...))

All 64 800 level-0 tiles, `LevelIndex(0, 0:64799)`, ascending, as a **lazy** vector.

Lazy because `PartialGrid` reads `first(rootcells(sys))` on every construction
(`src/fallbacks/partial_grid.jl`), and a materialised 64 800-element vector would
allocate about a megabyte per chunk built. The count is far above the small, cheap
collection the contract has in mind (12 for HEALPix), and it costs the generic tree
descent one cap evaluation per tile at the synthetic root — the price of a grid whose
base tessellation is the 1° graticule.
"""
DGG.rootcells(::CopernicusDEMSystem) = IdRange(Int32(0), Int64(0), NTILES)

# ===========================================================================
# The level grid: size, and positions <-> ids
# ===========================================================================

DGG.ncells(sys::CopernicusDEMSystem, l::Integer) =
    Int(l) == 0 ? NTILES :
    Int(l) == 1 ? Int(tables(sys).rowbase[end]) :
    throw(ArgumentError("level $l is outside $(DGG.levels(sys))"))

# The grid bounds-checks `i`, so this is the bijection and nothing else.
DGG.cellindex(::CopernicusDEMSystem, l::Integer, i::Int) = DGG.LevelIndex(l, i - 1)

"""
    cellposition(CopernicusDEMSystem(...), c) -> Union{Int,Nothing}

Closed form: `index + 1` for an in-range id, and `nothing` for one no cell has.
Never throws — a miss is an answer — so this is the one decoder that does not go
through the id guard.
"""
function DGG.cellposition(sys::CopernicusDEMSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    (l == 0 || l == 1) || return nothing
    0 <= c.index < DGG.ncells(sys, l) || return nothing
    return Int(c.index + 1)
end

# ===========================================================================
# The hierarchy
# ===========================================================================

"""
    parent(CopernicusDEMSystem(...), pixel) -> LevelIndex

The tile the pixel belongs to. Throws an `ArgumentError` on a level-0 cell, which is
a tile and is a root.
"""
function Base.parent(sys::CopernicusDEMSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    l == 1 || throw(ArgumentError(l == 0 ?
        "level-0 Copernicus DEM cell $c is a tile and has no parent" :
        "level $l is outside $(DGG.levels(sys))"))
    r, q, _, _ = decode(sys, c)
    return DGG.LevelIndex(0, tileordinal(r, q))
end

"""
    children(CopernicusDEMSystem(...), tile)

The tile's pixels in raster order — north row first, west to east — as a **lazy**
vector of `ncols * N` ids. Lazy because a GLO-30 tile in the 1x band has 12 960 000
of them and the tree cursor iterates without collecting.
"""
function DGG.children(sys::CopernicusDEMSystem{N}, c::DGG.LevelIndex) where {N}
    l = DGG.level(c)
    l == 0 || throw(ArgumentError(l == 1 ?
        "level-1 Copernicus DEM cell $c is a pixel, at max_level 1, and has no children" :
        "level $l is outside $(DGG.levels(sys))"))
    r, q, _, _ = decode(sys, c)
    return IdRange(Int32(1), tilebase(sys, r, q), Int(ncols(sys, r)) * N)
end

"""
    ancestor(CopernicusDEMSystem(...), c, l) -> LevelIndex

`c` itself at `l == level(c)`, and [`parent`](@ref) at `l == 0` for a pixel. Two
levels leave no third case.
"""
function DGG.ancestor(sys::CopernicusDEMSystem, c::DGG.LevelIndex, l::Integer)
    target = Int(l)
    lc = DGG.level(c)
    target <= lc || throw(ArgumentError(
        "ancestor level $target is deeper than the cell's own level $lc"))
    target >= 0 || throw(ArgumentError(
        "ancestor level $target is above the root level 0"))
    target == lc && return c
    return Base.parent(sys, c)
end

"""
    descendant_range(CopernicusDEMSystem(...), tile, 1) -> UnitRange{Int}

The tile's contiguous window of level-1 positions: `tilebase + 1 : tilebase + ncols*N`.

Exact and hole-free in both directions, which is what `has_sorted_subtrees == true`
asserts. It holds because the level-1 order is tile-major and raster-order within a
tile, so a tile's pixels are consecutive by construction and no id in the window
belongs to another tile.

`l == level(c)` is the cell's own one-element position range; `l < level(c)` throws an
`ArgumentError`.
"""
function DGG.descendant_range(sys::CopernicusDEMSystem{N}, c::DGG.LevelIndex,
        l::Integer) where {N}
    target = Int(l)
    lc = DGG.level(c)
    index = _checked_index(sys, c)          # also rejects a level outside 0:1
    target >= lc || throw(ArgumentError(
        "descendant level $target is above the cell's own level $lc"))
    target <= DGG.max_level(sys) || throw(ArgumentError(
        "descendant level $target is past max_level $(DGG.max_level(sys))"))
    if target == lc
        pos = Int(index + 1)
        return pos:pos
    end
    r, q, _, _ = decode(sys, c)
    base = tilebase(sys, r, q)
    return Int(base + 1):Int(base + ncols(sys, r) * Int64(N))
end

"""
    descendants(CopernicusDEMSystem(...), c, l)

Every level-`l` descendant of `c`, ascending, as the same **lazy** vector
[`children`](@ref) returns — the ids are [`descendant_range`](@ref) read off as
consecutive ordinals, and a GLO-30 tile's 12 960 000 of them are not worth
materialising.

!!! warning "This diverges from the interface"
    The interface docstring for [`descendants`](@ref) says the call *materializes*
    `O(subtree)` ids, and every other system hands back a freshly allocated `Vector` the
    caller owns. This method returns a lazy, read-only `AbstractVector` instead, which
    computes each id on indexing. Reading it is complete — `length`, `getindex`,
    iteration, `collect` — but nothing that writes to it works: no `setindex!`, no
    `push!`, no `sort!`, and no passing it to an API that mutates its argument. `collect`
    it first for any of those. The divergence buys the memory: one GLO-30 tile has
    12 960 000 level-1 descendants, at 16 bytes apiece in a `Vector{LevelIndex}`.
"""
function DGG.descendants(sys::CopernicusDEMSystem, c::DGG.LevelIndex, l::Integer)
    r = DGG.descendant_range(sys, c, l)     # validates `l` both ways
    return IdRange(Int32(l), Int64(first(r) - 1), length(r))
end

# The grid-level id guard, alongside `bands.jl`'s system-level one: it additionally
# pins the cell to the level of the grid it was handed to.
@inline function _checked_index(g::LevelGrid, c::DGG.LevelIndex)
    DGG.level(c) == g.level || throw(ArgumentError(
        "cell $c is at level $(DGG.level(c)), not the grid's level $(g.level)"))
    return _checked_index(g.system, c)
end

# ===========================================================================
# Geometry
# ===========================================================================

"""
    cell_box(sys, c) -> (west, east, south, north)

The cell's longitude/latitude box in **degrees**, closed on every side — the region
the DEM post samples. Every other geometric method here is a function of it.

# The registration, and why the box is offset

The AWS COGs are `RasterPixelIsPoint` (`AREA_OR_POINT=Point`): pixel CENTRES sit on
the integer-degree lattice. For a tile labelled `(lat_s, lon_w)` with `Δlat = 1/N` and
`Δlon = 1/ncols` degrees, centres are

    lon = lon_w + i*Δlon    i = 0 : ncols-1
    lat = lat_n - j*Δlat    j = 0 : N-1,  lat_n = lat_s + 1

so the first centre is exactly `(lon_w, lat_n)` and the GDAL geotransform origin is a
HALF PIXEL outside it: `(lon_w - Δlon/2, lat_n + Δlat/2)`. AWS deleted the east column
and south row of the original 3601-post DGED tile, so adjacent COG tiles abut with no
overlap and the pixel-centre lattice partitions the globe cleanly. That is the
convention this system indexes.

A pixel's box is its centre ± half a pixel. A tile's box is the union of its pixels'
boxes, which comes out as the nominal 1°x1° box shifted **half a pixel west and half a
pixel north**:

    lon in [lon_w - Δlon/2, lon_w + 1 - Δlon/2]
    lat in [lat_s + Δlat/2, lat_s + 1 + Δlat/2]

Tiles therefore abut exactly in latitude (Δlat is global) and in longitude within a
band (Δlon depends only on latitude), and the 64 800 tiles tile the sphere.

# The poles

Two corrections make that tiling exact rather than nearly so:

  - The top row of the `lat_s = 89` tiles would reach `90 + Δlat/2`. Its north edge is
    **clamped to +90**.
  - The bottom row of the `lat_s = -90` tiles stops at `-90 + Δlat/2`, leaving a
    half-pixel gap ring. Its south edge is **extended to -90**, making that row one and
    a half pixels tall.

Both give a cell one degenerate edge — every longitude at latitude ±90 is the same
point — so those cells are spherical TRIANGLES, and [`cell_boundary`](@ref) emits three
vertices for them, never four with a duplicate.
"""
function cell_box(sys::CopernicusDEMSystem{N}, c::DGG.LevelIndex) where {N}
    r, q, j, i = decode(sys, c)
    lat_s = _lat_s(r)
    lon_w = _lon_w(q)
    nc = ncols(sys, r)
    # A tile is its pixel grid's outer frame, so both levels are the same expression
    # over a column and a row interval.
    i_w, i_e, j_n, j_s = DGG.level(c) == 0 ? (0, Int(nc), 0, N) : (i, i + 1, j, j + 1)
    # `k / nc` and `k / N`, never `k * Δlon`: at the far edge `nc / nc` is exactly
    # `1.0`, so a cell's east edge is the same `Float64` as its neighbour's west edge
    # and the quads share a geodesic rather than nearly sharing one.
    half_dlon = (1 / nc) / 2
    half_dlat = (1 / N) / 2
    west = (lon_w + i_w / nc) - half_dlon
    east = (lon_w + i_e / nc) - half_dlon
    lat_n = lat_s + 1
    north = (lat_n - j_n / N) + half_dlat
    south = (lat_n - j_s / N) + half_dlat
    j_n == 0 && lat_s == 89 && (north = 90.0)
    j_s == N && lat_s == -90 && (south = -90.0)
    return (west, east, south, north)
end

# The exact ±90 vertices, as shared literals. `UnitSphereFromGeographic()((lon, ±90))`
# is not a substitute. It goes through `sincosd`, so `cosd(±90)` is exactly `0.0` and the
# image is exactly `(±0.0, ±0.0, ±1.0)`: `x = sind(0 or 180) * cosd(lon)` carries a sign
# bit wherever `cosd(lon) < 0`, and `y` wherever `sind(lon) < 0`. The separation from
# this point is therefore ZERO, not 1e-16 — no orientation or convexity test can tell
# the two apart, because `-0.0 == 0.0`.
#
# The hazard is identity, not magnitude. A signed zero is `===`-distinct and
# `isequal`-distinct, so building each pole cell's apex from its own longitude would
# hand out four bit patterns for one point: `Set`, `unique`, `Dict` keys and every
# `===`-based dedup or identity check would count up to four distinct pole vertices
# while `==` insists they are the same one. `atand` reads the sign bits too, which is
# how the same pole would land in different tiles — see [`cellat`](@ref).
const NORTH_POLE = GO.UnitSphericalPoint(0.0, 0.0, 1.0)
const SOUTH_POLE = GO.UnitSphericalPoint(0.0, 0.0, -1.0)

"""
    cell_boundary(grid, c) -> Vector{UnitSphericalPoint}

The cell's box as a **plain 4-corner great-circle quadrilateral**, implicitly closed and
counter-clockwise seen from outside the sphere, in the order
`(W,S) -> (E,S) -> (E,N) -> (W,N)`.

# Why there is no densification

A parallel is a small circle, not a great circle, so the ring's north and south edges
bow poleward of the true box edge by about `(Δλ²/8)·sin φ·cos φ` — 1.9e-5 rad (120 m)
for a 1° tile, about 10 μm for a 1-arcsec pixel. That small-angle expression is an
estimate rather than a bound, and it errs low.

The two bows nearly cancel, so ring and box differ in area only slightly. The
`"ring vs box"` testset in `test/systems/CopernicusDEM/runtests.jl` sweeps that
relative gap over every band and both pole rows, logs it as `worst_tile` /
`worst_pixel` / `worst_pole_pixel`, and bounds each: below 1e-8 for pixels outside the
±90 tile rows, below 1e-4 for tiles and for the pixels inside those rows, which are
half-pixel slivers and pole triangles.

Densifying would close that gap and cost this system its convexity: a densified
poleward edge reads as a chain of REFLEX vertices, which is exactly why HEALPix,
ISEA4R and A5 are non-conservative as regridding DESTINATIONS
(`test/systems/crosssystem/regridding_conservation.jl`, header). Undensified lat/lon
quads are convex, so this system is exact in both directions. [`cell_area`](@ref)
returns the exact box area rather than this ring's, and [`node_extent`](@ref) pads for
the bow analytically rather than by sampling.

Adjacent cells share their corner POINTS bit-identically within a band, so the shared
edge is literally the same geodesic and the quads tile the sphere with no gaps. Across
a band boundary (latitude 50/60/70/80/85) the two sides cut the shared parallel at
different longitudes, so their bows differ and the quads leave slivers of width
`|Δλ_below² − Δλ_above²|/8 · sin φ cos φ` — by that expression, 1.8e-12 rad (11.5 μm)
at latitude 50 in GLO-30, about 4e-7 of that pixel's own height, rising to 1.9e-11 rad
(122 μm) at latitude 85. The box tessellation is exact there; the quad tessellation is
exact to that order.

# Poles

A cell touching ±90 has two coincident corners and is emitted as a TRIANGLE, on the
same cyclic order with the duplicate dropped. The pole vertex is the exact literal
`UnitSphericalPoint(0.0, 0.0, ±1.0)`, not the image of `(lon, ±90)` under
`UnitSphereFromGeographic`, which would differ in x and y per longitude and turn a
triangle into a degenerate quad.
"""
function DGG.cell_boundary(sys::CopernicusDEMSystem, c::DGG.LevelIndex)
    west, east, south, north = cell_box(sys, c)
    north == 90.0 && return [TO_SPHERE((west, south)), TO_SPHERE((east, south)),
                             NORTH_POLE]
    south == -90.0 && return [SOUTH_POLE, TO_SPHERE((east, north)),
                              TO_SPHERE((west, north))]
    return [TO_SPHERE((west, south)), TO_SPHERE((east, south)),
            TO_SPHERE((east, north)), TO_SPHERE((west, north))]
end

"""
    cell_centroid(grid, c) -> UnitSphericalPoint

The midpoint of the cell's box, which is the DEM post itself for a pixel.

Strictly inside the published quad, as the contract requires: the quad's poleward edge
bows off the box's parallel by about `(Δλ²/8)·sin φ·cos φ` — 1.9e-5 rad for a 1° tile,
against a half-height of 8.7e-3 rad — so the midpoint clears it by a factor of 458.
"""
function DGG.cell_centroid(sys::CopernicusDEMSystem, c::DGG.LevelIndex)
    west, east, south, north = cell_box(sys, c)
    return TO_SPHERE(((west + east) / 2, (south + north) / 2))
end

"""
    cell_area(grid, c) -> Float64

The exact solid angle of the cell's lat/lon BOX, `Δλ · (sin φ_N − sin φ_S)` steradians,
in O(1).

The box, not the published ring: the ring is the box's undensified geodesic quad, and
the two differ in area by the bow [`cell_boundary`](@ref) describes and the
`"ring vs box"` testset measures.

Both regions tile the sphere, so both sum to 4π — but this closed form does **not**
telescope exactly: not every tile has `east - west === 1.0`, so consecutive terms do
not cancel to the last bit. The sum is accurate because those errors are tiny and
because pairwise summation stops them accumulating, not because anything cancels.
**Materialise** the 64 800 tile areas into a `Vector` so Julia's pairwise `sum`
reduces them; a generator argument — or any other sequential accumulation — lands
orders further from 4π and fails the `rtol = 1e-14` that
`"the boxes partition the sphere"` in `test/systems/CopernicusDEM/runtests.jl`
asserts. Do not loosen that tolerance to accommodate a generator.

`ConservativeRegridding` never reads this — it measures the ring itself
(`ConservativeRegridding/src/regridder/regridder.jl:102-103`) — so no conservation
assertion anywhere depends on the choice.
"""
function DGG.cell_area(g::LevelGrid, c::DGG.LevelIndex)
    _checked_index(g, c)
    west, east, south, north = cell_box(g.system, c)
    return deg2rad(east - west) * (sind(north) - sind(south))
end

"""
    cell_extent(grid, c) -> Extents.Extent{(:X, :Y)}

The cell's [`cell_box`](@ref), verbatim. Same rationale as [`cell_area`](@ref): the cell
IS the box, and the published ring is a slightly different region.

The generic fallback (`src/fallbacks/geometry.jl`) derives the extent from the ring, and
that is wrong here in a way that matters. The ring's poleward edge bows past the box's
parallel by the `(Δλ²/8)·sin φ·cos φ` [`cell_boundary`](@ref) quantifies — 0.00107° for
a 1° tile at mid latitudes, four GLO-30 pixel rows — so the derived `Y` upper bound
overshoots the box's north edge by that much and vertically adjacent tiles report
OVERLAPPING extents. Extents of a partition must abut, because callers use them to
decide which cells a query region can touch.
"""
function DGG.cell_extent(g::LevelGrid, c::DGG.LevelIndex)
    _checked_index(g, c)
    west, east, south, north = cell_box(g.system, c)
    return DGG.Extents.Extent(X = (west, east), Y = (south, north))
end

"""
    node_extent(CopernicusDEMSystem(...), c) -> SphericalCap

The cap centred on the cell centre, with radius the largest corner distance plus an
analytic pad for the great-circle bow, and one outward ULP.

Sound without sampling and without densification, in three steps.

 1. **The corners bound the box.** The farthest point of a lat/lon box from its midpoint
    centroid `(λ_c, φ_c)` is a corner, and all four corners are measured. Fix `φ`:
    `cos d = sin φ_c sin φ + cos φ_c cos φ cos(λ − λ_c)` falls as `|λ − λ_c|` grows, so
    the distance peaks at `λ ∈ {W, E}`. Now fix that `λ`: as a function of `φ` the same
    expression is `A sin φ + B cos φ` with `B = cos φ_c cos Δλ ≥ 0`, a sinusoid whose
    trough sits outside `[−90°, 90°]`; it is therefore unimodal on `[φ_S, φ_N]` and its
    minimum — the maximal distance — is at `φ ∈ {S, N}`. So `rmax` bounds every point of
    the box, hence every box vertex of every descendant.
 2. **The cell's own ring is already inside that cap.** Its vertices are the box corners,
    each within `rmax`, and a spherical cap is convex — it contains the geodesic between
    any two points it contains — so it contains the whole ring, bow and all.
 3. **The pad covers the descendants' bows.** The only geometry that leaves a descendant's
    box is that descendant's own ring, bowing poleward by about `(Δλ_child²/8)·sin φ·cos φ`.
    With `|sin φ cos φ| ≤ 1/2` the pad `Δλ²/16` — 1.9e-5 rad (120 m) for a 1° tile —
    already covers a bow at the CELL's own span; a child's span is at most `Δλ/120` on the
    shipped lattices (GLO-90's 10x band above latitude 85 is the coarsest, at 120 columns),
    so its bow is at least 14 400x inside the pad. The pad is belt-and-braces; steps 1 and
    2 are the proof.

Radii are far below 90°, so `require_convex_extents = true` holds and the harness's
vertex-sampling proxy is sound.
"""
function DGG.node_extent(sys::CopernicusDEMSystem, c::DGG.LevelIndex)
    centre = DGG.cell_centroid(sys, c)
    ring = DGG.cell_boundary(sys, c)
    rmax = maximum(US.spherical_distance(centre, p) for p in ring)
    west, east, _, _ = cell_box(sys, c)
    pad = deg2rad(east - west)^2 / 16
    return SphericalCap(centre, nextfloat(min(Float64(π), rmax + pad)))
end

# ===========================================================================
# Location
# ===========================================================================

"""
    cellat(grid, p::UnitSphericalPoint) -> LevelIndex

The cell containing `p`, by closed-form inversion: point -> (lon, lat) -> tile row ->
band -> tile column -> raster row and column. The steps are ordered because each needs
the previous one — the longitude spacing depends on the band, which depends on the
latitude row. A complete level grid covers the sphere, so this never returns `nothing`.

Ties on a shared boundary are decided by `floor` and are deterministic. A tile owns
`[west, east) x [south, north)`; within a tile the raster row is measured downward from
the north edge, so a point exactly on an interior row boundary goes to the row south of
it instead. Either way the level-0 answer is [`parent`](@ref) of the level-1 answer,
because the tile is chosen first and the raster indices are clamped into it.

The tie only decides anything when the point's latitude really is the edge, and the unit
sphere is a lossy carrier for that: `asind ∘ cosd(90 - ·)` moves a latitude by up to 10
ulps, and it is expansive near the equator, so a third of the 180 tile-row south edges —
60 of them, pinned as `found == 120` by `"cellat agrees with cell_box on south edges"` —
are not in its image at all and cannot be probed through a `UnitSphericalPoint`. What
this method guarantees for the ones that survive is exactness against
[`cell_box`](@ref): the south edge is `Float64(lat_s) + Δlat/2`, and
because `(x + h) - h` is not a `Float64` identity the `floor` below is repaired against
that expression rather than trusted.

At a pole the longitude is whatever `atan` makes of a signed zero, and that depends on
which longitude the point was BUILT from. `TO_SPHERE((lon, ±90))` has
`x = sind(0 or 180) * cosd(lon)`, which is `-0.0` whenever `cosd(lon) < 0`, and
`atand(±0.0, -0.0)` is `±180`. So a pole point built from `|lon| > 90` lands in the W180
tile of the pole row, and one built from `|lon| <= 90` — including the exact
`UnitSphericalPoint(0.0, 0.0, ±1.0)` literals that [`cell_boundary`](@ref) emits — lands
in the `lon_w = 0` tile. Deterministic either way, and the pole is on the boundary of all
360 tiles of the row regardless.
"""
function DGG.cellat(g::LevelGrid{N}, p::GO.UnitSphericalPoint) where {N}
    sys = g.system
    g.level == 0 || g.level == 1 || throw(ArgumentError(
        "level $(g.level) is outside $(DGG.levels(sys))"))
    lon, lat = FROM_SPHERE(p)
    lat = clamp(lat, -90.0, 90.0)
    half_dlat = (1 / N) / 2
    # A tile spans latitudes `[lat_s + Δlat/2, lat_s + 1 + Δlat/2)`. At `lat = 90` the
    # floor is 89; at `lat = -90` it is -91, and the clamp is the extended bottom row.
    lat_s = clamp(floor(Int, lat - half_dlat), -90, 89)
    # `cell_box` builds the south edge as `Float64(lat_s) + half_dlat`; adding then
    # subtracting `half_dlat` is not a Float64 identity, so repair against the edge itself.
    # The two branches are mutually exclusive — the first leaves `lat` at or above the new
    # tile's south edge, which is exactly what the second tests for — so each fires once.
    lat_s <  89 && lat >= Float64(lat_s + 1) + half_dlat && (lat_s += 1)
    lat_s > -90 && lat <  Float64(lat_s)     + half_dlat && (lat_s -= 1)
    r = _row(lat_s)
    nc = ncols(sys, r)
    half_dlon = (1 / nc) / 2
    # `[-180, 180)`, then the same half-pixel offset in longitude. `lon >= 180 - Δlon/2`
    # floors to 180, which is the W180 tile reached from the east; shifting `s` with it
    # keeps the within-tile fraction below in that tile's own frame.
    s = (lon - 360 * floor((lon + 180) / 360)) + half_dlon
    lon_w = floor(s)
    if lon_w >= 180
        lon_w = -180.0
        s -= 360
    end
    q = _col(Int(lon_w))
    g.level == 0 && return DGG.LevelIndex(0, tileordinal(r, q))
    i = clamp(floor(Int, (s - lon_w) * nc), 0, Int(nc) - 1)
    j = clamp(floor(Int, ((lat_s + 1 + half_dlat) - lat) * N), 0, N - 1)
    return DGG.LevelIndex(1, tilebase(sys, r, q) + Int64(j) * nc + Int64(i))
end

# ===========================================================================
# Topology
# ===========================================================================

# WHERE ADJACENCY IS DECIDED. A cell of either level IS its [`cell_box`](@ref), so two
# cells meet exactly where their boxes do, and each axis is settled in integers.
#
# LATITUDE. Level-1 rows are global: `J = r * N + j`, north to south over `180 * N` of
# them. Row `J`'s south edge and row `J + 1`'s north edge are the same rational — at a
# tile seam, `(lat_s(r) + 1) - N/N + 1/(2N)` against `(lat_s(r + 1) + 1) + 1/(2N)`,
# equal because `lat_s(r + 1) = lat_s(r) - 1` — so consecutive rows abut and rows two
# apart are separated by a whole row. Level-0 tile rows say the same one level up. A
# cell can therefore only meet cells in its own row and in the two beside it, and the
# pole clamps do not change that: they alter a row's HEIGHT, never which rows abut.
#
# LONGITUDE. A row with `nc` columns per degree puts its breakpoints at
# `(P k - 1) / (2 nc)` degrees east of -180: `P = 2` for a pixel row, whose cell `K`
# runs `(2K - 1)/(2nc) .. (2K + 1)/(2nc)`, and `P = 2 nc` for a tile row, whose tile `q`
# runs `(2 nc q - 1)/(2 nc) .. (2 nc (q + 1) - 1)/(2 nc)`. Both are the half-pixel
# registration `cell_box` describes — which is why a TILE's corners are not on integer
# degrees either, and why they move with the band — and both partition the circle.
#
# Comparing a row of `a` columns with one of `b` clears both denominators at once by
# scaling into units of `1/(2ab)` degrees, where every endpoint is an integer and the
# circle is `720 a b` long. `_facing` is that comparison and the only place adjacency is
# decided; nothing here reads a `Float64` or a tolerance.

# The two shapes of longitude lattice, as the breakpoint stride and the row length.
@inline _stride(nc::Int64, level::Int) = level == 0 ? 2 * nc : Int64(2)
@inline _rowlen(nc::Int64, level::Int) =
    level == 0 ? Int64(NCOLS_TILES) : Int64(NCOLS_TILES) * nc
@inline _nrows(::CopernicusDEMSystem{N}, level::Int) where {N} =
    level == 0 ? Int64(NROWS) : Int64(NROWS) * Int64(N)
@inline _gridrow(::CopernicusDEMSystem{N}, level::Int, J::Int64) where {N} =
    level == 0 ? Int(J) : Int(fld(J, Int64(N)))

"""
The cells of a row with `b` columns that meet cell `K` of a row with `a` columns:
`(lo, hi, touch_lo, touch_hi)`, where `lo:hi` overlap `K`'s longitude interval in
POSITIVE length and the flanking `lo - 1` / `hi + 1` meet it in a single POINT when the
matching flag is set. Indices are unreduced; callers take them modulo the row length.

`lo` is the facing cell holding `K`'s west edge and `hi` the last one starting before
its east edge, both by integer division on the cross-multiplied endpoints. The run
holds `1 + ceil(b/a)` cells at most — one more than the facing breakpoints that fit
strictly inside `K` — so it never exceeds three: the widest ratio the band table puts
side by side is 2:1, at latitude ±85.
"""
function _facing(a::Int64, b::Int64, level::Int, K::Int64)
    pa, pb = _stride(a, level), _stride(b, level)
    s_w = (pa * K - 1) * b              # K's west edge, in units of 1/(2ab) degrees
    s_e = (pa * (K + 1) - 1) * b        # and its east edge
    step = pb * a                       # one whole cell of the facing row
    lo = fld(s_w + a, step)
    hi = cld(s_e + a, step) - 1
    return (lo, hi, (pb * lo - 1) * a == s_w, (pb * (hi + 1) - 1) * a == s_e)
end

# `(row, column)` <-> id, at whichever level. Level 0 counts tile rows and tile columns,
# level 1 counts raster rows and columns globally, and `_facing` does not care which.
function _gridcell(sys::CopernicusDEMSystem{N}, level::Int, J::Int64, K::Int64) where {N}
    level == 0 && return DGG.LevelIndex(0, tileordinal(Int(J), Int(K)))
    r = fld(J, Int64(N))
    nc = ncols(sys, r)
    q, i = divrem(K, nc)
    return DGG.LevelIndex(1, tilebase(sys, Int(r), Int(q)) + (J - r * Int64(N)) * nc + i)
end

function _gridcoords(sys::CopernicusDEMSystem{N}, level::Int, c::DGG.LevelIndex) where {N}
    r, q, j, i = decode(sys, c)
    level == 0 && return (Int64(r), Int64(q))
    nc = ncols(sys, r)
    return (Int64(r) * Int64(N) + j, Int64(q) * nc + i)
end

# The neighbours across one of the cell's two parallels, WEST to EAST. Point contacts
# are dropped under `Edge()`, which is the whole of the difference between the two
# connectivities away from a pole.
function _across(sys::CopernicusDEMSystem, level::Int, J2::Int64, a::Int64, K::Int64,
        edge_only::Bool)
    b = ncols(sys, _gridrow(sys, level, J2))
    m = _rowlen(b, level)
    lo, hi, touch_lo, touch_hi = _facing(a, b, level, K)
    out = DGG.LevelIndex[]
    (touch_lo && !edge_only) && push!(out, _gridcell(sys, level, J2, mod(lo - 1, m)))
    for l in lo:hi
        push!(out, _gridcell(sys, level, J2, mod(l, m)))
    end
    (touch_hi && !edge_only) && push!(out, _gridcell(sys, level, J2, mod(hi + 1, m)))
    return out
end

# The 1-ring, in the documented order. The helpers above are called from here and from
# nowhere else.
function _ring1(sys::CopernicusDEMSystem{N}, level::Int, c::DGG.LevelIndex,
        connectivity::DGG.Connectivity) where {N}
    edge_only = connectivity isa DGG.Edge
    J, K = _gridcoords(sys, level, c)
    a = ncols(sys, _gridrow(sys, level, J))
    m = _rowlen(a, level)
    apex_n = J == 0                             # this cell's north edge is the +90 apex
    apex_s = J == _nrows(sys, level) - 1        # its south edge is the -90 apex
    # An apex is one point, so it carries the ring under `Vertex()` and nothing under
    # `Edge()`. When it carries the ring, the ring already holds the two laterals.
    apex_ring = !edge_only && (apex_n || apex_s)

    north = DGG.LevelIndex[]                    # the north side, EAST to WEST
    south = DGG.LevelIndex[]                    # the south side, WEST to EAST
    if apex_n
        # From a cell of the pole ring, the ring runs counter-clockwise in increasing
        # eastward offset: the eastern lateral first, over the pole, the western lateral
        # last. The cell half the ring away is the one that lies due north.
        apex_ring && (north = [_gridcell(sys, level, J, mod(K + t, m)) for t in 1:(m - 1)])
    else
        north = _across(sys, level, J - 1, a, K, edge_only)
        reverse!(north)
    end
    if apex_s
        apex_ring && (south = [_gridcell(sys, level, J, mod(K - t, m)) for t in 1:(m - 1)])
    else
        south = _across(sys, level, J + 1, a, K, edge_only)
    end

    out = DGG.LevelIndex[]
    sizehint!(out, length(north) + length(south) + 2)
    isempty(north) || push!(out, north[end])
    apex_ring || push!(out, _gridcell(sys, level, J, mod(K - 1, m)))
    append!(out, south)
    apex_ring || push!(out, _gridcell(sys, level, J, mod(K + 1, m)))
    isempty(north) || append!(out, view(north, 1:(length(north) - 1)))
    return out
end

# Outward shells, breadth-first over the 1-ring. Shell 1 IS the 1-ring and keeps its
# lattice order; the outer shells have no lattice order to keep and are wound by
# measured azimuth in the frame the generic walk fixes from ring 1, so `ring(c, k)` is
# the tail block of `neighbors(c, k)` at every `k`.
function _shells(g::LevelGrid, c::DGG.LevelIndex, steps::Int,
        connectivity::DGG.Connectivity)
    sys = g.system
    shells = Vector{DGG.LevelIndex}[]
    seen = Set{DGG.LevelIndex}((c,))
    frontier = DGG.LevelIndex[c]
    centre = DGG.cell_centroid(sys, c)
    frame = nothing
    for step in 1:steps
        next = DGG.LevelIndex[]
        for x in frontier, y in _ring1(sys, g.level, x, connectivity)
            y in seen && continue
            push!(seen, y)
            push!(next, y)
        end
        if step == 1
            isempty(next) || (frame = DGG.Fallbacks._ring_frame(g, centre, next))
        elseif frame !== nothing
            DGG.Fallbacks._wind!(next, g, centre, frame)
        end
        push!(shells, next)
        isempty(next) && break
        frontier = next
    end
    return shells
end

"""
    neighbors(grid, c, k = 1; connectivity = Vertex()) -> Vector{LevelIndex}
    ring(grid, c, k; connectivity = Vertex()) -> Vector{LevelIndex}

Cells within, or at exactly, `k` adjacency steps of `c`, excluding `c`, in the
interface's rotational order. Closed form at every cell of both levels: tile interiors,
tile edges and corners, the antimeridian, the band boundaries, the ±90 pole rows and
their half-pixel slivers, and the tile lattice itself.

# What adjacent means

Cells are neighbours under `Vertex()` when their closed [`cell_box`](@ref)es share at
least one point, and under `Edge()` when they share a segment of positive length. Both
are decided in exact integer arithmetic on the lattice — the rational endpoints of two
rows cross-multiplied to a common denominator — never by comparing coordinates within a
tolerance. Adjacency is therefore symmetric by construction, and a cell whose corner
falls strictly inside a longer edge of its neighbour is a `Vertex()` neighbour of it.

# The cases the lattice produces

  - **Interior to a tile**, and **tile edges and corners within a band**: eight under
    `Vertex()`, four under `Edge()`. The antimeridian is a tile edge like any other; the
    tiles either side of ±180 abut, at ids 359 tile columns apart.
  - **Band boundaries** (latitude 50/60/70/80/85) face two rows of different column
    count at one parallel. Every cell still meets the other side — the two rows tile the
    same parallel — but across a reduced ratio `p:q`, so a cell on the narrow side faces
    one or two of the wide side's and a cell on the wide side faces two or three of the
    narrow side's. Six to eight under `Vertex()`, four to six under `Edge()`. Corners
    coincide only where `p` and `q` are both odd, true of 5:3 alone among the five
    ratios (3:2, 4:3, 3:2, 5:3, 2:1) and so of latitude ±80 alone; elsewhere a corner
    lands strictly inside the facing edge, which is still a shared point and still a
    `Vertex()` neighbour.
  - **Pole rows** — raster row 0 of a `lat_s = 89` tile, row `N - 1` of a `lat_s = -90`
    one — are the triangles [`cell_boundary`](@ref) emits, and their apex is one point
    shared by the whole ring. So `Vertex()` gives the entire ring plus the three cells
    equatorward, which is what [`max_neighbors`](@ref) is sized for, while `Edge()` gives
    the two laterals and the one overlap equatorward: three.
  - **Level 0** is the same statement over tile rows. A tile's box is offset half a
    PIXEL, so tiles across a band boundary are offset by different amounts and their
    corners never coincide: eight neighbours within a band, seven across one, and 362 on
    a pole tile row.

# Order

Raster rows run north to south. The cycle is counter-clockwise seen from outside the
sphere — the north-side neighbours from east to west, the western lateral, the
south-side neighbours from west to east, the eastern lateral — and the list starts at
the north side's last member, the neighbour immediately west across the north edge. For
a cell interior to its tile that is the familiar

    Vertex()   NW, W, SW, S, SE, E, NE, N
    Edge()     N, W, S, E

and it stays that at a tile edge, where the same eight cells live in other tiles.
Where a side carries more or fewer cells the block grows or shrinks in place. On a pole
row the apex ring is that cell's north (or south) side and holds the laterals itself, so
it is enumerated over the pole from the eastern lateral to the western one; under
`Edge()`, where the apex carries nothing, the north side is empty and the list starts at
the western lateral instead.

`ring(c, k)` is the final block of `neighbors(c, k)`, and `neighbors(c, k)` is the rings
concatenated outward, as the interface requires. Rings past the first carry no lattice
order and are wound by measured azimuth about `cell_centroid(grid, c)`.
"""
function DGG.neighbors(g::LevelGrid, c::DGG.LevelIndex, k::Integer = 1;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    _checked_index(g, c)
    steps == 0 && return DGG.LevelIndex[]
    steps == 1 && return _ring1(g.system, g.level, c, connectivity)
    return reduce(vcat, _shells(g, c, steps, connectivity))
end

"""
    ring(grid, c, k; connectivity = Vertex()) -> Vector{LevelIndex}

The cells at adjacency distance exactly `k`, on the same closed-form adjacency and the
same rotational order [`neighbors`](@ref) documents; `k == 0` is `[c]`.
"""
function DGG.ring(g::LevelGrid, c::DGG.LevelIndex, k::Integer;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    _checked_index(g, c)
    steps == 0 && return DGG.LevelIndex[c]
    steps == 1 && return _ring1(g.system, g.level, c, connectivity)
    shells = _shells(g, c, steps, connectivity)
    return steps <= length(shells) ? shells[steps] : DGG.LevelIndex[]
end
