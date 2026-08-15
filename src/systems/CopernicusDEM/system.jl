# The two-level hierarchy: 64 800 tiles at level 0, one tile's raster at level 1.
# Every relation here is arithmetic on the dense 0-based ordinals of `bands.jl`.

# `levelgrid(CopernicusDEMSystem(...), l)` is the package's `HierarchicalLevelGrid`:
# all 64 800 tiles, or every pixel on Earth, in ordinal order. This system's fast
# paths hang off the alias, and the five primitives it forwards to are the
# `(sys, ...)` methods below.
const LevelGrid{N} = DGG.HierarchicalLevelGrid{CopernicusDEMSystem{N}}

# ===========================================================================
# System interface
# ===========================================================================

DGG.cellindextype(::CopernicusDEMSystem) = DGG.LevelIndex
DGG.levels(::CopernicusDEMSystem) = 0:1
DGG.has_sorted_subtrees(::CopernicusDEMSystem) = true

"""
    max_neighbors(CopernicusDEMSystem(...), connectivity) -> Int

`8` under `Vertex()` and `4` under `Edge()`: the Moore and von Neumann bounds of a
raster lattice.

!!! warning "No fast-path `neighbors`; the generic fallback answers"
    There is no `neighbors` method for this system, so a call resolves to the
    generic `AbstractGrid` walk in `src/fallbacks/locate.jl` — not a `MethodError`.
    That walk goes through `treeify` and does return an answer: at level 0 it gives
    an interior tile its eight Moore neighbours.

    The bound above is the *interior lattice* bound, and it is stated because
    [`max_neighbors`](@ref) has no default. It is **not** a claim about the whole
    sphere: a pixel in the top row of a pole tile meets every other pixel of the pole
    ring at the pole itself, and a pixel just below a band boundary meets up to ten
    coarser pixels above it. Anyone adding a fast-path [`neighbors`](@ref) to this
    system must revisit this number first.
"""
DGG.max_neighbors(::CopernicusDEMSystem, ::DGG.Vertex) = 8
DGG.max_neighbors(::CopernicusDEMSystem, ::DGG.Edge) = 4

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

Lazy on purpose: `PartialGrid` reads `first(rootcells(sys))` on every construction
(`src/fallbacks/partial_grid.jl`), and a materialised 64 800-element vector would
allocate about a megabyte per chunk built. The count is far above the "small, cheap
collection" the contract has in mind (12 for HEALPix), which costs the generic tree
descent one cap evaluation per tile at the synthetic root; that is the price of a grid
whose base tessellation is the 1° graticule.
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
    `O(subtree)` ids, and on every other system it hands back a freshly allocated
    `Vector` the caller owns. **This method does not.** It returns a **lazy, read-only
    `AbstractVector`** that computes each id on indexing. Reading is complete — `length`,
    `getindex`, iteration, `collect` — but nothing that writes works, because there is no
    array to write into: no `setindex!`, no `push!`, no `sort!`, and no passing it to an
    API that mutates its argument. **`collect` it first if you need any of those.** The
    divergence is deliberate: one GLO-30 tile has 12 960 000 level-1 descendants, and a
    `Vector{LevelIndex}` of them is 16 bytes apiece.
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
for a 1° tile and about 10 μm for a 1-arcsec pixel. That small-angle expression is an
estimate, not a bound, and it errs low: the exact bow of a 1° edge at latitude 45 is
1.9038830e-5 rad against the formula's 1.9038589e-5.

The two bows nearly cancel, so ring and box differ in area only slightly. The
`"ring vs box"` testset in `test/systems/CopernicusDEM/runtests.jl` measures that gap
over every band and both pole rows, `@info`s it as `worst_tile` / `worst_pixel` /
`worst_pole_pixel`, and bounds each: worst over tiles `< 1e-4` (5.1e-5, at the ±90
rows), worst over pixels outside the ±90 tile rows `< 1e-8` (8.1e-10 GLO-30, 3.8e-9
GLO-90), and worst over pixels inside them `< 1e-4` (1.4e-5 GLO-30, 1.2e-6 GLO-90) —
those cells being half-pixel slivers and pole triangles. Run it for the current
numbers rather than trusting these.

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
`|Δλ_below² − Δλ_above²|/8 · sin φ cos φ`: in GLO-30 that is 1.8e-12 rad (11.5 μm) at
latitude 50, about 4e-7 of that pixel's own height, rising to 1.9e-11 rad (122 μm) at
latitude 85. Those are values of the small-angle formula, not measurements of the
emitted rings. The box tessellation is exact there; the quad tessellation is exact to
that order.

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
telescope exactly. 890 of the 64 800 GLO-90 tiles have `east - west !== 1.0`, by up to
1.4e-14, so consecutive terms do not cancel to the last bit. The sum is accurate
because those errors are tiny and because pairwise summation stops them accumulating,
not because anything cancels: **materialise** the 64 800 tile areas into a `Vector` and
Julia's pairwise `sum` lands 3.6e-15 from 4π, while a generator — or any sequential
accumulation — lands 2.9e-12 (GLO-30) or 1.4e-12 (GLO-90) away and fails the
`rtol = 1e-14` the suite asserts. That assertion, and those figures, are
`"the boxes partition the sphere"` in `test/systems/CopernicusDEM/runtests.jl`; do not
loosen its tolerance to accommodate a generator.

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
