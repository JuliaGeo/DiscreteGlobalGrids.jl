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

`36N + 2` under `Vertex()` — a pole-row pixel touches the whole pole ring plus the
three pixels below it — and `6` under `Edge()`, attained on the wide side of a band
boundary. Both bound every cell of both levels, and both are attained.
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

All 64 800 level-0 tiles, `LevelIndex(0, 0:64799)`, ascending, as a **lazy** vector —
`PartialGrid` reads `first(rootcells(sys))` on every construction, and materialising
would allocate about a megabyte each time.
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

The AWS COGs are `RasterPixelIsPoint` (`AREA_OR_POINT=Point`): pixel CENTRES sit on
the integer-degree lattice, so a pixel's box is its centre ± half a pixel and a tile's
box is the nominal 1°x1° box shifted **half a pixel west and half a pixel north**.
Tiles abut exactly in latitude and, within a band, in longitude.

Two pole corrections keep the tiling exact: the top row of the `lat_s = 89` tiles is
clamped to +90, and the bottom row of the `lat_s = -90` tiles is extended to -90 (one
and a half pixels tall). Those cells have one degenerate edge and are spherical
TRIANGLES; [`cell_boundary`](@ref) emits three vertices for them.
"""
function cell_box(sys::CopernicusDEMSystem{N}, c::DGG.LevelIndex) where {N}
    r, q, j, i = decode(sys, c)
    lat_s = _lat_s(r)
    lon_w = _lon_w(q)
    nc = ncols(sys, r)
    # A tile is its pixel grid's outer frame, so both levels are the same expression
    # over a column and a row interval.
    i_w, i_e, j_n, j_s = DGG.level(c) == 0 ? (0, Int(nc), 0, N) : (i, i + 1, j, j + 1)
    # `k / nc`, never `k * Δlon`: a cell's east edge must be the same `Float64`
    # as its neighbour's west edge so the quads share a geodesic exactly.
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

# The exact ±90 vertices, as shared literals. `TO_SPHERE((lon, ±90))` carries signed
# zeros in x and y that vary with `lon`: `==`-equal but `===`- and `isequal`-distinct,
# and `atand` reads the sign bits — see [`cellat`](@ref).
const NORTH_POLE = GO.UnitSphericalPoint(0.0, 0.0, 1.0)
const SOUTH_POLE = GO.UnitSphericalPoint(0.0, 0.0, -1.0)

"""
    cell_boundary(grid, c) -> Vector{UnitSphericalPoint}

The cell's box as a **plain 4-corner great-circle quadrilateral**, implicitly closed and
counter-clockwise seen from outside the sphere, in the order
`(W,S) -> (E,S) -> (E,N) -> (W,N)`.

Not densified: undensified lat/lon quads are convex, which keeps this system exact as
a regridding destination, and the poleward bow off the true box edge is only about
`Δλ²/16` rad. Adjacent cells within a band share their corner points bit-identically,
so the quads tile the sphere with no gaps.

A cell touching ±90 is emitted as a TRIANGLE, same cyclic order with the duplicate
dropped, its apex the exact literal `UnitSphericalPoint(0.0, 0.0, ±1.0)`.
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

The midpoint of the cell's box — the DEM post itself for a pixel — strictly inside
the published quad.
"""
function DGG.cell_centroid(sys::CopernicusDEMSystem, c::DGG.LevelIndex)
    west, east, south, north = cell_box(sys, c)
    return TO_SPHERE(((west + east) / 2, (south + north) / 2))
end

"""
    cell_area(grid, c) -> Float64

The exact solid angle of the cell's lat/lon BOX, `Δλ · (sin φ_N − sin φ_S)` steradians,
in O(1) — the box, not the published ring, which differs by the bow
[`cell_boundary`](@ref) describes. Summing every cell's area to 4π needs pairwise
summation: materialise into a `Vector` before `sum`, not a generator.
"""
function DGG.cell_area(g::LevelGrid, c::DGG.LevelIndex)
    _checked_index(g, c)
    west, east, south, north = cell_box(g.system, c)
    return deg2rad(east - west) * (sind(north) - sind(south))
end

"""
    cell_extent(grid, c) -> Extents.Extent{(:X, :Y)}

The cell's [`cell_box`](@ref), verbatim — not derived from the ring as the generic
fallback would, whose poleward bow makes vertically adjacent tiles report overlapping
extents.
"""
function DGG.cell_extent(g::LevelGrid, c::DGG.LevelIndex)
    _checked_index(g, c)
    west, east, south, north = cell_box(g.system, c)
    return DGG.Extents.Extent(X = (west, east), Y = (south, north))
end

"""
    node_extent(CopernicusDEMSystem(...), c) -> SphericalCap

The cap centred on the cell centre, with radius the largest corner distance plus an
analytic `Δλ²/16` pad for the bow of descendants' rings, and one outward ULP. The
farthest point of a lat/lon box from its midpoint is a corner, so four distances
suffice. Radii stay far below 90°, so `require_convex_extents = true` holds.
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

The cell containing `p`, by closed-form inversion. A complete level grid covers the
sphere, so this never returns `nothing`.

Ties on a shared boundary are deterministic: a tile owns
`[west, east) x [south, north)`, and a point exactly on an interior raster-row
boundary goes to the row south of it. The level-0 answer is [`parent`](@ref) of the
level-1 answer. At a pole the longitude is whatever `atan` makes of the point's
signed zeros, so which pole-row tile answers depends on the longitude the point was
built from — deterministic either way, and the pole is on the boundary of all 360
tiles of the row regardless.
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
    # `(x + h) - h` is not a Float64 identity, so repair the floor against the south
    # edge exactly as `cell_box` builds it. The two branches are mutually exclusive.
    lat_s <  89 && lat >= Float64(lat_s + 1) + half_dlat && (lat_s += 1)
    lat_s > -90 && lat <  Float64(lat_s)     + half_dlat && (lat_s -= 1)
    r = _row(lat_s)
    nc = ncols(sys, r)
    half_dlon = (1 / nc) / 2
    # `[-180, 180)`, then the same half-pixel offset. `lon >= 180 - Δlon/2` floors to
    # 180 — the W180 tile reached from the east — so shift `s` into that tile's frame.
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

# WHERE ADJACENCY IS DECIDED. A cell IS its [`cell_box`](@ref), so two cells meet
# exactly where their boxes do, and each axis is settled in integers. Rows abut only
# with their immediate neighbours, so a cell can only meet cells in its own row and
# the two beside it. Longitude endpoints of two rows are cross-multiplied into units
# of `1/(2ab)` degrees, where every endpoint is an integer; `_facing` is that
# comparison and the only place adjacency is decided — nothing here reads a `Float64`
# or a tolerance.

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
its east edge. The run never exceeds three cells: the widest ratio the band table puts
side by side is 2:1.
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
        # The pole ring runs counter-clockwise in increasing eastward offset:
        # eastern lateral first, over the pole, western lateral last.
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
interface's rotational order. Closed form at every cell of both levels, including band
boundaries, the antimeridian, and the ±90 pole rows.

Cells are neighbours under `Vertex()` when their closed [`cell_box`](@ref)es share at
least one point, and under `Edge()` when they share a segment of positive length —
decided in exact integer arithmetic, never by tolerance. Across a band boundary a cell
faces up to three cells of the other side; a pole-row cell's apex is shared by the
whole pole ring, which is what [`max_neighbors`](@ref) is sized for.

The 1-ring is counter-clockwise seen from outside the sphere — north side east to
west, western lateral, south side west to east, eastern lateral — starting at the
neighbour immediately west across the north edge: `NW, W, SW, S, SE, E, NE, N` for an
interior cell, `N, W, S, E` under `Edge()`. On a pole row the apex ring is the cell's
north (or south) side and is enumerated over the pole from the eastern lateral to the
western; under `Edge()` the apex carries nothing and the list starts at the western
lateral. Rings past the first carry no lattice order and are wound by measured
azimuth about `cell_centroid(grid, c)`; `ring(c, k)` is the final block of
`neighbors(c, k)`.
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
