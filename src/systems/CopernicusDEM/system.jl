# Two levels: global 1° tiles, then each tile's raster.

# Complete level grids use global ordinal order.
const LevelGrid{N} = DGG.HierarchicalLevelGrid{CopernicusDEMSystem{N}}

# ===========================================================================
# System interface
# ===========================================================================

DGG.cellindextype(::CopernicusDEMSystem) = DGG.LevelIndex
DGG.levels(::CopernicusDEMSystem) = 0:1
DGG.has_sorted_subtrees(::CopernicusDEMSystem) = true

"""
    maxneighbors(CopernicusDEMSystem{N}(), connectivity) -> Int

`36N + 2` under `Vertex()` and `6` under `Edge()`. Both bounds are attained.
"""
DGG.maxneighbors(sys::CopernicusDEMSystem, ::DGG.Vertex) =
    360 * Int(max(ncols(sys, 0), ncols(sys, NROWS - 1))) + 2
DGG.maxneighbors(::CopernicusDEMSystem, ::DGG.Edge) = 6

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

All 64 800 level-0 tiles, `LevelIndex(0, 0:64799)`, as an ascending lazy vector.
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

Returns `index + 1` for an in-range id and `nothing` otherwise; it never throws.
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

The tile's `ncols * N` pixels as a lazy vector in north-to-south, west-to-east
raster order.
"""
function DGG.children(sys::CopernicusDEMSystem{N}, c::DGG.LevelIndex) where {N}
    l = DGG.level(c)
    l == 0 || throw(ArgumentError(l == 1 ?
        "level-1 Copernicus DEM cell $c is a pixel, at maxlevel 1, and has no children" :
        "level $l is outside $(DGG.levels(sys))"))
    r, q, _, _ = decode(sys, c)
    return IdRange(Int32(1), tilebase(sys, r, q), Int(ncols(sys, r)) * N)
end

"""
    ancestor(CopernicusDEMSystem(...), c, l) -> LevelIndex

Returns `c` at its own level or its level-0 [`parent`](@ref).
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

The tile's exact, contiguous level-1 position window:
`tilebase + 1 : tilebase + ncols*N`.

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
    target <= DGG.maxlevel(sys) || throw(ArgumentError(
        "descendant level $target is past maxlevel $(DGG.maxlevel(sys))"))
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

Every level-`l` descendant of `c`, ascending, as a lazy vector over
[`descendant_range`](@ref).

!!! warning "This diverges from the interface"
    Unlike the interface default, this returns a read-only `AbstractVector`, not
    an owned `Vector`. Call `collect` before mutation or passing it to mutating APIs.
"""
function DGG.descendants(sys::CopernicusDEMSystem, c::DGG.LevelIndex, l::Integer)
    r = DGG.descendant_range(sys, c, l)     # validates `l` both ways
    return IdRange(Int32(l), Int64(first(r) - 1), length(r))
end

# Also require the id to match the grid's level.
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

The cell's closed longitude/latitude box in **degrees**.

AWS COGs are `RasterPixelIsPoint` (`AREA_OR_POINT=Point`): pixel centres lie on
the lattice, and tile boxes extend half a pixel west and north of the nominal tile.

The top row of the `lat_s = 89` tiles is clamped to +90 and the bottom row of the
`lat_s = -90` tiles extended to -90; those cells are triangles, and
[`cell_boundary`](@ref) emits three vertices.
"""
function cell_box(sys::CopernicusDEMSystem{N}, c::DGG.LevelIndex) where {N}
    r, q, j, i = decode(sys, c)
    lat_s = _lat_s(r)
    lon_w = _lon_w(q)
    nc = ncols(sys, r)
    # A tile is its pixel grid's outer frame.
    i_w, i_e, j_n, j_s = DGG.level(c) == 0 ? (0, Int(nc), 0, N) : (i, i + 1, j, j + 1)
    # Division makes shared edges bit-identical between neighbours.
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

# Shared pole literals avoid longitude-dependent signed zeros.
const NORTH_POLE = GO.UnitSphericalPoint(0.0, 0.0, 1.0)
const SOUTH_POLE = GO.UnitSphericalPoint(0.0, 0.0, -1.0)

"""
    cell_boundary(grid, c) -> Helpers.SmallList{4,UnitSphericalPoint}

The closed box as a 4-corner great-circle quadrilateral, counter-clockwise from
outside the sphere, in the order
`(W,S) -> (E,S) -> (E,N) -> (W,N)`.

Edges are not densified, so the quads are convex — exact as clip windows — with a
poleward bow of about `Δλ²/16` radians. Adjacent cells within a band share
corners bit-identically.

A pole cell is a triangle with the duplicate corner dropped and an exact pole apex.

Storage is inline, as IGeo7's rings are: the ring, its [`closed_ring`](@ref), and
the polygon built on it are `isbits`, so a boundary read never reaches the heap.
"""
function DGG.cell_boundary(sys::CopernicusDEMSystem, c::DGG.LevelIndex)
    west, east, south, north = cell_box(sys, c)
    # A quad's worth of slots; a pole triangle leaves the fourth unused.
    ring = Helpers.empty_small_list(Val(4), NORTH_POLE)
    if north == 90.0
        ring = Helpers.small_push(ring, TO_SPHERE((west, south)))
        ring = Helpers.small_push(ring, TO_SPHERE((east, south)))
        return Helpers.small_push(ring, NORTH_POLE)
    elseif south == -90.0
        ring = Helpers.small_push(ring, SOUTH_POLE)
        ring = Helpers.small_push(ring, TO_SPHERE((east, north)))
        return Helpers.small_push(ring, TO_SPHERE((west, north)))
    end
    ring = Helpers.small_push(ring, TO_SPHERE((west, south)))
    ring = Helpers.small_push(ring, TO_SPHERE((east, south)))
    ring = Helpers.small_push(ring, TO_SPHERE((east, north)))
    return Helpers.small_push(ring, TO_SPHERE((west, north)))
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

The exact box solid angle, `Δλ · (sin φ_N − sin φ_S)` steradians. This is not the
published ring's area. Materialise before `sum` when accurate pairwise reduction matters.
"""
function DGG.cell_area(g::LevelGrid, c::DGG.LevelIndex)
    _checked_index(g, c)
    west, east, south, north = cell_box(g.system, c)
    return deg2rad(east - west) * (sind(north) - sind(south))
end

"""
    cell_extent(grid, c) -> Extents.Extent{(:X, :Y)}

The cell's [`cell_box`](@ref), not the bowed ring's extent.
"""
function DGG.cell_extent(g::LevelGrid, c::DGG.LevelIndex)
    _checked_index(g, c)
    west, east, south, north = cell_box(g.system, c)
    return DGG.Extents.Extent(X = (west, east), Y = (south, north))
end

"""
    node_extent(CopernicusDEMSystem(...), c) -> SphericalCap

A cap covering the cell and descendant rings, padded for edge bow and rounding.
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

The cell containing `p`. Complete level grids cover the sphere, so this never
returns `nothing`.

Boundary ownership is `[west, east) x [south, north)`; an interior raster-row
boundary belongs to the southern row. At a pole, signed zeros determine the
longitude and therefore which pole-row tile is returned.
"""
function DGG.cellat(g::LevelGrid{N}, p::GO.UnitSphericalPoint) where {N}
    sys = g.system
    g.level == 0 || g.level == 1 || throw(ArgumentError(
        "level $(g.level) is outside $(DGG.levels(sys))"))
    lon, lat = FROM_SPHERE(p)
    lat = clamp(lat, -90.0, 90.0)
    half_dlat = (1 / N) / 2
    # Clamp the shifted floor onto the extended pole rows.
    lat_s = clamp(floor(Int, lat - half_dlat), -90, 89)
    # Repair Float64 cancellation against the exact `cell_box` edge.
    lat_s <  89 && lat >= Float64(lat_s + 1) + half_dlat && (lat_s += 1)
    lat_s > -90 && lat <  Float64(lat_s)     + half_dlat && (lat_s -= 1)
    r = _row(lat_s)
    nc = ncols(sys, r)
    half_dlon = (1 / nc) / 2
    # Wrap the half-pixel-shifted longitude into the W180 tile.
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

# Adjacency follows closed `cell_box` intersections in exact integer arithmetic.

# The two shapes of longitude lattice, as the breakpoint stride and the row length.
@inline _stride(nc::Int64, level::Int) = level == 0 ? 2 * nc : Int64(2)
@inline _rowlen(nc::Int64, level::Int) =
    level == 0 ? Int64(NCOLS_TILES) : Int64(NCOLS_TILES) * nc
@inline _nrows(::CopernicusDEMSystem{N}, level::Int) where {N} =
    level == 0 ? Int64(NROWS) : Int64(NROWS) * Int64(N)
@inline _gridrow(::CopernicusDEMSystem{N}, level::Int, J::Int64) where {N} =
    level == 0 ? Int(J) : Int(fld(J, Int64(N)))

"""
Returns `(lo, hi, touch_lo, touch_hi)` for cells in a `b`-column row facing cell
`K` in an `a`-column row. `lo:hi` overlap by positive length; flagged flanking
cells touch at one point. Indices are unreduced and must be taken modulo row length.

The run contains at most three cells.
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

# Convert global row/column coordinates to ids at either level.
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

# Neighbours across a parallel, west to east; `Edge()` drops point contacts.
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

"""
    one_ring(grid, c, connectivity) -> Vector{LevelIndex}

The immediate neighbours of `c`, counter-clockwise seen from outside: `NW, W,
SW, S, SE, E, NE, N` for an interior `Vertex()` cell and `N, W, S, E` for
`Edge()`. A pole ring runs from the eastern lateral over the pole to the
western.
"""
function DGG.one_ring(g::LevelGrid, c::DGG.LevelIndex,
        connectivity::DGG.Connectivity)
    sys, level = g.system, g.level
    edge_only = connectivity isa DGG.Edge
    J, K = _gridcoords(sys, level, c)
    a = ncols(sys, _gridrow(sys, level, J))
    m = _rowlen(a, level)
    apex_n = J == 0                             # this cell's north edge is the +90 apex
    apex_s = J == _nrows(sys, level) - 1        # its south edge is the -90 apex
    # A pole apex carries its row only under `Vertex()`.
    apex_ring = !edge_only && (apex_n || apex_s)

    north = DGG.LevelIndex[]                    # the north side, EAST to WEST
    south = DGG.LevelIndex[]                    # the south side, WEST to EAST
    if apex_n
        # Counter-clockwise from the eastern to the western lateral.
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

"""
    neighbors(grid, c, k = 1; connectivity = Vertex()) -> Vector{LevelIndex}
    ring(grid, c, k; connectivity = Vertex()) -> Vector{LevelIndex}

Cells within, or exactly at, `k` adjacency steps of `c`, excluding `c`.

`Vertex()` requires any shared point between closed boxes; `Edge()` requires a
shared segment of positive length. Adjacency uses exact integer arithmetic. A
pole apex is shared by its whole row under `Vertex()`.

The 1-ring is counter-clockwise from outside: `NW, W, SW, S, SE, E, NE, N` for
an interior `Vertex()` cell and `N, W, S, E` for `Edge()`. Pole rings run from
the eastern lateral over the pole to the western. Later rings are ordered by
azimuth about the cell centre, from the spoke through the 1-ring's first entry;
`ring(c, k)` is the final block of `neighbors(c, k)`.
"""
function DGG.neighbors(g::LevelGrid, c::DGG.LevelIndex, k::Integer = 1;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = DGG.checked_steps(k)
    _checked_index(g, c)
    steps == 0 && return DGG.LevelIndex[]
    steps == 1 && return DGG.one_ring(g, c, connectivity)
    return reduce(vcat, DGG.adjacency_shells(g, c, steps, connectivity))
end

"""
    ring(grid, c, k; connectivity = Vertex()) -> Vector{LevelIndex}

The cells at adjacency distance exactly `k`, on the same closed-form adjacency and the
same rotational order [`neighbors`](@ref) documents; `k == 0` is `[c]`.
"""
function DGG.ring(g::LevelGrid, c::DGG.LevelIndex, k::Integer;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = DGG.checked_steps(k)
    _checked_index(g, c)
    steps == 0 && return DGG.LevelIndex[c]
    steps == 1 && return DGG.one_ring(g, c, connectivity)
    shells = DGG.adjacency_shells(g, c, steps, connectivity)
    return steps <= length(shells) ? shells[steps] : DGG.LevelIndex[]
end
