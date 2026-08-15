# The lattice: the system singleton, the DGED band table, the per-`N` column and
# prefix-sum tables, and the id codec that turns a tile label or a raster
# (row, column) into a `LevelIndex` and back.

# ===========================================================================
# The system
# ===========================================================================

"""
    CopernicusDEMSystem(30)   # GLO-30, 1 arcsec latitude spacing
    CopernicusDEMSystem(90)   # GLO-90, 3 arcsec latitude spacing

The Copernicus DEM lattice as a two-level `AbstractHierarchicalGridSystem`.

The type parameter is the number of **latitude intervals per degree** — 3600 for
GLO-30, 1200 for GLO-90 — not the nominal metre figure, because every formula in
this module is written in it. Any `N` divisible by 30 is structurally valid (all
six reduction factors then divide exactly); the test suite instantiates
`CopernicusDEMSystem{30}()` as a scaled twin of the shipped pair.

`levels(sys) == 0:1`: level 0 is a 1°x1° tile, level 1 is a pixel.
"""
struct CopernicusDEMSystem{N} <: DGG.AbstractHierarchicalGridSystem end

function CopernicusDEMSystem(spacing::Integer)
    spacing == 30 && return CopernicusDEMSystem{3600}()
    spacing == 90 && return CopernicusDEMSystem{1200}()
    throw(ArgumentError(
        "Copernicus DEM ships GLO-30 and GLO-90; got a nominal spacing of $spacing m. " *
        "For a scaled twin, name the lattice directly: CopernicusDEMSystem{N}()"))
end

Base.show(io::IO, ::CopernicusDEMSystem{3600}) = print(io, "CopernicusDEMSystem(30)")
Base.show(io::IO, ::CopernicusDEMSystem{1200}) = print(io, "CopernicusDEMSystem(90)")
Base.show(io::IO, ::CopernicusDEMSystem{N}) where {N} =
    print(io, "CopernicusDEMSystem{", N, "}()")

"Latitude intervals per degree: 3600 for GLO-30, 1200 for GLO-90."
lat_intervals(::CopernicusDEMSystem{N}) where {N} = N

# ===========================================================================
# The band table
# ===========================================================================

# The DGED longitude reduction factors, by the latitude band of a tile's
# EQUATOR-WARD edge. Stored doubled so `ncols = 2N / factor2` is exact integer
# division for every band: 1.5x and 5x were chosen by the product spec precisely
# so that 3600 and 1200 divide without remainder.
#
#   band (deg)   factor   GLO-30 cols   GLO-90 cols
#   [ 0, 50)     1x       3600          1200
#   [50, 60)     1.5x     2400           800
#   [60, 70)     2x       1800           600
#   [70, 80)     3x       1200           400
#   [80, 85)     5x        720           240
#   [85, 90)     10x       360           120
#
# Three independent confirmations: Product Handbook Table 3, the AWS bucket readme,
# and direct measurement of 44 real COGs. See the research document, section 1.2.
#
# This is the DGED table. DTED (`.dt1`/`.dt2`) has five bands with factors
# 1, 2, 3, 4, 6 and is not implemented here; the AWS buckets ship DGED.
const BAND_EDGES   = (0, 50, 60, 70, 80, 85)
const BAND_FACTOR2 = (2, 3, 4, 6, 10, 20)

"""
    band_factor2(lat_s) -> Int

Twice the longitude reduction factor for the tile whose lower-left corner latitude is
the integer `lat_s`, i.e. the tile spanning `lat_s` to `lat_s + 1`.

**The band is chosen by the tile's EQUATOR-WARD edge**, `min(|lat_s|, |lat_s + 1|)`,
with half-open intervals. In the northern hemisphere the label is the equator-ward
edge; in the southern hemisphere it is the pole-ward one. So `N50` (50 -> 51, edge 50)
is 1.5x and 2400 columns wide, while `S50` (-50 -> -49, edge **49**) is 1x and **3600**
columns wide. That asymmetry is not stated in either primary source; it was
established by measuring 44 real tiles and is what `test/systems/CopernicusDEM/`'s
fixture table pins.
"""
function band_factor2(lat_s::Integer)
    -90 <= lat_s <= 89 || throw(ArgumentError(
        "tile latitude $lat_s is outside -90:89"))
    e = min(abs(lat_s), abs(lat_s + 1))
    f = BAND_FACTOR2[1]
    for k in eachindex(BAND_EDGES)
        e >= BAND_EDGES[k] && (f = BAND_FACTOR2[k])
    end
    return f
end

# ===========================================================================
# The tile lattice
# ===========================================================================

"Per-`N` lattice tables: column count and cumulative pixel count, by tile row."
struct Tables
    ncols::Vector{Int64}     # 180 entries, index r+1, r = 89 - lat_s
    rowbase::Vector{Int64}   # 181 entries, pixels in tile rows 0:r-1 at index r+1
end

_lat_s(r::Integer) = 89 - Int(r)         # tile row 0 is lat_s = 89 (north-to-south)
_row(lat_s::Integer) = 89 - Int(lat_s)
_lon_w(q::Integer) = Int(q) - 180        # tile column 0 is lon_w = -180
_col(lon_w::Integer) = Int(lon_w) + 180

const NROWS = 180
const NCOLS_TILES = 360
const NTILES = NROWS * NCOLS_TILES       # 64 800

# Every product is widened to `Int64` before it is taken: `ncells(sys, 1)` is
# 620 524 800 000 at GLO-30, so a 32-bit `Int` build would wrap silently here and
# nowhere else.
function build_tables(N::Integer)
    N > 0 || throw(ArgumentError("N must be positive, got $N"))
    ncols = Vector{Int64}(undef, NROWS)
    rowbase = Vector{Int64}(undef, NROWS + 1)
    rowbase[1] = 0
    for r in 0:(NROWS - 1)
        f2 = band_factor2(_lat_s(r))
        (2 * N) % f2 == 0 || throw(ArgumentError(
            "N = $N does not divide evenly by the $(f2 / 2)x reduction factor; " *
            "a Copernicus DEM lattice needs 30 | N"))
        ncols[r + 1] = (2 * Int64(N)) ÷ f2
        rowbase[r + 2] = rowbase[r + 1] + Int64(NCOLS_TILES) * ncols[r + 1] * Int64(N)
    end
    return Tables(ncols, rowbase)
end

const GLO30_TABLES = build_tables(3600)
const GLO90_TABLES = build_tables(1200)
const OTHER_TABLES = Dict{Int,Tables}()
const OTHER_LOCK = ReentrantLock()

@inline tables(::CopernicusDEMSystem{3600}) = GLO30_TABLES
@inline tables(::CopernicusDEMSystem{1200}) = GLO90_TABLES
# Any other `N` is a scaled twin, used by the test suite so the conformance harness's
# `collect(children(...))` stays affordable. Looked up rather than dispatched, because
# there is no reason to compile a table into the method table for a one-off.
function tables(::CopernicusDEMSystem{N}) where {N}
    return lock(OTHER_LOCK) do
        get!(() -> build_tables(Int(N)), OTHER_TABLES, Int(N))
    end
end

"Columns in the tiles of tile row `r` (equivalently, of the tile at latitude `lat_s`)."
ncols(sys::CopernicusDEMSystem, r::Integer) = tables(sys).ncols[Int(r) + 1]
ncols_at(sys::CopernicusDEMSystem, lat_s::Integer) = ncols(sys, _row(lat_s))

# ===========================================================================
# The id codec
# ===========================================================================

"""
    tileordinal(r, q) -> Int

The 0-based level-0 id of the tile at tile row `r` (north to south, `r = 89 - lat_s`)
and tile column `q` (west to east, `q = lon_w + 180`): `r * 360 + q`.
"""
@inline tileordinal(r::Integer, q::Integer) = Int(r) * NCOLS_TILES + Int(q)

"""
    tilebase(sys, r, q) -> Int64

The 0-based level-1 id of the first (north-west) pixel of the tile at tile row `r`
and tile column `q`.

The prefix sum makes this O(1): `rowbase[r + 1]` counts every pixel in the tile rows
north of `r`, and within a row every tile has the same `ncols * N` pixels.
"""
tilebase(sys::CopernicusDEMSystem{N}, r::Integer, q::Integer) where {N} =
    tables(sys).rowbase[Int(r) + 1] + Int64(q) * ncols(sys, r) * Int64(N)

"""
    tilecell(sys, lat_s, lon_w) -> LevelIndex

The level-0 cell of the tile whose lower-left corner is the integer degree pair
`(lat_s, lon_w)` — the pair the AWS file name carries, e.g. `N50_00_E006_00` is
`(50, 6)` and `S90_00_W180_00` is `(-90, -180)`.
"""
function tilecell(::CopernicusDEMSystem, lat_s::Integer, lon_w::Integer)
    -90 <= lat_s <= 89 || throw(ArgumentError("tile latitude $lat_s is outside -90:89"))
    -180 <= lon_w <= 179 || throw(ArgumentError("tile longitude $lon_w is outside -180:179"))
    return DGG.LevelIndex(0, tileordinal(_row(lat_s), _col(lon_w)))
end

"""
    tilecorner(sys, c) -> (lat_s, lon_w)

The lower-left integer degree corner of the tile `c` (level 0) or of the tile
containing the pixel `c` (level 1). Inverse of [`tilecell`](@ref).
"""
function tilecorner(sys::CopernicusDEMSystem, c::DGG.LevelIndex)
    r, q, _, _ = decode(sys, c)
    return (_lat_s(r), _lon_w(q))
end

"""
    pixelcell(sys, tile, j, i) -> LevelIndex

The level-1 cell of the pixel at raster row `j` (0-based, **north row first**) and
raster column `i` (0-based, west to east) of `tile`.
"""
function pixelcell(sys::CopernicusDEMSystem{N}, tile::DGG.LevelIndex,
        j::Integer, i::Integer) where {N}
    DGG.level(tile) == 0 || throw(ArgumentError(
        "$tile is not a level-0 tile, so it has no raster of its own"))
    r, q, _, _ = decode(sys, tile)
    nc = ncols(sys, r)
    0 <= j < N || throw(ArgumentError("raster row $j is outside 0:$(N - 1)"))
    0 <= i < nc || throw(ArgumentError("raster column $i is outside 0:$(nc - 1)"))
    return DGG.LevelIndex(1, tilebase(sys, r, q) + Int64(j) * nc + Int64(i))
end

"""
    decode(sys, c) -> (r, q, j, i)

Tile row, tile column, and — for a level-1 cell — raster row and column within the
tile. A level-0 cell returns `(r, q, 0, 0)`.

The level-1 inverse is a binary search for the tile row followed by two `divrem`s,
because the pixels of a tile row are `360 * ncols(r) * N` consecutive ids and `ncols`
varies only by row.
"""
function decode(sys::CopernicusDEMSystem{N}, c::DGG.LevelIndex) where {N}
    index = _checked_index(sys, c)
    if DGG.level(c) == 0
        r, q = divrem(index, Int64(NCOLS_TILES))
        return (Int(r), Int(q), 0, 0)
    end
    rowbase = tables(sys).rowbase
    # `index < rowbase[end]` is guaranteed above, so this never lands past row 179.
    r = searchsortedlast(rowbase, index) - 1
    nc = ncols(sys, r)
    q, rest = divrem(index - rowbase[r + 1], nc * Int64(N))
    j, i = divrem(rest, nc)
    return (Int(r), Int(q), Int(j), Int(i))
end

# The id guard every decoder needs: `decode` would happily search the prefix table
# with an id no cell has and hand back a tile that does not exist rather than an
# error. (`cellposition` deliberately does NOT use this — there, a miss is
# `nothing`.) The `ncells` call is also what rejects a level outside `0:1`.
@inline function _checked_index(sys::CopernicusDEMSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    n = DGG.ncells(sys, l)
    0 <= c.index < n || throw(ArgumentError(
        "Copernicus DEM id $(c.index) is out of range 0:$(n - 1) at level $l"))
    return c.index
end
