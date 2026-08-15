# ---------------------------------------------------------------------------
# Measured fixtures for the Copernicus DEM suite. NO NETWORK: every number here
# was read off a real AWS Open Data COG once, by hand, and committed. The suite
# that reads them never touches the network, and neither does anything else in
# `test/`.
#
# Provenance: the research document's sections 9.3 (boundary probe, 36 rows — 28
# GLO-30 and 8 GLO-90), 9.4 (stratified random verification sweep, 33 GLO-30
# tiles) and 9.5 (GLO-90 southern-boundary confirmation, 11 rows), all produced
# by `ArchGDAL` against `https://copernicus-dem-{30,90}m.s3.amazonaws.com/`.
# 79 distinct COGs in all — 61 GLO-30 and 18 GLO-90, one fewer than the 80 rows
# because 9.5 re-measures 9.3's `S85_00_E000_00`, and both readings are kept
# below. The two arrays here are exactly those 61 and 18 tiles.
#
# These are the ONLY external evidence this system has. The band table in
# `src/systems/CopernicusDEM/bands.jl` agrees with the Product Handbook's Table 3
# and with the AWS bucket readme on the six reduction factors, but neither source
# states which side of a band boundary a SOUTHERN tile falls on — and the answer
# is not the one a reader of the tile label expects (`S50` is 1x and full width;
# `S51` is the 1.5x tile). That asymmetry is measurement, and only measurement,
# which is why these rows are here rather than derived.
# ---------------------------------------------------------------------------

"""
Measured GLO-30 (`N = 3600`) tiles: the lower-left integer degree corner the AWS
file name carries, and the tile's column count.

`S49_00` is `lat_s = -49`, `W002` is `lon_w = -2`. Grouped by band, and within
each band the northern rows come before the southern ones — the southern rows
straddling a band edge (`-50`/`-51`, `-60`/`-61`, `-70`/`-71`, `-80` … `-85`/`-86`)
are the ones that separate the equator-ward-edge rule from the tile label.
"""
const GLO30_TILES = [
    # 1x, 3600 columns: band [0, 50)
    (lat_s =   0, lon_w =   10, ncols = 3600), (lat_s =  49, lon_w =    0, ncols = 3600),
    (lat_s =  21, lon_w =   -2, ncols = 3600), (lat_s =  49, lon_w =  125, ncols = 3600),
    (lat_s =  10, lon_w =   -5, ncols = 3600), (lat_s = -14, lon_w =  -76, ncols = 3600),
    (lat_s = -23, lon_w = -177, ncols = 3600), (lat_s =  -2, lon_w =  149, ncols = 3600),
    (lat_s = -49, lon_w =   68, ncols = 3600), (lat_s = -50, lon_w =   68, ncols = 3600),
    # 1.5x, 2400 columns: band [50, 60)
    (lat_s =  50, lon_w =    0, ncols = 2400), (lat_s =  50, lon_w =    6, ncols = 2400),
    (lat_s =  59, lon_w =    4, ncols = 2400), (lat_s =  51, lon_w =   69, ncols = 2400),
    (lat_s =  59, lon_w =   84, ncols = 2400), (lat_s =  51, lon_w =  -10, ncols = 2400),
    (lat_s = -51, lon_w =   68, ncols = 2400), (lat_s = -54, lon_w =  -68, ncols = 2400),
    (lat_s = -56, lon_w =  -68, ncols = 2400), (lat_s = -51, lon_w =  -71, ncols = 2400),
    (lat_s = -59, lon_w =  -27, ncols = 2400), (lat_s = -60, lon_w =  -27, ncols = 2400),
    # 2x, 1800 columns: band [60, 70)
    (lat_s =  60, lon_w =    4, ncols = 1800), (lat_s =  69, lon_w =   14, ncols = 1800),
    (lat_s =  66, lon_w = -163, ncols = 1800), (lat_s =  63, lon_w =   47, ncols = 1800),
    (lat_s =  60, lon_w =  157, ncols = 1800), (lat_s = -61, lon_w =  -45, ncols = 1800),
    (lat_s = -67, lon_w =   54, ncols = 1800), (lat_s = -69, lon_w =   32, ncols = 1800),
    (lat_s = -70, lon_w =    0, ncols = 1800), (lat_s = -70, lon_w =   30, ncols = 1800),
    (lat_s = -70, lon_w =   12, ncols = 1800),
    # 3x, 1200 columns: band [70, 80)
    (lat_s =  70, lon_w =   18, ncols = 1200), (lat_s =  79, lon_w =   10, ncols = 1200),
    (lat_s =  70, lon_w =  -34, ncols = 1200), (lat_s =  72, lon_w =   84, ncols = 1200),
    (lat_s =  70, lon_w =  -86, ncols = 1200), (lat_s = -71, lon_w =    0, ncols = 1200),
    (lat_s = -79, lon_w =    0, ncols = 1200), (lat_s = -80, lon_w =    0, ncols = 1200),
    (lat_s = -79, lon_w =  129, ncols = 1200), (lat_s = -77, lon_w = -111, ncols = 1200),
    (lat_s = -72, lon_w =   44, ncols = 1200),
    # 5x, 720 columns: band [80, 85)
    (lat_s =  80, lon_w =   12, ncols =  720), (lat_s =  83, lon_w =  -25, ncols =  720),
    (lat_s =  82, lon_w =  -22, ncols =  720), (lat_s =  80, lon_w =   91, ncols =  720),
    (lat_s =  80, lon_w =  -55, ncols =  720), (lat_s = -84, lon_w =    0, ncols =  720),
    (lat_s = -85, lon_w =    0, ncols =  720), (lat_s = -84, lon_w =  -32, ncols =  720),
    (lat_s = -84, lon_w =  -25, ncols =  720), (lat_s = -84, lon_w =  113, ncols =  720),
    # 10x, 360 columns: band [85, 90)
    (lat_s = -86, lon_w =    0, ncols =  360), (lat_s = -89, lon_w =    0, ncols =  360),
    (lat_s = -90, lon_w =    0, ncols =  360), (lat_s = -90, lon_w = -180, ncols =  360),
    (lat_s = -86, lon_w =  126, ncols =  360), (lat_s = -89, lon_w =  101, ncols =  360),
    (lat_s = -89, lon_w =  119, ncols =  360),
]

"""
Measured GLO-90 (`N = 1200`) tiles, same shape as [`GLO30_TILES`](@ref).

`S85_00_E000_00` appears twice because it was measured twice, in section 9.3 and
again in section 9.5, and both readings are kept verbatim: 240 columns each.
"""
const GLO90_TILES = [
    (lat_s =   0, lon_w =   10, ncols = 1200), (lat_s =  50, lon_w =    6, ncols =  800),
    (lat_s =  49, lon_w =    0, ncols = 1200), (lat_s =  60, lon_w =    4, ncols =  600),
    (lat_s =  70, lon_w =   18, ncols =  400), (lat_s =  80, lon_w =   12, ncols =  240),
    (lat_s =  83, lon_w =  -25, ncols =  240), (lat_s = -49, lon_w =   68, ncols = 1200),
    (lat_s = -50, lon_w =   68, ncols = 1200), (lat_s = -51, lon_w =   68, ncols =  800),
    (lat_s = -60, lon_w =  -27, ncols =  800), (lat_s = -61, lon_w =  -45, ncols =  600),
    (lat_s = -70, lon_w =    0, ncols =  600), (lat_s = -71, lon_w =    0, ncols =  400),
    (lat_s = -85, lon_w =    0, ncols =  240), (lat_s = -86, lon_w =    0, ncols =  120),
    (lat_s = -90, lon_w =    0, ncols =  120), (lat_s = -90, lon_w = -180, ncols =  120),
    (lat_s = -85, lon_w =    0, ncols =  240),
]

"""
Six measured GDAL geotransforms, at the full precision `ArchGDAL` printed them.

A geotransform is `(origin_x, Δlon, 0, origin_y, 0, -Δlat)`; the origin is the
**outer corner** of the first pixel, which for these `AREA_OR_POINT=Point` rasters
is half a pixel WEST and half a pixel NORTH of the first pixel centre. `Δlon` is
recorded in arcseconds, exactly as the probe printed it, so the conversion to
degrees is visible in the test rather than pre-baked here; `Δlat` is `3600 / N`
arcseconds — 1″ at GLO-30, 3″ at GLO-90 — for every tile on Earth.
"""
const GEOTRANSFORMS = [
    (N = 3600, lat_s =   0, lon_w =   10,
     origin_x =    9.99986111111111, origin_y =   1.000138888888889, dlon_arcsec =  1.0),
    (N = 3600, lat_s =  50, lon_w =    6,
     origin_x =    5.999791666666667, origin_y = 51.00013888888889,  dlon_arcsec =  1.5),
    (N = 3600, lat_s =  80, lon_w =   12,
     origin_x =   11.999305555555555, origin_y = 81.00013888888888,  dlon_arcsec =  5.0),
    (N = 3600, lat_s = -90, lon_w = -180,
     origin_x = -180.0013888888889,  origin_y = -88.99986111111112,  dlon_arcsec = 10.0),
    (N = 1200, lat_s =  50, lon_w =    6,
     origin_x =    5.999375,         origin_y =  51.000416666666666, dlon_arcsec =  4.5),
    (N = 1200, lat_s = -90, lon_w = -180,
     origin_x = -180.00416666666666, origin_y = -88.99958333333333,  dlon_arcsec = 30.0),
]
