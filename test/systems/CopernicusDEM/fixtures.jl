# ---------------------------------------------------------------------------
# Committed ArchGDAL measurements from AWS Open Data COGs; tests use no network.
#
# The boundary samples establish the southern convention: `S50` is 1x and `S51`
# is 1.5x, a distinction absent from the handbooks.
# ---------------------------------------------------------------------------

"""
Measured GLO-30 tile corners and column counts, grouped by band. Southern boundary
rows pin the equatorward-edge rule.
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
Measured GLO-90 tiles in the same format. `S85_00_E000_00` appears twice because
both measurements are retained.
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
Measured GDAL geotransforms at ArchGDAL's printed precision. The origin is the
first pixel's outer corner; `Δlon` is in arcseconds.
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
