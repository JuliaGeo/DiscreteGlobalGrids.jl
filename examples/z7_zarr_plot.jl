# Demo: mapping a variable out of an IGEO7/Z7 DGGS-convention Zarr archive.
#
# The companion to `z7_zarr_read.jl`, which covers opening and selecting. This
# one takes the last step — draw it — and is separate because it is the only
# example that needs Makie.
#
# Every cell is rendered as the hexagon it actually is, from
# `DGGSZarr.cell_boundaries`, rather than as a scattered point. That matters
# beyond looks: the boundaries carry the archive's declared icosahedron
# placement (`dggs_vert0_lon = 11.2`) and latitude datum
# (`igeo7_wgs84_geodetic_conversion`), so the hexagons land where the data
# actually is. Drawn with `IGeo7.cell_boundary`'s defaults instead they would
# sit ~13 km north and 0.05° east, tessellating just as convincingly.
#
# Environment: unlike the other examples this one does **not** run in the
# package environment — CairoMakie and GeoMakie are not package dependencies.
# Use the docs environment, which already has them:
#
#     julia -t 4 --project=docs examples/z7_zarr_plot.jl
#
# Set DGGS_ZARR_TEST_DATA to the directory holding the archives if they are not
# in the default location, and DGGS_EXAMPLE_OUTPUT to choose where the PNG goes
# (default: the current directory).
using DiscreteGlobalGrids
using DiscreteGlobalGrids.DGGSZarr
using DiscreteGlobalGrids: IGeo7
using DimensionalData
using CairoMakie
using GeoMakie
import GeoInterface as GI

const ARCHIVES = get(ENV, "DGGS_ZARR_TEST_DATA",
    joinpath(homedir(), "dev", "build", "igeo7_z7_xarray_paper", "data", "working"))
const OUTPUT = get(ENV, "DGGS_EXAMPLE_OUTPUT", pwd())

# res 10 (~3.1k cells) rather than res 12 (~158k): both archives hold the same
# DEM, and at res 10 the hexagons are still individually visible and Cairo draws
# them in about a second. The script is resolution-agnostic otherwise — point it
# at the res-12 archive and it will draw that instead, more slowly.
# const ARCHIVE = joinpath(ARCHIVES, "pori_z7_r10_ranges.zarr")
const ARCHIVE = joinpath(ARCHIVES, "pori_z7_r12_ranges.zarr")

println("="^78)
println("z7_zarr_plot.jl — mapping an IGEO7/Z7 Zarr archive")
println("julia $(VERSION)  threads=$(Threads.nthreads())")
println("="^78)

if !isdir(ARCHIVE)
    println("\narchive not found at $ARCHIVE")
    println("set DGGS_ZARR_TEST_DATA to the directory holding pori_z7_r10_ranges.zarr")
    exit(0)
end

# --------------------------------------------------------------------------
# 1. Open, lazily, and take only the two variables we are going to draw.
# --------------------------------------------------------------------------

ds = open_dggs_dataset(ARCHIVE)
info = dggs_info(ds)
ids = dggs_cell_ids(ds)

println("\narchive   : ", basename(ARCHIVE))
println("grid      : ", info.name, " level ", info.level,
    " (compression \"", info.compression, "\")")
println("cells     : ", length(ids), " in ", IGeo7.z7_nranges(ids), " ranges")
println("variables : ", join(string.(keys(ds.cubes)), ", "))

"""
    layer(name) -> (values, label)

Materialize one variable of the archive as a plain `Vector`, with a label built
from its `long_name`/`units` attributes when the writer supplied them.

This read is the moment the archive stops being lazy — everything before it was
metadata and arithmetic. Note the `Array`: slicing a `YAXArray` gives back a
`YAXArray` (still carrying the DGGS dimension), and Makie wants the bare
numbers.
"""
function layer(name::Symbol)
    cube = ds.cubes[name]
    values = Array(cube[:])
    attrs = cube.properties
    long_name = get(attrs, "long_name", String(name))
    units = get(attrs, "units", "")
    return values, isempty(units) ? long_name : "$long_name ($units)"
end

# --------------------------------------------------------------------------
# 2. Cell geometry: one hexagon per cell, in the archive's own frame.
#
# `cell_boundaries` returns closed rings (the first corner repeated), which
# Makie tolerates; the repeat is dropped anyway so the polygons carry exactly
# their six distinct corners.
# --------------------------------------------------------------------------

println("\nbuilding ", length(ids), " cell polygons ...")
hexagons = map(cell_boundaries(ds)) do polygon
    ring = first(GI.coordinates(polygon))
    Makie.Polygon(Point2f.(ring[1:(end-1)]))
end

lons = [p[1] for hex in hexagons for p in hex.exterior]
lats = [p[2] for hex in hexagons for p in hex.exterior]
lon_range, lat_range = extrema(lons), extrema(lats)
pad = 0.005
extent = ((lon_range[1] - pad, lon_range[2] + pad),
          (lat_range[1] - pad, lat_range[2] + pad))
println("extent    : lon ", round.(lon_range; digits = 3),
    "  lat ", round.(lat_range; digits = 3))

# Ticks on a round 0.1° grid rather than wherever the auto-locator lands.
#
# 0.1° and not finer, deliberately: `GeoAxis` builds its own degree labels with
# `round(x; sigdigits = 3)` (GeoMakie `geoaxis.jl`) and ignores `xtickformat` /
# `ytickformat`, so on an extent this small a 0.05° step would print 58.15 and
# 58.25 both as "58.2°". Three distinct latitude labels beat five with three of
# them identical; don't reach for `ytickformat` to fix it, it is not wired up.
degree_ticks(lo, hi, step = 0.1) = round(lo / step) * step : step : hi

# --------------------------------------------------------------------------
# 3. Draw. One panel per variable, sharing the geographic frame.
# --------------------------------------------------------------------------

"""
    panel!(figure, column, values, label; colormap)

Draw one variable into `figure[1, column]` as coloured hexagons on a `GeoAxis`,
with its own colourbar underneath.
"""
function panel!(figure, column::Int, values, label::AbstractString; colormap)
    axis = GeoAxis(figure[1, column];
        title = label,
        limits = extent,
        xticks = degree_ticks(extent[1]...),
        yticks = degree_ticks(extent[2]...),
        xticklabelsize = 9,
        yticklabelsize = 9,
    )
    finite = filter(isfinite, values)
    plot = poly!(axis, hexagons;
        color = values,
        colormap,
        colorrange = extrema(finite),
        strokewidth = 0,
    )
    Colorbar(figure[2, column], plot; vertical = false, flipaxis = false,
        height = 10, ticklabelsize = 9)
    return axis
end

elevation, elevation_label = layer(:elevation)
slope, slope_label = layer(:slope_lookup)

# The maps are wider than tall (0.43° by 0.22°, and a degree of longitude is
# about half a degree of latitude at 58°N), and a GeoAxis holds that aspect —
# so the figure has to be sized to match or the panels float in whitespace.
figure = Figure(size = (1200, 470), fontsize = 12)
Label(figure[0, 1:2],
    "IGEO7 / Z7 level $(info.level) — $(length(ids)) cells, " *
    "vertex 0 at $(DGGSZarr.igeo7_vert0_lon(info))°E, WGS84";
    fontsize = 15, font = :bold)

panel!(figure, 1, elevation, elevation_label; colormap = :terrain)
panel!(figure, 2, slope, slope_label; colormap = :viridis)

rowsize!(figure.layout, 2, Fixed(24))       # the colourbar row, kept near its map
rowgap!(figure.layout, 1, 4)                # title -> maps
rowgap!(figure.layout, 2, 2)                # maps -> colourbars
colgap!(figure.layout, 30)

path = joinpath(OUTPUT, "z7_pori_r$(info.level).png")
save(path, figure)

println("\n", "="^78)
println("elevation : ", round(minimum(filter(isfinite, elevation)); digits = 1), " – ",
    round(maximum(filter(isfinite, elevation)); digits = 1), " m")
println("slope     : ", round(minimum(filter(isfinite, slope)); digits = 4), " – ",
    round(maximum(filter(isfinite, slope)); digits = 4), " m/m")
println("written   : ", path)
println("="^78)
