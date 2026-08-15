# # MOC storage: one value per cell, at whatever level the data needs
#
# Astronomers store a region as a **Multi-Order Coverage** map: HEALPix cells
# at mixed orders, coarse through the interior and fine along the edge, sorted
# so a point query is a binary search. The multi-order coverage page builds one
# as a *region* — the query side. This page gives the same structure *data*:
# one value per cell, at mixed levels.
#
# Read from the other direction it is adaptive mesh refinement. Keep the mesh
# coarse where the field is flat, refine it only where the field varies, and
# store one value per cell of the result. A global temperature field is mostly
# smooth: across the Sahara or the Siberian interior, whole neighbourhoods of
# cells agree to within a degree, and over the ocean this particular field is
# not there at all. So a mixed-level container can hold the field to a stated
# tolerance in a small fraction of the cells, and still answer "what is the
# value here?" for every leaf cell.
#
# The page goes leaf-first. Sample a real field onto a full HEALPix level,
# desample it adaptively, then read the compressed result back as if it were
# still full resolution. This first half builds the leaf data.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH",
    joinpath(tempdir(), "rasterdatasources")))

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GeometryOps as GO
import GeoInterface as GI
using Rasters, RasterDataSources
import ArchGDAL
using Statistics
using GLMakie, GeoMakie
GLMakie.activate!(inline = true)

# ## The source
#
# WorldClim's July mean temperature at 10 arc-minutes: one global raster in °C
# over land, `missing` over the ocean. RasterDataSources downloads it once and
# caches it under `RASTERDATASOURCES_PATH`.

raster = Raster(WorldClim{Climate}, :tavg; month = 7, res = "10m")
size(raster)

# ## Choosing the destination level
#
# A lon/lat pixel is widest at the equator and shrinks to nothing at the poles,
# so the equatorial pixel is the coarsest the source ever gets. A destination
# cell narrower than that resolves nothing the raster has; a cell much wider
# throws data away. HEALPix cells are exactly equal-area — `cell_area` returns
# `4π / (12 * 4^level)` steradians — so "across" below is √area, the convention
# on these pages.

R = 6371.0088                            # mean Earth radius, km
pixel_km = (10 / 60) * (2pi * R / 360)   # a 10 arc-minute pixel at the equator
println("raster: ", length(raster), " pixels, ", round(pixel_km; digits = 1),
        " km wide at the equator")
for l in 6:9
    g = DGG.levelgrid(DGG.HEALPixSystem(), l)
    across = sqrt(DGG.cell_area(g, DGG.CellVector(g)[1]) * R^2)
    println("level ", l, ": ", lpad(DGG.ncells(g), 8), " cells, ",
            lpad(round(across; digits = 1), 6), " km across, ",
            lpad(round(across / pixel_km; digits = 2), 5), " × that pixel")
end

# The last column decides it: level 8 is the finest level whose cells are still
# wider than the equatorial pixel, so sampling onto it cannot invent detail,
# and its leaf array is under a million values — small enough to hold densely,
# large enough that compressing it is worth doing. Level 7 would be defensible
# and four times cheaper to draw, at the price of discarding most of the
# raster; level 9 would be finer than the source along the equator.

L = 8
grid = DGG.levelgrid(DGG.HEALPixSystem(), L)
cells = DGG.CellVector(grid)

# `CellVector(grid)` is the whole level read as a lazy vector of cell ids, one
# per position. On HEALPix that position order is the nested order astronomy
# ships all-sky maps in; the HEALPix page checks that against Healpix.jl.

# ## Sampling onto the cells
#
# One value per cell, taken at the cell centroid with a nearest-neighbour
# lookup into the raster. This is the quick route — no geometry is clipped and
# the globe samples in a fraction of a second — and it costs little precisely
# because the level was chosen so a cell is about a pixel wide. The
# conservative alternative, which integrates every overlapping pixel and gets
# the areal mean exactly, is the subject of the regridding page. Everything
# below is indifferent to which of the two produced the values.

tolonlat = GO.UnitSpherical.GeographicFromUnitSphere()
centers = [tolonlat(DGG.cell_centroid(grid, c)) for c in cells]
tavg = [raster[X(Near(lon)), Y(Near(lat))] for (lon, lat) in centers]

# `Near` propagates the raster's `missing`s, so ocean cells arrive as `missing`
# rather than as a sentinel that later arithmetic would have to remember. That
# is the point: `missing` is the flattest field there is, and the ocean is most
# of the globe.

(; cells = length(tavg), land = count(!ismissing, tavg),
   sea = count(ismissing, tavg))

# ## The leaf array
#
# `CellLookup(grid)` reads a whole level as a cell axis, `Cells` is the
# DimensionalData dimension, and the result is an ordinary `DimArray` whose
# axis happens to be a discrete global grid.

A = DD.DimArray(Vector{Union{Float64, Missing}}(tavg), DGG.Cells(DGG.CellLookup(grid)))

# It indexes by position like any vector, and by point through the axis.
# (`DD.Contains` is DimensionalData's point selector, reached through `DD`
# because this package exports DE9IM's `Contains`, a different thing.)

A[DGG.Cells(DD.Contains((2.35, 48.86)))]  # July mean over Paris

# Equal-area cells mean the plain mean over the land cells is already the
# area-weighted land mean — no cosine weights anywhere.

(; land_mean = mean(skipmissing(A)), range = extrema(skipmissing(A)))

# ## The leaf data, drawn
#
# `cell_polygon` is the unit-sphere polygon of one cell and one `GO.transform`
# takes the whole vector to lon/lat. A thin column of cells straddles the
# antimeridian: in planar longitude their rings run the full width of the map
# and would paint bands across it, so they are dropped from the figure, as on
# the regridding page.

polys = GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygon.(Ref(grid), cells))
lonspan(p) = (ring = GI.coordinates(p)[1]; maximum(first, ring) - minimum(first, ring))
onmap = findall(p -> lonspan(p) < 180, polys)
length(cells) - length(onmap)

# Midwinter Antarctica reaches far below the rest of the range and would
# flatten every other contrast on a full-range colour scale, so the map clips
# at ±40 °C. Ocean cells are `missing`, drawn in gray. Makie reads a colour
# vector positionally, so this one is built over `parent(A)`: a comprehension
# over the `DimArray` itself would come back as a `DimArray`, and Makie would
# take that for a per-vertex colouring.

colors = [ismissing(v) ? NaN : v for v in parent(A)]

fig = Figure(size = (900, 470))
ax = GeoAxis(fig[1, 1]; dest = "+proj=moll",
    title = "WorldClim July mean temperature, HEALPix level $L")
plt = poly!(ax, polys[onmap]; color = colors[onmap], colormap = Reverse(:RdYlBu),
    colorrange = (-40, 40), highclip = :darkred, lowclip = :midnightblue,
    nan_color = :gray85, strokewidth = 0)
Colorbar(fig[1, 2], plt; label = "mean temperature (°C)")
fig

# That is the leaf array: one value per level-8 cell, every one of them stored,
# and most of them saying what a neighbour already said. The rest of the page
# asks how few of them the field actually needs.

# ## (continued in wave 3: pyramid, adaptive coarsening, compression)
