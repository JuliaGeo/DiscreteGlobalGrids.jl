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

# ## The obvious compression: a pyramid
#
# Fix a coarser level and average into it. `aggregate(f, A, l)` gives one value
# per level-`l` ancestor of the array's cells, `f` seeing the leaf values
# beneath it — the pyramid every tiled format builds, one call per level and no
# type of its own.
#
# The reducer has to say what an all-ocean group means. `mean` over a group
# holding one `missing` answers `missing`, and `mean(skipmissing(v))` over a
# group holding nothing but `missing`s answers `NaN`, an empty sum over an empty
# count. Neither reads as "no land here", so the reducer says it outright.

landmean(v) = count(!ismissing, v) == 0 ? missing : mean(skipmissing(v))

pyramid = [DGG.aggregate(landmean, A, l) for l in (6, 5, 4)]
[(; level = DGG.level(DD.lookup(P, DGG.Cells)), cells = length(P),
    land = count(!ismissing, parent(P))) for P in pyramid]

# Each level is drawn the way the leaf map was: polygons from the axis, then the
# same antimeridian drop. It has to be recomputed per cell set — a coarse cell
# straddles the line where none of its children did.

onmapof(ps) = findall(p -> lonspan(p) < 180, ps)
cellpolys(lk) = GO.transform(GO.GeographicFromUnitSphere(),
                             collect(DGG.getcell(DGG.PartialGrid(lk))))

# One cell's width, in kilometres, at whatever level it belongs to — the √area
# convention of the level table above, reused now that levels are mixed.
across_km(c) = sqrt(DGG.cell_area(DGG.levelgrid(DGG.HEALPixSystem(), DGG.level(c)), c) * R^2)

fig = Figure(size = (1000, 300))
for (k, P) in enumerate(pyramid)
    local lk = DD.lookup(P, DGG.Cells)
    local across = across_km(lk[1])
    local ps = cellpolys(lk)
    local keep = onmapof(ps)
    local panel = GeoAxis(fig[1, k]; dest = "+proj=moll",
        title = "level $(DGG.level(lk)), $(round(Int, across)) km across")
    poly!(panel, ps[keep];
        color = [ismissing(v) ? NaN : v for v in parent(P)][keep],
        colormap = Reverse(:RdYlBu), colorrange = (-40, 40),
        highclip = :darkred, lowclip = :midnightblue, nan_color = :gray85,
        strokewidth = 0)
end
fig

# Every panel stores a value for every cell of its level, ocean included, and
# the level is chosen for the whole globe at once — which is the trade. The same
# call, with a reducer that reports the spread instead of the mean, says what
# each fixed level costs where it costs most:

worst(v) = (m = collect(skipmissing(v)); isempty(m) ? missing : maximum(m) - minimum(m))
for l in 7:-1:4
    S = DGG.aggregate(worst, A, l)
    println("level ", l, ": worst cell replaces a spread of ",
            round(maximum(skipmissing(parent(S))); digits = 1), " °C")
end

# A level-6 cell on the Antarctic coast, 100 km across, stands in for tens of
# degrees; a level-6 cell in the middle of the Sahara stands in for a fraction
# of one. A fixed level cannot tell the two apart, and picking a level is
# picking which of them to be wrong about.

# ## Adaptive desampling
#
# `coarsen` asks the data instead. A cell replaces its subtree when the subtree
# is complete in the array and its leaf values span no more than `atol` — or
# when every one of them is `missing`, which is the ocean. A group mixing
# `missing` with data never merges, so a coastline is never averaged into the
# sea. The criterion is monotone, so what comes back is the coarsest cell that
# can speak for each leaf.
#
# Nothing in the field picks `atol`; the reader does. Here is what each choice
# costs:

for atol in (0.25, 0.5, 1.0, 2.0, 4.0)
    Mt = DGG.coarsen(A; atol)
    println("atol = ", rpad(atol, 4), " °C: ", lpad(length(Mt), 6), " cells stored, ",
            lpad(round(100 * length(Mt) / length(A); digits = 1), 4), "% of the leaf array")
end

# One degree is the granularity a monthly mean climatology is read at, and it is
# far below the field's own contrast — the land range spans over a hundred
# degrees. Take it.

M = DGG.coarsen(A; atol = 1.0)
mov = parent(DD.lookup(M, DGG.Cells))

# The result is a `DimArray` again, over a `MultiOrderLookup` this time: one
# value per stored cell, and the cells are at whatever level the data stopped
# at. Where they stopped is the whole story.

levs = DGG.level.(mov)
for l in sort(unique(levs))
    at = findall(==(l), levs)
    println("level ", l, ": ", lpad(length(at), 7), " cells, ",
            lpad(round(Int, across_km(mov[first(at)])), 5), " km across, ",
            lpad(count(ismissing, parent(M)[at]), 6), " ocean, ",
            lpad(count(!ismissing, parent(M)[at]), 6), " land")
end

# Nothing asked for level 2 — `minlevel` defaults to the shallowest level the
# system has, and the merging simply stopped there, on 34 cells 1600 km across
# in the open ocean. Almost everything coarser than level 6 is ocean; the single
# land cell that reaches level 4 sits on the equator in the western Amazon,
# where a July mean is flat over 165,000 km².
#
# Drawn twice: coloured by value, which should look like the leaf map, and
# coloured by cell level, which is the adaptive mesh itself.

mpolys = GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygons(mov))
mkeep = onmapof(mpolys)

fig = Figure(size = (1000, 340))
axv = GeoAxis(fig[1, 1]; dest = "+proj=moll", title = "value, atol = 1 °C")
pv = poly!(axv, mpolys[mkeep];
    color = [ismissing(v) ? NaN : v for v in parent(M)][mkeep],
    colormap = Reverse(:RdYlBu), colorrange = (-40, 40),
    highclip = :darkred, lowclip = :midnightblue, nan_color = :gray85, strokewidth = 0)
Colorbar(fig[2, 1], pv; vertical = false, label = "mean temperature (°C)")
axl = GeoAxis(fig[1, 2]; dest = "+proj=moll", title = "cell level")
pl = poly!(axl, mpolys[mkeep]; color = levs[mkeep],
    colormap = :isoluminant_cgo_70_c39_n256, colorrange = extrema(levs), strokewidth = 0)
Colorbar(fig[2, 2], pl; vertical = false, ticks = extrema(levs)[1]:extrema(levs)[2])
fig

# The right-hand map is the compression, stated: the ocean is a handful of huge
# cells, the flat continental interiors are large ones, and the leaf level
# survives exactly along the coastlines, the mountain ranges and the Antarctic
# margin — everywhere a hundred kilometres of ground is worth more than one
# number. The left-hand map is the field, unchanged to the eye.
#
# What that bought, in cells and in bytes:

(; stored = length(M), leaves = length(A),
   values_MiB = round(Base.summarysize(parent(M)) / 2^20; digits = 2),
   mesh_MiB = round(Base.summarysize(mov) / 2^20; digits = 2),
   dense_MiB = round(Base.summarysize(parent(A)) / 2^20; digits = 2))

# The values shrink by the cell ratio, a little over five times. The mesh that
# makes them addressable does not: it is a cell id and three interval bounds per
# stored cell, some forty bytes, and at this tolerance it gives most of that
# saving back. That is the honest accounting of a container built for lookup.
# The next section is where the arithmetic changes sign.

round(Int, Base.summarysize(mov) / length(M))  # bytes of index per stored cell

# ## Reading it back
#
# The mixed-level axis answers a point the way the leaf array did, and answers
# it with the cell that actually covers the point — which is what `Contains`
# means here. `At` is the other reading, the exact cell, and it raises rather
# than guessing.

for (name, lon, lat) in [("Lisbon", -9.14, 38.72), ("Denver", -104.99, 39.74),
                         ("central Sahara", 15.0, 24.0), ("Point Nemo", -123.4, -48.9)]
    v = M[DGG.Cells(DD.Contains((lon, lat)))]
    c = DGG.cellat(mov, lon, lat)
    across = sqrt(DGG.cell_area(DGG.levelgrid(DGG.HEALPixSystem(), DGG.level(c)), c) * R^2)
    println(rpad(name, 15), ismissing(v) ? "     missing" : lpad(round(v; digits = 1), 9) * " °C",
            "  from a level-", DGG.level(c), " cell, ", round(Int, across), " km across")
end

# Lisbon is on a coastline and gets a leaf cell; the Sahara and the Rockies are
# smooth enough to answer from a cell four times the area; the South Pacific
# answers `missing` from a cell 1600 km across. One array, four resolutions, one
# selector.
#
# `expand` goes the other way — it presents the whole thing at one level again:

E = DGG.expand(M, L)
DD.lookup(E, DGG.Cells) == DD.lookup(A, DGG.Cells)

# The same axis as `A`, cell for cell. The data is not the same, though: it is
# still one value per stored cell, plus one leaf count each, and the reindexing
# is a binary search. So the presentation level moves and the memory does not.

for l in (L, L + 1, L + 2, L + 4)
    El = DGG.expand(M, l)
    dense = Base.summarysize(parent(A)) / length(A) * length(El)
    println("presented at level ", l, ": ", lpad(length(El), 10), " cells, ",
            lpad(round(Base.summarysize(parent(El)) / 2^20; digits = 2), 6),
            " MiB stored, ", lpad(round(dense / 2^20; digits = 1), 8), " MiB if dense")
end

# That is where the sign changes. At the leaf level the expanded array already
# beats the dense one; four levels down it names two hundred million cells, a
# dense array of them would be well over a gigabyte, and the stored bytes have
# not moved at all. The mesh is paid for once and answers at every resolution.
#
# And it is within tolerance. `coarsen`'s criterion bounds the error outright —
# the stored value lies between the minimum and maximum of the leaves it
# replaces, and those differ by at most `atol` — so the leaf-by-leaf difference
# is a check on the implementation, not a hope:

err = abs.(collect(parent(E)) .- parent(A))
(; atol = 1.0, max_error = maximum(skipmissing(err)),
   mean_error = mean(skipmissing(err)),
   unchanged = count(iszero, skipmissing(err)),
   missing_agrees = all(ismissing.(collect(parent(E))) .== ismissing.(parent(A))))

# The ocean comes back `missing` in exactly the cells it went in as `missing`,
# because an all-`missing` group merges to `missing` and a mixed group does not
# merge at all. Two fifths of the land cells are bit-identical: those are the
# ones the mesh kept at the leaf level. Where the rest of the error lives is
# worth a map — aggregated to level 6 so it is one worst case per 100 km, using
# the same `aggregate` the pyramid did:

worsterr = DGG.aggregate(v -> (m = collect(skipmissing(v)); isempty(m) ? missing : maximum(m)),
                         DD.DimArray(err, DGG.Cells(DGG.CellLookup(grid))), 6)
epolys = cellpolys(DD.lookup(worsterr, DGG.Cells))
ekeep = onmapof(epolys)

fig = Figure(size = (900, 450))
ax = GeoAxis(fig[1, 1]; dest = "+proj=moll",
    title = "worst leaf error per 100 km cell, atol = 1 °C")
plt = poly!(ax, epolys[ekeep];
    color = [ismissing(v) ? NaN : v for v in parent(worsterr)][ekeep],
    colormap = :magma, colorrange = (0, 1.0), nan_color = :gray85, strokewidth = 0)
Colorbar(fig[1, 2], plt; label = "|expanded − leaf| (°C)")
fig

# The tolerance is spent where the field has structure at the merge scale —
# highest along the Andes, the Himalayan front, the Rift, the Antarctic margin
# — and moderately over the flat interiors their merged cells average across.
# The black zeros are the cells the mesh refused to merge at all: a coastline
# kept at the leaf level is not approximated, it is copied.

# ## What the container buys
#
# One array, mixed levels. `coarsen` chose a level per cell from the data and a
# tolerance rather than from a global guess; the axis answers points and regions
# on the mixed levels directly; `expand` presents the result at any level at all
# without materialising it. The cost is per *stored* cell and the resolution the
# array can answer at is per *leaf*, which is the whole trick: the same two
# megabytes speak for 786,432 cells or for two hundred million.
#
# The multi-order coverage page builds the same mixed-level structure as a
# *region* — cells chosen by an outline rather than by a field — and is where
# the query side lives. The regridding page replaces this page's
# nearest-neighbour sampling with the conservative, area-exact kind, which is
# what to reach for when the leaf values must be true areal means before
# anything compresses them.
