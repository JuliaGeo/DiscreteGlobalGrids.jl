# # Multi-order storage
#
# Astronomers store sky regions as **Multi-Order Coverage** maps: HEALPix
# cells at mixed orders, coarse inside, fine along the edge. The multi-order
# coverage page builds that as a *region*; this page attaches *data* — one
# value per cell, at mixed levels, as in adaptive mesh refinement.
#
# A smooth field needs leaf cells only where it varies, so a mixed-level
# container can hold it to a stated tolerance in a fraction of the leaf count
# and still answer every leaf query.
#
# The plan: sample a real field onto a full HEALPix level, coarsen it
# adaptively, then read the result back at full resolution.

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
# A lon/lat pixel is widest at the equator, so the equatorial pixel is the
# coarsest the source gets. HEALPix cells are equal-area (`cell_area` is
# `4π / (12 * 4^level)` steradians); "across" below is √area, as on the other
# pages.

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

# Level 8 is the finest level whose cells are still wider than the equatorial
# pixel, so sampling onto it cannot invent detail, and its leaf array is under
# a million values.

L = 8
grid = DGG.levelgrid(DGG.HEALPixSystem(), L)
cells = DGG.CellVector(grid)

# `CellVector(grid)` is the whole level as a lazy vector of cell ids, one per
# position — on HEALPix, the nested order astronomy ships all-sky maps in.

# ## Sampling onto the cells
#
# One value per cell, taken at the cell centroid by nearest-neighbour lookup
# into the raster. The regridding page has the conservative, area-exact
# alternative; everything below works the same either way.

tolonlat = GO.UnitSpherical.GeographicFromUnitSphere()
centers = [tolonlat(DGG.cell_centroid(grid, c)) for c in cells]
tavg = [raster[X(Near(lon)), Y(Near(lat))] for (lon, lat) in centers]

# `Near` propagates the raster's `missing`s, so ocean cells arrive as
# `missing` rather than as a sentinel.

(; cells = length(tavg), land = count(!ismissing, tavg),
   sea = count(ismissing, tavg))

# ## The leaf array
#
# `CellLookup(grid)` reads a whole level as a cell axis, `Cells` is the
# DimensionalData dimension, and the result is an ordinary `DimArray`.

A = DD.DimArray(Vector{Union{Float64, Missing}}(tavg), DGG.Cells(DGG.CellLookup(grid)))

# It indexes by position like any vector, and by point through the axis.
# (`DD.Contains` is DimensionalData's point selector, reached through `DD`
# because this package exports DE9IM's `Contains`, a different thing.)

A[DGG.Cells(DD.Contains((2.35, 48.86)))]  # July mean over Paris

# Equal-area cells make the plain mean over land cells the area-weighted land
# mean.

(; land_mean = mean(skipmissing(A)), range = extrema(skipmissing(A)))

# ## The leaf data, drawn
#
# `cell_polygon` gives each cell's unit-sphere polygon and `GO.transform`
# takes the vector to lon/lat. Cells straddling the antimeridian would paint
# full-width bands in planar longitude, so they are dropped from the figure.

polys = GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygon.(Ref(grid), cells))
lonspan(p) = (ring = GI.coordinates(p)[1]; maximum(first, ring) - minimum(first, ring))
onmap = findall(p -> lonspan(p) < 180, polys)
length(cells) - length(onmap)

# The colour scale clips at ±40 °C so midwinter Antarctica does not flatten
# the rest of the range; ocean cells are `missing`, drawn in gray. The colour
# vector is built over `parent(A)`: a comprehension over the `DimArray`
# returns a `DimArray`, which Makie would read as a per-vertex colouring.

colors = [ismissing(v) ? NaN : v for v in parent(A)]

fig = Figure(size = (900, 470))
ax = GeoAxis(fig[1, 1]; dest = "+proj=moll",
    title = "WorldClim July mean temperature, HEALPix level $L")
plt = poly!(ax, polys[onmap]; color = colors[onmap], colormap = Reverse(:RdYlBu),
    colorrange = (-40, 40), highclip = :darkred, lowclip = :midnightblue,
    nan_color = :gray85, strokewidth = 0)
Colorbar(fig[1, 2], plt; label = "mean temperature (°C)")
fig

# ## A fixed-level pyramid
#
# `aggregate(f, A, l)` gives one value per level-`l` ancestor of the array's
# cells, `f` seeing the leaf values beneath it; one call per level builds a
# pyramid.
#
# The reducer must define the all-ocean case: `mean` over a group holding any
# `missing` answers `missing`, and `mean(skipmissing(v))` over an all-`missing`
# group answers `NaN`. Neither means "no land here", so the reducer says it
# explicitly.

landmean(v) = count(!ismissing, v) == 0 ? missing : mean(skipmissing(v))

pyramid = [DGG.aggregate(landmean, A, l) for l in (6, 5, 4)]
[(; level = DGG.level(DD.lookup(P, DGG.Cells)), cells = length(P),
    land = count(!ismissing, parent(P))) for P in pyramid]

# Polygons and the antimeridian drop are recomputed per cell set — a coarse
# cell can straddle the line where none of its children did.

onmapof(ps) = findall(p -> lonspan(p) < 180, ps)
cellpolys(lk) = GO.transform(GO.GeographicFromUnitSphere(),
                             collect(DGG.getcell(DGG.PartialGrid(lk))))

# One cell's width in kilometres at its own level — the √area convention
# again, now that levels are mixed.
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

# A fixed level stores a value for every cell of that level, ocean included,
# and applies one resolution to the whole globe. A spread reducer shows what
# that costs:

worst(v) = (m = collect(skipmissing(v)); isempty(m) ? missing : maximum(m) - minimum(m))
for l in 7:-1:4
    S = DGG.aggregate(worst, A, l)
    println("level ", l, ": worst cell replaces a spread of ",
            round(maximum(skipmissing(parent(S))); digits = 1), " °C")
end

# A level-6 cell on the Antarctic coast, 100 km across, replaces a spread of
# tens of degrees; one in the Sahara replaces a fraction of a degree. A fixed
# level cannot treat the two differently.

# ## Adaptive desampling
#
# `coarsen` merges a sibling group when the group is complete in the array and
# its leaf values span at most `atol`, or when every value is `missing`. A
# group mixing `missing` with data never merges, so coastline cells are never
# averaged into the sea. The criterion is monotone, so each leaf ends up under
# the coarsest cell that satisfies it.
#
# `atol` is the caller's choice; a sweep shows what each costs:

for atol in (0.25, 0.5, 1.0, 2.0, 4.0)
    Mt = DGG.coarsen(A; atol)
    println("atol = ", rpad(atol, 4), " °C: ", lpad(length(Mt), 6), " cells stored, ",
            lpad(round(100 * length(Mt) / length(A); digits = 1), 4), "% of the leaf array")
end

# One degree is the granularity a monthly climatology is read at, and far
# below the field's range of over a hundred degrees; the rest of the page uses
# it.

M = DGG.coarsen(A; atol = 1.0)
mov = parent(DD.lookup(M, DGG.Cells))

# The result is a `DimArray` over a `MultiOrderLookup`: one value per stored
# cell, each at whatever level merging stopped. The lookup wraps a
# `MultiOrderVector` — `mov` above — a plain vector of cells usable with no
# datacube library; the lookup only makes it an axis. The rest of the page
# passes `mov` directly to `cell_polygons`, `cellat` and `summarysize`.

levs = DGG.level.(mov)
for l in sort(unique(levs))
    at = findall(==(l), levs)
    println("level ", l, ": ", lpad(length(at), 7), " cells, ",
            lpad(round(Int, across_km(mov[first(at)])), 5), " km across, ",
            lpad(count(ismissing, parent(M)[at]), 6), " ocean, ",
            lpad(count(!ismissing, parent(M)[at]), 6), " land")
end

# `minlevel` defaults to the shallowest level the system has; merging reached
# level 2 on its own, on 34 open-ocean cells 1600 km across. Almost everything
# coarser than level 6 is ocean; the one land cell at level 4 is in the
# western Amazon, where the July mean is flat over 165,000 km².
#
# Drawn twice: by value, which should match the leaf map, and by cell level,
# which is the mesh itself.

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

# The ocean is a few huge cells, flat interiors are large ones, and the leaf
# level survives along coastlines, mountain ranges and the Antarctic margin.
# By value, the field is unchanged to the eye.
#
# The cost in cells and bytes:

(; stored = length(M), leaves = length(A),
   values_MiB = round(Base.summarysize(parent(M)) / 2^20; digits = 2),
   mesh_MiB = round(Base.summarysize(mov) / 2^20; digits = 2),
   dense_MiB = round(Base.summarysize(parent(A)) / 2^20; digits = 2))

# The values shrink by the cell ratio, a little over five times. The index
# does not: a cell id and three interval bounds per stored cell — some forty
# bytes — gives most of that saving back at this tolerance.

round(Int, Base.summarysize(mov) / length(M))  # bytes of index per stored cell

# ## Reading it back
#
# `Contains` answers a point with the stored cell that covers it, at whatever
# level that cell is. `At` selects an exact cell and throws if it is absent.

for (name, lon, lat) in [("Lisbon", -9.14, 38.72), ("Denver", -104.99, 39.74),
                         ("central Sahara", 15.0, 24.0), ("Point Nemo", -123.4, -48.9)]
    v = M[DGG.Cells(DD.Contains((lon, lat)))]
    c = DGG.cellat(mov, lon, lat)
    across = sqrt(DGG.cell_area(DGG.levelgrid(DGG.HEALPixSystem(), DGG.level(c)), c) * R^2)
    println(rpad(name, 15), ismissing(v) ? "     missing" : lpad(round(v; digits = 1), 9) * " °C",
            "  from a level-", DGG.level(c), " cell, ", round(Int, across), " km across")
end

# Lisbon, on a coastline, answers from a leaf cell; the Sahara and the Rockies
# answer from a cell four times the area; the South Pacific answers `missing`
# from a cell 1600 km across.
#
# `expand` presents the whole array at one level again:

E = DGG.expand(M, L)
DD.lookup(E, DGG.Cells) == DD.lookup(A, DGG.Cells)

# The same axis as `A`, cell for cell. The storage is unchanged — one value
# per stored cell plus a leaf count each, indexed by binary search — so the
# presentation level moves while the memory does not.

for l in (L, L + 1, L + 2, L + 4)
    El = DGG.expand(M, l)
    dense = Base.summarysize(parent(A)) / length(A) * length(El)
    println("presented at level ", l, ": ", lpad(length(El), 10), " cells, ",
            lpad(round(Base.summarysize(parent(El)) / 2^20; digits = 2), 6),
            " MiB stored, ", lpad(round(dense / 2^20; digits = 1), 8), " MiB if dense")
end

# At the leaf level the expanded array is already smaller than the dense one;
# four levels down it presents two hundred million cells — over a gigabyte if
# dense — from the same stored bytes.
#
# `coarsen`'s criterion bounds the error: a stored value lies between the
# minimum and maximum of the leaves it replaces, which differ by at most
# `atol`. The leaf-by-leaf check:

err = abs.(collect(parent(E)) .- parent(A))
(; atol = 1.0, max_error = maximum(skipmissing(err)),
   mean_error = mean(skipmissing(err)),
   unchanged = count(iszero, skipmissing(err)),
   missing_agrees = all(ismissing.(collect(parent(E))) .== ismissing.(parent(A))))

# `missing` round-trips exactly: an all-`missing` group merges to `missing`
# and a mixed group never merges. Two fifths of the land cells are
# bit-identical — the ones kept at the leaf level. Mapping the rest, one worst
# case per level-6 cell:

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

# The error is largest where the field has structure at the merge scale — the
# Andes, the Himalayan front, the Rift, the Antarctic margin. The black zeros
# are cells kept at the leaf level: unmerged cells are copied, not
# approximated.

# ## Summary
#
# `coarsen` picks a level per cell from the data and a tolerance; the axis
# answers points and regions on the mixed levels directly; `expand` presents
# the result at any level without materialising it. Storage scales with stored
# cells and resolution with leaves: the same two megabytes present 786,432
# cells or two hundred million.
#
# The multi-order coverage page builds the same structure as a *region* — the
# query side. The regridding page replaces this page's nearest-neighbour
# sampling with conservative, area-exact regridding.
