# # The sky in HEALPix
#
# HEALPix is astronomy's grid: every CMB map and all-sky survey ships as a flat
# vector of `12 * nside^2` equal-area pixels in nested order. That order is
# exactly the position order of `levelgrid(HEALPixSystem(), level)`, so an
# astronomer's map and a DiscreteGlobalGrids data vector are the same vector.
#
# This page builds a synthetic all-sky map, checks that claim against
# Healpix.jl, and runs the two classic sky operations: a cone search and a
# galactic-plane cut.

import DiscreteGlobalGrids as DGG
import Healpix
import GeometryOps as GO
using Statistics, Random
using GLMakie, GeoMakie
GLMakie.activate!(inline = true)

# ## A synthetic sky
#
# Level 5 is `nside = 2^5 = 32`: 12288 pixels. `CellVector` reads the grid as a
# lazy vector of cell ids, one per position.

grid = DGG.levelgrid(DGG.HEALPixSystem(), 5)
cells = DGG.CellVector(grid)

# The fake sky is a diffuse band along the galactic plane, four point sources,
# and a little noise, evaluated at the cell centers read as galactic `(ℓ, b)`.

lonlat_tf = x -> GO.transform(GO.UnitSpherical.GeographicFromUnitSphere(), x)
centers = [lonlat_tf(DGG.cell_centroid(grid, c)) for c in cells]

separation(p, q) = acosd(clamp(sind(p[2]) * sind(q[2]) +
                               cosd(p[2]) * cosd(q[2]) * cosd(p[1] - q[1]), -1, 1))

sources = [(-80.0, 40.0, 4.0), (45.0, -30.0, 3.0), (150.0, 60.0, 2.5), (-160.0, -55.0, 2.0)]
Random.seed!(1234)
sky = map(centers) do c
    diffuse = exp(-(c[2] / 10)^2)
    blobs = sum(amp * exp(-(separation(c, (lon, lat)) / 2)^2) for (lon, lat, amp) in sources)
    diffuse + blobs + 0.05 * randn()
end

# ## The Healpix.jl correspondence
#
# The values vector *is* a nested-order `HealpixMap` — no reshuffle, no copy.
# Every cell center agrees with Healpix.jl's center for the same position:

m = Healpix.HealpixMap{Float64, Healpix.NestedOrder}(sky)
all(i -> collect(Healpix.pix2vecNest(m.resolution, i)) ≈ collect(DGG.cell_centroid(grid, cells[i])),
    1:length(cells))

# ## The all-sky map
#
# One `poly!` over the cell vector paints the whole sphere: Makie reads
# `cells` as one polygon per position, so `color = sky` lines up. Mollweide is
# the projection astronomers reach for; `+over` keeps the cells straddling
# ±180° from smearing across the map.

fig = Figure(size = (800, 420));
ax = GeoAxis(fig[1, 1]; dest = "+proj=moll +over",
    title = "Synthetic all-sky map, galactic coordinates")
plt = poly!(ax, cells; color = sky, colormap = :inferno, strokewidth = 0)
Colorbar(fig[1, 2], plt; label = "brightness")
fig

# ## Cone search
#
# The astronomer's spatial query: everything within 5° of a source. A
# `SphericalCap` is a first-class `query` target, handled exactly, and
# `cellposition` turns the returned ids into positions in `sky`. (`Within` in
# place of `Intersects` would keep only the cells wholly inside the cone.)

to_sphere = GO.UnitSpherical.UnitSphereFromGeographic()
lon0, lat0, _ = sources[1]
cone = GO.UnitSpherical.SphericalCap(to_sphere((lon0, lat0)), deg2rad(5))

idx = DGG.cellposition.(Ref(grid), DGG.query(grid, DGG.Intersects(cone)))
#
(; n = length(idx), cone_mean = mean(sky[idx]), sky_mean = mean(sky))

# The same query, drawn. `cells[idx]` is the returned cells as a cell vector,
# `sky[idx]` is their data, and both share the order of `idx`, so they pair up
# in `poly!` just like `cells` and `sky` did. An orthographic view centred on
# the source shows the neighbourhood faded and the selection at full strength:
# the query returned a disc of cells around the source (the cross).

near = findall(c -> separation(c, (lon0, lat0)) < 15, centers)

fig3 = Figure(size = (520, 540))
ax3 = GeoAxis(fig3[1, 1]; dest = "+proj=ortho +lon_0=$lon0 +lat_0=$lat0",
    limits = ((lon0 - 20, lon0 + 20), (lat0 - 16, lat0 + 16)),
    title = "Cells returned by the 5° cone search")
poly!(ax3, cells[near]; color = sky[near], colormap = :inferno,
    colorrange = extrema(sky), alpha = 0.25, strokewidth = 0)
plt3 = poly!(ax3, cells[idx]; color = sky[idx], colormap = :inferno,
    colorrange = extrema(sky), strokewidth = 0.4, strokecolor = :white)
scatter!(ax3, lon0, lat0; color = :cyan, marker = :xcross, markersize = 12)
Colorbar(fig3[1, 2], plt3; label = "brightness")
fig3

# ## Cutting the galactic plane
#
# To measure the extragalactic sky, astronomy masks the plane first. Cells are
# equal-area — `cell_area` is exactly `4π/npix` — so the masked mean needs no
# latitude weighting.

offplane = findall(c -> abs(c[2]) > 20, centers)
(; all_sky = mean(sky), off_plane = mean(sky[offplane]),
   area_is_exact = DGG.cell_area(grid, cells[1]) ≈ 4pi / DGG.ncells(grid))

# The same cut, drawn: plane cells in gray.

masked = [abs(c[2]) > 20 ? v : NaN for (c, v) in zip(centers, sky)]
fig2 = Figure(size = (800, 420))
ax2 = GeoAxis(fig2[1, 1]; dest = "+proj=moll +over",
    title = "Galactic-plane mask, |b| ≤ 20° in gray")
poly!(ax2, cells; color = masked, colormap = :inferno, nan_color = :gray70, strokewidth = 0)
fig2

# ## Any grid, same query
#
# Nothing above is HEALPix-specific except the nested order. The same cone
# search runs unchanged on any registered system:

length(DGG.query(DGG.levelgrid(DGG.IGeo7System(), 3), DGG.Intersects(cone)))
