# # The sky in HEALPix
#
# HEALPix is astronomy's grid: every CMB map and all-sky survey ships as a flat
# vector of `12 * nside^2` equal-area pixels, in nested order. That order is
# exactly `levelgrid(HEALPixSystem(), level)`'s dense position order, so an
# astronomer's map and a DiscreteGlobalGrids data vector are the same vector.
#
# This page builds a synthetic all-sky map in galactic coordinates, checks the
# Healpix.jl correspondence, and runs two classic sky operations: a cone search
# around a source and a galactic-plane cut.

import DiscreteGlobalGrids as DGG
import Healpix
import GeometryOps as GO
using Statistics, Random
using CairoMakie, GeoMakie
CairoMakie.activate!()

# Level 5 means `nside = 2^5 = 32`, i.e. 12288 pixels. The fake sky is a diffuse
# band of emission along the galactic plane, four bright point sources, and a
# little noise, evaluated at the cell centroids and read as galactic `(ℓ, b)`.

level = 5
grid = DGG.levelgrid(DGG.HEALPixSystem(), level)
cells = [DGG.cellindex(grid, i) for i in 1:DGG.ncells(grid)]

lonlat = GO.UnitSpherical.GeographicFromUnitSphere()
centers = [lonlat(DGG.cell_centroid(grid, c)) for c in cells]

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
# Position `i` of the level grid is Healpix.jl's nested pixel `i`, so the values
# vector *is* a `HealpixMap` — no reshuffle, no copy of anything but the values.

m = Healpix.HealpixMap{Float64, Healpix.NestedOrder}(sky)
theta, phi = Healpix.pix2angNest(m.resolution, 1)
(; same_pixels = m.pixels == sky,
   same_center = (rad2deg(phi), 90 - rad2deg(theta)) .≈ (mod(centers[1][1], 360), centers[1][2]))

# ## The all-sky map
#
# One `poly!` over the cell polygons paints the whole sphere; Mollweide is the
# projection astronomers reach for.

polys = GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygon.(Ref(grid), cells))

fig = Figure(size = (800, 420))
ax = GeoAxis(fig[1, 1]; dest = "+proj=moll +over",
    title = "Synthetic all-sky map, galactic coordinates")
plt = poly!(ax, polys; color = sky, colormap = :inferno, strokewidth = 0)
Colorbar(fig[1, 2], plt; label = "brightness")
fig

# ## Cone search
#
# The astronomer's spatial query: everything within 5° of a target. A
# `SphericalCap` is a first-class `query` target — handled exactly, without
# polygonising it — so the cone search is one call, and the returned ids become
# data positions through `cellposition`.

to_sphere = GO.UnitSpherical.UnitSphereFromGeographic()
lon0, lat0, _ = sources[1]
cone = GO.UnitSpherical.SphericalCap(to_sphere((lon0, lat0)), deg2rad(5))

hits = DGG.query(grid, DGG.Intersects(cone))
idx = [DGG.cellposition(grid, c) for c in hits]
(; n = length(idx), cone_mean = mean(sky[idx]), cone_max = maximum(sky[idx]),
   sky_mean = mean(sky))

# `Intersects` keeps every cell that meets the cone, including those only
# clipped by it. `Within` keeps the cells wholly inside, and the two bracket any
# centre-based rule an astronomer might prefer.

inner = DGG.query(grid, DGG.Within(cone))
centred = [i for i in idx
           if GO.UnitSpherical.spherical_distance(cone.point, to_sphere(centers[i])) <= cone.radius]
(; n_within = length(inner), n_centred = length(centred), n_intersecting = length(idx))

# ## Cutting the galactic plane
#
# To measure the extragalactic sky, astronomy masks the plane first. Cells are
# equal-area, so the masked mean needs no latitude weighting. `cell_area` on the
# level grid returns the exact `4π/npix`.

offplane = findall(c -> abs(c[2]) > 20, centers)
(; all_sky = mean(sky), off_plane = mean(sky[offplane]),
   cell_area = DGG.cell_area(grid, cells[1]), exact = 4pi / DGG.ncells(grid))

# The same mask, drawn: plane cells in gray, everything else as before.

masked = [abs(c[2]) > 20 ? v : NaN for (c, v) in zip(centers, sky)]
fig2 = Figure(size = (800, 420))
ax2 = GeoAxis(fig2[1, 1]; dest = "+proj=moll +over",
    title = "Galactic-plane mask, |b| ≤ 20° in gray")
poly!(ax2, polys; color = masked, colormap = :inferno, nan_color = :gray70, strokewidth = 0)
fig2

# ## The same cone on every system
#
# Nothing above is HEALPix-specific except the nested-order claim. A cap query
# is an interface method, so the cone search runs unchanged on every registered
# system, and on an `AuthalicSystem` wrap that re-reads the geometry at geodetic
# latitude.

for sys in (DGG.systems()..., DGG.AuthalicSystem(DGG.HEALPixSystem()))
    base = sys isa DGG.AuthalicSystem ? parent(sys) : sys
    l = base isa Union{DGG.IGeo7System, DGG.H3System} ? 3 : 4
    g = DGG.levelgrid(sys, l)
    name = sys isa DGG.AuthalicSystem ? "Authalic($(nameof(typeof(base))))" :
           string(nameof(typeof(sys)))
    println(rpad(name, 24), "level $l: ",
            lpad(length(DGG.query(g, DGG.Within(cone))), 4), " within, ",
            lpad(length(DGG.query(g, DGG.Intersects(cone))), 4), " intersecting")
end
