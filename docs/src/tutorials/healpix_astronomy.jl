# # The sky in HEALPix
#
# HEALPix is astronomy's grid: every CMB map and all-sky survey ships as a flat
# vector of `12 * nside^2` equal-area pixels. DiscreteGlobalGrids speaks the same
# nested scheme, so an astronomer's map is just a `DimArray` on a `HealpixLookup`.
# This tutorial builds a synthetic all-sky map in galactic coordinates, checks the
# Healpix.jl correspondence, and runs two classic sky operations: a cone search
# around a source and a galactic-plane cut.

using DiscreteGlobalGrids
using DiscreteGlobalGrids.HEALPix.HealpixLookups
import Healpix
import DimensionalData as DD
import GeometryOps as GO
import GeometryOps.SpatialTreeInterface as STI
using Statistics, Random
using CairoMakie, GeoMakie
CairoMakie.activate!()
nothing

# Level 5 means `nside = 2^5 = 32`, i.e. 12288 pixels. We evaluate a fake sky on
# the cell centers — a diffuse band of emission along the galactic plane, four
# bright point sources, and a little noise — reading lon/lat as galactic
# longitude and latitude.

level = 5                                  # nside = 2^level = 32
l = HealpixLookup(DGGSGlobeIds(HEALPixDGGS(), level))
centers = cell_centers(l)                  # (lon, lat) degrees, read as galactic (ℓ, b)

separation(p, q) = acosd(clamp(sind(p[2]) * sind(q[2]) +
    cosd(p[2]) * cosd(q[2]) * cosd(p[1] - q[1]), -1, 1))

sources = [(-80.0, 40.0, 4.0), (45.0, -30.0, 3.0), (150.0, 60.0, 2.5), (-160.0, -55.0, 2.0)]
Random.seed!(1234)
sky = map(centers) do c
    diffuse = exp(-(c[2] / 10)^2)
    blobs = sum(amp * exp(-(separation(c, (lon, lat)) / 2)^2) for (lon, lat, amp) in sources)
    diffuse + blobs + 0.05 * randn()
end
A = DD.DimArray(sky, Cells(l); name = :brightness)
nothing

# A full-globe lookup stores its sorted 0-based nested ids in exactly Healpix.jl's
# (1-based) nested pixel order, so the same vector *is* a `HealpixMap` — no
# reshuffle, no copy of anything but the values.

m = Healpix.HealpixMap{Float64, Healpix.NestedOrder}(parent(A))
theta, phi = Healpix.pix2angNest(m.resolution, 1)   # Healpix pixel 1 == lookup position 1
println("same pixel vector:  ", m.pixels == parent(A))
println("same pixel center:  ", (rad2deg(phi), 90 - rad2deg(theta)) == centers[1])
nothing

# ## The all-sky map
#
# One `poly!` over the cell polygons paints the whole sphere; Mollweide is the
# projection astronomers reach for.

fig = Figure(size = (800, 420))
ax = GeoMakie.GeoAxis(fig[1, 1]; dest = "+proj=moll +over",
    title = "Synthetic all-sky map, galactic coordinates")
plt = poly!(ax, cell_polygons(l); color = collect(parent(A)), colormap = :inferno)
Colorbar(fig[1, 2], plt; label = "brightness")
fig

# ## Cone search
#
# The astronomer's spatial query: everything within 5° of a target. A spherical
# cap around the first source, queried against the lookup's spatial tree, returns
# lookup positions — which index straight into `A`.

lon0, lat0, _ = sources[1]
center = GO.UnitSpherical.UnitSphereFromGeographic()((lon0, lat0))
cap = GO.UnitSpherical.SphericalCap(center, deg2rad(5))
tree = treeify(DGGSPartialGrid(l))
idx = STI.query(tree, intersects_cap(cap))
cone = A[Cells(idx)]
println(length(idx), " cells within 5° of (", lon0, ", ", lat0, ")")
println("cone mean = ", round(mean(cone); digits = 3),
    ",  cone max = ", round(maximum(cone); digits = 3),
    ",  sky mean = ", round(mean(A); digits = 3))
nothing

# The cap query is a superset — every cell whose bounding cap intersects the cone
# — so when exactness matters, refine by center distance.

to_sphere = GO.UnitSpherical.UnitSphereFromGeographic()
inside = filter(i -> GO.UnitSpherical.spherical_distance(center, to_sphere(centers[i])) <= deg2rad(5), idx)
println(length(inside), " of ", length(idx), " cell centers lie strictly within the cone")
nothing

# ## Cutting the galactic plane
#
# To measure the extragalactic sky, astronomy masks the plane first. Cells are
# equal-area, so the masked mean needs no latitude weighting.

offplane = findall(c -> abs(c[2]) > 20, centers)
println("all-sky mean:   ", round(mean(A); digits = 3))
println("|b| > 20° mean: ", round(mean(A[Cells(offplane)]); digits = 3))
nothing

# The same mask, drawn: plane cells in gray, everything else as before.

masked = [abs(c[2]) > 20 ? v : NaN for (c, v) in zip(centers, parent(A))]
fig2 = Figure(size = (800, 420))
ax2 = GeoMakie.GeoAxis(fig2[1, 1]; dest = "+proj=moll +over",
    title = "Galactic-plane mask, |b| ≤ 20° in gray")
poly!(ax2, cell_polygons(l); color = masked, colormap = :inferno, nan_color = :gray70)
fig2
