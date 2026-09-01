# # The sky in HEALPix
#
# HEALPix is astronomy's grid: every CMB map and all-sky survey ships as a flat
# vector of `12 * nside^2` equal-area pixels in nested order. That order is the
# index order of `levelgrid(HEALPixSystem(), level)`, so a HEALPix map is a
# `DimArray` over `Cells`, pixel for pixel.

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import Healpix
import GeometryOps as GO
using Statistics, Random, LinearAlgebra
using DiscreteGlobalGridsVisualization: dggpoly, dggpoly!
using GLMakie, GeoMakie
GLMakie.activate!(inline = true)

# ## Checking the pixel order against Healpix.jl
#
# Level `l` is `nside = 2^l`. Healpix.jl and the grid agree on the pixel count:

grid = DGG.levelgrid(DGG.HEALPixSystem(), 6)
lookup = DGG.CellLookup(grid)
resolution = Healpix.Resolution(64)
DGG.ncells(grid), 12 * 64^2

# and on where every pixel sits. The lookup indexes to cell ids, so `lookup[i]`
# is the cell Healpix.jl calls pixel `i`:

all(i -> collect(Healpix.pix2vecNest(resolution, i)) ≈
         collect(DGG.cell_centroid(grid, lookup[i])), 1:DGG.ncells(grid))

# ## Building a synthetic sky as a DimArray over Cells
#
# The field is a diffuse band along the galactic plane, four point sources and
# a little noise, read at the cell centres as galactic longitude and latitude
# `(ℓ, b)`. `synthetic_sky` takes a grid and returns a `DimArray` over a
# `Cells` dimension, which is the object every step below works on. A real map
# drops in at the same place: wrap a nested-order FITS vector at this `nside`
# as `DD.DimArray(pixels, DGG.Cells(lookup))`.

lonlat = GO.GeographicFromUnitSphere()
to_sphere = GO.UnitSphereFromGeographic()
sources = [(-80.0, 40.0, 4.0), (45.0, -30.0, 3.0), (150.0, 60.0, 2.5), (-160.0, -55.0, 2.0)]

function synthetic_sky(g)
    lk = DGG.CellLookup(g)
    Random.seed!(1234)
    brightness = map(DGG.cell_centroid.(g, lk)) do p
        ℓ, b = lonlat(p)
        diffuse = exp(-((b - 4 * sind(ℓ - 30)) / 8)^2)
        blobs = sum(amp * exp(-(rad2deg(GO.UnitSpherical.spherical_distance(
                        p, to_sphere((slon, slat)))) / 5)^2)
                    for (slon, slat, amp) in sources)
        diffuse + blobs + 0.05 * randn()
    end
    return DD.DimArray(brightness, DGG.Cells(lk); name = :brightness)
end

sky = synthetic_sky(grid)

# The values span the diffuse floor up to the brightest source:

extrema(sky)

# The same values are a Healpix.jl `HealpixMap` in nested order, and a
# `HealpixMap`'s pixels are a `DimArray` over `Cells`, sharing the one storage
# vector:

m = Healpix.HealpixMap{Float64, Healpix.NestedOrder}(parent(sky))
back = DD.DimArray(m.pixels, DGG.Cells(lookup); name = :brightness)
parent(back) === m.pixels

# ## Masking the galactic plane
#
# Astronomy masks the plane before measuring the extragalactic sky. Equal-area
# cells make the masked mean a plain `mean`: every cell weighs `4π/npix`.

DGG.cell_area(grid, first(lookup)) ≈ 4pi / DGG.ncells(grid)

# The mask is a `Bool` vector over the cells, keeping every centre more than
# 20° from the plane. Indexing by it returns a `DimArray` whose lookup carries
# the surviving cells:

b = last.(lonlat.(DGG.cell_centroid.(grid, lookup)))
offplane = abs.(b) .> 20
sky[offplane]

# The band holds most of the flux, so the off-plane mean is what the four point
# sources contribute:

mean(sky), mean(sky[offplane])

# `dggpoly!` takes the `DimArray` directly: cells from its `Cells` lookup,
# colours from its values. The masked cells are drawn in gray.

crange = extrema(sky)
fig = Figure(size = (960, 320))
ax1 = GeoAxis(fig[1, 1]; dest = "+proj=moll", xticks = -180:60:180, yticks = -60:30:60,
    xticklabelsvisible = false, title = "Synthetic all-sky map")
plt = dggpoly!(ax1, sky; color = sky, colormap = :inferno, strokewidth = 0)
ax2 = GeoAxis(fig[1, 2]; dest = "+proj=moll", xticks = -180:60:180, yticks = -60:30:60,
    xticklabelsvisible = false, title = "Galactic plane masked, |b| ≤ 20°")
dggpoly!(ax2, sky[offplane]; color = sky[offplane], colormap = :inferno,
    colorrange = crange, strokewidth = 0)
dggpoly!(ax2, sky[.!offplane]; color = "#d9d9d9", strokewidth = 0)
Colorbar(fig[1, 3], plt; label = "brightness")
colgap!(fig.layout, 8)
fig

# ## Cone search with a spherical cap
#
# A cone search selects every cell within 5° of a source. The cone is a
# `SphericalCap`, and a query predicate wrapped in `Cells` is a selector, so
# `sky[Cells(Intersects(cone))]` is the cone's piece of the map:

lon0, lat0, _ = sources[1]
cone = GO.UnitSpherical.SphericalCap(to_sphere((lon0, lat0)), deg2rad(5))
incone = sky[DGG.Cells(DGG.Intersects(cone))]

# The cone sits on the brightest source. Its mean against the all-sky mean:

mean(incone), mean(sky)

# `Within(cone)` keeps the cells lying wholly inside the cone. The difference
# between the two counts is the ring of cells straddling the rim:

inside = sky[DGG.Cells(DGG.Within(cone))]
length(incone), length(inside)

# ## Running the same cone search on H3
#
# The same selector runs on any grid. On the H3 level whose cells are nearest
# in size to this HEALPix level:

h3 = DGG.levelgrid(DGG.H3System(), DGG.levelfor(DGG.H3System(), DGG.cellsize(grid)))
sky_h3 = synthetic_sky(h3)
incone_h3 = sky_h3[DGG.Cells(DGG.Intersects(cone))]
inside_h3 = sky_h3[DGG.Cells(DGG.Within(cone))]
length(incone_h3), length(inside_h3)

# Both grids in one frame, azimuthal equidistant centred on the source, so the
# 5° cone is a circle on the page. Three tones:
#
# - dark purple: cells wholly inside the cone (`Within`);
# - light purple: cells crossing the rim, which `Intersects` adds;
# - gray: the surrounding sky.

rim_u = normalize(cross(cone.point, GO.UnitSpherical.UnitSphericalPoint(0.0, 0.0, 1.0))) #hide
rim_v = cross(cone.point, rim_u) #hide
rim = [lonlat(cos(cone.radius) .* cone.point .+ #hide
              sin(cone.radius) .* (cos(t) .* rim_u .+ sin(t) .* rim_v)) #hide
       for t in range(0, 2pi; length = 181)] #hide
context = GO.UnitSpherical.SphericalCap(cone.point, deg2rad(11))

fig2 = Figure(size = (860, 430))
for (k, (g, A, met, within)) in enumerate([(grid, sky, incone, inside),
                                           (h3, sky_h3, incone_h3, inside_h3)])
    ax = GeoAxis(fig2[1, k]; dest = "+proj=aeqd +lon_0=$lon0 +lat_0=$lat0",
        limits = ((lon0 - 16, lon0 + 16), (lat0 - 12, lat0 + 12)),
        xgridvisible = false, ygridvisible = false,
        xticklabelsvisible = false, yticklabelsvisible = false,
        title = "$(nameof(typeof(DGG.system(g)))) level $(DGG.level(g)): " *
                "$(length(met)) meet the cone, $(length(within)) inside")
    dggpoly!(ax, A[DGG.Cells(DGG.Intersects(context))]; color = "#e9ecef",
        strokecolor = ("#212529", 0.35), strokewidth = 0.4)
    dggpoly!(ax, met; color = "#cbb3dd", strokecolor = "#f8f9fa", strokewidth = 0.6)
    dggpoly!(ax, within; color = "#9558b2", strokecolor = "#f8f9fa", strokewidth = 0.6)
    lines!(ax, rim; color = :black, linewidth = 2)
    scatter!(ax, [lon0], [lat0]; color = :white, strokecolor = :black,
        strokewidth = 1.2, marker = :star5, markersize = 16)
end
fig2

# ## Doing it without the DimArray
#
# A FITS reader hands back a plain `Vector`, and two calls map a selection onto
# it: `query` returns the cell ids a predicate selects, and `globalindex` turns
# an id into its position in the level grid — the pixel number.

ids = DGG.query(grid, DGG.Intersects(cone))
idx = DGG.globalindex.(grid, ids)
mean(m.pixels[idx]) ≈ mean(incone)
