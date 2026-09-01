# # Choosing a grid

import DiscreteGlobalGrids as DGG                                              # hide
import Geodesy                                                                 # hide
using GLMakie, GeoMakie                                                        # hide
using DiscreteGlobalGridsVisualization: dggpoly, dggpoly!                      # hide
GLMakie.activate!(inline = true)                                               # hide
fig = Figure(size = (1500, 360), figure_padding = 4)                           # hide
for (k, sys) in enumerate((DGG.IGeo7System(), DGG.A5System(), DGG.H3System(),  # hide
                           DGG.HEALPixSystem(), DGG.ISEA4RSystem()))           # hide
    level = DGG.levelfor(sys, 800_000)                                         # hide
    ax = GlobeAxis(fig[1, k]; source = "+proj=longlat +R=1",                   # hide
        dest = Geodesy.Ellipsoid(; a = "1", b = "1"),                          # hide
        camera_longlat = (10, 25), camera_altitude = 1.9)                      # hide
    meshimage!(ax, -180 .. 180, -90 .. 90, fill("#f0faea", 1, 1);              # hide
        zlevel = -0.05, npoints = 300)                                         # hide
    dggpoly!(ax, DGG.levelgrid(sys, level);                                    # hide
        color = "#dcf5d7", strokecolor = "#2c7a1e", strokewidth = 0.9)         # hide
    lines!(ax, GeoMakie.coastlines(); color = ("#212529", 0.7),                # hide
        linewidth = 0.8, zlevel = 0.002)                                       # hide
    Label(fig[1, k, Top()], chopsuffix(String(nameof(typeof(sys))), "System"); # hide
        font = :bold, padding = (0, 0, 6, 0))                                  # hide
end                                                                            # hide
fig                                                                            # hide

# A discrete global grid system tiles the sphere with cells and refines them
# level by level. The five above cover most choices: IGeo7 and H3 lay hexagons
# on an icosahedron, A5 pentagons on a dodecahedron, HEALPix and ISEA4R
# quadrilaterals. Each globe is drawn at the level whose cells come nearest
# 800 km across. The [DGGS gallery](../all_dggs.md) draws every system the
# package ships.
#
# Three properties decide between them: the shape of a cell, how close the
# cells come to equal area, and how fine they go.

import DiscreteGlobalGrids as DGG
import GeometryOps as GO
using Statistics

# ## Compare cell shape, equal-area and fineness
#
# Neighbour tallies and area ratios are over every cell of the level nearest
# 500 km.
#
# | | IGeo7 | A5 | H3 | HEALPix | ISEA4R |
# |---|---|---|---|---|---|
# | **Cell shape** | hexagon, 12 pentagons | pentagon | hexagon, 12 pentagons | curvilinear quadrilateral | rhombus |
# | neighbours, edge / vertex | 6 / 6 (5 / 5 at a pentagon) | 5 / 6–8 | 6 / 6 (5 / 5 at a pentagon) | 4 / 7–8 | 4 / 7–9 |
# | suits | stencils and kernels: one class of neighbour | equal area with a single shape | joining against existing H3 ids | raster habits: 4 across an edge, 4 at a corner | raster habits |
# | **Equal-area**, max/min (middle 90 %) | 1.39 (1.02) | 1.01 (1.01) | 2.22 (1.58) | 1.00 (1.00) | 1.00 (1.00) |
# | suits | unweighted means | unweighted means | means weighted by `cell_area` | exact areal statistics | exact areal statistics |
# | **Fineness**: aperture | 7 | 4 | 7 | 4 | 4 |
# | levels | 0–19 | 0–29 | 0–15 | 0–29 | 0–29 |
# | coarsest → finest cell | 6520 km → 7 cm | 6524 km → 1 cm | 2026 km → 96 cm | 6520 km → 1 cm | 7142 km → 1 cm |
# | suits | steps of ÷2.65 | steps of ÷2 | steps of ÷2.65, 1 m floor | steps of ÷2 | steps of ÷2 |
#
# Compare systems by cell size. A level number counts refinements from a root
# that each system picks for itself, so the same number is a different size in
# every column:
#
# | level | IGeo7 | A5 | H3 | HEALPix | ISEA4R |
# |---|---|---|---|---|---|
# | 0 | 6520 km | 6524 km | 2026 km | 6520 km | 7142 km |
# | 1 | 2721 km | 2917 km | 776 km | 3260 km | 3571 km |
# | 2 | 1021 km | 1459 km | 298 km | 1630 km | 1786 km |
# | 3 | 386 km | 729 km | 113 km | 815 km | 893 km |
# | 4 | 146 km | 365 km | 43 km | 408 km | 446 km |
# | 5 | 55 km | 182 km | 16 km | 204 km | 223 km |
# | 6 | 21 km | 91 km | 6 km | 102 km | 112 km |
# | 7 | 8 km | 46 km | 2 km | 51 km | 56 km |
# | 8 | 3 km | 23 km | 0.9 km | 25 km | 28 km |
# | 9 | 1 km | 11 km | 0.3 km | 13 km | 14 km |
# | 10 | 0.4 km | 6 km | 0.1 km | 6 km | 7 km |
#
# ## Count a cell's neighbours
#
# A hexagon's six neighbours are alike: each shares a full edge, each centroid
# the same distance away. A quadrilateral has two classes — four across an
# edge, four more at a corner √2 further out — the same distinction every
# raster algorithm already makes. `connectivity` picks the class:

igeo7 = DGG.levelgrid(DGG.IGeo7System(), 4)
zurich = DGG.cellat(igeo7, 8.5, 47.4)
(DGG.neighborcount(igeo7, zurich; connectivity = DGG.Edge()),
 DGG.neighborcount(igeo7, zurich; connectivity = DGG.Vertex()))

#

healpix = DGG.levelgrid(DGG.HEALPixSystem(), 4)
diamond = DGG.cellat(healpix, 8.5, 47.4)
(DGG.neighborcount(healpix, diamond; connectivity = DGG.Edge()),
 DGG.neighborcount(healpix, diamond; connectivity = DGG.Vertex()))

# Twelve cells at every level of IGeo7 and H3 are pentagons, one at each vertex
# of the icosahedron. A kernel that assumes six neighbours has to handle those
# twelve; `neighborcount` returns 5 there.

count(c -> DGG.neighborcount(igeo7, c; connectivity = DGG.Edge()) == 5,
      DGG.CellVector(igeo7))

# ## Measure how equal the cell areas are
#
# An unweighted mean over cells is the areal mean over the ground they cover
# when the cells hold equal area. Largest cell area over smallest, across a
# whole level and across the middle ninety per cent of it:

map([DGG.IGeo7System(), DGG.A5System(), DGG.H3System(),
     DGG.HEALPixSystem(), DGG.ISEA4RSystem()]) do sys
    grid = DGG.levelgrid(sys, DGG.levelfor(sys, 500_000))
    areas = [DGG.cell_area(grid, c) for c in DGG.CellVector(grid)]
    lo, hi = quantile(areas, (0.05, 0.95))
    (; system = nameof(typeof(sys)),
       whole_level = round(maximum(areas) / minimum(areas); digits = 2),
       middle_90 = round(hi / lo; digits = 2))
end

# HEALPix and ISEA4R are exact by construction and A5 is within a per cent.
# IGeo7 is equal-area in its own projection and gives a little back to the
# sphere: its twelve pentagons hold about three quarters of a hexagon and
# account for the whole of its first number, which is why its middle ninety per
# cent is so much tighter. H3 puts its cells on gnomonic faces of the
# icosahedron, which buys a fast closed-form index at the price of a cell near a
# face corner holding nearly twice one at its centre — average over H3 cells
# with `cell_area` as the weight.
#
# ## Find the level for a cell size
#
# `cellsize(sys, level)` is the typical cell width in metres, the side of a
# square with the median cell's area:

DGG.cellsize(DGG.IGeo7System(), 5)

# `levelfor` inverts it. It takes metres, or anything carrying a resolution of
# its own — a raster, another grid:

DGG.levelfor(DGG.IGeo7System(), 25_000)

# Cell width shrinks by √7 per level on IGeo7 and H3 and by 2 on A5, HEALPix
# and ISEA4R, so `levelfor` returns the level nearest in ratio.
#
# ## Match the ellipsoid of the source
#
# A grid lays its cells on a sphere of equal area, and the package reads their
# geometry at *authalic* latitude. A GPS fix, a shapefile and a GeoTIFF carry
# *geodetic* latitude, on the WGS84 ellipsoid. The two differ by up to 0.128°,
# 14 km along a meridian at ±45°.
#
# `AuthalicSystem` reads a system's geometry at geodetic latitude, wrapping a
# system that computes at authalic latitude — four of the five:

DGG.AuthalicSystem(DGG.IGeo7System())

#

[DGG.AuthalicSystem(sys) for sys in (DGG.H3System(), DGG.HEALPixSystem(),
                                     DGG.ISEA4RSystem())]

# A5 converts to geodetic latitude inside its own projection, so its geometry
# is geodetic already — wrapping it would convert twice, and the constructor
# refuses:

try
    DGG.AuthalicSystem(DGG.A5System())
catch err
    err
end

# The wrapper moves the geometry and forwards the rest:
#
# | | verbs |
# |---|---|
# | read at geodetic latitude | `cell_boundary`, `cell_centroid`, `cell_area`, `cell_cap` |
# | takes its query point at geodetic latitude | `cellat` |
# | forwarded unchanged | cell ids, `level`, `neighbors`, `ring`, `parent`, `children`, ordering |
#
# Over every cell of IGeo7 level 4, the centroid moves along its meridian by:

sys = DGG.IGeo7System()
sphere = DGG.levelgrid(sys, 4)
ellipsoid = DGG.levelgrid(DGG.AuthalicSystem(sys), 4)
cells = DGG.CellVector(sphere)
lonlat = GO.GeographicFromUnitSphere()
shift = last.(lonlat.(DGG.cell_centroid.(ellipsoid, cells))) .-
        last.(lonlat.(DGG.cell_centroid.(sphere, cells)))
maximum(abs, shift)

#

maximum(abs, shift) * 111.2   # a degree of latitude is about 111.2 km

# A 14 km shift crosses a cell boundary once the cells are smaller than that.
# Zürich at 8 km cells lands in a different cell on each grid:

level = DGG.levelfor(sys, 10_000)
DGG.cellat(DGG.levelgrid(sys, level), 8.5, 47.4),
    DGG.cellat(DGG.levelgrid(DGG.AuthalicSystem(sys), level), 8.5, 47.4)

# So match the ellipsoid of the source. `regrid` reads source and destination
# coordinates in one frame and converts neither:
#
#   - A lon/lat raster, a shapefile or a GPS track carries geodetic latitude,
#     so regrid it onto `AuthalicSystem(sys)`.
#   - Data already on a DGGS carries authalic latitude, so it moves between
#     plain systems.
#   - Mixing the two misregisters by the 0.13° above, a cell and a half at
#     10 km.
#
# `AuthalicSystem` is also what reports a cell to somebody else's data:
# `cell_boundary` on a wrapped grid returns a ring that overlays a shapefile
# directly.
