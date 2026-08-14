# DGGS gallery

This page draws the six grid systems `systems()` returns: IGeo7, H3, HEALPix,
A5, S2, and ISEA4R. Every one of them answers the same interface, so one
function draws all six — a system and a level go in, and cell polygons come out.

The levels differ from panel to panel only so that the cells stay visible at the
size of the figure. Each globe is emitted as WGLMakie/Bonito HTML rather than as
a static raster image.

```@raw html
<style>
#VPContent .bonito-fragment canvas {
    max-width: 100%;
    height: auto !important;
}
</style>
```

```@example all-dggs
import DiscreteGlobalGrids as DGG
using GeoMakie
using WGLMakie
using Bonito
import Makie
import GeometryOps as GO

# Other pages use CairoMakie, so this one re-activates its own backend.
WGLMakie.activate!()

# `levelgrid` is the complete level; a position `i` names a cell through
# `cellindex`, and `cell_polygon` is that cell on the unit sphere.
function cell_polygons(sys, level)
    grid = DGG.levelgrid(sys, level)
    cells = [DGG.cellindex(grid, i) for i in 1:DGG.ncells(grid)]
    return GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygon.(Ref(grid), cells))
end

function globe_figure(sys, level)
    figure = Makie.Figure(size = (650, 580), figure_padding = 4)
    axis = GeoMakie.GlobeAxis(
        figure[1, 1];
        source = "+proj=longlat +R=1",
        dest = GeoMakie.Geodesy.Ellipsoid(; a = "1", b = "1"),
        camera_altitude = 2.0,
    )
    GeoMakie.meshimage!(axis, -180..180, -90..90, fill(colorant"white", 1, 1); zlevel = 0.0)
    Makie.lines!(axis, GeoMakie.coastlines(); color = (:gray35, 0.55), linewidth = 0.8, zlevel = 0.005)
    Makie.poly!(
        axis,
        cell_polygons(sys, level),
        color = :transparent,
        strokewidth = 0.9,
        strokecolor = :black,
        zlevel = 0.01,
    )
    figure
end

nothing
```

## IGeo7

Level 2. Hexagons with twelve pentagons, aperture 7, equal-area by construction.

```@example all-dggs
globe_figure(DGG.IGeo7System(), 2)
```

## H3

Level 1. The same hexagon-with-pentagon family, on libh3's gnomonic faces, so
not equal-area.

```@example all-dggs
globe_figure(DGG.H3System(), 1)
```

## HEALPix

Level 3. Curvilinear diamonds, exactly `4π/(12·4^l)` steradians each.

```@example all-dggs
globe_figure(DGG.HEALPixSystem(), 3)
```

## A5

Level 3. Cairo-style pentagons, equal-area; the one system whose canonical order
does not give contiguous `descendant_range`s.

```@example all-dggs
globe_figure(DGG.A5System(), 3)
```

## S2

Level 3. Geodesic quadrilaterals on the cube, with about a 2× area spread within
a level.

```@example all-dggs
globe_figure(DGG.S2System(), 3)
```

## ISEA4R

Level 3. Rhombi on ten icosahedral diamonds, exactly `4π/(10·4^l)` steradians
each.

```@example all-dggs
globe_figure(DGG.ISEA4RSystem(), 3)
```

## Ellipsoidal geometry

`AuthalicSystem` re-reads any of the six at geodetic latitude, leaving ids,
positions, hierarchy and ordering untouched. It is a system like any other, so
it goes through the very same call.

```@example all-dggs
globe_figure(DGG.AuthalicSystem(DGG.ISEA4RSystem()), 3)
```

## What differs between them

`systems()`'s own docstring is the comparison table — cell counts, cell shape,
equal-areaness, neighbour degree, and the traits that differ.

```@example all-dggs
for sys in DGG.systems()
    grid = DGG.levelgrid(sys, 3)
    println(rpad(nameof(typeof(sys)), 16),
            " levels ", DGG.levels(sys),
            ", ", lpad(DGG.ncells(grid), 6), " cells at level 3",
            ", sorted subtrees: ", DGG.has_sorted_subtrees(sys),
            ", id: ", nameof(DGG.cellindextype(sys)))
end
```
