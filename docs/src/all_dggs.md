# DGGS gallery

This page compares all six DGGS implementations that currently provide cell
geometry. The remaining eight entries returned by `all_systems()` are
registry-only: their metadata is available, but their cell-boundary mathematics
has not yet been ported and deliberately raises `NotPortedError`.

The levels below are chosen independently so that the cells remain visible at
the size of each figure. Every grid follows the same path from a 1-based dense
ordinal to a cell id and then to a unit-sphere polygon. Each globe is emitted as
WGLMakie/Bonito HTML rather than as a static raster image.

```@raw html
<style>
#VPContent .bonito-fragment canvas {
    max-width: 100%;
    height: auto !important;
}
</style>
```

```@example all-dggs
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
using GeoMakie
using WGLMakie
using Bonito
import Makie
import GeoInterface as GI, GeometryOps as GO

function cell_polygons(system, level)
    ordinals = 1:DGG.num_cells(system, level)
    cell_ids = ordinal_to_cell.((system,), (level,), ordinals)
    polygons = cell_polygon_unitsphere.((system,), (level,), cell_ids)
end

function globe_figure(system, level)
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
        cell_polygons(system, level) |> x -> GO.transform(GO.GeographicFromUnitSphere(), x),
        color = :transparent,
        strokewidth = 0.9,
        strokecolor = :black,
        zlevel = 0.01,
    )
    figure
end

nothing
```

## A5

Level 3.

```@example all-dggs
globe_figure(A5DGGS(), 3)
```

## H3

Level 1.

```@example all-dggs
globe_figure(H3DGGS(), 1)
```

## IGEO7

Level 2.

```@example all-dggs
globe_figure(IGEO7DGGS(), 2)
```

## HEALPix

Level 3.

```@example all-dggs
globe_figure(HEALPixDGGS(), 3)
```

## S2

Level 3.

```@example all-dggs
globe_figure(S2DGGS(), 3)
```

## ISEA4R

Level 3.

```@example all-dggs
globe_figure(ISEA4RDGGS(), 3)
```
