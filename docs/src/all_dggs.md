# DGGS gallery

This page compares all six DGGS implementations that currently provide cell
geometry. The remaining eight entries returned by `all_systems()` are
registry-only: their metadata is available, but their cell-boundary mathematics
has not yet been ported and deliberately raises `NotPortedError`.

The levels below are chosen independently so that the cells remain visible at
the size of the figure. Every grid follows the same path from a 1-based dense
ordinal to a cell id and then to a unit-sphere polygon. The resulting figure is
emitted as WGLMakie/Bonito HTML rather than as a static raster image.

```@example all-dggs
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
using GeoMakie
using WGLMakie
using Bonito
import Makie
import GeoInterface as GI, GeometryOps as GO

const gallery = (
    (; system = A5DGGS(),       level = 3),
    (; system = H3DGGS(),       level = 1),
    (; system = IGEO7DGGS(),    level = 2),
    (; system = HEALPixDGGS(),  level = 3),
    (; system = S2DGGS(),       level = 3),
    (; system = ISEA4RDGGS(),   level = 3),
)

function cell_polygons(system, level)
    ordinals = 1:DGG.num_cells(system, level)
    cell_ids = ordinal_to_cell.((system,), (level,), ordinals)
    polygons = cell_polygon_unitsphere.((system,), (level,), cell_ids)
end

colors = Makie.wong_colors()
figure = Makie.Figure(size = (1200, 700), fontsize = 18, figure_padding = 4)
Makie.rowgap!(figure.layout, 4)
Makie.colgap!(figure.layout, 4)

for (index, entry) in pairs(gallery)
    row, column = fldmod1(index, 3)
    row = row * 2 - 1 # for labels
    name = system_name(entry.system)
    label = Makie.Label(figure[row, column], text = "$(name) (level $(entry.level))", fontsize = 18, tellheight = true, tellwidth = false)
    axis = GeoMakie.GlobeAxis(
        figure[row+1, column];
        source = "+proj=longlat +R=1",
        dest = GeoMakie.Geodesy.Ellipsoid(; a = "1", b = "1"),
        camera_altitude = 2.0,
        # dest = "+proj=wintri",
        # title = "$(name) — level $(entry.level)",
    )
    GeoMakie.meshimage!(axis, -180..180, -90..90, fill(colorant"white", 1, 1); zlevel = 0.0)
    Makie.lines!(axis, GeoMakie.coastlines(); color = (:gray35, 0.55), linewidth = 0.8, zlevel = 0.005)
    Makie.poly!(
        axis,
        cell_polygons(entry.system, entry.level) |> x -> GO.transform(GO.GeographicFromUnitSphere(), x),
        color = :transparent,
        strokewidth = 0.9,
        strokecolor = :black,
        zlevel = 0.01,
    )
end

figure
```
