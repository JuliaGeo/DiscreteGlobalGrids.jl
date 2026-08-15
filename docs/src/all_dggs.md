# DGGS gallery

The six systems `systems()` returns, drawn the same way: `levelgrid(sys, level)`
goes straight into `Makie.poly!`. The levels differ per panel only because
apertures differ — each is chosen so the cells stay visible at figure size.

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
using WGLMakie   # tutorials use GLMakie; this page activates the web backend
using Bonito
import Makie

WGLMakie.activate!()

figure = Makie.Figure(size = (1200, 850), figure_padding = 4)
for (k, (sys, level)) in enumerate([
        DGG.IGeo7System() => 2, DGG.H3System() => 1, DGG.HEALPixSystem() => 3,
        DGG.A5System() => 3, DGG.S2System() => 3, DGG.ISEA4RSystem() => 3,
    ])
    row, col = fldmod1(k, 3)
    axis = GeoMakie.GlobeAxis(
        figure[row, col];
        source = "+proj=longlat +R=1",
        dest = GeoMakie.Geodesy.Ellipsoid(; a = "1", b = "1"),
        camera_altitude = 2.0,
    )
    GeoMakie.meshimage!(axis, -180..180, -90..90, fill(colorant"white", 1, 1); zlevel = 0.0)
    Makie.lines!(axis, GeoMakie.coastlines(); color = (:gray35, 0.55), linewidth = 0.8, zlevel = 0.005)
    Makie.poly!(axis, DGG.levelgrid(sys, level);
        color = :transparent, strokewidth = 0.9, strokecolor = :black, zlevel = 0.01)
    Makie.Label(figure[row, col, Makie.Top()], "$(nameof(typeof(sys))), level $level"; font = :bold)
end
figure
```

  - **IGeo7** — hexagons with twelve pentagons, aperture 7, equal-area by
    construction.
  - **H3** — the same hexagon family on libh3's gnomonic faces, so not
    equal-area.
  - **HEALPix** — curvilinear diamonds, exactly `4π/(12·4^l)` steradians each.
  - **A5** — Cairo-style pentagons, equal-area; the one system without
    contiguous `descendant_range`s.
  - **S2** — geodesic quadrilaterals on the cube, about a 2× area spread within
    a level.
  - **ISEA4R** — rhombi on ten icosahedral diamonds, exactly `4π/(10·4^l)`
    steradians each.

`AuthalicSystem` wraps any of the six to read geometry at geodetic latitude.
Ids, positions, hierarchy and ordering are untouched, so it draws the same
picture and is not in the sweep.

## What differs

At the *same* level the six disagree about almost everything — apertures 7 and 4
put their cell counts two orders of magnitude apart. The docstring of
`systems()` is the full comparison table.

```@example all-dggs
for sys in DGG.systems()
    grid = DGG.levelgrid(sys, 3)
    println(rpad(nameof(typeof(sys)), 14),
            " levels ", DGG.levels(sys),
            ", ", lpad(DGG.ncells(grid), 6), " cells at level 3",
            ", sorted subtrees: ", DGG.has_sorted_subtrees(sys),
            ", id: ", nameof(DGG.cellindextype(sys)))
end
```
