# DGGS gallery

The six systems, drawn the same way: `levelgrid(sys, level)` goes straight into
`dggpoly!`. The levels differ per panel only because apertures do.

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
using DiscreteGlobalGridsVisualization

WGLMakie.activate!()

panels = [
    DGG.IGeo7System() => 2, DGG.H3System() => 1, DGG.HEALPixSystem() => 3,
    DGG.A5System() => 3, DGG.S2System() => 3, DGG.ISEA4RSystem() => 3,
]

figure = Figure(size = (1200, 850), figure_padding = 4)
for (k, (sys, level)) in enumerate(panels)
    row, col = fldmod1(k, 3)
    axis = GlobeAxis(
        figure[row, col];
        source = "+proj=longlat +R=1",
        dest = GeoMakie.Geodesy.Ellipsoid(; a = "1", b = "1"),
        camera_altitude = 2.0,
    )
    meshimage!(axis, -180 .. 180, -90 .. 90, fill("#f0faea", 1, 1);
        zlevel = -0.05, npoints = 300)
    dggpoly!(axis, DGG.levelgrid(sys, level);
        color = "#dcf5d7", strokecolor = "#2c7a1e", strokewidth = 0.9)
    lines!(axis, GeoMakie.coastlines(); color = ("#212529", 0.7), linewidth = 0.8,
        zlevel = 0.002)
    Label(figure[row, col, Top()], "$(nameof(typeof(sys))), level $level"; font = :bold)
end
figure
```

## What differs

At the *same* level the six disagree about almost everything — apertures 7 and
4 put their cell counts two orders of magnitude apart:

```@example all-dggs
println(rpad("system", 10), lpad("levels", 8), lpad("cells at level 3", 18),
        "   id type")
for (sys, _) in panels
    println(rpad(chopsuffix(String(nameof(typeof(sys))), "System"), 10),
            lpad(string(DGG.levels(sys)), 8),
            lpad(DGG.ncells(DGG.levelgrid(sys, 3)), 18),
            "   ", nameof(DGG.cellindextype(sys)))
end
```

  - **IGeo7** — hexagons with twelve pentagons, aperture 7, equal-area by
    construction.
  - **H3** — the same hexagon family on gnomonic icosahedral faces, so not
    equal-area.
  - **HEALPix** — curvilinear diamonds, exactly `4π/(12·4^l)` steradians each.
  - **A5** — Cairo-style pentagons, equal-area.
  - **S2** — geodesic quadrilaterals on the cube, about a 2× area spread within
    a level.
  - **ISEA4R** — rhombi on ten icosahedral diamonds, exactly `4π/(10·4^l)`
    steradians each.

`AuthalicSystem` wraps a system to read its geometry at geodetic latitude. Ids,
indices, hierarchy and ordering are untouched, so it draws the same picture and
is not in the sweep. A5 is the exception: it converts to authalic latitude
inside its own projection, so its geometry is geodetic already and the wrapper
refuses it rather than converting twice.

[Choosing a grid](tutorials/choosing_a_grid.md) is the walk through the axes
that decide between them: cell shape, how close to equal-area, how fine, and
which sphere.
