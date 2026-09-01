# Writing a grid system

```@meta
CurrentModule = DiscreteGlobalGrids
```

DiscreteGlobalGrids describes a grid system as methods on two types you write: a
cell id (`<: AbstractCellIndex`) and a system (`<: AbstractHierarchicalGridSystem`).
Every algorithm in the package — point location, neighbourhoods, queries,
regridding, plotting, store I/O — is written once against those methods, so a
system that implements the required set gets all of them.

Thirteen methods are required. Three describe the cell id, an immutable
`isbits` struct subtyping `AbstractCellIndex`:

| method | contract |
|---|---|
| `level(c)` | the level the cell belongs to |
| `rawid(c)` | the integer the id stores; together with `level` it round-trips the id |
| `Base.isless(a, b)` | canonical cell order, which every sorted structure downstream uses |

Ten describe the system:

| method | contract |
|---|---|
| `cellindextype(sys)` | the id type the system hands out |
| `levels(sys)` | a `UnitRange`, coarsest level first |
| `ncells(sys, l)` | how many cells level `l` has |
| [`cellindex`](@ref)`(sys, l, i)` | dense index → id |
| [`globalindex`](@ref)`(sys, c)` | id → dense index, or `nothing` for a cell outside the level |
| `cell_boundary(sys, c)` | the cell's ring, closed and counter-clockwise, as unit-sphere points |
| `cell_centroid(sys, c)` | a point strictly inside the cell |
| `rootcells(sys)` | the coarsest level's cells, ascending |
| `Base.parent(sys, c)` | the cell one level up that contains `c` |
| `children(sys, c)` | the cells one level down that `c` contains, ascending |

Everything beyond those thirteen has a default or a generic implementation, and
falls into three groups:

| group | methods | what writing them buys |
|---|---|---|
| declarations | `has_congruent_refinement`, [`has_sorted_subtrees`](@ref), [`has_direct_location`](@ref), [`maxneighbors`](@ref), `maxring`, `winding`, `cap_inflation`, [`node_extent`](@ref) | the engines take a faster route once a system promises the property holds |
| store I/O | `idrank`, `idselect`, `idvalid`, `idcell`, `register_grid!` | `dggwrite` and `dggread` |
| fast paths | `cellat`, `localindex`, `neighbors`, `ring`, `cell_area`, `cell_extent`, `getcell`, `ancestor`, `descendants`, `subtree` | speed; the generic implementations already answer every one of them |

This page implements the thirteen for a longitude/latitude lattice, runs the
conformance suites over it, puts the package's high-level API through the
result, and then adds two of the optional groups.

## Define the cell id

Level `l` of the lattice is `8·2^l` columns by `4·2^l` rows of equal-angle
boxes, so a cell is named by its level, row and column. Ordering the ids by
that triple orders the cells west to east, then south to north, then coarse to
fine — the canonical order the interface asks for.

```@example extending
import DiscreteGlobalGrids as DGG
import GeometryOps as GO

struct LonLatCell <: DGG.AbstractCellIndex
    level::Int
    row::Int
    col::Int
end

ncolumns(l::Integer) = 8 << l
nrows(l::Integer) = 4 << l

DGG.level(c::LonLatCell) = c.level
DGG.rawid(c::LonLatCell) = (c.row - 1) * ncolumns(c.level) + c.col
Base.isless(a::LonLatCell, b::LonLatCell) =
    isless((a.level, a.row, a.col), (b.level, b.row, b.col))

LonLatCell(4, 49, 68)
```

`LonLatCell` is immutable and holds three `Int`s, so `==` and `hash` come from
Julia's defaults and already agree with `isless`.

## Define the system

The system itself carries no state — the lattice is fixed by the level — so it
is an empty struct. Five methods tie an id to a position in the level:

```@example extending
struct LonLatSystem <: DGG.AbstractHierarchicalGridSystem end

DGG.cellindextype(::LonLatSystem) = LonLatCell
DGG.levels(::LonLatSystem) = 0:12

DGG.ncells(::LonLatSystem, l::Integer) = ncolumns(l) * nrows(l)

function DGG.cellindex(::LonLatSystem, l::Integer, i::Integer)
    row, col = fldmod1(Int(i), ncolumns(l))
    return LonLatCell(Int(l), row, col)
end

DGG.globalindex(::LonLatSystem, c::LonLatCell) =
    (1 <= c.row <= nrows(c.level) && 1 <= c.col <= ncolumns(c.level)) ?
        DGG.rawid(c) : nothing

DGG.cellindex(LonLatSystem(), 4, 5000)
```

### Draw a cell on the unit sphere

`cell_boundary` returns the cell's ring on the unit sphere and `cell_centroid` a
point inside it. A cell's east and west edges are meridians, which are great
circles, so their endpoints describe them; its north and south edges are
parallels, which are small circles, so each is densified into `PARALLEL_STEPS`
segments.

```@example extending
const TO_SPHERE = GO.UnitSpherical.UnitSphereFromGeographic()
const PARALLEL_STEPS = 8

# (west, east, south, north), in degrees.
function cellbox(c::LonLatCell)
    dlon, dlat = 360 / ncolumns(c.level), 180 / nrows(c.level)
    return (-180 + (c.col - 1) * dlon, -180 + c.col * dlon,
            -90 + (c.row - 1) * dlat, -90 + c.row * dlat)
end

function DGG.cell_boundary(::LonLatSystem, c::LonLatCell)
    w, e, s, n = cellbox(c)
    ring = GO.UnitSpherical.UnitSphericalPoint{Float64}[]
    for k in 0:PARALLEL_STEPS                       # south parallel, west to east
        push!(ring, TO_SPHERE((w + (e - w) * k / PARALLEL_STEPS, s)))
    end
    for k in 0:PARALLEL_STEPS                       # north parallel, east to west
        push!(ring, TO_SPHERE((e - (e - w) * k / PARALLEL_STEPS, n)))
    end
    # A cell touching a pole collapses its polar edge to one point.
    return [p for (i, p) in enumerate(ring) if i == 1 || !(p ≈ ring[i - 1])]
end

function DGG.cell_centroid(::LonLatSystem, c::LonLatCell)
    w, e, s, n = cellbox(c)
    return TO_SPHERE(((w + e) / 2, (s + n) / 2))
end

DGG.cell_centroid(LonLatSystem(), LonLatCell(4, 49, 68))
```

### Link parents to children

Queries descend the hierarchy to a region without materialising the levels above
it, so `parent` and `children` are arithmetic on the id: no lookup table, no
geometry. Both ends of the level range throw an `ArgumentError`, and `children`
comes back ascending.

```@example extending
DGG.rootcells(::LonLatSystem) =
    [LonLatCell(0, r, k) for r in 1:nrows(0) for k in 1:ncolumns(0)]

function Base.parent(::LonLatSystem, c::LonLatCell)
    c.level == 0 && throw(ArgumentError("a root cell has no parent: $c"))
    return LonLatCell(c.level - 1, cld(c.row, 2), cld(c.col, 2))
end

function DGG.children(sys::LonLatSystem, c::LonLatCell)
    c.level == last(DGG.levels(sys)) &&
        throw(ArgumentError("a cell at the deepest level has no children: $c"))
    l = c.level + 1
    return [LonLatCell(l, 2c.row - 1, 2c.col - 1), LonLatCell(l, 2c.row - 1, 2c.col),
            LonLatCell(l, 2c.row,     2c.col - 1), LonLatCell(l, 2c.row,     2c.col)]
end

# The four children tile their parent exactly.
DGG.has_congruent_refinement(::LonLatSystem) = true

DGG.children(LonLatSystem(), LonLatCell(0, 3, 8))
```

## Run the conformance suites

[`DiscreteGlobalGridsConformanceTesting`](https://github.com/JuliaGeo/DiscreteGlobalGrids.jl/tree/main/lib/DiscreteGlobalGridsConformanceTesting)
checks two things.

**That the interface works.** `cellindex` and `localindex` are inverses; every
boundary ring is closed and winds counter-clockwise; every centroid lies inside
its own cell; `parent` and `children` are inverses and throw at both ends of the
hierarchy; every descendant's boundary lies inside every ancestor's extent.

**That what the system declares about itself is true.** `levels` ends where
`maxlevel` says it does; `cap_inflation` is at least 1; a system declaring
[`has_sorted_subtrees`](@ref) answers [`descendant_range`](@ref); a declared
[`maxneighbors`](@ref) is a real ceiling; a declared `winding` matches the order
the boundaries actually run in.

Each optional law group is gated on the method it tests. A group whose method
the system has not written is skipped, and the `@info` at the end names every
skip together with the signature that would open it. `test_grid_interface`
checks one grid:

```@example extending
using DiscreteGlobalGridsConformanceTesting
sys = LonLatSystem()
grid = DGG.levelgrid(sys, 4)
test_grid_interface(grid)
nothing # hide
```

`test_hierarchical_system` checks the hierarchy against every level the system
declares:

```@example extending
test_hierarchical_system(sys)
nothing # hide
```

Six methods are unwritten, and the skip list is the plan: `cellat`, `ancestor`,
`descendants`, `neighbors`, `ring` and `descendant_range`. The rest of this page
uses the grid as it stands, then writes the ones this lattice can answer in
closed form.

A third suite, `test_generic_fallbacks(sys)`, hides a system's fast paths from
dispatch and re-runs the adjacency and point-location laws against the generic
implementations underneath them.

## Locate a point, measure a cell, draw the grid

`cellat` is answered by the generic fallback, which descends the hierarchy from
the root cells:

```@example extending
DGG.cellat(grid, 10.5, 46.5)
```

`cell_area` is answered too, by taking the spherical area of the ring
`cell_boundary` returned. An equal-angle cell at the equator covers forty times
the area of one at the pole:

```@example extending
equator = DGG.cellat(grid, 0.0, 0.5)
polar = DGG.cellat(grid, 0.0, 89.5)
DGG.cell_area(grid, equator) / DGG.cell_area(grid, polar)
```

The same ratio over the whole level:

```@example extending
using GeoMakie, GLMakie
using DiscreteGlobalGridsVisualization: dggpoly, dggpoly!
GLMakie.activate!(inline = true)

cells = DGG.CellVector(grid)
areas = [DGG.cell_area(grid, c) * 6371.0^2 for c in cells]

fig = Figure(size = (760, 400))
ax = GeoAxis(fig[1, 1]; dest = "+proj=moll", xticklabelsvisible = false,
    title = "LonLatSystem level 4, cell area")
plt = dggpoly!(ax, cells; color = areas, colormap = :viridis,
    strokecolor = ("#ffffff", 0.35), strokewidth = 0.15)
Colorbar(fig[1, 2], plt; label = "km²")
fig
```

`levelfor` picks the level whose cells come nearest a size you name, and
`cellsize` answers the reverse question:

```@example extending
DGG.levelfor(sys, 100_000), round(Int, DGG.cellsize(sys, 4))
```

## Query and regrid on the lattice

`Intersects` selects the cells that meet a region; `Within`, the cells wholly
inside it:

```@example extending
import Extents
alps = Extents.Extent(X = (5.0, 16.0), Y = (44.0, 48.0))
length(DGG.query(grid, DGG.Intersects(alps))),
    length(DGG.query(grid, DGG.Within(alps)))
```

`regrid` puts a raster on the lattice with the same call every other page in
these docs uses, and returns a `Raster` whose one spatial dimension is `Cells`:

```@example extending
using Rasters
lon, lat = -179.5:1.0:179.5, -89.5:1.0:89.5
ras = Raster([20 - 0.4abs(y) + 3sind(x) for x in lon, y in lat],
    (X(lon; sampling = Rasters.Intervals(Rasters.Center())),
     Y(lat; sampling = Rasters.Intervals(Rasters.Center()))))
field = DGG.regrid(ras; to = grid)
```

`plan_regrid` is the same work stopped one step short, so the weights can be
reused across many rasters:

```@example extending
DGG.plan_regrid(ras; to = grid)
```

That `Cells` dimension carries the lookup, so the raster answers a lon/lat point
through the same generic `cellat`:

```@example extending
import DimensionalData as DD
field[DGG.Cells(DD.Contains((10.5, 46.5)))]
```

and a region selector:

```@example extending
field[DGG.Cells(DGG.Covering(alps))]
```

```@example extending
fig2 = Figure(size = (760, 400))
ax2 = GeoAxis(fig2[1, 1]; dest = "+proj=moll", xticklabelsvisible = false,
    title = "A one-degree raster, regridded")
plt2 = dggpoly!(ax2, field; color = Float32.(parent(field)),
    colormap = :thermal, strokewidth = 0)
Colorbar(fig2[1, 2], plt2; label = "field")
fig2
```

## Run a stencil, and see the pole

`neighbors` runs a geometric search here, from the generic fallback.
`ring(g, c, k)` is the `k`-th ring alone and `neighbors(g, c, k)` is everything
out to it:

```@example extending
c = DGG.cellat(grid, 10.5, 46.5)
length(DGG.neighbors(grid, c)),
    length(DGG.ring(grid, c, 2)),
    length(DGG.neighbors(grid, c, 2))
```

A stencil runs on that unchanged. `mapneighbors` calls a kernel once per cell
with `f(cell, value, neighbours)` and threads by default:

```@example extending
smoothed = DGG.mapneighbors((c, x, nbs) -> (x + sum(nbs)) / (1 + length(nbs)),
    cells, parent(field))
extrema(parent(field)), extrema(smoothed)
```

Adjacency is `Vertex()` by default: cells that touch at a corner are neighbours.
Every cell of a polar row touches the pole, so each one is a neighbour of every
other cell in that row.

```@example extending
table = DGG.adjacency(grid)
sort(unique(length.(table)))
```

Eight neighbours nearly everywhere, and 130 in the polar rows: the other 127
cells of the row, all touching at the pole, plus the three below. Vertex
adjacency means exactly this on a lattice with a polar row, and every
neighbourhood sweep sees it:

```@example extending
extrema(DGG.mapneighbors((c, x, nbs) -> length(nbs), cells, parent(field)))
```

A stencil written for this grid therefore needs a special case at the pole,
where a DGGS such as IGeo7 keeps every cell's neighbour count bounded.

## Answer point location in closed form

The lattice inverts analytically: longitude and latitude give the column and row
directly. Writing `cellat` on the level grid replaces the tree descent, and
[`has_direct_location`](@ref) tells the engines that the closed form exists.

```@example extending
const LonLatGrid = DGG.HierarchicalLevelGrid{LonLatSystem}
const FROM_SPHERE = GO.UnitSpherical.GeographicFromUnitSphere()

function DGG.cellat(g::LonLatGrid, p::DGG.UnitSphericalPoint)
    lon, lat = FROM_SPHERE(p)
    l = DGG.level(g)
    col = clamp(floor(Int, (lon + 180) * ncolumns(l) / 360) + 1, 1, ncolumns(l))
    row = clamp(floor(Int, (lat + 90) * nrows(l) / 180) + 1, 1, nrows(l))
    return LonLatCell(l, row, col)
end

DGG.has_direct_location(::LonLatSystem) = true

DGG.cellat(grid, 10.5, 46.5)
```

An optional method must return what the generic implementation returns, and the
conformance suite is where that is checked. The `cellat` law group opens, and
the grid run comes back with no skips at all:

```@example extending
test_grid_interface(grid)
nothing # hide
```

## Write and read a DGGS store

A DGGS store holds cell ids as integers a reader walks without the system in
hand, so `dggwrite` asks for four more methods and for a store name. Called
before they exist, it reports which:

```@example extending
using Zarr
try
    DGG.dggwrite(joinpath(mktempdir(), "lonlat.zarr"), DD.DimStack((; field)))
catch err
    showerror(stdout, err)
end
```

The four methods go on the level grid, and each is a one-liner here because
`rawid` is already the level's dense index:

| method | answers |
|---|---|
| `idrank(g, id)` | the zero-based rank of an integer id, total over the integer type |
| `idselect(g, r)` | the id at a rank |
| `idvalid(g, id)` | whether an integer names a cell of this grid |
| `idcell(g, id)` | the typed id an integer stands for |

```@example extending
DGG.idrank(g::LonLatGrid, id::Integer)  = clamp(Int(id) - 1, 0, DGG.ncells(g))
DGG.idselect(g::LonLatGrid, r::Integer) = Int(r) + 1
DGG.idvalid(g::LonLatGrid, id::Integer) = 1 <= id <= DGG.ncells(g)
DGG.idcell(g::LonLatGrid, id::Integer)  = DGG.cellindex(g, Int(id))

DGG.register_grid!("lonlat",
    DGG.GridReference("lonlat", LonLatSystem(), :linear, (:linear,)))
nothing # hide
```

With those, the round trip works and the system comes back out of the file:

```@example extending
path = DGG.dggwrite(joinpath(mktempdir(), "lonlat.zarr"), DD.DimStack((; field)))
store = DGG.dggread(path)
DD.lookup(store[:field], DGG.Cells)
```

The store answers the same region selector the in-memory raster did:

```@example extending
store[:field][DGG.Cells(DGG.Covering(alps))]
```

## Fast paths left to write

Four groups are still generic on this system, roughly in order of payoff:

  - **`neighbors` and `ring`** — lattice arithmetic, replacing a geometric
    search. [`maxneighbors`](@ref) declares the ceiling that lets the whole
    neighbourhood family use stack-allocated containers; the polar row's 130
    neighbours are what this lattice would have to declare, so a bound of 8 is
    out of reach.
  - **[`descendant_range`](@ref) with `has_sorted_subtrees(sys) = true`** — a
    range test for subtree membership, which the halo and grow engines use in
    place of scanning a level. Row-major order scatters a subtree across rows;
    a Morton or Hilbert index would make it contiguous, so this system keeps
    [`has_sorted_subtrees`](@ref) at `false`.
  - **[`node_extent`](@ref)** — a tight bound per cell, replacing the inflated
    bounding cap the default computes.
  - **`ancestor` and `descendants`** — arithmetic on the row and column,
    replacing a walk up `parent` or down `children`.

Re-run both suites after each: the skips turn into checks, as `cellat` did.
