# DiscreteGlobalGrids.jl

Discrete global grid systems (DGGS) for the Julia geo ecosystem: six systems —
IGEO7, H3, HEALPix, A5, S2, ISEA4R — behind one small interface, with every
algorithm written against the interface exactly once. A seventh,
`CopernicusDEMSystem`, wears the same interface so a DEM tile can be a regrid
*source* without an adapter; it is a raster lattice, not a DGGS, and is not in
`systems()`.

## The two tiers

`AbstractGrid` is one finite collection of cells on the sphere: a complete DGGS
level, a regional subset of one, or a standalone grid with no hierarchy at all.
Four required methods — `ncells`, `cellindex`, `cell_boundary`, `cell_centroid`
— buy the whole generic surface: `cell_polygon`, `cell_area`, `cell_extent`,
`cellat`, `neighbors`, `ring`, `adjacency`, `treeify`, `query`. On a SUBSET of
a level the topology verbs mean the complete level's answer clipped to
membership — omitted, not padded — so a stencil on a region is the same call it
is on the globe, and `halo` names the cells just outside it that the clipping
dropped.

`AbstractHierarchicalGridSystem` adds analytic parent/child structure. It is the
fast-path tier: tree pruning under the covering law of `node_extent`, contiguous
`descendant_range`s, sublinear `query`. Hierarchy is always an optimisation,
never a semantic — a system whose override disagrees with the generic answer is
wrong, not fast.

A bare `Int` argument is always an **index** in `1:ncells(grid)` — a
local index into that collection's own storage. A typed
`AbstractCellIndex` is always an **identity**, self-describing about its level,
so no call passes a level and an id side by side. All internal geometry is on
the unit sphere, as `UnitSphericalPoint` — `GeometryOps`', re-exported here so
an implementor writes it with no module path; longitude and latitude
appear only in explicitly named converting wrappers, in degrees.

## Setup

Julia ≥ 1.11. `ConservativeRegridding`, `GeometryOps` and `GeometryOpsCore`
resolve from git branches through `[sources]` in `Project.toml`, so the usual
`julia --project=.` followed by `Pkg.instantiate()` is all that is needed — no
sibling checkouts.

`Project.toml` declares a workspace: `test/`, `docs/` and
`lib/DiscreteGlobalGridsConformanceTesting/` share one manifest.

## Quick start

```julia
import DiscreteGlobalGrids as DGG

# A complete level is the entry point. `levelgrid` returns a
# `HierarchicalLevelGrid`, which holds the system and the level and nothing else.
grid = DGG.levelgrid(DGG.HEALPixSystem(), 4)
DGG.ncells(grid)                                  # 3072

# Indices <-> identities.
c = DGG.cellindex(grid, 1000)                     # a typed cell id
DGG.globalindex(grid, c)                          # 1000
DGG.level(c)                                      # 4

# Geometry, on the unit sphere.
DGG.cell_centroid(grid, c)
DGG.cell_area(grid, c)                            # steradians

# Location and topology. `(lon, lat)` wrappers take degrees.
DGG.cellat(grid, 8.5, 47.4)
DGG.neighbors(grid, c)                            # ring 1, CCW seen from outside
DGG.ring(grid, c, 2)                              # exactly distance 2
DGG.adjacency(grid)                               # every cell's stencil, as indices

# Spatial queries, with DE9IM predicate types and spherical semantics.
import Extents
DGG.query(grid, DGG.Intersects(Extents.Extent(X = (5, 12), Y = (45, 50))))

# The hierarchy, on the system rather than the grid — an id knows its level.
sys = DGG.HEALPixSystem()
DGG.children(sys, c)
parent(sys, c)
DGG.descendant_range(sys, c, 6)                   # indices in levelgrid(sys, 6)
region = DGG.subtree(sys, c, 6)                   # the subtree as an ordinary grid
DGG.border(region)                                # the border, O(border)
DGG.halo(region)                                  # the cells just OUTSIDE it
DGG.interior(region)                              # and the complement of the border

# All three are lazy: a halo can dwarf the border it wraps, so collecting is
# always the caller's call, and holes in the middle of a region count.
for p in DGG.halo(region)                         # O(depth) memory, resumable
    break
end
DGG.halo(region; cells = true)                    # ids rather than indices
```

Swapping `HEALPixSystem()` for `IGeo7System()`, `H3System()`, `A5System()`,
`S2System()` or `ISEA4RSystem()` changes nothing else, and `AuthalicSystem`
wraps any of them — except A5, whose geometry is geodetic already — to read
geometry at geodetic latitude. `DGG.systems()` lists
all six, and its docstring is the comparison table: cell counts, cell shape,
equal-areaness, and the traits that differ across them.

## Regions, kept compressed

A multi-order coverage is the coarsest cells that cover a target, at mixed
levels. Its two keywords are two modes, not a bound and a hint, and exactly one
of them is given:

```julia
igeo = DGG.IGeo7System()
ext = Extents.Extent(X = (6.0, 10.5), Y = (45.8, 47.8))   # a Switzerland-shaped box

accurate = DGG.query(igeo, DGG.MultiOrderCoverage(ext); level = 7)     # accuracy first
budget   = DGG.query(igeo, DGG.MultiOrderCoverage(ext); maxcells = 10) # cardinality first
```

`level` refines everything the target's boundary crosses down to a fixed depth
and lets the cell count fall where it may: 335 entries over levels 4 to 7 here.
`maxcells` refines the crossing cells breadth first, coarsest up, and stops when
the next replacement would not fit — "ten cells that cover California";
`maxlevel` bounds how deep the budget may descend.

Covering is a statement about the leaves: at the deepest level, every cell
meeting the target is a member of the set or a descendant of one. The union of
the *drawn* cells is that region only where refinement is congruent — HEALPix,
S2 and ISEA4R tile, IGEO7 and H3 leave slivers, A5 more.

`CellVector` reads either set as a lazy `AbstractVector` of ascending ids at one
level, stored as the leaf index windows they occupy rather than as the ids:

```julia
cv = DGG.CellVector(accurate)              # or of a grid, or of an explicit id vector
cv[3]                                      # the third id; nothing is materialised
DGG.localindex(cv, cv[3])                  # the inverse, or `nothing`
DGG.cellat(cv, 8.5, 47.4)                  # the cell of `cv` a point falls in
DGG.covering(cv, Extents.Extent(X = (7, 8), Y = (46, 47)))   # a sub-region
DGG.PartialGrid(cv)                        # read as a grid, O(1)
```

Memory is O(#windows), not O(#cells): 98 windows for the 1319 level-7 cells
those 335 entries expand to, and the *same* 98 for
`CellVector(accurate; level = 10)`, which names 452,417. `intersect` and
`issubset` are Base's, answered over the windows. `PartialGrid(cv)` is the
handshake with anything that wants a grid: index `k` of the grid is index
`k` of the vector, so data laid out against the vector needs no permutation.

`CellLookup` is that same type wearing a `DimensionalData.Lookup` hat, and the
`Cells` dimension puts it on a cube:

```julia
import DimensionalData as DD
using Statistics

lk = DGG.CellLookup(accurate)
A = DD.DimArray(rand(length(lk)), DGG.Cells(lk))

A[DGG.Cells(DD.At(lk[3]))]                    # a typed cell id
A[DGG.Cells(DD.Contains(8.0, 46.5))]          # a lon/lat point, through `cellat`
mean(A[DGG.Cells(DGG.Covering(ext))])         # a region; the view's axis is a CellLookup again
```

`At` and `Contains` stay qualified as DimensionalData's: this package exports
DE9IM's `Contains`, a predicate about two geometries rather than a selector, and
the two names must not collide. `Covering`, which DimensionalData has no
spelling for, is exported.

## Regridding

`regrid` moves a field from one cell collection onto another, first-order
conservative by default. A grid, a `CellVector`, a `CellLookup`, a
`MultiOrderCellSet` or a bare system is a destination as it stands:

```julia
grid = DGG.levelgrid(DGG.HEALPixSystem(), 6)
temps = DD.DimArray(rand(360, 180), (DD.X(-179.5:179.5), DD.Y(-89.5:89.5)))

tavg = DGG.regrid(temps; to = grid)                 # onto a whole level
sub  = DGG.regrid(temps; to = cv)                   # onto the region above
auto = DGG.regrid(temps; to = DGG.HEALPixSystem())  # level matched by cell area

# The other direction. `from` is required whenever the source is not a lon/lat
# raster: a grid, a `CellVector` or a `CellLookup` is a source as it stands, and
# a raster or a dimension tuple names the lon/lat destination.
back = DGG.regrid(tavg; to = temps, from = grid)

plan = DGG.plan_regrid(temps; to = grid)            # the operator alone, reusable
DGG.regrid(temps, plan)                             # ... applied, no keywords left
```

Weights are geometry: building the plan is the expensive half and reads no data,
and a plan carries the method, both spaces and the missing policy, so applying
one takes no keyword arguments. `missingpolicy` says what a partly covered
destination cell holds — `DGG.Weighted(t)` the coverage-normalised mean,
blanking cells covered less than `t`, `DGG.Extensive()` the raw conservative
sum. The destination's cells replace the source's spatial dimensions and every
other dimension passes through. Each space carries the manifold it computes on — the
unit sphere, everywhere here — and a pair that disagrees is an error rather than
a silent rescaling by `R^2`.

`treeify`, `ncells` and `getcell` are `ConservativeRegridding.Trees`' own
bindings, extended here, so the cell tree a regrid descends is the grid's own
with nothing wrapped around it.

A DGGS as the **source** conserves to `1e-13` on all six systems and on the
authalic wrap. A DGGS as the **destination** conserves only where the
destination cells' rings are convex — IGEO7 and S2 at every level the suite
sweeps, H3 at the even ones — because the clipper's Sutherland–Hodgman is an
intersection only against a convex clip window, and the destination is always
that window. Convexity is a property of the system *and* the level, so
`test/systems/crosssystem/regridding_conservation.jl` measures the rings rather
than carrying a list: it asserts the law where they come out convex, marks the
rest `@test_broken`, and names the upstream fix that closes them.

## Going further

Every script under `examples/` is an assertion-checked demo that exits non-zero
if a check fails — `julia -t 4 --project=. examples/regridding.jl`. The one
exception is `examples/copernicus_dem.jl`, which reads a COG and so needs the
docs environment: `julia -t auto --project=docs examples/copernicus_dem.jl`. The
seven tutorials under `docs/src/tutorials/` are Literate.jl sources run by the docs
build, each the shortest honest path to one result; `docs/src/index.md` lists
them and `docs/src/all_dggs.md` draws every system.

## Layout

| Path | Contents |
|:--|:--|
| `src/interface/` | the type vocabulary and every generic's contract — declarations and trait defaults, no algorithms |
| `src/fallbacks/` | the generic implementations: `HierarchicalLevelGrid`, `PartialGrid`, `AuthalicGrid`/`AuthalicSystem`, `HierarchicalGridCursor`, `MultiOrderCellSet`, `CellVector`, the subtree iterators, the clipped stencils, the query engine |
| `src/dimensionaldata.jl` | the cube face of `CellVector`: `CellLookup`, `Cells`, `Covering` |
| `src/systems/` | one directory per system, plus `src/systems/ISEA/` — the Snyder/icosahedron basis IGEO7 and ISEA4R share |
| `src/core/`, `src/Helpers/` | the authalic manifold pair, and shared allocation-free primitives |
| `lib/GlobalRegridding/` | the generic regridding engine — spaces, weights, plans, the lazy executor — consumed by the main package for `regrid`/`plan_regrid` |
| `lib/DiscreteGlobalGridsConformanceTesting/` | the test-only workspace package that makes the contracts executable |

## The systems

| System | Levels | Cells at level `l` | Cell shape | Equal-area | Canonical id |
|:--|:--|:--|:--|:--|:--|
| `IGeo7System` | `0:19` | `10·7^l + 2` | hexagons + 12 pentagons | yes | `Z7Cell` |
| `H3System` | `0:15` | `120·7^l + 2` | hexagons + 12 pentagons | no | `H3Cell` |
| `HEALPixSystem` | `0:29` | `12·4^l` | curvilinear diamonds | yes | `LevelIndex` (nested) |
| `A5System` | `0:29` | `12`, `60`, then `60·4^(l-1)` | pentagons | yes | `A5Cell` |
| `S2System` | `0:30` | `6·4^l` | geodesic quadrilaterals | no | `LevelIndex` |
| `ISEA4RSystem` | `0:29` | `10·4^l` | rhombi on ten diamonds | yes | `LevelIndex` |
| `CopernicusDEMSystem` | `0:1` | `64 800` tiles, then `360·N·Σ_r ncols(r)` pixels | lon/lat boxes | no | `LevelIndex` |

`CopernicusDEMSystem` is the odd row: the Copernicus DEM raster lattice —
1°×1° tiles over pixel-is-point rasters whose longitude spacing steps down at
latitude 50/60/70/80/85 — not a tessellation designed to be a DGGS, so it is
**not** in `systems()`. Its `neighbors` and `ring` are closed form at every
cell of both levels, decided in exact integer arithmetic rather than by
matching corners within a tolerance. It is here as the *source* side of a
DEM-to-DGGS regrid: `PartialGrid(CopernicusDEMSystem(90), tile, 1)` is one AWS
COG tile as an ordinary `AbstractGrid`, with no adapter and no corner matrix.
`examples/copernicus_dem.jl` moves one real tile onto IGEO7 and HEALPix and
checks every conservation law the move owes.

Native layers: H3 calls libh3 through `H3_jll`; every other system is pure Julia.
IGEO7 is a clean-room implementation but for one ported adjacency kernel (see
[Provenance](#provenance)); A5 ports upstream a5's arithmetic; HEALPix, S2 and
ISEA4R are closed-form charts with no external dependency.

No system defines a grid type. All seven return `HierarchicalLevelGrid` from
`levelgrid` and attach their fast paths to `HierarchicalLevelGrid{TheSystem}`:
`cellat`, `neighbors` and `ring` on all seven, and `cell_area` on the three
whose exact area is a closed form the published boundary only approximates
(HEALPix, ISEA4R, CopernicusDEM); the other four take the generic spherical area
of that boundary. Among the six in `systems()`, `border` over a rooted subtree is
an `O(border)` automaton on every one but A5, which walks the whole subtree;
`interior` shares that walk and emits the branches it prunes. Both are resumable
`EdgeCellIterator` / `InnerCellIterator` walks in `O(depth)` memory. A5 is also the one system without
`has_sorted_subtrees`, so `level_ranges` throws there and everything that would
use it takes the selection branch instead.

`halo` is the outside of that same boundary and is built the same way: a
resumable `SubtreeHaloIterator` in `O(depth)` memory, so a prefix of a deep halo
costs what the prefix costs and not what the ring would. HEALPix,
S2 and ISEA4R walk the band around their square block, one pruned quadtree
descent per face the halo touches; IGeo7 and H3 seed each neighbour's border
automaton with a calibrated arc and walk that; A5, again for want of
`descendant_range`, scans the target level. Only two of the seven engines count
in closed form — depth zero, which is the one-ring already in hand, and the
square in-face band, which is `4·side + 4` — so only those declare a `length`;
everywhere else `IteratorSize` is `SizeUnknown()` and there is no `length`
method at all, deliberately, because a `length` that walked the halo to answer
is the thing the design forbids. The same verb asks the same question of any
region — a `PartialGrid`, a `CellVector`, a `CellLookup` or a complete grid — so
a cell punched out of the middle of a subset joins that subset's halo.

The system submodules (`DiscreteGlobalGrids.H3` and friends) are deliberately
**not** exported: `H3`, `HEALPix`, `A5` and `S2` are also the names of
registered packages. Reach past the interface through the qualified module —
`DiscreteGlobalGrids.IGeo7.z7_string`, `DiscreteGlobalGrids.H3.H3Native`. Note
the capitalisation of `HEALPix`: it never shadows the registered Healpix.jl,
which this package's tests use as an independent oracle.

ISEA4R's diamond pairing, numbering and axis orientations are this package's own
convention with **no external oracle**; identifier compatibility with any
external ISEA4R product is not claimed. S2's native 64-bit `s2_cellid` is not
offered as an index scheme — the scaffold ordinal is canonical.

## Conformance

Every interface law has a property test a third-party implementor can run with
two calls. The suites live in a separate workspace package so `Test` never
becomes a dependency of `DiscreteGlobalGrids`:

```julia
using DiscreteGlobalGridsConformanceTesting

test_grid_interface(grid)          # the AbstractGrid contract
test_hierarchical_system(sys)      # the hierarchy contract, incl. the covering law
```

Each system's suite runs both against itself. The harness has its own tests
under `lib/DiscreteGlobalGridsConformanceTesting/test/`, which check that it
*catches* deliberately broken mock implementations rather than merely running.

## Tests

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

`test/runtests.jl` runs the interface suite, the fallback suites, one suite per
system, and a cross-system suite that sweeps `systems()` so registering a system
grows it automatically. Each is wrapped in its own module, because the systems
share generic vocabulary. The IGEO7 suite validates against recorded DGGRID
output in `test/systems/IGeo7/vectors/` and dominates the count.
**987,153 assertions, ~7m30s**, with 17 broken — all of them
destination-direction conservation arms in `regridding_conservation.jl`,
measured per system and level, waiting on the upstream clipper fix.

`lib/GlobalRegridding/` carries its own suite, which the root `Pkg.test()` does
not run: `julia --project=lib/GlobalRegridding -e 'using Pkg; Pkg.test()'`.

## Provenance

Migrated 2026-08-05 from the `dggs_lookup/` prototype tree in the
`vectordatacubes` workspace.

The IGEO7 implementation is a **clean-room** unit (`src/systems/ISEA/` +
`src/systems/IGeo7/`), with the one marked exception below. It replaces an
earlier implementation whose native layer was ported from an AGPL-licensed
reference, which is deliberately **excluded** here. That reference enters only
as an independent **black-box validation oracle**: the suite checks agreement
against dumps of its CLI output, never against its source. The full audit trail,
the 150 MB vector corpus and 24 MB of reference PDFs stay in `dggs_lookup/`; only
the ~9 MB of vectors the suite reads travel here.

### Exception: `src/systems/IGeo7/gbt.jl`

One file is source-level reuse rather than black-box validation, and says so in
its own header. The one-ring adjacency kernel is **ported from
[IGEO7.jl](https://github.com/allixender/IGEO7.jl)** (`src/IGEO7.jl`) by
**Alexander Kmoch** — specifically `get_neighbour`, `get_neighbours`,
`first_non_zero`, and the eight tables those read (`BASE_CELL_NEIGHBOURS`,
`EXCLUDE_NEIGHBOURS`, `ROTATIONS`, `POLE_0_ROTATIONS`, and the four GBT addition
tables). The reuse is covered by a **licence grant from Alexander Kmoch** to this
package.

> **TODO (Anshul):** record the grant's actual terms. This repository ships no
> `LICENSE` file yet; the terms belong next to it. Whether the grant is a
> relicence, a dual licence or something else is not recorded anywhere in this
> tree, so nothing beyond the fact of the grant is asserted here.

Attribution is owed whatever the terms are, so the file header names the upstream
repository, file, functions and author, and each ported table's docstring points
back at that header. What the port changed is shape, not arithmetic; the
counterclockwise ordering it emits in (`_encode_lattice_rot`, `z7grid.jl`) is
this package's own, because upstream defines no rotational order.

The port is **not trusted, it is checked**, against this package's own
independent implementation of the same question:
`IGeo7._cell_neighbors_ccw_geometric` (`src/systems/IGeo7/z7grid.jl`) derives
adjacency from the oracle-validated lattice and decoder instead of from digit
arithmetic. Two implementations that share no reasoning are the strongest
evidence available that either is right, so the geometric one stays in the tree
as the differential oracle even though nothing on the hot path calls it. Testset
`9b` of `test/systems/IGeo7/runtests.jl` pins the two together — same ids, same
order — on every cell of levels 0–3, on samples through level 19, across all
twelve pentagon chains and their two-ring neighbourhoods, and over the whole
`neighbors`/`ring` disc to `k = 3`.
