# DiscreteGlobalGrids.jl

Discrete global grid systems (DGGS) for the Julia geo ecosystem: six systems —
IGEO7, H3, HEALPix, A5, S2, ISEA4R — behind one small interface, with every
algorithm written against the interface exactly once.

## The two tiers

`AbstractGrid` is one finite collection of cells on the sphere: a complete DGGS
level, a regional subset of one, or a standalone grid with no hierarchy at all.
Four required methods — `ncells`, `cellindex`, `cell_boundary`, `cell_centroid`
— buy the whole generic surface: `cell_polygon`, `cell_area`, `cell_extent`,
`cellat`, `neighbors`, `ring`, `treeify`, `query`.

`AbstractHierarchicalGridSystem` adds analytic parent/child structure. It is the
fast-path tier: tree pruning under the covering law of `node_extent`, contiguous
`descendant_range`s, sublinear `query`. Hierarchy is always an optimisation,
never a semantic — a system whose override disagrees with the generic answer is
wrong, not fast.

A bare `Int` argument is always a **position** in `1:ncells(grid)`. A typed
`AbstractCellIndex` is always an **identity**, self-describing about its level,
so no call passes a level and an id side by side. All internal geometry is on
the unit sphere, as `GeometryOps.UnitSphericalPoint`; longitude and latitude
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

# Positions <-> identities.
c = DGG.cellindex(grid, 1000)                     # a typed cell id
DGG.cellposition(grid, c)                         # 1000
DGG.level(c)                                      # 4

# Geometry, on the unit sphere.
DGG.cell_centroid(grid, c)
DGG.cell_area(grid, c)                            # steradians

# Location and topology. `(lon, lat)` wrappers take degrees.
DGG.cellat(grid, 8.5, 47.4)
DGG.neighbors(grid, c)                            # ring 1, CCW seen from outside
DGG.ring(grid, c, 2)                              # exactly distance 2

# Spatial queries, with DE9IM predicate types and spherical semantics.
import Extents
DGG.query(grid, DGG.Intersects(Extents.Extent(X = (5, 12), Y = (45, 50))))

# The hierarchy, on the system rather than the grid — an id knows its level.
sys = DGG.HEALPixSystem()
DGG.children(sys, c)
parent(sys, c)
DGG.descendant_range(sys, c, 6)                   # positions in levelgrid(sys, 6)
DGG.subtree_border(sys, c, 6)                     # the rim, O(rim)
```

Swapping `HEALPixSystem()` for `IGeo7System()`, `H3System()`, `A5System()`,
`S2System()` or `ISEA4RSystem()` changes nothing else. `DGG.systems()` lists all
six, and its docstring is the comparison table: cell counts, cell shape,
equal-areaness, and the traits that differ across them.

`treeify`, `ncells` and `getcell` are `ConservativeRegridding.Trees`' own
bindings, extended here, so any grid is a regridding source with no wrapper:

```julia
import ConservativeRegridding as CR
regridder = CR.Regridder(destination, DGG.levelgrid(DGG.HEALPixSystem(), 4))
```

`examples/regridding.jl` is that claim as an assertion-checked script; run it
with `julia -t 4 --project=. examples/regridding.jl` and it exits non-zero if a
check fails. Every script under `examples/` is written that way and runs in this
project environment; the tutorials under `docs/src/tutorials/` are Literate.jl
sources, run by the docs build.

## Layout

| Path | Contents |
|:--|:--|
| `src/interface/` | the type vocabulary and every generic's contract — declarations and trait defaults, no algorithms |
| `src/fallbacks/` | the generic implementations: `HierarchicalLevelGrid`, `PartialGrid`, `AuthalicGrid`/`AuthalicSystem`, `HierarchicalGridCursor`, `MultiOrderCellSet`, the query engine |
| `src/systems/` | one directory per system, plus `src/systems/ISEA/` — the Snyder/icosahedron basis IGEO7 and ISEA4R share |
| `src/Helpers/` | shared allocation-free primitives |
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

Native layers: H3 calls libh3 through `H3_jll`; the other five are pure Julia.
IGEO7 is a clean-room implementation; A5 ports upstream a5's arithmetic; HEALPix,
S2 and ISEA4R are closed-form charts with no external dependency.

No system defines a grid type. All six return `HierarchicalLevelGrid` from
`levelgrid` and attach their fast paths — `cellat`, `neighbors`, `ring`,
`cell_area` — to `HierarchicalLevelGrid{TheSystem}`.

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

```julia
using Pkg
Pkg.activate("path/to/DiscreteGlobalGrids.jl")
Pkg.test()
```

`test/runtests.jl` runs the interface suite, the fallback suites, one suite per
system, and a cross-system suite that sweeps `systems()` so registering a system
grows it automatically. Each is wrapped in its own module, because the systems
share generic vocabulary. The IGEO7 suite validates against recorded DGGRID
output in `test/systems/IGeo7/vectors/` and dominates the count.
**918,836 assertions, ~90 s warm**, with 2 broken — A5's documented
`has_sorted_subtrees` trait skips.

## Provenance

Migrated 2026-08-05 from the `dggs_lookup/` prototype tree in the
`vectordatacubes` workspace.

The IGEO7 implementation is a **clean-room** unit (`src/systems/ISEA/` +
`src/systems/IGeo7/`). It replaces an earlier implementation whose native layer
was ported from an AGPL-licensed reference, which is deliberately **excluded**
here. That reference enters only as an independent **black-box validation
oracle**: the suite checks agreement against dumps of its CLI output, never
against its source. The full audit trail, the 150 MB vector corpus and 24 MB of
reference PDFs stay in `dggs_lookup/`; only the ~9 MB of vectors the suite reads
travel here.
