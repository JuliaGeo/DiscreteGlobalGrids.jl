# DiscreteGlobalGrids.jl

Discrete global grid systems (DGGS) for the Julia geo ecosystem — 14 registered
systems, six of which answer geometry, four of which are wired all the way into
one generic spatial-tree family shared with `ConservativeRegridding` and
`GeometryOps.SpatialTreeInterface`.

There is no per-system tree code. A system supplies operations — child
enumeration, parent, descendant interval, unit-sphere boundary — and the package
supplies the grid types, the traversal, the pruning and the caps. Adding a system
to the tree layer is one `<X>Kernel.jl` file and its include line. The kernel is
wired incrementally, and anything unwired throws `NotPortedError` rather than
guessing.

## Setup

`ConservativeRegridding` resolves via a `[sources]` path entry to the sibling
`../ConservativeRegridding.jl` checkout (its `Trees` submodule is not yet in a
registered release), so this package needs **Julia ≥ 1.11 and that sibling
checkout present**. Work in the package environment: `julia --project=.`.

## Quick start

```julia
using DiscreteGlobalGrids
using DiscreteGlobalGrids.HEALPix.HealpixLookups: HealpixLookup
using DiscreteGlobalGrids.H3.H3Lookups: H3Lookup
import DiscreteGlobalGrids.H3.H3Native
import ConservativeRegridding as CR
import GeometryOps as GO
import GeometryOps: SpatialTreeInterface as STI

to_sphere = GO.UnitSpherical.UnitSphereFromGeographic()
destination = [to_sphere((lon, lat)) for lon in range(0, 360; length=37),
                                         lat in range(-90, 90; length=19)]

# 1. A whole level as a conservative-regridding source. The DGGS singleton is the
#    only system-specific token in the call; swap it and nothing else moves.
regridder = CR.Regridder(destination, DGGSGrid(HEALPixDGGS(), 3))
size(regridder.intersections)                    # (648, 768) — (destination, source)

# 2. A stored cell axis becomes a spatial tree in one line. The grid holds
#    `lookup.data` itself — never copied, never reordered — so tree leaf i is
#    lookup position i, and a `DimArray` over that dimension lines up with a
#    `Regridder` column-for-column, with no permutation anywhere.
lookup = HealpixLookup(collect(Int64, 3000:3999); level=6)
tree = treeify(lookup)                           # via DGGSPartialGrid(lookup)
ncells(tree) == length(lookup.data)              # true

# 3. A globe-complete axis costs nothing to build: `DGGSGlobeIds` is lazy, so the
#    ids are computed on demand rather than stored. Indexing it degrades to an
#    ordinary explicit-id lookup.
globe = H3Lookup(DGGSGlobeIds(H3DGGS(), 9))
length(globe)                                    # 4842432842, allocated: none
ncells(treeify(globe))                           # same — a dense DGGSGrid cursor
globe[1:10]                                      # an H3Lookup backed by a Vector

# 4. One cell's subtree as a chunk — O(subtree) to build, never O(globe) — and a
#    cap query over it. Hits are positions in `chunk.ids`.
chunk = subtree_grid(H3DGGS(), H3Native.lonlat_to_cell(10.0, 45.0, 2);
                     root_level=2, leaf_level=5)
cap = GO.UnitSpherical.SphericalCap(to_sphere((10.0, 45.0)), 0.02)
hits = STI.query(treeify(chunk), intersects_cap(cap))
DiscreteGlobalGrids.cell_center(H3DGGS(), 5, chunk.ids[first(hits)])
```

Longer, assertion-checked versions live in `examples/` — `regridding.jl`,
`dimensionaldata.jl`, `chunked_h3.jl` — each running standalone under
`julia -t 4 --project=. examples/<name>.jl` and exiting non-zero if a check
fails. (`examples/healpix_demo.jl` additionally wants NaturalEarth and Rasters,
which are intentionally not dependencies.)

## Naming rules

Nothing system-level is re-exported at the top level: the systems intentionally
share generic vocabulary (`cell_center`, `cell_boundary`, `lonlat_to_cell`,
`Touching`, ...), so reach them through their submodule.

```julia
using DiscreteGlobalGrids
using DiscreteGlobalGrids.H3.H3Lookups   # cell_center, lonlat_to_cell, ...
```

The same rule runs the other way, which is why the kernel generics `num_cells`,
`cell_boundary` and `cell_center` are *not* exported even though every other one
is: a submodule already claims each of those names, and a `using
DiscreteGlobalGrids` next to a `using ...H3.H3Native` must not make either
ambiguous. Call them qualified — `DiscreteGlobalGrids.num_cells(system, level)`.

Note the capitalization of `HEALPix`: the submodule is deliberately *not* named
`Healpix`, so it never shadows the registered Healpix.jl package it imports.

## Grid systems

One submodule per system, all the same shape: a native layer, an `<X>Lookups`
`DimensionalData` integration nested inside it, and an `<X>Kernel.jl` wiring the
system into the operations kernel.

| Submodule | System | Native layer |
|:--|:--|:--|
| `A5` | A5 pentagonal DGGS | pure Julia |
| `H3` | H3 hexagonal DGGS | `H3_jll` (libh3) |
| `HEALPix` | HEALPix, as a nested id hierarchy (EOPF/GRID4EARTH conventions) and as a dense face grid | pure-Julia chart kernel; Healpix.jl in the lookup layer |
| `IGeo7` | IGEO7 (ISEA7H + Z7), clean-room implementation | pure Julia, stdlib-only core; geometry from `ISEA` |
| `ISEA4R` | ISEA4R, as a dense ten-diamond grid and as `ISEA4RDGGS` cell geometry | pure Julia; Snyder charts from `ISEA` |
| `ISEA9R` | ISEA9R, the same ten diamonds at aperture 9, and `ISEA9RDGGS` cell geometry | pure Julia; **imports `ISEA4R`'s chart unchanged** — the rhombus chart carries no aperture — plus base-9 index maps |
| `S2` | S2, as a dense cube-face grid and as `S2DGGS` cell geometry | pure-Julia closed forms; no s2geometry dependency |

`S2`, `ISEA4R` and `ISEA9R` have no `<X>Lookups` module: they are chart systems,
and a dimension of stored ids needs an id hierarchy none has yet. `ISEA9R` is the
thinnest instance of the shared face-grid layer: its `chart.jl` *imports*
`ISEA4R`'s chart rather than defining one, a claim
`test/ISEA9R/test_delegation.jl` checks with `===` on the function objects and
then bitwise on the two systems' grids and Regridders.

The ten-*root* layout both rhombic systems use is normative for ISEA9R — OGC
21-038r1 Annex B.2, *"The ten root rhombuses are formed by combining two
icosahedron triangles at their base"* — and DGGAL's
`RhombicIcosahedral9R::countZones(level)` returns `10 * 9^level`. The ten-diamond
*numbering*, the in-diamond axes and the in-diamond index are package
conventions with **no external oracle**, so no DGGRID/DGGAL/SST identifier
compatibility is claimed for either system — see
`docs/design/isea4r_diamond_layout.md`. `zonal` and `stencil` currently exist
only in `HealpixLookups`.

### Which systems answer what

| System | id hierarchy + grids | cell geometry | id the geometry takes |
|:--|:--|:--|:--|
| `HEALPixDGGS`, `H3DGGS`, `IGEO7DGGS`, `A5DGGS` | yes | yes | the system's canonical id |
| `S2DGGS` | no | yes | scaffold ordinal `face * 4^level + hilbert_position` |
| `ISEA4RDGGS` | no | yes | canonical `isea4r_ordinal` `diamond * 4^level + morton_position` |
| `ISEA9RDGGS` | no | yes | canonical `isea9r_ordinal` `diamond * 9^level + morton_position` (base-9 Morton) |
| the other seven, incl. `RHEALPixDGGS` | no | no | — |

How much of the kernel each wired system gets for free is decided by four
id-model traits (`cell_id_type`, `has_ordinal_ids`, `has_descendant_ranges`,
`has_exact_subtree_cap`) plus its `cell_cap_inflation` budget. Those values, what
each trait means, what the "registry-only" systems still record, and why the
cursor has three descent modes: `docs/reference/kernel_wiring.md`.

## Reference

| Page | Contents |
|:--|:--|
| `docs/reference/common_layer.md` | the abstract interface, system registry, operations kernel, grid/tree family, lookups and `DGGSGlobeIds`, shared submodules |
| `docs/reference/kernel_wiring.md` | id-model traits in full, the three cursor descent modes, the support matrix explained |
| `docs/reference/healpix_representations.md` | HEALPix as an id hierarchy and as a face grid over one chart kernel |
| `docs/design/` | design records: `full_globe_lookups.md`, `face_grid_generalization.md`, `isea4r_diamond_layout.md` |
| `docs/IGeo7/` | IGEO7 clean-room record — `CLEANROOM.md`, `PROVENANCE.md`, spec |

## Tests

```julia
using Pkg
Pkg.activate("path/to/DiscreteGlobalGrids.jl")
Pkg.test()
```

`test/runtests.jl` aggregates one suite per unit — `test/test_helpers.jl` and
`test/core/`, `test/A5/`, `test/H3/`, `test/HEALPix/`, `test/IGeo7/`,
`test/ISEA4R/`, `test/ISEA9R/`, `test/S2/` — each wrapped in its own module so the generic
vocabulary the systems share cannot collide across suites. The IGEO7 suite
validates against the oracle vectors in `test/IGeo7/vectors/` and dominates the
count. **513,337 assertions, ~80 s warm.**

## Provenance

Migrated 2026-08-05 from the `dggs_lookup/` prototype tree in the
`vectordatacubes` workspace.

The IGEO7 implementation shipped here is a **clean-room** unit (`ISEA` +
`IGeo7`). It replaces wholesale an earlier implementation whose native layer was
ported from an AGPL-licensed reference implementation, which is deliberately
**excluded** from this package. That reference implementation enters only as an
independent **black-box validation oracle**: the suite checks agreement against
dumps of its CLI output, never against its source. The clean-room records are
`docs/IGeo7/CLEANROOM.md` and `docs/IGeo7/PROVENANCE.md` (written under the
unit's former name, `IGeo7Clean`); the full oracle audit trail, the 150 MB
vector corpus and 24 MB of reference PDFs stay in `dggs_lookup/` — only the
~9 MB of vectors the suite reads travel here.
