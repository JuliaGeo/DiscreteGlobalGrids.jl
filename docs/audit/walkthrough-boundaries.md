# Newcomer walkthrough: domain workflows and package boundaries

This report covers UC-058 through UC-075 after the first documentation rewrite.
It starts from the public pages named in the inventory. Source and tests are
used only to check the answer a newcomer would need. Documentation gaps and API
gaps remain separate.

**Status words:** **Clear** means the revised public route answers the use-case
questions. **Partial** means useful pieces are documented but a material
question remains. **Unsupported** means the requested complete operation has no
public implementation.

## Findings

### UB-001 — The terrain tutorial does not define the full measurement contract

- **Kind / priority:** DOC GAP / P1
- **Use case:** UC-058
- **Affected files:** `docs/src/tutorials/hydrology.jl`;
  `docs/src/tutorials/out_of_core.jl`; `docs/src/tutorials/stencils.jl`
- **Evidence:** The hydrology tutorial identifies elevation and drop as metres,
  accumulation as square metres, and its custom downhill result as an elevation
  difference rather than slope. It explains NaN cells and pits. It does not say
  whether input DEM values are posts or pixel means, define the slope and
  roughness conventions, or describe flat and pentagon handling in the
  Geomorphometry calls.
- **A newcomer needs:** One table for each demonstrated terrain result: owning
  package, source sampling, units, neighborhood rule, missing-data rule, and
  behavior at pits, flats, boundaries, and exceptional cells.
- **Proposed correction:** Add the table beside the first terrain calculation.
  Link Geomorphometry for algorithms that this package only hosts on a cell
  axis. Keep the custom downhill kernel labelled as an example.

### UB-002 — The Copernicus DEM lattice lacks a public task route

- **Kind / priority:** DOC GAP / P1
- **Use case:** UC-059
- **Affected files:** CopernicusDEM module docstring;
  `docs/src/tutorials/choosing_a_grid.jl`;
  `docs/src/tutorials/hydrology.jl`
- **Evidence:** The module docstring clearly identifies GLO-30 and GLO-90,
  pixel-is-point sampling, two levels, latitude bands, ordering, WGS84, and the
  AWS DGED profile. The hydrology tutorial downloads a raster and transfers it
  to IGeo7; it does not construct `CopernicusDEMSystem` or walk a band
  transition, tile seam, rim, or polar limit. The production runner is under
  `scripts/`, outside the public site.
- **A newcomer needs:** A small local example that constructs both supported
  resolutions, locates a point, identifies its tile parent, and states the
  valid latitude and seam ownership rules.
- **Proposed correction:** Add a focused specialized-system section to the
  chooser or a short Copernicus DEM page. Link the module contract and label
  the distributed runner as experimental repository tooling.

### UB-003 — HEALPix ID conversion does not show how to reorder values

- **Kind / priority:** DOC GAP / P1
- **Use cases:** UC-060, UC-069
- **Affected files:** `docs/src/tutorials/healpix_astronomy.jl`; HEALPix
  system docstrings
- **Evidence:** The tutorial states `nside = 2^level`, `12*nside^2` cells,
  and nested storage. It tells readers to convert ring-order input before
  building a lookup but demonstrates only an already nested vector.
  `reindex` converts cell identifiers; it does not reorder a value vector.
  The tutorial assumes a galactic frame for its synthetic data and mentions a
  FITS reader without naming or testing one.
- **A newcomer needs:** An index-encoded ring vector whose values are explicitly
  permuted into nested cell order, plus a statement that coordinates carry no
  celestial-frame metadata and that FITS I/O belongs to another package.
- **Proposed correction:** Add a small ring-to-nested permutation example and
  name a tested FITS handoff. State the chosen coordinate frame beside the
  `Cells` lookup.

### UB-004 — Main-site plotting discovery stops before installation and limits

- **Kind / priority:** DOC GAP / P2
- **Use cases:** UC-061, UC-062
- **Affected files:** `docs/src/index.md`; `docs/src/all_dggs.md`;
  `lib/DiscreteGlobalGridsVisualization/README.md`
- **Evidence:** The companion README clearly distinguishes `dggpoly`,
  `dggsurface`, and `dggresample`; names backends, missing-color behavior,
  seams and poles, interpolation, and experimental status. The main site names
  the companion but does not link to its installation guide or summarize the
  interpolation and zoom-dependent aggregation choices.
- **A newcomer needs:** A clickable installation/reference route from the main
  site and one decision table for polygons, interpolated surfaces, and
  hierarchy-based resampling.
- **Proposed correction:** Link the companion README from Optional integrations
  and the gallery. Copy only the three-way choice and experimental-status
  sentence into the main site.

### UB-005 — The extension guide covers a hierarchy, not a general finite grid

- **Kind / priority:** DOC GAP / P1
- **Use cases:** UC-063, UC-064
- **Affected files:** `docs/src/extending.md`;
  `docs/src/internals/grid-contracts.md`; grid type docstrings
- **Evidence:** The guide opens by requiring a cell ID and an
  `AbstractHierarchicalGridSystem`, then gives thirteen system methods. That
  route is strong for hierarchical systems and links the conformance suites.
  It does not give the smaller contract for a standalone `AbstractGrid` or a
  task table for choosing a finite grid, hierarchical system, or quad-face
  system.
- **A newcomer needs:** Separate minimum method lists and conformance commands
  for the three extension shapes. The quad-face route must identify which face
  layout and projection hooks are required.
- **Proposed correction:** Add an extension-shape decision table before the
  worked hierarchy. Link each row to its exact contract and test suite.

### UB-006 — Optimized boundary hooks are discoverable and bounded

- **Kind / priority:** VERIFIED DOC COVERAGE
- **Use case:** UC-065
- **Evidence:** `docs/src/internals/boundary-engines.md` separates public
  results from internal engines, lists specialization hooks, says every engine
  returns exact results, and names fallback behavior. The system-capability page
  records order and memory bounds, seam filtering, exceptional hexagons, and
  A5's scan fallback.
- **Result:** No documentation or API gap found in this walkthrough.

### UB-007 — The custom partitioner contract is complete

- **Kind / priority:** VERIFIED DOC COVERAGE
- **Use case:** UC-066
- **Evidence:** `docs/src/api/partitioning.md` gives the
  `partitionlabels` signature, normalized-capacity input, original-row
  alignment, valid label range, common validation, empty-work behavior, and
  optional-backend loading boundary. It also says partitioning assigns work but
  does not execute it.
- **Result:** No documentation or API gap found in this walkthrough.

### UB-008 — Deprecations have no public migration index

- **Kind / priority:** DOC GAP / P2
- **Use case:** UC-067
- **Affected files:** `src/deprecated.jl`; `src/chunks.jl`; public API pages
- **Evidence:** `cellposition` forwards to `localindex`, and chunk
  `globalindices` forwards to `ownedindices`. Their docstrings explain the
  preserved meanings. No public page lists these migrations, so a reader must
  encounter a warning or search source to find the complete set.
- **A newcomer needs:** Old name, replacement, retained index space, and removal
  status in one short table.
- **Proposed correction:** Add a migration section to the landing or
  abstractions page. Include the noncontiguous-subset distinction for
  `cellposition`.

### UB-009 — The fast documentation check is not discoverable publicly

- **Kind / priority:** DOC GAP / P1
- **Use case:** UC-068
- **Affected files:** repository `README.md`; contributor guidance;
  `docs/make.jl`
- **Evidence:** `DGG_DOCS_FAST=true` now gives a strict, non-deploying build
  and `docs/audit/validation.md` records its exact exclusions. Neither the
  repository README nor a public contributor page gives that command. Running
  `docs/make.jl` without the flag keeps the graphics-heavy, warning-only,
  deployment-configured route.
- **A newcomer needs:** The fast command first, its skipped-example list, and a
  separate CI/full-site command with environment requirements.
- **Proposed correction:** Add a contributor validation section to the README
  or contributing guide and link the detailed validation record.

### UB-010 — Standard exchange formats need a documented scope decision

- **Kind / priority:** DOC GAP and FUTURE SCOPE DECISION / P1
- **Use case:** UC-069
- **Affected files:** store, geometry, regridding, and astronomy guides
- **Evidence:** The package writes its Zarr DGGS layouts. Geometry follows
  GeoInterface and raster regridding returns arrays, but the public surface has
  no GeoJSON, Shapefile, GeoTIFF, NetCDF, FITS, or standard MOC writer. The
  astronomy tutorial mentions FITS only as a possible source. A
  `MultiOrderCellSet` is not documented as a standards-compatible MOC.
- **A newcomer needs:** A format matrix: direct support, tested external
  package, required coordinate and identifier conversion, metadata retained,
  and unsupported formats.
- **Proposed correction:** Document tested handoffs without claiming that
  similar data structures implement a standard. Decide separately whether any
  direct writer belongs in this package before creating an API proposal.

### UB-011 — Physical-radius kernels are a composition with no task route

- **Kind / priority:** DOC GAP and FUTURE SCOPE DECISION / P1
- **Use case:** UC-070
- **Evidence:** The neighborhood API and stencil tutorial now state that
  `neighbors(grid, cell, k)` and `ring` can select topological neighborhoods,
  while `mapneighbors` supplies one ring only. The tutorial's gradient divides
  by explicit great-circle distance; its displayed Laplacian is an unscaled
  neighbor mean. There is no built-in Gaussian, physical Laplacian, or
  field-wide metric-radius sweep.
- **Result:** A caller can compose a cap or candidate selection with explicit
  centroid distances and user-defined weights. The documentation does not walk
  through that composition. A metric-radius operation is a separate scope
  decision; it is not the topological radius/ring sweep tracked as API-005.

### UB-012 — Adaptive mesh building blocks are not presented as a solver boundary

- **Kind / priority:** DOC GAP / P1
- **Use case:** UC-071
- **Affected files:** neighborhood, selection, boundary, and hierarchy pages
- **Evidence:** `member_neighbors` covers mixed-level adjacency, and the
  boundary guide says `grow` is single-level. The `grow` docstring says a
  `MultiOrderCellSet` has no method. No public page collects the absent
  finite-volume pieces: oriented faces, fluxes, restriction/prolongation,
  refinement criteria, time stepping, and a mixed-level halo.
- **A newcomer needs:** A boundary statement that lists usable geometry and
  adjacency primitives and separately lists solver operations the package does
  not supply.
- **Proposed correction:** Add a short “mesh algorithms supplied here” section
  to the hierarchy or neighborhood page. Do not present custom solver code as a
  package operation.

### UB-013 — Accelerator and execution boundaries are scattered

- **Kind / priority:** DOC GAP / P1
- **Use case:** UC-072
- **Affected files:** `docs/src/api/chunk-sweep.md`;
  `docs/src/api/partitioning.md`; `docs/src/index.md`
- **Evidence:** The docs clearly describe CPU threads, chunk assignment, and the
  executor's responsibility for retries, movement, and writes. They make no
  explicit support statement for GPU arrays, device kernels, automatic
  differentiation, or automatic cluster scheduling. An experimental Dagger
  runner exists only under `scripts/`.
- **A newcomer needs:** A support matrix that says which arrays and execution
  models are tested and labels the Dagger runner as experimental.
- **Proposed correction:** Add the matrix to the partitioning page. Keep
  assignment separate from execution.

### UB-014 — Remote and concurrent write limits need one explicit table

- **Kind / priority:** DOC GAP / P2
- **Use case:** UC-073
- **Affected files:** `docs/src/tutorials/store_io.jl`;
  `docs/src/api/subzone-layout.md`; Zarr extension docstrings
- **Evidence:** The tutorial says URL reads are supported, S3 needs AWSS3,
  writes are local or to an open group, and publication is an external upload.
  The subzone layout says whole columns can be written independently, while
  the implementation states that `SubzoneStore` holds no lock. The public
  route does not state the same-column race, atomicity, transaction, credential,
  append-coordinate, or retry boundary in one place.
- **A newcomer needs:** Read, create, reopen, append-column, concurrent-disjoint,
  conflicting-write, and remote-publication support in one table.
- **Proposed correction:** Add the table to the store API. State only the
  disjoint local-column behavior established by tests.

### UB-015 — Projected rasters and custom ellipsoids are implemented but hard to find

- **Kind / priority:** DOC GAP / P1
- **Use case:** UC-074
- **Affected files:** `docs/src/tutorials/choosing_a_grid.jl`;
  `docs/src/api/regridding.md`; `docs/src/index.md`;
  GlobalRegridding raster-space docstrings
- **Evidence:** `AuthalicSystem` accepts a custom `Geodesic` ellipsoid or a
  `Spherical(; radius)` manifold. GlobalRegridding accepts an explicit
  `Proj.Transformation` for projected raster coordinates. The chooser shows
  only the default WGS84 wrapper, and the main optional-integration table omits
  Proj because the extension lives in the workspace subpackage.
- **A newcomer needs:** One custom-ellipsoid example, one projected-raster
  example with `always_xy=true`, and a statement that DGGS geometry remains
  spherical rather than a general planar grid contract.
- **Proposed correction:** Add the two examples and link the GlobalRegridding
  coordinate contract from the regridding page.

### UB-016 — Scalability contracts are accurate but lack a planning worksheet

- **Kind / priority:** DOC GAP / P2
- **Use case:** UC-075
- **Affected files:** chooser, system-capability, boundary, chunk-sweep, and
  regridding pages
- **Evidence:** The docs expose cell counts, lazy boundary iterators,
  system-specific traversal costs, chunk ownership, lazy regridding budgets,
  and A5's scan fallback. These facts live on separate pages, and the chooser
  does not link to the system-capability cost page.
- **A newcomer needs:** A small checklist that estimates cell values, cell IDs,
  adjacency or halo storage, one active chunk and halo, regridding weights, and
  output storage before constructing a global level.
- **Proposed correction:** Link the capability page from the chooser and add a
  symbolic memory worksheet. Avoid performance numbers inferred from one run.

## Question-by-question discovery

This section records the answers found while following each inventory entry.
“Not stated” means the answer required source inspection or remains outside the
documented contract.

- **UC-058:** The tutorial identifies Copernicus elevations and downhill drops
  as metres and flow accumulation as square metres. Its custom downhill kernel
  skips `NaN`, leaves pits unrouted, and computes elevation difference rather
  than slope. The conventions used by Geomorphometry for roughness, slope,
  flats, and exceptional cells are not stated on the package site; those
  algorithms belong to Geomorphometry.
- **UC-059:** The public constructor accepts only `30` and `90`. Level 0 is a
  one-degree tile and level 1 is a sampled post. The module contract describes
  the latitude-band reductions, WGS84 coordinates, seam ownership, and polar
  limits. The site has no small workflow that puts those answers together or
  exercises a band edge and tile seam.
- **UC-060:** Level `l` has `nside = 2^l` and `12*nside^2` cells. Level grids
  and value storage use nested order. `reindex` converts an identifier; it does
  not move the value stored at that identifier. Celestial frame metadata is not
  represented, and this package has no FITS reader. `Intersects(cap)` selects
  cells whose footprints meet the cone; `Within(cap)` selects complete
  footprints inside it.
- **UC-061:** The plotting names and accepted grid, collection, and `Cells`
  inputs belong to `DiscreteGlobalGridsVisualization`. The companion README
  names CairoMakie, GLMakie, WGLMakie, and GeoMakie, requires one color per cell
  in collection order, omits `missing` and `NaN` colors, and describes seams,
  poles, projections, and outlines. The main site does not lead readers to that
  installation and limit information.
- **UC-062:** `dggsurface` interpolates between cell-center values and may use
  a height per cell; it does not change the stored field. `dggresample` chooses
  display resolution from the view and either samples or applies the supplied
  aggregation. Interactive use requires GLMakie or WGLMakie, and projected
  axes require GeoMakie. The companion labels all three plot APIs experimental;
  no numeric interactive-resolution guarantee is documented.
- **UC-063:** The public extension path does not give the minimum
  `AbstractGrid` method list, local-index rules, inherited queries, or fallback
  costs as one finite-grid contract. Those answers can be reconstructed from
  the interface source and conformance suite, so the use case is implemented
  but not newcomer-complete.
- **UC-064:** `HierarchicalLevelGrid` is the documented generic level-grid
  route. The revised trait pages state when sorted subtrees, congruent
  refinement, and direct location may be claimed and name their costs. The
  hierarchy conformance route is clear. Required quad-face layout and
  projection hooks are not collected into a public checklist.
- **UC-065:** The boundary-engine page names override hooks, order, memory,
  completeness, and fallback contracts. It also states seam filtering,
  exceptional-cell checks, and which helpers are internal. This use case is
  clear.
- **UC-066:** `partitionlabels` receives normalized positive capacities and
  must return one label in `1:nparts` per original row. Common code validates
  labels and restores task order. Empty work and excess partition behavior are
  documented, as is optional backend loading. Partitioning assigns work and
  does not execute it.
- **UC-067:** `cellposition` still means the collection-local position and
  forwards to `localindex`; `globalindex` is the complete-level position.
  `globalindices(::MapChunk/::ChunkCube)` forwards to `ownedindices`. These
  shims are actionable in their docstrings, but no public migration index
  collects them or the plotting-package moves.
- **UC-068:** The repository docs environment is Julia 1.13 with Documenter
  and Literate. Normal Literate conversion does not execute examples, but
  Documenter later executes generated `@example` blocks. The new fast mode is
  strict and does not deploy; it explicitly skips thirteen heavy pages. The
  public contributor route does not yet state its command, dependencies, or
  exclusions.
- **UC-069:** Core has no direct GeoJSON, Shapefile, GeoTIFF, NetCDF, FITS, or
  standard MOC writer. It writes its own Zarr DGGS layouts. GeoInterface
  geometry and regridded arrays are potential handoff values, but no external
  writer handoff is tested in the public guide. `MultiOrderCellSet` is not
  documented as standards-compatible MOC data, and the required CRS, frame,
  ID-order, and metadata conversions are therefore unspecified.
- **UC-070:** A topological ring is not a kilometre radius. `Intersects(cap)`
  selects footprints rather than centroids. The existing pieces allow a caller
  to select candidates, measure explicit spherical distances, and define area
  or distance weights. No built-in Gaussian or physical Laplacian exists, and
  the docs do not give the normalization and boundary rules for a composed
  physical kernel.
- **UC-071:** Public building blocks include cell geometry, one-level
  adjacency, and `member_neighbors` for mixed-level membership. No public
  mixed-level halo, oriented-face flux API, time stepper, refinement policy, or
  restriction/prolongation operator was found. Congruent refinement is true
  only for the systems named in the capability page. A complete adaptive PDE
  solver is outside the current public surface.
- **UC-072:** Documented execution is CPU threading plus logical chunk
  assignment. The executor remains responsible for movement, retries, and
  writes. No supported GPU array, device kernel, automatic-differentiation, or
  automatic cluster-scheduling contract was found. The repository Dagger
  runner is experimental script tooling.
- **UC-073:** URL and S3 reads are documented; S3 activation belongs to AWSS3.
  A URL destination is rejected, while a caller may write to a local path or
  already opened writable group and publish externally. The store does not
  manage credentials. `SubzoneStore` has no lock, and the docs establish only
  independent-column writing; they do not promise atomic conflicting writes,
  transactions, retries, or resumability.
- **UC-074:** `AuthalicSystem` accepts a `Geodesic` with semi-major axis and
  inverse flattening or a `Spherical` with radius. It maps geometry through an
  authalic sphere. GlobalRegridding accepts an explicit `Proj.Transformation`
  for projected raster coordinates. A planar grid is not part of the spherical
  DGGS contract. The main site does not collect the remaining unit and
  coordinate consequences or show either nondefault case.
- **UC-075:** Cell-count formulas, lazy boundary iterators, unknown-length
  cases, A5's full-level scan fallback, chunk ownership, and regridding storage
  budgets are documented. Complete grids and value vectors still scale with
  cell count, while the named boundary iterators keep depth-sized traversal
  state. No page combines data, identifiers, halo or adjacency, one active
  chunk, weight blocks, and output into a preflight memory estimate.

## Use-case coverage

| Use case | Result | Main evidence |
| --- | --- | --- |
| UC-058 | Partial | UB-001 |
| UC-059 | Partial | UB-002 |
| UC-060 | Partial | UB-003 |
| UC-061 | Partial | UB-004 |
| UC-062 | Partial | UB-004 |
| UC-063 | Partial | UB-005 |
| UC-064 | Clear for hierarchy; partial for quad-face discovery | UB-005 |
| UC-065 | Clear | UB-006 |
| UC-066 | Clear | UB-007 |
| UC-067 | Partial | UB-008 |
| UC-068 | Partial | UB-009 |
| UC-069 | Unsupported directly; external route undocumented | UB-010 |
| UC-070 | Composable from primitives; task route partial | UB-011 |
| UC-071 | Unsupported as a complete solver | UB-012 |
| UC-072 | Unsupported as a turnkey accelerator or scheduler route | UB-013 |
| UC-073 | Partial | UB-014 |
| UC-074 | Partial | UB-015 |
| UC-075 | Clear contracts; planning route partial | UB-016 |

## Evidence limits

The walkthrough did not download a DEM, open FITS data, render graphics, load a
native partitioner, contact a remote store, or run a distributed scheduler.
Those checks would test external systems rather than the documentation
questions above.

## Lightweight local checks

A Julia 1.13 scratch module loaded the repository environment through a private
`jld` daemon. One bounded test set made 41 assertions and printed
`boundary-walkthrough-probes-ok`. It checked:

- both public Copernicus DEM constructors, their two levels and 64,800 roots;
- point location and tile parents on both sides of the 50° band transition and
  the antimeridian, plus rejection of an unsupported nominal resolution;
- the complete level-2 HEALPix ring-to-nested permutation, value reordering,
  and identifier round trips;
- `cellposition` against `localindex` semantics on a noncontiguous cell vector;
- a minimal custom `partitionlabels` implementation and complete assignment;
- point location through `AuthalicSystem` on a nondefault ellipsoid; and
- documented HEALPix and IGeo7 cell-count formulas through level 5.

These probes establish that the small public primitives used as evidence still
work. They do not turn an absent task guide, external integration, or undecided
package boundary into verified support.

## Documentation-only source check

The 40 modified Julia files under `src/`, `ext/`, and `lib/` were parsed from
both the working tree and `HEAD`. After replacing only `Core.@doc` payloads and
removing line nodes, every syntax tree was equal. The result was
`non-doc-ast-differences=String[]`. This includes the GlobalRegridding files
whose reference spelling changed during strict documentation validation.
