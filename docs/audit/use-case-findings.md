# Documentation backlog from the 75 use-case walkthroughs

## Scope

This ledger consolidates the second audit stage. The source reports are:

- [UC-001–024](walkthrough-core.md)
- [UC-025–057](walkthrough-workflows.md)
- [UC-058–075](walkthrough-boundaries.md)

The walkthrough reports remain the question-level evidence. This ledger groups
questions only when one documentation change can answer them together. It does
not replace the reports. It does not mark a missing API as complete.

All entries remain open unless their status says otherwise. The A5 area entry
records the one narrow accuracy correction made after all source reports were
saved. No other entry was fixed during consolidation.

Priority follows the initial ledger:

- **P0:** A false contract can cause silent misuse or incorrect values.
- **P1:** A main workflow or support boundary is missing.
- **P2:** Navigation, planning, or an edge contract is incomplete.

## Deduplicated backlog

### UCDOC-001 — Publish reproducible setup and validation routes

- **Priority:** P1
- **Source findings:** UCOR-001, UB-009
- **Use cases:** UC-001, UC-068
- **Questions to answer:** Which package revision does the repository install command select? Do users need workspace projects? Which command runs the strict local documentation check? Which 13 pages keep their heavy examples disabled in fast mode?
- **Evidence:** `docs/src/index.md:74-88` and `README.md:24-31` use an unpinned repository URL. `Project.toml:6-7` shows that workspace projects are development inputs. `docs/audit/validation.md` records `DGG_DOCS_FAST=true`, but public contributor guidance does not.
- **Correction:** State the release or revision policy. State that normal package use does not need workspace projects. Publish the fast command, its exclusions, and the separate full-site requirements.
- **Related initial work:** DOC-018, DOC-020, DOC-021
- **Status:** Open.

### UCDOC-002 — Add a storage and scalability worksheet

- **Priority:** P2
- **Source findings:** UCOR-002, UB-016
- **Use cases:** UC-003, UC-075
- **Questions to answer:** How does a cell count become bytes for values, IDs, extra dimensions, adjacency, halos, weights, and output? Which objects materialize? Which objects stay lazy?
- **Evidence:** The chooser exposes `ncells` and physical resolution. The system, boundary, chunk, and regridding pages describe costs separately. No page combines them into a preflight estimate.
- **Correction:** Add a symbolic worksheet. Use element type, dimension lengths, active chunk size, halo width, and weight policy as inputs. Link system-specific traversal limits.
- **Related initial work:** DOC-002, DOC-003, DOC-016
- **Status:** Open.

### UCDOC-003 — Define raster sizing and coordinate inputs

- **Priority:** P1
- **Source findings:** UCOR-003, UB-015
- **Use cases:** UC-004, UC-074
- **Questions to answer:** Which X/Y lookup sampling forms work? Which units and coordinate frames apply? How does the fixed sample stride affect reproducibility? How does latitude affect the median? How do users supply a projected raster or a custom ellipsoid?
- **Evidence:** `src/sizing.jl:54-105` defines the deterministic sample and candidate-only `over` behavior. `test/interface/sizing.jl:42-69` demonstrates latitude dependence. GlobalRegridding accepts an explicit Proj transformation. The chooser shows only WGS84.
- **Correction:** Add accepted-axis and coordinate tables. Add a latitude-dependent size example, a custom-ellipsoid example, and a projected-raster example with `always_xy=true`.
- **Related initial work:** DOC-013, DOC-016, API-004
- **Status:** Open.

### UCDOC-004 — Qualify A5 equal-area claims by measurement surface

- **Priority:** P0
- **Source findings:** UCOR-012
- **Use case:** UC-002
- **Questions to answer:** Is A5 equal-area in ellipsoidal area, unit-sphere steradians, or both? Which quantity does `cell_area` return?
- **Evidence:** The gallery and capability table called A5 equal-area. The chooser called it nearly equal. `src/systems/A5/system.jl:4-14` says A5 is equal-area on its ellipsoid, while unit-sphere `cell_area` varies about 1%. A level-2 execution measured a maximum-to-minimum ratio of `1.0128715`.
- **Correction:** State “equal-area on its ellipsoid; unit-sphere steradians vary about 1%” in the two overview rows. Keep the longer system docstring as the authority.
- **Related initial work:** DOC-002
- **Status:** Corrected after the source reports were saved. This narrow correction fixes an inaccurate overview label. It does not close other second-stage work.

### UCDOC-005 — Route readers to each system's point-location tie rule

- **Priority:** P2
- **Source finding:** UCOR-004
- **Use case:** UC-006
- **Questions to answer:** Which incident cell owns an exact edge, corner, seam, or pole point for each system? Which rules can vary between floating-point platforms?
- **Evidence:** The generic `cellat` docstring promises that each system documents its rule. The system docstrings contain rules, but the user route does not link or summarize them.
- **Correction:** Add one tie-rule table or direct links beside `cellat`. State that a valid tie can select a different incident cell on another platform.
- **Related initial work:** DOC-002, DOC-007
- **Status:** Open.

### UCDOC-006 — Explain identifier context and value reordering

- **Priority:** P1
- **Source findings:** UCOR-005, UB-003
- **Use cases:** UC-007, UC-060, UC-069
- **Questions to answer:** Which typed IDs identify their system without a grid? When does `LevelIndex` still need system context? How does a reader permute a ring-order HEALPix value vector into nested cell order? Which celestial frame does a cell axis represent?
- **Evidence:** `LevelIndex` is shared by HEALPix, S2, and ISEA4R. The HEALPix tutorial tells readers to convert ring-order input but converts no value vector. `reindex` changes IDs; it does not reorder associated values.
- **Correction:** Add a system-context warning and a small ring-to-nested permutation. State that `Cells` does not carry celestial-frame metadata by itself.
- **Related initial work:** DOC-002, DOC-006
- **Status:** Open.

### UCDOC-007 — Document tree lifetime and node index space

- **Priority:** P2
- **Source finding:** UCOR-006
- **Use case:** UC-010
- **Questions to answer:** Does `query` accept or reuse a prepared tree? What object owns the tree? Are node indices local to the grid, local to the node, or complete-level positions?
- **Evidence:** The tree page names `treeify`. The optional multi-order example explains the synthetic root. `src/engine/cursor.jl:107-121,449-461` shows that node indices address the original grid. `query` has no prepared-tree argument.
- **Correction:** Add a short tree lifecycle contract and one index-space sentence. Do not imply a cache or reuse API that does not exist.
- **Related initial work:** DOC-003
- **Status:** Open.

### UCDOC-008 — Add a difficult global-geometry support matrix

- **Priority:** P1
- **Source finding:** UCOR-007
- **Use cases:** UC-012, UC-013
- **Questions to answer:** Which predicate results are impossible by dimensionality? How do antimeridian, pole, hole, multipart, empty, and large-polygon inputs behave? Which winding and validity rules apply? Must projected input be transformed first?
- **Evidence:** `docs/src/api/selecting-cells.md:67-96` defines direction and target types. Cross-system tests cover seams, holes, multipart shapes, and wide polygons. The public path has no shape and validity matrix.
- **Correction:** Add a matrix with supported input, coordinate frame, boundary rule, empty behavior, and tested systems. Give a projected-input preparation step.
- **Related initial work:** DOC-004
- **Status:** Open.

### UCDOC-009 — Complete the collection and cube construction contract

- **Priority:** P1
- **Source finding:** UCOR-008
- **Use cases:** UC-014, UC-015
- **Questions to answer:** Can a region be empty? Can `Cells` occupy any dimension? What happens when axis and data lengths differ? Which names, metadata, missing values, bands, and time dimensions survive construction and slicing?
- **Evidence:** The implementation accepts correctly typed empty IDs and non-leading cell dimensions. Tests raise `DimensionMismatch` for an axis-length error. The public examples use a leading or sole cell dimension.
- **Correction:** Add one non-leading cell-time cube. Include empty construction, dimension validation, and metadata preservation in the contract.
- **Related initial work:** DOC-006, DOC-015
- **Status:** Open.

### UCDOC-010 — State selector cardinality, errors, and laziness

- **Priority:** P1
- **Source finding:** UCOR-009
- **Use cases:** UC-007, UC-016
- **Questions to answer:** Which exact selectors raise when a cell is absent? Which spatial selectors return an empty result? Does an empty mask retain a cell lookup? Does each selection allocate, return a view, or preserve lazy storage?
- **Evidence:** Dimensional tests show `At` and point `Contains` errors. They show empty geometry and boolean results. The store tutorial shows a selected in-memory lookup but does not state the general rule.
- **Correction:** Add one selector outcome table. Include expected cardinality, absent behavior, lookup type, and storage behavior.
- **Related initial work:** DOC-001, DOC-006, DOC-019
- **Status:** Open.

### UCDOC-011 — Document region algebra as a task

- **Priority:** P1
- **Source finding:** UCOR-010
- **Use case:** UC-017
- **Questions to answer:** Which systems and levels can combine? How do `union`, `intersect`, `vcat`, and `issubset` differ? Which result order and duplicate rules apply? Do operations combine field values?
- **Evidence:** The selection page renders `expand` and `compact` only. `src/engine/region_algebra.jl` and collection tests define the other operations and compatibility errors.
- **Correction:** Add a region-algebra table and one overlapping-set example. State that these verbs combine membership, not values.
- **Related initial work:** DOC-006
- **Status:** Open.

### UCDOC-012 — Complete the zonal-statistics decision path

- **Priority:** P1
- **Source finding:** UCOR-011
- **Use case:** UC-018
- **Questions to answer:** Does the statistic use selected cells or overlap fractions? When must values use cell-area weights? What denominator applies with missing coastal data? How do overlapping zones count samples?
- **Evidence:** The zonal tutorial computes selection means and warns that they are not exact polygon statistics. It uses `skipmissing` without defining coverage semantics. It gives no overlap-weighted recipe.
- **Correction:** Add a method table and two small independent checks. Keep center-in-zone selection linked to API-001.
- **Related initial work:** DOC-011, DOC-013, DOC-016, API-001
- **Status:** Open.

### UCDOC-013 — Define terrain measurement and boundary conventions

- **Priority:** P1
- **Source finding:** UB-001
- **Use case:** UC-058
- **Questions to answer:** Are DEM inputs posts or pixel means? Which units and kernels define roughness and slope? How do flats, pits, missing cells, region edges, and pentagons behave? Which package owns each result?
- **Evidence:** The hydrology tutorial identifies several units and labels its custom downhill result. It does not state the full Geomorphometry input and edge contract.
- **Correction:** Add a result table beside the first terrain calculation. Link the owning package for composed algorithms.
- **Related initial work:** DOC-003, DOC-014
- **Status:** Open.

### UCDOC-014 — Add a Copernicus DEM task route

- **Priority:** P1
- **Source finding:** UB-002
- **Use case:** UC-059
- **Questions to answer:** How does a user construct GLO-30 or GLO-90 geometry, locate a point, find a tile parent, cross a latitude band, and handle seams and polar limits? Does construction download elevations?
- **Evidence:** The module docstring answers the lattice contract. The chooser only summarizes it. The hydrology tutorial downloads a raster and regrids it to IGeo7 instead of using `CopernicusDEMSystem`.
- **Correction:** Add a small local specialized-system example. Link the module contract and label repository runners as experimental.
- **Related initial work:** DOC-002, DOC-009
- **Status:** Open.

### UCDOC-015 — Link plotting installation and method limits

- **Priority:** P2
- **Source finding:** UB-004
- **Use cases:** UC-061, UC-062
- **Questions to answer:** How does a user install the companion package? When should the user choose polygons, interpolated surfaces, or hierarchy resampling? Which seam, pole, missing-color, aggregation, and experimental limits apply?
- **Evidence:** The companion README answers these questions. The main site names the package but does not link its installation and method guidance.
- **Correction:** Add a direct companion link and a short three-method decision table.
- **Related initial work:** DOC-018
- **Status:** Open.

### UCDOC-016 — Separate the three extension shapes

- **Priority:** P1
- **Source finding:** UB-005
- **Use cases:** UC-063, UC-064
- **Questions to answer:** When should an implementor add a finite grid, a hierarchical system, or a quad-face system? Which minimum methods and conformance checks apply to each?
- **Evidence:** The extension guide gives a strong hierarchical-system example. It does not give the smaller standalone-grid contract or the quad-face hook list.
- **Correction:** Add an extension-shape decision table. Link each shape to its exact methods and tests.
- **Related initial work:** DOC-003, DOC-008, DOC-009
- **Status:** Open.

### UCDOC-017 — Publish a migration index

- **Priority:** P2
- **Source finding:** UB-008
- **Use case:** UC-067
- **Questions to answer:** Which old names remain, what replaces each name, which index space stays unchanged, and when can a name be removed?
- **Evidence:** `src/deprecated.jl` maps `cellposition` to `localindex`. Chunk code maps `globalindices` to `ownedindices`. No public page lists the migrations.
- **Correction:** Add a table with old name, replacement, semantic note, and removal status.
- **Related initial work:** None
- **Status:** Open.

### UCDOC-018 — Document tested exchange-format handoffs

- **Priority:** P1
- **Source finding:** UB-010
- **Use case:** UC-069
- **Questions to answer:** Which standard formats have direct writers? Which external package handoffs are tested? Which coordinates, IDs, ordering, and metadata need conversion? Is a `MultiOrderCellSet` a standards-compatible MOC?
- **Evidence:** The package writes its own Zarr layouts. Geometry follows GeoInterface, and regridding returns dimensional arrays. No public route documents tested GeoJSON, Shapefile, GeoTIFF, NetCDF, FITS, or MOC handoffs.
- **Correction:** Add a format matrix. State “unsupported” where no direct writer or tested handoff exists. Do not infer standard compatibility from a similar structure.
- **Related initial work:** DOC-012, DOC-016, DOC-019
- **Status:** Open. Direct writers remain a future package-scope question, not a confirmed API obligation.

### UCDOC-019 — State the adaptive-solver boundary

- **Priority:** P1
- **Source finding:** UB-012
- **Use case:** UC-071
- **Questions to answer:** Which mixed-level geometry and adjacency pieces exist? Which oriented faces, fluxes, transfer operators, refinement rules, time stepping, and mixed-level halos are absent?
- **Evidence:** `member_neighbors` supplies mixed-level adjacency. `grow` is single-level. No public page separates those primitives from a complete adaptive solver.
- **Correction:** Add a supplied-primitives and absent-solver-pieces table. Do not present custom solver code as package support.
- **Related initial work:** DOC-005, DOC-014
- **Status:** Open. A full solver remains a future scope question.

### UCDOC-020 — Publish accelerator and execution support boundaries

- **Priority:** P1
- **Source finding:** UB-013
- **Use case:** UC-072
- **Questions to answer:** Which array types, CPU threads, processes, GPU kernels, automatic differentiation, and schedulers are tested? Does partitioning execute work? What status does the Dagger script have?
- **Evidence:** The public pages explain CPU threads and logical partitioning. They do not state GPU, AD, or automatic-cluster support. The Dagger runner exists only under `scripts/`.
- **Correction:** Add an execution support matrix. Label the script as experimental repository tooling.
- **Related initial work:** DOC-017, DOC-018
- **Status:** Open. Accelerator support remains a future scope question.

### UCDOC-021 — Consolidate remote and concurrent write limits

- **Priority:** P2
- **Source finding:** UB-014
- **Use case:** UC-073
- **Questions to answer:** Which remote reads work? Where can writes go? Are disjoint columns safe? What happens on same-column races? Which atomicity, retry, credential, transaction, and append guarantees exist?
- **Evidence:** The store tutorial documents URL reads and external upload. The subzone docs allow independent whole-column writes. `SubzoneStore` holds no lock. No page gives the complete boundary.
- **Correction:** Add one support table. State only the tested disjoint local-column behavior.
- **Related initial work:** DOC-019
- **Status:** Open. Remote transactions and general append remain future scope questions.

### UCDOC-022 — State when a DGGS source needs `from`

- **Priority:** P1
- **Source finding:** UW-001
- **Use case:** UC-037
- **Questions to answer:** Which dimensional sources infer their spatial space? Why does a DGGS `DimArray` over `Cells` still need `from`?
- **Evidence:** A walkthrough run omitted `from` and received `ArgumentError`. The same call passed with an explicit source grid. `lib/GlobalRegridding/src/api.jl:290-310` owns the rejection.
- **Correction:** Distinguish inferred raster X/Y axes from a DGGS value array. Add one DGGS-source example with `from`.
- **Related initial work:** DOC-016
- **Status:** Open.

### UCDOC-023 — Specify callback buffer ownership and lifetime

- **Priority:** P1
- **Source finding:** UW-002
- **Use cases:** UC-034, UC-053
- **Questions to answer:** May a callback retain or mutate neighbor sequences, requested rings, or chunk data? Does the contract change with threads?
- **Evidence:** Sweep and chunk pages define contents and index spaces. They do not define borrowing, ownership, retention, or mutation. Current allocation is an implementation detail.
- **Correction:** Mark each callback argument as borrowed or owned. State supported mutation and retention across calls and threads.
- **Related initial work:** DOC-003, DOC-014, DOC-015
- **Status:** Open.

### UCDOC-024 — Show how owned values map into a grown region

- **Priority:** P1
- **Source finding:** UW-003
- **Use case:** UC-033
- **Questions to answer:** How does a user add halo values, map original cells into the sorted grown region, compute owned outputs, and restore original order?
- **Evidence:** The boundary page explains grown membership. `grow` returns a new sorted `CellVector` without an ownership map.
- **Correction:** Add a small identity-based remapping example. Keep owned outputs separate from halo inputs.
- **Related initial work:** DOC-006, DOC-014
- **Status:** Open.

### UCDOC-025 — Define cost-distance missing and unreachable behavior

- **Priority:** P2
- **Source finding:** UW-004
- **Use case:** UC-035
- **Questions to answer:** Which cost values are valid? How should missing samples become barriers? What value marks an unreachable cell?
- **Evidence:** The stencil tutorial models barriers with `Inf`. Its example initializes all distances to `Inf`, so disconnected cells remain `Inf`. It does not define missing input.
- **Correction:** Label the function as example code. Define valid costs, missing conversion, and unreachable output.
- **Related initial work:** DOC-014
- **Status:** Open.

### UCDOC-026 — Explain mixed-level regridding output granularity

- **Priority:** P1
- **Source finding:** UW-005
- **Use case:** UC-041
- **Questions to answer:** Does a mixed-level target return one value per stored member or per expanded leaf? Which level and order apply? What happens for an empty target?
- **Evidence:** The regridding API lists `MultiOrderCellSet` as a target. `src/regridding.jl:289-294` converts it through `CellVector(set)` and `PartialGrid`, which produces a one-level leaf space.
- **Correction:** State the expansion level, order, and value meaning. Add a mixed-level example and validate empty-target behavior.
- **Related initial work:** DOC-005, DOC-016
- **Status:** Open.

### UCDOC-027 — Forbid or reject chunked source and destination aliasing

- **Priority:** P0
- **Source finding:** UW-006
- **Use case:** UC-052
- **Questions to answer:** May `dest` alias the swept data or a requested value field? If not, does the API reject it before writes begin?
- **Evidence:** A walkthrough run used the source parent as `dest`. Later chunk loads observed earlier writes, and the result silently differed from the eager one-ring result. A separate destination matched the eager result. `src/chunks.jl:686-707` does not state an alias rule.
- **Correction:** State that aliasing is unsupported unless the implementation adds snapshot semantics. Add a rejection or a focused non-alias validation.
- **Related initial work:** DOC-003, DOC-014
- **Status:** Open. This is a high-severity contract gap with an observed runtime mismatch. The walkthrough did not establish a promised in-place API.

### UCDOC-028 — Add a physical-radius neighborhood task route

- **Priority:** P1
- **Source finding:** UB-011
- **Use case:** UC-070
- **Questions to answer:** How does a user select candidates by physical radius, compute centroid distances, apply metric weights, and sweep a whole field? Which approximation and boundary rules apply?
- **Evidence:** The public neighborhood API supplies topological disks and rings. `mapneighbors` supplies one-ring callbacks. The stencil tutorial divides one gradient by an explicit great-circle distance, but its displayed Laplacian is an unscaled neighbor mean. No task route composes a metric-radius kernel.
- **Correction:** Add one explicit composition example and state its approximation limits. Keep a built-in metric kernel as a future package-scope decision.
- **Related initial work:** DOC-011, DOC-014, API-001
- **Status:** Open. This is separate from API-005, whose acceptance contract concerns topological radius and ring sweeps.

## Source-finding map

Every source finding from the three walkthrough reports appears here.

| Source finding | Consolidated result |
| --- | --- |
| UCOR-001 | UCDOC-001 |
| UCOR-002 | UCDOC-002 |
| UCOR-003 | UCDOC-003 |
| UCOR-004 | UCDOC-005 |
| UCOR-005 | UCDOC-006 |
| UCOR-006 | UCDOC-007 |
| UCOR-007 | UCDOC-008 |
| UCOR-008 | UCDOC-009 |
| UCOR-009 | UCDOC-010 |
| UCOR-010 | UCDOC-011 |
| UCOR-011 | UCDOC-012 |
| UCOR-012 | UCDOC-004 |
| UW-001 | UCDOC-022 |
| UW-002 | UCDOC-023 |
| UW-003 | UCDOC-024 |
| UW-004 | UCDOC-025 |
| UW-005 | UCDOC-026 |
| UW-006 | UCDOC-027 |
| UB-001 | UCDOC-013 |
| UB-002 | UCDOC-014 |
| UB-003 | UCDOC-006 |
| UB-004 | UCDOC-015 |
| UB-005 | UCDOC-016 |
| UB-006 | No backlog. The walkthrough verified the boundary-hook documentation. |
| UB-007 | No backlog. The walkthrough verified the custom partitioner contract. |
| UB-008 | UCDOC-017 |
| UB-009 | UCDOC-001 |
| UB-010 | UCDOC-018; future API-scope candidate |
| UB-011 | UCDOC-028; future API-scope candidate distinct from API-005 |
| UB-012 | UCDOC-019; future API-scope candidate |
| UB-013 | UCDOC-020; future API-scope candidate |
| UB-014 | UCDOC-021; future API-scope candidate |
| UB-015 | UCDOC-003 |
| UB-016 | UCDOC-002 |

## Status summary

- 28 deduplicated documentation entries were created.
- 27 entries remain open.
- UCDOC-004 received the one permitted narrow accuracy correction after the source reports were saved.
- Existing API-001 through API-005 remain open.
- The walkthroughs confirmed no new API obligation.
- Exchange writers, metric-radius kernels, adaptive solvers, accelerators, and remote transactions remain scope questions.
