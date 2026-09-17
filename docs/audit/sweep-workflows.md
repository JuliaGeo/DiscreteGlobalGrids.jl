# Workflow documentation sweep

This report records findings only. It does not change public documentation,
docstrings, or APIs. Local source and tests are the authority for behavior.

## Scope

The sweep covered neighborhood queries and sweeps, stencils, cell fields, stored
data, chunk execution, partitioning, regridding, optional extensions, tutorials,
and documentation navigation. A separate sweep covers grid systems, identifiers,
geometry, hierarchy, selection, and core collection types.

Severity has this meaning:

- **Critical:** The documentation states behavior that the API does not provide.
- **High:** A common workflow is blocked or can produce a wrong interpretation.
- **Medium:** A common task requires source reading or avoidable trial and error.
- **Low:** The issue affects clarity but does not usually block a task.

Each finding is either a **DOC GAP** or an **API GAP**. An API gap needs a later
design decision. This report does not prescribe a temporary workaround.

## Findings

### SW-001: Neighborhood sweeps cannot express a radius-k convolution

- **Kind:** API GAP
- **Severity:** High
- **Affected surface:** `mapneighbors`, `foreachneighbors`, `mapneighbors!`
- **Evidence:** Both sweep loops call `neighbors(..., 1)` directly at
  `src/engine/neighborhood.jl:267-279` and `src/engine/neighborhood.jl:284-297`.
  Their public signatures have no radius or neighborhood selector at
  `src/engine/neighborhood.jl:490-615`. The chunk runner delegates to the same
  one-ring operation at `src/chunks.jl:737-760`. In contrast, direct queries
  support `neighbors(grid, cell, k)` and `ring(grid, cell, k)` at
  `src/interface/grid.jl:316-464`.
- **A newcomer needs:** One global sweep whose callback receives either the disk
  through distance `k` or the exact ring at distance `k`. `neighbors(..., k)` is
  the disk, excludes the center, contains no duplicate cells, and concatenates
  exact rings from 1 through `k`. `ring(..., k)` is one exact shell;
  `ring(..., 0)` is the center. The callback should continue to receive the center
  separately.
- **Current API:** Only one-ring sweeps exist. Direct radius queries do not provide
  the global convolution operation or its chunked form.
- **Proposed correction:** Design a neighborhood selector or radius contract for
  all three sweep forms. Acceptance criteria should cover disk and exact-ring
  semantics, center inclusion, canonical ring-block order, no duplicates, subset
  clipping, eager and chunked equality, and rejection of a chunk plan whose halo
  is narrower than the requested neighborhood.

### SW-002: Halo width is documented as if it changes the sweep radius

- **Kind:** DOC GAP
- **Severity:** Critical
- **Affected surface:** `chunkplan`, `mapneighbors!`, stencil and out-of-core
  tutorials
- **Evidence:** `chunkplan` says that `halo = n` is what an n-ring stencil needs at
  `src/chunks.jl:152-167`. The chunk guide says any stencil within the plan width
  computes exactly at `docs/src/api/chunk-sweep.md:77-88`. The out-of-core tutorial
  repeats the rule at `docs/src/tutorials/out_of_core.jl:78-83`. However,
  `mapneighbors!` always calls one-ring `mapneighbors` at
  `src/chunks.jl:737-760`, and the eager sweep hard-codes radius 1 at
  `src/engine/neighborhood.jl:267-297`.
- **A newcomer needs:** A clear separation between data availability and callback
  reach. A wider halo supplies context to a custom `foreachchunk` callback. It
  does not change the neighbors passed by `mapneighbors!`.
- **Current documentation:** A reader can reasonably set `halo = 2` and expect a
  second-order convolution. They still receive one-ring values.
- **Proposed correction:** State the current one-ring limit beside every
  `mapneighbors` and `mapneighbors!` example. Describe a wider halo only as
  context available inside `foreachchunk`. In the stencil tutorial, state that
  repeated one-ring averaging is diffusion. It is not one radius-2 convolution:
  intermediate cells contribute with path multiplicity, and the center is
  introduced again on each pass. Link the future radius-k capability tracked in
  SW-001 without presenting a workaround as the supported workflow.

### SW-003: The three dimensional-array pass modes have no task-level home

- **Kind:** DOC GAP
- **Severity:** Medium
- **Affected surface:** `Neighbors`, `Values`, `NeighborSlices`, dimensional-array
  `mapneighbors`
- **Evidence:** The main docstring defines all three callback and output-shape
  contracts at `src/dimensionaldata.jl:564-654`. Behavioral tests cover
  multidimensional `Values` and `NeighborSlices` at
  `test/systems/crosssystem/mapneighbors.jl:187-248`. The neighborhood page lists
  only `Values` and `Neighbors` at `docs/src/api/neighbors.md:75-92`.
  `NeighborSlices` appears instead under spatial selection at
  `docs/src/api/selecting-cells.md:88-100`.
- **A newcomer needs:** One comparison beside `mapneighbors`: callback arguments,
  output dimensions, and the reason to choose each pass mode. A time-by-cells
  example should distinguish processing each time slice with `Values()` from
  giving one whole time series per cell with `NeighborSlices()`.
- **Current documentation:** The complete contract exists in a long docstring,
  but navigation places one mode on an unrelated page.
- **Proposed correction:** Put a compact three-row pass-mode table on the
  neighborhood page and render all three mode docstrings there. Keep the
  selection page focused on indexing and selectors.

### SW-004: The regridding reference omits the callable API and execution controls

- **Kind:** DOC GAP
- **Severity:** High
- **Affected surface:** `regrid`, `regrid!`, `plan_regrid`, `DGGSpace`, execution
  and storage policies
- **Evidence:** The method guide explains interpolation and missing-data behavior,
  but its only `@docs` block renders `Weighted`, `Extensive`, and `Conservative`
  at `docs/src/api/regridding-methods.md:67-85`. The public callable contracts are
  defined in `lib/GlobalRegridding/src/api.jl:11-242`, the execution storage
  types are defined in `lib/GlobalRegridding/src/plans.jl:110-444`, and the
  package-specific space is defined at `src/regridding.jl:11-20`.
  `docs/make.jl:26-30` registers
  the owner module so these re-exported docstrings can be rendered, but no API
  page renders them. Tutorials demonstrate plan reuse at
  `docs/src/tutorials/regridding.jl:112-121` and
  `docs/src/tutorials/between_grids.jl:154-164`.
- **A newcomer needs:** The full signatures, accepted source and destination
  spaces, return shapes, plan reuse, lazy execution, memory budget, storage
  policy, and errors. Method choice alone does not specify how to run a large job.
- **Current documentation:** Readers must inspect the subpackage source or infer
  keywords from tutorials. The README calls the grid-interface page the full
  reference at `README.md:100-112`, but that page does not cover regridding.
- **Proposed correction:** Add a regridding API section or page that renders the
  public callable and policy docstrings. Lead with the three normal workflows:
  one-shot, reusable plan, and in-place execution. Link to the existing method
  guide for interpolation details.

### SW-005: Partitioning examples are conceptual but look runnable

- **Kind:** DOC GAP
- **Severity:** Medium
- **Affected surface:** Regridding partition workflow
- **Evidence:** The first workflow imports `GlobalRegridding` directly and uses an
  undefined `regridplan` at `docs/src/api/partitioning.md:89-105`. The next uses
  undefined `dag`, `todochunks`, `tiles`, weights, and worker variables at
  `docs/src/api/partitioning.md:107-136`. These are plain `julia` fences, so the
  docs build does not execute them. `docs/make.jl:26-30` also states that
  `GlobalRegridding` is not a direct dependency of the docs environment.
- **A newcomer needs:** Either a minimal example they can run or an explicit
  schematic with every placeholder and its required type explained.
- **Current documentation:** The snippets use normal code-fence styling but cannot
  be copied into the documented environment.
- **Proposed correction:** Convert the smallest workflow to an executable
  `@example` built from a real public plan. Label the application-specific example
  as pseudocode and define each placeholder. Use the re-exported DGG surface unless
  the docs environment intentionally adds `GlobalRegridding` as a dependency.

### SW-006: Optional integrations lack one discovery point

- **Kind:** DOC GAP
- **Severity:** Medium
- **Affected surface:** Makie, Zarr, Metis, KaHyPar, and Scotch extensions
- **Evidence:** The package declares five extensions at `Project.toml:23-38`.
  Store I/O and partition backends have activation guidance at
  `docs/src/api/store-io.md:7-13` and `docs/src/api/partitioning.md:138-153`.
  Loading Makie also enables point and polygon conversion for grids, cell
  collections, mixed-level sets, coverages, and subtree iterators at
  `ext/MakieExt/cellsets.jl:1-68`; tests confirm these conversions at
  `test/plotting/runtests.jl:19-84`. No user documentation names that extension
  or distinguishes it from the visualization companion package.
- **A newcomer needs:** A short table that maps a capability to the package that
  activates it and to its reference page. For plotting, readers need to know what
  plain Makie conversion provides and what the visualization companion adds.
- **Current documentation:** Activation details are scattered. The core Makie
  integration is discoverable only from source and tests.
- **Proposed correction:** Add an optional-integrations table to installation or
  reference navigation. Include activation spelling and scope for all five
  extensions. Link plotting workflows to the companion package without hiding
  the base Makie conversions.

### SW-007: The store tutorial overstates URL support and understates encodings

- **Kind:** DOC GAP
- **Severity:** Medium
- **Affected surface:** `dggread`, store encodings, store tutorial
- **Evidence:** The tutorial says `dggread` opens public `gs://`, `s3://`, and
  `https://` stores at `docs/src/tutorials/store_io.jl:152-159`. The extension
  contract states that `s3://` needs `using AWSS3` at
  `ext/DiscreteGlobalGridsZarrExt/read.jl:74-78`. The tutorial introduction says
  that a store has two cell-ID encodings at
  `docs/src/tutorials/store_io.jl:1-8`, while the public reference lists dense,
  ranges, and implicit encodings at `docs/src/api/store-io.md:114-126`.
- **A newcomer needs:** The package needed for each remote URL scheme and the
  eligibility of each encoding. The tutorial's regional axis can demonstrate two
  encodings, while implicit encoding applies to a complete level.
- **Current documentation:** An S3 reader misses a required activation step, and
  the opening sentence makes a tutorial-local comparison sound exhaustive.
- **Proposed correction:** Add the AWSS3 requirement beside `s3://`. Change the
  introduction to “the two encodings applicable to this regional axis,” then link
  to implicit encoding for complete levels.

### SW-008: Workflow contracts are repeated in long docstrings and prose

- **Kind:** DOC GAP
- **Severity:** Medium
- **Affected surface:** Neighborhood sweeps, field requests, chunk sweeps, and
  store I/O
- **Evidence:** The `mapneighbors` docstring mixes callback contracts with cache,
  task scheduling, and chunk-number translation at
  `src/engine/neighborhood.jl:490-534`. The dimensional-array version repeats
  chunk-routing details at `src/dimensionaldata.jl:615-654`. The generic Zarr
  stubs carry large user contracts at `src/io/api.jl:26-106`, while the extension
  repeats and expands them at `ext/DiscreteGlobalGridsZarrExt/read.jl:64-122` and
  on the store API page at `docs/src/api/store-io.md:7-47`.
- **A newcomer needs:** A concise docstring with signatures, inputs, callback or
  return value, critical semantics, and important errors. Extended scheduling,
  caching, and storage rationale belongs on the workflow or architecture page.
- **Current documentation:** Repetition makes the stable contract harder to scan
  and creates several places that must stay synchronized.
- **Proposed correction:** Make the generic public function docstring the concise
  source of truth. Keep extension method docs limited to extension-specific
  behavior. Move operational detail to one linked workflow page. Apply this as a
  grouped cleanup instead of editing isolated sentences.

### SW-009: Documentation errors cannot fail the documentation job

- **Kind:** DOC GAP
- **Severity:** High
- **Affected surface:** Documentation validation and navigation
- **Evidence:** `docs/make.jl:92-93` sets `checkdocs = :none` and
  `warnonly = true`. Literate generation itself uses `execute = false` at
  `docs/make.jl:17-24`; Documenter executes only the generated `@example` blocks.
  The CI job runs the full site under Xvfb at `.github/workflows/CI.yml:111-168`
  and provisions 8 to 20 GB of swap because of the hydrology figures at
  `.github/workflows/CI.yml:120-146`.
- **A maintainer needs:** Broken references and executable-example failures to
  stop CI. A smaller validation path should catch links and light examples
  without requiring the full graphics-heavy build.
- **Current documentation:** The one full build is expensive and warning-only.
  Broken links, missing docstrings, or failed examples can be published.
- **Proposed correction:** Add a fast documentation validation mode that skips
  heavy figures but treats reference and example errors as failures. Keep
  `checkdocs = :none` only if undocumented-export coverage is intentionally out of
  scope. Remove `warnonly = true` after existing warnings are resolved.

### SW-010: Tutorial navigation has a dead landing route and an orphan page

- **Kind:** DOC GAP
- **Severity:** Medium
- **Affected surface:** README tutorial link and documentation page tree
- **Evidence:** The README links to `/dev/tutorials/` at `README.md:131-136`, but
  the page tree registers individual tutorial files and no tutorial index at
  `docs/make.jl:59-73`. `docs/src/abstractions.md:1-28` is absent from that page
  tree and has no incoming documentation link; its content overlaps the home and
  architecture introductions.
- **A newcomer needs:** Every landing link to resolve to a real chooser, and every
  maintained page to be reachable from navigation.
- **Current documentation:** The docs home already has a useful tutorial chooser
  at `docs/src/index.md:90-107`, but the README bypasses it. The abstractions page
  is invisible.
- **Proposed correction:** Point the README to the documentation home or add a real
  tutorials index. Remove or merge the orphan abstractions page, or add it to the
  page tree with a distinct purpose.

## Second-order neighbor discovery result

A newcomer can discover the direct queries `ring(grid, cell, 2)` and
`neighbors(grid, cell, 2)` from the neighborhood reference and the opening of the
stencil tutorial (`docs/src/api/neighbors.md:7-30` and
`docs/src/tutorials/stencils.jl:15-38`). They cannot discover a supported
second-order global convolution because the sweep API does not provide one. The
tutorial then moves directly to one-ring `mapneighbors` and repeated diffusion at
`docs/src/tutorials/stencils.jl:106-132`. SW-001 records the missing capability;
SW-002 records the misleading halo and tutorial wording. Documentation should
state this limit plainly until the API can express the operation.

## Validation commands and environment

The repository's documentation job uses Julia 1, Xvfb, a GL stack, and a large
swap file. Its effective commands are:

```sh
xvfb-run -s '-screen 0 1024x768x24' julia --project=docs -e \
  'using Pkg; Pkg.instantiate(; workspace = true)'
DATADEPS_ALWAYS_ACCEPT=true xvfb-run julia --project=docs docs/make.jl
```

On a workstation with a display and the graphics libraries already available,
`julia --project=docs docs/make.jl` is the direct form. This audit did not run
Julia or the full site build. It used static source, test, navigation, and link
inspection. A later fix pass should run the CI-equivalent command because the
tutorials contain executable `@example` blocks even though Literate generation
itself uses `execute = false`.

## Coverage ledger

| Area | Public source inspected | User documentation inspected | Behavioral evidence | Result |
|---|---|---|---|---|
| Neighborhood queries and sweeps | `src/interface/grid.jl`, `src/engine/neighborhood.jl`, neighborhood parts of `src/dimensionaldata.jl` | `docs/src/api/neighbors.md`, `docs/src/tutorials/stencils.jl`, neighborhood sections of `docs/src/tutorials/hydrology.jl` | `test/systems/crosssystem/neighborhood.jl`, `test/systems/crosssystem/mapneighbors.jl` | SW-001, SW-002, SW-003, SW-008 |
| Stencils, cell fields, and requested fields | `src/engine/stencil.jl`, `src/engine/cellfield.jl`, `src/engine/needs.jl` | `docs/src/api/neighbor-fields.md`, `docs/src/tutorials/stencils.jl`, `docs/src/tutorials/hydrology.jl` | Targeted stencil, needs, and cell-field tests | SW-001, SW-003, SW-008 |
| Chunk planning and execution | `src/chunks.jl` | `docs/src/api/chunk-sweep.md`, `docs/src/tutorials/out_of_core.jl` | Chunk-sweep and store-backed sweep tests | SW-001, SW-002, SW-008 |
| Store I/O and subzone storage | `src/io/api.jl`, public contracts in `src/io/description.jl`, `src/io/conventions.jl`, `src/io/encodings.jl`, `src/io/subzones.jl`, and `ext/DiscreteGlobalGridsZarrExt` | `docs/src/api/store-io.md`, `docs/src/api/subzone-layout.md`, `docs/src/tutorials/store_io.jl`, `docs/src/tutorials/out_of_core.jl` | Targeted read, write, convention, encoding, subzone, and round-trip tests | SW-007, SW-008 |
| Regridding | `src/regridding.jl`, public API in `lib/GlobalRegridding/src/api.jl` | `docs/src/api/regridding-methods.md`, `docs/src/tutorials/regridding.jl`, `docs/src/tutorials/between_grids.jl`, regridding sections of `docs/src/extending.md` | Regridding acceptance and plan-reuse tests | SW-004, SW-005 |
| Partitioning | `src/partitioning.jl`, `src/partitioning_backends.jl`, and the Metis, KaHyPar, and Scotch extensions | `docs/src/api/partitioning.md` | Backend and distributed partition tests | SW-005, SW-006 |
| Optional extensions and plotting | Extension declarations in `Project.toml`, `ext/MakieExt`, Zarr and partition extensions | Installation and extension references across `docs/src/index.md`, store, partition, tutorial, and visualization links | `test/plotting/runtests.jl` and targeted extension tests | SW-006, SW-007 |
| Tutorials and navigation | `docs/make.jl`, module exports in `src/DiscreteGlobalGrids.jl` | `README.md`, `docs/src/index.md`, `docs/src/all_dggs.md`, `docs/src/architecture.md`, `docs/src/abstractions.md`, `docs/src/extending.md`, all tutorial source files, and workflow API pages | Static page-tree and cross-reference search | SW-003, SW-004, SW-006, SW-009, SW-010 |

The audit read `docs/src/tutorials/between_grids.jl` as evidence only and did not
modify it. Its existing working-tree change is preserved. Detailed algorithm
internals were inspected only when needed to establish public behavior. Core grid
interfaces, systems, identifiers, geometry, hierarchy, selection, and region-type
contracts remain in the separate core sweep.
