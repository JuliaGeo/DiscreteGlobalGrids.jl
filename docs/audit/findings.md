# Documentation findings and correction plan

This ledger consolidates the initial core and workflow sweeps. The first documentation rewrite is complete.
Sixteen entries are **Verified** from the evidence named below. Five remain **Revised, awaiting validation** because a required exceptional or external path was not checked. The linked API work remains open.

The historical evidence remains in [sweep-core.md](sweep-core.md) and [sweep-workflows.md](sweep-workflows.md).
Use the source IDs to recover original line references. Those references describe the initial sweep and can shift after revisions.
Use-case IDs refer to [use-cases.md](use-cases.md). API work remains separate in [api-gaps.md](api-gaps.md).

## Priority and status

- **P0:** Correct a false contract or a statement likely to cause silent misuse first.
- **P1:** Correct a major discovery, interpretation, or validation gap next.
- **P2:** Add missing workflow guidance or improve navigation.
- **P3:** Improve wording consistency after the stronger issues.

A documentation entry can close while its linked API gap remains open.
For example, explaining the one-ring limit fixes misleading documentation but does not add radius-two sweeps.
Do not add a workaround and claim that the missing natural operation now exists.

After revision, use **Revised, awaiting validation**, then **Verified** only with recorded evidence.
Keep **Open** when a required correction or check remains unfinished.

## Deduplicated entries

### DOC-001: State the exact query result containers

- **Source findings:** SC-001
- **Priority:** P0
- **Use cases:** [UC-011](use-cases.md), [UC-014](use-cases.md), [UC-015](use-cases.md)
- **Affected files:** docs/src/api/selecting-cells.md; src/engine/query.jl; src/engine/multiorder.jl
- **Problem:** The guide calls a fixed-level query result a CellVector; the implementation returns a sorted Vector of cell IDs.
- **Planned correction:** Add a fixed-level versus mixed-level result table. Show supported conversions to a region or cell axis.
- **Validation:** Check result types and value alignment in one small fixed-level and one mixed-level query.
- **Linked API work:** None
- **Status:** Verified. UC-011 confirmed the fixed-level and mixed-level result containers, conversions, and value alignment with static and executed evidence.

### DOC-002: Provide one complete system-choice overview

- **Source findings:** SC-002, SC-003
- **Priority:** P1
- **Use cases:** [UC-002](use-cases.md), [UC-005](use-cases.md), [UC-008](use-cases.md), [UC-059](use-cases.md)
- **Affected files:** src/DiscreteGlobalGrids.jl; docs/src/tutorials/choosing_a_grid.jl; docs/src/all_dggs.md
- **Problem:** The registry description sounds exhaustive, while the choice table omits S2 and the user guides omit CopernicusDEM.
- **Planned correction:** Define the registry scope. Add S2 to the main comparison and CopernicusDEM as a specialized lattice. Summarize geometry, coordinates, identifier compatibility, and limits.
- **Validation:** Compare every exported system constructor with the overview. Check registry membership and documented compatibility restrictions.
- **Linked API work:** None
- **Status:** Verified. The system inventory and compatibility restrictions were checked across UC-002, UC-005, UC-008, and UC-059. UCDOC-004 records and corrects the later A5 surface-label inconsistency.

### DOC-003: Separate task reference from implementation contracts

- **Source findings:** SC-004, SW-008
- **Priority:** P1
- **Use cases:** [UC-007](use-cases.md), [UC-009](use-cases.md), [UC-010](use-cases.md), [UC-014](use-cases.md), [UC-063](use-cases.md), [UC-064](use-cases.md), [UC-065](use-cases.md), [UC-075](use-cases.md), [UC-027](use-cases.md), [UC-028](use-cases.md), [UC-034](use-cases.md), [UC-045](use-cases.md), [UC-052](use-cases.md), [UC-053](use-cases.md)
- **Affected files:** docs/src/api/grid-interface.md; docs/src/extending.md; docs/src/internals/boundary-engines.md; src/DiscreteGlobalGrids.jl; src/interface/system.jl; src/engine/cell_vector.jl; src/engine/neighborhood.jl; src/dimensionaldata.jl; src/io/api.jl; ext/DiscreteGlobalGridsZarrExt/read.jl; docs/src/api/neighbor-fields.md; docs/src/api/chunk-sweep.md; docs/src/api/store-io.md
- **Problem:** Usage details compete with long traversal contracts, representation details, and engine rationale. Workflow docstrings also repeat scheduling and storage contracts across generic functions, extensions, and guides.
- **Planned correction:** Lead public docstrings with inputs, output, units, limits, and one normal example. Move author-only contracts to extension or internals pages with links. Preserve technical guarantees. Keep the generic function contract authoritative. Limit extension docs to extension-specific behavior and place operational detail on one linked guide.
- **Validation:** Review representative user paths and contributor paths. Check that moved contract statements remain discoverable and accurate.
- **Linked API work:** None
- **Status:** Revised, awaiting validation. The strict build and all 75 walkthroughs checked representative paths, but UCDOC-007, UCDOC-016, UCDOC-023, and UCDOC-024 record remaining lifecycle, extension-shape, and callback-contract questions.

### DOC-004: Make predicate direction and target support explicit

- **Source findings:** SC-005
- **Priority:** P0
- **Use cases:** [UC-011](use-cases.md), [UC-012](use-cases.md), [UC-013](use-cases.md)
- **Affected files:** docs/src/api/selecting-cells.md; src/engine/query.jl
- **Problem:** Asymmetric relations can be read backward. The predicate list implies Crosses and all cap predicates are supported.
- **Planned correction:** State cell RELATION target explicitly. Add a supported-target matrix and asymmetric examples. Identify Crosses rejection and the actual cap predicate subset.
- **Validation:** Compare the matrix with predicate dispatch and rejection paths. Run one Contains/Within pair and an unsupported case.
- **Linked API work:** None
- **Status:** Verified. UC-011 and UC-012 checked relation direction, supported targets, and rejection behavior against the rendered contracts and runtime examples.

### DOC-005: Separate multi-order guarantees from refinement budgets

- **Source findings:** SC-007, SC-008
- **Priority:** P0
- **Use cases:** [UC-019](use-cases.md), [UC-021](use-cases.md), [UC-022](use-cases.md), [UC-023](use-cases.md), [UC-024](use-cases.md)
- **Affected files:** src/engine/multiorder.jl; src/engine/multiorder_budget.jl; docs/src/tutorials/multiorder.jl; docs/src/api/selecting-cells.md
- **Problem:** Coverage promises conflict across level and maxcells modes. The tutorial calls maxcells a hard cap although the seed can exceed it.
- **Planned correction:** Describe each mode separately. State guarantees only where established. Explain noncongruent refinement and seed-budget exceptions at first use. Label empirical measurements as measurements.
- **Validation:** Check descriptions against multi-order tests. Exercise a seed larger than the budget and compare congruent with noncongruent systems.
- **Linked API work:** None
- **Status:** Revised, awaiting validation. UC-021–023 confirm the written guarantees and a bounded budget run, but this pass did not execute an over-budget seed or compare congruent and noncongruent targets at runtime.

### DOC-006: Explain collection choices, conversions, and ownership

- **Source findings:** SC-009, SC-010
- **Priority:** P1
- **Use cases:** [UC-007](use-cases.md), [UC-014](use-cases.md), [UC-015](use-cases.md), [UC-017](use-cases.md), [UC-020](use-cases.md), [UC-021](use-cases.md)
- **Affected files:** docs/src/api/selecting-cells.md; src/engine/partial_grid.jl; src/engine/cell_vector.jl; src/dimensionaldata.jl
- **Problem:** There is no task-oriented conversion map. PartialGrid also requires canonical sorted unique IDs and retains its input storage.
- **Planned correction:** Add a conversion table by task and result type. State existing ordering, element-type, uniqueness, and aliasing requirements. Keep the missing owning constructor in the API ledger.
- **Validation:** Use a noncontiguous axis to check conversion ordering. Check the documented PartialGrid preconditions and ownership against its constructor.
- **Linked API work:** API-002
- **Status:** Verified. UC-007, UC-014, UC-020, and UC-021 checked conversion, position, ordering, and ownership contracts; the noncontiguous-position probe passed.

### DOC-007: Show geographic coordinates beside spherical geometry APIs

- **Source findings:** SC-011
- **Priority:** P2
- **Use cases:** [UC-005](use-cases.md), [UC-006](use-cases.md), [UC-009](use-cases.md)
- **Affected files:** docs/src/api/grid-interface.md; src/interface/grid.jl; docs/src/tutorials/choosing_a_grid.jl
- **Problem:** The core page returns unit-sphere coordinates but leaves longitude/latitude conversion in a specialized tutorial.
- **Planned correction:** Add a short centroid and boundary conversion example. State angle units, coordinate order, longitude convention, and the relevant frame guidance.
- **Validation:** Execute the example on a known non-equatorial cell and check geographic-to-sphere round trips.
- **Linked API work:** None
- **Status:** Revised, awaiting validation. UC-005 and UC-009 found the coordinate units and frame guidance discoverable, but the required non-equatorial geographic round trip was not rerun.

### DOC-008: Give hierarchy operations a runnable usage path

- **Source findings:** SC-013
- **Priority:** P2
- **Use cases:** [UC-019](use-cases.md), [UC-020](use-cases.md), [UC-021](use-cases.md)
- **Affected files:** docs/src/api/grid-interface.md; docs/src/api/boundaries.md; src/interface/system.jl
- **Problem:** Reference lists do not explain when to choose descendants, descendant_range, or subtree.
- **Planned correction:** Add a task table and one example from a parent through descendants to a region. Put the sorted-subtree availability limit beside descendant_range.
- **Validation:** Check the example on a supported system and identify A5 behavior explicitly.
- **Linked API work:** None
- **Status:** Verified. UC-019 and UC-020 followed the hierarchy route and executed children, descendants, descendant ranges, and subtree construction.

### DOC-009: Describe descendants without promising allocation or mutability

- **Source findings:** SC-014
- **Priority:** P2
- **Use cases:** [UC-019](use-cases.md), [UC-059](use-cases.md), [UC-064](use-cases.md), [UC-075](use-cases.md)
- **Affected files:** src/interface/system.jl; src/systems/CopernicusDEM/system.jl
- **Problem:** The generic description promises materialization while CopernicusDEM returns a lazy read-only vector.
- **Planned correction:** State only the common supported result operations, order, and mutability guarantees. Describe storage cost as implementation-specific. Preserve the actual interface contract.
- **Validation:** Compare all implementations before choosing an AbstractVector or iterable promise. Check CopernicusDEM and an ordinary system.
- **Linked API work:** None
- **Status:** Revised, awaiting validation. Static inspection confirmed the abstract result and CopernicusDEM contract, but the required ordinary and CopernicusDEM descendant comparison was not executed.

### DOC-010: Use explicit numeric level direction

- **Source findings:** SC-016
- **Priority:** P3
- **Use cases:** [UC-019](use-cases.md), [UC-020](use-cases.md), [UC-021](use-cases.md)
- **Affected files:** src/engine/region_algebra.jl; src/interface/system.jl; hierarchy examples
- **Problem:** Above and below can reverse the meaning of finer and coarser levels.
- **Planned correction:** Use finer/coarser and explicit inequalities. State that larger level numbers represent finer cells.
- **Validation:** Review every changed level condition against its argument checks.
- **Linked API work:** None
- **Status:** Verified. UC-019–021 checked the numeric level direction against argument contracts and executed hierarchy examples.

### DOC-011: Identify the absent centroid-selection operation

- **Source findings:** SC-006
- **Priority:** P1
- **Use cases:** [UC-012](use-cases.md), [UC-018](use-cases.md), [UC-070](use-cases.md)
- **Affected files:** docs/src/tutorials/zonal.jl; docs/src/api/selecting-cells.md
- **Problem:** The zonal comparison includes a center-in-zone rule but supplies no natural API spelling.
- **Planned correction:** State that the current query predicates test cell footprints and that a centroid selector is unavailable. Link the open API gap. Do not add a manual filter as its replacement.
- **Validation:** Check that the table distinguishes existing calls from proposed capability.
- **Linked API work:** API-001
- **Status:** Verified. UC-012 and UC-018 confirmed that footprint predicates do not provide centroid selection and that API-001 remains open.

### DOC-012: Narrow the collection polygon claim to supported inputs

- **Source findings:** SC-012
- **Priority:** P0
- **Use cases:** [UC-009](use-cases.md), [UC-014](use-cases.md), [UC-022](use-cases.md), [UC-061](use-cases.md), [UC-069](use-cases.md)
- **Affected files:** docs/src/api/boundaries.md; src/engine/multiorder_set.jl
- **Problem:** The guide recommends cell_polygons for collections although the only public method accepts MultiOrderCellSet.
- **Planned correction:** State the current accepted input precisely. Keep common collection support open in the API ledger. Do not imply that a documentation correction adds methods.
- **Validation:** Inspect methods and check the guide names only supported inputs.
- **Linked API work:** API-003
- **Status:** Verified. UC-009 and UC-022 checked the rendered claim against the public methods and confirmed the current MultiOrderCellSet limit; API-003 remains open.

### DOC-013: Correct the current levelfor area-of-interest contract

- **Source findings:** SC-015
- **Priority:** P0
- **Use cases:** [UC-003](use-cases.md), [UC-004](use-cases.md), [UC-074](use-cases.md)
- **Affected files:** src/sizing.jl; docs/src/api/grid-interface.md; docs/src/tutorials/choosing_a_grid.jl
- **Problem:** The over keyword is described as restricting both inputs, but it restricts candidate system levels only.
- **Planned correction:** Describe the current asymmetry and its effect on nonuniform targets. Keep symmetric regional comparison as open API work.
- **Validation:** Trace over through levelfor and target measurement. Check that the text does not promise regional target sampling.
- **Linked API work:** API-004
- **Status:** Verified. UC-003, UC-004, and UC-074 traced target and candidate measurement and confirmed that the revised text states the current asymmetry; API-004 remains open.

### DOC-014: Separate halo data from callback neighborhood reach

- **Source findings:** SW-001, SW-002
- **Priority:** P0
- **Use cases:** [UC-025](use-cases.md), [UC-026](use-cases.md), [UC-027](use-cases.md), [UC-030](use-cases.md), [UC-052](use-cases.md), [UC-054](use-cases.md)
- **Affected files:** src/chunks.jl; src/engine/neighborhood.jl; src/dimensionaldata.jl; docs/src/api/chunk-sweep.md; docs/src/api/neighbors.md; docs/src/tutorials/out_of_core.jl; docs/src/tutorials/stencils.jl
- **Problem:** Wider halo wording implies a wider built-in stencil, but all sweep callbacks receive one-ring neighbors.
- **Planned correction:** State the one-ring limit next to sweep entry points. Explain that halo width controls available chunk data. Distinguish radius-two disks, exact rings, and repeated one-ring diffusion. Keep the missing wider sweep open without prescribing a custom-loop substitute.
- **Validation:** Inspect callback membership with halo=2. Check that every second-order example states whether it is selection, a missing sweep, or repeated diffusion.
- **Linked API work:** API-005
- **Status:** Verified. UC-025–027, UC-030, UC-052, and UC-054 checked one-ring callback membership separately from halo width. API-005 remains open.

### DOC-015: Compare all dimensional-array pass modes beside the sweep

- **Source findings:** SW-003
- **Priority:** P2
- **Use cases:** [UC-015](use-cases.md), [UC-027](use-cases.md), [UC-028](use-cases.md), [UC-034](use-cases.md)
- **Affected files:** docs/src/api/neighbors.md; docs/src/api/selecting-cells.md; src/dimensionaldata.jl
- **Problem:** Neighbors, Values, and NeighborSlices lack one task-level comparison, and NeighborSlices is filed under selection.
- **Planned correction:** Add a table of callback inputs, output shape, and intended use. Place all mode reference entries with mapneighbors. Include a time-by-cells example.
- **Validation:** Compare Values processing per time slice with NeighborSlices receiving full time series. Check dimensions against existing tests.
- **Linked API work:** None
- **Status:** Verified. The focused Values and NeighborSlices examples passed, and UC-027, UC-028, and UC-034 confirmed the documented dimensions and callback modes.

### DOC-016: Document regridding calls and execution controls

- **Source findings:** SW-004
- **Priority:** P1
- **Use cases:** [UC-036](use-cases.md), [UC-037](use-cases.md), [UC-038](use-cases.md), [UC-039](use-cases.md), [UC-040](use-cases.md), [UC-041](use-cases.md), [UC-042](use-cases.md), [UC-043](use-cases.md), [UC-044](use-cases.md)
- **Affected files:** docs/src/api/regridding-methods.md or a dedicated regridding API page; docs/make.jl; README.md; src/regridding.jl; lib/GlobalRegridding/src/api.jl
- **Problem:** The method guide omits callable contracts and DGGSpace, leaving users to infer keywords and execution controls from source.
- **Planned correction:** Provide one-shot, reusable-plan, and in-place entry paths. Render the owned callable and policy docstrings. State accepted spaces, result shapes, memory/storage controls, and plan validity. Link to method choice separately.
- **Validation:** Build the reference locally with deployment disabled. Run small one-shot and reused-plan examples; check that all documented keywords exist.
- **Linked API work:** None
- **Status:** Verified. The strict reference build and focused one-shot, reusable-plan, and partitioning examples passed; UC-036–044 checked the public call routes.

### DOC-017: Make partitioning examples executable or explicitly schematic

- **Source findings:** SW-005
- **Priority:** P2
- **Use cases:** [UC-043](use-cases.md), [UC-055](use-cases.md), [UC-056](use-cases.md), [UC-057](use-cases.md)
- **Affected files:** docs/src/api/partitioning.md; docs/Project.toml if an explicit dependency is necessary
- **Problem:** Examples use undefined plans and application variables. One imports a package that is not a direct docs dependency.
- **Planned correction:** Make the smallest public-plan example executable. Label the application-specific example as schematic and define placeholders. Prefer the available re-exported surface where sufficient.
- **Validation:** Execute the small example in the docs environment and check that partitions cover every task exactly once.
- **Linked API work:** None
- **Status:** Verified. The focused synthetic and lazy-plan partition examples passed and assigned every task once; UC-055–057 checked the planning route.

### DOC-018: Provide one optional-integration discovery table

- **Source findings:** SW-006
- **Priority:** P2
- **Use cases:** [UC-001](use-cases.md), [UC-045](use-cases.md), [UC-056](use-cases.md), [UC-061](use-cases.md), [UC-062](use-cases.md)
- **Affected files:** docs/src/index.md; README.md; relevant plotting, store, and partition guides
- **Problem:** Activation instructions are scattered, and core Makie conversion is not distinguished from the visualization companion.
- **Planned correction:** List capability, activating package, activation spelling, and reference link for all declared extensions. Explain the scope of core Makie conversion and companion plot types.
- **Validation:** Compare the table with Project.toml extension declarations. Check each linked activation instruction and plotting ownership.
- **Linked API work:** None
- **Status:** Revised, awaiting validation. The extension table matches the declared extensions, but UC-061 and UC-062 found that plotting installation and method ownership still need a direct task route in UCDOC-015.

### DOC-019: Qualify remote store prerequisites and encoding scope

- **Source findings:** SW-007
- **Priority:** P2
- **Use cases:** [UC-045](use-cases.md), [UC-046](use-cases.md), [UC-047](use-cases.md), [UC-073](use-cases.md)
- **Affected files:** docs/src/tutorials/store_io.jl; docs/src/api/store-io.md; ext/DiscreteGlobalGridsZarrExt/read.jl
- **Problem:** The S3 example omits AWSS3 activation. The tutorial calls its two demonstrated encodings exhaustive despite implicit complete-level encoding.
- **Planned correction:** Name using AWSS3 beside s3://. Describe dense/ranges as this regional example and link to implicit encoding eligibility.
- **Validation:** Check protocol requirements against the reader contract and encoding eligibility against implementation. Test local encoding examples without requiring remote publication.
- **Linked API work:** None
- **Status:** Verified. UC-045–047 checked activation and encoding eligibility, and the local storage examples passed. Remote publication remains explicitly untested and is tracked in UCDOC-021.

### DOC-020: Add reliable local documentation validation

- **Source findings:** SW-009
- **Priority:** P1
- **Use cases:** [UC-001](use-cases.md), [UC-068](use-cases.md)
- **Affected files:** docs/make.jl; .github/workflows/CI.yml; contributor build guidance
- **Problem:** The full documentation build is expensive and warning-only. A successful process does not establish valid references or examples.
- **Planned correction:** Design a fast local validation path with deployment disabled and lightweight executable examples. Resolve existing warnings before making applicable errors fatal. Decide and state undocumented-export coverage separately.
- **Validation:** Run the chosen local mode. Confirm that controlled broken-reference and example failures are detected without deploying or requiring the heavy tutorial assets.
- **Linked API work:** None
- **Status:** Verified. The strict fast build, focused examples, and two fatal-failure fixtures passed on 2026-09-15; see [validation.md](validation.md).

### DOC-021: Repair the tutorial landing route and orphan page

- **Source findings:** SW-010
- **Priority:** P2
- **Use cases:** [UC-001](use-cases.md), [UC-014](use-cases.md), [UC-068](use-cases.md)
- **Affected files:** README.md; docs/make.jl; docs/src/index.md; docs/src/abstractions.md
- **Problem:** The README tutorial route has no registered landing page. The abstractions page is absent from navigation and duplicates other introductions.
- **Planned correction:** Point to an existing chooser or add a real tutorial index. Give the abstractions page a distinct linked purpose, or merge its useful content.
- **Validation:** Check all local landing targets and incoming links against the page tree.
- **Linked API work:** None
- **Status:** Verified. The strict build resolved the landing and navigation targets, and UC-001, UC-014, and UC-068 followed the linked routes.

## Source-to-ledger map

Every source finding appears below. Cross-listing is deliberate only when truthful current documentation and later API work require separate outcomes.

| Source finding | Documentation entry | API entry | Reason for shared tracking |
| --- | --- | --- | --- |
| SC-001 | DOC-001 | None | One documentation correction track. |
| SC-002 | DOC-002 | None | One documentation correction track. |
| SC-003 | DOC-002 | None | One documentation correction track. |
| SC-004 | DOC-003 | None | One documentation correction track. |
| SC-005 | DOC-004 | None | One documentation correction track. |
| SC-006 | DOC-011 | API-001 | State the existing limit now; implement the natural operation later. |
| SC-007 | DOC-005 | None | One documentation correction track. |
| SC-008 | DOC-005 | None | One documentation correction track. |
| SC-009 | DOC-006 | None | One documentation correction track. |
| SC-010 | DOC-006 | API-002 | State the existing limit now; implement the natural operation later. |
| SC-011 | DOC-007 | None | One documentation correction track. |
| SC-012 | DOC-012 | API-003 | State the existing limit now; implement the natural operation later. |
| SC-013 | DOC-008 | None | One documentation correction track. |
| SC-014 | DOC-009 | None | One documentation correction track. |
| SC-015 | DOC-013 | API-004 | State the existing limit now; implement the natural operation later. |
| SC-016 | DOC-010 | None | One documentation correction track. |
| SW-001 | DOC-014 | API-005 | State the existing limit now; implement the natural operation later. |
| SW-002 | DOC-014 | API-005 | State the existing limit now; implement the natural operation later. |
| SW-003 | DOC-015 | None | One documentation correction track. |
| SW-004 | DOC-016 | None | One documentation correction track. |
| SW-005 | DOC-017 | None | One documentation correction track. |
| SW-006 | DOC-018 | None | One documentation correction track. |
| SW-007 | DOC-019 | None | One documentation correction track. |
| SW-008 | DOC-003 | None | One documentation correction track. |
| SW-009 | DOC-020 | None | One documentation correction track. |
| SW-010 | DOC-021 | None | One documentation correction track. |

## Recommended rewrite scope

1. Correct false result, predicate, coverage, polygon, resolution, and halo claims: DOC-001, DOC-004, DOC-005, DOC-012, DOC-013, DOC-014.
2. Establish the newcomer route: DOC-002, DOC-006, DOC-007, DOC-008, DOC-011, DOC-015, DOC-016, DOC-017, DOC-018, DOC-019, DOC-021.
3. Rewrite long docstrings as one coordinated pass: DOC-003, DOC-009, DOC-010. Preserve extension contracts through linked author documentation.
4. Make validation reviewable and local: DOC-020. Validate focused examples before the full graphics-heavy site build.
5. Re-run assigned use cases from the published entry points. Record untested cases as pending, not passed.

The sequence is a correction plan, not an instruction to implement API gaps during the documentation pass.
Preserve the existing user edits in `docs/src/tutorials/between_grids.jl` throughout later work.
