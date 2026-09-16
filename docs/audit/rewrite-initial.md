# Initial documentation rewrite checkpoint

## Scope and status

DOC-001–019 and DOC-021 are **Revised, awaiting validation**.
DOC-020 has a separate validation-infrastructure owner.
All five API gaps remain open. This pass changes documentation, not API behavior.
The 75 newcomer walkthroughs have not been performed by the rewrite owner.

## Changes

- Corrected query result types, predicate direction, target support, and collection conversion guidance.
- Added system discovery for S2 and CopernicusDEM, plus identifier and coordinate guidance.
- Separated fixed-level coverage from budget-mode guarantees and stated the over-budget seed exception.
- Explained current centroid-selector, collection-polygon, regional-resolution, and wider-sweep limitations without adding workarounds.
- Distinguished one-ring callbacks, exact rings, radius-two disks, wider input halos, and repeated diffusion.
- Shortened key grid, collection, geometry-bound, neighborhood, dimensional-array, and IO docstrings.
- Preserved detailed implementation contracts in four linked internals pages.
- Added a regridding call reference and focused geometry, hierarchy, pass-mode, regridding, and partitioning examples.
- Added an optional-integration table, S3 activation guidance, and navigation repairs.

New reference pages:

- `docs/src/api/regridding.md`
- `docs/src/internals/grid-contracts.md`
- `docs/src/internals/collection-contracts.md`
- `docs/src/internals/system-capabilities.md`
- `docs/src/internals/workflow-contracts.md`

The validation owner controls `docs/make.jl`, including navigation additions and local validation mode.

## Validation evidence

`git diff --check` passed after the rewrite and reference cleanup.
A static comparison removed triple-quoted documentation from each changed source and extension file.
The remaining code matched HEAD in all 36 such files at this checkpoint.
The existing `between_grids.jl` user deletion remains unchanged: its diff contains one added and two removed lines.
The rewrite owner did not execute Julia or deploy documentation.

The first strict build run by the validation owner failed on cross-references.
Its log was `/Users/anshul/.cache/julia-daemon/docs-sweep_core_fastdocs-8fb08acb/daemon.log`.
That run used some source docstrings loaded before the final rewrite was saved.
The next run must load the revised package in a fresh Julia process.

The reference cleanup qualified actual owner-module bindings, added `StorageOrder` and `NeighborCallbackError` reference entries,
and replaced links to unrendered implementation-only names with code formatting.
This includes clear baseline reference failures in store, subzone, and boundary documentation.
A strict rerun is pending; no DOC entry is marked Verified.

## Handoff

The validation owner will report remaining reference or focused-example failures.
The next audit stage assigns all 75 use cases to newcomer walkthroughs.
Record findings first and consolidate them before another user-journey correction pass.
Do not close API gaps through documentation-only corrections or custom-loop substitutes.
