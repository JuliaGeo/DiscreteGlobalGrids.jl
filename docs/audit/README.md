# Newcomer documentation audit

The initial documentation rewrite and the question-level walkthrough of all 75
use cases are complete. The walkthroughs found follow-up documentation work and
originally confirmed five API gaps. Maintainer review subsequently dropped API-002
and clarified API-003, leaving four active proposals. See
[maintainer feedback](maintainer-feedback.md) for the decisions that supersede the audit.

## Read the audit

| File | Purpose | Current state |
| --- | --- | --- |
| [use-cases.md](use-cases.md) | Stable questions and expected outcomes for UC-001–075 | Inventory complete |
| [sweep-core.md](sweep-core.md) and [sweep-workflows.md](sweep-workflows.md) | Historical findings before the rewrite | Preserved as evidence |
| [findings.md](findings.md) | 21 initial documentation corrections | Status updated from the recorded checks |
| [rewrite-initial.md](rewrite-initial.md) | Files changed in the first documentation pass | Rewrite complete |
| [walkthrough-core.md](walkthrough-core.md) | UC-001–024 question-level evidence | 24 of 24 complete |
| [walkthrough-workflows.md](walkthrough-workflows.md) | UC-025–057 question-level evidence | 33 of 33 complete |
| [walkthrough-boundaries.md](walkthrough-boundaries.md) | UC-058–075 question-level evidence | 18 of 18 complete |
| [use-case-findings.md](use-case-findings.md) | Deduplicated second-stage documentation backlog | 28 entries: 27 open and one narrow correction recorded |
| [api-gaps.md](api-gaps.md) | Confirmed missing operations and separate future-scope questions | Four active proposals; API-002 withdrawn |
| [validation.md](validation.md) | Commands, passed checks, and explicit limits | Final local record |

## Walkthrough outcomes

Each report starts at the public documentation entry point and answers every
listed question. It marks source reading and runtime evidence separately.

- UC-001–024: 9 discoverable and 15 with a documentation gap. Existing
  API-001–004 also appear where a case needs a missing operation.
- UC-025–057: 24 discoverable, 7 with a documentation gap, and 2 blocked by
  existing API-005.
- UC-058–075: 2 clear, 13 partial or composable, and 3 unsupported as a direct
  or complete route.

The second-stage ledger maps UCOR-001–012, UW-001–006, and UB-001–016. UB-006
and UB-007 verified existing coverage. The remaining source findings map to 28
deduplicated documentation entries. UCDOC-027 records a high-severity runtime
mismatch: chunked source and destination aliasing can silently differ from an
eager sweep. The walkthrough did not establish a promised in-place API.

Exchange writers, physical-radius kernels, complete adaptive solvers,
accelerators, remote transactions, and broader planetary support remain future
scope questions. They are listed separately from confirmed API gaps.

## Validation

The following checks passed on Julia 1.13.0:

- package load and the package docstring doctest;
- a strict fresh-daemon fast documentation build;
- focused geometry, hierarchy, dimensional-array, regridding, storage,
  neighborhood, chunk, and partitioning examples;
- controlled broken-reference and throwing-example fixtures;
- 41 boundary-walkthrough assertions; and
- an AST comparison of all 40 modified Julia source files after removing only
  documentation payloads and line nodes.

The fast build creates every page but does not execute examples on `index.md`,
`all_dggs.md`, `extending.md`, or the ten tutorial pages. Those 13 pages need
graphics, downloads, optional development packages, or large data. The audit
did not test external downloads, graphics, deployment, remote stores, native
partitioners, or distributed schedulers. `checkdocs=:none` also remains an
explicit limit. See [validation.md](validation.md) for the exact command.

## Status rules

The initial rewrite changes documentation and docstrings only. The source AST
comparison found no other functionality change. Initial findings are marked
Verified only when their own validation evidence passed. Remaining second-stage
entries stay open for a future pass.

The one post-walkthrough correction qualifies A5's area surface in two overview
rows. A5 is equal-area on its ellipsoid; its unit-sphere `cell_area` varies by
about 1%. No other second-stage backlog item was fixed.

Four API proposals remain open after maintainer review. API-002 is withdrawn.
A custom loop, wider halo, or adapter does not close the active proposals. Preserve the existing user edit in
`docs/src/tutorials/between_grids.jl`.
