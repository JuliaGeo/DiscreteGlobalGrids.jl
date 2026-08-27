# Regridding successor handoff — 2026-08-22

This work was paused because the current machine ran out of disk space. Resume
on `claude/perf-ladder` from this file and the authoritative task-card plan:
`regrid-notes/2026-08-21-regridding-simplification-plan.md`.

## Repository state

- Main repository HEAD before this handoff commit: `ee586b0`.
- Phase 1A is in progress: A1 and A2 are complete; A3 is next but has an
  upstream ConservativeRegridding prerequisite. A4 has a second upstream CR
  prerequisite.
- A1: `47a7cf7 Simplify raster coordinate transformations`.
- A2 amended result: `41fb204 Add task-local Proj raster transformations`.
  The first A2 commit (`5b5bf1e`) had already reached the remote, so `ee586b0`
  non-destructively merges the two histories while retaining the amended tree.
- The main worktree was clean before this handoff file was added.
- Do not resume or reuse the old subagents; all active agents were interrupted.

## Completed behavior and gates

A1 removed the local longitude/latitude transform types, uses GeometryOps'
`UnitSphereFromGeographic` and `GeographicFromUnitSphere`, normalizes raster
geometry to one-coordinate calls, and adapts legacy two-argument closures once.
The full GlobalRegridding suite passed with 3,137 tests and one expected broken.

A2 added the Proj 1.9 weak extension. A shared locked state retains only the
clone template; every Julia task owns one context and two transformation clones
through a TLS owner. The owner has a registered, idempotent finalizer which
finalizes both transformations before destroying the context, plus complete
partial-construction cleanup through Proj.jl APIs. There is no direct `ccall` or
`PROJ_jll` use. The full suite passed with 3,241 tests and one expected broken.

Proj hot-path measurements over 500,000 warmed calls, best of seven, all zero
allocation:

- direct task-owned clone plus GeometryOps: 71.37 ns/call;
- task-prepared wrapper: 71.37 ns/call;
- safe one-key TLS wrapper: 115.23 ns/call.

Multi-point RasterGrid loops already hoist preparation. A3 must make hot raster
tree traversal use a short-lived task-prepared view; retained/shared indexes
must continue to contain the safe TLS view.

## A3 prerequisite and cutover

The pinned CR `TopDownQuadtreeCursor` is not yet a safe performance replacement
for `RasterCellTree`. In the audit, equivalent raster/raster traversal changed
from 0.234 s / 11.1 MB with the raw local tree to 0.386 s / 327.9 MB with stock
TopDown. Causes: `getchild(q)` rebuilt the complete child tuple for every yielded
child, and fixed 2x2 leaves create far more nodes than the current roughly
16-cell raster leaves. Restricted cursors preserve global IDs, but their local
`ncells` incorrectly sizes direct CR sparse outputs.

An **uncommitted and incompletely tested** prototype exists only in the old
temporary worktree `/tmp/cr-a3-uaZAWa`, branch
`codex/regrid-quadtree-prereq`, based on pinned CR `66ed54c`. It was not pushed
and should be recreated on the new machine. Its intended design is:

- direct one-child construction without rebuilding sibling tuples;
- configurable `leafsize`, default `(2,2)`, with GR requesting `(4,4)`;
- new `Trees.cell_index_count(tree)` for the dense global leaf-ID domain;
- default output sizing uses `cell_index_count`, while `ncells` and
  `split_weight` stay range-local;
- forwarding for tree wrappers and an IndexOffset override;
- focused restricted-cursor/global-ID/output-size/allocation tests.

Prepared-view lifetime law for the GR cutover:

- `LazyRegridArray.srcindex`, chunk indexes, `TileCells.tree`, and
  `CachedCellTree.tree` are shared and must never retain a task clone.
- Synchronous raster `candidatechunks!` may prepare a private cursor once.
- Serial conservative traversal may prepare private cursor copies once.
- Threaded CR frontier nodes cross tasks, so safe roots must remain unprepared.
  Without a CR task-local-tree hook, prepare an ephemeral grid view once per
  `cell_range_extent` call, not once per vertex.
- Make `RasterGridView` carry its callable in a field; the public constructor
  stores the safe callable and a private preparation helper substitutes the
  current task's direct wrapper only in nonescaping copies.

Measured projected `Trees.getvertex`: safe retained view 122.84 ns/call;
prepared local view 72.58 ns/call, zero allocations and identical values.

After the upstream CR commit is published and pinned, A3 should replace whole
and rectangular `RasterCellTree`s with restricted TopDown cursors, preserve
x-fast/y-fast global numbering, benchmark before deleting `MemoRasterTree`, and
leave scattered subsets on `RasterFlatTree` until B4.

## A4 prerequisite and cutover

Pinned CR's generic spherical `cell_range_extent` is unsound for wide/global
ranges: the audit found non-finite caps and finite caps missing descendant
geodesic edges. Current GR's generic fixed-quarter sampling is also unsound for
arbitrary transforms; a smooth injective warp produced 40 escaped vertices and
`candidatechunks! == []` for an owning chunk.

An **uncommitted and incompletely tested** prototype exists only in the old
temporary worktree `/tmp/cr-a4-1tw3kL`, branch
`codex/regrid-range-cap-prereq`, also based on `66ed54c`. It was not pushed and
should be recreated. Intended narrow fix:

- guard a non-finite or near-zero four-corner mean;
- walk the actual grid perimeter without materialization;
- outward-pad the maximum radius;
- return a canonical whole-sphere cap with radius `nextfloat(pi)` for any
  non-finite intermediate or padded radius greater than `pi/2`;
- rely on geodesic convexity only for caps no wider than a hemisphere;
- remove redundant explicit midpoint transforms;
- add exhaustive 8x4 all-rectangle coverage tests, dense edge samples, and
  global, stripe, polar, finite-miss, and injective-warp regressions.

After publishing/pinning that fix, A4 should delegate only non-geographic range
extents to corrected CR. Retain and precisely name the fast geographic analytic
specialization: on 2,000 rectangles of a 720x360 grid it took 0.229 ms versus
5.419 ms for raw CR perimeter traversal. Route chunk and node extents through
one seam, then delete the unsafe generic `_boxcap`/`_tighterof` path. Coverage
tests must use RasterGrid's declared geodesic cell polygons, include actual Proj
charts, check every cursor node and DiskArrays chunk, and assert zero data reads.

## Resume protocol

1. Check disk space, branch/worktree cleanliness, and remote state first.
2. Use Julia MCP for every Julia command; never run Julia through the shell.
3. Do not use Context7 for this refactoring/code-review work.
4. Recreate upstream CR work in isolated worktrees from `66ed54c`; commit and
   test each prerequisite separately. Do not pin unpublished commits or push
   upstream without the user's authorization.
5. Amend the authoritative plan with explicit upstream prerequisite task cards
   before A3/A4 integration.
6. Continue with one writer per worktree, subagents for bounded tasks, focused
   commits, and a phase checkpoint commit after A4.
7. Preserve the user's constraints: DiskArrays owns raster chunks; DGG uses its
   native hierarchy; no direct PROJ `ccall`; nearest/bilinear redesign remains
   deferred; lazy execution must ultimately consume the materialized graph.
