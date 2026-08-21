# A caching dual DFS in ConservativeRegridding, not in GeometryOps

2026-08-21. CR branch `claude/cached-dual-dfs`, pushed, no PR.

## The decision this implements

GeometryOps' `SpatialTreeInterface.dual_depth_first_search` is written to be as generic
as it can be, and the trees it is generic over — RTree, STRtree, `NaturalIndex`,
`FlatNoTree` — all *store* their extents. Caching child extents is pure overhead for
them. Regridding is the place where the other kind of tree turns up: a cursor over a
curvilinear grid, and downstream a DGGS cell cursor, has no extent to read, and derives
one by folding the perimeter of its index range into a bounding cap.

So the caching variant belongs to ConservativeRegridding, where those trees meet the dual
search, and GeometryOps' search stays as it is. GO PR #474 prototyped the caching in GO;
this is that design, ported into CR.

## What was ported, and what changed on the way

Source: GO branch `as/child-extents-alloc` (PR #474), three commits `cc8d59e8`,
`7dcc9b4a`, `e6c3a1e2`, against merge base `02750768`.

New file `src/utils/CachedDualDepthFirstSearch.jl` in CR — a module exporting
`cached_dual_depth_first_search` and defining the trait `children_extent_type`.

Ported verbatim in structure:

- the scratch stack: `_extent_stack` / `_as_stack` / `_stack_base` / `_carry` /
  `_fill_child_extents!` / `_child_extent`, and the descent that pushes node2's child
  extents once, reads them `nchild(node1)` times, and `resize!`s them off on the way out;
- the reuse rule — siblings share an extent type so one stack serves a node's children,
  and a level whose child type differs from the stack it inherited starts its own;
- `children_extent_type(node)`, asked of the *node* (falling back to the type), so both
  spellings fold to a compile-time constant;
- the explicit `return nothing` so `resize!`'s return value never escapes as if it were
  an `Action`.

Adapted:

- **The opt-in trait stays GeometryOps'.** `STI.node_extent_is_expensive` is already on
  GO main (it predates #474), it means exactly the right thing, and every DGG tree
  already declares it. CR reads it rather than introducing a parallel trait, so a tree
  that has opted in gets CR's caching for free and nothing needs re-declaring.
- **`children_extent_type` is CR's**, defined in
  `ConservativeRegridding.CachedDualDepthFirstSearch`, `@public` from the top-level
  module, defaulting to `nothing` = "same type as the node's own extent". Plain RTrees
  and STRtrees need no method.
- All `SpatialTreeInterface` calls are qualified `STI.…`; `@controlflow` comes from
  `GeometryOps.LoopStateMachine`.

Wired in (`src/regridder/intersection_areas.jl`, `src/utils/MultithreadedDualDepthFirstSearch.jl`):
both arms of `get_all_candidate_pairs` now run the caching search — the serial one
directly, the threaded one inside each frontier task, where each task descends its own
node pair and so gets its own stack. Nothing is shared across tasks.

Added beyond #474 (CR's own trees, which are the same case):
`Trees.extent_is_expensive(grid)` carries the grid-level fact, and
`QuadtreeCursor` / `TopDownQuadtreeCursor` / `IndexOffsetQuadtreeCursor` /
`ReorderedTopDownQuadtreeCursor` forward it to `STI.node_extent_is_expensive`. True for
every grid but `RegularGrid{<:Planar}`, where `cell_range_extent` is four array reads.
The wrappers that forward `node_extent` (`WithParallelizePolicy`,
`IndexLocalizerRewrapperTree`, `GeometryMaintainingTreeWrapper`) forward the trait with
it; `KnownFullSphereExtentWrapper`, which overrides `node_extent` with a constant, keeps
the `false` default.

Public API of CR is otherwise unchanged.

## Commits

```
4797607  Add a caching dual depth-first search for computed-extent trees
ad8a32b  Discover candidate pairs through the caching search
ff37e8c  Opt the quadtree cursors into cached child extents
66ed54c  Cover the caching search against the generic one
```

Base: `claude/parallel-assembly` @ `3b50a2f64ee824f2e48988c06ed6f74237fc6431` (PR #138).

## Tests

| | pass | fail | error |
|---|---|---|---|
| CR suite at `3b50a2f` (baseline, separate worktree) | 10234 | 0 | 0 |
| CR suite at `claude/cached-dual-dfs` | 10289 | 0 | 0 |

The 55 new assertions are `test/cached_dual_dfs.jl`. Most of them are the drop-in
contract stated directly: same pairs, **same order**, as `STI.dual_depth_first_search`,
over trees taking each arm of the caching decision — spherical `QuadtreeCursor` and
`TopDownQuadtreeCursor` (cached), planar regular cursors and a `FlatNoTree` of polygons
(not cached), mismatched resolutions, a cell-based grid against a regular one, and both
arms of `get_all_candidate_pairs`. The rest are the mechanics: that `_extent_stack`
folds to a concrete `Vector{SphericalCap{Float64}}` from the node type alone and reuses
an existing stack of the right type, that an opted-in tree derives strictly fewer
extents than the same tree opted out, that the descent allocates less than GeometryOps'
per-node-pair vectors and exactly as much as the generic search when the tree opts out,
that `children_extent_type` is honoured on the node and on the type (including a tree
whose extent type alternates by level), and that an `Action(:full_return)` still steers
the traversal.

Two pre-existing CR facts the tests had to route around, both unrelated to this change:

- a whole-sphere point set gives the all-NaN bounding cap, and a NaN cap intersects
  nothing, so a global-grid pair can prune to *no* candidates at all. The
  cross-resolution cases use regional grids so the comparison is not vacuous.
- `Trees.split_weight(::QuadtreeCursor)` goes through `ncells(cursor)`, and
  `QuadtreeCursor` has no no-argument `ncells` method, so the threaded frontier cannot
  take one. The threaded-arm test uses `TopDownQuadtreeCursor`.

## Benchmark

Scratch worktree `DGG-crcache-scratch` at DGG `0d0f069`, Manifest copied from
`DGG-perf-ladder`. The only worktree edit was the `[sources]` line for
ConservativeRegridding (git rev `6a4b997` ⇄ `path = /home/asinghvi17/geo/CR-cached-dfs`)
so the same script runs against old and new CR. **No DGG call-site adaptation was
needed** — consistent with what the parallel `claude/cr-pr138` job found. Nothing was
committed in the worktree.

Workload: IGeo7 level 5 (168,072 cells) against a 360x180 global raster (64,800 cells),
Julia 1.12.6, 8 threads, `nice -n 10`. Best of 3.

| discovery case | tree on side 2 | old wall | new wall | old alloc | new alloc |
|---|---|---|---|---|---|
| `get_all_candidate_pairs` serial | `CapCachedTree` (bucket 49) | 0.484 s | 0.452 s | 115.8 MiB | 114.7 MiB |
| `get_all_candidate_pairs` threaded | `CapCachedTree` | 0.087 s | 0.089 s | 126.7 MiB | 125.8 MiB |
| bare `HierarchicalGridCursor` dst, serial | bare cursor | 3.292 s | 3.300 s | 119.2 MiB | **80.9 MiB (-32.1%)** |
| bare `HierarchicalGridCursor` dst, threaded | bare cursor | 0.432 s | 0.444 s | 130.1 MiB | **92.4 MiB (-28.9%)** |
| raster destination (DGGS -> raster) | raw `RasterCellTree` | 3.866 s | 4.009 s | 114.9 MiB | 114.7 MiB |
| `intersection_areas`, whole space | `CapCachedTree` | 0.433 s | 0.564 s | 838.8 MiB | 840.3 MiB |

Candidate pairs are identical in count *and* hash in all five discovery
configurations (733,220 / 733,220 / 732,806 / 733,220 / 732,806). The weights from
`intersection_areas` over the whole space are identical **bit for bit** — `colptr`,
`rowval` equal, and every `nzval` compares `===` — 487,174 nonzeros, sum
`5.100658809728706e14`, and a spot-checked destination column matches entry for entry.

Reading the numbers honestly:

- **The win is allocation, not wall time.** GO main already avoids the *repeated
  `node_extent` calls* — `node_extent_is_expensive` landed before #474 and makes the
  generic search build a `[node_extent(child) for child in getchild(node)]` per visited
  node pair. What #474 (and this port) removes is that per-node-pair vector. Extent-call
  counts are unchanged, so wall time is flat to within run-to-run noise (+-3%; the
  `intersection_areas` row's -30% is one un-repeated measurement dominated by GC timing,
  not signal).
- **The win scales with visited internal node pairs on side 2.** Only the *second*
  tree's child extents are re-derived per opposing child, so only side 2 benefits. A
  bare 7-ary DGGS cursor gives 32%; `CapCachedTree`, whose `_bucketed` leaf size of 49
  collapses most of the hierarchy, gives 1%; a binary `RasterCellTree` gives 0.2%
  because its two-element vectors are noise against the `_rectcap` work itself.
- DGG's own conservative regrid puts the *destination* on side 2 and hands it over as
  `CapCachedTree`, so the default path sees the 1% row, not the 32% one. The 32% row is
  what a bare cursor destination costs today.

## Retargeting audit

| branch | PR | what is on it | genuinely GO? | action |
|---|---|---|---|---|
| `as/child-extents-alloc` | #474 | 3 commits, all one feature: the scratch stack, the `children_extent_type` trait, and their tests. Nothing else is entangled. | no — it is the regridding-specific caching | ported to CR; recommend closing #474 |
| `claude/perf-ladder-predicates` | #476 | `as/intersection_area`'s 7 commits plus two more: `spherical_orient`'s inlined stable cross product (exposing `RobustCrossProduct.min_stable_norm`), and `SphericalCap._intersects` settled on a squared-chord bracket instead of `atan2`. Both are geometry-kernel math with derivations and adjudication against 256-bit arithmetic. | yes | stays in GO |
| `as/intersection_area` | #473 | `intersection_area` for `ConvexConvexSutherlandHodgman`, `OverlayNG` and `FosterHormannClipping`; a formal ring-sink interface; `_OverlayInput` typed instead of `Any`; closed-ring trimming made explicit. | yes | stays in GO |
| `as/dualdfs-balanced-descent` | #463 | Based on an *older* main (`f21f45a`). Its first two commits — carrying node extents through the descent, and `node_extent_is_expensive` — have since landed on main under other SHAs. What is left is `dual_depth_first_search_balanced`, the `node_size` trait, and the `Lockstep`/`Balanced` descent-policy structs. | yes — it is a traversal policy over a generic trait, with a generic default (extent area) | stays in GO |

**Nothing needed extracting.** The expectation going in was that #474's three commits are
all the caching feature and that #476/#473 are genuinely GO; both held on inspection. No
new GO branch was created and nothing was pushed to GO.

One collision to note for whoever lands these: #463 and #474 both rewrite the
"neither node is a leaf" arm of `dual_depth_first_search`, and #463 still carries main's
old `_child_extents`. Whichever lands second needs its merge resolved by hand.

## Recommendation for GO PR #474

**Close it without merging.** The caching now lives in ConservativeRegridding as
`cached_dual_depth_first_search`, which is where the trees with computed extents are, and
GO's generic search keeps the shape it was designed for.

Two follow-ups for the GO side that are *not* done here and want an explicit call:

1. GO main's generic search still allocates a vector per visited node pair whenever
   `node_extent_is_expensive` is true. By the same argument that sends #474 to CR, GO
   could drop that caching entirely and let the generic search be a plain loop — the
   trees GO ships never set the trait, and the ones that do now have CR's search. That
   would be a small revert on main, not a merge of #474.
2. If instead GO keeps its caching, #474's scratch stack is a strict improvement to it
   and could be merged on its own merits — but then CR's copy is redundant and CR should
   just call GO's. These two are alternatives; pick one.

Either way `node_extent_is_expensive` itself stays: CR reads it as the opt-in.

## DGG-side wiring that `claude/cr-pr138` needs

**None.** Concretely:

- **No trait methods to add.** Every DGG tree already declares
  `STI.node_extent_is_expensive`: `HierarchicalGridCursor`, `RasterCellTree`,
  `CopernicusDEM.BlockCursor` and `CapCachedTree` (which forwards, correctly — it caches
  only *leaf* caps, so its internal-node extents really are derived) say true;
  `CellCapTree`, `MemoRasterTree`, `MemoBlockCursor`, `DGGChunkTree`, `PositionTreeNode`,
  `RasterFlatTree`, `CapQuery` say false. CR reads that trait directly, so all of them
  get the right behaviour with no edit.
- **No `children_extent_type` methods.** Every DGG tree's extents are
  `GO.UnitSpherical.SphericalCap{Float64}` at every level, node and children alike, so
  the `nothing` default is right everywhere. The trait exists for trees whose extent type
  changes with depth; DGG has none.
- **No call-site changes.** DGG reaches the search through
  `ConservativeRegridding.intersection_areas`, whose signature is untouched. Verified by
  running DGG at `0d0f069` unmodified against the new CR.

Worth knowing, though not required:

- `MemoRasterTree`'s docstring describes exactly the problem this change solves — "the
  dual-tree join re-derives a node's extent once per opposing node" — and its
  direct-mapped extent memo is now redundant for `node_extent`. Its *leaf* memo (the
  `(position, cap)` entry vectors, i.e. `child_indices_extents`) is a different
  concern and stays needed. Same split applies to `CapCachedTree`: its leaf cap vector
  still earns its keep; the `_CACHED_BUCKET_SIZE = 49` tuning, whose docstring records
  "+25% at leaf 50 and +441% at leaf 350" against a bare cursor, is worth re-measuring
  now that side-2 extents are cached on a stack.

## Open questions

- Wall time did not move. The remaining cost in the discovery phase is `node_extent`
  itself and the leaf `child_indices_extents` allocations, not the vectors this removes.
  If discovery is worth more attention, the leaf entries are the next target.
- Whether GO should keep any caching in its generic search at all (see the two
  alternatives above) is a call for the GO side, not made here.
- CR's `QuadtreeCursor` has no no-argument `ncells`, so `Trees.split_weight` throws for
  it and the threaded frontier cannot take one. Pre-existing; unrelated to this work but
  it surfaced while writing the tests.

## Artefacts

- CR clone/branch: `/home/asinghvi17/geo/CR-cached-dfs`, branch `claude/cached-dual-dfs`.
- CR baseline worktree: `/home/asinghvi17/geo/CR-baseline` (detached at `3b50a2f`).
- DGG scratch worktree: `/home/asinghvi17/geo/DGG-crcache-scratch` (detached at
  `0d0f069`), with `bench_discovery.jl`, `bench_pipeline.jl`, `bench_out/{old,new}.jls`
  and `Project.toml{.oldcr,.newcr}` — all uncommitted and throwaway.
