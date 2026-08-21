# 2026-08-19 — `_child_extents`: what it is, why it churns, and what the fix buys

The largest single allocation site left in the eager MOC regrid lives upstream,
in GeometryOps' dual-tree search. This note explains it end to end and records
the measured effect of fixing it. Upstream draft PR: <https://github.com/JuliaGeo/GeometryOps.jl/pull/474>.

## What the function does

`GO.SpatialTreeInterface.dual_depth_first_search` descends two spatial trees in
lockstep. In the branch where **neither** node is a leaf it runs a double loop —
every child of `node1` against every child of `node2` — and needs each child's
extent to test the predicate. Nothing stores those extents: a tree can opt in to
`node_extent_is_expensive(::Type{MyNode}) = true` to say "my `node_extent` is
*computed*, not read". For such a tree, deriving `node2`'s child extents inside
the inner loop would redo the same work once per child of `node1`, so GO hoists
them:

```julia
@inline function _child_extents(node::N) where {N}
    node_extent_is_expensive(N) || return nothing
    return [node_extent(child) for child in getchild(node)]   # ← one Vector per visited node pair
end
```

That is the whole function: a fresh `Vector` of the second node's child extents,
built once per *(node1, node2)* pair visited, thrown away when the pair is done.
The caching is a real win — it turns `nchild1 × nchild2` extent derivations into
`nchild1 + nchild2` — but it pays for it in heap traffic.

## How DGG gets there

```
DGG.regrid(dem; to = grid)
  → GlobalRegridding conservative operator: src RasterCellTree, dst CapCachedTree
    → CR.intersection_areas → multithreaded_dual_query
      → GO.SpatialTreeInterface.dual_depth_first_search   (one call per frontier pair, per task)
        → _child_extents(node2)                            ← here
```

Both DGG trees opt in. `RasterCellTree` (`rastergrid.jl:820`) derives a
`SphericalCap` for a pixel-index rectangle; `CapCachedTree` (`src/cap_cached_tree.jl:20`)
inherits the trait from `HierarchicalGridCursor` — its *leaf* caps are cached in
a vector, but an internal node's cap is still derived from the cursor. So every
internal-vs-internal node pair in the descent allocates one vector.

## Why that is gigabytes

The target workload — 3600² GLO-30 tile onto a level-13 IGEO7 MOC of 16,181,892
cells — generates **93,079,031 candidate pairs**, and IGEO7's fanout of 7 means
the descent spends most of its life in the both-internal branch. The allocation
profile (`regrid-notes/2026-08-19-post-stack-profile.md` §4) put
`dual_depth_first_search.jl:33` at **~7 GB, 42 % of candidate-pair-generation
bytes** — the biggest single site in the whole regrid.

That 7 GB was a *sampled share*, and it was slightly overstated: `Profile.Allocs`
folded DGG's own `child_indices_extents` buffer
(`Memory{Tuple{Int,SphericalCap}}`, `cap_cached_tree.jl:43`) into the same site
through inlining. Removing `_child_extents`' vector alone, measured
deterministically end to end, takes **38.58 → 33.61 GiB (−4.97 GiB, −12.9 %)**.
DGG's own leaf-entry buffer is the remainder and is a separate, DGG-side fix.

## The fix

Keep the caching, drop the per-pair vector: thread **one scratch stack** through
the descent. Each visited node writes its children's extents onto the end of that
stack, indexes them by offset for the duration of its double loop, and truncates
back on the way out. The stack is created lazily at the first node that needs one
and reused for the entire recursion, so the extent derivations are exactly the
same as before and the allocation count drops from *one per node pair* to *one
per traversal*. CR's multithreaded search enters `dual_depth_first_search` once
per task, so each task gets its own stack and no synchronization is involved.

A stack has one element type, so siblings must agree on their extent type. A new
opt-in trait, `children_extent_type`, says what that type is where it is not
simply the parent's; a level whose children type differs from the stack it
inherited starts its own. The traversal asks the node — `children_extent_type(node)`
— and a node method falls back to a type method, so a tree may answer either way
and both fold. DGG's trees are homogeneous and define nothing.

The public 4-arg and 6-arg entry points are unchanged; the stack rides on a new
7-arg internal method.

## What it buys, honestly

| warm rep, mean of 6 | GO main | with the fix | delta |
|:--|--:|--:|--:|
| wall | 52.07 s | **50.37 s** | −1.70 s (−3.3 %) |
| process CPU | 307.5 s | **303.4 s** | −4.1 CPU-s |
| GC | 5.33 s | **4.44 s** | −0.89 s (−17 %) |
| allocated | 38.58 GiB | **33.61 GiB** | −4.97 GiB (−12.9 %) |
| digest | `5b696475a3665634` | `5b696475a3665634` | identical |

Three A/B rounds, base and fix alternating in one window, 1-minute load 21–39
throughout, three reps each (reps 2–3 warm, all six above). The wall ranges do
not overlap — 51.20–53.58 s against 49.73–51.06 s.

**A trap worth remembering.** An intermediate version sized the fresh buffer with
`nchild(node2)` on every visited node and cost **+7 CPU-s** end to end:
`STI.nchild` on a `WindowCursor` loops over the child windows rather than reading
a stored count (`src/engine/cursor.jl:201`). Anything called once per node pair in
that descent is on a 27 M-times-per-regrid path. The reuse arm is now selected by
dispatch, so `nchild` is only reached on levels that genuinely allocate.

Isolated microbenchmark (fanout-7 k-d tree, two matched-resolution trees of
63,001 and 62,500 leaves, 250,000 pairs, 30,399 internal-vs-internal node pairs,
`node_extent_is_expensive = true`): **4.17 MiB in 30,399 allocs → 5.87 KiB in
7 allocs**, 21.99 → 21.56 ms min (24.42 → 21.67 ms median). The same tree shape
with the trait off is unchanged at 16 bytes / 1 alloc.

**The lever paid exactly what the profile predicted.** §6 lever 5 of the
post-stack note ceilinged it at −1.3 s, bounded by candidate-pair generation's
whole GC deficit; the measured −1.7 s wall is that ceiling plus the mutator work
the old code spent allocating and initialising 27 M short-lived vectors, which
shows up as −4.1 CPU-s beside the −0.9 s of GC. It stays a small lever: 3.3 % of
the call, against the two big ones (fusing `_fillcoo!` with the WeightBlock
`sparse()`, and destination-geometry reuse, worth −5.5 s and −6.4 s). Note the
figures carry shared-box noise of roughly ±1 s on any single rep — the claim
rests on the A/B being paired and on the two wall ranges not overlapping.

## The bump-allocator alternative, measured

The review proposed backing the per-node buffer with a Bumper.jl bump allocator
instead, so heterogeneous levels would cost nothing. I prototyped it (one
`SlabBuffer` per traversal, `alloc!` + checkpoint per node) and measured both:

| tree shape | GO main | typed stack | Bumper |
|:--|:--|:--|:--|
| homogeneous | 21.99 ms · 30,399 allocs | **21.6–21.7 ms · 6** | 21.66 ms · 1 |
| type changes every level | 18.24 ms · 30,399 | 18.75 ms · 96,629 | **17.62 ms · 1** |

A dead tie on the shape that exists; Bumper wins only on the shape that does not.
Against that it wants a new GO dependency, and its `UnsafeArray` cannot hold a
`Union` or non-isbits extent type at all (verified: `MethodError`) — a capability
both GO main and the stack have. `@no_escape` also rejects `return` anywhere in
its body, and macroexpands nested macros to find it, so `@controlflow` cannot
live inside the block; it needs a function barrier that changes how
`Action(:return, x)` propagates. **Recommendation: keep the typed stack**, with
the Bumper numbers in the PR so the GO maintainer can overrule with data.

## Risk

- **Extent-type homogeneity, now only among siblings.** One stack means one
  element type. The first version of this fix required it tree-wide; the review
  asked for that to be relaxed, so the child type is now a trait
  (`children_extent_type`) and a level that disagrees with its inherited stack
  starts its own. Cost of that generality: on a tree whose extent type changes at
  every level the traversal allocates ~3× the objects GO main does (96,629 vs
  30,399 on the isolated benchmark) and runs ~3 % slower — correct, but not free.
  No tree in GO, GR or DGG is heterogeneous, so DGG pays none of it.
- **Aliasing.** Correctness rests on push/truncate discipline: every frame records
  the stack length on entry, reads only what it pushed above that mark, and
  truncates back to it on the way out. A callee can only append above the
  caller's mark, so an early exit that skips a truncate (`Action(:return, …)`, or
  `Action(:full_return, …)` unwinding the whole traversal) leaves the caller's
  region intact and is cleaned up by the caller's own truncate.
- **Trees that do not opt in are untouched** — the trait folds at compile time,
  the `nothing` arm is the old code, and GO's "default trait allocates nothing"
  test still passes at exactly 0 bytes.
