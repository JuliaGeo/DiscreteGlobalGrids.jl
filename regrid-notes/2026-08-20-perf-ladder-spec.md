# Perf-ladder spec: source-side hot-path fixes for the CopDEM→IGeo7 production regrid

Written 2026-08-20 by the orchestrator, from the measured profile in
`regrid-notes/2026-08-20-production-profile.md` (+ `.ndjson`, 73 measurements).
Read that report first — every fix below cites its evidence.

## Context

- Live production dress rehearsal (PID in `/home/asinghvi17/geo/dggstores/copdem90-run.pid`)
  is running on worktree `/home/asinghvi17/geo/DGG-subzone-store`
  (branch `claude/copdem-production`). **Do not touch that worktree, the live
  store/log/conf under `/home/asinghvi17/geo/dggstores/`, or the run process.**
- Work happens in worktree **`/home/asinghvi17/geo/DGG-perf-ladder`**, branch
  **`claude/perf-ladder`**, based on `claude/copdem-production` @ `47961ea`.
  Root `Manifest.toml` already copied in with the CR git-tree-sha1 fix
  (`873cc732e77638eb44d7d92e06880585a31200a9` beside the repo-rev — do not
  regenerate it from scratch; `Pkg.instantiate` should just work).
  Bench env: `julia --project=bench` from the worktree root (workspace member;
  manifest at root). `ENV["RASTERDATASOURCES_PATH"]="/home/asinghvi17/geo/DiscreteGlobalGrids.jl/bench/data"`.
- Upstream pins: CR `6a4b997` (branch `claude/budget-frontier` of
  ConservativeRegridding.jl), GO `2825c47`. If a fix needs a CR-side change,
  branch the CR repo, push it, and update the pin SHA (precedent: the existing
  budget-frontier pin). Julia 1.12.6.
- **Box is shared** (64 cores, live run holds ~22, other tenants, load 47–68).
  Budget for this work: at most ~12 threads of benchmarking at a time, check
  `uptime` before heavy runs.
- Standing directives: principled library fixes, push branches, **never open PRs**.

## The measured baseline (must reproduce before changing anything)

One mid-latitude level-5 column (La=5, 823,543 L12 cells) costs ~28 core-seconds
at `-t 1`; ~28.3k cells/core-s. Hot path (flame-graph inclusive):
`node_extent` on the source `BlockCursor` 47.7% (of which `_box_cap` 41.9%),
`child_indices_extents` leaf caps 21.3%, `spherical_distance` 24.7%,
intersection/assembly arithmetic 24.7%. Alloc: 3.22 GiB/column, 71.6% of it
`Memory{Tuple{Int64, SphericalCap{Float64}}}` from `child_indices_extents`.
Reproduction recipe + experiment scripts table: report §7.

## Fixes, in order

### c1 — Memoize source-tree node extents per DEM tile  (target −26% CPU)
`STI.node_extent(::BlockCursor)` (`src/systems/CopernicusDEM/cursor.jl:245–248`)
recomputes `_node_box` → `_box_cap` (sincosd, UnitSphereFromGeographic,
spherical_distance) for every interior node the dual DFS visits, on every block
build, every column, every worker — for tile trees that never change.
Precedents to follow: destination side `CapCachedTree` (PR #44) and the raster
source's `MemoRasterTree` (`raster_tree_memo.jl`). Cache interior-node caps
keyed by source chunk; a tile's interior tree ≈ few thousand caps ≈ 100 KB;
26,475 tiles ≈ 2.6 GB total, so an LRU of the hot few hundred (workers process
8 contiguous columns per batch — strong tile locality). Must be safe under
concurrent readers (21 workers share the source DGGSpace).

### c2 — Stop per-leaf allocation in `child_indices_extents`  (target −9%)
The leaf-caps call allocates a fresh `Vector{Tuple{Int,Cap}}` per leaf — 71.6%
of all allocation. Use a per-task reusable buffer (STI accepts any indexable).
Find where this lives (STI interface implementation for BlockCursor — may be
DGG-side or CR-side; fix wherever it belongs, principled).

### c3 — Make the wave honest so CR's inner threading is reachable  (shape fix)
`_fillwave!` (lazy.jl:477) wraps every spawned block build in
`@with OUTER_PARALLEL => true`, suppressing CR's `_innerthreaded` /`_leafcaps`
threading — which measures **8.7× on 12 threads at 92% per-core efficiency**
(report §3.4) — whenever a column has >1 source chunk. Meanwhile the wave
itself delivers only 1.05–1.44× ideal (one block carries 44–100% of the work,
2–6 of the "connected" chunks yield zero pairs). Fix: `_wavesize`/`_fillwave!`
should estimate real work and NOT force OUTER_PARALLEL when the wave is narrow
(e.g. when one block dominates or nchunks is small), letting inner threading
take over. No CPU cut at fixed cores, but: W can drop from 21 to ~8 for the
same 24 cores, RSS ~3× lower, and it's the only path past 24 cores.
Pair with script change b2 (workers from a core budget, not thread count).
Careful: this file may be CR-side (the pinned repo) — if so, branch CR.

### c4 — Cache leaf caps per tile  (target another −10–15%)
Extension of c1 to the leaf level: 1.44M caps ≈ 46 MB/tile, small LRU, hot
within a batch of 8 contiguous columns. Do this only after c1/c2 land and
re-profiling confirms `child_indices_extents` is still the next sink (c2 may
change its cost profile).

### b2 — Script: `workers` from a core budget
In `scripts/copdem_production.jl`: compute `W` from a `--cores` budget
(`W = ceil(budget / measured cores-per-worker)`) instead of a thread count;
drop the `threadsper` fiction. After c3, cores-per-worker changes — recompute.

### Explicitly rejected (do not implement)
- c5 zero-pair chunk pruning: <0.5%, not worth it.
- Any change to `batch`/`cache`/`stripes`/`budget` defaults: tile path is
  0.02% of CPU; lowering `budget` splits the destination tile and is harmful.

## Acceptance criteria (every fix)

1. **Bit-identical output** on reference columns spanning regimes: mid-lat
   fed (e.g. columns 728–730 from the report), single-source (98241), polar
   multi-source (115424–115427), one 100%-NaN column. Compare written chunk
   bytes (or exact value equality incl. NaN positions) against the unmodified
   baseline built in the same session. c3 changes scheduling only — output
   must still be bit-identical.
2. **Measured speedup** per fix, same protocol as the report (`-t 1` for CPU
   cuts, `-t 12` W=1 for c3's shape), appended to a new
   `regrid-notes/2026-08-20-perf-ladder.ndjson` + summary md. State load at
   measurement time.
3. **Tests green**: full DGG suite once at the end (baseline: 987,456 pass /
   17 broken / 0 fail); GR/regridding suites after each fix. New caching gets
   its own tests (cache hit correctness, LRU eviction, concurrent access).
   Known: CI Julia 1.11-ubuntu allocation-law asserts fail on main already —
   not a signal. Allocation-law baselines may legitimately IMPROVE (c2);
   update those assertions with measured numbers, noting old → new.
4. **One commit per fix**, message explains the measured basis. Push
   `claude/perf-ladder` (and any CR branch) — **no PRs**.

## ADDENDUM (2026-08-20, after the improvement-hunter campaign)

`regrid-notes/2026-08-20-more-improvements.md` (+ 97-row ndjson) supersedes the
ordering above. **N1 lands first**: `PartialGrid.bucket_size` defaults to 0 →
the destination dual-DFS descends to 823,543 one-cell leaves; bucket_size=49
gives 2.35–2.71× per column, bit-identical, all regimes. All c1/c2/c4
measurements rebaseline against bucket_size=49 (their wins shrink; drop fixes
that no longer pay, recording the null result). Then N3 (isbits CopDEM
cell_boundary, −2.57%) and N2/N4 (GeometryOps spherical_orient strip +
normalization dedup, −6%+, needs GO branch off pin 2825c47 + repin). After c4,
re-sweep BlockCursor LEAF_CELLS=9. Hunter-measured dead ends (do not revisit):
leaf/leaf pre-pruning, per-tile knownempty, cell-polygon-per-tile memo, F3/D1
bisection, _fillcoo! fusion, batch ordering, Float32 clipper. c2's real win is
−3.4% wall, not −9%.

## Deliverables

- Pushed branch(es), commits c1..c4 + b2 (c4 optional if profiling says no).
- `regrid-notes/2026-08-20-perf-ladder.md` + `.ndjson` in the MAIN checkout
  (`/home/asinghvi17/geo/DiscreteGlobalGrids.jl/regrid-notes/`): per-fix
  measurements, final projected production ETA at 24 cores, any deviations
  from this spec with reasons.
- A short section: what the NEXT profile shows as the new top-3 sinks after
  all fixes (re-run the §7 hotpath experiment once at the end).
