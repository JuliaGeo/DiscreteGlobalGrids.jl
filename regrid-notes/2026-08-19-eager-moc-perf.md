# 2026-08-19 — Eager MOC regrid campaign (target: the maintainer's exact workload)

Target (verbatim): eager `DGG.regrid(dem; to = grid)` — Copernicus GLO-30 tile
N46_00_E010_00 (3600², in memory) → level-13 IGEO7 MultiOrderCoverage covering
of (10..11°E, 46..47°N), `PartialGrid(CellLookup(region))`, Conservative +
Weighted(0.5), `-t 8`, warm. ~76 s on the maintainer's machine at campaign start.

Environment: shared 64-core box (load noted per row), Julia 1.12.6,
DGG worktree pinned per item (baseline origin/main 0bc2b02), CR pinned
`claude/budget-frontier` 6a4b997, GeometryOps main c30ec910 (v0.1.44).
Raw rows in `2026-08-19-eager-moc-perf.ndjson` beside this file.
PRs are STACKED: item N+1 branches from item N's branch; bases retargeted to
main as the stack merges bottom-up.

## Baseline (worktree @ 0bc2b02, CR 6a4b997, GO c30ec910)

Warm reps, one `julia -t 8` session, `GC.gc()` between reps; digest = (sum,
sha256[1:16]) of the Float64-widened output vector (missing → sentinel).

| rep | wall | GC | churn | load | maxrss |
|--:|--:|--:|--:|:--|--:|
| 1 (first call) | 137.2 s | 17.7 s | 110.2 GiB | 31.7→28.5 | 8.4 GiB |
| 2 | 124.4 s | 16.3 s | 109.5 GiB | 28.5→26.3 | 10.0 GiB |
| 3 | 122.8 s | 16.6 s | 109.5 GiB | 26.3→35.8 | 10.0 GiB |

digest `sum=3.368822268175383e10  sha=5b696475a3665634`, n=16,181,892.

Stage decomposition (t8, load 30–42): candidate pairs 65.4 s / 15.5 GiB
(93,079,031 pairs); assemble (clip loop) 52.5 s / **83.2 GiB** (66,129,312
nnz); `sparse` 2.0 s; `_fillcoo!`+`WeightBlock` 5.4 s; apply 0.4 s.
Profile (t8, delay 10 ms, 176k samples, thread-summed): GC-flagged 51.3%,
scheduler idle 13.4%, IGEO7 dst geometry 13.2%, raster src geometry 9.3%,
SH clip 5.6%, descent 2.6%, dst tree walk 2.0%, area 1.7%, sparse 0.4%.
Inclusive: `cell_boundary` 10.5%, `node_extent` 7.5%, `_cornercap` 7.8%,
`spherical_distance` 7.5%, `cell_polygon` 5.2%, `index_to_cell` 2.7%,
`_child_window` 2.8%. So the ranked levers: (1) assemble churn → GC, (2) dst
per-visit geometry, (3) src per-visit caps, (4) polygon synthesis per pair.

## Per-item verdicts (stack order, each measured against the tip below it)

| # | branch / PR | what | warm wall (t8) | churn | digest |
|--:|---|---|---:|---:|---|
| — | baseline 0bc2b02 | | 122.8–124.4 s | 109.5 GiB | 5b696475a3665634 |
| 1 | claude/moc-dst-tree-cache #47 | whole-space dst tree: ids decoded once, leaf caps precomputed (CapCachedTree) | **103.0–103.4 s** (−16%) | 108.8 GiB | identical |
| 2 | claude/moc-area-only-clip #48 | spherical clip integrates ring in cache; no result polygon / area collect | **84.6–84.9 s** (−18%) | 52.3 GiB | identical |
| 3 | claude/moc-frozen-src-tree #49 | FrozenRasterTree: src node extents + leaf caps stored once | **63.2–64.8 s** (−25%) | 53.7 GiB (rss +2.4 GiB) | identical |
| 4 | claude/moc-clip-cell-reuse #50 | 64-slot per-task getcell memos (both sides) | **46.1–48.0 s** (−26%, load 36–40) | 38.7 GiB | identical |
| 5 | claude/moc-wholeblock-adopt | eager Conservative wholeblock adopts CR block (no COO round trip) | 53.5–54.6 s @ load 28–32 — SLOWER than item-4 row; churn +0.5 GiB; re-measuring interleaved | 39.2 GiB | identical |
| 6 | claude/moc-frozen-dst-tree | FrozenDGGTree: interior dst extents stored once, allocation-free leaves | (measuring, interleaved with 4 and 5) | | |

Items 4–6 measured across sessions on a loaded box; the interleaved re-run
(6→4→5 back-to-back) decides item 5's verdict and item 6's number.


## Interleaved re-measure (same window, back-to-back sessions)

| tip | warm wall | churn | load | note |
|---|---:|---:|:--|---|
| item 6 (with 5) | 54.4–55.4 s | 40.7 GiB | 22–27 | briefly overlapped a lib-suite rerun |
| item 4 | 48.6–50.1 s | 38.7 GiB | 30 | consistent with its first 46.1–48.0 window |
| item 5 | (aborted; wholeblock A/B below decides) | | | |

Items 5+6 measure as a net ~5 s regression against item 4 in the same window.

**Item 5 verdict: DROPPED.** Same-session alternating A/B of the plan build
(fast `wholeblock` vs the generic COO round trip, two rounds, load 25–27):
fast 52.5 / 53.1 s vs generic 46.0 / 46.6 s, fast +0.47 GiB churn — the
"remove the COO round trip" estimate not only fails to reproduce, the fast
path is reproducibly slower. Mechanism unresolved (post-stage should be
strictly less work; suspect allocator/locality interaction inside the
otherwise-identical assemble). Branch deleted; item 6 rebased onto item 4.

**Item 6 verdict: DROPPED.** Isolated on top of item 4 (commit 7dcf7c2):
47.9–50.5 s at load 32–37 vs item 4's 46.1–50.1 in adjacent windows — within
noise (<3%), with +1.5 GiB churn and +0.3 GiB rss from freezing 2.8M interior
nodes. GO's per-visit child-extent caching already amortizes interior extents
well enough that storing them buys nothing measurable. Branch deleted.

**Final stack: items 1–4** (#47 merged, #48, #49, #50).

## Review round (2026-08-19, asinghvi17) and re-audit

Three review comments, all addressed with commits folded into #48 (the open
stack's bottom) plus replies in-thread:

| PR | comment | disposition |
|---|---|---|
| #47 (merged) | "What is `_Cap`?? Better to give this a type param" | 07ad6cd: alias deleted; `DGGSpace`, `DGGChunkTree`, `CapCachedTree` carry the cap type as an inferred parameter |
| #47 (merged) | CapCachedTree "should be in a separate file" | 559ca75: moved verbatim to `src/cap_cached_tree.jl` |
| #48 | "new file please" (clip kernel) | 2530dec: moved verbatim to `lib/GlobalRegridding/src/spherical_clip_area.jl` |
| #49 | review: "wrong approach ... not viable in the presence of larger than memory data" | REWORKED, see below |

**Item 3 rework.** The frozen whole-rectangle tree materialized ~40 B/cell
(+2.4 GiB rss on the target); on the chunked out-of-core path `subtree` only
ever sees one source chunk, but the eager path froze the whole raster.
Replaced (fed9421) by `MemoRasterTree`: the lazy `RasterCellTree` wrapped with
a per-task direct-mapped memo (1024 extent slots + 256 leaf-entry slots,
~150 KB/task, keyed by rectangle + raster identity). Byte-identical caps;
memory O(slots) regardless of raster size.

Re-audit window (2026-08-19 15:29–16:10 UTC, serialized, digests all
`5b696475a3665634`):

| tip | warm wall (t8) | churn | maxrss | load |
|---|---:|---:|---:|:--|
| main tip 96d7f5d (=#47) | 101.8–110.4 s | 108.8 GiB | 11.0 GiB | 25–54 |
| #48 @ 07ad6cd (+review commits) | 85.7–87.2 s | 52.3 GiB | 10.2 GiB | 23–27 |
| (frozen item 3, dropped) | 61.5–62.4 s | 53.7 GiB | 12.4 GiB | 22–23 |
| (#50 on frozen, superseded) | 45.9–46.8 s | 38.7 GiB | 14.0 GiB | 21–23 |
| #49 @ fed9421 (memo rework) | 65.0–66.7 s | 53.6 GiB | 10.5 GiB | 20–35 |
| #50 @ abaf1d6 (on memo) | 50.6–50.8 s | 38.6 GiB | 11.0 GiB | 23–33 |

The memo keeps ~85% of the frozen win (−24% vs −28% on the stack) at a
quarter of the extra rss here and O(1) in general. Final stack: #48 → #49
(memo) → #50, cumulative 101.8→50.6 s in-window (~2×), and ~2.4× against the
campaign's 0bc2b02 baseline.

CI note: `Julia 1.11 - ubuntu` fails identically on main 96d7f5d and on #48
(same 8 failures, 987154/8/17 — H3/A5/IGeo7 allocation+neighbor tests,
`neighborhood.jl:152`, `subtree_halos.jl:1592`); pre-existing, not
stack-caused. 1.12 linux+mac and conformance green.

## Merges (2026-08-19, maintainer approved the memo rework: "Yes, free to merge.")

Bottom-up, each gated at a tree byte-identical to its post-merge main tree
(every merge was a content fast-forward; final `git diff main abaf1d6` empty):

| PR | head | gates (root / lib) | merge commit |
|---|---|---|---|
| #48 spherical clip in place (+3 review commits) | 07ad6cd | 987162/17/0 / 2170/1/0 | 0588703 |
| #49 per-task raster tree memos (rework) | fed9421 | 987162/17/0 / 2184/1/0 | 0ccab2a |
| #50 per-task cell-polygon memos | abaf1d6 | 987162/17/0 / 2190/1/0 | 1bb575d |

Root tally identical to the post-#47 baseline throughout; lib growth (486/1
pre-campaign → 2190/1) is the campaign's own testsets. The stray red CI job
(`Julia 1.11 - ubuntu`) fails identically on pre-campaign main (8 failures,
987154/8/17 — H3/A5/IGeo7 allocation+neighbor tests); the #48-run
Documentation job hung on a starved runner (main's identical-content build
passed in 8.5 min).

## Campaign close

Target workload, warm t8, digest `5b696475a3665634` (byte-identical
everywhere, n=16,181,892): campaign baseline 0bc2b02 122.8–124.4 s /
109.5 GiB churn → final main 1bb575d 50.6–50.8 s / 38.6 GiB in-window
(~2.4× wall, 2.8× churn). Against the post-#47 re-baseline in the same audit
window (101.8–110.4 s), the three merged PRs cut ~2.1×.

Remaining unimplemented levers (from the baseline profile, post-stack):
candidate-pair generation still dominates (dst-side IGEO7 geometry + tree
walk); scheduler idle ~13% at t8; items 5 (wholeblock adopt) and 6 (frozen
dst interior extents) measured as regressions/no-ops and were dropped —
records above. The `_chunklevel`/`_chunkwindows` defaults bug (perf-P3 D3/D4)
remains open.
