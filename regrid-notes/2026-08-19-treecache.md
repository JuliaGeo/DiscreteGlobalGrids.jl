# 2026-08-19 — TileCells restricted-tree memoization (PR #44)

Fixes lever 3 of `2026-08-19-thread-utilization.md` §6: budget-derived
destination tiles are not chunk-aligned, so every block build rebuilt a fresh
`CellCapTree` over the whole ~3.3M-cell tile — 320 rebuilds on the 3600² lazy
default path. `TileCells` now memoizes the tile's restricted tree under its
existing lock; one tree per tile serves all of its blocks.
`GR.cellcaptree_builds()` is the deterministic counter seam the tests use.

Merged as `af37e7a` (PR #44, branch from `9fc8514`). Measured on the shared
box at t8, CR pinned `508a637` (claude/budget-frontier), GeometryOps main
v0.1.44, Julia 1.12.6; raw rows beside this file in
`2026-08-19-treecache.ndjson` (before = `label:"base"`, after = `"cached"`).

MOC workload (whole GLO-30 tile N45_00_E010_00, lazy, 512×512 source chunks,
level-13 IGEO7 covering, defaults):

| metric | before (9fc8514+counter) | after (cache) |
|:--|--:|--:|
| read wall | 860.3 s (load ~22) | **179.3 s** (load ~43) |
| read util | 4.23 | 4.86 |
| `CellCapTree` builds | 320 | **5** (one per tile) |
| peak RSS (Sys.maxrss / time -v) | 5546 / 5679 MB | **3786 / 3877 MB** |
| digest / output hash | 6.904662615357636e9 / 0x05e3047a1370e472 | identical |

≥4.8× wall under worse load; the default path now matches the postfix
campaign's best chunk-aligned configuration (140–175 s). RSS drops ~1.8 GB —
a wave of 4 builds no longer holds 4 private tree copies. Small lazy cases
(single whole-space tile, zero builds either way): l512 warm 2.92→3.09 s,
l1024 warm 12.73→12.75 s, hashes identical — load noise only.

`_wavesize` left unchanged: its weight-budget bound models the per-block
sparse-matrix memory floor, which never charged for tree builds and is
unchanged by the cache. The `_chunklevel`/`_chunkwindows` defaults bug
(perf-P3 D3/D4, recommendations C2) remains open and separate.

Suites at the merge: root 987152 pass / 17 broken / 0 fail (unchanged);
lib/GlobalRegridding 475 / 1 / 0 (was 460 / 1 / 0; +15 = the new tests).
CI: all green except the pre-existing Julia 1.11 `@allocated` flake family
(identical 8 sites as PRs #40–#43).
