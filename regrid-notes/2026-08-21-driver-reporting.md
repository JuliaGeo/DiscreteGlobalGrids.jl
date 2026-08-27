# CopDEM driver: session vs store in the progress reporting

2026-08-21, branch `claude/driver-reporting` off `claude/perf-ladder` @ `0d0f069`.

Two reporting bugs in `scripts/copdem_production.jl`, both found in the completed
GLO-90 production run. Neither touched the data — every cell in the store was and
is correct. Both made a resumed run misreport itself, and both were believed.

## 1. The heartbeat ETA credited skipped chunks to this session

`heartbeat!` divided the session's elapsed time by `done + skipped`. A chunk
skipped at resume costs a set lookup, not a regrid, so folding the skips into the
session clock makes the rate arbitrarily fast. Session 3 of the production run
printed **"ETA 0.30 h"** with about three hours left: 54 k skips and 4.5 k real
chunks divided into one session's stopwatch as though all 58 k had been computed
in it.

`etaseconds(elapsed, computed, remaining)` now takes `computed` — chunks this
session actually regridded — and the heartbeat names both counts and labels every
rate `session`:

```
HEARTBEAT  162/400 chunks (3 computed this session, 159 skipped) | session 1.029e+03 cells, 132 cells/s | elapsed 0.00 h | ETA 0.17 h | ...
```

Measured on the kill/resume test below, against the old formula recomputed from
the same numbers. True remaining time from the 14:51:39 heartbeat was 0.024 h:

| heartbeat | chunks | this session | old ETA | new ETA |
|---|---|---|---|---|
| 14:51:18 | 162/400 | 3 | 0.0032 h | 0.17 h |
| 14:51:39 | 205/400 | 46 | 0.0077 h | 0.03 h |
| 14:52:26 | 305/400 | 146 | 0.0066 h | 0.01 h |

The old ETA is 3-4x optimistic through the whole run and never converges, because
its denominator is dominated by work it did not do. The new one is honest from the
first chunk and tightens as the sample grows.

## 2. The final banner was session-scoped without saying so

`RUN DONE` reported the `Progress` counters, which only ever describe the process
printing them. For a resumed run that reads as the whole dataset: session 3 signed
off with `9.7606e9 cells ... NaN 44.308%` when the store held 5.4500e10 cells at
26.892% NaN (see `2026-08-21-run-verification.md` §4).

The banner is now two lines, both labelled, and the store-wide one is read from
the same ledger the resume set comes from, so the two can never disagree:

```
RUN DONE  session: 241 chunks computed, 159 skipped, 8.2663e+04 cells in 0.03 h = 709 cells/s aggregate
RUN DONE  session: NaN 43.847% of the 8.2663e+04 cells this session wrote
RUN DONE  store total: 400 of 400 chunks written, 1.3720e+05 cells, NaN 35.729%
```

Two new functions in `copdem_store.jl` carry it:

- `doneledger(logpath)` — one entry per **chunk**, not per line. A chunk logged
  twice (recomputed because a crash lost its file but not its line) is one chunk,
  at its last line's values; counting lines would double its cells. `donechunks`
  now reads this instead of its own regex, so "which chunks are done" and "how
  many cells that is" come from one parse.
- `storetotals(logpath, written)` — cells and NaNs over the resume union.
  Chunks the union has from the file listing but the ledger has no line for come
  back as `unaccounted` rather than as zero, and the banner says so, so totals
  read off a truncated ledger are never mistaken for the whole store.

`donechunks` gained a `label` kwarg only so the post-run call does not print its
notes under the word "resume".

## The resume banner needed nothing

`resume: N chunks already written, M of T to do` already prints `length(done)`
from the script's own `donechunks` call. The misleading "~18760" in the production
log was a hand-written note in the transcript, not script output. Left alone.

## Testing

Kill/resume against a scratch store — level 8 over level-5 chunks, 400 chunks,
3 workers, synthetic elevations over the real GLO-90 tile list:

1. Session 1 hard-killed (`timeout -s KILL`, exit 137) after 159 chunks.
2. Session 2 resumed, skipped those 159, computed the remaining 241, finished.

Store totals in the banner (400 chunks / 1.3720e5 cells / 35.729% NaN) were
confirmed independently by summing the ndjson ledger outside Julia, and are
distinct from the session's (241 / 8.2663e4 / 43.847%) as they should be.

Six `reportchecks` assertions run on every `dryrun = true`, in the driver's
existing `check` style — two on the ETA arithmetic, two on synthetic ledgers
(duplicate line, unaccounted chunk), two against the run's own ledger. Dry run on
the populated store: **8 PASS, 0 FAIL** (6 new plus the 2 pre-existing).
