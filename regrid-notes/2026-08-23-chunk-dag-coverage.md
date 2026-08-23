# Predictive completeness of `chunk_dependency_graph`: mechanism and fix

Follow-on to `2026-08-23-full-run-attribution.md` §"Predicted-edge failure". That
report established the symptom — **[measured] 71** executor demand pairs absent
from the graph over the full CopDEM90 × IGeo7-L12 problem — and named the cause
only in general terms. This note pins the mechanism to a line of code, proposes
and implements a fix, and re-measures the same reconstruction at **[measured] 0**
missing pairs. Branch `claude/chunk-dag-coverage`, cut from `eb1d638`. Nothing
pushed.

## 1. Mechanism, proven

The report's stated cause — "a mismatch between the graph's direct chunk-cap
relation and the DGG hierarchical candidate relation" — is **correct in outline
but too weak to act on**. It is not a hierarchy, ancestor-rounding, pentagon or
dateline effect. It is this:

**The same level-0 Copernicus tile has two different covering caps, computed by
two different functions, and neither contains the other.**

| side | cap | built by |
|---|---|---|
| graph | `chunkextents(DGGSpace) = space.caps` | `DGG.node_extent(sys, id)` — centroid of `cell_box`, radius = max distance to the `cell_boundary` ring, pad `Δλ²/16` |
| executor | leaf entry of the level-0 frontier tree, via `_mappedfrontierchunks!` | `STI.node_extent(::BlockCursor)` → `_box_cap(_node_box(c), _leaf_pad(c))` — centre of the *node* box, radius = max distance to its four corners |

`_node_box` (`src/systems/CopernicusDEM/cursor.jl:78`) returns, for a tile block,

```julia
return (_lon_w(c.q0) - widest, Float64(_lon_w(c.q1) + 1), south, north)
```

The west edge carries the half-pixel band offset; **the east edge does not**. For
a single tile the correct east edge is `_lon_w(q1) + 1 - half_dlon` (that is
exactly what `cell_box` returns). Measured, for four of the 13 runtime
reproducers:

| tile | `ncols` | `ΔW` | `ΔE` | `ΔS` | `ΔN` |
|---|---:|---:|---:|---:|---:|
| N81 W062 | 240 | 0.0 | **+0.00208333** | 0.0 | 0.0 |
| S87 E043 | 120 | 0.0 | **+0.00416667** | 0.0 | 0.0 |
| N29 E121 | 1200 | 0.0 | **+0.00041667** | 0.0 | 0.0 |
| N18 W113 | 1200 | 0.0 | **+0.00041667** | 0.0 | 0.0 |

`ΔE` is exactly `half_dlon = (1/ncols)/2`. The pads are bit-identical. So the
executor's tile cap has its centre shifted a quarter-pixel east and a slightly
larger radius; the graph's cap reaches further west. **Neither is a superset of
the other**, verified directly:

```
tile N81 W062  ordinal=2998
  centre separation = 2.687120166430015e-6 rad;  r_exec - r_graph = 4.1598554865324155e-7
  exec cap contained in graph cap? false
  graph cap contained in exec cap? false
  dst root 11019: graph-side intersects = false ; executor-side intersects = true
```

That last line is the defect, reproduced deterministically. The same script
reproduces it for 5 of the 13 runtime source indices tried (147 → dst root
11019; 25259 → 123179; 23465 → 121539; 11497 → 77487; 13103 → 90638) — **5/5**.
No data, no threads, no scheduler: two caps and one predicate.

### It is not CopDEM-specific

The same class of defect exists on the shipped `RasterGrid` source path, which no
CopDEM fix would touch. `_rastercandidates!` answers from quadtree node extents
and, at a leaf that straddles a chunk boundary, pushes *every* chunk in the
leaf's rectangle with no per-chunk test at all. Measured on a 360×180 global
raster with 37×23 chunks and 2,485 destination caps: **[measured] 44** of
**[measured] 5,630** executor pairs are absent from the `chunkextents` cap
relation.

So the defect is structural: **`chunkextents(src_space)` is not the cap set
`chunkindex(src_space)` tests against, and nothing required it to be.** Both
relations are conservative supersets of the true overlap, but they are
conservative in *different directions*, so they cross. A graph that crosses the
executor's relation cannot be a refcount.

## 2. Options

| option | correctness | build cost (66,228 × 26,475, t8) | edges | notes |
|---|---|---:|---:|---|
| (a) build the graph from `candidatechunks!` on `chunkindex(src_space)` | correct by construction, every space, no per-space invariant | **[measured] 0.023 s → 0.121 s** | **[measured] 328,424 → 328,058** | one executor-equivalent query per destination chunk; loses the cap join and its latitude prefilter |
| (b) inflate the graph's caps to dominate the index | not provable — the required inflation is per-space and undocumented; an ε is a fudge | ~unchanged | grows | rejected |
| (c) fix `_node_box`'s east edge / return the canonical cell cap at level-0 leaves | fixes CopDEM only; leaves `RasterGrid` broken; still relies on two constructions agreeing bit-for-bit at the poles | unchanged | slightly fewer | good hygiene, not a fix |
| (d) executor-side re-pin of a retired source | masks the symptom, does not restore the refcount, and the reload it is meant to avoid has already happened | n/a | n/a | rejected |

**Recommendation: (a).** It costs **[measured] +0.098 s** on the production pair
— against an **[measured] 8.81 h** run — and it *removes* 366 spurious edges,
so refcounts get tighter as well as sound. (c) remains a worthwhile separate
tidy-up (`_node_box`'s east edge is genuinely loose for a single tile), but it is
not the fix and is not in this patch.

Memory: unchanged on the DGG path (`chunkindex(::DGGSpace) === space`). For
spaces with no specialised index, `chunkindex` now materialises one packed RTree
over `nsrc` caps inside the builder; that is the same object the lazy executor
builds anyway.

## 3. The patch

Four files, `lib/GlobalRegridding/src/chunkgraph.jl` being the only source
change:

- `chunk_dependency_graph` now calls `_chunkgraph(chunkextents(dst_space),
  chunkindex(src_space), nchunks(src_space), radius, refine)`.
- `_fillrow!` is one `candidatechunks!` call plus the `refine` filter. The
  latitude prefilter, `_caplatitude`, and the sorted-latitude machinery are gone.
- `prefilter` stays as an accepted, documented **no-op** so callers and the
  recorded POC script keep working.
- Docstrings restate the invariant as a correctness property and cite the
  measured 71/437 split.
- `test/scripts/copdem_policy.jl` and `scripts/chunk_dag_poc.jl` follow the
  private signature change; the POC's now-meaningless prefilter A/B is replaced
  by the demand-domination check on the real pair.

## 4. Proof

Two tests, each failing before and passing after.

- `lib/GlobalRegridding/test/test_chunkgraph.jl`, "the graph dominates what a
  lazy read demands": a `RasterGrid` source whose chunks do not align with the
  quadtree's splits. **[measured] 2 failures before, 0 after.**
- `test/systems/crosssystem/regrid.jl`, "the dependency graph holds every pair
  the chunk index answers": the real pair in miniature — CopDEM-90 level-0
  frontier (**[measured] 64,800** tiles) × a complete IGeo7 level-3 destination
  chunked at level 2 (**[measured] 492** chunks), asserted over *every*
  destination, with all twelve pentagons and both poles present and every source
  tile reached. **[measured] 6 unpredicted pairs before, 0 after**, in
  **[measured] 1.0 s**.

Suites: GlobalRegridding **[measured] 3708 pass / 1 broken / 0 fail** (the broken
one is pre-existing). DGG `regridding_conservation`, `regrid`,
`regrid_acceptance`, `scripts/copdem_source_mode`, `scripts/copdem_policy`: all
pass.

Full-scale reconstruction, the identical read-only script that produced the
report's numbers, re-run on the patched tree:

```
before: chunks=66228 tiles=26475 graph_edges=328424 actual_edges=328058 missing_pairs=71  extra_pairs=437
after:  chunks=66228 tiles=26475 graph_edges=328058 actual_edges=328058 missing_pairs=0   extra_pairs=0
```

**[measured] 0** missing pairs, and the two relations are now identical rather
than merely nested.

## 5. Residual risk

1. **Raster-source build cost is not measured at production scale.** On a
   4320×2160 raster with 162 chunks the index query costs **[measured] 35.1
   µs** per destination cap against **[measured] ~0.6 µs** for the brute cap
   join — a ~58× per-query ratio, though it also returns **[measured] 8,953**
   pairs where the cap join returns **[measured] 15,788**, so the graph gets
   much tighter. Extrapolated to 66,228 destinations that is ~2.3 s serial,
   ~0.3 s at t8; a far larger raster with a deeper cursor was not measured.
2. **`prefilter` is now inert.** Any caller that read a speedup out of it gets
   1.0×. The recorded `chunk_dag_poc` ndjson keys changed accordingly
   (`bruteforce_seconds`/`prefilter_speedup` → `demand_sweep_seconds`/
   `unpredicted_pairs`); prior records are not comparable on those keys.
3. **The graph now depends on `candidatechunks!` being thread-safe on a shared
   index.** It already had to be — the lazy executor calls it from every worker
   — but the graph builder is a new caller of that requirement.
4. **`_node_box`'s east edge is still loose** for a single tile (option (c)).
   Harmless now that nothing compares it against a second cap construction, but
   it makes the executor's tile candidates slightly wider than they need to be.
5. This changes *which* pairs the graph holds, so any stored adjacency
   (`chunk_dag_poc`'s `ADJACENCY` dump, and the driver's schedule derived from
   it) must be rebuilt, not reused.

STATUS: PROPOSED — not pushed, no PR.
