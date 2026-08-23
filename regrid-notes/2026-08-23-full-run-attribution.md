# Full-run authalic synthetic time and allocation attribution

## Decision summary

The accepted run completed its data product: **[measured] 66,228 / 66,228** destination columns are represented by the ledger/store union, and the store contains **[measured] 54,956** materialized chunk files plus **[measured] 11,282** all-NaN columns that intentionally need no file. The final nonzero exit was solely the predicted-edge invariant. It was not an incomplete computation.

The whole-run result was **[measured] 8.81 h** total wall, **[measured] 1,722,052 cells/s** over **[measured] 5.4541e10 cells**, **[estimated] 148.10 core-h**, **[estimated] 16.81 mean active cores**, **[measured] 25.19 GiB peak RSS**, and **[estimated] 17.8% GC-pause wall fraction**. The GC number is explicitly a model-based estimate: this production tip emitted no numeric GC-wall counter and none of the flat profile dumps contains a GC frame.

For the Phase-5 gate, W1/W2 has a real but bounded target. The structurally countable W1 redundancy is **[estimated] 6.079 GiB per global-average 100 columns**, or **[estimated] 3.931 TiB** over this run, from one copied denominator vector and one equally sized reference vector per source/destination block. W2 can additionally remove a COO copy and a second CSC construction, but the run recorded no `nnz` or allocation stacks, so its byte recovery cannot be quantified honestly. Those changes do not remove candidate discovery, clipping, per-block operators/memos, or the **[estimated] 3.068 GiB per 100 columns** destination cap-vector payload. The available evidence therefore supports implementing and measuring W1/W2; it does not support promising a several-fold allocation reduction or a sub-15% GC fraction.

The **[measured] 13** uncredited demands are genuine graph/executor coverage mismatches, not counter bookkeeping. A read-only reconstruction over every destination found **[measured] 71** executor candidate pairs absent from the graph and **[measured] 437** graph-only pairs. Every one of the **[measured] 13** runtime source indices is in the missing-pair set. Disabling the graph latitude prefilter produced the same **[measured] 328,424** edges and the same **[measured] 71** misses, so the fault is not the latitude-band optimization; it is a mismatch between the graph's direct chunk-cap relation and the DGG hierarchical candidate relation used by the lazy executor. No code was changed.

Throughout this report, every quantitative result is tagged **measured** when it comes directly from an artifact or exact reconstruction, and **estimated** when it is interpolated, scaled, modeled, or attributed from samples. Profile frame counts are measured inclusive flat counts; they overlap unless a disjoint subtraction is stated.

## Run completion and evidence inventory

The analyzed process is the accepted launch at commit `1654a703f5647e82affefaad98bfa43c6ac3e96d`, PID **[measured] 2,782,183**, launched **[measured] 2026-08-22 23:42 UTC** with outer scheduling, Julia `-t21`, **[measured] 40** worker tasks, **[measured] 4** GC threads, no concurrent sweeper, **[measured] 8**-column batches with taper, affinity ordering, refcount source caching, and no prefetch. The input was synthetic and performed **[measured] zero** network reads. Compute geometry was authalic; storage identity was bare IGeo7 with destination-geometry provenance.

Completion checks were independent:

- The final log banner landed at **[measured] 2026-08-23 08:30:32 UTC** and reports **[measured] 66,228 / 66,228** chunks computed this session, **[measured] 0** skipped, **[measured] 26.895%** NaN cells, and **[measured] 66,228 / 66,228** in the store.
- The done log has **[measured] 66,228** completion rows, and the columns cache lists **[measured] 66,228** columns.
- The store/ledger union is **[measured] 66,228 / 66,228**. The store occupies **[measured] 74 GB** on disk; the absent physical files correspond exactly to **[measured] 11,282** all-NaN columns.
- Final policy counters report **[measured] 26,475** retained source loads covering all **[measured] 26,475** source chunks, **[measured] 301,564** hits, **[measured] 6** joined loads, **[measured] 0** live tiles at exit, **[measured] 1,526** peak resident tiles / **[measured] 2.53 GiB**, and **[measured] 13** transient uncredited demands. All **[measured] 66,228** destinations retired.
- All **[measured] 4 / 4** requested profile dumps landed and are nonempty: manual early **[measured] 103,725** aggregate samples at about **[measured] 87%** sampled utilization, early-midlatitude **[measured] 96,775** at **[measured] 96%**, polar-heavy **[measured] 97,750** at **[measured] 96%**, and tail **[measured] 96,925** at **[measured] 96%**.

The accepted artifacts are the `glo90-synthetic-authalic-phase1` log, process watcher, profiles, ledger/columns cache, and Zarr store. Neither aborted-launch prefix is included anywhere in the calculations.

## Headline performance and epoch behavior

The watcher spans **[measured] 31,696.208 s** from its first to last usable sample and accumulates **[measured] 533,084.73 CPU-s**, or **[measured] 148.079 core-h**, over that covered interval. Integrating adjacent cumulative CPU-tick deltas is equivalent to integrating its instantaneous-core field. Adding the **[estimated] 22.58 CPU-s** before the first watcher interval and the final approximately **[estimated] 6 s** of process life gives **[estimated] 148.10 core-h** for the process. Dividing by total wall gives **[estimated] 16.81** integrated mean cores; this is not a whole-call average of a sampled gauge.

| Headline | Result |
|---|---:|
| Total wall | **[measured] 8.81 h** |
| Compute/session wall in final banner | **[measured] 8.80 h** |
| Aggregate throughput | **[measured] 1,722,052 cells/s** |
| Last heartbeat throughput | **[measured] 1,722,070 cells/s** |
| Integrated process CPU | **[estimated] 148.10 core-h** |
| Integrated mean cores | **[estimated] 16.81** |
| Peak RSS / exit RSS | **[measured] 25.19 GiB / 15.07 GiB** |
| GC-pause wall fraction | **[estimated] 17.8%**, plausible model range **[estimated] 16.8-19.9%** |
| Completed columns | **[measured] 66,228 / 66,228** |

There were no explicit phase markers, so “epoch” below means a temporal quartile of completed-column progress. Boundaries and rates are **estimated** by linear interpolation between measured 60-second watcher samples. This avoids assigning geographic labels from a handful of recent-latitude observations.

| Epoch | Completed columns | Wall | Mean cores | Throughput |
|---|---:|---:|---:|---:|
| First progress quartile | **[measured] 16,557** | **[estimated] 7,922.715 s** | **[estimated] 16.506** | **[estimated] 1,721,052 cells/s** |
| Second progress quartile | **[measured] 16,557** | **[estimated] 8,003.351 s** | **[estimated] 16.851** | **[estimated] 1,703,712 cells/s** |
| Third progress quartile | **[measured] 16,557** | **[estimated] 8,210.583 s** | **[estimated] 16.874** | **[estimated] 1,660,711 cells/s** |
| Final progress quartile | **[measured] 16,557** | **[estimated] 7,565.560 s** | **[estimated] 17.038** | **[estimated] 1,802,299 cells/s** |

The fastest final quartile and slowest third quartile show that progress was not monotonically degraded by polar work. The integrated result, rather than a single profile window, is the correct capacity number.

For a **[measured] 21-thread** recipient machine capable of this measured throughput, the compute-only ETA is `5.4541e10 / 1,722,052`, or **[estimated] 8.798 h (8 h 48 min)**. Adding the accepted run's measured setup/teardown scale gives **[estimated] 8.81 h (about 8 h 49 min)** end to end for the same synthetic, zero-network workload. This is not an ETA for real Copernicus downloads or a machine with different single-core/GC performance.

## Time attribution

### Method

The three scheduled profiles contain **[measured] 291,450** aggregate samples and are stable enough to form the whole-run mix. The manual health-gate capture is retained as a cross-check but excluded from proportions because its **[measured] 87%** sampled utilization and startup/early-wave mix differ from the scheduled captures' **[measured] 96%**.

The disjoint scheduled-profile partition is:

- Clipping: **[measured] 70,093** `IntersectionAreaOperator` samples.
- Non-clipping operator/memo work: **[measured] 26,703**, calculated as **[measured] 96,796** `BlockAreaOperator` minus the clipping frames.
- Sparse assembly outside block geometry: **[measured] 2,156**, calculated as **[measured] 98,952** `_assemble_sparse` minus `BlockAreaOperator`.
- Serial block finalization: **[measured] 3,974**, the disjoint sum of **[measured] 1,593** `_fillcoo!` and **[measured] 2,381** `WeightBlock` samples.
- Scheduler/idle: **[measured] 8,516** `poptask` samples.
- Store: two scheduled dumps censor `dggwrite!` below their **[measured] 20**-count printing threshold. Scaling the observed all-profile rate gives **[estimated] 61** scheduled-equivalent samples. All four dumps contain **[measured] 83** observed `dggwrite!` samples.
- Compute kernel: the exact residual after the buckets above, **[estimated] 206,650** samples. Named measured components inside it include **[measured] 26,554** cap-tree `_tiletree!` samples, **[measured] 66,709** `get_all_candidate_pairs!` samples, and the **[measured] 26,703** non-clipping operator/memo samples; the remaining **[estimated] 86,684** samples are other tree traversal, geometry, memo, and kernel work.

The production heartbeat has no numeric GC-wall counter. Across all four flat dumps there are **[measured] 0** `jl_gc`, sweep, or collector frames. `collect` occurrences are ordinary array/data collection, not GC. Julia's stop-the-world intervals are consequently missing from the conditional “workers running” profile mix. GC is estimated from the integrated-core deficit with one caller parked: `C = 20(1-g) + k*g`, where **[measured] 20** threads can do worker work and **[estimated] k = 2.1** effective cores remain busy during a pause. The `k` estimate is calibrated from the matched pre-fix and t21-remeasure W40 runs: **[measured] 14.737** mean cores with **[measured] 29.54%** GC gives **[estimated] k = 2.18**, while **[measured] 14.115** cores with **[measured] 32.91%** GC gives **[estimated] k = 2.12**. With **[estimated] C = 16.81**, this yields **[estimated] g = 17.8%**; allowing **[estimated] k = 1-4** gives **[estimated] 16.8-19.9%**. This is a model, not a hidden production counter.

### Whole-run buckets

The core-hour column first assigns **[estimated] 3.293 core-h** to GC using that model, then distributes the remaining **[estimated] 144.807 core-h** in the scheduled-profile proportions. “Throughput-equivalent wall” is core-hours divided by **[estimated] 16.81** integrated mean cores; these values sum to total wall but are not literal mutually exclusive timers because **[measured] 40** worker tasks overlap. The GC row separately reports its estimated stop-the-world wall.

| Phase | Profile basis | Share of conditional scheduled samples | Estimated core-hours | Estimated throughput-equivalent wall | Literal/observable wall evidence |
|---|---:|---:|---:|---:|---:|
| Compute kernel, including non-clipping operator/memo | **[estimated] 206,650 samples** | **[estimated] 70.904%** | **[estimated] 102.674 core-h** | **[estimated] 6.108 h (366.5 min)** | No exclusive timer |
| Clipping | **[measured] 70,093 samples** | **[measured] 24.050%** | **[estimated] 34.826 core-h** | **[estimated] 2.072 h (124.3 min)** | No exclusive timer |
| Sparse assembly outside block geometry | **[measured] 2,156 samples** | **[measured] 0.740%** | **[estimated] 1.071 core-h** | **[estimated] 0.064 h (3.8 min)** | No exclusive timer |
| Store write | **[estimated] 61 samples** | **[estimated] 0.021%** | **[estimated] 0.030 core-h** | **[estimated] 0.0018 h (about 6 s)** | Blocking I/O below sampling resolution cannot be excluded |
| GC | **[measured] 0 visible GC frames** | Not present in conditional profiles | **[estimated] 3.293 core-h** | **[estimated] 0.196 h (11.8 min)** | **[estimated] 1.568 h (94.1 min)** stop-the-world wall |
| Scheduler + idle | **[measured] 8,516 samples** | **[measured] 2.922%** | **[estimated] 4.231 core-h** | **[estimated] 0.252 h (15.1 min)** | No exclusive timer |
| Serial block finalization tail | **[measured] 3,974 samples** | **[measured] 1.364%** | **[estimated] 1.974 core-h** | **[estimated] 0.117 h (7.0 min)** | Per-block serial work, overlapped across outer workers |
| Total | **[measured] 291,450 conditional samples** | **[measured] 100.000%** | **[estimated] 148.10 core-h** | **[measured] 8.81 h** | **[measured] 8.81 h** process wall |

### Serial tail, called out explicitly

There are two different tails and they should not be conflated:

1. The repeated block-local serial finalization (`_fillcoo!` followed by `WeightBlock`) is **[measured] 3,974 / 291,450 samples (1.364%)**, equivalent to **[estimated] 7.0 min** of whole-run throughput wall. It is serial within one block but overlaps among outer workers.
2. The global queue drain from the first worker's final completion at **[measured] 08:30:07 UTC** to the last worker's final completion at **[measured] 08:30:28 UTC** is only **[measured] 21 s**. Shutdown from that point to the final banner is **[measured] 4 s**. This is the literal end-of-run serial/underfilled tail.

### Comparison with the 2026-08-19 post-stack profile

That prior report is a spherical, eager, warm, **[measured] 8-thread** MOC case with explicit phase timers and allocation profiling; the present run is authalic, chunked, outer-W40, and has CPU flat profiles. Comparisons are directional, not byte- or phase-identical.

| Comparable signal | Prior spherical eager result | Current authalic full-run profile |
|---|---:|---:|
| Cap-tree construction | **[measured] 8.1% CPU** | **[measured] 9.11%** of scheduled samples (`_tiletree!`) |
| Candidate-pair discovery | **[measured] 41.5% CPU** | **[measured] 22.89%** of scheduled samples |
| Clipping/COO family | **[measured] 47.0% CPU** | **[measured] 33.95%** for clipping + non-clipping operator/memo + sparse overhead |
| Sparse/serial block finalization | **[measured] about 3.3% CPU** and **[measured] 7.77 s / 50.82 s** wall tail | **[measured] 1.364%** block-finalization samples; **[estimated] 7.0 min** throughput-equivalent wall; **[measured] 21 s** global drain |
| GC | **[measured] 5.43 s / 50.82 s (10.7% wall)** | **[estimated] 17.8%** pause wall, with no production counter |
| Peak RSS | **[measured] 38.58 GiB** | **[measured] 25.19 GiB** |

The apparent candidate-share improvement is partly a denominator and execution-shape change; it is not evidence of a directly transferable speedup. Likewise, the prior `Profile.Allocs` ownership is materially stronger allocation evidence than the current CPU-only dumps.

## Allocation attribution

### What is measured, and what is not

The matched **[measured] 100-column** assembly-cache A/B is the only byte-total anchor: cache-less **[measured] 285.1588 GiB**, pooled **[measured] 242.2104 GiB**, saving **[measured] 42.9484 GiB (15.06%)**. GC fell **[measured] 31.40% -> 26.25%**, wall fell **[measured] 100.40 s -> 93.34 s (7.03%)**, and mean cores rose **[measured] 11.78 -> 12.35**. Commit `de6d6aca` is included in the overnight run, so **[measured] 242.2104 GiB per that 100-column workload** is the relevant unexplained remainder.

Scaling that Antarctic/sparse-subset rate uniformly to **[measured] 66,228** columns would imply **[estimated] 156.65 TiB** of cumulative allocations and **[estimated] 27.78 TiB** saved by the pool. These are deliberately labeled estimates, not full-run measurements: block degree, clipping `nnz`, latitude, and NaN behavior vary geographically. The production shutdown emitted no total allocated-byte counter, and its profile files are CPU flat profiles rather than `Profile.Allocs`. An exact byte split of the **[measured] 242.2104 GiB** remainder is therefore impossible from these artifacts.

The primary profile evidence, aggregated across all **[measured] 4** dumps, is:

| Inclusive frame | Measured flat count | Allocation implication |
|---|---:|---|
| `BlockAreaOperator` | **[measured] 118,366** | Operator/memo/clipping family is pervasive |
| `IntersectionAreaOperator` | **[measured] 85,865** | Nested clipping work |
| Non-clipping operator/memo difference | **[measured] 32,501** | Operator, maps, cells, and memo setup outside the clip frame |
| `get_all_candidate_pairs!` | **[measured] 86,882** | Candidate buffers/tree traversal |
| `_tiletree!` / cap-tree root | **[measured] 33,999** | Destination geometry/cap-tree creation |
| Base `Array` / `GenericMemory` | **[measured] 29,866 / 30,974** | Broad array allocation/copy pressure, not owner-separable |
| `_growend!` / `resize!` | **[measured] 7,248 / 1,133** | Dynamic-vector growth in candidate/COO paths |
| SparseArrays `sparse(...)` | **[measured] 4,535** | CSC construction; overlaps callers |
| `_fillcoo!` | **[measured] 2,082** | First CSC copied into COO vectors |
| `WeightBlock` | **[measured] 3,370** | Second CSC plus denominator copy |
| `copy` | **[measured] 852** | Copies exist but the flat name does not identify all owners |
| cursor-memo `node_extent` / `_taskmemo` | **[measured] 4,467 / 380** | Cell/cap memo and cursor lifecycle |
| `dggwrite!` / DiskArrays `collect` | **[measured] 83 / 2,135** | Store frame is tiny; output collection is visible |

These counts are CPU samples, not allocation counts. They identify active ownership paths but cannot apportion bytes proportional to the counts.

### Owner 1: denominator and reference-vector lifecycle

Every source/destination block constructs a `WeightCOO` denominator, `WeightBlock` copies it, and `getblock!` then constructs a same-length `CachedBlock.ref`. A destination column has **[measured] 823,543** cells, so each `Vector{Float64}` payload is **[derived/estimated] 6,588,344 bytes**. Production counters and the reconstructed candidate relation agree on **[measured] 328,058** demands/blocks: **[measured] 26,475** cache-load leaders + **[measured] 301,564** hits + **[measured] 6** joined demands + **[measured] 13** uncredited demands.

One such vector over all blocks is **[estimated] 2,012.922 GiB (1.966 TiB)** of cumulative payload. Retaining one denominator is necessary; the copied denominator plus reference vector are the W1-recoverable redundancy, **[estimated] 4,025.845 GiB (3.931 TiB)** over the run. At the run-wide **[estimated] 4.9535 blocks/column** mean, that is **[estimated] 6.079 GiB per 100 columns**, or **[estimated] 2.51%** of the matched **[measured] 242.2104 GiB** anchor. This is a conservative structural payload calculation, not a measured claim that the geographically specific A/B had exactly the global-average block degree.

### Owner 2: `_fillcoo!` CSC -> COO -> CSC round trip

The conservative block path first returns a CSC from `_intersectionareas`; `_fillcoo!` then appends its entries to `rows::Vector{Int}`, `cols::Vector{Int}`, and `vals::Vector{Float64}`; `WeightBlock` calls `sparse(...)` to build a second CSC. The redundant payload is at least **[derived/estimated] 24 bytes per `nnz`** for the three COO arrays plus **[derived/estimated] 16 bytes per `nnz`** for the second CSC's row/value arrays, before vector capacity slack, CSC column pointers, sort/workspace, and copies. The **[measured] 7,248** `_growend!`, **[measured] 2,082** `_fillcoo!`, **[measured] 4,535** `sparse(...)`, and **[measured] 3,370** `WeightBlock` frames make this a credible churn owner.

There is no recorded full-run or A/B `nnz` count, so converting that formula to GiB would be fabrication. Its **[measured] 1.364%** serial CPU-sample share does not imply a similarly small allocation share: bulk array construction can allocate many bytes in few samples.

### Owner 3: per-block operators, memos, and candidate structures

`BlockAreaOperator`, index maps, source/destination `CellMemo`s, boundary/cap memoization, intersection trees, and candidate vectors are created for each block. The all-dump subtraction leaves **[measured] 32,501** non-clipping `BlockAreaOperator` samples, supplemented by **[measured] 4,467** cursor-memo `node_extent`, **[measured] 380** `_taskmemo`, **[measured] 86,882** candidate-pair, **[measured] 7,248** growth, and **[measured] 1,133** resize frames.

The assembly pool only reuses `SparseMatrixAssemblyCache` buffers. It does not own these operators, `WeightCOO`, `CellMemo`, candidate arrays, or tree/cursor structures. Their combined bytes cannot be separated with CPU stacks, and no split is assigned here.

### Owner 4: destination geometry and cap-tree construction

The cached destination tree materializes one `Vector{SphericalCap{Float64}}` for each destination column, not once for each source block. Each cap is **[measured from the runtime type] 40 bytes**, so the raw vector payload is **[derived/estimated] 32,941,720 bytes per column**, **[estimated] 3.068 GiB per 100 columns**, and **[estimated] 1.984 TiB** over **[measured] 66,228** columns. This is **[estimated] 1.27%** of the **[measured] 242.2104 GiB** anchor on a 100-column basis, before vector/object overhead. The **[measured] 33,999** cap-tree samples support this code path. W1/W2 does not directly remove it.

### Owner 5: store collection, serialization, and payload

One dense Float32 output vector is only **[derived/estimated] 3,294,172 bytes per column**, or **[estimated] 0.3068 GiB per 100 columns**. The entire logical dense output is **[estimated] 203.18 GiB**, while compression/omitted NaN chunks produce **[measured] 74 GB** on disk. The profiles show **[measured] 83** `dggwrite!` samples and **[measured] 2,135** DiskArrays `collect` samples across all dumps. That evidence does not support store serialization as the dominant owner of a **[measured] 242.2104 GiB** allocation total. It also cannot exclude short-lived compression buffers or unsampled blocking I/O, so the store is not assigned a zero-byte share.

### What remains inseparable

The countable two-vector W1 redundancy plus destination cap payload is only **[estimated] 9.147 GiB per global-average 100 columns (3.78%)** against the matched **[measured] 242.2104 GiB** remainder. A necessary denominator, CSC/COO/CSC entries and workspace, per-block operators/memos, candidate vectors, cell boundaries, intersection geometry, and Julia/runtime overhead make up the unresolved majority. The CPU flat profiles prove those paths are active but provide no sound byte ratios. No further percentage split is claimed.

### W1/W2 gate: plausible recovery and non-recovery

Plausible recovery:

- W1 can eliminate the redundant denominator copy and `CachedBlock.ref` lifecycle: **[estimated] 6.079 GiB per global-average 100 columns**, **[estimated] 3.931 TiB** full-run-equivalent payload.
- W2 can adopt the first CSC and eliminate the intermediate COO arrays and second CSC. Recovery scales with unrecorded `nnz`; the formula above is evidence-based, but a GiB figure is not available.
- Both can shorten allocation-driven GC even though their direct CPU frames are small. The amount must be measured at driver scale.

Not recovered:

- **[measured] 70,093 / 291,450** scheduled clipping samples, candidate discovery, and the non-clipping operator/memo family.
- The **[estimated] 3.068 GiB per 100 columns** destination cap vector.
- Source tile payload, store output/compression, cell-boundary geometry, and other memo/cursor arrays.
- GC caused by all of those remaining owners. Nothing here substantiates a target below **[estimated] 15%** GC or a several-fold total-allocation reduction from W1/W2 alone.

## Predicted-edge failure: graph miss, not accounting artifact

I reconstructed both relations read-only, using the accepted ledger's **[measured] 66,228** destination chunk roots and the production source space of **[measured] 26,475** tiles:

1. `chunk_dependency_graph(dstspace, srcspace; radius)` produced **[measured] 328,424** edges.
2. For each destination, the exact `candidatechunks!` query used by `LazyRegridArray._connectedsource!` produced **[measured] 328,058** demand pairs.
3. Set comparison found **[measured] 71** demand pairs absent from the graph and **[measured] 437** conservative graph-only pairs.
4. Rebuilding with `prefilter=false` still produced **[measured] 328,424** edges and left all **[measured] 71** misses. The latitude-band prefilter is exonerated.

The runtime's **[measured] 13** uncredited sources all occur in those **[measured] 71** missing relations:

| Runtime source index | One missing destination position / chunk root | Source tile |
|---:|---:|---|
| **[measured] 147** | **[measured] 7,070 / 11,019** | `N81 W062` |
| **[measured] 8,571** | **[measured] 16,887 / 30,576** | `N44 W117` |
| **[measured] 17,945** | **[measured] 30,354 / 60,031** | `S23 E016` |
| **[measured] 13,478** | **[measured] 38,657 / 73,334** | `N16 E075` |
| **[measured] 11,497** | **[measured] 41,937 / 77,487** | `N29 E121` |
| **[measured] 3,555** | **[measured] 45,580 / 82,625** | `N65 E100` |
| **[measured] 6,682** | **[measured] 45,627 / 82,672** | `N53 E089` |
| **[measured] 13,103** | **[measured] 47,578 / 90,638** | `N18 W113` |
| **[measured] 23,465** | **[measured] 54,272 / 121,539** | `S82 E049` |
| **[measured] 25,259** | **[measured] 55,643 / 123,179** | `S87 E043` |
| **[measured] 21,780** | **[measured] 63,600 / 157,452** | `S77 E136` |
| **[measured] 20,647** | **[measured] 63,695 / 157,547** | `S73 E143` |
| **[measured] 25,079** | **[measured] 65,926 / 165,197** | `S87 W137` |

The counter reports only a missing demand made after a source's predicted consumer refcount has already reached **[measured] 0**. The other **[measured] 58** missing pairs happened while that source still had at least one predicted consumer, so the cache served them without labeling them uncredited. Thus **[measured] 13** is scheduling-dependent evidence of a larger **[measured] 71-pair** relation mismatch, not evidence that the counter randomly miscounted.

This is a genuine predictive-completeness defect in the graph relative to executor reads, although candidate pairs are conservative and do not prove a positive cell-level intersection. It is independent of synthetic values and can recur on the real-data path with the same grids, potentially causing transient extra tile loads after graph retirement. Correctness was preserved because the cache serves uncredited demands; the accepted store is complete. No graph or executor fix was attempted in this analysis-only task.

## Mandatory caveats

1. **Sparse-subset uncredited-demand artifact.** Both **[measured] 100-column** mini-runs completed their stores but reported **[measured] 1** uncredited demand for source chunk **[measured] 25,955** and exited **[measured] 1**. That known sparse-subset harness accounting artifact was not investigated. The full run has **[measured] 13** different source indices, does not include **[measured] 25,955**, and the full relation reconstruction proves its incidents are real graph/executor misses. They must not be waved off as the mini-run artifact.
2. **Spherical versus authalic.** This run is authalic. Every prior store, byte-identity baseline, and the t21-remeasure campaign is spherical. This is a fresh measurement baseline and is not byte-comparable. Authalic point kernels add **[measured in microbenchmarks] 17-46%**; only a **[estimated] low-single-digit wall percentage** was expected after whole-driver dilution.
3. **Assembly pool already included.** The assembly-cache pool's **[measured] 15.06%** allocation reduction and **[measured] 7.03%** wall reduction are already included through commit `de6d6aca`. The remaining churn is after that improvement, not before it. The task-local wiring experiment was null at **[measured] +0.009% allocation** because the driver creates ephemeral block tasks.
4. **Phase 1 supplied no driver-scale speedup.** The t21 remeasure found W40 **[measured] 4.81% slower** than pre-fix, with GC rising **[measured] 29.54% -> 32.91%** and allocation rising **[measured] 2.53%**. Benchmark-scale ConservativeRegridding results of **[measured in that benchmark] -95% allocation** and **[measured in that benchmark] 5.2% GC** did not materialize at driver scale and must not be used as expectations.
5. **Aborted launches excluded.** The SIGUSR1-deferral launch under saturated outer waves was aborted and fixed in `f1945e7`; the W=38 heuristic launch was aborted and fixed in `1654a70`. Their artifacts at `glo90-synthetic-authalic-phase1-aborted-profile-defer-20260822T2321Z` and `glo90-synthetic-authalic-phase1-aborted-W38-20260822T2340Z` are preserved but contribute **[measured] 0** samples, watcher intervals, ledger rows, or counters to this report.

## Gate conclusion

Proceeding to a narrowly scoped W1/W2 driver-scale experiment is evidence-supported. The gate should require matched allocation totals, GC wall from a harness that actually emits it, `nnz`/block counts, correctness parity, and the graph/executor predictive check. The present report establishes a conservative, countable W1 floor and strong qualitative W2 evidence, while explicitly leaving the unmeasurable majority unsplit. It also elevates the **[measured] 71-pair** graph/executor coverage mismatch as a separate real-data readiness finding; it should not be “fixed” by interpreting the **[measured] 13** runtime counter as mere accounting.

STATUS: DONE
