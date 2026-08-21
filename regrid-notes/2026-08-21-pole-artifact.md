# FINAL CopDEM synthetic-store pole artifact diagnosis

**Date:** 2026-08-21  
**Production store (read-only):**
`/home/asinghvi17/geo/dggstores/copdem90-igeo7-l12-synthetic.zarr`  
**Verdict:** **not persisted corruption and not a historical numerical-code
bug.** The store is a deterministic **hybrid-source** product. The production
driver ran with `real=auto`, so four cached real Copernicus GeoTIFFs replaced the
analytic synthetic source at S85/S89/S90. The later shakedown ran with
`real=:none`, so it correctly produced the analytic field instead. The apparent
3.3 km residual is a comparison between two different source fields.

The range scan found **11 columns / 466,323 finite cells** outside the analytic
field's rigorous `[-1400, 1600] m` envelope. All 11 were written in physical run
session 3 at `e2f90d1`, with the concurrent GC sweeper off. Era-matched
`e2f90d1 + real=auto` recomputation is bit-identical to production in all four
tested columns, including 123203 and 123204. The same commit and Manifest with
only `real=none` changed produces in-range pole columns and agrees with the
later shakedown to at most `1.22e-4 m`.

## 1. Analytic bound

At production-era commit `e2f90d1`, `scripts/copdem_production.jl` defines

```text
f(lambda, phi) = 1000 sin(3 lambda) cos(2 phi)
               +  500 cos(7 lambda) sin(5 phi) + 100.
```

Each product of sine/cosine factors lies in `[-1, 1]`. Therefore

```text
-1000 - 500 + 100 <= f <= 1000 + 500 + 100,
```

so every analytic source post lies in the rigorous global envelope
`[-1400, 1600] m`. A conservative mean with non-negative weights, including a
renormalised mean over valid posts, cannot leave the convex hull of its finite
source values. This is the bound used for the scan.

The production driver contains an important exception to that argument: its
documented default is `real=auto`, and its source implementation says cached
GeoTIFFs decode for real while every other listed tile is synthesised. Each
production launch logged these four files:

- `Copernicus_DSM_COG_30_S85_00_E000_00_DEM.tif`
- `Copernicus_DSM_COG_30_S89_00_E000_00_DEM.tif`
- `Copernicus_DSM_COG_30_S90_00_E000_00_DEM.tif`
- `Copernicus_DSM_COG_30_S90_00_E001_00_DEM.tif`

Thus the analytic envelope is a valid detector of non-analytic source data, but
it is not a valid physical range bound for those real elevations.

## 2. Full-store range scan

The scanner opened the Zarr in mode `"r"`, enumerated the 54,917 physically
written `elevation` chunks, and streamed one decompressed column per thread. It
ran under `nice -n 10`, Julia 1.12.6, `-t 8`, and the safe single-field
`--gcthreads=4`. It read all 54,917 chunks in 181.94 seconds. All-NaN columns
have no chunk file and cannot contain an out-of-range finite value.

All violations are above `1600 m`; there are none below `-1400 m`.

| column | centroid latitude | outside cells | finite min (m) | finite max (m) | smallest outside (m) | real source footprint |
|---:|---:|---:|---:|---:|---:|---|
| 122975 | -84.527199 | 49,637 | -347.77924 | 2546.4668 | 1600.5034 | S85 E000 |
| 122977 | -84.232301 | 141,539 | -514.04596 | 2502.6875 | 1600.4283 | S85 E000 |
| 122979 | -85.032677 | 36,542 | -442.64798 | 2579.2683 | 1600.5371 | S85 E000 |
| 122981 | -84.730039 | 88,680 | -543.63880 | 2552.2010 | 1601.1693 | S85 E000 |
| 123147 | -83.732502 | 4,491 | -477.03638 | 2435.7900 | 1608.7020 | S85 E000 north-post fringe |
| 123172 | -87.839876 | 8,554 | -556.62960 | 2680.9856 | 1615.0426 | S89 E000 |
| 123176 | -88.345348 | 50,259 | -584.06964 | 2695.1165 | 1600.2650 | S89 E000 |
| 123199 | -89.336266 | 33,614 | -885.49260 | 2796.6907 | 1602.7502 | S89/S90 E000/E001 |
| 123202 | -88.846935 | 44,020 | -585.84564 | 2742.0122 | 1600.8232 | S89/S90 E000/E001 |
| 123203 | -89.735117 | 7,754 | -963.19745 | 2830.4058 | 1601.4515 | S90 E000/E001 |
| 123204 | -89.533430 | 1,233 | -396.74210 | 2792.6772 | 1603.5569 | S90 E000/E001 |
| **total** | **-89.735 to -83.733** | **466,323** | | **2830.4058** | | **11 columns** |

Column 123147 looks just north of the nominal S85 tile by its ancestor centroid,
but its 4,491 violations occupy longitude `0.0233..0.9967 E` and latitude
`-84.0331..-83.99968`, exactly the northern pixel-support fringe of the S85 E000
pixel-is-point tile. It is not a detached second cluster.

The refined tile/column dependency graph produces 15 possible cached-real-tile
candidates. Comparing all 15 against shakedown A (`real=none`) finds material
differences above `1e-3 m` in exactly the 11 range-failing columns above. The
other four candidates (122969, 122978, 123171, 165221) are equal or differ by at
most `1.22e-4 m`, consistent with later geometric arithmetic drift. Across the
11 material columns, 475,415 cells differ from the pure-synthetic shakedown by
more than `1e-3 m`; 466,323 of those differences are independently detectable by
the analytic range test.

## 3. True run-session boundaries and attribution

The done-log timestamps and run log reconstruct three physical writer sessions:

| physical session | launch / done-log span (UTC) | code | GC sweeper | columns written | end |
|---|---|---|---:|---:|---|
| S1 | 2026-08-20 18:32:30 / 18:33:05 to 2026-08-21 04:13:35 | `47961ea` (`claude/copdem-production`) | on (`8,1`) | 28,011 | deliberate SIGTERM for upgrade |
| S2 | 2026-08-21 04:15:02 / 04:15:37 to 07:04:39 | `7627ee3` (`claude/perf-ladder`) | on (`4,1`) | 26,315 | SIGSEGV near columns 121514-121700 |
| S3 | 2026-08-21 08:13:45 / 08:14:45 to 13:53:35 | `e2f90d1` (`claude/perf-ladder`) | **off** (`4`) | 11,852 | clean completion |

All 11 affected columns are S3 output:

| column | done timestamp (UTC) | worker | time after S2's last completed row | near segfault window? |
|---:|---|---:|---:|---|
| 122975 | 08:53:33 | 1 | 1:48:54 | no |
| 122977 | 08:55:08 | 1 | 1:50:29 | no |
| 122979 | 08:57:04 | 1 | 1:52:25 | no |
| 122981 | 08:51:00 | 23 | 1:46:21 | no |
| 123147 | 09:02:30 | 23 | 1:57:51 | no |
| 123172 | 09:11:51 | 16 | 2:07:12 | no |
| 123176 | 09:04:08 | 1 | 1:59:29 | no |
| 123199 | 09:09:54 | 22 | 2:05:15 | no |
| 123202 | 09:21:42 | 22 | 2:17:03 | no |
| 123203 | 09:25:48 | 22 | 2:21:09 | no |
| 123204 | 09:29:07 | 22 | 2:24:28 | no |

Pattern: **all S3, all Antarctic cached-real-tile footprints, time-clustered over
38 minutes, and none in a sweeper-on session or the segfault frontier.** This is
the opposite of the persisted page-corruption hypothesis's prediction.

## 4. Era-matched recomputation

A detached worktree was created at
`/home/asinghvi17/geo/DGG-pole-diag-s3`, exact commit
`e2f90d15547446793e781ca0495f0a0fdca0da0b`. Its Manifest was copied from
`/home/asinghvi17/geo/DGG-subzone-store/Manifest.toml`; both copies have SHA-256
`ad9dab81405e8eaa01d864f3bee2ce59cde5ef80dd2f9c040b7ad6a921ec66a6`.

Four columns were recomputed to
`/home/asinghvi17/geo/scratch-stores/pole-diag-s3.zarr` with one worker, four
Julia threads, `--gcthreads=2` (single field; banner `gcsweep=0`), the production
data root, and era-default `real=auto`. The run logged four real tiles decoded.

| column | production vs era `real=auto` bit differences | NaN-mask differences | max absolute difference | result |
|---:|---:|---:|---:|---|
| 122975 | 0 | 0 | 0 | bit-identical |
| 123176 | 0 | 0 | 0 | bit-identical |
| 123203 | 0 | 0 | 0 | bit-identical |
| 123204 | 0 | 0 | 0 | bit-identical |

A controlled second run used the same worktree, Manifest, worker/thread shape,
and inputs except `real=none`, writing only to
`pole-diag-s3-synthetic-only.zarr`:

| column | era `real=auto` outside range | era `real=none` outside range | production vs era `real=none` max residual | era `real=none` vs later shakedown A |
|---:|---:|---:|---:|---:|
| 123203 | 7,754 | 0 | 3301.082123 m | 23 bits, max 1.22e-4 m |
| 123204 | 1,233 | 0 | 3270.526215 m | 14 bits, max 6.10e-5 m |

This same-commit source-mode A/B is stronger than a landmark code bisect for
this question: no code change is required to turn the result from production-
exact real elevation into the analytic field. Consequently the requested
landmark bisect over `47961ea`, `e2f90d1`, `0d0f069`, and `22cf2cf` is not
applicable. Snyder closed form and the later kernels account only for the shown
sub-millimetre Float32 drift; they cannot account for the 3.3 km source-field
difference.

## 5. Verdict and verification-status amendment

### Classification

The values are **persisted, but they are not corruption**. They are exactly what
the era code deterministically computes when `real=auto` sees the four cached
GeoTIFFs. Nor are they a historical conservative-weight or pole-geometry bug.
Relative to an intended all-synthetic-store contract, this is a **historical
configuration/provenance bug**: `real=auto` admitted real data into a store whose
title and root `source` attribute describe it as synthetic.

What made the shakedown “good” was its explicit `real=:none`, not the GeometryOps
repin, perf-ladder kernels, Snyder closed form, or DAG driver.

### Earlier `COMPLETE + integrity PASS`

The physical claims remain valid:

- 66,178/66,178 columns are complete;
- all chunks decompress and their lengths/NaN counts agree with the ledger;
- sampled columns reproduce the era code bit-for-bit.

The semantic claim needs amendment. The store is not globally a pure analytic-
synthetic oracle. Record it as:

> **COMPLETE; physical integrity PASS; source provenance HYBRID; pure-synthetic
> semantic purity FAIL in 11 columns / 466,323 out-of-envelope cells.**

The earlier verification did not merely miss column 123203 by random sampling.
It sampled three affected columns—122975, 123176, and 123199—explicitly labelled
them `real tile`, and reproduced them using the same `real=auto` configuration.
Its analytic oracle deliberately skipped real-backed tiles and separately
checked that all four real-tile centres were finite. That correctly verified the
hybrid product it was told to verify, but it could not detect a mismatch between
the store's “synthetic” label and the desired pure-synthetic interpretation.

### Comparison-oracle blast radius

For comparisons to `real=none` pure-synthetic runs, exclude these columns:

```text
122975, 122977, 122979, 122981, 123147, 123172,
123176, 123199, 123202, 123203, 123204
```

All other refined cached-real-tile candidates were tested and have no material
`>1e-3 m` source-mode difference. The store remains a valid era-matched hybrid
oracle outside the 11-column exclusion set, subject to the already documented
`<=1.22e-4 m` later arithmetic drift when comparing different code landmarks.

## 6. Read-only assurance and artifacts

The production Zarr and ledger were only opened for reading. The ledger MD5 is
still `84e11e22cd263a8692fe034f1948abca`, matching the earlier verification, and
the production metadata mtimes remain the original 2026-08-20 creation times.
No `.zarray` or `.zattrs` was changed.

Scratch artifacts:

- `/home/asinghvi17/geo/scratch-stores/pole-diag-tools/range-scan.ndjson`
- `/home/asinghvi17/geo/scratch-stores/pole-diag-tools/era-comparison.ndjson`
- `/home/asinghvi17/geo/scratch-stores/pole-diag-tools/source-mode-comparison.ndjson`
- `/home/asinghvi17/geo/scratch-stores/pole-diag-tools/real-influence-comparison.ndjson`
- `/home/asinghvi17/geo/scratch-stores/pole-diag-s3.zarr`
- `/home/asinghvi17/geo/scratch-stores/pole-diag-s3-synthetic-only.zarr`

## 7. Open questions

1. Was the four-tile hybrid exception intended for this FINAL store, or was the
   intended deliverable pure synthetic? The era script and old verification
   explicitly knew about the exception, but the store title/root source metadata
   do not disclose it clearly.
2. Because the store is immutable, the metadata cannot be corrected in place.
   The external catalog/verification verdict should carry the hybrid-source
   amendment, or a new pure store should be generated with `real=none`.
3. If a bit-for-bit pure-synthetic oracle is required rather than a `1e-3 m`
   numerical oracle, regenerate the 11 columns (or the whole store) at one chosen
   code landmark; cross-landmark geometric changes already cause isolated
   `1e-4 m` Float32 drift independent of this source-mode issue.
