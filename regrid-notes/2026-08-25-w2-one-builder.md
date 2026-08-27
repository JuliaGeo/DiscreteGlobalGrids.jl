# W2 — one build path, and Conservative adopts what it assembles

- Date: 2026-08-25
- Card: `regrid-notes/2026-08-21-regridding-simplification-plan.md`,
  "Phase 5 — one final weight block representation", Task W2
- Commit: `Unify final weight block construction`

## 1. The seam

Every weight block, eager or chunked, is assembled by `weightblock`
(`lib/GlobalRegridding/src/plans.jl`), which dispatches on `outputsampling` and
whose area branch calls

```julia
pairblock(method, dst_space, dst_inds, src_space, src_inds) -> WeightBlock
```

The generic `pairblock` fills one `WeightCOO` through `buildweights!` and
assembles a matrix from it. That remains the only hook a method has to supply,
and every point method and every third-party area method reaches it unchanged.

`weightblock` dispatches on `outputsampling` so a sampling *may* specialise the
assembly, and one does not: the `::Points` method it shipped with was
byte-identical to the `::Sampling` method it shadowed, so it decided nothing and
is gone. Every sampling now assembles a pair through `pairblock`, and the
whole-tile route a point method with a `sampler` takes is selected once per
plan, by `tilesampler` in `blockfor`, not per block here. The docstrings on
`weightblock` and `outputsampling` say that rather than describing two paths
through `weightblock` that were one.

`Conservative` specializes `pairblock` (`conservative.jl`). It measures the
intersection areas exactly as before, then **adopts the matrix
ConservativeRegridding assembles** as the block's weights and reads each
destination's denominator off that matrix once, with `_blockdenom`. Handing that
vector to `WeightBlock(weights, denom)` makes it the block's reference as well,
so one pass over the stored entries produces the denominator, and the block
holds one copy of it. No conservative block is copied into a coordinate list and
out of it again, on any route.

A degenerate side — an empty `dst_inds` or `src_inds` — `invoke`s the generic
route, which keeps the empty-side contract exactly: no destination returns
before denominators are declared and reports none, no source returns after and
reports a denominator of zeros, and neither compares the two manifolds.

`wholeblock(::Conservative, …)` is deleted. The eager whole domain is the
generic `wholeblock` over `1:ncells(dst_space)` against `1:ncells(src_space)`,
which is `weightblock` and therefore `pairblock(::Conservative, …)`; a chunked
plan's pair is `buildblock`, which is the same two calls. One builder path, one
final block representation, and the specialization that used to serve the eager
domain alone now serves both.

`buildweights!(::Conservative, …)` survives as the documented generic assembly.
It is what the degenerate sides `invoke`, and what a method wrapping
`Conservative` reaches when it forwards that hook alone.

**The forwarding rule**, stated on `pairblock` and in `buildweights!`'s
docstring: a method that wraps another and forwards its build takes the inner
method's build by forwarding `pairblock` — and `sampler` too, for a point
method. Forwarding `buildweights!` alone reaches the generic route whatever the
inner method assembles for itself. Both are legitimate; the rule says which one
a wrapper gets, so an instrumented method measures the path it means to.

Nothing about point methods, `TileWeights`, the tile route, missing policies or
values changes.

## 2. What it verifies

`lib/GlobalRegridding/test/test_conservative.jl`, 73 assertions net:

- the seam's block and the generic route's are identical *stored structure*:
  equal `colptr` and `rowval`, `===` `nzval` and `===` denominators, over the
  whole domain and over each of four chunk pairs of a non-contiguous source
  chunking; the eager `wholeblock` is that same block, entry for entry;
- a conservative pair build allocates no coordinate list: on a 288×4608 pair
  with 5,568 nonzeros, serially, the adopted route allocates 3,093,744 bytes
  against the generic route's 3,734,576 — 640,832 bytes less, 17.2 %, where the
  coordinate list's own three vectors are 133,632 bytes at minimum;
- partition invariance through the seam: two non-contiguous source chunks
  reassemble the whole block's columns exactly, and equal the coordinate-list
  reference;
- the empty sides keep their asymmetry — no destination reports no denominator
  and an empty reference, no source a denominator of zeros the block references
  — and both agree with the generic route on size, weights and denominator type;
- a disjoint pair reports zero denominators rather than the row sums of no
  weights;
- a wrapper forwarding `pairblock` gets the inner method's assembly, entry for
  entry, and its `buildweights!` is never called; a wrapper forwarding
  `buildweights!` alone is called exactly once and produces the same values; a
  chunked plan makes the same choice for both;
- one block carries one denominator, computed where it was built: a reused eager
  plan builds once, applies the stored denominator, and answers a second
  application identically and equal to plain `Conservative`.

Unchanged and still passing: the third-party emitters `T4SplitMethod`,
`T4TileMethod`, `T4PlaceCount`, `T5PlaceCount`, `UnimplementedMethod`,
`CountingMethod`, `T6CountingMethod(Conservative())` and `T7CountingMethod`,
all of which supply or forward `buildweights!` and nothing else.

The mutants these kill: a specialization whose adopted matrix differs from the
generic assembly in layout, value or denominator; a chunk-local column map taken
from the space's indices rather than the chunk's, which partition invariance
catches; an adoption that keeps the coordinate list anyway, which the allocation
floor catches; an empty side that reports the wrong denominator, which would
change what `Weighted` blanks; a seam that silently drops a wrapper back to the
generic route, or one that refuses a wrapper that only forwards `buildweights!`.

## 3. Measurements

Same session, interleaved A/B, three repetitions of each state, Julia 1.12.6,
8 threads, `powermode 2`, inner-threaded shape (one top-level build, so
`OUTER_PARALLEL` is false). Peak memory is `Sys.maxrss()` sampled after the
timed build in a fresh process per mode. Both benchmarks produce identical
nonzero counts, weight sums and block sizes in every run.

`benchmark/conservative_block_baseline.jl` — IGeo7 level 5, 168,072 destination
cells from a 360×180 `RasterGrid`, 487,174 nonzeros:

| mode | measure | before | after |
|---|---|---|---|
| direct (eager) | seconds, median (min–max) | 5.258 (5.222–5.716) | 5.211 (4.957–5.217) |
| direct (eager) | allocated, mean | 856,487,328 | 854,319,088 |
| direct (eager) | max RSS MiB, mean | 846.0 | 870.9 |
| chunked | seconds, median (min–max) | 2.752 (2.709–2.835) | 2.685 (2.616–2.815) |
| chunked | allocated, mean | 722,540,875 | 674,627,061 |
| chunked | max RSS MiB, mean | 872.7 | 847.5 |

`benchmark/conservative_roundtrip_baseline.jl` — 3600×1800 to 360×180
`RasterGrid`, 7,120,800 nonzeros:

| mode | measure | before | after |
|---|---|---|---|
| direct | seconds, mean | 13.026 | 12.305 |
| direct | allocated, mean | 8,288,806,165 | 8,296,036,779 |
| direct | max RSS MiB, mean | 2,733.7 | 2,790.7 |
| chunked | seconds, mean | 12.396 | 12.278 |
| chunked | allocated, mean | 9,111,076,891 | 8,298,191,840 |
| chunked | max RSS MiB, mean | 2,751.4 | 2,723.5 |

The chunked route's allocations fall by 812,885,051 bytes on the round-trip
workload, which is the round trip Phase 0 measured at 810,045,520 bytes, and by
47,913,814 bytes (6.6 %) on the production tree. After the change the chunked
route allocates what the eager route allocates, to within 0.03 %. Eager time is
unchanged: the medians differ by 0.9 % against a 9.5 % spread within the before
state alone.

Peak RSS is the weaker instrument here. On the round-trip workload the
run-to-run spread is 200–320 MiB, so Phase 0's 144 MiB round-trip RSS delta is
not resolvable on this machine; the allocation delta is, and it is the one that
reproduces exactly.

On the production benchmark, over five runs of each state, the chunked mode's
peak drops 27 MiB (3.1 %) and the eager mode's *rises* 24.5 MiB (2.9 %), with
ranges that in both cases barely touch. The eager rise is not live data — the
block is byte-identical and the run allocates 0.25 % *less* — and it is not
resident code: compiling the same eager chain on a 12×6-from-36×18 pair and
stopping there peaks at 677.7 and 650.4 MiB before, 656.1 and 649.5 MiB after.
It is unexplained, it is smaller than the 43 MiB the eager mode's own peak
varies by across runs, and it moves the opposite way from the chunked mode
measured in the same session.

## 4. What Phase 6 needs

- The build seam is method-side and space-agnostic. Nothing on it names a
  destination's array shape or output type, and no `RegridSpace` specializes
  `pairblock`, so O1's generic output application meets one `WeightBlock`
  whatever the destination space is, and adding or generalizing a space adds no
  build path.
- `wholeblock` and `buildblock` are the only builders, and each is two calls
  deep to the seam. A keyword or default that moves into `plan_regrid` changes
  neither.
- A method type is still the only thing that may specialize a build. Anything
  wrapping a method — including instrumentation a phase gate installs — must
  forward `pairblock` to be measured on the path production takes.
