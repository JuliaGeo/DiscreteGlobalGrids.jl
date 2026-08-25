# G3 — dependency graph identity and row views

- Date: 2026-08-23
- Card: `regrid-notes/2026-08-21-regridding-simplification-plan.md`, Phase 2,
  Task G3.
- Branch: `claude/g3-graph-identity`, cut from `claude/g1-graph-oracles`
  @ `c13eca9` (PR #70's head). Stacked PR into `claude/g1-graph-oracles`.
- Commit: `Add reusable dependency graph identities`.

The card asked for four things: stamp the graph with chunk counts, support
radius, space identity and a serializable narrow-phase tag; add validated
`restrict(graph, destination)` row views that share the parent relation; retain
every existing accessor and the bidirectional CSR; and keep Graphs.jl
compatibility over the same CSR without a second graph type. All four landed.

---

## 1. The identity

Three new types in `lib/GlobalRegridding/src/chunkgraph.jl`, all serializable and
none of them holding a space, a closure, or a type parameter.

```julia
struct SpaceStamp                    struct DependencyIdentity
    tag::Symbol                          dst::SpaceStamp
    ncells::Int64                        src::SpaceStamp
    nchunks::Int64                       radius::Float64
    digest::UInt64                       narrow::Symbol
end                                  end
```

`ChunkDependencyGraph` gained one `id::DependencyIdentity` field and lost its
loose `radius::Float64`; `dependency_radius(g)` now reads `g.id.radius` and is
otherwise unchanged. The whole record is **120 bytes**, constant, on every graph
in the matrix — including the 3.35 MB production one.

`spacestamp(space)` digests the *name* of the space's type, `ncells`, `nchunks`,
and every covering cap `chunkextents` reports, in chunk order.

### What invalid reuse this catches

`validate_dependencies(graph, dst_space, src_space; radius, narrow,
destinations)` returns the graph or throws an `ArgumentError` naming the
mismatch. Five checks, and the failure is at the point of adoption, not at read
time:

| what moves | caught by | why it must fail |
|---|---|---|
| a different destination grid, resolution, extent or chunking | `dst` stamp | different rows entirely |
| the same for the source | `src` stamp | different columns entirely |
| source and destination swapped | stamps compared in role | the relation is not symmetric |
| a space of a different type | stamp `tag` | different index, different relation |
| a method with wider support | `radius <= g.id.radius` | the relation is monotone only the other way |
| a narrow phase the caller did not ask for | `narrow` equality | the caller would silently get a thinned relation |
| the full relation where a narrow phase was expected | `narrow` equality | the caller would silently get a wider one |
| a graph narrowed by an anonymous closure | `UNNAMED_NARROW` rejected unconditionally | nothing records which pairs it dropped |
| a row view offered where a whole-space graph is needed | `isrestricted` + row count | a plan would compute a fraction of its destination |
| a row view over the wrong rows | `globaldestinations(g) == destinations` | the column's rows are not these rows |

The narrow-phase tag is a *claim about the relation*, so both halves must agree
and both mismatches are errors at construction: `narrow` without `refine`, and
`refine` with `narrow = :none`. The unsupplied value is `nothing`, not `:none`,
precisely so those are distinguishable. A `refine` given with no tag is stamped
`:unnamed`, which never validates in either direction — an anonymous closure has
no identity, so a relation it narrowed cannot be certified as the relation any
later caller wants. Production's own narrow phase is now named
`:copdem_tile_lonlat_box` at `scripts/copdem_production.jl:764`, which is what
makes a refined production graph reusable at all.

### What it does NOT catch — stated plainly

A fingerprint is not a proof, and the card asked for this to be said rather than
papered over. Equal stamps still permit:

1. **A 64-bit digest collision.** Unlikely; silent when it happens.
2. **Different type parameters behind the same type name.** Only `nameof` is
   stamped. Rendering the full parametric type of a `RasterGrid` costs
   **[measured] 0.79 ms** — more than building the whole relation on a small
   case — because `show(::Type)` searches every module for an alias. In practice
   a parameter that changes the relation changes a cap or a count too, but the
   parameters themselves are unchecked.
3. **Identical caps over different cell geometry.** This one is a real hole and
   it is worse on the source side than the destination side. The destination
   half of the relation is a function of the caps alone, so equal destination
   caps really do give equal rows. The source half is not: it comes from
   `chunkindex`, a native hierarchy that need not test the caps `chunkextents`
   reports — the exact divergence PR #69 exists to fix. Two source spaces with
   equal chunk caps and different hierarchies can produce different relations and
   identical stamps.
4. **In-place mutation of a space after the graph was built.** A stamp is a
   snapshot, not a live binding.

Nothing cheaper is sound and nothing sound is cheap: the only exact check is to
rebuild the relation and compare it, which is the work the identity exists to
avoid. What was rejected, and why:

- **Object identity (`===` on the stored spaces).** Exact, but not serializable,
  and it makes the graph hold the spaces alive — the thing the type deliberately
  avoids. A `WeakRef` fallback would make validation depend on when the GC ran.
- **`objectid(typeof(space))` folded into the digest.** Cheap and more
  discriminating than `nameof`, but not guaranteed stable across sessions, which
  would quietly break "write the identity down beside the relation".

## 2. `restrict`, and how it shares the parent CSR

`restrict(g, destinations)` lives in `chunkgraph.jl` beside the graph and returns
a `ChunkDependencyGraph` — the same type, no second graph. The struct grew one
field:

```julia
dstrows::Vector{Int}   # local row -> global destination chunk; EMPTY == identity
```

Empty means the identity map, which is every whole-space graph, so a graph built
over a whole destination space allocates nothing for it and pays one predictable
branch in `_row`. Restriction composes: because a whole-space graph's row `d`
*is* destination chunk `d`, one vector serves as both the index into the shared
offsets and the global chunk number, at any nesting depth.

- **Shared, by reference, never copied:** `dstoff` (the destination-major
  offsets) and `srcof` (the entire edge array). A test asserts `view.dstoff ===
  g.dstoff` and `view.srcof === g.srcof`; if that ever becomes a copy the view
  has stopped being a view.
- **Rebuilt, and it has to be:** `srcoff`/`dstof`, the source-major direction. A
  view's `consumersof` must count only the rows the view holds, or a refcount
  taken from it retires nothing. That is a counting-sort over the selected rows'
  edges, `O(selected edges + nsourcechunks)`, with no spatial index, no cap test
  and no query of either space. `_transpose` was generalized to take a row
  selection; the whole-space call passes `OneTo(ndst)` and is the same code.

Global destination identity is retained through `globaldestinations(g)`,
`globaldestination(g, d)` and `localdestination(g, chunk)`, and through the
identity record, which still stamps the **whole** destination space — a view
knows what it is a view *of*. Rows must be strictly ascending; a walk order is a
separate permutation applied by whoever walks the rows, exactly as it already is
for a whole-space graph.

Everything else is retained: `sourcesof`, `consumersof`, `sourcedegree`,
`consumerdegree`, the vertex/chunk converters, and the full Graphs.jl interface
over the same CSR. `Graphs.ne` now reads the private source-major array rather
than `length(srcof)`, because on a view `srcof` is the parent's whole edge array.
The suite runs `is_bipartite`, `connected_components`, `edges`, `has_edge` and
the neighbour sum against a row view.

## 3. Measured: `restrict` against a rebuild

`benchmark/chunk_graph_gates.jl` gained a third question and a
`restrict_measurement` block. The rebuild arm does every piece of work
`chunk_dependency_graph` does for a column-sized destination — stamp both spaces,
take the destination caps, build the source chunk index, one `candidatechunks!`
per row, assemble the CSR, transpose it — inside the timed region, and hands back
a real `ChunkDependencyGraph`. Both sides therefore produce the same kind of
object and pay the same `O(nsourcechunks)` transpose. It asserts they agree row
for row **and** in both CSR directions before reporting a ratio.

It understates the rebuild: it does not construct the per-column destination
space, which a real per-column rebuild also pays.

Julia 1.12.6, `-t 8 --gcthreads=4`, 5 samples, medians. Block = ⌈ndst/16⌉
contiguous destinations.

| case | ndst | nsrc | restrict 1 | rebuild 1 | × | block | restrict blk | rebuild blk | × |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| raster-small | 72 | 24 | 2 µs | 143 µs | **91×** | 5 | 2 µs | 163 µs | **77×** |
| raster-support | 72 | 24 | 2 µs | 134 µs | 79× | 5 | 2 µs | 151 µs | 66× |
| raster-nonuniform | 72 | 49 | 2 µs | 201 µs | **104×** | 5 | 2 µs | 201 µs | 89× |
| dgg-complete | 72 | 492 | 5 µs | 197 µs | 40× | 5 | 5 µs | 330 µs | 67× |
| dgg-crosssystem | 72 | 24 | 2 µs | 81 µs | 46× | 5 | 2 µs | 77 µs | 33× |
| dgg-rooted | 7 | 72 | 2 µs | 99 µs | 56× | 1 | 2 µs | 110 µs | 54× |
| dgg-sparse | 175 | 24 | 2 µs | 83 µs | 43× | 11 | 2 µs | 131 µs | 55× |
| polar-source | 72 | 15 | 2 µs | 101 µs | 54× | 5 | 2 µs | 209 µs | 111× |
| antimeridian-source | 72 | 12 | 1 µs | 73 µs | 53× | 5 | 2 µs | 74 µs | 42× |
| raster-4320-162chunks | 3 432 | 162 | 2 µs | 243 µs | 102× | 215 | 6 µs | 2 229 µs | **370×** |
| raster-4320-1800chunks | 3 432 | 1 800 | 32 µs | 1 101 µs | 34× | 215 | 39 µs | 6 275 µs | 159× |
| dgg-large | 3 432 | 492 | 9 µs | 181 µs | 19× | 215 | 18 µs | 2 491 µs | 141× |
| **copdem90-igeo7-l12** | 66 175 | 26 475 | **36 µs** | **959 µs** | **26×** | 4 136 | **140 µs** | **11 244 µs** | **80×** |

Allocation, production pair: 424 KB against 541 KB for one destination, and
**518 KB against 8.52 MB** — 16× less — for the 4 136-chunk column.

### Two corrections to the card's premise

The card said a per-column rebuild "now costs a full spatial index build". That
is true only on the **generic** path, and **neither shipped native space takes
it**: `chunkindex(::DGGSpace) = space` (a field read) and
`chunkindex(::RasterGrid)` is a cursor. `source_index_seconds` measures below
0.5 µs in every case in the matrix. The saving is real and large, but it comes
from the destination caps and the per-row queries, not from an R-tree build.

Second, the `O(nsourcechunks)` transpose term is not free. Restricting to a
**single** destination on a 26 475-source problem allocates a 26 476-entry offset
array; against a bare re-query of one row (24 µs, no identity, no CSR) `restrict`
at 39 µs is *slower*. Restrict a column, not a row. Documented on `restrict`.

## 4. Cost of the identity

Stamping is inside `chunk_dependency_graph`, so it is a component of the graph
timing rather than something beside it. `identity_seconds` is a new ndjson
column.

| case | graph before | graph after | Δ | identity | identity share | alloc Δ | graph object Δ |
|---|---:|---:|---:|---:|---:|---:|---:|
| raster-small | 0.34 ms | 0.42 ms | +23% | 20 µs | 4.8% | +1 072 B | +120 B |
| raster-nonuniform | 0.71 ms | 0.70 ms | −1% | 32 µs | 4.5% | +2 096 B | +120 B |
| dgg-complete | 1.37 ms | 1.35 ms | −1% | 9 µs | 0.6% | +32 B | +120 B |
| dgg-sparse | 0.27 ms | 0.36 ms | +34% | 4 µs | 1.1% | +32 B | +120 B |
| polar-source | 0.19 ms | 0.17 ms | −9% | 13 µs | 7.8% | +688 B | +120 B |
| raster-4320-162chunks | 34.1 ms | 34.3 ms | +1% | 118 µs | 0.3% | +6 632 B | +120 B |
| raster-4320-1800chunks | 62.7 ms | 63.7 ms | +2% | 824 µs | 1.3% | +72 104 B | +120 B |
| dgg-large | 32.7 ms | 30.6 ms | −6% | 37 µs | 0.1% | +32 B | +120 B |
| **copdem90-igeo7-l12** | **115.3 ms** | **116.7 ms** | **+1.2%** | **758 µs** | **0.6%** | **+32 B** | **+120 B** |

The identity is 0.1–9.6% of build time, median ~1.3%; the large percentage deltas
are on cases whose absolute time is a few hundred microseconds and whose
run-to-run spread already exceeds the change (G1 recorded 0.1219 s for the
production case on this same commit; both runs here are 0.115/0.117). The graph
object grows by a constant 120 bytes, not proportionally.

Allocation grows by **32 bytes** wherever the source is a `DGGSpace`, because
`chunkextents` there is a stored field. It grows proportionally to the source
chunk count wherever the source is a `RasterGrid` — up to 72 KB at 1 800 chunks —
because stamping materializes `chunkextents(src)`, which the indexed builder
otherwise never needs. This is the one place stamping does real work; it is
1.7% of the build on that case.

## 5. Gates: before and after

`benchmark/chunk_graph_gates.jl`, 13 cases including the production
GLO-90 × IGeo7-L12 pair with a local tile list, both arms, 26 ndjson rows.

**The relation is identical, row for row.** Across all 26 rows, `edges`,
`demanded_pairs`, `demand_missing`, `oracle_pairs`, `oracle_missing`,
`only_here` and `missing_here` are unchanged — including the production pair's
326 064 indexed edges, 326 386 latjoin edges, and the 72 / 394 crossing.

Both runs print:

```
cases run: 13.  Oracle-checked: 9.  Oracle skipped: 4.
verdict: PASS on 9 oracle-checked case(s); 4 case(s) unchecked
```

That is a real geometric verdict on nine cases, not a `NOT CHECKED` run.

The `:latjoin` arm now stamps the same identity the production builder does, so
the two arms stay comparable as objects and not only as relations; the stamping
cost is inside the timed region for both.

## 6. Suites

| suite | before (`c13eca9`) | after | delta |
|---|---|---|---|
| `lib/GlobalRegridding/test` | 3 737 pass / 1 broken / 0 fail | **3 828 pass / 1 broken / 0 fail** | +91 pass |
| `test/systems/crosssystem/regrid.jl` | 202 pass / 0 fail | **202 pass / 0 fail** | none; output byte-identical |
| `test/scripts/copdem_policy.jl` | not separately baselined | **87 pass / 0 fail** | see below |
| `test/scripts/copdem_source_mode.jl` | not separately baselined | **7 pass / 0 fail** | none |

Every delta explained: the +91 are the six new testsets in
`lib/GlobalRegridding/test/test_chunkgraph.jl` — space stamps, graph identity,
the narrow-phase tag, invalid reuse, row views sharing the parent relation, and a
row view checked against the cell-geometry oracle. Nothing existing changed
count, and the one pre-existing broken test is still broken.

`test/scripts/copdem_policy.jl` needed a one-line edit: it builds a graph from
raw caps through the private `_chunkgraph`, which now takes an identity. It
passes `DependencyIdentity()`, the empty identity, which matches no space and
which `validate_dependencies` refuses to certify — the honest record for a
relation built by hand. Its assertion count is unchanged (the seam testset's six
assertions are the ones that were already there); it was not run before the
change, so it is reported as an after-state, not a delta.

## 7. Residual uncertainty

- **The source-side fingerprint hole (§1.3) is real and unfixed.** Equal source
  chunk caps do not imply equal source relations, because the relation comes from
  the hierarchy rather than the caps. Closing it needs an identity for the
  *index*, which no space exposes today. Adding a `chunkindexstamp` hook to the
  qualified space interface would close it; that is an interface change and
  outside this card.
- Only `nameof` is stamped, so type parameters are unchecked (§1.2).
- The production case's graph timings move ±5% run to run at t8 on a shared box.
  The +1.2% recorded for the identity is inside that band; the identity's own
  758 µs (0.6%) is the number to trust.
- `restrict` costs `O(nsourcechunks)` even for one row (§3). No caller in tree
  restricts to a single row today, and G4's per-column plans will not.
- Nothing yet *calls* `validate_dependencies` or `restrict` in production. G4
  owns that: `ChunkedPlan` becomes the sole graph owner, `refine` is confined to
  lazy `plan_regrid`, and per-column plans take validated row views instead of
  rebuilding. This card deliberately stopped at the mechanism.
- `benchmark/regridding_plan_baseline.jl` was not re-measured; it times
  `chunk_dependency_graph` on the production pair and will read ~1% higher.
