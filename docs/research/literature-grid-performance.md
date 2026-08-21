# Performance audit: the literature grid families

The nine systems added by this branch — ISEA3H, ISEA4H, ISEA4T, rHEALPix,
AusPIX, and the four IVEA/RTEA rhombic systems — measured against the six
already on `main`, per operation, with the fixes that came out of it.

## Method

Each system is measured at the level whose cell count is nearest 200 000, so the
numbers compare work per cell rather than work per level index. Two hundred
sampled cells per operation, best of three sweeps, every result stored into a
sink so that no call is eliminated as dead — without that, the cheap integer
paths measure as 0ns and the picture is wrong in the direction that flatters.
`neighbors` is `Vertex()` unless the column says `Edge()`. The numbers are from
one machine (Apple M-series, Julia 1.12.6) and are meaningful as ratios, not as
absolutes.

The six systems on `main` establish the band a new system should land in.

Every fix below is under the full suite: 1 205 462 assertions pass, none fail,
62 are the `@test_skip`s the branch already carried, 7m00.9s.

## Where the families stand now

| system | level | `cell_centroid` | `cell_boundary` | `node_extent` | `cellat` | `neighbors` V | `neighbors` E | `neighbors` k=2 |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| `IGeo7System` | 5 | 215ns | 1.0µs | 1.1µs | 210ns | 79ns | 75ns | 4.9µs |
| `H3System` | 4 | 137ns | 594ns | 657ns | 179ns | 89ns | 90ns | 369ns |
| `HEALPixSystem` | 7 | 25ns | 362ns | 697ns | 28ns | 103ns | 92ns | 7.6µs |
| `A5System` | 7 | 485ns | 1.7µs | 1.8µs | 1.1µs | 9.6µs | 7.3µs | 42µs |
| `S2System` | 8 | 21ns | 38ns | 65ns | 32ns | 152ns | 92ns | 7.8µs |
| `ISEA4RSystem` | 7 | 149ns | 4.2µs | 4.7µs | 88ns | 336ns | 193ns | 22µs |
| `ISEA3HSystem` | 9 | 270ns | 8.8µs | 9.0µs | **33.5µs** | **205µs** | **198µs** | **1.3ms** |
| `ISEA4HSystem` | 7 | 150ns | 8.7µs | 9.0µs | **29.4µs** | **205µs** | **195µs** | **1.3ms** |
| `ISEA4TSystem` | 7 | 158ns | 3.1µs | 3.3µs | 128ns | 33.3µs | 4.0µs | 382µs |
| `RHEALPixSystem` | 5 | 80ns | 1.0µs | 1.4µs | 54ns | 6.2µs | 89ns | 64µs |
| `AusPIXSystem` | 5 | 115ns | 2.2µs | 1.5µs | 114ns | 6.4µs | 88ns | 66µs |
| `IVEA4RSystem` | 7 | 547ns | **69µs** | 9.5µs | 3.0µs | 6.7µs | 4.0µs | 69µs |
| `IVEA9RSystem` | 5 | 551ns | **71µs** | 9.6µs | 2.9µs | 6.8µs | 4.1µs | 70µs |
| `RTEA4RSystem` | 7 | 617ns | **77µs** | 10.5µs | 3.1µs | 7.2µs | 4.3µs | 74µs |
| `RTEA9RSystem` | 5 | 611ns | **77µs** | 10.6µs | 3.1µs | 7.2µs | 4.3µs | 74µs |

Bold marks a column where a system is more than ten times the slowest of the
six references. ISEA3H and ISEA4H are outside on lookup and on adjacency, which
reads from lookup. The four rhombic systems are outside on `cell_boundary`
alone, because they draw 128 points per cell where every other congruent system
draws 32 — an accuracy choice, not a defect; see the dead end below. ISEA4T,
rHEALPix and AusPIX are inside the band on every operation. Before the fixes,
ISEA4T's `Vertex()` adjacency was 785µs, eighty times the band.

## What was wrong, and what it cost

Six defects, each one repeated work rather than an algorithm that had to be
slow. Before-and-after at the same level, same harness.

### 1. The rhombic inverse projection bisected 58 times

`inverse_triangle` solved `spherical_area(A, B, slerp(B, C, t)) == target` by 58
halvings, and called `slerp` inside the loop — so the arc `BC`'s length and its
sine were recomputed on all 58 of them. The map is smooth and strictly
increasing and `[0, 1]` brackets the root from the start, so Illinois (regula
falsi with the stale endpoint halved) keeps the bracket and converges
superlinearly. The arc is hoisted out.

This is the inner loop of everything geometric on those four systems: one call
per centroid, per boundary vertex, per extent sample.

`cell_centroid`: 2.9µs → 547ns. The new solver is also *more accurate* — at the
level-16 apex probe the area residual is 1.0e-25 against the old 1.6e-18.

### 2. The rhombic node extent sampled twice the perimeter it needed

`node_extent` is a *cap*, not a polygon: it has to contain the subtree, and a
cap that is slightly too wide is conservative rather than wrong. It sampled
eight segments per chart edge. Measured across levels 0–7, halving that to four
widens the cap by 6.7% — a 14% area increase — against halving the inverse
projections at every node of every tree descent. Two segments would widen it by
20% and one by 47%, which is where the lost pruning costs more than the
projections it saved.

`node_extent`: 97µs → 9.6µs, of which 5.2× is the solver above and 2× is this.

### 3. Every counter-clockwise ring sort recomputed its keys per comparison

Four systems wrote `sort!(cells; by = c -> azimuth(cell_centroid(g, c)))`. `by`
runs on every *comparison*, not once per element, so a nine-cell ring built
about twenty-five centroids — each an inverse projection — to order nine. Keys
are now built once and `sortperm` orders the permutation.

### 4. rHEALPix re-derived its seam topology on every crossing

`edge_neighbors` is exact integer arithmetic on the ternary lattice until a move
leaves its root square, and then it called `_seam_metadata`, which derives the
cut's target square, edge and orientation *from the projection* — four forward
projections and two inverse ones. `vertex_neighbors` makes nine such moves.

The topology depends on the two polar-square placements alone: the derivation
round-trips plane → sphere → plane and the longitude origin cancels. That is 16
placements × 6 roots × 4 directions = 384 entries, now built once at load and
checked against the derivation for all 16 placements at three longitude origins.

`edge_neighbors`: 146ns and allocation-free.

### 5. ISEA4T rebuilt each candidate's boundary once per vertex of the subject's

The vertex star filtered a three-step edge-BFS disc by chart-vertex incidence:

```julia
out = [n for n in candidates if any(_same_vertex(a, b)
    for a in ring for b in DGG.cell_boundary(g, n))]
```

In a *flattened* generator (`for a in … for b in …`) the inner iterable is a
closure over `a` and is rebuilt for every `a`. `ring` has 24 points, so each
candidate's densified 24-point boundary was constructed 24 times — and vertex
incidence is a statement about the three *corners*, which is 3 inverse
projections rather than 24. Corners now, once per candidate.

`neighbors` V: 785µs → 33.3µs. `neighbors` k=2: 10.3ms → 382µs.

### 6. The ISEA3H/4H inverse re-walked every prefix it extended

`_locate_in_root` descends the central-place prefix tree keeping the 24 closest
candidates per digit, and scored each by rebuilding its whole development-plane
offset with `_hex_dev` — which allocated a digit vector and summed from digit
one. The offset is a *prefix sum*, so extending a candidate is one term. The
frontier now carries the running sum, `_hex_dev` is allocation-free, the
direction and radial-scale tables are precomputed, and the two `atan2`s per
candidate are down to about half of one (the query's own angle is hoisted out of
the descent, and the wedge fold settles the upper half-plane on the sign of the
imaginary part).

`cellat`: 150µs → 33.5µs (ISEA3H level 9), 114µs → 29.4µs (ISEA4H level 7).
`neighbors` follows it, being six `cellat` calls: 832µs → 205µs.

## The dead end: level-aware boundary densification

The rhombic `cell_boundary` chords every chart edge with 32 great-circle
segments at every level — 128 points per cell, where HEALPix, ISEA4R and
rHEALPix all draw 32. The obvious economy is that the deviation of a chord from
the edge it approximates is quadratic in the chart *step* `1/(nseg · nside)`, so
accuracy is a property of the product: hold `nseg · nside >= 32` and a level-7
cell needs one segment per edge rather than 32, at the same 2e-5 rad bound the
fixed count buys at level 0. That is a **169× speedup**, 372µs to 2.2µs, and it
is wrong.

It fails on MIXED-LEVEL sets, which is what `MultiOrderCoverage` produces. A
parent's edge and the child edges along it are then chorded at different steps,
and the lens between the two polylines is in neither polygon. On a coverage of
California it left five sampled interior points in no cell at all — two on
IVEA4R, three on RTEA4R — and every one of the five sat where a level-6 cell
abuts level-7 ones, not at the coverage rim. These systems are congruent, four
children tile their parent, and `multiorder_polygons.jl` pins that tiling at
zero slivers. It is a broken law, not a tolerance.

Nesting the polylines is the fix and cannot be afforded. A parent's samples
contain its children's exactly when `nseg(l) == s · nseg(l+1)`, which is the
same condition as `nseg · nside` constant — so the scheme above *is* nested,
right up until `nseg` hits its floor of 1 and stops halving. Below that floor
the parent samples a midpoint its children do not, and pushing the floor to the
bottom of the tree means `s^(maxlevel − l)` segments at level `l`: 2^25 of them
at a root.

So the count has to be level-independent, and the only question is what it
should be. At a constant `nseg` the parent/child mismatch is
`~0.02/(nseg · nside)²` and shrinks quadratically with depth on its own; 32
puts it at 4.7e-9 rad by level 6. Eight — ISEA4R's density, and ISEA4R passes
the same law — also leaves zero slivers on this fixture, and takes
`cell_boundary` from 70µs to 17.7µs. It costs an order of magnitude of edge
accuracy at level 0, from 2e-5 rad to 3.1e-4, which is the bound the
implementer chose deliberately. That is a call about the geometry contract
rather than an optimisation, so 32 stands and the number is recorded here.

### Effect on the suite

Each suite `include`d on its own in a fresh process, before and after measured
the same way, so the numbers carry that process's compilation with them.

| suite | before | after |
|:--|--:|--:|
| `test/systems/ISEAGrids/` | 44.8s | 18.3s |
| `test/systems/IVEARTEA/` | 16.2s | 9.0s |
| `test/systems/RHEALPix/` | 4.6s | 5.0s |
| `test/systems/crosssystem/stencils.jl` | 120.9s | 68.7s |

`stencils.jl` does strictly more work after than before: the member-adjacency
sweep grew from six systems to all fifteen.

rHEALPix is the one that does not move, and the reason is worth stating rather
than rounding away. Its suite is dominated by the projection oracle and by the
SUID codec, neither of which touches `edge_neighbors`; the seam table shows up
in the operation (146ns, allocation-free) and not in the file. A suite timing
measures the suite.

Whole suite after: 1 205 462 assertions, no failures, 7m00.9s.

## What is still slow, and why

Three residues, none of them repeated work. Each is an architecture, and each
would need a different implementation rather than a fix.

**ISEA3H/4H `cellat` is a search, not an inverse.** The central-place hierarchy
has no closed-form point-to-cell, so `_locate_in_root` descends the prefix tree
keeping 24 candidates per digit, over three roots per face. That is
`O(level · 24 · aperture)` scored candidates for one lookup, and `neighbors`
costs six of them because a neighbour is found by locating the point a lattice
step away.

Closing it means working in digit space instead of in the plane, and the two
systems are not in the same position for that. ISEA3H already carries the
modified-balanced-ternary direction digits of Sahr (2008) — `hex.jl` says so at
the forward construction, and `isea-family.md` cites the paper as the primary
indexing source — so the digits are there and only the *inverse* and the
*step* are missing. Whether MBT admits a closed-form neighbour step the way
Generalized Balanced Ternary does for the aperture-7 hierarchy in `gbt.jl` is a
question for that paper and for DGGAL's `RI3H.ec`, which the branch already
treats as a BSD reference; it is not something this audit established. ISEA4H
has no such digit algebra on record here at all — `hex.jl` calls A4 the ordinary
class-I triangular lattice, which points at lattice coordinates rather than at a
named arithmetic.

**ISEA3H/4H `node_extent` is 9µs and `cap_inflation` was 3.5–4.5.**
Neither system overrides `node_extent`, so the 9µs is `cell_boundary`
(8.7–8.8µs above) reduced to a cap and scaled by the constant. Both are
central-place refinements and a cell's children genuinely lie outside its
polygon, so the cap does have to inflate past the cell — that much the geometry
forces. How far it has to inflate is a measurement: the smallest `f` with
`spherical_distance(centres) + r_desc <= f · r_anc` over every
ancestor/descendant pair, exhaustive at the low levels, is ~2.0 for ISEA3H and
~1.68 for ISEA4H. The declared 4.5 and 3.5 are more than double that. They are a
magnitude-sum limit of the per-level overhang series — a triangle inequality
that assumes every level's displacement points the same way — where the
central-place digit directions rotate between levels and mostly cancel. Tree
pruning is correspondingly poor, which is why these two are the slowest systems
in the multi-order coverage tests even after the fixes, and more than half of
that cap radius buys nothing. The oversize is also the only reason the ISEA3H
and ISEA4H conformance calls pass `require_convex_extents=false`: at 4.5 and 3.5
an extent exceeds a hemisphere at ISEA3H levels 0–1 and ISEA4H level 0, and
nowhere else. Non-convexity is the declared constant's doing, not the
refinement's.

The bound turned out to be exactly derivable, so the constants are now **2.3 and
2.0** rather than a measurement plus a guess. The covering factor of a descendant
`n` levels down is `sqrt(3)·|Σ q^j u_j| + q^n`; bounding `|Σ|` by `Σ|·|` is where
`3.3661` and `2.7321` came from, and it assumes every step points the same way.
Aperture 4's digit directions are constant, so its worst case really is a single
direction repeated — but that gives `sqrt(3) = 1.7321`, and the old `1 + sqrt(3)`
had added a whole ancestor circumradius for a descendant that shrinks to nothing.
Aperture 3's directions rotate: odd depths use a fixed gauge, even depths are
gauged by the preceding digit. Taking the support function of the reachable
displacement set — a DP over `(depth, previous digit)`, swept over direction with
a `sec(π/M)` enclosure and a geometric tail bound — gives a supremum of exactly
**2**, approached by the alternating sequence `1,2,1,2,…` and only from an
even-depth ancestor; odd-depth ancestors reach only `sqrt(3)`. Brute force over
every admissible sequence to depth 13 gives 1.99942 with 0.0023 of tail left, and
the model is checked digit-by-digit against `_dev_step` itself.

The sphere adds about 1%, not the 59% the chart's worst-case anisotropy (1.5864)
would suggest: an ancestor and its subtree share one local Jacobian, so the
distortion very nearly cancels in a ratio instead of multiplying it. Exhaustively
over every ancestor at levels 0–4 with seven-level subtrees, the sphere-side
worst is **2.0211** for ISEA3H — on the odd ancestor levels, as predicted — and
**1.721** for ISEA4H over levels 0–3. The new factors carry 14% and 16% over
those, put the widest extent at 86° and 75°, and let both systems drop
`require_convex_extents=false` and take the harness's convexity assertions like
every other system. `test/systems/ISEAGrids/runtests.jl` reproduces the
measurement, which nothing did before — that is how 4.5 went unchallenged.

**ISEA4T's vertex star is geometric.** 33µs against 4µs for the edge star: the
edge star crosses each chart edge at its inverse-projected midpoint and asks
`cellat`, and the vertex star then expands a three-step edge-BFS and filters it
by corner incidence. Both are correct and neither needs a face-adjacency table,
which is what they were written to avoid. The analytic alternative is
path-arithmetic neighbours on the triangle quadtree plus an icosahedral
face-adjacency table, and it would make the vertex star as cheap as the edge one.

**The rhombic systems have no `descendant_range`.** Their canonical order is
row-major within a root, so a subtree is `s` rows of `s` cells and never a
contiguous run — `has_sorted_subtrees` is false and stays false. The consequence
is not speed of `neighbors` but of everything range-based: `subtree_border` walks
the whole subtree instead of its rim, `member_neighbors` builds a per-call member
dictionary instead of binary-searching curve keys, and `CellVector` stores
positions instead of windows.

This one is available: an aperture-4 rhombic system in **Morton order** has
contiguous descendant ranges and can reuse `SquareRimEngine`/
`SquareInteriorEngine` under `MortonCurve()` exactly as `ISEA4RSystem` does. The
sealed DGGAL oracle vectors pin only level 0, where `nside == 1` and row-major
and Morton coincide, so the reindexing is oracle-safe. It would change canonical
cell identity for IVEA4R and RTEA4R, which is a decision about the published
index rather than an optimisation, and the aperture-9 pair would need a nonary
curve and a 3×3 generalisation of the two square engines before the family could
be consistent.

**None of the nine has a specialized halo engine.** Found when main's subtree
halo work was merged in, which is why it is not in the table above: the shipped
walk on these systems *is* the generic one, and the generic one costs what the
full-grid reference scan costs. At depth 7 from a level-0 root, collecting the
halo takes 4.46s on rHEALPix against 4.46s for the brute-force scan it is
supposed to beat — while S2, which has an engine, does the same job in 0.006s
against its own 0.72s scan. The cost is structural rather than per-call, so it
lands on every halo law at once: the cross-system halo suite went from minutes
to 222 of the 230-minute merge run before its depth ladders were sized by
subtree size rather than by level number.

The four aperture-9 systems are worst, and for a compounding reason. Level 7 is
a 16k-cell subtree at aperture 4 and a 4.8M-cell one at aperture 9, and IVEA9R
and RTEA9R also lack `descendant_range`, so the subset walk cannot prune by
range *and* has no engine. One complete subset walk asks 535,818 questions on
IVEA9R where rHEALPix — same aperture, same subset size, but with sorted
subtrees — asks 9,612, and S2 asks 1,572. Those counts are now pinned in
`test/systems/crosssystem/subtree_halos.jl`, so closing the gap shows up as a
failing pin rather than as an unnoticed improvement.

The two square families already have the engines to copy: `SquareRimEngine` and
`SquareInteriorEngine` for the aperture-4 rhombic systems under the Morton
reindexing described above, and rHEALPix's own square blocks for the aperture-9
ones once a nonary curve exists. ISEA3H/4H would need a hex-arc engine of the
kind `HexArcHaloEngine` already provides at aperture 7.

## Recommended order of further work

1. **Digit-space adjacency for ISEA3H/4H.** Largest single win: it removes
   `cellat` from the neighbour path entirely and takes `neighbors` from 205µs
   toward the sub-microsecond band. Start by reading Sahr (2008) and DGGAL's
   `RI3H.ec` to settle whether ISEA3H's MBT digits support a closed-form step —
   that is the open question, not a settled result — and expect ISEA4H to need
   its own answer.
2. **Path-arithmetic vertex star for ISEA4T.** Removes the BFS-and-filter, ~8×
   on `neighbors` V.
3. **A decision on the rhombic boundary density.** Eight segments per edge
   instead of 32 is 4× on `cell_boundary` — 70µs to 17.7µs — and still tiles
   the mixed-level coverage exactly. It is ISEA4R's density, and it gives up
   the 2e-5 rad level-0 bound. A call for whoever owns the geometry contract;
   see the dead end above for the measurement.
4. **Morton reindexing for IVEA4R/RTEA4R**, with a nonary curve and a generalised
   square engine for the aperture-9 pair. Unlocks `descendant_range` and
   everything built on it. Needs a decision on canonical identity first.
5. **An interior fast path for rHEALPix `vertex_neighbors`.** A cell whose row
   and column are both interior to its root has its eight neighbours at fixed
   lattice offsets in a fixed rotational order, so the geometric sort — eight
   centroids — is only needed on the ~1.6% of cells that touch a seam.
