# The M2 baseline: what `Centroid()` costs a sweep today, and what a perfect
# cache could give back.
#
#     julia -t auto --project=benchmark benchmark/needs_centroid.jl
#
# `needs = (Value(dem), Centroid())` computed `cell_centroid` once per touch
# before the M2a working set landed: once for the visited cell and once for
# each of its clipped neighbours, so `1 + degree` computations per cell where a
# per-sweep working set would want one. Three requests bracket that on the same
# cells, with the same kernel wherever the kernel can be the same:
#
#   1. `(Value(dem),)`                  the value-only floor: the sweep, the
#                                       clip and the rings, no geometry.
#   2. `(Value(dem), Centroid())`       today's behaviour, the thing to beat.
#   3. `(Value(dem), Value(table))`     the same centroids read from a full
#                                       precomputed table -- the ceiling a
#                                       perfect cache reaches, since a cache
#                                       cannot be cheaper than an array read.
#
# Requests 2 and 3 run the *same* kernel over the *same* numbers (asserted
# equal, elementwise), so `2 - 3` is centroid production alone. `2 - 1` also
# carries the geometry arithmetic the value-only kernel does not do, so it is
# the upper bound on what a cache could address; read the two together.
#
# A tight `cell_centroid` loop over the same cells prices one call in
# isolation, which is what the `1 + degree` arithmetic multiplies.
#
# IGeo7's aperture is 7, so a rooted subtree holds exactly 7^depth cells and
# the reachable sizes step 117,649 -> 823,543 -> 5,764,801: nothing lands
# inside 1e6..2e6, and 823,543 is the closest rung to it.
#
# Measured 2026-08-24, M-series macOS, Julia 1.12.6, 8 threads, best of 3:
#
#   request                       cells      seq s     thr s      alloc B
#   (Value,)                    117,649    0.02071   0.00319      950,336
#   (Value, Centroid)           117,649    0.10814   0.01508      950,336
#   (Value, Value(table))       117,649    0.02803   0.00530      950,336
#   (Value,)                    823,543    0.14990   0.02219    6,602,816
#   (Value, Centroid)           823,543    0.78296   0.10977    6,602,816
#   (Value, Value(table))       823,543    0.20446   0.03025    6,602,816
#
# All three requests allocate the same bytes -- the output vector and the
# sweep's own overhead -- so `needs` itself allocates nothing and the centroid
# ring stays on the stack. What `Centroid()` costs is time, not memory:
# 80.8% / 80.9% of the whole sweep, 743 / 769 ns per cell, against 115.5 /
# 119.2 ns for one isolated `cell_centroid` (x 7 = 808 / 834 ns, so the sweep
# pays the isolated call per touch and, in place, a little less).
#
# A prebuilt table removes 91.6% / 91.4% of that -- the whole sweep runs
# 3.86x / 3.83x faster, saving 97.3 / 100.4 ns per touch, which is again one
# `cell_centroid` and nothing on top of it. Pay the table's own build out of
# the same budget -- which a per-sweep working set must, since it still
# computes one centroid per cell -- and it removes 76.4% / 76.0%, just under
# the 6/7 = 85.7% the `1 + degree -> 1` model predicts; the shortfall is the
# ring read that replaces the call. The table itself is 2.8 MB / 19.8 MB, and
# that O(ncells) footprint is the whole reason M2 wants a bounded window
# instead.
#
# Read ratios, not absolutes: the table above is one run, and across five runs
# of this file every absolute moved by up to ~25% together (the threaded rows
# most) while the derived ratios stayed put -- centroid share of the
# `Centroid()` sweep 77.7-85.2%, the prebuilt table's share of that 89.9-94.0%,
# one isolated `cell_centroid` 115.5-145.0 ns.
#
# After the M2a centroid working set, measured 2026-08-24. `Centroid()` is now
# answered from a per-task direct-mapped window of 16,384 slots keyed by local
# index, so the sweep computes 1.026 / 1.037 centroids per cell instead of the
# 6.96 / 6.98 it touches. Both columns below are single runs of this file in
# one session, best of 3 per cell, before and after on the same machine. Rows 1
# and 3 are unchanged code, and they are what says how that session's machine
# compared to the baseline capture above: the floor and the ceiling had moved
# +20.7% and +30.2% at 117,649 cells, +4.4% and +2.7% at 823,543. Read the
# pairs, not the absolutes:
#
#   request                       cells   before ms    after ms        alloc B
#   (Value,)                    117,649       24.99       21.78         950,336
#   (Value, Centroid)           117,649      132.42       48.58       1,474,752
#   (Value, Value(table))       117,649       36.49       29.52         950,336
#   (Value,)                    823,543      156.51      153.52       6,602,816
#   (Value, Centroid)           823,543      823.39      348.25       7,127,232
#   (Value, Value(table))       823,543      210.05      210.81       6,602,816
#
# The surcharge over the value-only floor falls 107.43 -> 26.80 ms and
# 666.88 -> 194.73 ms: 75.1% and 70.8% of it removed, against the 6/7 = 85.7%
# a cache that never missed and cost nothing to consult would reach. Threaded,
# reported not gated: `(Value, Centroid)` 33.79 -> 9.07 ms and 124.07 ->
# 51.13 ms. The 524,416 B the `Centroid()` rows now carry is the window itself
# -- 16,384 tags and 16,384 points -- allocated once per task, not per cell;
# the other rows still allocate the output alone.
#
# The M2a gate was two absolute thresholds -- <= 47 ms at 117,649 cells and
# <= 340 ms at 823,543 -- and this run missed both: 48.58 ms is 3.4% over, and
# 348.25 ms is 2.4% over. Those thresholds encoded "at least 70% of the
# surcharge removed" at the baseline machine's speed, and the same-session
# before/after pair above -- whose unchanged floor and ceiling rows had moved
# +20.7%/+30.2% and +4.4%/+2.7% from that capture -- removes 75.1% and 70.8%.
# The criterion is met; the absolute numbers it was written as are not, and
# the ruling accepted it on the criterion (design note, "M2a outcome").
#
# The window is keyed by the local index, so what it is worth is a property of
# the visit order, and the two shuffled rows at the small size price that: a
# fixed random permutation costs `(Value, Centroid)` 108.79 ms against a
# same-order `(Value, Value(table))` control at 30.50 ms, a surcharge of
# 78.29 ms where storage order pays 18.67 ms -- 4.19x, and 4.43x on the run
# before it (both 2026-08-24, this file). Storage order and any permutation
# that preserves locality keep the hit rate; a random one does not.
#
# Where the remaining 26.80 ms goes, per touch of the 6.963 the small case
# makes: 9.4 ns builds the ring and reads a value into it, which the prebuilt
# table pays too; 18.1 ns is the 14.7% of touches that miss the window, times
# a 123.1 ns `cell_centroid`; ~5 ns is the probe itself, over what reading a
# plain array costs. A wider window did not buy that back -- W = 1024, 4096,
# 16384 and 65536 all landed inside run-to-run noise, measured during
# implementation and not reproducible from this file -- because the misses a
# bigger one removes (1.080 -> 1.004 centroids per cell) are worth under a
# millisecond here. W is a hit-rate choice rather than a time one.
#
# After `Centroid()` became a CELL FIELD, measured 2026-08-24. `Centroid()` now
# resolves to `Value(cellfield(cell_centroid, cv))`; the window is the field's
# per-task reader, and it is handed the cell the sweep already decoded, so a
# miss computes without decoding it again. Two rows price the two spellings of
# a field the caller builds itself:
#
#   4. `(Value, Value(field, full))`    a field whose `known` covers the whole
#                                       collection -- read straight through,
#                                       with no window at all, so this is row
#                                       3 with one wrapper on it.
#   5. `(Value, Value(field, border))`  a field whose `known` is the border
#                                       alone (2,184 of 117,649 cells; 6,558
#                                       of 823,543), so nearly every window
#                                       miss pays a membership probe against
#                                       that subset before computing.
#
# Three states of this package's source, run interleaved in ONE session on one
# machine (IGeo7, Julia 1.12.6, 8 threads, best of 3 per cell): HEAD = the M2a
# `_CentroidWindow`; SPIKE = the first cell-field version, whose window miss
# re-decoded `cv[k]`; FIELD = the adopted one, whose miss is handed the cell.
# Medians over 5 / 4 / 5 runs, [min, max] beside the small ones, milliseconds:
#
#   row                           HEAD (M2a)          SPIKE               FIELD
#   small (Value,)             20.45 [20.1, 21.0] 18.66 [17.7, 20.7] 17.84 [17.4, 20.6]
#   small (Value, Centroid)    46.91 [45.4, 48.4] 45.31 [43.6, 47.8] 44.33 [42.9, 48.0]
#   small (Value, table)       28.38 [27.3, 28.9] 27.13 [26.0, 28.4] 26.61 [25.8, 28.6]
#   small (Value, field full)                   -              27.11 27.03 [25.9, 28.6]
#   small (Value, field bord)                   -              45.82 45.43 [44.2, 47.5]
#   small (Value, Cent) shuf  111.3 [108, 115]   122.0 [118, 126]   110.7 [108, 118]
#   large (Value,)            147.77             137.79             130.32
#   large (Value, Centroid)   333.66             325.39             317.47
#   large (Value, table)      201.42             194.92             193.14
#   large (Value, field full)      -             195.07             191.62
#   large (Value, field bord)      -             331.70             328.26
#   small (Value, Cent) thr     7.05               6.76               6.54
#   large (Value, Cent) thr    48.33              47.27              46.61
#
# Read the differences, not the rows. The centroid surcharge over the prebuilt
# table -- what the window costs where it hits -- is 19.06 / 18.18 / 17.72 ms
# small and 133.16 / 129.85 / 124.32 ms large, so the field is 7.0% and 6.6%
# under HEAD and the whole `(Value, Centroid)` sweep 5.5% and 4.9% under it:
# nothing was paid for making the centroid an ordinary vector read.
#
# The shuffled row is where a re-decode shows, because the window misses on
# nearly every read there. Its surcharge is 80.59 / 91.33 / 81.33 ms: the spike
# ran 13.3% over HEAD, and handing the miss the cell puts it back (+0.9% of
# HEAD, inside this file's run-to-run spread). The shuffled/storage ratio this
# file prints reads 4.28 [4.08, 4.33] at HEAD, 5.03 [4.73, 5.26] for the spike
# and 4.59 [4.30, 5.22] for the field -- over HEAD's not because the numerator
# grew but because the denominator shrank 7%, which is the paragraph above.
#
# What a partial `known` costs, from rows 2 and 5: the border field computes
# the centroids the window would have, plus one membership probe against the
# border on each miss, for +1.25 ms small (+2.8%) and +14.24 ms large (+4.5%)
# -- 10 to 17 ns a probe over the ~1.03 misses per cell. That is the bar a
# precomputed subset has to clear: it has to save more `cell_centroid` calls
# than its probes cost. A complete `known` (row 4) is row 3 to within noise
# (+0.4 ms small, -1.5 ms large), so the wrapper itself is free.
#
# The gate is a fraction of a same-session pair, never an absolute: the field
# must not be slower than HEAD on any sequential row, and its shuffled
# surcharge must land inside HEAD's run-to-run spread. Both hold above.
# Absolutes from an earlier session do not transfer -- the value-only floor
# alone moved between 17.4 and 20.6 ms inside this one.

using Printf
using Random: shuffle, Xoshiro
import DiscreteGlobalGrids as DGG
import DimensionalData as DD
using DiscreteGlobalGrids: Value, Centroid
const F = DGG.Engine

const RESULTS = Any[]
# Nothing here may be dead code: every measured call's result is folded in.
const SINK = Ref(0.0)

function measure(label, f; reps = 3, alloc = true)
    f()                                   # warm-up: compile, touch every page
    a = alloc ? (@allocated f()) : -1
    best = Inf
    for _ in 1:reps
        GC.gc(false)
        best = min(best, @elapsed f())
    end
    @printf("%-46s %10.3f ms %14s B\n", label, best * 1e3,
        alloc ? string(a) : "-")
    push!(RESULTS, (label, best * 1e9, a))
    return (time = best, alloc = a)
end

# A subtree rooted at one cell, exactly the fixture the `needs` tests use.
rooted(sys, base, depth) =
    DGG.subtree(sys, DGG.cellindex(DGG.levelgrid(sys, base), 3), base + depth)

# A smooth synthetic elevation with a deterministic rough component, laid out
# on the collection's cell axis. No RNG, so the numbers are reproducible.
synthdem(n) = [500.0 * sin(0.00137 * k) + 120.0 * cos(0.0119 * k) +
               3.0 * mod(k * 2654435761, 97) for k in 1:n]

# --- kernels ---------------------------------------------------------------

# Value-only: the steepest *drop*, no geometry. Reads the centre and every
# slot of the one ring it is given.
@inline function dropkernel(center, rings)
    v0 = center[1]
    vs = rings[1]
    best = 0.0
    @inbounds for i in eachindex(vs)
        d = v0 - vs[i]
        d > best && (best = d)
    end
    return best
end

# Steepest descent: drop per unit chord length on the unit sphere. Reads the
# centre value, the centre centroid, and every slot of both rings, so nothing
# it is handed can be elided. Requests 2 and 3 both run THIS function; only
# where `rings[2]` came from differs.
@inline function slopekernel(center, rings)
    v0 = center[1]
    p0 = center[2]
    vs = rings[1]
    ps = rings[2]
    best = 0.0
    @inbounds for i in eachindex(vs)
        q = ps[i]
        dx = p0[1] - q[1]
        dy = p0[2] - q[2]
        dz = p0[3] - q[3]
        s = sqrt(dx * dx + dy * dy + dz * dz)
        g = s > 0.0 ? (v0 - vs[i]) / s : 0.0
        g > best && (best = g)
    end
    return best
end

# One `cell_centroid` per cell and nothing else: the per-call cost the
# `1 + degree` arithmetic multiplies. All three components are summed so no
# part of the transform is dead.
function centroidloop(g, cv)
    s = 0.0
    @inbounds for i in eachindex(cv)
        p = DGG.cell_centroid(g, cv[i])
        s += p[1] + p[2] + p[3]
    end
    return s
end

# --- one size --------------------------------------------------------------

const ROWS = Any[]

function runsize(sys, tag, base, depth; permuted::Bool = false)
    cv = F.CellVector(rooted(sys, base, depth))
    n = length(cv)
    dem = synthdem(n)
    @printf("\n--- %s: IGeo7 subtree(L%d cell 3) to L%d, %d cells, degree %d ---\n",
        tag, base, base + depth, n, DGG.maxneighbors(sys, DGG.Vertex()))

    # The full precomputed table: the ceiling's input, priced on its own.
    build() = [DGG.cell_centroid(cv.grid, c) for c in cv]
    tb = measure("$tag  centroid table build", build)
    table = build()
    SINK[] += table[1][1]

    v_only(thr) = DGG.mapneighbors(dropkernel, cv; needs = (Value(dem),),
        threaded = thr)
    v_cent(thr) = DGG.mapneighbors(slopekernel, cv;
        needs = (Value(dem), Centroid()), threaded = thr)
    v_tab(thr) = DGG.mapneighbors(slopekernel, cv;
        needs = (Value(dem), Value(table)), threaded = thr)

    # The ceiling must answer the same numbers as the thing it bounds, or the
    # comparison is between two different computations.
    @assert v_cent(false) == v_tab(false)
    @assert v_cent(false) == v_cent(true)
    @assert length(v_only(false)) == n
    SINK[] += sum(v_cent(false)) + sum(v_only(false))

    reqs = (("(Value,)", v_only), ("(Value, Centroid)", v_cent),
        ("(Value, Value(table))", v_tab))
    times = Float64[]
    for (label, f) in reqs
        s = measure("$tag  $label  sequential", () -> f(false))
        t = measure("$tag  $label  threaded", () -> f(true); alloc = false)
        push!(times, s.time)
        push!(ROWS, (tag = tag, req = label, n = n, seq = s.time,
            thr = t.time, alloc = s.alloc))
    end

    # Two spellings of the same field, priced against the rows above.
    #
    #   - a COMPLETE field: `known` covers the collection, so the sweep reads
    #     it straight through with no window at all. This row is the
    #     `Value(table)` row with one layer of wrapper on it, and any gap
    #     between them is that wrapper.
    #   - a PARTIAL field: `known` is the border alone, a few percent of the
    #     cells, so nearly every window miss pays a membership probe against
    #     the subset before computing. This row is what the probe costs.
    bd = cv[collect(DGG.border(cv))]
    bcube = DD.DimArray(
        [DGG.cell_centroid(cv.grid, c) for c in bd],
        (DGG.Cells(DGG.CellLookup(bd)),))
    full = DGG.cellfield(DGG.cell_centroid, cv; known = table)
    part = DGG.cellfield(DGG.cell_centroid, cv; known = bcube)
    f_full(thr) = DGG.mapneighbors(slopekernel, cv;
        needs = (Value(dem), Value(full)), threaded = thr)
    f_part(thr) = DGG.mapneighbors(slopekernel, cv;
        needs = (Value(dem), Value(part)), threaded = thr)
    @assert f_full(false) == v_cent(false)
    @assert f_part(false) == v_cent(false)
    @printf("%-46s %10d of %d\n", "$tag  border cells known", length(bd), n)
    for (label, f) in (("(Value, Value(field, full))", f_full),
        ("(Value, Value(field, border))", f_part))
        s = measure("$tag  $label  sequential", () -> f(false))
        t = measure("$tag  $label  threaded", () -> f(true); alloc = false)
        push!(ROWS, (tag = tag, req = label, n = n, seq = s.time,
            thr = t.time, alloc = s.alloc))
    end

    # The window is keyed by the local index, so what it is worth is a
    # property of the visit order. A shuffled `order` has no locality in that
    # index and the window stops hitting; the same-order table row is the
    # control, since reading an array does not care about locality, so the
    # ratio below is the window's loss and not the permuted cursor's.
    permratio = NaN
    if permuted
        perm = shuffle(Xoshiro(42), 1:n)
        p_cent() = DGG.mapneighbors(slopekernel, cv;
            needs = (Value(dem), Centroid()), order = perm, threaded = false)
        p_tab() = DGG.mapneighbors(slopekernel, cv;
            needs = (Value(dem), Value(table)), order = perm, threaded = false)
        # Results are stored in subset index order, so the order may not change
        # a single number.
        @assert p_cent() == v_cent(false)
        @assert p_tab() == v_tab(false)
        pc = measure("$tag  (Value, Centroid)  shuffled order", p_cent)
        pt = measure("$tag  (Value, Value(table))  shuffled order", p_tab)
        permratio = (pc.time - pt.time) / (times[2] - times[3])
        @printf("%-46s %10.2fx\n",
            "$tag  centroid surcharge shuffled / storage", permratio)
    end

    # One `cell_centroid`, in isolation.
    cl = measure("$tag  cell_centroid loop ($n calls)",
        () -> centroidloop(cv.grid, cv))
    SINK[] += centroidloop(cv.grid, cv)
    @printf("%-46s %10.1f ns/call\n", "$tag  cell_centroid per call",
        cl.time / n * 1e9)

    return (tag = tag, n = n, t1 = times[1], t2 = times[2], t3 = times[3],
        percall = cl.time / n * 1e9, build = tb.time, buildbytes = tb.alloc,
        permratio = permratio)
end

# --- run -------------------------------------------------------------------

println("julia ", VERSION, "   threads = ", Threads.nthreads())

const SYS = DGG.IGeo7System()
const DEGREE = DGG.maxneighbors(SYS, DGG.Vertex())

# 7^6 = 117,649 cells and 7^7 = 823,543; see the aperture note at the top.
const SUMMARY = [runsize(SYS, "small", 1, 6; permuted = true),
    runsize(SYS, "large", 1, 7)]

# --- table -----------------------------------------------------------------

println()
println("=== baseline: needs requests on IGeo7 rooted subtrees ===")
@printf("%-24s %10s %10s %10s %14s %14s %14s\n", "request", "cells", "seq s",
    "thr s", "alloc B", "cent/cell 2-1", "gain 2-3 s")
for tag in ("small", "large")
    s = only(filter(x -> x.tag == tag, SUMMARY))
    for r in filter(x -> x.tag == tag, ROWS)
        # Both derived columns are properties of request 2 against its two
        # bracketing rows, so they are stated on that row and nowhere else.
        wanted = r.req == "(Value, Centroid)"
        cent = wanted ? @sprintf("%11.1f ns", (s.t2 - s.t1) / s.n * 1e9) :
               "             -"
        gain = wanted ? @sprintf("%14.5f", s.t2 - s.t3) : "             -"
        @printf("%-24s %10d %10.5f %10.5f %14d %14s %14s\n", r.req, r.n,
            r.seq, r.thr, r.alloc, cent, gain)
    end
end

println()
println("=== reading the table ===")
for s in SUMMARY
    @printf("--- %s (%d cells) ---\n", s.tag, s.n)
    @printf("  %-52s %10.5f s\n", "centroid table build", s.build)
    @printf("  %-52s %10d B\n", "centroid table bytes", s.buildbytes)
    @printf("  %-52s %10.1f ns\n", "one cell_centroid, isolated", s.percall)
    @printf("  %-52s %10.1f ns\n",
        "  x (1 + degree) = $(1 + DEGREE) per cell", (1 + DEGREE) * s.percall)
    @printf("  %-52s %10.1f ns\n", "observed centroid cost per cell (2 - 1)",
        (s.t2 - s.t1) / s.n * 1e9)
    @printf("  %-52s %9.1f %%\n", "centroid share of the Centroid() sweep",
        100 * (s.t2 - s.t1) / s.t2)
    @printf("  %-52s %10.1f ns\n",
        "one in-sweep centroid touch ((2-3)/($(1 + DEGREE)n))",
        (s.t2 - s.t3) / ((1 + DEGREE) * s.n) * 1e9)
    @printf("  %-52s %9.1f %%\n", "removed by the prebuilt table ((2-3)/(2-1))",
        100 * (s.t2 - s.t3) / (s.t2 - s.t1))
    # A per-sweep working set still computes one centroid per cell, so the
    # table's own build is part of what it would have to pay: this is the
    # figure the `1 + degree -> 1` model actually predicts.
    @printf("  %-52s %9.1f %%\n",
        "  with the build paid too ((2-3-build)/(2-1))",
        100 * (s.t2 - s.t3 - s.build) / (s.t2 - s.t1))
    @printf("  %-52s %9.1f %%\n",
        "predicted by $(1 + DEGREE) -> 1 calls ($DEGREE/$(1 + DEGREE))",
        100 * DEGREE / (1 + DEGREE))
    @printf("  %-52s %10.2fx\n", "speedup of the ceiling over today (2 / 3)",
        s.t2 / s.t3)
end

println()
println("sink = ", SINK[])
println("RESULTS_BEGIN")
for (l, t, a) in RESULTS
    @printf("%s|%.1f|%d\n", l, t, a)
end
println("RESULTS_END")
