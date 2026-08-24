# The M2 baseline: what `Centroid()` costs a sweep today, and what a perfect
# cache could give back.
#
#     julia -t auto --project=benchmark benchmark/needs_centroid.jl
#
# `needs = (Value(dem), Centroid())` recomputes `cell_centroid` once per touch:
# once for the visited cell and once for each of its clipped neighbours, so
# `1 + degree` computations per cell where a per-sweep working set would want
# one. Three requests bracket that on the same cells, with the same kernel
# wherever the kernel can be the same:
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

using Printf
import DiscreteGlobalGrids as DGG
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

function runsize(sys, tag, base, depth)
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

    # One `cell_centroid`, in isolation.
    cl = measure("$tag  cell_centroid loop ($n calls)",
        () -> centroidloop(cv.grid, cv))
    SINK[] += centroidloop(cv.grid, cv)
    @printf("%-46s %10.1f ns/call\n", "$tag  cell_centroid per call",
        cl.time / n * 1e9)

    return (tag = tag, n = n, t1 = times[1], t2 = times[2], t3 = times[3],
        percall = cl.time / n * 1e9, build = tb.time, buildbytes = tb.alloc)
end

# --- run -------------------------------------------------------------------

println("julia ", VERSION, "   threads = ", Threads.nthreads())

const SYS = DGG.IGeo7System()
const DEGREE = DGG.maxneighbors(SYS, DGG.Vertex())

# 7^6 = 117,649 cells and 7^7 = 823,543; see the aperture note at the top.
const SUMMARY = [runsize(SYS, "small", 1, 6), runsize(SYS, "large", 1, 7)]

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
