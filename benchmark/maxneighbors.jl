# Measured 2026-08-17 on M-series macOS, Julia 1.12.6, 4 threads (min of 9):
# declared systems are byte-identical in serial allocations and within ±4% in
# time vs the SmallVector build. The fallback path itself (same rings, Vector
# buffer forced): HEALPix _sweep! +7.9% / +10.2 MB, mapneighbors threaded
# +9.9%; IGeo7 _sweep! +4.4% / +6.4 MB, mapneighbors threaded +5.9%; per-cell
# transient 256 B (HEALPix, M=8) / 160 B (IGeo7, M=6). undef+resize! is the
# load-bearing detail: push!-from-empty costs 4-6x the allocation and 6x the
# slowdown.
#
# The cost of the fallback itself.
#
# Part 1 isolates the buffer choice on a REAL system: the same HEALPix/IGeo7
# CellVector, the same automaton one-rings, driven through `_sweep!` with an
# explicit capacity of `Val(M)` (declared) or `nothing` (undeclared). Nothing
# but the ring container differs.
#
# Part 2 is the end-to-end public API on the audit's toy, in a declared and an
# undeclared copy. That toy has no fast `neighbors`, so it also shows what the
# buffer costs RELATIVE to the geometric fallback an undeclared system is
# likely to be using anyway.

using Printf
import InteractiveUtils
import DiscreteGlobalGrids as DGG
const F = DGG.Engine
include(joinpath(@__DIR__, "toys.jl"))

subset(sys, lvl, n) = DGG.PartialGrid(sys, lvl,
    [DGG.cellindex(DGG.levelgrid(sys, lvl), i) for i in 1:n])

const RESULTS = Any[]

function measure(label, f; reps = 9)
    f()
    a = @allocated f()
    best = Inf
    for _ in 1:reps
        GC.gc(false)
        best = min(best, @elapsed f())
    end
    @printf("%-52s %10.3f ms %14d B\n", label, best * 1e3, a)
    push!(RESULTS, (label, best * 1e9, a))
    return nothing
end

println("julia ", VERSION, "   threads = ", Threads.nthreads())

# --- Part 1: same system, capacity forced ----------------------------------

# `_sweep!` with a counting closure: the ring is built, read, and dropped.
function sweepsum(cv, conn, cap)
    s = Ref(0)
    F._sweep!(cv, conn, 1:length(cv), cap) do k, c, nbrs
        s[] += length(nbrs)
    end
    return s[]
end

# `_adjacency_chunk` is the CSR row builder; it takes the capacity explicitly.
function chunkbuild(cv, conn, cap)
    offsets = Vector{Int}(undef, length(cv) + 1)
    offsets[1] = 1
    return F._adjacency_chunk(cv, conn, 1:length(cv), F.ClippedRows(), offsets, cap)
end

println("\n=== Part 1: declared vs forced-Vector, same system ===")
for (sysname, sys, lvl, n) in (("HEALPix L6", DGG.HEALPixSystem(), 6, 40_000),
                               ("IGeo7 L5", DGG.IGeo7System(), 5, 40_000))
    pg = subset(sys, lvl, n)
    cv = F.CellVector(pg)
    conn = DGG.Vertex()
    M = DGG.maxneighbors(sys, conn)
    @assert sweepsum(cv, conn, Val(M)) == sweepsum(cv, conn, nothing)
    @assert chunkbuild(cv, conn, Val(M)) == chunkbuild(cv, conn, nothing)
    println("--- ", sysname, "  (", length(cv), " cells, M=", M, ") ---")
    measure("$sysname  _sweep! SmallVector{$M}", () -> sweepsum(cv, conn, Val(M)))
    measure("$sysname  _sweep! Vector (fallback)", () -> sweepsum(cv, conn, nothing))
    measure("$sysname  _adjacency_chunk SmallVector{$M}", () -> chunkbuild(cv, conn, Val(M)))
    measure("$sysname  _adjacency_chunk Vector (fallback)", () -> chunkbuild(cv, conn, nothing))
    # Threaded: allocation in parallel tasks is where GC pressure bites.
    thr = F.GOCore.booltype(true)
    mapdeg(cap) = F._mapstore!((c, nb) -> length(nb), Vector{Int}(undef, length(cv)),
        cv, conn, F.StorageOrder(), thr, cap)
    @assert mapdeg(Val(M)) == mapdeg(nothing)
    measure("$sysname  mapneighbors thr SmallVector{$M}", () -> mapdeg(Val(M)))
    measure("$sysname  mapneighbors thr Vector (fallback)", () -> mapdeg(nothing))
    println()
end

# --- Part 2: end-to-end on the toy -----------------------------------------

println("=== Part 2: toy system, declared vs undeclared, public API ===")
lvl = 4
n = 1500
for (name, sys) in (("Octant-declared", OctantDeclared()),
                    ("Octant-bare    ", OctantBare()))
    pg = subset(sys, lvl, n)
    cv = F.CellVector(pg)
    data = collect(1.0:length(cv))
    println("--- ", name, "  (", length(cv), " cells, maxneighbors=",
        repr(DGG.maxneighbors(sys, DGG.Vertex())), ") ---")
    measure("$name  neighbors(cv) sweep",
        () -> (s = 0; for (_, nb) in DGG.neighbors(cv); s += length(nb); end; s))
    measure("$name  mapneighbors deg (serial)",
        () -> DGG.mapneighbors((c, nb) -> length(nb), cv; threaded = false))
    measure("$name  mapneighbors data-mean (serial)",
        () -> DGG.mapneighbors((c, v, vs) -> (v + sum(vs; init = 0.0)) /
                                             (1 + length(vs)), cv, data;
            threaded = false))
    measure("$name  adjacency (serial)", () -> DGG.adjacency(cv; threaded = false))
    measure("$name  adjacency (threaded)", () -> DGG.adjacency(cv))
    println()
end

# --- Part 3: per-cell allocation of a single clipped ring ------------------

println("=== Part 3: one clipped ring, @allocated ===")
for (sysname, sys, lvl, n) in (("HEALPix L6", DGG.HEALPixSystem(), 6, 40_000),
                               ("IGeo7 L5", DGG.IGeo7System(), 5, 40_000))
    pg = subset(sys, lvl, n)
    cv = F.CellVector(pg)
    conn = DGG.Vertex()
    M = DGG.maxneighbors(sys, conn)
    c = cv[1000]
    ring = DGG.neighbors(cv.grid, c, 1; connectivity = conn)
    f1() = F._indexed(cv, 1, 1, ring, Val(M))
    f2() = F._indexed(cv, 1, 1, ring, nothing)
    f1(); f2()
    @printf("%-52s %14d B\n", "$sysname  _indexed SmallVector{$M}", @allocated f1())
    @printf("%-52s %14d B\n", "$sysname  _indexed Vector (fallback)", @allocated f2())
end

# --- Part 4: where a stack container stops paying (STATIC_RING_CAP) --------
#
# `Engine.STATIC_RING_CAP` (64 elements) and `STATIC_RING_BYTES` (512) are
# COMPILE-time limits, so what is measured is emitted code size per element,
# not run time. `SmallVector{N,T}` is a distinct type per `N`; each `N`
# respecialises the neighbourhood stack, and past the cliff each one costs more
# to compile and emits more code for the same work.
#
# Measured 2026-08-23, M-series macOS, Julia 1.12, SmallCollections 0.6.3:
# 8-byte ids step at N=65 (315 -> 492 B/elem, +56%) with no recovery above
# (752 at 96) -- that cliff is what STATIC_RING_CAP catches. 4-byte ids show no
# cliff through 160, so the element cap is merely conservative there. 16-byte
# ids show no cliff either but cost ~3.6x per element everywhere, so
# STATIC_RING_BYTES bounds the total rather than catching a jump.

const SC = DGG.SmallCollections

println("=== Part 4: SmallVector code size per element ===")

struct Pair16Bench
    a::UInt64
    b::UInt64
end
Base.isless(x::Pair16Bench, y::Pair16Bench) = x.a < y.a

# A clip-shaped loop: the shape `_clip` and `_indexed` actually compile.
for (T, tag, ns) in ((UInt64, "8-byte id", (32, 48, 56, 62, 64, 65, 66, 72, 80, 96)),
                     (UInt32, "4-byte id", (64, 96, 128, 144, 160)),
                     (Pair16Bench, "16-byte id", (16, 24, 32, 40, 48, 64)))
    println("--- $tag (sizeof = $(sizeof(T))) ---")
    @printf("%6s %10s %12s %10s  %s\n", "N", "payloadB", "nativeB", "B/elem", "capacity")
    for N in ns
        fn = Symbol("_ladder_", tag[1], "_", N)
        @eval function $fn(src::Vector{$T}, lo::$T)
            v = SC.SmallVector{$N,$T}()
            @inbounds for i in eachindex(src)
                src[i] > lo && (v = SC.push(v, src[i]))
            end
            return v
        end
        f = getfield(@__MODULE__, fn)
        na = length(sprint(io -> InteractiveUtils.code_native(io, f, (Vector{T}, T))))
        cap = F._capacity(N, T)
        @printf("%6d %10d %12d %10.1f  %s\n", N, N * sizeof(T), na, na / N,
            cap === nothing ? "heap" : "Val{$N}")
    end
    println()
end

# --- Part 5: the carried walk vs the azimuth sort, same system -------------

# Both paths get the SAME grid, the same cells and the same one-ring automaton,
# so the only difference is how a shell is ordered: `_shells_azimuth` sorts every
# candidate by bearing around the centre (a centroid per candidate), while the
# carried walk inherits the order from the neighbour that produced it and only
# pins the shell's starting phase. IGeo7 declares `CounterClockwise` and a
# `maxring` of `6k`, so the public API here is the carried walk with stack
# buffers; Part 5 forces the azimuth path through the internal entry to show
# what that declaration bought.
#
# Read the `k = 1` rows differently from the rest. `ring`/`neighbors` have
# always answered `k <= 1` straight from the automaton, so that row is not two
# versions of the public API but the short-circuit against the shell machinery
# it skips — it prices the short-circuit, and the `0 B` is the property the
# `@allocated` tests in each system's suite pin down. The `k >= 2` rows are the
# real before-and-after: same rings, same automaton, different ordering.
#
# Measured 2026-08-23, IGeo7 L8, 4 cells per call, min of 9 (machine in low
# power mode, so read the ratios and the byte counts, not the absolute times):
#
#     ring k=1  3.95x   0 B vs 2880 B      neighbors k=1  5.30x   0 B vs 3008 B
#     ring k=2  1.85x  1728 B vs 9536 B    neighbors k=2  1.86x  1920 B vs 9920 B
#     ring k=3  1.95x  2304 B vs 14400 B   neighbors k=3  1.97x  2944 B vs 15040 B
#     ring k=4  1.89x  2880 B vs 30144 B   neighbors k=4  1.81x  4160 B vs 31040 B
println("\n=== Part 5: carried walk vs azimuth sort (IGeo7 L8) ===")
let sys = DGG.IGeo7System(), g = DGG.levelgrid(sys, 8), conn = DGG.Vertex(),
        FB = DGG.Fallbacks
    cells = [DGG.cellindex(g, i) for i in (3_000_000, 17, 250_000, 41_237_881)]
    # Consume the shell into a scalar: returning a container to an untyped
    # caller boxes it, which would price the measurement rather than the walk.
    drain(v) = (s = UInt(0); for x in v; s ⊻= hash(x); end; s)
    wound_ring(k)  = (s = UInt(0); for c in cells; s ⊻= drain(DGG.ring(g, c, k)); end; s)
    wound_disc(k)  = (s = UInt(0); for c in cells; s ⊻= drain(DGG.neighbors(g, c, k)); end; s)
    function azim(k, discp)
        s = UInt(0)
        for c in cells
            sh = FB._shells_azimuth(x -> DGG.one_ring(g, x, conn), g, c, k)
            if discp
                for r in sh; s ⊻= drain(r); end
            else
                k <= length(sh) && (s ⊻= drain(@inbounds sh[k]))
            end
        end
        return s
    end
    @printf("%-52s %13s %14s\n", "", "time", "alloc")
    for k in 1:4
        measure(@sprintf("ring k=%d   carried walk", k), () -> wound_ring(k))
        measure(@sprintf("ring k=%d   azimuth sort", k), () -> azim(k, false))
        measure(@sprintf("neighbors k=%d   carried walk", k), () -> wound_disc(k))
        measure(@sprintf("neighbors k=%d   azimuth sort", k), () -> azim(k, true))
    end
end

# --- Part 6: `Val(k)` against `k` ---------------------------------------

# The two forms answer the same cells in the same order; all that differs is
# whether the capacity a system declares is a compile-time constant. It only is
# when `k` is in the type, which is why every call below writes `Val(2)` rather
# than passing a variable — `Val(k)` built from a run-time `k` is a dynamic
# dispatch and gives back everything the form was for.
#
# `k <= 1` is answered from the automaton either way, so it is here only to show
# that the `Val` form does not cost anything to reach it.
#
# Measured 2026-08-23, IGeo7 L8, 4 cells per call, min of 9. Read the bytes, not
# the times: at these sizes the timer's resolution is coarser than the
# difference, and the form removes allocations rather than work.
#
#     ring k=1       0 B ->    0 B     neighbors k=2   480 B -> 112 B
#     ring k=2     432 B ->  224 B     neighbors k=4  1040 B -> 208 B
#     ring k=3     576 B ->  320 B
#     ring k=4     720 B ->  416 B
#
# and `ring(g, c, Val(2))` infers `SmallVector{12,Z7Cell}` where
# `ring(g, c, 2)` infers `Vector{Z7Cell}` — same cells, different container.
println("\n=== Part 6: Val(k) vs k (IGeo7 L8) ===")
let g = DGG.levelgrid(DGG.IGeo7System(), 8), c = DGG.cellindex(g, 3_000_000)
    drain(v) = (s = UInt(0); for x in v; s ⊻= hash(x); end; s)
    iring1(g, c) = drain(DGG.ring(g, c, 1));        vring1(g, c) = drain(DGG.ring(g, c, Val(1)))
    iring2(g, c) = drain(DGG.ring(g, c, 2));        vring2(g, c) = drain(DGG.ring(g, c, Val(2)))
    iring3(g, c) = drain(DGG.ring(g, c, 3));        vring3(g, c) = drain(DGG.ring(g, c, Val(3)))
    iring4(g, c) = drain(DGG.ring(g, c, 4));        vring4(g, c) = drain(DGG.ring(g, c, Val(4)))
    idisc2(g, c) = drain(DGG.neighbors(g, c, 2));   vdisc2(g, c) = drain(DGG.neighbors(g, c, Val(2)))
    idisc4(g, c) = drain(DGG.neighbors(g, c, 4));   vdisc4(g, c) = drain(DGG.neighbors(g, c, Val(4)))
    for (label, fi, fv) in (("ring k=1", iring1, vring1), ("ring k=2", iring2, vring2),
            ("ring k=3", iring3, vring3), ("ring k=4", iring4, vring4),
            ("neighbors k=2", idisc2, vdisc2), ("neighbors k=4", idisc4, vdisc4))
        measure("$label   k::Integer", () -> fi(g, c))
        measure("$label   Val(k)", () -> fv(g, c))
    end
    # The point of the form: a concrete, stack-sized return type. Infer the
    # neighbourhood call itself — the timed helpers above drain into a `UInt64`,
    # so inferring those would report the accumulator and not the container.
    tring2(g, c) = DGG.ring(g, c, Val(2))
    tdisc2(g, c) = DGG.neighbors(g, c, Val(2))
    iring2t(g, c) = DGG.ring(g, c, 2)
    for (label, f) in (("ring(g, c, Val(2))", tring2),
            ("neighbors(g, c, Val(2))", tdisc2), ("ring(g, c, 2)", iring2t))
        t = only(Base.return_types(f, (typeof(g), typeof(c))))
        @printf("  %-30s %-8s %s\n", "infers $label",
            isconcretetype(t) ? "concrete" : "ABSTRACT", t)
    end
end

println()
println("RESULTS_BEGIN")
for (l, t, a) in RESULTS
    @printf("%s|%.1f|%d\n", l, t, a)
end
println("RESULTS_END")
