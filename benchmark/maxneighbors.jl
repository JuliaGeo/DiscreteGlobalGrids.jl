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
    f1() = F._positioned(cv, 1, 1, ring, Val(M))
    f2() = F._positioned(cv, 1, 1, ring, nothing)
    f1(); f2()
    @printf("%-52s %14d B\n", "$sysname  _positioned SmallVector{$M}", @allocated f1())
    @printf("%-52s %14d B\n", "$sysname  _positioned Vector (fallback)", @allocated f2())
end

println()
println("RESULTS_BEGIN")
for (l, t, a) in RESULTS
    @printf("%s|%.1f|%d\n", l, t, a)
end
println("RESULTS_END")
