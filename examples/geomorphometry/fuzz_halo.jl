#!/usr/bin/env julia
# Randomised stress test of the subtree/subset halo API.
#
#   julia --project=test examples/geomorphometry/fuzz_halo.jl [ncases] [seed]
#
# Randomizes system, root level, root cell, depth, connectivity, and elevation
# field across three case families:
#
#   A. chunked terrain    — compare whole-grid and chunk-plus-halo metrics and
#                           count unavailable reads. Full grids are capped at
#                           `MAXCELLS`.
#   B. deep subtree       — compare the halo with a one-ring oracle built only
#                           from subtree cells, allowing larger target levels.
#   C. random subset      — halo(PartialGrid|CellVector|CellLookup) over a
#                           random member set, holes and all.
#
# The script also reports heap growth and repeated-call allocation stability.

include("Harness.jl")
using .Harness
using .Harness: SphericalTerrain, whole_field
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: Vertex, Edge
using Random, Printf

const MAXCELLS = 60_000
const FIELDS = (:noise, :harmonic, :step, :ramp, :spike)
const CONNS = (Vertex(), Edge())

"""
    deep_subtree_failures(sys, root, level, conn)

Halo vs a brute-force halo built only from the subtree's own cells. `O(subtree)`
rather than `O(grid)`, so this reaches levels where a whole level grid is out
of the question — which is where the halo API earns its keep.
"""
function deep_subtree_failures(sys, root, level, conn)
    fails = Failure[]
    tag = describe(sys, root, level, conn)
    g = DGG.levelgrid(sys, level)
    r = chunk_range(sys, root, level)
    lo, hi = first(r), last(r)
    seen = Set{Int}()
    for p in r
        for m in DGG.neighbors(g, DGG.cellindex(g, p), 1; connectivity = conn)
            q = DGG.cellposition(g, m)
            (lo <= q <= hi) || push!(seen, q)
        end
    end
    want = sort!(collect(seen))
    got = collect(DGG.halo(DGG.subtree(sys, root, level); connectivity = conn))
    issorted(got) || push!(fails, Failure(:unsorted, tag, "halo not ascending"))
    length(unique(got)) == length(got) ||
        push!(fails, Failure(:duplicate, tag, "halo has duplicates"))
    got == want || push!(fails, Failure(:halo_mismatch, tag,
        "deep halo != brute force: missing=$(setdiff(want, got)) extra=$(setdiff(got, want))"))
    return fails, length(want), length(r)
end

function random_subset(rng, n, style)
    if style === :blob
        lo = rand(rng, 1:n)
        hi = min(n, lo + rand(rng, 1:max(1, n ÷ 4)))
        v = collect(lo:hi)
    elseif style === :scatter
        v = sort!(unique!(rand(rng, 1:n, max(2, n ÷ 8))))
    elseif style === :holed
        lo = rand(rng, 1:n)
        hi = min(n, lo + rand(rng, 4:max(5, n ÷ 3)))
        v = collect(lo:hi)
        for _ in 1:max(1, length(v) ÷ 6)
            length(v) > 2 && deleteat!(v, rand(rng, 2:length(v) - 1))
        end
    else # :bands
        v = Int[]
        for p in 1:n
            (p % 7) < 3 && push!(v, p)
        end
    end
    return isempty(v) ? [1] : v
end

function main(ncases::Int, seed::UInt)
    rng = MersenneTwister(seed)
    syslist = collect(DGG.systems())
    fails = Failure[]
    stats = Dict(:A => 0, :B => 0, :C => 0)
    biggest = (0, "")
    halocells = 0
    heap0 = Base.gc_live_bytes()
    t0 = time()

    @printf("fuzz: %d cases, seed = 0x%016x\n", ncases, seed)
    println("-"^78)

    for i in 1:ncases
        family = rand(rng, (:A, :A, :B, :C))
        sys = rand(rng, syslist)
        conn = rand(rng, CONNS)
        ml = DGG.maxlevel(sys)

        if family === :A
            rootlevel = rand(rng, 0:3)
            depth = rand(rng, 0:4)
            target = min(rootlevel + depth, ml)
            target < rootlevel && continue
            (sys isa DGG.A5System && target > 6) && continue
            DGG.ncells(sys, target) > MAXCELLS && continue
            DGG.ncells(sys, rootlevel) > 40_000 && continue
            groot = DGG.levelgrid(sys, rootlevel)
            root = DGG.cellindex(groot, rand(rng, 1:DGG.ncells(groot)))
            f, s = check_subtree_case(sys, root, target, conn, rand(rng, FIELDS))
            append!(fails, f); stats[:A] += 1; halocells += s.nhalo
            s.nhalo > biggest[1] && (biggest = (s.nhalo, describe(sys, root, target, conn)))
        elseif family === :B
            rootlevel = rand(rng, 0:4)
            depth = rand(rng, 3:6)
            target = min(rootlevel + depth, ml)
            target <= rootlevel && continue
            # A5 uses the generic geometry walk, so keep randomized cases shallow.
            (sys isa DGG.A5System && target > 6) && continue
            # Limit the number of materialized subtree cells.
            sub = DGG.ncells(sys, target) ÷ max(DGG.ncells(sys, rootlevel), 1)
            sub > 120_000 && continue
            groot = DGG.levelgrid(sys, rootlevel)
            DGG.ncells(sys, rootlevel) > 40_000 && continue
            root = DGG.cellindex(groot, rand(rng, 1:DGG.ncells(groot)))
            f, nh, nc = deep_subtree_failures(sys, root, target, conn)
            append!(fails, f); stats[:B] += 1; halocells += nh
            nh > biggest[1] && (biggest = (nh, describe(sys, root, target, conn)))
        else
            level = rand(rng, 1:4)
            DGG.ncells(sys, level) > 12_000 && continue
            k = gridctx(sys, level, conn)
            members = random_subset(rng, length(k), rand(rng, (:blob, :scatter, :holed, :bands)))
            append!(fails, check_subset_case(sys, level, conn, members;
                label = "fuzz case $i"))
            stats[:C] += 1
        end

        if i % 200 == 0
            @printf("  %5d/%d  A=%d B=%d C=%d  fails=%d  heap=%+.1f MB  %.0fs\n",
                i, ncases, stats[:A], stats[:B], stats[:C], length(fails),
                (Base.gc_live_bytes() - heap0) / 2^20, time() - t0)
        end
    end

    println("-"^78)
    println("### heap under a fixed configuration (no new harness cache keys)")
    # Cache keys are fixed here, so repeated calls do not add grid or metric entries.
    fixsys = DGG.HEALPixSystem()
    gridctx(fixsys, 5, Vertex()); whole_field(fixsys, 5, Vertex(), :noise, 1)
    groot = DGG.levelgrid(fixsys, 2)
    for p in 1:DGG.ncells(groot)
        check_subtree_case(fixsys, DGG.cellindex(groot, p), 5, Vertex(), :noise;
            spike_in_halo = false, control = false)
    end
    GC.gc(); GC.gc()
    h1 = Base.gc_live_bytes()
    reps = 0
    for _ in 1:6, p in 1:DGG.ncells(groot)
        f, _ = check_subtree_case(fixsys, DGG.cellindex(groot, p), 5, Vertex(), :noise;
            spike_in_halo = false, control = false)
        append!(fails, f); reps += 1
    end
    GC.gc(); GC.gc()
    @printf("  %d repeats of %d warmed configurations: heap %+.2f MB\n",
        reps, DGG.ncells(groot), (Base.gc_live_bytes() - h1) / 2^20)

    println("-"^78)
    println("### allocation stability (same case, 200 repeats)")
    sys = DGG.HEALPixSystem()
    root = DGG.cellindex(DGG.levelgrid(sys, 2), 7)
    collect(DGG.halo(DGG.subtree(sys, root, 8)))      # warm
    a1 = @allocated collect(DGG.halo(DGG.subtree(sys, root, 8)))
    for _ in 1:200; collect(DGG.halo(DGG.subtree(sys, root, 8))); end
    a2 = @allocated collect(DGG.halo(DGG.subtree(sys, root, 8)))
    @printf("  halo(subtree) allocations: first %d B, after 200 repeats %d B\n", a1, a2)
    a1 == a2 || push!(fails, Failure(:alloc_drift, "HEALPix L2 root -> L8",
        "per-call allocation drifted from $a1 B to $a2 B over 200 identical calls"))

    GC.gc(); GC.gc()
    @printf("  heap delta over whole run: %+.1f MB\n",
        (Base.gc_live_bytes() - heap0) / 2^20)

    println("-"^78)
    @printf("cases run: A(chunked terrain)=%d  B(deep subtree)=%d  C(random subset)=%d\n",
        stats[:A], stats[:B], stats[:C])
    @printf("halo cells generated: %d;  largest single halo: %d  (%s)\n",
        halocells, biggest[1], biggest[2])
    @printf("elapsed %.0fs\n", time() - t0)
    println("FAILURES: ", length(fails))
    for f in fails[1:min(end, 40)]; println(f); end
    return length(fails)
end

const NCASES = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2000
const SEED = length(ARGS) >= 2 ? parse(UInt, ARGS[2]) : 0x5EED_0FA1_0BAD_0001 % UInt
exit(main(NCASES, SEED) == 0 ? 0 : 1)
