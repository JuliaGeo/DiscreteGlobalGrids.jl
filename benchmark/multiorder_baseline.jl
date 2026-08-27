# Multi-order coverage (MOC) container baseline.
#
#     julia --project=benchmark --threads=1 benchmark/multiorder_baseline.jl
#
# Seven arms, all on synthetic HEALPix fixtures, no downloads:
#
#   1. INFERENCE   — the return type of every public entry point. The package
#                    carries a juliac ambition, so an answer that is neither
#                    concrete nor the small union the docstring promises is a
#                    defect rather than a tuning knob.
#   2. LOOKUP      — `localindex` / `covering_index` against container size.
#                    The interval index is a binary search, so cost grows like
#                    log n, not n, and nothing allocates.
#   3. SET ALGEBRA — `union` / `intersect` / `setdiff` / `symdiff` and the two
#                    predicates against size: sorted-interval merges, linear in
#                    the operands.
#   4. HIERARCHY   — `aggregate`, `coarsen` and container construction, in leaf
#                    cells per second on a level-8 field.
#   5. EXPAND      — the lazy presentation: per-element access is one binary
#                    search and allocates nothing, `collect` is a run fill.
#   6. STORE IO    — a compacted write+read against the expanded write of the
#                    same field, same session, same data.
#   7. REGRID      — `BarycentricPoint` off the stored cells against off their
#                    reference-level expansion, same session, same field: build
#                    and apply, and what the reference level costs.
#
# Plus a region note: `covering_indices` and the `Covering` selector are bounded
# not by the container but by `query(sys, MultiOrderCoverage(...); level)`,
# which is measured on its own so the two costs are not confused.
#
# The IO arm needs Zarr, which is a weak dependency: it self-skips when absent.
#
# Checkpoint 2026-08-27, Julia 1.12.7, 1 thread, M-series laptop, 140 s and
# 1 563 MiB peak of which most is loading and compiling this script:
#
#   * every entry point concrete, `Union{T,Nothing}`, or the documented
#     two-way shape union;
#   * `localindex` / `covering_index` 8.6 -> 42 ns and ZERO bytes across
#     1e3 -> 1e6 cells, so 1000x the container for 5x the lookup;
#   * set algebra 9-20 ns per output cell, flat in n;
#   * `coarsen` on 786 432 leaves: 0.91 ms at atol 0.05 (10 896 stored cells),
#     12.0 ms at atol 0.001 (778 992); `aggregate` 1.2 ms one level up;
#   * `expand` access 9.5 ns and zero bytes; `collect` 0.07 ns/cell;
#   * a compacted store is 0.095x the expanded one on disk and 2.7x slower to
#     read back, because reading validates every (level, id) pair;
#   * `BarycentricPoint` off 624 stored cells onto 12 288: 819 ns and 0 bytes
#     per destination point, 3.9 MiB of plan. Re-keying those SAME cells two
#     reference levels deeper (16x the leaves) leaves the native plan at 1.02x
#     time and 1.00x bytes; the expansion's goes to 3.52x bytes, so at that
#     reference level the native plan is 0.22x the expanded one, and applying
#     it 0.07x. Build TIME is the destination sweep either way — the source
#     presentation moves the bytes, not the clock.
#
# The one cost that is NOT the container's is the region arm: see its note.

import BenchmarkTools
import DimensionalData as DD
import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import GlobalRegridding as GR
import Printf: @printf, @sprintf

# Loading Zarr is what makes the store extension exist; arm 6 checks and skips.
try
    @eval import Zarr
catch
end

const EN = DGG.Engine
const CL = DGG.CellLookups
const LK = DD.Lookups

const SYS = DGG.HEALPixSystem()

const REGION = GI.Polygon([GI.LinearRing([(-20.0, 20.0), (40.0, 20.0),
    (40.0, 60.0), (-20.0, 60.0), (-20.0, 20.0)])])

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# `n` disjoint level-`l` cells at every `step`-th index, so no two are adjacent
# and the container really holds `n` intervals rather than one merged run.
function spread_container(n::Int; level::Int = 9, first::Int = 1, step::Int = 2)
    g = DGG.levelgrid(SYS, level)
    first + step * (n - 1) <= DGG.ncells(g) ||
        error("level $level is too shallow for $n cells at step $step")
    cells = [DGG.cellindex(g, first + step * (k - 1)) for k in 1:n]
    return EN.MultiOrderVector(SYS, cells; reference_level = level)
end

# Deterministic query cells at the container's own level, half of them stored
# and half not, so the hit and the miss path are both exercised.
query_cells(n::Int; level::Int = 9, q::Int = 10_000) =
    (g = DGG.levelgrid(SYS, level);
     [DGG.cellindex(g, 1 + (7919 * (k - 1)) % (2n)) for k in 1:q])

# A smooth field over a whole level — the input `coarsen` is designed for,
# where large regions agree to within `atol` and merge.
function level_field(l::Int)
    g = DGG.levelgrid(SYS, l)
    n = DGG.ncells(g)
    vals = Vector{Float64}(undef, n)
    @inbounds for k in 1:n
        vals[k] = DGG.cell_centroid(g, DGG.cellindex(g, k))[3]
    end
    return DGG.CellVector(g), vals
end

# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------

# `gcsample` collects before every sample. It costs wall time and buys the
# only thing that matters on a laptop: the heap never carries an arm's garbage
# into the next one.
function estimate(f; seconds = 1.0, samples = 12, gc = true)
    trial = BenchmarkTools.run(BenchmarkTools.@benchmarkable($f());
        samples, evals = 1, seconds, gcsample = gc)
    e = BenchmarkTools.minimum(trial)
    return (; ns = e.time, bytes = e.memory, allocs = e.allocs)
end

rule(title) = println("\n", title, "\n", "="^length(title))

row(cols, widths) =
    println(join((rpad(string(c), w) for (c, w) in zip(cols, widths)), "  "))

fmt_ns(x) = x < 1e3 ? @sprintf("%.1f ns", x) :
            x < 1e6 ? @sprintf("%.2f us", x / 1e3) :
            x < 1e9 ? @sprintf("%.2f ms", x / 1e6) : @sprintf("%.2f s", x / 1e9)

fmt_bytes(b) = b < 1024 ? "$b B" :
               b < 2^20 ? @sprintf("%.1f KiB", b / 1024) :
               @sprintf("%.1f MiB", b / 2^20)

# ---------------------------------------------------------------------------
# 1. Inference
# ---------------------------------------------------------------------------
#
# Keyword calls go through a typed wrapper on purpose: inference through a
# closure over a script global sees `Any` and reports a false instability.

w_coarsen(cv, values, atol) = DGG.coarsen(cv, values; atol)
w_coarsen_dim(A, atol) = DGG.coarsen(A; atol)
w_cellvector(mov, l) = DGG.CellVector(mov; level = l)
w_mov(sys, cells, ref) = EN.MultiOrderVector(sys, cells; reference_level = ref)

# What each answer is allowed to be:
#   :concrete — one type, no exceptions;
#   :optional — `Union{T,Nothing}`, the "no such cell" contract;
#   :shape    — a documented two-way union over the answer's SHAPE (the window
#               form a `CellVector` picks, or container-or-plain-vector).
function verdict(rt, expected)
    isconcretetype(rt) && return ("concrete", true)
    rt === Any && return ("Any", false)
    if rt isa Union
        parts = Base.uniontypes(rt)
        if expected === :optional
            ok = length(parts) == 2 && Nothing in parts &&
                 all(isconcretetype, parts)
            return (ok ? "Union{T,Nothing}" : "$(length(parts))-way Union", ok)
        end
        ok = expected === :shape && length(parts) <= 2 &&
             all(isconcretetype, parts)
        return ("$(length(parts))-way Union", ok)
    end
    # A `(cells, values)` pair whose only looseness is the same shape union.
    if expected === :shape && rt isa DataType && rt <: Tuple
        ok = all(_shape_ok, rt.parameters)
        return ("Tuple of $(length(rt.parameters))", ok)
    end
    return ("abstract $(nameof(rt))", false)
end

_shape_ok(T) = isconcretetype(T) ||
               (T isa Union && length(Base.uniontypes(T)) <= 2 &&
                all(isconcretetype, Base.uniontypes(T)))

function inference_arm(cv, vals)
    rule("1. INFERENCE — return type of each entry point")
    coarse, cvals = DGG.coarsen(cv, vals; atol = 0.05)
    lk = DGG.MultiOrderLookup(coarse)
    A = DD.DimArray(cvals, (DGG.Cells(lk),))
    Acv = DD.DimArray(vals, (DGG.Cells(DGG.CellLookup(cv)),))
    _, edata = DGG.expand(coarse, cvals, DGG.level(cv))

    ID = DGG.LevelIndex
    TS, TM, TCV = typeof(SYS), typeof(coarse), typeof(cv)
    TLK, TA, TACV, TED = typeof(lk), typeof(A), typeof(Acv), typeof(edata)
    TP = typeof(DGG.cell_centroid(DGG.levelgrid(SYS, DGG.level(cv)), cv[1]))
    TR = typeof(REGION)

    probes = [
        ("MultiOrderVector(sys, cells)", EN.MultiOrderVector, Tuple{TS,Vector{ID}}, :concrete),
        ("MultiOrderVector(...; reference_level)", w_mov, Tuple{TS,Vector{ID},Int}, :concrete),
        ("MultiOrderVector(::MultiOrderCellSet)", EN.MultiOrderVector,
            Tuple{EN.MultiOrderCellSet{TS,ID}}, :concrete),
        ("getindex(mov, ::Int)", getindex, Tuple{TM,Int}, :concrete),
        ("getindex(mov, ::Vector{Int})", getindex, Tuple{TM,Vector{Int}}, :shape),
        ("iterate(mov)", iterate, Tuple{TM}, :optional),
        ("localindex(mov, cell)", DGG.localindex, Tuple{TM,ID}, :optional),
        ("localindex(mov, point)", DGG.localindex, Tuple{TM,TP}, :optional),
        ("covering_index(mov, cell)", DGG.covering_index, Tuple{TM,ID}, :optional),
        ("cellat(mov, point)", DGG.cellat, Tuple{TM,TP}, :optional),
        ("in(cell, mov)", in, Tuple{ID,TM}, :concrete),
        ("union(a, b)", union, Tuple{TM,TM}, :concrete),
        ("intersect(a, b)", intersect, Tuple{TM,TM}, :concrete),
        ("setdiff(a, b)", setdiff, Tuple{TM,TM}, :concrete),
        ("symdiff(a, b)", symdiff, Tuple{TM,TM}, :concrete),
        ("complement(mov)", EN.complement, Tuple{TM}, :concrete),
        ("issubset(a, b)", issubset, Tuple{TM,TM}, :concrete),
        ("isdisjoint(a, b)", isdisjoint, Tuple{TM,TM}, :concrete),
        ("covering(mov, region)", EN.covering, Tuple{TM,TR}, :concrete),
        ("covering_indices(mov, region)", DGG.covering_indices, Tuple{TM,TR}, :concrete),
        ("CellVector(mov; level)", w_cellvector, Tuple{TM,Int}, :concrete),
        ("aggregate(f, cv, vals, l)", DGG.aggregate,
            Tuple{typeof(sum),TCV,Vector{Float64},Int}, :shape),
        ("aggregate(f, A, l)", DGG.aggregate, Tuple{typeof(sum),TACV,Int}, :shape),
        ("coarsen(cv, vals; atol)", w_coarsen, Tuple{TCV,Vector{Float64},Float64}, :concrete),
        ("coarsen(A; atol)", w_coarsen_dim, Tuple{TACV,Float64}, :concrete),
        ("expand(mov, vals, l)", DGG.expand, Tuple{TM,Vector{Float64},Int}, :concrete),
        ("expand(A, l)", DGG.expand, Tuple{TA,Int}, :concrete),
        ("getindex(expanded, ::Int)", getindex, Tuple{TED,Int}, :concrete),
        ("collect(expanded)", collect, Tuple{TED}, :concrete),
        ("MultiOrderLookup(mov)", DGG.MultiOrderLookup, Tuple{TM}, :concrete),
        ("getindex(lookup, ::Int)", getindex, Tuple{TLK,Int}, :concrete),
        ("getindex(lookup, ::Vector{Int})", getindex, Tuple{TLK,Vector{Int}}, :shape),
        ("selectindices(lk, At(cell))", LK.selectindices, Tuple{TLK,LK.At{ID}}, :concrete),
        ("selectindices(lk, Contains(cell))", LK.selectindices,
            Tuple{TLK,LK.Contains{ID}}, :concrete),
        ("selectindices(lk, Contains(lon,lat))", LK.selectindices,
            Tuple{TLK,LK.Contains{Tuple{Float64,Float64}}}, :concrete),
        ("selectindices(lk, Covering(region))", LK.selectindices,
            Tuple{TLK,DGG.Covering{TR}}, :concrete),
        ("GR._asspace(mov, name)", GR._asspace, Tuple{TM,String}, :concrete),
        ("GR._asspace(lookup, name)", GR._asspace, Tuple{TLK,String}, :concrete),
    ]

    widths = (40, 20, 10)
    row(("entry point", "inferred", "ok?"), widths)
    println("-"^74)
    bad = String[]
    for (name, f, T, expected) in probes
        rt = Base.infer_return_type(f, T)
        (label, ok) = verdict(rt, expected)
        ok || push!(bad, "$name -> $rt")
        row((name, label, ok ? "yes" : "NO"), widths)
    end
    println()
    if isempty(bad)
        println("every entry point concrete, `Union{T,Nothing}`, or the " *
                "documented two-way shape union")
    else
        println("FAILURES:")
        foreach(b -> println("  ", b), bad)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# 2. Lookup scaling
# ---------------------------------------------------------------------------

function lookup_batch(f, mov, qs)
    t = 0
    for c in qs
        k = f(mov, c)
        k === nothing || (t += k)
    end
    return t
end

function lookup_arm(sizes)
    rule("2. LOOKUP — covering-ancestor cost against container size")
    widths = (12, 16, 12, 12, 12)
    row(("cells", "op", "ns/op", "bytes/op", "ns / log2 n"), widths)
    println("-"^70)
    for n in sizes
        mov = spread_container(n)
        qs = query_cells(n)
        q = length(qs)
        for (label, f) in (("localindex", DGG.localindex),
                           ("covering_index", DGG.covering_index))
            g = () -> lookup_batch(f, mov, qs)
            g()
            e = estimate(g; gc = false)
            per = e.ns / q
            row((n, label, @sprintf("%.1f", per), @sprintf("%.2f", e.bytes / q),
                @sprintf("%.2f", per / log2(n))), widths)
        end
        mov = qs = nothing
        GC.gc(true)
    end
    println("\n1000x the cells for ~5x the lookup. The right-hand column " *
            "rises rather than\nstaying flat because a longer binary search " *
            "misses cache more often, not because\nthe search grows with n.")
end

# ---------------------------------------------------------------------------
# 3. Set algebra
# ---------------------------------------------------------------------------

function algebra_arm(sizes)
    rule("3. SET ALGEBRA — merge cost against operand size")
    println("`a` holds every 2nd cell of level 9, `b` every 4th, so b < a and " *
            "every answer\nis proportional to n rather than collapsing to a " *
            "handful of coarse cells.")
    widths = (10, 12, 11, 12, 12, 12)
    row(("cells", "op", "out cells", "time", "bytes", "ns/out"), widths)
    println("-"^76)
    for n in sizes
        a = spread_container(n; step = 2)
        b = spread_container(fld(n, 2); step = 4)
        for (label, f) in (("union", union), ("intersect", intersect),
                           ("setdiff", setdiff), ("symdiff", symdiff))
            out = f(a, b)
            g = () -> f(a, b)
            e = estimate(g)
            m = max(length(out), 1)
            row((n, label, length(out), fmt_ns(e.ns), fmt_bytes(e.bytes),
                @sprintf("%.1f", e.ns / m)), widths)
        end
        for (label, f) in (("issubset", issubset), ("isdisjoint", isdisjoint),
                           ("complement", EN.complement))
            g = label == "complement" ? (() -> f(a)) : (() -> f(a, b))
            g()
            e = estimate(g)
            row((n, label, "-", fmt_ns(e.ns), fmt_bytes(e.bytes), "-"), widths)
        end
        a = b = nothing
        GC.gc(true)
    end
end

# ---------------------------------------------------------------------------
# 4. Hierarchy verbs
# ---------------------------------------------------------------------------

function hierarchy_arm(cv, vals)
    leaf = DGG.level(cv)
    rule("4. HIERARCHY — aggregate, coarsen and construction on level $leaf")
    n = length(cv)
    widths = (26, 12, 12, 14, 12)
    row(("op", "leaf cells", "time", "leaf cells/s", "bytes"), widths)
    println("-"^80)

    for l in (leaf - 1, leaf - 3)
        g = () -> DGG.aggregate(sum, cv, vals, l)
        g()
        e = estimate(g; seconds = 3.0, samples = 5)
        row(("aggregate -> level $l", n, fmt_ns(e.ns),
            @sprintf("%.0f M", n / (e.ns / 1e9) / 1e6), fmt_bytes(e.bytes)), widths)
    end

    for atol in (0.001, 0.05)
        g = () -> DGG.coarsen(cv, vals; atol)
        mov, _ = g()
        e = estimate(g; seconds = 3.0, samples = 5)
        row(("coarsen atol=$atol", n, fmt_ns(e.ns),
            @sprintf("%.0f M", n / (e.ns / 1e9) / 1e6), fmt_bytes(e.bytes)), widths)
        println("    -> $(length(mov)) stored cells, levels " *
                "$(minimum(DGG.level, mov)):$(maximum(DGG.level, mov)), " *
                "$(round(100 * length(mov) / n; digits = 2))% of the leaves")
    end

    # Construction, the other superlinear-looking step. Ascending input is the
    # case every caller inside the package hands it; shuffled input pays a sort.
    g9 = DGG.levelgrid(SYS, 9)
    for n2 in (10_000, 100_000)
        asc = [DGG.cellindex(g9, 2k - 1) for k in 1:n2]
        shuffled = asc[vcat(2:2:n2, 1:2:n2)]
        for (label, cells) in (("construct sorted", asc),
                               ("construct shuffled", shuffled))
            g = () -> EN.MultiOrderVector(SYS, cells; reference_level = 9)
            g()
            e = estimate(g; seconds = 2.0, samples = 8)
            row(("$label $n2", n2, fmt_ns(e.ns),
                @sprintf("%.0f M", n2 / (e.ns / 1e9) / 1e6), fmt_bytes(e.bytes)),
                widths)
        end
    end
    GC.gc()
end

# ---------------------------------------------------------------------------
# 5. Expand
# ---------------------------------------------------------------------------

function access_batch(v, n)
    t = zero(eltype(v))
    for k in 1:n
        t += v[k]
    end
    return t
end

function expand_arm(cv, vals)
    leaf = DGG.level(cv)
    rule("5. EXPAND — the lazy presentation at level $leaf")
    mov, mvals = DGG.coarsen(cv, vals; atol = 0.05)
    ecv, edata = DGG.expand(mov, mvals, leaf)
    stored, leaves = length(mov), length(edata)
    println("$stored stored values presented as $leaves cells " *
            "($(round(leaves / stored; digits = 1))x); summarysize ",
        fmt_bytes(Base.summarysize(mvals)), " stored vs ",
        fmt_bytes(Base.summarysize(edata)), " presented")
    widths = (26, 12, 12, 14)
    row(("op", "time", "bytes", "per element"), widths)
    println("-"^68)

    q = min(100_000, leaves)
    for (label, g, per) in (
        ("getindex x $q", () -> access_batch(edata, q), q),
        ("_leaf_offsets ($stored)", () -> CL._leaf_offsets(mov, leaf), stored),
        ("collect ($leaves)", () -> collect(edata), leaves),
        ("expand (build)", () -> DGG.expand(mov, mvals, leaf), 0))
        g()
        e = estimate(g; seconds = 2.0, samples = 10)
        row((label, fmt_ns(e.ns), fmt_bytes(e.bytes),
            per == 0 ? "-" : @sprintf("%.2f ns", e.ns / per)), widths)
    end
    GC.gc()
    return mov, mvals, ecv, edata
end

# ---------------------------------------------------------------------------
# The region bound, which is not the container's
# ---------------------------------------------------------------------------

function region_arm()
    rule("R. REGION — `covering_indices` is bounded by the coverage query")
    println("`covering_indices(mov, target)` is one binary search per coverage " *
            "interval over\nthe container, but it first builds the coverage at " *
            "the container's reference level.\nThat second cost is " *
            "`MultiOrderCoverage`'s, not the container's, and it dominates:")
    widths = (16, 14, 14, 14)
    row(("query level", "set cells", "time", "bytes"), widths)
    println("-"^62)
    for l in 5:7
        g = () -> DGG.query(SYS, DGG.MultiOrderCoverage(REGION); level = l)
        set = g()
        e = estimate(g; seconds = 2.0, samples = 3)
        row((l, length(set.cells), fmt_ns(e.ns), fmt_bytes(e.bytes)), widths)
        set = nothing
        GC.gc(true)
    end
    mov = spread_container(20_000; level = 6, step = 2)
    g = () -> DGG.covering_indices(mov, REGION)
    ks = g()
    e = estimate(g; seconds = 2.0, samples = 3)
    println("\ncovering_indices(20 000 cells, ref 6) -> $(length(ks)) cells, ",
        fmt_ns(e.ns), ", ", fmt_bytes(e.bytes),
        "\n  — the same order as the level-6 query row above, so the container " *
        "adds\n    essentially nothing to it.")
    mov = nothing
    GC.gc(true)
end

# ---------------------------------------------------------------------------
# 6. Store IO
# ---------------------------------------------------------------------------

dirsize(path) = sum(filesize(joinpath(root, f))
                    for (root, _, files) in walkdir(path) for f in files; init = 0)

function io_arm(leaf::Int)
    rule("6. STORE IO — compacted against the expanded write of the same field")
    if Base.get_extension(DGG, :DiscreteGlobalGridsZarrExt) === nothing
        println("Zarr is unavailable (it is a weak dependency); arm skipped. " *
                "Run under --project=benchmark.")
        return
    end
    cv, vals = level_field(leaf)
    mov, mvals = DGG.coarsen(cv, vals; atol = 0.05)
    M = DD.DimArray(mvals, (DGG.Cells(DGG.MultiOrderLookup(mov)),); name = :field)
    ecv, edata = DGG.expand(mov, mvals, leaf)
    E = DD.DimArray(collect(edata), (DGG.Cells(DGG.CellLookup(ecv)),); name = :field)
    println("the same level-$leaf field: $(length(mov)) compacted cells against " *
            "$(length(ecv)) expanded ones")

    root = mktempdir()
    widths = (14, 12, 12, 12, 14)
    row(("store", "cells", "write", "read", "on disk"), widths)
    println("-"^68)
    results = Tuple{String,Float64,Float64,Int}[]
    for (label, A) in (("compacted", M), ("expanded", E))
        # `dggwrite` refuses a non-empty store, so every sample gets its own
        # directory. Both arms pay that same fixed cost.
        k = Ref(0)
        gw = () -> DGG.dggwrite(joinpath(root, "$label-$(k[] += 1).zarr"), A)
        gw()
        ew = estimate(gw; seconds = 3.0, samples = 3)
        path = joinpath(root, "$label-1.zarr")
        gr = () -> DGG.dggread(path)
        S = gr()
        er = estimate(gr; seconds = 3.0, samples = 3)
        sz = dirsize(path)
        row((label, length(DD.lookup(A, DGG.Cells)), fmt_ns(ew.ns),
            fmt_ns(er.ns), fmt_bytes(sz)), widths)
        push!(results, (label, ew.ns, er.ns, sz))
        # The round trip has to agree, or the comparison means nothing.
        back = S[:field]
        length(back) == length(A) ||
            error("$label round trip changed the cell count")
        S = back = nothing
        GC.gc(true)
    end
    c, e = results
    @printf("\ncompacted / expanded: write %.2fx, read %.2fx, on disk %.3fx\n",
        c[2] / e[2], c[3] / e[3], c[4] / e[4])
    rm(root; recursive = true, force = true)
    GC.gc()
end

# ---------------------------------------------------------------------------
# 7. Regrid
# ---------------------------------------------------------------------------
#
# An interpolating method reads sample sites, so it takes the stored cells as
# they are. The comparison is against the only other way to interpolate a
# mixed-level field — expand to the reference level first, where every leaf
# under a stored cell repeats one value — measured in the same session on the
# same numbers, so only the source presentation differs.

# Typed wrappers, as arm 1: a keyword call through a closure over a script
# global reports a false answer for both inference and allocation.
w_plan(A, dst, m) = GR.plan_regrid(A; to = dst, method = m, lazy = false)
w_apply(A, plan) = GR.regrid(A, plan)
w_cube(mov, values) = DD.DimArray(values, (DGG.Cells(DGG.MultiOrderLookup(mov)),))

function w_sweep(row, smp, pts)
    n = 0
    for p in pts
        GR.ismapped(GR.weightsat!(row, smp, p)) && (n += 1)
    end
    return n
end

function regrid_arm(leaf::Int, dstlevel::Int; atol = 0.25)
    rule("7. REGRID — BarycentricPoint off the stored cells, and off their expansion")
    cv, vals = level_field(leaf)
    mov, mvals = DGG.coarsen(cv, vals; atol)
    cv = vals = nothing
    GC.gc(true)
    A = w_cube(mov, mvals)
    E = DGG.expand(A, leaf)
    dst = DGG.levelgrid(SYS, dstlevel)
    m = GR.BarycentricPoint()
    stored, leaves = length(mov), DGG.ncells(DGG.levelgrid(SYS, leaf))
    println("a level-$leaf field coarsened at atol $atol to $stored stored " *
            "cells ($(round(leaves / stored; digits = 1))x fewer than its " *
            "$leaves leaves),\nonto level $dstlevel ($(DGG.ncells(dst)) cells)")

    widths = (14, 9, 9, 11, 12, 11)
    row(("route", "source", "columns", "build", "build bytes", "apply"), widths)
    println("-"^76)
    results = Dict{String,NTuple{3,Float64}}()
    for (label, data, apply) in (("stored", A, true), ("expanded", E, true))
        gb = () -> w_plan(data, dst, m)
        plan = gb()
        eb = estimate(gb; seconds = 3.0, samples = 5)
        ns = 0.0
        if apply
            ga = () -> w_apply(data, plan)
            ga()
            ns = estimate(ga; seconds = 3.0, samples = 5).ns
        end
        row((label, Int(GR.ncells(plan.src_space)),
                size(plan.block.weights, 2), fmt_ns(eb.ns), fmt_bytes(eb.bytes),
                apply ? fmt_ns(ns) : "-"), widths)
        results[label] = (eb.ns, Float64(eb.bytes), ns)
        plan = nothing
        GC.gc(true)
    end

    # The reference level names the leaves and nothing else. Re-keying the SAME
    # stored cells two levels deeper multiplies the leaf count by sixteen: the
    # native plan should not move at all, and the expansion's should follow the
    # leaves. Only the builds are measured here — applying the deep expansion
    # would materialise a leaf per value and say nothing new.
    deep = EN.MultiOrderVector(SYS, collect(mov); reference_level = leaf + 2)
    D = w_cube(deep, mvals)
    deepleaves = DGG.ncells(DGG.levelgrid(SYS, leaf + 2))
    println("\nthe same $stored stored cells re-keyed to reference level " *
            "$(leaf + 2) — $deepleaves leaves, $(round(Int, deepleaves / leaves))x:")
    for (label, data) in (("stored/deep", D), ("expanded/deep", DGG.expand(D, leaf + 2)))
        gd = () -> w_plan(data, dst, m)
        plan = gd()
        e = estimate(gd; seconds = 3.0, samples = 3)
        row((label, Int(GR.ncells(plan.src_space)),
                size(plan.block.weights, 2), fmt_ns(e.ns), fmt_bytes(e.bytes),
                "-"), widths)
        results[label] = (e.ns, Float64(e.bytes), 0.0)
        plan = nothing
        GC.gc(true)
    end
    D = deep = nothing
    GC.gc(true)

    s, e = results["stored"], results["expanded"]
    sd, ed = results["stored/deep"], results["expanded/deep"]
    @printf("\nstored / expanded, at reference level %d: build %.2fx time, %.2fx bytes; apply %.3fx\n",
        leaf, s[1] / e[1], s[2] / e[2], s[3] / e[3])
    @printf("stored / expanded, at reference level %d: build %.2fx time, %.3fx bytes\n",
        leaf + 2, sd[1] / ed[1], sd[2] / ed[2])
    @printf("deep / shallow, same stored cells:  stored %.2fx time %.2fx bytes, expanded %.2fx time %.2fx bytes\n",
        sd[1] / s[1], sd[2] / s[2], ed[1] / e[1], ed[2] / e[2])

    # And one query: the fan is fixed capacity, so a located stencil reaches no
    # heap however big the container is.
    smp = GR.sampler(m, GR.sourcespacefor(mov, m))
    pts = [GR.cellcentroid(DGG.DGGSpace(dst), i) for i in 1:DGG.ncells(dst)]
    grow = () -> w_sweep(GR.WeightRow(), smp, pts)
    mapped = grow()
    eq = estimate(grow; seconds = 3.0, samples = 5)
    @printf("weightsat! over %d points: %s each, %s total, %d mapped\n",
        length(pts), fmt_ns(eq.ns / length(pts)), fmt_bytes(eq.bytes), mapped)
    A = E = mov = mvals = smp = pts = nothing
    GC.gc(true)
end

# ---------------------------------------------------------------------------

rss() = round(Sys.maxrss() / 2^20; digits = 1)
mark(name) = println("\n[peak RSS after $name: $(rss()) MiB]")

function main()
    println("DiscreteGlobalGrids multi-order baseline — Julia ", Base.VERSION,
        ", ", Threads.nthreads(), " thread(s)")
    mark("load and compile")
    small_cv, small_vals = level_field(6)
    inference_arm(small_cv, small_vals)
    small_cv = small_vals = nothing
    GC.gc(true)
    mark("arm 1")

    lookup_arm((1_000, 10_000, 100_000, 1_000_000))
    GC.gc(true)
    mark("arm 2")
    algebra_arm((1_000, 10_000, 100_000))
    GC.gc(true)
    mark("arm 3")

    cv, vals = level_field(8)
    hierarchy_arm(cv, vals)
    mark("arm 4")
    expand_arm(cv, vals)
    cv = vals = nothing
    GC.gc(true)
    mark("arm 5")

    region_arm()
    mark("region")
    regrid_arm(6, 5)
    GC.gc(true)
    mark("arm 7")
    io_arm(6)
    println("\npeak RSS ", rss(), " MiB")
end

main()
