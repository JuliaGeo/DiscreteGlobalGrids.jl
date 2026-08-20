# Benchmark `halo(subset)` against a stable rooted-subtree control.
#
#     julia --project=test scripts/bench/halo_subset_scaling.jl
#
# The script uses only public APIs so results remain comparable across commits.
# The `xctl` column normalizes timings by a HEALPix square-band walk and is more
# stable across machines than the absolute millisecond columns.
#
# The cases cover the former 2,048-member batching boundary, a subtree with an
# interior hole, and a subset scattered across the sphere.
#
# For the holed subtree, members grow quadratically with the perimeter-sized
# halo. Fitted exponents therefore distinguish input-sized from halo-sized work.

using DiscreteGlobalGrids
using Printf

best(f, reps = 3) = minimum(begin
    t0 = time_ns()
    f()
    (time_ns() - t0) / 1e6
end for _ in 1:reps)

sys = HEALPixSystem()
root = cellindex(levelgrid(sys, 0), 1)

# Stable control that does not use the subset path.
whole = subtree(sys, root, 10)
control = best(() -> collect(halo(whole)))
@printf("control  halo(subtree(HEALPix, level-0 root, l = 10))  %8.3f ms\n\n", control)

# Former batching boundary.

println("the batch cliff — a contiguous block of the level-8 grid")
@printf("%10s %10s %10s %8s\n", "members", "halo", "ms", "xctl")
grid8 = levelgrid(sys, 8)
for n in (2047, 2048, 2049, 4096)
    pg = PartialGrid(sys, 8, [cellindex(grid8, p) for p in 1:n])
    h = length(collect(halo(pg)))
    ms = best(() -> collect(halo(pg)))
    @printf("%10d %10d %10.3f %8.1f\n", n, h, ms, ms / control)
end

# Scaling for a subtree with one interior cell removed.

function holed(sys, c, l)
    r = descendant_range(sys, c, l)
    g = levelgrid(sys, l)
    ids = [cellindex(g, p) for p in r]
    deleteat!(ids, length(ids) ÷ 2)          # create an interior halo component
    return PartialGrid(sys, l, ids; root = c)
end

println("\nthe scaling law — a rooted level-0 subtree with one cell punched out")
@printf("%4s %10s %8s %10s %8s %10s %8s\n",
    "l", "members", "halo", "grid ms", "xctl", "vector ms", "xctl")
ns = Int[]; halos = Int[]; tg = Float64[]; tv = Float64[]
for l in 5:10
    pg = holed(sys, root, l)
    cv = CellVector(pg)
    h = length(collect(halo(pg)))
    mg = best(() -> collect(halo(pg)))
    mv = best(() -> collect(halo(cv)))
    push!(ns, length(pg.ids)); push!(halos, h); push!(tg, mg); push!(tv, mv)
    @printf("%4d %10d %8d %10.3f %8.1f %10.3f %8.1f\n",
        l, length(pg.ids), h, mg, mg / control, mv, mv / control)
end

slope(x, y) = ((n = length(x); mx = sum(log, x) / n; my = sum(log, y) / n);
sum((log(x[i]) - mx) * (log(y[i]) - my) for i in 1:n) /
sum((log(x[i]) - mx)^2 for i in 1:n))

@printf("\nfitted exponents   %-8s %8s %8s\n", "", "members", "halo")
@printf("  PartialGrid      %-8s %8.3f %8.3f\n", "", slope(ns, tg), slope(halos, tg))
@printf("  CellVector       %-8s %8.3f %8.3f\n", "", slope(ns, tv), slope(halos, tv))
println("(a walk sized by the input reads ~1.0 and ~2.0; one sized by the")
println(" answer reads ~0.5 and ~1.0 — the halo of a block is its perimeter.")
println(" The two containers differ by their membership search: a `CellVector`")
println(" keeps a holed subtree as two windows, so its search is O(1) here,")
println(" while a `PartialGrid` binary-searches an id vector that is growing.)")

# Scattered subsets have no compact boundary to follow.
# A subset scattered uniformly over the whole sphere has a halo of the same
# order as the level grid, making `O(halo)` and `O(ncells)` equivalent here.

using Random

println("\nthe scattered subset — uniformly random cells of the level-7 grid")
@printf("%10s %10s %10s %10s %8s\n", "members", "halo", "grid ms", "vector ms", "xctl")
grid7 = levelgrid(sys, 7)
rng = MersenneTwister(7)
for frac in (0.02, 0.1, 0.5)
    k = round(Int, ncells(grid7) * frac)
    ps = sort(randperm(rng, ncells(grid7))[1:k])
    pg = PartialGrid(sys, 7, [cellindex(grid7, p) for p in ps])
    cv = CellVector(pg)
    h = length(collect(halo(cv)))
    mg = best(() -> collect(halo(pg)))
    mv = best(() -> collect(halo(cv)))
    @printf("%10d %10d %10.3f %10.3f %8.1f\n", k, h, mg, mv, mv / control)
end
