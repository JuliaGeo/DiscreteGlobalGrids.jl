# ---------------------------------------------------------------------------
# `halo(subset)` — what the walk costs, measured against an untouched control.
#
#     julia --project=test scripts/bench/halo_subset_scaling.jl
#
# PUBLIC API ONLY, so the same file runs unchanged on any commit and the two
# runs are comparable. That is the whole point of it: every number this package
# states about `halo`'s cost is a number this script prints.
#
# WHY THERE IS A CONTROL. `subtree_halo` on a HEALPix block takes the square
# band walk, which no change to the subset path touches, so its time is this
# machine's clock rather than the measurement. Absolute milliseconds move with
# the machine's power state and with whatever else is running on it; the RATIO
# to the control does not. Read the `xctl` column, not the `ms` column.
#
# THREE CASES:
#
#   * THE BATCH CLIFF. A subset one cell over a bounding-cap batch limit used to
#     fall back to a full-sphere cap, which prunes nothing, and the walk then
#     descended the whole hierarchy. 2048 against 2049 members is the same
#     question asked twice with the same answer either side.
#   * THE SCALING. A rooted subtree with one interior cell punched out is the
#     irregular chunk `halo` exists for — the hole law is exactly what
#     `subtree_halo` cannot express. Its member count quadruples per level and
#     its halo doubles, so the two exponents below separate cleanly: a walk
#     sized by the input reads ~2.0 in the halo, one sized by the answer ~1.0.
#   * THE SCATTERED SUBSET, which has no boundary to follow and is therefore the
#     case a boundary-following walk could plausibly lose on. It is here so that
#     the report is the whole cost and not the favourable half of it.
# ---------------------------------------------------------------------------

using DiscreteGlobalGrids
using Printf

best(f, reps = 3) = minimum(begin
    t0 = time_ns()
    f()
    (time_ns() - t0) / 1e6
end for _ in 1:reps)

sys = HEALPixSystem()
root = cellindex(levelgrid(sys, 0), 1)

# The control: untouched by anything the subset path does.
control = best(() -> subtree_halo(sys, root, 10))
@printf("control  subtree_halo(HEALPix, level-0 root, l = 10)  %8.3f ms\n\n", control)

# --- the batch cliff -------------------------------------------------------

println("the batch cliff — a contiguous block of the level-8 grid")
@printf("%10s %10s %10s %8s\n", "members", "halo", "ms", "xctl")
grid8 = levelgrid(sys, 8)
for n in (2047, 2048, 2049, 4096)
    pg = PartialGrid(sys, 8, [cellindex(grid8, p) for p in 1:n])
    h = length(collect(halo(pg)))
    ms = best(() -> collect(halo(pg)))
    @printf("%10d %10d %10.3f %8.1f\n", n, h, ms, ms / control)
end

# --- the scaling law -------------------------------------------------------

function holed(sys, c, l)
    r = descendant_range(sys, c, l)
    g = levelgrid(sys, l)
    ids = [cellindex(g, p) for p in r]
    deleteat!(ids, length(ids) ÷ 2)          # one interior cell, the hole law
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

# --- the case with no boundary to follow -----------------------------------
#
# A subset scattered uniformly over the whole sphere has a halo of the same
# order as the level, so `O(halo)` and `O(ncells)` are the same statement and
# there is nothing to win.

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
