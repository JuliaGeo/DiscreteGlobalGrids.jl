# Demo: a DimensionalData cell axis backed by a multi-order coverage.
#
# A regional cube wants ONE dimension naming a region's cells at a leaf level.
# Materialising that dimension is the naive move: a Switzerland-sized region at
# IGEO7 level 12 is millions of ids. The compact form already exists —
#
#     set    = DGG.query(sys, DGG.MultiOrderCoverage(region); level = leaf)
#     ranges = DGG.level_ranges(set, leaf)
#
# `set` is a few hundred mixed-level cells; `ranges` is their expansion to
# sorted, disjoint POSITION ranges at the leaf level. Their concatenation is
# exactly the leaf id vector, so a DD lookup can be O(#ranges) in memory and
# still answer `length`, `getindex` and every selector.
#
# The first half of this file is that arithmetic, done by hand, with the laws
# the lookup has to satisfy asserted. The second half is the API those laws
# should be hidden behind, written out as the acceptance test for T16 — it is
# COMMENTED OUT because nothing implements it yet.
#
# Environment: needs nothing beyond DiscreteGlobalGrids and its dependencies:
#
#     julia -t 4 --project=. examples/dimensionaldata.jl
#
# It ends in PASS/FAIL assertions and exits non-zero if any of them fail.

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GeoInterface as GI

const FAILURES = Ref(0)
function check(name, ok; detail="")
    ok || (FAILURES[] += 1)
    println(ok ? "PASS  " : "FAIL  ", rpad(name, 56), detail)
    return ok
end
note(text) = println("      ", text)

const SYS = DGG.IGeo7System()
const LEAF = 9
# A Switzerland-shaped box. Any GeoInterface polygon works; a rectangle keeps
# the demo's output stable.
const REGION = GI.Polygon([GI.LinearRing([(6.0, 45.8), (10.5, 45.8), (10.5, 47.8),
    (6.0, 47.8), (6.0, 45.8)])])

println("="^78)
println("dimensionaldata.jl — a leaf-level cell axis, held as a multi-order coverage")
println("julia $(VERSION)  threads=$(Threads.nthreads())")
println("="^78)

# --------------------------------------------------------------------------
# 1. The compact backing.
# --------------------------------------------------------------------------

set = DGG.query(SYS, DGG.MultiOrderCoverage(REGION); level=LEAF)
ranges = DGG.level_ranges(set, LEAF)
nleaf = sum(length, ranges)

check("the coverage is mixed-level", length(unique(DGG.level.(collect(set)))) > 1;
    detail="$(length(set)) cells over levels " *
           "$(minimum(DGG.level, set)):$(maximum(DGG.level, set))")
check("ranges are sorted and disjoint",
    issorted(first.(ranges)) && all(last(ranges[i]) < first(ranges[i+1])
                                    for i in 1:length(ranges)-1);
    detail="$(length(ranges)) ranges, $nleaf leaf cells")
note("compression: $(length(ranges)) stored ranges for $nleaf leaf cells " *
     "($(round(nleaf / length(ranges); digits=1))x)")

# --------------------------------------------------------------------------
# 2. The arithmetic a lookup would hide.
#
# The logical content is the concatenation of the ranges. Position `k` in that
# concatenation maps to a leaf POSITION by a binary search over the cumulative
# range lengths, and only then to a typed id — nothing is materialised.
# --------------------------------------------------------------------------

const LEAFGRID = DGG.levelgrid(SYS, LEAF)
const OFFSETS = cumsum(length.(ranges))          # offsets[j] = cells up to range j

function leaf_position(k::Int)
    1 <= k <= nleaf || throw(BoundsError(ranges, k))
    j = searchsortedfirst(OFFSETS, k)
    return first(ranges[j]) + (k - (j == 1 ? 0 : OFFSETS[j-1])) - 1
end

leaf_id(k::Int) = DGG.cellindex(LEAFGRID, leaf_position(k))

function concat_position(c::DGG.AbstractCellIndex)
    p = DGG.cellposition(LEAFGRID, c)
    p === nothing && return nothing
    j = searchsortedfirst(last.(ranges), p)
    (j <= length(ranges) && p in ranges[j]) || return nothing
    return (j == 1 ? 0 : OFFSETS[j-1]) + (p - first(ranges[j])) + 1
end

# LAW 1 — position <-> id round trips.
check("position -> id -> position round trips",
    all(concat_position(leaf_id(k)) == k
        for k in (1, 2, nleaf ÷ 3, nleaf ÷ 2, nleaf)))
# LAW 2 — the lazy form equals the materialised leaf vector.
materialised = DGG.cellindex.(Ref(LEAFGRID), reduce(vcat, collect.(ranges)))
check("lazy ids == the materialised leaf vector",
    length(materialised) == nleaf && all(leaf_id(k) == materialised[k] for k in 1:nleaf))
# LAW 3 — a cell outside the region has no position.
outside = DGG.cellat(LEAFGRID, -60.0, 0.0)
check("a cell outside the coverage has no position",
    concat_position(outside) === nothing)

# --------------------------------------------------------------------------
# 3. Selectors, also by hand.
#
# A point selector is `cellat` on the leaf grid followed by the same search; a
# polygon selector is a second coverage intersected with the first.
# --------------------------------------------------------------------------

k = concat_position(DGG.cellat(LEAFGRID, 8.0, 46.5))
check("a lon/lat point selects one position", k isa Int;
    detail="(8.0, 46.5) -> position $k -> $(leaf_id(k))")

zurich = GI.Polygon([GI.LinearRing([(8.4, 47.3), (8.7, 47.3), (8.7, 47.5),
    (8.4, 47.5), (8.4, 47.3)])])
sub = DGG.level_ranges(DGG.query(SYS, DGG.MultiOrderCoverage(zurich); level=LEAF), LEAF)
selected = Int[k for r in sub for p in r
               for k in (concat_position(DGG.cellindex(LEAFGRID, p)),) if k !== nothing]
check("a polygon selects a contiguous run of positions",
    !isempty(selected) && issorted(selected);
    detail="$(length(selected)) of $nleaf positions")

# --------------------------------------------------------------------------
# 4. What DimensionalData can do with this TODAY.
#
# A plain id vector is a valid DD lookup, so `At` works — and that is the whole
# of it. The vector is O(leaf cells), which is exactly the cost the multi-order
# backing exists to avoid, and no other selector knows what a cell is.
# --------------------------------------------------------------------------

A = DD.DimArray(Float64.(1:nleaf), (DD.Dim{:cell}(materialised),); name=:demo)
check("DD.At(id) resolves through a materialised id vector",
    A[cell=DD.At(materialised[k])] == Float64(k))
note("that dimension holds $nleaf ids; the coverage behind it holds $(length(set))")

# ==========================================================================
# ASPIRATIONAL — the API this file wants, as the acceptance test for T16.
# NOTHING BELOW THIS LINE RUNS. Every line above is the same computation done
# by hand, so the semantics are pinned; what is missing is the type.
# ==========================================================================
#
# ```julia
# # -- construction ------------------------------------------------------
# # MOC-backed: memory is O(#coverage entries), never O(#leaf cells).
# lk = DGG.CellLookup(set)                          # leaf level = set's reference level
# lk = DGG.CellLookup(set; level = LEAF)            # or expand to a deeper leaf
#
# # The two degenerate cases fall out of the same type:
# lk = DGG.CellLookup(DGG.levelgrid(SYS, LEAF))     # one full-level range
# lk = DGG.CellLookup(DGG.PartialGrid(SYS, LEAF, ids))  # an explicit id list
#
# # -- the collection surface, all lazy -----------------------------------
# length(lk) == nleaf                               # sum of the range lengths
# lk[k]                                             # position -> typed leaf id
# DGG.cellposition(lk, c)                           # typed leaf id -> position, or nothing
# collect(lk) == materialised                       # the equivalence law
# DGG.level(lk) == LEAF                             # every cell in it is a leaf
# parent(lk) === set                                # the backing, for a second coverage op
#
# # -- as a DimensionalData dimension --------------------------------------
# A = DD.DimArray(zeros(length(lk)), DGG.Cells(lk); name = :dem)
#
# A[DGG.Cells(DD.At(c))]                            # typed id      -> one value
# A[DGG.Cells(DD.Contains(8.0, 46.5))]              # lon/lat point -> `cellat` -> one value
# A[DGG.Cells(DGG.Covering(zurich))]                # polygon       -> T15 coverage ∩ backing
# A[DGG.Cells(DGG.Covering(cap))]                   # any `query` target, same path
#
# # A selector returns a view over a NEW `CellLookup` whose backing is the
# # intersected coverage — still O(#entries), so subsetting never materialises.
# sub = A[DGG.Cells(DGG.Covering(zurich))]
# DD.lookup(sub, DGG.Cells) isa DGG.CellLookup
#
# # -- the grid a lookup describes ----------------------------------------
# # `PartialGrid` and `CellLookup` are the same set seen from two sides, so the
# # conversion is free and the regridder can consume the cube's own axis:
# grid = DGG.PartialGrid(lk)
# regridder = CR.Regridder(MANIFOLD, destination, grid)
# CR.regrid!(out, regridder, parent(A))             # parent(A)[i] IS grid position i
#
# # -- A5, which has no descendant ranges ----------------------------------
# # `level_ranges` throws there, so `CellLookup` must fall back to explicit
# # sorted leaf positions per entry (selection mode) rather than ranges. Same
# # type, same laws, different backing — or a documented refusal.
# DGG.CellLookup(DGG.query(DGG.A5System(), DGG.MultiOrderCoverage(REGION); level = 7))
# ```

println()
note("the block at the bottom of this file is the T16 acceptance test")
note("today's spelling: level_ranges(set, leaf) + searchsortedfirst, by hand")

println()
println(FAILURES[] == 0 ? "ALL CHECKS PASSED" : "$(FAILURES[]) CHECK(S) FAILED")
exit(FAILURES[] == 0 ? 0 : 1)
