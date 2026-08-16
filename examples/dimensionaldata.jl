# Demo: a DimensionalData cell axis backed by a multi-order coverage.
#
# A regional cube wants ONE dimension naming a region's cells at a leaf level.
# Materialising that dimension is the naive move: a Switzerland-sized region at
# IGEO7 level 12 is millions of ids. The compact form already exists —
#
#     set    = DGG.query(sys, DGG.MultiOrderCoverage(region); level = leaf)
#     ranges = DGG.level_ranges(set, leaf)
#
# `set` is a few thousand mixed-level cells; `ranges` is their expansion to
# sorted, disjoint POSITION ranges at the leaf level. Their concatenation is
# exactly the leaf id vector, so a DD lookup can be O(#ranges) in memory and
# still answer `length`, `getindex` and every selector.
#
# `DGG.CellLookup` is that lookup and this file is its acceptance test. It was
# written before the type existed, as twenty-five lines of `cumsum` and
# `searchsortedfirst` over `level_ranges`, with the laws asserted so that the
# implementation had a pinned specification; every one of those laws is still
# below, now asked of the type instead of of the arithmetic.
#
# Environment: needs nothing beyond DiscreteGlobalGrids and its dependencies:
#
#     julia -t 4 --project=. examples/dimensionaldata.jl
#
# It ends in PASS/FAIL assertions and exits non-zero if any of them fail.

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GlobalRegridding as GR
import GeometryOps as GO
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
# 2. The lookup.
#
# `CellLookup(set)` is semantically the leaf id vector and structurally the
# ranges. Position `k` maps to a leaf POSITION by a binary search over the
# cumulative range lengths, and only then to a typed id — nothing is
# materialised, which is what the memory section below measures.
# --------------------------------------------------------------------------

const LEAFGRID = DGG.levelgrid(SYS, LEAF)

lk = DGG.CellLookup(set)                        # leaf level = set's reference level
check("the lookup is the leaf id vector's length", length(lk) == nleaf;
    detail="$(length(lk)) cells")
check("every cell in it is a leaf", DGG.level(lk) == LEAF &&
                                    all(DGG.level(lk[k]) == LEAF for k in (1, nleaf ÷ 2, nleaf)))
check("the backing is the set itself", DGG.cellset(lk) === set;
    detail="cellset(lk) === set, for a second coverage op")
# `parent` is the VALUES, which is what DimensionalData reads: some thirty of
# its `Lookup` methods derive their behaviour from it, `Where` among them.
check("parent is the lazy id vector",
    parent(lk) isa AbstractVector{eltype(lk)} && length(parent(lk)) == nleaf)

# LAW 1 — position <-> id round trips.
check("position -> id -> position round trips",
    all(DGG.cellposition(lk, lk[k]) == k
        for k in (1, 2, nleaf ÷ 3, nleaf ÷ 2, nleaf)))
# LAW 2 — the lazy form equals the materialised leaf vector.
materialised = DGG.cellindex.(Ref(LEAFGRID), reduce(vcat, collect.(ranges)))
check("lazy ids == the materialised leaf vector",
    length(materialised) == nleaf && collect(lk) == materialised)
# LAW 3 — a cell outside the region has no position.
outside = DGG.cellat(LEAFGRID, -60.0, 0.0)
check("a cell outside the coverage has no position",
    DGG.cellposition(lk, outside) === nothing)

# The same type, re-expanded. The set's reference level is 9; asking for 12
# names 343 times as many cells and stores the same 666 windows.
deep = DGG.CellLookup(set; level=12)
check("one set, any leaf level",
    length(deep) == 343 * nleaf;
    detail="level 12: $(length(deep)) cells")

# --------------------------------------------------------------------------
# 3. The cube, and the three questions a cell axis is asked.
#
# `At` and `Contains` stay spelled `DD.`-qualified: they are DimensionalData's
# selectors and this package exports its own `Contains`, the DE9IM predicate,
# which would collide with the selector under a plain `using`.
# --------------------------------------------------------------------------

A = DD.DimArray(Float64.(1:nleaf), DGG.Cells(lk); name=:dem)
check("a DimArray over the lookup", DD.lookup(A, DGG.Cells) === lk;
    detail=sprint(show, DD.lookup(A, DGG.Cells)))

c = lk[nleaf÷2]
check("At(id) selects one value", A[DGG.Cells(DD.At(c))] == Float64(nleaf ÷ 2);
    detail="$c -> position $(nleaf ÷ 2)")

k = DGG.cellposition(lk, DGG.cellat(LEAFGRID, 8.0, 46.5))
check("Contains(lon, lat) selects the cell the point is in",
    A[DGG.Cells(DD.Contains(8.0, 46.5))] == Float64(k);
    detail="(8.0, 46.5) -> position $k -> $(lk[k])")

zurich = GI.Polygon([GI.LinearRing([(8.4, 47.3), (8.7, 47.3), (8.7, 47.5),
    (8.4, 47.5), (8.4, 47.3)])])
sub = A[DGG.Cells(DGG.Covering(zurich))]
selected = Int[k for r in DGG.level_ranges(
                   DGG.query(SYS, DGG.MultiOrderCoverage(zurich); level=LEAF), LEAF)
               for p in r
               for k in (DGG.cellposition(lk, DGG.cellindex(LEAFGRID, p)),) if k !== nothing]
check("Covering(polygon) == the coverage expansion",
    parent(sub) == Float64.(selected);
    detail="$(length(sub)) of $nleaf positions")

# Any `query` target takes the same path, so a cap is a selector too.
cap = GO.UnitSpherical.SphericalCap(GO.UnitSpherical.UnitSphereFromGeographic()((8.5, 46.8)),
    0.005)
check("Covering(cap) too", !isempty(A[DGG.Cells(DGG.Covering(cap))]);
    detail="$(length(A[DGG.Cells(DGG.Covering(cap))])) cells within 0.005 rad of (8.5, 46.8)")

# A selection is a view over a NEW CellLookup whose backing is the intersected
# region, so subsetting a cube never materialises the axis either.
sublk = DD.lookup(sub, DGG.Cells)
check("a subset is a cell axis again", sublk isa DGG.CellLookup &&
                                       collect(sublk) == [lk[k] for k in selected];
    detail=sprint(show, sublk))

# --------------------------------------------------------------------------
# 4. The degenerate cases, which are the same type.
# --------------------------------------------------------------------------

whole = DGG.CellLookup(DGG.levelgrid(SYS, 5))
check("a whole level is one window",
    length(whole) == DGG.ncells(DGG.levelgrid(SYS, 5)) &&
    collect(whole)[1:4] == [DGG.cellindex(DGG.levelgrid(SYS, 5), i) for i in 1:4];
    detail=sprint(show, whole))

partial = DGG.CellLookup(DGG.PartialGrid(SYS, LEAF, materialised))
check("an explicit id list is the same axis", partial == lk && collect(partial) == materialised;
    detail=sprint(show, partial))

# --------------------------------------------------------------------------
# 5. What the compression is worth.
#
# `PartialGrid` and `CellLookup` are the same set seen from two sides, so the
# conversion is free and the regridder consumes the cube's own axis: grid
# position `i` IS `parent(A)[i]`, with no permutation anywhere.
# --------------------------------------------------------------------------

check("the lookup as a grid is aligned with the cube",
    all(DGG.cellindex(DGG.PartialGrid(lk), i) == lk[i] for i in (1, 2, nleaf ÷ 2, nleaf)))

check("memory is the windows, not the cells",
    Base.summarysize(deep) == Base.summarysize(lk);
    detail="$(Base.summarysize(lk)) bytes for $nleaf cells and for $(length(deep))")
note("the materialised level-12 vector would be " *
     "$(round(Int, 343 * nleaf * sizeof(eltype(lk)) / 1024)) KiB")

# Regridding off the cube's own axis. `DGGSpace` is how a grid names itself as a
# regridding SOURCE, and a lon/lat destination is a `RasterGrid` over two axes of
# cell centres — here a 9x5 box on the region. `Weighted(0)` divides each
# destination by the source area it actually saw, so the answer is a mean
# whatever the coverage was: this checks alignment rather than conservativity,
# which is the property the axis is responsible for.
const DST = GR.RasterGrid(DD.X(DD.Sampled(range(6.25, 10.25; length=9))),
    DD.Y(DD.Sampled(range(46.0, 47.6; length=5))))
elevation = [40.0 + 20 * sinpi(k / nleaf) for k in 1:nleaf]
zonal = DGG.regrid(elevation; to=DST, from=DGG.DGGSpace(DGG.PartialGrid(lk)),
    missingpolicy=GR.Weighted(0))
check("regridding the cube's data through its own axis",
    all(isfinite, zonal) && minimum(zonal) >= minimum(elevation) &&
    maximum(zonal) <= maximum(elevation);
    detail="$(length(zonal)) destination cells, means in " *
           "$(round(minimum(zonal); digits=2))..$(round(maximum(zonal); digits=2))")

# --------------------------------------------------------------------------
# 6. A5, which has no descendant ranges.
#
# `level_ranges` throws there, so the lookup is built by SELECTION instead:
# `descendants` names the leaves, they are resolved to positions, sorted, and
# compressed like any other position list. Same type, same laws.
# --------------------------------------------------------------------------

const A5 = DGG.A5System()
a5set = DGG.query(A5, DGG.MultiOrderCoverage(REGION); level=7)
check("A5 has no descendant ranges to expand",
    !DGG.has_sorted_subtrees(A5) &&
    (try
        DGG.level_ranges(a5set, 7)
        false
    catch e
        e isa ArgumentError
    end))

a5lk = DGG.CellLookup(a5set)
a5ids = sort!(reduce(vcat, [DGG.descendants(A5, c, 7) for c in a5set]))
check("and a cell axis all the same",
    collect(a5lk) == a5ids &&
    all(DGG.cellposition(a5lk, a5lk[k]) == k for k in eachindex(a5ids));
    detail=sprint(show, a5lk))
note("selection mode's cost is the construction, which walks the leaves;")
note("an A5 leaf can also sit outside its own ancestor's footprint, so a")
note("Covering selection over-covers by whatever the refinement does")

println()
note("what this replaces: twenty-five lines of cumsum + searchsortedfirst over")
note("level_ranges(set, leaf), written out once per page before T16")

println()
println(FAILURES[] == 0 ? "ALL CHECKS PASSED" : "$(FAILURES[]) CHECK(S) FAILED")
exit(FAILURES[] == 0 ? 0 : 1)
