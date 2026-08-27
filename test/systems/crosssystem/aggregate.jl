# ---------------------------------------------------------------------------
# The two aggregation verbs, without DimensionalData.
#
# `aggregate` reduces a `(CellVector, values)` pair to one fixed coarser level;
# `coarsen` merges complete sibling groups within a tolerance and returns the
# mixed-level container.
#
# Swept on a radix-4 system (HEALPix) and a radix-7 one (IGeo7): IGeo7's
# sibling subtrees have unequal sizes, so a mean of child means differs from a
# mean of leaves there — an implementation that summarised children instead of
# leaves would pass on HEALPix alone. The oracles group with `ancestor` and a
# `Dict`, not the implementations' descendant-range walks.
# ---------------------------------------------------------------------------

module AggregateTests

using Test
import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import GeometryOps as GO
import Statistics: mean

const EN = DGG.Engine
const LONLAT = GO.UnitSpherical.GeographicFromUnitSphere()

# The same Switzerland box the coverage suites use.
const REGION = GI.Polygon([GI.LinearRing([(6.0, 45.8), (10.5, 45.8), (10.5, 47.8),
    (6.0, 47.8), (6.0, 45.8)])])

# `dense`/`coarse` are the whole-level pair, `deep`/`shallow` the coverage pair,
# and `sub` the level a rooted subtree is expanded to for `coarsen`.
const SWEEP = [
    (system=DGG.HEALPixSystem(), dense=3, coarse=1, deep=8, shallow=6, sub=3),
    (system=DGG.IGeo7System(), dense=2, coarse=1, deep=6, shallow=4, sub=3),
]

# Swept, not tuned; fixtures pick the first that yields a mixed-level answer.
const TOLERANCES = (0.0, 2.0, 5.0, 8.0, 10.0, 25.0, 45.0)

sysname(sys) = string(nameof(typeof(sys)))

# --- fixtures --------------------------------------------------------------

# A whole rooted subtree — complete by construction; completeness is broken by
# removing cells.
rooted_subtree(sys, l) = DGG.CellVector(DGG.subtree(sys, first(DGG.rootcells(sys)), l))

# Centroid latitude in degrees: a smooth, grid-derived field.
centroid_lat(cv) = (g = DGG.levelgrid(DGG.system(cv), DGG.level(cv));
[LONLAT(DGG.cell_centroid(g, c))[2] for c in cv])

# Leaf index -> index into a data array laid out against `cv`.
dataindex(cv) = (g = DGG.levelgrid(DGG.system(cv), DGG.level(cv));
Dict(DGG.localindex(g, c) => k for (k, c) in enumerate(cv)))

# --- oracles ---------------------------------------------------------------

# Dict-grouping oracle for `aggregate`.
function brute_aggregate(f, sys, ids, values, l)
    grid = DGG.levelgrid(sys, l)
    groups = Dict{Any,Vector{Int}}()
    order = Any[]
    for (k, c) in enumerate(ids)
        a = DGG.ancestor(sys, c, l)
        haskey(groups, a) || (groups[a] = Int[]; push!(order, a))
        push!(groups[a], k)
    end
    sort!(order; by=a -> DGG.localindex(grid, a))
    return order, [f(values[groups[a]]) for a in order], groups
end

# Per-leaf oracle for `coarsen`: climb from `minlevel`, take the first
# complete-and-within-tolerance ancestor (the criterion is monotone, so the
# first hit is the coarsest).
function brute_coarsen(sys, cv, values, atol, by, minlevel)
    L = DGG.level(cv)
    index = dataindex(cv)
    leafvals(r) = [values[index[p]] for p in r]
    mergeable(vs) = (nm = count(ismissing, vs);
    nm == length(vs) || (nm == 0 && maximum(vs) - minimum(vs) <= atol))
    cells, vals = Any[], Any[]
    for (k, c) in enumerate(cv)
        stored = c
        for l in minlevel:(L-1)
            a = DGG.ancestor(sys, c, l)
            r = DGG.descendant_range(sys, a, L)
            all(p -> haskey(index, p), r) || continue
            mergeable(leafvals(r)) || continue
            stored = a
            break
        end
        isempty(cells) || cells[end] != stored || continue
        push!(cells, stored)
        if stored == c
            push!(vals, values[k])
        else
            vs = leafvals(DGG.descendant_range(sys, stored, L))
            push!(vals, count(ismissing, vs) == length(vs) ? missing : by(vs))
        end
    end
    return cells, vals
end

# Per leaf, the value stored for its covering cell. Requires every stored cell
# complete.
function covering_values(sys, cells, vals, cv)
    L = DGG.level(cv)
    index = dataindex(cv)
    out = Vector{eltype(vals)}(undef, length(cv))
    for (j, c) in enumerate(cells), p in DGG.descendant_range(sys, c, L)
        out[index[p]] = vals[j]
    end
    return out
end

# ---------------------------------------------------------------------------
# aggregate
# ---------------------------------------------------------------------------

@testset "aggregate groups by ancestor: $(sysname(f.system))" for f in SWEEP
    sys, dense, coarse = f.system, f.dense, f.coarse
    deep, shallow = f.deep, f.shallow

    @testset "a whole level" begin
        cv = DGG.CellVector(DGG.levelgrid(sys, dense))
        values = Float64.(eachindex(cv))
        cells, out = DGG.aggregate(sum, cv, values, coarse)
        want, wantvals, _ = brute_aggregate(sum, sys, collect(cv), values, coarse)

        @test cells isa DGG.CellVector
        @test DGG.level(cells) == coarse
        @test DGG.system(cells) == sys
        # A whole level aggregates to a whole level — no coarse cell dropped.
        @test length(cells) == DGG.ncells(DGG.levelgrid(sys, coarse))
        @test collect(cells) == want
        @test out == wantvals
        # `sum` totals agree only if every leaf is counted exactly once.
        @test sum(out) == sum(values)
    end

    @testset "a coverage subset, with partial groups" begin
        set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=deep)
        cv = DGG.CellVector(set)
        values = Float64.(eachindex(cv))
        cells, out = DGG.aggregate(sum, cv, values, shallow)
        want, wantvals, groups = brute_aggregate(sum, sys, collect(cv), values, shallow)

        @test collect(cells) == want
        @test out == wantvals
        @test sum(out) == sum(values)
        @test length(cells) < length(cv)
        # A coarse cell with no present descendant is ABSENT, not reduced over
        # nothing.
        @test length(cells) < DGG.ncells(DGG.levelgrid(sys, shallow))
        # The fixture must contain partial groups.
        @test count(a -> length(groups[a]) <
                         length(DGG.descendant_range(sys, a, deep)), want) > 0
    end

    @testset "the reducer sees one contiguous view per group" begin
        set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=deep)
        cv = DGG.CellVector(set)
        values = Float64.(eachindex(cv))
        seen = Any[]
        DGG.aggregate(v -> (push!(seen, v); sum(v)), cv, values, shallow)

        # Views into the caller's array, tiling it contiguously in order. A
        # group straddling a window gap is still one range, since the data
        # array indexes only the cells `cv` holds.
        @test all(v -> v isa SubArray && parent(v) === values, seen)
        spans = [only(parentindices(v)) for v in seen]
        @test all(s -> s isa AbstractUnitRange, spans)
        @test first(first(spans)) == 1
        @test last(last(spans)) == length(values)
        @test all(first(spans[i+1]) == last(spans[i]) + 1 for i in 1:length(spans)-1)
    end

    @testset "aggregating twice is aggregating once" begin
        # Associative `f` plus refining groups: two-step equals one-step.
        set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=deep)
        cv = DGG.CellVector(set)
        values = Float64.(eachindex(cv))
        mid, midvals = DGG.aggregate(sum, cv, values, shallow)
        twice, twicevals = DGG.aggregate(sum, mid, midvals, coarse)
        once, oncevals = DGG.aggregate(sum, cv, values, coarse)
        @test collect(twice) == collect(once)
        @test twicevals ≈ oncevals
    end
end

@testset "aggregate refuses what it cannot answer" begin
    sys, leaf = DGG.HEALPixSystem(), 3
    cv = DGG.CellVector(DGG.levelgrid(sys, leaf))
    values = Float64.(eachindex(cv))

    # Same-level and deeper targets throw.
    @test_throws ArgumentError DGG.aggregate(sum, cv, values, leaf)
    @test_throws ArgumentError DGG.aggregate(sum, cv, values, leaf + 1)
    @test_throws ArgumentError DGG.aggregate(sum, cv, values, -1)
    # A values length that does not match the cells throws.
    @test_throws ArgumentError DGG.aggregate(sum, cv, values[1:end-1], leaf - 1)

    # A5 has no descendant ranges, so `aggregate` refuses it.
    a5 = DGG.A5System()
    @test !DGG.has_sorted_subtrees(a5)
    a5cv = DGG.CellVector(DGG.levelgrid(a5, 1))
    @test_throws ArgumentError DGG.aggregate(sum, a5cv, Float64.(eachindex(a5cv)), 0)
end

# ---------------------------------------------------------------------------
# coarsen: the core, checked without the container
# ---------------------------------------------------------------------------

@testset "coarsen merges what the criterion allows: $(sysname(f.system))" for f in SWEEP
    sys, L = f.system, f.sub
    cv = rooted_subtree(sys, L)
    n = length(cv)
    root = first(DGG.rootcells(sys))
    top = first(DGG.levels(sys))
    parents = DGG.descendants(sys, root, L - 1)

    # Deterministic fields: constant per level-`L-1` sibling group, all
    # distinct, and banded.
    parentgroup = [Float64(DGG.localindex(DGG.levelgrid(sys, L - 1),
        DGG.ancestor(sys, c, L - 1))) for c in cv]
    distinct = Float64.(eachindex(cv))
    banded = [Float64((k - 1) ÷ 3) for k in 1:n]

    @testset "against the per-leaf climb" begin
        for (values, atol) in ((parentgroup, 0.0), (distinct, 0.0),
            (banded, 2.0), (fill(7.0, n), 0.0))
            cells, vals = EN._coarsen(cv, values; atol, by=mean, minlevel=top)
            want, wantvals = brute_coarsen(sys, cv, values, atol, mean, top)
            @test cells == want
            @test isequal(vals, wantvals)
        end
        # The cells come out ready to be keyed by interval: ascending, disjoint.
        cells, _ = EN._coarsen(cv, banded; atol=2.0)
        starts = [first(DGG.descendant_range(sys, c, L)) for c in cells]
        @test issorted(starts) && allunique(starts)
    end

    @testset "piecewise-constant data collapses to its pieces" begin
        cells, vals = EN._coarsen(cv, parentgroup; atol=0.0)
        # Exactly the level-`L-1` cells: constant within each, differing
        # between them.
        @test cells == parents
        @test vals ≈ unique(parentgroup)
        @test length(cells) < n
    end

    @testset "atol = 0 on distinct data is the identity" begin
        cells, vals = EN._coarsen(cv, distinct; atol=0.0)
        @test cells == collect(cv)
        @test vals == distinct
        # An unmerged leaf keeps its own value, so an integer field stays `Int`
        # (routing through `by = mean` would make it `Float64`).
        ints = collect(1:n)
        _, intvals = EN._coarsen(cv, ints; atol=0)
        @test intvals == ints
        @test eltype(intvals) === Int
    end

    @testset "a span wider than its own integer type never merges" begin
        # `typemax - typemin` wraps negative, which a bare `<= atol` reads as a
        # flat group.
        wide = fill(0, n)
        wide[1] = typemin(Int)
        wide[2] = typemax(Int)
        cells, vals = EN._coarsen(cv, wide; atol=0)
        @test cells[1] == cv[1] && cells[2] == cv[2]
        @test vals[1] == typemin(Int) && vals[2] == typemax(Int)
    end

    @testset "the flat field collapses to one cell, and `minlevel` stops it" begin
        flat = fill(7.0, n)
        cells, vals = EN._coarsen(cv, flat; atol=0.0)
        @test cells == [root]
        @test vals == [7.0]
        # `minlevel` stops the climb at that level.
        for stop in (L - 1, L)
            cells, vals = EN._coarsen(cv, flat; atol=0.0, minlevel=stop)
            @test cells == DGG.descendants(sys, root, stop)
            @test all(==(7.0), vals)
        end
    end

    @testset "`by` summarises the LEAVES, not the children" begin
        # The stored value is the mean of every leaf. On IGeo7 sibling subtree
        # sizes differ, so a mean of child means is a different number.
        cells, vals = EN._coarsen(cv, distinct; atol=Float64(n))
        @test cells == [root]
        @test only(vals) ≈ mean(distinct)

        sizes = [length(DGG.descendant_range(sys, k, L)) for k in DGG.children(sys, root)]
        childmeans = Float64[]
        off = 0
        for s in sizes
            push!(childmeans, mean(distinct[off+1:off+s]))
            off += s
        end
        if allequal(sizes)
            @test mean(childmeans) ≈ mean(distinct)     # HEALPix: the two agree
        else
            @test !(mean(childmeans) ≈ mean(distinct))  # IGeo7: they do not
            @test !(only(vals) ≈ mean(childmeans))
        end

        # `by` is honoured; the default is `mean`.
        _, maxvals = EN._coarsen(cv, distinct; atol=Float64(n), by=maximum)
        @test only(maxvals) == maximum(distinct)
        @test isequal(EN._coarsen(cv, banded; atol=2.0, by=mean),
            EN._coarsen(cv, banded; atol=2.0))
    end

    @testset "the error bound the default `by` buys" begin
        lat = centroid_lat(cv)
        lengths = Int[]
        for atol in TOLERANCES
            cells, vals = EN._coarsen(cv, lat; atol)
            push!(lengths, length(cells))
            stored = covering_values(sys, cells, vals, cv)
            # With `by = mean` and extremes within `atol`, every leaf is within
            # `atol` of its stored value.
            @test all(abs(lat[k] - stored[k]) <= atol for k in eachindex(lat))
        end
        # A larger tolerance never stores more; over this sweep, strictly
        # fewer — so the bound is not met by refusing to merge.
        @test issorted(lengths; rev=true)
        @test first(lengths) == n
        @test last(lengths) < n
    end

    @testset "completeness is what makes the cell set recoverable" begin
        # Holes are punched per sibling group (first leaf of every third), not
        # by index stride: a small stride hits every radix-7 group and would
        # leave the mixed-level assertions below vacuous.
        grid = DGG.levelgrid(sys, L)
        dropped = Set(first(DGG.descendant_range(sys, a, L))
                      for (j, a) in enumerate(parents) if j % 3 == 0)
        keep = [k for (k, c) in enumerate(cv) if DGG.localindex(grid, c) ∉ dropped]
        sub = cv[keep]
        flat = fill(2.0, length(sub))
        cells, vals = EN._coarsen(sub, flat; atol=0.0)
        want, wantvals = brute_coarsen(sys, sub, flat, 0.0, mean, top)
        @test cells == want
        @test isequal(vals, wantvals)
        @test root ∉ cells

        # The stored subtrees name exactly the input's indices.
        named = sort!(reduce(vcat,
            [collect(DGG.descendant_range(sys, c, L)) for c in cells]))
        @test named == sort!([DGG.localindex(grid, c) for c in sub])
        # The holes left some groups whole and broke others.
        @test any(c -> DGG.level(c) < L, cells)
        @test any(c -> DGG.level(c) == L, cells)
    end

    @testset "the index-list window shape answers the same" begin
        # Every fixture above stores `RangeWindows`; this one forces
        # `IndexWindows` — first sibling group kept whole, the rest thinned
        # to every other cell — so the merge path runs on that shape too.
        g = length(DGG.descendant_range(sys, parents[1], L))
        thin = cv[[k for k in 1:n if k <= g || isodd(k)]]
        @test EN.windows(thin) isa EN.IndexWindows
        @test length(thin) < n

        # Flat on the whole group, distinct after it.
        values = [k <= g ? 1.0 : Float64(k) for k in eachindex(thin)]
        cells, vals = EN._coarsen(thin, values; atol=0.0)
        want, wantvals = brute_coarsen(sys, thin, values, 0.0, mean, top)
        @test cells == want
        @test isequal(vals, wantvals)
        # The whole group merged, the thinned cells did not, and nothing
        # outside `thin` was named.
        @test parents[1] in cells
        @test any(c -> DGG.level(c) == L, cells)
        @test all(c -> DGG.level(c) == L ? c in thin : true, cells)

        # The fixed-level verb on the same window shape.
        acells, avals = DGG.aggregate(sum, thin, values, L - 1)
        awant, awantvals, _ = brute_aggregate(sum, sys, collect(thin), values, L - 1)
        @test collect(acells) == awant
        @test avals == awantvals
        @test sum(avals) == sum(values)
    end

    @testset "`missing` is read three ways" begin
        # First level-`L-1` group all `missing`, second mixed, rest data.
        index = dataindex(cv)
        values = Vector{Union{Float64,Missing}}(fill(4.0, n))
        for p in DGG.descendant_range(sys, parents[1], L)
            values[index[p]] = missing
        end
        values[index[first(DGG.descendant_range(sys, parents[2], L))]] = missing

        cells, vals = EN._coarsen(cv, values; atol=0.0)
        want, wantvals = brute_coarsen(sys, cv, values, 0.0, mean, top)
        @test cells == want
        @test isequal(vals, wantvals)

        # An all-`missing` group merges to `missing`.
        @test parents[1] in cells
        @test ismissing(vals[findfirst(==(parents[1]), cells)])
        # A mixed group never merges (ocean is never averaged into land).
        @test parents[2] ∉ cells
        @test any(c -> DGG.level(c) == L && DGG.ancestor(sys, c, L - 1) == parents[2],
            cells)
        @test eltype(vals) === Union{Float64,Missing}
    end
end

@testset "coarsen refuses what it cannot answer" begin
    sys, L = DGG.HEALPixSystem(), 3
    cv = rooted_subtree(sys, L)
    values = Float64.(eachindex(cv))

    @test_throws ArgumentError EN._coarsen(cv, values[1:end-1]; atol=0.0)
    @test_throws ArgumentError EN._coarsen(cv, values; atol=0.0, minlevel=L + 1)
    @test_throws ArgumentError EN._coarsen(cv, values; atol=0.0, minlevel=-1)

    a5 = DGG.A5System()
    a5cv = DGG.CellVector(DGG.levelgrid(a5, 1))
    @test_throws ArgumentError DGG.coarsen(a5cv, Float64.(eachindex(a5cv)); atol=0.0)
end

# ---------------------------------------------------------------------------
# The container the core is wrapped in
# ---------------------------------------------------------------------------

@testset "coarsen hands back a MultiOrderVector: $(sysname(f.system))" for f in SWEEP
    sys, L = f.system, f.sub
    cv = rooted_subtree(sys, L)
    n = length(cv)
    lat = centroid_lat(cv)

    # First tolerance that yields a genuinely mixed-level container.
    atol = nothing
    for a in TOLERANCES
        cells, _ = EN._coarsen(cv, lat; atol=a)
        if length(unique(DGG.level.(cells))) > 1
            atol = a
            break
        end
    end
    @test atol !== nothing

    mov, vals = DGG.coarsen(cv, lat; atol)
    cells, corevals = EN._coarsen(cv, lat; atol)

    @test mov isa DGG.MultiOrderVector
    @test DGG.system(mov) == sys
    @test EN.reference_level(mov) == L
    @test length(vals) == length(mov)
    @test length(unique(DGG.level.(collect(mov)))) > 1
    # Container order equals the core's, so values still name their cells.
    @test collect(mov) == cells
    @test vals == corevals

    # Expanded back to the leaf level, the container is the input.
    @test DGG.CellVector(mov; level=L) == cv

    # Every leaf is covered by one stored cell whose value is within `atol`.
    for k in (1, n ÷ 3, n ÷ 2, n)
        j = EN.covering_index(mov, cv[k])
        @test j !== nothing
        @test abs(vals[j] - lat[k]) <= atol
    end
end

end # module AggregateTests
