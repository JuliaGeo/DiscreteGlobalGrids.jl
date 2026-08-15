# ---------------------------------------------------------------------------
# MOC storage — the two aggregation verbs, WITHOUT DimensionalData.
#
# `aggregate` reduces a `(CellVector, values)` pair to one fixed coarser level;
# `coarsen` merges complete sibling groups within a tolerance and hands back the
# mixed-level container. Contract: docs/design/moc-storage.md §2.
#
# Both are swept on one radix-4 system (HEALPix, whose four children tile their
# parent) and one radix-7 system (IGeo7, whose root cells have six children and
# whose sibling subtrees are of UNEQUAL size). The second is not decoration: a
# mean of child means equals a mean of leaves exactly when the siblings are
# equinumerous, so an implementation that summarised the children instead of the
# leaves would pass on HEALPix alone.
#
# The oracles are written against `ancestor` and a `Dict` — the route the
# implementations deliberately do not take, which walks descendant ranges and
# window lookups instead. An oracle sharing the implementation's idiom would be
# a restatement rather than a check.
#
# The laws, in order:
#
#   * AGGREGATE — the grouping is `ancestor` at the target level over the cells
#     that are PRESENT; partial groups reduce over what is there; the value
#     segments tile the data array once each, contiguously; and reducing twice
#     with an associative `f` is reducing once.
#   * COARSEN — the merge criterion (complete AND within tolerance, `missing`
#     read three ways), the value stored (`by` of the LEAF values), the error
#     bound the default `by` buys, and the cell-set recovery that completeness
#     exists for.
#   * THE CONTAINER — the wrapper is thin: cells and values still line up, and
#     the recovery law holds through `CellVector(mov; level = ...)`.
# ---------------------------------------------------------------------------

module AggregateTests

using Test
import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import GeometryOps as GO
import Statistics: mean

const FB = DGG.Fallbacks
const LONLAT = GO.UnitSpherical.GeographicFromUnitSphere()

# The Switzerland box the coverage suites use, so the fixtures read against each
# other.
const REGION = GI.Polygon([GI.LinearRing([(6.0, 45.8), (10.5, 45.8), (10.5, 47.8),
    (6.0, 47.8), (6.0, 45.8)])])

# `dense`/`coarse` are the whole-level pair, `deep`/`shallow` the coverage pair,
# and `sub` the level a rooted subtree is expanded to for `coarsen`.
const SWEEP = [
    (system=DGG.HEALPixSystem(), dense=3, coarse=1, deep=8, shallow=6, sub=3),
    (system=DGG.IGeo7System(), dense=2, coarse=1, deep=6, shallow=4, sub=3),
]

# Tolerances are INPUTS, so they are swept rather than tuned: the laws below
# hold at every one of them, and the fixture picks whichever produces a
# genuinely mixed-level answer instead of this file pinning a number.
const TOLERANCES = (0.0, 2.0, 5.0, 8.0, 10.0, 25.0, 45.0)

sysname(sys) = string(nameof(typeof(sys)))

# --- fixtures --------------------------------------------------------------

# A whole rooted subtree: complete by construction, so the completeness rule is
# exercised by REMOVING cells from it rather than by whatever a coverage
# happened to leave ragged.
subtree(sys, l) = DGG.CellVector(DGG.PartialGrid(sys, first(DGG.rootcells(sys)), l))

# The latitude of each cell's centroid, in degrees — a smooth field that is a
# fact about the grid rather than a number typed into this file.
centroid_lat(cv) = (g = DGG.levelgrid(DGG.system(cv), DGG.level(cv));
[LONLAT(DGG.cell_centroid(g, c))[2] for c in cv])

# Leaf position -> index into a data array laid out against `cv`.
dataindex(cv) = (g = DGG.levelgrid(DGG.system(cv), DGG.level(cv));
Dict(DGG.cellposition(g, c) => k for (k, c) in enumerate(cv)))

# --- oracles ---------------------------------------------------------------

# `aggregate`, spelled as the `Dict` grouping the implementation avoids.
function brute_aggregate(f, sys, ids, values, l)
    grid = DGG.levelgrid(sys, l)
    groups = Dict{Any,Vector{Int}}()
    order = Any[]
    for (k, c) in enumerate(ids)
        a = DGG.ancestor(sys, c, l)
        haskey(groups, a) || (groups[a] = Int[]; push!(order, a))
        push!(groups[a], k)
    end
    sort!(order; by=a -> DGG.cellposition(grid, a))
    return order, [f(values[groups[a]]) for a in order], groups
end

# `coarsen`, spelled bottom-up per leaf: climb from `minlevel` and take the
# first ancestor that is both complete in the present set and within tolerance.
# The criterion is monotone, so the first hit on the way up is the coarsest one.
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

# For every leaf of `cv`, the value stored for the cell covering it. Only
# meaningful where every stored cell is complete, which a whole subtree is.
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
        # A whole level aggregates to a whole level: every coarse cell has a
        # present descendant, so nothing is dropped.
        @test length(cells) == DGG.ncells(DGG.levelgrid(sys, coarse))
        @test collect(cells) == want
        @test out == wantvals
        # Every leaf value is counted exactly once — which `sum` reports and a
        # skipped or double-counted group boundary would not.
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
        # And the fixture must really be ragged, or the law above is the whole
        # level's law wearing a subset's name.
        @test count(a -> length(groups[a]) <
                         length(DGG.descendant_range(sys, a, deep)), want) > 0
    end

    @testset "the reducer sees one contiguous view per group" begin
        set = DGG.query(sys, DGG.MultiOrderCoverage(REGION); level=deep)
        cv = DGG.CellVector(set)
        values = Float64.(eachindex(cv))
        seen = Any[]
        DGG.aggregate(v -> (push!(seen, v); sum(v)), cv, values, shallow)

        # A view into the caller's array, never a copy and never a gather, and
        # the segments tile it in order. A group whose leaves straddle a window
        # gap is still ONE range here, because the data array counts only the
        # cells `cv` holds.
        @test all(v -> v isa SubArray && parent(v) === values, seen)
        spans = [only(parentindices(v)) for v in seen]
        @test all(s -> s isa AbstractUnitRange, spans)
        @test first(first(spans)) == 1
        @test last(last(spans)) == length(values)
        @test all(first(spans[i+1]) == last(spans[i]) + 1 for i in 1:length(spans)-1)
    end

    @testset "aggregating twice is aggregating once" begin
        # `sum` is associative over a partition and the level-`shallow` groups
        # refine the level-`coarse` ones, so the two routes agree exactly. A
        # segmentation that lost or duplicated a boundary would not.
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

    # Aggregating to a vector's own level is the identity; saying so is more
    # useful than a slow copy, and aggregating DOWN is a mistake either way.
    @test_throws ArgumentError DGG.aggregate(sum, cv, values, leaf)
    @test_throws ArgumentError DGG.aggregate(sum, cv, values, leaf + 1)
    @test_throws ArgumentError DGG.aggregate(sum, cv, values, -1)
    # Values that do not line up with the cells are an error, not a silently
    # truncated answer.
    @test_throws ArgumentError DGG.aggregate(sum, cv, values[1:end-1], leaf - 1)

    # A5 has no descendant ranges, so a sibling group is not a slice of
    # anything — the same refusal `level_ranges` makes.
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
    cv = subtree(sys, L)
    n = length(cv)
    root = first(DGG.rootcells(sys))
    top = first(DGG.levels(sys))
    parents = DGG.descendants(sys, root, L - 1)

    # Deterministic, no RNG: a field constant on each level-`L-1` sibling group
    # and different between them, a field of distinct values, and a banded one.
    parentgroup = [Float64(DGG.cellposition(DGG.levelgrid(sys, L - 1),
        DGG.ancestor(sys, c, L - 1))) for c in cv]
    distinct = Float64.(eachindex(cv))
    banded = [Float64((k - 1) ÷ 3) for k in 1:n]

    @testset "against the per-leaf climb" begin
        for (values, atol) in ((parentgroup, 0.0), (distinct, 0.0),
            (banded, 2.0), (fill(7.0, n), 0.0))
            cells, vals = FB._coarsen(cv, values; atol, by=mean, minlevel=top)
            want, wantvals = brute_coarsen(sys, cv, values, atol, mean, top)
            @test cells == want
            @test isequal(vals, wantvals)
        end
        # The cells come out ready to be keyed by interval: ascending, disjoint.
        cells, _ = FB._coarsen(cv, banded; atol=2.0)
        starts = [first(DGG.descendant_range(sys, c, L)) for c in cells]
        @test issorted(starts) && allunique(starts)
    end

    @testset "piecewise-constant data collapses to its pieces" begin
        cells, vals = FB._coarsen(cv, parentgroup; atol=0.0)
        # Exactly the level-`L-1` cells: each is constant, and no level-`L-2`
        # cell is, because the pieces were built to differ.
        @test cells == parents
        @test vals ≈ unique(parentgroup)
        @test length(cells) < n
    end

    @testset "atol = 0 on distinct data is the identity" begin
        cells, vals = FB._coarsen(cv, distinct; atol=0.0)
        @test cells == collect(cv)
        @test vals == distinct
        # A leaf keeps its OWN value rather than a one-element summary of it,
        # so an integer field survives as one: a leaf routed through the default
        # `by` would come back `Float64`.
        ints = collect(1:n)
        _, intvals = FB._coarsen(cv, ints; atol=0)
        @test intvals == ints
        @test eltype(intvals) === Int
    end

    @testset "the flat field collapses to one cell, and `minlevel` stops it" begin
        flat = fill(7.0, n)
        cells, vals = FB._coarsen(cv, flat; atol=0.0)
        @test cells == [root]
        @test vals == [7.0]
        # The climb stops where it is told to, and what stops there is the whole
        # level under the root.
        for stop in (L - 1, L)
            cells, vals = FB._coarsen(cv, flat; atol=0.0, minlevel=stop)
            @test cells == DGG.descendants(sys, root, stop)
            @test all(==(7.0), vals)
        end
    end

    @testset "`by` summarises the LEAVES, not the children" begin
        # Merging the whole subtree at the root: the stored value is the mean of
        # every leaf. On IGeo7 the root's children hold unequal numbers of
        # leaves, so a mean of child means is a DIFFERENT number — which is the
        # assertion that tells the two implementations apart.
        cells, vals = FB._coarsen(cv, distinct; atol=Float64(n))
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

        # `by` is honoured, and the package's default is the arithmetic mean it
        # says it is.
        _, maxvals = FB._coarsen(cv, distinct; atol=Float64(n), by=maximum)
        @test only(maxvals) == maximum(distinct)
        @test isequal(FB._coarsen(cv, banded; atol=2.0, by=mean),
            FB._coarsen(cv, banded; atol=2.0))
    end

    @testset "the error bound the default `by` buys" begin
        lat = centroid_lat(cv)
        lengths = Int[]
        for atol in TOLERANCES
            cells, vals = FB._coarsen(cv, lat; atol)
            push!(lengths, length(cells))
            stored = covering_values(sys, cells, vals, cv)
            # The pinned bound: with `by = mean` between the extremes, and the
            # extremes no more than `atol` apart, no leaf is further than `atol`
            # from what is stored for it.
            @test all(abs(lat[k] - stored[k]) <= atol for k in eachindex(lat))
        end
        # A larger tolerance never stores more, and over this sweep it stores
        # strictly fewer — so the bound above is not being met by refusing to
        # merge anything.
        @test issorted(lengths; rev=true)
        @test first(lengths) == n
        @test last(lengths) < n
    end

    @testset "completeness is what makes the cell set recoverable" begin
        # Punch holes by GROUP, not by stride: the first leaf of every third
        # sibling group is gone. A position stride reads the same but is vacuous
        # on radix 7 — five or more consecutive positions always contain a
        # multiple of five, so a `% 5` stride breaks EVERY IGeo7 group and the
        # mixed-level guard below would have nothing left to guard.
        grid = DGG.levelgrid(sys, L)
        dropped = Set(first(DGG.descendant_range(sys, a, L))
                      for (j, a) in enumerate(parents) if j % 3 == 0)
        keep = [k for (k, c) in enumerate(cv) if DGG.cellposition(grid, c) ∉ dropped]
        sub = cv[keep]
        flat = fill(2.0, length(sub))
        cells, vals = FB._coarsen(sub, flat; atol=0.0)
        want, wantvals = brute_coarsen(sys, sub, flat, 0.0, mean, top)
        @test cells == want
        @test isequal(vals, wantvals)
        @test root ∉ cells

        # The recovery law, stated where it lives: the stored subtrees name
        # exactly the positions the input held, no more and no fewer.
        named = sort!(reduce(vcat,
            [collect(DGG.descendant_range(sys, c, L)) for c in cells]))
        @test named == sort!([DGG.cellposition(grid, c) for c in sub])
        # Not vacuous in either direction: the holes left some groups whole and
        # broke others.
        @test any(c -> DGG.level(c) < L, cells)
        @test any(c -> DGG.level(c) == L, cells)
    end

    @testset "the position-list window shape answers the same" begin
        # Every fixture above stores `RangeWindows`: the compression picks that
        # shape whenever the cells run in long blocks, which a whole subtree and
        # a group-punched one both do. So the OTHER shape — a bare sorted
        # position list, and the window lookups written for it — went unread by
        # this file, and `_next_position(::PositionWindows, ...)` could search
        # from the wrong end without a single assertion noticing.
        #
        # This keeps the first sibling group whole and thins everything after it
        # to every other cell: gaps enough that the heuristic stores positions,
        # and one complete group so the MERGE path is reached on that shape too.
        g = length(DGG.descendant_range(sys, parents[1], L))
        thin = cv[[k for k in 1:n if k <= g || isodd(k)]]
        @test FB.windows(thin) isa FB.PositionWindows
        @test length(thin) < n

        # Flat on the whole group, distinct after it.
        values = [k <= g ? 1.0 : Float64(k) for k in eachindex(thin)]
        cells, vals = FB._coarsen(thin, values; atol=0.0)
        want, wantvals = brute_coarsen(sys, thin, values, 0.0, mean, top)
        @test cells == want
        @test isequal(vals, wantvals)
        # Not vacuous in either direction, which is what makes the comparison
        # above worth making: the whole group merged, the thinned ones did not,
        # and nothing outside `thin` was named.
        @test parents[1] in cells
        @test any(c -> DGG.level(c) == L, cells)
        @test all(c -> DGG.level(c) == L ? c in thin : true, cells)

        # And the fixed-level verb, which reads the same windows as intervals.
        acells, avals = DGG.aggregate(sum, thin, values, L - 1)
        awant, awantvals, _ = brute_aggregate(sum, sys, collect(thin), values, L - 1)
        @test collect(acells) == awant
        @test avals == awantvals
        @test sum(avals) == sum(values)
    end

    @testset "`missing` is read three ways" begin
        # The first level-`L-1` group all `missing`, the second mixed, the rest
        # data.
        index = dataindex(cv)
        values = Vector{Union{Float64,Missing}}(fill(4.0, n))
        for p in DGG.descendant_range(sys, parents[1], L)
            values[index[p]] = missing
        end
        values[index[first(DGG.descendant_range(sys, parents[2], L))]] = missing

        cells, vals = FB._coarsen(cv, values; atol=0.0)
        want, wantvals = brute_coarsen(sys, cv, values, 0.0, mean, top)
        @test cells == want
        @test isequal(vals, wantvals)

        # An all-`missing` region is perfectly flat, and merges to `missing`.
        @test parents[1] in cells
        @test ismissing(vals[findfirst(==(parents[1]), cells)])
        # A mixed one never merges, whatever its data half looks like — which is
        # what keeps a coastline from averaging ocean into land.
        @test parents[2] ∉ cells
        @test any(c -> DGG.level(c) == L && DGG.ancestor(sys, c, L - 1) == parents[2],
            cells)
        @test eltype(vals) === Union{Float64,Missing}
    end
end

@testset "coarsen refuses what it cannot answer" begin
    sys, L = DGG.HEALPixSystem(), 3
    cv = subtree(sys, L)
    values = Float64.(eachindex(cv))

    @test_throws ArgumentError FB._coarsen(cv, values[1:end-1]; atol=0.0)
    @test_throws ArgumentError FB._coarsen(cv, values; atol=0.0, minlevel=L + 1)
    @test_throws ArgumentError FB._coarsen(cv, values; atol=0.0, minlevel=-1)

    a5 = DGG.A5System()
    a5cv = DGG.CellVector(DGG.levelgrid(a5, 1))
    @test_throws ArgumentError DGG.coarsen(a5cv, Float64.(eachindex(a5cv)); atol=0.0)
end

# ---------------------------------------------------------------------------
# The container the core is wrapped in
# ---------------------------------------------------------------------------

@testset "coarsen hands back a MultiOrderVector: $(sysname(f.system))" for f in SWEEP
    sys, L = f.system, f.sub
    cv = subtree(sys, L)
    n = length(cv)
    lat = centroid_lat(cv)

    # The tolerance the fixture chooses rather than one this file pins: the
    # first that leaves a genuinely MIXED-level container, which is the only
    # shape that can catch a wrapper flattening one.
    atol = nothing
    for a in TOLERANCES
        cells, _ = FB._coarsen(cv, lat; atol=a)
        if length(unique(DGG.level.(cells))) > 1
            atol = a
            break
        end
    end
    @test atol !== nothing

    mov, vals = DGG.coarsen(cv, lat; atol)
    cells, corevals = FB._coarsen(cv, lat; atol)

    @test mov isa DGG.MultiOrderVector
    @test DGG.system(mov) == sys
    @test FB.reference_level(mov) == L
    @test length(vals) == length(mov)
    @test length(unique(DGG.level.(collect(mov)))) > 1
    # The wrapper is thin: the container's order IS the core's, so the values
    # still name the cells they were computed from. A constructor that permuted
    # would break exactly here.
    @test collect(mov) == cells
    @test vals == corevals

    # The recovery law through the bridge: the container expanded back to the
    # leaf level is the input, window for window.
    @test DGG.CellVector(mov; level=L) == cv

    # And the compression story end to end: every leaf is covered by exactly one
    # stored cell, and what is stored there is within `atol` of the leaf's own
    # value.
    for k in (1, n ÷ 3, n ÷ 2, n)
        j = FB.covering_position(mov, cv[k])
        @test j !== nothing
        @test abs(vals[j] - lat[k]) <= atol
    end
end

end # module AggregateTests
