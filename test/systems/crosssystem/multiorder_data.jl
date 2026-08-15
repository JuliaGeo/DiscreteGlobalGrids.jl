# ---------------------------------------------------------------------------
# MOC storage — the DimensionalData face of the mixed-level container.
#
# `MultiOrderLookup` claims to BE the container's id vector while answering a
# cube's questions about it, and `expand` claims to present every leaf while
# storing only the mesh. Everything below is one of those two claims:
#
#   * TWO KINDS OF MEMBERSHIP, AS SELECTORS — `At` is EXACT and `Contains` is
#     COVERING. A leaf strictly under a stored coarse cell must therefore be a
#     `SelectorError` for the first and that cell's value for the second, and a
#     cell DEEPER than the reference level must answer the second too, through
#     its reference-level ancestor. Confusing the two is the mistake the
#     container exists to make impossible; here it is the mistake a selector
#     table can quietly reintroduce.
#   * THE SUBSET FORK — an ascending subset is a sorted disjoint interval list
#     again and stays a `MultiOrderLookup`; a reordered one is not, and takes
#     the `Categorical` fallback `CellLookup` takes. Stated through `getindex`,
#     through `vcat` (which arrives by `rebuild`, not by `getindex`) and through
#     a DISORDERED `vcat`, which is the one that tells a real check from a
#     rubber stamp.
#   * EXPAND IS LAZY — the deliverable. Presenting the same mesh one level
#     deeper multiplies the elements and must not move `Base.summarysize` at
#     all, because what is stored is the values and one leaf count per stored
#     cell. Materialising it fills runs, and must agree element for element with
#     the elementwise reads it replaces.
#   * THE ROUND-TRIP BOUND — `expand(coarsen(A; atol), L)` is within `atol` of
#     `A` everywhere, and EXACTLY `A` where the data was piecewise constant.
#   * THE VERBS' RESTRICTION — one dimension over `Cells` in v1, and the right
#     lookup for the verb. Each of the three refuses the other's axis.
#
# Swept on one radix-4 system (HEALPix) and one radix-7 system (IGeo7), for the
# reason `aggregate.jl` gives: sibling subtrees of unequal size are what tell a
# per-cell leaf count apart from a constant one. A5 has no container at all —
# `multiorder_vector.jl` owns that exclusion.
#
# No tolerance is pinned here. The fixture sweeps candidates and keeps the first
# that produces a genuinely MIXED-level axis, because an axis that happened to
# be flat would make every law above vacuous — and it asserts that it found one.
# ---------------------------------------------------------------------------

module MultiOrderDataTests

using Test
import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GeometryOps as GO
import SmallCollections

const FB = DGG.Fallbacks
const CL = DGG.CellLookups
const LONLAT = GO.UnitSpherical.GeographicFromUnitSphere()

const SWEEP = [
    (system=DGG.HEALPixSystem(), leaf=3),
    (system=DGG.IGeo7System(), leaf=3),
]

# Inputs, not tuning: the fixture keeps the first that leaves a mixed-level
# container, so a change to the merge criterion moves which one is chosen rather
# than breaking a number typed into this file.
const TOLERANCES = (2.0, 5.0, 8.0, 12.0, 20.0, 45.0)

sysname(sys) = string(nameof(typeof(sys)))

# A whole rooted subtree: complete by construction, so every stored cell of the
# coarsened container covers a full sibling family and the covering law below
# has something to say about every leaf.
subtree(sys, l) = DGG.CellVector(DGG.PartialGrid(sys, first(DGG.rootcells(sys)), l))

# The latitude of each cell's centroid, in degrees — a smooth field that is a
# fact about the grid rather than a number typed into this file.
centroid_lat(cv) = (g = DGG.levelgrid(DGG.system(cv), DGG.level(cv));
[LONLAT(DGG.cell_centroid(g, c))[2] for c in cv])

# The whole fixture: a leaf axis, its data, and the coarsest tolerance-driven
# mesh over it that is genuinely mixed-level.
function fixture(sys, L)
    cv = subtree(sys, L)
    lat = centroid_lat(cv)
    A = DD.DimArray(lat, DGG.Cells(DGG.CellLookup(cv)); name=:lat)
    for atol in TOLERANCES
        M = DGG.coarsen(A; atol)
        mov = parent(DD.lookup(M, DGG.Cells))
        length(unique(DGG.level, collect(mov))) > 1 &&
            return (; cv, lat, A, M, mov, atol)
    end
    return nothing
end

# `Covering`'s oracle, spelled leaf by leaf: every leaf the coverage names,
# resolved through the compression verb. The implementation compares INTERVALS
# instead, so this is a different sentence about the same set.
function covering_byhand(mov, sys, L, target)
    grid = DGG.levelgrid(sys, L)
    set = DGG.query(sys, DGG.MultiOrderCoverage(target); level=L)
    out = Int[]
    for r in DGG.level_ranges(set, L), p in r
        k = FB.covering_position(mov, DGG.cellindex(grid, p))
        k === nothing || push!(out, k)
    end
    return unique!(sort!(out))
end

# ---------------------------------------------------------------------------
# The axis, and the two questions it answers
# ---------------------------------------------------------------------------

@testset "a mixed-level cell axis: $(sysname(f.system))" for f in SWEEP
    sys, L = f.system, f.leaf
    fx = fixture(sys, L)
    # The premise of everything below, asserted rather than assumed.
    @test fx !== nothing
    cv, lat, A, M, mov, atol = fx.cv, fx.lat, fx.A, fx.M, fx.mov, fx.atol
    lk = DD.lookup(M, DGG.Cells)
    vals = parent(M)
    ids = collect(mov)

    @testset "what it is" begin
        @test lk isa DD.Lookups.Lookup
        @test lk isa DGG.MultiOrderLookup
        # `parent` is the VALUES — the container — for the reason the
        # `CellLookup` file spells out at its own `parent`.
        @test parent(lk) === mov
        @test collect(lk) == ids
        @test length(lk) == length(mov)
        @test eltype(lk) === DGG.cellindextype(sys)
        @test DGG.system(lk) == sys
        @test FB.reference_level(lk) == L
        @test length(unique(DGG.level, ids)) > 1
        # "Forward" in the container's own key — the reference-level interval
        # STARTS — which is what every verb here searches in. Stated as the
        # starts rather than as `issorted(ids)`, because `isless` on an id
        # compares its level first and a mixed-level container is unsorted in
        # THAT sense on most systems.
        @test DD.Lookups.order(lk) === DD.Lookups.ForwardOrdered()
        starts = [first(DGG.descendant_range(sys, c, L)) for c in ids]
        @test issorted(starts) && allunique(starts)
        @test DD.Lookups.val(lk) === parent(lk)
        @test DD.Lookups.metadata(lk) === DD.Lookups.NoMetadata()
        @test DD.Lookups.bounds(lk) == (first(ids), last(ids))
        @test DD.name(M) === :lat
        @test length(M) == length(mov)
        # It really is a compression of the leaf axis, or the mesh is decoration.
        @test length(M) < length(A)
    end

    # The two verbs, on the axis itself, before any selector plumbing can
    # confuse them.
    @testset "At is exact and Contains is covering" begin
        @test all(DGG.cellposition(lk, ids[k]) == k for k in eachindex(ids))
        @test all(FB.covering_position(lk, ids[k]) == k for k in eachindex(ids))

        j = findfirst(c -> DGG.level(c) < L, ids)
        @test j !== nothing                     # a mixed axis has a coarse cell
        coarse = ids[j]
        leaf = first(DGG.descendants(sys, coarse, L))
        @test leaf != coarse
        @test DGG.cellposition(lk, leaf) === nothing
        @test FB.covering_position(lk, leaf) == j
        # Deeper than the reference level: keyed through its reference-level
        # ancestor, which is the deliberate half of `covering_position`.
        deeper = first(DGG.descendants(sys, coarse, L + 1))
        @test DGG.cellposition(lk, deeper) === nothing
        @test FB.covering_position(lk, deeper) == j
        # And an ANCESTOR of a stored cell is held by neither: a container is
        # not its own coarsening.
        up = DGG.ancestor(sys, coarse, DGG.level(coarse) - 1)
        @test DGG.cellposition(lk, up) === nothing
        @test FB.covering_position(lk, up) === nothing
    end

    @testset "the selectors are those two verbs" begin
        for k in (1, length(ids) ÷ 2, length(ids))
            @test M[DGG.Cells(DD.At(ids[k]))] == vals[k]
            @test M[DGG.Cells(DD.Contains(ids[k]))] == vals[k]
        end

        j = findfirst(c -> DGG.level(c) < L, ids)
        coarse = ids[j]
        leaf = first(DGG.descendants(sys, coarse, L))
        deeper = first(DGG.descendants(sys, coarse, L + 1))
        # The pairing, as selectors: `At` refuses a cell the axis does not
        # store, `Contains` resolves it to the cell that speaks for it.
        @test_throws DD.Lookups.SelectorError M[DGG.Cells(DD.At(leaf))]
        @test M[DGG.Cells(DD.Contains(leaf))] == vals[j]
        @test M[DGG.Cells(DD.Contains(deeper))] == vals[j]
        @test DD.Lookups.hasselection(lk, DD.At(ids[j]))
        @test !DD.Lookups.hasselection(lk, DD.At(leaf))
        @test DD.Lookups.hasselection(lk, DD.Contains(leaf))

        # A point falls INSIDE a cell rather than naming one, so both point
        # selectors are the covering question — read at whatever level the cell
        # it lands in happens to sit.
        for k in (1, j, length(ids))
            grid = DGG.levelgrid(sys, DGG.level(ids[k]))
            lon, y = LONLAT(DGG.cell_centroid(grid, ids[k]))
            @test M[DGG.Cells(DD.Contains(lon, y))] == vals[k]
            @test M[DGG.Cells(DD.At(lon, y))] == vals[k]
        end

        # A point in a root subtree the axis does not hold names nothing.
        other = DGG.rootcells(sys)[end]
        lon, y = LONLAT(DGG.cell_centroid(DGG.levelgrid(sys, DGG.level(other)), other))
        @test DGG.cellat(mov, lon, y) === nothing
        @test_throws DD.Lookups.SelectorError M[DGG.Cells(DD.Contains(lon, y))]
    end

    @testset "Covering keeps a stored cell whole" begin
        cap = GO.UnitSpherical.SphericalCap(
            DGG.cell_centroid(DGG.levelgrid(sys, L), cv[1]), 0.15)
        byhand = covering_byhand(mov, sys, L, cap)
        # Non-vacuous in both directions: a proper, non-empty subset.
        @test !isempty(byhand)
        @test length(byhand) < length(mov)

        sub = M[DGG.Cells(DGG.Covering(cap))]
        @test parent(sub) == vals[byhand]
        sublk = DD.lookup(sub, DGG.Cells)
        @test sublk isa DGG.MultiOrderLookup
        @test collect(sublk) == ids[byhand]
        @test FB.reference_level(sublk) == L
        @test all(DGG.cellposition(sublk, sublk[k]) == k for k in eachindex(byhand))

        # A region no cell of the axis meets selects nothing at all.
        away = GO.UnitSpherical.SphericalCap(
            -DGG.cell_centroid(DGG.levelgrid(sys, L), cv[1]), 0.05)
        @test isempty(M[DGG.Cells(DGG.Covering(away))])
    end

    # `vcat` arrives at the lookup through `rebuild`, not through `getindex`,
    # and two ascending disjoint halves are one container again. Swept rather
    # than left to the mechanics testset below because DimensionalData decides
    # whether to even try by comparing the axes' VALUES, and mixed-level ids
    # compare differently per system.
    @testset "vcat of two disjoint ascending halves" begin
        n = length(M)
        joined = vcat(M[1:(n÷2)], M[(n÷2+1):n])
        jlk = DD.lookup(joined, DGG.Cells)
        @test jlk isa DGG.MultiOrderLookup
        @test jlk == lk
        @test collect(jlk) == ids
        @test parent(joined) == vals
    end
end

# ---------------------------------------------------------------------------
# expand: the compression, presented at one level
# ---------------------------------------------------------------------------

@testset "expand presents the leaves and stores the mesh: $(sysname(f.system))" for
    f in SWEEP

    sys, L = f.system, f.leaf
    fx = fixture(sys, L)
    @test fx !== nothing
    cv, lat, A, M, mov, atol = fx.cv, fx.lat, fx.A, fx.M, fx.mov, fx.atol
    vals = parent(M)

    E = DGG.expand(M, L)
    data = parent(E)

    @testset "the axis is the leaf axis it came from" begin
        elk = DD.lookup(E, DGG.Cells)
        @test elk isa DGG.CellLookup
        @test DGG.level(elk) == L
        @test parent(elk) == cv
        @test length(E) == length(A)
        @test DD.name(E) === :lat
    end

    @testset "every leaf reads the cell covering it" begin
        want = [vals[FB.covering_position(mov, c)] for c in cv]
        @test data == want
        @test all(data[k] == want[k] for k in eachindex(want))
        # The bound `coarsen` promises, read end to end through the cube.
        @test all(abs(data[k] - lat[k]) <= atol for k in eachindex(lat))
        @test E[DGG.Cells(DD.At(cv[3]))] == want[3]
    end

    @testset "materialising fills runs, and agrees with reading" begin
        @test collect(data) == [data[k] for k in eachindex(data)]
        dest = Vector{eltype(data)}(undef, length(data))
        @test copyto!(dest, data) == collect(data)
        @test collect(data) isa Vector{eltype(data)}
    end

    # The deliverable. One mesh, two presentation levels: the deeper names many
    # more cells and stores the same values and the same one leaf count per
    # stored cell, so nothing moves.
    @testset "memory is O(#stored cells), not O(#leaves)" begin
        deep = DGG.expand(M, L + 1)
        @test length(deep) > length(E)
        @test Base.summarysize(parent(deep)) == Base.summarysize(data)
        # Against the thing it replaces: one word per presented element is the
        # floor for a materialised vector, and this is under it.
        @test Base.summarysize(parent(deep)) < 8 * length(deep)
        # The values it presents are still the covering cell's, one level down.
        deepcv = parent(DD.lookup(deep, DGG.Cells))
        @test all(parent(deep)[k] == vals[FB.covering_position(mov, deepcv[k])]
                  for k in (1, length(deep) ÷ 2, length(deep)))
    end

    @testset "where the data was flat, the round trip is exact" begin
        # Constant on each level-`L-1` sibling group and different between them,
        # so `atol = 0` merges exactly those groups and nothing above them.
        coarse = DGG.levelgrid(sys, L - 1)
        piece = [Float64(DGG.cellposition(coarse, DGG.ancestor(sys, c, L - 1)))
                 for c in cv]
        P = DD.DimArray(piece, DGG.Cells(DGG.CellLookup(cv)))
        C = DGG.coarsen(P; atol=0.0)
        @test length(C) < length(P)
        @test collect(parent(DGG.expand(C, L))) == piece
    end

    @testset "expand refuses a level it cannot name" begin
        # The container holds level-`L` cells, which have no interval above `L`.
        @test_throws ArgumentError DGG.expand(M, L - 1)
    end
end

# ---------------------------------------------------------------------------
# The three cube verbs, against the cores they wrap
# ---------------------------------------------------------------------------

@testset "the DimArray verbs wrap the cores: $(sysname(f.system))" for f in SWEEP
    sys, L = f.system, f.leaf
    fx = fixture(sys, L)
    @test fx !== nothing
    cv, lat, A, atol = fx.cv, fx.lat, fx.A, fx.atol

    @testset "aggregate" begin
        Ag = DGG.aggregate(sum, A, L - 1)
        cells, want = DGG.aggregate(sum, cv, lat, L - 1)
        alk = DD.lookup(Ag, DGG.Cells)
        @test alk isa DGG.CellLookup
        @test DGG.level(alk) == L - 1
        @test parent(alk) == cells
        @test parent(Ag) == want
        @test DD.name(Ag) === :lat
    end

    @testset "coarsen" begin
        M = DGG.coarsen(A; atol)
        mov, want = DGG.coarsen(cv, lat; atol)
        @test parent(DD.lookup(M, DGG.Cells)) == mov
        @test parent(M) == want
        # The keywords reach the core rather than being dropped on the way:
        # stopping the climb at the leaf level is the identity, and a `by` that
        # is not the mean says so in the values.
        flat = DGG.coarsen(A; atol, minlevel=L)
        @test length(flat) == length(A)
        @test parent(flat) == lat
        @test parent(DGG.coarsen(A; atol, by=maximum)) !=
              parent(DGG.coarsen(A; atol, by=minimum))
    end
end

@testset "the verbs refuse the axis they cannot read" begin
    sys, L = DGG.HEALPixSystem(), 3
    fx = fixture(sys, L)
    @test fx !== nothing
    A, M = fx.A, fx.M

    # Each verb has ONE axis shape it reads, and says which when handed the
    # other: a mesh cannot be aggregated to a fixed level, a leaf axis cannot be
    # expanded.
    @test_throws ArgumentError DGG.aggregate(sum, M, L - 1)
    @test_throws ArgumentError DGG.coarsen(M; atol=1.0)
    @test_throws ArgumentError DGG.expand(A, L)

    # One dimension in v1, and a `Cells` dimension at that.
    lat = parent(A)
    two = DD.DimArray(hcat(lat, lat), (DGG.Cells(DD.lookup(A, DGG.Cells)), DD.X(1:2)))
    @test_throws ArgumentError DGG.aggregate(sum, two, L - 1)
    @test_throws ArgumentError DGG.coarsen(two; atol=1.0)
    plain = DD.DimArray(lat, DD.X(eachindex(lat)))
    @test_throws ArgumentError DGG.aggregate(sum, plain, L - 1)
    @test_throws ArgumentError DGG.expand(plain, L)
end

# ---------------------------------------------------------------------------
# The DimensionalData plumbing, on one system
# ---------------------------------------------------------------------------

@testset "lookup mechanics" begin
    sys, L = DGG.IGeo7System(), 3
    fx = fixture(sys, L)
    @test fx !== nothing
    M, mov = fx.M, fx.mov
    lk = DD.lookup(M, DGG.Cells)
    ids = collect(lk)
    vals = parent(M)

    @testset "the subset fork" begin
        @test lk[:] === lk
        @test lk[[1, 3, 5]] isa DGG.MultiOrderLookup
        @test collect(lk[[1, 3, 5]]) == ids[[1, 3, 5]]
        # A reordered subset is not a sorted disjoint interval list, and says so
        # by wearing the ordinary lookup instead of lying about it.
        @test !(lk[[3, 1]] isa DGG.MultiOrderLookup)
        @test lk[[3, 1]] isa DD.Lookups.Categorical
        @test collect(lk[[3, 1]]) == ids[[3, 1]]
        mask = falses(length(lk))
        mask[2] = mask[4] = true
        @test collect(lk[mask]) == ids[[2, 4]]
        @test_throws BoundsError lk[falses(length(lk) + 1)]
        # A neighbour list is a `SmallVector`; indexing an axis by one is
        # ambiguous unless the tie against SmallCollections' own method is
        # broken, exactly as it is for `CellLookup`.
        @test collect(lk[SmallCollections.SmallVector{8,Int}([3, 4, 5])]) == ids[3:5]
        @test collect(reverse(lk)) == reverse(ids)
        @test reverse(lk) isa DD.Lookups.Lookup
    end

    @testset "rebuild, which is how the cube concatenates" begin
        @test DD.Lookups.rebuild(lk) === lk
        @test DD.Lookups.rebuild(lk; data=parent(lk)) === lk
        @test DD.Lookups.rebuild(lk; data=ids) isa DGG.MultiOrderLookup
        @test DD.Lookups.rebuild(lk; data=ids) == lk
        @test DD.Lookups.rebuild(lk; data=mov) isa DGG.MultiOrderLookup
        @test DD.Lookups.rebuild(lk; data=reverse(ids)) isa DD.Lookups.Categorical
        @test_throws ArgumentError DD.Lookups.rebuild(lk; data=[1, 2, 3])
        @test occursin("cat", sprint(showerror,
            try
                DD.Lookups.rebuild(lk; data=[1, 2, 3])
            catch err
                err
            end))
        # A cell deeper than the reference level has no interval to be keyed by,
        # and that is an error rather than a silently deeper container.
        deep = first(DGG.descendants(sys, ids[1], L + 1))
        @test_throws ArgumentError DD.Lookups.rebuild(lk; data=[deep])
    end

    @testset "concatenation and reduction" begin
        n = length(M)
        lo, hi = M[1:(n÷2)], M[(n÷2+1):n]
        for joined in (vcat(lo, hi), cat(lo, hi; dims=DGG.Cells))
            jlk = DD.lookup(joined, DGG.Cells)
            @test jlk isa DGG.MultiOrderLookup
            @test jlk == lk
            @test collect(jlk) == ids
            @test parent(joined) == vals
        end
        # Disjoint with a GAP between them, which is the case a rebuild that
        # only handled adjacency would get wrong: the two halves no longer tile
        # anything, and the join is still a container.
        gapped = vcat(M[1:(n÷2)], M[(n÷2+2):n])
        glk = DD.lookup(gapped, DGG.Cells)
        @test glk isa DGG.MultiOrderLookup
        @test collect(glk) == ids[[1:(n÷2); (n÷2+2):n]]
        @test glk != lk

        # A disordered join is not an ascending interval list, and `rebuild` is
        # where that is decided — `vcat` never gets there, because
        # DimensionalData refuses misaligned ordered lookups one layer above.
        @test DD.Lookups.rebuild(lk; data=vcat(ids[(n÷2+1):n], ids[1:(n÷2)])) isa
              DD.Lookups.Categorical

        r = sum(M; dims=DGG.Cells)
        @test size(r) == (1,)
        @test r[1] == sum(vals)
        @test DD.lookup(r, DGG.Cells) isa DD.Lookups.NoLookup
        @test DD.Lookups.reducelookup(lk) isa DD.Lookups.NoLookup
    end

    @testset "equality, bounds and show" begin
        same = DGG.MultiOrderLookup(DGG.MultiOrderVector(sys, ids; reference_level=L))
        @test lk == same
        @test lk != DGG.MultiOrderLookup(DGG.MultiOrderVector(sys, ids[1:end-1];
            reference_level=L))
        @test DGG.MultiOrderLookup(lk) === lk

        @test DD.Lookups.bounds(lk[Int[]]) == (nothing, nothing)

        s = sprint(show, lk)
        @test occursin("MultiOrderLookup", s)
        @test occursin("IGeo7System", s)
        @test occursin(string(length(ids)), s)
        @test occursin("MultiOrderLookup",
            sprint(show, MIME"text/plain"(), DD.dims(M, DGG.Cells)))

        @test_throws DimensionMismatch DD.DimArray(zeros(3), DGG.Cells(lk))
    end

    # The other way in: a coverage read as an axis directly, with no data verb
    # in between. A multi-order query is already the shape this lookup wants.
    @testset "a coverage is an axis in its own right" begin
        set = DGG.query(sys, DGG.MultiOrderCoverage(
                GO.UnitSpherical.SphericalCap(FB.unit_point(8.0, 46.5), 0.05));
            level=6)
        setlk = DGG.MultiOrderLookup(set)
        @test setlk isa DGG.MultiOrderLookup
        @test collect(setlk) == collect(set)
        @test length(unique(DGG.level, collect(setlk))) > 1
    end
end

end # module MultiOrderDataTests
