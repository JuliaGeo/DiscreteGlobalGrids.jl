# ---------------------------------------------------------------------------
# The DimensionalData face of the mixed-level container: `MultiOrderLookup`
# and the cube methods of `aggregate`, `coarsen` and `expand`.
#
#   * selectors: `At` is exact membership, `Contains` is covering — including
#     cells deeper than the reference level, via their ancestor there.
#   * subsets: an ascending subset stays a `MultiOrderLookup`; a reordered one
#     falls back to `Categorical`. `vcat` is swept over every split — an
#     ordered-lookup claim would let DimensionalData's cat pre-check silently
#     drop the axis at level-inverting splits.
#   * `expand` is lazy: a deeper presentation level grows the element count
#     but not `Base.summarysize`. `expand(coarsen(A; atol), L)` is within
#     `atol` of `A`, and exact where the data was piecewise constant.
#   * the verbs are 1-D over `Cells`, and each refuses the other's axis shape.
#
# Swept on HEALPix (radix 4) and IGeo7 (radix 7); see aggregate.jl. Tolerances
# are swept and the first mixed-level fixture kept, asserted to exist.
# ---------------------------------------------------------------------------

module MultiOrderDataTests

using Test
import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GeometryOps as GO
import SmallCollections

const FB = DGG.Fallbacks
const EN = DGG.Engine
const CL = DGG.CellLookups
const LONLAT = GO.UnitSpherical.GeographicFromUnitSphere()

const SWEEP = [
    (system=DGG.HEALPixSystem(), leaf=3),
    (system=DGG.IGeo7System(), leaf=3),
]

# The fixture keeps the first tolerance that yields a mixed-level container.
const TOLERANCES = (2.0, 5.0, 8.0, 12.0, 20.0, 45.0)

sysname(sys) = string(nameof(typeof(sys)))

# A whole rooted subtree — complete, so every stored cell covers a full
# sibling family.
rooted_subtree(sys, l) = DGG.CellVector(DGG.subtree(sys, first(DGG.rootcells(sys)), l))

# Centroid latitude in degrees: a smooth, grid-derived field.
centroid_lat(cv) = (g = DGG.levelgrid(DGG.system(cv), DGG.level(cv));
[LONLAT(DGG.cell_centroid(g, c))[2] for c in cv])

# Leaf axis, its data, and the first tolerance whose mesh is mixed-level.
function fixture(sys, L)
    cv = rooted_subtree(sys, L)
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

# Leaf-by-leaf oracle for `Covering`; the implementation compares intervals
# instead.
function covering_byhand(mov, sys, L, target)
    grid = DGG.levelgrid(sys, L)
    set = DGG.query(sys, DGG.MultiOrderCoverage(target); level=L)
    out = Int[]
    for r in DGG.level_ranges(set, L), p in r
        k = EN.covering_position(mov, DGG.cellindex(grid, p))
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
    @test fx !== nothing
    cv, lat, A, M, mov, atol = fx.cv, fx.lat, fx.A, fx.M, fx.mov, fx.atol
    lk = DD.lookup(M, DGG.Cells)
    vals = parent(M)
    ids = collect(mov)

    @testset "what it is" begin
        @test lk isa DD.Lookups.Lookup
        @test lk isa DGG.MultiOrderLookup
        # `parent` is the container, as for `CellLookup`.
        @test parent(lk) === mov
        @test collect(lk) == ids
        @test length(lk) == length(mov)
        @test eltype(lk) === DGG.cellindextype(sys)
        @test DGG.system(lk) == sys
        @test EN.reference_level(lk) == L
        @test length(unique(DGG.level, ids)) > 1
        # `Unordered`: ids compare level-first under `isless`, so a mixed-level
        # axis is unsorted in the sense DimensionalData reads — an ordered
        # claim makes cat silently drop the axis. The container's own
        # interval-start order is asserted beside it.
        @test DD.Lookups.order(lk) === DD.Lookups.Unordered()
        starts = [first(DGG.descendant_range(sys, c, L)) for c in ids]
        @test issorted(starts) && allunique(starts)
        @test DD.Lookups.val(lk) === parent(lk)
        @test DD.Lookups.metadata(lk) === DD.Lookups.NoMetadata()
        # The generic unordered answer; a `(first, last)` pair here would not
        # be an interval.
        @test DD.Lookups.bounds(lk) === (nothing, nothing)
        @test DD.name(M) === :lat
        @test length(M) == length(mov)
        # Really a compression: fewer stored cells than leaves.
        @test length(M) < length(A)
    end

    @testset "At is exact and Contains is covering" begin
        @test all(DGG.cellposition(lk, ids[k]) == k for k in eachindex(ids))
        @test all(EN.covering_position(lk, ids[k]) == k for k in eachindex(ids))

        j = findfirst(c -> DGG.level(c) < L, ids)
        @test j !== nothing                     # a mixed axis has a coarse cell
        coarse = ids[j]
        leaf = first(DGG.descendants(sys, coarse, L))
        @test leaf != coarse
        @test DGG.cellposition(lk, leaf) === nothing
        @test EN.covering_position(lk, leaf) == j
        # Deeper than the reference level: resolved through its ancestor there.
        deeper = first(DGG.descendants(sys, coarse, L + 1))
        @test DGG.cellposition(lk, deeper) === nothing
        @test EN.covering_position(lk, deeper) == j
        # An ancestor of a stored cell is neither stored nor covered.
        up = DGG.ancestor(sys, coarse, DGG.level(coarse) - 1)
        @test DGG.cellposition(lk, up) === nothing
        @test EN.covering_position(lk, up) === nothing
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
        # `At` refuses an unstored cell; `Contains` resolves it to its
        # covering cell.
        @test_throws DD.Lookups.SelectorError M[DGG.Cells(DD.At(leaf))]
        @test M[DGG.Cells(DD.Contains(leaf))] == vals[j]
        @test M[DGG.Cells(DD.Contains(deeper))] == vals[j]
        @test DD.Lookups.hasselection(lk, DD.At(ids[j]))
        @test !DD.Lookups.hasselection(lk, DD.At(leaf))
        @test DD.Lookups.hasselection(lk, DD.Contains(leaf))

        # Both point selectors ask the covering question.
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
        @test EN.reference_level(sublk) == L
        @test all(DGG.cellposition(sublk, sublk[k]) == k for k in eachindex(byhand))

        # A region no cell of the axis meets selects nothing at all.
        away = GO.UnitSpherical.SphericalCap(
            -DGG.cell_centroid(DGG.levelgrid(sys, L), cv[1]), 0.05)
        @test isempty(M[DGG.Cells(DGG.Covering(away))])
    end

    # Swept over every split: DimensionalData's cat pre-check reads the axis
    # values, and an ordered-lookup claim would silently drop the axis exactly
    # at splits where a coarse id follows a deep one.
    @testset "vcat of two disjoint ascending halves, at every split" begin
        n = length(M)
        for s in 1:(n-1)
            joined = vcat(M[1:s], M[(s+1):n])
            @test joined isa DD.AbstractDimArray
            jlk = DD.lookup(joined, DGG.Cells)
            @test jlk isa DGG.MultiOrderLookup
            if !(jlk isa DGG.MultiOrderLookup)
                break   # one failing split is diagnosis enough
            end
            @test jlk == lk && parent(joined) == vals
        end
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
        want = [vals[EN.covering_position(mov, c)] for c in cv]
        @test data == want
        @test all(data[k] == want[k] for k in eachindex(want))
        # The bound `coarsen` promises, read end to end through the cube.
        @test all(abs(data[k] - lat[k]) <= atol for k in eachindex(lat))
        @test E[DGG.Cells(DD.At(cv[3]))] == want[3]
    end

    @testset "the pair form is the DimArray form without the wrapper" begin
        pcv, pdata = DGG.expand(mov, vals, L)
        @test pcv == cv
        @test collect(pdata) == collect(data)
        @test_throws ArgumentError DGG.expand(mov, vals[1:(end-1)], L)
        # And the container's own leaf-level expansion, as `expand` reads it
        # on the other region types.
        @test DGG.expand(mov, L) == cv
    end

    @testset "materialising fills runs, and agrees with reading" begin
        @test collect(data) == [data[k] for k in eachindex(data)]
        dest = Vector{eltype(data)}(undef, length(data))
        @test copyto!(dest, data) == collect(data)
        @test collect(data) isa Vector{eltype(data)}
    end

    # One mesh at two presentation levels: the deeper names more cells and
    # stores the same bytes.
    @testset "memory is O(#stored cells), not O(#leaves)" begin
        deep = DGG.expand(M, L + 1)
        @test length(deep) > length(E)
        @test Base.summarysize(parent(deep)) == Base.summarysize(data)
        # Under the one-word-per-element floor of a materialised vector.
        @test Base.summarysize(parent(deep)) < 8 * length(deep)
        # The values it presents are still the covering cell's, one level down.
        deepcv = parent(DD.lookup(deep, DGG.Cells))
        @test all(parent(deep)[k] == vals[EN.covering_position(mov, deepcv[k])]
                  for k in (1, length(deep) ÷ 2, length(deep)))
    end

    @testset "where the data was flat, the round trip is exact" begin
        # Constant per level-`L-1` sibling group, so `atol = 0` merges exactly
        # those groups.
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
        # Keywords reach the core: `minlevel = L` is the identity, and `by`
        # changes the values.
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

    # Each verb reads one axis shape and refuses the other.
    @test_throws ArgumentError DGG.aggregate(sum, M, L - 1)
    @test_throws ArgumentError DGG.coarsen(M; atol=1.0)
    @test_throws ArgumentError DGG.expand(A, L)

    # Only one dimension, and it must be `Cells`.
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
        # A reordered subset falls back to `Categorical`.
        @test !(lk[[3, 1]] isa DGG.MultiOrderLookup)
        @test lk[[3, 1]] isa DD.Lookups.Categorical
        @test collect(lk[[3, 1]]) == ids[[3, 1]]
        mask = falses(length(lk))
        mask[2] = mask[4] = true
        @test collect(lk[mask]) == ids[[2, 4]]
        @test_throws BoundsError lk[falses(length(lk) + 1)]
        # Indexing by a `SmallVector` (a neighbour list) needs the
        # SmallCollections ambiguity tie-break, as for `CellLookup`.
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
        # A cell deeper than the reference level re-keys the axis rather than
        # refusing it, so a concatenation of unequal depths still joins.
        deep = first(DGG.descendants(sys, ids[1], L + 1))
        dlk = DD.Lookups.rebuild(lk; data=[deep])
        @test dlk isa DGG.MultiOrderLookup
        @test collect(dlk) == [deep]
        @test DGG.reference_level(dlk) == L + 1
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
        # Halves with a gap between them still join to a container.
        gapped = vcat(M[1:(n÷2)], M[(n÷2+2):n])
        glk = DD.lookup(gapped, DGG.Cells)
        @test glk isa DGG.MultiOrderLookup
        @test collect(glk) == ids[[1:(n÷2); (n÷2+2):n]]
        @test glk != lk

        # A disordered join is decided at `rebuild` and falls back to
        # `Categorical`.
        @test DD.Lookups.rebuild(lk; data=vcat(ids[(n÷2+1):n], ids[1:(n÷2)])) isa
              DD.Lookups.Categorical

        # Two axes at different reference levels join at the deeper one.
        deeper = first(DGG.descendants(sys, ids[end], L + 1))
        tail = DD.DimArray([one(eltype(vals))],
            DGG.Cells(DGG.MultiOrderLookup(
                DGG.MultiOrderVector(sys, [deeper]; reference_level=L + 1))))
        mixed = vcat(M[1:(n-1)], tail)
        mlk = DD.lookup(mixed, DGG.Cells)
        @test mlk isa DGG.MultiOrderLookup
        @test DGG.reference_level(mlk) == L + 1
        @test collect(mlk) == vcat(ids[1:(n-1)], [deeper])

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

    # A coverage becomes an axis directly, with no data verb in between.
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
