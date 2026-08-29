# Eager plans, execution, and public API.

import DimensionalData as DD
import SparseArrays
# Load and test `GlobalRegriddingRastersExt`.
import Rasters

# Weight-build counter.
mutable struct CountingMethod <: AbstractRegriddingMethod
    inner::ToyDiagonalMethod
    builds::Int
end

CountingMethod(; kw...) = CountingMethod(ToyDiagonalMethod(; kw...), 0)

function buildweights!(coo::WeightCOO, method::CountingMethod,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    countbuild!(method)
    return buildweights!(coo, method.inner, dst_space, dst_inds, src_space, src_inds)
end

@testset "Executor" begin
    space = ToyLonLatSpace(6, 3)
    n = ncells(space)

    @testset "identity" begin
        field = rand(6, 3)
        method = ToyDiagonalMethod(; scale = 3.0)

        # Weighted returns means; Extensive returns scaled sums.
        mean = regrid(field; to = space, from = space, method,
            missingpolicy = Weighted())
        @test mean ≈ vec(field)
        sums = regrid(field; to = space, from = space, method,
            missingpolicy = Extensive())
        @test sums ≈ 3 .* vec(field)

        # A 2-D source flattens to a plain vector over the destination's cells.
        @test mean isa Vector{Float64}
        @test length(mean) == n

        dest = fill(NaN, n)
        regrid!(dest, field; to = space, from = space, method,
            missingpolicy = Extensive())
        @test dest ≈ sums
    end

    @testset "missing data, block with a denominator" begin
        field = rand(6, 3) .+ 1
        method = ToyDiagonalMethod(; scale = 3.0)
        kw = (; to = space, from = space, method, missingpolicy = Weighted(0.5))

        whole = regrid(field; kw...)
        holed = copy(field)
        holes = [localindex(space, 2, 2), localindex(space, 5, 1)]
        holed[2, 2] = NaN
        holed[5, 1] = NaN
        punched = regrid(holed; kw...)

        # Unrelated NaNs do not affect valid destinations.
        kept = setdiff(1:n, holes)
        @test punched[kept] ≈ whole[kept]
        # All-invalid destinations are blanked.
        @test all(isnan, punched[holes])

        # `missing` in, `missing` out.
        withmissing = Array{Union{Missing,Float64}}(field)
        withmissing[2, 2] = missing
        mixed = regrid(withmissing; kw...)
        @test ismissing(mixed[holes[1]])
        @test Vector{Float64}(mixed[kept]) ≈ whole[kept]

        # `Extensive` blanks nothing: the lost mass is a raw sum of zero.
        raw = regrid(holed; to = space, from = space, method,
            missingpolicy = Extensive())
        @test raw[holes] == [0.0, 0.0]
    end

    @testset "missing data, block without a denominator" begin
        field = rand(6, 3) .+ 1
        method = ToyDiagonalMethod(; scale = 2.0, withdenom = false)

        # Blocks without denominators still normalize by valid row weight.
        weighted = regrid(field; to = space, from = space, method,
            missingpolicy = Weighted(0.5))
        extensive = regrid(field; to = space, from = space, method,
            missingpolicy = Extensive())
        @test weighted ≈ vec(field)
        @test extensive ≈ 2 .* weighted

        # Row sums provide coverage thresholds without denominators.
        holed = copy(field)
        holed[3, 2] = NaN
        h = localindex(space, 3, 2)
        blanked = regrid(holed; to = space, from = space, method,
            missingpolicy = Weighted(0.5))
        rest = setdiff(1:n, (h,))
        @test isnan(blanked[h])
        @test blanked[rest] ≈ vec(field)[rest]
        @test regrid(holed; to = space, from = space, method,
            missingpolicy = Extensive())[h] == 0.0
    end

    @testset "missingval sentinel" begin
        # Declared sentinels behave exactly like NaN under both policies.
        field = rand(6, 3) .+ 1
        holes = [localindex(space, 2, 2), localindex(space, 5, 1)]
        sentinel = copy(field)
        sentinel[2, 2] = -9999.0
        sentinel[5, 1] = -9999.0
        nanned = copy(field)
        nanned[2, 2] = NaN
        nanned[5, 1] = NaN

        for policy in (Weighted(0.5), Extensive())
            kw = (; to = space, from = space,
                method = ToyDiagonalMethod(; scale = 3.0), missingpolicy = policy)
            @test all(isequal.(regrid(sentinel; kw..., missingval = -9999.0),
                regrid(nanned; kw...)))
        end

        # Undeclared sentinels remain ordinary data.
        plain = regrid(sentinel; to = space, from = space,
            method = ToyDiagonalMethod(), missingpolicy = Extensive())
        @test plain[holes] == [-9999.0, -9999.0]

        # Metadata detects normalized and CF nodata keys as strings or symbols.
        withmeta(key) = DD.DimArray(sentinel, (DD.X(1:6), DD.Y(1:3));
            metadata = DD.Metadata(Dict(key => -9999.0)))
        @test all(GR.sourcemissingval(withmeta(k)) == -9999.0
                  for k in ("missingval", "_FillValue", "missing_value",
            :missingval, :_FillValue, :missing_value))
        @test GR.sourcemissingval(sentinel) === nothing
        @test GR.sourcemissingval(DD.DimArray(sentinel, (DD.X(1:6), DD.Y(1:3)))) === nothing
        @test GR.sourcemissingval(withmeta("units")) === nothing

        declared = withmeta("_FillValue")
        flat(x) = vec(x isa DD.AbstractDimArray ? parent(x) : x)
        dkw = (; to = space, from = space, method = ToyDiagonalMethod(; scale = 3.0),
            missingpolicy = Extensive())
        @test all(isequal.(flat(regrid(declared; dkw...)), flat(regrid(nanned; dkw...))))

        # Explicit `missingval` overrides metadata.
        @test flat(regrid(declared; to = space, from = space,
            method = ToyDiagonalMethod(), missingpolicy = Extensive(),
            missingval = nothing)) == plain

        # The Rasters extension reads field-based sentinels and normalizes absence.
        let dims = (DD.X(1:6), DD.Y(1:3))
            @test GR.sourcemissingval(Rasters.Raster(sentinel, dims;
                missingval = -9999.0)) == -9999.0
            @test GR.sourcemissingval(Rasters.Raster(nanned, dims;
                missingval = NaN)) === NaN
            @test GR.sourcemissingval(Rasters.Raster(sentinel, dims;
                missingval = nothing)) === nothing
            @test GR.sourcemissingval(Rasters.Raster(
                convert(Matrix{Union{Missing,Float64}}, sentinel), dims;
                missingval = missing)) === nothing
            # Raster fields take precedence over metadata.
            @test GR.sourcemissingval(Rasters.Raster(sentinel, dims;
                missingval = -9999.0, metadata = DD.Metadata(Dict("_FillValue" => 0.0)))) ==
                  -9999.0
            # Raster sentinels reach plans without an explicit keyword.
            rkw = (; to = space, from = space,
                method = ToyDiagonalMethod(; scale = 3.0), missingpolicy = Extensive())
            @test all(isequal.(flat(regrid(Rasters.Raster(sentinel, dims;
                    missingval = -9999.0); rkw...)), flat(regrid(nanned; rkw...))))
        end

        # Declared integer sentinels force a validity scan.
        ints = rand(1:100, 6, 3)
        ints[2, 2] = -9999
        out = regrid(ints; to = space, from = space, method = ToyDiagonalMethod(),
            missingpolicy = Weighted(0.5), missingval = -9999)
        @test isnan(out[holes[1]])
        @test out[setdiff(1:n, holes[1])] ≈ vec(ints)[setdiff(1:n, holes[1])]
    end

    @testset "N-D pass-through" begin
        cube = DD.DimArray(rand(6, 3, 12), (DD.X(1:6), DD.Y(1:3), DD.Ti(1:12)))
        plan = plan_regrid(cube; to = space, from = space,
            method = ToyDiagonalMethod(; scale = 2.0), missingpolicy = Extensive())
        out = regrid(cube, plan)

        # Destination cells first, the non-spatial dimensions after, unchanged.
        @test out isa DD.AbstractDimArray
        @test size(out) == (n, 12)
        @test DD.hasdim(out, DD.Ti)
        @test collect(DD.lookup(out, DD.Ti)) == 1:12
        @test !DD.hasdim(out, DD.X) && !DD.hasdim(out, DD.Y)

        # N-D slices match independent 2-D applications.
        for k in 1:12
            @test Array(out)[:, k] ≈ regrid(parent(cube)[:, :, k], plan)
        end
    end

    @testset "plan reuse" begin
        field = rand(6, 3)
        method = CountingMethod(; scale = 1.5)
        plan = plan_regrid(field; to = space, from = space, method,
            missingpolicy = Extensive())
        @test method.builds == 1

        once = regrid(field, plan)
        twice = regrid(field, plan)
        @test once == twice
        @test once ≈ 1.5 .* vec(field)
        # Reusing a plan does not rebuild weights.
        @test method.builds == 1
    end

    @testset "accumulation is allocation-free" begin
        field = rand(6, 3)
        plan = plan_regrid(field; to = space, from = space,
            method = ToyDiagonalMethod(), missingpolicy = Weighted())
        block = plan.block
        num = zeros(n)
        cover = zeros(n)
        ref = block.reference
        x = vec(field)

        GR.applyblock!(num, cover, block, x, nothing, ref)
        @test (@allocated GR.applyblock!(num, cover, block, x, nothing, ref)) <= 128
        GR.applyblock!(num, cover, block, x, x, ref)
        @test (@allocated GR.applyblock!(num, cover, block, x, x, ref)) <= 128
    end

    @testset "one reference vector lives in the final block" begin
        inds = ownedindices(space, 1)
        m = length(inds)

        # A denominated block references its denominator itself: one vector, not
        # a copy of one, so nothing can hold a stale second opinion of it.
        coo = WeightCOO(m)
        buildweights!(coo, ToyDiagonalMethod(; scale = 2.0), space, inds, space, inds)
        denominated = WeightBlock(coo, m, m)
        @test denominated.denom == fill(2.0, m)
        @test denominated.reference === denominated.denom

        # A method reporting no denominator allocates none at any point, and its
        # block references the weights' row sums instead of zeros.
        bare = WeightCOO(m)
        @test bare.denom === nothing
        buildweights!(bare, ToyDiagonalMethod(; scale = 2.0, withdenom = false),
            space, inds, space, inds)
        @test bare.denom === nothing
        point = WeightBlock(bare, m, m)
        @test point.denom === nothing
        @test point.reference == fill(2.0, m)
        @test point.weights == denominated.weights

        # Dense and sparse weights of the same operator agree on the reference.
        dense = WeightBlock(Matrix(point.weights), nothing)
        @test dense.reference == point.reference
        @test dense.reference !== point.reference

        # A row no source reaches has reference zero, and the policy decides
        # what that destination becomes — the block states the fact only.
        holed = WeightBlock(SparseArrays.sparse([1, 3], [1, 3], [2.0, 2.0], 3, 3), nothing)
        @test holed.reference == [2.0, 0.0, 2.0]
        num = zeros(3)
        cover = zeros(3)
        GR.applyblock!(num, cover, holed, ones(3))
        @test cover == holed.reference
        blanked = Vector{Union{Missing,Float64}}(undef, 3)
        GR.finalize!(blanked, num, cover, holed.reference, Weighted(0.5))
        @test ismissing(blanked[2])
        GR.finalize!(blanked, num, cover, holed.reference, Extensive())
        @test blanked[2] == 0.0

        # Empty sides build and apply. A block with no destinations has an empty
        # reference; one with no sources has a zero reference of full length.
        nodst = WeightBlock(SparseArrays.sparse(Int[], Int[], Float64[], 0, 4), nothing)
        @test isempty(nodst.reference)
        @test GR.applyblock!(Float64[], Float64[], nodst, zeros(4)) == Float64[]
        nosrc = WeightBlock(SparseArrays.sparse(Int[], Int[], Float64[], 3, 0), nothing)
        @test nosrc.reference == zeros(3)
        n0, c0 = zeros(3), zeros(3)
        GR.applyblock!(n0, c0, nosrc, Float64[])
        @test n0 == zeros(3) && c0 == zeros(3)
        # A builder that declares denominators for an empty destination still
        # produces a denominated block, and it references its own empty vector.
        empty = WeightBlock(GR.markdenominated!(WeightCOO(0)), 0, 4)
        @test empty.denom == Float64[]
        @test empty.reference === empty.denom

        # Resident bytes are the weights plus exactly one vector, whichever kind
        # of reference the block carries.
        vecbytes = 8 * m
        wbytes = 16 * SparseArrays.nnz(denominated.weights) +
                 8 * (size(denominated.weights, 2) + 1)
        @test GR._blockbytes(denominated) == wbytes + vecbytes + 64
        @test GR._blockbytes(point) == wbytes + vecbytes + 64
    end

    @testset "signed weights, coverage of their own" begin
        # A second-order stencil: non-negative overlap areas plus a zero-row-sum
        # gradient correction. The values are signed; the coverage is the areas.
        W = SparseArrays.sparse([1, 1, 2, 2], [1, 2, 2, 3],
            [1.2, -0.2, 0.8, 0.2], 2, 3)
        C = SparseArrays.sparse([1, 1, 2, 2], [1, 2, 2, 3],
            [0.8, 0.2, 0.5, 0.5], 2, 3)
        denom = [1.0, 1.0]
        block = WeightBlock(W, denom, C)
        @test block.coverage === C
        @test block.reference === block.denom
        @test contains(repr(block), "denom, coverage")

        # Clean sources take the reference, exactly as an unsigned block does.
        src = [4.0, 10.0, 20.0]
        num = zeros(2)
        cover = zeros(2)
        GR.applyblock!(num, cover, block, src)
        @test num ≈ Matrix(W) * src
        @test cover == denom

        # A hole at source 2 leaves each destination covered by its remaining
        # area — never by the signed weights, which here would report 1.2 and
        # 0.2 and misnormalize both destinations.
        holed = [4.0, NaN, 20.0]
        fill!(num, 0.0)
        fill!(cover, 0.0)
        GR.applyblock!(num, cover, block, holed, holed)
        @test num ≈ [1.2 * 4.0, 0.2 * 20.0]
        @test cover ≈ [0.8, 0.5]
        @test cover ≉ [1.2, 0.2]

        # `Weighted` divides by that coverage and blanks against the reference:
        # destination 2 keeps half its area and falls below a 0.6 threshold.
        out = Vector{Union{Missing,Float64}}(undef, 2)
        GR.finalize!(out, num, cover, block.reference, Weighted(0.5))
        @test out ≈ [1.2 * 4.0 / 0.8, 0.2 * 20.0 / 0.5]
        GR.finalize!(out, num, cover, block.reference, Weighted(0.6))
        @test out[1] ≈ 1.2 * 4.0 / 0.8
        @test ismissing(out[2])
        # `Extensive` reports the undivided sums, coverage or none.
        GR.finalize!(out, num, cover, block.reference, Extensive())
        @test out ≈ num

        # An all-valid mask reaches the coverage the clean path took.
        fill!(num, 0.0)
        fill!(cover, 0.0)
        GR.applyblock!(num, cover, block, src, src)
        @test num ≈ Matrix(W) * src
        @test cover ≈ denom
        # And accumulating coverage separately allocates nothing either.
        @test (@allocated GR.applyblock!(num, cover, block, src, src)) <= 128

        # Dense weights and dense coverage answer the same, and both sparse
        # column walks agree with the dense loop over an empty-column stretch.
        dense = WeightBlock(Matrix(W), denom, Matrix(C))
        n2, c2 = zeros(2), zeros(2)
        GR.applyblock!(n2, c2, dense, holed, holed)
        @test n2 ≈ [1.2 * 4.0, 0.2 * 20.0]
        @test c2 ≈ [0.8, 0.5]
        for ncols in (20, 20_000)
            rows, cs = [1, 1, 2], [4, 9, 12]
            Wn = SparseArrays.sparse(rows, cs, [1.5, -0.5, 2.0], 2, ncols)
            Cn = SparseArrays.sparse(rows, cs, [0.7, 0.3, 2.0], 2, ncols)
            @test GR._walknonzeros(Cn) == (ncols == 20_000)
            wide = WeightBlock(Wn, [1.0, 2.0], Cn)
            x = collect(1.0:ncols)
            x[9] = NaN
            n3, c3 = zeros(2), zeros(2)
            GR.applyblock!(n3, c3, wide, x, x)
            @test n3 ≈ Matrix(Wn) * map(v -> isnan(v) ? 0.0 : v, x)
            @test c3 ≈ Matrix(Cn) * map(v -> isnan(v) ? 0.0 : 1.0, x)
        end

        # Coverage is indexed like the values it stands in for.
        @test_throws DimensionMismatch WeightBlock(W, denom,
            SparseArrays.sparse(Int[], Int[], Float64[], 2, 4))
    end

    @testset "coverage through the COO" begin
        coo = WeightCOO(2)
        @test !GR.hascoverage(coo)
        for (j, k, w, a) in ((1, 1, 1.2, 0.8), (1, 2, -0.2, 0.2),
            (2, 2, 0.8, 0.5), (2, 3, 0.2, 0.5))
            addweight!(coo, j, k, w)
            GR.addcoverage!(coo, j, k, a)
        end
        # Duplicate entries are summed on both lists alike.
        addweight!(coo, 2, 3, 0.1)
        GR.addcoverage!(coo, 2, 3, 0.1)
        adddenom!(coo, 1, 1.0)
        adddenom!(coo, 2, 1.0)
        @test GR.hascoverage(coo)
        @test contains(repr(coo), "denom, coverage")

        block = WeightBlock(coo, 2, 3)
        @test Matrix(block.weights) ≈ [1.2 -0.2 0.0; 0.0 0.8 0.3]
        @test Matrix(block.coverage) ≈ [0.8 0.2 0.0; 0.0 0.5 0.6]
        @test block.reference === block.denom

        # A builder that declares coverage but reports none covers nothing: the
        # block reads a zero operator rather than falling back to its weights.
        empty = WeightCOO(2)
        addweight!(empty, 1, 1, 1.2)
        GR.markcovered!(empty)
        @test GR.hascoverage(empty)
        zeroed = WeightBlock(empty, 2, 3)
        @test SparseArrays.nnz(zeroed.coverage) == 0
        num, cover = zeros(2), zeros(2)
        GR.applyblock!(num, cover, zeroed, ones(3), ones(3))
        @test num ≈ [1.2, 0.0]
        @test cover == zeros(2)

        # A COO that reports no coverage builds the block it always did, and its
        # weights go on serving as their own coverage.
        plain = WeightCOO(2)
        addweight!(plain, 1, 1, 1.0)
        addweight!(plain, 2, 3, 3.0)
        @test !GR.hascoverage(plain)
        bare = WeightBlock(plain, 2, 3)
        @test bare.coverage === nothing
        @test bare.denom === nothing
        @test bare.reference == [1.0, 3.0]
        @test !contains(repr(bare), "coverage")
        n4, c4 = zeros(2), zeros(2)
        GR.applyblock!(n4, c4, bare, [2.0, 0.0, NaN], [2.0, 0.0, NaN])
        @test n4 ≈ [2.0, 0.0]
        @test c4 ≈ [1.0, 0.0]
    end

    @testset "repeated eager application allocates no reference" begin
        wide = ToyLonLatSpace(40, 20)
        nwide = Int(ncells(wide))
        plan = plan_regrid(zeros(40, 20); to = wide, from = wide,
            method = ToyDiagonalMethod(), missingpolicy = Weighted())
        @test plan.block.reference === plan.block.denom
        src = reshape(collect(1.0:nwide), nwide, 1)
        dst = zeros(Float64, nwide, 1)

        GR.applyplan!(dst, plan, src)
        first = copy(dst)
        # Two accumulators and no third vector: a per-application reference would
        # cost another `8 * nwide` bytes on top.
        @test (@allocated GR.applyplan!(dst, plan, src)) < 3 * 8 * nwide skip = VERSION < v"1.12"
        @test dst == first
    end

    @testset "sparse column walks agree" begin
        # Sparse and dense traversal paths agree across empty columns.
        rows = [1, 3, 1, 2]
        colvals = [4, 4, 9, 12]
        vals = [2.0, 5.0, 1.5, 0.5]
        for ncols in (20, 20_000)
            W = SparseArrays.sparse(rows, colvals, vals, 3, ncols)
            @test GR._walknonzeros(W) == (ncols == 20_000)
            block = WeightBlock(W, nothing)
            src = collect(1.0:ncols)
            src[9] = NaN
            num = zeros(3)
            cover = zeros(3)
            GR.applyblock!(num, cover, block, src, src)
            dense = Matrix(W)
            @test num ≈ dense * map(x -> isnan(x) ? 0.0 : x, src)
            @test cover ≈ dense * map(x -> isnan(x) ? 0.0 : 1.0, src)

            # And the no-mask path, whose coverage comes off the reference.
            n2 = zeros(3)
            c2 = zeros(3)
            clean = collect(1.0:ncols)
            GR.applyblock!(n2, c2, block, clean)
            @test n2 ≈ dense * clean
            @test c2 ≈ vec(sum(dense; dims = 2))
        end
    end

    @testset "API surface" begin
        field = rand(6, 3)
        method = ToyDiagonalMethod()

        # `to` is a space at this layer, and says so.
        @test_throws ArgumentError plan_regrid(field; to = (6, 3), from = space, method)
        # `lazy = true` plans to a chunked plan and builds no weights doing so.
        @test plan_regrid(field; to = space, from = space,
            method, lazy = true) isa ChunkedPlan
        # A source that does not flatten to the space's cells is caught before
        # any weight is applied.
        @test_throws DimensionMismatch regrid(rand(5, 3); to = space, from = space, method)
        # The eager path names the lazy-only knobs it was handed rather than
        # dropping them silently.
        @test_throws "`chunks` or `budget`" plan_regrid(field; to = space,
            from = space, method, chunks = (3, 3), budget = 2^10)
        @test_throws "`storage`" plan_regrid(field; to = space, from = space,
            method, storage = PerChunk())
        @test_throws ArgumentError plan_regrid(field; to = space, from = space,
            method, lazy = true, budget = 0)
        # Default lazy storage is bounded by the weight budget.
        bounded = plan_regrid(field; to = space, from = space, method, lazy = true,
            budget = 2^16)
        @test bounded.storage.maxbytes == GR.weightbudget(2^16)
        @test GR.weightbudget(2^16) + GR.databudget(2^16) == 2^16
    end
end
