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

function build_weights!(coo::WeightCOO, method::CountingMethod,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    countbuild!(method)
    return build_weights!(coo, method.inner, dst_space, dst_inds, src_space, src_inds)
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
        holes = [cellposition(space, 2, 2), cellposition(space, 5, 1)]
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
        h = cellposition(space, 3, 2)
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
        holes = [cellposition(space, 2, 2), cellposition(space, 5, 1)]
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
        ref = GR.blockreference!(Vector{Float64}(undef, n), block)
        x = vec(field)

        GR.applyblock!(num, cover, block, x, nothing, ref)
        @test (@allocated GR.applyblock!(num, cover, block, x, nothing, ref)) <= 128
        GR.applyblock!(num, cover, block, x, x, ref)
        @test (@allocated GR.applyblock!(num, cover, block, x, x, ref)) <= 128
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
        # Lazy-path knobs are accepted and ignored by the eager path.
        @test regrid(field; to = space, from = space, method, chunks = (3, 3),
            budget = 2^10) ≈ vec(field)
        @test_throws ArgumentError plan_regrid(field; to = space, from = space,
            method, budget = 0)
        # Weight storage is a lazy-path knob and says so rather than being
        # silently dropped by a plan that holds one whole-domain block.
        @test_throws ArgumentError plan_regrid(field; to = space, from = space,
            method, storage = PerChunk())
        # Default lazy storage is bounded by the weight budget.
        bounded = plan_regrid(field; to = space, from = space, method, lazy = true,
            budget = 2^16)
        @test bounded.storage.maxbytes == GR.weightbudget(2^16)
        @test GR.weightbudget(2^16) + GR.databudget(2^16) == 2^16
    end
end
