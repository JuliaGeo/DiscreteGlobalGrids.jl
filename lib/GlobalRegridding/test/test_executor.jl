# Plans, the eager executor, and the user API. Owned by task T5.
#
# Everything here runs against `toyspaces.jl`: the toy lon–lat space and the toy
# diagonal method, whose weights are `scale * I` and hand-checkable, so nothing
# in this file depends on a weight-construction task.

import DimensionalData as DD

# A method that counts the weight builds it is asked for.
mutable struct CountingMethod <: AbstractRegriddingMethod
    inner::ToyDiagonalMethod
    builds::Int
end

CountingMethod(; kw...) = CountingMethod(ToyDiagonalMethod(; kw...), 0)

function build_weights!(coo::WeightCOO, method::CountingMethod,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    method.builds += 1
    return build_weights!(coo, method.inner, dst_space, dst_inds, src_space, src_inds)
end

@testset "Executor" begin
    space = ToyLonLatSpace(6, 3)
    n = ncells(space)

    @testset "identity" begin
        field = rand(6, 3)
        method = ToyDiagonalMethod(; scale = 3.0)

        # `scale * I` with a denominator: the coverage-normalized mean is the
        # field itself, the raw sum is the field times the weight. Swapping the
        # numerator and the denominator, or the two policies, breaks one of
        # these two.
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

        # A destination fed only by valid sources is untouched by NaNs
        # elsewhere: the divisor is the valid weight, not the total weight.
        kept = setdiff(1:n, holes)
        @test punched[kept] ≈ whole[kept]
        # A destination fed only by invalid sources is blanked, not zero and not
        # a division by zero.
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

        # `Weighted` divides by the accumulated valid weight whether or not the
        # block declared a denominator, so a block with none still returns the
        # mean; `Extensive` returns the undivided sum.
        weighted = regrid(field; to = space, from = space, method,
            missingpolicy = Weighted(0.5))
        extensive = regrid(field; to = space, from = space, method,
            missingpolicy = Extensive())
        @test weighted ≈ vec(field)
        @test extensive ≈ 2 .* weighted

        # `Weighted` still blanks a destination whose stencil lost its mass to
        # invalid data, measured against the row sum rather than a denominator.
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

        # Every slice is what the 2-D call on that slice gives, through the very
        # same plan: the slice loop neither transposes nor reorders.
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
        # Applying a plan builds no weights, however often it is applied.
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
    end
end
