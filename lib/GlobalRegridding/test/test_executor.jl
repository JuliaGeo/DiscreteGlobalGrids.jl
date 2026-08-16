# Plans, the eager executor, and the user API. Owned by task T5.
#
# Everything here runs against `toyspaces.jl`: the toy lon–lat space and the toy
# diagonal method, whose weights are `scale * I` and hand-checkable, so nothing
# in this file depends on a weight-construction task.

import DimensionalData as DD
import SparseArrays

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

    @testset "missingval sentinel" begin
        # A sentinel is nothing but a third spelling of invalid, so a field
        # holding one must give exactly what the same field gives with its
        # sentinels replaced by NaN — bit for bit, under either policy. A
        # sentinel folded into the weights, or into the coverage but not the
        # numerator, breaks one of these.
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

        # Without the declaration the sentinel is just data — the guard against
        # a `missingval` that leaks into a plan that was never given one.
        plain = regrid(sentinel; to = space, from = space,
            method = ToyDiagonalMethod(), missingpolicy = Extensive())
        @test plain[holes] == [-9999.0, -9999.0]

        # A source that declares its own sentinel needs no keyword: the default
        # of `missingval` is what the source says, so a reader that carries
        # nodata in its metadata is handled the same as one whose caller
        # repeats it. Both CF spellings and the normalized one are read, under
        # a String key or a Symbol one; a source that declares nothing answers
        # nothing, which is what keeps the sentinel-as-data case above honest.
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

        # …and the caller still overrides it, back to the sentinel-as-data
        # answer. Without this the default could be an unconditional read of the
        # metadata and nothing would notice.
        @test flat(regrid(declared; to = space, from = space,
            method = ToyDiagonalMethod(), missingpolicy = Extensive(),
            missingval = nothing)) == plain

        # An integer field cannot hold NaN or `missing`, so a sentinel is the
        # only way it can carry nodata at all: the element-type shortcut that
        # skips the validity scan must not fire when one is declared.
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

    @testset "sparse column walks agree" begin
        # The accumulation walks a near-empty block nonzero-first and a denser
        # one column-by-column. Both must give the dense answer, and the
        # nonzero-first walk has to recover a column from a nonzero's position
        # through repeated `colptr` entries — leading, interior and trailing
        # empty columns are exactly where that goes wrong.
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
        # A lazy plan bounds its weights by the budget's weight share unless the
        # caller names a storage. This is the L3 failure P1 measured: an
        # unbounded default accumulates the whole operator.
        bounded = plan_regrid(field; to = space, from = space, method, lazy = true,
            budget = 2^16)
        @test bounded.storage.maxbytes == GR.weightbudget(2^16)
        @test GR.weightbudget(2^16) + GR.databudget(2^16) == 2^16
    end
end
