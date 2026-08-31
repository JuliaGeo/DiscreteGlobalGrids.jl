# `DirectNearest` answers exactly what `NearestCell` answers.
#
# The method exists to skip the weight matrix, so the tests below assert two
# things and nothing else: the two methods agree element-wise on every route and
# every policy, and the direct plan holds nothing sized by a cell.

@testset "DirectNearest" begin

    # The destination is wider than the source, so its outer cells fall outside
    # coverage and exercise the "no source cell here" branch; its chunks span
    # the full longitude range, so each destination tile reads more than one
    # source chunk.
    dn_src = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
        chunks = (4, 2))
    dn_dst = ToyLonLatSpace(16, 8; lon = (-60.0, 60.0), lat = (-30.0, 30.0),
        chunks = (16, 4))
    dn_field = collect(reshape(1.0:32.0, 8, 4))

    # A source carrying both invalid kinds: NaN and a nodata sentinel.
    dn_poisoned = copy(dn_field)
    dn_poisoned[1:7:end] .= NaN
    dn_poisoned[3:11:end] .= -9999.0

    dn_eager(method, data; policy = Weighted(0.5), missingval = nothing) =
        regrid(data; to = dn_dst, from = dn_src, method = method,
            missingpolicy = policy, missingval = missingval, lazy = false)

    dn_lazyplan(method; policy = Weighted(0.5), missingval = nothing) =
        ChunkedPlan(method, policy, dn_dst, dn_src, PerChunk(),
            GR.DEFAULT_BUDGET, nothing, missingval)

    dn_lazy(method, data; policy = Weighted(0.5), missingval = nothing) =
        LazyRegridArray(data, dn_lazyplan(method; policy, missingval))[:]

    # `isequal`, so a NaN matches a NaN and a `missing` matches a `missing`.
    dn_same(a, b) = length(a) == length(b) && all(isequal(a[i], b[i]) for i in eachindex(a, b))

    @testset "method traits" begin
        @test GR.outputsampling(DirectNearest()) === GR.outputsampling(NearestCell())
        @test supportradius(DirectNearest(), dn_src) == supportradius(NearestCell(), dn_src)

        # The un-specialized route is `NearestCell`'s own, entry for entry.
        inds = 1:ncells(dn_dst)
        sinds = 1:ncells(dn_src)
        a = WeightCOO(length(inds))
        b = WeightCOO(length(inds))
        buildweights!(a, DirectNearest(), dn_dst, inds, dn_src, sinds)
        buildweights!(b, NearestCell(), dn_dst, inds, dn_src, sinds)
        @test (a.rows, a.cols, a.vals) == (b.rows, b.cols, b.vals)
    end

    @testset "eager plan and apply" begin
        direct = plan_regrid(dn_field; to = dn_dst, from = dn_src,
            method = DirectNearest(), lazy = false)
        weighted = plan_regrid(dn_field; to = dn_dst, from = dn_src,
            method = NearestCell(), lazy = false)
        @test direct isa GR.NearestDirectPlan
        @test weighted isa DirectPlan

        # `regrid`, the one-shot form, and `regrid!` all agree.
        @test dn_same(regrid(dn_field, direct), regrid(dn_field, weighted))
        @test dn_same(dn_eager(DirectNearest(), dn_field), dn_eager(NearestCell(), dn_field))
        outa = Array{Float64}(undef, ncells(dn_dst))
        outb = Array{Float64}(undef, ncells(dn_dst))
        regrid!(outa, dn_field, direct)
        regrid!(outb, dn_field, weighted)
        @test dn_same(outa, outb)

        # Destinations outside the source are blanked, and something is.
        @test count(isnan, regrid(dn_field, direct)) > 0

        # `applyplan!` writes every slice of a multi-slice source.
        stack = hcat(vec(dn_field), 2 .* vec(dn_field))
        @test dn_same(vec(regrid(reshape(stack, 8, 4, 2), direct)),
            vec(regrid(reshape(stack, 8, 4, 2), weighted)))
    end

    @testset "the plan holds nothing sized by a cell" begin
        # Two destinations 16x apart in cell count. The direct plan is the two
        # spaces and four scalars, so its size does not move; the weighted plan
        # carries a row per destination cell, so its size does.
        big = ToyLonLatSpace(64, 32; lon = (-60.0, 60.0), lat = (-30.0, 30.0),
            chunks = (64, 8))
        smalldirect = plan_regrid(dn_field; to = dn_dst, from = dn_src,
            method = DirectNearest(), lazy = false)
        bigdirect = plan_regrid(dn_field; to = big, from = dn_src,
            method = DirectNearest(), lazy = false)
        smallweighted = plan_regrid(dn_field; to = dn_dst, from = dn_src,
            method = NearestCell(), lazy = false)
        bigweighted = plan_regrid(dn_field; to = big, from = dn_src,
            method = NearestCell(), lazy = false)

        @test ncells(big) == 16 * ncells(dn_dst)
        @test Base.summarysize(bigdirect) == Base.summarysize(smalldirect)
        @test Base.summarysize(bigdirect) < 512
        @test Base.summarysize(bigweighted) > 8 * Base.summarysize(smallweighted)
    end

    @testset "lazy route" begin
        plan = dn_lazyplan(DirectNearest())
        # The point of this destination: a tile that reads more than one chunk.
        @test length(GR.sourcesof(GR.dependencies(plan), 1)) > 1
        @test dn_same(dn_lazy(DirectNearest(), dn_field),
            dn_lazy(NearestCell(), dn_field))
        # The lazy answer is the eager one.
        @test dn_same(dn_lazy(DirectNearest(), dn_field),
            vec(regrid(dn_field; to = dn_dst, from = dn_src,
                method = NearestCell(), lazy = false)))
        # A partial read agrees with the same window of the whole.
        whole = dn_lazy(DirectNearest(), dn_field)
        A = LazyRegridArray(dn_field, dn_lazyplan(DirectNearest()))
        @test dn_same(A[20:100], whole[20:100])
    end

    @testset "policies and sentinels" begin
        cases = ((Weighted(0.5), nothing), (Weighted(0.0), nothing),
            (Extensive(), nothing), (Weighted(0.5), -9999.0), (Extensive(), -9999.0))
        for (policy, mv) in cases, data in (dn_field, dn_poisoned)
            @test dn_same(dn_eager(DirectNearest(), data; policy, missingval = mv),
                dn_eager(NearestCell(), data; policy, missingval = mv))
            @test dn_same(dn_lazy(DirectNearest(), data; policy, missingval = mv),
                dn_lazy(NearestCell(), data; policy, missingval = mv))
        end
        # The poisoned source really does blank destinations the clean one does not.
        @test count(isnan, dn_eager(DirectNearest(), dn_poisoned)) >
              count(isnan, dn_eager(DirectNearest(), dn_field))

        # A declared sentinel is what the unmapped destinations come back
        # holding, on both routes.
        blanks = count(isnan, dn_eager(DirectNearest(), dn_field))
        @test blanks > 0
        @test count(==(-9999.0),
            dn_eager(DirectNearest(), dn_field; missingval = -9999.0)) == blanks
        @test count(==(-9999.0), LazyRegridArray(dn_field,
            dn_lazyplan(DirectNearest()); missingval = -9999.0)[:]) == blanks
    end

    @testset "missing-carrying sources" begin
        withmissing = Array{Union{Missing,Float64}}(dn_field)
        withmissing[2:5:end] .= missing
        @test dn_same(dn_eager(DirectNearest(), withmissing),
            dn_eager(NearestCell(), withmissing))
        @test dn_same(dn_lazy(DirectNearest(), withmissing),
            dn_lazy(NearestCell(), withmissing))
        @test count(ismissing, dn_eager(DirectNearest(), withmissing)) > 0
    end
end
