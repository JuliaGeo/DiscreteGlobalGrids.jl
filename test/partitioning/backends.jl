module PartitionBackendTests

using Test
import Random
import Serialization
import DiscreteGlobalGrids as DGG

@testset "optional partitioner constructors and load hints" begin
    for (constructor, backend, extension) in (
            (DGG.KaHyParPartition, :KaHyPar_jll, :DiscreteGlobalGridsKaHyParExt),
            (DGG.ScotchPartition, :Scotch, :DiscreteGlobalGridsScotchExt))
        @test_throws ArgumentError constructor(seed=-1)
        @test_throws ArgumentError constructor(seed=big(typemax(Cint)) + 1)
        @test_throws ArgumentError constructor(imbalance=NaN)
        @test_throws ArgumentError constructor(imbalance=-0.1)
        @test Base.get_extension(DGG, extension) === nothing
        error = try
            DGG.partition(DGG.PartitionProblem(Vector{Vector{Int}}()), 2;
                algorithm=constructor())
        catch caught
            caught
        end
        @test error isa DGG.PartitionBackendUnavailable
        @test error.backend === backend
        @test occursin("using $backend", sprint(showerror, error))
        @test occursin("Pkg.add(\"$backend\")", sprint(showerror, error))
    end
    @test_throws ArgumentError DGG.ScotchPartition(maxedges=0)
end

import KaHyPar_jll
import Scotch

function checkassignment(problem, result, k)
    @test DGG.npartitions(result) == k
    @test sort(reduce(vcat, result.parts; init=Int[])) == eachindex(problem.ids)
    @test result.ids == problem.ids
    for p in 1:k
        rows = filter(r -> result.assignment[r] == p, problem.order)
        @test DGG.partindices(result, p) == rows
        @test DGG.partchunks(result, p) == problem.ids[rows]
        sources = sort!(unique!(reduce(vcat, problem.reads[rows]; init=Int[])))
        @test DGG.partsources(result, p) == problem.sourceids[sources]
        @test DGG.partweights(result)[p] ≈ sum(problem.weights[rows]; init=0.0)
    end
end

@testset "native backend assignments" begin
    for algorithm in (DGG.KaHyParPartition(seed=17), DGG.ScotchPartition(seed=17))
        @testset "$(nameof(typeof(algorithm)))" begin
            problem = DGG.PartitionProblem([[1], [1], [1], [2], [2], [2]])
            result = DGG.partition(problem, 2; algorithm)
            @test length(unique(result.assignment[1:3])) == 1
            @test length(unique(result.assignment[4:6])) == 1
            @test result.assignment[1] != result.assignment[4]
            @test DGG.partition(problem, 2; algorithm).assignment == result.assignment

            uneven = DGG.PartitionProblem([[1], [1], [1], [2], [2], [2], [2], [2], [2]])
            capacity = DGG.partition(uneven, 2; algorithm, capacities=[1, 2])
            @test DGG.partweights(capacity) == [3, 6]
            checkassignment(uneven, capacity, 2)
            weighted = DGG.PartitionProblem(problem.reads; weights=[2, 2, 2, 1, 1, 1])
            @test DGG.partweights(DGG.partition(weighted, 2;
                algorithm, capacities=[1, 2])) == [3, 6]
            sharing = DGG.PartitionProblem([[1, 3], [1, 4], [2, 3], [2, 4]];
                sourceweights=[1, 1, 100, 100])
            affinity = DGG.partition(sharing, 2; algorithm).assignment
            @test affinity[1] == affinity[3]
            @test affinity[2] == affinity[4]
            @test affinity[1] != affinity[2]
            checkassignment(problem, DGG.partition(problem, 2;
                algorithm, capacities=[nextfloat(0.0), 1.0]), 2)


            for reads in (Vector{Vector{Int}}(), [Int[]], fill(Int[], 6), fill([1], 6))
                p = DGG.PartitionProblem(reads; weights=zeros(length(reads)))
                for k in (1, 2, 8)
                    checkassignment(p, DGG.partition(p, k; algorithm), k)
                end
            end

            rng = Random.MersenneTwister(42)
            for _ in 1:24
                n = rand(rng, 2:20)
                m = rand(rng, 1:8)
                reads = [findall(rand(rng, m) .< 0.35) for _ in 1:n]
                p = DGG.PartitionProblem(reads;
                    weights=Float64.(rand(rng, 0:5, n)),
                    sourceweights=Float64.(rand(rng, 0:5, m)),
                    ids=collect(-n:-1), sourceids=collect(-m:-1),
                    order=Random.randperm(rng, n))
                k = rand(rng, 2:4)
                result = DGG.partition(p, k; algorithm, capacities=rand(rng, k) .+ 0.2)
                checkassignment(p, result, k)
            end
            disparate = DGG.PartitionProblem([[1], [1], [2], [2]];
                weights=[floatmax(Float64)/4, 1.0, 0.0, nextfloat(0.0)],
                sourceweights=[nextfloat(0.0), 1.0, floatmax(Float64)])
            checkassignment(disparate, DGG.partition(disparate, 2; algorithm), 2)
            tasks = [Threads.@spawn DGG.partition(problem, 2; algorithm).assignment for _ in 1:4]
            @test all(==(DGG.partition(problem, 2; algorithm).assignment), fetch.(tasks))
            io = IOBuffer()
            Serialization.serialize(io, capacity)
            seekstart(io)
            @test Serialization.deserialize(io).assignment == capacity.assignment
        end
    end
end

@testset "Scotch private random context" begin
    Scotch.random_seed(53)
    Scotch.random_reset()
    expected = Scotch.random_val()
    Scotch.random_reset()
    DGG.partition(DGG.PartitionProblem([[1], [1], [2], [2]]), 2;
        algorithm=DGG.ScotchPartition())
    @test Scotch.random_val() == expected
end

@testset "shared-source representation" begin
    problem = DGG.PartitionProblem([[1, 2], [1], [1, 3], [4]];
        sourceweights=[7, 3, 0, 5, 9])
    ext = Base.get_extension(DGG, :DiscreteGlobalGridsKaHyParExt)
    offsets, pins, weights = ext._hypergraph(problem)
    @test offsets == [0, 3]
    @test pins == [0, 1, 2]
    @test weights == [7]
    dense = DGG.PartitionProblem(fill([1], 1001))
    @test_throws ArgumentError DGG.partition(dense, 2;
        algorithm=DGG.ScotchPartition(maxedges=1000))
    checkassignment(dense, DGG.partition(dense, 2; algorithm=DGG.KaHyParPartition()), 2)
end

end
