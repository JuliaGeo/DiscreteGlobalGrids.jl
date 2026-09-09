module PartitioningTests

using Test
import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import DimensionalData as DD
import Serialization

struct Alternating <: DGG.AbstractPartitioningAlgorithm end
DGG.partitionlabels(::Alternating, problem::DGG.PartitionProblem, nparts::Int,
    capacities::Vector{Float64}) = [mod1(i, nparts) for i in eachindex(problem.ids)]

struct BadLabels <: DGG.AbstractPartitioningAlgorithm
    labels::Vector{Int}
end
DGG.partitionlabels(algorithm::BadLabels, problem::DGG.PartitionProblem, nparts::Int,
    capacities::Vector{Float64}) = algorithm.labels

struct ToyInput
    rows::Vector{Vector{Int}}
end
DGG.partitionproblem(input::ToyInput; kwargs...) = DGG.PartitionProblem(input.rows; kwargs...)

@testset "partition problem" begin
    reads = [[3, 1, 3], [2]]
    weights = [2, 4]
    problem = DGG.PartitionProblem(reads; weights, sourceweights=[1, 2, 3],
        ids=[20, 10], sourceids=[300, 100, 200], order=[2, 1])
    reads[1][1] = 2
    weights[1] = 99
    @test problem.reads == [[1, 3], [2]]
    @test problem.weights == [2.0, 4.0]
    @test problem.sourceweights == [1.0, 2.0, 3.0]
    @test problem.ids == [20, 10]
    @test problem.sourceids == [300, 100, 200]
    @test problem.order == [2, 1]

    @test_throws ArgumentError DGG.PartitionProblem([[0]])
    @test_throws ArgumentError DGG.PartitionProblem([[2]]; sourceids=[1])
    @test_throws ArgumentError DGG.PartitionProblem([Int[]]; weights=[])
    @test_throws ArgumentError DGG.PartitionProblem([Int[]]; weights=[NaN])
    @test_throws ArgumentError DGG.PartitionProblem([Int[]]; weights=[-1])
    externalids = DGG.PartitionProblem([[1]]; ids=[0], sourceids=[-7])
    @test externalids.ids == [0]
    @test externalids.sourceids == [-7]
    @test_throws ArgumentError DGG.PartitionProblem([Int[], Int[]]; ids=[1, 1])
    @test_throws ArgumentError DGG.PartitionProblem([Int[], Int[]]; order=[2, 2])
    @test_throws ArgumentError DGG.PartitionProblem([Int[], Int[]]; order=[1])
    @test_throws ArgumentError DGG.PartitionProblem([Int[]];
        sourceweights=[1], sourceids=[1, 2])

    @test_throws ArgumentError DGG.MetisPartition(seed=-1)
    @test_throws ArgumentError DGG.MetisPartition(seed=big(typemax(Cint)) + 1)
    @test_throws ArgumentError DGG.MetisPartition(imbalance=Inf)
    @test_throws ArgumentError DGG.MetisPartition(imbalance=1.01)
    @test_throws ArgumentError DGG.MetisPartition(maxedges=0)
end

@testset "logical assignments" begin
    emptyresult = DGG.partition(DGG.PartitionProblem(Vector{Vector{Int}}()), 3)
    @test DGG.npartitions(emptyresult) == 3
    @test all(isempty, emptyresult.parts)

    one = DGG.partition(DGG.PartitionProblem([[1], [1]]), 1)
    @test one.assignment == [1, 1]
    @test DGG.partindices(one, 1) == [1, 2]

    sparse = DGG.partition(DGG.PartitionProblem([[1], [2]]), 4)
    @test DGG.npartitions(sparse) == 4
    @test sort(reduce(vcat, sparse.parts)) == [1, 2]
    @test count(isempty, sparse.parts) == 2

    zeroresult = DGG.partition(DGG.PartitionProblem([[1], [2], [3], [4]];
        weights=zeros(4)), 2)
    @test DGG.partweights(zeroresult) == [0.0, 0.0]
    @test length.(zeroresult.parts) == [2, 2]

    capacity = DGG.partition(DGG.PartitionProblem(fill(Int[], 6)), 2;
        capacities=[1, 2])
    @test length.(capacity.parts) == [2, 4]

    uneven = DGG.partition(DGG.PartitionProblem(fill(Int[], 4);
        weights=[5, 1, 1, 1]), 2)
    @test DGG.partweights(uneven) == [5.0, 3.0]

    orderedproblem = DGG.PartitionProblem([[2, 1], [2], [3], [1]];
        ids=[40, 30, 20, 10], sourceids=[900, 700, 800], order=[3, 1, 4, 2])
    ordered = DGG.partition(orderedproblem, 2)
    @test DGG.partindices(ordered, 1) == [3, 1]
    @test DGG.partindices(ordered, 2) == [4, 2]
    @test DGG.partchunks(ordered, 1) == [20, 40]
    @test DGG.partsources(ordered, 1) == [900, 700, 800]
    @test_throws BoundsError DGG.partchunks(ordered, 0)
    @test_throws BoundsError DGG.partsources(ordered, 3)

    alternating = DGG.partition(ToyInput([[1], [2], [3]]), 2;
        algorithm=Alternating(), ids=[10, 20, 30])
    @test alternating.assignment == [1, 2, 1]
    @test_throws ArgumentError DGG.partition(orderedproblem, 2;
        algorithm=BadLabels([1]))
    @test_throws ArgumentError DGG.partition(orderedproblem, 2;
        algorithm=BadLabels([1, 2, 3, 1]))
    @test_throws ArgumentError DGG.partition(orderedproblem, 0)
    @test_throws ArgumentError DGG.partition(orderedproblem, 2; capacities=[1])
    @test_throws ArgumentError DGG.partition(orderedproblem, 2; capacities=[1, 0])
    @test_throws ArgumentError DGG.partition(orderedproblem, 2;
        capacities=[floatmin(Float64), floatmax(Float64)])
    @test_throws ArgumentError DGG.partition(
        DGG.PartitionProblem([Int[], Int[]]; weights=fill(floatmax(Float64), 2)), 1)

    io = IOBuffer()
    Serialization.serialize(io, ordered)
    seekstart(io)
    restored = Serialization.deserialize(io)
    @test restored.ids == ordered.ids
    @test restored.assignment == ordered.assignment
    @test restored.parts == ordered.parts
    @test restored.chunks == ordered.chunks
    @test restored.sources == ordered.sources
    @test restored.weights == ordered.weights
    @test fieldtypes(typeof(restored)) == (Vector{Int}, Vector{Int},
        Vector{Vector{Int}}, Vector{Vector{Int}}, Vector{Vector{Int}}, Vector{Float64})
end

@testset "missing Metis backend" begin
    @test Base.get_extension(DGG, :DiscreteGlobalGridsMetisExt) === nothing
    error = try
        DGG.partition(DGG.PartitionProblem(Vector{Vector{Int}}()), 1;
            algorithm=DGG.MetisPartition())
        nothing
    catch caught
        caught
    end
    @test error isa DGG.PartitionBackendUnavailable
    rendered = sprint(showerror, error)
    @test occursin("using Metis", rendered)
    @test occursin("Pkg.add(\"Metis\")", rendered)
end

import Metis

@testset "MapChunkPlan adapter" begin
    grid = DGG.levelgrid(DGG.S2System(), 2)
    lookup = DGG.CellLookup(DGG.CellVector(grid))
    data = DD.DimArray(Float64.(1:DGG.ncells(grid)), DGG.Cells(lookup))
    plan = DGG.chunkplan(data; chunks=8)
    permuted = DGG.MapChunkPlan(plan.lookup, reverse(plan.chunks),
        plan.width, plan.connectivity)
    problem = DGG.partitionproblem(permuted)
    @test problem.ids == [chunk.index for chunk in reverse(plan.chunks)]
    @test problem.sourceids == problem.ids
    @test problem.weights == Float64[
        length(chunk.range) + length(chunk.halo) for chunk in reverse(plan.chunks)]
    @test problem.sourceweights == Float64[length(chunk.range) for chunk in reverse(plan.chunks)]
    @test all(row -> issorted(row) && allunique(row), problem.reads)

    assignment = DGG.partition(plan, 3; order=reverse(1:length(plan)))
    expected = zeros(length(data))
    actual = zeros(length(data))
    kernel(cell, value, neighbors) = value + sum(neighbors; init=0.0)
    DGG.mapneighbors!(expected, kernel, data, plan; threaded=false)
    visits = zeros(Int, length(data))
    for part in 1:DGG.npartitions(assignment)
        piece = plan[DGG.partindices(assignment, part)]
        DGG.mapneighbors!(actual, kernel, data, piece; threaded=false)
        for chunk in piece
            visits[DGG.ownedindices(chunk)] .+= 1
        end
    end
    @test actual == expected
    @test all(==(1), visits)
    @test sort(reduce(vcat, assignment.chunks)) == sort([chunk.index for chunk in plan])
    @test plan[[3, 1]][1].index == plan[3].index
    @test_throws BoundsError plan[[length(plan) + 1]]

    missing = DGG.MapChunkPlan(plan.lookup, plan.chunks[2:end],
        plan.width, plan.connectivity)
    @test_throws ArgumentError DGG.partitionproblem(missing)
    duplicate = DGG.MapChunkPlan(plan.lookup, [plan.chunks[1]; plan.chunks],
        plan.width, plan.connectivity)
    @test_throws ArgumentError DGG.partitionproblem(duplicate)
end

@testset "GlobalRegridding adapters" begin
    system = DGG.S2System()
    source = DGG.DGGSpace(DGG.levelgrid(system, 1); chunklevel=0)
    destination = DGG.DGGSpace(DGG.levelgrid(system, 2); chunklevel=1)
    graph = GR.chunk_dependency_graph(destination, source)
    problem = DGG.partitionproblem(graph; sourceids=collect(101:(100 + GR.nsourcechunks(graph))))
    @test problem.ids == collect(1:GR.ndestinationchunks(graph))
    @test length(problem.sourceids) == GR.nsourcechunks(graph)
    @test problem.reads == [collect(Int, GR.sourcesof(graph, row))
                            for row in 1:GR.ndestinationchunks(graph)]

    selected = [2, 4]
    restricted = GR.restrict(graph, selected)
    restrictedproblem = DGG.partitionproblem(restricted)
    @test restrictedproblem.ids == selected
    @test length(restrictedproblem.sourceids) == GR.nsourcechunks(graph)

    plan = GR.ChunkedPlan(GR.NearestCell(), GR.Weighted(1.0), destination, source;
        dependencies=graph)
    planproblem = DGG.partitionproblem(plan)
    @test planproblem.weights == Float64[
        length(GR.ownedindices(destination, row)) for row in 1:GR.ndestinationchunks(graph)]
    @test planproblem.sourceweights == Float64[
        length(GR.ownedindices(source, row)) for row in 1:GR.nsourcechunks(graph)]
    bare = GR.ChunkedPlan(GR.NearestCell(), GR.Weighted(1.0), destination, source;
        dependencies=false)
    @test_throws ArgumentError DGG.partitionproblem(bare)
end

@testset "Metis extension" begin
    @test Base.get_extension(DGG, :DiscreteGlobalGridsMetisExt) !== nothing

    clusters = DGG.PartitionProblem([[1], [1], [1], [2], [2], [2]])
    algorithm = DGG.MetisPartition(seed=7)
    labels = DGG.partition(clusters, 2; algorithm).assignment
    @test length(unique(labels[1:3])) == 1
    @test length(unique(labels[4:6])) == 1
    @test labels[1] != labels[4]
    @test labels == DGG.partition(clusters, 2; algorithm).assignment

    @test isempty(DGG.partition(DGG.PartitionProblem(Vector{Vector{Int}}()), 3;
        algorithm).assignment)
    noedges = DGG.partition(DGG.PartitionProblem([Int[], Int[], Int[]]), 2;
        algorithm)
    @test sort(reduce(vcat, noedges.parts)) == [1, 2, 3]
    @test DGG.partition(clusters, 1; algorithm).assignment == ones(Int, 6)
    @test DGG.npartitions(DGG.partition(clusters, 8; algorithm)) == 8

    huge = DGG.PartitionProblem([[1], [1], [2], [2]];
        weights=[floatmax(Float64) / 2, floatmax(Float64) / 4, 1.0, 0.0],
        sourceweights=[floatmax(Float64), floatmax(Float64) / 2])
    @test length(DGG.partition(huge, 2; algorithm).assignment) == 4
    disparate = DGG.PartitionProblem([[1], [1], [2], [2]];
        sourceweights=[nextfloat(0.0), 1.0, floatmax(Float64)])
    @test length(DGG.partition(disparate, 2; algorithm).assignment) == 4
    @test_throws ArgumentError DGG.partition(
        DGG.PartitionProblem(fill([1], 4)), 2;
        algorithm=DGG.MetisPartition(maxedges=2))

    savedoptions = copy(Metis.options)
    tasks = [Threads.@spawn DGG.partition(clusters, 2; algorithm).assignment for _ in 1:4]
    @test all(labels -> length(labels) == 6, fetch.(tasks))
    @test Metis.options == savedoptions
end

end
