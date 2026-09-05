module PartitionDistributedTests

using Test
import Distributed
import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import Metis

const FIXTURE = joinpath(@__DIR__, "worker_fixture.jl")
if !isdefined(Main, :PartitionWorkerFixtures)
    Base.include(Main, FIXTURE)
end
const Fixtures = Main.PartitionWorkerFixtures

@testset "Copernicus regridding partition metadata" begin
    system = DGG.CopernicusDEMSystem(90)
    tile = DGG.cellat(DGG.levelgrid(system, 0), 0.5, 0.5)
    source = DGG.DGGSpace(DGG.subtree(system, tile, 1); chunklevel=0)
    targetsystem = DGG.IGeo7System()
    root = DGG.cellat(DGG.levelgrid(targetsystem, 4), 0.5, 0.5)
    destination = DGG.DGGSpace(DGG.subtree(targetsystem, root, 6); chunklevel=5)
    plan = GR.ChunkedPlan(GR.NearestCell(), GR.Weighted(1.0), destination, source;
        dependencies=true)
    relation = GR.dependencies(plan)
    problem = DGG.partitionproblem(plan)
    assignment = DGG.partition(problem, 2; algorithm=DGG.MetisPartition())
    @test GR.ncells(source) == 1_440_000
    @test problem.sourceweights == [GR.ncells(source)]
    @test length(problem.ids) == GR.nchunks(destination) == 7
    @test problem.reads == [collect(GR.sourcesof(relation, row)) for row in 1:7]
    @test GR.dependencies(plan) === relation
    @test sort(reduce(vcat, [DGG.partchunks(assignment, p) for p in 1:2])) == 1:7
    @test all(sources -> all(==(1), sources),
        [DGG.partsources(assignment, p) for p in 1:2])
end

@testset "partition assignments across processes" begin
    project = dirname(Base.active_project())
    processes = Distributed.addprocs(2; exeflags=`--project=$project --threads=1`)
    try
        @sync for process in processes
            @async Distributed.remotecall_fetch(Base.include, process, Main, FIXTURE)
        end

        @testset "stencil parts preserve owned outputs" begin
            data, plan = Fixtures.stencil_fixture()
            expected = zeros(size(data))
            DGG.mapneighbors!(expected, Fixtures.stencil, data, plan; threaded=false)
            assignment = DGG.partition(plan, 3;
                algorithm=DGG.MetisPartition(), capacities=[1.0, 2.0, 1.0])
            futures = [Distributed.remotecall(Fixtures.stencil_part,
                processes[mod1(part, length(processes))], assignment, part)
                for part in 1:DGG.npartitions(assignment)]
            received = fetch.(futures)
            result = fill(NaN, size(expected))
            visits = zeros(Int, length(expected))
            for part in received
                result[part.indices] = part.values
                visits[part.indices] .+= 1
                @test !part.metis_loaded
            end
            @test all(==(1), visits)
            @test result == expected
            @test sort(reduce(vcat, getproperty.(received, :chunks))) ==
                  sort([chunk.index for chunk in plan])
        end

        @testset "regridding assigns destinations and retains source dependencies" begin
            data, plan = Fixtures.regrid_fixture()
            expected = collect(GR.regrid(data, plan))
            assignment = DGG.partition(plan, 3; algorithm=DGG.MetisPartition())
            futures = [Distributed.remotecall(Fixtures.regrid_part,
                processes[mod1(part + 1, length(processes))], assignment, part)
                for part in 1:DGG.npartitions(assignment)]
            received = fetch.(futures)
            result = fill(NaN, size(expected))
            visits = zeros(Int, length(expected))
            for part in received
                result[part.indices] = part.values
                visits[part.indices] .+= 1
                needed = sort!(unique!(reduce(vcat,
                    [collect(GR.sourcesof(GR.dependencies(plan), c)) for c in part.chunks];
                    init=Int[])))
                @test part.sources == needed
                @test !part.metis_loaded
            end
            @test all(==(1), visits)
            @test isequal(result, expected)
            @test all(isfinite, result)
        end
    finally
        Distributed.rmprocs(processes)
    end
end

end
