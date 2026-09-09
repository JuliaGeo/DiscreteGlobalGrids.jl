module PartitionWorkerFixtures

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GlobalRegridding as GR

partitioner_loaded() = any(name -> Base.get_extension(DGG, name) !== nothing,
    (:DiscreteGlobalGridsMetisExt, :DiscreteGlobalGridsKaHyParExt, :DiscreteGlobalGridsScotchExt))

function stencil_fixture()
    grid = DGG.levelgrid(DGG.S2System(), 2)
    lookup = DGG.CellLookup(DGG.CellVector(grid))
    data = DD.DimArray(Float64.(1:DGG.ncells(grid)), DGG.Cells(lookup))
    return data, DGG.chunkplan(data; chunks=8)
end

stencil(cell, value, neighbors) = value + sum(neighbors; init=0.0)

function stencil_part(assignment, part)
    data, plan = stencil_fixture()
    result = fill(NaN, size(data))
    selected = plan[DGG.partindices(assignment, part)]
    DGG.mapneighbors!(result, stencil, data, selected; threaded=false)
    indices = reduce(vcat, [collect(DGG.ownedindices(c)) for c in selected]; init=Int[])
    return (; indices, values=result[indices], chunks=DGG.partchunks(assignment, part),
        partitioner_loaded=partitioner_loaded())
end

function regrid_fixture()
    system = DGG.S2System()
    source = DGG.DGGSpace(DGG.levelgrid(system, 1); chunklevel=0)
    destination = DGG.DGGSpace(DGG.levelgrid(system, 2); chunklevel=1)
    data = Float64.(1:GR.ncells(source))
    plan = GR.ChunkedPlan(GR.NearestCell(), GR.Weighted(1.0), destination, source;
        dependencies=true)
    return data, plan
end

function regrid_part(assignment, part)
    data, plan = regrid_fixture()
    output = GR.regrid(data, plan)
    chunks = DGG.partchunks(assignment, part)
    indices = reduce(vcat,
        [collect(GR.ownedindices(plan.dst_space, c)) for c in chunks]; init=Int[])
    values = Float64[]
    for chunk in chunks
        append!(values, output[GR.ownedindices(plan.dst_space, chunk)])
    end
    return (; indices, values, chunks,
        sources=DGG.partsources(assignment, part),
        partitioner_loaded=partitioner_loaded())
end

end
