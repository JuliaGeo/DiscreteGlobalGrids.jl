module DiscreteGlobalGridsScotchExt

import DiscreteGlobalGrids as DGG
import Scotch

const Idx = Scotch.SCOTCH_Num
const SCOTCH_LOCK = ReentrantLock()

function DGG.partitionlabels(algorithm::DGG.ScotchPartition,
        problem::DGG.PartitionProblem, nparts::Int, capacities::Vector{Float64})
    nwork = length(problem.ids)
    if nwork == 0 || nparts == 1 || nwork <= nparts
        return DGG.partitionlabels(DGG.WeightedContiguous(), problem, nparts, capacities)
    end
    nwork <= typemax(Idx) || throw(ArgumentError("too many work rows for Scotch $Idx"))
    xadj, adjncy, affinities = DGG._partition_graph(problem, algorithm.maxedges, Idx)
    isempty(adjncy) && return DGG.partitionlabels(
        DGG.WeightedContiguous(), problem, nparts, capacities)
    vwgt = max.(one(Idx), DGG._partition_quantize(problem.weights, Idx; budget=max(nwork, 1_000_000)))
    adjwgt = DGG._partition_quantize(affinities, Idx; budget=max(length(affinities), 1_000_000))
    targets = DGG._partition_quantize(capacities, Idx; budget=max(nparts, 1_000_000))

    return lock(SCOTCH_LOCK) do
        graph = Scotch.graph_build(xadj, adjncy; index_start=0,
            v_weights=vwgt, e_weights=adjwgt)
        context = Scotch.context_alloc()
        bound = nothing
        architecture = nothing
        strategy = nothing
        try
            Scotch.random_clone(context)
            Scotch.random_seed(context, algorithm.seed)
            Scotch.random_reset(context)
            Scotch.context_option!(context, :deterministic, true)
            Scotch.context_option!(context, :fixed_seed, true)
            bound = Scotch.bind_graph(context, graph)
            architecture = Scotch.arch_complete_graph(nparts; weights=targets)
            strategy = Scotch.strat_build(:graph_map; parts=nparts,
                imbalance_ratio=algorithm.imbalance)
            # The bound graph borrows the original graph's native storage.
            labels = GC.@preserve graph Scotch.graph_map(bound, architecture, strategy)
            return Int.(labels) .+ 1
        finally
            strategy === nothing || finalize(strategy)
            architecture === nothing || finalize(architecture)
            bound === nothing || finalize(bound)
            finalize(context)
            finalize(graph)
        end
    end
end

end
