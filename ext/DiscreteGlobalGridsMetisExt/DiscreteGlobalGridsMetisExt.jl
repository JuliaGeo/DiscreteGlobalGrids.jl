module DiscreteGlobalGridsMetisExt

import DiscreteGlobalGrids as DGG
import Metis

const LM = Metis.LibMetis
const Idx = LM.idx_t
const RealT = LM.real_t

# METIS installs process-wide signal traps during a partition call.
const METIS_LOCK = ReentrantLock()

function DGG.partitionlabels(algorithm::DGG.MetisPartition,
        problem::DGG.PartitionProblem, nparts::Int, capacities::Vector{Float64})
    nwork = length(problem.ids)
    if nwork == 0 || nparts == 1 || nwork <= nparts
        return DGG.partitionlabels(DGG.WeightedContiguous(), problem, nparts, capacities)
    end
    nwork <= typemax(Idx) || throw(ArgumentError(
        "$nwork work rows exceed Metis $(Idx) vertex capacity"))

    xadj, adjncy, affinities = DGG._partition_graph(problem, algorithm.maxedges, Idx)
    isempty(adjncy) && return DGG.partitionlabels(
        DGG.WeightedContiguous(), problem, nparts, capacities)

    workweights = all(iszero, problem.weights) ? ones(nwork) : problem.weights
    vwgt = DGG._partition_quantize(workweights, Idx)
    adjwgt = DGG._partition_quantize(affinities, Idx)
    tpwgts = _targetweights(capacities)
    options = Vector{Idx}(undef, Int(LM.METIS_NOPTIONS))
    _check(:METIS_SetDefaultOptions, LM.METIS_SetDefaultOptions(options))
    options[Int(LM.METIS_OPTION_OBJTYPE) + 1] = Idx(LM.METIS_OBJTYPE_CUT)
    options[Int(LM.METIS_OPTION_SEED) + 1] = Idx(algorithm.seed)
    options[Int(LM.METIS_OPTION_UFACTOR) + 1] =
        Idx(max(1, round(Int, 1000 * algorithm.imbalance)))
    options[Int(LM.METIS_OPTION_NUMBERING) + 1] = 0

    nvtxs = Ref(Idx(nwork))
    ncon = Ref(Idx(1))
    count = Ref(Idx(nparts))
    edgecut = Ref(Idx(0))
    labels = Vector{Idx}(undef, nwork)
    status = lock(METIS_LOCK) do
        LM.METIS_PartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, C_NULL,
            adjwgt, count, tpwgts, C_NULL, options, edgecut, labels)
    end
    _check(:METIS_PartGraphKway, status)
    all(label -> 0 <= label < nparts, labels) || throw(ErrorException(
        "METIS returned a partition label outside 0:$(nparts - 1)"))
    return Int.(labels .+ one(Idx))
end

function _targetweights(capacities::Vector{Float64})
    targets = max.(RealT.(capacities), floatmin(RealT))
    targets ./= RealT(sum(Float64, targets))
    largest = argmax(targets)
    targets[largest] += RealT(1 - sum(Float64, targets))
    targets[largest] > 0 || throw(ArgumentError(
        "capacities cannot be represented as positive Metis real_t targets"))
    return targets
end

function _check(functionname::Symbol, status)
    status == LM.METIS_OK || throw(LM.MetisError(functionname, Cint(status)))
    return nothing
end

end
