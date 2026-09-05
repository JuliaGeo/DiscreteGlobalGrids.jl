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

    xadj, adjncy, affinities = _project(problem, algorithm.maxedges)
    isempty(adjncy) && return DGG.partitionlabels(
        DGG.WeightedContiguous(), problem, nparts, capacities)

    workweights = all(iszero, problem.weights) ? ones(nwork) : problem.weights
    vwgt = _quantize(workweights)
    adjwgt = _quantize(affinities)
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

function _project(problem::DGG.PartitionProblem, maxedges::Int)
    nwork = length(problem.ids)
    consumers = [Int[] for _ in problem.sourceids]
    for row in 1:nwork, source in problem.reads[row]
        push!(consumers[source], row)
    end

    contributing = [problem.sourceweights[source] for source in eachindex(consumers)
                    if length(consumers[source]) > 1 && problem.sourceweights[source] > 0]
    maxsource = isempty(contributing) ? 0.0 : maximum(contributing)
    edges = Dict{Tuple{Int,Int},Float64}()
    maxmetisedges = Int(typemax(Idx)) ÷ 2
    if maxsource > 0
        for source in eachindex(consumers)
            sourceweight = problem.sourceweights[source]
            sourceweight == 0 && continue
            rows = consumers[source]
            degree = length(rows)
            degree <= 1 && continue
            if isempty(edges)
                cliqueedges = degree * (degree - 1) ÷ 2
                cliqueedges > maxedges && _edge_budget_error(maxedges)
                cliqueedges <= maxmetisedges || throw(ArgumentError(
                    "projected graph exceeds Metis $(Idx) CSR capacity"))
            end
            contribution = max((sourceweight / maxsource) / (degree - 1),
                nextfloat(0.0))
            for j in 2:degree, i in 1:(j - 1)
                edge = (rows[i], rows[j])
                if haskey(edges, edge)
                    combined = edges[edge] + contribution
                    isfinite(combined) || throw(ArgumentError(
                        "projected affinity weights exceed Float64 range"))
                    edges[edge] = combined
                else
                    length(edges) < maxedges || _edge_budget_error(maxedges)
                    length(edges) < maxmetisedges || throw(ArgumentError(
                        "projected graph exceeds Metis $(Idx) CSR capacity"))
                    edges[edge] = contribution
                end
            end
        end
    end

    adjacency = [Tuple{Int,Float64}[] for _ in 1:nwork]
    for ((left, right), weight) in edges
        push!(adjacency[left], (right, weight))
        push!(adjacency[right], (left, weight))
    end
    xadj = Vector{Idx}(undef, nwork + 1)
    xadj[1] = 0
    adjncy = Vector{Idx}(undef, 2 * length(edges))
    affinities = Vector{Float64}(undef, length(adjncy))
    cursor = 1
    for row in 1:nwork
        sort!(adjacency[row]; by=first)
        for (neighbor, weight) in adjacency[row]
            adjncy[cursor] = Idx(neighbor - 1)
            affinities[cursor] = weight
            cursor += 1
        end
        xadj[row + 1] = Idx(cursor - 1)
    end
    return xadj, adjncy, affinities
end

function _edge_budget_error(maxedges::Int)
    throw(ArgumentError(
        "projected graph exceeds maxedges=$maxedges; use WeightedContiguous " *
        "or a hypergraph backend suited to dense resource sharing"))
end

# One unit preserves each positive weight; the remaining bounded budget carries ratios.
function _quantize(values::AbstractVector{Float64})
    result = zeros(Idx, length(values))
    positive = findall(>(0), values)
    isempty(positive) && return result
    for i in positive
        result[i] = one(Idx)
    end
    maxvalue = maximum(values)
    scaledtotal = sum(value / maxvalue for value in values)
    integermax = Int(typemax(Idx))
    target = max(length(positive), min(integermax ÷ 4, 1_000_000_000))
    remaining = target - length(positive)
    used = 0
    for i in positive
        increment = floor(Int, (values[i] / maxvalue) / scaledtotal * remaining)
        increment = min(increment, remaining - used)
        result[i] += Idx(increment)
        used += increment
    end
    return result
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
