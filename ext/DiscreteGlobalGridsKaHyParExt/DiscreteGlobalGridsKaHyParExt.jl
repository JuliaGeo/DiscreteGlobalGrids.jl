module DiscreteGlobalGridsKaHyParExt

import DiscreteGlobalGrids as DGG
using KaHyPar_jll: libkahypar

const KAHYPAR_LOCK = ReentrantLock()
const CONFIG = joinpath(@__DIR__, "km1.ini")

function _hypergraph(problem::DGG.PartitionProblem)
    consumers = [Cuint[] for _ in problem.sourceids]
    for row in eachindex(problem.reads), source in problem.reads[row]
        push!(consumers[source], Cuint(row - 1))
    end
    offsets = Csize_t[0]
    pins = Cuint[]
    weights = Float64[]
    for source in eachindex(consumers)
        rows = consumers[source]
        if length(rows) >= 2 && problem.sourceweights[source] > 0
            append!(pins, rows)
            push!(offsets, length(pins))
            push!(weights, problem.sourceweights[source])
        end
    end
    length(pins) <= typemax(Cuint) || throw(ArgumentError(
        "too many source dependencies for KaHyPar"))
    return offsets, pins, weights
end

function DGG.partitionlabels(algorithm::DGG.KaHyParPartition,
        problem::DGG.PartitionProblem, nparts::Int, capacities::Vector{Float64})
    nwork = length(problem.ids)
    if nwork == 0 || nparts == 1 || nwork <= nparts
        return DGG.partitionlabels(DGG.WeightedContiguous(), problem, nparts, capacities)
    end
    nwork <= typemax(Cint) ÷ 8 || throw(ArgumentError("too many work rows for KaHyPar"))
    offsets, pins, sourceweights = _hypergraph(problem)
    isempty(sourceweights) && return DGG.partitionlabels(
        DGG.WeightedContiguous(), problem, nparts, capacities)
    length(sourceweights) <= typemax(Cuint) || throw(ArgumentError(
        "too many shared sources for KaHyPar"))

    # KaHyPar requires positive vertex weights, including for zero-cost work.
    vertexweights = max.(Cint(1), DGG._partition_quantize(problem.weights, Cint;
        budget=max(nwork, 1_000_000)))
    # Connectivity can charge an edge up to k-1 times; its objective is a Cint.
    edgebudget = Int(typemax(Cint)) ÷ 4 ÷ (nparts - 1)
    edgeweights = DGG._partition_quantize(sourceweights, Cint; budget=edgebudget)
    total = sum(Int64, vertexweights)
    targets = [max(1, ceil(Int64, (1 + algorithm.imbalance) * capacity * total))
               for capacity in capacities]
    targets[argmax(capacities)] += max(0, total - sum(targets))
    sum(targets) <= typemax(Cint) || throw(ArgumentError(
        "KaHyPar capacity limits exceed its integer range"))
    blockweights = Cint.(targets)
    labels = fill(Cint(-1), nwork)

    # KaHyPar's C API shares process-global random state.
    lock(KAHYPAR_LOCK) do
        context = ccall((:kahypar_context_new, libkahypar), Ptr{Cvoid}, ())
        context == C_NULL && throw(OutOfMemoryError())
        try
            ccall((:kahypar_configure_context_from_file, libkahypar), Cvoid,
                (Ptr{Cvoid}, Cstring), context, CONFIG)
            ccall((:kahypar_set_seed, libkahypar), Cvoid,
                (Ptr{Cvoid}, Cint), context, algorithm.seed)
            ccall((:kahypar_supress_output, libkahypar), Cvoid,
                (Ptr{Cvoid}, Bool), context, true)
            ccall((:kahypar_set_custom_target_block_weights, libkahypar), Cvoid,
                (Cint, Ptr{Cint}, Ptr{Cvoid}), nparts, blockweights, context)
            valid = ccall((:kahypar_validate_input, libkahypar), Bool,
                (Cuint, Cuint, Ptr{Csize_t}, Ptr{Cuint}, Ptr{Cint}, Ptr{Cint}, Bool),
                nwork, length(edgeweights), offsets, pins, edgeweights, vertexweights, false)
            valid || throw(ArgumentError("KaHyPar rejected the weighted hypergraph"))
            objective = Ref{Cint}(0)
            ccall((:kahypar_partition, libkahypar), Cvoid,
                (Cuint, Cuint, Cdouble, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Csize_t},
                 Ptr{Cuint}, Ref{Cint}, Ptr{Cvoid}, Ptr{Cint}),
                nwork, length(edgeweights), algorithm.imbalance, nparts,
                vertexweights, edgeweights, offsets, pins, objective, context, labels)
        finally
            ccall((:kahypar_context_free, libkahypar), Cvoid, (Ptr{Cvoid},), context)
        end
    end
    return Int.(labels) .+ 1
end

end
