function _partition_graph(problem::PartitionProblem, maxedges::Int, ::Type{Idx}) where {Idx<:Integer}
    nwork = length(problem.ids)
    consumers = [Int[] for _ in problem.sourceids]
    for row in 1:nwork, source in problem.reads[row]
        push!(consumers[source], row)
    end

    contributing = [problem.sourceweights[source] for source in eachindex(consumers)
                    if length(consumers[source]) > 1 && problem.sourceweights[source] > 0]
    maxsource = isempty(contributing) ? 0.0 : maximum(contributing)
    edges = Dict{Tuple{Int,Int},Float64}()
    maxnativeedges = Int(typemax(Idx)) ÷ 2
    if maxsource > 0
        for source in eachindex(consumers)
            sourceweight = problem.sourceweights[source]
            sourceweight == 0 && continue
            rows = consumers[source]
            degree = length(rows)
            degree <= 1 && continue
            if isempty(edges)
                cliqueedges = degree * (degree - 1) ÷ 2
                cliqueedges > maxedges && _partition_edge_budget_error(maxedges)
                cliqueedges <= maxnativeedges || throw(ArgumentError(
                    "projected graph exceeds $(Idx) CSR capacity"))
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
                    length(edges) < maxedges || _partition_edge_budget_error(maxedges)
                    length(edges) < maxnativeedges || throw(ArgumentError(
                        "projected graph exceeds $(Idx) CSR capacity"))
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

function _partition_edge_budget_error(maxedges::Int)
    throw(ArgumentError(
        "projected graph exceeds maxedges=$maxedges; use WeightedContiguous " *
        "or KaHyParPartition for dense resource sharing"))
end

# One unit preserves each positive weight; the remaining bounded budget carries ratios.
function _partition_quantize(values::AbstractVector{Float64}, ::Type{Idx};
        budget::Int=min(Int(typemax(Idx)) ÷ 4, 1_000_000_000)) where {Idx<:Integer}
    1 <= budget <= typemax(Idx) ÷ 4 || throw(ArgumentError(
        "weight budget exceeds the native integer range"))
    result = zeros(Idx, length(values))
    positive = findall(>(0), values)
    isempty(positive) && return result
    for i in positive
        result[i] = one(Idx)
    end
    maxvalue = maximum(values)
    scaledtotal = sum(value / maxvalue for value in values)
    length(positive) <= budget || throw(ArgumentError("too many positive weights for the native integer budget"))
    target = budget
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
