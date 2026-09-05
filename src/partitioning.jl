"""
    AbstractPartitioningAlgorithm

Supertype for chunk partitioning algorithms. Implement
[`partitionlabels`](@ref) to add an algorithm.
"""
abstract type AbstractPartitioningAlgorithm end

"""
    WeightedContiguous()

Cut work into capacity-weighted contiguous segments in the problem's traversal
order.
"""
struct WeightedContiguous <: AbstractPartitioningAlgorithm end

"""
    MetisPartition(; seed=0, imbalance=0.03, maxedges=1_000_000)

Partition work by shared-source affinity through the optional Metis backend.
`maxedges` bounds the projected graph's unique undirected edges.
METIS represents `imbalance` in thousandths and uses `0.001` for a requested
zero.
"""
struct MetisPartition <: AbstractPartitioningAlgorithm
    seed::Cint
    imbalance::Float64
    maxedges::Int

    function MetisPartition(; seed::Integer=0, imbalance::Real=0.03,
            maxedges::Integer=1_000_000)
        seed, u = _partition_seed_imbalance(seed, imbalance)
        m = _partition_maxedges(maxedges)
        return new(seed, u, m)
    end
end

"""
    KaHyParPartition(; seed=0, imbalance=0.03)

Group chunks by shared source data with KaHyPar's connectivity objective.
Each source costs its weight for every additional partition that reads it.
Load `KaHyPar_jll` to enable this optional backend. Capacity targets become native
integer upper bounds, including `imbalance`; indivisible chunks may exceed them.
"""
struct KaHyParPartition <: AbstractPartitioningAlgorithm
    seed::Cint
    imbalance::Float64
    function KaHyParPartition(; seed::Integer=0, imbalance::Real=0.03)
        seed, u = _partition_seed_imbalance(seed, imbalance)
        new(seed, u)
    end
end

"""
    ScotchPartition(; seed=0, imbalance=0.03, maxedges=1_000_000)

Map the shared-source graph onto workers with relative compute capacities.
Load `Scotch` to enable this optional backend. `maxedges` bounds the projected
graph's unique undirected edges.
"""
struct ScotchPartition <: AbstractPartitioningAlgorithm
    seed::Cint
    imbalance::Float64
    maxedges::Int
    function ScotchPartition(; seed::Integer=0, imbalance::Real=0.03,
            maxedges::Integer=1_000_000)
        seed, u = _partition_seed_imbalance(seed, imbalance)
        new(seed, u, _partition_maxedges(maxedges))
    end
end

function _partition_seed_imbalance(seed::Integer, imbalance::Real)
    0 <= seed <= typemax(Cint) || throw(ArgumentError(
        "seed must fit a non-negative Cint, got $seed"))
    u = Float64(imbalance)
    isfinite(u) && 0 <= u <= 1 || throw(ArgumentError(
        "imbalance must be finite and between 0 and 1, got $imbalance"))
    return Cint(seed), u
end

function _partition_maxedges(maxedges::Integer)
    m = try
        Int(maxedges)
    catch
        throw(ArgumentError("maxedges must fit Int, got $maxedges"))
    end
    m > 0 || throw(ArgumentError("maxedges must be positive, got $maxedges"))
    return m
end

"""
    PartitionProblem(reads; weights=nothing, sourceweights=nothing, ids=nothing,
                     sourceids=nothing, order=nothing)

A bipartite chunk-partitioning problem. Each row of `reads` contains source
resource positions used by one indivisible work item.

The fields `reads`, `weights`, `sourceweights`, `ids`, `sourceids`, and `order`
are owned canonical vectors. Treat them as read-only. `order` contains work row
positions, while `ids` and `sourceids` carry stable external identities.
"""
struct PartitionProblem
    reads::Vector{Vector{Int}}
    weights::Vector{Float64}
    sourceweights::Vector{Float64}
    ids::Vector{Int}
    sourceids::Vector{Int}
    order::Vector{Int}
end

function PartitionProblem(reads; weights=nothing, sourceweights=nothing,
        ids=nothing, sourceids=nothing, order=nothing)
    rawrows = collect(reads)
    canonical = Vector{Vector{Int}}(undef, length(rawrows))
    largest = 0
    for (i, row) in pairs(rawrows)
        values = try
            collect(row)
        catch
            throw(ArgumentError("reads row $i must be a collection of resource positions"))
        end
        parsed = Vector{Int}(undef, length(values))
        for (j, value) in pairs(values)
            parsed[j] = _partition_integer(value, "reads[$i][$j]")
            parsed[j] >= 1 || throw(ArgumentError(
                "reads[$i][$j] must be positive, got $(parsed[j])"))
        end
        sort!(parsed)
        unique!(parsed)
        canonical[i] = parsed
        isempty(parsed) || (largest = max(largest, parsed[end]))
    end

    rawsourceweights = sourceweights === nothing ? nothing : collect(sourceweights)
    rawsourceids = sourceids === nothing ? nothing : collect(sourceids)
    nresources = if rawsourceweights !== nothing
        rawsourceids === nothing || length(rawsourceids) == length(rawsourceweights) ||
            throw(ArgumentError("sourceweights and sourceids must have equal lengths"))
        length(rawsourceweights)
    elseif rawsourceids !== nothing
        length(rawsourceids)
    else
        largest
    end
    largest <= nresources || throw(ArgumentError(
        "reads reference resource $largest, but only $nresources resources were supplied"))

    nwork = length(canonical)
    workweights = weights === nothing ? ones(nwork) : _partition_weights(weights, nwork, "weights")
    resourceweights = rawsourceweights === nothing ? ones(nresources) :
        _partition_weights(rawsourceweights, nresources, "sourceweights")
    workids = ids === nothing ? collect(1:nwork) : _partition_ids(ids, nwork, "ids")
    resourceids = rawsourceids === nothing ? collect(1:nresources) :
        _partition_ids(rawsourceids, nresources, "sourceids")
    traversal = order === nothing ? collect(1:nwork) : _partition_order(order, nwork)
    return PartitionProblem(canonical, workweights, resourceweights, workids,
        resourceids, traversal)
end

function _partition_integer(value, name::AbstractString)
    value isa Integer || throw(ArgumentError("$name must be an integer, got $(repr(value))"))
    return try
        Int(value)
    catch
        throw(ArgumentError("$name must fit Int, got $value"))
    end
end

function _partition_weights(values, n::Int, name::AbstractString)
    raw = collect(values)
    length(raw) == n || throw(ArgumentError(
        "$name must contain $n entries, got $(length(raw))"))
    out = Vector{Float64}(undef, n)
    for i in eachindex(raw)
        out[i] = try
            Float64(raw[i])
        catch
            throw(ArgumentError("$name[$i] must convert to Float64, got $(repr(raw[i]))"))
        end
        isfinite(out[i]) && out[i] >= 0 || throw(ArgumentError(
            "$name[$i] must be finite and non-negative, got $(repr(raw[i]))"))
    end
    return out
end

function _partition_ids(values, n::Int, name::AbstractString)
    raw = collect(values)
    length(raw) == n || throw(ArgumentError(
        "$name must contain $n entries, got $(length(raw))"))
    out = Vector{Int}(undef, n)
    seen = Set{Int}()
    for i in eachindex(raw)
        id = _partition_integer(raw[i], "$name[$i]")
        id in seen && throw(ArgumentError("$name must contain unique IDs; $id is repeated"))
        out[i] = id
        push!(seen, id)
    end
    return out
end

function _partition_order(values, n::Int)
    raw = collect(values)
    length(raw) == n || throw(ArgumentError(
        "order must contain $n row positions, got $(length(raw))"))
    out = Vector{Int}(undef, n)
    seen = falses(n)
    for i in eachindex(raw)
        row = _partition_integer(raw[i], "order[$i]")
        1 <= row <= n || throw(ArgumentError(
            "order[$i] must be a row position in 1:$n, got $row"))
        seen[row] && throw(ArgumentError("order must be a permutation; row $row is repeated"))
        out[i] = row
        seen[row] = true
    end
    return out
end

"""
    partitionproblem(input; kwargs...) -> PartitionProblem

Adapt a chunk plan or dependency graph to a partitioning problem. Packages may
add methods for their own work descriptions.
"""
function partitionproblem end

function partitionproblem(plan::MapChunkPlan; weights=nothing, sourceweights=nothing,
        ids=nothing, sourceids=nothing, order=nothing)
    nwork = length(plan)
    ncells = length(plan.lookup)
    if ncells == 0
        nwork == 0 || throw(ArgumentError(
            "a MapChunkPlan over an empty lookup must contain no chunks"))
    elseif nwork == 0
        throw(ArgumentError("MapChunkPlan owns no chunks of its $ncells-cell lookup"))
    end

    positions = sortperm(1:nwork; by=i -> first(plan.chunks[i].range))
    starts = Int[]
    expected = 1
    for position in positions
        chunk = plan.chunks[position]
        isempty(chunk.range) && throw(ArgumentError(
            "MapChunk $(chunk.index) has an empty owned range"))
        first(chunk.range) == expected || throw(ArgumentError(
            "MapChunkPlan owned ranges must form a disjoint complete cover of 1:$ncells"))
        last(chunk.range) <= ncells || throw(ArgumentError(
            "MapChunk $(chunk.index) owns indices beyond its $ncells-cell lookup"))
        push!(starts, first(chunk.range))
        expected = last(chunk.range) + 1
    end
    expected == ncells + 1 || throw(ArgumentError(
        "MapChunkPlan owned ranges must form a disjoint complete cover of 1:$ncells"))

    planids = [chunk.index for chunk in plan.chunks]
    _partition_ids(planids, nwork, "MapChunk indices")
    all(>(0), planids) || throw(ArgumentError("MapChunk indices must be positive"))
    readrows = Vector{Vector{Int}}(undef, nwork)
    defaultweights = Vector{Float64}(undef, nwork)
    defaultsources = Vector{Float64}(undef, nwork)
    for row in 1:nwork
        chunk = plan.chunks[row]
        issorted(chunk.halo) && allunique(chunk.halo) || throw(ArgumentError(
            "MapChunk $(chunk.index) halo must contain ascending unique indices"))
        resources = Int[row]
        for cell in chunk.halo
            1 <= cell <= ncells || throw(ArgumentError(
                "MapChunk $(chunk.index) halo index $cell is outside 1:$ncells"))
            cell in chunk.range && throw(ArgumentError(
                "MapChunk $(chunk.index) halo contains owned index $cell"))
            sortedposition = searchsortedlast(starts, cell)
            sortedposition >= 1 || throw(ArgumentError(
                "MapChunk $(chunk.index) halo index $cell belongs to no owned range"))
            resource = positions[sortedposition]
            cell in plan.chunks[resource].range || throw(ArgumentError(
                "MapChunk $(chunk.index) halo index $cell belongs to no owned range"))
            push!(resources, resource)
        end
        readrows[row] = resources
        defaultweights[row] = length(chunk.range) + length(chunk.halo)
        defaultsources[row] = length(chunk.range)
    end
    return PartitionProblem(readrows;
        weights=weights === nothing ? defaultweights : weights,
        sourceweights=sourceweights === nothing ? defaultsources : sourceweights,
        ids=ids === nothing ? planids : ids,
        sourceids=sourceids === nothing ? planids : sourceids,
        order)
end

function partitionproblem(graph::GR.ChunkDependencyGraph;
        weights=nothing, sourceweights=nothing, ids=nothing, sourceids=nothing,
        order=nothing)
    nwork = GR.ndestinationchunks(graph)
    nresources = GR.nsourcechunks(graph)
    reads = [collect(Int, GR.sourcesof(graph, row)) for row in 1:nwork]
    defaultids = [GR.destinationchunk(graph, row) for row in 1:nwork]
    return PartitionProblem(reads;
        weights=weights === nothing ? ones(nwork) : weights,
        sourceweights=sourceweights === nothing ? ones(nresources) : sourceweights,
        ids=ids === nothing ? defaultids : ids,
        sourceids=sourceids === nothing ? collect(1:nresources) : sourceids,
        order)
end

function partitionproblem(plan::GR.ChunkedPlan;
        weights=nothing, sourceweights=nothing, ids=nothing, sourceids=nothing,
        order=nothing)
    graph = GR.dependencies(plan)
    graph === nothing && throw(ArgumentError(
        "ChunkedPlan has no dependency relation; construct it with dependencies enabled"))
    nwork = GR.ndestinationchunks(graph)
    nresources = GR.nsourcechunks(graph)
    defaultweights = [Float64(length(GR.ownedindices(plan.dst_space, row)))
                      for row in 1:nwork]
    defaultsources = [Float64(length(GR.ownedindices(plan.src_space, source)))
                      for source in 1:nresources]
    return partitionproblem(graph;
        weights=weights === nothing ? defaultweights : weights,
        sourceweights=sourceweights === nothing ? defaultsources : sourceweights,
        ids, sourceids, order)
end

"""
    partitionlabels(algorithm, problem, nparts, capacities) -> Vector{Int}

Return one label in `1:nparts` for each problem row. `capacities` contains
validated normalized targets. Extensions define this hook for custom
`AbstractPartitioningAlgorithm` subtypes.
"""
function partitionlabels end

function partitionlabels(::WeightedContiguous, problem::PartitionProblem,
        nparts::Int, capacities::Vector{Float64})
    nwork = length(problem.ids)
    nwork == 0 && return Int[]
    maximumweight = maximum(problem.weights)
    costs = maximumweight == 0 ? ones(nwork) : problem.weights ./ maximumweight
    ordered = costs[problem.order]
    prefix = vcat(0.0, cumsum(ordered))
    total = prefix[end]
    cuts = Vector{Int}(undef, nparts - 1)
    previous = 0
    target = 0.0
    for part in 1:(nparts - 1)
        target += capacities[part] * total
        upper = searchsortedfirst(prefix, target) - 1
        upper = clamp(upper, previous, nwork)
        lower = max(previous, upper - 1)
        cut = abs(prefix[lower + 1] - target) <= abs(prefix[upper + 1] - target) ?
            lower : upper
        cuts[part] = cut
        previous = cut
    end
    labels = Vector{Int}(undef, nwork)
    firstposition = 1
    for part in 1:nparts
        lastposition = part == nparts ? nwork : cuts[part]
        for position in firstposition:lastposition
            labels[problem.order[position]] = part
        end
        firstposition = lastposition + 1
    end
    return labels
end

"""
    PartitionBackendUnavailable

The requested optional partitioning backend has not been loaded.
"""
struct PartitionBackendUnavailable <: Exception
    backend::Symbol
    algorithm::AbstractPartitioningAlgorithm
end

function Base.showerror(io::IO, err::PartitionBackendUnavailable)
    print(io, "partitioning backend ", err.backend, " is unavailable for ",
        nameof(typeof(err.algorithm)))
    Base.Experimental.show_error_hints(io, err)
end

# The extension's concrete Vector signature is more specific than this fallback.
function partitionlabels(algorithm::MetisPartition, problem::PartitionProblem,
        nparts::Integer, capacities::AbstractVector{<:Real})
    throw(PartitionBackendUnavailable(:Metis, algorithm))
end

function partitionlabels(algorithm::KaHyParPartition, problem::PartitionProblem,
        nparts::Integer, capacities::AbstractVector{<:Real})
    throw(PartitionBackendUnavailable(:KaHyPar_jll, algorithm))
end

function partitionlabels(algorithm::ScotchPartition, problem::PartitionProblem,
        nparts::Integer, capacities::AbstractVector{<:Real})
    throw(PartitionBackendUnavailable(:Scotch, algorithm))
end

"""
    ChunkPartition

A serializable logical chunk assignment. Its vector fields are owned by the
result and exposed read-only through the accessors below.
"""
struct ChunkPartition
    ids::Vector{Int}
    assignment::Vector{Int}
    parts::Vector{Vector{Int}}
    chunks::Vector{Vector{Int}}
    sources::Vector{Vector{Int}}
    weights::Vector{Float64}
end

"""
    partition(problem, nparts; algorithm=WeightedContiguous(), capacities=ones(nparts))
        -> ChunkPartition

Assign indivisible work rows to exactly `nparts` logical partitions. Empty
partitions are retained.
"""
function partition(problem::PartitionProblem, nparts::Integer;
        algorithm::AbstractPartitioningAlgorithm=WeightedContiguous(), capacities=nothing)
    count = _partition_integer(nparts, "nparts")
    count > 0 || throw(ArgumentError("nparts must be positive, got $nparts"))
    rawcapacities = capacities === nothing ? ones(count) : capacities
    targets = _partition_weights(rawcapacities, count, "capacities")
    all(>(0), targets) || throw(ArgumentError("capacities must be positive"))
    scale = maximum(targets)
    targets ./= scale
    targets ./= sum(targets)
    all(>(0), targets) || throw(ArgumentError(
        "capacities are too different in magnitude to normalize to positive Float64 targets"))

    rawlabels = partitionlabels(algorithm, problem, count, targets)
    length(rawlabels) == length(problem.ids) || throw(ArgumentError(
        "$(typeof(algorithm)) returned $(length(rawlabels)) labels for " *
        "$(length(problem.ids)) work rows"))
    labels = Vector{Int}(undef, length(rawlabels))
    for i in eachindex(rawlabels)
        labels[i] = _partition_integer(rawlabels[i], "partition label $i")
        1 <= labels[i] <= count || throw(ArgumentError(
            "$(typeof(algorithm)) returned label $(labels[i]) outside 1:$count at row $i"))
    end

    rows = [Int[] for _ in 1:count]
    chunks = [Int[] for _ in 1:count]
    sourcepositions = [BitSet() for _ in 1:count]
    loads = zeros(count)
    for row in problem.order
        part = labels[row]
        push!(rows[part], row)
        push!(chunks[part], problem.ids[row])
        union!(sourcepositions[part], problem.reads[row])
        load = loads[part] + problem.weights[row]
        isfinite(load) || throw(ArgumentError(
            "work weights in partition $part have a sum outside Float64 range"))
        loads[part] = load
    end
    sources = [[problem.sourceids[position] for position in positions]
               for positions in sourcepositions]
    return ChunkPartition(copy(problem.ids), labels, rows, chunks, sources, loads)
end

function partition(input, nparts::Integer;
        algorithm::AbstractPartitioningAlgorithm=WeightedContiguous(), capacities=nothing,
        kwargs...)
    return partition(partitionproblem(input; kwargs...), nparts; algorithm, capacities)
end

"""Return the exact logical partition count, including empty partitions."""
npartitions(result::ChunkPartition) = length(result.parts)

function _partition_number(result::ChunkPartition, part::Integer)
    p = _partition_integer(part, "partition")
    1 <= p <= npartitions(result) || throw(BoundsError(result, part))
    return p
end

"""Return input row positions in traversal order. Treat the vector as read-only."""
partindices(result::ChunkPartition, part::Integer) =
    result.parts[_partition_number(result, part)]

"""Return stable work IDs in traversal order. Treat the vector as read-only."""
partchunks(result::ChunkPartition, part::Integer) =
    result.chunks[_partition_number(result, part)]

"""Return unique source IDs in source-axis order. Treat the vector as read-only."""
partsources(result::ChunkPartition, part::Integer) =
    result.sources[_partition_number(result, part)]

"""Return original work-weight sums by partition. Treat the vector as read-only."""
partweights(result::ChunkPartition) = result.weights

function Base.getindex(plan::MapChunkPlan, rows::AbstractVector{<:Integer})
    selected = Vector{MapChunk}(undef, length(rows))
    for (i, row) in pairs(rows)
        position = _partition_integer(row, "MapChunkPlan row")
        1 <= position <= length(plan) || throw(BoundsError(plan, row))
        selected[i] = plan.chunks[position]
    end
    return MapChunkPlan(plan.lookup, selected, plan.width, plan.connectivity)
end
