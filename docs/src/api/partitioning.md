# Assigning chunks to workers

```@meta
CurrentModule = DiscreteGlobalGrids
```

Partitioning assigns whole chunks to logical groups that an executor can place
on threads or processes. Use it to balance work and, with a suitable algorithm,
keep computations that read the same source data together.

The same API accepts a neighbourhood [`MapChunkPlan`](@ref), a
`GlobalRegridding.ChunkedPlan`, or its destination-to-source dependency graph.
Each result records the assigned chunks and the sources each partition needs.
The result contains ordinary numeric arrays and can be serialized without an
open store, a regridding plan, or the partitioning library.

## Describe the work

A [`PartitionProblem`](@ref) has one row per work chunk. Each row lists the
source resources that chunk reads. Work weights estimate computation; source
weights estimate the cost of reading or transferring shared data.

```@example partitioning
import DiscreteGlobalGrids as DGG

problem = DGG.PartitionProblem([[1], [1], [2], [2]];
    ids=[11, 12, 21, 22], sourceids=[101, 102],
    weights=[1.0, 1.0, 1.0, 1.0])
assignment = DGG.partition(problem, 2)

[(chunks=DGG.partchunks(assignment, p),
  sources=DGG.partsources(assignment, p)) for p in 1:2]
```

Chunks 11 and 12 read source 101; chunks 21 and 22 read source 102. The
default [`WeightedContiguous`](@ref) algorithm cuts the supplied traversal
order into groups with approximately equal work.

`capacities=[1, 2]` asks the second partition to carry twice the work of the
first. These are relative compute capacities. Chunks remain indivisible, so
their sizes limit the balance any algorithm can achieve. Capacity targets are
not memory limits.

```@example partitioning
unequal = DGG.partition(problem, 2; capacities=[1, 2])
DGG.partweights(unequal)
```

The API keeps row positions separate from application IDs:

| value | meaning |
|---|---|
| `order` | a permutation of input row positions |
| `partindices(assignment, p)` | assigned input rows, in traversal order |
| `partchunks(assignment, p)` | the corresponding application chunk IDs |
| `partsources(assignment, p)` | unique application source IDs needed by the partition |

The requested number of logical partitions is preserved. Some can be empty
when there is little work. Partition numbers are independent of process IDs;
an executor chooses where each partition runs.

## Partition a neighbourhood sweep

Build the full chunk plan, then partition it. [`partitionproblem`](@ref)
identifies which storage chunks supply each chunk's owned cells and halo.
Its default work estimate counts the owned and halo cells; override `weights`
when measured kernel costs or other array dimensions matter.

```julia
plan = DGG.chunkplan(A; halo=1)
assignment = DGG.partition(plan, Threads.nthreads())
out = zeros(Float64, size(A))

@sync for p in 1:DGG.npartitions(assignment)
    piece = plan[DGG.partindices(assignment, p)]
    Threads.@spawn DGG.mapneighbors!(out, kernel, A, piece; threaded=false)
end
```

Each piece writes its owned indices. For a stored destination, also align work
ownership with its physical write chunks so concurrent writes do not update
the same storage block.

The current chunk runner reads each chunk's halo separately. Grouping related
chunks creates an opportunity for reuse; an executor needs a cache or halo
exchange to realize it. See [chunk sweeps](chunk-sweep.md) for the loading and
kernel contracts.

## Partition a regridding run

A regridding partition assigns **destination chunks**. Source chunks remain
read dependencies: several partitions may need the same source.

```julia
import GlobalRegridding as GR

assignment = DGG.partition(regridplan, 4)
destination_chunks = DGG.partchunks(assignment, 1)
source_chunks = DGG.partsources(assignment, 1)
```

This reads `GR.dependencies(regridplan)`, the relation the plan already owns.
It builds neither regridding weights nor a replacement dependency graph.
The relation is conservative: a listed source may be ruled out when the
actual interpolation or overlap weights are built.

For a run whose graph rows correspond to application-specific chunk IDs,
provide the mapping explicitly. For example, a Copernicus DEM run can use its
destination store chunk numbers and source tile numbers:

```julia
problem = DGG.partitionproblem(dag.graph;
    ids=todochunks, sourceids=tiles, order=dag.order,
    weights=estimated_work, sourceweights=tile_bytes)
assignment = DGG.partition(problem, length(worker_ids);
    capacities=worker_capacities)

jobs = [(chunks=DGG.partchunks(assignment, p),
         tiles=DGG.partsources(assignment, p))
        for p in 1:DGG.npartitions(assignment)]
```

Here `estimated_work` follows graph row order, `tile_bytes` follows source
order, and `dag.order` is the existing traversal permutation. On a restricted
dependency graph, default destination IDs retain their original chunk numbers;
`partindices` still refers to the restricted graph's local rows.

Send each job and the input configuration to its worker. The worker opens its
own data handles and schedules the assigned destination chunks. Tasks within
that worker can pull work dynamically while sharing a tile cache. Recompute
cache consumer counts for the worker's assigned destinations; global consumer
counts include work that other workers will finish.

For resumed runs, describe only the pending destinations and retain their
application IDs. An assignment describes ownership; the executor still owns
task completion, retries, data exchange, and output writes.

## Use METIS

Load Metis.jl to enable [`MetisPartition`](@ref):

```julia
using Metis

assignment = DGG.partition(problem, 4;
    algorithm=DGG.MetisPartition(seed=42, imbalance=0.03))
```

The algorithm type is available from DiscreteGlobalGrids even before Metis.jl
is loaded. Calling it then raises [`PartitionBackendUnavailable`](@ref) with a
hint to load `Metis`. Only the process building the assignment needs Metis.jl.

The extension builds a graph whose vertices are work chunks. Two chunks are
connected when they read a common source. A source with `d` consumers adds
`sourceweight / (d - 1)` to each pair's connection weight. METIS balances work
weights while reducing connections that cross partitions.

This graph approximates shared-data costs. It does not count exact replicated
bytes or enforce a cache budget. A source with many consumers creates many
pairwise connections; `maxedges` bounds the graph expansion. Native weights
are quantized to METIS integers. Use consistent inputs and a fixed seed for
repeatable calls; partition labels are not a persistent identity across
library versions.

## Add a partitioner

Implement [`partitionproblem`](@ref) to adapt another kind of plan to work rows
and source dependencies. The problem retains the full dependency relation so
an algorithm can use it directly or build its own graph representation.

For another algorithm, subtype [`AbstractPartitioningAlgorithm`](@ref) and
implement [`partitionlabels`](@ref):

```julia
function DGG.partitionlabels(algorithm::MyPartitioner,
        problem::DGG.PartitionProblem, nparts::Int,
        capacities::Vector{Float64})
    # Return one logical partition number per original problem row.
end
```

The caller supplies validated, normalized capacities. Return integer labels
in `1:nparts`, aligned with the original rows even when `problem.order` is a
permutation. The common API validates those labels and constructs the result's
ordered work lists and unique source lists.

```@docs
PartitionProblem
partitionproblem
partition
ChunkPartition
npartitions
partindices
partchunks
partsources
partweights
AbstractPartitioningAlgorithm
WeightedContiguous
MetisPartition
partitionlabels
PartitionBackendUnavailable
```
