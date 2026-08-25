"""
    GlobalRegridding

Regrid spherical cell collections eagerly or in chunks.

The source and destination spaces both implement [`RegridSpace`](@ref).
Regridding methods build geometry-only sparse [`WeightBlock`](@ref)s through
[`buildweights!`](@ref).
Plans contain the method, spaces, missing-data policy, storage, and memory
budget, so applying a plan takes no keywords:

    regrid(data; to, method = Conservative())   # build a plan, apply it, drop it
    plan = plan_regrid(data; to, method)        # keep it
    regrid(data, plan)                          # reuse across slices and reads

[`Weighted`](@ref) returns coverage-normalized means; [`Extensive`](@ref)
returns raw sums.
"""
module GlobalRegridding

import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
import GeometryOps: SpatialTreeInterface as STI

import ConservativeRegridding
import ConservativeRegridding: Trees
# Extend the original bindings so spaces work with `ConservativeRegridding`.
import ConservativeRegridding.Trees: ncells, getcell
import GeometryOpsCore: manifold

import DimensionalData as DD
import DiskArrays
using Base.ScopedValues: ScopedValue, @with
import Graphs
import SparseArrays
using SparseArrays: SparseMatrixCSC, sparse

# Keep `SphericalCap` generic for its two-argument constructor.
const US = GO.UnitSpherical
const USPoint = GO.UnitSphericalPoint{Float64}
const SphericalCap = GO.UnitSpherical.SphericalCap
const Cap = GO.UnitSpherical.SphericalCap{Float64}
const Extents = GO.Extents
const FlexibleRTrees = GO.FlexibleRTrees

include("shared.jl")
include("spaces.jl")
include("rastergrid.jl")
include("methods.jl")
include("conservative.jl")
include("intersection_area.jl")
include("interpolation.jl")
include("barycentric.jl")
include("discovery.jl")
# `chunkgraph.jl` precedes `plans.jl`: a `ChunkedPlan` owns the one
# `ChunkDependencyGraph` it exposes, so the plan's field type names the graph's.
include("chunkgraph.jl")
include("plans.jl")
include("executor.jl")
include("lazy.jl")
include("api.jl")

# Space interface
export RegridSpace
export celltree, nchunks, ownedindices, ncells, getcell
export cellcentroid, cellat, hascellchart, manifold

# `cellindices` stays exported for the deprecation shim in `spaces.jl`.
export cellindices

# Included spaces
export RasterGrid

# Methods
export AbstractRegriddingMethod
export Conservative, NearestCell, BilinearPoint, BarycentricPoint
export buildweights!, supportradius

# `build_weights!` and `support_radius` stay exported for the deprecation shims
# in `methods.jl`.
export build_weights!, support_radius

export WeightCOO, addweight!, adddenom!

# Missing-data policies
export AbstractMissingPolicy, Weighted, Extensive

# Plans
export AbstractRegriddingPlan, WeightBlock
export DirectPlan, ChunkedPlan, PerChunk, Spilled

# The one chunk dependency relation a plan owns. Not exported: `dependencies`
# is too generic a name to put in a user's namespace unqualified.
public dependencies

# User API
export regrid, regrid!, plan_regrid
export LazyRegridArray

# Qualified extension and observability APIs.
public knownempty, sourcemissingval, chunkat, cellarea
public residency, LazyStats, ShapedRegridArray
public spilledfiles, usesreference
public outputsampling, destinationdims, dimsource

# Qualified `RegridSpace` extension hooks. Their declarations and contracts are
# grouped by responsibility in spaces.jl; they stay unexported to avoid generic
# names in user namespaces.
public subtree
public chunkextents, chunkextent, chunkindex, candidatechunks!
public chunkranges
public chartaxes, chartcoords, chartlocalindex, chartperiod, chartspacing
public _asspace

# `chartposition` stays public for the deprecation shim in `spaces.jl`.
public chartposition

# Other qualified extension hooks used by package integrations.
public resolvespatialdims
public _prepare_raster_transform_pair, _task_prepared_raster_transform

# The chunk dependency graph. Public but not exported: these names are generic
# enough that exporting them into a user's namespace would be presumptuous.
public ChunkDependencyGraph, chunk_dependency_graph
public sourcesof, consumersof, sourcedegree, consumerdegree
public srcvertex, dstvertex, srcchunk, dstchunk
public issrcvertex, isdstvertex, srcvertices, dstvertices
public nsourcechunks, ndestinationchunks, dependency_radius
# Graph identity and row views: what makes one relation reusable by a plan that
# did not build it, and what a per-column plan restricts it to.
public SpaceStamp, spacestamp, DependencyIdentity, dependency_identity
public narrowphase, UNNAMED_NARROW, validate_dependencies
public restrict, isrestricted, subspace_dependencies
public globaldestinations, globaldestination, localdestination
# The relation's own inputs, kept: where per-chunk cap metadata lives.
public hasextents, destinationextents, sourceextents
public destinationextent, sourceextent

# DiscreteGlobalGrids extends the qualified space contract in every
# responsibility it customizes: `subtree`; `chunkextents`, `chunkindex`, and
# `candidatechunks!`; `chunkranges`; and `dimsource`/`_asspace`. In particular,
# its native DGG and CopernicusDEM chunk paths do not use private discovery
# hooks. The Proj extension separately specializes the two public preparation
# hooks above for task-owned native transforms.

end # module GlobalRegridding
