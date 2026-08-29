"""
    GlobalRegridding

Regrid spherical cell data eagerly or in chunks.

[`RegridSpace`](@ref) supplies source and destination geometry. Regridding
methods build geometry-only [`WeightBlock`](@ref)s through
[`buildweights!`](@ref). Plans collect the method, spaces, missing-data policy,
storage, and memory budget for reuse:

    regrid(data; to, method = Conservative())
    plan = plan_regrid(data; to, method)
    regrid(data, plan)

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
using StableTasks: StableTasks, StableTask

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
# `ChunkedPlan` needs `ChunkDependencyGraph` defined first.
include("chunkgraph.jl")
include("plans.jl")
include("executor.jl")
include("lazy.jl")
include("api.jl")
# Its specializations require `api.jl` and `lazy.jl` bindings.
include("directnearest.jl")

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
export Conservative, NearestCell, BarycentricPoint
export DirectNearest
export buildweights!, supportradius

# Deprecation shims require the legacy names to remain exported.
export build_weights!, support_radius

export WeightCOO, addweight!, adddenom!

# Missing-data policies
export AbstractMissingPolicy, Weighted, Extensive

# Plans
export AbstractRegriddingPlan, WeightBlock
export DirectPlan, ChunkedPlan, PerChunk, Spilled
export NearestDirectPlan

# Qualified access keeps the generic name out of user namespaces.
public dependencies

# User API
export regrid, regrid!, plan_regrid
export LazyRegridArray

# Qualified extension and observability APIs.
public knownempty, sourcemissingval, chunkat, cellarea
public residency, LazyStats, ShapedRegridArray
public spilledfiles, usesreference
public outputsampling, destinationdims, dimsource
# Source-presentation hooks and the method traits that select a presentation.
public sourceview, sourcespacefor, checksource
public sourcesampling, refinementinvariant

# Qualified access keeps generic extension-hook names out of user namespaces.
public subtree, expensivecellgeometry
public chunkextents, chunkextent, chunkindex, candidatechunks!
# Point sampling separates topology preparation from per-point lookup.
public hasdualcells, dualcellat, samplerstate
public Sampler, DualCell, BasisKind, Bilinear, MeanValue, chartat
public chunkranges
public chartaxes, chartcoords, chartlocalindex, chartperiod, chartspacing
public _asspace

# `chartposition` stays public for the deprecation shim in `spaces.jl`.
public chartposition

public resolvespatialdims
public _prepare_raster_transform_pair, _task_prepared_raster_transform

# Qualified access keeps generic graph names out of user namespaces.
public ChunkDependencyGraph, chunk_dependency_graph
public sourcesof, consumersof, sourcedegree, consumerdegree
public srcvertex, dstvertex, srcchunk, dstchunk
public issrcvertex, isdstvertex, srcvertices, dstvertices
public nsourcechunks, ndestinationchunks, dependency_radius
# Identity validates relation reuse; row views support per-column plans.
public SpaceStamp, spacestamp, DependencyIdentity, dependency_identity
public narrowphase, UNNAMED_NARROW, validate_dependencies
public restrict, isrestricted, subspace_dependencies
public destinationchunks, destinationchunk, destinationrow
# Relation-owned caps preserve the geometry used at construction.
public hasextents, destinationextents, sourceextents
public destinationextent, sourceextent

end # module GlobalRegridding
