"""
    GlobalRegridding

Regrid spherical cell collections eagerly or in chunks.

The source and destination are GridInterface grids (`GridInterface.isgrid`):
cells, trees, chunks, point location and neighbourhoods are read through
GridInterface's verbs. Regridding methods build geometry-only sparse
[`WeightBlock`](@ref)s through [`buildweights!`](@ref).
Plans contain the method, grids, missing-data policy, storage, and memory
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
# The grid contract. `ncells` is GridInterface's binding (what
# `ConservativeRegridding` reads); every other verb is spelled qualified.
import GridInterface
import GridInterface: ncells
using GridInterface.Structured: AbstractStructuredGrid, TopDownQuadtreeCursor, TileIndices
import GridInterface.Structured: isstructured
import GridInterface.Structured: nfaces, facesize, facevertex, faceposition, facecell, blockextent
import GridInterface.Structured: hasfaceaxes, faceaxes, facecoords, faceperiod, facelocalindex
import GridInterface.Structured: facecross, facefan
# Dual-cell answers are fixed-capacity, as GridInterface's rings are.
import SmallCollections
using SmallCollections: SmallVector

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
using GeometryOps.FlexibleRTrees: RTree, HPR

include("shared.jl")
include("spaces.jl")
include("locators.jl")
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
# The weightless nearest-cell spike: it specializes `eagerplan` from `api.jl`
# and `_readdestination!` from `lazy.jl`, so it comes after both.
include("directnearest.jl")

# Included grids. `ncells` is GridInterface's, re-exported so a session that
# regrids can count a grid's cells without a second import.
export RasterGrid
export ncells

# Methods
export AbstractRegriddingMethod
export Conservative, NearestCell, BarycentricPoint
export DirectNearest
export buildweights!, supportradius

# Candidate locators
export CandidateLocator, TreeLocator, AnalyticLocator
public overlappairs, locatecell

# `build_weights!` and `support_radius` stay exported for the deprecation shims
# in `methods.jl`.
export build_weights!, support_radius

export WeightCOO, addweight!, adddenom!

# Missing-data policies
export AbstractMissingPolicy, Weighted, Extensive

# Plans
export AbstractRegriddingPlan, WeightBlock
export DirectPlan, ChunkedPlan, PerChunk, Spilled
export NearestDirectPlan

# The one chunk dependency relation a plan owns. Not exported: `dependencies`
# is too generic a name to put in a user's namespace unqualified.
public dependencies

# User API
export regrid, regrid!, plan_regrid
export LazyRegridArray

# Qualified extension and observability APIs.
public knownempty, sourcemissingval, cellarea
# The destination's nodata convention: what a regrid blanks with and declares.
public outputmissingval, destinationmissingval, rebuildoutput
public residency, LazyStats, ShapedRegridArray
public spilledfiles, usesreference
public outputsampling, destinationdims, dimsource

# The execution hooks a package extends for its grids. Their contracts are in
# `spaces.jl`, `methods.jl` and `barycentric.jl`; they stay unexported to avoid
# generic names in user namespaces. Everything grid-shaped (cells, trees,
# chunks, location, neighbours, charts) is read through GridInterface's verbs.
public chunkranges, samplerstate
# The source-chunk index a dependency relation is discovered through, and its
# one query; a grid whose chunks a packed R-tree indexes poorly may answer both.
public chunkindex, candidatechunks!
# Making a source grid point-samplable with a state of its own: `samplerstate`
# prepares whatever the lookup reads, `weightsat!` answers one point on the
# `Sampler` carrying it, and `dualcellat` names the `DualCell` its nodes and
# `BasisKind` describe. `chartat` is where a source without a cell chart says
# which plane its nodes are written in.
public weightsat!, dualcellat, Sampler, DualCell, BasisKind, Bilinear, MeanValue, chartat
public DualState, ChartState, WeightRow, meanvalueweights!, addentry!, snapunit
# What `weightsat!` answers.
public WeightStatus, WeightsMapped, WeightsRim, WeightsOutside, WeightsDegenerate, ismapped
public _asgrid

# Other qualified extension hooks used by package integrations.
public resolvespatialdims
public _prepare_raster_transform_pair, _task_prepared_raster_transform
public _crs_native_to_unit_sphere, _projected_crs_native_to_unit_sphere

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
public destinationchunks, destinationchunk, destinationrow
# The relation's own inputs, kept: where per-chunk cap metadata lives.
public hasextents, destinationextents, sourceextents
public destinationextent, sourceextent

# DiscreteGlobalGrids answers the grid verbs for its grids (`treeify`,
# `subtree`, `chunks`, `cellat`, `neighbors`, ...) and extends only the
# execution hooks above (`samplerstate`, `weightsat!`, `supportradius`,
# `destinationdims`, `dimsource`, `_asgrid`). The Proj extension separately
# specializes the two public preparation hooks above for task-owned native
# transforms.

end # module GlobalRegridding
