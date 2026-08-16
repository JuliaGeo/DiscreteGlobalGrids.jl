"""
    GlobalRegridding

Regrid spherical cell collections eagerly or in chunks.

Both inputs implement [`RegridSpace`](@ref). Regridding methods build
geometry-only sparse [`WeightBlock`](@ref)s through [`build_weights!`](@ref).
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
import LinearAlgebra
import SparseArrays
using SparseArrays: SparseMatrixCSC, sparse

# Keep `SphericalCap` generic for its two-argument constructor.
const US = GO.UnitSpherical
const USPoint = GO.UnitSphericalPoint{Float64}
const SphericalCap = GO.UnitSpherical.SphericalCap
const Cap = GO.UnitSpherical.SphericalCap{Float64}

include("spaces.jl")
include("rastergrid.jl")
include("methods.jl")
include("conservative.jl")
include("interpolation.jl")
include("plans.jl")
include("discovery.jl")
include("executor.jl")
include("lazy.jl")
include("api.jl")

# Space interface
export RegridSpace
export celltree, chunktree, nchunks, cellindices, ncells, getcell
export cellcentroid, cellat, hascellchart, manifold

# Included spaces
export RasterGrid

# Methods
export AbstractRegriddingMethod
export Conservative, NearestCell, BilinearPoint
export build_weights!, support_radius
export WeightCOO, addweight!, adddenom!

# Missing-data policies
export AbstractMissingPolicy, Weighted, Extensive

# Plans
export AbstractRegriddingPlan, WeightBlock
export DirectPlan, ChunkedPlan, PerChunk, Spilled

# User API
export regrid, regrid!, plan_regrid
export LazyRegridArray

# Qualified extension and observability APIs.
public knownempty, isvalidvalue, sourcemissingval, chunkat
public residency, LazyStats
public weightbudget, databudget, spilledfiles, usesreference

end # module GlobalRegridding
