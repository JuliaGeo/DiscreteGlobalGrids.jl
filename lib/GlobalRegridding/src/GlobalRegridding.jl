"""
    GlobalRegridding

Regridding between two cell collections on the sphere, for grids too large to
hold in memory.

A participant on either side is a [`RegridSpace`](@ref): a cell-level spatial
tree, a chunk-level spatial tree, lazy cell polygons, and a chunk-to-cell index
map. The package has no notion of what a space *is* — a raster, a discrete
global grid, an unstructured mesh all enter through the same contract — and
ships one implementation, [`RasterGrid`](@ref), for dimensional rasters.

A method ([`Conservative`](@ref), [`NearestCell`](@ref),
[`BilinearPoint`](@ref)) supplies weights and nothing else, through
[`build_weights!`](@ref). Every method here is linear in the source data, so
every method compiles to the same object — a sparse weight matrix with an
optional per-destination denominator ([`WeightBlock`](@ref)) — and the executor
sees only that. Weights are geometry-only: a builder never sees data, IO, or
laziness, which is what makes the same weights serve eager, chunked, and
streaming execution unchanged.

A plan ([`AbstractRegriddingPlan`](@ref)) is self-contained — method, spaces,
missing policy, storage and budget — so applying one takes no keyword
arguments:

    regrid(data; to, method = Conservative())   # build a plan, apply it, drop it
    plan = plan_regrid(data; to, method)        # keep it
    regrid(data, plan)                          # reuse across slices and reads

Missing data is handled once, at finalize: [`Weighted`](@ref) returns a
coverage-normalized mean and blanks destinations below its coverage threshold;
[`Extensive`](@ref) returns raw conservative sums.

# Layout

  - `src/spaces.jl` — the [`RegridSpace`](@ref) contract.
  - `src/rastergrid.jl` — the raster space.
  - `src/methods.jl` — methods, the [`build_weights!`](@ref) hook, missing
    policies.
  - `src/conservative.jl`, `src/interpolation.jl` — weight construction.
  - `src/plans.jl` — [`WeightBlock`](@ref) and the plan types.
  - `src/discovery.jl` — connected chunk pairs from the two chunk trees.
  - `src/executor.jl` — matvec, accumulate, finalize.
  - `src/lazy.jl` — the lazy destination array.
  - `src/api.jl` — [`regrid`](@ref), [`regrid!`](@ref), [`plan_regrid`](@ref).
"""
module GlobalRegridding

import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
import GeometryOps: SpatialTreeInterface as STI

import ConservativeRegridding
import ConservativeRegridding: Trees
# Extend the original `Trees` bindings rather than shadow them: a space is then
# addressable by `ConservativeRegridding` without a wrapper, and a package that
# re-exports both surfaces introduces no binding ambiguity.
import ConservativeRegridding.Trees: ncells, getcell
# The space contract's manifold accessor is GeometryOpsCore's, extended here.
import GeometryOpsCore: manifold

import DimensionalData as DD
import DiskArrays
import LinearAlgebra
import SparseArrays
using SparseArrays: SparseMatrixCSC, sparse

# Unit-sphere aliases. `SphericalCap` stays the UnionAll so its two-argument
# constructor remains available; `Cap` is the concrete node-extent type.
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

# --- The space contract ----------------------------------------------------
export RegridSpace
export celltree, chunktree, nchunks, cellindices, ncells, getcell
export cellcentroid, cellat, hascellchart, manifold

# --- Shipped spaces --------------------------------------------------------
export RasterGrid

# --- Methods and their hook ------------------------------------------------
export AbstractRegriddingMethod
export Conservative, NearestCell, BilinearPoint
export build_weights!, support_radius
export WeightCOO, addweight!, adddenom!

# --- Missing-data policies -------------------------------------------------
export AbstractMissingPolicy, Weighted, Extensive

# --- Plans -----------------------------------------------------------------
export AbstractRegriddingPlan, WeightBlock
export DirectPlan, ChunkedPlan, PerChunk, Spilled

# --- User API --------------------------------------------------------------
export regrid, regrid!, plan_regrid
export LazyRegridArray

end # module GlobalRegridding
