# Regridding methods build geometry-only weights for a pair of chunks.

"""
    AbstractRegriddingMethod

Base type for algorithms that combine source-cell values into destination-cell
values.

Methods implement [`buildweights!`](@ref). Stencils that extend beyond
overlapping cells also implement [`supportradius`](@ref). Every method must
produce weights that are linear and independent of field data.
"""
abstract type AbstractRegriddingMethod end

"""
    refinementinvariant(method) -> Bool

Return whether refining each source cell into children with replicated values
preserves the method's result.

Area methods preserve the covered value and total overlap. Nearest-cell methods
preserve the containing cell's value. Methods that blend by site position
generally change their result because refinement moves sample sites. The default
is `false`; [`sourceview`](@ref) substitutes refined geometry only for methods
that opt in.
"""
refinementinvariant(::AbstractRegriddingMethod) = false

"""
    Conservative()

Weight source cells by their spherical intersection area with each destination.
With [`Extensive`](@ref), this preserves the covered integral; with
[`Weighted`](@ref), it returns coverage-normalized means. Requires cell polygons.
"""
struct Conservative <: AbstractRegriddingMethod end

refinementinvariant(::Conservative) = true

"""
    NearestCell()

Give weight 1 to the source cell containing each destination centroid.

Requires [`cellcentroid`](@ref) for the destination and [`cellat`](@ref) for the
source. A centroid outside source coverage emits an empty row for the missing
policy to resolve. Use [`Conservative`](@ref) when integral preservation is
required.
"""
struct NearestCell <: AbstractRegriddingMethod end

refinementinvariant(::NearestCell) = true

"""
    BarycentricPoint(; poles = NearestCell())

Interpolate between source sample sites at each destination sample site.

The stencil uses the source dual cell containing the destination site. Its basis
selects tensor Q1 coordinates on quadrilaterals and mean-value coordinates on
convex polygons; the triangular case gives barycentric coordinates. The
nonnegative, unit-sum weights keep results within the source-value range.

The method requires [`cellcentroid`](@ref) for the destination and point queries
on the source. Sites outside the source dual complex emit an empty row.

`poles` handles sources whose sample rows end before a pole:

  - [`NearestCell`](@ref) uses the nearest site in the polemost row.
  - `nothing` leaves the site unmapped.

Sources whose sample sites reach both poles ignore this policy.
"""
struct BarycentricPoint{P} <: AbstractRegriddingMethod
    poles::P
    function BarycentricPoint(poles::P) where {P}
        (poles isa NearestCell || poles === nothing) || throw(ArgumentError(
            "BarycentricPoint(poles = $(repr(poles))) is not a polar policy; " *
            "pass NearestCell() to take the nearest polemost sample site, or " *
            "nothing to leave points beyond the polemost row unmapped"))
        return new{P}(poles)
    end
end

BarycentricPoint(; poles = NearestCell()) = BarycentricPoint(poles)

# The shortest call that reconstructs it, so the opt-out is what stands out.
Base.show(io::IO, m::BarycentricPoint) = print(io, "BarycentricPoint(",
    m.poles isa NearestCell ? "" : "poles = $(repr(m.poles))", ")")

"""
    outputsampling(method::AbstractRegriddingMethod) -> DimensionalData.Lookups.Sampling

Return the sampling metadata written on destination axes. Area-based methods
use `Intervals(Center())`, the default; point methods use `Points()`.

This trait also selects weight construction. A point method with a
[`sampler`](@ref) builds whole destination tiles, while [`weightblock`](@ref)
may specialize pair assembly by sampling type.
"""
outputsampling(::AbstractRegriddingMethod) = DD.Lookups.Intervals(DD.Lookups.Center())
outputsampling(::NearestCell) = DD.Lookups.Points()
outputsampling(::BarycentricPoint) = DD.Lookups.Points()

"""
    sourcesampling(method::AbstractRegriddingMethod) -> DimensionalData.Lookups.Sampling

Return the sampling a method reads from its source. `Intervals(Center())`, the
default, reads source cell areas and requires a gap-free polygon cover.
`Points()` reads source sample sites.

This trait is independent of [`outputsampling`](@ref): a method may read areas
and write points. Sources with multiple presentations use it through
[`sourcespacefor`](@ref) and [`sourceview`](@ref).
"""
sourcesampling(::AbstractRegriddingMethod) = DD.Lookups.Intervals(DD.Lookups.Center())
sourcesampling(::NearestCell) = DD.Lookups.Points()
sourcesampling(::BarycentricPoint) = DD.Lookups.Points()

# Weight construction

"""
    WeightCOO(ndst::Int)

Create a chunk-local coordinate-list accumulator for `ndst` destination cells.
`rows` and `cols` index the builder's `dst_inds` and `src_inds`; block assembly
sums duplicate entries.

`denom` stores optional destination denominators. [`markdenominated!`](@ref) or
the first [`adddenom!`](@ref) allocates the vector. Point methods leave it as
`nothing`.
"""
mutable struct WeightCOO
    const ndst::Int
    const rows::Vector{Int}
    const cols::Vector{Int}
    const vals::Vector{Float64}
    denom::Union{Nothing,Vector{Float64}}
end

WeightCOO(ndst::Integer) = WeightCOO(Int(ndst), Int[], Int[], Float64[], nothing)

Base.length(coo::WeightCOO) = length(coo.vals)

Base.show(io::IO, coo::WeightCOO) =
    print(io, "WeightCOO(ndst=", coo.ndst, ", entries=", length(coo),
        coo.denom === nothing ? "" : ", denom", ")")

"""
    addweight!(coo::WeightCOO, dst_local::Int, src_local::Int, w::Real)

Add `w` to the weight of source `src_local` in destination `dst_local`. Both
are chunk-local indices within the builder's `dst_inds` and `src_inds`, not the
spaces' local indices.
"""
function addweight!(coo::WeightCOO, dst_local::Int, src_local::Int, w::Real)
    push!(coo.rows, dst_local)
    push!(coo.cols, src_local)
    push!(coo.vals, Float64(w))
    return coo
end

"""
    adddenom!(coo::WeightCOO, dst_local::Int, d::Real)

Add `d` to the denominator of chunk-local destination `dst_local`. Report only
the share from the current source chunk.
"""
function adddenom!(coo::WeightCOO, dst_local::Int, d::Real)
    v = coo.denom
    if v === nothing
        v = zeros(Float64, coo.ndst)
        coo.denom = v
    end
    v[dst_local] += Float64(d)
    return coo
end

"""
    markdenominated!(coo::WeightCOO)

Declare denominator support and return `coo`. This allocates a zero-filled
vector, allowing a builder with zero destination coverage to produce a
denominated block of zeros.
"""
function markdenominated!(coo::WeightCOO)
    coo.denom === nothing && (coo.denom = zeros(Float64, coo.ndst))
    return coo
end

"""
    WeightRow()

Create one reusable row of source-cell indices and weights for point sampling.

The parallel vectors assign `weights[k]` to source cell `indices[k]`.
[`weightsat!`](@ref) clears the row before reuse, so a sweep grows it only to
its widest stencil. Each task must own its row.
"""
struct WeightRow
    indices::Vector{Int}
    weights::Vector{Float64}
end

WeightRow() = WeightRow(Int[], Float64[])

Base.length(row::WeightRow) = length(row.indices)
Base.isempty(row::WeightRow) = isempty(row.indices)

function Base.empty!(row::WeightRow)
    empty!(row.indices)
    empty!(row.weights)
    return row
end

Base.show(io::IO, row::WeightRow) = print(io, "WeightRow(", length(row), " entries)")

# Omitting zero weights keeps boundary stencils limited to contributing nodes.
@inline function _addentry!(row::WeightRow, i::Integer, w::Float64)
    w > 0 || return row
    push!(row.indices, Int(i))
    push!(row.weights, w)
    return row
end

"""
    WeightStatus

Describe the result of [`weightsat!`](@ref) for one destination point.
`WeightsMapped` carries a complete row; diagnostic states leave the row empty.
Execution distinguishes them through [`ismapped`](@ref).

Diagnostic states:

  - `WeightsOutside`: the point lies outside the supplied cell, dual complex, or
    source coverage;
  - `WeightsRim`: the collection lacks a sample site required by the dual cell;
  - `WeightsDegenerate`: repeated nodes, zero area, a fold, a reflex corner, or
    inverse-map failure makes the cell unusable.
"""
@enum WeightStatus::UInt8 begin
    WeightsMapped
    WeightsOutside
    WeightsRim
    WeightsDegenerate
end

"""
    ismapped(status::WeightStatus) -> Bool

Return whether `status` says the row carries a stencil.
"""
@inline ismapped(status::WeightStatus) = status === WeightsMapped

"""
    buildweights!(coo, method, dst_space, dst_inds, src_space, src_inds)

Append chunk-local weights for `dst_inds` and `src_inds`, then return `coo`.

Builder requirements:

  - Emit weights only for cells in `src_inds`; geometry inspection may extend
    beyond that set.
  - Keep weights independent of field data and execution order.

This required method hook is also the generic fallback for [`pairblock`](@ref).
Wrapper methods must forward `pairblock` and, for point methods,
[`sampler`](@ref) to preserve specialized construction.
"""
function buildweights!(coo::WeightCOO, method::AbstractRegriddingMethod,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    throw(ArgumentError(
        "buildweights! is not implemented for $(typeof(method)) from " *
        "$(typeof(src_space)) to $(typeof(dst_space))"))
end

"""
    supportradius(method, src_space::RegridSpace) -> Float64

Return a conservative upper bound, in radians, on stencil reach beyond a source
chunk. The default is `0.0`. Overestimates add discovery work; underestimates
can omit required weights.

Point stencils may use sample sites whose cells never overlap the destination.
Chunk discovery expands source caps by this bound so the dependency relation
contains every chunk a tile may read.
"""
supportradius(::AbstractRegriddingMethod, ::RegridSpace) = 0.0

# Deprecations forward calls only; extensions must implement the current hooks.

"""
    build_weights!(coo, method, dst_space, dst_inds, src_space, src_inds)

Deprecated. Use [`buildweights!`](@ref); this method forwards directly.
"""
function build_weights! end

@deprecate build_weights!(coo::WeightCOO, method::AbstractRegriddingMethod,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds) buildweights!(
    coo, method, dst_space, dst_inds, src_space, src_inds) false

"""
    support_radius(method, src_space::RegridSpace) -> Float64

Deprecated. Use [`supportradius`](@ref); this method forwards directly.
"""
function support_radius end

@deprecate support_radius(method::AbstractRegriddingMethod, src_space::RegridSpace) supportradius(
    method, src_space) false

# Missing-data policies

"""
    AbstractMissingPolicy

Define how accumulated destinations handle missing source values.
"""
abstract type AbstractMissingPolicy end

"""
    Weighted(threshold::Real = 0.5)

Return each weighted sum divided by its valid source weight. Destinations below
the valid-coverage `threshold` become missing. `threshold` must lie in `[0, 1]`.
The denominator is the source-covered portion of each destination.
"""
struct Weighted <: AbstractMissingPolicy
    threshold::Float64
    function Weighted(threshold::Real = 0.5)
        0 <= threshold <= 1 ||
            throw(ArgumentError("Weighted threshold must lie in [0, 1], got $threshold"))
        return new(Float64(threshold))
    end
end

"""
    Extensive()

Return undivided weighted sums. With [`Conservative`](@ref), this preserves the
integral over valid source cells; missing cells contribute no mass.
"""
struct Extensive <: AbstractMissingPolicy end
