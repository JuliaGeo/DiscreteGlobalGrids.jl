# Regridding methods build geometry-only weights for a pair of chunks.

"""
    AbstractRegriddingMethod

How source cell values are combined into destination cell values.

Methods implement [`buildweights!`](@ref), plus [`supportradius`](@ref) when
their stencil extends beyond overlapping cells. Weights must be linear and
independent of field data.
"""
abstract type AbstractRegriddingMethod end

"""
    Conservative()

Weight source cells by their spherical intersection area with each destination.
With [`Extensive`](@ref), this preserves the covered integral; with
[`Weighted`](@ref), it returns coverage-normalized means. Requires cell polygons.
"""
struct Conservative <: AbstractRegriddingMethod end

"""
    NearestCell()

Give weight 1 to the source cell containing each destination centroid.

Requires [`cellcentroid`](@ref) of the destination space and [`cellat`](@ref)
of the source space. A destination centroid outside the source's coverage emits
no entry at all; the missing policy decides what that destination cell becomes.

This method does not preserve integrals.
"""
struct NearestCell <: AbstractRegriddingMethod end

"""
    BilinearPoint()

Bilinearly interpolate the source chart at each destination centroid.

The stencil is written on the source space's chart, so the source must answer
`true` to [`hascellchart`](@ref); the destination must provide
[`cellcentroid`](@ref).

This point sample does not preserve integrals.
"""
struct BilinearPoint <: AbstractRegriddingMethod end

"""
    outputsampling(method::AbstractRegriddingMethod) -> DimensionalData.Lookups.Sampling

Return the sampling a method gives the destination it writes. Area-based methods
report `Intervals(Center())`, the default; point samples report `Points()`.

This trait also selects the method's weight build path in [`weightblock`](@ref):
`Points()` takes the point path, every other sampling the area path.
"""
outputsampling(::AbstractRegriddingMethod) = DD.Lookups.Intervals(DD.Lookups.Center())
outputsampling(::NearestCell) = DD.Lookups.Points()
outputsampling(::BilinearPoint) = DD.Lookups.Points()

# Weight construction

"""
    WeightCOO(ndst::Int)

A chunk-local coordinate-list accumulator. `rows` and `cols` are chunk-local
indices within the builder's `dst_inds` and `src_inds`. `denom` stores optional
per-destination denominators. Duplicate entries are summed when the block is
assembled.
"""
mutable struct WeightCOO
    const rows::Vector{Int}
    const cols::Vector{Int}
    const vals::Vector{Float64}
    const denom::Vector{Float64}
    hasdenom::Bool
end

WeightCOO(ndst::Integer) =
    WeightCOO(Int[], Int[], Float64[], zeros(Float64, ndst), false)

Base.length(coo::WeightCOO) = length(coo.vals)

Base.show(io::IO, coo::WeightCOO) =
    print(io, "WeightCOO(ndst=", length(coo.denom), ", entries=", length(coo),
        coo.hasdenom ? ", denom" : "", ")")

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
    coo.denom[dst_local] += Float64(d)
    coo.hasdenom = true
    return coo
end

"""
    markdenominated!(coo::WeightCOO)

Declare that `coo` carries denominators, without adding to any of them.
"""
function markdenominated!(coo::WeightCOO)
    coo.hasdenom = true
    return coo
end

"""
    buildweights!(coo, method, dst_space, dst_inds, src_space, src_inds)

Append chunk-local weights for `dst_inds` and `src_inds`, then return `coo`.
Builders may inspect geometry outside `src_inds`, but must emit weights only for
sources inside it. Otherwise weights are duplicated across chunk blocks.

Weight construction must not depend on field data or execution order.
"""
function buildweights!(coo::WeightCOO, method::AbstractRegriddingMethod,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    throw(ArgumentError(
        "buildweights! is not implemented for $(typeof(method)) from " *
        "$(typeof(src_space)) to $(typeof(dst_space))"))
end

"""
    supportradius(method, src_space::RegridSpace) -> Float64

Return the maximum angular distance, in radians, that the method's stencil
extends beyond a source chunk. The default is `0.0`. Overestimates add discovery
work; underestimates can omit required weights.
"""
supportradius(::AbstractRegriddingMethod, ::RegridSpace) = 0.0

# `build_weights!` and `support_radius` are the old names of `buildweights!`
# and `supportradius` and forward to them, so a call of an old name answers the
# same with a deprecation warning. Only callers are carried: a method that
# defines an old name supplies neither hook, and the generics dispatch on the
# new ones.

"""
    build_weights!(coo, method, dst_space, dst_inds, src_space, src_inds)

Deprecated. Use [`buildweights!`](@ref), which this forwards to, so existing
calls keep their old behaviour exactly.
"""
function build_weights! end

@deprecate build_weights!(coo::WeightCOO, method::AbstractRegriddingMethod,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds) buildweights!(
    coo, method, dst_space, dst_inds, src_space, src_inds) false

"""
    support_radius(method, src_space::RegridSpace) -> Float64

Deprecated. Use [`supportradius`](@ref), which this forwards to, so existing
calls keep their old behaviour exactly.
"""
function support_radius end

@deprecate support_radius(method::AbstractRegriddingMethod, src_space::RegridSpace) supportradius(
    method, src_space) false

# Missing-data policies

"""
    AbstractMissingPolicy

How accumulated destinations handle missing source values.
"""
abstract type AbstractMissingPolicy end

"""
    Weighted(threshold::Real = 0.5)

Return each weighted sum divided by its valid source weight. Destinations below
the valid-coverage `threshold` become missing. `threshold` must be in `[0, 1]`.
Coverage is measured against the source-covered portion, not the full
destination area.
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
