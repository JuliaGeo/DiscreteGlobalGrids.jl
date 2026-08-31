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
    BarycentricPoint(; poles = NearestCell())

Interpolate between source sample sites at each destination sample site.

  - The stencil is the dual cell of source sample sites containing the point,
    weighted by the coordinates that cell's basis names: tensor Q1 on a
    quadrilateral, or mean-value coordinates on a convex polygon, which on a
    triangle are that triangle's barycentric coordinates.
  - Weights are nonnegative and sum to one, so the result lies between the
    source values it came from. Integrals are not preserved.
  - Requires [`cellcentroid`](@ref) of the destination space and a source space
    that answers point queries; a destination outside the source's dual complex
    emits no entry at all and the missing policy decides what it becomes.
  - `poles` is the policy where a source's own sample sites stop short of a
    pole: [`NearestCell`](@ref) takes the nearest site of the polemost row with
    weight one, `nothing` leaves those points unmapped. A source whose sites
    reach the poles — every conforming grid — has no such region and ignores it.
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

Return the sampling a method gives the destination it writes. Area-based methods
report `Intervals(Center())`, the default; point samples report `Points()`.

This trait also selects the build route. A chunked plan selects the whole-tile
route with it, through [`tilesampler`](@ref): a method reporting `Points()` and
supplying a [`sampler`](@ref) builds one destination tile at a time.
[`weightblock`](@ref) dispatches on it as well, so a sampling may specialise the
per-pair assembly; today every sampling assembles a pair through
[`pairblock`](@ref).
"""
outputsampling(::AbstractRegriddingMethod) = DD.Lookups.Intervals(DD.Lookups.Center())
outputsampling(::NearestCell) = DD.Lookups.Points()
outputsampling(::BarycentricPoint) = DD.Lookups.Points()

# Weight construction

"""
    CoverageCOO()

A chunk-local coordinate list of coverage weights: the non-negative share of a
destination one source cell covers, in the weight list's own convention.
[`WeightCOO`](@ref) holds one only when its value weights cannot serve as
coverage themselves.
"""
struct CoverageCOO
    rows::Vector{Int}
    cols::Vector{Int}
    vals::Vector{Float64}
end

CoverageCOO() = CoverageCOO(Int[], Int[], Float64[])

Base.length(cov::CoverageCOO) = length(cov.vals)

"""
    WeightCOO(ndst::Int)

A chunk-local coordinate-list accumulator over `ndst` destination cells. `rows`
and `cols` are chunk-local indices within the builder's `dst_inds` and
`src_inds`. Duplicate entries are summed when the block is assembled.

`denom` holds optional per-destination denominators and is `nothing` until a
builder declares them, through [`markdenominated!`](@ref) or the first
[`adddenom!`](@ref); a method that reports none — every point sample — leaves it
`nothing` and allocates no denominator vector.

`coverage` holds an optional [`CoverageCOO`](@ref), `nothing` until
[`markcovered!`](@ref) or the first [`addcoverage!`](@ref) declares it. Weights
that are all non-negative double as coverage and need no list. Signed weights do
— summed over the valid sources they measure nothing, being off by the
correction and possibly negative.
"""
mutable struct WeightCOO
    const ndst::Int
    const rows::Vector{Int}
    const cols::Vector{Int}
    const vals::Vector{Float64}
    denom::Union{Nothing,Vector{Float64}}
    coverage::Union{Nothing,CoverageCOO}
end

WeightCOO(ndst::Integer) =
    WeightCOO(Int(ndst), Int[], Int[], Float64[], nothing, nothing)

Base.length(coo::WeightCOO) = length(coo.vals)

Base.show(io::IO, coo::WeightCOO) =
    print(io, "WeightCOO(ndst=", coo.ndst, ", entries=", length(coo),
        coo.denom === nothing ? "" : ", denom",
        coo.coverage === nothing ? "" : ", coverage", ")")

"""
    addweight!(coo::WeightCOO, dst_local::Int, src_local::Int, w::Real)

Add `w` to the value weight of source `src_local` in destination `dst_local`,
both chunk-local indices within the builder's `dst_inds` and `src_inds`.

`w` may be negative, in which case the method must report coverage separately
through [`addcoverage!`](@ref): such weights no longer measure how much of a
destination a source covers.
"""
function addweight!(coo::WeightCOO, dst_local::Int, src_local::Int, w::Real)
    push!(coo.rows, dst_local)
    push!(coo.cols, src_local)
    push!(coo.vals, Float64(w))
    return coo
end

"""
    addcoverage!(coo::WeightCOO, dst_local::Int, src_local::Int, a::Real)

Add `a` to the coverage of destination `dst_local` by source `src_local`, in
[`addweight!`](@ref)'s chunk-local indices, and declare that `coo` carries
coverage. `a` is the non-negative weight — an overlap area, conservatively —
`src_local` contributes when valid.

Only signed methods need this; non-negative weights serve as their own coverage.
"""
function addcoverage!(coo::WeightCOO, dst_local::Int, src_local::Int, a::Real)
    cov = coo.coverage
    if cov === nothing
        cov = CoverageCOO()
        coo.coverage = cov
    end
    push!(cov.rows, dst_local)
    push!(cov.cols, src_local)
    push!(cov.vals, Float64(a))
    return coo
end

"""
    signedweights(method::AbstractRegriddingMethod) -> Bool

Return whether `method` emits negative value weights, and so reports coverage
separately ([`addcoverage!`](@ref)). Defaults to `false`.

The lazy path reads this to size accumulators before it has a block: signed
weights let a hole reach a destination without reaching its coverage, so the
executor carries the two extra sums [`degradetainted!`](@ref) needs. A method
forwarding [`buildweights!`](@ref) to another should forward this too.
"""
signedweights(::AbstractRegriddingMethod) = false

"""
    markcovered!(coo::WeightCOO)

Declare that `coo` carries a coverage list, adding nothing to it, and return
`coo`. Allocating the empty list here gives a builder that covers no destination
a block whose coverage is a zero operator, not its signed value weights.
"""
function markcovered!(coo::WeightCOO)
    coo.coverage === nothing && (coo.coverage = CoverageCOO())
    return coo
end

"""
    hascoverage(coo::WeightCOO) -> Bool

Return whether `coo` carries a coverage list of its own. When it does not, its
value weights are the coverage.
"""
hascoverage(coo::WeightCOO) = coo.coverage !== nothing

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

Declare that `coo` carries denominators, without adding to any of them, and
return `coo`. This is where the zero-filled denominator vector is allocated, so
a builder that reports coverage for no destination still produces a denominated
block of zeros.
"""
function markdenominated!(coo::WeightCOO)
    coo.denom === nothing && (coo.denom = zeros(Float64, coo.ndst))
    return coo
end

# Per-point weights

"""
    WeightRow()

The source cells one destination point takes its value from, and their weights.

`indices` and `weights` are parallel: entry `k` gives weight `weights[k]` to the
source cell the source space calls `indices[k]`. [`weightsat!`](@ref) clears the
row on entry, so one row serves a whole sweep and grows only to the widest
stencil it has met. A row belongs to one task.
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

# A zero weight is dropped rather than stored, so a point on an edge or at a
# node emits only the nodes that carry it.
@inline function _addentry!(row::WeightRow, i::Integer, w::Float64)
    w > 0 || return row
    push!(row.indices, Int(i))
    push!(row.weights, w)
    return row
end

"""
    WeightStatus

What [`weightsat!`](@ref) made of one destination point.

`WeightsMapped` says the row holds a complete stencil. Every other value says
the row is empty and the destination takes no weights at all; they differ only
in the reason, which is diagnostic. Execution reads
[`ismapped`](@ref) and nothing else.

  - `WeightsOutside` — the source covers no cell the point belongs to: it lies
    outside the dual complex, outside the cell that was offered for it, or
    outside the source's coverage altogether;
  - `WeightsRim` — a required sample site is not in the source collection, so no
    dual cell exists there. A space's own `weightsat!` answers this;
  - `WeightsDegenerate` — the cell is unusable: repeated nodes, no area, a fold,
    a reflex corner, or an inverse map that did not converge.
"""
@enum WeightStatus::UInt8 begin
    WeightsMapped
    WeightsOutside
    WeightsRim
    WeightsDegenerate
end

"""
    ismapped(status::WeightStatus) -> Bool

Whether `status` says the row carries a stencil.
"""
@inline ismapped(status::WeightStatus) = status === WeightsMapped

"""
    buildweights!(coo, method, dst_space, dst_inds, src_space, src_inds)

Append chunk-local weights for `dst_inds` and `src_inds`, then return `coo`.
Builders may inspect geometry outside `src_inds`, but must emit weights only for
sources inside it. Otherwise weights are duplicated across chunk blocks.

Weight construction must not depend on field data or execution order.

This is the one hook a method must supply, and the generic assembly every
[`pairblock`](@ref) falls back to. A method that wraps another and forwards
`buildweights!` therefore builds through this route even where the inner method
assembles a block of its own; forward [`pairblock`](@ref) as well — and
[`sampler`](@ref) too, for a point method — to take the inner method's own
build.
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

This is a *bound* on the stencil, not an estimate of it. Chunk discovery is cap
overlap plus this radius, and cap overlap alone is not a superset of a point
stencil's reach: a destination cell finer than a source cell can be bracketed by
sample sites whose own cells it never touches. A method that supplies a
[`sampler`](@ref) therefore declares here how far its stencils reach, so that
the dependency relation stays a superset of the chunks its tiles read; a chunked
read refuses a tile whose weights name a chunk the relation does not.
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

The valid source weight is accumulated from the block's coverage weights, or
from its value weights when it reports none ([`WeightBlock`](@ref)); a method
with signed value weights is therefore normalized by the coverage it declared,
not by its own weights.
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
