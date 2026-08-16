# Methods and the one hook they implement. A method is asked for weights over a
# pair of chunks and given nothing else — no data, no IO, no array — so the same
# weights serve eager, chunked, and streaming execution unchanged.

"""
    AbstractRegriddingMethod

How source cell values are combined into destination cell values.

Every method here is **linear in the source data**, so every method compiles to
the same object: a sparse weight matrix with an optional per-destination
denominator. A method therefore implements exactly one thing,
[`build_weights!`](@ref), plus [`support_radius`](@ref) if its stencil reaches
beyond the cells it overlaps.

Data-dependent rules — nearest *valid* neighbour, gap-filling — are not methods
in this sense and are out of scope: weights are geometry-only, and geometry
does not know where the missing values are.
"""
abstract type AbstractRegriddingMethod end

"""
    Conservative()

Area-weighted regridding: each weight is the spherical area of the intersection
of a source cell with a destination cell.

The default, and the only exactly conservative method here — with
[`Extensive`](@ref) it preserves the global integral, and with
[`Weighted`](@ref) it returns the coverage-normalized cell mean.

Requires cell polygons of both spaces and nothing else, so it works between any
two spaces.
"""
struct Conservative <: AbstractRegriddingMethod end

"""
    NearestCell()

One weight-1 entry per destination cell, at the source cell containing that
destination cell's centroid.

Requires [`cellcentroid`](@ref) of the destination space and [`cellat`](@ref)
of the source space. A destination centroid outside the source's coverage emits
no entry at all; the missing policy decides what that destination cell becomes.

Not conservative: it neither preserves integrals nor respects cell areas.
"""
struct NearestCell <: AbstractRegriddingMethod end

"""
    BilinearPoint()

Bilinear interpolation of the source field, evaluated at the destination cell's
centroid.

The stencil is written on the source space's chart, so the source must answer
`true` to [`hascellchart`](@ref); the destination must provide
[`cellcentroid`](@ref).

Not conservative, and not integral-consistent either — it samples a point
rather than averaging over the destination cell.
"""
struct BilinearPoint <: AbstractRegriddingMethod end

# ===========================================================================
# The weight-building hook
# ===========================================================================

"""
    WeightCOO(ndst::Int)

The coordinate-list accumulator a method appends weights to.

`rows` and `cols` are **chunk-local**: `rows[k]` indexes into the `dst_inds`
vector the builder was handed, `cols[k]` into `src_inds`. Cell positions never
enter — that is what lets one block be built once and reused wherever those two
chunks meet.

`denom[j]` accumulates the per-destination denominator for local destination
`j`, which for [`Conservative`](@ref) is the destination cell's covered area.
It is carried only if a builder ever calls [`adddenom!`](@ref); a block built
without one finalizes as a raw sum.

Duplicate `(row, col)` entries are summed when the block is assembled, so a
builder may emit a destination's shares one at a time.
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

Add `w` to the weight of local source `src_local` in local destination
`dst_local`. Indices are positions within the builder's `dst_inds` and
`src_inds`, not cell positions.
"""
function addweight!(coo::WeightCOO, dst_local::Int, src_local::Int, w::Real)
    push!(coo.rows, dst_local)
    push!(coo.cols, src_local)
    push!(coo.vals, Float64(w))
    return coo
end

"""
    adddenom!(coo::WeightCOO, dst_local::Int, d::Real)

Add `d` to local destination `dst_local`'s denominator, and mark the block as
carrying one. Blocks over the same destination chunk sum their denominators, so
a builder reports only the share its own source chunk contributes.
"""
function adddenom!(coo::WeightCOO, dst_local::Int, d::Real)
    coo.denom[dst_local] += Float64(d)
    coo.hasdenom = true
    return coo
end

"""
    build_weights!(coo, method, dst_space, dst_inds, src_space, src_inds)

Append `method`'s weights for the cell pairs between `dst_inds` of `dst_space`
and `src_inds` of `src_space` to `coo`, and return `coo`.

`dst_inds` and `src_inds` are the cell positions of one chunk on each side, as
[`cellindices`](@ref) returns them. Entries are addressed by position within
those vectors; see [`WeightCOO`](@ref).

The **partition invariant**: a builder may read geometry anywhere — a stencil
straddling a source chunk boundary must look at cell centres outside
`src_inds` — but may emit entries only for source cells that are *in*
`src_inds`. Every entry then belongs to exactly one block, and the executor's
accumulation over source chunks sums a straddling stencil's shares back to the
whole. A builder that emits an out-of-chunk entry silently double-counts.

Nothing about data, chunk order, or execution reaches here. The same call must
produce the same weights whether the caller is regridding one field or a
thousand, eagerly or one chunk at a time.

Implementations that need to locate points or read a chart require
[`cellat`](@ref), [`cellcentroid`](@ref), or [`hascellchart`](@ref) of the
relevant space and should say so when they error.
"""
function build_weights!(coo::WeightCOO, method::AbstractRegriddingMethod,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    throw(ArgumentError(
        "build_weights! is not implemented for $(typeof(method)) from " *
        "$(typeof(src_space)) to $(typeof(dst_space))"))
end

"""
    support_radius(method, src_space::RegridSpace) -> Float64

The angular radius, in radians, by which a source chunk's extent must be
dilated before it is tested for connection to a destination chunk. Defaults to
`0.0`.

Chunk-connectivity discovery compares chunk extents. That is exact only for a
method whose weights vanish outside the geometric overlap — `0.0` for
[`Conservative`](@ref) and [`NearestCell`](@ref). A method whose stencil reaches
past the cells it overlaps must report that reach here, or a destination chunk
lying wholly inside one source chunk's footprint will never be paired with the
neighbouring chunk its stencil needs, and the missing policy will quietly
renormalize a truncated stencil instead of reporting anything.

Report an upper bound: half a source cell for a 4-point stencil, the patch
radius for a least-squares patch. Over-reporting costs discovery work; under-
reporting is wrong.
"""
support_radius(::AbstractRegriddingMethod, ::RegridSpace) = 0.0

# ===========================================================================
# Missing-data policies
# ===========================================================================

"""
    AbstractMissingPolicy

What the executor does with a destination cell once its weighted source
contributions have been accumulated — the single place missing data is spelled.

[`Weighted`](@ref) or [`Extensive`](@ref).
"""
abstract type AbstractMissingPolicy end

"""
    Weighted(threshold::Real = 0.5)

Coverage-normalized mean: each destination value is its accumulated weighted
sum divided by the accumulated weight of the source cells that were **not**
missing, and destinations whose valid coverage falls below `threshold` — as a
fraction of the block's denominator — are themselves missing.

The default policy, and the one that makes a regrid invariant to how the source
is padded: adding missing cells around a source field changes neither the
divisor nor the result.

`threshold` is a fraction in `[0, 1]`. `0` keeps every destination that saw any
valid source at all; `1` keeps only fully covered destinations.
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

Raw sums: each destination value is its accumulated weighted sum, undivided.

With [`Conservative`](@ref) weights this is the exactly conservative answer —
the global integral of the destination field equals that of the source over the
covered region — and it is the policy a conservation test asserts against.
Missing source cells contribute nothing, which is a real loss of mass rather
than something to normalize away.
"""
struct Extensive <: AbstractMissingPolicy end
