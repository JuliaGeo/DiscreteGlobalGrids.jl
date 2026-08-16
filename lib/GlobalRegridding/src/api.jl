# The user-facing verbs. Bare `regrid` builds an ephemeral plan and drops it;
# `plan_regrid` hands one back. Applying a plan takes no keyword arguments
# because the plan already answers every question they would.

# Whether a source is chunked storage rather than a plain in-memory array.
# `DiskArrays.haschunks` is total, so this is total.
isdiskbacked(data) = DiskArrays.haschunks(data) isa DiskArrays.Chunked

"""
    regrid(data; to, method = Conservative(), missingpolicy = Weighted(0.5),
           lazy = isdiskbacked(data), chunks = nothing, budget = 2^30)
    regrid(data, plan::AbstractRegriddingPlan)

Regrid `data` onto the cells named by `to`.

Spatial dimensions are replaced by the destination's; every other dimension
passes through unchanged, and the output dimensions are
`(destination spatial dim, other dims...)`. One spatial plan is built and
reused across the non-spatial slices, so an N-D call costs one weight
construction, not one per field.

# Keyword arguments

  - `to`: the destination. The generic core takes a [`RegridSpace`](@ref);
    packages that supply spaces may accept their own target spellings and
    resolve them here.
  - `method`: an [`AbstractRegriddingMethod`](@ref). [`Conservative`](@ref) by
    default, the only exactly conservative choice.
  - `missingpolicy`: [`Weighted`](@ref) — coverage-normalized mean, blanking
    destinations below the threshold — or [`Extensive`](@ref) — raw sums.
  - `lazy`: return a [`LazyRegridArray`](@ref) instead of a materialized array.
    Defaults to whether `data` is chunked storage. Constructing one reads no
    source data at all.
  - `chunks`: the destination's chunking, as a `DiskArrays.GridChunks` or a
    tuple of chunk sizes. `nothing` derives it from the destination space.
  - `budget`: approximate bytes of source data that may be resident while one
    destination chunk is produced. Purely a performance knob — connected source
    chunks are held together when they fit and streamed one at a time when they
    do not, to the same answer either way.

The one-argument-plus-plan form takes no keyword arguments: a plan already
carries the method, both spaces, the missing policy, and the storage and budget
policy. Bare `regrid` is exactly `plan_regrid`, apply, drop — a one-shot call
retains nothing.
"""
function regrid end

regrid(data; to, method::AbstractRegriddingMethod = Conservative(),
    missingpolicy::AbstractMissingPolicy = Weighted(0.5),
    lazy::Bool = isdiskbacked(data), chunks = nothing, budget::Integer = 2^30) =
    error("implemented in T5")

regrid(data, plan::AbstractRegriddingPlan) = error("implemented in T5")

"""
    regrid!(dest, data; to, method = Conservative(), missingpolicy = Weighted(0.5),
            lazy = isdiskbacked(data), chunks = nothing, budget = 2^30)
    regrid!(dest, data, plan::AbstractRegriddingPlan)

Regrid `data` into the preallocated `dest` and return `dest`.

Destination chunks are filled one at a time over `eachchunk(dest)`, each with
its own scratch, so `dest` need not fit alongside the source. `dest`'s
non-spatial dimensions must match `data`'s and its spatial dimension must match
the destination space.

Keyword arguments are [`regrid`](@ref)'s. The plan form takes none.
"""
function regrid! end

regrid!(dest, data; to, method::AbstractRegriddingMethod = Conservative(),
    missingpolicy::AbstractMissingPolicy = Weighted(0.5),
    lazy::Bool = isdiskbacked(data), chunks = nothing, budget::Integer = 2^30) =
    error("implemented in T5")

regrid!(dest, data, plan::AbstractRegriddingPlan) = error("implemented in T5")

"""
    plan_regrid(data; to, method = Conservative(), missingpolicy = Weighted(0.5),
                lazy = isdiskbacked(data), chunks = nothing, budget = 2^30)
        -> AbstractRegriddingPlan

Build the regridding operator from `data`'s spaces onto `to` without applying
it, for reuse across fields, slices, and repeated reads.

Keyword arguments are [`regrid`](@ref)'s, and the plan retains all of them:
`regrid(data, plan)` and `regrid!(dest, data, plan)` accept none.

Building a plan reads **no source data**. Weights are geometry-only, and a
chunked plan's blocks are built on first touch rather than up front, so
planning a regrid of a terabyte source is as cheap as planning a megabyte one.
Reusing the plan is what makes the second and later applications free of weight
construction — the reason to reach for `plan_regrid` over bare
[`regrid`](@ref).
"""
plan_regrid(data; to, method::AbstractRegriddingMethod = Conservative(),
    missingpolicy::AbstractMissingPolicy = Weighted(0.5),
    lazy::Bool = isdiskbacked(data), chunks = nothing, budget::Integer = 2^30) =
    error("implemented in T5")
