# A weightless nearest-cell regrid: no COO, no sparse block, no per-cell weight.
#
# `NearestCell` answers "which source cell holds this destination's sample site"
# by *building* that answer into a sparse matrix of ones and then applying the
# matrix. Every entry is 1.0, every row has exactly one entry, and the apply is
# a gather dressed as a matvec. `DirectNearest` computes the same answer and
# skips both halves: preparation holds nothing but the two spaces, and the
# application asks `cellat` for the source cell and copies the value across.
#
# The result is defined to be identical to `NearestCell` under the same policy.
# Under `Weighted(t)` a mapped destination has cover 1 and reference 1, so
# `num/cover` is the source value itself and no threshold in `[0, 1]` can blank
# it; an unmapped destination, or one whose source value is invalid, has cover 0
# and is blanked. Under `Extensive` the same two cases sum to zero.

"""
    DirectNearest()

Take each destination's value from the source cell containing its sample site,
without building weights for it.

The stencil is [`NearestCell`](@ref)'s exactly — one source cell, weight one —
and so are the results, element for element, under either missing policy and
with any nodata sentinel. What differs is that no [`WeightCOO`](@ref),
[`WeightBlock`](@ref) or sparse matrix is ever assembled: a plan holds only the
two spaces, and the apply loop calls [`cellat`](@ref) per destination cell and
copies the value across.

# Which of the two to use

Prefer `DirectNearest` when the operator is applied about as many times as it is
built — a one-shot [`regrid`](@ref), or a chunked run whose plan serves one
destination column and is then dropped, which is every lazy read. It allocates
far less: the plan is `O(1)` rather than a row per destination cell, and a
chunked read carries one `Int` per cell of the tile instead of a sparse block
per source chunk. That is what it buys, and under many concurrent workers
sharing one heap it is most of what matters.

Prefer `NearestCell` when one plan is applied to *many different sources*. Its
matrix locates every destination once, for good, and each later application is a
gather over stored indices; `DirectNearest` re-locates every destination on
every application, and locating is the expensive half. Both methods locate once
for all the non-spatial slices of a single source, so a multi-slice source is
not the case this distinguishes. `NearestCell` is also the one to reach for when
the operator itself is wanted — to inspect, to store, or to apply outside this
package — since `DirectNearest` never materializes one.

[`buildweights!`](@ref) is still supplied, so any route that has not been
specialized for this method (a chunk pair built through [`pairblock`](@ref), a
weight file, a test) falls back to `NearestCell`'s own assembly and answers the
same.
"""
struct DirectNearest <: AbstractRegriddingMethod end

outputsampling(::DirectNearest) = DD.Lookups.Points()

supportradius(::DirectNearest, ::RegridSpace) = 0.0

# The fallback route, so nothing that has not been specialized breaks.
buildweights!(coo::WeightCOO, ::DirectNearest, dst_space::RegridSpace, dst_inds,
    src_space::RegridSpace, src_inds) =
    buildweights!(coo, NearestCell(), dst_space, dst_inds, src_space, src_inds)

# --------------------------------------------------------------------------
# The eager plan: preparation is the two spaces and nothing else
# --------------------------------------------------------------------------

"""
    NearestDirectPlan(method, missingpolicy, dst_space, src_space, missingval, sampling)

An eager [`DirectNearest`](@ref) plan. It holds no weights at all — the fields
are the two spaces, the policy, the nodata sentinel and the sampling override —
so building one costs nothing and allocates nothing per destination cell.
"""
struct NearestDirectPlan{M<:DirectNearest,P<:AbstractMissingPolicy,
                         D<:RegridSpace,S<:RegridSpace,V,
                         L<:Union{Nothing,DD.Lookups.Sampling}} <: AbstractRegriddingPlan
    method::M
    missingpolicy::P
    dst_space::D
    src_space::S
    missingval::V
    sampling::L
end

Base.show(io::IO, plan::NearestDirectPlan) =
    print(io, "NearestDirectPlan(", ncells(plan.dst_space), " ← ",
        ncells(plan.src_space), " cells, no weights)")

"""
    eagerplan(method, missingpolicy, dst_space, src_space, missingval, sampling)

The eager plan `method` wants. The default is one whole-domain
[`DirectPlan`](@ref) built through [`wholeblock`](@ref); a method whose apply
needs no weights returns a plan of its own here instead.
"""
eagerplan(method::AbstractRegriddingMethod, missingpolicy::AbstractMissingPolicy,
    dst_space::RegridSpace, src_space::RegridSpace, missingval, sampling) =
    DirectPlan(method, missingpolicy, dst_space, src_space,
        wholeblock(method, dst_space, src_space), missingval, sampling)

eagerplan(method::DirectNearest, missingpolicy::AbstractMissingPolicy,
    dst_space::RegridSpace, src_space::RegridSpace, missingval, sampling) =
    NearestDirectPlan(method, missingpolicy, dst_space, src_space, missingval, sampling)

destinationdims(plan::NearestDirectPlan) = destinationdims(plan.dst_space,
    something(plan.sampling, outputsampling(plan.method)))

# What an unmapped destination, or one whose source value is invalid, becomes.
_directblank(::Weighted, ::Type{T}, missingval) where {T} =
    _blankvalue(T, missingval)
_directblank(::Extensive, ::Type{T}, missingval) where {T} = convert(T, 0.0)

function regrid(data, plan::NearestDirectPlan; missingval = outputmissingval(data))
    sd, othersizes, src = _flatten(data, plan)
    ndst = Int(ncells(plan.dst_space))
    out = Array{outputeltype(eltype(data), missingval)}(undef, ndst, othersizes...)
    applyplan!(reshape(out, ndst, prod(othersizes)), plan, src, missingval)
    return wrapoutput(out, data, sd, destinationdims(plan), missingval)
end

function regrid!(dest, data, plan::NearestDirectPlan;
    missingval = destinationmissingval(dest))
    _, othersizes, src = _flatten(data, plan)
    ndst = Int(ncells(plan.dst_space))
    dstdims = destinationdims(plan)
    shaped = dstdims === nothing ? (ndst, othersizes...) :
             (map(length, dstdims)..., othersizes...)
    size(dest) == shaped || size(dest) == (ndst, othersizes...) ||
        throw(DimensionMismatch(
            "destination of size $(size(dest)) cannot hold a regrid of size $shaped"))
    raw = dest isa DD.AbstractDimArray ? parent(dest) : dest
    applyplan!(reshape(raw, ndst, prod(othersizes)), plan, src, missingval)
    return dest
end

"""
    applyplan!(out, plan::NearestDirectPlan, src, missingval = _maskedvalue(eltype(out))) -> out

Sample the source directly, one destination cell at a time. No weights are read
because none exist: the loop locates the source cell at the destination's sample
site and copies its value, blanking the destination with `missingval` where
there is no such cell or its value is invalid.

The loop is over destination cells and every iteration writes only its own row,
so it is threaded when more than one thread is available.
"""
function applyplan!(out::AbstractMatrix, plan::NearestDirectPlan, src::AbstractMatrix,
    missingval = _maskedvalue(eltype(out)))
    dst_space, src_space = plan.dst_space, plan.src_space
    ndst = Int(ncells(dst_space))
    size(src, 1) == Int(ncells(src_space)) || throw(DimensionMismatch(
        "plan expects $(ncells(src_space)) source cells, got $(size(src, 1))"))
    size(out, 1) == ndst || throw(DimensionMismatch(
        "plan produces $ndst destination cells, output holds $(size(out, 1))"))
    size(out, 2) == size(src, 2) || throw(DimensionMismatch(
        "$(size(src, 2)) source slices into $(size(out, 2)) output slices"))
    sites = samplesites(dst_space)
    blank = _directblank(plan.missingpolicy, eltype(out), missingval)
    mv = plan.missingval
    if Threads.nthreads() > 1 && !OUTER_PARALLEL[]
        Threads.@threads for j in 1:ndst
            _directcell!(out, src, src_space, sites, j, blank, mv)
        end
    else
        for j in 1:ndst
            _directcell!(out, src, src_space, sites, j, blank, mv)
        end
    end
    return out
end

@inline function _directcell!(out::AbstractMatrix, src::AbstractMatrix,
    src_space::RegridSpace, sites, j::Int, blank, mv)
    T = eltype(out)
    i = cellat(src_space, @inbounds sites[j])
    if i === nothing
        @inbounds for s in axes(out, 2)
            out[j, s] = blank
        end
        return out
    end
    k = Int(i)
    @inbounds for s in axes(out, 2)
        x = src[k, s]
        out[j, s] = _isvalid(x, mv) ? convert(T, x) : blank
    end
    return out
end

# --------------------------------------------------------------------------
# The lazy path: the executor's tiles and source reads, sampled directly
# --------------------------------------------------------------------------

"""
    _readdestination!(out, A::LazyRegridArray{...,<:ChunkedPlan{DirectNearest}}, ...)

The chunked read for [`DirectNearest`](@ref). It reuses the executor's tiling,
its dependency-relation reach check, its source hold and its chunk reads, and
replaces the weight build and the matvec with two passes over the tile:

 1. locate: one [`cellat`](@ref) per destination cell of the tile, giving the
    source cell it takes its value from and the source chunk that owns it. The
    ascending unique chunks are the tile's manifest — exactly what
    [`tileweights`](@ref) would have produced — and are checked against the
    plan's relation the same way. Destination cells are then ordered by chunk
    with a counting sort, so each chunk is read once and touched once.
 2. gather: for each source chunk in ascending order, read its buffer through
    the same [`SourceHold`](@ref) and copy each destination's value out of it.

Per-tile state is one `Int` per destination cell plus the output tile itself.
No `WeightCOO`, no `sparse`, no `WeightBlock`, no accumulator pair.
"""
function _readdestination!(out::AbstractMatrix,
    A::LazyRegridArray{T,N,NS,NO,AA,P}, cellr::UnitRange{Int},
    others::NTuple{NO,UnitRange{Int}}, nslices::Int) where {T,N,NS,NO,AA,
    P<:ChunkedPlan{DirectNearest}}
    plan = A.plan
    src_space, dst_space = plan.src_space, plan.dst_space
    mv = plan.missingval
    blank = _directblank(plan.missingpolicy, T, A.missingval)
    sites = samplesites(dst_space)
    groups = _slicegroups(A, others)
    strides = _slicestrides(others)
    hold = SourceHold(_emptybuffer(A), databudget(plan.budget), A.stats)
    srcchunks = Int[]
    srcranges = NTuple{NS,UnitRange{Int}}[]
    slots = zeros(Int, Int(nchunks(src_space)))
    srcidx = Int[]
    chunkof = Int[]
    perm = Int[]
    starts = Int[]
    cursor = Int[]
    vals = Matrix{T}(undef, 0, nslices)
    for t in _coveringtiles(A, cellr)
        dinds = _tileindices(A, t)
        nd = length(dinds)
        length(srcidx) == nd || (srcidx = Vector{Int}(undef, nd))
        length(chunkof) == nd || (chunkof = Vector{Int}(undef, nd))
        length(perm) == nd || (perm = Vector{Int}(undef, nd))
        size(vals) == (nd, nslices) || (vals = Matrix{T}(undef, nd, nslices))
        fill!(vals, blank)
        # --- pass 1: locate ------------------------------------------------
        empty!(srcchunks)
        @inbounds for j in 1:nd
            s = cellat(src_space, sites[Int(dinds[j])])
            if s === nothing
                srcidx[j] = 0
                chunkof[j] = 0
            else
                si = Int(s)
                srcidx[j] = si
                c = Int(chunkat(src_space, si))
                chunkof[j] = c
                if slots[c] == 0
                    slots[c] = 1
                    push!(srcchunks, c)
                end
            end
        end
        sort!(srcchunks)
        nc = length(srcchunks)
        @inbounds for k in 1:nc
            slots[srcchunks[k]] = k
        end
        # The manifest must sit inside the relation's rows for this tile, for
        # the same reason `_selectsource!` checks a `TileWeights` manifest.
        bad = _outsiderows(A, t, srcchunks)
        bad == 0 || throw(ArgumentError(_reachmessage(A, t, bad)))
        _sourceranges!(srcranges, A, srcchunks)
        # Counting sort of the tile's cells by source chunk.
        length(starts) == nc + 1 || (starts = Vector{Int}(undef, nc + 1))
        length(cursor) == nc + 1 || (cursor = Vector{Int}(undef, nc + 1))
        fill!(starts, 0)
        @inbounds for j in 1:nd
            c = chunkof[j]
            c == 0 && continue
            starts[slots[c]] += 1
        end
        acc = 1
        @inbounds for k in 1:nc
            n = starts[k]
            starts[k] = acc
            cursor[k] = acc
            acc += n
        end
        starts[nc+1] = acc
        @inbounds for j in 1:nd
            c = chunkof[j]
            c == 0 && continue
            k = slots[c]
            perm[cursor[k]] = j
            cursor[k] += 1
        end
        @inbounds for k in 1:nc
            slots[srcchunks[k]] = 0
        end
        # --- pass 2: gather ------------------------------------------------
        keep = _canhold(hold, srcranges, groups, sizeof(eltype(hold.scratch)))
        for k in 1:nc
            s = srcchunks[k]
            sr = srcranges[k]
            ncell = prod(map(length, sr))
            cmap = indexmap(ownedindices(src_space, s))
            rng = starts[k]:(starts[k+1]-1)
            for (gi, gpos) in enumerate(groups)
                gr = _grouprange(others, gpos)
                if knownempty(A.source, (sr..., gr...))
                    A.stats.skipped += 1
                    continue
                end
                buf = _sourcefor!(hold, A, (s, gi), sr, gr, gpos, keep)
                _scattergroup!(vals, buf, ncell, gpos, strides, srcidx, perm, rng,
                    cmap, mv, blank)
            end
        end
        _writedirect!(out, vals, dinds, cellr)
    end
    return out
end

# Copy one source chunk's values into the destination cells that named it.
function _scattergroup!(vals::Matrix{T}, buf::AbstractArray, ncell::Int,
    pos::NTuple{NO,UnitRange{Int}}, strides::NTuple{NO,Int}, srcidx::Vector{Int},
    perm::Vector{Int}, rng::UnitRange{Int}, cmap, mv, blank) where {T,NO}
    src = reshape(buf, ncell, :)
    kk = 0
    for I in CartesianIndices(pos)
        kk += 1
        col = 1
        @inbounds for d in 1:NO
            col += (I[d] - 1) * strides[d]
        end
        @inbounds for q in rng
            j = perm[q]
            c = localindex(cmap, srcidx[j])
            x = src[c, kk]
            vals[j, col] = _isvalid(x, mv) ? convert(T, x) : blank
        end
    end
    return vals
end

# Scatter the tile's finished values into the requested cell window.
function _writedirect!(out::AbstractMatrix, vals::Matrix, dinds,
    cellr::AbstractUnitRange)
    lo, hi = first(cellr), last(cellr)
    off = lo - 1
    for col in axes(vals, 2)
        @inbounds for (j, p) in enumerate(dinds)
            lo <= p <= hi || continue
            out[p-off, col] = vals[j, col]
        end
    end
    return out
end
