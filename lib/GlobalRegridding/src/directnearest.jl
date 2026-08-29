"""
    DirectNearest()

Copy each destination value from the source cell containing its sample site.
The apply loop calls [`cellat`](@ref) directly and materializes no weights.

`DirectNearest` produces the same one-cell, unit-weight stencil as
[`NearestCell`](@ref), including missing-data results. Choose between them by
plan use:

# Selection guide

  - Prefer `DirectNearest` for one-shot and lazy reads. Its eager plan is `O(1)`,
    and each lazy tile stores one source index per destination cell.
  - Prefer `NearestCell` when one plan serves many sources or when callers need
    to inspect, store, or apply the weight operator. Its matrix locates every
    destination once and reuses those indices.

The fallback [`buildweights!`](@ref) delegates to `NearestCell`, preserving the
same result on generic assembly routes.
"""
struct DirectNearest <: AbstractRegriddingMethod end

outputsampling(::DirectNearest) = DD.Lookups.Points()

sourcesampling(::DirectNearest) = DD.Lookups.Points()

supportradius(::DirectNearest, ::RegridSpace) = 0.0

refinementinvariant(::DirectNearest) = true

buildweights!(coo::WeightCOO, ::DirectNearest, dst_space::RegridSpace, dst_inds,
    src_space::RegridSpace, src_inds) =
    buildweights!(coo, NearestCell(), dst_space, dst_inds, src_space, src_inds)

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
_directblank(::Weighted, ::Type{T}) where {T} = _maskedvalue(T)
_directblank(::Extensive, ::Type{T}) where {T} = convert(T, 0.0)

function regrid(data, plan::NearestDirectPlan)
    sd, othersizes, src = _flatten(data, plan)
    ndst = Int(ncells(plan.dst_space))
    out = Array{outputeltype(eltype(data))}(undef, ndst, othersizes...)
    applyplan!(reshape(out, ndst, prod(othersizes)), plan, src)
    return wrapoutput(out, data, sd, destinationdims(plan))
end

function regrid!(dest, data, plan::NearestDirectPlan)
    _, othersizes, src = _flatten(data, plan)
    ndst = Int(ncells(plan.dst_space))
    dstdims = destinationdims(plan)
    shaped = dstdims === nothing ? (ndst, othersizes...) :
             (map(length, dstdims)..., othersizes...)
    size(dest) == shaped || size(dest) == (ndst, othersizes...) ||
        throw(DimensionMismatch(
            "destination of size $(size(dest)) cannot hold a regrid of size $shaped"))
    raw = dest isa DD.AbstractDimArray ? parent(dest) : dest
    applyplan!(reshape(raw, ndst, prod(othersizes)), plan, src)
    return dest
end

"""
    applyplan!(out, plan::NearestDirectPlan, src) -> out

Sample the source directly into `out`, one destination cell at a time. Each
iteration locates the containing source cell, copies valid values, and writes
the policy's blank value for unmapped or invalid inputs.

Independent destination rows permit threading when multiple threads are
available.
"""
function applyplan!(out::AbstractMatrix, plan::NearestDirectPlan, src::AbstractMatrix)
    dst_space, src_space = plan.dst_space, plan.src_space
    ndst = Int(ncells(dst_space))
    size(src, 1) == Int(ncells(src_space)) || throw(DimensionMismatch(
        "plan expects $(ncells(src_space)) source cells, got $(size(src, 1))"))
    size(out, 1) == ndst || throw(DimensionMismatch(
        "plan produces $ndst destination cells, output holds $(size(out, 1))"))
    size(out, 2) == size(src, 2) || throw(DimensionMismatch(
        "$(size(src, 2)) source slices into $(size(out, 2)) output slices"))
    sites = samplesites(dst_space)
    blank = _directblank(plan.missingpolicy, eltype(out))
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

"""
    _readdestination!(out, A::LazyRegridArray{...,<:ChunkedPlan{DirectNearest}}, ...)

Read one [`DirectNearest`](@ref) destination tile through the existing executor
and source cache. Two passes replace weight construction and matrix application:

 1. Locate each destination's source cell and owning chunk, validate the unique
    chunk manifest against the dependency relation, then group destinations by
    chunk.
 2. Read each source chunk once through [`SourceHold`](@ref) and gather its
    destination values.

Per-tile state stores one source index per destination cell plus the output
tile.
"""
function _readdestination!(out::AbstractMatrix,
    A::LazyRegridArray{T,N,NS,NO,AA,P}, cellr::UnitRange{Int},
    others::NTuple{NO,UnitRange{Int}}, nslices::Int) where {T,N,NS,NO,AA,
    P<:ChunkedPlan{DirectNearest}}
    plan = A.plan
    src_space, dst_space = plan.src_space, plan.dst_space
    mv = plan.missingval
    blank = _directblank(plan.missingpolicy, T)
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
