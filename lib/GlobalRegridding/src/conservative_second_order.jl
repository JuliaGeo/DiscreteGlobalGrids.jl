# Second-order conservative weights: overlap areas plus a first-moment gradient
# correction, `W = A + M G`, with the coverage `A` kept apart from the signed `W`.

"""
    ConservativeSecondOrder()

Weight source cells by their spherical intersection area with each destination,
and correct each overlap by the source cell's recovered gradient, so that a
field varying linearly across a source cell is carried into the destinations
that split it rather than flattened to the cell's mean.

  - Each source cell `i` is reconstructed as `f̄ᵢ + gᵢ ⋅ (u − ūᵢ)` in the tangent
    plane at its own mean position, where `gᵢ` is the unweighted least-squares
    fit of the neighbouring cells' means ([`cellneighbors`](@ref),
    [`gradientstencil!`](@ref)). The correction has zero mean over the cell, so
    the integral of every source cell the destination covers entirely is
    preserved exactly whatever the gradient; that, not positivity, is what
    makes the method conservative. A cell the destination covers only in part
    — at the rim of a partial destination — hands over the reconstruction's
    integral over the covered part, which differs from the flat share
    [`Conservative`](@ref) hands over by the gradient's moment there.
  - A destination's value weights are `A + M G`: the overlap areas
    ([`Conservative`](@ref) exactly, bit for bit) plus each overlap's first
    moment about the source mean folded through the gradient stencil. Weights
    are therefore signed and can reach one cell past the overlapping sources;
    a destination can overshoot the source values it draws on.
  - The error in a smooth field is second order in the source cell width,
    against first order for [`Conservative`](@ref): each cell's reconstruction
    is exact for a field linear in that cell's tangent coordinates, and the
    chart's curvature is the order below. A cell with fewer than two
    non-collinear neighbours keeps a zero gradient, which is first order there.
  - Coverage is the non-negative overlap area, reported apart from the signed
    weights, so [`Weighted`](@ref) normalizes and thresholds by area covered
    and [`Extensive`](@ref) preserves the covered integral. Missing source
    values inside a gradient stencil are read as zero by the fixed operator,
    which biases that cell's gradient; where a source carries missing values,
    fill it first or use [`Conservative`](@ref).
  - Requires cell polygons on both sides, and on the source
    [`cellneighbors`](@ref) and [`celldiameter`](@ref); both have geometric
    fallbacks. Source and destination manifolds must match.
"""
struct ConservativeSecondOrder <: AbstractRegriddingMethod end

# A block's weights for the sources of a chunk read the overlaps of the cells
# around them, one ring out, so discovery must pair a destination with every
# source chunk within one cell of it.
supportradius(::ConservativeSecondOrder, src_space::RegridSpace) =
    celldiameter(src_space)

preparesdestination(::ConservativeSecondOrder, dst_space::RegridSpace) =
    expensivecellgeometry(dst_space)

# `BlockAreaOperator` keeps an overlap the way the pair operator does.
@inline _covers(pm::PolygonMoments) = pm.area > 0

# Gradient stencils

"""
    GradientStencils

The reconstruction of every cell in a block's extended source set — a chunk's
cells and their one ring — ready to fold overlaps through, in one flat layout:
cell `p` has tangent frame `(e1[p], e2[p])` and mean position `mean[p]`, and
its value-weight terms are `cols[k] => coeffs[k]` for `k in ptr[p]:ptr[p+1]-1`,
such that the cell's gradient is `Σ coeffs[k] f[cols[k]]` — the self term
already merged, so a constant field has zero gradient exactly. `cols` are
chunk-local source columns only; a neighbour outside the chunk belongs to
another block's stencil, which folds its own share of this cell. A cell with
no terms has no gradient: no area, or too few neighbours to fix one.

Flat rather than one vector per cell because a whole-space build sweeps
millions of cells, and the sweep is otherwise its allocations.
"""
struct GradientStencils
    e1::Vector{NTuple{3,Float64}}
    e2::Vector{NTuple{3,Float64}}
    mean::Vector{NTuple{3,Float64}}
    ptr::Vector{Int}
    cols::Vector{Int}
    coeffs::Vector{NTuple{2,Float64}}
end

@inline _meanposition(pm::PolygonMoments) =
    (pm.moment[1] / pm.area, pm.moment[2] / pm.area, pm.moment[3] / pm.area)

const _NOWHERE = (0.0, 0.0, 0.0)

# Mean positions of the extended set, by position in it, and of the ring beyond
# it — a neighbour of a halo cell — by index, looked up as they are met. A
# whole space has no ring beyond it and never touches the table.
struct _MeanPositions{M}
    means::Vector{NTuple{3,Float64}}
    extmap::M
    beyond::Dict{Int,NTuple{3,Float64}}
    space::RegridSpace
end

function _MeanPositions(space::RegridSpace, ext_inds, extmap)
    means = Vector{NTuple{3,Float64}}(undef, length(ext_inds))
    for (p, i) in enumerate(ext_inds)
        pm = cellmoments(space, i)
        means[p] = pm.area > 0 ? _meanposition(pm) : _NOWHERE
    end
    return _MeanPositions(means, extmap, Dict{Int,NTuple{3,Float64}}(), space)
end

@inline function (mp::_MeanPositions)(i::Int)
    p = localindex(mp.extmap, i)
    p == 0 || return @inbounds mp.means[p]
    return get!(mp.beyond, i) do
        pm = cellmoments(mp.space, i)
        pm.area > 0 ? _meanposition(pm) : _NOWHERE
    end
end

@inline _unit(c::NTuple{3,Float64}) = (s = sqrt(c[1]^2 + c[2]^2 + c[3]^2);
    USPoint(c[1] / s, c[2] / s, c[3] / s))

"""
    _gradientstencils(src_space, ext_inds, extmap, srcmap) -> GradientStencils

The stencils of every cell in `ext_inds`, restricted to the chunk `srcmap`
names. A cell without area, or whose neighbours cannot fix a gradient, gets
no terms and keeps first order.
"""
function _gradientstencils(src_space::RegridSpace, ext_inds, extmap, srcmap)
    n = length(ext_inds)
    positions = _MeanPositions(src_space, ext_inds, extmap)
    e1s = Vector{NTuple{3,Float64}}(undef, n)
    e2s = Vector{NTuple{3,Float64}}(undef, n)
    ptr = Vector{Int}(undef, n + 1)
    cols = Int[]
    terms = NTuple{2,Float64}[]
    coeffs = NTuple{2,Float64}[]
    neighbourmeans = NTuple{3,Float64}[]
    ptr[1] = 1
    for (p, i) in enumerate(ext_inds)
        c = positions.means[p]
        e1s[p], e2s[p] = c === _NOWHERE ? (_NOWHERE, _NOWHERE) : tangentframe(_unit(c))
        ptr[p+1] = ptr[p]
        c === _NOWHERE && continue
        neighbours = cellneighbors(src_space, i)
        empty!(neighbourmeans)
        for k in neighbours
            push!(neighbourmeans, positions(k))
        end
        gradientstencil!(coeffs, e1s[p], e2s[p], c, neighbourmeans) || continue
        sx = sy = 0.0
        for (k, coeff) in zip(neighbours, coeffs)
            sx += coeff[1]
            sy += coeff[2]
            col = localindex(srcmap, k)
            col == 0 && continue
            push!(cols, col)
            push!(terms, coeff)
        end
        col = localindex(srcmap, i)
        if col != 0
            push!(cols, col)
            push!(terms, (-sx, -sy))
        end
        ptr[p+1] = length(cols) + 1
    end
    return GradientStencils(e1s, e2s, positions.means, ptr, cols, terms)
end

# The chunk's cells and their one ring: every cell whose overlaps a weight of
# this block reads. A whole space has no ring beyond it and keeps its range.
function _haloindices(src_space::RegridSpace, src_inds)
    _iswholespace(src_space, src_inds) && return src_inds
    inside = indexmap(src_inds)
    halo = Int[]
    for k in src_inds, i in cellneighbors(src_space, k)
        localindex(inside, i) == 0 && push!(halo, i)
    end
    isempty(halo) && return src_inds
    ext = vcat(collect(Int, src_inds), halo)
    return unique!(sort!(ext))
end

# Assembly

"""
    buildweights!(coo, ::ConservativeSecondOrder, dst_space, dst_inds, src_space, src_inds)

Append the second-order weights for the two chunks: each overlap's area for
its own source, plus its first moment folded through the gradient stencils of
every source cell it lies in. Coverage and denominators carry the areas alone.

Overlaps are measured for the chunk's cells and their one ring, since a
weight on a chunk cell reads the moments of the cells around it, but weights
are emitted for the chunk's own cells only.
"""
function buildweights!(coo::WeightCOO, ::ConservativeSecondOrder,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    isempty(dst_inds) && return coo
    markdenominated!(coo)
    markcovered!(coo)
    isempty(src_inds) && return coo
    _secondorderweights!(coo, dst_space, indexmap(dst_inds),
        _destinationtree(dst_space, dst_inds), _cellmemo(dst_space, dst_inds),
        src_space, src_inds)
    return coo
end

function pairblock(method::ConservativeSecondOrder, dst_space::RegridSpace,
    dst_cache::DestinationCache, src_space::RegridSpace, src_inds)
    ndst = length(dst_cache.inds)
    (ndst == 0 || isempty(src_inds)) &&
        return pairblock(method, dst_space, dst_cache.inds, src_space, src_inds)
    coo = markcovered!(markdenominated!(WeightCOO(ndst)))
    _secondorderweights!(coo, dst_space, dst_cache.map, _destinationtree(dst_cache),
        dst_cache, src_space, src_inds)
    return WeightBlock(coo, ndst, length(src_inds))
end

function pairblock(method::ConservativeSecondOrder, dst_space::RegridSpace,
    prepared::DestinationTree, src_space::RegridSpace, src_inds)
    ndst = length(prepared.inds)
    (ndst == 0 || isempty(src_inds)) &&
        return pairblock(method, dst_space, prepared.inds, src_space, src_inds)
    coo = markcovered!(markdenominated!(WeightCOO(ndst)))
    _secondorderweights!(coo, dst_space, prepared.map, _destinationtree(prepared),
        _cellmemo(dst_space, prepared.inds), src_space, src_inds)
    return WeightBlock(coo, ndst, length(src_inds))
end

function _secondorderweights!(coo::WeightCOO, dst_space::RegridSpace, dstmap,
    dst_tree, dstcells, src_space::RegridSpace, src_inds)
    m = _sharedmanifold(dst_space, src_space)
    m isa GO.Spherical || throw(ArgumentError(
        "ConservativeSecondOrder measures spherical overlaps; got manifold $m"))
    srcmap = indexmap(src_inds)
    ext_inds = _haloindices(src_space, src_inds)
    extmap = indexmap(ext_inds)
    op = BlockAreaOperator(IntersectionMomentOperator(m), dstmap, extmap,
        _cellmemo(src_space, ext_inds), dstcells)
    overlaps = _intersectionareas(m, dst_tree, subtree(src_space, ext_inds), op)
    stencils = _gradientstencils(src_space, ext_inds, extmap, srcmap)
    return _foldoverlaps!(coo, overlaps, stencils, ext_inds, srcmap)
end

# `overlaps` is destinations by extended sources. A column's own share goes to
# its source when the chunk owns it; its moments fold into every chunk column
# its stencil names, wherever the column's cell lives.
function _foldoverlaps!(coo::WeightCOO, overlaps::SparseArrays.AbstractSparseMatrixCSC,
    stencils::GradientStencils, ext_inds, srcmap)
    rows = SparseArrays.rowvals(overlaps)
    vals = SparseArrays.nonzeros(overlaps)
    @inbounds for (p, i) in enumerate(ext_inds)
        own = localindex(srcmap, i)
        t1, t2 = stencils.ptr[p], stencils.ptr[p+1] - 1
        e1, e2, c = stencils.e1[p], stencils.e2[p], stencils.mean[p]
        for t in SparseArrays.nzrange(overlaps, p)
            pm = vals[t]
            _covers(pm) || continue
            j = rows[t]
            if own != 0
                addweight!(coo, j, own, pm.area)
                addcoverage!(coo, j, own, pm.area)
                adddenom!(coo, j, pm.area)
            end
            t1 > t2 && continue
            d = (pm.moment[1] - pm.area * c[1], pm.moment[2] - pm.area * c[2],
                pm.moment[3] - pm.area * c[3])
            mx, my = tangentcoords(e1, e2, d)
            for k in t1:t2
                coeff = stencils.coeffs[k]
                addweight!(coo, j, stencils.cols[k], mx * coeff[1] + my * coeff[2])
            end
        end
    end
    return coo
end
