# Cell adjacency and least-squares gradient recovery over cell means.

# --------------------------------------------------------------------------
# Generic cell adjacency
# --------------------------------------------------------------------------

"""
    cellneighbors(space::RegridSpace, i::Int) -> Vector{Int}

The geometric fallback: the cells whose polygon shares a vertex with cell `i`'s,
which on a lattice is the eight-cell ring, diagonals included.

For a once-per-plan sweep only; a space that knows its topology should answer
from that. Each call queries [`celltree`](@ref) with cell `i`'s cap — one
[`getcell`](@ref) per candidate, `O(candidates × vertices²)` comparisons, and a
fresh cell tree if the space packs one per query.

Vertices match within `1e-6` of cell `i`'s cap radius, floored so a degenerate
cap still compares something. The test measures chord, which at that scale is
the shorter of chord and arc, so it never loosens.
"""
function cellneighbors(space::RegridSpace, i::Int)
    1 <= i <= ncells(space) || throw(BoundsError(space, i))
    cap = _packedcellcap(space, i)
    ring = _cellvertices(space, i)
    tol = max(1e-6 * Float64(cap.radius), 1e-12)
    tol2 = tol * tol
    out = Int[]
    STI.depth_first_search(Base.Fix1(Extents.intersects, cap), celltree(space)) do j
        j != i && _sharesvertex(ring, _cellvertices(space, j), tol2) && push!(out, j)
        return nothing
    end
    # Ascending and once each, whatever order the tree reported its leaves in.
    sort!(out)
    unique!(out)
    return out
end

# One cell's ring as plain 3-vectors, the closing repeat included: carrying the
# duplicate costs one comparison and dropping it costs a branch per vertex.
function _cellvertices(space::RegridSpace, i::Int)
    out = NTuple{3,Float64}[]
    for p in GI.getpoint(getcell(space, i))
        push!(out, (Float64(GI.x(p)), Float64(GI.y(p)), Float64(GI.z(p))))
    end
    return out
end

function _sharesvertex(a::Vector{NTuple{3,Float64}}, b::Vector{NTuple{3,Float64}},
        tol2::Float64)
    for p in a, q in b
        d1, d2, d3 = p[1] - q[1], p[2] - q[2], p[3] - q[3]
        d1 * d1 + d2 * d2 + d3 * d3 <= tol2 && return true
    end
    return false
end

"""
    celldiameter(space::RegridSpace) -> Float64

The geometric fallback: twice the widest leaf cap of [`celltree`](@ref), capped
at `pi`. A leaf cap covers its cell, so twice its radius covers the diameter.

For plan construction only — it visits every leaf, `O(ncells(space))`.
"""
celldiameter(space::RegridSpace) =
    min(Float64(pi), 2 * _widestleafcap(celltree(space)))

function _widestleafcap(node)
    r = 0.0
    if STI.isleaf(node)
        for (_, extent) in STI.child_indices_extents(node)
            r = max(r, _capradius(extent))
        end
    else
        for child in STI.getchild(node)
            r = max(r, _widestleafcap(child))
        end
    end
    return r
end

# The cell-tree contract: every node extent is a spherical cap.
_capradius(cap::SphericalCap) = Float64(cap.radius)
_capradius(extent) = throw(ArgumentError(
    "a cell tree's node extents must be SphericalCaps; got $(typeof(extent))"))

# --------------------------------------------------------------------------
# Tangent frames
# --------------------------------------------------------------------------

# `_point3`, `_dot3`, `_cross3` and `_framefirst` come from barycentric.jl, so a
# gradient is fitted in the same frame a chart is drawn in.

"""
    tangentframe(n::USPoint) -> (e1, e2)

A right-handed orthonormal basis of the tangent plane at the unit vector `n`,
with `e1 × e2 == n`.

`n` alone fixes the frame — the axis least parallel to `n` names `e1`, never
close to degenerate — and must be a unit vector, so normalize a mean position
first.
"""
@inline function tangentframe(n::USPoint)
    u = _point3(n)
    e1 = _framefirst(u)
    return (e1, _cross3(u, e1))
end

"""
    tangentcoords(e1, e2, v) -> (Float64, Float64)

`v`'s projection onto the tangent plane, written in the frame `(e1, e2)`. `v` is
any 3-vector; the difference of two mean positions is the usual argument.
"""
@inline function tangentcoords(e1::NTuple{3,Float64}, e2::NTuple{3,Float64}, v)
    w = _point3(v)
    return (_dot3(e1, w), _dot3(e2, w))
end

# --------------------------------------------------------------------------
# Least-squares gradient recovery
# --------------------------------------------------------------------------

"""
    gradientstencil!(coeffs, e1, e2, c, neighbours) -> Bool

Fill `coeffs` with the coefficients that recover a cell's tangent gradient from
its neighbours' cell means, and return whether the stencil holds one.

`c` and `neighbours` are mean positions as 3-vectors of any length; `(e1, e2)`
is the frame the gradient comes back in ([`tangentframe`](@ref)); `coeffs` is
resized to `length(neighbours)`, so one buffer serves a whole sweep.

    g = sum(coeffs[k] .* (f[k] - f_c) for k in 1:length(neighbours))

is the unweighted least-squares fit of `Δf_k ≈ g ⋅ d_k` over the tangent offsets
`d_k = tangentcoords(e1, e2, neighbours[k] - c)`, whose normal equations
`S g = Σ d_k Δf_k` with `S = Σ d_k d_kᵀ` give `coeffs[k] = S⁻¹ d_k`.

The caller applies the differences and the self coefficient `-sum(coeffs)`,
which is what recovers a zero gradient from a constant field exactly.

Returns `false` with `coeffs` empty, leaving the caller a zero gradient, when
there are fewer than two neighbours or `det(S) <= 1e-10 * (S₁₁ + S₂₂)^2` —
neighbours strung out along one line through the cell.
"""
function gradientstencil!(coeffs::Vector{NTuple{2,Float64}},
        e1::NTuple{3,Float64}, e2::NTuple{3,Float64}, c,
        neighbours::AbstractVector)
    n = length(neighbours)
    if n < 2
        empty!(coeffs)
        return false
    end
    resize!(coeffs, n)
    c3 = _point3(c)
    s11 = s12 = s22 = 0.0
    for k in 1:n
        p = _point3(neighbours[k])
        d = tangentcoords(e1, e2, (p[1] - c3[1], p[2] - c3[2], p[3] - c3[3]))
        @inbounds coeffs[k] = d
        s11 += d[1] * d[1]
        s12 += d[1] * d[2]
        s22 += d[2] * d[2]
    end
    det = s11 * s22 - s12 * s12
    trace = s11 + s22
    # Negated so a NaN determinant fails here rather than inverting a matrix
    # nothing in the geometry produced.
    if !(det > 1e-10 * trace * trace)
        empty!(coeffs)
        return false
    end
    for k in 1:n
        @inbounds d1, d2 = coeffs[k]
        @inbounds coeffs[k] = ((s22 * d1 - s12 * d2) / det,
            (s11 * d2 - s12 * d1) / det)
    end
    return true
end

"""
    gradientstencil(e1, e2, c, neighbours) -> Union{Nothing,Vector{NTuple{2,Float64}}}

The allocating form of [`gradientstencil!`](@ref): the coefficients, or
`nothing` where that returns `false`. A sweep should reuse one buffer through
the mutating form.
"""
function gradientstencil(e1::NTuple{3,Float64}, e2::NTuple{3,Float64}, c,
        neighbours::AbstractVector)
    coeffs = Vector{NTuple{2,Float64}}(undef, length(neighbours))
    return gradientstencil!(coeffs, e1, e2, c, neighbours) ? coeffs : nothing
end
