# ---------------------------------------------------------------------------
# S2 lattice neighbours and cube-face seam crossings. `SEAM` is derived from
# the chart face frames. Across an edge, coordinates undergo an exact signed
# permutation; `st_to_uv`'s exact oddness makes the same transform valid on the
# integer lattice without reprojection.
#
# `wrap_xyf` crosses one face edge. A diagonal step across a cube corner returns
# `nothing`: three cells meet there, so no fourth diagonal neighbour exists.
# Thus the 24 face-corner cells have seven Moore neighbours above level 0.
# Counter-clockwise lattice offsets remain counter-clockwise on the sphere
# because the gnomonic chart preserves orientation.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Face edges
#
# `d ∈ 1:4` names an edge of a face by the lattice direction that leaves through
# it: 1 = `+s` (`u = +1`), 2 = `+t` (`v = +1`), 3 = `-s` (`u = -1`),
# 4 = `-t` (`v = -1`). The "along-edge" axis is `v̂` for the `±s` edges and `û`
# for the `±t` edges, and the centred coordinate along it is `b = 2i + 1 - nside`
# where `i` is whichever of `ix`, `iy` is still in range.
# ---------------------------------------------------------------------------

@inline _edge_outward(face::Int, d::Int) =
    d == 1 ? FACE_U_AXIS[face + 1] :
    d == 2 ? FACE_V_AXIS[face + 1] :
    d == 3 ? map(-, FACE_U_AXIS[face + 1]) : map(-, FACE_V_AXIS[face + 1])

@inline _edge_along(face::Int, d::Int) =
    (d == 1 || d == 3) ? FACE_V_AXIS[face + 1] : FACE_U_AXIS[face + 1]

_dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]

"""
    _seam_entry(face, d) -> (g, d', σ)

The cube-edge crossing that leaves `face` through edge `d`: the face `g` on the
other side, the edge `d'` of `g` the crossing arrives at, and the sign `σ`
relating the two faces' centred along-edge lattice coordinates (`b' = σ b`).

Derived exactly from signed-axis face frames. [`SEAM`](@ref) caches all 24
entries.
"""
function _seam_entry(face::Int, d::Int)
    outward = Tuple(_edge_outward(face, d))
    along = Tuple(_edge_along(face, d))
    g = findfirst(h -> FACE_NORMAL[h + 1] == outward, 0:5) - 1
    u_g = FACE_U_AXIS[g + 1]
    v_g = FACE_V_AXIS[g + 1]
    w_f = FACE_NORMAL[face + 1]
    # `ŵ_f` decides WHICH edge of `g` the crossing lands on; the along-edge axis
    # decides the sign of the coordinate that runs along it.
    on_u = _dot3(w_f, u_g)
    if on_u != 0
        return (g, on_u == 1 ? 1 : 3, _dot3(along, v_g))
    else
        return (g, _dot3(w_f, v_g) == 1 ? 2 : 4, _dot3(along, u_g))
    end
end

"""
    SEAM[face + 1][d] -> (g, d', σ)

The cube-edge adjacency table: leaving 0-based `face` through edge `d` (1 =
`+s`, 2 = `+t`, 3 = `-s`, 4 = `-t`) arrives on face `g` at its edge `d'`, with
the centred along-edge lattice coordinate multiplied by `σ ∈ {+1, -1}`.

Built at load time by [`_seam_entry`](@ref) from the chart face frames.
"""
const SEAM = ntuple(f -> ntuple(d -> _seam_entry(f - 1, d), 4), 6)

"""
    wrap_xyf(ix, iy, face, nside) -> Union{NTuple{3,Int},Nothing}

One lattice step, with cube-edge crossing: the `(ix, iy, face)` of the cell at
face-local lattice position `(ix, iy)` of `face`, where `ix` and `iy` may each
be out of range by one (`-1` or `nside`).

  - Both in range: returned unchanged.
  - Exactly one out of range: the crossing of that face edge, through
    [`SEAM`](@ref).
  - Both out of range: `nothing`, because the step crosses a cube corner where
    no fourth cell exists.

Pure integer arithmetic; no geometry is evaluated and no tolerance is involved.
Inputs more than one step out of range are unsupported and unchecked.
"""
function wrap_xyf(ix::Integer, iy::Integer, face::Integer, nside::Integer)
    n = Int64(nside)
    x = Int64(ix)
    y = Int64(iy)
    (0 <= x < n && 0 <= y < n) && return (Int(x), Int(y), Int(face))
    outx = x < 0 ? -1 : (x >= n ? 1 : 0)
    outy = y < 0 ? -1 : (y >= n ? 1 : 0)
    (outx != 0 && outy != 0) && return nothing
    d = outx == 1 ? 1 : outx == -1 ? 3 : outy == 1 ? 2 : 4
    b = (d == 1 || d == 3) ? 2y + 1 - n : 2x + 1 - n
    g, dprime, sigma = SEAM[Int(face) + 1][d]
    k = (sigma * b + n - 1) >> 1          # the along-edge lattice index on `g`
    dprime == 1 && return (Int(n - 1), Int(k), g)
    dprime == 3 && return (0, Int(k), g)
    dprime == 2 && return (Int(k), Int(n - 1), g)
    return (Int(k), 0, g)
end

"""
    NEIGHBOR_OFFSETS

The eight lattice offsets `(dx, dy)` of the Moore neighbourhood, in
**counter-clockwise order seen from outside the sphere, starting at `+s`** —
the direction of increasing `ix`, i.e. of increasing chart coordinate `s`.

The orientation-preserving chart makes the plain `0°, 45°, …` lattice sweep
counter-clockwise on the sphere.
"""
const NEIGHBOR_OFFSETS = ((1, 0), (1, 1), (0, 1), (-1, 1),
                          (-1, 0), (-1, -1), (0, -1), (1, -1))

"""
    _neighbor_cycle(connectivity) -> Tuple

The indices of [`NEIGHBOR_OFFSETS`](@ref) to visit, in rotational order.

`Vertex()` uses all eight offsets. `Edge()` keeps the four axis offsets
`+s, +t, -s, -t`; one-component offsets share an edge and diagonal offsets
share only a vertex.
"""
_neighbor_cycle(::DGG.Vertex) = (1, 2, 3, 4, 5, 6, 7, 8)
_neighbor_cycle(::DGG.Edge) = (1, 3, 5, 7)

"""
    lattice_neighbors(ordinal, level, connectivity) -> SmallVector{8,Int64}

The scaffold ordinals adjacent to `ordinal` at `level`, in the rotational order
of [`_neighbor_cycle`](@ref): counter-clockwise seen from outside the sphere,
starting at the `+s` lattice direction.

Invalid cube-corner diagonals are omitted. Duplicate cells are removed while
retaining the first occurrence, preserving rotational order.
"""
@inline _neighbor_value(::Type{Int64}, ::Integer, h::Int64) = h
@inline _neighbor_value(::Type{DGG.LevelIndex}, level::Integer, h::Int64) =
    DGG.LevelIndex(level, h)

function _lattice_neighbors(::Type{T}, ordinal::Integer, level::Integer,
        connectivity::DGG.Connectivity) where {T}
    nside = Int64(1) << Int(level)
    ix, iy, face = hilbert_to_xyf(ordinal, nside)
    out = DGG.Helpers.empty_small_list(Val(8),
        _neighbor_value(T, level, Int64(0)))
    for m in _neighbor_cycle(connectivity)
        dx, dy = NEIGHBOR_OFFSETS[m]
        w = wrap_xyf(ix + dx, iy + dy, face, nside)
        w === nothing && continue
        h = xyf_to_hilbert(w[1], w[2], w[3], nside)
        h == ordinal && continue
        item = _neighbor_value(T, level, h)
        item in out && continue
        out = DGG.Helpers.small_push(out, item)
    end
    return SmallVector{8,T}(out)
end

lattice_neighbors(ordinal::Integer, level::Integer,
    connectivity::DGG.Connectivity) =
    _lattice_neighbors(Int64, ordinal, level, connectivity)
