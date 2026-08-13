# ---------------------------------------------------------------------------
# The 3x3 lattice neighbourhood on the cube, and what happens at a seam
#
# Within one face the neighbourhood of a cell is the eight lattice offsets
# around it. Across a face boundary it is a cube-edge crossing, and this file is
# the arithmetic of that crossing: a small table [`SEAM`](@ref), DERIVED from
# the face frames of `chart.jl` at load time rather than transcribed, and the
# `wrap_xyf` lattice step that reads it.
#
# ## Why the crossing is a signed permutation, and why it is exact in `(s, t)`
#
# Face `f` carries the orthonormal frame `(û, v̂, ŵ)`, and a point of its chart
# is the direction `ŵ + u û + v v̂`. Consider the `u = +1` edge. Every point on
# it is `ŵ_f + û_f + β v̂_f` for `β ∈ [-1, 1]`. The direction `û_f` is one of the
# six signed axes, so it IS the outward normal `ŵ_g` of exactly one other face
# `g` — that is the face across the seam, and it is found by table lookup on
# `FACE_NORMAL` with no geometry at all. Rewriting the same point in `g`'s frame,
#
#     ŵ_f + û_f + β v̂_f  =  ŵ_g + (ŵ_f + β v̂_f)
#
# and both `ŵ_f` and `v̂_f` are orthogonal to `ŵ_g` (`= û_f`), hence lie in
# `span(û_g, v̂_g)`. Each is a signed axis, so each is `±û_g` or `±v̂_g` — and
# they cannot be the same one, being orthogonal to each other. So exactly one of
# two things happens:
#
#   * `ŵ_f = ±û_g`  ⟹  `u' = ±1` (the crossing lands on `g`'s `u = ±1` edge)
#     and `v' = σ β` with `σ = ±1`;
#   * `ŵ_f = ±v̂_g`  ⟹  `v' = ±1` and `u' = σ β`.
#
# That is: **(target face, target edge, sign)**, which is exactly the three
# integers `SEAM[f + 1][d]` holds.
#
# The last step is why the crossing is integer lattice arithmetic rather than a
# re-projection. The correspondence above is stated in `(u, v)`, but the lattice
# lives in `(s, t)`, and `st_to_uv` is odd about `s = 1/2` **exactly** in
# floating point (`st_to_uv(1 - s) == -st_to_uv(s)`). So a sign flip in `u` is a
# reflection `s ↦ 1 - s` in `s`, and the same signed permutation carries the
# lattice: writing the centred lattice coordinate `b = 2i + 1 - nside` (which
# runs over `|b| < nside` with the parity of `nside + 1` — odd at every level
# above 0, and just `b = 0` at `nside = 1`), the crossing sends `b ↦ σ b` and
# puts the result in the row or column of `g` adjacent to the target edge. No
# `uv_to_st`, no rounding, no `nside`-dependent tolerance.
#
# ## What this buys, and what it costs
#
# `wrap_xyf` handles a crossing of ONE face edge. A step that would cross two at
# once — the diagonal offset at a cell in a face CORNER — returns `nothing`, and
# that is correct rather than a limitation: three faces meet at a cube corner,
# each contributing one cell, so the corner cell's diagonal "neighbour" across
# the corner is not a fourth cell. It is one of the two cells already reached by
# the two single-edge crossings, and dropping it is what leaves those 24 cells
# (four per face, at every level ≥ 1) with **seven** Moore neighbours instead of
# eight. `test/systems/S2/runtests.jl` checks the whole answer against shared
# boundary geometry rather than against this argument.
#
# ## Rotational order
#
# The interface's `neighbors` order is ROTATIONAL: counter-clockwise seen from
# OUTSIDE the sphere, from a documented start. Unlike HEALPix — whose compass
# tuple turned out to run the wrong way — the S2 lattice needs no reversal, and
# the reason is the same one that makes `cell_corners` counter-clockwise: the
# gnomonic chart is orientation-PRESERVING from `(s, t)` onto the sphere seen
# from outside (see the argument in `chart.jl`). So the offsets taken
# counter-clockwise in the `(dx, dy)` plane are counter-clockwise on the sphere.
#
# Measured, not merely argued: the azimuths of the neighbour centres about the
# cell centre wrap exactly once — the conformance harness's own winding test —
# for every cell of every level swept in the suite, on face interiors, across
# seams, and at cube corners, under both connectivities.
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

Derived from the face frames by the argument in this file's header — three dot
products of signed axes, so the result is exact integer data with no geometry
evaluated. [`SEAM`](@ref) caches all 24 of them.
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

Built once at load time by [`_seam_entry`](@ref) from `chart.jl`'s face frames,
so it cannot drift out of step with the charts the geometry is evaluated
through. It is a lookup table only in the sense that it is *precomputed*: it is
24 entries of integer data derived from 6 orthonormal frames, not a transcribed
constant anyone has to trust.
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
  - **Both** out of range: `nothing`. That step would cross a cube CORNER, where
    three faces meet and there is no fourth cell to name — see this file's
    header. It is why the 24 cells in face corners have seven Moore neighbours.

Pure integer arithmetic; no geometry is evaluated and no tolerance is involved.
Inputs further out of range than one step are not supported and are not checked:
this is the one-ring step, and multi-ring neighbourhoods are built by iterating
it.
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

Counter-clockwise in `(dx, dy)` *is* counter-clockwise on the sphere here,
because the gnomonic chart preserves orientation; see this file's header and
the corner-winding argument in `chart.jl`. The cycle is therefore the plain
`0°, 45°, 90°, …` sweep of the lattice plane, with no reversal.
"""
const NEIGHBOR_OFFSETS = ((1, 0), (1, 1), (0, 1), (-1, 1),
                          (-1, 0), (-1, -1), (0, -1), (1, -1))

"""
    _neighbor_cycle(connectivity) -> Tuple

The positions of [`NEIGHBOR_OFFSETS`](@ref) to visit, in rotational order.

`Vertex()` is all eight — `+s, +s+t, +t, -s+t, -s, -s-t, -t, +s-t`. `Edge()`
keeps the four axis offsets `+s, +t, -s, -t` in the same cyclic order, which
stays counter-clockwise for free: restricting a cycle cannot change its winding.

An S2 cell is an axis-aligned rectangle in its face chart, so — unlike HEALPix,
whose pixel is a diamond ROTATED 45° against its lattice — the labels mean what
they look like: a one-component offset shares a whole cell edge, a two-component
offset shares only the single lattice corner between them. That holds across a
cube seam too, and the suite checks it geometrically rather than by assertion.
"""
_neighbor_cycle(::DGG.Vertex) = (1, 2, 3, 4, 5, 6, 7, 8)
_neighbor_cycle(::DGG.Edge) = (1, 3, 5, 7)

"""
    lattice_neighbors(ordinal, level, connectivity) -> Vector{Int64}

The scaffold ordinals adjacent to `ordinal` at `level`, in the rotational order
of [`_neighbor_cycle`](@ref): counter-clockwise seen from outside the sphere,
starting at the `+s` lattice direction.

Steps that leave the cube (the corner diagonals at the 24 face-corner cells) are
dropped from the cycle rather than leaving a hole, and repeats are dropped
keeping the FIRST occurrence — which is what preserves the cycle. Repeats
happen only at level 0, where `nside == 1` makes a face a single cell and all
four in-cycle steps land on four distinct faces anyway; the guard is there
because a cheap `in` test on a length-≤8 vector costs nothing and a silent
duplicate would violate the interface contract.
"""
function lattice_neighbors(ordinal::Integer, level::Integer, connectivity::DGG.Connectivity)
    nside = Int64(1) << Int(level)
    ix, iy, face = hilbert_to_xyf(ordinal, nside)
    out = Int64[]
    sizehint!(out, 8)
    for m in _neighbor_cycle(connectivity)
        dx, dy = NEIGHBOR_OFFSETS[m]
        w = wrap_xyf(ix + dx, iy + dy, face, nside)
        w === nothing && continue
        h = xyf_to_hilbert(w[1], w[2], w[3], nside)
        (h == ordinal || h in out) && continue
        push!(out, h)
    end
    return out
end
