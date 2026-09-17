# `gridRingUnsafe` supplies counter-clockwise shells. `gridDisk` is not a
# substitute: its order is deterministic but not rotational — at res 0, 33 of
# the 122 k = 1 sequences wind 0 or -1 turns. At pentagon distortion,
# rings are recovered from distance-bucketed disks and sorted by azimuth, then
# id. The starting direction is libh3's deterministic per-cell choice.

const MAX_NEIGHBORS = 6

"""
    neighbors(grid::H3LevelGrid, c::H3Cell, k = 1; connectivity = Vertex())

Cells within `k` grid steps, excluding `c`, as counter-clockwise shells
concatenated outward. [`ring`](@ref) is the final shell. Hexagons have six
immediate neighbours and pentagons five, without padding.

`gridRingUnsafe` supplies the order when available; pentagon-distorted rings use
azimuth order with id tie-breaking.

`k <= 1` returns a `SmallVector{6,H3Cell}` and **allocates nothing at all**,
including at pentagons.

`k >= 2` returns a `Vector{H3Cell}`, since `3k(k+1)` outgrows any static bound;
the `Val(k)` form answers the same cells in the same order from a stack buffer
while `3k(k+1)` fits [`static_capacity`](@ref).
"""
Base.@constprop :aggressive function neighbors(grid::H3LevelGrid, c::H3Cell, k::Integer=1;
        connectivity::Connectivity=Vertex())
    steps = DGG.checked_steps(k)
    steps == 0 && return SmallVector{MAX_NEIGHBORS,H3Cell}()
    steps == 1 && return one_ring(grid, c, connectivity)
    out = H3Cell[]
    append!(out, one_ring(grid, c, connectivity))
    for j in 2:steps
        append!(out, _ring_vector(c, j))
    end
    return out
end

"""
    neighborcount(grid::H3LevelGrid, c::H3Cell; connectivity = Vertex()) -> Int

Return 5 for pentagons and 6 for other cells, for either connectivity. The
count uses one libh3 pentagon test and does not construct the ring.
"""
DGG.neighborcount(grid::H3LevelGrid, c::H3Cell;
    connectivity::Connectivity=Vertex()) = ispentagon(c) ? 5 : 6

"""
    ring(grid::H3LevelGrid, c::H3Cell, k; connectivity = Vertex())

The cells at grid distance **exactly** `k` from `c`, counter-clockwise seen
from outside. `ring(grid, c, 0)` is `[c]`.

The result is the final shell returned by [`neighbors`](@ref)`(grid, c, k)`.

Uses libh3's O(k) shell walk, or an O(k²) pentagon-safe disk fallback ordered by
azimuth.
"""
Base.@constprop :aggressive function ring(grid::H3LevelGrid, c::H3Cell, k::Integer;
        connectivity::Connectivity=Vertex())
    steps = DGG.checked_steps(k)
    steps == 0 && return H3Cell[c]
    steps == 1 && return collect(one_ring(grid, c, connectivity))
    return _ring_vector(c, steps)
end

# ===========================================================================
# The shells
#
# libh3 walks its own shells, so this system implements the `one_ring` hook and
# keeps its native `neighbors`/`ring` rather than the shared breadth-first walk.
# ===========================================================================

"""
    one_ring(grid, c, connectivity) -> SmallVector{6,H3Cell}

The immediate neighbours of `c`, counter-clockwise seen from outside, starting
at libh3's deterministic per-cell direction. Allocation-free on both paths:
`gridRingUnsafe` where it applies, and an azimuth-sorted disk at a pentagon
seam.
"""
function one_ring(::H3LevelGrid, c::H3Cell, ::Connectivity)
    shell = H3Native.grid_ring_unsafe_1(c.id)
    out = SmallVector{MAX_NEIGHBORS,H3Cell}()
    if shell !== nothing
        for id in shell
            id == 0 && continue
            out = SmallCollections.push(out, H3Cell(id))
        end
        return out
    end
    # At a pentagon seam, recover the ring from the disk and sort it by azimuth.
    frame = _tangent_frame(c.id)
    keyed = SmallVector{MAX_NEIGHBORS,Tuple{Float64,H3Cell}}()
    for id in H3Native.grid_disk_1(c.id)
        (id == 0 || id == c.id) && continue
        keyed = _insert_sorted(keyed, (_azimuth(frame, id), H3Cell(id)))
    end
    for (_, cell) in keyed
        out = SmallCollections.push(out, cell)
    end
    return out
end

# The `Val` forms: `K` is a type parameter, so the ring bound `6K` and the
# disc bound `3K(K + 1)` are constants and libh3 writes each shell straight
# into a stack buffer. The order is the native one the `Integer` forms above
# answer, cell for cell; only the container changes.
function neighbors(grid::H3LevelGrid, c::H3Cell, ::Val{K};
        connectivity::Connectivity=Vertex()) where {K}
    DGG.checked_steps(K)
    K == 0 && return SmallVector{MAX_NEIGHBORS,H3Cell}()
    K == 1 && return one_ring(grid, c, connectivity)
    bound = DGG.static_capacity(DGG.Fallbacks._disc_bound(grid, Val(K),
        connectivity), H3Cell)
    bound === nothing && return neighbors(grid, c, K; connectivity)
    out = _shellbuf(bound)
    for x in one_ring(grid, c, connectivity)
        push!(out, x)
    end
    for j in 2:K
        _ring_into!(out, c, j, Val(K))
    end
    return SmallVector(out)
end

function ring(grid::H3LevelGrid, c::H3Cell, ::Val{K};
        connectivity::Connectivity=Vertex()) where {K}
    DGG.checked_steps(K)
    K == 0 && return H3Cell[c]
    K == 1 && return one_ring(grid, c, connectivity)
    cap = DGG.static_capacity(DGG.maxring(grid, K, connectivity), H3Cell)
    cap === nothing && return ring(grid, c, K; connectivity)
    out = _shellbuf(cap)
    _ring_into!(out, c, K, Val(K))
    return SmallVector(out)
end

@inline _shellbuf(::Val{N}) where {N} = SmallCollections.MutableSmallVector{N,H3Cell}()

# Append the shell at distance `k <= K` to `out`: libh3's ring walk, or at a
# pentagon distortion the distance-bucketed disk sorted by azimuth, both into
# buffers sized by `K` so they stay on the stack.
@inline function _ring_into!(out, c::H3Cell, k::Int, ::Val{K}) where {K}
    shell = H3Native.grid_ring_unsafe_static(c.id, k, Val(6K))
    if shell !== nothing
        for i in 1:6k
            id = @inbounds shell[i]
            id == 0 && continue
            push!(out, H3Cell(id))
        end
        return nothing
    end
    cells, dists = H3Native.grid_disk_distances_static(c.id, k,
        Val(3K * (K + 1) + 1))
    frame = _tangent_frame(c.id)
    keyed = SmallVector{6K,Tuple{Float64,H3Cell}}()
    for i in eachindex(cells)
        id = @inbounds cells[i]
        (id == 0 || Int(@inbounds dists[i]) != k) && continue
        keyed = _insert_sorted(keyed, (_azimuth(frame, id), H3Cell(id)))
    end
    for (_, cell) in keyed
        push!(out, cell)
    end
    return nothing
end

# Rings at k >= 2 with a run-time `k`, where the answer outgrows any static bound.
function _ring_vector(c::H3Cell, k::Int)
    shell = H3Native.grid_ring_unsafe(c.id, k)
    shell !== nothing && return [H3Cell(id) for id in shell if id != 0]
    cells, dists = H3Native.grid_disk_distances(c.id, k)
    frame = _tangent_frame(c.id)
    keyed = Tuple{Float64,H3Cell}[]
    for (id, d) in zip(cells, dists)
        (id == 0 || Int(d) != k) && continue
        push!(keyed, (_azimuth(frame, id), H3Cell(id)))
    end
    sort!(keyed)
    return [cell for (_, cell) in keyed]
end

# ===========================================================================
# The tangent frame the fallback orders by
# ===========================================================================

# A right-handed tangent basis at a cell centre. Increasing
# `atan(d.north, d.east)` is counter-clockwise when viewed from outside.
@inline function _tangent_frame(id::UInt64)
    c = H3Native.cell_center_cartesian(id)
    ex, ey = -c[2], c[1]
    n = sqrt(ex * ex + ey * ey)
    # Use a fixed east vector when longitude is undefined at a pole.
    east = n <= 1e-12 ? (1.0, 0.0, 0.0) : (ex / n, ey / n, 0.0)
    north = (c[2] * east[3] - c[3] * east[2],
             c[3] * east[1] - c[1] * east[3],
             c[1] * east[2] - c[2] * east[1])
    return c, east, north
end

@inline function _azimuth(frame, id::UInt64)
    c, east, north = frame
    p = H3Native.cell_center_cartesian(id)
    # The neighbour direction projected into the tangent plane at `c`.
    radial = p[1] * c[1] + p[2] * c[2] + p[3] * c[3]
    dx = p[1] - radial * c[1]
    dy = p[2] - radial * c[2]
    dz = p[3] - radial * c[3]
    return atan(dx * north[1] + dy * north[2] + dz * north[3],
                dx * east[1] + dy * east[2] + dz * east[3])
end

# Insert into an immutable sorted `SmallVector` without allocating.
@inline function _insert_sorted(v::SmallVector{N,T}, x::T) where {N,T}
    i = length(v)
    v = SmallCollections.push(v, x)
    while i >= 1 && x < @inbounds v[i]
        v = SmallCollections.setindex(v, @inbounds(v[i]), i + 1)
        i -= 1
    end
    return SmallCollections.setindex(v, x, i + 1)
end
