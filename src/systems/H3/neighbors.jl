# `gridRingUnsafe` supplies counter-clockwise shells. `gridDisk` is not a
# substitute: its order is deterministic but not rotational — at res 0, 33 of
# the 122 k = 1 sequences wind 0 or -1 turns. At pentagon distortion,
# rings are recovered from distance-bucketed disks and sorted by azimuth, then
# id. The starting direction is libh3's deterministic per-cell choice.

const MAX_NEIGHBORS = 6

"""
    neighbors(grid::LevelGrid, c::H3Cell, k = 1; connectivity = Vertex())

Cells within `k` grid steps, excluding `c`, as counter-clockwise shells
concatenated outward. [`ring`](@ref) is the final shell. Hexagons have six
immediate neighbours and pentagons five, without padding.

`gridRingUnsafe` supplies the order when available; pentagon-distorted rings use
azimuth order with id tie-breaking.

`k <= 1` returns a `SmallVector{6,H3Cell}` and **allocates nothing at all**,
including at pentagons.

`k >= 2` returns a `Vector{H3Cell}`, since `3k(k+1)` outgrows any static bound.
"""
function neighbors(grid::LevelGrid, c::H3Cell, k::Integer=1;
        connectivity::Connectivity=Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return SmallVector{MAX_NEIGHBORS,H3Cell}()
    steps == 1 && return _ring1(c)
    out = H3Cell[]
    append!(out, _ring1(c))
    for j in 2:steps
        append!(out, _ring_vector(c, j))
    end
    return out
end

"""
    ring(grid::LevelGrid, c::H3Cell, k; connectivity = Vertex())

The cells at grid distance **exactly** `k` from `c`, counter-clockwise seen
from outside. `ring(grid, c, 0)` is `[c]`.

Identical, element for element, to the last `length` entries of
[`neighbors`](@ref)`(grid, c, k)` — the two share one implementation, so the
disc really is its shells concatenated rather than merely agreeing as a set.

Uses libh3's O(k) shell walk, or an O(k²) pentagon-safe disk fallback ordered by
azimuth.
"""
function ring(grid::LevelGrid, c::H3Cell, k::Integer;
        connectivity::Connectivity=Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return H3Cell[c]
    steps == 1 && return collect(_ring1(c))
    return _ring_vector(c, steps)
end

# ===========================================================================
# The shells
# ===========================================================================

# The k = 1 ring, allocation-free on both paths.
function _ring1(c::H3Cell)
    shell = H3Native.grid_ring_unsafe_1(c.id)
    out = SmallVector{MAX_NEIGHBORS,H3Cell}()
    if shell !== nothing
        for id in shell
            id == 0 && continue
            out = SmallCollections.push(out, H3Cell(id))
        end
        return out
    end
    # Pentagon seam: recover the ring from the stack-buffer disk and impose the
    # rotational order ourselves.
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

# Rings at k >= 2, where the answer outgrows any static bound.
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

# A right-handed tangent basis at a cell centre: `east x north` is the outward
# normal, so increasing `atan(d.north, d.east)` is counter-clockwise seen from
# OUTSIDE the sphere — the same handedness `cell_boundary` winds in.
@inline function _tangent_frame(id::UInt64)
    c = H3Native.cell_center_cartesian(id)
    ex, ey = -c[2], c[1]
    n = sqrt(ex * ex + ey * ey)
    # No H3 cell centre lands exactly on a pole, but a frame that silently
    # produced NaNs there would corrupt an order rather than fail.
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

# Insert into a sorted `SmallVector`, keeping it immutable and allocation-free.
#
# `SmallCollections.sort` is deliberately not used anywhere in this file: it
# returns a `MutableSmallVector`, which is both a different concrete type from
# the `k = 0` answer and a heap allocation.
@inline function _insert_sorted(v::SmallVector{N,T}, x::T) where {N,T}
    i = length(v)
    v = SmallCollections.push(v, x)
    while i >= 1 && x < @inbounds v[i]
        v = SmallCollections.setindex(v, @inbounds(v[i]), i + 1)
        i -= 1
    end
    return SmallCollections.setindex(v, x, i + 1)
end
