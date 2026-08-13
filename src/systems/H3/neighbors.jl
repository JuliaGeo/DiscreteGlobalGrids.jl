# ---------------------------------------------------------------------------
# Adjacency
#
# `gridDisk` is the pentagon-safe enumeration: it walks the hex lattice and
# leaves a zero in any slot the walk could not fill, which is what happens
# around the twelve pentagons where a disk has fewer cells than `3k(k+1) + 1`.
# `gridRingUnsafe` is the O(k) hollow-ring walk and is much cheaper, but it
# refuses outright near a pentagon — hence the two paths in `ring`.
#
# Connectivity is accepted and has no effect: on a hexagonal/pentagonal
# tessellation, sharing a vertex and sharing an edge are the same relation
# (three cells meet at every vertex, and any two of those three already share an
# edge), so `Vertex()` and `Edge()` name the same neighbours. See
# `max_neighbors`.
# ---------------------------------------------------------------------------

const MAX_NEIGHBORS = 6

"""
    neighbors(grid::H3Grid, c::H3Cell, k = 1; connectivity = Vertex())

The cells within `k` grid steps of `c`, excluding `c`, **ascending by
[`H3Cell`](@ref) order**.

libh3 promises no order from `gridDisk`, so this imposes one — the canonical id
order the rest of the package sorts and binary-searches in. Deterministic order
is part of the interface contract, and "whatever the C library happened to
emit" is not one.

A hexagon has six neighbours and a pentagon five; around a pentagon the disk is
genuinely smaller rather than zero-padded, which is what the "absent, never
padded" rule means here.

# Container

`k <= 1` returns a fixed-capacity `SmallVector{6,H3Cell}` and allocates
nothing, which is what makes a whole-grid neighbour sweep free of garbage.
Larger `k` returns a `Vector{H3Cell}`, since `3k(k+1)` outgrows any static
bound.
"""
function neighbors(grid::H3Grid, c::H3Cell, k::Integer=1;
        connectivity::Connectivity=Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return SmallVector{MAX_NEIGHBORS,H3Cell}()
    steps == 1 && return _neighbors1(c)
    out = H3Cell[]
    for id in H3Native.grid_disk(c.id, steps)
        (id == 0 || id == c.id) && continue
        push!(out, H3Cell(id))
    end
    return sort!(out)
end

# The k = 1 disk: at most seven slots, so the whole answer fits in the static
# bound and never reaches the heap.
function _neighbors1(c::H3Cell)
    out = SmallVector{MAX_NEIGHBORS,H3Cell}()
    for id in H3Native.grid_disk(c.id, 1)
        (id == 0 || id == c.id) && continue
        out = SmallCollections.push(out, H3Cell(id))
    end
    return sort(out)
end

"""
    ring(grid::H3Grid, c::H3Cell, k; connectivity = Vertex())

The cells at grid distance **exactly** `k` from `c`, ascending.
`ring(grid, c, 0)` is `[c]`.

Two paths, same answer:

  - `gridRingUnsafe`, an O(k) walk around the shell, whenever libh3 can
    complete it;
  - `gridDiskDistances` — O(k²), the whole disk with each cell's distance —
    when it cannot, which is precisely when a pentagon lies within `k` steps and
    the hex-lattice walk has nowhere to turn.

The fallback is not an approximation of the fast path; it is the same set,
computed the expensive way, which is why the pentagon seams do not need a
special case anywhere above this function.
"""
function ring(grid::H3Grid, c::H3Cell, k::Integer;
        connectivity::Connectivity=Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return H3Cell[c]
    shell = H3Native.grid_ring_unsafe(c.id, steps)
    if shell === nothing
        # Pentagon distortion: fall back to the disk, keeping the cells libh3
        # reports at exactly this distance.
        cells, dists = H3Native.grid_disk_distances(c.id, steps)
        out = H3Cell[]
        for (id, d) in zip(cells, dists)
            (id == 0 || Int(d) != steps) && continue
            push!(out, H3Cell(id))
        end
        return sort!(out)
    end
    out = H3Cell[]
    for id in shell
        id == 0 && continue
        push!(out, H3Cell(id))
    end
    return sort!(out)
end
