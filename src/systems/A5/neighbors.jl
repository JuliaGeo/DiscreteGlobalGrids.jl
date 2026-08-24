# A5 supplies unordered adjacency sets. Shells are therefore sorted by azimuth
# counter-clockwise in a tangent frame. Ring 1 starts at the smallest id; outer
# rings use the same spoke. `edge_only=true` implements `Edge()` connectivity.

# Container capacity for one-step results under either connectivity.
const MAX_NEIGHBORS = 11

_edge_only(::Vertex) = false
_edge_only(::Edge) = true

"""
    one_ring(grid, c, connectivity) -> SmallVector{11,A5Cell}

The immediate neighbours of `c` counter-clockwise seen from outside, starting at
the smallest [`A5Cell`](@ref) id: the native adjacency set, which arrives
unordered, wound about `c`'s centroid.

Throws `ArgumentError` unless `c` is valid at the grid resolution.
"""
function one_ring(grid::LevelGrid, c::A5Cell, connectivity::Connectivity)
    level(c) == grid.level || throw(ArgumentError(
        "A5 cell $c is at resolution $(level(c)), not this grid's $(grid.level)"))
    isvalid(c) || throw(ArgumentError("A5 cell $c is not a valid cell"))
    shell = map(A5Cell,
        A5Native._get_global_cell_neighbors(c.id; edge_only=_edge_only(connectivity)))
    if length(shell) > 1
        centre = cell_centroid(grid, c)
        e1, e2, zero = DGG._ring_frame(grid, centre, minimum(shell))
        keyed = DGG.Helpers.empty_small_list(Val(MAX_NEIGHBORS), (0.0, c))
        for d in shell
            phase = DGG.Fallbacks._phase(DGG.Fallbacks._azimuth(
                centre, e1, e2, cell_centroid(grid, d)) - zero)
            keyed = DGG.Helpers.small_push(keyed, (phase, d))
        end
        shell = map(last, DGG.Helpers.small_sort(keyed))
    end
    return SmallVector{MAX_NEIGHBORS,A5Cell}(shell)
end

"""
    neighbors(grid::LevelGrid, c::A5Cell, k = 1; connectivity = Vertex())

Cells within `k` grid steps, excluding `c`, as counter-clockwise shells
concatenated outward. Ring 1 starts at the smallest [`A5Cell`](@ref) id and all
outer rings use the same spoke. [`ring`](@ref) is the final shell.

[`Vertex()`](@ref Vertex) includes corner-only neighbours; [`Edge()`](@ref Edge)
does not. See [`maxneighbors`](@ref).

`k <= 1` returns a `SmallVector{11,A5Cell}` — sized by the `Vertex()` bound
under both connectivities. `k >= 2` returns a `Vector{A5Cell}`.

Throws `ArgumentError` unless `c` is valid at the grid resolution.
"""
Base.@constprop :aggressive function neighbors(grid::LevelGrid, c::A5Cell, k::Integer=1;
        connectivity::Connectivity=Vertex())
    steps = DGG.checked_steps(k)
    steps == 0 && return SmallVector{MAX_NEIGHBORS,A5Cell}()
    steps == 1 && return one_ring(grid, c, connectivity)
    return DGG.shell_disc(grid, c, steps, connectivity)
end

"""
    ring(grid::LevelGrid, c::A5Cell, k; connectivity = Vertex())

The cells at grid distance **exactly** `k` from `c`, counter-clockwise seen from
outside. `ring(grid, c, 0)` is `[c]`.

The result is the final shell returned by [`neighbors`](@ref)`(grid, c, k)`.
"""
Base.@constprop :aggressive function ring(grid::LevelGrid, c::A5Cell, k::Integer;
        connectivity::Connectivity=Vertex())
    steps = DGG.checked_steps(k)
    steps == 0 && return A5Cell[c]
    steps == 1 && return one_ring(grid, c, connectivity)
    # An exhausted component yields an empty shell, not a missing one.
    return DGG.shell_ring(grid, c, steps, connectivity)
end
