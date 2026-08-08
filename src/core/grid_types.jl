# ---------------------------------------------------------------------------
# Generic DGGS grid types
#
# One grid family for every kernel-wired system, replacing the per-system
# full-grid / partial-grid clones this package used to carry. A grid is a value
# describing *which* cells exist; `Trees.treeify` turns it into the spatial-tree
# cursor defined in `generic_cursor.jl`.
#
# Construction validates eagerly exactly what it can do without a wired kernel:
# the `max_level` bound on `level`, and — for `DGGSPartialGrid` — element type,
# strict ascending id order and root membership. `DGGSGrid` deliberately does
# NOT call `num_cells`: registered-but-unwired systems
# (`DGGSGrid(RHEALPixDGGS(), 40)`) must still construct, so a level the id
# encoding can hold but no uniform grid occupies (`DGGSGrid(A5DGGS(), 30)`,
# A5's res-30 encoding gap) constructs here and throws on the first operation
# that needs the cell count.
#
# `DGGSPartialGrid` calls it in one place only, and only for
# `has_ordinal_ids` systems: there the count is `root_count * radix^level`,
# pure trait arithmetic that needs no wiring, and it is what turns the sorted
# ids' two endpoints into a complete range check. Anything else would have to
# read the ids, which is the lookup layer's job, not this one's.
# ---------------------------------------------------------------------------

"""
    DGGSGrid(system, level; manifold=GO.Spherical())

The complete grid of `system` at `level` (`num_cells(system, level)` cells),
addressed by dense ordinal — see [`cell_to_ordinal`](@ref).

`manifold` declares the surface the grid is referenced to; see
[`grid_manifold`](@ref) for what it does and does not change.
"""
struct DGGSGrid{S<:AbstractDGGS,M<:GOCore.Manifold}
    system::S
    level::Int
    manifold::M
    function DGGSGrid(system::AbstractDGGS, level::Integer;
            manifold::GOCore.Manifold=GO.Spherical())
        lvl = Int(level)
        lvl >= 0 || throw(ArgumentError("level must be non-negative"))
        limit = max_level(system)
        limit === nothing || lvl <= limit ||
            throw(ArgumentError("$(system_name(system)) level must be in 0:$limit"))
        # Rejects `Planar` / `AutoManifold` here rather than at the first
        # traversal: a grid whose frame nobody declared is exactly what the
        # manifold field exists to make impossible.
        authalic_sphere(manifold)
        return new{typeof(system),typeof(manifold)}(system, lvl, manifold)
    end
end

"""
    DGGSPartialGrid(system, level, ids; bucket_size=0)

A subset of one level's cells, addressed by position in `ids` — which must be
strictly ascending canonical ids (the order every `<X>Lookup` stores, so
`Trees.getcell` indices line up with lookup/data order). `bucket_size > 0`
stops tree descent once a node covers that few stored leaves.

The `root_level`/`root_id` fields carry an optional subtree rooting (see
[`subtree_grid`](@ref)); `root_level == -1` means the whole sphere.

# What is checked

Always: `level` against [`max_level`](@ref), `eltype(ids)` against
[`cell_id_type`](@ref), strict ascending order, and — for a hand-rooted grid —
that the ids really are descendants of `(root_level, root_id)` (complete and
O(1) where [`has_descendant_ranges`](@ref) holds, an endpoint-only guard
otherwise).

Cell-id *validity* splits by id model, and the split is worth knowing because
the two halves are not equally covered:

  * ordinal-id systems (`has_ordinal_ids`, i.e. HEALPix here) are checked
    completely: the ids at `level` are exactly `0:num_cells - 1`, and since the
    vector is already known sorted, `first(ids) >= 0 && last(ids) < num_cells`
    bounds every id in two comparisons. `DGGSPartialGrid(HEALPixDGGS(), 0,
    Int64[100, 200])` is an `ArgumentError`, not a grid that explodes inside
    chart math halfway through a traversal.
  * structural-id systems (H3, IGEO7, A5) are not checked here at all — an id
    is an opaque bit encoding whose validity only the system can judge, and
    doing it per id is the expensive check their lookups keep behind
    `validate=true`. Their `<X>Lookup` constructors do check that the ids are
    sorted and that the *endpoints* encode the declared resolution; interior
    ids are trusted unless the lookup was built with `validate=true`. An
    invalid interior id therefore surfaces at the first kernel call that
    decodes it, as that system's own error type.

# From a lookup

Every kernel-wired system also defines a one-argument constructor over its
`DimensionalData` lookup — `DGGSPartialGrid(l)` for a `HealpixLookup`,
`H3Lookup`, `IGeo7Lookup` or `A5Lookup` — which reads the system and the level
off the lookup and forwards its `kwargs` here (`bucket_size` / `root_level` /
`root_id`). Everything above still runs on that path; nothing is skipped
because the ids arrived from a lookup.

The lookup's id vector is passed through, never copied or reordered, so
`grid.ids === l.data`: leaf index `i` of the resulting tree is position `i` of
the lookup, hence of any `DimArray` built on that dimension. That is what makes
a `ConservativeRegridding.Regridder` line up with the data array without a
permutation. `Trees.treeify(GO.Spherical(), l)` is the same path in one step —
except over a globe-complete lookup (`l.data isa DGGSGlobeIds`), which it sends
to the dense [`DGGSGrid`](@ref) cursor instead. This constructor stays open on
one, for a caller who explicitly wants the partial cursor over a globe: every
check above is O(1) on a `DGGSGlobeIds`, since strict ascent is answered from
the type and the range and root-membership checks read only the two endpoints.
"""
struct DGGSPartialGrid{S<:AbstractDGGS,V<:AbstractVector,ID,M<:GOCore.Manifold}
    system::S
    level::Int
    ids::V
    bucket_size::Int
    root_level::Int
    root_id::ID
    manifold::M
    # Validation lives in the inner constructor so no construction path — the
    # cursor's internals included — can skip the invariants the cursor assumes.
    function DGGSPartialGrid(system::AbstractDGGS, level::Integer, ids::AbstractVector,
            bucket_size::Integer, root_level::Integer, root_id,
            manifold::GOCore.Manifold)
        lvl = Int(level)
        lvl >= 0 || throw(ArgumentError("level must be non-negative"))
        limit = max_level(system)
        limit === nothing || lvl <= limit ||
            throw(ArgumentError("$(system_name(system)) level must be in 0:$limit"))
        bucket_size >= 0 || throw(ArgumentError("bucket_size must be non-negative"))
        -1 <= Int(root_level) <= lvl || throw(ArgumentError("expected -1 <= root_level <= level"))
        eltype(ids) === cell_id_type(system) || throw(ArgumentError(
            "ids must be a $(cell_id_type(system)) vector for $(system_name(system)), got eltype $(eltype(ids))"))
        # `Helpers.strictly_increasing`, which this used to spell as
        # `issorted(ids; lt=(<=))` and which is the same predicate, because it
        # is the one a `DGGSGlobeIds` answers from its type
        # (`globe_ids.jl`): the ordinal contract already guarantees the ascent,
        # and walking it would be `num_cells` kernel calls — 4.8e9 on a res-9 H3
        # globe. One function, so the lookups' identical check and this one are
        # short-circuited together.
        Helpers.strictly_increasing(ids) ||
            throw(ArgumentError("ids must be strictly ascending"))
        if has_ordinal_ids(system) && !isempty(ids)
            # The ids at `lvl` are exactly `0:total - 1` by the definition of
            # the trait, and the vector is sorted by the line above, so its two
            # endpoints bound every entry: this is O(1) *and* complete. Without
            # it an out-of-range id builds a structurally valid grid and throws
            # only later, from inside the chart math of a traversal.
            total = num_cells(system, lvl)
            (first(ids) >= 0 && last(ids) < total) || throw(ArgumentError(
                "$(system_name(system)) level-$lvl cell ids must be in 0:$(total - 1), got $(first(ids)):$(last(ids))"))
        end
        rid = cell_id_type(system)(root_id)
        if Int(root_level) >= 0 && !isempty(ids)
            # A hand-rooted grid whose ids stray outside the root's subtree
            # would silently hide those leaves from traversal (still addressable
            # through `Trees.getcell`, never reached under the root), so
            # membership is checked here. With a descendant interval the sorted
            # endpoints bound every id, making the check complete and O(1);
            # otherwise the endpoint parent checks are a partial guard.
            if has_descendant_ranges(system)
                lo, hi = descendant_range(system, Int(root_level), rid, lvl)
                (lo <= first(ids) && last(ids) <= hi) || throw(ArgumentError(
                    "ids must all be descendants of the root cell ($(Int(root_level)), $rid)"))
            else
                (cell_parent(system, lvl, first(ids), Int(root_level)) == rid &&
                 cell_parent(system, lvl, last(ids), Int(root_level)) == rid) || throw(ArgumentError(
                    "ids must all be descendants of the root cell ($(Int(root_level)), $rid)"))
            end
        end
        authalic_sphere(manifold)  # see the note in `DGGSGrid`
        return new{typeof(system),typeof(ids),typeof(rid),typeof(manifold)}(
            system, lvl, ids, Int(bucket_size), Int(root_level), rid, manifold)
    end
end

function DGGSPartialGrid(system::AbstractDGGS, level::Integer, ids::AbstractVector;
        bucket_size::Integer=0, root_level::Integer=-1, root_id=zero(cell_id_type(system)),
        manifold::GOCore.Manifold=GO.Spherical())
    return DGGSPartialGrid(system, level, ids, bucket_size, root_level, root_id, manifold)
end

"""
    subtree_grid(system, root_id; root_level, leaf_level, bucket_size=0) -> DGGSPartialGrid

The subtree of cell `(root_level, root_id)` at `leaf_level`, as a partial grid
over its materialized descendants — O(subtree), never O(globe). Leaves are
numbered 1:subtree_leaf_count in ascending id order, so a `Regridder` built on
the resulting tree lines up with per-chunk data arrays.

```julia
chunk = subtree_grid(H3DGGS(), H3.H3Native.lonlat_to_cell(10.0, 45.0, 2);
                     root_level=2, leaf_level=6)
tree = treeify(chunk)   # 2401 leaves, indexed 1:2401
```

Both levels are keywords on purpose. The two of them are bare `Integer`s that
would otherwise straddle the id positionally, and a system whose ids are also
integers (every one of them) cannot tell a swapped call from a correct one —
it would just build a tree of the wrong subtree. There is no positional form.
"""
function subtree_grid(system::AbstractDGGS, root_id;
        root_level::Integer, leaf_level::Integer, bucket_size::Integer=0,
        manifold::GOCore.Manifold=GO.Spherical())
    0 <= Int(root_level) <= Int(leaf_level) ||
        throw(ArgumentError("expected 0 <= root_level <= leaf_level"))
    ids = cell_descendants(system, root_level, root_id, leaf_level)
    return DGGSPartialGrid(system, leaf_level, ids; bucket_size,
        root_level=Int(root_level), root_id, manifold)
end

"""
    grid_manifold(grid) -> GeometryOpsCore.Manifold

The manifold a grid declares — `Spherical` for a plain spherical grid,
`Geodesic` for one referenced to an ellipsoid.

This is the *declared* frame, not the compute frame. `GOCore.best_manifold`
returns [`authalic_sphere(grid_manifold(grid))`](@ref) instead, because the
tessellation lives on the authalic sphere and because `ConservativeRegridding`
dispatches only on `Planar`/`Spherical` — a `Geodesic` compute manifold is a
`MethodError`.

What the declared manifold changes today: it is carried into the grid's type,
so grids on different ellipsoids are different types; it fixes the authalic
radius `best_manifold` reports; and it is the ellipsoid any serialization
(CF `earth_radius` / `semi_major_axis`, the Zarr DGGS `ellipsoid` object,
`dggh.parameters.ellipsoid`) must declare. What it does *not* change: the unit
sphere the charts, boundaries and caps are computed on. Nothing in the tree
layer is warped by it.

!!! note "Value mismatch is not a type mismatch"
    Two grids on different ellipsoids are different types, so dispatch sees
    them. But two `Geodesic`s with *different parameters* share one type, so a
    WGS84 grid and a GRS80 grid are type-identical. Anything that combines two
    grids has to compare the manifolds by value.
"""
grid_manifold(grid::DGGSGrid) = grid.manifold
grid_manifold(grid::DGGSPartialGrid) = grid.manifold
grid_manifold(cursor) = grid_manifold(cursor.grid)
