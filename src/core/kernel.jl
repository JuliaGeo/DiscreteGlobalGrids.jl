# ---------------------------------------------------------------------------
# Per-system operations kernel
#
# The trait functions in `interface.jl` describe *facts* about a grid system;
# the kernel generics below describe *operations*: hierarchy navigation, dense
# ordinal addressing, unit-sphere geometry, and the sorted-id pruning hook that
# the generic spatial trees (`grid_types.jl` / `generic_cursor.jl`) are built
# on. Methods are wired per system next to the native layers
# (`src/<System>/<System>Kernel.jl`); every generic here either derives from
# another kernel operation or throws `NotPortedError`.
#
# Argument order is uniformly `(system, level, id)`: `level` is the refinement
# level of the cell that `id` refers to, and any target level comes last.
# Systems with structural ids (H3, Z7) encode the level inside the id and may
# ignore the argument, but the uniform order keeps generic call sites free of
# per-system special cases.
#
# Two id models exist:
# - Ordinal ids (`has_ordinal_ids(system) == true`): the canonical ids at
#   `level` are exactly `0:num_cells(system, level) - 1` and each cell's
#   children occupy one radix-sized block (nested HEALPix). Every hierarchy
#   and ordinal operation below derives from `radix` arithmetic for free.
# - Structural ids: opaque bit encodings (H3 index, Z7). Systems wire native
#   methods for the hierarchy/ordinal group.
#
# `num_cells`, `cell_boundary`, and `cell_center` are deliberately not
# exported from the package namespace: the system submodules already export
# same-vocabulary names (`IGeo7.num_cells`, `IGeo7.cell_boundary`, ...), and
# the main module promises never to collide with those. Use them qualified.
# The full claimant list is in `src/DiscreteGlobalGrids.jl`'s export block.
# Every other generic here — `cell_polygon_unitsphere` included — is exported.
# ---------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Id-model traits
# --------------------------------------------------------------------------

"""
    cell_id_type(system::AbstractDGGS) -> Type{<:Integer}

Concrete integer type of the canonical cell id. Defaults to `Int64` (the
ordinal-id model); structural-id systems override (`UInt64` for H3 and Z7).
"""
cell_id_type(system::AbstractDGGS) = Int64

"""
    has_ordinal_ids(system::AbstractDGGS) -> Bool

`true` when the canonical ids at every `level` are exactly
`0:num_cells(system, level) - 1` *and* each cell's children occupy one
radix-sized contiguous block (nested HEALPix). Such systems get the whole
hierarchy/ordinal kernel group derived from `radix` arithmetic for free.
"""
has_ordinal_ids(system::AbstractDGGS) = false

"""
    has_descendant_ranges(system::AbstractDGGS) -> Bool

`true` when [`descendant_range`](@ref) is wired: a node's leaf-level
descendants occupy one contiguous interval in ascending canonical-id order.
Systems without it fall back to per-id `cell_parent` membership filtering in
the generic partial trees. Defaults to [`has_ordinal_ids`](@ref).
"""
has_descendant_ranges(system::AbstractDGGS) = has_ordinal_ids(system)

"""
    has_exact_subtree_cap(system::AbstractDGGS) -> Bool

`true` when [`subtree_cap`](@ref) is both O(1) and geographically tight —
the system's parent cells contain their descendants, so the cap of a node's
own cell bounds any stored subset as well as a union cap would (HEALPix).
The generic partial trees then use `subtree_cap` for internal-node extents
instead of the O(stored) [`cells_cap`](@ref), matching the old per-system
trees' O(1) node extents. Defaults to `false` — aperture-7 systems want the
tighter stored-id union caps.
"""
has_exact_subtree_cap(system::AbstractDGGS) = false

"""
    has_congruent_geometry(system::AbstractDGGS) -> Bool

`true` when every parent cell geographically contains its descendants (a
congruent refinement — nested HEALPix). That containment is what makes
parent-level geometry sound for whole subtrees: it is the fact behind
HEALPix's O(1) exact [`subtree_cap`](@ref) override, and it makes
parent-outline polygon predicates (`pred_disjoint` to prune,
`pred_covers` to bulk-accept) valid for any traversal that wants them —
wire [`subtree_polygon_unitsphere`](@ref) alongside. (The shipped query
descent measures that classification as a net loss under the current
predicate engine and does not consult it; see the design note in
`core/lookup_ops.jl`.)

Defaults to `false` — in the aperture-7 systems (H3, IGEO7, A5) children
overhang their parent, so no parent polygon bounds the subtree and every
surviving leaf must be tested exactly.
"""
has_congruent_geometry(system::AbstractDGGS) = false

# --------------------------------------------------------------------------
# Hierarchy
# --------------------------------------------------------------------------

"""
    root_ids(system::AbstractDGGS) -> AbstractVector

Canonical ids of the level-0 cells, ascending.
"""
function root_ids(system::AbstractDGGS)
    has_ordinal_ids(system) || throw(NotPortedError(system_name(system), :root_ids,
        "Wire the native root-cell enumeration."))
    return collect(cell_id_type(system), 0:(root_count(system) - 1))
end

"""
    cell_children(system, level, id) -> AbstractVector

Canonical ids of the immediate (level + 1) children of cell `(level, id)`,
ascending. Pentagon cells in aperture-7 systems return fewer than
`radix(system)` children.

The ordinal default range-checks the id, `0 <= id < num_cells(system, level)`,
as in [`descendant_range`](@ref): the radix block is arithmetically fine for an
id that is not a cell, and a caller handed `id * radix .+ (0:radix-1)` has no
way to tell that its parent never existed.
"""
function cell_children(system::AbstractDGGS, level::Integer, id)
    has_ordinal_ids(system) || throw(NotPortedError(system_name(system), :cell_children,
        "Wire the native child enumeration."))
    # Free on traversals: the generic cursor's dense descent never comes
    # through here — `_child_ids` builds this exact block inline as a
    # `UnitRange` for every `has_ordinal_ids` system, and the structural
    # systems it does route through `cell_children` answer from their own
    # wirings (`core/generic_cursor.jl`). What is left is direct callers and
    # `cell_descendants`' level-by-level expansion, which is O(subtree) anyway.
    total = num_cells(system, level)
    0 <= id < total || throw(ArgumentError(
        "$(system_name(system)) level-$(Int(level)) cell id $id is out of range 0:$(total - 1)"))
    b = radix(system)
    base = cell_id_type(system)(id) * b
    return collect(base:(base + b - 1))
end

"""
    cell_parent(system, level, id, parent_level) -> id

Canonical id of the ancestor of cell `(level, id)` at `parent_level <= level`.

`parent_level` outside `0:level` is an `ArgumentError`, and the ordinal default
range-checks the id against its *own* level too,
`0 <= id < num_cells(system, level)`, as in [`descendant_range`](@ref): the
integer division answers for any id at all, so a nonexistent one divides down
to a perfectly ordinary — and wrong — ancestor.
"""
function cell_parent(system::AbstractDGGS, level::Integer, id, parent_level::Integer)
    has_ordinal_ids(system) || throw(NotPortedError(system_name(system), :cell_parent,
        "Wire the native parent computation."))
    0 <= parent_level <= level || throw(ArgumentError("expected 0 <= parent_level <= level"))
    # Two integer comparisons, `descendant_range`'s guard verbatim. The one
    # traversal that calls `cell_parent` per stored id — the partial cursor's
    # membership filter — already runs a `searchsortedfirst` beside it, and
    # gets here only for a system that turned `has_descendant_ranges` off.
    total = num_cells(system, level)
    0 <= id < total || throw(ArgumentError(
        "$(system_name(system)) level-$(Int(level)) cell id $id is out of range 0:$(total - 1)"))
    return cell_id_type(system)(id) ÷ radix(system)^(Int(level) - Int(parent_level))
end

"""
    cell_descendants(system, level, id, leaf_level) -> AbstractVector

Canonical ids of all `leaf_level` descendants of cell `(level, id)`,
ascending. The generic fallback expands [`cell_children`](@ref) level by
level (O(subtree)); systems with a native multi-level enumeration should
override it.

`leaf_level < level` is an `ArgumentError` in every wiring — a kernel-level
argument mistake gets one error type across systems, so generic code can catch
it. Deeper invalidity (an unusable id, a `leaf_level` past the system's own
maximum) stays the system's own error type, since only the system can describe
it; `IGEO7` raises `IGeo7.InvalidZ7Error` there, for instance. On the ordinal
default the id itself is range-checked, as in [`descendant_range`](@ref).
"""
function cell_descendants(system::AbstractDGGS, level::Integer, id, leaf_level::Integer)
    level <= leaf_level || throw(ArgumentError("expected level <= leaf_level"))
    if has_descendant_ranges(system) && has_ordinal_ids(system)
        # `descendant_range` runs the id guard for this branch.
        lo, hi = descendant_range(system, level, id, leaf_level)
        return collect(lo:hi)
    end
    if has_ordinal_ids(system)
        # The same guard as `descendant_range`, for the ordinal system that
        # turned the range trait off: the radix expansion below would otherwise
        # walk a nonexistent id down into a whole subtree of nonexistent ids.
        total = num_cells(system, level)
        0 <= id < total || throw(ArgumentError(
            "$(system_name(system)) level-$(Int(level)) cell id $id is out of range 0:$(total - 1)"))
    end
    current = cell_id_type(system)[cell_id_type(system)(id)]
    for l in Int(level):(Int(leaf_level) - 1)
        next = cell_id_type(system)[]
        for cell in current
            append!(next, cell_children(system, l, cell))
        end
        current = next
    end
    return current
end

"""
    num_cells(system, level) -> Int64

Number of *valid* cells at `level`. Distinct from [`leaf_count`](@ref), which
counts logical id slots: aperture-7 systems have pentagon gaps, so e.g. H3 has
`2 + 120 * 7^level` cells in `122 * 7^level` slots. Ordinal-id systems have no
gaps, so the two coincide and the trait arithmetic is used directly.
"""
function num_cells(system::AbstractDGGS, level::Integer)
    has_ordinal_ids(system) || throw(NotPortedError(system_name(system), :num_cells,
        "Wire the native cell count."))
    return Int64(leaf_count(system, level))
end

"""
    subtree_leaf_count(system, level, id, leaf_level) -> Int64

Number of `leaf_level` descendants of cell `(level, id)`. Used by the generic
trees for `ncells`, dense leaf indexing, and the parallelization policy, so a
wired O(1) method matters; the generic fallback enumerates the subtree.

The ordinal default range-checks the id, `0 <= id < num_cells(system, level)`,
as in [`descendant_range`](@ref), and cannot return a wrapped count: its
`radix^(leaf_level - level)` is the same product [`leaf_interval`](@ref)
guards, so a `leaf_level` past [`max_level`](@ref) is an `ArgumentError` and a
system that pins no maximum gets an `OverflowError`.
"""
function subtree_leaf_count(system::AbstractDGGS, level::Integer, id, leaf_level::Integer)
    if has_ordinal_ids(system)
        # The same two comparisons as `descendant_range`, for the same reason:
        # below is pure radix arithmetic, happy to size the subtree of a cell
        # that does not exist. `Trees.ncells` asks per node, and `num_cells` is
        # `radix^level`.
        total = num_cells(system, level)
        0 <= id < total || throw(ArgumentError(
            "$(system_name(system)) level-$(Int(level)) cell id $id is out of range 0:$(total - 1)"))
        delta = Int(leaf_level) - Int(level)
        # A reversed level pair keeps the `DomainError` `^` has always raised
        # here — which error type each wiring reports for that is a separate
        # question — and, crucially, never reaches `_checked_pow`, whose empty
        # product would answer 1.
        delta >= 0 || return Int64(radix(system))^delta
        _checked_level(system, leaf_level)
        return Int64(_checked_pow(radix(system), delta))
    end
    return Int64(length(cell_descendants(system, level, id, leaf_level)))
end

# --------------------------------------------------------------------------
# Neighbors
#
# Same-level edge adjacency, the operation the lookup-level halo table and
# `stencil` (`core/lookup_ops.jl`) are built on. The neighbor count of a cell
# is a small compile-time-boundable number — 6 for the hexagonal aperture-7
# systems (5 at the 12 pentagons), 8 for HEALPix (7 at its 24 degree-3-vertex
# pixels) — so the container is a `SmallCollections.SmallVector` sized by the
# `max_neighbors` trait: fixed capacity, variable length, no heap allocation.
# --------------------------------------------------------------------------

"""
    max_neighbors(system::AbstractDGGS) -> Int

Capacity bound for [`cell_neighbors`](@ref)' container: the largest number of
edge neighbors any cell of the system has. 6 for the hexagonal aperture-7
systems (pentagons have 5), 8 for HEALPix (24 pixels per grid have 7). Wired
next to `cell_neighbors`; the fallback throws [`NotPortedError`](@ref).
"""
function max_neighbors(system::AbstractDGGS)
    throw(NotPortedError(system_name(system), :max_neighbors,
        "Wire the static neighbor-count bound next to cell_neighbors."))
end

"""
    cell_neighbors(system, level, id) -> SmallVector{max_neighbors(system),cell_id_type(system)}

Canonical ids of the cells sharing an edge with cell `(level, id)`, in
ascending canonical-id order. Hexagon cells in the aperture-7 systems have 6,
the 12 pentagons 5; HEALPix pixels have 8 (24 per grid have 7 — corner
neighbors included, following the HEALPix convention where the stencil
neighborhood is the 3×3 lattice block).

The relation is symmetric and never includes the cell itself. Directional
(ring-)ordered neighborhoods stay native — `HealpixLookups.nested_neighbors`
keeps the SW..S compass order for consumers that key on direction; this
kernel operation trades that for an order every system can promise.
"""
function cell_neighbors(system::AbstractDGGS, level::Integer, id)
    throw(NotPortedError(system_name(system), :cell_neighbors,
        "Wire the native edge-neighbor enumeration."))
end

# Ascending insertion for the wirings that collect neighbors from an unordered
# native enumeration (libh3's gridDisk, the HEALPix compass tuple): capacity is
# at most 8, so a binary-search insert into the immutable SmallVector beats
# materializing and sorting a heap vector.
@inline _insert_sorted(v::SmallVector{N,T}, x::T) where {N,T} =
    SmallCollections.insert(v, searchsortedfirst(v, x), x)

# --------------------------------------------------------------------------
# Dense ordinals
#
# The 1-based position of a cell among all valid cells of its level, in
# ascending canonical-id order. This is the leaf numbering the generic dense
# trees expose through `Trees.getcell`/`child_indices_extents`, so
# `cell_to_ordinal` must be strictly monotone in the id.
# --------------------------------------------------------------------------

"""
    cell_to_ordinal(system, level, id) -> Int

1-based position of cell `(level, id)` among all valid cells at `level` in
ascending canonical-id order. Inverse of [`ordinal_to_cell`](@ref).

The ordinal default range-checks the id, `0 <= id < num_cells(system, level)`,
as in [`descendant_range`](@ref): `id + 1` is otherwise a position in a level
that has no such position, and every consumer of this numbering (the dense
trees' leaf windows) treats it as one that does.
"""
function cell_to_ordinal(system::AbstractDGGS, level::Integer, id)
    has_ordinal_ids(system) || throw(NotPortedError(system_name(system), :cell_to_ordinal,
        "Wire the native dense-index computation."))
    # Two integer comparisons on a traversal hot path — `Trees.getcell`
    # resolves every leaf through this pair — which is what the guard has to
    # stay, and does: `num_cells` is `radix^level`.
    total = num_cells(system, level)
    0 <= id < total || throw(ArgumentError(
        "$(system_name(system)) level-$(Int(level)) cell id $id is out of range 0:$(total - 1)"))
    return Int(id) + 1
end

"""
    ordinal_to_cell(system, level, ordinal) -> id

Canonical id of the `ordinal`-th (1-based) cell at `level` in ascending
canonical-id order. Inverse of [`cell_to_ordinal`](@ref).

An `ordinal` outside `1:num_cells(system, level)` is an
[`OrdinalRangeError`](@ref) here, as it is in every wired system (A5, H3,
IGEO7): the ordinals of a level are an index space, this is the one operation
that indexes into it, and a caller who lands outside it needs to be told which
level and which range — which is exactly what a `BoundsError` against
`1:num_cells` cannot say. Unguarded, `ordinal_to_cell(HEALPixDGGS(), 0, 10^9)`
answered id `999999999` on a level whose ids are `0:11`.
"""
function ordinal_to_cell(system::AbstractDGGS, level::Integer, ordinal::Integer)
    has_ordinal_ids(system) || throw(NotPortedError(system_name(system), :ordinal_to_cell,
        "Wire the native dense-index computation."))
    # The same two comparisons as `cell_to_ordinal`, on the same hot path; the
    # message costs nothing until something prints it.
    total = num_cells(system, level)
    1 <= ordinal <= total || throw(OrdinalRangeError(
        system_name(system), Int(level), Int(ordinal), total))
    return cell_id_type(system)(ordinal - 1)
end

# --------------------------------------------------------------------------
# Pruning hook
# --------------------------------------------------------------------------

"""
    descendant_range(system, level, id, leaf_level) -> (lo, hi)

Inclusive canonical-id bounds of the `leaf_level` descendants of cell
`(level, id)`. The contract is two-sided:

- every descendant's id lies in `[lo, hi]`, and
- every *valid* `leaf_level` id in `[lo, hi]` is a descendant.

So a binary search of any ascending vector of valid leaf ids against
`[lo, hi]` yields exactly the stored descendants — invalid ids inside the
range (pentagon gaps) can never be stored. Only meaningful when
[`has_descendant_ranges`](@ref) is `true`; the generic trees check the trait
and fall back to [`cell_parent`](@ref) membership filtering otherwise.

As with [`cell_descendants`](@ref), `leaf_level < level` is an `ArgumentError`
in every wiring; system-specific invalidity keeps the system's error type.

Every wiring also rejects an `id` that is not a cell *at* `level`, because the
interval it would otherwise return is arithmetically fine and semantically
nonsense — `subtree_grid(HEALPixDGGS(), 50; root_level=0, leaf_level=2)` would
build a well-formed grid of 16 cells that do not exist. The ordinal default
checks `0 <= id < num_cells(system, level)`; the structural wirings check the
resolution the id encodes against `level` (H3, IGEO7).
"""
function descendant_range(system::AbstractDGGS, level::Integer, id, leaf_level::Integer)
    has_ordinal_ids(system) || throw(NotPortedError(system_name(system), :descendant_range,
        "Wire the sorted-id descendant interval, or leave has_descendant_ranges false."))
    level <= leaf_level || throw(ArgumentError("expected level <= leaf_level"))
    # Two integer comparisons, the ordinal counterpart of H3's one-bit-op
    # resolution guard: an ordinal id carries no level of its own, so the range
    # of the level is all there is to check it against. Negligible beside the
    # binary-search descent this sits under, and `num_cells` is `radix^level`.
    total = num_cells(system, level)
    0 <= id < total || throw(ArgumentError(
        "$(system_name(system)) level-$(Int(level)) cell id $id is out of range 0:$(total - 1)"))
    span = cell_id_type(system)(radix(system))^(Int(leaf_level) - Int(level))
    lo = cell_id_type(system)(id) * span
    return (lo, lo + span - one(span))
end

# --------------------------------------------------------------------------
# Geometry
#
# All geometry is on the unit sphere — the form `GeometryOps.Spherical()` and
# `ConservativeRegridding` consume. The cap constants and formulas below are
# lifted verbatim from the per-system tree modules they replace: leaf caps
# inflate the max center-to-vertex distance by 1.2 (covers child overhang in
# aperture-7 hierarchies, where a parent does not geographically contain its
# children), exact union caps inflate by 1.0001. The leaf-cap factor is a
# per-system trait (`cell_cap_inflation`) because the overhang is a property of
# the refinement, not of the sphere.
# --------------------------------------------------------------------------

const CELL_CAP_INFLATION = 1.2
const SUBTREE_CAP_EXACT_LIMIT = 2048

"""
    cell_cap_inflation(system::AbstractDGGS) -> Float64

Factor [`cell_cap`](@ref) inflates a cell's max center-to-vertex distance by,
so that the cell's cap also bounds its descendants where a parent does not
geographically contain its children. Defaults to `CELL_CAP_INFLATION`; a system
whose measured subtree overhang leaves too little headroom under that default
raises it (A5, whose descendants reach 1.452 of the parent radius — 1.469
extrapolated, because its pentagon lattice only *approximately* nests — against
H3's 1.052 and IGEO7's 1.048). Every value is validated empirically in the
per-system kernel test suites (`test/<System>/test_*_kernel.jl`).
"""
cell_cap_inflation(system::AbstractDGGS) = CELL_CAP_INFLATION

"""
    cell_boundary(system, level, id; closed=false) -> Vector{UnitSphericalPoint}

Boundary vertices of cell `(level, id)` on the unit sphere, in ring order.
`closed = true` repeats the first point at the end (polygon-ring form).
Systems whose raw boundaries carry projection-distortion vertices (H3) return
the cleaned ring here so every downstream consumer sees one geometry.
"""
function cell_boundary(system::AbstractDGGS, level::Integer, id; closed::Bool=false)
    throw(NotPortedError(system_name(system), :cell_boundary,
        "Wire the native cell boundary as unit-sphere points."))
end

"""
    cell_center(system, level, id) -> UnitSphericalPoint

Center of cell `(level, id)` on the unit sphere. The generic fallback
normalizes the mean of the boundary vertices; systems with a native center
should override it.
"""
function cell_center(system::AbstractDGGS, level::Integer, id)
    points = cell_boundary(system, level, id)
    center = reduce(+, points) / length(points)
    norm = sqrt(sum(abs2, center))
    norm <= eps(Float64) && throw(ArgumentError(
        "boundary vertices of cell ($level, $id) average to the origin; wire a native cell_center"))
    return GO.UnitSphericalPoint(center[1] / norm, center[2] / norm, center[3] / norm)
end

"""
    cell_polygon_unitsphere(system, level, id) -> GI.Polygon

The cell's closed boundary ring as a unit-sphere polygon. Derived from
[`cell_boundary`](@ref); this is what the generic trees hand to
`ConservativeRegridding` through `Trees.getcell`.
"""
function cell_polygon_unitsphere(system::AbstractDGGS, level::Integer, id)
    return GI.Polygon([GI.LinearRing(cell_boundary(system, level, id; closed=true))])
end

"""
    subtree_polygon_unitsphere(system, level, id, leaf_level) -> Union{Nothing, GI.Polygon}

A unit-sphere polygon that **exactly bounds** the union of the `leaf_level`
descendant cell polygons of cell `(level, id)`, or `nothing` when no such
polygon is available — the default, and the only sound answer where a parent
does not geographically contain its children (the aperture-7 systems).

The contract is geometric, not approximate: every point of every descendant's
[`cell_polygon_unitsphere`](@ref) must lie inside (or on) the returned
polygon, so a traversal may prune a subtree whose outline is disjoint from a
geometry and bulk-accept one it covers. A wiring may also answer `nothing`
selectively — HEALPix does above a densification cutoff, where the outline
stops being worth building. Only meaningful where
[`has_congruent_geometry`](@ref) holds; note the shipped query descent does
not consult it (see the design note in `core/lookup_ops.jl`).
"""
subtree_polygon_unitsphere(system::AbstractDGGS, level::Integer, id, leaf_level::Integer) =
    nothing

"""
    cell_cap(system, level, id) -> SphericalCap

Bounding cap of cell `(level, id)`: centered on [`cell_center`](@ref) with the
max center-to-vertex distance inflated by [`cell_cap_inflation`](@ref). The
inflation keeps single-cell caps usable as subtree caps in hierarchies where
children overhang their parent.
"""
function cell_cap(system::AbstractDGGS, level::Integer, id)
    center = cell_center(system, level, id)
    points = cell_boundary(system, level, id)
    radius = maximum(GO.UnitSpherical.spherical_distance(center, point) for point in points)
    return GO.UnitSpherical.SphericalCap(center,
        nextfloat(min(Float64(pi), radius * cell_cap_inflation(system) + 1e-9)))
end

"""
    cells_cap(system, level, ids) -> SphericalCap

Bounding cap of an arbitrary batch of same-level cells: the exact union cap of
their boundary vertices, or the full sphere once the batch exceeds
`SUBTREE_CAP_EXACT_LIMIT` (the exact cap stops paying for itself). Used by the
generic partial trees for internal-node extents.
"""
function cells_cap(system::AbstractDGGS, level::Integer, ids)
    isempty(ids) && return full_sphere_extent()
    length(ids) == 1 && return cell_cap(system, level, first(ids))
    length(ids) > SUBTREE_CAP_EXACT_LIMIT && return full_sphere_extent()
    points = GO.UnitSphericalPoint{Float64}[]
    sizehint!(points, 6 * length(ids))
    for id in ids
        append!(points, cell_boundary(system, level, id))
    end
    isempty(points) && return full_sphere_extent()
    center_vec = reduce(+, points) / length(points)
    center_norm = sqrt(sum(abs2, center_vec))
    center_norm <= eps(Float64) && return full_sphere_extent()
    center = GO.UnitSphericalPoint(ntuple(i -> center_vec[i] / center_norm, 3))
    radius = 0.0
    for point in points
        radius = max(radius, GO.UnitSpherical.spherical_distance(center, point))
    end
    return GO.UnitSpherical.SphericalCap(center,
        nextfloat(min(Float64(pi), radius * 1.0001 + 1e-12)))
end

"""
    intersects_cap(cap, extent) -> Bool
    intersects_cap(cap) -> predicate

Whether the spherical `cap` meets `extent`, and — curried — the predicate the
`SpatialTreeInterface` traversals take, which is what makes a cap query over a
DGGS tree a one-liner:

```julia
tree = treeify(DGGSGrid(HEALPixDGGS(), 5))
cap = cell_cap(HEALPixDGGS(), 2, 7)
indices = STI.query(tree, intersects_cap(cap))
```

Node extents on these trees are always `SphericalCap`s (see
[`cell_cap`](@ref) / [`cells_cap`](@ref)), so this is the whole predicate
vocabulary a query needs. It exists so callers stop reaching into
`GeometryOps.UnitSpherical._intersects`, which is private.
"""
intersects_cap(cap::GO.UnitSpherical.SphericalCap,
    extent::GO.UnitSpherical.SphericalCap) = GO.UnitSpherical._intersects(cap, extent)

intersects_cap(cap::GO.UnitSpherical.SphericalCap) = Base.Fix1(intersects_cap, cap)

"""
    subtree_cap(system, level, id, leaf_level) -> SphericalCap

Bounding cap of the full `leaf_level` subtree of cell `(level, id)`: an exact
union cap while the subtree has at most `SUBTREE_CAP_EXACT_LIMIT` leaves, then
the inflated [`cell_cap`](@ref) of the cell itself. Systems where a parent
geographically contains its children (HEALPix) should override with the exact
parent cap.
"""
function subtree_cap(system::AbstractDGGS, level::Integer, id, leaf_level::Integer)
    level == leaf_level && return cell_cap(system, level, id)
    if subtree_leaf_count(system, level, id, leaf_level) <= SUBTREE_CAP_EXACT_LIMIT
        return cells_cap(system, leaf_level, cell_descendants(system, level, id, leaf_level))
    end
    return cell_cap(system, level, id)
end
