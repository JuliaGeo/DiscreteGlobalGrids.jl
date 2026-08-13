# IGeo7Kernel.jl — wiring of the package-level operations kernel for the
# `IGEO7DGGS` singleton (see the operation contracts in `src/core/kernel.jl`).
#
# Every method here is a thin adapter over the public native API listed in
# `IGeo7.jl`'s export block: the hierarchy is Z7 prefix arithmetic, the dense
# ordinal is `cell_to_index`/`index_to_cell`, and the geometry is the Snyder
# chart's cartesian boundary ring. Two operations the kernel needs are not
# native API of their own and are built here from this module's integer layer:
#
# - `subtree_leaf_count`, O(1) from the resolution delta (`7^Δ` for a hexagon,
#   `p(Δ) = (5·7^Δ + 1)/6` for a pentagon — `grid.jl`'s `PENT_COUNT`), where
#   the generic fallback would enumerate the whole subtree;
# - `descendant_range`, the sorted-id pruning interval, by masking digit slots
#   (derivation in the block comment above it).
#
# Unqualified calls below are this module's own names; everything from the
# package namespace is reached through `DGG`. `num_cells`, `cell_boundary` and
# `cell_center` exist in both namespaces and are deliberately never conflated:
# `DGG.cell_boundary` is the unit-sphere kernel generic, `cell_boundary` the
# native `(lon, lat)` one.
#
# CLEANROOM: no sealed source was opened; only the public native API and this
# module's own clean-room layers (`z7.jl`, `grid.jl`) are used.

import ..DiscreteGlobalGrids as DGG
import GeometryOps as GO
import GeometryOpsCore as GOCore
import SmallCollections
using SmallCollections: SmallVector

# --------------------------------------------------------------------------
# Id model
#
# The Z7 `UInt64` *is* the cell id (design.md Section 3) — a structural id, so
# `has_ordinal_ids` stays `false` and the whole hierarchy/ordinal/pruning group
# is wired natively below rather than derived from radix arithmetic.
# --------------------------------------------------------------------------

DGG.cell_id_type(::DGG.IGEO7DGGS) = UInt64

# `root_count` is a verified fact now, wired next to the other IGEO7 traits in
# `src/core/systems/igeo7.jl` (`== Z7_NUM_BASES == 12`, pinned against
# `test/IGeo7/vectors/res0_cells.csv`). Nothing in this kernel reads it:
# `root_ids` is wired natively below, and the trait-derived
# `leaf_count = root_count · radix^level` counts id *slots*, `12 · 7^r`, not
# the `10 · 7^r + 2` valid cells `num_cells` reports.

# --------------------------------------------------------------------------
# Hierarchy
#
# `level` is redundant for a structural id and is only checked where a
# mismatch would silently corrupt a result (`descendant_range`).
# --------------------------------------------------------------------------

DGG.root_ids(::DGG.IGEO7DGGS) = res0_cells()

DGG.cell_children(::DGG.IGEO7DGGS, level::Integer, id) = cell_to_children(UInt64(id))

DGG.cell_parent(::DGG.IGEO7DGGS, level::Integer, id, parent_level::Integer) =
    cell_to_parent(UInt64(id), parent_level)

# The `level <= leaf_level` guard is the kernel's, not Z7's: a negative delta is
# an argument mistake generic code makes, and the kernel contract promises one
# error type for it across systems. Left to `cell_to_children` it would surface
# as `InvalidZ7Error(:descendant_res)` — right about the Z7 layer, wrong about
# whose mistake it is, and unique to this system. Everything past the guard is
# genuinely Z7 validity and keeps the native error.
function DGG.cell_descendants(::DGG.IGEO7DGGS, level::Integer, id, leaf_level::Integer)
    Int(level) <= Int(leaf_level) || throw(ArgumentError("expected level <= leaf_level"))
    return cell_to_children(UInt64(id), leaf_level)
end

DGG.num_cells(::DGG.IGEO7DGGS, level::Integer) = num_cells(level)

# Subtree size in O(1): a hexagon prefix has all seven child digits at every
# level below it (`7^Δ`), a pentagon prefix is missing its base's deleted digit
# for as long as it stays a pentagon (`PENT_COUNT[Δ+1] = (5·7^Δ + 1)/6`) —
# `grid.jl`'s tables, the same two closed forms `_subtree_count` uses to size
# `cell_to_children`'s output.
function DGG.subtree_leaf_count(::DGG.IGEO7DGGS, level::Integer, id, leaf_level::Integer)
    z = UInt64(id)
    res = get_resolution(z)                      # validates the id
    res <= leaf_level <= MAX_RESOLUTION ||
        throw(InvalidZ7Error(:descendant_res, z, _z7_int(leaf_level), res))
    delta = Int(leaf_level) - res
    return is_pentagon(z) ? (@inbounds PENT_COUNT[delta+1]) : (@inbounds POW7[delta+1])
end

# --------------------------------------------------------------------------
# Neighbors
#
# `_cell_neighbors` (grid.jl) already answers the kernel's contract — edge
# neighbors, ascending, 6 for a hexagon and 5 for a pentagon — in the native
# container (`Helpers.SmallList`, since the native core stays stdlib +
# `Helpers` only); the wiring re-seats it in the kernel's `SmallVector`.
# --------------------------------------------------------------------------

DGG.max_neighbors(::DGG.IGEO7DGGS) = 6

function DGG.cell_neighbors(::DGG.IGEO7DGGS, level::Integer, id)
    out = SmallVector{6,UInt64}()
    for neighbor in _cell_neighbors(UInt64(id))
        out = DGG._insert_sorted(out, neighbor)
    end
    return out
end

function neighbors(::DGG.IGEO7DGGS, index::IGEO7Index)
    out = SmallVector{6,IGEO7Index}()
    for id in _cell_neighbors(index.id)
        out = SmallCollections.push(out, IGEO7Index(id))
    end
    return out
end

neighbors(index::IGEO7Index) = neighbors(DGG.IGEO7DGGS(), index)

@inline function Base.in(
        index::IGEO7Index,
        indices::IGEO7Indices{<:DGG.DGGSSubtreeIds{DGG.IGEO7DGGS}},
    )
    get_resolution(index) == indices.resolution || return false
    return index.id in indices.ids
end

@inline function Base.in(
        index::IGEO7Index,
        indices::IGEO7Indices{<:DGG.DGGSGlobeIds{DGG.IGEO7DGGS}},
    )
    return get_resolution(index) == indices.resolution
end

function edges(
        indices::IGEO7Indices{<:DGG.DGGSSubtreeIds{DGG.IGEO7DGGS}},
    )
    tile = indices.ids
    return IGEO7Index.(border_descendants(tile.root_id, tile.level))
end

edges(::IGEO7Indices{<:DGG.DGGSGlobeIds{DGG.IGEO7DGGS}}) = IGEO7Index[]

function _complete_subtree(indices::IGEO7Indices)
    isempty(indices) && return nothing
    if indices.contiguous && indices.first_ordinal == 1 &&
            length(indices) == num_cells(indices.resolution)
        return :globe
    end
    lo, hi = first(indices).id, last(indices).id
    z7_base_cell(lo) == z7_base_cell(hi) || return nothing
    root_level = 0
    for level in 1:indices.resolution
        _z7_digit(lo, level) == _z7_digit(hi, level) || break
        root_level = level
    end
    root = z7_parent(lo, root_level)
    DGG.subtree_leaf_count(
        DGG.IGEO7DGGS(), root_level, root, indices.resolution) == length(indices) ||
        return nothing
    range = DGG.descendant_range(
        DGG.IGEO7DGGS(), root_level, root, indices.resolution)
    return range == (lo, hi) ? (root_level, root) : nothing
end

function edges(indices::IGEO7Indices)
    isempty(indices) && return IGEO7Index[]
    subtree = _complete_subtree(indices)
    subtree === :globe && return IGEO7Index[]
    if subtree !== nothing
        _, root = subtree
        return IGEO7Index.(border_descendants(root, indices.resolution))
    end
    out = IGEO7Index[]
    for index in indices
        any(neighbor -> neighbor ∉ indices, neighbors(index)) && push!(out, index)
    end
    return out
end

# The subtree rim without a single neighbor query: `border_descendants`
# (grid.jl) decides membership from the Z7 digits. The fallback's pruned
# frontier is already O(result) in cells, so what this drops is the constant —
# the `grid_disk`-shaped neighbor sweep and the ancestor walk it runs per
# candidate per level — not an exponent. The
# `level <= leaf_level` guard is the kernel's, for the reason `cell_descendants`
# gives above; everything past it is Z7 validity and keeps the native error.
function DGG.subtree_border(::DGG.IGEO7DGGS, level::Integer, id, leaf_level::Integer)
    Int(level) <= Int(leaf_level) || throw(ArgumentError("expected level <= leaf_level"))
    return border_descendants(UInt64(id), leaf_level)
end

# --------------------------------------------------------------------------
# Dense ordinals
#
# `cell_to_index` is already the kernel's contract: 1-based rank among the
# cells of the id's own resolution in ascending id order, strictly monotone in
# the id (ascending `UInt64` order is (base, digit-string) lexicographic order,
# and the rank walk sums whole sibling subtrees in that order).
# --------------------------------------------------------------------------

DGG.cell_to_ordinal(::DGG.IGEO7DGGS, level::Integer, id) = cell_to_index(UInt64(id))

# The range check is deliberately written twice. `index_to_cell` runs it
# natively and throws a `BoundsError`, which its docstring pins as part of the
# native API shape and which direct native callers still get. The kernel
# promises the opposite thing — one error type for an out-of-range ordinal
# across every system, `DGG.OrdinalRangeError`, naming the system, the level
# and the level's ordinal range (`src/core/kernel.jl`) — so the wiring checks
# the range itself rather than relaying the native error. Two integer
# comparisons and a `num_cells` table lookup ahead of an O(res) walk.
function DGG.ordinal_to_cell(::DGG.IGEO7DGGS, level::Integer, ordinal::Integer)
    total = num_cells(level)                     # native; also validates `level`
    1 <= ordinal <= total || throw(DGG.OrdinalRangeError(
        DGG.system_name(DGG.IGEO7DGGS()), Int(level), Int(ordinal), Int(total)))
    return index_to_cell(ordinal, level)
end

# --------------------------------------------------------------------------
# Pruning hook
#
# A res-`R` id's slots `R+1:20` are all the padding sentinel `7`, identical for
# every cell of that resolution, so comparing two res-`R` ids as `UInt64`s
# compares `(base, d_1 … d_R)` lexicographically. The res-`R` descendants of a
# res-`r` cell are exactly the ids sharing its `(base, d_1 … d_r)` prefix, so
# they occupy one contiguous interval whose endpoints are the digit-extremal
# members: all-`0` in slots `r+1:R` for `lo`, all-`6` for `hi` (digit 7 is
# padding, never a child; digit 0 and digit 6 are never the deleted pentagon
# digit — that is 2 or 5 — so both endpoints are always *valid* cells).
#
# Both are one mask of the id itself, which already carries 7s in every slot
# below `r`:
#
#   slots  = tail(r) ⊻ tail(R)          bits of digit slots r+1 … R
#   lo     = z & ~slots                 those slots -> 0, slots R+1… stay 7
#   hi     = z & ~(slots & SLOT_LSB)    those slots 7 -> 6, slots R+1… stay 7
#
# TRAP (design review): `hi` must keep the 7-padding below `R`. Writing 6s into
# `r+1:R` and then zero-padding below `R` yields an id *smaller* than the true
# all-6 descendant, i.e. one-sided under-coverage that no random containment
# test would notice — the tight-endpoint check in `test_igeo7_kernel.jl`
# (`extrema(cell_descendants) == descendant_range`) is what pins it.
#
# The two-sided contract follows: any res-`R` id in `[lo, hi]` must share the
# prefix (a differing prefix sorts outside on one side or the other), hence is
# a descendant; pentagon gaps inside the interval are simply invalid ids that
# no lookup can hold.
# --------------------------------------------------------------------------

# Verified by the checklist in `test/IGeo7/test_igeo7_kernel.jl` (tight
# endpoints, ordered/disjoint sibling ranges, subtree sums == num_cells, deep
# endpoint validity up to res 19) before this was switched on.
DGG.has_descendant_ranges(::DGG.IGEO7DGGS) = true

function DGG.descendant_range(::DGG.IGEO7DGGS, level::Integer, id, leaf_level::Integer)
    # Same kernel-level guard as `cell_descendants` above, and for the same
    # reason: a negative delta is one error type on every system.
    Int(level) <= Int(leaf_level) || throw(ArgumentError("expected level <= leaf_level"))
    z = UInt64(id)
    res = get_resolution(z)                      # validates the id
    # Cheap cursor-bookkeeping check: a mismatched `level` means the caller
    # lost track of which level the id belongs to, and every bound below would
    # be silently wrong.
    res == Int(level) || throw(ArgumentError(
        "IGEO7 cell $(z7_to_hex(z; prefix=true)) is at resolution $res, not level $level"))
    res <= leaf_level <= MAX_RESOLUTION ||
        throw(InvalidZ7Error(:descendant_res, z, _z7_int(leaf_level), res))
    slots = _z7_tail_mask(res) ⊻ _z7_tail_mask(leaf_level)
    return (z & ~slots, z & ~(slots & Z7_SLOT_LSB))
end

# --------------------------------------------------------------------------
# Geometry
#
# The native ring is already the unit sphere in cartesian form, so the wiring
# is a point-type change; `cell_center` goes through the native `(lon, lat)`
# center (exactly the icosahedron vertex for a pentagon) rather than the
# kernel's boundary-mean fallback.
# --------------------------------------------------------------------------

function DGG.cell_boundary(::DGG.IGEO7DGGS, level::Integer, id; closed::Bool=false)
    ring = cell_boundary_cartesian(UInt64(id); closed_ring=closed)
    return [GO.UnitSphericalPoint(p[1], p[2], p[3]) for p in ring]
end

function DGG.cell_center(::DGG.IGEO7DGGS, level::Integer, id)
    lon, lat = cell_center(UInt64(id))
    return GO.UnitSphericalPoint(lonlat_to_xyz(lon, lat))
end

# `cell_cap` / `cells_cap` / `subtree_cap` keep the kernel definitions: the
# `CELL_CAP_INFLATION = 1.2` leaf cap does cover a whole ISEA7H subtree. The
# CAP-VALIDATION test measured the union ratio (max descendant-vertex distance
# from the cap center over the cell's own max vertex distance) exhaustively for
# deltas 1-5 at res 0-1 and over pentagon neighborhoods at res 4/6: worst
# 1.0482, i.e. every descendant sits at 87% of the wired cap radius, and the
# ratio converges geometrically. Convergence is two-step here, not per-delta:
# the Eisenstein lattice's chirality alternates with resolution parity
# (`_encode_lattice`'s odd/even branch), so the subtree bulges in alternating
# directions and the increments alternate large/small while each *same-parity*
# increment shrinks by 4x or more. See `test/IGeo7/test_igeo7_kernel.jl`.

# --------------------------------------------------------------------------
# Convenience constructors
# --------------------------------------------------------------------------

# Partial grid over a lookup's cells. Id validity at `l.resolution` is trusted
# from the lookup (`IGeo7Lookup` checks it on construction when built with
# `validate=true`); the generic constructor re-checks only eltype and strict
# ordering, never cell validity. `kwargs` reach `DGGSPartialGrid`'s
# `bucket_size` / `root_level` / `root_id`.
DGG.DGGSPartialGrid(l::IGeo7Lookups.IGeo7Lookup; kwargs...) =
    DGG.DGGSPartialGrid(DGG.IGEO7DGGS(), l.resolution, l.data; kwargs...)

# What the generic lookup operations (`neighbor_indices`, `stencil`, `zonal`)
# ask of a lookup: which system, which level.
DGG.dggs_system(::IGeo7Lookups.IGeo7Lookup) = DGG.IGEO7DGGS()
DGG.dggs_level(l::IGeo7Lookups.IGeo7Lookup) = l.resolution

# Treeifying a lookup directly is the shortest path from a `DimensionalData`
# dimension to a `Regridder`, and it needs nothing from this file: the method
# is generic over `AbstractDGGSLookup` (`core/lookups.jl`), routing a stored id
# vector through the constructor above and a `DGGSGlobeIds` through the dense
# `DGGSGrid` instead. All that is per-system is the manifold, which is what
# makes the one-argument `treeify(l)` resolve at all.
GOCore.best_manifold(::IGeo7Lookups.IGeo7Lookup) = GO.Spherical()

# --------------------------------------------------------------------------
# Tile stencils
#
# The subtree neighbor stepper: edge adjacency inside one tile from the Z7
# digits alone, replacing the geometric `_cell_neighbors` on the path that
# needs it most. Included last because it needs both `DGG` and this module's
# engine layer.
# --------------------------------------------------------------------------

include("tile_neighbors.jl")
