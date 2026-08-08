# ---------------------------------------------------------------------------
# H3 operations-kernel wiring
#
# H3 is the reference *structural-id* system of `core/kernel.jl`: the level
# lives inside the 64-bit index, ids are not dense, and pentagon subtrees are
# missing one child each — so every operation is wired natively here (through
# `H3Native`, i.e. the H3 C library) rather than derived from `radix`.
#
# Three pieces are not plain `H3Native` calls but algorithms carried over from
# the old per-system H3 tree module, which this kernel replaced and which no
# longer exists — the kernel is the only tree layer now:
# - the cumulative base-cell table behind `cell_to_ordinal`/`ordinal_to_cell`,
# - the distortion-vertex removal pipeline behind `cell_boundary`, and
# - the native center used by `cell_cap`.
#
# `descendant_range` is digit arithmetic on the index layout; see the comment
# above it. It makes H3 a `has_descendant_ranges` system, which is what lets
# the generic partial cursor prune with `searchsorted` instead of per-id
# `cell_parent` filtering.
# ---------------------------------------------------------------------------

import ..DiscreteGlobalGrids as DGG
import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI

# --------------------------------------------------------------------------
# Id-model traits
# --------------------------------------------------------------------------

DGG.cell_id_type(::DGG.H3DGGS) = UInt64

# Verified in `test/H3/test_h3_kernel.jl`: tight `descendant_range`
# endpoints over complete parent levels, disjoint ordered sibling ranges, and
# subtree sizes summing to `num_cells` — jointly the two-sided range contract.
DGG.has_descendant_ranges(::DGG.H3DGGS) = true

# --------------------------------------------------------------------------
# Hierarchy
# --------------------------------------------------------------------------

# `getRes0Cells` already returns base cells 0:121 in order, and base-cell bits
# (45-51) outrank every digit, so the vector is ascending as `root_ids` wants —
# asserted in the tests, never sorted here. Cached because `ordinal_to_cell`
# indexes it per call; `root_ids` hands out copies so no caller can corrupt it.
const _H3_ROOT_IDS = H3Native.res0_cells()

DGG.root_ids(::DGG.H3DGGS) = copy(_H3_ROOT_IDS)

DGG.cell_children(::DGG.H3DGGS, level::Integer, id) =
    H3Native.cell_to_children(id, Int(level) + 1)

DGG.cell_parent(::DGG.H3DGGS, level::Integer, id, parent_level::Integer) =
    H3Native.cell_to_parent(id, parent_level)

# `cellToChildren` enumerates any level gap directly (and in ascending index
# order), so the kernel's level-by-level expansion is never used for H3.
function DGG.cell_descendants(::DGG.H3DGGS, level::Integer, id, leaf_level::Integer)
    Int(level) <= Int(leaf_level) || throw(ArgumentError("expected level <= leaf_level"))
    return H3Native.cell_to_children(id, leaf_level)
end

DGG.num_cells(::DGG.H3DGGS, level::Integer) = H3Native.num_cells(level)

# O(1): `cellToChildrenSize` is closed-form digit arithmetic, pentagon gaps
# included. This is what `ncells` and the parallelize policy call per node.
DGG.subtree_leaf_count(::DGG.H3DGGS, level::Integer, id, leaf_level::Integer) =
    H3Native.cell_to_children_size(id, leaf_level)

# --------------------------------------------------------------------------
# Dense ordinals
#
# Cells are numbered base cell by base cell, and within a base cell by child
# position (the old per-system H3 tree's dense indexing), so
# one cumulative table per level turns the ordinal into a single binary search
# plus an O(1) H3 position lookup. `cell_to_child_pos` is monotone in the id
# (it is the digit path read as a number, pentagon gaps skipped), so the
# resulting numbering is monotone in the canonical id as the kernel requires.
# --------------------------------------------------------------------------

const _H3_ROOT_ENDS = ntuple(H3Native.MAX_RESOLUTION + 1) do resolution_index
    counts = [H3Native.cell_to_children_size(root, resolution_index - 1) for root in _H3_ROOT_IDS]
    return cumsum(counts)
end

function _checked_int(n::Integer)
    n <= typemax(Int) || throw(OverflowError("H3 count $n does not fit in Int"))
    return Int(n)
end

# `level` (not the id's own resolution field) picks the cumulative table, so a
# caller that mislabels a cell gets a wrong ordinal rather than a silent
# cross-level answer; `descendant_range` is where that mismatch is caught.
function DGG.cell_to_ordinal(::DGG.H3DGGS, level::Integer, id)
    cell = UInt64(id)
    root_index = H3Native.get_base_cell(cell) + 1
    root_ends = _H3_ROOT_ENDS[Int(level) + 1]
    previous_end = root_index == 1 ? Int64(0) : root_ends[root_index - 1]
    return _checked_int(previous_end + H3Native.cell_to_child_pos(cell, 0) + 1)
end

# An ordinal outside `1:num_cells` is a `DGG.OrdinalRangeError`, the one type
# every wiring of this operation raises (`src/core/kernel.jl`): the caller is
# told which system's level it overran and what that level's ordinals are,
# which the `BoundsError` this used to throw could not say.
function DGG.ordinal_to_cell(::DGG.H3DGGS, level::Integer, ordinal::Integer)
    total = H3Native.num_cells(level)
    1 <= ordinal <= total || throw(DGG.OrdinalRangeError(
        DGG.system_name(DGG.H3DGGS()), Int(level), Int(ordinal), _checked_int(total)))
    root_ends = _H3_ROOT_ENDS[Int(level) + 1]
    root_index = searchsortedfirst(root_ends, ordinal)
    previous_end = root_index == 1 ? Int64(0) : root_ends[root_index - 1]
    return H3Native.child_pos_to_cell(ordinal - previous_end - 1, _H3_ROOT_IDS[root_index], level)
end

# --------------------------------------------------------------------------
# Pruning hook
#
# H3 index layout: the resolution field is bits 52-55, and digit `k`
# (`k in 1:15`) is the three bits whose low bit sits at `45 - 3k`. Every valid
# index 7-pads the digit slots below its own resolution, and digits are bounded
# by 6, so among the `leaf_level` descendants of a cell the smallest is the
# all-0 digit path and the largest the all-6 one. Base-cell and higher digits
# are shared by the whole subtree and outrank the varying slots, so the two
# endpoints bracket exactly the descendants: any id between them carries the
# same prefix, and any digit 7 inside the varying window belongs to no valid
# `leaf_level` cell. That is the two-sided contract.
# --------------------------------------------------------------------------

const _H3_RESOLUTION_MASK = UInt64(0x0f) << 52

_h3_digit_shift(k::Int) = 45 - 3k
_h3_resolution(id::UInt64) = Int((id >> 52) & 0x0f)

function DGG.descendant_range(::DGG.H3DGGS, level::Integer, id, leaf_level::Integer)
    cell = UInt64(id)
    lvl = Int(level)
    leaf = Int(leaf_level)
    lvl <= leaf || throw(ArgumentError("expected level <= leaf_level"))
    0 <= leaf <= H3Native.MAX_RESOLUTION ||
        throw(ArgumentError("H3 resolution must be in 0:$(H3Native.MAX_RESOLUTION)"))
    # One bit-op guard: a cursor that hands over a cell whose real resolution is
    # not `level` would otherwise get a silently wrong (and unsound) interval.
    _h3_resolution(cell) == lvl || throw(ArgumentError(
        "H3 cell 0x$(string(cell; base=16)) has resolution $(_h3_resolution(cell)), not level $lvl"))
    digits_mask = UInt64(0)
    low_bit_mask = UInt64(0)
    for k in (lvl + 1):leaf
        shift = _h3_digit_shift(k)
        digits_mask |= UInt64(7) << shift
        low_bit_mask |= UInt64(1) << shift
    end
    base = (cell & ~_H3_RESOLUTION_MASK) | (UInt64(leaf) << 52)
    # `lo` clears the new digit slots (the all-center descendant). `hi` clears
    # only their low bit (7 -> 6, the largest valid digit) and leaves the slots
    # below `leaf_level` at their 7 padding: writing 6s and then zero-padding
    # below `leaf_level` would put `hi` *under* the true maximum descendant and
    # silently drop the tail of every subtree.
    return (base & ~digits_mask, base & ~low_bit_mask)
end

# --------------------------------------------------------------------------
# Geometry
#
# `cellToBoundary` inserts extra vertices where a cell crosses an icosahedron
# face boundary (up to 10 for a hexagon); they are projection distortion, not
# cell corners, and they make the ring non-convex, which the spherical clipper
# cannot handle. The cleanup drops
# vertices one at a time, preferring the removal that restores geodesic
# convexity and perturbs the spherical area least, until 6 (5 for a pentagon)
# remain. Everything downstream — `cell_polygon_unitsphere`, `cell_cap`,
# `cells_cap` — sees only the cleaned ring.
# --------------------------------------------------------------------------

_unit_point(xyz) = GO.UnitSphericalPoint(xyz[1], xyz[2], xyz[3])

_polygon_from_points(points) = GI.Polygon([GI.LinearRing(points)])

function _closed_ring(points)
    isempty(points) && return copy(points)
    isclosed = points[end] == points[1]
    ring = Vector{eltype(points)}(undef, length(points) + !isclosed)
    copyto!(ring, points)
    isclosed || (ring[end] = ring[1])
    return ring
end

function _is_geodesically_convex(points)
    length(points) < 3 && return false
    for i in eachindex(points)
        edge_start = points[i]
        edge_end = points[mod1(i + 1, length(points))]
        for point in points
            GO.UnitSpherical.spherical_orient(edge_start, edge_end, point) < 0 && return false
        end
    end
    return true
end

function _convexity_violations(points)
    count = 0
    for i in eachindex(points)
        edge_start = points[i]
        edge_end = points[mod1(i + 1, length(points))]
        for point in points
            count += GO.UnitSpherical.spherical_orient(edge_start, edge_end, point) < 0
        end
    end
    return count
end

_spherical_area(points) = GO.area(GO.Spherical(), _polygon_from_points(_closed_ring(points)))

_without_index(points, i::Integer) = [points[j] for j in eachindex(points) if j != i]

function _remove_distortion_vertices(cell_id::UInt64, points)
    target = H3Native.is_pentagon(cell_id) ? 5 : 6
    length(points) <= target && return points

    result = copy(points)
    target_area = _spherical_area(result)
    while length(result) > target
        candidates = map(eachindex(result)) do i
            candidate = _without_index(result, i)
            convex = _is_geodesically_convex(candidate)
            violations = _convexity_violations(candidate)
            area_error = abs(_spherical_area(candidate) - target_area)
            return (; i, candidate, convex, violations, area_error)
        end
        sort!(candidates; by=c -> (c.convex ? 0 : 1, c.violations, c.area_error))
        result = candidates[1].candidate
    end
    return result
end

function DGG.cell_boundary(::DGG.H3DGGS, level::Integer, id; closed::Bool=false)
    cell = UInt64(id)
    points = _remove_distortion_vertices(
        cell, _unit_point.(H3Native.cell_boundary_cartesian(cell; closed_ring=false)))
    return closed ? _closed_ring(points) : points
end

# Native cell center (`cellToLatLng`), not the kernel's boundary mean: it is the
# center the H3 hierarchy is built around, and the one `cell_cap` measures from.
function DGG.cell_center(::DGG.H3DGGS, level::Integer, id)
    lon, lat = H3Native.cell_center(id)
    λ = deg2rad(lon)
    φ = deg2rad(lat)
    cosφ = cos(φ)
    return GO.UnitSphericalPoint(cosφ * cos(λ), cosφ * sin(λ), sin(φ))
end

# `cell_cap` is the kernel default (native center + max vertex distance inflated
# by `CELL_CAP_INFLATION`), which reproduces the old per-system H3 tree's cell
# extent exactly. The 1.2 inflation is what covers descendant overhang — H3
# parents do not contain their children geographically — and
# `test/H3/test_h3_kernel.jl` measures the actual overhang against it.

# --------------------------------------------------------------------------
# Lookup convenience
# --------------------------------------------------------------------------

"""
    DGGSPartialGrid(l::H3Lookup; kwargs...)

The lookup's stored cells as a generic partial grid. Cell-id validity at
`l.resolution` is trusted from the lookup (which checks it when constructed
with `validate=true`); the generic constructor only re-checks ordering and
element type. `kwargs` reach `DGGSPartialGrid`'s `bucket_size` / `root_level` /
`root_id`.
"""
DGG.DGGSPartialGrid(l::H3Lookups.H3Lookup; kwargs...) =
    DGG.DGGSPartialGrid(DGG.H3DGGS(), l.resolution, l.data; kwargs...)

# Treeifying a lookup directly is the shortest path from a `DimensionalData`
# dimension to a `Regridder`, and it needs nothing from this file: the method
# is generic over `AbstractDGGSLookup` (`core/lookups.jl`), routing a stored id
# vector through the constructor above and a `DGGSGlobeIds` through the dense
# `DGGSGrid` instead. All that is per-system is the manifold, which is what
# makes the one-argument `treeify(l)` resolve at all.
GOCore.best_manifold(::H3Lookups.H3Lookup) = GO.Spherical()
