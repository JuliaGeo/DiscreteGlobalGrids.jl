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
# Two more are digit arithmetic on the index layout rather than libh3 calls:
# `descendant_range` (see the comment above it), which makes H3 a
# `has_descendant_ranges` system and so lets the generic partial cursor prune
# with `searchsorted` instead of per-id `cell_parent` filtering; and
# `subtree_border`, which reads the subtree's rim off the digits instead of
# sweeping the subtree with `gridDisk`.
# ---------------------------------------------------------------------------

import ..DiscreteGlobalGrids as DGG
import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
using SmallCollections: SmallVector

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
# Neighbors
#
# `gridDisk(k = 1)` is the origin plus its edge neighbors, zero-padded where a
# pentagon truncates the disk — the pentagon-safe enumeration (`gridRingUnsafe`
# fails outright near pentagons). Dropping the origin and the padding leaves
# exactly the 6 neighbors (5 for a pentagon); libh3 promises no order, so the
# wiring sorts into the kernel's ascending contract.
# --------------------------------------------------------------------------

DGG.max_neighbors(::DGG.H3DGGS) = 6

function DGG.cell_neighbors(::DGG.H3DGGS, level::Integer, id)
    cell = UInt64(id)
    out = SmallVector{6,UInt64}()
    for neighbor in H3Native.grid_disk(cell, 1)
        (neighbor == 0 || neighbor == cell) && continue
        out = DGG._insert_sorted(out, neighbor)
    end
    return out
end

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
# Subtree border
#
# Which `leaf_level` descendants of a cell share an edge with a cell outside
# its subtree, from the index digits alone — no `gridDisk`, no geometry, no
# pass over the subtree. The rim of a depth-`d` hexagon subtree holds
# `3^(d+1) - 3` cells against the `7^d` the subtree holds, and this enumerates
# only the rim, so what the generic fallback (`src/core/kernel.jl`) spends as a
# `cell_neighbors` sweep per candidate per level becomes O(result).
#
# The mechanism is the one `src/IGeo7/grid.jl` derives for IGEO7 under its own
# `--- subtree border ---`, and it transfers because H3 refines the same way:
# aperture 7, a center child plus six around it, addressed by the seven IJK
# digits. What that buys, in the same two steps:
#
# - *Adjacency lemma.* A child of `x` and any edge neighbor of that child lie in
#   `x` or in one of `x`'s six edge neighbors — the displacement between two
#   cells of adjacent parents is at most one aperture-7 step. So the rim is
#   decided level by level, locally, and a rim cell descends only from a rim
#   cell: the search prunes.
# - *Center child.* Digit 0 sits at the middle of its parent's seven, enclosed
#   by its own siblings, so it is never on the rim — and no rim suffix contains
#   a zero digit at *any* position, which is why the digit loop below starts at
#   1 rather than 0.
#
# What a rim cell carries is *which* of the six directions it is exposed in, and
# those sets are always contiguous arcs. A state is therefore the arc
# `{s, s+1, ..., s+L-1}` (mod 6) of `_H3_DIGIT_DIR` positions, with `L` the
# number of rim children the cell has — 2, 3 or 4 — plus `L = 6` for the subtree
# root, which is exposed all round. `_h3_border_step` is that arc's transition.
#
# The direction ring is the geometric cyclic order of the IJK unit vectors,
# `K(1), JK(3), J(2), IJ(6), I(4), IK(5)`, inverted into digit order as
# `_H3_DIGIT_DIR`. It is the digit *alphabet*, fixed for every resolution; the
# per-resolution rotation of the grid lives in H3's Class II / Class III
# alternation, not in the digits.
#
# That alternation is why the table is level-parity dependent: consecutive H3
# resolutions are rotated against each other by `acos(sqrt(3/7))` ~ 19.1 deg,
# in alternating senses, exactly as IGEO7's refinement alternates chirality
# (`chi_is_c`, `src/IGeo7/engine.jl`). Hence the two mirror-image branches
# below, selected on the *child's* absolute resolution.
#
# **The parity roles are swapped relative to IGEO7.** The branch H3 takes at an
# even child level is the one IGEO7 takes at an odd one, and vice versa. That is
# not a transcription slip: the two systems anchor their alternation at
# different ends (H3 res 0 is Class II; IGEO7's `chi_is_c` turns `+ALPHA` at
# even levels), and the offsets `t + 1` / `t - 1` / `t - 2` follow whichever
# sense the refinement turns in. A reader comparing the two files should expect
# them to differ exactly here and nowhere else.
#
# The census follows from the table as it does for IGEO7: an arc of length `L`
# has `L` rim children, of which one is a 2-arc, one a 3-arc and the rest
# 4-arcs, so `n_4 - n_2` is invariant (6 for a hexagon root, 5 for a pentagon)
# and `B(d) = 3 B(d-1) + (n_4 - n_2)` — `3^(d+1) - 3` and `5 (3^d - 1) / 2`,
# which is what `sizehint!` uses below.
#
# Pentagons. All 12 pentagon base cells (4, 14, 24, 38, 49, 58, 63, 72, 83, 97,
# 107, 117) delete the same child, digit 1, the K axis — H3 has no per-cell
# deleted digit like Z7's two-valued `Z7_DELETED_DIGIT`, so the arc table needs
# no pentagon variant, only the skip. And the flag drops to `false` one level
# down for good: H3's pentagon descendants are exactly the all-digit-0 chain,
# and rim suffixes carry no zero digit, so only the subtree root can ever be a
# pentagon on a rim path.
#
# STATUS — this table was *fitted from observation, not derived.* Unlike
# IGEO7's, where the arc transitions fall out of the Eisenstein-lattice
# arithmetic in `engine.jl`, there is no proof here that these offsets are
# forced; they are the assignment that reproduces libh3's answers everywhere it
# has been checked. What has been checked, and what a reader may therefore rely
# on:
#   - exhaustive agreement with the definition (enumerate the subtree, keep the
#     cells with a `gridDisk` neighbor outside it) at depths 0-4 over 63 roots
#     spanning both parities, all 12 pentagon base cells, higher-resolution
#     pentagons, base-cell-crossing roots and random roots at res 1-8; and at
#     depths 5-6 over a five-root subset;
#   - agreement with the generic kernel fallback — an independent
#     implementation, which asks the grid at every level — at depths 7-9;
#   - sampled soundness and classification agreement at depth 10;
#   - an exhaustive search over all 720 digit-to-direction assignments, holding
#     the transitions fixed, which leaves exactly the six rotations of
#     `_H3_DIGIT_DIR` — the ring is pinned up to the rotation that the root
#     state `(6, 0)` makes immaterial.
# `test/H3/test_border.jl` re-runs the first three at reduced breadth, so the
# evidence survives as a regression net. Nothing above rules out a divergence at
# some depth or configuration not covered; treat a mismatch as a bug in this
# table, not in libh3.
# --------------------------------------------------------------------------

"""
    _H3_DIGIT_DIR[digit] -> Int

Position `0:5` on the direction ring of the IJK digit `digit in 1:6`, i.e. the
inverse of the geometric cyclic order `K(1), JK(3), J(2), IJ(6), I(4), IK(5)`.
"""
const _H3_DIGIT_DIR = (0, 2, 1, 4, 5, 3)

"""
    _h3_border_step(state, digit, level) -> NTuple{2,Int}

State a `digit` child at absolute resolution `level` inherits from a cell in
`state`, or `(0, 0)` when that child is interior to the subtree. A state is
`(L, s)`: the arc of exposed directions `s, s+1, ..., s+L-1` (mod 6) in
`_H3_DIGIT_DIR` positions. `L == 6` is the subtree root, the one state with no
arc ends.
"""
@inline function _h3_border_step(state::NTuple{2,Int}, digit::Int, level::Int)
    L, s = state
    (L == 0 || digit == 0) && return (0, 0)
    t = @inbounds _H3_DIGIT_DIR[digit]
    o = mod(t - s, 6)
    o < L || return (0, 0)
    if iseven(level)                         # child level is Class II
        L < 6 && o == 0 && return (2, mod(t + 1, 6))
        L < 6 && o == L - 1 && return (3, mod(t - 1, 6))
        return (4, mod(t - 1, 6))
    else                                     # child level is Class III (mirror)
        L < 6 && o == 0 && return (3, mod(t - 1, 6))
        L < 6 && o == L - 1 && return (2, mod(t - 2, 6))
        return (4, mod(t - 2, 6))
    end
end

"""
    _h3_fill_border!(out, z, res, target, state, pentagon) -> out

Digit-lexicographic DFS over the border automaton, appending every res-`target`
rim descendant of `z` to `out`. `z` already carries `target` in its resolution
field, so each step only overwrites one digit slot and the leaves come out
fully formed.

Digits run `1:6` because digit 0 is never on the rim (see above), and `pentagon`
is passed `false` to every recursive call rather than retested: the first digit
of a rim suffix is already nonzero, so no node below the root is a pentagon.
"""
function _h3_fill_border!(out::Vector{UInt64}, z::UInt64, res::Int, target::Int,
    state::NTuple{2,Int}, pentagon::Bool)
    shift = _h3_digit_shift(res + 1)
    cleared = z & ~(UInt64(7) << shift)
    for digit in 1:6
        pentagon && digit == 1 && continue    # the K axis, deleted under a pentagon
        child_state = _h3_border_step(state, digit, res + 1)
        child_state[1] == 0 && continue
        child = cleared | (UInt64(digit) << shift)
        if res + 1 == target
            push!(out, child)
        else
            _h3_fill_border!(out, child, res + 1, target, child_state, false)
        end
    end
    return out
end

"""
    _h3_border_descendants(cell, level, leaf_level) -> Vector{UInt64}

The rim of `cell`'s subtree at `leaf_level`, ascending. Validates `cell` and the
levels first; see `DGG.subtree_border` below for the error contract.
"""
function _h3_border_descendants(cell::UInt64, level::Int, leaf::Int)
    0 <= leaf <= H3Native.MAX_RESOLUTION ||
        throw(ArgumentError("H3 resolution must be in 0:$(H3Native.MAX_RESOLUTION)"))
    # libh3 validates neither `cellToChildren` nor `cellToChildrenSize`, so the
    # `cell_descendants` route the generic fallback takes for its degenerate
    # case buys H3 no id checking at all (`cellToChildren(garbage, res)` happily
    # returns a subtree of garbage). Hence the explicit `isValidCell`: without
    # it a malformed index — a cleared padding slot, a digit 7 in an active
    # slot, base cell 122+, a leading K digit under a pentagon — would come back
    # as a confidently enumerated rim of cells that do not exist.
    H3Native.is_valid_cell(cell) ||
        throw(ArgumentError("H3 cell 0x$(string(cell; base=16)) is not a valid cell"))
    # As in `descendant_range`: a mislabelled `level` is a caller bug, not a
    # cheaper spelling of "read it off the index". The generic fallback resolves
    # membership with `cell_parent(..., level)`, so it would answer a different
    # question here rather than the same one faster.
    _h3_resolution(cell) == level || throw(ArgumentError(
        "H3 cell 0x$(string(cell; base=16)) has resolution $(_h3_resolution(cell)), not level $level"))
    # A depth-0 subtree is the cell itself, whose whole neighborhood is outside
    # it. Spelled directly because the validation above has already done what
    # the fallback delegates to `cell_descendants` for.
    leaf == level && return UInt64[cell]
    # The resolution field moves to `leaf` once; the digit slots between `level`
    # and `leaf` are filled in on the way down and the ones below `leaf` keep
    # the 7 padding they already carry.
    z = (cell & ~_H3_RESOLUTION_MASK) | (UInt64(leaf) << 52)
    p3 = 3^(leaf - level)
    out = Vector{UInt64}()
    sizehint!(out, H3Native.is_pentagon(cell) ? (5 * (p3 - 1)) ÷ 2 : 3 * p3 - 3)
    return _h3_fill_border!(out, z, level, leaf, (6, 0), H3Native.is_pentagon(cell))
end

# The subtree rim without a single neighbor query, in O(result) rather than the
# fallback's O(subtree) — the block comment above says how, and says what the
# transition table rests on.
#
# Enumeration is digit-lexicographic depth-first, which *is* ascending id order
# here: the emitted ids agree on every bit except the digit slots between
# `level` and `leaf_level`, higher slots sit in higher bits, and each slot is
# visited in increasing digit order — so the result is sorted by construction
# and never sorted again. `test/H3/test_border.jl` asserts that rather than
# trusting it.
#
# `leaf_level < level` is the kernel's own argument mistake and is an
# `ArgumentError` as in `cell_descendants`. Past that, H3 has no error type of
# its own — `descendant_range` reports a bad resolution and a level/index
# mismatch as `ArgumentError` too — so an out-of-range `leaf_level`, an index
# that is not a valid H3 cell, and an index whose resolution is not `level` are
# all `ArgumentError`s naming the offending value.
function DGG.subtree_border(::DGG.H3DGGS, level::Integer, id, leaf_level::Integer)
    Int(level) <= Int(leaf_level) || throw(ArgumentError("expected level <= leaf_level"))
    return _h3_border_descendants(UInt64(id), Int(level), Int(leaf_level))
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

# What the generic lookup operations (`neighbor_indices`, `stencil`, `zonal`)
# ask of a lookup: which system, which level.
DGG.dggs_system(::H3Lookups.H3Lookup) = DGG.H3DGGS()
DGG.dggs_level(l::H3Lookups.H3Lookup) = l.resolution

# Treeifying a lookup directly is the shortest path from a `DimensionalData`
# dimension to a `Regridder`, and it needs nothing from this file: the method
# is generic over `AbstractDGGSLookup` (`core/lookups.jl`), routing a stored id
# vector through the constructor above and a `DGGSGlobeIds` through the dense
# `DGGSGrid` instead. All that is per-system is the manifold, which is what
# makes the one-argument `treeify(l)` resolve at all.
GOCore.best_manifold(::H3Lookups.H3Lookup) = GO.Spherical()
