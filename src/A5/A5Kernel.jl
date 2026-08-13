# ---------------------------------------------------------------------------
# A5 operations-kernel wiring
#
# A5 is the reference *non-uniform-radix* system of `core/kernel.jl`. Its
# hierarchy has three regimes, so nothing here may be derived from a `radix`
# (`src/core/systems/a5.jl` deliberately wires none):
#
#   level -1  the world cell            — outside the kernel's level space
#   level  0  12 dodecahedron pentagons — `res0_cells()`
#   level  1  5 quintants per pentagon  — the pentagon cut into triangles
#   level >1  4 Hilbert children each   — aperture 4 from here down
#
# so a level-`r >= 1` grid has `60 * 4^(r - 1)` cells, and `subtree_leaf_count`
# is `5 * 4^(leaf - 1)` from a root and `4^(leaf - level)` below it. Every
# operation delegates to `A5Native` — the per-system A5 tree module this kernel
# replaced is gone.
#
# Trait-wise A5 is the most conservative system in the package: ids are
# structural (`has_ordinal_ids` false), no `descendant_range` is wired
# (`has_descendant_ranges` false — the ids of a subtree *are* contiguous, but
# pinning the two-sided contract across the res-0/res-1 regime changes is out
# of scope), and parents do not contain their children geographically
# (`has_exact_subtree_cap` false). All three are the kernel defaults, so no
# method is written for them; the generic partial cursor therefore takes its
# `cell_parent` membership-filter path here — A5 is that path's first real
# system.
# ---------------------------------------------------------------------------

import ..DiscreteGlobalGrids as DGG
import GeometryOps as GO
import GeometryOpsCore as GOCore
using SmallCollections: SmallVector

# --------------------------------------------------------------------------
# Id model
# --------------------------------------------------------------------------

DGG.cell_id_type(::DGG.A5DGGS) = UInt64

# --------------------------------------------------------------------------
# Hierarchy
#
# `level` is redundant for a structural id — `A5Native` reads the resolution
# out of the id itself — and is only used to name the target resolution.
# --------------------------------------------------------------------------

# The native enumeration is ascending *below* res 0: a subtree there keeps its
# quintant bits and varies only the Hilbert digits, which sit above the
# resolution marker bit, so digit order is id order. The res-0 fan-out is not:
# `cell_to_children` walks `segment = 0:4` while `serialize` orders by
# `quintant = 5 * origin + (segment - first_quintant) mod 5`, a rotation. One
# `issorted` pass covers that case — and any caller that mislabels a cell's
# level — for O(n) rather than the O(n log n) of sorting unconditionally.
_ascending(ids::Vector{UInt64}) = issorted(ids; lt=(<=)) ? ids : sort!(ids)

DGG.root_ids(::DGG.A5DGGS) = collect(UInt64, A5Native.res0_cells())

DGG.cell_children(::DGG.A5DGGS, level::Integer, id) =
    _ascending(collect(UInt64, A5Native.cell_to_children(UInt64(id), Int(level) + 1)))

DGG.cell_parent(::DGG.A5DGGS, level::Integer, id, parent_level::Integer) =
    A5Native.cell_to_parent(UInt64(id), parent_level)

# `cell_to_children` enumerates any resolution gap directly, so the kernel's
# level-by-level expansion is never used for A5.
function DGG.cell_descendants(::DGG.A5DGGS, level::Integer, id, leaf_level::Integer)
    Int(level) <= Int(leaf_level) || throw(ArgumentError("expected level <= leaf_level"))
    return _ascending(collect(UInt64, A5Native.cell_to_children(UInt64(id), leaf_level)))
end

# `12` at level 0, `60 * 4^(level - 1)` below — and an `ArgumentError` past
# `MAX_GRID_RESOLUTION`, since res 30 exists only for the 42 quintants that fit
# the 64-bit encoding and so has no uniform full-world grid.
DGG.num_cells(::DGG.A5DGGS, level::Integer) = Int64(A5Native.num_cells(level))

# O(1) from the three regimes above; the generic fallback would enumerate the
# whole subtree. This is what `ncells` and the parallelize policy call per node.
function DGG.subtree_leaf_count(::DGG.A5DGGS, level::Integer, id, leaf_level::Integer)
    lvl = Int(level)
    leaf = Int(leaf_level)
    lvl <= leaf || throw(ArgumentError("expected level <= leaf_level"))
    0 <= leaf <= A5Native.MAX_GRID_RESOLUTION || throw(ArgumentError(
        "A5 full-grid resolution must be in 0:$(A5Native.MAX_GRID_RESOLUTION)"))
    lvl == leaf && return Int64(1)
    lvl == 0 && return Int64(5) * Int64(4)^(leaf - 1)
    return Int64(4)^(leaf - lvl)
end

# --------------------------------------------------------------------------
# Neighbors
#
# `A5Native._get_global_cell_neighbors` is A5's own adjacency walk and already
# returns `sort!(collect(::Set))`, so this wiring is a container change, a mode
# choice and two guards — not a computation.
#
# WHY `edge_only = true`. The native keyword defaults to `false`, which is a
# *vertex* neighborhood (everything sharing a corner), and that default is what
# A5 uses internally: `_spherical_to_cell` widens its point-location search
# with it, where extra candidates only cost time. The kernel's contract is edge
# adjacency — it is what `stencil` and the lookup halo table mean — and the two
# sets genuinely differ. Measured over complete levels: res 1 gives 3 edge
# against 11 vertex neighbors; res 2-4 give a flat 5 edge against 6, 7 or 8
# vertex neighbors (60 cells of each level report 6 and 60 report 8, the rest
# 7). Res 0 is the one level where the keyword changes nothing — a dodecahedron
# face shares a vertex with no face it does not also share an edge with, and
# `_get_res0_neighbors` ignores the keyword outright. Everywhere else passing
# the default would type-check and quietly answer a different question, so the
# keyword is spelled out.
#
# FIRST EXERCISE. Nothing in the package called `edge_only = true` before this
# wiring — `_spherical_to_cell` was the function's only caller and it takes the
# default — so the res-1 branch and the `edge_only` filters in
# `_find_quintant_neighbor_s` / `_get_boundary_neighbors` were, until now,
# dead code. The neighbor testsets in `test/A5/test_a5_kernel.jl` are the whole
# of their coverage, and they assert the kernel contract directly rather than
# trusting the native code: ascending, unique, self-excluded, same-resolution,
# and above all *symmetric* — the res >= 2 answer is a `Set` union of two
# independently computed families (the within-quintant lattice step and the
# quintant/face-seam special cases), which is exactly the shape that produces a
# one-directional edge. It does not: symmetry holds over every cell of res 0-8.
#
# THE BOUND IS MEASURED, not inferred from the cell shape. Sweeping complete
# levels 0-8 (1_310_712 cells) and 400-cell ordinal samples with their two-hop
# neighborhoods at res 9, 12, 15, 20, 25 and 29, the edge degree is not merely
# bounded by 5 but *constant* per regime: exactly 5 at res 0 (a dodecahedron
# face has 5 face-adjacent faces), exactly 3 at res 1 (a quintant is a triangle
# and its three edges reach the two sibling quintants and one across a face
# seam), exactly 5 at res >= 2 (the pentagon lattice is edge-to-edge, at the
# apex and along the seams as much as in a quintant's interior). No cell
# anywhere in that sweep reported 4, 6 or more edge neighbors, so 5 is the
# bound with no slack rather than a padded guess — the same discipline as the
# `cell_cap_inflation` measurement above, and the testset is the evidence.
#
# RES 30 IS REFUSED. `MAX_RESOLUTION` is 30 but `MAX_GRID_RESOLUTION` is 29:
# only 42 of the 60 quintants fit the res-30 encoding, so there is no complete
# res-30 grid (which is why `num_cells` and `subtree_leaf_count` already throw
# there). Adjacency is worse than incomplete at res 30 — it is wrong. A res-30
# neighbor landing in an unsupported quintant hits `serialize`'s
# `S >> 2, MAX_RESOLUTION - 1` fallback and comes back as a *res-29* id, so the
# native answer mixes resolutions, drops to 4 entries and is not symmetric.
# Rather than hand that to the kernel, this method range-checks the level the
# way `cell_to_ordinal` does and refuses.
# --------------------------------------------------------------------------

DGG.max_neighbors(::DGG.A5DGGS) = 5

function DGG.cell_neighbors(::DGG.A5DGGS, level::Integer, id)
    lvl = Int(level)
    0 <= lvl <= A5Native.MAX_GRID_RESOLUTION || throw(ArgumentError(
        "A5 full-grid resolution must be in 0:$(A5Native.MAX_GRID_RESOLUTION)"))
    cell = UInt64(id)
    # The id guard `cell_to_ordinal` uses, for the same reason: a structural id
    # carries its own resolution, and the native walk reads it out of the id
    # rather than out of `level`, so a mislabelled level would otherwise return
    # a different level's neighbors under this level's name. `get_resolution`
    # also rejects the world cell (resolution -1), which `deserialize` would
    # hand on to an `InexactError` deep in the lattice arithmetic. An id whose
    # quintant bits name no origin still reaches the native `BoundsError`, as
    # it does through `cell_children` and `cell_to_ordinal`.
    A5Native.get_resolution(cell) == lvl || throw(ArgumentError(
        "cell $(repr(cell)) is at A5 resolution $(A5Native.get_resolution(cell)), not $lvl"))
    out = SmallVector{5,UInt64}()
    # The native return is already ascending and duplicate-free (it is a `Set`
    # collected through `sort!`), so this is a straight refill; `_insert_sorted`
    # keeps the kernel's ordering contract stated here rather than borrowed
    # from A5Native's internals, at the cost of a binary search over <= 5 ids.
    for neighbor in A5Native._get_global_cell_neighbors(cell; edge_only=true)
        out = DGG._insert_sorted(out, neighbor)
    end
    return out
end

# `subtree_border` keeps the kernel's generic fallback. A5 has no index-only
# shortcut of the kind `IGEO7` reads off its Z7 digits: rim membership there is
# a digit predicate because the tiling is congruent and rep-7, while an A5
# subtree is four Hilbert children that cover the parent's area but not its
# footprint (`has_exact_subtree_cap` is false for the same reason). The
# fallback's level-by-level expansion is therefore the end state here, and that
# is fine: it pays a `cell_neighbors` sweep and a `cell_parent` ascent per rim
# child per level where IGEO7 reads the answer off its digits, but the frontier
# it carries is the rim, not the subtree.
#
# One caveat, because the fallback prunes rather than enumerating: it drops the
# children of any cell that touches nothing outside the subtree, on the premise
# that a cell's children touch nothing outside the children of that cell and of
# its edge neighbors. Measured over complete levels, that closure holds at res 0
# and from res 2 down (0 violations over all 12 + 240 cells and a res-3 sample)
# but *fails at res 1*: 120 violations over the 60 quintants, two per quintant.
# They are the res-2 cells at a pentagon apex, whose edge neighbors reach into
# pentagons that are only *vertex* neighbors of the parent quintant — a
# consequence of the res-1 regime being triangles with 3 edges where the res-2
# lattice below it is pentagons with 5.
#
# The gap is unreachable, which is why the fallback is still exact here: a
# quintant is never interior to any subtree the pruning sees. Its cross-face
# neighbor always lies in a different pentagon, so for a res-0 root all five
# quintants stay in the frontier and nothing is pruned at that step; for a res-1
# root the root *is* the frontier and the first expansion tests res-2 cells
# directly; a deeper root never traverses res 1 at all. The testset pins both
# halves — the closure counts per level, the "no quintant is ever interior"
# fact, and then the whole operation against a brute-force pass over
# `cell_descendants` at root levels 0-3 and depths 0-6.

# --------------------------------------------------------------------------
# Dense ordinals
#
# Ascending id order at a fixed resolution is `(quintant, S)` lexicographic:
# `serialize` lays an id out as `quintant | S | marker`, with `quintant` in the
# top bits, the `2 * (level - 1)` Hilbert digits below it and the resolution
# marker below those. `quintant = 5 * origin.id + (segment - first_quintant)
# mod 5` runs over `0:59` and `S` over `0:4^(level - 1) - 1`, so
#
#     ordinal = quintant * 4^(level - 1) + S + 1
#
# is both a bijection onto `1:num_cells` and strictly monotone in the id, which
# is exactly the kernel's contract. Level 0 is the 12-entry root table.
# --------------------------------------------------------------------------

_a5_quintant(cell::A5Native.A5Cell) =
    5 * cell.origin.id + mod(cell.segment - cell.origin.first_quintant, 5)

function DGG.cell_to_ordinal(::DGG.A5DGGS, level::Integer, id)
    cell = UInt64(id)
    lvl = Int(level)
    0 <= lvl <= A5Native.MAX_GRID_RESOLUTION || throw(ArgumentError(
        "A5 full-grid resolution must be in 0:$(A5Native.MAX_GRID_RESOLUTION)"))
    deserialized = A5Native.deserialize(cell)
    # A structural id carries its own resolution, so a mislabelled `level` would
    # otherwise silently return an ordinal of a different level's grid — the
    # `get_resolution(id) == level` discipline the range systems already have.
    deserialized.resolution == lvl || throw(ArgumentError(
        "cell $(repr(cell)) is at A5 resolution $(deserialized.resolution), not $lvl"))
    lvl == 0 && return searchsortedfirst(A5Native.res0_cells(), cell)
    return Int(_a5_quintant(deserialized) * Int64(4)^(lvl - 1) + Int64(deserialized.S) + 1)
end

# An ordinal outside `1:num_cells` is a `DGG.OrdinalRangeError`, the one type
# every wiring of this operation raises (`src/core/kernel.jl`): the caller is
# told which system's level it overran and what that level's ordinals are,
# which the `BoundsError` this used to throw could not say.
function DGG.ordinal_to_cell(::DGG.A5DGGS, level::Integer, ordinal::Integer)
    lvl = Int(level)
    total = DGG.num_cells(DGG.A5DGGS(), lvl)          # also range-checks `lvl`
    1 <= ordinal <= total || throw(DGG.OrdinalRangeError(
        DGG.system_name(DGG.A5DGGS()), lvl, Int(ordinal), total))
    lvl == 0 && return A5Native.res0_cells()[Int(ordinal)]
    span = Int64(4)^(lvl - 1)
    quintant, S = divrem(Int64(ordinal) - 1, span)
    origin = A5Native.ORIGINS[Int(quintant ÷ 5) + 1]
    segment = mod(Int(quintant % 5) + origin.first_quintant, 5)
    return A5Native.serialize(A5Native.A5Cell(origin, segment, UInt64(S), lvl))
end

# --------------------------------------------------------------------------
# Pruning hook
#
# Not wired: `has_descendant_ranges` stays at its `false` default, so the
# generic trees fall back to `cell_parent` membership filtering and no
# `descendant_range` method is needed.
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# Geometry
#
# A5 cell edges are straight in the dodecahedral face plane, not on the sphere,
# so the native ring subdivides them (`segments = :auto`: `2^(6 - resolution)`
# pieces per edge down to res 6, one above) before inverse-projecting. That
# subdivided ring is the cell's actual shape, and it is what every downstream
# consumer sees — the regridding polygon and the caps alike. Level 0 rings have
# 5 vertices, level 1 rings 3 (a quintant is a triangular slice of a pentagon)
# and level >= 2 rings 5 again, times the subdivision.
# --------------------------------------------------------------------------

# FRAME: the ring must come from the native *lonlat* boundary, not from
# `cell_boundary_cartesian`. The latter stops at `A5Native._to_cartesian`, which
# is A5's internal spherical frame: it applies neither the `LONGITUDE_OFFSET`
# un-rotation nor the authalic -> geodetic latitude conversion that
# `A5Native._to_lonlat` does. Mixing it with the `cell_to_lonlat` center below
# put boundary and center 93 degrees apart in longitude and blew leaf caps up to
# ~2.3 rad. Both ends of this file now land in the same geographic frame as
# `A5Native.lonlat_to_cell`, which is what the trees, the caps and every
# cross-system regrid actually address.
#
# `segments = :auto` is the native default and stays: it is what makes the ring
# the cell's real spherical shape rather than a chord polygon.
#
# One consequence worth naming: A5 is equal-area on the ELLIPSOID, so in the
# geodetic frame its unit-sphere areas carry that latitude conversion — ~1% peak
# to peak across a level. Undoing the conversion on the wired ring restores
# equal area to round-off, which is what the geometry testset asserts. Alignment
# with `lonlat_to_cell` and with every other system is worth that 1%: the
# alternative frame is exactly-equal-area but 0.19 degrees off the coordinates
# A5 itself reports.
function _a5_unit_point(lon, lat)
    λ = deg2rad(lon)
    φ = deg2rad(lat)
    cosφ = cos(φ)
    return GO.UnitSphericalPoint(cosφ * cos(λ), cosφ * sin(λ), sin(φ))
end

# The native `closed_ring = true` is NOT the kernel's `closed = true`: it seals
# the ring *before* reversing it, so the result is the open ring rotated by one
# vertex rather than the open ring with its first point repeated. The kernel
# contract is the latter, so the closing point is appended here.
function DGG.cell_boundary(::DGG.A5DGGS, level::Integer, id; closed::Bool=false)
    ring = A5Native.cell_boundary(UInt64(id); closed_ring=false, segments=:auto)
    points = [_a5_unit_point(p[1], p[2]) for p in ring]
    closed && !isempty(points) && push!(points, points[1])
    return points
end

# Native cell center (`cell_to_lonlat`, the face-plane centroid pulled back
# through the equal-area projection), not the kernel's boundary mean: it is the
# center the A5 lattice is built around, and one call instead of a ring.
function DGG.cell_center(::DGG.A5DGGS, level::Integer, id)
    lon, lat = A5Native.cell_to_lonlat(UInt64(id))
    return _a5_unit_point(lon, lat)
end

# A5 overhangs its ancestors far further than the aperture-7 systems do, so the
# package default is nowhere near enough. The cause is structural: an A5 cell is
# a fixed pentagon placed at a lattice point and scaled by `2^-resolution`
# (`A5Native._get_pentagon_vertices`), and that pentagon is not a rep-4 tile —
# the four cells a Hilbert digit names as "children" cover the parent's area but
# not its footprint. Measured on random points, ~37% of them land in a res-`r`
# cell whose `cell_to_parent` is not the res-`r-1` cell containing them (0% for
# res 0 -> 1, where the quintant cut *is* exact). `has_exact_subtree_cap` being
# false is that fact; the inflation below is its size.
#
# CAP-VALIDATION (`test/A5/test_a5_kernel.jl`) measures the union ratio — max
# descendant-vertex distance from the cell's cap center over the cell's own max
# vertex distance — exhaustively over res 0-2 and on samples at res 3, 5 and 8,
# deltas 1-6, plus a delta 1-8 convergence probe. Worst measured 1.45159 (res 8,
# delta 6); the increments halve cleanly from delta 4 on, so the geometric tail
# puts the supremum at 1.46872 (H3 measures 1.052, IGEO7 1.048). The ratio is
# essentially level-independent from res 3 down — the drift is a fixed fraction
# of the current cell size at every level, and the series converges. A5
# therefore pins its inflation at 1.75: 16% headroom over the extrapolated
# supremum, at the price of caps ~2.4x the default's area.
#
# The overhang is real geometry, not a bad center: the boundary-mean center the
# kernel would fall back to measures 1.3283 at res 2 / 1.3937 at res 3 against
# the native center's 1.3547 / 1.4082 over the same deltas — a ~2% difference
# that changes no budget, for a full ring per center instead of one call.
DGG.cell_cap_inflation(::DGG.A5DGGS) = 1.75

# `cells_cap` / `subtree_cap` keep the kernel definitions.

# --------------------------------------------------------------------------
# Lookup convenience
# --------------------------------------------------------------------------

"""
    DGGSPartialGrid(l::A5Lookup; kwargs...)

The lookup's stored cells as a generic partial grid. Cell-id validity at
`l.resolution` is trusted from the lookup (which checks it when constructed
with `validate=true`); the generic constructor only re-checks ordering and
element type. `kwargs` reach `DGGSPartialGrid`'s `bucket_size` / `root_level` /
`root_id`.
"""
DGG.DGGSPartialGrid(l::A5Lookups.A5Lookup; kwargs...) =
    DGG.DGGSPartialGrid(DGG.A5DGGS(), l.resolution, l.data; kwargs...)

# What the generic lookup operations (`neighbor_indices`, `stencil`, `zonal`)
# ask of a lookup: which system, which level. With `cell_neighbors` wired above
# those three now answer for A5 through their generic definitions; this pair is
# what points them at the right system and resolution.
DGG.dggs_system(::A5Lookups.A5Lookup) = DGG.A5DGGS()
DGG.dggs_level(l::A5Lookups.A5Lookup) = l.resolution

# Treeifying a lookup directly is the shortest path from a `DimensionalData`
# dimension to a `Regridder`, and it needs nothing from this file: the method
# is generic over `AbstractDGGSLookup` (`core/lookups.jl`), routing a stored id
# vector through the constructor above and a `DGGSGlobeIds` through the dense
# `DGGSGrid` instead. All that is per-system is the manifold, which is what
# makes the one-argument `treeify(l)` resolve at all.
GOCore.best_manifold(::A5Lookups.A5Lookup) = GO.Spherical()
