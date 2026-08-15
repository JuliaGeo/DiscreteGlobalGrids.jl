# ---------------------------------------------------------------------------
# Hierarchical-grid fast-path contracts. Arguments use `(sys, c)` with any
# target level last; cell ids encode their own level. Overrides must preserve
# base-interface semantics.
# ---------------------------------------------------------------------------

# ===========================================================================
# Required
# ===========================================================================

"""
    cellindextype(sys::AbstractHierarchicalGridSystem) -> Type{<:AbstractCellIndex}

The canonical cell index type of `sys`: the scheme its grids return from
[`cellindex`](@ref), its ids sort in, and its neighbour containers have as
`eltype`.

**Required.** One canonical scheme per system; alternates are declared with
[`cellindextypes`](@ref) and converted with [`reindex`](@ref).
"""
function cellindextype end

"""
    levels(sys::AbstractHierarchicalGridSystem) -> AbstractUnitRange{Int}

The valid refinement levels of `sys`, coarsest first: `first(levels(sys))` is
the level of [`rootcells`](@ref) and `last(levels(sys))` is
[`max_level`](@ref).

**Required.**

A system without an intrinsic depth limit must still bound the range at the
deepest level supported by its id encoding and `Int` cell counts.
"""
function levels end

"""
    levelgrid(sys::AbstractHierarchicalGridSystem, l::Integer) -> AbstractGrid

The **complete** grid of `sys` at level `l`: every cell the system has at that
level, in the system's canonical dense order.

**Derived, with a default**, and the entry point every consumer uses. The
default is `HierarchicalLevelGrid(sys, l)`, checked against
[`levels`](@ref) — so a system implements the five level-grid primitives below
and never writes a grid type. Overriding this is the escape hatch for a system
whose grid genuinely carries state beyond `(sys, l)`; none of the six shipped
here does.

`system(levelgrid(sys, l)) === sys` and `level(levelgrid(sys, l)) == l`. The
grid is a lightweight descriptor, not a materialised cell list — constructing
one must not be O(cells).

`l` outside [`levels(sys)`](@ref levels) throws an `ArgumentError`. For a
subset of a level, see `PartialGrid`.
"""
function levelgrid end

"""
    ncells(sys::AbstractHierarchicalGridSystem, l::Integer) -> Int
    cellindex(sys::AbstractHierarchicalGridSystem, l::Integer, i::Int) -> AbstractCellIndex
    cellposition(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex) -> Union{Int,Nothing}
    cell_boundary(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex)
    cell_centroid(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex)

The complete level of `sys`, described cell by cell: its size, its dense order
in both directions, and the geometry of one cell.

**Required**, unless the system overrides [`levelgrid`](@ref) with a grid type
of its own and answers the [`AbstractGrid`](@ref) contract there.

These are the five methods [`HierarchicalLevelGrid`](@ref) forwards to. Each
carries its grid-level contract — [`ncells`](@ref), [`cellindex`](@ref),
[`cellposition`](@ref), [`cell_boundary`](@ref), [`cell_centroid`](@ref) —
minus the two things the grid method has already settled:

  - `cellindex` may assume `i in 1:ncells(sys, l)`. The grid bounds-checks.
  - `cellposition` may assume `c` is in [`cellindextype(sys)`](@ref
    cellindextype) and at the level being asked about: the grid answers
    `nothing` for a cell at another level, and reindexes the alternate schemes
    first. It still returns `nothing` for a canonical id that names no cell at
    all — a malformed encoding, or an index past the end of its level.
  - the geometry pair may assume only the **level**, which the grid enforces by
    throwing an `ArgumentError` rather than by answering. It does *not* reindex:
    an id in an alternate scheme at the right level is forwarded as it stands,
    so a system that wants to accept one must say so with its own method. Left
    alone, the fallthrough is a `MethodError`, exactly as it is on the system
    method itself.

The geometry pair takes **no level argument**: an [`AbstractCellIndex`](@ref)
is self-describing, so a level beside it would only be a second source of truth
that could disagree with the first.
"""
ncells(::AbstractHierarchicalGridSystem, ::Integer)

"""
    rootcells(sys::AbstractHierarchicalGridSystem)

The top-level cells of `sys` — the base tessellation everything else refines —
in ascending canonical order, all at level `first(levels(sys))`.

**Required.**

These are the roots the generic tree descent starts from, so this must be a
small, cheap collection (12 for HEALPix, 122 for H3, 12 for IGeo7).
"""
function rootcells end

"""
    Base.parent(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex) -> AbstractCellIndex

The parent of `c`: the cell at `level(c) - 1` whose [`children`](@ref) contain
`c`.

**Required**, and required to be **analytic** — index arithmetic on the id,
with no lookup table, no allocation and no geometry. Tree descent and ancestry
walks call this in inner loops.

Calling it on a root cell (`level(c) == first(levels(sys))`) throws an
`ArgumentError`; there is no `nothing` sentinel here, because a caller that
descends the hierarchy always knows the level it is at.

`parent(sys, c)` and `children(sys, c)` are inverses:
`c in children(sys, parent(sys, c))`.
"""
Base.parent(::AbstractHierarchicalGridSystem, ::AbstractCellIndex)

"""
    children(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex)

The immediate children of `c` — its refinement at `level(c) + 1` — in ascending
canonical order.

**Required**, and required to be **analytic**, for the same reason as
[`parent`](@ref).

Child count may vary by cell; generic code must not assume it equals the nominal
aperture.

Calling it on a cell at `max_level(sys)` throws an `ArgumentError`.
"""
function children end

"""
    node_extent(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex) -> GO.UnitSpherical.SphericalCap

The covering region of the entire subtree rooted at `c`.

**Required.**

> `node_extent(sys, c)` contains the geometry of **every descendant of `c`, at
> every depth** — every point of every cell boundary in the subtree, all the
> way down to `max_level(sys)`.

Tree pruning depends on this covering law. Over-coverage only reduces pruning;
under-coverage can omit valid results. For a convex cap of angular radius at
most 90°, containing all boundary vertices also contains their great-circle
arcs. Non-convex extents must establish containment of the full geometry.

A cell's own boundary need not cover its descendants; aperture-7 children, for
example, can extend beyond the parent boundary.

The generic implementation inflates the cell's bounding cap by
[`cap_inflation(sys)`](@ref cap_inflation). Systems may provide a tighter
covering cap.

The result must cover through `max_level(sys)`, independent of a caller's
planned traversal depth.

Extents are `SphericalCap`s throughout this package, at every node of every
tree, which is what lets one predicate vocabulary serve all of them.
"""
function node_extent end

# ===========================================================================
# Traits
# ===========================================================================

"""
    has_sorted_subtrees(sys::AbstractHierarchicalGridSystem) -> Bool

Whether the descendants of any cell, at any fixed deeper level, occupy a
**contiguous interval** of that level's canonical dense order.

`false` by default. A system opting in must implement
[`descendant_range`](@ref). The property enables range-based subtree membership
and traversal.

Systems whose canonical order is a space-filling curve (nested HEALPix, Z7,
H3's resolution-major order) have this property; it is a fact about the
ordering, and declaring it falsely produces silently wrong subtree answers.
"""
has_sorted_subtrees(::AbstractHierarchicalGridSystem) = false

"""
    cap_inflation(sys::AbstractHierarchicalGridSystem) -> Float64

The factor by which the default [`node_extent`](@ref) inflates a cell's own
bounding-cap radius so that the cap covers the cell's entire subtree.

Defaults to `1.2`.

The value bounds how far descendants extend beyond a cell's bounding cap and
must be validated for the system's refinement geometry. Systems overriding
`node_extent` ignore it.

Raising it costs query time (looser pruning). Setting it too low is a
correctness bug — see the covering law in [`node_extent`](@ref).
"""
cap_inflation(::AbstractHierarchicalGridSystem) = 1.2

"""
    max_neighbors(sys::AbstractHierarchicalGridSystem, connectivity::Connectivity = Vertex()) -> Int

A **static** upper bound on the number of `connectivity`-neighbours of any cell
of `sys`, at any level.

The static bound permits fixed-capacity neighbour containers. Individual cells
may have fewer neighbours.

There is no default: a system that has not thought about the bound gets a
`MethodError` rather than a silently wrong capacity.
"""
function max_neighbors end

max_neighbors(sys::AbstractHierarchicalGridSystem) = max_neighbors(sys, Vertex())

"""
    max_level(sys::AbstractHierarchicalGridSystem) -> Int

The deepest valid level of `sys`: `last(levels(sys))`, which is the default
implementation.
"""
max_level(sys::AbstractHierarchicalGridSystem) = last(levels(sys))

cellindextypes(sys::AbstractHierarchicalGridSystem) = (cellindextype(sys),)

# ===========================================================================
# Derived hierarchy operations
#
# Declared here with their contracts; implemented generically in
# `src/fallbacks/`, and overridden by systems that can answer them in closed
# form (which, for `descendant_range`, is the only way they are answered).
# ===========================================================================

"""
    ancestor(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex, l::Integer) -> AbstractCellIndex

The ancestor of `c` at level `l`, for `first(levels(sys)) <= l <= level(c)`.

`ancestor(sys, c, level(c))` is `c` itself, and `ancestor(sys, c, level(c) - 1)`
is [`parent(sys, c)`](@ref parent). `l > level(c)` throws an `ArgumentError` —
an ancestor is never deeper than the cell, and a system with a closed-form
answer (drop `level(c) - l` digits) should override rather than iterate.
"""
function ancestor end

"""
    descendants(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex, l::Integer)

Every descendant of `c` at level `l`, in ascending canonical order, for
`level(c) <= l <= max_level(sys)`.

`descendants(sys, c, level(c))` is `[c]`. `l < level(c)` throws an
`ArgumentError` (uniformly across systems, so generic code can catch it).

This materializes `O(subtree)` ids. Use [`descendant_range`](@ref) when
available, or [`subtree_border`](@ref) when only the rim is needed.
"""
function descendants end

"""
    subtree_border(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex, l::Integer; connectivity::Connectivity = Vertex())

The **rim** of `c`'s subtree at level `l`: every level-`l` descendant of `c`
that has a neighbour which is *not* a descendant of `c`.

`subtree_border(sys, c, level(c))` is `[c]` — a depth-0 subtree is the cell
itself, and its entire neighbourhood lies outside it. `l < level(c)` throws an
`ArgumentError`, as it does for [`descendants`](@ref).

`connectivity` selects which adjacency defines "has a neighbour outside", with
the same meaning as in [`neighbors`](@ref). It changes nothing where the two
adjacencies coincide (H3 and IGeo7, whose vertices are all 3-valent) or where
the rim comes out the same set either way (HEALPix); A5's 4-valent corners make
them genuinely different relations, so assume nothing.

The generic fallback walks [`descendant_range`](@ref) and tests each cell's
[`neighbors`](@ref). Systems may override with an `O(rim)` algorithm. Order is
ascending canonical order unless documented otherwise.

This is `collect` of `EdgeCellIterator(sys, c, l; connectivity)`, which is the
same walk resumable and in `O(depth)` memory — reach for the iterator when the
rim is large or a prefix of it will do.

See also [`subtree_interior`](@ref), the complement.
"""
function subtree_border end

"""
    subtree_interior(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex, l::Integer; connectivity::Connectivity = Vertex())

The level-`l` descendants of `c` that are **not** on the rim — the complement of
[`subtree_border`](@ref) within [`descendants`](@ref), in ascending canonical
order.

    subtree_border(sys, c, l) ∪ subtree_interior(sys, c, l) == descendants(sys, c, l)

with the two disjoint. `subtree_interior(sys, c, level(c))` is empty: the cell
itself is its own rim.

This materializes most of the subtree. It is `collect` of
`InnerCellIterator(sys, c, l; connectivity)`, which generates the interior from
the rim walk's pruned branches — never a border set, never the rim materialized
— in `O(depth)` memory. For large subtrees, use the iterator.
"""
function subtree_interior end

"""
    rim_engine(sys, c, target::Int, connectivity)
    interior_engine(sys, c, target::Int, connectivity)
    halo_engine(sys, c, target::Int, connectivity)

The iteration engine `EdgeCellIterator` / `InnerCellIterator` /
`SubtreeHaloIterator` forwards the whole iteration protocol to — the single
place a system overrides to ship an `O(rim)` subtree walk or an `O(halo)` halo
walk, and the single place both the lazy and the eager
([`subtree_border`](@ref) / [`subtree_interior`](@ref) / `subtree_halo`) faces
of it read.

An engine is any iterator over `cellindextype(sys)`. All three methods own their
level validation, so their `ArgumentError`s are the ones the eager verbs raise.

Engine selection is **private multiple dispatch on the system type**: a system
ships a fast path by adding one method here, not by setting a trait anyone can
read and not by anything inspecting the method table at runtime. That keeps each
engine's protocol monomorphic, and it keeps "which walk did I get?" out of the
public surface, where it would become a compatibility promise.

The generic implementations walk [`descendant_range`](@ref) with one
[`ancestor`](@ref) test per cell, and materialize where
[`has_sorted_subtrees`](@ref) is `false`. `halo_engine`'s generic implementation
walks the hierarchy from OUTSIDE the subtree instead — see `SubtreeHaloIterator`
— because the halo is not a sub-interval of any one subtree's position range.
"""
function rim_engine end

@doc (@doc rim_engine)
function interior_engine end

@doc (@doc rim_engine)
function halo_engine end

"""
    lattice_decode(sys, c) -> (ix, iy, face)
    lattice_cell(sys, level::Int, ix, iy, face) -> cell
    face_orientation(sys, face) -> UInt8

The three lines an aperture-4 system writes to let the shared square halo walk
cross its seams. Only the systems whose cells ARE an aligned square lattice per
face — HEALPix, S2, ISEA4R — implement them; there is no generic fallback,
because there is no generic lattice.

`lattice_decode` and `lattice_cell` are the system's existing face-lattice codec
under one name (`nested_to_xyf`/`xyf_to_nested`, `hilbert_to_xyf`/
`xyf_to_hilbert`, `morton_to_xyd`/`xyd_to_morton`), with `face` 0-based and
`(ix, iy)` in `0:2^level - 1`. `face_orientation` is the curve state a face's
ROOT is read under, before any position bits are consumed — `0x0` for the two
Morton systems, and S2's odd-face swap for the Hilbert one.

The walk needs no seam table of its own: it asks the system for the neighbours
of a few rim cells and reads the answers back through `lattice_decode`. That is
the whole per-system surface, which is why a fourth square system would need
nothing else.
"""
function lattice_decode end

@doc (@doc lattice_decode)
function lattice_cell end

@doc (@doc lattice_decode)
function face_orientation end

"""
    hex_child_direction(sys, c) -> Int
    seeded_rim_engine(sys, c, target::Int, arclen::Int, start::Int)

The two lines an aperture-7 system writes to let the shared calibrated halo walk
(`hex_halo_engine`) approach a subtree from its neighbours. Only the systems
that already own a subtree-rim automaton over an arc of exposed lattice
directions — H3 and IGeo7 — implement them; there is no generic fallback,
because there is no generic automaton.

`hex_child_direction` is the position `0:5` on the six-direction ring of the step
from `c`'s parent to `c`, and `-1` for the centre child (which has no direction)
or for a root cell (which has no parent). It is the system's existing digit →
direction table (`_H3_DIGIT_DIR`, `SIGMA_J`) read through the cell's own last
digit — no child list is searched.

`seeded_rim_engine` is the system's rim automaton entered at an ARBITRARY arc
`(arclen, start)` rather than at the fully exposed `(6, 0)` a subtree root gets:
`c`'s level-`target` descendants reachable along the arc of exposed directions
`start, start+1, …, start+arclen-1 (mod 6)`, ascending, in `O(depth)` memory. It
must carry `c`'s own pentagon deletion on the root frame, since a calibrated arc
is seeded at a cell that may be a pentagon. It declares `SizeUnknown()`: the
closed-form rim census counts the `(6, 0)` walk and does not describe a seeded
one.

Neither validates `c` — both are called from `hex_halo_engine` on cells that
came out of `neighbors` and `children`, and the driver owns the level guard.
"""
function hex_child_direction end

@doc (@doc hex_child_direction)
function seeded_rim_engine end

"""
    descendant_range(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex, l::Integer) -> UnitRange{Int}

The contiguous interval of **positions** in `levelgrid(sys, l)`'s canonical
dense order occupied by the descendants of `c` at level `l`.

Available only when [`has_sorted_subtrees(sys)`](@ref has_sorted_subtrees) is
`true`; otherwise there is no method and the call is a `MethodError`.

Both directions are required:

 1. every level-`l` descendant of `c` has its position in the range, and
 2. every position in the range is a level-`l` descendant of `c`.

The range is over valid dense positions, not raw ids, and therefore has no id
encoding gaps. It can be intersected with sorted position vectors by binary
search.

Sibling ranges are disjoint and partition the parent's range in canonical order.

`l < level(c)` throws an `ArgumentError`.
"""
function descendant_range end
