# ---------------------------------------------------------------------------
# Hierarchical-grid fast-path contracts. Arguments use `(sys, c)` with any
# target level last; cell ids encode their own level. Overrides must preserve
# base-interface semantics.
# ---------------------------------------------------------------------------

# ===========================================================================
# The implementor surface. Each entry states whether it is required or
# defaulted; `AbstractHierarchicalGridSystem` tabulates the split.
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
[`maxlevel`](@ref).

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
whose grid genuinely carries state beyond `(sys, l)`; none of the seven shipped
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
    method itself — except for an id in a scheme the system does not claim at
    all, which is an `ArgumentError` naming both systems.

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

Calling it on a cell at `maxlevel(sys)` throws an `ArgumentError`.
"""
function children end

"""
    node_extent(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex) -> GO.UnitSpherical.SphericalCap

The covering region of the entire subtree rooted at `c`.

**Defaulted**: the generic implementation inflates the cell's own bounding cap
by [`cap_inflation(sys)`](@ref cap_inflation), and a system able to compute a
tighter covering cap overrides it. The covering law below holds either way, and
validating it is the system's responsibility.

> `node_extent(sys, c)` contains the geometry of **every descendant of `c`, at
> every depth** — every point of every cell boundary in the subtree, all the
> way down to `maxlevel(sys)`.

Tree pruning depends on this covering law. Over-coverage only reduces pruning;
under-coverage can omit valid results. For a convex cap of angular radius at
most 90°, containing all boundary vertices also contains their great-circle
arcs. Non-convex extents must establish containment of the full geometry.

A cell's own boundary need not cover its descendants; aperture-7 children, for
example, can extend beyond the parent boundary.

The result must cover through `maxlevel(sys)`, independent of a caller's
planned traversal depth.

Extents are `SphericalCap`s throughout this package, at every node of every
tree, which is what lets one predicate vocabulary serve all of them.
"""
function node_extent end

"""
    maxneighbors(sys::AbstractHierarchicalGridSystem, connectivity::Connectivity = Vertex()) -> Union{Int,Nothing}

A **static** upper bound on the number of `connectivity`-neighbours of any cell
of `sys`, at any level, or `nothing` when the system declares no bound.

**Sizes the neighbourhood family.** An `Int` bound permits the fixed-capacity
stack containers behind [`neighbors`](@ref) and [`ring`](@ref) on a subset and
behind [`adjacency`](@ref); the complete-level verbs and the subtree family
never ask for it. Defaults to `nothing` — never a guessed
capacity — and the same machinery then buffers one-rings in a heap `Vector`,
allocating once per cell: identical answers, the slow path. Declaring the bound
is a speed decision, not a correctness requirement. Individual cells may have
fewer neighbours than the bound.
"""
maxneighbors(::AbstractHierarchicalGridSystem, ::Connectivity) = nothing

maxneighbors(sys::AbstractHierarchicalGridSystem) = maxneighbors(sys, Vertex())

# ===========================================================================
# Traits with defaults
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
    maxlevel(sys::AbstractHierarchicalGridSystem) -> Int

The deepest valid level of `sys`: `last(levels(sys))`, which is the default
implementation.
"""
maxlevel(sys::AbstractHierarchicalGridSystem) = last(levels(sys))

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
`level(c) <= l <= maxlevel(sys)`.

`descendants(sys, c, level(c))` is `[c]`. `l < level(c)` throws an
`ArgumentError` (uniformly across systems, so generic code can catch it).

This materializes `O(subtree)` ids. Use [`descendant_range`](@ref) when
available, or [`border`](@ref)`(subtree(sys, c, l))` when only the border is
needed.
"""
function descendants end

"""
    subtree(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex, l::Integer; bucket_size = 0) -> PartialGrid

`c`'s subtree at level `l`, as a region: the [`PartialGrid`](@ref) holding every
level-`l` descendant of `c`, rooted at `c`.

This is the one spelling of a subtree. It is what gives the subtree family the
same currency every other region has — `halo(subtree(sys, c, l))`,
`border(subtree(sys, c, l))`, `adjacency(subtree(sys, c, l); halo = 1)` — rather
than a parallel set of verbs taking `(sys, c, l)` argument tuples.

`l < level(c)` throws an `ArgumentError`. Construction is `O(1)` where
[`has_sorted_subtrees`](@ref) holds, since the ids are then the level grid's own
over a known position range; elsewhere it materialises
[`descendants`](@ref). `bucket_size` is [`PartialGrid`](@ref)'s.
"""
function subtree end

"""
    border_engine(sys, c, target::Int, connectivity)
    interior_engine(sys, c, target::Int, connectivity)
    halo_engine(sys, c, target::Int, connectivity)

Return the iteration engine used by `EdgeCellIterator`, `InnerCellIterator`, or
`SubtreeHaloIterator`. A system may override these methods with an `O(border)` or
`O(halo)` traversal; eager operations collect the same engine.

An engine is any iterator over `cellindextype(sys)`. All three methods own their
level validation, so their `ArgumentError`s are the ones the eager verbs raise.

Engine selection uses private dispatch on the system type and is not part of
the public compatibility surface.

The generic border and interior engines scan [`descendant_range`](@ref), or
materialize descendants when [`has_sorted_subtrees`](@ref) is `false`. The
generic halo engine walks cells outside the subtree because a halo is not a
single descendant interval.

Specializations may enumerate conservative candidates, but must filter them by
the requested adjacency. If their preconditions cannot be verified, they must
return `generic_halo_engine(sys, c, target, connectivity)`.
"""
function border_engine end

@doc (@doc border_engine)
function interior_engine end

@doc (@doc border_engine)
function halo_engine end

"""
    lattice_decode(sys, c) -> (ix, iy, face)
    lattice_cell(sys, level::Int, ix, iy, face) -> cell
    face_orientation(sys, face) -> UInt8

Define the face-lattice operations used by the shared square halo traversal.
Only systems with an aligned square lattice per face implement these methods.

`lattice_decode` and `lattice_cell` convert between cell ids and face-local
coordinates, with `face` 0-based and
`(ix, iy)` in `0:2^level - 1`. `face_orientation` is the curve state a face's
root uses before consuming position bits.

The traversal derives seam rectangles by decoding neighbours of border cells.

`SquareBandEngine` also requires these invariants:

 1. `cellindextype(sys) === LevelIndex`. The engine emits `LevelIndex` and
    declares it as its `eltype`, unconditionally.
 2. Ids follow `face * faceside^2 + curvecode`, with 0-based face and curve
    code. The emit step constructs ids with this arithmetic.
 3. Interior face adjacency is the 3×3 lattice, so an in-face band requires no
    additional adjacency check.

A fourth square system holding all three subtypes
[`AbstractQuadFaceGridSystem`](@ref) and writes only the three methods above. One
that does not writes its own [`halo_engine`](@ref border_engine) instead.
"""
function lattice_decode end

@doc (@doc lattice_decode)
function lattice_cell end

@doc (@doc lattice_decode)
function face_orientation end

"""
    hex_child_direction(sys, c) -> Int
    seeded_border_engine(sys, c, target::Int, arclen::Int, start::Int)

Define the operations used by the calibrated aperture-7 halo traversal. H3 and
IGeo7 implement them using their subtree-border automata.

`hex_child_direction` returns the position `0:5` of the parent-to-child step on
the direction ring, or `-1` for a centre child or root cell.

`seeded_border_engine` enters the system's border automaton at an arbitrary arc
`(arclen, start)` rather than at the fully exposed `(6, 0)` a subtree root gets:
`c`'s level-`target` descendants reachable along the arc of exposed directions
`start, start+1, …, start+arclen-1 (mod 6)`, ascending, in `O(depth)` memory. It
must carry `c`'s own pentagon deletion on the root frame, since a calibrated arc
is seeded at a cell that may be a pentagon. It declares `SizeUnknown()`: the
closed-form border census counts the `(6, 0)` walk and does not describe a seeded
one.

Neither method validates `c`; `hex_halo_engine` supplies cells returned by
`neighbors` and `children` and performs level validation.
"""
function hex_child_direction end

@doc (@doc hex_child_direction)
function seeded_border_engine end

"""
    descendant_range(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex, l::Integer) -> UnitRange{Int}

The contiguous interval of **positions** in `levelgrid(sys, l)`'s canonical
dense order occupied by the descendants of `c` at level `l`.

Available only when [`has_sorted_subtrees(sys)`](@ref has_sorted_subtrees) is
`true`; otherwise there is no method and the call is a `MethodError`. A system
that declares the trait and implements nothing gets an `ArgumentError` naming
both at the first call, rather than that `MethodError`.

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

# The trait's obligation, diagnosed where it is first broken. Without this the
# consumer sees a `MethodError` whose candidate list is other systems' methods.
function descendant_range(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        l::Integer)
    has_sorted_subtrees(sys) && throw(ArgumentError(
        "$(nameof(typeof(sys))) declares has_sorted_subtrees = true, which obliges " *
        "descendant_range(::$(typeof(sys)), ::$(typeof(c)), ::Integer); it is not implemented"))
    throw(MethodError(descendant_range, (sys, c, l)))
end
