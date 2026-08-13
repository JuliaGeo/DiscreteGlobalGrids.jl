# ---------------------------------------------------------------------------
# The hierarchical grid system interface.
#
# Everything here is a fast path. A system's methods let generic code prune
# trees, walk subtrees, and answer queries in sublinear time; none of them may
# change what an answer *is*. If a system's override and the base-interface
# fallback disagree, the system is wrong.
#
# The argument convention is `(sys, c)` — a system and a typed cell id, and the
# target level last where there is one (`ancestor(sys, c, l)`). A cell id knows
# its own level (`level(c)` is total), so the old `(system, level, id)` triple
# is gone, and with it every way for a level and an id to disagree.
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

A system with no intrinsic depth limit still returns a bounded range — the
deepest level its id encoding and its `Int` cell counts remain valid at. That
bound is a fact about the implementation, and stating it is what keeps
`max_level` total and keeps `radix^level` arithmetic one comparison away from
silent overflow rather than one multiply away from it.
"""
function levels end

"""
    levelgrid(sys::AbstractHierarchicalGridSystem, l::Integer) -> AbstractGrid

The **complete** grid of `sys` at level `l`: every cell the system has at that
level, in the system's canonical dense order.

**Required.**

`system(levelgrid(sys, l)) === sys` and `level(levelgrid(sys, l)) == l`. The
grid is a lightweight descriptor, not a materialised cell list — constructing
one must not be O(cells).

`l` outside [`levels(sys)`](@ref levels) throws an `ArgumentError`. For a
subset of a level, see `PartialGrid`.
"""
function levelgrid end

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

The count is the system's aperture for most cells, and *less* for the
exceptional ones: a pentagon in an aperture-7 icosahedral system has six
children, not seven. Generic code must never assume a fixed child count; that
is what makes the pentagon cases work by construction rather than by patch.

Calling it on a cell at `max_level(sys)` throws an `ArgumentError`.
"""
function children end

"""
    node_extent(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex) -> GO.UnitSpherical.SphericalCap

The covering region of the whole subtree rooted at `c`.

**Required.**

# The covering law

> `node_extent(sys, c)` contains the geometry of **every descendant of `c`, at
> every depth** — every point of every cell boundary in the subtree, all the
> way down to `max_level(sys)`.

(For a geodesically convex extent — a spherical cap with angular radius at most
90° — containing a cell's boundary *vertices* implies containing the great-arc
edges between them, which is why the conformance suite may sample vertices as
its proxy. An implementation whose extents are not convex owes the full law,
not the proxy.)

This is the contract that makes generic tree pruning **correct**. Every
traversal in this package — [`query`](@ref), [`cellat`](@ref), the
`ConservativeRegridding` dual tree walk — discards a node without looking
inside it the moment its extent misses the target. If the law is violated by
one cell at one depth, that cell is silently dropped from an answer, and
nothing downstream can detect it. It is therefore a **property-tested
contract**: `DiscreteGlobalGridsConformanceTesting.test_hierarchical_system(sys)`
samples cells, walks their subtrees several levels down, and asserts
containment.

Over-covering is always safe and only ever costs time. Under-covering is a
correctness bug.

# Implementing it

A cell's own exact boundary is **not** in general a valid node extent. Under
aperture 7 the children of a hexagon poke out past their parent's edges, so a
tight parent polygon violates the law at depth 1 — this is why the old
"congruent geometry" trait was a trap, and why there is no way to opt out of
covering here.

The generic default is the cell's bounding cap inflated by
[`cap_inflation(sys)`](@ref cap_inflation), which is sound for every system in
scope and O(1) at every depth. Systems that can do better override:

  - **HEALPix** returns the exact subtree cap — nested children are contained
    in their parent, so the tight cap is already legal, and it is still O(1).
  - Systems whose overhang is well characterised lower `cap_inflation` rather
    than overriding this.

Note that the answer must **not** depend on how deep the caller intends to
descend. There is no leaf-level argument: an extent that covers the subtree to
depth `d` but not to `d + 1` is not a node extent.

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

`false` by default — a system opts in by declaring it and implementing
[`descendant_range`](@ref).

This is the single trait that used to be spelled two ways (`supports_prefix_ranges`
for the declaration, `has_descendant_ranges` for the wiring). It is worth a lot
when true: subtree membership becomes two `searchsorted` calls against a sorted
id vector, a tree cursor can carry a position *window* instead of a
materialised selection, and a multi-order cell set can be expanded to any level
as sorted disconnected ranges.

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

The number is a property of the system's refinement geometry: how far past its
parent's bounding cap the deepest descendant can reach, in the limit. It is
measured, not guessed, and it belongs in the system's own test suite as a
sampled bound. Systems whose children overhang more must raise it; systems that
override `node_extent` outright ignore it.

Raising it costs query time (looser pruning). Setting it too low is a
correctness bug — see the covering law in [`node_extent`](@ref).
"""
cap_inflation(::AbstractHierarchicalGridSystem) = 1.2

"""
    max_neighbors(sys::AbstractHierarchicalGridSystem, connectivity::Connectivity = Vertex()) -> Int

A **static** upper bound on the number of `connectivity`-neighbours of any cell
of `sys`, at any level.

Static is the point: it is a property of the types alone, so
[`neighbors`](@ref) can size a `SmallCollections.SmallVector` at compile time
and sweep a grid without allocating. Individual cells may have fewer — a
pentagon has five neighbours where hexagons have six, and 24 HEALPix pixels
have seven where the rest have eight; the bound covers the maximum, and the
container is variable-length below it.

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

This is O(subtree) and materialises: the result of asking for level 20
descendants of a root cell is not a thing that fits in memory. Reach for
[`descendant_range`](@ref) when the system supports it, and for the subtree
border/interior iterators when only the rim is wanted.
"""
function descendants end

"""
    descendant_range(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex, l::Integer) -> UnitRange{Int}

The contiguous interval of **positions** in `levelgrid(sys, l)`'s canonical
dense order occupied by the descendants of `c` at level `l`.

Available only when [`has_sorted_subtrees(sys)`](@ref has_sorted_subtrees) is
`true`; otherwise there is no method and the call is a `MethodError`.

# The two-sided contract

Both directions hold, and generic code depends on both:

 1. every level-`l` descendant of `c` has its position in the range, and
 2. every position in the range is a level-`l` descendant of `c`.

So intersecting the range with any *sorted* vector of level-`l` positions
yields exactly the descendants present, by two binary searches and no
per-cell parent walk. Because these are positions among valid cells and not raw
id bounds, the interval has no holes to skip: pentagon id gaps are never in it.

Sibling ranges are disjoint and, taken in order, partition the parent's range —
which is what makes ordering a multi-order cell set by `first(descendant_range)`
at a fixed reference depth equal to depth-first curve order.

`l < level(c)` throws an `ArgumentError`.
"""
function descendant_range end
