# ---------------------------------------------------------------------------
# Core grid, system, and cell-index types. `level` and `rawid` are declared here
# for `LevelIndex`; their contracts are documented in `interface/grid.jl`.
# ---------------------------------------------------------------------------

function level end
function rawid end

"""
    abstract type AbstractGrid

One finite collection of unit-sphere cells: a complete DGGS level, a regional
subset, or a standalone structured grid. A grid need not cover the sphere;
[`treeify`](@ref)'s root extent defines its coverage.

# Index types and index spaces

A grid defines the canonical dense order `1:ncells(grid)`. The same cell can be
named three ways, differing in what kind of index it is and which space that
index counts in:

  - a **cell index** — an [`AbstractCellIndex`](@ref): a typed identity, and the
    only one of the three that names the cell with no collection at all;
  - a **local index** (a **storage index**) — a bare `Int` into this
    collection's own order, `1:ncells(grid)`;
  - a **global index** — a bare `Int` into the complete grid at that level.

On a complete grid the local and global indices of a cell coincide, because its
storage IS the level. On a subset they do not, and code that confuses them reads
one cell's data for another. Never store a bare `Int` without knowing which of
the two it is.

The canonical order is the grid's own choice, but it must be stable for the
lifetime of the grid object and consistent with [`cellindex`](@ref) /
[`localindex`](@ref), which are inverses of each other over it.

# Required interface

A grid type writes exactly four methods:

| method | contract |
|---|---|
| [`ncells(grid)`](@ref ncells) | number of cells |
| [`cellindex(grid, i)`](@ref cellindex) | local index → canonical cell index |
| [`cell_boundary(grid, c)`](@ref cell_boundary) | exact boundary ring, unit-sphere points |
| [`cell_centroid(grid, c)`](@ref cell_centroid) | representative interior point |

Everything else—[`localindex`](@ref), [`cell_polygon`](@ref),
[`cell_area`](@ref), [`cell_extent`](@ref), [`getcell`](@ref),
[`cellat`](@ref), [`neighbors`](@ref), [`ring`](@ref), [`treeify`](@ref),
and [`query`](@ref)—is provided generically and may be optimized without
changing semantics. The generic [`localindex`](@ref) scans `1:ncells(grid)`
linearly, so a grid that can search should override it.

A hierarchical system needs no grid type: it answers this interface with the
five level-grid primitives listed under
[`AbstractHierarchicalGridSystem`](@ref), which [`HierarchicalLevelGrid`](@ref)
forwards to. That list is five rather than four because a linear scan is not an
acceptable [`globalindex`](@ref) for a complete level.

A grid produced by a hierarchical system reports it through [`system`](@ref)
and [`level`](@ref); a standalone grid returns `nothing` from both and stops at
the base interface.

See also [`AbstractHierarchicalGridSystem`](@ref), [`AbstractCellIndex`](@ref).
"""
abstract type AbstractGrid end

# A grid is one argument to every cell verb, never a container to iterate over:
# `cell_centroid.(grid, cells)` broadcasts the cells, not the grid.
Base.broadcastable(g::AbstractGrid) = Ref(g)

"""
    abstract type AbstractHierarchicalGridSystem

A family of grid levels related by analytic parent/child structure. The system
names cells; [`levelgrid`](@ref) returns a complete level grid. Hierarchical
methods accelerate base-grid operations without changing their results.

# Required interface

Identity and hierarchy:

| method | contract |
|---|---|
| [`cellindextype(sys)`](@ref cellindextype) | canonical typed id type |
| [`levels(sys)`](@ref levels) | valid level range |
| [`rootcells(sys)`](@ref rootcells) | top-level cells, ascending |
| [`Base.parent(sys, c)`](@ref parent) | analytic parent, no lookup tables |
| [`children(sys, c)`](@ref children) | analytic children, no lookup tables |

`parent` is a method on `Base.parent`, so it is spelled `Base.parent(sys, c)`
at the definition site and reached unqualified from any session.

The five level-grid primitives, which answer the [`AbstractGrid`](@ref)
interface for the complete level:

| method | contract |
|---|---|
| [`ncells(sys, l)`](@ref ncells) | number of cells at level `l` |
| [`cellindex(sys, l, i)`](@ref cellindex) | index → canonical typed id |
| [`globalindex(sys, c)`](@ref globalindex) | id → global index, or `nothing` |
| [`cell_boundary(sys, c)`](@ref cell_boundary) | exact boundary ring, unit-sphere points |
| [`cell_centroid(sys, c)`](@ref cell_centroid) | representative interior point |

[`maxneighbors(sys, connectivity)`](@ref maxneighbors) **sizes the
neighbourhood family**: [`neighbors`](@ref), [`ring`](@ref) on a subset and
[`adjacency`](@ref) use it for their fixed-capacity containers. It defaults to `nothing` — no bound declared — and
the same verbs then buffer in heap `Vector`s: identical answers, one
allocation per cell. Declaring the bound is the fast path.

# Defaults an implementor may override

| method | default |
|---|---|
| [`levelgrid(sys, l)`](@ref levelgrid) | `HierarchicalLevelGrid(sys, l)`, checked against [`levels`](@ref) |
| [`node_extent(sys, c)`](@ref node_extent) | the cell's bounding cap, inflated — covers descendant *geometry*, not descendant caps; see the covering law |
| [`cap_inflation(sys)`](@ref cap_inflation) | `1.2` |
| [`maxlevel(sys)`](@ref maxlevel) | `last(levels(sys))` |
| [`has_sorted_subtrees(sys)`](@ref has_sorted_subtrees) | `false`; declaring it `true` obliges [`descendant_range`](@ref) |
| [`has_congruent_refinement(sys)`](@ref has_congruent_refinement) | `false`; `true` asserts that children tile their parent |
| [`has_direct_location(sys)`](@ref has_direct_location) | `false`; declaring it `true` obliges [`cellat`](@ref) on the level grid |

# Grid methods and system methods

Identity and geometry dispatch on the **system**: the tables above are all
`(sys, ...)` methods, and the level-grid primitives are what
[`HierarchicalLevelGrid`](@ref) forwards the base interface to.

Everything else dispatches on the **grid**. A system's fast paths —
[`cellat`](@ref), [`neighbors`](@ref), [`ring`](@ref), [`cell_area`](@ref),
[`treeify`](@ref), the subtree engines — attach to
`HierarchicalLevelGrid{typeof(sys)}`, for which each system here keeps a local
alias:

    import DiscreteGlobalGrids as DGG
    const LevelGrid = DGG.HierarchicalLevelGrid{MySystem}

    DGG.cellat(g::LevelGrid, p::DGG.UnitSphericalPoint) = ...

Every fast path is optional, and must return what the generic implementation
would have returned.

See also [`AbstractGrid`](@ref), [`AbstractQuadFaceGridSystem`](@ref),
[`node_extent`](@ref).
"""
abstract type AbstractHierarchicalGridSystem end

# As for a grid: a system broadcasts as one value, so `ancestor.(sys, cells, 3)`
# needs no `Ref`.
Base.broadcastable(sys::AbstractHierarchicalGridSystem) = Ref(sys)

"""
    abstract type AbstractQuadFaceGridSystem <: AbstractHierarchicalGridSystem

A system whose cells are an aligned `2^level × 2^level` lattice on each of
`nbasefaces(sys)` congruent faces, named by the dense 0-based id
`face * 4^level + curvecode` and indexed at `id + 1`. S2, HEALPix and ISEA4R
are that one family.

A subtype inherits every method that identity alone determines: the hierarchy
block ([`rootcells`](@ref), [`parent`](@ref), [`children`](@ref),
[`ancestor`](@ref), [`descendant_range`](@ref), [`descendants`](@ref)), the
level-grid arithmetic ([`ncells`](@ref), [`cellindex`](@ref),
[`globalindex`](@ref), [`cellindextype`](@ref),
[`has_sorted_subtrees`](@ref)), and the subtree engines
([`border_engine`](@ref)/`interior_engine`/`halo_engine`), which read a subtree as
the square lattice block it is.

It must still provide its own face layout and projection: [`levels`](@ref),
[`maxneighbors`](@ref), [`cell_boundary`](@ref), [`cell_centroid`](@ref),
[`node_extent`](@ref), [`cellat`](@ref), [`one_ring`](@ref), the lattice codec
hooks [`lattice_decode`](@ref)/`lattice_cell`/`face_orientation`, and the three
declarations `nbasefaces`, `systemname`, `idname`. A system whose curve carries
orientation state also overrides `subtree_curve` and `subtree_orientation`.

See also [`AbstractHierarchicalGridSystem`](@ref).
"""
abstract type AbstractQuadFaceGridSystem <: AbstractHierarchicalGridSystem end

"""
    abstract type AbstractCellIndex

A typed, self-describing name for one cell of one system.

A cell index is an identity, not an offset into any collection, and names the
same cell in complete grids, subsets, or without a grid at all.

# Required of every subtype

  - **isbits.** Indices must be immutable and allocation-free.
  - **[`level(c)`](@ref level) is total.** Every id encodes its level.
  - **[`rawid(c)`](@ref rawid)** returns the encoded integer.
  - **A total order.** `Base.isless` must implement canonical cell order, with
    consistent `==` and `hash`.

One scheme per system is canonical ([`cellindextype`](@ref)); alternates are
reached through [`reindex`](@ref) and listed by [`cellindextypes`](@ref).
"""
abstract type AbstractCellIndex end

"""
    LevelIndex(level, index) <: AbstractCellIndex

An id consisting of an explicit `level` and a system-defined linear `index`.

The system documents whether `index` is zero- or one-based. It must increase in
canonical cell order; `LevelIndex` orders lexicographically by `(level, index)`.

Here `index` is an identity component, not an offset into a collection. Use
[`globalindex`](@ref) for the index in the complete grid and
[`localindex`](@ref) for the index in a subset's own storage.

Construction does not validate level or index ranges. Validation occurs when an
id is used with a system or grid.

```jldoctest
julia> c = LevelIndex(3, 17); (level(c), rawid(c))
(3, 17)
```
"""
struct LevelIndex <: AbstractCellIndex
    # `Int32` rather than `Int8`: a level is compared and incremented against
    # `Int`s everywhere, and an 8-bit field would silently wrap at 127 for a
    # system with no `maxlevel`. The struct is 8-byte aligned either way, so
    # the two narrower choices cost exactly the same 16 bytes.
    level::Int32
    index::Int64
end

level(c::LevelIndex) = Int(c.level)
rawid(c::LevelIndex) = c.index

Base.isless(a::LevelIndex, b::LevelIndex) =
    isless((a.level, a.index), (b.level, b.index))

Base.show(io::IO, c::LevelIndex) = print(io, "LevelIndex(", c.level, ", ", c.index, ")")

"""
    abstract type AbstractCellVector{ID} <: AbstractVector{ID}

An ascending collection of cells from one system and one level, addressed by
local index. This is the **region** contract in vector form: the four region verbs
([`halo`](@ref), [`border`](@ref), [`interior`](@ref), [`adjacency`](@ref)),
the neighbourhood sweeps, regridding and plotting are all written against it,
so a new backing gets the whole surface by subtyping rather than by
reimplementing it.

Two backings ship. [`CellVector`](@ref) COMPUTES its ids from compressed
index windows; [`ChunkedCellVector`](@ref) reads the ids a store WROTE, from
its chunk manifest. What differs is where element `k` comes from and what it
costs, not what it means.

# Required interface

| method | contract |
|---|---|
| [`system(cv)`](@ref system) | the grid system the cells belong to |
| [`level(cv)`](@ref level) | the single level they are all at |
| `Base.size(cv)` | `(n,)` — the number of cells |
| `Base.getindex(cv, k::Int)` | local index `k` → typed cell id, `1`-based |
| [`localindex(cv, c)`](@ref localindex) | cell id → local index, or `nothing` |

`getindex` and `localindex` are inverses over `1:length(cv)`, and the ids
they range over are strictly ascending. A subtype that cannot promise ascent is
not a cell vector; it is an unordered list of cells.

# Cost is a property of the backing, not of the contract

`localindex` is closed-form arithmetic on one backing and may decode a stored
chunk on another. Generic code that resolves indices in an order the backing
did not choose is correct on both and cheap on only one, which is what
[`chunkplan`](@ref) exists to fix: it names the traversal order that keeps a
chunk-backed vector reading each chunk once.

See also [`AbstractCellLookup`](@ref), the `DimensionalData` face of the same
contract, and [`PartialGrid`](@ref), the grid-shaped sibling.
"""
abstract type AbstractCellVector{ID} <: AbstractVector{ID} end

"""
    abstract type Connectivity

How two cells must meet to count as adjacent. See [`neighbors`](@ref).

Concrete singletons: [`Vertex`](@ref) (Moore, the default) and [`Edge`](@ref)
(von Neumann).
"""
abstract type Connectivity end

"""
    Vertex() <: Connectivity

Moore connectivity: cells are adjacent when they share at least a vertex. This
is the default connectivity.

`Vertex()` and [`Edge()`](@ref Edge) coincide where exactly three cells meet at
each vertex. Higher-valence vertices add corner-only neighbours. This includes
A5's 4-valent corners and ISEA4R's 4-valent lattice corners and 5-valent
icosahedral vertices.
"""
struct Vertex <: Connectivity end

"""
    Edge() <: Connectivity

von Neumann connectivity: two cells are adjacent only if they share a whole
**edge**. The opt-in restriction of the [`Vertex()`](@ref Vertex) default.
"""
struct Edge <: Connectivity end

"""
    abstract type Winding

The order a system's [`one_ring`](@ref) arrives in, declared rather than
promised in prose. See [`winding`](@ref).

**What reads it.** The shell walk behind `neighbors(grid, c, k)` and
`ring(grid, c, k)` for `k >= 2` has to know the rotational order of each ring.
A declared turn ([`CounterClockwise`](@ref), [`Clockwise`](@ref)) lets it
propagate that order outward from the one-rings it already has. Otherwise it
measures the order geometrically instead — a [`cell_centroid`](@ref) per cell of
every ring, plus a sort.

Concrete singletons: [`CounterClockwise`](@ref), [`Clockwise`](@ref),
[`CustomOrder`](@ref) and [`Unordered`](@ref).
"""
abstract type Winding end

"""
    CounterClockwise() <: Winding

`one_ring` is one counter-clockwise turn seen from **outside** the sphere,
starting at the system's own start direction. The order [`neighbors`](@ref)
states, and the one every system in this package declares.
"""
struct CounterClockwise <: Winding end

"""
    Clockwise() <: Winding

`one_ring` is one clockwise turn seen from outside the sphere. A rotational
winding like [`CounterClockwise`](@ref), read in the other direction: the shell
walk reverses it and is otherwise unchanged, so declaring this costs nothing
against declaring the counter-clockwise turn.
"""
struct Clockwise <: Winding end

"""
    CustomOrder() <: Winding

`one_ring` has a deterministic order the shell walk may **not** carry outward.
Callers may rely on it being stable between calls; `k >= 2` is measured by
azimuth as under [`Unordered`](@ref).

Weaker than [`CounterClockwise`](@ref), and not the same as having no turn.
A5 is the case: its one-rings *are* counter-clockwise when measured, but its
shells are not rotational copies of them — its rings grow `8, 18, 29, 39` rather
than linearly — so there is no outward order to propagate and the geometric sort
is the answer rather than a fallback. A system whose one-ring is genuinely
unsorted wants [`Unordered`](@ref) instead.
"""
struct CustomOrder <: Winding end

"""
    Unordered() <: Winding

`one_ring` promises no order at all, not even stability between calls. The
weakest declaration, and the safe default for a system that has not stated
otherwise. The shell walk measures azimuth, as it does for
[`CustomOrder`](@ref).
"""
struct Unordered <: Winding end
