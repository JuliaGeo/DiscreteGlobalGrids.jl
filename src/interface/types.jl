# ---------------------------------------------------------------------------
# The type vocabulary of the grid interface.
#
# Three abstract types carry the whole design: a grid is one finite collection
# of cells, a system is a hierarchy that produces grids, and a cell index is a
# typed name for a cell. Everything else in `interface/` is a generic function
# written against these three.
#
# `level` and `rawid` are declared (bodiless) here because `LevelIndex`, the one
# concrete cell index this layer ships, implements them immediately below. Their
# contracts are documented in `interface/grid.jl`, next to the rest of the
# identity plumbing.
# ---------------------------------------------------------------------------

function level end
function rawid end

"""
    abstract type AbstractGrid

One finite collection of cells on the unit sphere: a complete DGGS level, a
regional subset of one, or a standalone structured grid (tripolar, gaussian)
with no hierarchy at all.

# Coverage

A grid is **not** assumed to partition the sphere. Partial coverage is a
property of the base interface, not an exception to it: nothing generic may
assume that every point of the sphere lies in some cell, and there is no
`coverage_extent` primitive — coverage is whatever [`treeify`](@ref)'s root
extent says it is.

# Position vs identity

A grid imposes a **canonical dense order** on its cells, `1:ncells(grid)`. A
bare `Int` passed to any function in this package is always a *position* in
that order — the storage coordinate that data arrays, lookups and regridding
matrices are laid out against. A typed [`AbstractCellIndex`](@ref) is always an
*identity*: a name relative to a system, meaningful with no grid in hand. Ids
are never bare integers, so `f(grid, i::Int)` and `f(grid, c::AbstractCellIndex)`
coexist unambiguously and mean different things.

The canonical order is the grid's own choice, but it must be stable for the
lifetime of the grid object and consistent with [`cellindex`](@ref) /
[`cellposition`](@ref), which are inverses of each other over it.

# Required interface

An implementor writes exactly four methods:

| method | contract |
|---|---|
| [`ncells(grid)`](@ref ncells) | number of cells |
| [`cellindex(grid, i)`](@ref cellindex) | position → canonical typed id |
| [`cell_boundary(grid, c)`](@ref cell_boundary) | exact boundary ring, unit-sphere points |
| [`cell_centroid(grid, c)`](@ref cell_centroid) | representative interior point |

Everything else — [`cellposition`](@ref), [`cell_polygon`](@ref),
[`cell_area`](@ref), [`cell_extent`](@ref), [`getcell`](@ref),
[`cellat`](@ref), [`neighbors`](@ref), [`ring`](@ref), [`treeify`](@ref),
[`query`](@ref) — is provided generically and may be overridden for speed,
never for different semantics.

A grid produced by a hierarchical system reports it through [`system`](@ref)
and [`level`](@ref); a standalone grid returns `nothing` from both and stops at
the base interface.

See also [`AbstractHierarchicalGridSystem`](@ref), [`AbstractCellIndex`](@ref).
"""
abstract type AbstractGrid end

"""
    abstract type AbstractHierarchicalGridSystem

A hierarchy of grids: a family of levels related by analytic parent/child
structure. The system is the object that *names* cells; a grid is one level of
it (or a subset of one level), obtained with [`levelgrid`](@ref).

Hierarchy is the fast path, never the contract. Every algorithm in this package
is written against [`AbstractGrid`](@ref) first; a system's methods let the
generic code prune trees, walk subtrees and answer queries in sublinear time,
but they may not change what the answers *are*.

# Required interface

| method | contract |
|---|---|
| [`cellindextype(sys)`](@ref cellindextype) | canonical typed id type |
| [`levels(sys)`](@ref levels) | valid level range |
| [`levelgrid(sys, l)`](@ref levelgrid) | complete `AbstractGrid` at level `l` |
| [`rootcells(sys)`](@ref rootcells) | top-level cells |
| [`parent(sys, c)`](@ref parent) | analytic parent, no lookup tables |
| [`children(sys, c)`](@ref children) | analytic children, no lookup tables |
| [`node_extent(sys, c)`](@ref node_extent) | covering region — see the covering law |

# Traits

[`has_sorted_subtrees`](@ref) (default `false`), [`cap_inflation`](@ref)
(default `1.2`), [`max_neighbors`](@ref), [`max_level`](@ref). Descriptive
metadata (cell shape, aperture, equal-areaness) is never load-bearing for an
algorithm and lives outside this interface.

See also [`AbstractGrid`](@ref), [`node_extent`](@ref).
"""
abstract type AbstractHierarchicalGridSystem end

"""
    abstract type AbstractCellIndex

A typed, self-describing name for one cell of one system.

A cell index is an *identity*, not a position: it means the same cell whether
it is read against a complete level grid, a regional subset, or no grid at all.
Bare integers are never cell indices — see [`AbstractGrid`](@ref) on position
vs identity.

# Required of every subtype

  - **isbits.** Cell indices are stored in dense arrays and small stack
    containers by the million; they must be immutable and allocation-free.
    Canonical schemes are an 8-byte wrapper with the level encoded in-band
    (`Z7Cell`, `H3Cell`) or a level field plus a linear index
    ([`LevelIndex`](@ref)).
  - **[`level(c)`](@ref level) is total.** Every id knows its own level with no
    system, grid or table in hand. This is why the interface never needs the
    `(system, level, id)` argument triple.
  - **[`rawid(c)`](@ref rawid)** returns the encoded integer.
  - **A total order.** `Base.isless` must implement the system's canonical cell
    order, and `==`/`hash` must agree with it. Generic code sorts results by id
    and binary-searches sorted id vectors; "unspecified order" is never allowed
    anywhere in this package.

One scheme per system is canonical ([`cellindextype`](@ref)); alternates are
reached through [`reindex`](@ref) and listed by [`cellindextypes`](@ref).
"""
abstract type AbstractCellIndex end

"""
    LevelIndex(level, index) <: AbstractCellIndex

The canonical id of a system whose cells at each level are a dense linear
range: an explicit `level` field plus the linear `index` of the cell within
that level. Nested HEALPix is the archetype.

`index` is the system's own linear numbering, and the system documents whether
it is 0- or 1-based (a `LevelIndex` is agnostic; it only requires that the
numbering be strictly increasing in canonical cell order, so that
`isless(::LevelIndex, ::LevelIndex)` — lexicographic in `(level, index)` — *is*
that order).

Note that `index` is an identity, not a grid position: for a complete level
grid the two coincide up to the base offset, but for a
subset grid they do not, and [`cellposition`](@ref) is the only way to go from
one to the other.

```jldoctest
julia> c = LevelIndex(3, 17);

julia> level(c), rawid(c)
(3, 17)

julia> LevelIndex(3, 17) < LevelIndex(3, 18) < LevelIndex(4, 0)
true
```
"""
struct LevelIndex <: AbstractCellIndex
    # `Int32` rather than `Int8`: a level is compared and incremented against
    # `Int`s everywhere, and an 8-bit field would silently wrap at 127 for a
    # system with no `max_level`. The struct is 8-byte aligned either way, so
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
    abstract type Connectivity

How two cells must meet to count as adjacent. See [`neighbors`](@ref).

Concrete singletons: [`Vertex`](@ref) (Moore, the default) and [`Edge`](@ref)
(von Neumann).
"""
abstract type Connectivity end

"""
    Vertex() <: Connectivity

Moore connectivity: two cells are adjacent if they share **at least a vertex**.

This is the default everywhere in this package, because it is the superset a
consumer cannot reconstruct from an edge-only answer — dropping the corner
neighbours of a `Vertex()` result is one filter, while recovering them from an
`Edge()` result needs the grid again.

On hexagonal and pentagonal grids `Vertex()` and [`Edge()`](@ref Edge)
coincide; on quadrilateral grids (HEALPix, S2) `Vertex()` adds the four corner
neighbours.
"""
struct Vertex <: Connectivity end

"""
    Edge() <: Connectivity

von Neumann connectivity: two cells are adjacent only if they share a whole
**edge**. The opt-in restriction of the [`Vertex()`](@ref Vertex) default.
"""
struct Edge <: Connectivity end
