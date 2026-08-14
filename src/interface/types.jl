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

# Position vs identity

A grid defines the canonical dense order `1:ncells(grid)`. A bare `Int` is a
position in that order; an [`AbstractCellIndex`](@ref) is a typed cell identity.

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

Everything else—[`cellposition`](@ref), [`cell_polygon`](@ref),
[`cell_area`](@ref), [`cell_extent`](@ref), [`getcell`](@ref),
[`cellat`](@ref), [`neighbors`](@ref), [`ring`](@ref), [`treeify`](@ref),
and [`query`](@ref)—is provided generically and may be optimized without
changing semantics.

A grid produced by a hierarchical system reports it through [`system`](@ref)
and [`level`](@ref); a standalone grid returns `nothing` from both and stops at
the base interface.

See also [`AbstractHierarchicalGridSystem`](@ref), [`AbstractCellIndex`](@ref).
"""
abstract type AbstractGrid end

"""
    abstract type AbstractHierarchicalGridSystem

A family of grid levels related by analytic parent/child structure. The system
names cells; [`levelgrid`](@ref) returns a complete level grid. Hierarchical
methods accelerate base-grid operations without changing their results.

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

[`has_sorted_subtrees`](@ref), [`cap_inflation`](@ref),
[`max_neighbors`](@ref), and [`max_level`](@ref).

See also [`AbstractGrid`](@ref), [`node_extent`](@ref).
"""
abstract type AbstractHierarchicalGridSystem end

"""
    abstract type AbstractCellIndex

A typed, self-describing name for one cell of one system.

A cell index is an identity, not a grid position, and names the same cell in
complete grids, subsets, or without a grid.

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

`index` is an identity component, not a grid position. Use
[`cellposition`](@ref) for the position in a complete or partial grid.

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
