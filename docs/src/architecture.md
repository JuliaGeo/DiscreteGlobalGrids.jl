# Architecture

This page describes the abstractions in DiscreteGlobalGrids: what each one is, what it
demands of an implementor, and why it exists. It is written for someone who has not read
the source.

A discrete global grid system tessellates the sphere into cells, usually recursively, so
that every cell has an identifier and a geometry. The package implements several of them and
makes them answerable through one set of verbs while preserving the closed-form
arithmetic each one has of its own.

## Grids and systems

The two root types are disjoint: neither is a subtype of the other.

An `AbstractGrid` is a finite, concrete collection of cells. It has a length, and index
`i` means a local index — a place in `1:ncells(grid)`. An `AbstractHierarchicalGridSystem` is a rule for
tessellating the sphere, with analytic parent and child structure. It is not a collection
and has no length.

This split lets the package accept grids that no hierarchy produced — a tripolar
ocean grid, a cubed sphere, anything whose neighbours must be derived geometrically. Such a
grid supplies four methods and gets tree building, point location, adjacency, and geometry
from generic fallbacks. A grid that *is* backed by a system keeps those fallbacks as a
floor and overrides whatever it can do faster.

A grid must implement exactly four things: `ncells`, `cellindex`, `cell_boundary`, and
`cell_centroid`. Everything else in the grid interface has a working default. A system must
implement `cellindextype`, `levels`, `rootcells`, `parent`, `children`, and five level-grid
primitives described below.

## Level grids

No shipped system defines a grid type of its own. [`levelgrid`](@ref) returns a
`HierarchicalLevelGrid`, which holds nothing but the system and the level, and forwards the
four grid primitives to five system-level ones: `ncells(sys, l)`, `cellindex(sys, l, i)`,
`globalindex(sys, c)`, `cell_boundary(sys, c)`, and `cell_centroid(sys, c)`.

The additional `globalindex` method gives a complete level an efficient id-to-position
operation. Its generic grid version scans the cells, which is suitable for small or exotic
grids; hierarchical systems usually provide a faster implementation.

Everything beyond those four — `cellat`, `neighbors`, `ring`, `cell_area` — a system
attaches to `HierarchicalLevelGrid{ItsOwnSystem}`. Fast paths dispatch on the type parameter,
while grid construction remains independent of the system's implementation.

## Cell indexing

Many grid systems define their own identifier scheme; some define several; some define
none. Computing an identifier, or going from an identifier back to an index in an array,
can be expensive.

So there are two ways to name a cell, and both are always available. A bare `Int` is
always an **index** in `1:ncells(grid)` — a local index into that collection's own
storage. An `AbstractCellIndex` is a typed identity that records its
own level. There is no wrapper type for indices — the two are separated by method
signature alone.

A cell index must be an immutable isbits struct with a total `Base.isless` giving canonical
order, plus `level` and `rawid`. Equality and hashing come from Julia's structural equality
on isbits values, so immutability is part of the contract.

The single gate between the two worlds is internal: a cell of the wrong level, or of a
scheme the system does not claim, resolves to `nothing` rather than throwing. That is what
lets `globalindex` be asked safely about any cell from anywhere.

`reindex` converts a cell between alternate schemes of one system. Only HEALPix has more
than one — nested and ring — so every other system gets identity-or-throw.

## The quad-face family

`AbstractQuadFaceGridSystem` covers systems whose identifier is `face * 4^level + curvecode`
over some number of congruent square lattices. HEALPix, S2 and ISEA4R belong to it.

A member declares three things — how many base faces it has, and two names used in error
messages — and inherits parent, children, ancestor, descendant ranges, descendants, cell
counts, identifier and index arithmetic, sorted-subtree status, the Morton codec, and the
border and interior walkers. All of it is bit arithmetic on the identifier. In practice a
new member of this family needs six methods: the three declarations plus `levels`,
`cell_boundary` and `cell_centroid`.

This is the largest piece of reuse in the package. IGeo7, H3 and A5 are not in the family;
each provides its own identity and hierarchy layer.

There is no chart abstraction. The three chart files share a naming convention and two
helper calls, but no abstract type and no shared interface. They are independent
projections that happen to be organised the same way.

## Neighbours

The hinge of the adjacency layer is `one_ring`, which returns a cell's immediate
neighbours. It is not exported, but every system implements it, and `neighbors` and `ring`
are built on top. A grid without a system gets a geometric fallback that finds neighbours by
matching shared vertices against a spatial tree.

Adjacency is not one relation, so `Connectivity` is an explicit argument rather than a
property of the grid. `Vertex()` and `Edge()` can differ sharply: a system may have many
more corner-touching neighbours than edge-sharing neighbours.

`Winding` declares whether a system's ring order is a rotation that can be carried outward
to shells at distance two and beyond, or whether those shells must be re-sorted by measured
azimuth. Declaring a winding is a speed decision, not a correctness one — the geometric sort
is always available and always correct.

`maxneighbors` is an upper bound, not a count. Declaring one lets the neighbourhood machinery
use fixed-capacity storage and avoid a heap allocation when the bound is known.

## Regions

A region is a subset of one complete level, or a complete level itself. `PartialGrid`, any
`AbstractCellVector` — computed or read from a store — and a complete grid are all regions;
`subtree` is the verb that turns a subtree into one.

On a region, `neighbors` and `ring` answer the complete level's question clipped to
membership. Distance is measured in the system, not in the subset, so a hole does not
lengthen a path around itself.

Four verbs answer about the region as a whole. `halo` walks what is outside it, including
cells punched out of its middle. `border` and `interior` split what is inside. `adjacency`
tables every one-ring at once.

The first three are lazy iterators, serial, and use memory proportional to tree depth rather
than to the answer. `adjacency` is materialised and threaded, and it is the one
that keeps its halo. Its CSR rows support three layouts: out-of-region members can be dropped,
marked in place with a zero, or addressed into a `[region; halo]` buffer. The two
full-width shapes preserve slot index against the canonical ring, so a slot number is a
direction; the clipped shape preserves order only.

`halo` selects among several engines: a width-one band walk for the quad-face systems, a
seeded boundary automaton for the hexagonal ones, and a generic outside-first hierarchy walk
that every specialisation falls back to and that a newly registered system inherits.

## Cell vectors

`AbstractCellVector` is the region contract in vector form, and what the four verbs, the
neighbourhood sweeps, regridding and plotting are all written against. Five methods are
required of it: `system`, `level`, `size`, `getindex`, and `localindex`. Two backings
ship, and they differ in where element `k` comes from, not in what it means.

`CellVector` is the computed one, and the container the rest normalise to. It stores runs of
consecutive *indices* in the complete level grid rather than lists of identifiers, in one
of two forms chosen by whichever is smaller: a run-length form, or explicit indices.

A complete level is a single run. A subtree of a system with sorted subtrees is also a
single run, computed without touching a cell. Set operations walk intervals rather than
cells. Indexing is a binary search over runs and allocates nothing.

`ChunkedCellVector` is the stored one — see [Store IO](tutorials/store_io.md). `region` is the verb that turns any
cell vector into the compressed one; it is the identity on `CellVector` and a memoised
conversion on the stored twin, so the region machinery has exactly one implementation of
each walk.

`PartialGrid` is the grid-shaped sibling: a sorted identifier vector plus the complete level
grid it forwards geometry to. It is a grid, so data can be laid out against its local
indices, and it is a region, so the four verbs accept it.

## Mixed-level sets

`MultiOrderCellSet` holds cells at different levels — the result of an adaptive query that
stopped refining where a coarse cell already sufficed.

It needs its own type for two reasons. There is no single index space across levels, so
it stores sort keys instead of indices. And it deliberately has no `halo`, `border`,
`interior` or `adjacency`, because those questions have no meaning without one level to ask
them at; `member_neighbors` is the cross-level substitute.

It also carries a per-cell flag recording whether containment was *proven* or merely not
asked. The flag is asymmetric because exact containment is more expensive than intersection;
the traversal skips it wherever the answer would not change the result.

## The cube layer

`AbstractCellLookup` makes a cell vector usable as a DimensionalData lookup, so a DGGS axis
behaves like any other dimension. `Base.parent` returns the vector, and every cell verb a
cube supports is defined once against that one method. Code meaning "the cell dimension of
this cube" dispatches on this type and so accepts a cube from a store as readily as one
built in memory.

`CellLookup` is the lookup over a `CellVector`, `ChunkedCellLookup` the one over a stored
axis. Slicing and selection work unmodified, and — the point of the exercise — a slice stays
compressed whenever the result is still ascending. What a subset *is* stays each lookup's
own: a computed window set stays a window set, and a stored axis stops being stored.

`Cells` is the dimension. `Covering` is a selector that takes a geometry and returns the
indices covering it.

## Store IO

Reading a DGGS store involves three concerns, each represented by an abstraction.

`CellEncoding` describes how cell identifiers are laid out on disk — densely, as ranges, or
implicitly — together with the rank and select arithmetic needed to resolve a selector
without scanning the identifier array.

`DGGSConvention` describes how to recognise and interpret a particular community's metadata.
Detection runs every registered convention, sorts
what they return so an explicit declaration beats a fingerprint, and merges the results
field by field. Two conventions disagreeing about the same field is an error rather than a
silent precedence choice.

`StoreDescription` is the merged, format-neutral result. None of these three knows anything
about Zarr; the Zarr dependency is weak and gates only the actual reading and writing.

`ChunkedCellVector` is the axis those three produce, and it is an `AbstractCellVector`, so a
cube from `dggread` is a region like any other. `region` converts it to the compressed form,
and what that costs is the encoding's rather than the axis's length: a stored interval is a
run of consecutive *ranks*, and a rank plus one is an index in the complete level grid, so
the ranges and implicit encodings convert by arithmetic over what they already hold and read
nothing. Only a dense axis has no such structure to borrow; its ids are read once, in the
order that touches each chunk once, and the result is kept.

## Following the chunk lines

A cell-at-a-time pass over a lazy cube can decode the same storage chunk repeatedly when a
ring crosses chunk boundaries. The chunk traversal follows the store's existing chunk grid
and exposes the work as a plan and a runner.

[`chunkplan`](@ref) returns a `MapChunkPlan`: per chunk, the axis indices it owns and the axis
indices outside it that its cells' rings reach. Boundaries come from the data's own chunk
grid, so an irregularly chunked store — one chunk per ancestor subtree, say — is planned on
its real boundaries. Building a plan walks each chunk's boundary through `halo`, which is
CPU and no IO.

The plan is a value because three decisions are made on it: the order to read in, how to cut
it up — `split(plan, n)` is what parallelises a sweep, so there is no `threaded` keyword —
and what it will cost before anything is read.

`foreachchunk` runs it, handing the callback a `ChunkCube`: the chunk's cells *and its halo*,
in memory, as an ordinary cube over a `CellLookup`. Every verb in the package then works on
it unmodified. The owned cells are contiguous within it, because every halo cell is outside
the chunk's own run, so `localindices` is a range and `globalindices` says where those results
belong in the full axis.

Because the halo holds every axis neighbour of every owned cell, a stencil reaching no
further than the plan's halo width computes on a chunk exactly what it computes on the whole
axis. `mapneighbors!` is the streaming form built on that, and `mapneighbors` with
`pass = Values()` takes the same route by itself whenever the data is chunked. `Neighbors()`
deliberately does not: its callback closes over the original array and indexes it by the
handles it is given, so a sweep over blocks would hand it block indices to read from the
whole cube.

## Regridding

Regridding lives in a separate package, `GlobalRegridding`, which owns the verbs and does
not depend on this one. The main package implements its contract for DGGS grids by adding
methods to imported bindings, so a session that loads both sees one function per name.

Its central abstraction is `RegridSpace`: anything with cells, chunks, and a manifold.
A raster is one; a DGGS grid is another. Six methods are required of a space.

What a DGGS space contributes is the idea that a chunk is a non-empty ancestor subtree,
chosen at whichever level lands closest to a target cell count. Chunk discovery is then a
descent of the hierarchy the grid already has, rather than a second spatial index built
alongside it.

Method, missing-data policy and weight storage are three orthogonal axes, each spelled as a
singleton tag: `Conservative`, `NearestCell` or `BilinearPoint`; `Weighted` or `Extensive`;
`PerChunk` or `Spilled`.

## Conformance testing

The executable form of these contracts is
`lib/DiscreteGlobalGridsConformanceTesting`, whose two entry points assert the laws a grid
and a system must satisfy: index round-trips, boundary rings that are closed implicitly and
wound counter-clockwise seen from outside, centroids strictly inside their own cells,
neighbour symmetry and two-hop closure, ring and disc consistency, and the covering law —
that a descendant's boundary lies inside every ancestor's node extent, at every depth.

A system opts out by capability rather than by name: the suite checks whether a method is
actually specialised and skips the corresponding laws if it is not. A new grid or system is
expected to pass both suites.
