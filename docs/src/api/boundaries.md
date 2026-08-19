# Subtree and subset boundaries

```@meta
CurrentModule = DiscreteGlobalGrids
```

A boundary has two faces, and this package gives each one a verb.

The **inside** is [`subtree_border`](@ref) — the cells of a subtree that have a
neighbour outside it — and its complement [`subtree_interior`](@ref). The
**outside** is [`subtree_halo`](@ref): the cells at the same level that are *not*
in the subtree but touch it. Neither derives the other, and a chunked
computation needs both — the border is what a chunk sends, the halo is what it
fetches.

[`halo`](@ref) asks the outside question of an arbitrary same-level subset
instead of a subtree, so a cell punched out of the middle of a subset is part of
its halo. All of these are `collect` of a lazy, resumable, `O(depth)`-memory
iterator; materialising is the caller's decision, because a halo can be far
larger than the rim it wraps.

Two properties of the walk are contracts rather than accidents, and both are
worth knowing before writing any code against these verbs:

  - **The order is ascending `cellposition` on the target grid**, strictly
    increasing, on every system and through every engine. A consumer indexing
    position-ordered storage needs no `sort` and no `sortperm`.
    [`halo_positions`](@ref) hands that stream over directly, without
    materialising the ids on the way. These two verbs are fetch lists, and they
    are the *only* neighbourhood verbs here that are ascending: [`neighbors`](@ref),
    [`ring`](@ref), [`halo_table`](@ref) and [`member_neighbors`](@ref) all
    answer counter-clockwise seen from outside the sphere, in ids and in
    positions alike.
  - **`length` is refused wherever the count is unproved.** Face seams, poles
    and pentagons break perimeter formulas, so most engines declare
    `SizeUnknown()` and define no `length` at all. [`halo_sizehint`](@ref) is the
    separate, explicitly approximate answer for a caller who only wants to
    `sizehint!` a vector.

## The outside face

```@docs
subtree_halo
SubtreeHaloIterator
halo
DiscreteGlobalGrids.Fallbacks.SubsetHaloIterator
```

## Position space, and sizing

```@docs
halo_positions
DiscreteGlobalGrids.Fallbacks.HaloPositionIterator
halo_sizehint
```

## The inside face

```@docs
subtree_border
subtree_interior
EdgeCellIterator
InnerCellIterator
```

## Adjacency on a subset

The two verbs a halo is most often confused with. [`halo_table`](@ref) is the
*in-set* stencil — which of a subset's own cells each of its cells touches —
and [`member_neighbors`](@ref) is adjacency across the levels of a
[`MultiOrderCellSet`](@ref), which has no `halo` because it has no single level
to answer at.

```@docs
halo_table
member_neighbors
neighbors
ring
```

## The vocabulary these are stated in

Listed here because the docstrings above cross-reference them, not because this
is where they belong: a general API reference for the package does not exist
yet, and these entries are what make the links on this page resolve.

That reference stops here rather than going one step further on purpose. Every
`@ref` written *by the boundary family* now lands on this page, with one
exception: `Base.IteratorSize`, which is Base's docstring and not this site's to
render. What does *not* resolve is the next ring out — the container and
interface docstrings pulled in below cross-reference the rest of the package,
and following those links to closure means rendering 530 entries, which is the
whole package and then some. Those stay dead until a full reference exists to
catch them.

!!! note "Three names print someone else's docstring first"

    `ncells`, `treeify` and `getcell` are `ConservativeRegridding.Trees`
    bindings this package extends rather than owns. `?ncells` in the REPL prints
    the upstream docstring above this package's; both are there, and the second
    is the one that describes a grid.

```@docs
Connectivity
Vertex
Edge
levelgrid
cellindex
cellposition
descendant_range
descendants
has_sorted_subtrees
node_extent
max_neighbors
LevelIndex
PartialGrid
CellVector
CellLookup
MultiOrderCellSet
MultiOrderCoverage
```

## The engines

Internal, and documented because a halo verb's behaviour *is* its engine's: the
verb chooses one at construction by private dispatch, and the choice is what
decides the cost, the count contract and the emit rule. A system adds a fast
path by adding a `halo_engine` method, so this section is the reference for
someone writing one — not public API, and not stable.

Every engine here is exact. A specialization may enumerate a conservative
candidate band, but it runs the adjacency test on every candidate before
yielding it.

```@docs
DiscreteGlobalGrids.Fallbacks.halo_engine
DiscreteGlobalGrids.Fallbacks.generic_halo_engine
DiscreteGlobalGrids.Fallbacks.RingHaloEngine
DiscreteGlobalGrids.Fallbacks.OutsideWalkEngine
DiscreteGlobalGrids.Fallbacks.ScanHaloEngine
DiscreteGlobalGrids.Fallbacks.subset_halo_engine
DiscreteGlobalGrids.Fallbacks.geometry_halo_engine
```

### The aperture-4 band

HEALPix, S2 and ISEA4R: a subtree is an aligned square block in one face's
lattice, and its halo is the width-one band around it.

```@docs
DiscreteGlobalGrids.Fallbacks.SquareBandEngine
DiscreteGlobalGrids.Fallbacks.square_halo_engine
DiscreteGlobalGrids.Fallbacks.FaceRect
DiscreteGlobalGrids.Fallbacks.SquareBandWalk
DiscreteGlobalGrids.Fallbacks.NoCheck
DiscreteGlobalGrids.Fallbacks.NativeCheck
```

### The aperture-7 directed walk

IGeo7 and H3: a subtree's halo lies under the root's own same-level neighbours,
and is reached by seeding each neighbour's rim automaton with the arc that faces
the root.

```@docs
DiscreteGlobalGrids.Fallbacks.hex_halo_engine
DiscreteGlobalGrids.Fallbacks.HexChildHaloEngine
DiscreteGlobalGrids.Fallbacks.HexChildWalk
DiscreteGlobalGrids.Fallbacks.HexArcHaloEngine
DiscreteGlobalGrids.Fallbacks.HexArcWalk
DiscreteGlobalGrids.Fallbacks.HexNeighbour
DiscreteGlobalGrids.Fallbacks._hex_calibrate
DiscreteGlobalGrids.Fallbacks._hex_validate
DiscreteGlobalGrids.Fallbacks._minimal_arc
```

### Adjacency providers

Candidate enumeration and adjacency testing are kept apart: the walk decides
which cells to consider, a provider decides whether one of them touches the
subject. The geometry provider shares no topology with the indexed one, which is
what makes it an oracle rather than a second opinion.

```@docs
DiscreteGlobalGrids.Fallbacks.IndexedNeighbors
DiscreteGlobalGrids.Fallbacks.ForcedGeometry
DiscreteGlobalGrids.Fallbacks.SubsetMembership
```

### The rim automata the halo walks borrow

```@docs
DiscreteGlobalGrids.Fallbacks.SquareRimEngine
DiscreteGlobalGrids.Fallbacks.SquareInteriorEngine
DiscreteGlobalGrids.Fallbacks.ScanRimEngine
DiscreteGlobalGrids.Fallbacks.MortonCurve
DiscreteGlobalGrids.Fallbacks.quadrant_step
DiscreteGlobalGrids.Fallbacks.face_orientation
DiscreteGlobalGrids.Fallbacks.cells_cap
```

## Index

```@index
Pages = ["api/boundaries.md"]
```
