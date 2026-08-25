# Region boundaries

```@meta
CurrentModule = DiscreteGlobalGrids
```

A boundary has two faces, and this package gives each one a verb. All three take
a *region* — a [`PartialGrid`](@ref), a [`CellVector`](@ref), a
[`CellLookup`](@ref) or a complete grid — and [`subtree`](@ref) is how a rooted
subtree becomes one.

The **inside** is [`border`](@ref) — the cells of the region that have a
neighbour outside it — and its complement [`interior`](@ref). The **outside** is
[`halo`](@ref): the cells at the same level that are *not* in the region but
touch it. Neither derives the other, and a chunked computation needs both — the
border is what a chunk sends, the halo is what it fetches.

Because the region is arbitrary, a cell punched out of the middle of a subset is
part of that subset's halo. All three verbs return a lazy, resumable,
`O(depth)`-memory iterator; materialising is the caller's decision, because a
halo can be far larger than the border it wraps.

Three properties are contracts rather than accidents, and all three are worth
knowing before writing any code against these verbs:

  - **Indices by default, ids on request.** Each verb yields indices and
    takes `cells = true` to yield cell ids instead. A halo's indices are
    global indices, into the complete level grid; a border's and an
    interior's are local indices, into the region.
  - **The order is ascending index on the target grid**, strictly increasing,
    on every system and through every engine. A consumer indexing
    index-ordered storage needs no `sort` and no `sortperm`. These are fetch
    lists, and they are the *only* neighbourhood verbs here that are ascending:
    [`neighbors`](@ref), [`ring`](@ref), [`adjacency`](@ref) and
    [`member_neighbors`](@ref) all answer counter-clockwise seen from outside the
    sphere, in ids and in indices alike.
  - **`length` is refused wherever the count is unproved.** Face seams, poles
    and pentagons break perimeter formulas, so most engines declare
    `SizeUnknown()` and define no `length` at all.
    [`sizehint`](@ref DiscreteGlobalGrids.sizehint) is the separate, explicitly
    approximate answer for a caller who only wants to `sizehint!` a vector.

Neither `collect` nor `Set` is special-cased: `collect` fills a vector in the
walk's own order, sized by the hint when there is one, and `Set` gives
membership. Both consume the walk once.

## The region, and its two faces

```@docs
subtree
halo
border
interior
```

## The walks

```@docs
SubtreeHaloIterator
DiscreteGlobalGrids.Engine.SubsetHaloIterator
DiscreteGlobalGrids.RegionSide
EdgeCellIterator
InnerCellIterator
```

## Index space, and sizing

```@docs
halo_indices
DiscreteGlobalGrids.Engine.HaloIndexIterator
DiscreteGlobalGrids.sizehint
```

## Adjacency on a region

The verbs a halo is most often confused with. [`adjacency`](@ref) is the cached
one-ring table over a region — which cells each of its cells touches, in three
row shapes — and [`member_neighbors`](@ref) is adjacency across the levels of a
[`MultiOrderCellSet`](@ref), which has no halo because it has no single level to
answer at.

```@docs
adjacency
AdjacencyTable
halocells
haloindices
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
localindex
globalindex
descendant_range
descendants
has_sorted_subtrees
has_direct_location
node_extent
maxneighbors
LevelIndex
PartialGrid
AbstractCellVector
CellVector
AbstractCellLookup
CellLookup
region
MultiOrderCellSet
MultiOrderCoverage
grow
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
DiscreteGlobalGrids.halo_engine
DiscreteGlobalGrids.Engine.generic_halo_engine
DiscreteGlobalGrids.Engine.RingHaloEngine
DiscreteGlobalGrids.Engine.OutsideWalkEngine
DiscreteGlobalGrids.Engine.ScanHaloEngine
DiscreteGlobalGrids.Engine.subset_halo_engine
DiscreteGlobalGrids.Engine.geometry_halo_engine
```

### The aperture-4 band

HEALPix, S2 and ISEA4R: a subtree is an aligned square block in one face's
lattice, and its halo is the width-one band around it.

```@docs
DiscreteGlobalGrids.Engine.SquareBandEngine
DiscreteGlobalGrids.Engine.square_halo_engine
DiscreteGlobalGrids.Engine.FaceRect
DiscreteGlobalGrids.Engine.SquareBandWalk
DiscreteGlobalGrids.Engine.NoCheck
DiscreteGlobalGrids.Engine.NativeCheck
```

### The aperture-7 directed walk

IGeo7 and H3: a subtree's halo lies under the root's own same-level neighbours,
and is reached by seeding each neighbour's border automaton with the arc that faces
the root.

```@docs
DiscreteGlobalGrids.Engine.hex_halo_engine
DiscreteGlobalGrids.Engine.HexChildHaloEngine
DiscreteGlobalGrids.Engine.HexChildWalk
DiscreteGlobalGrids.Engine.HexArcHaloEngine
DiscreteGlobalGrids.Engine.HexArcWalk
DiscreteGlobalGrids.Engine.HexNeighbour
DiscreteGlobalGrids.Engine._hex_calibrate
DiscreteGlobalGrids.Engine._hex_validate
DiscreteGlobalGrids.Engine._minimal_arc
```

### Adjacency providers

Candidate enumeration and adjacency testing are kept apart: the walk decides
which cells to consider, a provider decides whether one of them touches the
subject. The geometry provider shares no topology with the indexed one, which is
what makes it an oracle rather than a second opinion.

```@docs
DiscreteGlobalGrids.Engine.IndexedNeighbors
DiscreteGlobalGrids.Engine.ForcedGeometry
DiscreteGlobalGrids.Engine.SubsetMembership
```

### The border automata the halo walks borrow

```@docs
DiscreteGlobalGrids.Fallbacks.SquareBorderEngine
DiscreteGlobalGrids.Fallbacks.SquareInteriorEngine
DiscreteGlobalGrids.Fallbacks.ScanBorderEngine
DiscreteGlobalGrids.Fallbacks.MortonCurve
DiscreteGlobalGrids.Fallbacks.quadrant_step
DiscreteGlobalGrids.face_orientation
DiscreteGlobalGrids.Fallbacks.cells_cap
```

## Index

```@index
Pages = ["api/boundaries.md"]
```
