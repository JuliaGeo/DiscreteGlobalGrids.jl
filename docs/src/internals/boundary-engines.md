# Boundary traversal engines

```@meta
CurrentModule = DiscreteGlobalGrids
```

This reference is for contributors implementing or tuning boundary traversal.
For the public operations and their index conventions, see
[Region boundaries](../api/boundaries.md).

## Boundary iterators

```@docs
SubtreeHaloIterator
DiscreteGlobalGrids.Engine.SubsetHaloIterator
DiscreteGlobalGrids.RegionSide
EdgeCellIterator
InnerCellIterator
```

## Iterator indices and sizing

```@docs
halo_indices
DiscreteGlobalGrids.Engine.HaloIndexIterator
DiscreteGlobalGrids.sizehint
```

## Engine selection

A boundary operation selects an engine to enumerate candidates and test their
adjacency. A grid system can provide a specialised `halo_engine` method to use
its topology efficiently. These interfaces are internal and may change.

Every engine returns exact results. A candidate band can be larger than the
halo; adjacency tests filter it before the iterator yields cells.

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

The traversal chooses candidate cells; an adjacency provider tests whether they
touch the subject. The geometry provider gives an independent check of results
from the topology-based provider.

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
Pages = ["internals/boundary-engines.md"]
```
