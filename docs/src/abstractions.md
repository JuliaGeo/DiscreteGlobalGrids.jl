# Abstractions

This file attempts to go through most of the abstractions in this package, and explain what they are and why they are needed.

## Grids and systems

The idea is to separate the concept of a "grid", a collection of cells, from a "grid system", describing a way to tessellate the sphere.
For example, a grid of H3 at level 12 is a grid; H3 is itself a grid system.

But a grid need not be backed by a hierarchical system.  Consider the tripolar grid of ocean simulation, or a cubed sphere grid, or other
similarly exotic grids.  These are innately non hierarchical, and some of them might even require neighbors to be derived geometrically.

In DiscreteGlobalGrids, we want to support any case you can throw at us, _and_ have it be performant.  But when you have information about e.g.
hierarchy, to derive a tree from, then we need to hook into that too.

Therefore we have these two abstractions.  This also allows partial grids to take advantage of accelerations - they are not full systems, nor fully
instantiated levels of a system, but still backed by some system (if they are).

### AbstractGrid interface

### AbstractHierarchicalGridSystem interface

## Cell indexing

Many DGGSes define their own indexing schemes.  Some may not even have an indexing scheme at all.
And in some cases, computing the cell index or going from a local index in the array to the cell index may be expensive.

There's a difference between "cell index" (the semantic identity) and "local index" (an index into an array).  Some systems may even
have multiple - see HEALPix with its nested, ring, morton and zuniq orderings.

We made the choice that a local index in the array is always a supported way to address
a cell.  