# Collection and traversal contracts

```@meta
CurrentModule = DiscreteGlobalGrids
```

This page preserves storage and pruning details for implementors. Start with
[Selecting cells](../api/selecting-cells.md) for collection choices and conversions.

## CellVector storage and expansion

CellVector(set::MultiOrderCellSet; level = set's reference level)
    CellVector(grid::AbstractGrid)
    CellVector(sys, level, ids::AbstractVector)

An immutable, lazy `AbstractVector` of strictly ascending cell ids from one
level of one hierarchical system. It stores their leaf-grid index windows
instead of the ids.

Semantically `cv` **is** the id vector: `length(cv)` is the number of cells,
`cv[k]` is the `k`th of them, `collect(cv)` is the vector itself. What is
*stored* is sorted, disjoint intervals or a sorted index list in
`levelgrid(system(cv), level(cv))`, so memory is O(number of windows) rather
than O(number of cells). `cv[k]` searches the windows' cumulative lengths and
resolves one `cellindex`; [`localindex`](@ref DiscreteGlobalGrids.localindex) runs the inverse. Nothing is
materialised.

[`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup) provides the DimensionalData wrapper and delegates its
collection operations to this type.

# The three ways in

  - a [`MultiOrderCellSet`](@ref), the compressed coverage, optionally
    re-expanded to a deeper `level` than the set's own reference level;
  - an [`AbstractGrid`](@ref) — `levelgrid(sys, l)`, a whole level, which is one
    window; or a [`PartialGrid`](@ref), which is that subset's indices (one
    window when the subset is a subtree, the explicit list when it is
    scattered);
  - `sys, level, ids` — an explicit **strictly ascending** id vector, validated
    and run-compressed.

`CellVector(cv)` returns `cv` unchanged.

# The verbs

```julia
cv[k]                          # the kth cell id
localindex(cv, c)              # its inverse: Int, or `nothing`
c in cv                        # membership, O(log #windows)
localindex(cv, lon, lat)       # the index of the cell a point falls in
cellat(cv, lon, lat)           # that cell's id instead
covering(cv, polygon)          # the sub-vector a region's coverage names
intersect(cv, other)           # two vectors at the same level, O(#windows)
PartialGrid(cv)                # read as a grid, O(1) — the regridding handshake
cellset(cv)                    # what it was built from
```

`Base.summarysize(cv)` is O(number of windows) plus the object referenced by
[`cellset`](@ref). Indexing is O(log(number of windows)) and allocation-free.

Where [`has_sorted_subtrees`](@ref) holds, [`level_ranges`](@ref) constructs the
windows in
O(#entries). Where it does not (A5), a cell's descendants are not one interval
of their level, so the vector is built by explicit selection: `descendants` names the
leaves, they are resolved to indices and sorted, and the result is
run-compressed like any other index list. Every method above is unchanged
as any other index list. This construction visits every leaf. The stored
form ranges from a small set of windows to one index per cell.

For A5, descendants need not lie inside their parent's footprint, so expanding a coverage
names leaves the target does not touch — most visibly inside a hole. A
[`covering`](@ref) subset on A5 is therefore a superset of the cells that meet
the region, by the same margin the refinement itself is; see
[`MultiOrderCoverage`](@ref).

## Subtree extent contract

node_extent(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex) -> GO.UnitSpherical.SphericalCap

The covering region of the entire subtree rooted at `c`.

**Defaulted**: the generic implementation inflates the cell's own bounding cap
by [`cap_inflation(sys)`](@ref cap_inflation), and a system able to compute a
tighter covering cap overrides it. The covering law below holds either way, and
validating it is the system's responsibility.

> **Guaranteed.** `node_extent(sys, c)` contains the *geometry* of every
> descendant of `c`, at every depth — every point of every cell boundary in the
> subtree, all the way down to `maxlevel(sys)`, independent of a caller's
> planned traversal depth.
>
> **Not guaranteed.** It need not contain the descendants' `cell_cap`s or their
> own `node_extent`s, and for real systems it does not. A cap is an inflated or
> sampled bound *around* geometry, not geometry, so a child's cap can stand
> partly outside its parent's extent while every child polygon lies inside it.

Both halves are load-bearing. A cell's own boundary does not cover its
subtree — aperture-7 children extend beyond the parent boundary — so `cell_cap`
is never a substitute for `node_extent`. Equally, the extents down a branch are
not a nested chain of caps: a consumer that compares a parent's `node_extent`
against a child's `node_extent` or `cell_cap`, or that expects radii to shrink
with depth, is testing a property this function does not have and never
promised.

Pruning is sound against queries derived from geometry, which is what tree
descent asks: a query cap disjoint from `node_extent(sys, c)` meets no
descendant cell polygon, so discarding the subtree can drop no result. Every
node is tested against the query, never against its parent. Over-coverage only
reduces pruning; under-coverage can omit valid results. For a convex cap of
angular radius at most 90°, containing all boundary vertices also contains their
great-circle arcs. Non-convex extents must establish containment of the full
geometry.

A merge of the children's caps (`Extents.union` over `cell_cap`) is a bound over
caps — a different and looser object, and one that says nothing about levels
below those children. It is not `node_extent` and does not substitute for it.

Extents are `SphericalCap`s throughout this package, at every node of every
tree, which is what lets one predicate vocabulary serve all of them.
