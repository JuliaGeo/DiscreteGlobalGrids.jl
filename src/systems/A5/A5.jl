# ---------------------------------------------------------------------------
# T10 — A5.
#
# felixpalmer's A5 grid on the redesigned interface. The projection, the id
# codec, the hierarchy arithmetic, the adjacency and the boundary all come from
# `native.jl`, a verbatim carry-over of the pure-Julia port of upstream a5, so
# what this module contributes is wiring, a dense order, and a rotational
# neighbour order the upstream does not have (it sorts).
#
# A5 is the package's first system with `has_sorted_subtrees == false`: its
# grids treeify to a *selection-mode* cursor and its multi-order sets take the
# `(level, ordinal)` fallback sort. That is deliberate — see the trait's
# docstring in `system.jl` — and `test/systems/A5/` drives treeify, query and
# coverage through those paths on purpose.
# ---------------------------------------------------------------------------

"""
    DiscreteGlobalGrids.A5

The [A5](https://a5geo.org) discrete global grid system: an equal-area
pentagonal tiling of a dodecahedron, resolutions `0:29`.

  - [`A5System`](@ref) — the system singleton.
  - [`A5Cell`](@ref) — the canonical id, a `UInt64` in upstream a5's own
    serialization, with the resolution in-band.
  - [`A5Grid`](@ref) — one complete resolution, from `levelgrid(A5System(), l)`.
  - [`A5Native`](@ref) — the ported upstream arithmetic, if you want it directly.

# The hierarchy has three regimes

A5 is **not** a fixed-radix system, and nothing here may be derived from an
aperture:

| level | cells | refinement |
|---|---|---|
| 0 | 12 | the dodecahedron's faces |
| 1 | 60 | each face cut into 5 triangular quintants |
| ≥ 2 | `60·4^(l-1)` | 4 Hilbert children per cell |

# The canonical order

Cells at one level are numbered **quintant by quintant**, and within a quintant
by the **Hilbert state `S`**:

    position = quintant · 4^(level-1) + S + 1,     quintant = 5·origin + segment

which is exactly ascending order of the raw `UInt64` at a fixed level (see
[`A5Cell`](@ref)). Level 0 is the twelve faces in origin order.

[`has_sorted_subtrees`](@ref) is nonetheless **`false`** here — the conservative
default, kept because the two-sided `descendant_range` contract has not been
verified across the res-0 → res-1 fan-out, whose quintant numbering is a
*rotation* of a5's own segment walk. See that trait's docstring in `system.jl`
for what it costs and what to do instead.

# Fast paths over the generic fallbacks

| operation | how |
|---|---|
| [`cellat`](@ref) | a5's own `lonlat_to_cell` point location |
| [`cellindex`](@ref) / [`cellposition`](@ref) | the closed-form ordinal above |
| [`neighbors`](@ref) / [`ring`](@ref) | a5's adjacency walk, wound counter-clockwise |
| [`ancestor`](@ref) / [`descendants`](@ref) | `cell_to_parent` / `cell_to_children` across any level gap |

[`node_extent`](@ref) is the generic default — the cell's cap inflated by
[`cap_inflation`](@ref) — with the factor raised to `1.75`; see that docstring
for the measurement. `subtree_border` keeps the generic fallback: an A5 subtree
is four Hilbert children that cover the parent's area but not its footprint, so
there is no digit predicate to read a rim off.
"""
module A5

# Exactly the generics this module defines methods on, plus the types those
# methods dispatch on. `node_extent` is deliberately absent: A5 takes the
# generic default (an inflated cell cap), and declaring `cap_inflation` is the
# whole of its say in the matter. So is `descendant_range`: A5 does not declare
# `has_sorted_subtrees`, so it owes no method.
import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, Connectivity, Vertex, Edge,
    ncells, cellindex, cell_boundary, cell_centroid,
    cellposition, rawid,
    cellat, neighbors, ring, system, level,
    cellindextype, levels, levelgrid, rootcells, children,
    cap_inflation, max_neighbors, has_sorted_subtrees,
    ancestor, descendants

import GeometryOps as GO
import SmallCollections
using SmallCollections: SmallVector

# The unit-sphere vocabulary, spelled the same way the fallback substrate
# spells it.
const USPoint = GO.UnitSphericalPoint{Float64}

# The ported upstream arithmetic first: everything below is a wiring of it.
include("native.jl")

include("cell.jl")
include("system.jl")
include("geometry.jl")
include("neighbors.jl")

end # module A5
