# Grid-family implementation decisions

This is the cross-family answer to the implementation questions in
`2026-08-13-implementation-tasks.md`. The supporting derivations, source audit,
and edge cases are in [the ISEA note](isea-family.md), [the rHEALPix/AusPIX
note](rhealpix-auspix.md), and [the IVEA/RTEA note](ivea-rtea.md).

## Decisions

| Registered system | Kernel to implement | Canonical package identity and hierarchy | External oracle |
|---|---|---|---|
| ISEA3H | Existing Snyder ISEA face map plus aperture-3 triangular-lattice Voronoi cells | DGGRID-compatible `Z3Cell` prefix tree in the package's existing 11.25-degree spherical frame | DGGRID 9.0b black-box geometry, neighbor, indexing-parent, and Z3 crosswalk vectors; a smaller, separate 11.20-degree orientation corpus; point-location seam probes remain a fixture addition |
| ISEA4H | Existing Snyder map plus aperture-4 lattice Voronoi cells | Root-major, base-4 `LevelIndex`; DGGRID ZORDER is an alternate codec | DGGRID 9.0b black-box vectors |
| ISEA4T | Existing Snyder map plus recursive four-way face-triangle subdivision | Zero-based face plus base-4 path in `LevelIndex` | DGGRID 9.0b geometry/centre vectors; face/path crosswalk must be fixed geometrically |
| rHEALPix | HEALPix equal-area projection, polar-triangle rearrangement, and a `3 x 3` square refinement | Published SUID: roots `N O P Q R S`, then row-major digits `0...8`; prefix parentage | MIT option of rHEALPixDGGS 0.6.0, independently projection-checked against MIT PROJ |
| AusPIX | No new geometric kernel: WGS84/authalic wrapper around rHEALPix with Greenwich meridian, polar squares `(0,0)`, and `N_side=3` | The same rHEALPix SUID; no evidence supports the old separate-ordinal claim | The rHEALPix corpus includes a sealed AusPIX/WGS84 profile |
| IVEA variants | Great-circle slice-and-dice equal-area projection with the icosahedron vertex as radial vertex; shared 5-by-6 ten-rhombus atlas | OGC/DGGAL ZIRS for IVEA9R/3H/7H and corresponding DGGAL IDs for 4R; `_Z7` is an alternate encoding | Generated BSD-3-Clause DGGAL 0.0.6 level-zero reconnaissance; a pinned post-fix source build is still required for implementation-gating refinement vectors |
| RTEA variants | The DGGAL RT(S)EA meaning: the same slice-and-dice construction with an icosahedron-edge midpoint/rhombic-triacontahedron face centre as radial vertex | DGGAL ZIRS; `_Z7` is an alternate encoding | Generated BSD-3-Clause DGGAL 0.0.6 level-zero reconnaissance plus published invariants; deeper pinned-source vectors remain an acquisition gate |

The newer Wang et al. projection also called RTEA is not the registered DGGAL
family and must not silently replace it.

## Corrections to the registry assumptions

1. ISEA hex refinements and IVEA/RTEA `3H`/`7H` refinements are central-place
   grids, not nested polygon partitions. Their full geometric parent relation
   is multi-valued. A singular `parent` is therefore a chosen primary tree,
   not the complete covering relationship.
2. DGGAL RTEA rhombic grids have ten hierarchy roots in the common atlas, not
   thirty. Thirty is the number of rhombic-triacontahedron faces. Hex variants
   have twelve level-zero pentagons.
3. AusPIX is a profile of rHEALPix and uses its SUID. It is not a second cell
   ordinal scheme.
4. OGC ISEA3H ZIRS and DGGRID Z3 are different indexing systems and use
   different standard orientations. OGC ZIRS is not a drop-in serialization
   of the package's chosen Z3 tree.
5. A finite corner ring does not exactly represent a straight chart edge after
   inverse ISEA, rHEALPix, IVEA, or RTEA projection. The initial API should
   explicitly call it a corner ring and offer densified edge samples for
   rendering, intersection, and extent validation.

## Shared implementation structure

The implementations divide cleanly into three reusable layers:

1. **Projection charts.** Reuse the existing Snyder ISEA chart for ISEA3H/4H/4T.
   Add one generic slice-and-dice fundamental-triangle kernel parameterized by
   radial vertex for IVEA and RTEA. Add one HEALPix plus polar-rearrangement
   chart for rHEALPix/AusPIX.
2. **Planar refinements.** Implement aperture-3 and aperture-4 triangular
   lattices, the face-triangle quadtree, rHEALPix's base-3 square subdivision,
   and DGGAL-compatible `4R/9R/3H/7H` atlas grids.
3. **Topology and codecs.** Keep geometric covering parents distinct from the
   package's singular primary parent. Implement Z3/SUID/ZIRS codecs outside the
   projection kernels and derive public CCW neighbor cycles from chart topology,
   not from an oracle library's iteration order.

## Implementation order

1. Implement and exhaustively test ISEA4T; it has exact chart containment and a
   conventional four-child tree.
2. Implement rHEALPix and expose AusPIX as a profile. The mathematics and
   hierarchy are completely pinned, and the corpus covers seams and all polar
   placements.
3. Implement ISEA4H, then ISEA3H. Reuse the lattice/seam machinery, but test
   primary prefix parentage separately from geometric covering relationships.
4. Implement the generic slice-and-dice projection and validate one fundamental
   triangle before adding the 5-by-6 atlas. Then add IVEA and RTEA variants in
   the order `4R`, `9R`, `3H`, `7H`; the hex systems require the most care around
   multi-parent hierarchy semantics.

## Acceptance gates

For every system, implementation completion requires: ID/position bijection;
analytic global counts and area; centre, corner-ring, point-location and seam
agreement with the sealed oracle; reciprocal neighbor sets remapped to the
documented CCW cycle; explicit tie ownership; parent/child inverse tests for the
chosen primary relation; and independently sampled validation that every
reported extent contains all descendant geometry. Oracle packages remain test
fixture producers only and are never runtime dependencies.
