# System capabilities and traversal costs

```@meta
CurrentModule = DiscreteGlobalGrids
```

For choosing a system, start with [Choosing a grid](../tutorials/choosing_a_grid.md).
This page records system traits, boundary algorithms, and interoperability limits for implementors.

## Registered global systems

The size column is [`cellsize`](@ref) in km at levels 0, 4 and 8 — the side of a
square with the median cell's area, which is the only size a system whose cells
differ in shape and area has. Ask [`levelfor`](@ref) for the level a given
resolution wants rather than reading a level off this table.

| system | cells at level `l` | ≈ size (km) at `l` = 0 / 4 / 8 | cell shape | equal-area |
|---|---|---|---|---|
| `IGeo7System` | `10·7^l + 2` | 6520 / 146 / 3.0 | hexagons + 12 pentagons | approximately; pentagons differ; see `IGeo7.equal_area_steradians` |
| `H3System` | `120·7^l + 2` | 2026 / 42.6 / 0.87 | hexagons + 12 pentagons | no (libh3's gnomonic faces) |
| `HEALPixSystem` | `12·4^l` | 6520 / 408 / 25.5 | curvilinear diamonds | yes, exactly `4π/(12·4^l)` |
| `A5System` | `12`, `60`, then `60·4^(l-1)` | 6524 / 365 / 22.8 | pentagons (Cairo-style) | on its ellipsoid; unit-sphere `cell_area` varies by about 1% |
| `S2System` | `6·4^l` | 9220 / 584 / 36.1 | geodesic quadrilaterals | no; ~2.08× within-level spread |
| `ISEA4RSystem` | `10·4^l` | 7142 / 446 / 27.9 | rhombi on ten diamonds | yes, exactly `4π/(10·4^l)` |

Important cross-system traits:

  - **Neighbour degree varies, and [`Vertex`](@ref)/[`Edge`](@ref) need not
    coincide.** [`maxneighbors`](@ref) is an upper bound. IGeo7/H3 have degree
    6 except for 12 pentagons of degree 5, with identical connectivities. A5's
    4-valent corners give
    `maxneighbors(A5System(), Vertex()) == 11` against
    `maxneighbors(A5System(), Edge()) == 5`; at resolution 1, some cells have
    11 vertex-neighbours and 3 edge-neighbours. Above level 1, ISEA4R has ten
    degree-9 cells at icosahedral vertices 0 and 11, thirty degree-7 cells, and
    degree 8 elsewhere; level 0 is 6-regular.
  - **[`node_extent`](@ref).** S2, ISEA4R, and HEALPix provide exact uninflated
    subtree caps. IGeo7, H3, and A5 use inflated caps; A5 sets
    [`cap_inflation`](@ref) to `1.75`.
  - **[`has_sorted_subtrees`](@ref).** True except for A5, whose canonical order
    has not established the two-sided [`descendant_range`](@ref DiscreteGlobalGrids.descendant_range) contract.
  - **[`has_congruent_refinement`](@ref).** True for HEALPix, S2 and ISEA4R,
    whose four children tile their parent exactly; false for IGeo7, H3 and A5.
    A `maxcells` [`MultiOrderCellSet`](@ref) descends through meeting cells
    alone where it holds, and through cells that miss as well where it does not.
  - **[`has_direct_location`](@ref).** True for every system here: each names
    the cell containing a point from the point's coordinates, so a
    [`PartialGrid`](@ref) over any of them locates through its complete level
    rather than by searching its own cells.
  - **[`border`](@ref) on a subtree region.** IGeo7, H3, HEALPix, ISEA4R, and
    S2 provide `O(border)` walkers; A5 uses the `O(subtree)` fallback. Each is a
    resumable [`EdgeCellIterator`](@ref) / [`InnerCellIterator`](@ref) in
    `O(depth)` memory, which is what `border` and [`interior`](@ref) read.
  - **[`halo`](@ref) on a subtree region.** The same boundary from outside, as a resumable
    [`SubtreeHaloIterator`](@ref) in `O(depth)` memory. `l == level(c)` is the
    cell's own one-ring on every system, emitted ascending without a sort.
    Deeper, five of the six take a specialization and A5 does not:

      + **HEALPix, S2 and ISEA4R** walk the width-one band around their square
        block, one pruned quadtree descent per face the halo touches, taken in
        face order — which is canonical order, so the merge is concatenation.
        `O(halo + depth)` time. A block nowhere flush with its face edge has a
        closed-form count (`4·side + 4` under [`Vertex`](@ref), `4·side` under
        [`Edge`](@ref)) and therefore a `length`; a block that crosses a seam
        has a conservative rectangle band that every candidate is filtered
        through the native one-ring before yielding, and declares
        `SizeUnknown()`.
      + **IGeo7 and H3** approach the halo from the root's own same-level
        neighbours. One level down the calibration is already the answer; deeper,
        each neighbour's subtree-border automaton is seeded with an exposed-direction
        arc calibrated by observation and walked to the target. Neighbouring
        subtrees occupy disjoint [`descendant_range`](@ref DiscreteGlobalGrids.descendant_range)s, so walking them in
        range order is canonical without a heap or a seen-set, and every
        candidate goes through the native one-ring. Neither hexagonal engine
        declares a `length`: the observed counts (`3^(d+1) + 3` around a
        hexagon, `5(3^d + 1)/2` around a pentagon) are validated by enumeration,
        not derived from the transition recurrence.
      + **A5** has no [`descendant_range`](@ref DiscreteGlobalGrids.descendant_range) to prune a descent by and no
        validated boundary automaton, so it scans the target level in `O(1)`
        memory and `O(ncells)` time. Its aperture and Hilbert-like indexing are
        not evidence of a square fast path, and none is inferred from them.

    The generic outside-first hierarchy walk — canonical-order candidates with
    [`node_extent`](@ref) cap pruning — is still what every specialization's
    guards return to and what a newly registered system inherits.
  - **[`halo`](@ref) on any other region.** The same question about an
    arbitrary subset, and always an iterator. A rooted grid holding a complete
    subtree delegates to
    [`SubtreeHaloIterator`](@ref) and keeps its system's specialization;
    everything else — a hole, a forgotten root, an arbitrary id list — takes an
    outside-first walk against membership, pruned by the subset's own index
    spans rather than by geometry: a block the subset holds entire is retired by
    one lookup, and a block no NEIGHBOUR of which it touches is retired by the
    coarse-containment law. The walk follows the subset's boundary, so its cost
    is the halo's and not the subset's. A cell punched
    out of the middle of a subset is outside it and touches it, so it joins the
    halo. A5 is again the exception to the delegation: without
    [`has_sorted_subtrees`](@ref) there is no way to recognise a held subtree, so
    even its rooted complete grid takes the subset walk — to the same answer.
  - **[`adjacency`](@ref).** Every one-ring of a region at once, CSR and
    counter-clockwise, with the out-of-region members dropped (`halo = 0`),
    addressed into a `[region; halo]` buffer (`halo = 1`), or marked with `0`
    in place (`halo = :mark`). The two complete-width shapes preserve slot
    indices against the canonical `one_ring`, so a direction code is a property
    of the cell; the clipped shape preserves order only. Rows exist for
    in-region indices alone, so `halo` above 1 throws and points at
    [`grow`](@ref).
  - **Cross-level adjacency ([`member_neighbors`](@ref)).** Boundary sharing in
    the geometric sense on HEALPix, S2 and ISEA4R, whose four children tile
    their parent exactly; the hierarchy's own relation on IGEO7, H3 and A5,
    where they do not and a member's footprint is not its descendants' union.

## Interoperability caveats

  - ISEA4R's face pairing, diamond numbering, and `(x, y)` orientations are
    package-defined and anchored on vertices `(0, 11)`. Compatibility with
    external ISEA4R identifiers, including DGGAL, is not claimed; establish a
    fixture-derived permutation before interchange.
  - S2's native 64-bit `s2_cellid` is not an available [`reindex`](@ref) scheme.
    The scaffold ordinal is the canonical identifier.

Use [`levels`](@ref) and [`levelgrid`](@ref) to construct queryable grids. Each
system module documents its identifier codec and optimized operations.

## Multi-order budget observations

Budget-mode coverage is exact on congruent systems. Tests also measure misses
on noncongruent systems; those measurements are not guarantees for arbitrary targets.
The initial implementation notes recorded these fractions of a state-sized target missed:

| Measurement | IGeo7 | Authalic IGeo7 | H3 | A5 |
| --- | --- | --- | --- | --- |
| Target outside emitted polygons | 2% | 3% | 15% | 30% |
| Intersecting reference-level cells not represented | 1% | 2% | 2% | 18% |

Consult the multi-order budget tests for their actual targets and thresholds.
Fixed-level queries can continue overhang descent to a fixed depth. Budget
queries cannot rely on that same stopping rule, so their geometric and leaf
coverage must be assessed separately. Expanding a budget set reports its
represented descendants; it does not repair cells the budget traversal missed.
