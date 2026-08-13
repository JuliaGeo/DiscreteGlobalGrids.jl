# ISEA3H, ISEA4H, and ISEA4T: implementation mathematics and oracle plan

Status: research/design note, 2026-08-13. This note is intentionally more
specific than the retired registry stubs. It distinguishes published grid
mathematics, a package indexing choice, and facts observed from a black-box
DGGRID executable.

## Conclusions

1. All three grids can reuse the package's existing `ISEA` icosahedron and
   Snyder equal-area face chart. The missing kernels are planar lattice or
   triangle-subdivision kernels plus global seam canonicalization; the Snyder
   projection does not need to be written again.
2. Use level 0 for the coarsest grid. ISEA3H and ISEA4H have 12 level-0
   pentagons centered on icosahedron vertices. ISEA4T has the 20 icosahedron
   faces as level-0 triangles. Some older ISEA3H publications call the
   32-zone grid “resolution 1” and simply omit the 12-zone dual icosahedron;
   DGGRID and this package include level 0.
3. The hexagonal hierarchies are central-place, not nested polygon
   partitions. A fine zone at a coarse vertex participates in three
   geometric parents. DGGRID's `.chd` output therefore has six children for
   a pentagon and seven for an ordinary hexagon, even when the aperture is 3
   or 4. That relation cannot be used directly for this package's required
   single-valued `parent`/`children` interface.
4. For ISEA3H, use a DGGRID-Z3-compatible prefix spanning tree as the
   package's canonical hierarchy and a packed `Z3Cell` as the canonical ID.
   It gives analytic parent/children, a dense order, and contiguous subtree
   ranges. The OGC ISEA3H ZIRS is valuable but is a different, non-prefix
   identifier and sub-zone ordering; implement it later as an interoperability
   codec/system, not as the first strict-tree kernel.
5. For ISEA4H, use the analogous ZORDER prefix tree and a 0-based
   `LevelIndex` dense ordinal as canonical. Preserve DGGRID's ZORDER digit
   string as an alternate codec. For ISEA4T, use a 20-face, four-child
   triangle path and a 0-based `LevelIndex` ordinal.
6. Keep the exhaustive DGGRID oracle at the repository's existing standard
   ISEA placement (spherical vertex 0 at 11.25 degrees E,
   58.282525588539 degrees N). Add a smaller second orientation corpus at
   11.20 degrees E to pin the OGC orientation. Do not replace the 11.25-degree
   corpus.

## Sources and license screen

Only papers/standards and permissively licensed software were used as
implementation references.

| Source | What it establishes | License/use decision |
|---|---|---|
| Snyder, *An Equal-Area Map Projection for Polyhedral Globes* (1992), [DOI 10.3138/27H7-8K88-4882-1752](https://doi.org/10.3138/27H7-8K88-4882-1752) | Forward/inverse equal-area face map | Paper; primary mathematical source |
| Sahr, White, and Kimerling, *Geodesic Discrete Global Grid Systems* (2003), [DOI 10.1559/152304003100011090](https://doi.org/10.1559/152304003100011090) | Five DGGS design choices; ISEA3H construction; aperture/class terminology | Paper; primary grid source |
| Carr, Kahn, Sahr, and Olsen, *ISEA Discrete Global Grids* (1997), [newsletter PDF](https://mason.gmu.edu/~dcarr/lib/v8n2.pdf) | Level-1 construction, alternating ISEA3H orientation, `1/sqrt(3)` scale, non-nesting | Paper; primary early specification |
| Sahr, *Location coding on icosahedral aperture 3 hexagon discrete global grids* (2008), [DOI 10.1016/j.compenvurbsys.2007.11.005](https://doi.org/10.1016/j.compenvurbsys.2007.11.005) | Q2DI, modified balanced ternary, and aperture-3 path addressing | Paper; primary indexing source |
| OGC API - DGGS Part 1, Annex B.4 (2025), [OGC 21-038r1](https://docs.ogc.org/is/21-038r1/21-038r1.html#_isea3h_dggrs_definition) | Current OGC ISEA3H orientation, text/uint64 ZIRS, scanline sub-zone order | Open standard; normative interoperability source |
| [PROJ ISEA documentation](https://proj.org/en/stable/operations/projections/isea.html), [license](https://proj.org/en/stable/about.html#license) | Independent projection behavior and orientation cross-check | X/MIT. Current individual polyhedral files also carry permissive Apache-2.0/BSD provenance; preserve their file-level notices if code is copied |
| [geogrid](https://github.com/mocnik-science/geogrid), especially `ISEA3H.java` | Independent planar aperture-3 center lattice and nearest-center quantization | MIT, copyright Heidelberg University |
| [DGGAL](https://github.com/ecere/dggal), especially `RI3H.ec`, `I3HSubZones.ec`, and `ISEA3H.ec` | Independent ISEA3H topology, OGC ZIRS, seams, neighbors, and sub-zone order | BSD-3-Clause, copyright Ecere Corporation |

Explicit exclusions:

- DGGRID is AGPL. Its source is not an implementation reference. Its released
  executable is used only as a black-box oracle, which the task explicitly
  permits. The command documentation is used only to configure that oracle.
- `dggridR` and `dggrid4py` wrap DGGRID and are not implementation references.
- OpenEAGGR is LGPL-3.0. It was evaluated and rejected under the requested
  non-GPL rule, even though it contains ISEA3H/ISEA4T code. Published reviews
  also report known ISEA3H hierarchy-traversal problems in that library.

## Common projection and orientation

Let each planar icosahedron face be an equilateral triangle of side

```text
R_EA = sqrt(4*pi / (15*sqrt(3)))
L    = sqrt(3) * R_EA = sqrt(4*pi / (5*sqrt(3))).
```

Its area is `sqrt(3)*L^2/4 = pi/5`, equal to one twentieth of the
unit sphere's area. Snyder's map is area preserving between the spherical
face and this triangle. The local implementation already exposes the needed
operations as `ISEA.snyder_fwd(p) -> (face, complex_position)` and
`ISEA.snyder_inv_xyz(face, position)`. It also owns face assignment, inverse
iteration, and the standard 12-vertex frame. Reuse those functions rather
than importing another projection kernel.

The repository's identity orientation is the historical spherical ISEA/DGGRID
placement:

```text
vertex-0 longitude = 11.25 degrees E
vertex-0 latitude  = atan(golden_ratio) = 58.282525588538995 degrees N
azimuth to the specified adjacent vertex = 0 degrees (due north)
```

`ISEA.Orientation` rotates world coordinates into this fixed grid frame on
input and back to world coordinates on output. Therefore all three systems
should retain the fixed chart and store an `Orientation` in the system object.

### OGC orientation is deliberately different

OGC 21-038r1 Annex B.4 defines ISEA3H with vertex longitude 11.20 degrees E
and vertex latitude 58.397145907431 degrees **geodetic WGS84**. That geodetic
latitude becomes 58.282525588539 degrees on the WGS84 authalic sphere. Thus:

- a unit-sphere oracle uses `(11.20, 58.282525588539, 0)`;
- an ellipsoidal OGC-facing wrapper accepts geodetic latitude and applies the
  package's authalic transform first;
- passing 58.397145907431 directly to DGGRID's spherical grid would create the
  wrong orientation.

Orientation is part of the DGGRS, not merely a serialization choice. An OGC
zone ID and a DGGRID Z3 string cannot be called alternate IDs for the same
cell unless both the indexing conversion and the orientation are explicit.

## Planar grid mathematics

Use a triangular lattice with basis

```text
e1 = (1, 0)
e2 = (1/2, sqrt(3)/2).
```

An aligned hex refinement with integer parameters `(h,k)` has aperture

```text
A = h^2 + h*k + k^2
```

and lattice transform

```text
M(h,k) = [ h  -k
           k h+k ],              det(M) = A.
```

The finer center lattice uses the inverse transform. Aperture 3 is `(h,k) =
(1,1)`: scale by `1/sqrt(3)` and rotate by 30 degrees (the sign depends on
the chosen face basis). Aperture 4 is `(2,0)`: scale by `1/2` with no
rotation. Implementing the matrix, rather than parity-specific floating
constants, makes the class-I/class-II change visible and testable.

The hexagons are planar Voronoi cells of these center lattices. If the nearest
center spacing is `d`, a regular planar hexagon has side/circumradius
`s = d/sqrt(3)`. At an icosahedron vertex five face wedges glue into a
pentagon; everywhere else six wedges form a hexagon. Project the planar
Voronoi vertices through `snyder_inv_xyz` and canonicalize duplicates on face
seams.

### ISEA3H

The equivalent Goldberg frequency parameters at level `r` are

```text
r = 2j:     (m,n) = (3^j, 0)       class I
r = 2j + 1: (m,n) = (3^j, 3^j)     class II
T = m^2 + m*n + n^2 = 3^r.
```

The sequence begins `GP(1,0)` (12-zone dodecahedral dual), `GP(1,1)`
(truncated icosahedron, 32 zones), `GP(3,0)` (92 zones), and `GP(3,3)`
(272 zones). Carr et al. describe the same rule geometrically: odd published
resolutions have a central hexagon with a base parallel to the face base;
even published resolutions have a vertex pointing toward the base; each step
scales an edge by `1/sqrt(3)`.

For a face-local implementation matching geogrid's convention at DGGRID level
`r >= 1`, let `l = L / sqrt(3)^(r-1)`. In the parity-normalized face frame,
center coordinates can be written

```text
x = nx * l/2
y = (ny + (isodd(nx) ? 1/2 : 0)) * l/sqrt(3).
```

Swap the face-local x/y roles on the alternate refinement parity. The planar
hex vertex offsets in that normalized frame are

```text
(+l/3, 0), (+l/6, +l/(2sqrt(3))), (-l/6, +l/(2sqrt(3))),
(-l/3, 0), (-l/6, -l/(2sqrt(3))), (+l/6, -l/(2sqrt(3))).
```

The matrix formulation should be the production implementation; these
coordinates are useful as a direct independent test.

### ISEA4H

ISEA4H is class I at every level:

```text
(m,n) = (2^r, 0),   T = 4^r.
```

On each face, use all barycentric triangular-mesh vertices with denominator
`2^r`, glued across common icosahedron edges. The nearest-center spacing is
`d = L/2^r`, and each interior planar Voronoi hexagon has side
`L/(sqrt(3)*2^r)`. Level 1 therefore places centers at the 12 icosahedron
vertices and 30 edge midpoints, giving 42 zones.

### ISEA4T

At level `r`, divide every face edge into `2^r` parts. Each root face contains
`4^r` equal-area triangles. For a CCW parent `(a,b,c)`, define

```text
ab = (a+b)/2, bc = (b+c)/2, ca = (c+a)/2

digit 0: (a,  ab, ca)
digit 1: (ab, b,  bc)
digit 2: (ca, bc, c)
digit 3: (ab, bc, ca)       # central, inverted relative to the parent outline
```

All tuples above are CCW in a conventional upright parent. Recursing this
rule gives unique parentage, exact containment in projection space, and four
children per cell. Map the three corners and the planar barycenter through
the existing inverse Snyder chart for the boundary vertices and region point.

## Counts and areas

For either hex grid with aperture `A` (`A=3` or `4`):

```text
ncells(r)       = 10*A^r + 2
pentagons(r)    = 12
hexagons(r)     = 10*(A^r - 1)
hex area        = 4*pi / (10*A^r) steradians
pentagon area   = (5/6) * hex area.
```

The “hex area” at level 0 is the hypothetical regular-hex unit; each of the
12 actual pentagons has area `4*pi/12`.

For ISEA4T:

```text
ncells(r) = 20*4^r
cell area = 4*pi/ncells(r) = pi/(5*4^r) steradians.
```

Useful exhaustive-oracle sizes are:

| level | ISEA3H | ISEA4H | ISEA4T |
|---:|---:|---:|---:|
| 0 | 12 | 12 | 20 |
| 1 | 32 | 42 | 80 |
| 2 | 92 | 162 | 320 |
| 3 | 272 | 642 | 1,280 |
| 4 | 812 | 2,562 | 5,120 |
| 5 | 2,432 | 10,242 | 20,480 |
| 6 | 7,292 | 40,962 | 81,920 |
| 7 | 21,872 | 163,842 | 327,680 |
| 8 | 65,612 | 655,362 | 1,310,720 |
| 9 | 196,832 | 2,621,442 | 5,242,880 |

## Canonical identifiers and hierarchy semantics

### Recommended ISEA3H canonical ID: `Z3Cell`

DGGRID's Z3 digit string consists of a two-decimal-digit root `00` through
`11`, followed by one ternary digit per level. Roots `00` through `09` each
have three prefix children; the two polar roots `10` and `11` have only the
all-zero chain. This produces exactly `10*3^r+2` IDs.

Use an isbits `Z3Cell(::UInt64)` patterned after the package's `Z7Cell`:

```text
bits 63:60       root 0:11
then 30 x 2-bit digits, most significant digit first
active digit     0, 1, or 2
unused digit     3
```

The first padding digit gives the level; no padding means level 30. Structural
validity additionally requires roots 10 and 11 to contain only active zero
digits. This bounds `levels(sys)` to `0:30`, well inside `Int64` cell counts.

For a root `q < 10` and base-3 path value `p` at level `r`, the 1-based dense
position is

```text
q*3^r + p + 1.
```

The two polar positions are `10*3^r+1` and `10*3^r+2`. Parent drops the final
digit. Children append `0:2`, except that a polar cell appends only `0`.
At a fixed target level, all descendants of a prefix form one contiguous
dense interval, so `has_sorted_subtrees = true` and `descendant_range` is
closed-form.

This `children` is the chosen **primary prefix tree**, not DGGRID's geometric
`.chd` relation. Every prefix child is present among the parent's geometric
children, but a geometric vertex child is present in multiple `.chd` rows.
Tests must keep the two concepts separate.

### OGC ISEA3H ZIRS comparison

OGC Annex B.4 instead identifies a zone via the dual ISEA9R rhombus grid.
Its uint64 layout is:

```text
bits  1:0   subhex: 0 even level; 1 center odd; 2 top-right odd; 3 bottom-right odd
bits 52:2   51-bit sub-rhombus row-major index (zero for polar zones)
bits 56:53  root rhombus 0:9, A north pole, B south pole
bits 61:57  associated ISEA9R level
bits 63:62  zero
```

The ISEA3H level is `2*isea9r_level + (subhex != 0)`. DGGAL bounds its
implementation at ISEA3H level 33 because a `3^16 by 3^16` rhombus is the
largest row-major index that fits the 51-bit field. Text IDs look like
`C0-2B-A`. The standard orders a parent's geometric sub-zones in tightly
packed scanlines, alternating scanline direction by parity.

Advantages are an official OGC identifier, a self-describing 64-bit value,
and a deterministic order for every geometric sub-zone set. Disadvantages
for this package's current interface are that it is not a prefix tree,
subtrees are not dense intervals in the global scanline order, and the OGC
geometric child relation is multi-parent. A package-specific “primary parent”
would still have to be chosen. Therefore:

- ship Z3 first for the strict-tree API and the existing 11.25-degree frame;
- later expose the OGC layout as `I3HZoneIndex`/text conversion, backed by the
  BSD DGGAL mathematics;
- only advertise `OGC/1.0/ISEA3H` when the 11.20-degree orientation and WGS84
  authalic conversion are also selected;
- declare `has_sorted_subtrees = false` for the OGC scanline ordering.

### ISEA4H

DGGRID's ZORDER digit string has the same two-digit roots, followed by one
base-4 digit per level. Roots `00:09` have four prefix children; polar roots
`10:11` have only the zero child. Use that prefix tree, but make the canonical
identity a 0-based `LevelIndex(level, ordinal)` in root-major/base-4 path
order:

```text
q < 10: ordinal = q*4^r + path_value
q = 10: ordinal = 10*4^r
q = 11: ordinal = 10*4^r + 1.
```

This permits levels `0:29` before `ncells` exceeds `Int64`. Parent and
children are digit arithmetic on the ordinal, and descendants are contiguous.
Expose `ISEA4HZOrderIndex` or a string codec through `reindex` if DGGRID
address compatibility is desired. As with ISEA3H, `.chd` is a geometric
overlap oracle, not the four-child primary tree.

### ISEA4T

Use `LevelIndex(level, face*4^level + path)` with zero-based face `0:19` and
the four digits defined above. The canonical dense position is `index+1`.
Parent divides the path by 4; children append digits `0:3`. Descendants are
contiguous and `levels = 0:29` is the largest range whose global cell count
fits `Int64`.

DGGRID `SEQNUM` is not this canonical order. Build an oracle crosswalk by
matching DGGRID level-0 polygon centroids to the package's 20 face records,
then classify each reported child geometrically against the four published
subtriangles. Never infer the canonical digit from the order of a `.chd` row.

## Point location, seams, neighbors, and boundaries

The implementation sequence should be:

1. rotate a world point into the fixed ISEA grid frame;
2. call `snyder_fwd` to get a face and planar coordinate;
3. quantize to the nearest triangular-lattice center (hex grids), or select a
   barycentric subtriangle (ISEA4T);
4. canonicalize face-edge duplicates into the Z3/ZORDER/face-path address;
5. on exact ties, return the smallest canonical ID among all incident-face
   candidates. Pin this rule with explicit seam vectors.

For hex grids, six planar lattice offsets give neighbors. Crossing an
icosahedron edge applies the corresponding face-to-face lattice transform;
at an icosahedron vertex one direction is deleted, leaving five neighbors.
DGGRID's neighbor output is counter-clockwise, but it does not specify the
package's required start direction. Rotate the verified cycle to a documented
start (for example the face-local +e1 edge); do not sort by ID.

For ISEA4T, edge neighbors can be derived directly from the face/path triangle
and the three icosahedron seam maps. Vertex connectivity is larger and
level-dependent at mesh vertices; enumerate incident triangles and order them
by spherical azimuth.

Snyder maps straight planar edges to curved spherical arcs. A ring containing
only inverse-projected planar corners is exact as a list of corners but a
consumer that connects them by great-circle arcs has an approximation to the
true ISEA boundary. The package must state one convention consistently:

- recommended initial convention: return the inverse-projected planar
  corners, matching existing IGeo7 behavior and DGGRID `densification 0`;
- compute analytic equal areas from the formulas above, not the great-arc
  approximation;
- add an internal edge-densification helper for rendering/intersection and an
  oracle generated with `densification 16`;
- do not describe the three/six-corner great-arc polygon as the exact Snyder
  curved boundary.

Hex children protrude beyond a parent's polygon. Measure the limiting subtree
cap overhang from levels 0 through at least 8 and set `cap_inflation` from the
maximum plus a floating-point margin. ISEA4T children are contained exactly in
the parent in chart space, but a spherical cap around only three mapped
corners still needs validation against densified Snyder edges.

## DGGRID oracle acquisition

The committed corpus is black-box output from a locally built DGGRID 9.0b
executable at commit `04bf5ed1372b174b9349faca8c265d112f6d8587`.
Neither the generator nor this research consulted DGGRID source. The exact
generator is [`scripts/oracles/dggrid/generate_isea_oracles.py`](../../scripts/oracles/dggrid/generate_isea_oracles.py),
and its complete CLI metafile is
[`scripts/oracles/dggrid/isea.meta.in`](../../scripts/oracles/dggrid/isea.meta.in).
Regenerate from the repository root with:

```sh
python3 scripts/oracles/dggrid/generate_isea_oracles.py \
  --dggrid /absolute/path/to/the/pinned/dggrid
```

The script refuses a binary whose stdout does not report `DGGRID version
9.0b`, checks every whole-earth record count against the formulas above, and
writes a byte count and SHA-256 for every output in its manifest. It fixes a
unit sphere, 15 fractional digits, `densification 0`, longitude wrapping, and
the 11.25-degree standard orientation. Boundaries are GeoJSON rather than the
older exploratory AIGEN output. Z3 padding is explicitly digit 3, avoiding the
DGGRID 8-to-9 default change.

The sealed suites are:

- [`test/oracles/ISEA3H/dggrid-9.0b`](../../test/oracles/ISEA3H/dggrid-9.0b):
  exhaustive levels 0:5; Z3 digit-string centers, GeoJSON boundaries, CCW
  neighbors, packed-INT64 Z3 center labels, and indexing-parent rows for
  levels 1:5;
- [`test/oracles/ISEA4H/dggrid-9.0b`](../../test/oracles/ISEA4H/dggrid-9.0b):
  exhaustive levels 0:4; ZORDER digit-string centers, GeoJSON boundaries,
  CCW neighbors, packed-INT64 ZORDER center labels, and indexing-parent rows
  for levels 1:4;
- [`test/oracles/ISEA4T/dggrid-9.0b`](../../test/oracles/ISEA4T/dggrid-9.0b):
  exhaustive levels 0:4 with SEQNUM centers and GeoJSON boundaries; DGGRID
  does not emit triangular neighbor or hierarchy files here; and
- [`test/oracles/ISEA3H/ogc-annex-b-orientation-dggrid-9.0b`](../../test/oracles/ISEA3H/ogc-annex-b-orientation-dggrid-9.0b):
  levels 0:3 centers and boundaries at longitude 11.20 degrees and authalic
  spherical latitude 58.282525588538995 degrees.

The last suite isolates the 0.05-degree orientation change. It is explicitly
not an OGC ZIRS oracle: the OGC latitude 58.397145907431 degrees is WGS84
geodetic, its sphere equivalent is the authalic latitude used above, and
DGGRID's Z3 labels are a different codec. Obtain OGC text/uint64 address
vectors from the OGC examples or DGGAL (BSD-3-Clause) and cross-check them by
centroid.

The committed `*-parents.txt` files are DGGRID **indexing-parent** output for
the Z3/ZORDER hierarchies. They are the appropriate oracle for the chosen
single-parent package trees. They are not the full geometric overlap relation
and must not be interpreted as proof that the central-place polygons are
nested. DGGRID geometric-child output was useful during exploration but is
not part of the sealed corpus.

Possible future additions, kept separate from the sealed baseline, are a
compact `densification 16` boundary corpus and deeper center-only levels. The
uncommitted exploratory conda-forge 8.44 runs established the CLI shape only;
they are not provenance for any committed fixture.

Finally create a deterministic `TRANSFORM_POINTS` corpus containing centers,
all level-0 face edges and vertices, points `2^-40` radians on each side of
those seams, selected cell-edge midpoints with both-sided perturbations, both
poles, and antimeridian cases. A locator metafile is:

```text
dggrid_operation TRANSFORM_POINTS
dggs_type @TYPE@
dggs_res_spec @LEVEL@
dggs_orient_specify_type SPECIFIED
dggs_vert0_lon 11.25
dggs_vert0_lat 58.282525588539
dggs_vert0_azimuth 0.0
precision 17
point_input_file_type TEXT
input_file_name /absolute/path/to/probes.txt
input_address_type GEO
input_delimiter ","
output_file_type TEXT
output_file_name /absolute/path/to/probe_addresses.txt
output_delimiter ","
@ADDRESS_BLOCK@
```

Record alongside every fixture set:

- executable version and SHA-256;
- the complete metafile and stdout parameter echo;
- SHA-256 for each output;
- line count checked against the formulas above;
- orientation name and whether coordinates are spherical or WGS84 geodetic;
- DGGRID's `z3_invalid_digit` value, even though digit strings avoid the
  packed-INT64 padding change between DGGRID 8 and 9.

## Required validation before implementation is called complete

1. Exhaustive ID/position bijection at every sealed level.
2. Every center agrees with the DGGRID oracle within a tolerance derived from
   its printed precision; inverse decode returns the same ID.
3. Boundary corners agree as cyclic rings, allowing only rotation of the
   starting vertex, never reflection.
4. Neighbor sets agree and the final public cycle is CCW from its documented
   start; pentagons have five, hexagons six.
5. Prefix parent/children are inverses and agree exactly with the DGGRID
   indexing-parent rows. Test the separate multi-parent geometric relation
   only if a future corpus intentionally records it.
6. ISEA4T canonical child geometry matches the four explicit midpoint
   triangles and their union tiles the parent in face space.
7. Counts, pentagon counts, and analytic areas match the formulas exactly.
8. The 11.20-degree corpus equals the 11.25-degree implementation after the
   expected rigid grid-frame rotation.
9. `node_extent` contains sampled descendants through the deepest practical
   oracle level, with the maximum overhang and safety margin recorded.
10. All seam/tie probes are deterministic on repeated runs and on both sides
    of every icosahedron edge.

## Remaining uncertainties/design work

- The package's exact 10-diamond seam tables for Z3/ZORDER still have to be
  written. DGGAL's BSD implementation and Sahr's Q2DI paper are acceptable
  references; DGGRID source is not.
- A Z3 prefix primary parent is a deliberate package hierarchy layered on the
  central-place geometry. It is not the full OGC sub-zone relation. Public
  documentation must say this plainly.
- The proposed ISEA4H `LevelIndex` order is package-defined even though its
  path is DGGRID-ZORDER-compatible. Seal the root and digit crosswalk before
  exposing IDs.
- Face numbering and child-digit gauge for ISEA4T are package choices. Pin
  them once against level-0/1 geometry; changing them later is an identifier
  breaking change.
- The finite boundary-ring convention cannot make a Snyder-curved edge exact.
  The interface documentation and `cell_area` behavior must be reconciled
  with the convention chosen above.
- PROJ's ISEA implementation has recently changed: current documentation and
  source branches do not always expose the same legacy `mode/aperture` grid
  parameters. Use PROJ only for projection cross-checks, not as a grid-address
  oracle.
