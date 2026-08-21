# IVEA and RTEA: mathematics, DGGS semantics, and oracle acquisition

Status: implementation research, 2026-08-13. This note deliberately uses only
papers/standards and permissively licensed software as implementation references.
DGGRID source code is excluded. DGGAL is usable because it is BSD-3-Clause, not
GPL. If DGGRID output is ever added as a second black-box oracle, preserve only
the generated facts and provenance, not code-derived logic.

## Recommended interpretation

The names in the old registry refer to DGGAL's families:

- **IVEA** is the *Icosahedral Vertex-oriented great-circle Equal Area*
  slice-and-dice projection of van Leeuwen and Strebe.
- **RTEA**, in this repository, should mean DGGAL's *Rhombic Triacontahedron
  (Snyder) Equal Area* projection, also written RT(S)EA in some material. It is
  the same 120-fundamental-triangle construction with the radial vertex moved
  from an icosahedron vertex to an icosahedron-edge midpoint, the centre of a
  rhombic-triacontahedron face.
- The suffix is the refinement aperture and cell type: `4R` and `9R` are
  rhombic cells; `3H` and `7H` are hexagonal/pentagonal cells. `_Z7` is an
  alternate zone-index representation of the same aperture-7 grid, not another
  projection.

The most important correction to the old metadata is that **RTEA does not have
30 canonical roots in DGGAL's DGGS hierarchy**. Thirty is the number of faces of
the rhombic triacontahedron. DGGAL maps the 120 fundamental triangles into the
same ten-rhombus 5x6 atlas used by its other rhombic-icosahedral grids. Thus:

- `IVEA4R`, `IVEA9R`, `RTEA4R`, and `RTEA9R` have 10 level-zero root cells;
- the `3H` and `7H` systems have 12 level-zero pentagons, not 10 or 30.

There is also a newer projection called RTEA in Wang et al. (2025), based on
two vertex-oriented triangles per rhombic-triacontahedron face. It must not be
silently substituted for the old registry's DGGAL family. The implementation
owner should record this choice explicitly; the recommendation here is DGGAL
RT(S)EA for compatibility with the registered names and identifiers.

## Primary and permissible sources

| Source | What it establishes | Software licence/version |
|---|---|---|
| van Leeuwen & Strebe, [“A Slice-and-Dice Approach to Area Equivalence in Polyhedral Map Projections”](https://doi.org/10.1559/152304006779500687) (2006), with [author/project abstract](https://www.tableau.com/research/publications/slice-and-dice-approach-area-equivalence-polyhedral-map-projections) | General great-circle slice-and-dice equal-area construction and IVEA | Paper, not a code dependency |
| Snyder, [“An Equal-Area Map Projection for Polyhedral Globes”](https://doi.org/10.3138/27H7-8K88-4882-1752) (1992) | Original equal-area radial construction and the icosahedral face-centred case | Paper |
| Hall et al., [“Disdyakis Triacontahedron DGGS”](https://doi.org/10.3390/ijgi9050315) (2020) | 120 spherical fundamental triangles, slice-and-dice formulae, inverse construction, and hierarchy context | Paper is CC BY 4.0 |
| Liang et al., [“Construction of a rhombic triacontahedron discrete global grid system”](https://doi.org/10.1080/17538947.2022.2130459) (2022) | Rhombic-triacontahedron DGGS construction and terminology | Open-access paper |
| Wang et al., [“Equal-Area Projection for Construction of Rhombic Triacontahedron Grid Systems”](https://doi.org/10.13203/j.whugis20220231) (2025), [publisher PDF](https://ch.whu.edu.cn/fileWHDXXBXXKXB/journal/article/whdxxbxxkxb/2025/3/0c85784c-3b41-4a62-8c9b-e48942e37eb9.pdf) | The newer, different vertex-oriented RTEA; useful for disambiguation, not the default compatibility target | Paper |
| OGC API — DGGS — Part 1, [OGC 21-038r1, Annex B](https://docs.ogc.org/is/21-038r1/21-038r1.html#annex-b) | Informative specifications for IVEA9R, IVEA3H and IVEA7H, including orientation, topology and ZIRS | OGC standard; IVEA annexes are informative |
| [DGGAL](https://github.com/ecere/dggal), [source at `e16cea7`](https://github.com/ecere/dggal/tree/e16cea7d930e603e09a8310edcd8f58218016e8f), [licence](https://github.com/ecere/dggal/blob/e16cea7d930e603e09a8310edcd8f58218016e8f/LICENSE) | Complete forward/inverse projection, atlas, identifiers, hierarchy, geometry and relationships for all requested IVEA/RTEA variants | BSD-3-Clause; pin commit `e16cea7d930e603e09a8310edcd8f58218016e8f` (2026-07-20) |
| [DGGAL 0.0.6 on PyPI](https://pypi.org/project/dggal/0.0.6/) | Stable packaged fallback | BSD-3-Clause; version 0.0.6, tag/commit `c323c4c444522f16fd6eba0c56ec65714a147c8c`, released 2025-11-13 |
| [A5](https://github.com/felixpalmer/a5), [source at `b4bc14c`](https://github.com/felixpalmer/a5/tree/b4bc14cb320c83c2a9f3f42410eeec1814f3e607), [licence](https://github.com/felixpalmer/a5/blob/b4bc14cb320c83c2a9f3f42410eeec1814f3e607/LICENSE) | Independent generic slice-and-dice equal-area implementation; useful for an IVEA local-triangle cross-check after matching vertex roles and affine output triangle | Apache-2.0; package version 0.9.0 at commit `b4bc14cb320c83c2a9f3f42410eeec1814f3e607` (2026-08-03) |
| [PROJ ISEA documentation](https://proj.org/en/stable/operations/projections/isea.html), [PROJ download/version page](https://proj.org/en/stable/download.html), [licence](https://proj.org/en/stable/about.html) | A useful negative result: PROJ implements Snyder ISEA, not IVEA or RTEA | MIT/X; current documentation and local inspection were for PROJ 9.8.1 (2026-04-10) |

DGGAL's most relevant source units at the pinned commit are
[`icoVertexGreatCircle.ec`](https://github.com/ecere/dggal/blob/e16cea7d930e603e09a8310edcd8f58218016e8f/src/projections/icoVertexGreatCircle.ec),
[`ri5x6.ec`](https://github.com/ecere/dggal/blob/e16cea7d930e603e09a8310edcd8f58218016e8f/src/projections/ri5x6.ec),
and the `RI3H`, `RI4R`, `RI7H`, and `RI9R` files under
[`src/dggrs`](https://github.com/ecere/dggal/tree/e16cea7d930e603e09a8310edcd8f58218016e8f/src/dggrs).
These are permissible implementation references under BSD-3-Clause.

The pinned DGGAL commit is one commit newer than 0.0.6 and specifically says it
fixes aperture-7 odd-level point-to-zone conversion. Oracle generation should
therefore use the pinned source commit, not silently use the older wheel. Pin
the associated eC source commit in the oracle manifest too.

## The common spherical construction

### Fundamental triangles and radial-vertex variants

An icosahedron face is divided by its three edge midpoints and face centre into
six spherical triangles. Across twenty faces this gives 120 congruent
disdyakis-triacontahedron fundamental triangles. Each has spherical excess

```text
Omega = 4*pi / 120 = pi/30 = 6 degrees.
```

Let the unpermuted fundamental-triangle vertices be:

```text
E = midpoint of an icosahedron edge
V = an icosahedron vertex
F = centre of an icosahedron face
```

Their spherical angles are respectively `90°`, `36°`, and `60°`. With
`phi = (1 + sqrt(5))/2`, the three side lengths used by DGGAL are

```text
|EV| = atan(1/phi)
|VF| = acos(sqrt((phi + 1)/3))
|FE| = atan(2/phi^2).
```

The projections differ only in which vertex is the radial point `A`:

| Projection | `A` (radial) | Other vertices | Angle at `A` |
|---|---|---|---:|
| IVEA | `V` | `E`, `F` | 36° |
| ISEA | `F` | `E`, `V` | 60° |
| DGGAL RTEA / RT(S)EA | `E` | `V`, `F` | 90° |

The description “RTEA swaps vertices” is therefore literal: the edge midpoint
becomes radial. That point is the centre of one of the 30 rhombic faces of the
dual rhombic triacontahedron.

### Forward projection on one oriented triangle

Let `(A,B,C)` be a consistently oriented unit-sphere triangle and
`(A',B',C')` its corresponding planar triangle. For input unit vector `P`:

1. Intersect great circle `AP` with great circle `BC`. Select the antipode `D`
   that lies on the minor arc `BC`.
2. Compute the area fraction
   `q = area_spherical(A,B,D) / Omega`.
3. Compute the radial fraction
   `h = sin(|AP|/2) / sin(|AD|/2)`.
4. Map by planar barycentric weights

```text
(b_A, b_B, b_C) = (1-h, h*(1-q), h*q)
P' = b_A*A' + b_B*B' + b_C*C'.
```

This is the useful implementation form of the slice-and-dice construction.
It preserves area because `q` partitions equal spherical area in the angular
direction while the half-chord ratio supplies the correct radial area scale.

Do not compute the radial fraction as a difference of nearly equal cosines.
DGGAL evaluates `sin(angle/2)` as the norm of a cross product with the
normalized midpoint, falling back to half the chord length for extremely small
angles.

### Inverse projection on one oriented triangle

First compute the planar barycentric weights `(b_A,b_B,b_C)`. Return the exact
corresponding vertex before doing divisions if a weight is within the chosen
vertex tolerance of one. Otherwise:

```text
h     = 1 - b_A
theta = (b_C/h) * Omega.
```

`theta` is `area_spherical(A,B,D)`. One stable closed form for recovering `D`
from this prescribed area, as used by DGGAL's vector inverse, is:

```text
c01 = dot(A,B)
c12 = dot(B,C)
c20 = dot(C,A)
s12 = sqrt(max(0, 1-c12^2))
T   = dot(A, cross(B,C))             # signed; keep the chosen winding
S   = sin(theta)
K   = 1-cos(theta) = 2*sin(theta/2)^2
f   = S*T + K*(c01*c12-c20)
g   = K*s12*(1+c01)
delta = 2*atan2(g,f)
D   = slerp(B,C, delta/|BC|).
```

Then recover the radial distance and point:

```text
x = 2*asin(clamp(h*sin(|AD|/2), -1, 1))
P = slerp(A,D, x/|AD|).
```

DGGAL further reduces the two SLERPs to vector operations on the normal path,
but retains the two-SLERP construction as the degeneracy fallback. The simpler
formula above is the better first Julia implementation: it is auditable and can
be checked directly against oracle vectors. Optimise only after correctness.

## Atlas, orientation, and coordinates

DGGAL unfolds the twenty projected triangular faces into ten unit rhombi in a
staircase embedded in a 5-by-6 planar chart. Face vertices and exact chart
incidences are the tables `vertices5x6` and `icoIndices` in permissively
licensed `ri5x6.ec`; reproduce them with attribution or derive an equivalent
table and check all seams against the oracle.

The OGC/DGGAL orientation places the first icosahedron vertex at authalic
latitude

```text
atan(phi) = 58.282525588538994... degrees
```

and longitude `11.20° E`, with the adjacent designated vertex due north. On the
WGS84 ellipsoid the first value corresponds to approximately
`58.397145907431°` geodetic latitude. This distinction is easy to miss:

- the mathematical projection operates on a sphere/authalic latitude;
- DGGAL's default geographic CLI input/output is geodetic WGS84 and performs
  authalic conversion internally;
- this repository's geometric core is a unit sphere, so it should not embed a
  WGS84 conversion in the projection kernel.

Oracle records should contain both DGGAL's geodetic coordinates and converted
authalic/unit-vector coordinates. Tests of the spherical kernel should compare
unit vectors; a separate wrapper can test WGS84 conversion.

## Refinement, cell counts, topology, and identifiers

For aperture `a`, the global counts are:

```text
rhombic aR: N(level) = 10*a^level
hexagonal aH: N(level) = 10*a^level + 2
```

Every `3H` or `7H` level has twelve pentagons. The other cells are hexagons.
On the unit sphere, an equal-area hexagon has area `4*pi/(10*a^level)` and a
pentagon has `5/6` of that area. This gives exactly `4*pi` globally.

| Variant | Scale/refinement | Global topology | Canonical DGGAL/OGC facts | DGGAL packed-ID maximum level |
|---|---|---|---|---:|
| `4R` | planar rows/columns double; 4 children geometrically | 10 root rhombi; four edge-neighbours | ID is level letter, root `0..9`, and uppercase hexadecimal row-major index | 25 |
| `9R` | planar rows/columns triple; 9 children geometrically | 10 root rhombi; four edge-neighbours | Same structure with scale `3^level`; IVEA9R is specified in OGC Annex B | 16 |
| `3H` | aperture 3 Goldberg refinement | 12 pentagons, all others hexagons; 5 or 6 neighbours | even `l=2k`: `GP(3^k,0)`; odd `l=2k+1`: `GP(3^k,3^k)`; IVEA3H ZIRS in OGC Annex B | 33 |
| `7H` | aperture 7 Goldberg refinement | 12 pentagons, all others hexagons; 5 or 6 neighbours | even `l=2k`: `GP(7^k,0)`; odd `l=2k+1`: `GP(2*7^k,7^k)`; normal and `_Z7` representations; IVEA7H in OGC Annex B | 19 |

The 64-bit limits above are DGGAL representation limits, not mathematical
limits. The repository can initially expose the same limits for interoperability.

### Hierarchy is not always a tree

For `4R` and `9R`, the geometric hierarchy has one parent and respectively four
or nine children. For the hexagonal refinements, boundary overlap means a cell
can have more than one geometric parent: DGGAL reports up to three parents for
`3H` and two for `7H`, and up to seven or thirteen children respectively.
DGGAL also defines a primary/centroid parent.

Consequences for this repository's interface:

1. The singular `parent(cell)` is a design decision, not a fact that can be
   inferred from `radix`. The recommended interoperable choice is DGGAL's
   primary/centroid parent.
2. If possible, expose all covering parents separately. Tests should distinguish
   `primary_parent` from `covering_parents`.
3. `children(parent)` must be defined consistently with the chosen parent
   relation; DGGAL's full geometric children are not an ordinary radix tree.
4. Set `has_sorted_subtrees = false` for the canonical scanline/row-major ZIRS.
   Even rhombic row-major descendants are generally not one contiguous ordinal
   interval. A Morton-style alternate index could be added later.

### Canonical identifiers

Use DGGAL/OGC ZIRS as the compatibility identifier rather than inventing a new
one. In particular, the level-zero IDs observed in DGGAL are:

- rhombic variants: `A0-0` through `A9-0`;
- hexagonal variants: twelve pentagons `A0-0-A` through `AB-0-A`.

Use the OGC-defined representations for IVEA3H, IVEA7H and IVEA9R. Use DGGAL's
corresponding representation for IVEA4R and the RTEA variants. Treat `7H_Z7` as
an alternate encoding reachable through `reindex`; it must share cell geometry
and relationships with ordinary `7H`.

## Cell boundaries and neighbours

Cells are straight-sided polygons in the equal-area 5x6 plane, inverse-projected
to the sphere. Their spherical edges are generally **not great-circle arcs**.
Therefore a list of inverse-projected corners is sufficient for topology but
not necessarily for accurate spherical polygon area, extent, cap or rendering.
DGGAL distinguishes ordinary vertices from refined WGS84 vertices for this
reason.

The implementation must choose and document one of these contracts:

- `cell_boundary` returns only topological corners, while a separate refined
  boundary method samples the actual projected edge; or
- `cell_boundary` returns a tolerance-controlled densified curve.

The first is simpler and preserves stable vertex counts. In either case, oracle
files should store both corners and a refined ring.

DGGAL's rhombic relationship code enumerates the four planar directions as
top, left, right and bottom. That is not the repository's required cyclic CCW
neighbour ordering. For all variants, capture the oracle's direction code and
geometry, then sort/remap around the cell centre into the repository's chosen
CCW rotational cycle. Test the order independently at seams and pentagons; do
not assume the library's iteration order is already CCW.

## PROJ and other permissive implementations

PROJ 9.8.1's projection registry contains `+proj=isea`, a spherical Snyder ISEA
implementation. It does not contain `ivea`, `rtea`, or a rhombic-
triacontahedron equal-area projection. Its ISEA inverse is limited to plane
output mode. Thus PROJ cannot generate IVEA/RTEA oracle vectors and should not
be used as a name-compatible substitute.

A5's Apache-2.0 equal-area module is an independent implementation of the same
generic slice-and-dice construction. Its dodecahedral/Snyder-equivalent framing
can validate the IVEA fundamental-triangle formula after explicitly matching
the three spherical vertex roles and applying an affine transformation to the
same planar triangle. A5 does not provide DGGAL's IVEA/RTEA 5x6 atlas or ZIRS,
so it is only a projection-kernel cross-check.

No second permissively licensed, complete RTEA implementation was found. For
RTEA, DGGAL should be the external software oracle, with independence supplied
by the published formula and strong invariants: forward/inverse round trip,
constant Jacobian/area, all 120 triangle vertices and centroids, icosahedral
symmetries, global cell count/area, and seam incidence. Do not claim two
independent implementations where only one exists.

## Oracle acquisition

### Reproducible source build

Prefer a pinned source build on Linux in a disposable directory. DGGAL's build
expects its eC checkout as a sibling. Record both resolved commits; the eC line
below intentionally records rather than assumes its moving default revision.

```sh
mkdir oracle-src && cd oracle-src
git clone https://github.com/ecere/eC.git
git clone https://github.com/ecere/dggal.git
git -C dggal checkout e16cea7d930e603e09a8310edcd8f58218016e8f
git -C eC rev-parse HEAD > eC.commit
git -C dggal rev-parse HEAD > dggal.commit
make -C eC -j4
make -C dggal -j4
```

Use the executable produced by that build and record its path, compiler, target,
flags and shared-library hashes. The packaged fallback is:

```sh
python3.13 -m venv .oracle-venv
.oracle-venv/bin/python -m pip install 'dggal==0.0.6'
.oracle-venv/bin/dgg ivea3h info
```

The fallback is unsuitable as the only aperture-7 point-to-zone oracle because
it precedes the pinned fix. Also, the 0.0.6 macOS arm64 wheel inspected during
this research bundled x86_64 `libecrt.dylib`/`libdggal.dylib` with an arm64
extension, and Python 3.14 had no compatible wheel. A Rosetta x86_64 Python 3.9
environment ran it, but a Linux source build is the reproducible route.

### CLI reconnaissance commands

Run these for every concrete name `ivea4r`, `ivea9r`, `ivea3h`, `ivea7h`,
`ivea7h_z7`, `rtea4r`, `rtea9r`, `rtea3h`, `rtea7h`, and `rtea7h_z7` where the
CLI exposes it:

```sh
dgg ivea3h info
dgg ivea3h level 3
dgg ivea3h list 0
dgg ivea3h grid 3 > ivea3h-level3-crs84.geojson
dgg ivea3h -crs ico grid 3 > ivea3h-level3-ico.geojson
dgg ivea3h geom A4-0-A
dgg ivea3h -crs ico geom A4-0-A
dgg ivea3h zone 34,-70
dgg ivea3h sub A4-0-A -depth 3
dgg ivea3h index A4-0-A B2-5-C
```

The `zone` CLI defaults its level from the DGGRS setting; set it explicitly in
the scripted API or via the documented CLI option rather than relying on a
session default.

For structured files, use DGGAL's Python or C binding rather than parsing the
human-readable `info` command. The needed API surface is:

```text
getZoneTextID, getZoneLevel
getZoneWGS84Centroid, getZoneCRSCentroid
getZoneWGS84Vertices, getZoneRefinedWGS84Vertices
getZoneNeighbors, getZoneParents, getZoneChildren
countSubZones, getSubZoneAtIndex
getZoneFromWGS84Centroid / point-to-zone equivalent
```

The binding example
[`bindings_examples/py/info.py`](https://github.com/ecere/dggal/blob/e16cea7d930e603e09a8310edcd8f58218016e8f/bindings_examples/py/info.py)
is the permissible starting point for the generator.

### Proposed committed oracle layout

Use one directory per concrete variant and a shared versioned schema:

```text
test/oracles/dggal/<variant>/
  manifest.toml
  levels.csv
  cells.jsonl
  projection.csv
  lookup.csv
  hierarchy.jsonl
```

`manifest.toml` should contain:

```text
schema version; UTC generation time; generator commit and SHA-256
DGGAL commit; eC commit; DGGAL licence identifier
compiler/OS/architecture; exact commands and environment
ellipsoid and authalic-conversion constants
icosahedron orientation and longitude convention
variant, maximum sampled level, random seed
SHA-256 for every data file
```

`levels.csv` should contain `variant, level, global_cell_count,
hex_area_unit_sphere, pent_area_unit_sphere`, including exact formula strings as
well as evaluated values.

Each `cells.jsonl` record should contain:

```text
variant, id, packed_u64_if_applicable, level, enumeration_ordinal, cell_kind
centre_geodetic_lonlat, centre_authalic_lonlat, centre_unit_xyz, centre_5x6_xy
corner_5x6_xy, corner_geodetic_lonlat, corner_authalic_xyz
refined_boundary_geodetic_lonlat, refined_boundary_authalic_xyz
neighbours [{direction_code, id}], covering_parents, primary_parent, children
reported_area, expected_area
```

`projection.csv` is specifically a projection oracle, independent of the cell
schema. Store geographic input in both geodetic and authalic form, unit XYZ,
5x6 output, selected face and fundamental triangle if available, inverse result,
and error. Include all 12 icosahedron vertices, 20 face centres, 30 edge
midpoints, 120 fundamental-triangle centroids, interior seeded points, all seams,
and points on both sides of seams.

`lookup.csv` should contain fixed point-to-cell cases: poles; equator/cardinal
longitudes; every special vertex/edge/face centre; a deterministic global
lattice; cell centres; and seam/vertex probes perturbed by approximately
`1e-12` radians and `2^-40` chart units. Store the tie-case classification so a
future implementation does not accidentally treat an arbitrary boundary result
as universally canonical.

`hierarchy.jsonl` should preserve both primary and all covering relationships,
plus DGGAL subzone counts and indices. This catches the most likely false-tree
implementation error in `3H`/`7H`.

Generate complete low levels rather than only random examples. A reasonable
first corpus is all cells through level 3 for `4R`, `9R`, and `3H`, all cells
through level 2 or 3 for `7H`, plus deterministic deep-level samples up to the
packed-ID limits. Always include every pentagon and every seam-adjacent cell at
each sampled level.

### Oracle validation before commit

The generator should fail unless all of the following hold:

1. observed counts equal `10*a^level` or `10*a^level+2`;
2. exactly twelve hex-family cells are pentagons;
3. reported areas sum to the WGS84/authalic globe area and normalized areas sum
   to `4*pi`, within a declared tolerance;
4. neighbour relationships are reciprocal and every neighbour shares the
   expected projected boundary;
5. all projected and geographic centres round-trip;
6. ordinary and `_Z7` identifiers reindex to identical geometry;
7. primary parents occur in covering parents, and hierarchy counts agree with
   DGGAL's subzone APIs;
8. every chart seam has matching inverse geometry from each incident face;
9. A5 agrees with matched IVEA local-triangle cases; and
10. two independent runs in clean environments produce byte-identical sorted
    facts, apart from explicitly excluded timestamp/build metadata.

## Numerical and semantic pitfalls

- **Face ties:** select faces/fundamental triangles with a deterministic rule.
  DGGAL's face logic uses explicit `<`/`<=` choices and an approximately
  `1e-11` nudge. Record exact boundary probes rather than hiding this rule.
- **Great-circle antipodes:** `cross(cross(A,P), cross(B,C))` returns a line,
  hence two antipodal candidates. Choose the point on minor arc `BC`.
- **Vertex division:** inverse barycentrics require `b_C/(1-b_A)`; return exact
  vertices before division.
- **Domain drift:** clamp dot products and `asin` inputs and renormalize output.
- **Cancellation:** use `1-cos(theta)=2*sin(theta/2)^2` and a stable half-chord.
- **Quadrants:** use `atan2`. DGGAL notes that a published Hall/Lee–Mortari
  rearrangement needs both numerator and denominator signs negated to obtain the
  correct `atan2` quadrant for non-right RTEA subtriangles.
- **Optimiser sensitivity:** DGGAL compiles the projection paths with unsafe
  math optimisations disabled and retains a two-SLERP degeneracy fallback.
- **Pole longitude:** longitude is undefined at a pole. Compare unit vectors.
- **Coordinate model:** never compare DGGAL geodetic latitude directly with a
  spherical/authalic implementation.
- **Curved edges:** inverse-projected straight chart edges are not generally
  great circles; densify before testing spherical extents or polygon areas.
- **Neighbour order:** library direction enumeration is not proof of the
  repository's required CCW rotational order.
- **RTEA naming:** keep the DGGAL RT(S)EA and 2025 Wang RTEA definitions separate.

## Decisions still required

1. Confirm that all registered RTEA variants mean DGGAL RT(S)EA, not Wang et
   al.'s 2025 projection. Recommended: DGGAL for registry compatibility.
2. Confirm DGGAL/OGC ZIRS as canonical for the concrete IVEA variants and DGGAL
   ZIRS for RTEA. Recommended: yes; `_Z7` is an alternate codec.
3. Define the singular `parent` contract for `3H`/`7H`. Recommended: DGGAL's
   primary/centroid parent, with a separate all-covering-parents relation.
4. Decide whether `cell_boundary` means corners or a densified curved boundary.
   Recommended: corners, plus an explicit refined-boundary API.
5. Decide the spherical-core versus WGS84-wrapper split. Recommended: authalic
   unit sphere in the kernel, explicit ellipsoidal conversion outside it.
6. Choose the start direction for canonical CCW neighbour cycles, then derive
   and test the remapping from DGGAL direction codes at ordinary cells, seams
   and all twelve pentagons.
7. Pin an eC commit alongside DGGAL after proving the clean source build; do not
   leave the oracle dependent on eC's moving default branch.
8. Decide whether all DGGAL maximum levels are public limits or merely oracle
   sampling limits. Recommended initially: match them for 64-bit ZIRS interop.

With those decisions recorded, these grids are implementable from the common
projection kernel plus four refinement/indexing kernels. IVEA and DGGAL RTEA
share all atlas and DGGS logic; only the radial-vertex permutation changes in
the projection. The external oracle should be acquired before implementation,
with the pinned BSD DGGAL build and the schema above.
