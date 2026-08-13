# rHEALPix and AusPIX: implementation math and oracle plan

## Conclusions

The implementation target is well specified. rHEALPix is not a new spherical
polyhedral kernel: it is HEALPix's equal-area map projection followed by a
piecewise rigid rearrangement of its eight polar triangles, then an ordinary
square quadtree-like hierarchy with `3 × 3` refinement. The exact hierarchy
used here has six roots `N O P Q R S`, nine row-major children `0…8`, and
prefix identifiers.

AusPIX does **not** define another cell codec in the available authoritative
implementation. It is a profile of rHEALPix:

- WGS84 ellipsoid;
- Greenwich/zero central meridian;
- `north_square = 0`, `south_square = 0`;
- `N_side = 3` (aperture 9);
- the same identifiers, such as `R7751215231`.

The old registry claim `:auspix_id_with_rhealpix_ordinal` should therefore be
replaced by “rHEALPix SUID under the AusPIX profile” unless a separate GA
standard defining another ordinal is found. The OGC Testbed-16 report explicitly
says that AusPIX uses the default hard-coded WGS84 system based on the zero
meridian and instantiates related profiles with `(0,0), N_side=3`.

For this repository's unit-sphere interface, implement rHEALPix on the unit
authalic sphere. AusPIX should be a profile/I/O wrapper that converts WGS84
geodetic latitude to authalic latitude before using that spherical kernel. It
should not pretend that geodetic latitude itself is equal-area on a unit sphere.

Primary sources are the complete 2013 definition and algorithms in [Gibb,
Raichev, and Speth's preprint](https://raichev.net/files/rhealpix_dggs_preprint.pdf)
([archived DOI record](https://doi.org/10.7931/J2D21VHM)), the shorter published
[2016 rHEALPix paper](https://doi.org/10.1088/1755-1315/34/1/012012), and the
underlying HEALPix projection derivation by [Calabretta and Roukema
(2007)](https://doi.org/10.1111/j.1365-2966.2007.12297.x). AusPIX profile
evidence comes from the [Geoscience Australia repository](https://github.com/GeoscienceAustralia/AusPIX_DGGS)
and the [OGC Testbed-16 engineering report](https://docs.ogc.org/per/20-039r2.html).

## Reference implementations and licensing

| source | use | license finding |
|---|---|---|
| [`manaakiwhenua/rhealpixdggs-py` v0.6.0, commit `5929a73`](https://github.com/manaakiwhenua/rhealpixdggs-py/tree/5929a73b427a33d66a64800051d077ce36bbf901) | pinned black-box oracle and a readable implementation of the paper | Its [top-level license](https://github.com/manaakiwhenua/rhealpixdggs-py/blob/5929a73b427a33d66a64800051d077ce36bbf901/LICENSE) expressly offers LGPL-3.0 **or MIT**, at the user's choice; the [MIT text](https://github.com/manaakiwhenua/rhealpixdggs-py/blob/5929a73b427a33d66a64800051d077ce36bbf901/LICENSE-MIT) is present. We elect MIT. |
| [PROJ rHEALPix](https://proj.org/en/stable/operations/projections/rhealpix.html) | independent forward/inverse projection cross-check, not the grid hierarchy | [MIT](https://proj.org/en/stable/about.html). |
| [`GeoscienceAustralia/AusPIX_DGGS` commit `4120e2d`](https://github.com/GeoscienceAustralia/AusPIX_DGGS/tree/4120e2da5b61dafadc17aa719e71dd0c4bf19fbf) | evidence for the AusPIX profile and published IDs | Repository [license is Apache-2.0](https://github.com/GeoscienceAustralia/AusPIX_DGGS/blob/4120e2da5b61dafadc17aa719e71dd0c4bf19fbf/LICENSE). |

There are two license hygiene warnings. First, rHEALPixDGGS 0.6.0's
`pyproject.toml` still says GPL even though the root license, MIT deed, and
README explicitly grant the MIT alternative. The root grant is unambiguous,
but this metadata discrepancy should be recorded in any legal inventory.
Second, the old AusPIX repository bundles upstream engine files carrying old
LGPL headers despite its Apache top-level license. Do not copy that bundled
engine. Use GA's own wrapper/profile material as evidence and use the modern
dual-MIT upstream or the papers for math. The oracle generator imports only
the pinned modern upstream and elects MIT.

As an independent check, the generator compares the pure-Python spherical
projection against MIT-licensed PROJ for every polar placement on a fixed
global grid. With pyproj 3.7.2 / PROJ 9.5.1, the largest discrepancy is
recorded in the manifest and is at floating-point roundoff scale (below
`2e-14` authalic radii).

## Projection mathematics

All formulas below use radians. Let geodetic longitude and latitude be
`(λ, φ)`, ellipsoid eccentricity `e`, semi-major axis `a`, central meridian
`λ₀`, authalic latitude `β`, and authalic radius `R_q`. Normalize longitude
first: `λ′ = wrap[-π,π)(λ - λ₀)`.

For an ellipsoid,

```text
q(φ) = (1-e²) sinφ / (1-e² sin²φ)
       - (1-e²)/(2e) log((1-e sinφ)/(1+e sinφ))

q_p  = 1 - (1-e²)/(2e) log((1-e)/(1+e))
β    = asin(q(φ)/q_p)

R_q  = a sqrt( 1/2 * [1 - (1-e²)/(2e) log((1-e)/(1+e))] ).
```

For a sphere, `β = φ` and `R_q = a`. The inverse `β → φ` has no simple closed
form; use a well-tested authalic-latitude inverse series or a monotone root
solve. WGS84's transition in geodetic latitude is the inverse-authalic image
of `asin(2/3)`.

### HEALPix stage

Let `β₀ = asin(2/3)`. Work first in authalic-radius units.

For `|β| ≤ β₀`:

```text
X = λ′
Y = (3π/8) sinβ.
```

For `|β| > β₀`:

```text
σ   = sqrt(3(1 - |sinβ|))
c   = min(3, floor(2λ′/π + 2))
λ_c = -3π/4 + cπ/2
X   = λ_c + (λ′ - λ_c)σ
Y   = sign(β) (π/4)(2 - σ).
```

The physical projected point is `(x,y)=R_q(X,Y)`. The inverse is:

```text
|Y| ≤ π/4:
    λ′ = X
    β  = asin(8Y/(3π))

π/4 < |Y| < π/2:
    τ   = 2 - 4|Y|/π
    c   = min(3, floor(2X/π + 2))
    λ_c = -3π/4 + cπ/2
    λ′  = λ_c + (X - λ_c)/τ
    β   = sign(Y) asin(1 - τ²/3).
```

At a pole the reference canonicalizes longitude to `-π`. Finally apply the
inverse authalic transform and add `λ₀`.

### Polar-square rearrangement

Let `M = [0 -1; 1 0]`, an anticlockwise quarter-turn. The forward HEALPix
polar triangle number is selected by `X` intervals
`[-π,-π/2), [-π/2,0), [0,π/2), [π/2,π]`, giving `c=0,1,2,3`.
The transition lines belong to the equatorial region; north is `Y>π/4` and
south is `Y<-π/4`.

For north-square position `n` and south-square position `s`, define triangle
tips

```text
t_c^+ = (-3π/4 + cπ/2,  π/2)
t_c^- = (-3π/4 + cπ/2, -π/2)
u_n   = (-3π/4 + nπ/2,  π/2)
u_s   = (-3π/4 + sπ/2, -π/2).
```

The rearrangement is the rigid motion

```text
north: p_r = M^(c-n) (p_h - t_c^+) + u_n
south: p_r = M^(-(c-s)) (p_h - t_c^-) + u_s
equator: p_r = p_h.
```

The inverse applies the opposite rotations. Inside a polar square, select
`c` by its two diagonals. The paper gives the exact case table in Appendix B.
The pinned implementation uses that table with `1e-15` fuzz at diagonals;
`projection_unit_all_polar_placements.csv` deliberately seals seam behavior
for all 16 `(n,s)` placements. A native implementation should use an explicit
half-open diagonal ownership rule, then test that it agrees with these vectors,
rather than spreading epsilons throughout the kernel.

The rearrangement consists only of rotations and translations, so it preserves
the HEALPix equal-area factor and orientation.

## Planar grid, identifiers, and hierarchy

For aperture 9 (`N_side=3`), the level-0 square width is
`w₀ = R_q π/2`. Root upper-left corners in authalic-radius units are:

| root | upper-left `(X,Y)` |
|---|---|
| `N` | `(-π + nπ/2, 3π/4)` |
| `O` | `(-π, π/4)` |
| `P` | `(-π/2, π/4)` |
| `Q` | `(0, π/4)` |
| `R` | `(π/2, π/4)` |
| `S` | `(-π + sπ/2, -π/4)` |

Each square is split into row-major children:

```text
0 1 2
3 4 5
6 7 8
```

For ID `root d₁…d_l`, write `row_i = d_i ÷ 3` and `col_i = d_i % 3`.
Then

```text
w_l = w₀ 3^-l
ul_x = root_ul_x + w₀ Σ(col_i 3^-i)
ul_y = root_ul_y - w₀ Σ(row_i 3^-i)
nucleus = (ul_x + w_l/2, ul_y - w_l/2).
```

The exact cell count is `6·9^l`; its equal area is
`R_q² (2π/3) 9^-l` on the ellipsoid, or `2π/(3·9^l)` steradians on the unit
sphere. Parent truncates the last digit; children append `0…8`. Prefix order is
subtree order. At a fixed level, the zero-based dense ordinal is

```text
root_number·9^l + Σ(d_i·9^(l-i)), root_number(N…S)=0…5.
```

Thus `cellposition = ordinal0 + 1` in Julia. A `UInt64` dense ordinal supports
only levels `0:19`; the mathematics is unbounded, but `levels(sys)` must report
the representable range of the chosen cell type. A structural root-plus-digits
cell type or `LevelIndex(level, ordinal)` both work; the former preserves the
published text ID without repeated base conversion.

The same SUID is the AusPIX ID. There is no evidence for a second
`auspix_id_with_rhealpix_ordinal` scheme.

## Geometry primitives

### Boundary

The true boundary is defined exactly by parameterizing the four straight edges
of the planar square and applying the inverse rHEALPix projection pointwise.
This yields the published quad, skew-quad, dart, and cap shapes. The finite
vertex list alone is not an exact boundary for curved edges. Darts have only
three geometric vertices; cap cells have no ordinary polygon corners, even
though four planar corner images are useful topological seam markers.

This exposes an interface mismatch: `cell_boundary` currently asks for an
“exact boundary ring” made only of points. No finite point ring exactly
represents general rHEALPix edges. The implementation must either:

1. return a documented fixed/adaptive densification (the initial practical
   choice, consistent with existing HEALPix handling), or
2. extend the geometry vocabulary with parameterized boundary arcs.

The generated oracle records five samples per planar edge, 16 perimeter points
with implicit closure. It starts at the planar upper-left corner and proceeds
counterclockwise as viewed from outside the authalic sphere. Each sample records
separate geodetic and authalic longitude/latitude coordinates; its 3-D unit
direction is formed from the authalic coordinates. These are exact samples of
the true edge parameterization, not a claim that the connecting great-circle
chords are exact.

### Nucleus versus centroid

The inverse image of the planar square center is the published **nucleus**. It
is interior, deterministic, recursively aligned for odd `N_side`, and is the
right implementation of this repository's `cell_centroid` contract, which only
requires a representative interior point. Document the terminology because it
is generally not an area centroid.

The paper's `centroid` is an average of longitude and latitude over the cell,
not the normalized 3-D spherical area centroid. Cap centroids equal their
nuclei; quad longitude equals nucleus longitude and latitude is the mean of
latitudinal extremes; darts preserve nucleus longitude; skew-quad/dart latitude
requires numerical integration because inverse authalic latitude has no closed
form. `centroids.jsonl` records this reference quantity only as a diagnostic and
must be compared with tolerance.

### Cell-at and boundary ownership

Project the input to the rHEALPix plane, select a root square, then compute
normalized offsets `dx,dy ∈ [0,1]`. The first `l` base-3 digits of `dy` and
`dx` give row and column at every depth, with `digit=3·row+col`. This is O(level)
and needs no search.

The paper's ownership convention is half-open: `N` and `S` contain none of
their outer edges; `O…R` contain left, top, and bottom edges; recursively a child
contains left and top edges unless it shares a parent boundary, when it inherits
the parent's ownership. The Python oracle implements this mostly through strict
root comparisons and `floor`, but also nudges an exactly-one normalized offset
inward by half a finest-cell width. That nudge depends on the reference object's
artificial `max_resolution` and should **not** become production semantics.
Choose the paper's recursive half-open rule as canonical and use
`boundary_ties.csv` to identify any deliberate mismatches with the reference.
For every root and level-1 source cell, that corpus crosses the outward normal
of each of the four planar edges (vertical perturbations on horizontal edges
and horizontal perturbations on vertical edges), and probes all four 2-D
quadrants plus the exact point at every corner. Probes that leave the projection
image are retained explicitly with `in_projection_image=false`.

### Neighbors

Every cell has four edge neighbors. The root edge topology for `(n,s)=(0,0)` is:

```text
N: down O, right P, up Q, left R
O: up N, right P, down S, left R
P: up N, right Q, down S, left O
Q: up N, right R, down S, left P
R: up N, right O, down S, left Q
S: up O, right P, down Q, left R
```

Within a face, move row/column normally. Crossing a root edge maps to the
adjacent root and rotates every remaining base-9 digit by 0, 1, 2, or 3 quarter
turns. A digit `(row,col)` rotated anticlockwise becomes
`(N-1-col,row)`; applying this recursively gives the reference neighbor
automaton. The sealed `edge_neighbors` direction maps cover all cells through
level 2, including every polar seam.

For the interface's default `Vertex()` connectivity, roots have four neighbors;
other cells have seven at cube/dart singularities or eight ordinarily. The
generated `vertex_neighbors.jsonl` derives these sets independently by boundary
corner incidence. It intentionally sorts IDs and is not an order oracle.

For the required rotational API order, use the local projected-up direction as
the documented start and sort neighbor nuclei by tangent-plane azimuth CCW.
For `Edge()`, the orientation-preserving map gives the simple CCW cycle
`up, left, down, right`. For `Vertex()`, azimuth sorting correctly interleaves
corner and edge neighbors and handles seven-neighbor singularities.

## Mapping to the package interface

| interface operation | implementation |
|---|---|
| `ncells(levelgrid(sys,l))` | `6·9^l` with checked integer arithmetic |
| `cellindex` / `cellposition` | root/base-9 ordinal formula above |
| `parent` / `children` | truncate / append digits |
| `has_sorted_subtrees` | `true` |
| `descendant_range` | prefix interval `[ord·9^k, (ord+1)·9^k-1]` at depth `k` |
| `cell_centroid` | inverse projection of planar nucleus |
| `cell_boundary` | inverse-project sampled planar square perimeter; document densification |
| `cellat` | forward project, root selection, simultaneous base-3 expansion |
| `neighbors(..., Edge())` | face-grid step plus cross-face digit rotation |
| `neighbors(..., Vertex())` | edge+corner topology, azimuth ordered |
| `node_extent` | descendants lie inside the parent cell exactly; a sound cap over a densified perimeter plus conservative margin is sufficient initially |
| `max_neighbors` | 4 for `Edge()`, 8 for `Vertex()` |

## Generated oracle artifacts

The generator is [`scripts/oracles/rhealpix/generate_rhealpix_oracles.py`](../../scripts/oracles/rhealpix/generate_rhealpix_oracles.py),
with its exact environment in [`requirements.txt`](../../scripts/oracles/rhealpix/requirements.txt).
It produced [`test/oracles/rhealpix`](../../test/oracles/rhealpix):

| file | coverage |
|---|---|
| `projection_unit_all_polar_placements.csv` | forward/inverse unit-sphere projection, all 16 `(n,s)` placements, equatorial/polar transitions, poles, seams |
| `projection_auspix_wgs84.csv` | WGS84/AusPIX forward/inverse points in metres/degrees, including Australian cities |
| `hierarchy.csv` | every cell through level 3: ID, ordinal, parent, children, upper-left corner, width |
| `cells.jsonl` | both profiles, every cell through level 2: area, shape, geodetic/authalic nucleus, authalic unit direction, 16 CCW boundary samples in all three coordinate views, four edge neighbors |
| `lookup.csv` | deterministic geodetic point-to-cell probes at levels 0,1,2,3,5,10 |
| `boundary_ties.csv` | all four edge-normal crossings and all four corner quadrants, plus exact points, at levels 0,1,2,5; projection-image exits are explicit |
| `vertex_neighbors.jsonl` | complete unit-sphere vertex-neighbor sets through level 2 |
| `centroids.jsonl` | roots plus four examples per cell shape for nucleus/coordinate-centroid distinction |
| `manifest.json`, `SHA256SUMS` | schema version 2, coordinate semantics, CCW winding, source commit, elected license, versions, profiles, checksums |

Generation and checksum verification commands are in the adjacent
[`README.md`](../../scripts/oracles/rhealpix/README.md). The generator was run
twice in the isolated worktree and reproduced byte-identical checksums.

## Remaining decisions and uncertainties

1. **Finite ring versus exact curve:** decide package-wide whether a densified
   point ring is allowed to stand for a curved exact boundary. This is the only
   mathematical/interface blocker.
2. **Canonical system parameters:** use `(0,0), N_side=3, lon_0=0` as the
   no-argument rHEALPix system because that is the paper's running example and
   AusPIX profile. If arbitrary `n,s` are exposed, they are part of system
   identity and change neighbor/ID geometry even though the ID grammar stays
   the same.
3. **Centroid naming:** implement the nucleus for `cell_centroid`; do not silently
   use the reference package's longitude/latitude average.
   On ellipsoidal profiles, never form its unit direction directly from geodetic
   latitude: convert geodetic latitude to authalic latitude first.
4. **Boundary ties:** adopt the paper's recursive half-open rule. Matching the
   Python `max_resolution`-dependent nudge exactly would be a design regression.
5. **AusPIX scope:** the available GA material establishes a profile and use
   conventions, not a separate mathematical DGGS or codec. Treat AusPIX as an
   rHEALPix wrapper until a normative source proves otherwise.
6. **License metadata discrepancy:** the selected upstream commit gives an
   explicit MIT choice, but retain the license evidence with the fixture and do
   not copy the older LGPL-header AusPIX engine.
