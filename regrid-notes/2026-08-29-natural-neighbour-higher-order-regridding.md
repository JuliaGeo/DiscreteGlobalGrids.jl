# Natural-neighbour and higher-order regridding

**Status:** research and design note, not an implementation plan.

**Package decision:** implement these algorithms directly in
`GlobalRegridding.jl`. A separate `NaturalRegridding.jl` package is not needed
for the first implementation. The natural-neighbour coordinate calculation can
live as an internal geometry kernel. It should be extracted only if it later
becomes useful independently of global regridding.

This note extends
[`generic-barycentric-patch-regridding.md`](generic-barycentric-patch-regridding.md)
and [`esmf-regridding-methods.md`](esmf-regridding-methods.md).

## Conclusion

The central idea holds in a useful but weaker form.

- On a perfect planar triangular lattice, Sibson coordinates inside a reference
  hexagonal Voronoi cell are universal, compact, non-negative functions. A
  generic query has four or five positive source weights, all drawn from the
  containing site and its first ring.
- Farin, Sibson-1, and Hiyoshi interpolants can also become one fixed sparse
  matrix when derivative recovery is a fixed linear operator. Their effective
  support is the composition of the natural-neighbour basis with the derivative
  estimator, not an intrinsic property of the interpolant.
- The proposed rule `C0 -> one ring, C1 -> two rings, C2 -> three rings` is
  therefore a useful design target, not a theorem. It can be achieved for some
  derivative estimators. The defaults in `NaturalNeighbours.jl` can reach four
  rings.
- One exact universal stencil does not exist across a spherical DGGS. Local
  geometry varies with position and resolution; pentagons, face seams, and
  Delaunay flips are genuine exceptional cases. Reference stencils can still be
  an accelerated approximation with a geometry-exact fallback.
- Smooth point interpolation, point-to-area reconstruction, and conservative
  cell-average remapping are three distinct operators. They can share topology,
  local coordinates, quadrature, derivative recovery, and sparse storage, but
  they should not share one semantic method name.
- Higher-order folded stencils are signed. `WeightRow` and the sparse assembly
  path should accept negative values. Coverage must be computed separately from
  the signed value operator.

The most promising sequence is:

1. implement geometry-exact planar/local-chart Sibson-0 point weights;
2. integrate those basis functions over spherical destination cells;
3. add fixed linear derivative recovery and Farin/Hiyoshi point methods;
4. implement second-order conservative remapping from spherical overlap first
   moments;
5. only then investigate tabulated regular-hex fast paths.

## What was checked

The investigation covered:

- `NaturalNeighbours.jl` v1.3.6 and its Sibson, Laplace, Sibson-1, Farin, and
  Hiyoshi implementations;
- `DelaunayTriangulation.jl` v1.6.6, including Bowyer-Watson insertion,
  Voronoi geometry, and weighted triangulations;
- ESMF patch and second-order conservative source code and documentation;
- the current `GlobalRegridding` sparse weight and tiled execution paths;
- exact lattice calculations and Julia experiments run in a persistent Julia
  daemon.

`NaturalNeighbours.jl` computes Sibson and Laplace coordinates by inserting the
query through a non-mutating Bowyer-Watson operation and inspecting the insertion
envelope. Sibson-1, Farin, and Hiyoshi use the same zeroth-order Sibson
coordinates and add source gradients or Hessians. Their larger data stencils
come from derivative generation rather than a larger natural-neighbour cavity.

`DelaunayTriangulation.jl` supports weighted triangulations and power diagrams,
but `NaturalNeighbours.jl` currently uses ordinary circumcentres in its
coordinate formulas. Weighted natural-neighbour coordinates are therefore not
already available by switching a triangulation type.

## Exact result on the planar hexagonal lattice

Use the unit triangular lattice

\[
p_{ij}=\left(i+\frac{j}{2},\frac{\sqrt 3}{2}j\right).
\]

The Voronoi cell of the origin is

\[
|x|\leq\frac12,\qquad
|x+\sqrt3y|\leq1,\qquad
|x-\sqrt3y|\leq1.
\]

By dihedral symmetry it is enough to study the wedge

\[
0\leq y\leq x/\sqrt3,\qquad 0<x<1/2.
\]

In that wedge, the ordinary positive support is the origin plus the east,
northeast, and southeast sites. A fifth northwest site enters on one side of
the Bowyer-Watson event curve

\[
x^2+y^2-\frac{2y}{\sqrt3}=0.
\]

Rotations and reflections produce twelve generic support regions: six with
four positive sites and six with five. Thus “one universal function over one
hexagon” needs a piecewise representation with explicit symmetry and event
regions. It is still small, but not one polynomial over the cell.

On the symmetry line `q=(t,0)`, `0 <= t <= 1/2`, the four nonzero Sibson
coordinates reduce to

\[
\lambda_O=\frac{(1-t)(3-t)}{3},\qquad
\lambda_E=\frac{t(2+t)}{3},\qquad
\lambda_U=\lambda_D=\frac{t(1-t)}{3}.
\]

These are non-negative, sum to one, and reproduce the query position. More
generally, the weights are piecewise rational functions: circumcentres are
rational functions of the query coordinates, and stolen Voronoi areas are
polygon-area expressions formed from those circumcentres. The event boundaries
are circumcircle arcs.

A structured scan over 6,561 query points and a random scan over 10,000 queries
found no positive second-ring weight. At four exact Voronoi vertices, numerical
bookkeeping can name a second-ring site with zero or roundoff-sized weight. That
is a degeneracy, not genuine support. The maximum generic support contained five
sites.

This validates the strong version of the reference-cell claim for an infinite,
perfect planar lattice. It does not validate the same claim for a distorted or
spherical lattice. An affine stretch preserves graph rings while changing
natural-neighbour weights; natural-neighbour coordinates are invariant under
similarities, not arbitrary affine transformations.

## Smoothness, exactness, and effective rings

The interpolation family separates into two layers:

1. query-local Sibson geometry;
2. derivatives estimated at the source sites and folded into the value weights.

If `B0`, `B1`, and `B2` evaluate the value, gradient, and Hessian parts of a
basis, and `G` and `H` are fixed linear derivative estimators, the final point
operator has the form

\[
W_{\mathrm{point}}=B_0+B_1G+B_2H.
\]

It is therefore one geometry-only sparse matrix. This stops being true when a
limiter, mask-dependent fit, or data-dependent stencil selection is introduced.

| Method | Nominal continuity | Polynomial reproduction with exact jets | Natural-neighbour geometry | Example effective support after derivative recovery |
|---|---|---|---|---|
| Sibson-0 | globally `C0`; `C1` away from data sites | affine fields | first ring | first ring |
| Sibson-1 | `C1` | affine and radial/spherical quadratics, not all `P2` | first ring | second ring with first-ring gradients |
| Farin | `C1` | all quadratics | first ring | second ring with first-ring gradients |
| Hiyoshi-2 | `C2` | all cubics | first ring | third ring with a compact quadratic jet fit; fourth with a direct cubic fit |

The exactness distinctions matter. Sibson-1 should not be advertised as a
general quadratic-reproducing method. Farin is the appropriate `C1`, `P2`
candidate; Hiyoshi is the `C2`, `P3` candidate.

Impulse experiments on a regular lattice gave representative folded stencils:

| Method and derivative estimator | Nonzeros | Furthest ring | Most negative weight |
|---|---:|---:|---:|
| Sibson-0 | 4 | 1 | 0 |
| Sibson-1, first-ring gradient | 14 | 2 | -0.0467 |
| Farin, first-ring gradient | 14 | 2 | -0.0435 |
| Hiyoshi, quadratic fit over two rings | 30 | 3 | -0.0149 |
| Hiyoshi, cubic fit over three rings | 52 | 4 | -0.0343 |

Counts vary with the query and accidental cancellation, but the reach and sign
patterns are stable. With `NaturalNeighbours.jl`'s direct cubic derivative
generation, even Sibson-1 and Farin can reach four rings.

Negative lobes are expected once derivatives are eliminated. They are what
allow higher-order polynomial reproduction from value-only degrees of freedom.
They also permit overshoot near fronts. The implementation should expose this
honestly and add a limiter later as a separate nonlinear execution mode.

## All area integration should be spherical

The plane has two possible roles, which must not be conflated:

1. it can define or approximate an interpolation basis in local coordinates;
2. it can provide coordinates in which a spherical integral is evaluated.

It should not replace the sphere as the production integration measure.

Given a nodal reconstruction

\[
f_h(x)=\sum_i\phi_i(x)f_i,
\]

the point-to-area weights for the actual spherical destination polygon `D_j`
are

\[
W_{ji}=\frac{1}{A_j}
       \int_{D_j\subset S^2}\phi_i(x)\,dA_{S^2}.
\]

There are two valid ways to define `phi_i`.

### Intrinsic spherical natural neighbours

Construct the spherical Voronoi diagram from unit-vector sites, insert the
query, and measure stolen **spherical** areas. The corresponding Delaunay
triangulation is obtained from the three-dimensional convex hull. This basis is
non-negative, nodal, and a partition of unity.

It does not generally satisfy either

\[
x=\sum_i\lambda_i p_i
\]

in ambient three-space or

\[
\sum_i\lambda_i\log_x(p_i)=0
\]

in the query tangent plane. Constants remain exact; planar affine precision is
only recovered asymptotically as cells shrink. The tangent moment residual is a
useful convergence diagnostic.

### Local-chart natural neighbours

Project the nearby source sites and query into a tangent chart, compute the
planar basis there, and integrate its pullback with spherical area:

\[
W_{ji}=\frac{1}{A_j}
       \int_{\Pi(D_j)}
       \phi_i^{\mathrm{plane}}(u)
       J_{\Pi^{-1}}(u)\,du.
\]

The basis is then projection-dependent even though the measure is spherical.
A Lambert azimuthal equal-area chart removes the area Jacobian from this
formula, provided the mapped domain is represented accurately. A gnomonic
chart maps great-circle edges to straight lines but requires a varying area
Jacobian. An azimuthal-equidistant chart makes radial distance exact but is not
area-preserving.

For angular radius `c=h/R`, the leading projection distortions are:

| Projection | Leading local distortion |
|---|---|
| Gnomonic | radial scale `1+c^2+...`; transverse scale `1+c^2/2+...`; area `1+3c^2/2+...` |
| Lambert azimuthal equal-area | exact area; directional scales `1 +/- c^2/8+...` |
| Azimuthal equidistant | exact radial scale; transverse scale `1+c^2/6+...` |

Thus a local projected basis can be a controlled `O(h^2/R^2)` approximation,
but spherical integration remains authoritative.

### Spherical fan quadrature

Triangulate each spherical destination polygon into a fan. For a fan triangle
with unit-vector vertices `a`, `b`, and `c`, one convenient radial
parameterization is

\[
r=(1-\xi-\eta)a+\xi b+\eta c,\qquad x=r/\|r\|,
\]

with surface Jacobian

\[
J(\xi,\eta)=R^2
\frac{|\det(a,b,c)|}{\|r\|^3}.
\]

Adaptive quadrature is appropriate because natural-neighbour support changes
across circumcircle events and the basis is piecewise rational. Positive
quadrature applied to Sibson-0 preserves non-negative integrated weights.

Planar integration remains useful for exact lattice analysis, reference-table
construction, and convergence tests. It should not define production DGGS cell
averages.

## Point-to-area reconstruction is not automatically conservative

Integrating a partition-of-unity basis gives destination row sums of one and
therefore exact constants. That is not the same as conserving source DGGS cell
mass.

The basis assigns each source node an induced mass

\[
m_i=\int_{S^2}\phi_i(x)\,dA.
\]

There is no reason for `m_i` to equal the area of the published DGGS cell around
that node, particularly when the reported centre is not the generator of that
cell. Point-to-area natural-neighbour reconstruction should therefore be named
and documented as a reconstruction operator, not a conservative cell-average
operator.

## Second-order conservative regridding

Second-order conservative regridding starts from source **cell means**, not
point samples. In a tangent coordinate system attached to source cell `S_i`,
reconstruct

\[
f_i(x)=\bar f_i+g_i\mathbin{\cdot}
       \left(u_i(x)-\bar u_i\right),
\qquad
\bar u_i=\frac{1}{A_i}\int_{S_i}u_i(x)\,dA_{S^2}.
\]

For overlap `O_ji = D_j intersect S_i`, compute spherical area and first moment

\[
A_{ji}=\int_{O_{ji}}dA_{S^2},\qquad
M_{ji}=\int_{O_{ji}}
       \left(u_i(x)-\bar u_i\right)dA_{S^2}.
\]

The destination mean is

\[
\bar f_j=\frac{1}{A_j}\sum_i
\left(A_{ji}\bar f_i+M_{ji}\mathbin{\cdot}g_i\right).
\]

If gradients are recovered by a fixed linear operator `g=G f`, the assembled
matrix is

\[
W_{\mathrm{FV}}=A_D^{-1}(A+M G).
\]

This is the direct extension of the current first-order overlap matrix. It is
conservative because the linear correction has zero mean over each source cell,
not because its final weights are positive. A third-order finite-volume method
adds Hessian recovery and spherical second moments.

ESMF follows this pattern: it recovers source gradients from vertex-sharing
source cells, computes overlap areas and centroids/first moments, emits signed
weights, and falls back to first order for degenerate gradient stencils. Its
second-order method has no monotonicity limiter and can overshoot. ESMF stores
non-negative source and destination coverage fractions separately from the
signed high-order value weights.

This method is more accurate and conservative, but it is not a globally `C1`
interpolant. The piecewise-linear reconstruction can jump at source-cell
boundaries. It solves a different problem from Farin or Hiyoshi interpolation.

## Sphere, non-Voronoi cells, and pentagons

### Published cells need not be centre Voronoi cells

Natural-neighbour point interpolation is defined by the point set. It does not
require the published DGGS polygons to equal the Voronoi cells of their reported
centres. Conservative remapping is defined by the published source and
destination polygons. Keeping these two geometries distinct is cleaner than
forcing one to approximate the other.

An additively weighted power diagram cannot generally make arbitrary published
polygons into the cells of fixed generators. A scalar site weight moves a
bisector but does not provide enough freedom to reproduce arbitrary edge
orientations and vertices. Power diagrams may improve a specially constructed
grid, but they are not a general repair for H3 or ISEA centre/polygon mismatch.

### The twelve pentagons

For a closed mostly hexagonal triangulation, the twelve valence-five defects are
topologically required:

\[
\sum_v(6-\operatorname{valence}(v))=12.
\]

They defeat a single translation-invariant stencil. They do not defeat natural
neighbours: the geometry-exact construction handles variable valence directly.

Loop subdivision provides a useful analogy, not a direct theorem about natural
neighbours. On the regular valence-six triangular mesh, Loop's limit surface is
the quartic three-direction box spline and is `C2`. At an extraordinary
valence-five vertex, tuned subdivision masks generally provide only `C1`
regularity. A Loop-like globally coupled reconstruction should therefore expect
special pentagon stencils and reduced regularity. Sibson/Farin/Hiyoshi schemes
must be analysed on their own; they do not automatically inherit Loop's `C1`
limit.

For acceleration, cache or tabulate by a geometry class such as

```text
(method, resolution, topology class, normalized local geometry, event region)
```

rather than only by normalized position in one regular hexagon. Regular
hexagons can use a reference fast path. Pentagons, seams, large distortion, and
near-Delaunay-flip configurations should use the exact kernel.

## Non-triangular and arbitrary grids

Natural-neighbour interpolation is a construction on sample sites, not on
source cell shapes. Any nondegenerate planar or spherical point set has a
Delaunay/Voronoi complex, so raster, pentagonal, and irregular grids can use it.
The supplied grid adjacency is an acceleration hint; the Delaunay insertion
envelope determines the semantic support.

A square lattice illustrates why topology alone is insufficient. Sibson weights
equal raster Q1 weights only in regions where the four positive natural
neighbours are exactly the rectangle corners. Elsewhere adjacent-square sites
enter; a query such as `(0.1, 0.8)` in a unit square has six natural neighbours.
On a stretched rectangular lattice the discrepancy is larger. Natural-neighbour
interpolation should therefore not silently replace the existing raster Q1
fast path.

There is also no positive, vertex-value-only “second-order barycentric” basis
that reproduces all quadratics at a generic interior point. If non-negative
weights satisfy constant and linear precision, strict convexity and Jensen's
inequality prevent exact reproduction of every quadratic. Higher exactness
requires at least one of:

- signed value-only weights;
- derivative degrees of freedom;
- additional edge/interior control values;
- a nonlinear or data-dependent reconstruction.

The general alternatives are:

| Construction | Strength | Cost or limitation |
|---|---|---|
| Q1 or generalized barycentric coordinates | local, positive, linearly exact | no general `P2` exactness |
| Quadratic polygonal finite elements | `P2` exact | edge/interior degrees of freedom; commonly signed basis functions |
| MLS/GMLS local polynomial recovery | topology-generic, fixed sparse `P2/P3` operator | signed weights; conditioning and boundary policy |
| Farin/Hiyoshi natural neighbours | interpolatory with `C1/C2` continuity | derivative recovery, larger signed stencil |
| Compact or global RBF | flexible smoothness and geometry | shape parameter, conditioning, potentially broad support |

For arbitrary grids, compact MLS/GMLS is the most direct interpretation of a
“second-order barycentric” value-only method. Natural neighbours are most
compelling where interpolation, locality, positivity at order zero, and smooth
Voronoi-based geometry are specifically desired.

## Comparison with existing methods

| Method | Source semantics | Exactness / smoothness | Positive | Conservative | Typical support |
|---|---|---|---|---|---|
| Nearest | point | constants; discontinuous | yes | no | 1 |
| Triangle barycentric | point | affine; elementwise `C0` | inside triangle | no | 3 |
| Raster Q1 | point | bilinear in raster chart; `C0` | inside cell | no | 4 |
| Sibson-0 | point | affine planar precision; `C0` globally | yes | no | natural neighbours |
| Farin | point plus recovered gradient | `C1`, `P2` | no after folding | no | natural neighbours composed with gradient stencil |
| Hiyoshi | point plus recovered jets | `C2`, `P3` | no after folding | no | larger composed stencil |
| ESMF patch | point/element recovery | quadratic least-squares patches blended by low-order shape functions | no | no | union of corner stars |
| RBF | point | kernel-dependent smoothness and reproduction | generally no | no | compact or global |
| First-order conservative | cell mean | piecewise constant | yes | yes | geometric overlaps |
| Second-order conservative | cell mean plus recovered gradient | higher-order finite volume; not globally smooth | no | yes | overlaps composed with gradient stencil |

ESMF patch is conceptually closer to MLS recovery than to natural-neighbour
coordinates. Both produce a larger signed pointwise matrix. ESMF patch fits
quadratic polynomials around element corners and blends them with low-order
finite-element shape functions; Farin instead blends value and derivative jets
using natural-neighbour geometry.

## Integration into `GlobalRegridding.jl`

The algorithms belong beside the existing nearest, barycentric, and
conservative methods because they need the same grid traits, lazy tile planning,
masks, and sparse execution machinery.

Suggested public method families are:

```julia
NaturalNeighborPoint(; scheme = :sibson,
                       geometry = :spherical,
                       derivatives = :compact,
                       boundary = :missing)

NaturalNeighborAverage(; scheme = :sibson,
                         geometry = :spherical,
                         quadrature = :adaptive,
                         boundary = :missing)

ConservativeSecondOrder(; gradient = :least_squares,
                          limiter = nothing,
                          fallback = :first_order)
```

Names should continue the conventions of the actual method API when implemented;
the important distinction is `Point`, reconstructed destination `Average`, and
conservative source-cell remapping.

An internal layout could be:

```text
GlobalRegridding/src/
  methods/natural_neighbor.jl
  methods/derivative_recovery.jl
  methods/spherical_basis_quadrature.jl
  methods/conservative_second_order.jl
```

No package boundary is required. If the first file eventually becomes a stable,
geometry-only implementation of spherical and planar natural-neighbour
coordinates with independent users, it can be extracted without changing the
`GlobalRegridding` method API.

### Sparse-weight changes

The current sparse builder rejects non-positive entries in
`lib/GlobalRegridding/src/methods.jl`. Higher-order methods require:

1. signed `WeightRow`/`WeightCOO` entries;
2. pruning by `abs(w) > tolerance`, not `w > 0`;
3. row-sum validation that tolerates cancellation;
4. a non-negative coverage operator stored separately from the signed value
   operator;
5. explicit missing-stencil semantics.

Coverage cannot be inferred by summing signed value weights. Renormalizing the
surviving part of a masked signed stencil is also unsafe: it destroys polynomial
reproduction and can amplify cancellation. Reasonable policies are:

- require the complete stencil selected at plan construction;
- freeze masks into the plan and rebuild weights when masks change;
- fall back to a lower-order complete stencil;
- use a nonlinear, mask-aware recovery path that is explicitly not one fixed
  matrix.

### Geometry and topology hooks

The implementation needs reusable internal operations for:

- ordered `k`-ring neighborhoods;
- local tangent bases and logarithmic/projection maps;
- spherical Delaunay/Voronoi insertion or a robust bounded local equivalent;
- spherical polygon intersection area, centroid, first moment, and later second
  moment;
- derivative recovery with rank and conditioning diagnostics;
- topology/geometry class identification for optional cached fast paths.

## Proposed implementation sequence

### Phase 0: contracts and reference tests

- Permit signed sparse value weights.
- Separate coverage from value weights throughout planning and execution.
- Add exact polynomial-reproduction tests in planar coordinates.
- Add spherical constant, symmetry, seam, and conservation tests.
- Record deterministic tie rules for cocircular Delaunay configurations.

### Phase 1: Sibson-0 point interpolation

- Prototype from `NaturalNeighbours.jl` on a local planar chart.
- Implement an intrinsic spherical oracle using spherical Voronoi insertion.
- Compare projection error on regular cells, pentagons, face seams, and distorted
  cells over resolution.
- Retain barycentric P1 and raster Q1 as existing specialized methods.

### Phase 2: spherical point-to-area reconstruction

- Integrate Sibson-0 basis functions over actual spherical destination polygons.
- Use positive adaptive spherical fan quadrature.
- Validate row sums, positivity, rotational invariance, and quadrature
  convergence.
- Document explicitly that this is not source-cell conservation.

### Phase 3: higher-order point methods

- Add a reusable fixed linear derivative-recovery operator.
- Implement Farin before Sibson-1 because Farin has full quadratic precision.
- Add Hiyoshi with selectable compact quadratic and cubic jet recovery.
- Measure effective ring reach, negative lobes, condition numbers, and error near
  fronts.
- Keep limiting as a later nonlinear mode.

### Phase 4: second-order conservative

- Extend spherical intersections to return first moments in each source tangent
  frame.
- Assemble `A_D^-1 (A + M G)`.
- Add first-order fallback for rank-deficient gradient stencils.
- Test global integral conservation independently of value row sums.
- Compare directly with ESMF on matched meshes.

### Phase 5: fast paths

- Tabulate regular planar-hex event regions and basis functions.
- Fit resolution- and geometry-class corrections for spherical regular cells.
- Use pentagon, seam, high-distortion, and near-flip cases as exact-kernel
  fallbacks.
- Benchmark build cost, matrix density, application cost, and accuracy separately.

## Acceptance and falsification tests

The universal-stencil optimization should be rejected or narrowed if any of the
following fail:

1. **Planar symmetry:** rotations and reflections do not map weights to the
   corresponding support permutation.
2. **Event continuity:** values or promised derivatives jump while crossing a
   Delaunay insertion event.
3. **Support containment:** a positive ideal-lattice Sibson-0 weight occurs
   beyond the first ring away from exact degeneracies.
4. **Spherical convergence:** reference-table error does not decay under
   refinement for regular hexagons.
5. **Defect behaviour:** pentagon and seam errors do not remain localized or do
   not converge under exact fallback.
6. **Polynomial reproduction:** folded Farin and Hiyoshi matrices fail `P2` and
   `P3` tests with exact or advertised recovered jets.
7. **Conservation:** the second-order conservative matrix fails the source-to-
   destination global integral test on complete spherical coverage.
8. **Chunk invariance:** eager and differently chunked plans do not produce the
   same stencil and result.

Near a cocircular Delaunay flip, two triangulations may be equally valid. The
coordinate result should be invariant in the limit or use a deterministic tie
policy. A topological first ring supplied by the DGGS must not silently override
the geometric natural-neighbour envelope.

## Sources

- R. Sibson, [“A vector identity for the Dirichlet tessellation”](https://doi.org/10.1017/S0305004100056589), 1980.
- G. Farin, [“Surfaces over Dirichlet tessellations”](https://doi.org/10.1016/0167-8396(90)90036-Q), 1990.
- H. Hiyoshi and K. Sugihara, [“Two generalizations of an interpolant based on Voronoi diagrams”](https://doi.org/10.11540/bjsiam.12.2_176), 2002.
- H. Hiyoshi, [“Stable computation of a natural neighbour interpolation”](https://doi.org/10.1504/IJCSE.2007.014460), 2007.
- J.-D. Boissonnat and F. Cazals, [“Smooth surface reconstruction via natural neighbour interpolation of distance functions”](https://doi.org/10.1016/S0925-7721(01)00018-9), 2002.
- J.-D. Boissonnat and J. Flötotto, [“A local coordinate system on a surface”](https://doi.org/10.1016/S0010-4485(03)00059-9), 2004.
- M. S. Floater, A. Gillette, and N. Sukumar, [“Gradient bounds for Wachspress coordinates on polytopes” and generalized barycentric survey](https://doi.org/10.1017/S0962492914000129), 2014.
- M. Rand, A. Gillette, and C. Bajaj, [“Quadratic serendipity finite elements on polygons using generalized barycentric coordinates”](https://doi.org/10.1090/S0025-5718-2014-02807-X), 2014.
- P. Lancaster and K. Salkauskas, [“Surfaces generated by moving least squares methods”](https://doi.org/10.1090/S0025-5718-1981-0616367-1), 1981.
- D. Mirzaei, R. Schaback, and M. Dehghan, [“On generalized moving least squares and diffuse derivatives”](https://doi.org/10.1093/imanum/drr030), 2012.
- C. T. Loop, [*Smooth subdivision surfaces based on triangles*](https://www.microsoft.com/en-us/research/publication/smooth-subdivision-surfaces-based-on-triangles/), 1987.
- U. Reif, [“A unified approach to subdivision algorithms near extraordinary vertices”](https://doi.org/10.1016/0167-8396(94)00007-F), 1995.
- P. W. Jones, [“First- and second-order conservative remapping schemes for grids in spherical coordinates”](https://doi.org/10.1175/1520-0493(1999)127%3C2204%3AFASOCR%3E2.0.CO%3B2), 1999.
- [`NaturalNeighbours.jl` interpolation mathematics](https://danielvandh.github.io/NaturalNeighbours.jl/stable/interpolation_math/).
- [`NaturalNeighbours.jl` Sibson implementation](https://github.com/DanielVandH/NaturalNeighbours.jl/blob/24be1a0a6d92f32dac5fed779cd2017eebd26eca/src/interpolation/coordinates/sibson.jl#L2-L29).
- [`NaturalNeighbours.jl` Bowyer-Watson envelope handling](https://github.com/DanielVandH/NaturalNeighbours.jl/blob/24be1a0a6d92f32dac5fed779cd2017eebd26eca/src/interpolation/utils.jl#L1-L58).
- [ESMF second-order conservative documentation](https://earthsystemmodeling.org/docs/nightly/develop/ESMF_refdoc/node5.html#sec:interpolation:conserve_2ndorder).
- [ESMF second-order conservative implementation](https://github.com/esmf-org/esmf/blob/1082792c4b565616f4ce056d3507d76ce0661df2/src/Infrastructure/Mesh/src/Regridding/ESMCI_Conserve2ndInterp.C#L1268-L1417).

