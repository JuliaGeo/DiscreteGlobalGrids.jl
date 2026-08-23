# Generic barycentric and patch-like point regridding

> **Implementation status:** This is a mathematical design note. The
> authoritative API, chunk-planning, Copernicus DEM, scope, and task decisions
> are in `regrid-notes/2026-08-23-barycentric-regridding-plan.md`. In
> particular, that plan supersedes the tentative `sampleelement`, `patchsites`,
> and `stencilreach` proposals below.

## Recommendation

Implement `BarycentricPoint` on the **dual interpolation complex** of the
source grid, and keep `PatchPoint` as a distinct higher-order method. A source
cell value is attached to that cell's centroid, so source centroids—not the
corners of source polygons—are the interpolation nodes. Each primal grid vertex
induces one dual element whose vertices are the centroids of its incident source
cells.

| primal source grid | typical vertex valence | dual interpolation element | basis |
|---|---:|---|---|
| quadrilateral raster | 4 | quadrilateral | tensor-product Q1 (bilinear) |
| hexagons/pentagons | 3 away from defects | triangle | P1 barycentric |
| triangles | usually 6 | polygon | generalized barycentric |
| mixed polygons | variable | polygon | generalized barycentric or a fixed triangulation |

This is the important meaning of “generic”: topology selects the local element;
geometry computes its weights. A query-local selection of arbitrary nearby
centroids is not enough, because independently selected triangles can disagree
across cell boundaries and make the interpolant discontinuous.

## Mathematical contract

For destination point `p` and source sample sites `xᵢ`, point regridding emits

```math
    \hat f(p) = \sum_i w_i(p) f_i, \qquad
    \sum_i w_i(p)=1, \qquad
    \sum_i w_i(p)x_i=p.
```

The last equality is interpreted in the chosen local 2-D coordinates. The first
condition reproduces constants; both together give linear precision. Inside a
convex element, non-negative weights additionally give a local maximum
principle. This is a **point sample at the destination centroid**, not a
destination-cell average, and it is not conservative.

- A triangle has unique P1 barycentric weights and reproduces affine fields.
- A polygon needs a named generalized-coordinate family. Mean-value
  coordinates are a reasonable default for convex, irregular dual polygons;
  the choice is not mathematically unique.
- A quadrilateral raster should use Q1 shape functions
  `((1-u)(1-v), u(1-v), uv, (1-u)v)`. On a rectangular native chart these
  exactly reproduce `a + bx + cy + dxy`: this is the promised bilinear case.
  On a curvilinear quadrilateral, solve the inverse isoparametric map for
  `(u,v)`; do not substitute distances to four points.

The current `BilinearPoint` instead requires global separable `chartaxes` and
therefore only describes `RasterGrid`. Its raster kernel is useful as the Q1
fast path, but the global chart trait is not the generic abstraction.

## Coordinates and the sphere

Search and containment should continue to use the package's unit-sphere
geometry. Weight construction needs a nonsingular local 2-D representation:

- For a raster, use its declared native chart. “Bilinear” is chart-dependent,
  so longitude periodicity and the branch across the antimeridian remain part
  of the raster contract.
- For an unstructured spherical element, map its sample sites to a tangent
  plane with the logarithmic map about `p` (or a fixed element centre), then
  compute P1/generalized barycentric weights there. Linear precision then means
  `Σwᵢ logₚ(xᵢ)=0`; there is no global affine coordinate system on a sphere.

Raw longitude/latitude must not be used for generic stencils: it is singular at
the poles, discontinuous at the antimeridian, and metrically anisotropic.
ESMF likewise performs spherical regridding in 3-D Cartesian space to avoid pole
and periodicity singularities. Local elements must be small enough to lie in a
common hemisphere; a folded tangent projection or non-unique inverse map is a
degeneracy, not a stencil.

## Discovering the correct neighbourhood

1. Locate the source cell containing `p` with `cellat`.
2. Locate the containing **dual element** among that cell's corners.
3. Return its incident source-cell positions in cyclic order and its element
   kind (`TriangleP1`, `QuadQ1`, or `PolygonP1`).
4. Compute weights and emit only entries present in the current `src_inds`, as
   required by `build_weights!`; other source chunks emit their shares.

`DGGSpace` already has centroids and counter-clockwise `Edge()`/`Vertex()`
neighbourhoods, but `RegridSpace` exposes no topology and an ordered one-ring
does not identify which cells are incident on each primal vertex when vertex
valence exceeds three. A truly generic implementation therefore needs a new
high-level hook rather than coordinate matching of polygon corners:

```julia
sampleelement(space, p) -> SampleElement(source_positions, kind, geometry)
patchsites(space, element, corner, degree) -> source_positions
stencilreach(space, method) -> conservative angular bound
```

`RasterGrid` can construct its four-node dual element arithmetically. Hexagonal
DGGS implementations can derive three-node elements from adjacent pairs in the
ordered edge ring. General spaces should provide native primal-vertex incidence
behind `sampleelement`. A deterministic, globally cached spherical Delaunay
triangulation of centroids is a possible fallback, but a different local
triangulation per query is not.

The method-facing API can remain small:

```julia
BarycentricPoint(; polygon=:meanvalue, boundary=:missing)
PatchPoint(; degree=2, rings=:adaptive, fallback=:linear,
             boundary=:missing, rcond=sqrt(eps(Float64)))

pointweights(method, src_space, p) -> (source_positions, weights)
```

`build_weights!` becomes generic over `pointweights`. `support_radius` should
delegate to a space-supplied, safe bound for the required dual element or
topological rings; underestimating it can silently omit cross-chunk weights.

## What “patch” means

ESMF patch recovery is not generalized barycentric interpolation. In 2-D it
constructs one quadratic polynomial for each corner of the source cell that
contains the destination point. Each polynomial is least-squares fitted to
source data in cells surrounding that corner; the evaluated corner patches are
then blended at the destination. ESMF supports this on structured grids and
unstructured polygons, reports a substantially larger sparse matrix than
bilinear, and explicitly does not promise range preservation.

For sites `qᵢ` in local coordinates, a quadratic patch uses

```math
 \phi(q)=[1,x,y,x^2,xy,y^2]^T,\qquad
 \min_a \sum_i \rho_i\left(\phi(q_i)^Ta-f_i\right)^2.
```

It is still a geometry-only sparse operator: compute

```math
 w^T=\phi(p)^T(A^TWA)^{-1}A^TW
```

with scaled coordinates and QR/SVD, never normal equations in implementation.
With a full-rank fit, it reproduces quadratics. At least six geometrically
independent sites are necessary; more are desirable for conditioning.

There are two honest APIs:

- **ESMF-style `PatchPoint`**: requires corner incidence, one cloud per corner,
  and the host element's blending basis. This can reuse `sampleelement` and
  `patchsites`.
- **Moving-least-squares patch**: fit one adaptive `Vertex()` k-ring around the
  host cell. This is easier and genuinely topology-generic, but should be called
  patch-like/MLS rather than claimed to reproduce ESMF's algorithm.

Expand rings until the design matrix has acceptable numerical rank. If that
fails, apply an explicit `quadratic → linear → nearest/missing` fallback. Patch
weights are generally signed, so the current missing-value renormalization is
not a meaningful coverage measure and destroys polynomial reproduction. The
safe default is to require the complete patch stencil; mask-aware patch weight
generation would need a separately planned geometry/mask operation.

## Boundaries, defects, and determinism

- The union of interior dual elements covers only the convex hull of sample
  sites. A regional grid therefore has a half-cell boundary strip with no
  barycentric stencil even though `cellat` succeeds. Choose explicitly among
  `:missing`, one-sided signed extrapolation, ghost nodes, or nearest clamping;
  do not silently clamp under the name barycentric.
- Reject or diagnose repeated sites, zero-area triangles, self-intersecting or
  non-convex dual elements unsupported by the selected coordinate family,
  failed Q1 inversion, tangent-plane folding, and rank-deficient patches.
- On a shared dual edge, use a deterministic element tie rule. Conforming P1/Q1
  elements give the same edge value, up to roundoff.
- Pentagon/icosahedral defects are ordinary variable-valence topology, not pole
  exceptions. Missing neighbours in partial grids are boundary cases.

## Vectors and tensors

The scalar sparse weights may be applied componentwise only when components use
one common fixed basis. East/north components live in different tangent frames.
Before interpolation, rotate or parallel-transport every source vector into the
destination tangent frame (or interpolate embedded 3-D vectors and project back).
Transport both slots of a rank-2 tensor, `T ↦ R T Rᵀ`. Otherwise even a constant
physical vector changes spuriously across the sphere. Transport is unique only
for local, non-antipodal stencils, which is another reason to bound support.

## Minimal acceptance laws

1. Identity at source sample sites and `Σw=1` on every valid stencil.
2. Exact local affine reproduction on triangles/polygons; exact
   `a+bx+cy+dxy` reproduction for raster Q1.
3. Hex/pentagon tests across face seams and raster tests across the periodic
   longitude seam.
4. Eager, chunked, and differently chunked plans produce the same weights;
   a stencil crossing a chunk boundary loses no term.
5. Explicit tests for each boundary/fallback policy and every degeneracy above.
6. Quadratic reproduction and rank-triggered fallback for `PatchPoint`, plus a
   test demonstrating that overshoot is permitted.

## Sources

- [ESMF Reference Manual: bilinear, higher-order patch, spherical and extrapolation behaviour](https://earthsystemmodeling.org/docs/release/ESMF_8_7_0/ESMF_refdoc/node5.html)
- [Current ESMPy overview: supported grid/mesh geometries and regridding in 3-D Cartesian space](https://earthsystemmodeling.org/esmpy/)
- [Floater, *Mean value coordinates* (2003)](https://doi.org/10.1016/S0167-8396(03)00002-5)
- [Floater, Gillette and Sukumar, *Generalized barycentric coordinates and applications* (2014)](https://doi.org/10.1017/S0962492914000129)
- [Zienkiewicz and Zhu, *The superconvergent patch recovery... Part 1* (1992)](https://doi.org/10.1002/nme.1620330702)
