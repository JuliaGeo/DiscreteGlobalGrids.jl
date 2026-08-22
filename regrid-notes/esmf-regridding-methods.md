# ESMF regridding methods: concise geometry note

Scope: ESMF/ESMPy 8.9.1, the current release as of 2026-08-21. ESMF first builds a sparse operator

\[
u^{d}_j=\sum_i W_{ji}u^{s}_i,
\]

then applies it by sparse matrix multiplication. The distinctions below are therefore distinctions in how a containing cell or neighborhood is found and how each row of `W` is formed. See the [ESMF 8.9.1 release table](https://earthsystemmodeling.org/static/releases.html), [regridding overview/status](https://earthsystemmodeling.org/regrid/), and [Reference Manual, Regrid section](https://earthsystemmodeling.org/docs/release/latest/ESMF_refdoc/node5.html).

## Bilinear: local finite-element interpolation

**Documented/source facts.** ESMF locates the source element containing each destination point, maps the point to the element's reference coordinates, and evaluates the element's nodal shape functions. Its source implements

- triangle: \(N=(1-\xi-\eta,\xi,\eta)\), i.e. barycentric interpolation;
- quadrilateral: \(N=\tfrac14((1-\xi)(1-\eta),(1+\xi)(1-\eta),(1+\xi)(1+\eta),(1-\xi)(1+\eta))\), i.e. isoparametric bilinear interpolation.

These formulas are explicit in [`ESMCI_ShapeFunc.C`](https://github.com/esmf-org/esmf/blob/v8.9.1/src/Infrastructure/Mesh/src/Regridding/ESMCI_ShapeFunc.C#L179-L196) and [`ESMCI_ShapeFunc.C`](https://github.com/esmf-org/esmf/blob/v8.9.1/src/Infrastructure/Mesh/src/Regridding/ESMCI_ShapeFunc.C#L263-L309); weight extraction evaluates the containing element's field and records sensitivities to its local degrees of freedom ([`ESMCI_Interp.C`](https://github.com/esmf-org/esmf/blob/v8.9.1/src/Infrastructure/Mesh/src/Regridding/ESMCI_Interp.C#L2701-L2782)). Thus a rectilinear raster gives the usual four weights \((1-t)(1-s),t(1-s),ts,(1-t)s\).

In 2-D ESMF accepts polygons with any number of sides, but polygons with more than four sides are represented internally by an ear-clipped triangle set ([ESMF status page](https://earthsystemmodeling.org/regrid/), [`ESMCI_Mesh_Glue.C`](https://github.com/esmf-org/esmf/blob/v8.9.1/src/Infrastructure/Mesh/src/ESMCI_Mesh_Glue.C#L1132-L1191), [`ESMCI_MathUtil.C`](https://github.com/esmf-org/esmf/blob/v8.9.1/src/Infrastructure/Mesh/src/ESMCI_MathUtil.C#L1202-L1320)). For an element-centred `Mesh` field, integrated regridding first creates a dual mesh: original cell centres become dual nodes, and the ordered cells incident on an original vertex form a dual element ([`ESMF_FieldRegrid.F90`](https://github.com/esmf-org/esmf/blob/v8.9.1/src/Infrastructure/Field/src/ESMF_FieldRegrid.F90#L2653-L2682), [`ESMCI_MeshDual.C`](https://github.com/esmf-org/esmf/blob/v8.9.1/src/Infrastructure/Mesh/src/ESMCI_MeshDual.C#L394-L492)).

**Inference for a generic method.** ESMF's “bilinear” is best understood as *lowest-order, element-local Lagrange interpolation*, not as one four-point raster formula. A six-sided source polygon with vertex data becomes piecewise barycentric over its chosen triangulation; cell-centred hexagon data is instead interpolated on the dual mesh of neighbouring cell centres. The former is triangulation-dependent and is not a generalized n-gon barycentric scheme such as mean-value or Wachspress coordinates. A new generic API should make that policy explicit.

## Patch: quadratic recovery over a topological star

**Documented/source facts.** After locating the containing source element, ESMF constructs one polynomial patch for each of its corner nodes. For corner \(a\), its neighborhood is the topological star of elements incident on \(a\). A least-squares polynomial is fitted to field samples from those elements; masked incident elements are omitted. In 2-D the documented patch degree is 2. The source uses the nine-term tensor-product basis \(x^m y^n\), \(0\le m,n\le2\), samples the finite-element field at quadrature points, and solves the least-squares system with `DGELSD` ([manual algorithm](https://earthsystemmodeling.org/docs/release/latest/ESMF_refdoc/node3.html), [`ESMCI_PatchRecovery.C`](https://github.com/esmf-org/esmf/blob/v8.9.1/src/Infrastructure/Mesh/src/Regridding/ESMCI_PatchRecovery.C#L59-L99), [`ESMCI_PatchRecovery.C`](https://github.com/esmf-org/esmf/blob/v8.9.1/src/Infrastructure/Mesh/src/Regridding/ESMCI_PatchRecovery.C#L677-L874)).

If \(q_a\) is the recovered patch and \(N_a\) is the same lowest-order element shape function used by bilinear interpolation, ESMF evaluates

\[
u^d(x)=\sum_{a\in\text{corners}(K)}N_a(x)q_a(x).
\]

The source calls these blending functions a “bilinear partition of unity” ([`ESMCI_PatchRecovery.C`](https://github.com/esmf-org/esmf/blob/v8.9.1/src/Infrastructure/Mesh/src/Regridding/ESMCI_PatchRecovery.C#L1483-L1518)). The stencil therefore contains nodes from all elements incident on the containing element's corners. ESMF reports roughly four times the bilinear matrix size for quadrilateral grids; patch is smoother but neither monotone nor range-preserving ([Reference Manual](https://earthsystemmodeling.org/docs/release/latest/ESMF_refdoc/node5.html)).

**Inference for a generic method.** Patch recovery is a second layer on top of a valid local-coordinate/shape-function method. It requires topology (incident-element stars), a well-conditioned local tangent coordinate system on the sphere, enough independent samples, and an explicit boundary/degeneracy fallback. It is not merely “use more nearest neighbours”.

## Conservative and nearest methods

**Documented facts.** First-order conservative regridding treats each source cell value as constant and uses intersection areas

\[
W_{ji}=\frac{A(S_i\cap D_j)}{A(D_j)}
\]

for destination-area (`DSTAREA`) normalization. With destination covered fraction \(f_j\), fraction-area normalization uses \(A(S_i\cap D_j)/(A(D_j)f_j)\). The former preserves the covered integral but produces fraction-diluted values in partially covered destination cells; the latter produces local averages, with \(f_j\) still needed in integral calculations. Second-order conservative regridding adds a reconstructed source-gradient contribution and may overshoot. See [first/second-order and normalization sections](https://earthsystemmodeling.org/docs/release/latest/ESMF_refdoc/node5.html) and [ESMPy `NormType`](https://earthsystemmodeling.org/esmpy_doc/release/latest/html/NormType.html).

`NEAREST_STOD` gives each destination row one unit-weight entry for its nearest source. `NEAREST_DTOS` assigns every source to its nearest destination; multiple sources assigned to one destination are summed, and some destinations may remain empty. Equal-distance ties use the smallest sequence index ([Reference Manual](https://earthsystemmodeling.org/docs/release/latest/ESMF_refdoc/node5.html)).

## Geometry and supported inputs

| Method | 2-D source geometry | Neighborhood/support | Key limits |
|---|---|---|---|
| Bilinear | multi-tile `Grid`; polygonal `Mesh`; `XGrid` | containing triangle/quad after any internal triangulation | no self-intersecting cells; `LocStream` is destination-only |
| Patch | single-tile `Grid`; polygonal `Mesh`; `XGrid` | union of incident-element stars at containing-cell corners | 2-D only; larger stencil; may overshoot |
| Conservative | multi-tile `Grid`; polygonal `Mesh`; `XGrid` | all cells with non-zero geometric overlap | needs cell boundaries and cell-centred/element data; no `LocStream` |
| Nearest | `Grid`, `Mesh`, `LocStream`, `XGrid` | nearest point(s), no containing cell needed | topology is irrelevant except for later creep fill |

These are the current documented capabilities ([bilinear/patch/conservative sections](https://earthsystemmodeling.org/docs/release/latest/ESMF_refdoc/node5.html)). “Hexagon” here means a 2-D six-sided polygon. A 3-D *hexahedron* is a different supported element: bilinear/trilinear works on it, while 3-D patch is unsupported.

On the sphere ESMF works in 3-D Cartesian space to avoid longitude periodicity and the pole singularity. Bilinear and patch may use Cartesian chord or great-circle line types (Cartesian is default); conservative methods use great-circle cell edges. Edges spanning half a sphere or more are unsupported/ambiguous, and long edges may not represent straight longitude-latitude boundaries ([spherical line-type discussion](https://earthsystemmodeling.org/docs/release/latest/ESMF_refdoc/node5.html)).

## Masks, unmapped points, and extrapolation

**Documented facts.** For bilinear/patch, masks live on points: masking one source corner invalidates every source cell using it, while a masked destination point gets no row. Conservative masks are cell masks. Extrapolation weights are appended only for unmapped destinations and are not supported for conservative regridding because they would break conservation. Options are nearest source, inverse-distance average over the nearest \(N\) sources (weights proportional to \(d^{-P}\) and row-normalized; ESMPy defaults \(N=8,P=2\)), and destination-topology creep fill. Creep fill averages already-filled immediate neighbours level by level and may leave points unmapped. The default unmapped action is an error; `IGNORE` omits those rows. Detection is unavailable for `NEAREST_DTOS`, which behaves as ignore ([mask/extrapolation sections](https://earthsystemmodeling.org/docs/release/latest/ESMF_refdoc/node5.html), [ESMPy `Regrid`](https://earthsystemmodeling.org/esmpy_doc/release/latest/html/regrid.html), [ESMPy `ExtrapMethod`](https://earthsystemmodeling.org/esmpy_doc/release/latest/html/ExtrapMethod.html)).

## Practical implication

**Inference.** A generic neighborhood regridder should separate four contracts: (1) containment and local coordinates, (2) source data location and any primal-to-dual conversion, (3) the interpolation kernel—triangle barycentric, quad bilinear, explicit n-gon coordinates, or patch recovery—and (4) unmapped/extrapolation policy. Calling all of these “barycentric” would hide materially different geometry, especially for cell-centred hexagonal grids.
