# K1 — the barycentric kernels

- Date: 2026-08-25
- Plan: `regrid-notes/2026-08-23-barycentric-regridding-plan.md`
- Scope: `lib/GlobalRegridding/src/barycentric.jl`, the `BarycentricPoint`
  docstring in `methods.jl`, `lib/GlobalRegridding/test/test_barycentric.jl`

## 1. Two bases, because there are two interpolants

`BasisKind` names `Bilinear` and `MeanValue`. They are different interpolants,
not two spellings of one: `Bilinear` inverts the isoparametric map of the unit
square onto a quadrilateral and takes the tensor Q1 weights of the recovered
coordinates, which reproduce `a + bx + cy + dxy` on an axis-aligned rectangle;
`MeanValue` weights a convex polygon of any node count by mean-value
coordinates, which reproduce affine fields and nothing more.

The third kind, `Linear`, is gone. Mean-value coordinates on a triangle *are*
that triangle's barycentric coordinates — the identity is exact, not an
approximation — so a triangle is a `MeanValue` cell and the separate P1 kernel
was one implementation of a formula the polygon kernel already carried.
`NO_DUALCELL` says `MeanValue`; `dualweights!` sends `Bilinear` to
`bilinearweights!` and everything else to `meanvalueweights!`; the kind is still
the cell's own statement and never a reading of its node count.

## 2. Mean value from distances

Write `sᵢ` for the vector from the query point to node `i`, `rᵢ` for its length
and `ŝᵢ` for the unit vector along it. Each edge `(i, j)` contributes

```math
t = \operatorname{sign}\det(sᵢ, sⱼ)\,\frac{\|ŝᵢ - ŝⱼ\|}{\|ŝᵢ + ŝⱼ\|},
```

the tangent of half the angle the edge subtends at the point, and node `i` takes
`(t₋ + t₊) / rᵢ`, normalized. The two norms are `2 sin(θ/2)` and `2 cos(θ/2)`,
so the ratio is the half-angle tangent with no `atan`, `tan`, `sin` or `cos`
evaluated anywhere: distances and one determinant sign. A ring's turning
direction multiplies the determinant, so a clockwise cell reads the same as a
counter-clockwise one and the positivity guard still means what it says.

Every rule around the formula is unchanged: a point at a node takes that node
alone; a point on an edge takes the edge's two nodes by their distances along
it; a point outside is `WeightsOutside` and is never clamped in; a ring with a
repeated node, a straight corner or a reflex corner is rejected by
`_orientation` before any weight is formed; a non-finite or non-positive weight
is `WeightsDegenerate`. The kernel allocates nothing once the row is warm.

Steady-state cost of one call, one session, one benchmark script with the
kernel source swapped between arms:

| half-angle tangent | triangle | pentagon | hexagon |
|---|---:|---:|---:|
| `tan(atan(det, dot)/2)`, before this card | 100.0 ns | 126.2 ns | 145.0 ns |
| `det / (rᵢrⱼ + sᵢ·sⱼ)`, GeometryOps' `t_value` | 23.8 ns | 34.9 ns | 40.1 ns |
| `±‖ŝᵢ - ŝⱼ‖ / ‖ŝᵢ + ŝⱼ‖`, shipped | 31.1 ns | 46.3 ns | 55.7 ns |

## 3. What it verifies

`lib/GlobalRegridding/test/test_barycentric.jl`:

- a triangle read as a mean-value cell against the closed form the test
  computes from signed areas, at interior points, at each node, on an edge and
  outside, to `1e-12` — the assertion that lets the P1 kernel go;
- a sweep of convex rings of three to nine nodes, each at three interior
  points: every query mapped with one entry per node, every weight positive,
  partition of unity and linear reproduction to `1e-12`, and the same weights
  from the reversed ring — which is what fails if the turning direction is
  dropped from the determinant;
- the mean-value cases already pinned: a regular hexagon weighting its nodes
  equally at the centre, an irregular pentagon reproducing affine fields, nodes
  and edge points emitting only the nodes that carry them, and every degeneracy
  class rejected;
- the `Bilinear` cases, untouched but for the kind list.

Against the kernel as it stood before the rewrite, on 3,920 generated cases —
convex polygons of three to nine nodes, half of them clockwise, queried inside,
at a node, on an edge, outside, and at `1e-2` and `1e-6` of the way from an edge
to the centroid — every status and every index row is identical and the largest
weight difference is `6.1e-16`. On a further 45 cases chosen to be hostile —
every degeneracy class, corners straight to within `1e-14` … `1e-3`, nodes
separated by `0` … `1e-6`, and points `0` … `1e-4` either side of an edge —
every status and index row is identical again; weights agree to `5.6e-17`
except on razor triangles, where a cell of aspect ratio `2.5e-12` differs by
`6.1e-6`, `2.5e-10` by `8.2e-8` and `2.5e-7` by `9.6e-11`.

## 4. Why not the determinant-over-dot form

GeometryOps computes the same tangent as
`det(sᵢ, sⱼ) / (rᵢrⱼ + sᵢ·sⱼ)`. It is algebraically equal and about 28 % faster
per call, and it was measured and rejected: as a point approaches an edge the
subtended angle approaches `π`, and that denominator is then a difference of
nearly equal numbers.

On the same fixtures, against 256-bit truth on a triangle:

| perpendicular offset, as a fraction of the cell's size | `det`/dot error | shipped error |
|---|---:|---:|
| `1e-2` | `5.0e-15` | `5.6e-17` |
| `1e-4` | `3.0e-13` | `5.6e-17` |
| `1e-6` | `1.1e-11` | `5.6e-17` |
| `1e-8` | `1.5e-9` | `5.6e-17` |
| `1e-10` | denominator `0`, `NaN` | `5.6e-17` |

Over the 3,920-case capture the `det`/dot form reaches `7.6e-9`; in all 83
triangle cases where the two forms differed by more than `1e-13` it was the
further from truth. Below an offset of about `1e-9` the denominator cancels to
exactly zero and the kernel answers `WeightsDegenerate`: four of the 45 hostile
cases changed status that way, where the shipped form changes none. The
containment test absorbs offsets below `1e-12`, so the band that would have
changed is roughly `1e-12` to `1e-9` of a cell's size — reachable when P4 builds
dual triangles for a hexagonal source, where the cost is a destination dropped
by the missing policy rather than interpolated.
