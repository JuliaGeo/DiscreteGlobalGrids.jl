# Choosing the production regridding semantics

`scripts/copdem_production.jl`, `scripts/copdem_store.jl`, `scripts/launch_copdem_bench.sh`,
`benchmark/copdem_semantics.jl`, `docs/src/api/regridding-methods.md`.

## The configuration

`CONFIG.method`, from `COPDEM_METHOD`, names what the run computes. `regridmethod` gives the
method and `regridpolicy` the missing policy its weights need; the global plan in `dagplan`
and every column's plan in `regrid_chunk` take both, so the relation the schedule and the
cache read is built at the same `supportradius` the columns regrid with.
`scripts/launch_copdem_bench.sh point` launches it as the other three do.

| `method` | weights | policy | destination cell holds |
|---|---|---|---|
| `:conservative` (default) | area overlap | `Weighted(0.5)` | the coverage-normalised mean of the pixels it overlaps |
| `:point` | Q1 / barycentric between posts | `Weighted(1)` | a sample at its centroid |
| `:nearest`, `:nearest-direct` | one post | `Weighted(0.5)` | the nearest post's value |

An area row is a spectrum a threshold cuts. A point row is complete or empty — four posts or
none — so the threshold becomes a yes/no switch on whether an absent post may be
interpolated across, and `Weighted(1)` says no. It blanks nothing the loose setting keeps on
a source with no missing values: 24,210 cells placed under it and 24,210 under
`Weighted(0)`, on every twin arm. A new store records the name as `regridding_method` and
refuses to reopen under a different one, as a spherical store already refuses authalic
columns.

## What the two make of the same data

`benchmark/copdem_semantics.jl`, one process per table, `-t 1`, idle machine, cold.

### Real GLO-90 N46 E010 to IGeo7 level 12, 2,313,802 destination cells

| method | cold sweep | source chunk reads | placed | peak RSS |
|---|---:|---:|---:|---:|
| `Conservative()` | 36.049 s | 8 | 2,311,288 (99.891 %) | 1.69 GiB |
| `BarycentricPoint()` | 5.333 s | 8 | 2,305,690 (99.649 %) | 1.88 GiB |

Under `--heap-size-hint=1200M budget=134217728`, which keeps the arm inside the 2 GiB cap
and makes the point column 67.8 % collection; at the default heap the two are 36.087 s and
2.206 s, peaking at 2.33 GiB. So the ratio runs 6.8x to 16.4x with the heap the point path
is given, the values the same either way.

Both read the same eight chunks and agree to 3.02 m on average, 35.7 m at worst, over the
2,305,690 cells both placed, on Alpine relief. The 5,598 the area method placed alone,
0.24 %, are the half-pixel strip outside the outermost posts; none went the other way.

### Analytic fields on the scaled lattice twin

Thirty latitude intervals to the degree, every band ratio of a real product; the holding
10-12 E, 49-51 N puts the 50-degree band edge and its transitions inside it. The destination
is IGeo7 level 9, 1,124 m cells against 2,997 m posts, finer than the source, which is where
the two readings part company. The `rim` arm pushes it a degree past the holding.

| arm | field | area error | point error | both placed | area only | mean \|diff\| |
|---|---|---:|---:|---:|---:|---:|
| linear | `1000 + 37λ − 61φ` m | 0.423 m mean, 2.123 max | 4.2e-13 m mean, 1.8e-12 max | 24,210 of 25,553 | 903 | 0.4168 m |
| quadratic | `1000 + 90(λ−11)² + 140(φ−50)²` m | 1.094 m mean, 8.393 max | 0.0528 m mean, 0.0951 max | 24,210 of 25,553 | 903 | 1.0392 m |
| rim | linear | | | 24,210 of 101,381 | 1,334 | 0.4168 m |

Q1 and barycentric coordinates reproduce an affine field exactly, so the point column is
roundoff. The area column is not wrong: an area mean of an affine field is that field at the
*area's* centroid — of the posts covering a cell, not the cell's own — half a source pixel
away, 1.0 m at 61 m per degree. The 1,334 rim cells it places alone, 5.22 % of its own, are
stencils naming a post the holding does not hold; nothing is substituted and nothing is
renormalised.

### The poles

The twin's polemost post row stands at 89.991667 N, 926 m short of the pole. The source is
that whole band of tiles, all 360 degrees of longitude, since the policy takes the nearest
site of that row and a holding short in longitude would put its own rim in the way; the
destination is the 792 IGeo7 level-12 cells whose centroids fall north of it.

| method | cold sweep | reads | placed | error vs the field |
|---|---:|---:|---:|---:|
| `Conservative()` | 0.457 s | 1,892 | 792 | mean 23.04 m, max 98.12 m |
| `BarycentricPoint()` | 0.004 s | 416 | 792 | mean 23.03 m, max 107.57 m |
| `poles = nothing` | | | 54 | mean 0.00022 m, max 0.00046 m |

The field is a ramp in the plane tangent at the pole; a field affine in longitude would
report its own discontinuity there. Over the 54 cells that have a real dual cell the point
method reproduces that ramp to half a millimetre; over the other 738 both methods
effectively hold the nearest post of a row 926 m away and land within 0.31 m of each other
on average, 21.7 m at worst. The policy is 93.2 % of that strip.

## No regression to the conservative run

`scripts/copdem_production.jl` on a bounded subset — four level-5 chunks of a 2x2-degree
synthetic tile block at level 12, `maxchunks = 4`, one worker, `budget = 2^26`,
`schedule = :canonical` — against the same script without the option; fresh alternating
processes, `-t 1`, idle machine.

| rep | without, cells/s | with, cells/s | without, total | with, total |
|---|---:|---:|---:|---:|
| 1 | 179,916 | 207,043 | 31.614 s | 29.326 s |
| 2 | 212,696 | 209,721 | 29.785 s | 28.903 s |
| 3 | 211,557 | 207,603 | 28.546 s | 29.042 s |
| median | 211,557 | 207,603 | 29.785 s | 29.042 s |

**0.98x on the median regrid rate and on the median total wall**, inside the 15 % spread the
unmodified script shows across its own three reps, and the values are identical: 3,294,172
cells compared between the two stores, 2,221,171 nodata in both, **0 differing**. Peak RSS
1.75-2.20 GiB against a 1.56 GiB single-chunk floor, most of it the package and the store,
not the regrid, which the byte budget bounds.

## Promoted to `public`

A source space can now be made point-samplable from outside the package, nothing exported to
do it: `hasdualcells` declares that the space answers with dual cells, `samplerstate`
prepares what the lookup wants ready once, `dualcellat` answers one point, `Sampler` and
`DualCell` are what those dispatch on and return, `BasisKind`/`Bilinear`/`MeanValue` are a
`DualCell`'s third argument, and `chartat` names the plane the nodes are in for a space with
no cell chart of its own. `WeightRow`, `weightsat!`, `sampler`, `samplesites` and the
coordinate kernels stay internal: replacing the whole query is the Copernicus DEM route, not
yet a contract to rely on outside. `BarycentricPoint` joins the method names
`DiscreteGlobalGrids` re-exports, beside `Conservative`, `NearestCell` and `DirectNearest`.
