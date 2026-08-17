# GlobalRegridding.jl

Regrid spherical cell collections eagerly or in chunks. The source and
destination are both `RegridSpace`s, so a regridding method is written once and
runs between any pair of them.

```julia
regrid(data; to, method = Conservative())   # build a plan, apply it, drop it
plan = plan_regrid(data; to, method)        # keep it
regrid(data, plan)                          # reuse across slices and reads
```

A dimensional source comes back on the destination's own axes — `RasterGrid`
gives `(X, Y)` — sampled as the method implies: `Intervals` for area-based
methods, `Points` for point samples, or whatever `sampling` says.

`RasterGrid` is the built-in space. `DiscreteGlobalGrids.jl` supplies
`DGGSpace` for its grid systems, and is the reference for what a space package
has to provide.

## Extension surface

A package that supplies its own space implements the `RegridSpace` interface —
`celltree`, `chunktree`, `nchunks`, `cellindices`, `ncells`, `getcell`,
`cellcentroid`, `cellat`, `hascellchart`, `manifold` — all exported.

Five further names are unexported but load-bearing from outside, and their
signatures are as fixed as the exported ones:

| Name | Where | Role |
|:--|:--|:--|
| `_asspace(target, name[, src_space])` | `src/api.jl` | resolve a `to`/`from` argument spelling into a `RegridSpace` |
| `subtree(space, inds)` | `src/conservative.jl` | cell tree restricted to one chunk |
| `chunkextents(space)` | `src/discovery.jl` | per-chunk spherical caps |
| `resolvespatialdims(data, nsrc)` | `src/executor.jl` | which array dimensions a regrid replaces |
| `dimsource(lookup)` | `src/spaces.jl` | the `from` a lookup already names |

Observability and finer extension points are marked `public` in
`src/GlobalRegridding.jl`.

## Tests

```
julia --project=lib/GlobalRegridding -e 'using Pkg; Pkg.test()'
```
