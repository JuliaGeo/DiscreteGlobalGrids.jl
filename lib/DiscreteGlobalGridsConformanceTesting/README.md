# DiscreteGlobalGridsConformanceTesting.jl

Check that a grid implementation follows the DiscreteGlobalGrids interfaces.

This companion to [DiscreteGlobalGrids.jl](../../README.md) provides reusable
property tests for grid authors. It checks relationships between operations,
such as whether a child's parent is the original cell, alongside geometry and
indexing rules. Keeping these tests in a separate package avoids adding Julia's
`Test` library to the main package's dependencies.

Use the suites when adding a grid system or changing its implementation:

- `test_grid_interface` checks a grid's cell identifiers, boundaries, centroids,
  and point location.
- `test_hierarchical_system` checks relationships across levels, spatial bounds,
  neighbors, and rings.
- `test_generic_fallbacks` checks the generic algorithms on your system's geometry,
  with its specialized methods hidden.

## Quick start

Use Julia 1.11 or later. Install the packages in the environment you use for tests:

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaGeo/DiscreteGlobalGrids.jl")
Pkg.add(url = "https://github.com/JuliaGeo/DiscreteGlobalGrids.jl",
        subdir = "lib/DiscreteGlobalGridsConformanceTesting")
```

This example checks an S2 grid and its hierarchy:

```julia
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGridsConformanceTesting

sys = DGG.S2System()
grid = DGG.levelgrid(sys, 2)

test_grid_interface(grid; n_samples = 16)
test_hierarchical_system(sys; levels = 0:3, n_samples = 8)
```

Replace `sys` with your own system. For a standalone grid without a hierarchy,
call only `test_grid_interface`. The example restricts the sampled levels to
`0:3`; omit `levels` to sample across the system's full supported range.

Run the generic algorithms separately on a small grid:

```julia
test_generic_fallbacks(sys; levels = 0:1, n_samples = 4)
```

Generic neighbor searches can cost more than specialized implementations, so
start with coarse levels. Each suite produces ordinary Julia test sets, with
failures grouped by the rule they violate.

## API

| API | Purpose |
| --- | --- |
| `test_grid_interface(grid; ...)` | Check the base grid interface. |
| `test_hierarchical_system(sys; ...)` | Check hierarchy, spatial bounds, and neighborhood rules. |
| `test_generic_fallbacks(sys; ...)` | Check generic implementations using the system's geometry. |

`n_samples` controls cell sampling. The system suites also accept `levels`
and `n_levels` to control level sampling. All three accept `rng` for reproducible
sampling and `label` for the test-set name. See the
[suite docstrings](src/DiscreteGlobalGridsConformanceTesting.jl) for tolerances
and other options.

## How it works

Each run samples cells with a seeded random generator. Checks in that run use
the same sampled cells, and repeating a default call reproduces the sample.
For a custom seed, load `Random` and pass `rng = Random.MersenneTwister(seed)`.
Use a fresh generator each time you want to reproduce a run.

The grid suite checks identifier round trips, boundary orientation, points on
the unit sphere, and centroids inside their cells. The hierarchy suite checks
parent–child relationships and whether spatial bounds contain descendant cells.
It also checks neighbor symmetry and ring ordering under the selected
connectivity rules.

By default, the suites skip optional operations without specialized methods.
They report the reasons at the end of the run; these skips appear as `Broken`
results in Julia's test summary. `test_generic_fallbacks` wraps the system to
hide its specialized methods and checks the generic implementations explicitly.

These suites check interface rules. Algorithms such as spatial queries and
regridding have their own tests in DiscreteGlobalGrids. See
[Writing a grid system](../../docs/src/extending.md) for the implementation guide.

## License

This package uses the repository's [MIT license](LICENSE.md).

## AI disclosure

This package was created with the help of AI agents, including Claude and Codex,
and will continue to be developed with these agents.
