# DiscreteGlobalGridsConformanceTesting

Property-test suites for implementations of the
[`DiscreteGlobalGrids`](../..) interfaces. This is a separate test-only package
so `DiscreteGlobalGrids` does not depend on Julia's `Test` standard library.

```julia
using DiscreteGlobalGridsConformanceTesting

test_grid_interface(grid)
test_hierarchical_system(system)
```

Those two test what a system *implements*: a law group whose method a system
has not written is skipped rather than failed, and every skip states its reason
at the end of the run.

```julia
test_generic_fallbacks(system)
```

runs the same adjacency and point-location laws against the **generic
fallbacks** — the implementations that answer for a system which wrote none of
them — by presenting the system to dispatch with its fast paths hidden. It is
the expensive one, and it answers the other question: whether what will
actually answer obeys the contract on this system's geometry.

The package is a member of the repository's Julia workspace. Repository
tests resolve it through the path source declared in `test/Project.toml`.
