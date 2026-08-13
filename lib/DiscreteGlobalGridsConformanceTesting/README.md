# DiscreteGlobalGridsConformanceTesting

Property-test suites for implementations of the
[`DiscreteGlobalGrids`](../..) interfaces. This is a separate test-only package
so `DiscreteGlobalGrids` does not depend on Julia's `Test` standard library.

```julia
using DiscreteGlobalGridsConformanceTesting

test_grid_interface(grid)
test_hierarchical_system(system)
```

The package is a member of the repository's Julia workspace. Repository
tests resolve it through the path source declared in `test/Project.toml`.
