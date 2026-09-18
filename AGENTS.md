# Working in this repository

## Environments

Manifests are disposable. Every environment in this monorepo is fully described by its `Project.toml`: workspace members under the root `[workspace]`, path dependencies under `[sources]`, extensions under `[extensions]` and `[weakdeps]`. `Manifest.toml` files are gitignored and must never be copied between checkouts or worktrees.

- A fresh checkout or worktree starts with `julia --project=<env> -e 'using Pkg; Pkg.instantiate()'`.
- When a package or extension fails to load, or a path dependency is missing, delete the environment's `Manifest.toml` and instantiate again before debugging anything else. A stale manifest silently hides extensions added to `Project.toml` after it was written.
- To add a dependency, edit `Project.toml` (`[deps]`, and `[sources]` for a path) and re-resolve. Never hand-edit a manifest.

## Tests and docs

- One test file: `julia --project=test test/<path>.jl` from the repo root. Each suite is its own module.
- GlobalRegridding: `julia --project=lib/GlobalRegridding -e 'using Pkg; Pkg.test()'`, which is what CI runs. `lib/CopernicusUtils` runs the same way.
- Full suite: `julia --project=. -e 'using Pkg; Pkg.test()'`, about 18 minutes.
- Strict docs build: `DGG_DOCS_FAST=true julia --project=docs docs/make.jl`.

## Public names

A new public name is exported or declared `public` in `src/DiscreteGlobalGrids.jl`, added to the literal list in the module-surface testset of `test/interface/runtests.jl`, given a docstring, and listed in a `@docs` block under `docs/src/api/`. An unlisted docstring is silently unpublished.
