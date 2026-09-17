# Local documentation validation

This file records the validation added for DOC-020 and its limits. The fast
build is local, uses the existing documentation environment, and never deploys.

## Fast build

Run:

```sh
DGG_DOCS_FAST=true julia --project=docs docs/make.jl
```

The fast mode:

- copies `docs/src` to a temporary directory;
- generates all ten Literate tutorial pages in that directory;
- builds every page in the navigation tree;
- executes examples on the API and internals pages;
- treats Documenter warnings as fatal with `warnonly=false`;
- disables remote source links and deployment;
- uses Documenter's plain HTML renderer, so it does not load Vitepress, Bonito,
  GeoMakie, Makie, or WGLMakie.

The mode logs every page whose examples it does not execute. The list is
`index.md`, `all_dggs.md`, `extending.md`, and the ten generated tutorial
pages. Those pages need graphics, downloads, optional development packages, or
large datasets. The staged Markdown keeps their code as ordinary `julia`
fences, so the published source is unchanged.

Literate's `execute=false` controls conversion only. With
`Literate.DocumenterFlavor()`, conversion still creates `@example` blocks,
and Documenter executes them later. The fast mode changes those blocks only in
its temporary copy.

## Checks run

The local docs environment loaded DiscreteGlobalGrids 0.1.0, Documenter 1.19.0,
and Literate 2.21.0 under Julia 1.13.0. These checks passed on 2026-09-15:

- package load in the repository environment;
- the package docstring doctest, 1 of 1;
- strict fast build from a fresh documentation daemon;
- direct execution of the focused examples for grid geometry and hierarchy,
  `Values` and `NeighborSlices`, eager regridding and plan reuse, and
  partitioning a synthetic problem and a lazy regridding plan;
- controlled failure fixtures for a broken `@ref` and a throwing `@example`.

Both fixtures made `makedocs` fail with the expected
`:cross_references` or `:example_block` category. The fixture harness caught
those expected failures and printed `fatal-docs-failure-fixtures-ok`.

The strict build initially exposed invalid references in rewritten pages and in
rendered docstrings. The corrections qualify public bindings where Documenter
can resolve them and leave unrendered implementation helpers as code literals.
The final fresh-daemon run exited successfully.

For agents following the Julia daemon workflow, the equivalent command used
for the final build was:

```sh
jld --project=docs --name=sweep_core_fastdocs --idle-timeout=2h eval \
  'ENV["DGG_DOCS_FAST"] = "true"; include(joinpath(pwd(), "docs", "make.jl"))'
```

`Documenter.doctest(DiscreteGlobalGrids)` must run in the daemon's `Main`
module with the current manual pages. A scratch module cannot resolve the
manual pages' `CurrentModule = DiscreteGlobalGrids` metadata even when the
package binding exists in `Main`.

## Limits

The fast mode does not execute the landing page, gallery, extension guide, or
tutorial examples. It does not test external downloads, large datasets,
graphics, Vitepress or Bonito output, GL setup, or deployment. The existing CI
documentation job remains the full-site check and provides Xvfb, graphics
libraries, automatic DataDeps acceptance, and extra swap for the hydrology
tutorial.

`checkdocs=:none` remains explicit in both modes. The strict fast build checks
the pages, references, rendered docstrings, doctests, and active examples, but
does not require every exported binding to appear in an `@docs` block.
