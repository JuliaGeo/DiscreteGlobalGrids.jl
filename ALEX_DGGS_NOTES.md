# AGENTS.md

Discrete global grid systems (DGGS) for Julia Geo — 14 registered systems, six
with cell geometry. This file records the non-obvious facts; the README and
in-file docstrings are the detailed reference.

## Setup (fails without this)

- Requires **Julia ≥ 1.11**.
- `ConservativeRegridding` resolves via a `[sources]` path entry
  (`{rev=..., url=...}`) pointing at a **sibling checkout**
  `../ConservativeRegridding.jl` (its `Trees` submodule is not in a release).
  That checkout must exist before `instantiate`/`test` — CI clones it explicitly.
- Work in the package environment: `julia --project=.` (a `[workspace]` spans
  `test/` and `docs/`, so those are separate environments).

## Commands

- Full test suite: `julia --project=. -e 'using Pkg; Pkg.test()'` — ~528k
  assertions, ~2m20s warm.
- Single suite (fast iteration): `julia --project=. test/IGeo7/runtests.jl`
  (same for `test/core`, `test/H3`, `test/DGGSZarr`, …). Each suite is its own
  module.
- `test/DGGSZarr` has two tiers: `test_monotonic.jl` is self-contained, while
  `test_archives.jl` reads archives that live **outside this repo** and skips
  itself (loudly) when they are absent. Point it elsewhere with
  `DGGS_ZARR_TEST_DATA=<dir>`; the default is
  `~/dev/build/igeo7_z7_xarray_paper/data/working`.
- Examples (assertion-checked, exit non-zero on failure):
  `julia -t 4 --project=. examples/<name>.jl`. One exception:
  `examples/z7_zarr_plot.jl` needs CairoMakie/GeoMakie, which are **not**
  package deps — run it with `--project=docs`, which already has them.
  (`examples/healpix_demo.jl` is broken on `main`, unrelated to any of this.)
- Docs: `julia --project=docs docs/make.jl` (needs Xvfb/GL on CI — Makie).

## Architecture

- **No per-system tree code.** A system supplies operations only — child
  enumeration, parent, descendant interval, unit-sphere boundary — and the
  package supplies grid types, traversal, pruning, caps. Adding a system is one
  `<X>Kernel.jl` file plus its `include` line; unwired ops throw `NotPortedError`.
- `src/DiscreteGlobalGrids.jl`'s `include` order *is* the dependency order:
  `Helpers` → `core/*` (interface, kernel, globe_ids, lookups, lookup_ops,
  grid_types, generic_cursor, face_grid, systems) → `ISEA` → per-system modules
  → `io/*` (I/O integrations, which read the systems and so come last).
- Every system has the same shape: a native layer, an `<X>Lookups`
  (DimensionalData) module nested in it, and an `<X>Kernel.jl` wiring file.
  `S2`/`ISEA4R` are chart systems: no `<X>Lookups` module, geometry only.

## Export / naming rules (easy to get wrong)

- Nothing system-level is re-exported at top level. Systems share generic
  vocabulary (`cell_center`, `cell_boundary`, `lonlat_to_cell`, `Touching`), so
  reach them through their submodule: `using DiscreteGlobalGrids.H3.H3Lookups`.
- Kernel generics `num_cells`, `cell_boundary`, `cell_center` are **not**
  exported (a submodule already claims each name); call them qualified —
  `DiscreteGlobalGrids.num_cells(system, level)`.
- Submodule capitalization is `HEALPix` (deliberately not `Healpix`, so it never
  shadows the registered Healpix.jl package).
- The `IGEO7` registered package in `Project.toml` deps is **not imported in
  `src`** — it is the black-box oracle only. The clean-room code lives in the
  `IGeo7` submodule.

## IGEO7 (ISEA7H + Z7) — cells, geometry, indexing

Files in `src/IGeo7/`, include order = dependency order:
`z7.jl` (pure integer) → `z7_ranges.jl` (monotonic line + lazy id vectors, also
pure integer) → `engine.jl` (Eisenstein lattice) → `grid.jl` (geometry + dense
indexing) → `IGeo7Lookups.jl` → `IGeo7Kernel.jl`.

- **The Z7 `UInt64` is the cell id.** `cell_to_z7` / `z7_to_cell` are
  validation-only identities. Bit layout: base cell (0:11) in bits `[63:60]`;
  digit `k` (1..20) in bits `[62-3k : 60-3k]`; digit `7` is the padding
  sentinel, not a child. Canonical full-world order = ascending `UInt64`.
- Two resolution caps: `MAX_RESOLUTION = 19` (has geometry), `Z7_MAX_RESOLUTION
  = 20` (valid for prefix/string/hex ops only — geometry constructors reject it).
- **String codec**: `z7_to_string` / `z7_from_string` → `"00"`, `"0800433"`.
  **Hex codec**: `z7_to_hex` / `z7_from_hex` → fixed 16-char, order-preserving.
- **Dense index**: `cell_to_index` / `index_to_cell` (1-based rank at own
  resolution). `num_cells(r) = 10·7^r + 2` (12 pentagons + `10·7^r − 10`
  hexagons). Pentagons have no deleted-digit child, so their subtree size is
  `p(d) = (5·7^d + 1)/6`.
- **Geometry** (`grid.jl`): `cell_center(id)` → `(lon, lat)`; `cell_boundary(id)`
  → `Vector{(lon,lat)}` (6/5 corners, `closed_ring=true` repeats first);
  `cell_boundary_cartesian(id)` → unit-sphere xyz tuples; `cell_area(id)` is a
  closed form (Snyder is exactly equal-area). `lonlat_to_cell`/`lonlat_to_z7`
  decode via Snyder forward → face's 3 corner bases → strict re-encode accept.
- **Errors**: `InvalidZ7Error(reason::Symbol, value, got, limit[, input])` — the
  human-readable message is built lazily in `Base.showerror`, never on the throw
  path. `index_to_cell` out-of-range throws **`BoundsError`** (deliberate API
  shape); the kernel layer translates to `DGG.OrdinalRangeError`.
- Lookups: `IGeo7Lookups.IGeo7Lookup` wraps sorted `UInt64` ids (sorted ascending
  = canonical order, so lookup position == `cell_to_index`). Globe-complete dims
  use `DGGSGlobeIds` (lazy, two words). `Touching` selects by polygon intersection,
  `Contains` by cell center.

## DGGS-Zarr reader (`src/io/DGGSZarr.jl`)

Reads [zarr-conventions/dggs
v1](https://github.com/zarr-conventions/dggs/blob/v1/README.md) archives as lazy
`YAXArrays.Dataset`s whose spatial dim is a real `IGeo7Lookup`. Entry point
`open_dggs_dataset`; demo `examples/z7_zarr_read.jl`. Five things bite here:

- **`z7_to_monotonic` is not `cell_to_index`.** The `compression: "ranges"`
  form stores `[start, end]` pairs on the **base-7 positional line**
  (`base·7^r + Σ dₖ·7^(r−k)`, `12·7^r` slots, deleted-digit paths included).
  `cell_to_index` is the pentagon-aware dense rank over `10·7^r + 2` real
  cells. Both are strictly increasing in the packed `UInt64`, so they agree on
  *order* and differ on *spacing* — expanding a range table with the wrong one
  silently yields the wrong `N`. `src/IGeo7/z7_ranges.jl`; pinned by
  `test_monotonic.jl` "distinct from cell_to_index".
- **Zarr.jl reverses dimension order.** On-disk `(R, 2)` arrives as `(2, R)`,
  so `_ARRAY_DIMENSIONS` must be reversed (`DGGSZarr.julia_dimensions`). 1-D
  variables are unaffected, which is why the mistake survives a casual test.
- **`YAXArrays.open_dataset` cannot open a ranges archive** — it demands a
  coordinate array per dimension and dies on `ranges`/`bounds` with
  `KeyError: key "bounds" not found`. Hence the `Dataset` is assembled from a
  raw `zopen`.
- **Placement and datum are archive-declared and non-default.** The
  Python-written archives set `dggs_vert0_lon = 11.2` (package default
  `ISEA_LON0 = 11.25`) and `igeo7_wgs84_geodetic_conversion = true`, while
  `IGeo7.cell_center` returns *authalic* latitude. Ignoring the first shifts
  longitude 0.05°; ignoring the second shifts latitude ~0.115° (~13 km).
  Neither fails — it just returns wrong coordinates. `DGGSZarr.cell_centers` /
  `cell_boundaries` / `sel_latlon` honour both; the same-named `IGeo7Lookups`
  functions deliberately do not.
- **Lazy id vectors** (`Z7RangeIds`, `Z7CachedIds`, both `<: Z7LazyIds`) are the
  `DGGSGlobeIds` pattern again: `IGeo7Lookup{<:Z7RangeIds}` is an ordinary
  lookup. Both need the `IGeo7Lookup(::Z7LazyIds)` constructor fast path or
  `DD.rebuild` materializes the dimension. `Z7CachedIds.getindex` serves the
  **first and last** elements from the source and materializes for any other
  index — deliberate, not an oversight: the lookup constructor and DD's compact
  display probe exactly those two, while walking a chunked compressed store one
  element at a time cost 33 s / 49 GiB for 1.3 MB of ids.

`examples/z7_zarr_plot.jl` maps a variable as real hexagons. Note that
`GeoAxis` builds degree tick labels itself with `round(x; sigdigits = 3)` and
**ignores `xtickformat`/`ytickformat`** (GeoMakie `src/geoaxis.jl`), so on a
sub-degree extent any tick step finer than 0.1° prints duplicate labels.

## Provenance constraint (do not violate)

IGEO7 is **clean-room**. `z7_ranges.jl` is unaffected: the base-7 place-value
reading of a digit string is first-principles arithmetic, and it matches the
`z7py` writer (the user's own code), not the reference. The AGPL reference
implementation is excluded from the
package and enters only as an independent black-box validation oracle: the suite
checks against dumps of its CLI output in `test/IGeo7/vectors/*.csv` (100% exact
z7 decode on all 196,080 oracle cell centers, res 1–5). Never port from or open
the reference source; every fitted constant is recorded against the vectors it
was fit to (see the `Provenance`/`[fitted]` citations in `z7.jl`/`grid.jl`).

## Stale docs note

The README references `docs/IGeo7/`, `docs/design/`, and `docs/reference/`, but
those directories are **not present in this checkout** — only `docs/src/` exists
(`index.md`, `all_dggs.md`, `tutorials/`). Don't spend time hunting for them.
