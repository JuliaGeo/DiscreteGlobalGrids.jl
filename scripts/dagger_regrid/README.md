# Dagger regrid

An experimental distributed CopDEM runner. Our graph and scheduler choose
complete destination-chunk batches; Dagger places those batches on worker
processes. Workers keep their own source cache and write disjoint Zarr chunks.
`main.jl` stays top-level so its setup and coordinator loop can be evaluated
section by section in a Julia REPL.

Run it from the repository root:

```sh
julia --project=scripts/dagger_regrid -p 2 -t 1 scripts/dagger_regrid/main.jl
```

`DAGGER_REGRID_MODE` selects `run` (default), `canary`, or `smoke`. Set
`COPDEM_MAXCHUNKS=1` for a one-chunk canary. The shared CopDEM configuration is
in `copdem_helpers.jl`; Dagger-specific admission and cache settings are near
the top of `main.jl`. `RASTERDATASOURCES_PATH` selects the data directory; the
public Copernicus tile list is downloaded and cached there when missing.

The isolated project tracks Dagger's `master` branch. Instantiate it once with
`julia --project=scripts/dagger_regrid -e 'using Pkg; Pkg.instantiate()'`.
