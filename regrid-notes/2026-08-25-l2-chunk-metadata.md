# L2 — one reading of a source's declared chunking

- Date: 2026-08-25
- Card: `regrid-notes/2026-08-21-regridding-simplification-plan.md`,
  "Phase 8 — private lazy state", Task L2
- Commit: `Normalize lazy DiskArrays chunk metadata`

## 1. The readings that existed

Four places asked what a source declares about its own chunking, each spelling
the test out again:

| where | what it asked | what it answered |
|---|---|---|
| `isdiskbacked` (`api.jl:4`) | `haschunks(data) isa Chunked` | the `lazy` default, and whether `flatsource` materializes before reshaping |
| `_passthroughchunks` (`lazy.jl:419`) | `haschunks`, then `eachchunk`, then `isa GridChunks`, then a dimension count | the output's chunks on the pass-through dimensions |
| `_sourceotherchunks` (`lazy.jl:431`) | the same four, again | the non-spatial chunk ranges a read splits and `knownempty` is asked about |
| `_spatialchunks` (`rastergrid.jl:410`) | `haschunks`, then `eachchunk`, then `isa GridChunks` | a derived `RasterGrid`'s `xchunks` and `ychunks` |

The middle two ran on every lazy application, one after the other, so a source
that declares chunked storage was asked `eachchunk` twice. A plan declaring a
chunking of its own brought that back to once: the output took the plan's, and
only the read groups still came off the source.

## 2. The reading that remains

`SourceChunking` (`lib/GlobalRegridding/src/lazy.jl`) is that one reading. It
holds what the source declares about its pass-through dimensions in the three
shapes the array needs:

```julia
struct SourceChunking{NO,P<:Tuple}
    passthrough::P                                 # as declared
    splits::NTuple{NO,Vector{UnitRange{Int}}}      # the same chunks as ranges
    groups::Vector{NTuple{NO,UnitRange{Int}}}      # their combinations
end
```

A `LazyRegridArray` builds one when it is applied to a source and holds it as
`A.chunking`. The output's chunk grid takes `passthrough` where the plan
declares no chunking of its own; `_slicegroups` splits a requested slice range
along `splits`; `_allempty` asks `knownempty` about `groups`. Nothing else asks
the source, and the array keeps no second description: the two fields that held
one, `otherchunks` and `othergroups`, are gone.

The declaration is usable when `eachchunk` answers a `GridChunks` over the
regrid's own dimensions — the spatial dimensions the source flattens over,
followed by the pass-through ones. Unchunked storage, an answer that is not a
chunk grid, and a grid describing a different number of dimensions are all "no
usable declaration", and each pass-through dimension is then one whole chunk.
None of the three is an error, which is what the source-side cases pin; a
malformed `chunks` on the *plan* still throws, with its existing messages.

Spatial chunks are not carried. A lazy array never needs them: a source's
spatial chunking is the source space's to describe, and a read addresses it
through `chunkranges`.

`eachchunk` is now consulted exactly once per lazy application over a source
that declares chunked storage, and not at all over one that does not.

## 3. `isdiskbacked` is `declareschunks`

The predicate never tested where data lives. It is `declareschunks(data)`
(`api.jl`), and its callers are the `lazy` default in `plan_regrid` and its
three documented signatures, `flatsource`'s materialize-before-reshape test,
`SourceChunking`, and `_spatialchunks` — which now calls it rather than
spelling the same test a fourth time.

No deprecation shim: the name was neither exported nor `public`, and the
repository had no caller of it outside `lib/GlobalRegridding/src`.

The residence test stays separate and now says so where it is defined:
`_isdisksource` (`lazy.jl`) is what decides whether a source is read through
`DiskArrays.readblock!` or copied from, and an in-memory array that declares a
chunking is exactly the case that separates the two.

## 4. Raster spatial chunks stay where they are

`_spatialchunks` still reads the declaration itself. Deriving it from the lazy
array's value would take one of two couplings: a `RasterGrid` is constructed
from an array long before — and often without — any lazy array, so either the
space constructor would have to be handed the lazy array's reading, or the
reading would have to be parked on the plan, which is reusable across
compatible arrays while a chunk description is not. What the two do share is
the predicate, and a raster's chunks are its own geometry, fixed at
construction: it is the same question asked of a different object at a
different time, not a second answer to the same one.

## 5. Types

The spatial count reaches `SourceChunking` as a `Val`, so each dimension's
chunk description keeps its own type. `Base.return_types` for a concrete source
gives one concrete `SourceChunking` in every case, and `@inferred` passes for a
regular grid, a grid mixing regular and irregular dimensions, and a plain array
that declares nothing. The mixed grid is what the old `_passthroughchunks` could
not type: reading `chunks[nspatial+i]` with a run-time `nspatial` gave a tuple
of `Union{RegularChunks,IrregularChunks}`, which the array's `chunks` field then
carried.

No field is abstract, no container holds `Any`, and no function here returns a
`Union`: an unusable declaration is a `SourceChunking` over whole dimensions,
never `nothing`.

## 6. What was checked against the previous behaviour

The same fixtures were run against the untouched checkout and against this one:
sources declaring regular chunks, irregular chunks, whole-dimension chunks,
nothing at all, a non-grid answer, and a grid of the wrong rank, each with and
without a plan chunking, plus a source with no pass-through dimensions. The
output chunk grid (ranges *and* chunk-description types), the read groups, the
slice splitting of a partial request, the destination tiling, the values, the
`chunks`-keyword error messages and the `lazy` default were identical in every
case. The one difference is the count: `eachchunk` calls per application fell
from 2 to 1 wherever the source declares chunked storage, and stayed at 1 where
the plan declared a chunking and at 0 where the source declares none.

No test was dropped: nothing tested `_passthroughchunks` or `_sourceotherchunks`
apart from the results they fed, which the new cases pin directly.

## 7. Suites

- `lib/GlobalRegridding`, `--depwarn=yes`: 4757 passed, 0 failed, 0 errored, 1
  broken, no deprecation warning. The baseline is 4713 / 0 / 0 / 1; the 44 new
  assertions are the source-chunking cases.
- Bridge cross-check (the three DGG↔regridding suites): `passed=383 failed=0
  errored=0 broken=12`, unchanged.
- Weight and value comparator against the conservative, nearest and bilinear
  fixtures: 54 entries, all `==` and all `isequal`.
