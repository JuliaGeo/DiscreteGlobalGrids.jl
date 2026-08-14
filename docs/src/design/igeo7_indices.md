# Fixed-resolution IGEO7 indexing

## Status

Design specification for a first implementation. The initial scope is a
one-dimensional `DimArray`/`Raster` whose only dimension has an `IGeo7Lookup`,
and geometric translation within one non-pentagon IGEO7 chart at one
resolution.

Crossing an icosahedron edge or a pentagon distortion is deliberately out of
scope for the first implementation. See [Deferred topology work](#deferred-topology-work).

## Problem

The hydrology tutorial produces and reloads a level-13 IGEO7 raster:

```julia
using JLD2
rl = load("dem_igeo7_ras.bin", "dem_igeo7_ras")
```

The example raster has 5,764,801 cells. Its lookup contains the complete
level-13 descendants of one non-pentagon level-5 cell, in ascending Z7 order.
Callers should be able to:

1. iterate with canonical IGEO7 cell indices;
2. index the raster with such an index;
3. enumerate stored or global edge neighbors;
4. express one-cell and larger geometric displacements in readable hex
   coordinates; and
5. compose these operations into future stencil algorithms.

Three different concepts must not be conflated:

- a **cell identity** is the canonical Z7 `UInt64`;
- an **array position** is a 1-based `Int` into one particular lookup; and
- a **geometric displacement** is a signed vector on the resolution's
  Eisenstein lattice.

The Z7 ID contains generalized balanced ternary digits, but its `UInt64`
layout is a packed hierarchical code, not an integer coordinate. Ordinary
subtraction of two Z7 values includes base-cell bits, padding sentinels,
per-level carries, alternating aperture-7 chirality, and pentagon deletion.
It therefore does not produce a geometric displacement. The existing raw `horner` and `decode_step` machinery is the correct seam for
that conversion. `_encode_lattice` is not: it applies per-cell cone wrapping,
so subtracting two of its outputs is not translation-invariant when the cells
fall on opposite sides of that cut.

## Design principles

This is one in-process module. No Adapter is warranted.

The Interface is intentionally smaller than the Implementation:

- callers see identity, position, displacement, iteration, indexing, and
  neighbors;
- the Implementation hides Z7 digits, ordinal offsets, GBT carry handling,
  chirality, Eisenstein conversion, lookup search, and neighbor steppers.

This Depth gives callers Leverage while keeping topology rules local to the
IGEO7 implementation.

## Interface

### Canonical cell identity

```julia
struct IGEO7Index
    id::UInt64

    function IGEO7Index(id::Integer)
        # Convert without truncation, validate as an IGEO7 geometry cell,
        # and reject resolution 20.
    end
end

IGEO7Index(text::AbstractString) = IGEO7Index(z7_from_string(text))
```

`IGEO7Index` has value semantics. Equality, ordering, and hashing use `id`.
Ascending `IGEO7Index` order is ascending canonical Z7 order.

Required methods:

```julia
Base.convert(::Type{UInt64}, index::IGEO7Index)
Base.show(io::IO, index::IGEO7Index)
Base.isless(a::IGEO7Index, b::IGEO7Index)
get_resolution(index::IGEO7Index) -> Int
```

`show` should use the existing readable Z7 string, for example
`IGEO7Index("0c4d9...")`, rather than print only a decimal integer.

`IGEO7Index` is not an `Integer`: arithmetic on a cell identity is meaningful
only with a geometric displacement.

### Geometric displacement

```julia
struct RelativeIGEO7Index
    a::Int64
    b::Int64
    resolution::UInt8
end
```

`a + b*omega` is an Eisenstein integer in the raw resolution lattice produced
by the GBT Horner evaluation, before pentagon collapse and cone wrapping.
Within a non-pentagon subtree, the omitted transformation is constant and
this raw frame is affine. Components are signed because every useful
displacement has an inverse. `resolution` is part of the value because the
same lattice-unit count has a different physical scale at another resolution,
and applying a difference across resolutions is an error.

Do not encode this as the arithmetic difference between packed Z7 `UInt64`
values. A future implementation may pack the two signed components and the
resolution after proving bounds, but that is an internal representation
change and must not alter the Interface.

Required vector arithmetic:

```julia
Base.:-(d::RelativeIGEO7Index)
Base.:+(a::RelativeIGEO7Index, b::RelativeIGEO7Index)
Base.:-(a::RelativeIGEO7Index, b::RelativeIGEO7Index)
Base.:*(n::Integer, d::RelativeIGEO7Index)
directioncode(d::RelativeIGEO7Index) -> UInt8
```

Binary operations require equal resolutions and throw `DimensionMismatch`
otherwise. Provide `RelativeIGEO7Index(resolution)` as the zero displacement;
`zero(displacement)` preserves its resolution. `zero(RelativeIGEO7Index)`
remains undefined because a resolution is required.

For unit displacements, `directioncode` packs the shifted axial components
into a four-bit lookup key and returns `1:6`; zero maps to `0`. Larger
displacements are rejected because they do not name one edge direction.

### Human-readable hex coordinates

```julia
struct HexIndex
    i::Int64
    j::Int64
    k::Int64

    function HexIndex(i::Integer, j::Integer, k::Integer)
        i + j + k == 0 ||
            throw(ArgumentError("hex coordinates must satisfy i + j + k == 0"))
        # Checked Int64 conversion.
    end
end
```

The fields must be signed, not `UInt8`. The six unit directions contain
negative components, and larger jumps can exceed 255 cells.

The conversion follows the cube coordinates already documented by
`IGeo7.hex_round`:

```julia
Base.convert(::Type{HexIndex}, d::RelativeIGEO7Index) =
    HexIndex(d.a, d.b - d.a, -d.b)

RelativeIGEO7Index(h::HexIndex, resolution::Integer)
```

The six edge moves are:

```julia
HexIndex( 1, -1,  0)
HexIndex( 1,  0, -1)
HexIndex( 0,  1, -1)
HexIndex(-1,  1,  0)
HexIndex(-1,  0,  1)
HexIndex( 0, -1,  1)
```

These are local chart directions, not global east/northeast/etc. Compass
names would become misleading after face rotation.

### Cell arithmetic

```julia
Base.:-(a::IGEO7Index, b::IGEO7Index) -> RelativeIGEO7Index
Base.:+(index::IGEO7Index, d::RelativeIGEO7Index) -> IGEO7Index
Base.:-(index::IGEO7Index, d::RelativeIGEO7Index) -> IGEO7Index

trytranslate(index::IGEO7Index, d::RelativeIGEO7Index) ->
    Union{IGEO7Index,Nothing}
```

For the first implementation:

- both cells in `a - b` must have the same resolution;
- their lowest common Z7 ancestor must be non-pentagonal, establishing one
  affine lattice frame;
- the displacement resolution must match the cell resolution;
- translation must decode to a valid cell in that same frame; and
- translation must round-trip exactly:
  `result - index == d`.

`a - b` throws `ArgumentError` when no supported common frame exists.
`trytranslate` returns `nothing` when translation leaves the supported chart.
`index + d` and `index - d` call `trytranslate` and throw `DomainError` when
it returns `nothing`.

The successful arithmetic laws are:

```julia
a - a == RelativeIGEO7Index(get_resolution(a))
b + (a - b) == a
a - (a - b) == b
(a + d) - a == d
a + RelativeIGEO7Index(get_resolution(a)) == a
```

`a - a` returns the zero displacement even when `a` is a pentagon. Every
other law above is scoped to cells and translations satisfying the supported
non-pentagon-frame rules.

Large jumps use the same operations as one-cell moves. They must not be
implemented as repeated neighbor stepping.

### Array iteration and indexing

The direct integration is restricted by dispatch to a one-dimensional
`DimensionalData.AbstractDimArray` whose sole dimension contains an
`IGeo7Lookup`. Ordinary one- and two-dimensional rasters must retain Base's
existing indexing.

Conceptually, the narrow array type is:

```julia
const IGEO7DimVector =
    DD.AbstractDimArray{T,1,<:Tuple{<:DD.Dimension{<:IGeo7Lookup}}} where T
```

The implementer must use a valid non-ambiguous signature for the installed
DimensionalData version rather than exporting this alias.

Required methods:

```julia
Base.eachindex(A::IGEO7DimVector) -> IGEO7Indices
Base.getindex(A::IGEO7DimVector, index::IGEO7Index)
Base.setindex!(A::IGEO7DimVector, value, index::IGEO7Index)
Base.checkbounds(::Type{Bool}, A::IGEO7DimVector, index::IGEO7Index)
```

`IGEO7Indices` is a lazy `AbstractVector{IGEO7Index}` over the lookup's ID
vector:

```julia
struct IGEO7Indices{V<:AbstractVector{UInt64}} <: AbstractVector{IGEO7Index}
    ids::V
    resolution::Int
    first_ordinal::Int
    contiguous::Bool
end
```

Iteration order is lookup order, which `IGeo7Lookup` already requires to be
strictly ascending canonical-ID order. `IGEO7Indices` adds no per-cell
storage. Membership is O(1) for a lazy subtree or globe, O(resolution) for a
contiguous explicit lookup such as the saved hydrology raster, and O(log n)
for a genuinely sparse explicit lookup.

`A[index]` resolves the ID through the lookup:

- a `DGGSGlobeIds` lookup uses `cell_to_ordinal`;
- a `DGGSSubtreeIds` lookup uses `subtree_position`;
- an explicit sorted vector uses `cell_position`/binary search; and
- an explicit vector whose endpoint ordinals and length prove that it is a
  contiguous ordinal interval may use ordinal subtraction as a fast path.

An absent cell throws `BoundsError(A, index)`. `checkbounds(Bool, ...)`
returns `false`. Integer indexing (`A[1]`) remains ordinary positional
indexing.

Do not add a broad `eachindex(::AbstractDimArray)` method followed by a runtime
DGGS check. The dispatch itself must distinguish a real IGEO7 dimension so
ordinary rasters are unaffected.

Because changing `eachindex(A)` affects generic array code, acceptance is
conditional on compatibility tests for display, copying, broadcasting,
mapping, and reductions in both DimensionalData and Rasters. If narrow
dispatch still breaks those operations, implement the wrapper described under
"A mandatory raster wrapper" rather than ship a surprising `AbstractArray`.

For arrays with additional dimensions, do not change `eachindex` in the first
implementation. A later design can combine `IGEO7Index` with the other axes,
or provide a lightweight wrapper if DimensionalData's type parameters do not
permit safe narrow dispatch.

### Neighbors

Provide two scopes through dispatch:

```julia
neighbors(index::IGEO7Index) -> SmallVector{6,IGEO7Index}
neighbors(::IGEO7DGGS, index::IGEO7Index) -> SmallVector{6,IGEO7Index}
neighbors(A::IGEO7DimVector, index::IGEO7Index) -> SmallVector{6,IGEO7Index}
cellbearing(A::IGEO7DimVector, from::IGEO7Index, to::IGEO7Index) -> Float64
celldistance(A::IGEO7DimVector, from::IGEO7Index, to::IGEO7Index) -> Float64
cellarea(A::IGEO7DimVector, index::IGEO7Index) -> Float64
edges(A::IGEO7DimVector) -> Vector{IGEO7Index}
```

The system form returns every global edge neighbor: six for a hexagon and
five for a pentagon. It delegates to `cell_neighbors`.

The array form returns only neighbors stored in `A`, preserving ascending
canonical-ID order. It delegates to the existing neighbor machinery and
filters with lookup membership. Off-coverage neighbors are absent, matching
the current `stencil` and `subtree_stencil` semantics.

Neighbor order is the counterclockwise order of the six Eisenstein units,
constructed directly before the system kernel's canonical-ID sort.

`celldistance` requires both cells to be stored and returns the great-circle
distance between their centers in metres on IGEO7's WGS84 authalic sphere.
`cellarea` requires the cell to be stored and returns its closed-form area in
square metres.

`edges` returns every stored cell with at least one global neighbor outside
the stored coverage. Complete subtrees use `border_descendants` in O(result);
arbitrary sparse lookups use neighbor membership as the correct fallback.

## Example

```julia
using JLD2
using Rasters
using DiscreteGlobalGrids

rl = load("dem_igeo7_ras.bin", "dem_igeo7_ras")

for index in eachindex(rl)
    center = rl[index]
    for neighbor in neighbors(rl, index)
        drop = center - rl[neighbor]
        # Hydrology operation...
    end
end

index = first(eachindex(rl))
resolution = get_resolution(index)

step = RelativeIGEO7Index(HexIndex(1, -1, 0), resolution)
next = trytranslate(index, step)
next === nothing || @assert rl[next] isa eltype(rl)

farther = trytranslate(index, 100 * step)
```

For a pair of cells:

```julia
a, b = eachindex(rl)[100], eachindex(rl)[10_000]
d = a - b

@assert convert(HexIndex, d).i +
        convert(HexIndex, d).j +
        convert(HexIndex, d).k == 0
@assert b + d == a
```

There is intentionally no `convert(LinearIndex, d)`. Julia has
`LinearIndices`, not a scalar `LinearIndex`, and a geometric displacement has
no anchor-independent space-filling-curve offset. The explicit conversions
are:

```julia
cell_to_position(A, index::IGEO7Index) -> Int
position_to_cell(A, position::Integer) -> IGEO7Index
```

The difference between two returned positions is a storage-order difference,
not a geometric displacement.

## Implementation seam

Place the public types and generic-facing methods in a focused IGEO7 indexing
file included by `IGeo7.jl`. Keep DimensionalData-specific methods beside
`IGeo7Lookups`.

Export `IGEO7Index`, `RelativeIGEO7Index`, `HexIndex`, `trytranslate`,
`neighbors`, `celldistance`, `get_resolution`, `cell_to_position`, and
`position_to_cell` from `IGeo7`. Import and re-export those same bindings from
`DiscreteGlobalGrids` so the direct raster example works with
`using DiscreteGlobalGrids`. These names do not collide with the
system-specific names that the package intentionally leaves unexported.

The Implementation should reuse:

- `horner` or an allocation-free equivalent over a Z7 ID's active digits for
  canonical ID to raw axial coordinates;
- `decode_step`/`digit_decode!` or the proven `_suffix_lattice` and
  `_lattice_suffix` logic for exact inverse conversion;
- the existing exact Eisenstein arithmetic and `hex_round` cube-coordinate
  convention;
- `cell_neighbors` for global adjacency;
- `DGGSGlobeIds`, `DGGSSubtreeIds`, `cell_position`, and
  `subtree_position` for array lookup; and
- `neighbor_stepper`/`step_neighbors` for a complete subtree fast path.

Do not create a second stencil engine. Future stencil syntax should compile
to the existing `stencil` or `subtree_stencil` implementation.

Private helpers should own:

1. finding and validating a shared non-pentagon frame;
2. checked conversion between Z7 and the raw Horner lattice;
3. exact translated-coordinate digit decoding, same-base Z7 reconstruction,
   cell validation, and round-trip validation; and
4. detection of the narrow IGEO7 `DimArray` type.

Do not route cell arithmetic through `_encode_lattice`, `dev_to_xyz`, or
`_xyz_to_z7`. Those functions correctly canonicalize spherical geometry, but
their cell-dependent cone-wrap branch is precisely what the supported local
affine frame excludes.

## Errors

| Condition | Result |
|:--|:--|
| Invalid Z7 ID or resolution 20 | existing `InvalidZ7Error` |
| Cell/difference resolution mismatch | `DimensionMismatch` |
| No shared supported non-pentagon frame for subtraction | `ArgumentError` |
| Translation leaves the supported frame | `nothing` from `trytranslate` |
| `+`/`-` translation leaves the supported frame | `DomainError` |
| Cell is not stored in an indexed array | `BoundsError` |
| Invalid cube coordinates | `ArgumentError` |
| Checked integer conversion overflows | `InexactError` |

Do not silently reinterpret an unsupported geometric operation as ordinal
arithmetic.

## Performance requirements

- All three scalar index types are immutable and isbits.
- Constructing and iterating `IGEO7Indices` allocates no vector of wrappers.
- `IGEO7Index +/- RelativeIGEO7Index` performs work proportional to the
  resolution, not to displacement magnitude.
- `IGEO7Index - IGEO7Index` performs work proportional to the resolution.
- Global neighbor enumeration remains bounded and uses a fixed-capacity
  container.
- A complete non-pentagon subtree uses the existing integer neighbor stepper
  for stored-neighbor iteration.
- Direct indexing of globe and subtree-backed lookups uses arithmetic, not
  binary search.
- No operation materializes the 5,764,801 lookup indices merely to iterate.

## Acceptance tests

### Type and validation tests

- `IGEO7Index` accepts valid resolution 0:19 IDs and rejects invalid/resolution
  20 IDs with `InvalidZ7Error`.
- `HexIndex` accepts exactly the signed triples summing to zero.
- all scalar index types are isbits.
- ordinary 1-D and 2-D arrays and non-DGGS rasters retain their existing
  `eachindex`.

### Iteration and array tests

- loading `dem_igeo7_ras.bin` with both `Rasters` and `JLD2` in scope
  reconstructs a real raster;
- `length(eachindex(rl)) == length(rl)`;
- every yielded index is an `IGEO7Index`;
- iteration is ascending by ID;
- `rl[index] == parent(rl)[position]` for sampled and boundary positions;
- absent and wrong-resolution indices fail as specified;
- a `DGGSSubtreeIds`-backed lookup follows the arithmetic fast path; and
- display, `copy`, broadcast, `map`, and reductions retain normal
  DimensionalData/Rasters behavior after the `eachindex` specialization.

### Arithmetic tests

- the algebraic laws above hold exhaustively for manageable non-pentagon
  subtrees and on a seeded level-13 sample from the tutorial raster;
- all six unit `HexIndex` values agree with `cell_neighbors` for interior
  cells;
- larger random displacements round-trip exactly when their target remains in
  the supported frame;
- resolution mismatches fail;
- pentagon and face-edge cases fail explicitly rather than returning a wrong
  cell; and
- arithmetic allocates zero bytes after compilation.

### Neighbor tests

- `neighbors(IGEO7DGGS(), index)` agrees with `cell_neighbors`;
- `neighbors(rl, index)` is the stored subset of global neighbors;
- interior tutorial-raster cells have six stored neighbors;
- border cells omit off-raster neighbors; and
- pentagons return five distinct global neighbors.

## Rejected alternatives

### Packed Z7 subtraction

Rejected because numeric `UInt64` subtraction is subtraction of a serialized
prefix code, not of GBT/Eisenstein coordinates. Bit shifts can extract digits,
but geometric subtraction still requires mixed-radix carry handling and the
alternating chirality already implemented by the lattice engine.

### Relative index as array position

Rejected because it would make `+1` mean "next along canonical Z7 order", not
"adjacent cell". Array positions remain plain `Int`.

### Global cube coordinates

Rejected for the first implementation because an icosahedral hex/pentagon
grid has no single translation frame over the sphere. Exposing one would hide
face rotations and the topological defects at the twelve pentagons.

### A mandatory raster wrapper

Fallback if narrow DimensionalData dispatch proves ambiguous, brittle, or
breaks generic `AbstractArray` behavior. The loaded raster's dimension type
carries `IGeo7Lookup`, so direct integration should be tried and tested first.
A wrapper is preferable to a broad runtime-checked `eachindex` override or an
array that violates Base expectations.

## Deferred topology work

- [ ] Define frame transport across an icosahedron edge, including the
  rotation/reflection applied to a `RelativeIGEO7Index`.
- [ ] Specify whether multi-edge translation is path-independent. If it is
  not, introduce an explicit cursor carrying transported frame state rather
  than changing the meaning of `RelativeIGEO7Index`.
- [ ] Define movement at and around pentagons, including the missing direction
  and holonomy around a pentagon.
- [ ] Extend translation tests to every face edge, vertex neighborhood, and
  both hemispheres.
- [ ] Generalize identity, position, and displacement vocabulary to other DGGS
  systems only after IGEO7 semantics are stable.
- [ ] Design multi-dimensional `DimArray` iteration without changing ordinary
  Cartesian indexing.
