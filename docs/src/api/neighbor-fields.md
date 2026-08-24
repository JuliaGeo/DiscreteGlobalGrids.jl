# Requesting neighbour fields

```@meta
CurrentModule = DiscreteGlobalGrids
```

A neighbourhood kernel almost never wants only cells. It wants an elevation, an
index into its own storage, a centroid — and with [`Neighbors`](@ref) it is
handed cell handles and has to go and get each of those itself, once per
neighbour, from inside the innermost loop: one array read per neighbour, issued
by the callback for every cell in the sweep, or a `cell_centroid` call whose trig
is recomputed for every edge the cell shares. A **request** turns that around.
`needs = (Value(dem), Centroid())` states up front which per-neighbour
quantities the kernel reads, and the sweep streams exactly those out of the one
membership clip it already performs per cell. The callback reaches back into
nothing.

Each entry of `needs` is an [`AbstractNeed`](@ref). There are four of them, and
one — [`Index`](@ref DiscreteGlobalGrids.Index) — also names the index space it
answers in:

| need | the visited cell | a ring slot |
|---|---|---|
| [`Cell`](@ref)`()` | its cell id | the neighbour's cell id |
| [`Index`](@ref DiscreteGlobalGrids.Index)`(`[`Local`](@ref)`())` | its index in the collection | the neighbour's |
| [`Index`](@ref DiscreteGlobalGrids.Index)`(`[`Global`](@ref)`())` | its [`globalindex`](@ref) | the neighbour's |
| [`Index`](@ref DiscreteGlobalGrids.Index)`(T)` | its id `reindex`ed to scheme `T` | the neighbour's |
| [`Value`](@ref)`(a)` | `a` at its index | `a` at the neighbour's index |
| [`Centroid`](@ref)`()` | its centroid, on the unit sphere | the neighbour's |

`Value` takes any vector laid out against the collection, and `needs` may carry
as many of them as the kernel reads — two arrays of different element types give
two rings of those element types. A single one is what [`Values`](@ref) already
passes; a request is that idea without the one-array limit and with the geometry
in it. Everything else about the call is unchanged: `needs` is a keyword on
[`mapneighbors`](@ref) and [`foreachneighbors`](@ref), alongside `order`,
`threaded` and `connectivity`, over a [`CellVector`](@ref), a
[`PartialGrid`](@ref), a [`CellLookup`](@ref) or a dimensional array.

## The callback

`f(center, rings)`. Both arguments are tuples with one entry per need, in the
order `needs` states them; there is no special case for a single need.
`center` holds the visited cell's answers and `rings` holds the neighbours'.
Three things about them are the contract:

  - **Rings are field-major.** `rings[j]` is need `j`'s value for *every*
    clipped neighbour, so a two-need request gives two rings, not one ring of
    pairs. Slot `i` of every ring names the same neighbour, in the order
    [`neighbors`](@ref) states — counter-clockwise seen from outside the sphere,
    clipped to membership. A kernel that wants neighbour-major records writes
    `zip(rings...)`, on its own side; the sweep does not build them, because
    most kernels reduce one field at a time and would pay to take the tuples
    apart again. Each ring is a `SmallVector` where the system declares a
    [`maxneighbors`](@ref) and a `Vector` where it does not, which is the same
    rule the [`Values`](@ref) pass follows.

  - **`Index(Local())` is the index in the collection the caller passed** — the
    [`CellVector`](@ref), or the cube's cell axis — never an index into some
    block the sweep chose internally. On a complete grid that equals the global
    index; on a subset it does not, and the subset's numbering is the one you
    get. That pin holds today, and it binds whatever comes next: there is no
    chunked form of `needs` yet, and when there is, it must translate each
    chunk's numbering back to the caller's so its result equals the whole-axis
    sweep's, cell for cell. In this milestone a request does not take the chunk
    route at all; [`mapneighbors!`](@ref) over a
    [plan](@ref "Sweeping a cube along its chunk lines") is what does, in
    [`Values`](@ref)' form.

  - **`Centroid()` is `cell_centroid` on the unit sphere, answered from a
    working set.** A cell's centroid is touched once as a centre and once from
    every ring that names it, so each sweep task keeps a bounded cache of
    recent centroids and computes each one on its first touch inside it —
    automatic, with no knob, and the same values `cell_centroid` gives either
    way. The key is the local index rather than the visit position, which is
    what keeps a locality-preserving `order` hitting; a random permutation
    scatters the keys, the window essentially stops hitting, and the centroid
    surcharge measured 4.2–4.4× storage order's on the same fixture
    (`benchmark/needs_centroid.jl`). In storage order it settles at roughly
    one computation per cell instead of one per touch, and that one is most
    of what is left: at 117,649 cells the sweep is 48.6 ms against 29.5 ms
    with the whole grid's centroids prebuilt as `Value(table)`, over a
    21.8 ms value-only floor. The table is materially faster — some 70% of
    the residual surcharge — at 24 bytes per cell and a build pass, so the
    bounded window is the default and the table the opt-in.
    The point is a direction, not a place: distance and bearing want a radius,
    and stay in the kernel.

## Steepest descent, in one pass

A drainage kernel is the case the request form was built for: it reads one data
array and the geometry, and it reads both for every neighbour. On the unit
sphere the fall between two cells is an elevation difference over an angle, and
the direction is a unit tangent vector — no radius enters, so nothing here
commits to a datum.

```@example needs
using DiscreteGlobalGrids
using DiscreteGlobalGrids: Value, Centroid     # public, not exported
using LinearAlgebra: dot, normalize

sys   = IGeo7System()
cells = CellVector(subtree(sys, cellindex(levelgrid(sys, 1), 3), 5))

# A cone rising away from the first cell, so every cell drains towards it.
sink = cell_centroid(sys, cells[1])
dem  = [1000 * acos(clamp(dot(sink, cell_centroid(sys, c)), -1.0, 1.0))
        for c in cells]
length(cells), extrema(dem)
```

```@example needs
function steepest(center, rings)
    z, p   = center                   # this cell's elevation and centroid
    zs, ps = rings                    # one ring per need, field-major
    best, dir = 0.0, zero(p)
    for (zn, pn) in zip(zs, ps)       # neighbour-major records, on our side
        c = clamp(dot(p, pn), -1.0, 1.0)
        drop = (z - zn) / acos(c)     # metres per radian: no radius involved
        drop > best && ((best, dir) = (drop, normalize(pn - c * p)))
    end
    return (best, dir)
end

fall, direction = mapneighbors(steepest, cells; needs = (Value(dem), Centroid()))
nothing # hide
```

The kernel returned a concrete tuple, so the sweep split it the way it splits
any other: one vector per component, both in collection index order.

```@example needs
typeof(fall), typeof(direction)
```

Cell 1 sits at the bottom of the cone and has nowhere to go: no neighbour is
lower, the loop never fires, and it keeps the zero fall and zero direction it
started with. Every other cell drains, and its direction is a unit vector in the
tangent plane at its own centroid.

```@example needs
(fall[1], direction[1]), (fall[end], direction[end])
```

## Any per-cell quantity

`Centroid()` is not a special case inside the sweep. It resolves to a **cell
field** — a vector over the collection whose entries are computed by a per-cell
function instead of stored — and the sweep then reads that field like any other
`Value`. [`cellfield`](@ref) builds one, so these two requests are the same
request:

```@example needs
using DiscreteGlobalGrids: cellfield          # public, not exported

byname  = mapneighbors(steepest, cells; needs = (Value(dem), Centroid()))
spelled = mapneighbors(steepest, cells;
    needs = (Value(dem), Value(cellfield(cell_centroid, cells))))
byname == spelled
```

Any function of `(grid, cell)` can be a field — `cell_area` as readily
as `cell_centroid`. A field is pure: it never mutates and never
remembers, so one field is safely read by every task of a threaded sweep, and
what remembers is the bounded window the sweep gives each task.

`known` hands the field what you have already computed. A vector on the
collection's cell axis is the complete case, and a complete field is read
straight through with no window at all:

```@example needs
table = [cell_centroid(sys, c) for c in cells]
whole = cellfield(cell_centroid, cells; known = table)
mapneighbors(steepest, cells; needs = (Value(dem), Value(whole))) == byname
```

A one-dimensional cube on a `Cells` dimension is the partial case: the
cells it carries are read from it, and every other cell is computed. Nothing
requires the whole table, so precompute only the part that pays — the
[`border`](@ref), say, whose neighbours lie outside the collection's own index
range and so miss the window most often:

```@example needs
using DimensionalData: DimArray

edge = cells[collect(border(cells))]
part = cellfield(cell_centroid, cells;
    known = DimArray([cell_centroid(sys, c) for c in edge],
                     (Cells(CellLookup(edge)),)))
length(edge), mapneighbors(steepest, cells;
    needs = (Value(dem), Value(part))) == byname
```

A field is read by local index, so it must be over the collection being swept.
One built over other cells is an `ArgumentError` raised before the sweep
starts, not a silently wrong answer.

## The needs

```@docs
AbstractNeed
Cell
Index
Local
Global
Value
Centroid
cellfield
```

## Index

```@index
Pages = ["api/neighbor-fields.md"]
```
