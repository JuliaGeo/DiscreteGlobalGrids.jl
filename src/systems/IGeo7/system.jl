# ---------------------------------------------------------------------------
# The new-interface wiring: `Z7Cell`, `IGeo7System`, `IGeo7Grid`.
#
# Everything below is an adapter over the native layers included ahead of it
# (`z7.jl`, `engine.jl`, `z7grid.jl`), which are ported verbatim from the
# verified clean-room implementation and are the oracle-validated arithmetic.
# No projection maths is rederived here; what changes is only the shape of the
# calls — the old `(system, level, id)` triple becomes a typed cell that knows
# its own level, and dense ordinals become grid *positions*.
#
# Namespace note: `z7grid.jl` defines native `cell_boundary`, `cell_area` and
# `cell_center` that answer in `(lon, lat)` degrees. Those names are NOT
# imported from the package, so inside this module they stay the native ones and
# the interface generics are reached through `DGG.`; the two are never
# conflated. This is the same discipline the old `IGeo7Kernel.jl` used.
# ---------------------------------------------------------------------------

"""
    Z7Cell(id::UInt64) <: AbstractCellIndex
    Z7Cell(s::AbstractString)

The canonical cell id of [`IGeo7System`](@ref): the Z7 `UInt64` itself, whose
level is carried in-band by the digit slots (4 bits of base cell `0:11`, then
twenty 3-bit digits, `7` marking a padded slot past the resolution).

`level(c)` is the count of leading active digit slots, so it needs no system and
no table. `rawid(c)` is the `UInt64`, which is what `z7_to_hex` /
`z7_to_string` encode and what a file or a C library should carry.

Ordering is ascending `UInt64`, which **is** the system's canonical cell order:
the base cell occupies the high bits and the padding sentinel `7` sorts after
every active digit, so comparing two same-level ids compares
`(base, d_1 … d_r)` lexicographically — the space-filling curve order that makes
[`has_sorted_subtrees`](@ref) true and every subtree one contiguous interval of
grid positions.

Construction does **not** validate, matching `LevelIndex`: an id is a cheap
name, and validation happens where a name meets a system or a grid
(`cellposition` answers `nothing`, the geometry entry points throw
[`InvalidZ7Error`](@ref)). The string constructor is the exception — it goes
through `z7_from_string`, which validates as it parses.

```julia
julia> c = Z7Cell("0941054");

julia> level(c)
5
```
"""
struct Z7Cell <: AbstractCellIndex
    id::UInt64
end

Z7Cell(c::Z7Cell) = c
Z7Cell(id::Unsigned) = Z7Cell(UInt64(id))
Z7Cell(s::AbstractString) = Z7Cell(z7_from_string(s))

# `_z7_leading_resolution` rather than `z7_resolution`: `level` is total on
# `AbstractCellIndex` by contract, so it reports the leading-digit count of
# whatever it is handed rather than validating the padding tail. Ids that are
# not cells are rejected where they meet a grid or a geometry call.
level(c::Z7Cell) = _z7_leading_resolution(c.id)
rawid(c::Z7Cell) = c.id

Base.isless(a::Z7Cell, b::Z7Cell) = isless(a.id, b.id)

Base.show(io::IO, c::Z7Cell) =
    print(io, "Z7Cell(\"", is_valid_cell(c.id) ? z7_to_string(c.id) : z7_to_hex(c.id), "\")")

"""
    z7_string(c::Z7Cell) -> String
    z7_hex(c::Z7Cell; prefix=false) -> String

The two text codecs of a cell id, on the typed wrapper: the Z7 string form
(two base characters then one character per active digit) and the fixed
16-character hexadecimal form. Both round-trip through [`Z7Cell`](@ref).
"""
z7_string(c::Z7Cell) = z7_to_string(c.id)

# The docstring above documents both codecs, but Julia attaches it to the
# binding on the line that follows it — so `z7_hex` needs the attachment made
# explicitly, or the exported name ships undocumented.
@doc (@doc z7_string)
z7_hex(c::Z7Cell; prefix::Bool=false) = z7_to_hex(c.id; prefix)

"""
    is_pentagon(c::Z7Cell) -> Bool

Whether `c` is one of the twelve pentagons of its level — every active digit
zero, one per icosahedron vertex. Pentagons are the cells with five neighbours
and six children instead of six and seven.
"""
is_pentagon(c::Z7Cell) = is_pentagon(c.id)

"""
    is_valid_cell(c::Z7Cell) -> Bool

Whether `c` names a cell that exists: a well-formed Z7 id at a level in
`0:$(MAX_RESOLUTION)` whose digit chain does not take one of the twelve
pentagons' deleted branches.

This is the total, non-throwing test — the same one [`cellposition`](@ref)
answers `nothing` from. The geometry entry points instead throw
[`InvalidZ7Error`](@ref), because a caller asking for the boundary of a cell
that does not exist has a bug rather than a miss.
"""
is_valid_cell(c::Z7Cell) = is_valid_cell(c.id)

# ---------------------------------------------------------------------------
# The system
# ---------------------------------------------------------------------------

"""
    IGeo7System() <: AbstractHierarchicalGridSystem

The IGEO7 discrete global grid system: an aperture-7 hexagonal hierarchy on the
icosahedron (ISEA7H) with Z7 indexing, in the standard ISEA placement.

Twelve pentagons — one per icosahedron vertex — and `10·7^r − 10` hexagons tile
level `r`, for `10·7^r + 2` cells in all. Refinement is aperture 7: a hexagon
has seven children, a pentagon six, and children **overhang** their parent's
boundary, which is why [`node_extent`](@ref) is an inflated cap rather than the
cell polygon (see the covering law).

| trait | value |
|:--|:--|
| [`cellindextype`](@ref) | [`Z7Cell`](@ref) |
| [`levels`](@ref) | `0:19` |
| [`has_sorted_subtrees`](@ref) | `true` |
| [`max_neighbors`](@ref) | `6` (pentagons have 5) |
| [`cap_inflation`](@ref) | `1.2` (the interface default) |

Geometry is the Snyder equal-area chart on the icosahedron, shared with the rest
of the ISEA family through the [`ISEA`](@ref) module.

The projection is exactly equal-area, but the published cell is **not** the
chart's equal-area region: [`cell_boundary`](@ref) reports the corner ring, and
those rings tile the sphere exactly while carrying slightly unequal areas
(+1.6% on hexagons and −9.9% on pentagons at level 1, narrowing as cells
shrink). So [`cell_area`](@ref) here is the ring's area — the area of the true
cell, which for IGEO7 is the ring — and the chart's closed form
`4π/(10·7^r)` steradians for a hexagon, `5/6` of that for a pentagon, is a
*different quantity* available separately as
[`equal_area_steradians`](@ref). Reach for that one when you want the
system's nominal equal-area figure, and for `cell_area` when you want the area
of the polygon this package will actually intersect, regrid and draw.

Agreement with DGGRID is pinned by the sealed oracle vectors in
`test/IGeo7/vectors/`: all 196,080 published cell centres at levels 1–5 decode
to their exact Z7 string.
"""
struct IGeo7System <: AbstractHierarchicalGridSystem end

DGG.cellindextype(::IGeo7System) = Z7Cell

DGG.levels(::IGeo7System) = 0:MAX_RESOLUTION

DGG.has_sorted_subtrees(::IGeo7System) = true

# Hexagons and pentagons: vertex adjacency and edge adjacency coincide, so the
# bound is 6 under either connectivity (a pentagon reaches 5).
DGG.max_neighbors(::IGeo7System, ::Connectivity) = 6

# `cap_inflation` keeps the interface default of 1.2. It is not a guess here:
# the covering ratio (a subtree's farthest descendant vertex over the cell's own
# cap radius) was measured exhaustively for depths 1-5 at levels 0-1 and over
# pentagon neighbourhoods at levels 4 and 6, worst case 1.0482 — every
# descendant sits inside 87% of the wired radius, and the ratio converges
# geometrically in two-step (chirality alternates with level parity). The
# conformance covering-law suite re-checks it by sampling.

"""
    rootcells(::IGeo7System) -> SmallVector{12,Z7Cell}

The twelve level-0 cells, one pentagon per icosahedron vertex, ascending. The
level-0 id of base `b` is `(b << 60) | 0x0fff…f`, so ascending id order is
ascending base order and these are exactly positions `1:12` of
`levelgrid(sys, 0)`.
"""
function DGG.rootcells(::IGeo7System)
    out = SmallVector{Z7_NUM_BASES,Z7Cell}()
    for base in 0:(Z7_NUM_BASES-1)
        out = SmallCollections.push(out, Z7Cell((UInt64(base) << Z7_BASE_SHIFT) | Z7_PAD_MASK))
    end
    return out
end

"""
    parent(sys::IGeo7System, c::Z7Cell) -> Z7Cell

The cell one level coarser, by dropping `c`'s last digit — filling its slot, and
every slot below it, with the padding sentinel. Pure bit arithmetic.

Throws an `ArgumentError` for a level-0 root, which has no parent.
"""
function Base.parent(::IGeo7System, c::Z7Cell)
    res = _geometry_checked(c.id)
    res > 0 || throw(ArgumentError(
        "IGeo7 cell $(z7_to_string(c.id)) is a level-0 root and has no parent"))
    return Z7Cell(c.id | _z7_tail_mask(res - 1))
end

"""
    children(sys::IGeo7System, c::Z7Cell) -> SmallVector{7,Z7Cell}

The immediate children of `c`, ascending: **seven** for a hexagon, **six** for a
pentagon, whose base's deleted digit (2 for bases 0–5, 5 for bases 6–11) has no
subtree while the digit prefix is still all zero. Appending a digit is one shift
and one or, so ascending digit order is ascending id order and the result needs
no sort.

Throws an `ArgumentError` at `max_level`, where a child would need a
twenty-first digit slot that has no geometry.
"""
function DGG.children(::IGeo7System, c::Z7Cell)
    res = _geometry_checked(c.id)
    res < MAX_RESOLUTION || throw(ArgumentError(
        "IGeo7 cell $(z7_to_string(c.id)) is at max_level $MAX_RESOLUTION and has no children"))
    out = SmallVector{7,Z7Cell}()
    for z in z7_children(c.id)
        out = SmallCollections.push(out, Z7Cell(z))
    end
    return out
end

"""
    ancestor(sys::IGeo7System, c::Z7Cell, l::Integer) -> Z7Cell

The ancestor of `c` at level `l`, in closed form: one mask that pads every digit
slot below `l`. `ancestor(sys, c, level(c))` is `c` itself.

Throws an `ArgumentError` for `l` outside `0:level(c)`.
"""
function DGG.ancestor(::IGeo7System, c::Z7Cell, l::Integer)
    res = _geometry_checked(c.id)
    target = Int(l)
    0 <= target <= res || throw(ArgumentError(
        "IGeo7 ancestor level must be in 0:$res, got $target"))
    return Z7Cell(c.id | _z7_tail_mask(target))
end

"""
    descendants(sys::IGeo7System, c::Z7Cell, l::Integer) -> Vector{Z7Cell}

Every descendant of `c` at level `l`, ascending. Enumeration is
digit-lexicographic depth-first, which *is* ascending id order, so the result
needs no sort. The count is `7^d` for a hexagon and `(5·7^d + 1)/6` for a
pentagon, `d = l - level(c)`.

Throws an `ArgumentError` for `l` outside `level(c):max_level`. This
materialises — reach for [`descendant_range`](@ref) when positions will do.
"""
function DGG.descendants(::IGeo7System, c::Z7Cell, l::Integer)
    res = _geometry_checked(c.id)
    target = Int(l)
    res <= target <= MAX_RESOLUTION || throw(ArgumentError(
        "IGeo7 descendant level must be in $res:$MAX_RESOLUTION, got $target"))
    out = Vector{Z7Cell}(undef, _subtree_count(c.id, res, target))
    _fill_descendant_cells!(out, 1, c.id, res, target)
    return out
end

function _fill_descendant_cells!(out::Vector{Z7Cell}, pos::Int, z::UInt64, res::Int, target::Int)
    if res == target
        @inbounds out[pos] = Z7Cell(z)
        return pos + 1
    end
    deleted = @inbounds Z7_DELETED_DIGIT[z7_base_cell(z)+1]
    active = Z7_PAD_MASK ⊻ _z7_tail_mask(res)
    pentagon = (z & active) == zero(UInt64)
    shift = _z7_shift(res + 1)
    cleared = z & ~(UInt64(7) << shift)
    for digit in 0:6
        pentagon && digit == deleted && continue
        pos = _fill_descendant_cells!(out, pos, cleared | (UInt64(digit) << shift), res + 1, target)
    end
    return pos
end

"""
    descendant_range(sys::IGeo7System, c::Z7Cell, l::Integer) -> UnitRange{Int}

The contiguous interval of **positions** in `levelgrid(sys, l)` occupied by the
descendants of `c` at level `l`.

Two facts make this O(level) integer work rather than a subtree walk. Ascending
id order is the level's canonical position order, and the level-`l` descendants
of `c` are exactly the ids sharing `c`'s digit prefix — so they are contiguous,
starting at the all-zero-suffix descendant (digit 0 is never a deleted pentagon
digit, so that id is always a valid cell) and running for the subtree's size,
which is `7^d` for a hexagon and `(5·7^d + 1)/6` for a pentagon.

Throws an `ArgumentError` for `l` outside `level(c):max_level`.
"""
function DGG.descendant_range(::IGeo7System, c::Z7Cell, l::Integer)
    res = _geometry_checked(c.id)
    target = Int(l)
    res <= target <= MAX_RESOLUTION || throw(ArgumentError(
        "IGeo7 descendant level must be in $res:$MAX_RESOLUTION, got $target"))
    slots = _z7_tail_mask(res) ⊻ _z7_tail_mask(target)
    first_position = cell_to_index(c.id & ~slots)
    return first_position:(first_position+_subtree_count(c.id, res, target)-1)
end

# ---------------------------------------------------------------------------
# The level grid
# ---------------------------------------------------------------------------

"""
    IGeo7Grid(level) <: AbstractGrid

One complete level of [`IGeo7System`](@ref): all `10·7^level + 2` cells, in the
canonical dense order (ascending Z7 id, which is the space-filling curve order).

A lightweight descriptor — it stores the level and nothing else, so building one
is O(1) and `cellindex` / `cellposition` are O(level) digit walks rather than
table lookups. Obtain one with `levelgrid(IGeo7System(), l)`.
"""
struct IGeo7Grid <: AbstractGrid
    level::Int
end

"""
    levelgrid(sys::IGeo7System, l::Integer) -> IGeo7Grid

The complete IGEO7 grid at level `l`. Throws an `ArgumentError` for `l` outside
`levels(sys)`.
"""
function DGG.levelgrid(::IGeo7System, l::Integer)
    target = Int(l)
    0 <= target <= MAX_RESOLUTION || throw(ArgumentError(
        "IGeo7 level must be in 0:$MAX_RESOLUTION, got $target"))
    return IGeo7Grid(target)
end

DGG.system(::IGeo7Grid) = IGeo7System()
level(g::IGeo7Grid) = g.level

"""
    ncells(g::IGeo7Grid) -> Int

`10·7^level + 2`: twelve pentagons and `10·7^level − 10` hexagons. Fits `Int`
through level 19 (1.14e17).
"""
DGG.ncells(g::IGeo7Grid) = Int(num_cells(g.level))

"""
    cellindex(g::IGeo7Grid, i::Int) -> Z7Cell

The cell at position `i`, by inverting the positional rank walk: peel the base
cell's block, then at each level take the pentagon child while the remainder
fits its subtree and otherwise divide by `7^depth` to pick the hexagon sibling,
re-inserting the deleted digit's gap. O(level), allocation-free.

`i` outside `1:ncells(g)` throws a `BoundsError`.
"""
DGG.cellindex(g::IGeo7Grid, i::Int) = Z7Cell(index_to_cell(i, g.level))

"""
    cellposition(g::IGeo7Grid, c::Z7Cell) -> Union{Int,Nothing}

The position of `c` in the level's dense order, or `nothing` when `c` is not a
cell of this grid — a different level, or an id that is not a valid cell at all.
The walk adds the subtree size of every earlier sibling at each digit, which is
O(level) and needs no table.
"""
function DGG.cellposition(g::IGeo7Grid, c::Z7Cell)
    is_valid_cell(c.id) || return nothing
    _z7_leading_resolution(c.id) == g.level || return nothing
    return cell_to_index(c.id)
end

"""
    cell_boundary(g::IGeo7Grid, c::Z7Cell) -> Vector{UnitSphericalPoint}

The exact boundary ring of `c` on the unit sphere: six corners for a hexagon,
five for a pentagon, **implicitly closed** (the first vertex is not repeated)
and counter-clockwise seen from outside the sphere.

Hexagon corners are the lattice-space midpoint constructions in the cell's
base's dev frame, canonicalised across the cone cut by the cell centre's branch;
pentagon corners are the bisectors between the pentagon's five ring slots. Edges
are straight in the Snyder chart and are reported as their endpoints, which the
package then reads as great-circle arcs.
"""
function DGG.cell_boundary(::IGeo7Grid, c::Z7Cell)
    ring = cell_boundary_cartesian(c.id; closed_ring=false)
    out = Vector{USPoint}(undef, length(ring))
    @inbounds for i in eachindex(ring)
        p = ring[i]
        out[i] = USPoint(p[1], p[2], p[3])
    end
    return out
end

"""
    cell_centroid(g::IGeo7Grid, c::Z7Cell) -> UnitSphericalPoint

The centre of `c` on the unit sphere, strictly interior to the cell: the cell's
physical lattice point pulled back through the Snyder chart, and for a pentagon
exactly the icosahedron vertex it surrounds.

This is the true centre of the cell in the equal-area chart, which is what the
oracle centre dumps publish.
"""
DGG.cell_centroid(::IGeo7Grid, c::Z7Cell) = USPoint(_cell_center_xyz(c.id, _geometry_checked(c.id)))

"""
    cellat(g::IGeo7Grid, p::UnitSphericalPoint) -> Z7Cell

The cell containing `p`, in closed form — no tree descent. The Snyder forward
map picks the containing face, its three corner bases are tried nearest-first
through the dev-frame slot maps, and each candidate is accepted only if a strict
re-encode reproduces the rounded lattice point. Exactly one owner accepts, by
global-lattice consistency, so there is no nearest-centre arbitration.

Never `nothing`: a complete level covers the sphere.

**Ties.** A point exactly on a shared boundary is resolved by the decoder's
rounding-tie fallback, which takes the equally near owner in ascending Voronoi
margin. That is the documented rule, and it is deterministic in the sense the
interface asks for — see [`cellat`](@ref)'s contract, which is per-platform
determinism, not cross-platform bit-identity. IGEO7 does not strengthen it:
the margin comparison is in floating point, so a different CPU or libm may
resolve an exactly-equidistant pair the other way. What holds everywhere is
that the winner is one of the cells genuinely incident to the point.
"""
DGG.cellat(g::IGeo7Grid, p::GO.UnitSphericalPoint) =
    Z7Cell(_xyz_to_z7((Float64(p[1]), Float64(p[2]), Float64(p[3])), g.level))

# `cell_area` is deliberately NOT overridden — see `equal_area_steradians` for
# why the closed form is not the same quantity.

"""
    equal_area_steradians(c::Z7Cell) -> Float64

The **ideal equal-area** solid angle of `c` in steradians, in closed form:
`4π/(10·7^r)` for a hexagon at level `r` and `5/6` of that for a pentagon. These
sum to exactly `4π` over a level, by construction.

This is *not* `cell_area(grid, c)`, and the difference is real rather than
numerical. The Snyder chart is exactly equal-area, so the true IGEO7 cell — the
preimage of a chart hexagon — has exactly this area. But the boundary this
package publishes is the cell's **corner ring**, and the great-circle arcs
between those corners are not the true cell edges: they cut a slightly different
region out of the sphere.

Both regions are honest tessellations. The corner rings tile the sphere exactly
(adjacent cells share corners, so they share whole edges, and their areas sum to
`4π` to full double precision at every level), but they are *not* equal-area:
measured against the closed form, level-1 hexagons run +1.6% and pentagons
−9.9%, and the spread narrows as cells shrink (±2.3% at level 3).

[`cell_area`](@ref) is contractually the area of the **true cell** — and for
IGEO7 the true cell is the published ring, since it is the rings that tile the
sphere and the rings that this package intersects, regrids and draws. So
`cell_area` is left to the generic ring-derived implementation and always
agrees with [`cell_boundary`](@ref) and [`cell_polygon`](@ref); the closed form
here is the *other* quantity, the chart's nominal equal-area figure. (Contrast
HEALPix, where the true cell is an analytic diamond and the published ring is a
densified approximation of it, so its `cell_area` overrides to the closed form
instead. Same contract, opposite conclusion, because the systems differ in
which region is the cell.)

Reach for this function when the equal-area *property*
is what matters — area-weighted statistics, sanity checks against DGGRID's
published areas — and for the ellipsoidal area multiply by the authalic radius
squared (`ISEA.R_AUTHALIC^2`).
"""
function equal_area_steradians(c::Z7Cell)
    res = _geometry_checked(c.id)
    hex = 4pi / (10 * 7.0^res)
    return z7_is_pentagon(c.id) ? 5 * hex / 6 : hex
end

# ---------------------------------------------------------------------------
# Topology
# ---------------------------------------------------------------------------

"""
    neighbors(g::IGeo7Grid, c::Z7Cell, k = 1; connectivity = Vertex()) -> SmallVector{6,Z7Cell}

The cells within `k` adjacency steps of `c`, excluding `c`.

On a hexagonal grid vertex and edge adjacency coincide, so `Vertex()` and
`Edge()` return the same six cells (five for a pentagon) and `max_neighbors` is
6 under either.

# Order

**Rotational**, per the interface contract: rings 1..`k` concatenated outward,
each ring counter-clockwise seen from outside the sphere. So
[`ring`](@ref)`(g, c, k)` is exactly the trailing block of `neighbors(g, c, k)`,
element for element, and the two are computed by one walk so they cannot
disagree.

Ring 1 starts at the dev frame's `+1` reference direction — the six Eisenstein
unit steps in their lattice order, which is the same winding `cell_boundary`
reports its ring in. A pentagon yields five: at the cone apex two of the six
unit directions fold onto the same physical slot, and the duplicate id drops
out.

Rings 2 and outward have no lattice cycle of their own to read — the shell of a
hex disc is not a unit-step orbit — so they are ordered by measurement, the way
the interface docstring prescribes: by azimuth about the cell centroid,
counter-clockwise seen from outside, taking the direction of ring 1's first
entry as the zero. Every ring therefore starts on the same spoke, and exact
azimuth ties break by ascending id.

`k == 0` is empty; `k < 0` throws an `ArgumentError`.

Neighbours are computed by exact lattice arithmetic — one Eisenstein unit step
on the cell's physical lattice point — and the position is turned back into an
id by the same decoder `cellat` uses, whose strict re-encode rejects anything
that is not exactly the cell standing there. Pentagon seams need no special
case.

The container is the static-capacity `SmallVector{6,Z7Cell}` at `k <= 1`, where
the bound is [`max_neighbors`](@ref) and the call does not allocate, and a plain
`Vector{Z7Cell}` above it, where the disc has no static bound.
"""
function DGG.neighbors(g::IGeo7Grid, c::Z7Cell, k::Integer=1;
    connectivity::Connectivity=Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    _level_checked(g, c)
    steps == 0 && return SmallVector{6,Z7Cell}()
    steps == 1 && return _neighbors1(c)
    shells = _shells(g, c, steps)
    isempty(shells) && return Z7Cell[]
    return reduce(vcat, shells)
end

# The k == 1 primitive, in CCW order and in the static-capacity container.
function _neighbors1(c::Z7Cell)
    out = SmallVector{6,Z7Cell}()
    for z in _cell_neighbors_ccw(c.id)
        out = SmallCollections.push(out, Z7Cell(z))
    end
    return out
end

# The breadth-first walk BOTH `neighbors` and `ring` read, so that the disc is
# the shells concatenated by construction rather than by coincidence. Shell `j`
# is the cells at adjacency distance exactly `j`, in the contract's rotational
# order.
function _shells(g::IGeo7Grid, c::Z7Cell, steps::Int)
    shells = Vector{Z7Cell}[]
    steps >= 1 || return shells

    # Ring 1 comes out of the lattice already in CCW order, and its first entry
    # is the zero direction every outer ring is measured against.
    first_ring = collect(_neighbors1(c))
    isempty(first_ring) && return shells
    push!(shells, first_ring)
    reference = DGG.cell_centroid(g, first(first_ring))

    seen = Set{Z7Cell}(first_ring)
    push!(seen, c)
    frontier = first_ring
    for _ in 2:steps
        next = Z7Cell[]
        for x in frontier
            for y in _neighbors1(x)
                y in seen && continue
                push!(seen, y)
                push!(next, y)
            end
        end
        isempty(next) && break
        _sort_ccw!(next, g, c, reference)
        push!(shells, next)
        frontier = next
    end
    return shells
end

"""
    ring(g::IGeo7Grid, c::Z7Cell, k; connectivity = Vertex())

The cells at adjacency distance **exactly** `k`. `ring(g, c, 0)` is `[c]`, and
`ring(g, c, 1)` is [`neighbors`](@ref) at `k == 1`.

Shares [`neighbors`](@ref)' walk, so this is that function's trailing block:
`neighbors(g, c, k)` is `vcat(ring(g, c, 1), ..., ring(g, c, k))`, and the
order contract is the one stated there.
"""
function DGG.ring(g::IGeo7Grid, c::Z7Cell, k::Integer;
    connectivity::Connectivity=Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    _level_checked(g, c)
    steps == 0 && return Z7Cell[c]
    steps == 1 && return _neighbors1(c)
    shells = _shells(g, c, steps)
    # A walk that ran out of cells before reaching `steps` has no shell there:
    # the ring is genuinely empty, not missing.
    steps <= length(shells) || return Z7Cell[]
    return shells[steps]
end

# ---------------------------------------------------------------------------
# Rotational ordering of the outer shells
#
# An orthonormal frame in the tangent plane at the subject cell's centroid, with
# `u` pointing at the reference direction (ring 1's first entry) and `v` chosen
# so that the rotation u -> v is counter-clockwise SEEN FROM OUTSIDE. That is
# `v = p x u` and not `u x p`: for a point `p` on the unit sphere, `p x u`
# leads `u` by a quarter turn in the right-handed sense about the outward
# normal, which is what "counter-clockwise from outside" means.
# ---------------------------------------------------------------------------

function _tangent_frame(centre, toward)
    d = (toward[1] - centre[1], toward[2] - centre[2], toward[3] - centre[3])
    radial = d[1] * centre[1] + d[2] * centre[2] + d[3] * centre[3]
    t = (d[1] - radial * centre[1], d[2] - radial * centre[2],
         d[3] - radial * centre[3])
    n = sqrt(t[1]^2 + t[2]^2 + t[3]^2)
    # A degenerate reference means the zero direction coincides with the centre
    # or its antipode, which a distinct neighbouring cell centre cannot do. The
    # fallback only keeps this total.
    if n <= eps(Float64)
        t = abs(centre[3]) < 0.9 ? (0.0, 0.0, 1.0) : (1.0, 0.0, 0.0)
        n = 1.0
    end
    e1 = (t[1] / n, t[2] / n, t[3] / n)
    e2 = (centre[2] * e1[3] - centre[3] * e1[2],
          centre[3] * e1[1] - centre[1] * e1[3],
          centre[1] * e1[2] - centre[2] * e1[1])
    return e1, e2
end

function _azimuth(centre, e1, e2, p)
    d = (p[1] - centre[1], p[2] - centre[2], p[3] - centre[3])
    a = atan(d[1] * e2[1] + d[2] * e2[2] + d[3] * e2[3],
             d[1] * e1[1] + d[2] * e1[2] + d[3] * e1[3])
    return a < 0 ? a + 2 * Float64(π) : a
end

function _sort_ccw!(shell::Vector{Z7Cell}, g::IGeo7Grid, c::Z7Cell,
        reference)
    length(shell) <= 1 && return shell
    centre = DGG.cell_centroid(g, c)
    e1, e2 = _tangent_frame(centre, reference)
    # Ties by ascending id, per the contract: azimuth is the key, the id is the
    # tiebreak, so the order is total and reproducible.
    sort!(shell; by = z -> (_azimuth(centre, e1, e2, DGG.cell_centroid(g, z)), z))
    return shell
end

# A cell handed to a grid operation must belong to that grid's level; otherwise
# every id below would be silently at the wrong resolution.
@inline function _level_checked(g::IGeo7Grid, c::Z7Cell)
    res = _geometry_checked(c.id)
    res == g.level || throw(ArgumentError(
        "IGeo7 cell $(z7_to_string(c.id)) is at level $res, not this grid's level $(g.level)"))
    return res
end

# ---------------------------------------------------------------------------
# Subtree border
#
# Which descendants of a cell touch a cell outside its subtree, from the Z7
# digits alone — no neighbour query and no geometry. T7 added the generic
# `subtree_border` hook, so this is now IGeo7's method on it rather than a
# module-local name.
# ---------------------------------------------------------------------------

"""
    subtree_border(sys::IGeo7System, c::Z7Cell, l::Integer; connectivity = Vertex()) -> Vector{Z7Cell}

IGeo7's method on the package's [`subtree_border`](@ref) generic.

The descendants of `c` at level `l` that share an edge with a cell outside `c`'s
subtree — the subtree's rim — ascending. `l == level(c)` returns `[c]`, whose
whole neighbourhood is outside its own subtree.

`connectivity` is accepted for the generic's signature and does not change the
answer: on a hexagonal grid sharing a vertex and sharing an edge are the same
relation.

Decided by a six-state automaton over the digit string (a border cell's exposed
directions always form a contiguous arc of the six unit steps, and the arc's
length is exactly its number of border children), so the cost is `O(result)`:
the rim of a hexagon subtree holds `3^(d+1) − 3` cells at depth `d` against the
`7^d` the subtree holds, and a pentagon's holds `5·(3^d − 1)/2`.

Throws an `ArgumentError` for `l` outside `level(c):max_level`.
"""
function subtree_border(::IGeo7System, c::Z7Cell, l::Integer;
        connectivity::Connectivity=Vertex())
    res = _geometry_checked(c.id)
    target = Int(l)
    res <= target <= MAX_RESOLUTION || throw(ArgumentError(
        "IGeo7 descendant level must be in $res:$MAX_RESOLUTION, got $target"))
    return [Z7Cell(z) for z in border_descendants(c.id, target)]
end

"""
    subtree_border_count(sys::IGeo7System, c::Z7Cell, l::Integer) -> Int

The size of [`subtree_border`](@ref) without enumerating it: `3^(d+1) − 3` for a
hexagon subtree and `5·(3^d − 1)/2` for a pentagon one, `d = l - level(c)`.
"""
function subtree_border_count(::IGeo7System, c::Z7Cell, l::Integer)
    res = _geometry_checked(c.id)
    target = Int(l)
    res <= target <= MAX_RESOLUTION || throw(ArgumentError(
        "IGeo7 descendant level must be in $res:$MAX_RESOLUTION, got $target"))
    return Int(_border_count(c.id, res, target))
end
