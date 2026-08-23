# Typed interface over the native Z7 and ISEA geometry layers. Native geometry
# names remain local; package interface methods are qualified with `DGG.`.

"""
    Z7Cell(id::UInt64) <: AbstractCellIndex
    Z7Cell(s::AbstractString)

Canonical Z7 cell identifier. The `UInt64` stores a 4-bit base followed by
twenty 3-bit digit slots; `7` pads slots below the cell's level. Ascending raw
identifiers give canonical lexicographic order and contiguous subtrees.

Unsigned construction does not validate. String construction validates while
parsing; grid and geometry operations validate other identifiers when used.

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

# Attach the shared codec docstring to the second method.
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

This is a non-throwing validity check. Geometry methods instead throw
[`InvalidZ7Error`](@ref) for invalid ids.
"""
is_valid_cell(c::Z7Cell) = is_valid_cell(c.id)

# ---------------------------------------------------------------------------
# The system
# ---------------------------------------------------------------------------

"""
    IGeo7System() <: AbstractHierarchicalGridSystem

An aperture-7 hexagonal hierarchy on the icosahedron using the Snyder ISEA
chart and [`Z7Cell`](@ref) identifiers. Level `r` contains `10·7^r + 2` cells:
twelve pentagons and `10·7^r - 10` hexagons. Hexagons have seven children and
pentagons six; child cells can overhang the parent boundary.

[`cell_area`](@ref) measures the published great-circle corner ring.
[`equal_area_steradians`](@ref) returns the distinct nominal Snyder-chart area.
"""
struct IGeo7System <: AbstractHierarchicalGridSystem end

DGG.cellindextype(::IGeo7System) = Z7Cell

DGG.levels(::IGeo7System) = 0:MAX_RESOLUTION

DGG.has_sorted_subtrees(::IGeo7System) = true

# Hexagons and pentagons: vertex adjacency and edge adjacency coincide, so the
# bound is 6 under either connectivity (a pentagon reaches 5).
DGG.maxneighbors(::IGeo7System, ::Connectivity) = 6

# Z7's one-rings come out of the GBT automaton counter-clockwise, so the shell
# walk carries that order outward instead of measuring azimuth per cell.
DGG.winding(::IGeo7System, ::Connectivity) = DGG.CounterClockwise()

# The default `cap_inflation == 1.2` covers the observed maximum descendant
# overhang ratio of `1.0482` — descendant boundary vertices against an ancestor's
# cell-cap radius. Descendant caps are separate bounds and are not covered.

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

The immediate children of `c`, ascending: seven for a hexagon and six for a
pentagon, whose base's deleted digit (2 for bases 0–5, 5 for bases 6–11) has no
subtree while the digit prefix is still all zero. Appending a digit is one shift
and one or, so ascending digit order is ascending id order and the result needs
no sort.

Throws an `ArgumentError` at `maxlevel`, where a child would need a
twenty-first digit slot that has no geometry.
"""
function DGG.children(::IGeo7System, c::Z7Cell)
    res = _geometry_checked(c.id)
    res < MAX_RESOLUTION || throw(ArgumentError(
        "IGeo7 cell $(z7_to_string(c.id)) is at maxlevel $MAX_RESOLUTION and has no children"))
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

Throws an `ArgumentError` for `l` outside `level(c):maxlevel`. This
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

The contiguous position interval of `c`'s level-`l` descendants. Prefix order
makes this an `O(level)` calculation; the size is `7^d` for a hexagon and
`(5·7^d + 1)/6` for a pentagon.

Throws an `ArgumentError` for `l` outside `level(c):maxlevel`.
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

# Grid descriptor for all `10·7^l + 2` cells in ascending Z7 order.
const LevelGrid = DGG.HierarchicalLevelGrid{IGeo7System}

"""
    ncells(::IGeo7System, l::Integer) -> Int

`10·7^l + 2`: twelve pentagons and `10·7^l − 10` hexagons. Fits `Int` through
level 19 (1.14e17).
"""
DGG.ncells(::IGeo7System, l::Integer) = Int(num_cells(Int(l)))

"""
    cellindex(::IGeo7System, l::Integer, i::Int) -> Z7Cell

The cell at position `i` of level `l`, by inverting the positional rank walk:
peel the base cell's block, then at each level take the pentagon child while the
remainder fits its subtree and otherwise divide by `7^depth` to pick the hexagon
sibling, re-inserting the deleted digit's gap. O(level), allocation-free.
"""
DGG.cellindex(::IGeo7System, l::Integer, i::Int) = Z7Cell(index_to_cell(i, Int(l)))

"""
    cellposition(::IGeo7System, c::Z7Cell) -> Union{Int,Nothing}

The position of `c` in its own level's dense order, or `nothing` when `c` is not
a valid cell at all. The walk adds the subtree size of every earlier sibling at
each digit, which is O(level) and needs no table. The grid has already rejected
a cell from another level.
"""
function DGG.cellposition(::IGeo7System, c::Z7Cell)
    is_valid_cell(c.id) || return nothing
    return cell_to_index(c.id)
end

"""
    cell_boundary(::IGeo7System, c::Z7Cell) -> Helpers.SmallList{6,UnitSphericalPoint}

The exact boundary ring of `c` on the unit sphere: six corners for a hexagon,
five for a pentagon, **implicitly closed** (the first vertex is not repeated)
and counter-clockwise seen from outside the sphere.

Hexagon corners are the lattice-space midpoint constructions in the cell's
base's dev frame, canonicalised across the cone cut by the cell centre's branch;
pentagon corners are the bisectors between the pentagon's five ring slots. Edges
are straight in the Snyder chart and are reported as their endpoints, which the
package then reads as great-circle arcs.
"""
function DGG.cell_boundary(::IGeo7System, c::Z7Cell)
    out = Helpers.empty_small_list(Val(6), USPoint(1.0, 0.0, 0.0))
    for p in cell_boundary_cartesian(c.id; closed_ring=false)
        out = Helpers.small_push(out, USPoint(p[1], p[2], p[3]))
    end
    return out
end

"""
    cell_centroid(::IGeo7System, c::Z7Cell) -> UnitSphericalPoint

The centre of `c` on the unit sphere, strictly interior to the cell: the cell's
physical lattice point pulled back through the Snyder chart, and for a pentagon
exactly the icosahedron vertex it surrounds.
"""
DGG.cell_centroid(::IGeo7System, c::Z7Cell) = USPoint(_cell_center_xyz(c.id, _geometry_checked(c.id)))

"""
    cellat(g::LevelGrid, p::UnitSphericalPoint) -> Z7Cell

Return the cell containing `p` by Snyder projection and strict lattice
re-encoding. Complete levels never return `nothing`. Boundary ties choose an
incident cell by ascending Voronoi margin; exact ties may differ across
floating-point platforms.
"""
DGG.cellat(g::LevelGrid, p::GO.UnitSphericalPoint) =
    Z7Cell(_xyz_to_z7((Float64(p[1]), Float64(p[2]), Float64(p[3])), g.level))

# `cell_area` uses the boundary ring; `equal_area_steradians` reports the
# distinct Snyder-chart area.

"""
    equal_area_steradians(c::Z7Cell) -> Float64

Nominal Snyder-chart solid angle: `4π/(10·7^r)` for a level-`r` hexagon and
`5/6` of that for a pentagon. This differs from [`cell_area`](@ref), which
measures the published great-circle corner ring — a structural gap, not a
numerical one: the rings tile exactly but run +1.6% on level-1 hexagons and
−9.9% on pentagons, narrowing as cells shrink. Multiply by
`ISEA.R_AUTHALIC^2` for ellipsoidal area.
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
    neighbors(g::LevelGrid, c::Z7Cell, k = 1; connectivity = Vertex()) -> SmallVector{6,Z7Cell}

Cells within `k` adjacency steps of `c`, excluding `c`. Vertex and edge
connectivity coincide. Rings are concatenated outward and ordered
counterclockwise; ring 1 starts at the development frame's `+1` direction and
outer rings use the same azimuth reference. Exact azimuth ties use ascending
identifier order.

`k == 0` is empty and `k < 0` throws. For `k <= 1` the result is a
`SmallVector{6,Z7Cell}`; larger discs return `Vector{Z7Cell}`.
"""
function DGG.neighbors(g::LevelGrid, c::Z7Cell, k::Integer=1;
    connectivity::Connectivity=Vertex())
    steps = DGG.checked_steps(k)
    _level_checked(g, c)
    steps == 0 && return SmallVector{6,Z7Cell}()
    steps == 1 && return DGG.one_ring(g, c, connectivity)
    shells = DGG.adjacency_shells(g, c, steps, connectivity)
    isempty(shells) && return Z7Cell[]
    return reduce(vcat, shells)
end

"""
    neighborcount(g::LevelGrid, c::Z7Cell; connectivity = Vertex()) -> Int

Return 5 for pentagons and 6 for other cells, for either connectivity, without
constructing the ring.
"""
function DGG.neighborcount(g::LevelGrid, c::Z7Cell;
        connectivity::Connectivity=Vertex())
    _level_checked(g, c)
    return is_pentagon(c) ? 5 : 6
end

"""
    one_ring(grid, c, connectivity) -> SmallVector{6,Z7Cell}

The immediate neighbours of `c`, counter-clockwise from the development frame's
`+1` direction. Five entries at a pentagon, six elsewhere.
"""
function DGG.one_ring(::LevelGrid, c::Z7Cell, ::Connectivity)
    out = SmallVector{6,Z7Cell}()
    for z in _cell_neighbors_ccw(c.id)
        out = SmallCollections.push(out, Z7Cell(z))
    end
    return out
end

"""
    ring(g::LevelGrid, c::Z7Cell, k; connectivity = Vertex())

The cells at adjacency distance **exactly** `k`. `ring(g, c, 0)` is `[c]`, and
`ring(g, c, 1)` is [`neighbors`](@ref) at `k == 1`.

Shares [`neighbors`](@ref)' walk, so this is that function's trailing block:
`neighbors(g, c, k)` is `vcat(ring(g, c, 1), ..., ring(g, c, k))`, and the
order contract is the one stated there.
"""
function DGG.ring(g::LevelGrid, c::Z7Cell, k::Integer;
    connectivity::Connectivity=Vertex())
    steps = DGG.checked_steps(k)
    _level_checked(g, c)
    steps == 0 && return Z7Cell[c]
    steps == 1 && return DGG.one_ring(g, c, connectivity)
    shells = DGG.adjacency_shells(g, c, steps, connectivity)
    # Return an empty ring after the traversal exhausts the component.
    steps <= length(shells) || return Z7Cell[]
    return shells[steps]
end

# A cell handed to a grid operation must belong to that grid's level; otherwise
# every id below would be silently at the wrong resolution.
@inline function _level_checked(g::LevelGrid, c::Z7Cell)
    res = _geometry_checked(c.id)
    res == g.level || throw(ArgumentError(
        "IGeo7 cell $(z7_to_string(c.id)) is at level $res, not this grid's level $(g.level)"))
    return res
end

# Subtree borders are derived directly from Z7 digits.

function DGG.border_engine(::IGeo7System, c::Z7Cell, target::Int,
        connectivity::Connectivity)
    return Z7BorderEngine(c.id, _z7_subtree_checked(c, target), target)
end

function DGG.interior_engine(::IGeo7System, c::Z7Cell, target::Int,
        connectivity::Connectivity)
    return Z7InteriorEngine(c.id, _z7_subtree_checked(c, target), target)
end

# ---------------------------------------------------------------------------
# Hexagonal halo support
# ---------------------------------------------------------------------------

# The ring position of the step from a cell's parent to the cell, read off the
# cell's own last digit through the same table `_border_step` uses. Digit 0 is
# the centre child, which has no direction, and a base cell has no parent.
#
# RAW `SIGMA_J`, carrying none of the encode rotation `g` that
# `_encode_lattice_rot` documents, because this number is never geometry: the arc
# `_hex_calibrate` builds out of it is handed straight back to `_border_step`,
# which re-reads `SIGMA_J[digit]` for the very same children one call later. The
# raw digit frame is the only frame either side names, so neither owes a
# rotation. Putting one on one side alone would not cancel out on the other —
# `g` follows a cell's own first nonzero digit and its angle against the cone
# cut, so even siblings can disagree on it — and the seeded arc would face where
# none of that neighbour's children lie: a SHORT halo, which is the one way this
# walk answers wrong instead of falling back. `g` belongs where dev-frame order
# is asked for, which is `_cell_neighbors_ccw` and nowhere on this path.
function DGG.hex_child_direction(::IGeo7System, c::Z7Cell)
    res = z7_resolution(c.id)
    res == 0 && return -1
    digit = _z7_digit(c.id, res)
    digit == 0 && return -1
    return @inbounds SIGMA_J[digit]
end

# Unvalidated on purpose: `hex_halo_engine` owns the level guard and only ever
# passes cells that came out of `neighbors`. `_z7_subtree_checked` is the entry
# point for the public verbs, which do not know that.
DGG.seeded_border_engine(::IGeo7System, c::Z7Cell, target::Int, arclen::Int,
        start::Int) = Z7ArcEngine(c.id, z7_resolution(c.id), target,
    Int8(arclen), Int8(start))

# The halo is approached from the neighbouring subtrees rather than from the
# root's own, because a subtree's halo is not an interval of anything Z7 can
# name. See `hex_halo_engine` for the calibration, the containment argument, and
# the guards that send a case back to the generic walk.
DGG.halo_engine(sys::IGeo7System, c::Z7Cell, target::Int,
    connectivity::Connectivity) =
    DGG.hex_halo_engine(sys, c, target, connectivity)

# The one level guard both walks share, and the cell's own resolution.
function _z7_subtree_checked(c::Z7Cell, target::Int)
    res = _geometry_checked(c.id)
    res <= target <= MAX_RESOLUTION || throw(ArgumentError(
        "IGeo7 descendant level must be in $res:$MAX_RESOLUTION, got $target"))
    return res
end

"""
    subtree_border_count(sys::IGeo7System, c::Z7Cell, l::Integer) -> Int

The size of `border(subtree(sys, c, l))` without enumerating it: `3^(d+1) − 3` for a
hexagon subtree and `5·(3^d − 1)/2` for a pentagon one, `d = l - level(c)`.
"""
function subtree_border_count(::IGeo7System, c::Z7Cell, l::Integer)
    res = _geometry_checked(c.id)
    target = Int(l)
    res <= target <= MAX_RESOLUTION || throw(ArgumentError(
        "IGeo7 descendant level must be in $res:$MAX_RESOLUTION, got $target"))
    return Int(_border_count(c.id, res, target))
end
