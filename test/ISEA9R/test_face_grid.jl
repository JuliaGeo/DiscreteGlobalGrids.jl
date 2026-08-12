module Isea9rFaceGridTestSuite

# Tests for `src/ISEA9R/face_grid.jl`: dense ISEA9R grids (all `10 nside²` cells
# of one resolution) as spatial trees, under a swappable data ordering.
#
# This file is the sibling of `test/ISEA4R/test_face_grid.jl`, section for
# section, with the resolutions moved onto powers of three. The three things it
# pins down are the same three, in rising order of load-bearingness:
#
# 1. *Orderings are bijections.* `data_index` / `lattice_index` must be exact
#    mutual inverses over `1:10nside²`, because that pair is the whole alignment
#    contract: column `j` of a `Regridder` is data position `j`. For ISEA9R the
#    Morton half is base-9 (base-3 digit interleave) and exists only at
#    `nside = 3^k`, which `validate_ordering` must enforce at construction.
#
# 2. *Node extents contain their subtree's geometry.* This is what makes tree
#    pruning sound, and for ISEA9R the O(1) four-corner extent `face_grid.jl`
#    ships is justified by *measurement*, not by proof. The measurement is this
#    system's own — `nside ∈ (3, 9, 27)`, the aperture-9 block shapes — not the
#    ISEA4R sibling's result borrowed: the chart is shared but the blocks a
#    cursor builds over a `3^k × 3^k` lattice are not the blocks it builds over
#    a `2^k × 2^k` one. If that assertion ever fires, the documented fix is to
#    switch `cap_policy(::Isea9rFaceSystem)` in `face_grid.jl` to
#    `PerimeterWalkCap()`; the stock default is sound and merely slower.
#
# 3. *The grid regrids conservatively, and the Morton ordering nests.* There is
#    no id-hierarchy *tree* to compare against — `DGGSGrid(ISEA9RDGGS(), level)`
#    needs `has_ordinal_ids`, which is still false — so its place is taken by
#    self-, cross-resolution and cross-system `Regridder` checks. The DGGS
#    *geometry* path (`cell_polygon(ISEA9RDGGS(), level, id)`, wired in
#    `src/ISEA9R/Isea9rKernel.jl`) does exist; a test below smoke-checks that it
#    agrees with this grid, and `test/ISEA9R/test_isea9r_kernel.jl` sweeps it
#    bitwise.
#
# What is NOT here, deliberately: anything about the chart itself. ISEA9R does
# not have one of its own — it imports `ISEA4R`'s, which
# `test/ISEA4R/test_diamonds.jl` pins. `test/ISEA9R/test_delegation.jl` checks
# that the import really is an import and that the two systems agree bitwise.
#
# Naming: `ISEA9R.cell_polygon` is a function in the `ISEA9R` namespace and is
# NOT a method of the top-level `cell_polygon(::AbstractDGGS, level, id)` that
# `using DiscreteGlobalGrids` brings into scope here. Both are used below, so
# the ISEA9R one is always written qualified and the DGGS one bare.
#
# Regridders here are built on the *unit* sphere (`GO.Spherical(radius = 1)`)
# rather than the default Earth-radius manifold, so that conservation reads as
# `4π` and per-entry tolerances are absolute numbers rather than fractions of
# 5.1e14 m².

using Test
using Printf
using Random
import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
import GeometryOps: SpatialTreeInterface as STI
import ConservativeRegridding as CR
import ConservativeRegridding: Trees
import SparseArrays

using DiscreteGlobalGrids
using DiscreteGlobalGrids.ISEA9R
using DiscreteGlobalGrids.HEALPix                 # for the cross-system case
import DiscreteGlobalGrids.ISEA9R as ISEA9R
import DiscreteGlobalGrids.ISEA4R as ISEA4R
using DiscreteGlobalGrids.ISEA9R: DiamondChartGrid, Isea9rFaceRoot,
    num_cells, data_index, lattice_index, validate_ordering, ispow3,
    xyd_to_point, cell_corners, xyd_to_rowmajor, xyd_to_morton, morton_to_xyd

const US = GO.UnitSpherical
# Unit-radius sphere: areas come out in steradians, so conservation is `4π`.
const UNIT = GO.Spherical(radius=1.0)

# Recorded so the numbers land in the test log (and in the milestone report).
# `-Inf` rather than `0.0` as the neutral element: several of the quantities
# below are *negative* by design (a cap overhang inside the cap), and a `0.0`
# floor would quietly report them as zero.
const MEASURED = Dict{String,Float64}()
record!(key, value) = (MEASURED[key] = max(get(MEASURED, key, -Inf), value))

# --------------------------------------------------------------------------
# Helpers — the ISEA4R sibling's, verbatim
# --------------------------------------------------------------------------

cross3(a, b) = (a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])

# Signed area of the 4-gon seen from *outside* the sphere; positive ⇔ CCW.
function ccw_measure(corners)
    acc = (0.0, 0.0, 0.0)
    n = length(corners)
    for i in 1:n
        acc = acc .+ cross3(corners[i], corners[i%n+1])
    end
    outward = reduce((a, b) -> a .+ Tuple(b), corners; init=(0.0, 0.0, 0.0))
    return sum(acc .* outward)
end

ring_points(poly) = collect(GI.getpoint(GI.getexterior(poly)))
# `getcell` closes the ring, so the four distinct corners are points 1:4.
open_ring(poly) = ring_points(poly)[1:4]

"Depth-first walk over every node of a diamond cursor (the node itself included)."
function walk_nodes(f::Function, node)
    f(node)
    STI.isleaf(node) && return nothing
    for child in STI.getchild(node)
        walk_nodes(f, child)
    end
    return nothing
end

"""
Every `(data index, leaf extent)` pair the traversal can ever test, harvested
from `STI.child_indices_extents` — i.e. the literal leaf-level predicate input,
not a reconstruction of it.
"""
function leaf_extents(root)
    exts = Dict{Int,US.SphericalCap{Float64}}()
    for d in 1:STI.nchild(root)
        walk_nodes(STI.getchild(root, d)) do node
            STI.isleaf(node) || return nothing
            for (index, extent) in STI.child_indices_extents(node)
                exts[index] = extent
            end
            return nothing
        end
    end
    return exts
end

"Random unit vector / random cap, from a seeded RNG."
function random_point(rng)
    v = randn(rng, 3)
    v ./= sqrt(sum(abs2, v))
    return GO.UnitSphericalPoint(v[1], v[2], v[3])
end
random_cap(rng) = US.SphericalCap(random_point(rng), 0.01 + 1.2 * rand(rng))

"Is any point of cell `(ix, iy, diamond)` inside `cap`? (dense chart sampling)"
function cell_meets_cap(cap, ix, iy, diamond, nside, samples=4)
    for a in 0:samples, b in 0:samples
        p = xyd_to_point((ix + a / samples) / nside, (iy + b / samples) / nside, diamond)
        US.spherical_distance(cap.point, p) <= cap.radius && return true
    end
    return false
end

"""
Every index-range block a `TopDownQuadtreeCursor` over an `nside × nside` grid
can ever hold: the range-bisection recursion from `(1:n, 1:n)` down, mirroring
`STI.getchild`, plus every `1×1` leaf range (the ranges
`STI.child_indices_extents` builds per-cell caps from).
"""
function cursor_blocks(n)
    out = Tuple{UnitRange{Int},UnitRange{Int}}[]
    function rec(ir, jr)
        push!(out, (ir, jr))
        all(length.((ir, jr)) .<= 2) && return nothing
        li, lj = length(ir), length(jr)
        if li == 1
            s = lj ÷ 2
            rec(ir, jr[1:s])
            rec(ir, jr[(s + 1):end])
        elseif lj == 1
            s = li ÷ 2
            rec(ir[1:s], jr)
            rec(ir[(s + 1):end], jr)
        else
            si, sj = li ÷ 2, lj ÷ 2
            rec(ir[1:si], jr[1:sj])
            rec(ir[1:si], jr[(sj + 1):end])
            rec(ir[(si + 1):end], jr[1:sj])
            rec(ir[(si + 1):end], jr[(sj + 1):end])
        end
        return nothing
    end
    rec(1:n, 1:n)
    for i in 1:n, j in 1:n
        push!(out, (i:i, j:j))
    end
    return unique(out)
end

"""
The radius `circle_from_four_corners` would use *before* its `1.0001` slack: the
farthest of the four corners and their four slerp midpoints from the cap centre.
Recomputed rather than obtained as `cap.radius / 1.0001`, which would carry the
ulp of the multiply-then-divide round trip into a quantity compared against `0`.
"""
function preslack_radius(center, bl, tl, br, tr)
    p1, p2, p3, p4 = bl, br, tr, tl            # the CCW cycle `_spherical_cap` uses
    mids = (US.slerp(p1, p2, 0.5), US.slerp(p2, p3, 0.5),
        US.slerp(p3, p4, 0.5), US.slerp(p4, p1, 0.5))
    return maximum(US.spherical_distance(center, p) for p in (p1, p2, p3, p4, mids...))
end

# --------------------------------------------------------------------------
# 1. Construction and validation
# --------------------------------------------------------------------------

@testset "Isea9rFaceSpace / Isea9rFaceGrid construction" begin
    @test_throws ArgumentError Isea9rFaceSpace(0)
    @test_throws ArgumentError Isea9rFaceSpace(-3)
    @test Isea9rFaceSpace(1).nside == 1
    @test num_cells(Isea9rFaceSpace(3)) == 90
    @test num_cells(Isea9rFaceSpace(9)) == 810

    # Upper bound: `3^18`, the largest power of three that keeps `10 * nside^2`
    # inside `Int64` *and* inside the domain the borrowed ISEA4R chart's
    # seam-exactness argument covers (`nside <= 2^29`). Deliberately NOT the
    # sibling's `2^29`: a delegating system may not outrun the argument it
    # delegates to, and `3^k` is the only resolution at which ISEA9R has an id
    # space at all. See the derivation at `max_nside(::Isea9rFaceSystem)`.
    @test Isea9rFaceSpace(3^18).nside == 3^18
    @test_throws ArgumentError Isea9rFaceSpace(3^18 + 1)
    @test_throws ArgumentError Isea9rFaceSpace(3^19)
    @test_throws ArgumentError Isea9rFaceSpace(2^29)     # the ISEA4R bound is bigger
    @test 10 * Int64(3^18)^2 < typemax(Int64)
    @test 10 * Int128(3^19)^2 > typemax(Int64)
    @test 3^18 < 2^29
    # ...and the error message names the bound in the form it is documented in.
    @test occursin("3^18", sprint(showerror, try
        Isea9rFaceSpace(3^19)
    catch e
        e
    end))

    # The base-9 Morton code interleaves one base-3 digit per level, so it
    # exists only on a 3^k x 3^k diamond. It must be refused at *construction*,
    # not from inside a traversal.
    for nside in (2, 4, 5, 6, 8, 10, 16)
        @test !ispow3(nside)
        @test_throws ArgumentError Isea9rFaceGrid(nside; ordering=MortonOrder())
    end
    @test_throws ArgumentError Isea9rFaceGrid(Isea9rFaceSpace(6), MortonOrder())
    @test_throws ArgumentError Isea9rFaceGrid(0)
    @test_throws ArgumentError Isea9rFaceGrid(0; ordering=MortonOrder())
    @test ispow3(1) && ispow3(3) && ispow3(9) && ispow3(27) && ispow3(3^18)
    @test !ispow3(0) && !ispow3(-3)

    # Row-major carries no such restriction — that is the point of this layer.
    for nside in (1, 2, 3, 4, 5, 7, 9)
        grid = Isea9rFaceGrid(nside)
        @test grid.ordering isa RowMajorOrder          # the default
        @test grid.space.nside == nside
        @test num_cells(grid) == 10nside^2
    end
    @test Isea9rFaceGrid(9; ordering=MortonOrder()).ordering isa MortonOrder
    @test RowMajorOrder() isa AbstractIsea9rOrdering
    @test MortonOrder() isa AbstractIsea9rOrdering
    # The optional third contract method defaults to "accept everything". It
    # takes the *space*, not a loose `nside`: the resolution reaches this layer
    # only through the type that has already checked it.
    @test validate_ordering(RowMajorOrder(), Isea9rFaceSpace(4)) === nothing
    @test validate_ordering(MortonOrder(), Isea9rFaceSpace(9)) === nothing
    @test_throws ArgumentError validate_ordering(MortonOrder(), Isea9rFaceSpace(4))

    @test sprint(show, Isea9rFaceSpace(9)) == "Isea9rFaceSpace(9)"
    @test occursin("nside=9", sprint(show, Isea9rFaceGrid(9)))
    @test occursin("810 cells", sprint(show, Isea9rFaceGrid(9)))
end

# The ISEA9R orderings are *distinct types* from the ISEA4R ones even where the
# arithmetic behind them is shared, so that a grid of one system cannot be
# built with the other's ordering — the check that keeps "radix 9" from
# silently meaning "radix 4".
@testset "ISEA9R and ISEA4R orderings do not cross" begin
    @test !(ISEA9R.RowMajorOrder() isa ISEA4R.AbstractIsea4rOrdering)
    @test !(ISEA4R.RowMajorOrder() isa AbstractIsea9rOrdering)
    @test ISEA9R.RowMajorOrder() !== ISEA4R.RowMajorOrder()
    @test typeof(ISEA9R.MortonOrder()) !== typeof(ISEA4R.MortonOrder())

    @test_throws ArgumentError Isea9rFaceGrid(9; ordering=ISEA4R.MortonOrder())
    @test_throws ArgumentError Isea9rFaceGrid(9; ordering=ISEA4R.RowMajorOrder())
    @test_throws ArgumentError ISEA4R.Isea4rFaceGrid(8; ordering=ISEA9R.MortonOrder())
    @test_throws ArgumentError Isea9rFaceGrid(9; ordering=DiscreteGlobalGrids.S2.RowMajorOrder())
    @test_throws ArgumentError Isea9rFaceGrid(9; ordering=DiscreteGlobalGrids.HEALPix.NestedOrder())
end

# `Isea9rFaceRoot` is independently constructible — it is the REPL/test
# convenience form — so it must not be a way *around* the checks
# `Isea9rFaceGrid` runs. Its resolution is an `Isea9rFaceSpace` (so `nside` is
# checked by construction) and it re-runs both ordering checks itself.
@testset "Isea9rFaceRoot construction validates like the grid" begin
    # (a) `nside`, through the space the convenience form builds.
    @test_throws ArgumentError Isea9rFaceRoot(0)
    @test_throws ArgumentError Isea9rFaceRoot(-5, RowMajorOrder())
    @test_throws ArgumentError Isea9rFaceRoot(3^19, RowMajorOrder())
    @test Isea9rFaceRoot(3^18, RowMajorOrder()).space.nside == 3^18

    # (b) another system's ordering: an `ArgumentError` here rather than a
    # `MethodError` from inside a traversal.
    @test_throws ArgumentError Isea9rFaceRoot(9, ISEA4R.RowMajorOrder())
    @test_throws ArgumentError Isea9rFaceRoot(
        GO.Spherical(), Isea9rFaceSpace(9), DiscreteGlobalGrids.HEALPix.NestedOrder())

    # (c) `validate_ordering`, which now runs at root construction too.
    @test_throws ArgumentError Isea9rFaceRoot(4, MortonOrder())
    @test_throws ArgumentError Isea9rFaceRoot(GO.Spherical(), Isea9rFaceSpace(5), MortonOrder())
    @test Isea9rFaceRoot(9, MortonOrder()).ordering isa MortonOrder
end

@testset "treeify and tree toplevel" begin
    for (nside, ordering) in ((9, MortonOrder()), (5, RowMajorOrder()),
                              (3, MortonOrder()), (4, RowMajorOrder()))
        grid = Isea9rFaceGrid(nside; ordering)
        # One-argument `treeify` resolves through `best_manifold`.
        @test GOCore.best_manifold(grid) == GO.Spherical()
        root = treeify(grid)
        @test root isa Isea9rFaceRoot
        @test root.space == Isea9rFaceSpace(nside)
        @test root.space.nside == nside
        @test root.ordering === ordering
        @test GOCore.best_manifold(root) == GO.Spherical()
        # Explicit-manifold form, and idempotent passthrough on the tree.
        @test treeify(GO.Spherical(), grid) isa Isea9rFaceRoot
        @test treeify(root) === root
        @test treeify(GO.Spherical(), root) === root
        @test Trees.treeify(UNIT, root) === root
        # The direct constructor defaults the manifold for REPL use.
        @test Isea9rFaceRoot(nside, ordering) == root

        @test STI.isspatialtree(typeof(root))
        @test !STI.isleaf(root)
        @test STI.nchild(root) == 10
        @test Trees.ncells(root) == 10nside^2
        @test length(collect(STI.getchild(root))) == 10
        @test occursin("nside=$nside", sprint(show, root))

        # The root bounds nothing tighter than the whole sphere: the ten
        # diamonds tile it.
        @test STI.node_extent(root).radius >= Float64(pi)

        # Children are stock quadtree cursors over per-diamond chart grids.
        for d in 1:10
            child = STI.getchild(root, d)
            @test child isa Trees.TopDownQuadtreeCursor{<:DiamondChartGrid}
            @test child.grid.face == d - 1
            @test child.grid.space == Isea9rFaceSpace(nside)
            @test Trees.ncells(child.grid, 1) == nside
            @test Trees.ncells(child.grid, 2) == nside
            @test GOCore.manifold(child.grid) == GO.Spherical()
        end

        @test_throws BoundsError Trees.getcell(root, 0)
        @test_throws BoundsError Trees.getcell(root, 10nside^2 + 1)
        @test length(collect(Trees.getcell(root))) == 10nside^2
    end
end

# --------------------------------------------------------------------------
# 2. The ordering contract
# --------------------------------------------------------------------------

@testset "ordering bijections (nside = $nside)" for nside in (1, 2, 3, 4, 5, 9)
    ncell = 10nside^2
    space = Isea9rFaceSpace(nside)
    orderings = ispow3(nside) ? (RowMajorOrder(), MortonOrder()) : (RowMajorOrder(),)
    for ordering in orderings
        # `lattice_index` is onto the whole lattice and `data_index` inverts it.
        seen = Set{NTuple{3,Int}}()
        for j in 1:ncell
            ix, iy, diamond = lattice_index(ordering, space, j)
            @test 0 <= ix < nside && 0 <= iy < nside && 0 <= diamond <= 9
            push!(seen, (ix, iy, diamond))
            @test data_index(ordering, space, ix, iy, diamond) == j
        end
        @test length(seen) == ncell

        # ... and the other direction, over the lattice.
        for diamond in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
            j = data_index(ordering, space, ix, iy, diamond)
            @test 1 <= j <= ncell
            @test lattice_index(ordering, space, j) == (ix, iy, diamond)
        end
    end

    # The two shipped orderings are exactly the chart's two index maps, with the
    # single 0-based/1-based reconciliation both need.
    ix, iy = min(1, nside - 1), min(2, nside - 1)
    @test data_index(RowMajorOrder(), space, ix, iy, 3) ==
          xyd_to_rowmajor(ix, iy, 3, nside) + 1
    if ispow3(nside)
        @test data_index(MortonOrder(), space, ix, iy, 3) ==
              xyd_to_morton(ix, iy, 3, nside) + 1
        @test lattice_index(MortonOrder(), space, 7) == morton_to_xyd(6, nside)
    end
end

# The codecs' own guard rails, in the direction a user-supplied id enters
# through. `nside` that is not `3^k` is refused by BOTH codecs, which is what
# keeps a base-9 code from being read at a resolution where it means nothing.
@testset "base-9 codec range and resolution guards" begin
    for nside in (2, 4, 5, 8)
        @test_throws ArgumentError xyd_to_morton(0, 0, 0, nside)
        @test_throws ArgumentError morton_to_xyd(0, nside)
    end
    for nside in (1, 3, 9)
        ncell = 10nside^2
        @test_throws ArgumentError morton_to_xyd(-1, nside)
        @test_throws ArgumentError morton_to_xyd(ncell, nside)
        @test morton_to_xyd(ncell - 1, nside) == (nside - 1, nside - 1, 9)
        @test_throws ArgumentError xyd_to_morton(nside, 0, 0, nside)
        @test_throws ArgumentError xyd_to_morton(0, -1, 0, nside)
        @test_throws ArgumentError xyd_to_morton(0, 0, 10, nside)
        @test_throws ArgumentError xyd_to_morton(0, 0, -1, nside)
    end
end

# `MortonOrder` is written against the ordinal the `ISEA9RDGGS` registry entry
# records: `diamond * 9^level + position` at `nside = 3^level`, with `position`
# pinned here to the base-9 Morton (Z-order) code. What is pinned here is the
# ordinal arithmetic alone; the geometry comparison against the DGGS kernel over
# the same ordinals is `test/ISEA9R/test_isea9r_kernel.jl`.
@testset "MortonOrder position j is isea9r_ordinal j - 1" begin
    # An independent base-3 digit interleave, so this is not the source checking
    # itself: base-9 digit `k` of the code is `ix_k + 3 * iy_k`.
    interleave(ix, iy) = sum((((ix ÷ 3^k) % 3) + 3 * ((iy ÷ 3^k) % 3)) * 9^k for k in 0:11)
    for level in 0:2
        nside = 3^level
        for diamond in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
            j = data_index(MortonOrder(), Isea9rFaceSpace(nside), ix, iy, diamond)
            ordinal = j - 1
            position = ordinal - diamond * 9^level
            @test 0 <= position < 9^level
            @test position == interleave(ix, iy)
            @test ordinal == diamond * 9^level + position
            @test morton_to_xyd(ordinal, nside) == (ix, iy, diamond)
        end
        @test root_count(ISEA9RDGGS()) == 10            # the 10 in `10 * 9^level`
        @test radix(ISEA9RDGGS()) == 9
        @test supports_prefix_ranges(ISEA9RDGGS())
    end
end

# The radix-9 prefix property, at the level of the ordinal itself: the nine
# children of ordinal `p` are `9p:9p+8`, and they are the nine lattice cells
# that subdivide `p`'s. This is the arithmetic
# `supports_prefix_ranges(ISEA9RDGGS()) == true` is a claim about.
@testset "base-9 Morton ordinals nest: children of p are 9p:9p+8" begin
    for level in 0:1
        n = 3^level
        fine = 3n
        for p in 0:(10 * n^2 - 1)
            ix, iy, diamond = morton_to_xyd(p, n)
            kids = [morton_to_xyd(9p + k, fine) for k in 0:8]
            @test all(k -> k[3] == diamond, kids)
            # Exactly the 3x3 block of the parent's lattice cell, once each.
            @test sort([(k[1], k[2]) for k in kids]) ==
                  sort([(3ix + a, 3iy + b) for a in 0:2, b in 0:2] |> vec)
            # ...and the interval is what `leaf_interval` would compute.
            @test leaf_interval(ISEA9RDGGS(), level, p, level + 1) == (9p):(9p + 8)
        end
    end
    @test leaf_count(ISEA9RDGGS(), 2) == 10 * 9^2
    @test leaf_interval(ISEA9RDGGS(), 0, 9, 2) == (9 * 81):(10 * 81 - 1)
end

@testset "getcell round-trips through the ordering (nside = $nside, $ordering)" for
        (nside, ordering) in ((9, MortonOrder()), (9, RowMajorOrder()),
                              (3, MortonOrder()), (5, RowMajorOrder()))

    space = Isea9rFaceSpace(nside)
    root = treeify(Isea9rFaceGrid(nside; ordering))
    ncell = 10nside^2
    for j in 1:ncell
        ix, iy, diamond = lattice_index(ordering, space, j)
        # `Trees.getcell(root, j)` is the geometry side of the alignment rule.
        @test ring_points(Trees.getcell(root, j)) ==
              ring_points(ISEA9R.cell_polygon(ix, iy, diamond, nside))
        # Closed CCW 4-gon, with the corners the chart emits, in its order.
        pts = ring_points(Trees.getcell(root, j))
        @test length(pts) == 5 && pts[1] == pts[5]
        @test Tuple(pts[1:4]) == cell_corners(ix, iy, diamond, nside)
    end
    # The iterator form must agree with the indexed one, position by position.
    @test collect(Trees.getcell(root)) == [Trees.getcell(root, j) for j in 1:ncell]

    # The per-diamond grid's own index maps target the *global* data layout.
    for diamond in 0:9
        g = DiamondChartGrid(GO.Spherical(), space, diamond, ordering)
        for ix in 0:(nside - 1), iy in 0:(nside - 1)
            j = Trees.cartesian_to_linear_idx(g, CartesianIndex(ix + 1, iy + 1))
            @test j == data_index(ordering, space, ix, iy, diamond)
            @test Trees.linear_to_cartesian_idx(g, j) == CartesianIndex(ix + 1, iy + 1)
            @test ring_points(Trees.getcell(g, ix + 1, iy + 1)) ==
                  ring_points(Trees.getcell(root, j))
            @test Trees.getvertex(g, ix + 1, iy + 1) ==
                  xyd_to_point(ix / nside, iy / nside, diamond)
        end
    end
end

# --------------------------------------------------------------------------
# 3. Node extents contain their subtree's geometry
#
# The soundness check for pruning, by exhaustion, plus the pre-registered
# measurement the O(1) override rests on.
# --------------------------------------------------------------------------

@testset "node extents contain their block geometry (nside = $nside)" for nside in (3, 4, 5, 9)
    root = treeify(Isea9rFaceGrid(nside; ordering=RowMajorOrder()))
    worst_vertex = -Inf
    worst_sample = -Inf
    nodes = 0
    for d in 1:10
        cursor = STI.getchild(root, d)
        grid = cursor.grid
        walk_nodes(cursor) do node
            nodes += 1
            extent = STI.node_extent(node)
            imin, imax = extrema(node.leafranges[1]); imax += 1
            jmin, jmax = extrema(node.leafranges[2]); jmax += 1
            for i in imin:imax, j in jmin:jmax
                p = Trees.getvertex(grid, i, j)
                worst_vertex = max(worst_vertex,
                    US.spherical_distance(extent.point, p) - extent.radius)
            end
            for a in 0:8, b in 0:8
                x = ((imin - 1) + (imax - imin) * a / 8) / nside
                y = ((jmin - 1) + (jmax - jmin) * b / 8) / nside
                p = xyd_to_point(x, y, grid.face)
                worst_sample = max(worst_sample,
                    US.spherical_distance(extent.point, p) - extent.radius)
            end
            return nothing
        end
    end
    @test nodes > 10                                   # the walk actually descended
    # Non-positive overhang: every probe sits inside the cap.
    @test worst_vertex <= 0
    @test worst_sample <= 0
    record!("node-extent overhang, lattice vertices (rad)", worst_vertex)
    record!("node-extent overhang, dense chart samples (rad)", worst_sample)

    # Leaf extents likewise contain their own cell, and the leaf index space is
    # exactly `1:10nside²` with no index emitted twice.
    exts = leaf_extents(root)
    @test sort!(collect(keys(exts))) == collect(1:(10nside^2))
    worst_leaf = -Inf
    for (j, extent) in exts
        for p in open_ring(Trees.getcell(root, j))
            worst_leaf = max(worst_leaf,
                US.spherical_distance(extent.point, p) - extent.radius)
        end
    end
    @test worst_leaf <= 0
    record!("leaf-extent overhang, cell corners (rad)", worst_leaf)
end

# The standing justification for the O(1) four-corner `STI.node_extent` override
# in `face_grid.jl`. Pre-registered decision rule, recorded before the numbers
# were in and identical to the ISEA4R sibling's except for the resolutions:
# adopt the override if the worst PRE-SLACK overhang over every block the cursor
# can create is `<= 0` (or if the worst ratio is a stable `<= 5%`). The
# resolutions are this system's own — `3, 9, 27` are the aperture-9 block shapes
# — because the sibling's measurement is about `2^k × 2^k` blocks and says
# nothing about these.
#
# Measured: exactly `0.0` at every `nside` (13, 117 and 1045 blocks per diamond)
# — the farthest sampled point of every block is one of its own four corners,
# seam-straddling blocks and whole-diamond roots included. The assertion below
# is that measurement, kept as a test so a regression in the chart or in
# `circle_from_four_corners` cannot silently invalidate the override.
@testset "four corners are extremal on every cursor block (pre-slack)" begin
    worst_overhang = -Inf
    worst_ratio = -Inf
    for nside in (3, 9, 27)
        blocks = cursor_blocks(nside)
        for diamond in 0:9, (ir, jr) in blocks
            imin, imax = extrema(ir); imax += 1
            jmin, jmax = extrema(jr); jmax += 1
            corner(i, j) = xyd_to_point((i - 1) / nside, (j - 1) / nside, diamond)
            bl, tl = corner(imin, jmin), corner(imin, jmax)
            br, tr = corner(imax, jmin), corner(imax, jmax)
            cap = Trees.circle_from_four_corners((bl, tl, br, tr), ())
            pre = preslack_radius(cap.point, bl, tl, br, tr)
            worst = -Inf
            for a in 0:16, b in 0:16
                x = ((imin - 1) + (imax - imin) * a / 16) / nside
                y = ((jmin - 1) + (jmax - jmin) * b / 16) / nside
                worst = max(worst, US.spherical_distance(cap.point, xyd_to_point(x, y, diamond)))
            end
            worst_overhang = max(worst_overhang, worst - pre)
            worst_ratio = max(worst_ratio, (worst - pre) / pre)
        end
    end
    @test worst_overhang <= 0
    @test worst_ratio <= 0
    record!("pre-slack 4-corner overhang over all cursor blocks (rad)", worst_overhang)
    record!("pre-slack 4-corner overhang over all cursor blocks (ratio)", worst_ratio)
end

# --------------------------------------------------------------------------
# 4. Cap queries against brute force
# --------------------------------------------------------------------------

@testset "STI queries vs brute force (nside = $nside, $ordering)" for
        (nside, ordering) in ((9, MortonOrder()), (3, MortonOrder()), (5, RowMajorOrder()))

    root = treeify(Isea9rFaceGrid(nside; ordering))
    ncell = 10nside^2
    exts = leaf_extents(root)
    lattice = [lattice_index(ordering, Isea9rFaceSpace(nside), j) for j in 1:ncell]

    rng = MersenneTwister(20260806 + nside)
    pruned_total = 0
    leaf_positive_total = 0
    tested_total = 0
    for _ in 1:30
        cap = random_cap(rng)
        # Instrumented predicate: `tested_total` counts every extent the
        # traversal actually evaluates, internal blocks included.
        answer = sort!(STI.query(root, function (extent)
            tested_total += 1
            return intersects_cap(cap, extent)
        end))

        @test allunique(answer)
        @test all(j -> 1 <= j <= ncell, answer)

        # Upper bound: the traversal only ever *returns* leaves whose own leaf
        # extent passes the predicate, so no result can be outside this set.
        leaf_positives = [j for j in 1:ncell if intersects_cap(cap, exts[j])]
        @test issubset(answer, leaf_positives)

        # Lower bound — the load-bearing half. Any cell that genuinely meets the
        # cap must survive every ancestor's pruning test. (Ground truth by dense
        # chart sampling of the cell: a sampled point inside the cap is a
        # witness that the cell really intersects it.)
        truth = [j for j in 1:ncell if cell_meets_cap(cap, lattice[j]..., nside)]
        @test issubset(truth, answer)

        leaf_positive_total += length(leaf_positives)
        pruned_total += length(leaf_positives) - length(answer)
    end
    # Internal pruning has to actually happen: a regression that inflated the
    # block caps to trivial full-sphere extents would still satisfy both subset
    # assertions above (it only ever *adds* leaves to the answer) and would
    # silently pass the rest of this testset. Under that regression every
    # internal node passes, so the traversal degenerates to testing every node
    # of every diamond tree — strictly MORE extents than there are leaves
    # (`30 * ncell`). Real pruning comes in comfortably under that.
    @test tested_total < 30 * ncell

    # Informational: how much of the leaf-extent-positive set internal pruning
    # removes. NOT asserted positive — a leaf whose own cap meets the query cap
    # while an ancestor's does not is a rare accident of cap slack, not an
    # invariant. `tested_total` above is the load-bearing one.
    record!("leaf-extent positives pruned by ancestors (fraction)",
        pruned_total / max(leaf_positive_total, 1))
end

# --------------------------------------------------------------------------
# 5. Regridder structure
#
# The place `test/HEALPix/test_face_grid.jl` compares the face-grid path against
# the id-hierarchy path. ISEA9R has no id-hierarchy *grid* (that needs
# `has_ordinal_ids`), so (a) checks what it does have — the DGGS geometry path
# over the same ordinals — and records what is still missing, while (b)-(e)
# stand in for the matrix comparison.
# --------------------------------------------------------------------------

@testset "the DGGS geometry path agrees; the ISEA9R id hierarchy is not wired" begin
    # `cell_polygon` here is the top-level `cell_polygon(::AbstractDGGS, level,
    # id)` — NOT `ISEA9R.cell_polygon(ix, iy, diamond, nside)`, which is what
    # the rest of this file exercises. It used to throw `NotPortedError`; since
    # `src/ISEA9R/Isea9rKernel.jl` it answers over the canonical
    # `isea9r_ordinal`, which is exactly this grid's `MortonOrder` data position
    # minus one. A smoke check of that correspondence lives here; the full
    # bitwise sweep, the caps and the still-unwired hierarchy group are
    # `test/ISEA9R/test_isea9r_kernel.jl`'s business.
    root = treeify(Isea9rFaceGrid(9; ordering=MortonOrder()))
    for id in (0, 17, 10 * 81 - 1)
        @test ring_points(cell_polygon(ISEA9RDGGS(), 2, id)) ==
              ring_points(Trees.getcell(root, id + 1))
    end
    # What is still NOT wired: the hierarchy/ordinal/pruning group (deferred,
    # not blocked — the radix-9 arithmetic over these ordinals is exact, which
    # is what `supports_prefix_ranges` asserts at the interface level).
    @test_throws NotPortedError cell_children(ISEA9RDGGS(), 0, 0)
    @test_throws NotPortedError descendant_range(ISEA9RDGGS(), 0, 0, 1)
    @test max_level(ISEA9RDGGS()) === nothing          # unbounded level
    @test is_equal_area(ISEA9RDGGS())
    @test aperture(ISEA9RDGGS()) == 9
    @test canonical_index_name(ISEA9RDGGS()) === :isea9r_ordinal
    # The registry facts this milestone flipped, next to the geometry they
    # describe: ten root rhombuses (OGC 21-038r1 Annex B.2; DGGAL RI9R.ec
    # `countZones(0) == 10`), radix 9, prefix ranges over the package ordinal.
    @test root_count(ISEA9RDGGS()) == 10
    @test supports_prefix_ranges(ISEA9RDGGS())
end

@testset "self-regridding is diagonal (nside = 3)" begin
    nside = 3
    ncell = 10nside^2
    tree = treeify(Isea9rFaceGrid(nside; ordering=RowMajorOrder()))
    R = CR.Regridder(UNIT, tree, tree; threaded=false, normalize=false)

    @test size(R.intersections) == (ncell, ncell)
    # Neighbouring cells share bit-identical chart corners, so their
    # intersection is a zero-area sliver and is dropped: the matrix is exactly
    # diagonal. (The same holds for HEALPix, S2 and ISEA4R.)
    @test SparseArrays.nnz(R.intersections) == ncell
    offdiag = sum(abs, R.intersections) - sum(abs, SparseArrays.diag(R.intersections))
    @test offdiag <= 1e-12
    record!("self-regrid off-diagonal mass (sr)", offdiag)
    @test all(i -> R.intersections[i, i] > 0, 1:ncell)
    @test SparseArrays.diag(R.intersections) ≈ R.dst_areas
    @test isapprox(sum(R.intersections), 4π; rtol=1e-10)
    record!("|sum(intersections) - 4π|, self nside=3", abs(sum(R.intersections) - 4π))
end

# A found-in-the-building defect, recorded rather than hidden. It is NOT an
# ISEA9R defect: it reproduces exactly on the committed ISEA4R grid at
# `Isea4rFaceGrid(9)` and `Isea4rFaceGrid(16)`, with no ISEA9R code loaded, and
# the minimal reproduction bypasses this package's tree layer entirely —
#
#   p1 = ISEA4R.cell_polygon(8, 7, 9, 9); p2 = ISEA4R.cell_polygon(0, 0, 5, 9)
#   GO.area(UNIT, GO.intersection(GO.ConvexConvexSutherlandHodgman(UNIT), p1, p2;
#                                 target = GI.PolygonTrait()))   # 0.0, correct
#   GO.area(UNIT, GO.intersection(GO.ConvexConvexSutherlandHodgman(UNIT), p2, p1;
#                                 target = GI.PolygonTrait()))   # 0.0166..., wrong
#
# — i.e. the spherical convex clipper is argument-order dependent on one
# edge-adjacent pair of convex cells and returns a whole cell as the
# intersection of two disjoint ones. Both cells are convex (all four spherical
# turns positive) and they meet only along a cross-diamond border, so nothing in
# this package's contract is violated; the fix belongs upstream in
# `GeometryOps.ConvexConvexSutherlandHodgman`. See
# `docs/design/isea9r_layout.md` §8.
#
# What this testset asserts is therefore the invariants that DO hold at
# `nside = 9`, plus a ceiling on the spurious mass loose enough to keep passing
# once upstream fixes it. `nside = 9` is level 2, the first ISEA9R resolution a
# user is likely to reach for, so silence here would be the wrong trade.
@testset "self-regrid at nside 9 hits an upstream convex-clip false positive" begin
    nside = 9
    ncell = 10nside^2
    tree = treeify(Isea9rFaceGrid(nside; ordering=RowMajorOrder()))
    R = CR.Regridder(UNIT, tree, tree; threaded=false, normalize=false)

    # The grid itself is sound: every cell overlaps itself, the diagonal is the
    # cell areas, and the areas sum to the sphere.
    @test all(i -> R.intersections[i, i] > 0, 1:ncell)
    @test SparseArrays.diag(R.intersections) ≈ R.dst_areas
    @test isapprox(sum(R.dst_areas), 4π; rtol=1e-12)
    @test isapprox(sum(R.src_areas), 4π; rtol=1e-12)

    # The clipper is not: a handful of off-diagonal entries may appear, each
    # bounded by a cell area. The assertion is a ceiling, not the current value.
    offdiag = sum(abs, R.intersections) - sum(abs, SparseArrays.diag(R.intersections))
    @test SparseArrays.nnz(R.intersections) >= ncell
    @test offdiag <= 5 * (4π / ncell)
    record!("upstream convex-clip spurious self-overlap at nside 9 (sr)", offdiag)
    record!("...as a fraction of 4π", offdiag / 4π)
end

@testset "cross-resolution Morton nesting (nside 9 <- 3)" begin
    fine = treeify(Isea9rFaceGrid(9; ordering=MortonOrder()))
    coarse = treeify(Isea9rFaceGrid(3; ordering=MortonOrder()))
    R = CR.Regridder(UNIT, fine, coarse; threaded=false, normalize=false)
    M = R.intersections
    @test size(M) == (10 * 81, 10 * 9)

    # The base-9 Morton prefix property (`position ÷ 9` is the parent's
    # position) made observable through the clipper. It shows up as *dominance*,
    # not as an exact 9-nonzero pattern: unlike S2, ISEA9R cell edges are not
    # great circles, so a coarse cell's 4-gon is not the union of its nine
    # children's 4-gons and a few percent of its area leaks into geometric
    # neighbours. That is a property of the chord approximation, not of the
    # ordering — ISEA4R's nested 4 <- 2 leaks 5.4% by the same mechanism.
    worst_spurious = 0.0
    for j in 1:size(M, 2)
        span = M.colptr[j]:(M.colptr[j + 1] - 1)
        rows = M.rowval[span]
        vals = M.nzval[span]
        children = collect((9j - 8):(9j))
        # All nine children are present ...
        @test issubset(children, rows)
        # ... and they are the nine largest entries of the column.
        order = sortperm(vals; rev=true)
        @test sort(rows[order[1:9]]) == children
        spurious = sum(v for (r, v) in zip(rows, vals) if !(r in children); init=0.0)
        worst_spurious = max(worst_spurious, spurious / sum(vals))
    end
    @test worst_spurious <= 0.05
    record!("cross-resolution mass outside the nine Morton children (fraction)", worst_spurious)

    # Conservation is unaffected by that redistribution: it is a reshuffle
    # between neighbours, not a loss.
    @test isapprox(vec(sum(M, dims=2)), R.dst_areas; rtol=1e-12)
    @test isapprox(vec(sum(M, dims=1)), R.src_areas; rtol=1e-12)
    @test isapprox(sum(M), 4π; rtol=1e-10)
    record!("|sum(intersections) - 4π|, Morton 9<-3", abs(sum(M) - 4π))
end

@testset "non-nested cross-resolution (nside 5 <- 4)" begin
    R = CR.Regridder(UNIT, treeify(Isea9rFaceGrid(5)), treeify(Isea9rFaceGrid(4));
        threaded=false, normalize=false)
    @test size(R.intersections) == (250, 160)
    @test isapprox(sum(R.intersections), 4π; rtol=1e-10)
    record!("|sum(intersections) - 4π|, 5<-4", abs(sum(R.intersections) - 4π))
    dst = CR.regrid!(zeros(250), R, fill(3.5, 160))
    @test maximum(abs, dst .- 3.5) <= 1e-12
end

@testset "cross-system: ISEA9R against HEALPix (nside 9 <- 8)" begin
    i9_tree = treeify(Isea9rFaceGrid(9; ordering=MortonOrder()))
    hp_tree = treeify(HealpixFaceGrid(8; ordering=RingOrder()))
    R = CR.Regridder(UNIT, i9_tree, hp_tree; threaded=false, normalize=false)

    @test size(R.intersections) == (10 * 81, 12 * 64)
    # Two completely different tessellations of the same sphere: the shared
    # total is the real check that both sets of polygons are wound the same way
    # and cover everything exactly once.
    @test isapprox(sum(R.intersections), 4π; rtol=1e-10)
    @test isapprox(sum(R.src_areas), 4π; rtol=1e-10)
    @test isapprox(vec(sum(R.intersections, dims=2)), R.dst_areas; rtol=1e-12)
    record!("|sum(intersections) - 4π|, ISEA9R <- HEALPix", abs(sum(R.intersections) - 4π))

    # A constant field regrids to itself — the defining property of a
    # conservative (mean-preserving) operator, here across systems.
    dst = CR.regrid!(zeros(10 * 81), R, fill(3.5, 12 * 64))
    @test maximum(abs, dst .- 3.5) <= 1e-12
    record!("constant-field regrid deviation, ISEA9R <- HEALPix", maximum(abs, dst .- 3.5))
end

# --------------------------------------------------------------------------
# 6. Row-major vs Morton: one grid, two orderings, one permutation
# --------------------------------------------------------------------------

@testset "row-major/Morton Regridders differ by exactly the index permutation" begin
    nside = 9
    ncell = 10nside^2
    morton_tree = treeify(Isea9rFaceGrid(nside; ordering=MortonOrder()))
    rowmajor_tree = treeify(Isea9rFaceGrid(nside; ordering=RowMajorOrder()))
    M = CR.Regridder(UNIT, morton_tree, morton_tree; threaded=false, normalize=false)
    R = CR.Regridder(UNIT, rowmajor_tree, rowmajor_tree; threaded=false, normalize=false)

    # σ maps a Morton data position to the row-major data position of the same
    # cell: position i holds Morton id i - 1, whose lattice cell has row-major
    # id `xyd_to_rowmajor(...)`, hence data position `+ 1`.
    sigma = [Int(xyd_to_rowmajor(morton_to_xyd(i - 1, nside)..., nside)) + 1 for i in 1:ncell]
    @test sort(sigma) == collect(1:ncell)              # it is a permutation
    @test sigma != collect(1:ncell)                    # and not the identity

    # Both orderings hand the clipper *bit-identical* polygons — same chart
    # kernel, same corner order, only the data slot differs — and assembly with
    # `threaded = false` is deterministic. So this is exact equality, not
    # `isapprox`: any drift here would mean the ordering leaked into the
    # geometry, which is precisely what this layer is designed to prevent.
    # (It holds through the upstream clip anomaly above too: the same wrong
    # value lands in the permuted slot.)
    @test M.intersections == R.intersections[sigma, sigma]
    @test M.dst_areas == R.dst_areas[sigma]
    @test M.src_areas == R.src_areas[sigma]
    # And the geometry itself, cell by cell.
    for i in 1:ncell
        @test ring_points(Trees.getcell(morton_tree, i)) ==
              ring_points(Trees.getcell(rowmajor_tree, sigma[i]))
    end
end

# --------------------------------------------------------------------------
# 7. Non-power-of-three conservation and the equal-area claim
#
# The whole reason this layer is separate from the (future) id hierarchy:
# `nside = 4` and `nside = 5` have no aperture-9 id space at all, but they are
# perfectly good diamond grids and must regrid conservatively.
# --------------------------------------------------------------------------

@testset "each diamond carries exactly 4π/10" begin
    @test is_equal_area(ISEA9RDGGS())
    grid = Isea9rFaceGrid(1)
    R = CR.Regridder(UNIT, treeify(grid), treeify(grid); threaded=false, normalize=false)
    @test length(R.dst_areas) == 10
    # At `nside = 1` the cell boundary IS the diamond boundary, four icosahedron
    # edges, which are great-circle arcs — so the 4-gon is exact here and the
    # equal-area statement is directly visible. This is also the level-0 zone
    # area OGC 21-038r1 Annex B.2 states for ISEA9R: `4π / (10 * 9^0)`.
    for a in R.dst_areas
        @test isapprox(a, 4π / 10; rtol=1e-12)
    end
    record!("max |diamond area - 4π/10| / (4π/10)",
        maximum(abs.(R.dst_areas .- 4π / 10) ./ (4π / 10)))
end

@testset "non-power-of-three conservation (nside = $nside)" for nside in (4, 5)
    @test !ispow3(nside)
    ncell = 10nside^2
    grid = Isea9rFaceGrid(nside; ordering=RowMajorOrder())
    R = CR.Regridder(UNIT, treeify(grid), treeify(grid); threaded=false, normalize=false)

    @test size(R.intersections) == (ncell, ncell)
    @test isapprox(sum(R.intersections), 4π; rtol=1e-10)
    @test isapprox(sum(R.dst_areas), 4π; rtol=1e-10)
    @test isapprox(sum(R.src_areas), 4π; rtol=1e-10)
    record!("|sum(intersections) - 4π|, nside=$nside", abs(sum(R.intersections) - 4π))

    # Every cell overlaps itself: no cell may be missed by the traversal.
    @test all(i -> R.intersections[i, i] > 0, 1:ncell)

    # A constant field regrids to itself — the defining property of a
    # conservative (mean-preserving) operator.
    dst = CR.regrid!(zeros(ncell), R, fill(3.5, ncell))
    @test maximum(abs, dst .- 3.5) <= 1e-12
    record!("constant-field regrid deviation, nside=$nside", maximum(abs, dst .- 3.5))

    # The chart is exactly equal-area (`test/ISEA4R/test_diamonds.jl` pins that
    # on densified boundaries), but the 4-gon handed to the clipper is a chord
    # approximation of a curved cell, so per-cell areas deviate. The deviation
    # is edge bulge plus corner-cell skew at the twelve icosahedron vertices.
    # The SUM stays exact to ~1e-14 because within-diamond bulges cancel
    # pairwise on bit-identical shared edges.
    target = 4π / (10nside^2)
    worst = maximum(abs.(R.dst_areas .- target) ./ target)
    @test worst < 0.20
    record!("per-cell 4-gon area deviation, nside=$nside (relative)", worst)
end

# --------------------------------------------------------------------------
# 8. CCW discipline
#
# The convex-clip kernel clips a clockwise ring to EMPTY, so a reversed ring
# yields silent zero intersection areas rather than an error. Every polygon this
# layer emits must wind CCW as seen from outside the sphere.
# --------------------------------------------------------------------------

@testset "getcell polygons are CCW from outside (nside = $nside)" for nside in (3, 4, 9)
    ordering = RowMajorOrder()
    root = treeify(Isea9rFaceGrid(nside; ordering))
    worst = Inf
    for j in 1:(10nside^2)
        worst = min(worst, ccw_measure(open_ring(Trees.getcell(root, j))))
    end
    @test worst > 0
    # Same for the per-diamond grid's `(i, j)` accessor, which is the interface
    # method the cursor machinery calls.
    for diamond in 0:9
        g = DiamondChartGrid(GO.Spherical(), Isea9rFaceSpace(nside), diamond, ordering)
        for i in 1:nside, j in 1:nside
            @test ccw_measure(open_ring(Trees.getcell(g, i, j))) > 0
        end
    end
    record!("min CCW measure (nside=$nside)", worst)
end

# --------------------------------------------------------------------------
# 9. Seams, at grid level
#
# `test/ISEA4R/test_diamonds.jl` pins the seam tolerances on the chart; what
# matters here is that they survive into the polygons the clipper sees at the
# aperture-9 resolutions.
# --------------------------------------------------------------------------

@testset "seam and diamond borders in the emitted polygons (nside = $nside)" for nside in (3, 9)
    root = treeify(Isea9rFaceGrid(nside; ordering=RowMajorOrder()))
    corners = Dict{NTuple{3,Int},NTuple{4,Any}}()
    for diamond in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
        corners[(ix, iy, diamond)] =
            Tuple(open_ring(Trees.getcell(root,
                data_index(RowMajorOrder(), Isea9rFaceSpace(nside), ix, iy, diamond))))
    end
    # Cells straddling the seam (`ix == iy`) share bit-identical corners with
    # their neighbours in both directions — the seam-ownership rule surviving
    # into the emitted geometry.
    for diamond in 0:9, k in 0:(nside - 2)
        here = corners[(k, k, diamond)]
        east = corners[(k + 1, k, diamond)]
        north = corners[(k, k + 1, diamond)]
        @test here[1] === north[4] && here[2] === north[3]
        @test here[1] === east[2] && here[4] === east[3]
    end

    # Across diamonds the polygons cannot be bit-identical — the two sides come
    # from different Snyder faces — but they must agree to the Newton floor.
    worst = 0.0
    border = Dict{Int,Vector{Any}}()
    for diamond in 0:9
        g = DiamondChartGrid(GO.Spherical(), Isea9rFaceSpace(nside), diamond, RowMajorOrder())
        ps = Any[]
        for i in 1:(nside + 1), j in 1:(nside + 1)
            (i == 1 || i == nside + 1 || j == 1 || j == nside + 1) || continue
            push!(ps, Trees.getvertex(g, i, j))
        end
        border[diamond] = ps
    end
    for d in 0:9, p in border[d]
        best = Inf
        for e in 0:9
            e == d && continue
            for q in border[e]
                best = min(best, US.spherical_distance(p, q))
            end
        end
        worst = max(worst, best)
    end
    @test worst < 1e-13
    record!("cross-diamond border mismatch in getvertex (rad)", worst)
end

@printf("[ISEA9R face grid] measured:\n")
for key in sort!(collect(keys(MEASURED)))
    @printf("[ISEA9R face grid]   %-58s %+.3e\n", key, MEASURED[key])
end

end # module Isea9rFaceGridTestSuite
