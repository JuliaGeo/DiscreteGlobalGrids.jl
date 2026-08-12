module Isea9rKernelTestSuite

# Tests for `src/ISEA9R/Isea9rKernel.jl`: the operations-kernel wiring of
# `ISEA9RDGGS`. The sibling of `test/ISEA4R/test_isea4r_kernel.jl`, section for
# section, at radix 9.
#
# The ten-diamond *layout* is OGC 21-038r1 Annex B.2's ("The ten root rhombuses
# are formed by combining two icosahedron triangles at their base") and DGGAL's
# `countZones(level) = 10 * 9^level`; the *numbering* of those ten, the
# in-diamond axes and the in-diamond index are this package's own conventions
# with no external oracle behind them. So — exactly as in
# `test/ISEA4R/test_isea4r_kernel.jl` — nothing here checks against a reference
# implementation, and nothing here may be read as a DGGRID / DGGAL /
# SphericalSpatialTrees compatibility claim. See
# `docs/design/isea9r_layout.md`.
#
# What is checked is the two claims the wiring makes:
#
# 1. *The kernel IS the chart, bitwise.* `cell_boundary(ISEA9RDGGS(), level,
#    id)` and the dense diamond grid at `nside = 3^level` under `MortonOrder`
#    must emit LITERALLY the same `Float64`s for the same `isea9r_ordinal`,
#    because they share one evaluation of `cell_corners ∘ morton_to_xyd`.
#
# 2. *The milestone boundary is where the docs say it is.* Geometry answers; the
#    hierarchy / ordinal / pruning group still throws `NotPortedError` — the
#    same line the ISEA4R sibling holds, deliberately, so the two systems stay
#    one decision rather than two. The last testsets pin that so the trait
#    cannot flip silently.
#
# Naming: `ISEA9R.cell_polygon` is a function in the `ISEA9R` namespace, NOT a
# method of the kernel generic of the same name. The ISEA9R ones are always
# written qualified here and the DGGS ones through `DGG.` (or bare, for the
# exported `cell_polygon`). `cell_corners` / `cell_center` reach this module
# from `ISEA4R` by import — that is the delegation `test_delegation.jl` pins.

using Test
using Printf
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids
using DiscreteGlobalGrids.ISEA9R
import DiscreteGlobalGrids.ISEA9R as ISEA9R
import DiscreteGlobalGrids.ISEA4R as ISEA4R
using DiscreteGlobalGrids.ISEA9R: DiamondChartGrid, cell_corners, cell_center,
    morton_to_xyd, xyd_to_morton, xyd_to_point
import ConservativeRegridding: Trees
import GeometryOps as GO
import GeoInterface as GI

const S = ISEA9RDGGS()
const US = GO.UnitSpherical

# Bit-level equality, not `==`: `==` would accept `0.0` for `-0.0`, and a signed
# zero that differs between two paths is already a fork.
identical(a, b) = all(k -> a[k] === b[k], 1:3)

ring_points(poly) = collect(GI.getpoint(GI.getexterior(poly)))

cross3(a, b) = (a[2] * b[3] - a[3] * b[2],
                a[3] * b[1] - a[1] * b[3],
                a[1] * b[2] - a[2] * b[1])

# Signed area of the 4-gon seen from *outside* the sphere; positive ⇔ CCW.
function ccw_measure(corners)
    acc = (0.0, 0.0, 0.0)
    n = length(corners)
    for i in 1:n
        acc = acc .+ cross3(corners[i], corners[i % n + 1])
    end
    outward = reduce((a, b) -> a .+ Tuple(b), corners; init=(0.0, 0.0, 0.0))
    return sum(acc .* outward)
end

# Levels 0:2 are swept exhaustively (810 ordinals at level 2); level 3 (7290) is
# sampled deterministically — first/middle/last of every diamond plus a fixed
# stride — so a failure is reproducible verbatim without an RNG.
sweep_ids(level) = collect(Int64, 0:(10 * 9^level - 1))
function sample_ids(level)
    per_diamond = 9^level
    offsets = unique(filter(o -> 0 <= o < per_diamond,
        [0, 1, per_diamond ÷ 3, per_diamond ÷ 2, per_diamond - 1]))
    ids = Int64[diamond * per_diamond + o for diamond in 0:9 for o in offsets]
    append!(ids, Int64.(0:(10 * per_diamond ÷ 37):(10 * per_diamond - 1)))
    return sort!(unique!(ids))
end

const MEASURED = Dict{String,Float64}()
record!(key, value) = (MEASURED[key] = max(get(MEASURED, key, -Inf), value))

# --------------------------------------------------------------------------
# 1. The id model: what is wired and what deliberately is not
# --------------------------------------------------------------------------

@testset "ISEA9R kernel id model" begin
    @test DGG.cell_id_type(S) === Int64

    # The hierarchy/ordinal/pruning group is OFF — a scope line, not a
    # limitation of the id space: `isea9r_ordinal` IS the canonical index and
    # `9p:9p+8` really are cell `p`'s children (see `Isea9rKernel.jl`'s header).
    @test !DGG.has_ordinal_ids(S)
    @test !DGG.has_descendant_ranges(S)
    @test !DGG.has_exact_subtree_cap(S)
    @test canonical_index_name(S) === :isea9r_ordinal

    # The registry facts the ordinal arithmetic runs on. `root_count` and
    # `supports_prefix_ranges` are the two this milestone flipped: ten root
    # rhombuses per OGC 21-038r1 Annex B.2 and DGGAL's `countZones`, and exact
    # radix-9 prefix arithmetic over the PACKAGE ordinal (not over DGGAL's zone
    # id, whose descendants are rectangular blocks).
    @test root_count(S) == 10
    @test radix(S) == 9
    @test aperture(S) == 9
    @test max_level(S) === nothing
    @test supports_prefix_ranges(S)
    @test is_equal_area(S)
    @test cell_shape(S) === :rhomb
    @test base_solid(S) === :icosahedron

    # `cell_cap_inflation` is still the shared default, and it is now dead code
    # for ISEA9R: `cell_cap` is overridden with the exact four-corner cap below.
    @test DGG.cell_cap_inflation(S) == 1.2
end

# --------------------------------------------------------------------------
# 2. Geometry is the chart kernel, bitwise — and so is the diamond grid's
# --------------------------------------------------------------------------

@testset "ISEA9R kernel geometry is the chart kernel, bitwise (level $level)" for level in 0:2
    nside = 3^level
    # `MortonOrder` puts ordinal `p` at data position `p + 1`, so no permutation
    # sits between the two sides.
    root = treeify(Isea9rFaceGrid(nside; ordering=MortonOrder()))
    worst_ccw = Inf
    for id in sweep_ids(level)
        ix, iy, diamond = morton_to_xyd(id, nside)
        chart = cell_corners(ix, iy, diamond, nside)

        boundary = DGG.cell_boundary(S, level, id)
        @test length(boundary) == 4
        @test all(i -> identical(boundary[i], chart[i]), 1:4)

        closed = DGG.cell_boundary(S, level, id; closed=true)
        @test length(closed) == 5
        @test all(i -> identical(closed[i], chart[i]), 1:4)
        @test identical(closed[5], chart[1])
        @test closed[end] == closed[1]

        # CCW as seen from outside — the hard contract every polygon emitter in
        # this package carries (a CW ring clips to EMPTY, silently). Worth
        # sweeping per level here rather than trusting the chart tests: the
        # seam-straddling cells are the ones at risk, and the Morton sweep hits
        # every one of them.
        worst_ccw = min(worst_ccw, ccw_measure(boundary))

        center = DGG.cell_center(S, level, id)
        @test identical(center, cell_center(ix, iy, diamond, nside))
        @test sum(abs2, center) ≈ 1.0 atol = 1e-15

        # The polygon rings — what `ConservativeRegridding` actually clips — on
        # the kernel path and on the diamond-grid path, for the same ordinal.
        kernel_ring = ring_points(DGG.cell_polygon_unitsphere(S, level, id))
        face_ring = ring_points(Trees.getcell(root, Int(id) + 1))
        @test length(kernel_ring) == length(face_ring) == 5
        @test all(i -> identical(kernel_ring[i], face_ring[i]), eachindex(kernel_ring))
    end
    @test worst_ccw > 0
    # Recorded negated, because `record!` keeps the max and this is a minimum:
    # a value below zero means every ring on every level swept CCW.
    record!("least CCW measure, negated (< 0 ⇒ all CCW)", -worst_ccw)
end

@testset "ISEA9R kernel geometry is the chart kernel, bitwise (sampled, level 3)" begin
    level = 3
    nside = 3^level
    root = treeify(Isea9rFaceGrid(nside; ordering=MortonOrder()))
    for id in sample_ids(level)
        ix, iy, diamond = morton_to_xyd(id, nside)
        chart = cell_corners(ix, iy, diamond, nside)
        boundary = DGG.cell_boundary(S, level, id)
        @test all(i -> identical(boundary[i], chart[i]), 1:4)
        @test identical(DGG.cell_center(S, level, id),
                        cell_center(ix, iy, diamond, nside))
        kernel_ring = ring_points(DGG.cell_polygon_unitsphere(S, level, id))
        face_ring = ring_points(Trees.getcell(root, Int(id) + 1))
        @test all(i -> identical(kernel_ring[i], face_ring[i]), eachindex(kernel_ring))
        # The round trip through the ordinal decomposition the registry records:
        # `id == diamond * 9^level + morton9_position`, with the position in range.
        @test xyd_to_morton(ix, iy, diamond, nside) == id
        @test 0 <= id - diamond * 9^level < 9^level
    end
end

# --------------------------------------------------------------------------
# 3. `cell_polygon` at the interface level
# --------------------------------------------------------------------------

@testset "cell_polygon(ISEA9RDGGS(), level, id) answers (level $level)" for level in 0:2
    nside = 3^level
    for id in sweep_ids(level)
        polygon = cell_polygon(S, level, id)          # the interface generic
        @test polygon isa GI.Polygon
        # Same object shape as the kernel generic it delegates to, point for
        # point: the interface method is a pure bridge, not a second geometry.
        unit = DGG.cell_polygon_unitsphere(S, level, id)
        @test all(i -> identical(ring_points(polygon)[i], ring_points(unit)[i]), 1:5)
        # ...and the same polygon the chart-side `ISEA9R.cell_polygon` builds.
        ix, iy, diamond = morton_to_xyd(id, nside)
        @test all(i -> identical(ring_points(polygon)[i],
                                 ring_points(ISEA9R.cell_polygon(ix, iy, diamond, nside))[i]), 1:5)
        # `cell_extent` is derived from `cell_polygon`, so it answers too — and
        # it really bounds the ring it was derived from.
        extent = cell_extent(S, level, id)
        @test all(p -> extent.X[1] <= p[1] <= extent.X[2] &&
                       extent.Y[1] <= p[2] <= extent.Y[2] &&
                       extent.Z[1] <= p[3] <= extent.Z[2],
                  DGG.cell_boundary(S, level, id))
    end
end

# --------------------------------------------------------------------------
# 4. Caps
#
# `cell_cap` is overridden with the exact four-corner cap. Snyder cell edges are
# NOT great circles, so unlike S2 this is a measured property rather than a
# proved one: the pre-registered measurement at `cap_policy(::Isea9rFaceSystem)`
# — this system's own, at `nside ∈ (3, 9, 27)` — found a worst pre-slack
# overhang of exactly 0.0 over a dense sampling of every block, leaf blocks
# included. This testset re-runs the leaf case against these caps specifically.
# --------------------------------------------------------------------------

@testset "ISEA9R cell caps contain the cell (level $level)" for level in 0:2
    nside = 3^level
    worst = -Inf
    for id in sweep_ids(level)
        cap = DGG.cell_cap(S, level, id)
        ix, iy, diamond = morton_to_xyd(id, nside)
        for point in DGG.cell_boundary(S, level, id)
            worst = max(worst, US.spherical_distance(cap.point, point) - cap.radius)
        end
        @test US.spherical_distance(cap.point, DGG.cell_center(S, level, id)) <= cap.radius
        # Dense sampling of the cell's chart rectangle — the part that matters
        # here, since the curved Snyder edges bulge off the corner-to-corner
        # chords and a corners-only check would not see it.
        for a in 0:8, b in 0:8
            p = xyd_to_point((ix + a / 8) / nside, (iy + b / 8) / nside, diamond)
            worst = max(worst, US.spherical_distance(cap.point, p) - cap.radius)
        end
        # `subtree_cap` at `leaf_level == level` degenerates to the cell cap;
        # deeper it needs the unwired hierarchy (asserted in section 6).
        @test DGG.subtree_cap(S, level, id, level) === cap
    end
    @test worst <= 0
    record!("cell-cap overhang, corners + dense chart samples (rad)", worst)
end

@testset "ISEA9R cell caps are the diamond-grid four-corner cap, bitwise" begin
    for level in 0:2
        nside = 3^level
        for id in sweep_ids(level)
            ix, iy, diamond = morton_to_xyd(id, nside)
            # The four block corners exactly as the face-grid layer harvests
            # them (`Trees.getvertex` on the per-diamond chart grid), in
            # `circle_from_four_corners`' documented `(BL, TL, BR, TR)` slots —
            # which is the permutation `Isea9rKernel.jl` applies to the chart
            # ring.
            g = DiamondChartGrid(GO.Spherical(), Isea9rFaceSpace(nside), diamond, MortonOrder())
            bl = Trees.getvertex(g, ix + 1, iy + 1)
            tl = Trees.getvertex(g, ix + 1, iy + 2)
            br = Trees.getvertex(g, ix + 2, iy + 1)
            tr = Trees.getvertex(g, ix + 2, iy + 2)
            @test DGG.cell_cap(S, level, id) === Trees.circle_from_four_corners((bl, tl, br, tr), ())
        end
    end
end

@testset "ISEA9R cell caps beat the generic inflated cap" begin
    # What the generic `cell_cap` would have returned: max center-to-vertex
    # distance × `cell_cap_inflation`. The exact cap must be strictly tighter
    # everywhere — that is the whole reason for the override.
    worst_ratio = -Inf
    for level in 0:2
        for id in sweep_ids(level)
            cap = DGG.cell_cap(S, level, id)
            center = DGG.cell_center(S, level, id)
            radius = maximum(US.spherical_distance(center, p)
                             for p in DGG.cell_boundary(S, level, id))
            generic = nextfloat(min(Float64(pi), radius * DGG.cell_cap_inflation(S) + 1e-9))
            @test cap.radius < generic
            worst_ratio = max(worst_ratio, cap.radius / generic)
        end
    end
    record!("exact cap / generic inflated cap radius, worst", worst_ratio)
end

# --------------------------------------------------------------------------
# 5. Id and level validation — inherited from the codecs, not reinvented
# --------------------------------------------------------------------------

@testset "ISEA9R kernel id/level validation" begin
    # Out-of-range ordinals are `morton_to_xyd`'s `ArgumentError`; the kernel
    # adds no validation layer of its own. (No `OrdinalRangeError` here: that
    # type belongs to `ordinal_to_cell`, the 1-based dense-position generic,
    # which for ISEA9R still throws `NotPortedError` — see section 6. The
    # sibling `Isea4rKernel.jl` draws the same line.)
    for level in 0:2
        count = 10 * 9^level
        @test_throws ArgumentError DGG.cell_boundary(S, level, count)
        @test_throws ArgumentError DGG.cell_boundary(S, level, -1)
        @test_throws ArgumentError DGG.cell_center(S, level, count)
        @test_throws ArgumentError DGG.cell_cap(S, level, count)
        @test_throws ArgumentError cell_polygon(S, level, count)
        # The last valid id does answer, so the bound is exactly `count - 1`.
        @test DGG.cell_boundary(S, level, count - 1) isa Vector
    end
    # A negative level is `3^Int(level)`'s `DomainError`, exactly as `2^level`
    # is for HEALPix, S2 and ISEA4R — matched deliberately rather than
    # invented, so the kernel's failure modes do not differ per system.
    @test_throws DomainError DGG.cell_boundary(S, -1, 0)
    @test_throws DomainError DGG.cell_boundary(ISEA4RDGGS(), -1, 0)
end

# --------------------------------------------------------------------------
# 6. THE BOUNDARY: geometry answers, the hierarchy does not
#
# `has_ordinal_ids(::ISEA9RDGGS)` is false, so every hierarchy/ordinal/pruning
# generic falls through to its `NotPortedError`. Flipping that trait is a
# one-line change that would make all of these pass — which is exactly why it
# must come with its own kernel-test battery rather than slip in unnoticed.
# --------------------------------------------------------------------------

@testset "the ISEA9R id hierarchy is still not wired" begin
    @test_throws NotPortedError DGG.root_ids(S)
    @test_throws NotPortedError DGG.cell_children(S, 0, 0)
    @test_throws NotPortedError DGG.cell_parent(S, 1, 9, 0)
    @test_throws NotPortedError DGG.cell_descendants(S, 0, 0, 1)
    @test_throws NotPortedError DGG.cell_to_ordinal(S, 0, 0)
    @test_throws NotPortedError DGG.ordinal_to_cell(S, 0, 1)
    @test_throws NotPortedError DGG.descendant_range(S, 0, 0, 1)
    @test_throws NotPortedError DGG.num_cells(S, 0)
    @test_throws NotPortedError DGG.subtree_leaf_count(S, 0, 0, 1)
    @test_throws NotPortedError DGG.subtree_cap(S, 0, 0, 1)

    # The *interface*-level prefix arithmetic is a different function, and for
    # ISEA9R it started answering with this milestone: it derives from
    # `supports_prefix_ranges` / `radix` alone and claims nothing about
    # geometry. Worth pinning next to the kernel's silence so the two are not
    # confused for one another — and worth pinning at all, because the trait it
    # rests on was `false` before.
    @test leaf_count(S, 2) == 10 * 9^2
    @test leaf_interval(S, 0, 9, 2) == (9 * 81):(10 * 81 - 1)
    @test child_ids(S, 0, 3) == collect(27:35)
    @test leaf_interval(S, 1, 5, 3) == (5 * 81):(6 * 81 - 1)

    # `cells_cap` needs only boundaries, so it answers too.
    cap = DGG.cells_cap(S, 1, Int64[0, 1, 2, 3])
    @test cap isa US.SphericalCap
    @test all(id -> all(p -> US.spherical_distance(cap.point, p) <= cap.radius,
                        DGG.cell_boundary(S, 1, id)), 0:3)
end

# --------------------------------------------------------------------------
# 7. The sibling still answers over ITS ordinal
#
# The two systems share a chart but not an id space: the same `(level, id)` pair
# means different cells, because `4^level` and `9^level` differ and so do the
# two Morton radices. Pinned here so a future hoist cannot quietly merge them.
# --------------------------------------------------------------------------

@testset "ISEA4R and ISEA9R ordinals are different id spaces" begin
    # Level 2: ISEA4R has 160 cells, ISEA9R has 810.
    @test DGG.cell_boundary(ISEA4RDGGS(), 2, 159) isa Vector
    @test_throws ArgumentError DGG.cell_boundary(ISEA4RDGGS(), 2, 160)
    @test DGG.cell_boundary(S, 2, 809) isa Vector
    @test_throws ArgumentError DGG.cell_boundary(S, 2, 810)
    # Level 0 is the one level where the two id spaces coincide — ten diamonds,
    # one cell each — and there the geometry is identical, bit for bit.
    for id in 0:9
        @test all(i -> identical(DGG.cell_boundary(S, 0, id)[i],
                                 DGG.cell_boundary(ISEA4RDGGS(), 0, id)[i]), 1:4)
    end
end

@printf("[ISEA9R kernel] measured:\n")
for key in sort!(collect(keys(MEASURED)))
    @printf("[ISEA9R kernel]   %-52s %+.3e\n", key, MEASURED[key])
end

end # module Isea9rKernelTestSuite
