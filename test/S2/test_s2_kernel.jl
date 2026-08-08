module S2KernelTestSuite

# Tests for `src/S2/S2Kernel.jl`: the operations-kernel wiring of `S2DGGS`.
#
# There is no S2 oracle in this repository (no s2geometry fixtures, no
# s2geometry dependency), so — exactly as in `test/S2/test_chart.jl` — nothing
# here checks against a reference implementation. What it checks instead is the
# two claims the wiring actually makes:
#
# 1. *The kernel IS the chart, bitwise.* `cell_boundary(S2DGGS(), level, id)`
#    and the dense face grid at `nside = 2^level` under `HilbertOrder` must emit
#    LITERALLY the same `Float64`s for the same scaffold ordinal, because they
#    share one evaluation of `cell_corners ∘ hilbert_to_xyf`. This is the honest
#    analogue of `test/HEALPix/test_healpix_kernel.jl`'s "geometry is the chart
#    kernel, bitwise" testset — HEALPix compares two *grid* paths through a
#    `DGGSGrid`, which S2 has no access to (`has_ordinal_ids` is false), so the
#    comparison runs polygon ring against polygon ring.
#
# 2. *The milestone boundary is where the docs say it is.* Dense geometry
#    enumeration answers; the hierarchy / pruning group still throws `NotPortedError`,
#    because S2's canonical index is the native 64-bit `s2_cellid` and this
#    package indexes cells by the registry's scaffold ordinal
#    `face * 4^level + hilbert_position` instead. The last testset pins that,
#    so a future `has_ordinal_ids(::S2DGGS) = true` cannot land silently.
#
# Naming: `S2.cell_corners` / `S2.cell_center` / `S2.cell_polygon` are functions
# in the `S2` namespace, NOT methods of the kernel generics of the same name.
# The S2 ones are always written qualified here and the DGGS ones through `DGG.`
# (or bare, for the exported `cell_polygon`).

using Test
using Printf
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids
using DiscreteGlobalGrids.S2
using DiscreteGlobalGrids.S2: FaceChartGrid, cell_corners, hilbert_to_xyf, xyf_to_hilbert,
    stf_to_point
import ConservativeRegridding: Trees
import GeometryOps as GO
import GeoInterface as GI

const S = S2DGGS()
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

# Levels 0:3 are swept exhaustively (384 ordinals at level 3); level 5 is
# sampled deterministically — first/middle/last of every face plus a fixed
# stride — so a failure is reproducible verbatim without an RNG.
sweep_ids(level) = collect(Int64, 0:(6 * 4^level - 1))
function sample_ids(level)
    per_face = 4^level
    offsets = unique(filter(o -> 0 <= o < per_face, [0, 1, per_face ÷ 3, per_face ÷ 2, per_face - 1]))
    ids = Int64[face * per_face + o for face in 0:5 for o in offsets]
    append!(ids, Int64.(0:(6 * per_face ÷ 37):(6 * per_face - 1)))
    return sort!(unique!(ids))
end

const MEASURED = Dict{String,Float64}()
record!(key, value) = (MEASURED[key] = max(get(MEASURED, key, -Inf), value))

# --------------------------------------------------------------------------
# 1. The id model: what is wired and what deliberately is not
# --------------------------------------------------------------------------

@testset "S2 kernel id model" begin
    # The kernel's default id type. Scaffold ordinals are plain dense `Int64`s;
    # a native `s2_cellid` port would make this `UInt64`, which is one of the
    # reasons the two id spaces are not conflated.
    @test DGG.cell_id_type(S) === Int64

    # The whole hierarchy/ordinal/pruning group is OFF, on purpose: the
    # canonical index is `s2_cellid`, not the scaffold ordinal these geometry
    # methods take, so deriving the hierarchy from `root_count`/`radix` would
    # answer in the wrong coordinate system.
    @test !DGG.has_ordinal_ids(S)
    @test !DGG.has_descendant_ranges(S)
    @test !DGG.has_exact_subtree_cap(S)
    @test canonical_index_name(S) === :s2_cellid

    # Unchanged registry facts the ordinal arithmetic *would* use if the native
    # id decision ever went the other way (see `S2Kernel.jl`'s header).
    @test root_count(S) == 6
    @test radix(S) == 4
    @test max_level(S) == 30
    @test supports_prefix_ranges(S)

    # `cell_cap_inflation` is still the shared default, and it is now dead code
    # for S2: `cell_cap` is overridden with the exact four-corner cap below.
    @test DGG.cell_cap_inflation(S) == 1.2
end

# --------------------------------------------------------------------------
# 2. Geometry is the chart kernel, bitwise — and so is the face grid's
# --------------------------------------------------------------------------

@testset "S2 kernel geometry is the chart kernel, bitwise (level $level)" for level in 0:3
    nside = 2^level
    ids = sweep_ids(level)
    # `HilbertOrder` puts scaffold ordinal `p` at data position `p + 1`, so no
    # permutation sits between the two sides.
    root = treeify(S2FaceGrid(nside; ordering=HilbertOrder()))
    worst_ccw = Inf
    for id in ids
        ix, iy, face = hilbert_to_xyf(id, nside)
        chart = cell_corners(ix, iy, face, nside)

        boundary = DGG.cell_boundary(S, level, id)
        @test length(boundary) == 4
        @test all(i -> identical(boundary[i], chart[i]), 1:4)

        closed = DGG.cell_boundary(S, level, id; closed=true)
        @test length(closed) == 5
        @test all(i -> identical(closed[i], chart[i]), 1:4)
        @test identical(closed[5], chart[1])
        @test closed[end] == closed[1]

        # CCW as seen from outside — the hard contract every polygon emitter in
        # this package carries (a CW ring clips to EMPTY, silently).
        worst_ccw = min(worst_ccw, ccw_measure(boundary))

        # Centers: S2's own definition, and a unit vector.
        center = DGG.cell_center(S, level, id)
        @test identical(center, S2.cell_center(ix, iy, face, nside))
        @test sum(abs2, center) ≈ 1.0 atol = 1e-15

        # The polygon rings — what `ConservativeRegridding` actually clips — on
        # the kernel path and on the face-grid path, for the same ordinal.
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

@testset "S2 kernel geometry is the chart kernel, bitwise (sampled, level 5)" begin
    level = 5
    nside = 2^level
    root = treeify(S2FaceGrid(nside; ordering=HilbertOrder()))
    for id in sample_ids(level)
        ix, iy, face = hilbert_to_xyf(id, nside)
        chart = cell_corners(ix, iy, face, nside)
        boundary = DGG.cell_boundary(S, level, id)
        @test all(i -> identical(boundary[i], chart[i]), 1:4)
        @test identical(DGG.cell_center(S, level, id), S2.cell_center(ix, iy, face, nside))
        kernel_ring = ring_points(DGG.cell_polygon_unitsphere(S, level, id))
        face_ring = ring_points(Trees.getcell(root, Int(id) + 1))
        @test all(i -> identical(kernel_ring[i], face_ring[i]), eachindex(kernel_ring))
        # The round trip through the ordinal decomposition the registry records:
        # `id == face * 4^level + hilbert_position`, with the position in range.
        @test xyf_to_hilbert(ix, iy, face, nside) == id
        @test 0 <= id - face * 4^level < 4^level
    end
end

# --------------------------------------------------------------------------
# 3. `cell_polygon` at the interface level — the milestone-5 connection
# --------------------------------------------------------------------------

@testset "cell_polygon(S2DGGS(), level, id) answers (level $level)" for level in 0:2
    nside = 2^level
    for id in sweep_ids(level)
        polygon = cell_polygon(S, level, id)          # the interface generic
        @test polygon isa GI.Polygon
        # Same object shape as the kernel generic it delegates to, point for
        # point: the interface method is a pure bridge, not a second geometry.
        unit = DGG.cell_polygon_unitsphere(S, level, id)
        @test all(i -> identical(ring_points(polygon)[i], ring_points(unit)[i]), 1:5)
        # ...and the same polygon the chart-side `S2.cell_polygon` builds.
        ix, iy, face = hilbert_to_xyf(id, nside)
        @test all(i -> identical(ring_points(polygon)[i],
                                 ring_points(S2.cell_polygon(ix, iy, face, nside))[i]), 1:5)
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
# `cell_cap` is overridden with the exact four-corner cap. Three things to pin:
# it contains the cell (corners AND a dense sampling of the cell's chart
# rectangle, since a cap that bounds only the corners of a curved patch would
# be unsound in general — for S2 it is sound because the edges are geodesics,
# which is the proof recorded at `cap_policy(::S2FaceSystem)`), it is the same
# cap the face-grid layer's `FourCornerCap` computes from the same four lattice
# vertices, and it is tighter than the generic 1.2-inflated formula it replaces.
# --------------------------------------------------------------------------

@testset "S2 cell caps contain the cell (level $level)" for level in 0:3
    nside = 2^level
    worst = -Inf
    for id in sweep_ids(level)
        cap = DGG.cell_cap(S, level, id)
        ix, iy, face = hilbert_to_xyf(id, nside)
        for point in DGG.cell_boundary(S, level, id)
            worst = max(worst, US.spherical_distance(cap.point, point) - cap.radius)
        end
        @test US.spherical_distance(cap.point, DGG.cell_center(S, level, id)) <= cap.radius
        # Dense interior/boundary sampling of the cell's chart rectangle.
        for a in 0:8, b in 0:8
            p = stf_to_point((ix + a / 8) / nside, (iy + b / 8) / nside, face)
            worst = max(worst, US.spherical_distance(cap.point, p) - cap.radius)
        end
        # `subtree_cap` at `leaf_level == level` degenerates to the cell cap;
        # deeper it needs the unported hierarchy (asserted in section 6).
        @test DGG.subtree_cap(S, level, id, level) === cap
    end
    @test worst <= 0
    record!("cell-cap overhang, corners + dense chart samples (rad)", worst)
end

@testset "S2 cell caps are the face-grid four-corner cap, bitwise" begin
    for level in 0:2
        nside = 2^level
        for id in sweep_ids(level)
            ix, iy, face = hilbert_to_xyf(id, nside)
            # The four block corners exactly as the face-grid layer harvests
            # them (`Trees.getvertex` on the per-face chart grid), in
            # `circle_from_four_corners`' documented `(BL, TL, BR, TR)` slots —
            # which is the permutation `S2Kernel.jl` applies to the chart ring.
            g = FaceChartGrid(GO.Spherical(), S2FaceSpace(nside), face, HilbertOrder())
            bl = Trees.getvertex(g, ix + 1, iy + 1)
            tl = Trees.getvertex(g, ix + 1, iy + 2)
            br = Trees.getvertex(g, ix + 2, iy + 1)
            tr = Trees.getvertex(g, ix + 2, iy + 2)
            @test DGG.cell_cap(S, level, id) === Trees.circle_from_four_corners((bl, tl, br, tr), ())
        end
    end
end

@testset "S2 cell caps beat the generic inflated cap" begin
    # What the generic `cell_cap` would have returned: max center-to-vertex
    # distance × `cell_cap_inflation`. The exact cap must be strictly tighter
    # everywhere — that is the whole reason for the override.
    worst_ratio = -Inf
    for level in 0:3
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

@testset "S2 kernel id/level validation" begin
    # Out-of-range scaffold ordinals are `hilbert_to_xyf`'s `ArgumentError`; the
    # kernel adds no validation layer of its own.
    for level in 0:3
        count = 6 * 4^level
        @test_throws ArgumentError DGG.cell_boundary(S, level, count)
        @test_throws ArgumentError DGG.cell_boundary(S, level, -1)
        @test_throws ArgumentError DGG.cell_center(S, level, count)
        @test_throws ArgumentError DGG.cell_cap(S, level, count)
        @test_throws ArgumentError cell_polygon(S, level, count)
        # The last valid id does answer, so the bound is exactly `count - 1`.
        @test DGG.cell_boundary(S, level, count - 1) isa Vector
    end
    # A negative level is `2^Int(level)`'s `DomainError`, exactly as for
    # HEALPix — matched deliberately rather than invented, so the kernel's
    # failure modes do not differ per system.
    @test_throws DomainError DGG.cell_boundary(S, -1, 0)
    @test_throws DomainError DGG.cell_boundary(HEALPixDGGS(), -1, 0)
end

# --------------------------------------------------------------------------
# 6. THE BOUNDARY: dense geometry enumeration answers, the hierarchy does not
#
# `has_ordinal_ids(::S2DGGS)` is false, so hierarchy/pruning generics fall
# through to `NotPortedError`. `num_cells` and `ordinal_to_cell` are narrower:
# they enumerate the scaffold ids already accepted by geometry. This testset
# keeps the
# milestone honest: flipping that trait without wiring (and testing) the native
# `s2_cellid` semantics must FAIL here, not pass quietly.
# --------------------------------------------------------------------------

@testset "the S2 id hierarchy is still not ported" begin
    for level in 0:5
        count = 6 * 4^level
        @test DGG.num_cells(S, level) == count
        @test DGG.ordinal_to_cell(S, level, 1) == 0
        @test DGG.ordinal_to_cell(S, level, count) == count - 1
        @test_throws OrdinalRangeError DGG.ordinal_to_cell(S, level, 0)
        @test_throws OrdinalRangeError DGG.ordinal_to_cell(S, level, count + 1)
    end

    @test_throws NotPortedError DGG.root_ids(S)
    @test_throws NotPortedError DGG.cell_children(S, 0, 0)
    @test_throws NotPortedError DGG.cell_parent(S, 1, 4, 0)
    @test_throws NotPortedError DGG.cell_descendants(S, 0, 0, 1)
    @test_throws NotPortedError DGG.cell_to_ordinal(S, 0, 0)
    @test_throws NotPortedError DGG.descendant_range(S, 0, 0, 1)
    @test_throws NotPortedError DGG.subtree_leaf_count(S, 0, 0, 1)
    # `subtree_cap` above level equality needs `cell_descendants`, so it throws
    # too — the honest state of a system with geometry but no hierarchy.
    @test_throws NotPortedError DGG.subtree_cap(S, 0, 0, 1)

    # `cells_cap`, by contrast, needs only boundaries, so it now answers.
    cap = DGG.cells_cap(S, 1, Int64[0, 1, 2, 3])
    @test cap isa US.SphericalCap
    @test all(id -> all(p -> US.spherical_distance(cap.point, p) <= cap.radius,
                        DGG.cell_boundary(S, 1, id)), 0:3)
end

@printf("[S2 kernel] measured:\n")
for key in sort!(collect(keys(MEASURED)))
    @printf("[S2 kernel]   %-52s %+.3e\n", key, MEASURED[key])
end

end # module S2KernelTestSuite
