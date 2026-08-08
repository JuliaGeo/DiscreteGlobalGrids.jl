module HealpixKernelTestSuite

# Tests for `src/HEALPix/HealpixKernel.jl`: the operations-kernel wiring of
# `HEALPixDGGS`. Ground truth is Healpix.jl itself (`boundariesRing`,
# `pix2angNest`, `nside2npix`) plus the nested-HEALPix id arithmetic written
# out locally below, so every cell is checked against the convention, not
# against another implementation in this package.
#
# Since milestone 3 the kernel's geometry is evaluated by `src/HEALPix/chart.jl`
# rather than by Healpix.jl, which splits the geometry tests in two:
#
#   * "geometry vs Healpix.jl" is now a *cross-check* — two independent
#     implementations of the same closed forms, so agreement is to round-off
#     (atol 1e-14) and the measured deviations are printed. It stays for as
#     long as Healpix.jl is a dependency, because it is the only assertion here
#     that reaches outside the package for its truth.
#   * "geometry is the chart kernel, bitwise" is the exactness claim, and it is
#     stronger than the `==`-against-Healpix.jl it replaced: the id hierarchy
#     and the dense face grid (`face_grid.jl`) must emit *literally the same*
#     `Float64`s, because they now share one evaluation.

using Test
using Printf
using Random
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: HEALPixDGGS, H3DGGS, DGGSPartialGrid, DGGSGrid, subtree_grid
import DiscreteGlobalGrids.HEALPix
using DiscreteGlobalGrids.HEALPix.HealpixLookups: HealpixLookup
import ConservativeRegridding as CR
import ConservativeRegridding: Trees
import GeometryOps as GO
import GeometryOps: SpatialTreeInterface as STI
import GeoInterface as GI
import Healpix
import SparseArrays

const S = HEALPixDGGS()

resolution(level) = Healpix.Resolution(2^level)

# Ground-truth corners straight from Healpix.jl (0-based nested -> 1-based ring).
function truth_corners(level, pixel)
    res = resolution(level)
    cartesian = Healpix.boundariesRing(res, Healpix.nest2ring(res, pixel + 1), 1, Float64)
    return [GO.UnitSphericalPoint(cartesian[i, 1], cartesian[i, 2], cartesian[i, 3]) for i in 1:4]
end

# The pixel outline sampled far more densely than the four corners: what the
# caps have to bound, and what descendant boundaries converge to.
function truth_outline(level, pixel; step::Int=64)
    res = resolution(level)
    cartesian = Healpix.boundariesRing(res, Healpix.nest2ring(res, pixel + 1), step, Float64)
    return [GO.UnitSphericalPoint(cartesian[i, 1], cartesian[i, 2], cartesian[i, 3])
            for i in 1:size(cartesian, 1)]
end

# Per-face first/middle/last pixels (face corners and the base-tiling seams),
# plus the two polar pixels and a ring of equator pixels.
function sample_pixels(level)
    per_face = 4^level
    offsets = unique(filter(o -> 0 <= o < per_face, [0, 1, per_face ÷ 2, per_face - 1]))
    pixels = Int64[]
    for face in 0:11, offset in offsets
        push!(pixels, face * per_face + offset)
    end
    res = resolution(level)
    push!(pixels, Healpix.ang2pixNest(res, 0.0, 0.0) - 1)
    push!(pixels, Healpix.ang2pixNest(res, Float64(pi), 0.0) - 1)
    for k in 0:7
        push!(pixels, Healpix.ang2pixNest(res, Float64(pi) / 2, k * Float64(pi) / 4) - 1)
    end
    return sort!(unique!(pixels))
end

@testset "HEALPix kernel traits" begin
    @test DGG.cell_id_type(S) === Int64
    @test DGG.has_ordinal_ids(S)
    @test DGG.has_descendant_ranges(S)
    # Parents contain their children here and the 4-corner cap is O(1), so
    # partial internal nodes want `subtree_cap`, not a stored-id union cap.
    # Aperture-7 systems keep the default: their parent caps are 1.2-inflated.
    @test DGG.has_exact_subtree_cap(S)
    @test !DGG.has_exact_subtree_cap(H3DGGS())
    # `num_cells` is the kernel default (root_count * radix^level); the point of
    # the check is that the default already equals Healpix's own count.
    for level in 0:10
        @test DGG.num_cells(S, level) == Healpix.nside2npix(Healpix.order2nside(level))
        @test DGG.num_cells(S, level) == 12 * 4^level
    end
    @test DGG.root_ids(S) == collect(Int64, 0:11)
end

# Cross-check against the external implementation. These were bitwise `==`
# while `HealpixKernel.jl` *called* `boundariesRing` / `pix2vecNest`; they are
# tolerances now that the kernel evaluates `chart.jl` instead, because the two
# implementations do the same arithmetic in a different order. The corner ORDER
# is still asserted exactly — a permuted ring would blow past atol by ~1e0
# (measured: the three non-identity cyclic shifts deviate by 0.75, 1.41, 0.75),
# so this doubles as the "same ring `boundariesRing` emits" check that
# order-sensitive downstream consumers depend on.
@testset "HEALPix kernel geometry vs Healpix.jl" begin
    max_corner = 0.0
    max_center = 0.0
    max_cap_point = 0.0
    max_cap_radius = 0.0
    for level in 0:8, pixel in sample_pixels(level)
        corners = truth_corners(level, pixel)
        boundary = DGG.cell_boundary(S, level, pixel)
        @test length(boundary) == 4
        for i in 1:4
            max_corner = max(max_corner, maximum(abs, boundary[i] - corners[i]))
        end
        @test all(i -> isapprox(boundary[i], corners[i]; atol=1e-14), 1:4)
        closed = DGG.cell_boundary(S, level, pixel; closed=true)
        @test length(closed) == 5
        @test closed[1:4] == boundary                # same call, so still exact
        @test closed[end] == closed[1]

        theta, phi = Healpix.pix2angNest(resolution(level), pixel + 1)
        center = DGG.cell_center(S, level, pixel)
        max_center = max(max_center, maximum(abs,
            center - GO.UnitSphericalPoint(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta))))
        @test center[1] ≈ sin(theta) * cos(phi) atol = 1e-14
        @test center[2] ≈ sin(theta) * sin(phi) atol = 1e-14
        @test center[3] ≈ cos(theta) atol = 1e-14
        @test sum(abs2, center) ≈ 1.0 atol = 1e-12

        # Caps are the exact 4-corner caps.
        cap = DGG.cell_cap(S, level, pixel)
        truth_cap = Trees.circle_from_four_corners(corners, ())
        max_cap_point = max(max_cap_point, maximum(abs, cap.point - truth_cap.point))
        max_cap_radius = max(max_cap_radius, abs(cap.radius - truth_cap.radius))
        @test isapprox(cap.point, truth_cap.point; atol=1e-14)
        @test cap.radius ≈ truth_cap.radius atol = 1e-14
        @test DGG.subtree_cap(S, level, pixel, level + 4) == cap
        # Ground truth for "exact": the densely sampled pixel outline is inside.
        # `truth_outline` is Healpix.jl's own densification, so this keeps
        # bounding the *external* notion of the pixel, not the chart's.
        @test all(point -> GO.UnitSpherical.spherical_distance(cap.point, point) <= cap.radius,
            truth_outline(level, pixel))
        @test GO.UnitSpherical.spherical_distance(cap.point, center) <= cap.radius

        polygon = DGG.cell_polygon_unitsphere(S, level, pixel)
        @test collect(GI.getpoint(polygon)) == [boundary; [boundary[1]]]
    end
    @printf("[HEALPix chart vs Healpix.jl] max |Δ| per coordinate: corners = %.3e, center = %.3e, cap.point = %.3e, cap.radius = %.3e\n",
        max_corner, max_center, max_cap_point, max_cap_radius)
end

# ---------------------------------------------------------------------------
# The milestone-3 unification claim: the id-hierarchy path and the dense face
# grid are ONE evaluation, not two that agree. Every assertion here is bitwise
# (`===` per coordinate) — an ulp of drift means the two paths have forked and
# `Regridder`s built on either would stop being entry-for-entry interchangeable
# (which is exactly what `test_face_grid.jl` asserts with `==` on the sparse
# matrices).
# ---------------------------------------------------------------------------

# Bit-level equality, not `==`: `==` would accept `0.0` for `-0.0`, and a
# signed zero that differs between the two paths is already a fork.
identical(a, b) = all(k -> a[k] === b[k], 1:3)

@testset "HEALPix kernel geometry is the chart kernel, bitwise" begin
    for level in (0, 1, 2, 3, 5)
        nside = 2^level
        count = 12 * 4^level
        # Levels 0-2 exhaustively; above that the per-face/seam/pole sample plus
        # a seeded random draw (seeded so a failure is reproducible verbatim).
        pixels = level <= 2 ? collect(Int64, 0:(count - 1)) :
            sort!(unique!(vcat(sample_pixels(level),
                Int64.(rand(MersenneTwister(20260805 + level), 0:(count - 1), 256)))))
        # The face-grid representation of the same resolution. `NestedOrder`
        # puts nested id `p` at data position `p + 1`, which is also the DGGS
        # ordinal, so no permutation sits between the two sides.
        root = Trees.treeify(GO.Spherical(), HEALPix.HealpixFaceGrid(nside; ordering=HEALPix.NestedOrder()))
        for pixel in pixels
            chart = HEALPix.pixel_corners(HEALPix.nested_to_xyf(pixel, nside)..., nside)
            boundary = DGG.cell_boundary(S, level, pixel)
            @test length(boundary) == 4
            @test all(i -> identical(boundary[i], chart[i]), 1:4)
            closed = DGG.cell_boundary(S, level, pixel; closed=true)
            @test all(i -> identical(closed[i], chart[i]), 1:4)
            @test identical(closed[5], chart[1])

            @test identical(DGG.cell_center(S, level, pixel),
                HEALPix.pixel_center(HEALPix.nested_to_xyf(pixel, nside)..., nside))

            # Caps: same `circle_from_four_corners` call on the same corners.
            @test DGG.cell_cap(S, level, pixel) ===
                  Trees.circle_from_four_corners(chart, ())

            # The polygon rings — what `ConservativeRegridding` actually clips.
            kernel_ring = collect(GI.getpoint(DGG.cell_polygon_unitsphere(S, level, pixel)))
            face_ring = collect(GI.getpoint(Trees.getcell(root, Int(pixel) + 1)))
            @test length(kernel_ring) == length(face_ring) == 5
            @test all(i -> identical(kernel_ring[i], face_ring[i]), eachindex(kernel_ring))

            # The interface-level `cell_polygon(::AbstractDGGS, level, id)` is a
            # pure bridge to `cell_polygon_unitsphere` (wired alongside the S2
            # and ISEA4R ones so that no system with fully wired geometry still
            # throws `NotPortedError` at the interface).
            @test collect(GI.getpoint(DGG.cell_polygon(S, level, pixel))) == kernel_ring
        end
    end
end

@testset "HEALPix kernel hierarchy vs nested-pixel arithmetic" begin
    for level in 0:6, pixel in sample_pixels(level)
        @test DGG.cell_children(S, level, pixel) == collect(4pixel:(4pixel + 3))
        @test DGG.cell_to_ordinal(S, level, pixel) == pixel + 1
        @test DGG.ordinal_to_cell(S, level, pixel + 1) == pixel
        for delta in 1:4
            leaf_level = level + delta
            # The nested convention's half-open descendant interval [lo, hi).
            lo = Int64(pixel) * Int64(4)^delta
            hi = lo + Int64(4)^delta
            @test DGG.descendant_range(S, level, pixel, leaf_level) == (lo, hi - 1)
            @test DGG.subtree_leaf_count(S, level, pixel, leaf_level) == hi - lo
            @test DGG.subtree_leaf_count(S, level, pixel, leaf_level) == 4^delta
            descendants = DGG.cell_descendants(S, level, pixel, leaf_level)
            @test descendants == collect(lo:(hi - 1))
            @test issorted(descendants; lt=(<=))
            @test all(id -> DGG.cell_parent(S, leaf_level, id, level) == pixel, descendants)
        end
    end
    # Ordinals are dense and monotone at every level: ids *are* the ordinals.
    for level in 0:4
        ids = collect(Int64, 0:(DGG.num_cells(S, level) - 1))
        @test [DGG.cell_to_ordinal(S, level, id) for id in ids] == collect(1:length(ids))
        @test [DGG.ordinal_to_cell(S, level, o) for o in 1:length(ids)] == ids
    end
end

@testset "HEALPix subtree caps contain descendants" begin
    # Exhaustive over levels 0:3, then sampled at 4:6: every vertex of every
    # descendant must sit inside the parent's exact 4-corner cap. This is what
    # justifies overriding `subtree_cap` with `cell_cap`.
    worst = 0.0
    for level in 0:3, pixel in 0:(DGG.num_cells(S, level) - 1), delta in 1:3
        cap = DGG.subtree_cap(S, level, pixel, level + delta)
        for descendant in DGG.cell_descendants(S, level, pixel, level + delta),
            point in DGG.cell_boundary(S, level + delta, descendant)

            worst = max(worst, GO.UnitSpherical.spherical_distance(cap.point, point) / cap.radius)
        end
    end
    @test worst <= 1.0
    for level in 4:6, pixel in sample_pixels(level), delta in 1:3
        cap = DGG.subtree_cap(S, level, pixel, level + delta)
        for descendant in DGG.cell_descendants(S, level, pixel, level + delta),
            point in DGG.cell_boundary(S, level + delta, descendant)

            worst = max(worst, GO.UnitSpherical.spherical_distance(cap.point, point) / cap.radius)
        end
    end
    @test worst <= 1.0
    @printf("[HEALPix cap-validation] worst descendant-vertex / subtree_cap radius ratio = %.6f\n", worst)

    # The exact 4-corner cap is barely looser than the generic union cap over
    # the same subtree (which recenters on the descendant vertices), so the
    # override costs no meaningful pruning power while saving the O(4^Δ)
    # boundary evaluations the union cap needs.
    looseness = 0.0
    for level in 0:2, pixel in 0:(DGG.num_cells(S, level) - 1)
        cap = DGG.subtree_cap(S, level, pixel, level + 3)
        union_cap = DGG.cells_cap(S, level + 3, DGG.cell_descendants(S, level, pixel, level + 3))
        looseness = max(looseness, cap.radius / union_cap.radius)
    end
    @test looseness <= 1.10
    @printf("[HEALPix cap-validation] exact 4-corner cap / union cap radius, worst = %.6f\n", looseness)
end

@testset "HEALPix grid convenience constructors" begin
    level = 5
    ids = Int64[0, 1, 2, 3, 100, 4095, 4096, 12 * 4^level - 1]
    lookup = HealpixLookup(ids; level)
    grid = DGGSPartialGrid(lookup)
    @test grid.system === S
    @test grid.level == level
    @test grid.ids === lookup.data          # stored by reference, never copied
    @test grid.bucket_size == 0
    @test grid.root_level == -1
    bucketed = DGGSPartialGrid(lookup; bucket_size=4)
    @test bucketed.bucket_size == 4
    @test bucketed.ids === lookup.data

    # Round trip: grid ids are exactly the lookup's data order, so leaf index i
    # of a tree over this grid is lookup position i.
    @test [grid.ids[i] for i in eachindex(lookup.data)] == lookup.data
    @test DGG.cell_polygon_unitsphere(S, grid.level, grid.ids[2]) isa GI.Polygon

    dense = DGGSGrid(S, 3)
    @test DGG.num_cells(dense.system, dense.level) == 768

    lo, hi = DGG.descendant_range(S, 1, 5, 4)
    chunk = subtree_grid(S, 5; root_level=1, leaf_level=4)
    @test chunk.ids == collect(lo:hi)
    @test length(chunk.ids) == 4^3
    @test chunk.root_level == 1 && chunk.root_id == 5
end

@testset "HEALPix partial node extents use the exact subtree cap" begin
    # One face's worth of sparsely stored ids: the internal node above them is
    # exactly where the old per-system HEALPix node used its O(1) pixel cap and
    # the generic cursor used to fall back to an O(stored) union cap.
    leaf_level = 5
    face = 3
    lo, hi = DGG.descendant_range(S, 0, face, leaf_level)
    ids = collect(lo:7:hi)
    grid = DGGSPartialGrid(S, leaf_level, ids)
    tree = Trees.treeify(GO.Spherical(), grid)
    @test STI.isspatialtree(typeof(tree))
    children = collect(STI.getchild(tree))
    @test length(children) == 1                     # only face 3 stores anything
    node = only(children)                           # level 0, internal
    @test !STI.isleaf(node)
    extent = STI.node_extent(node)
    exact = DGG.subtree_cap(S, 0, face, leaf_level)
    # Sound on either path: a node extent must bound every cell it stores.
    @test all(ids) do id
        all(vertex -> GO.UnitSpherical.spherical_distance(extent.point, vertex) <= extent.radius,
            DGG.cell_boundary(S, leaf_level, id))
    end
    # The trait branch is in: internal partial nodes report the exact O(1)
    # parent-pixel cap, nothing looser (a regression here must FAIL, not warn).
    @test extent == exact
end

# ---------------------------------------------------------------------------
# `ConservativeRegridding.Regridder` over the generic grids. The `HealpixLookup`
# route (`Regridder(l, l)`) is covered in `test/HEALPix/runtests.jl`; these are
# the same round trips entered through `DGGSGrid` / `DGGSPartialGrid`.
# ---------------------------------------------------------------------------
@testset "HEALPix Regridder round trips" begin
    dense = CR.Regridder(DGGSGrid(S, 1), DGGSGrid(S, 1); threaded=false, normalize=false)
    @test size(dense.intersections) == (48, 48)
    @test SparseArrays.nnz(dense.intersections) == 48            # diagonal only
    @test isapprox(sum(dense.intersections), sum(dense.dst_areas); rtol=1e-10)
    @test isapprox(sum(dense.dst_areas), 4pi * GO.Spherical().radius^2; rtol=1e-9)

    ids = Int64[0, 1, 2, 3, 16, 17, 18, 19]
    grid = DGGSPartialGrid(S, 2, ids)
    partial = CR.Regridder(grid, grid; threaded=false, normalize=false)
    @test size(partial.intersections) == (length(ids), length(ids))
    @test SparseArrays.nnz(partial.intersections) == length(ids)
    @test isapprox(sum(partial.intersections), sum(partial.dst_areas); rtol=1e-10)
end

end # module HealpixKernelTestSuite
