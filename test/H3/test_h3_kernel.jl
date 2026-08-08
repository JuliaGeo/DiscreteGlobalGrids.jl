module H3KernelTestSuite

# Tests for `src/H3/H3Kernel.jl`: the operations-kernel wiring of `H3DGGS`.
#
# The bulk of the file is the mandatory `descendant_range` / ordinal
# verification checklist (see the operation contracts in `src/core/kernel.jl`)
# — H3 only gets `has_descendant_ranges = true`, and the generic cursor only
# gets its `searchsorted` pruning path, if all five items hold — plus the
# cap-validation sweep that measures how far descendants overhang their
# ancestor's cap (the `CELL_CAP_INFLATION = 1.2` factor the kernel applies).

using Test
using Printf
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: H3DGGS, DGGSPartialGrid, DGGSGrid, subtree_grid
using DiscreteGlobalGrids.H3.H3Lookups: H3Lookup
import ConservativeRegridding as CR
import GeometryOps as GO
import GeoInterface as GI
import SparseArrays

const H3N = DGG.H3.H3Native
const S = H3DGGS()
const ROOTS = DGG.root_ids(S)
const RES1 = reduce(vcat, [DGG.cell_children(S, 0, root) for root in ROOTS])

# Cells around each pentagon: the children of the pentagon's parent, i.e. the
# pentagon itself plus its 5 same-parent siblings — the corner of the H3
# hierarchy where digit sequences are deleted and distortion is largest.
function pentagon_neighborhood(res)
    parents = unique(H3N.cell_to_parent(pentagon, res - 1) for pentagon in H3N.get_pentagons(res))
    return sort!(reduce(vcat, [DGG.cell_children(S, res - 1, parent) for parent in parents]))
end

# ---------------------------------------------------------------------------
# Checklist items 3 and 4: the two-sided `descendant_range` contract.
#
# (a) tight endpoints        `extrema(cell_descendants) == descendant_range`
# (b) ordered, disjoint      consecutive sorted parents never overlap
# (c) complete count         Σ subtree_leaf_count == num_cells
# and the strongest single statement, which needs no enumeration at all:
# the number of VALID cells inside `[lo, hi]` (an ordinal difference, since
# ordinals count valid cells in ascending id order) equals the subtree size —
# so no non-descendant valid cell can hide inside the interval.
# ---------------------------------------------------------------------------
function range_report(parents, level, leaf_level; enumerate_descendants::Bool=true)
    tight = true
    ascending = true
    ordered = true
    exact_ordinal_width = true
    total = Int64(0)
    previous_hi = nothing
    for parent in parents
        lo, hi = DGG.descendant_range(S, level, parent, leaf_level)
        count = DGG.subtree_leaf_count(S, level, parent, leaf_level)
        total += count
        if enumerate_descendants
            descendants = DGG.cell_descendants(S, level, parent, leaf_level)
            ascending &= issorted(descendants; lt=(<=))
            tight &= (minimum(descendants) == lo && maximum(descendants) == hi)
            tight &= (length(descendants) == count)
        end
        exact_ordinal_width &=
            (DGG.cell_to_ordinal(S, leaf_level, hi) - DGG.cell_to_ordinal(S, leaf_level, lo) + 1 == count)
        previous_hi === nothing || (ordered &= previous_hi < lo)
        previous_hi = hi
    end
    return (; tight, ascending, ordered, exact_ordinal_width, total)
end

# ---------------------------------------------------------------------------
# CAP-VALIDATION: how far do a cell's `leaf_level` descendants stick out past
# the cell's own boundary? Measured as `max distance from cell_center to any
# descendant boundary vertex / max distance to the cell's own vertices`, both
# taken from the geometry exactly as wired (`cell_center`, cleaned
# `cell_boundary`). A ratio below `CELL_CAP_INFLATION` means the wired
# `cell_cap` already covers the whole subtree.
# ---------------------------------------------------------------------------
function union_ratios(parents, level, deltas)
    ratios = zeros(length(deltas))
    for parent in parents
        center = DGG.cell_center(S, level, parent)
        own = maximum(GO.UnitSpherical.spherical_distance(center, vertex)
                      for vertex in DGG.cell_boundary(S, level, parent))
        for (k, delta) in enumerate(deltas)
            # Angular distance is monotone in the negated dot product, so the
            # inner loop stays free of `acos`.
            smallest_dot = 1.0
            for descendant in DGG.cell_descendants(S, level, parent, level + delta),
                vertex in DGG.cell_boundary(S, level + delta, descendant)

                projection = center[1] * vertex[1] + center[2] * vertex[2] + center[3] * vertex[3]
                projection < smallest_dot && (smallest_dot = projection)
            end
            ratios[k] = max(ratios[k], acos(clamp(smallest_dot, -1.0, 1.0)) / own)
        end
    end
    return ratios
end

@testset "H3 kernel traits and hierarchy" begin
    @test DGG.cell_id_type(S) === UInt64
    @test !DGG.has_ordinal_ids(S)
    @test DGG.has_descendant_ranges(S)
    for res in 0:H3N.MAX_RESOLUTION
        @test DGG.num_cells(S, res) == 2 + 120 * 7^res
        @test DGG.num_cells(S, res) == H3N.num_cells(res)
    end

    # Checklist 1: root ids ascending as returned, base cells 0:121 in order.
    @test length(ROOTS) == 122
    @test ROOTS == H3N.res0_cells()          # asserted, never sorted
    @test issorted(ROOTS; lt=(<=))
    @test [H3N.get_base_cell(root) for root in ROOTS] == collect(0:121)
    @test all(root -> H3N.is_valid_cell(root) && H3N.get_resolution(root) == 0, ROOTS)
    @test DGG.root_ids(S) !== DGG.root_ids(S)     # callers get their own copy

    @test length(RES1) == DGG.num_cells(S, 1)
    @test issorted(RES1; lt=(<=))
    @test all(id -> DGG.cell_parent(S, 1, id, 0) == H3N.cell_to_parent(id, 0), RES1)
    @test DGG.cell_descendants(S, 0, ROOTS[1], 0) == [ROOTS[1]]
    @test_throws ArgumentError DGG.cell_descendants(S, 2, RES1[1], 1)

    # Pentagon subtree sizes: (5 * 7^Δ + 1) / 6 rather than 7^Δ.
    for pentagon in H3N.get_pentagons(1), delta in 1:5
        @test DGG.subtree_leaf_count(S, 1, pentagon, 1 + delta) == (5 * 7^delta + 1) ÷ 6
    end
    for hexagon in RES1[1:20], delta in 1:5
        H3N.is_pentagon(hexagon) && continue
        @test DGG.subtree_leaf_count(S, 1, hexagon, 1 + delta) == 7^delta
    end
end

@testset "H3 checklist 2: children and descendants ascending" begin
    for root in ROOTS, delta in 1:3
        @test issorted(DGG.cell_descendants(S, 0, root, delta); lt=(<=))
    end
    for res in (1, 2, 5)
        pentagons = H3N.get_pentagons(res)
        hexagons = [child for child in DGG.cell_children(S, res - 1, H3N.cell_to_parent(first(pentagons), res - 1))
                    if !H3N.is_pentagon(child)]
        samples = vcat(pentagons, hexagons, DGG.cell_descendants(S, 0, ROOTS[end], res)[1:5])
        for sample in samples, delta in 1:3
            children = DGG.cell_children(S, res, sample)
            @test issorted(children; lt=(<=))
            @test length(children) == (H3N.is_pentagon(sample) ? 6 : 7)
            @test issorted(DGG.cell_descendants(S, res, sample, res + delta); lt=(<=))
        end
    end
end

@testset "H3 checklist 3: two-sided descendant_range contract" begin
    # Complete parent level res 0, exhaustively.
    for leaf_res in 1:3
        report = range_report(ROOTS, 0, leaf_res)
        @test report.tight
        @test report.ascending
        @test report.ordered
        @test report.exact_ordinal_width
        @test report.total == DGG.num_cells(S, leaf_res)
    end
    # Complete parent level res 1, exhaustively.
    for leaf_res in 2:4
        report = range_report(RES1, 1, leaf_res)
        @test report.tight
        @test report.ascending
        @test report.ordered
        @test report.exact_ordinal_width
        @test report.total == DGG.num_cells(S, leaf_res)
    end
    # Pentagon neighborhoods, where the digit deletion makes subtrees smaller
    # than 7^Δ while the range endpoints must stay exactly on the extremes.
    for res in (4, 6), delta in 1:3
        report = range_report(pentagon_neighborhood(res), res, res + delta)
        @test report.tight
        @test report.ascending
        @test report.ordered
        @test report.exact_ordinal_width
    end
    # Every id in the range that is *not* a descendant must be invalid: sample
    # the interval and check directly (the ordinal-width test proves this
    # globally, this makes the failure mode legible).
    parent = H3N.get_pentagons(2)[1]
    lo, hi = DGG.descendant_range(S, 2, parent, 4)
    descendants = Set(DGG.cell_descendants(S, 2, parent, 4))
    for id in lo:UInt64(max(1, (hi - lo) ÷ 500)):hi
        @test (id in descendants) == H3N.is_valid_cell(id)
    end
end

@testset "H3 checklist 4: deep-delta endpoints without enumeration" begin
    for leaf_res in (6, 9, 12, 15)
        previous_hi = nothing
        for root in ROOTS
            lo, hi = DGG.descendant_range(S, 0, root, leaf_res)
            @test H3N.is_valid_cell(lo) && H3N.is_valid_cell(hi)
            @test H3N.get_resolution(lo) == leaf_res && H3N.get_resolution(hi) == leaf_res
            @test H3N.cell_to_parent(lo, 0) == root && H3N.cell_to_parent(hi, 0) == root
            @test lo < hi
            # Next sibling starts strictly after this one ends.
            previous_hi === nothing || @test previous_hi < lo
            previous_hi = hi
        end
    end
    # Same, from deep pentagon/hexagon parents down to the maximum resolution.
    # Deep parents come from point lookups, never from enumerating a subtree.
    points = [(-122.41795, 37.77594), (0.0, 0.0), (139.6917, 35.6895),
        (-58.3816, -34.6037), (18.4241, -33.9249)]
    for res in (5, 9), leaf_res in (res + 1, 12, 15)
        parents = sort!(unique!(vcat(H3N.get_pentagons(res),
            [H3N.lonlat_to_cell(lon, lat, res) for (lon, lat) in points])))
        previous_hi = nothing
        for parent in parents
            lo, hi = DGG.descendant_range(S, res, parent, leaf_res)
            @test H3N.is_valid_cell(lo) && H3N.is_valid_cell(hi)
            @test H3N.cell_to_parent(lo, res) == parent && H3N.cell_to_parent(hi, res) == parent
            @test H3N.get_resolution(hi) == leaf_res
            previous_hi === nothing || @test previous_hi < lo
            previous_hi = hi
        end
    end
    # Degenerate delta: a cell is its own only `leaf_level == level` descendant.
    self = DGG.cell_descendants(S, 0, ROOTS[1], 3)[1]
    @test DGG.descendant_range(S, 3, self, 3) == (self, self)
    # The resolution guard: a cell handed over at the wrong level is rejected
    # rather than silently given an unsound interval.
    @test_throws ArgumentError DGG.descendant_range(S, 1, ROOTS[1], 4)
    @test_throws ArgumentError DGG.descendant_range(S, 4, RES1[1], 2)
    @test_throws ArgumentError DGG.descendant_range(S, 0, ROOTS[1], 16)
end

@testset "H3 checklist 5: ordinal monotonicity" begin
    for res in 0:2
        ids = res == 0 ? ROOTS :
              reduce(vcat, [DGG.cell_descendants(S, 0, root, res) for root in ROOTS])
        @test length(ids) == DGG.num_cells(S, res)
        @test issorted(ids; lt=(<=))
        @test [DGG.cell_to_ordinal(S, res, id) for id in ids] == collect(1:length(ids))
        @test [DGG.ordinal_to_cell(S, res, ordinal) for ordinal in 1:length(ids)] == ids
    end
    # Samples straddling base-cell boundaries and pentagon subtrees at res 5-9.
    for res in 5:9
        boundaries = cumsum([DGG.subtree_leaf_count(S, 0, root, res) for root in ROOTS])
        @test last(boundaries) == DGG.num_cells(S, res)
        ordinals = Int[]
        for boundary in boundaries
            append!(ordinals, (boundary - 1):(boundary + 1))
        end
        for pentagon in H3N.get_pentagons(res)
            ordinal = DGG.cell_to_ordinal(S, res, pentagon)
            append!(ordinals, (ordinal - 1):(ordinal + 1))
            lo, hi = DGG.descendant_range(S, res - 2, H3N.cell_to_parent(pentagon, res - 2), res)
            append!(ordinals, [DGG.cell_to_ordinal(S, res, lo), DGG.cell_to_ordinal(S, res, hi)])
        end
        filter!(ordinal -> 1 <= ordinal <= DGG.num_cells(S, res), ordinals)
        sort!(unique!(ordinals))
        cells = [DGG.ordinal_to_cell(S, res, ordinal) for ordinal in ordinals]
        @test all(cell -> H3N.is_valid_cell(cell) && H3N.get_resolution(cell) == res, cells)
        @test issorted(cells; lt=(<=))                       # monotone in the id
        @test [DGG.cell_to_ordinal(S, res, cell) for cell in cells] == ordinals
    end
    # The kernel-uniform error for an ordinal that names no cell of its level.
    @test_throws DGG.OrdinalRangeError DGG.ordinal_to_cell(S, 2, 0)
    @test_throws DGG.OrdinalRangeError DGG.ordinal_to_cell(S, 2, DGG.num_cells(S, 2) + 1)
    err = try
        DGG.ordinal_to_cell(S, 2, DGG.num_cells(S, 2) + 1)
    catch e
        e
    end
    @test err.system === :H3
    @test err.level == 2
    @test err.total == DGG.num_cells(S, 2)
    @test occursin("H3", sprint(showerror, err))
end

@testset "H3 cap validation" begin
    scopes = [
        ("res 0, all 122 cells", ROOTS, 0, 1:5),
        ("res 1, all 842 cells", RES1, 1, 1:4),
        ("res 1, 24-cell sample", sort!(vcat(H3N.get_pentagons(1), RES1[1:70:842])), 1, 1:5),
        ("res 4, pentagon nbhd", pentagon_neighborhood(4), 4, 1:5),
        ("res 6, pentagon nbhd", pentagon_neighborhood(6), 6, 1:5),
    ]
    worst = 0.0
    for (label, parents, level, deltas) in scopes
        ratios = union_ratios(parents, level, deltas)
        increments = diff(ratios)
        @printf("[H3 cap-validation] %-22s (%3d cells) deltas %d:%d ratios %s increments %s\n",
            label, length(parents), first(deltas), last(deltas),
            join((@sprintf("%.5f", r) for r in ratios), " "),
            join((@sprintf("%.5f", d) for d in increments), " "))
        worst = max(worst, maximum(ratios))
        # (ii) the bound itself, and the geometric envelope that justifies
        # extrapolating past the measured deltas. H3 alternates Class II and
        # Class III between consecutive resolutions, so the increments come as a
        # two-step staircase (on pentagon neighborhoods every other increment is
        # exactly zero); the convergence statement is therefore per *two*
        # deltas, which still gives a geometric — ratio 1/sqrt(2) — tail.
        @test maximum(ratios) <= 1.10
        @test all(>=(-1e-12), increments)                    # ratios never fall back
        @test all(k -> increments[k + 2] <= increments[k] / 2 + 1e-12,
            1:(length(increments) - 2))
        # Tail of that envelope beyond the last measured delta, summed to
        # infinity: Σ_{j>0} inc[K+j] <= (inc[K-1] + inc[K]) * (1/2 + 1/4 + ...).
        extrapolated = ratios[end] + increments[end - 1] + increments[end]
        @test extrapolated < DGG.CELL_CAP_INFLATION * 0.95
        @printf("[H3 cap-validation] %-22s extrapolated supremum %.5f\n", label, extrapolated)
        # And the invariant that actually matters: descendants are inside the
        # wired cap, whose radius is `own * CELL_CAP_INFLATION`.
        @test maximum(ratios) < DGG.CELL_CAP_INFLATION
    end
    @printf("[H3 cap-validation] worst union ratio over all scopes = %.5f (cap inflation %.2f)\n",
        worst, DGG.CELL_CAP_INFLATION)
    @test worst <= 1.10

    # Direct containment against the caps as the trees consume them.
    for level in (0, 4), parent in (level == 0 ? ROOTS : pentagon_neighborhood(4)), delta in 1:2
        cap = DGG.cell_cap(S, level, parent)
        @test cap.radius <= Float64(pi)
        for descendant in DGG.cell_descendants(S, level, parent, level + delta),
            vertex in DGG.cell_boundary(S, level + delta, descendant)

            @test GO.UnitSpherical.spherical_distance(cap.point, vertex) <= cap.radius
        end
    end
end

@testset "H3 cleaned cell boundaries" begin
    samples = UInt64[]
    append!(samples, ROOTS)
    append!(samples, H3N.get_pentagons(0))
    for res in (1, 2, 4), pentagon in H3N.get_pentagons(res)
        push!(samples, pentagon)
        append!(samples, DGG.cell_children(S, res, pentagon))
    end
    append!(samples, DGG.cell_descendants(S, 0, ROOTS[end], 2)[1:40])
    # A known face-crossing cell: 7 raw vertices, 6 after cleaning.
    face_crossing = parse(UInt64, "8101bffffffffff"; base=16)
    push!(samples, face_crossing)
    push!(samples, H3N.lonlat_to_cell(-122.41795063018799, 37.775938728915946, 9))
    push!(samples, H3N.lonlat_to_cell(0.0, 0.0, 15))
    unique!(samples)

    n_cleaned = 0
    worst_intact_area = 0.0
    worst_cleaned_area = 0.0
    for id in samples
        res = H3N.get_resolution(id)
        raw = H3N.cell_boundary_cartesian(id; closed_ring=false)
        boundary = DGG.cell_boundary(S, res, id)
        @test length(boundary) == (H3N.is_pentagon(id) ? 5 : 6)
        cleaned = length(raw) > length(boundary)
        n_cleaned += cleaned
        # Cleaning only ever drops vertices — it never invents them.
        raw_points = [GO.UnitSphericalPoint(p[1], p[2], p[3]) for p in raw]
        @test all(point -> point in raw_points, boundary)
        closed = DGG.cell_boundary(S, res, id; closed=true)
        @test length(closed) == length(boundary) + 1
        @test closed[1:(end - 1)] == boundary
        @test closed[end] == closed[1]
        @test DGG.cell_boundary(S, res, id) == boundary   # not mutated by `closed`
        polygon = DGG.cell_polygon_unitsphere(S, res, id)
        @test collect(GI.getpoint(polygon)) == closed
        # Area against H3's own `cellAreaRads2`. Rings that needed no cleaning
        # must match to round-off; rings that lost distortion vertices lose the
        # slivers those vertices carried (worst case: pentagons, whose raw ring
        # is 10 vertices for 5 corners). That is inherited from the boundary
        # pipeline the kernel took over from the old per-system H3 tree, and is
        # reported, not asserted away.
        error = abs(abs(GO.area(GO.Spherical(; radius=1.0), polygon)) / H3N.cell_area(id) - 1)
        if cleaned
            worst_cleaned_area = max(worst_cleaned_area, error)
        else
            worst_intact_area = max(worst_intact_area, error)
        end
    end
    @test worst_intact_area < 1e-9
    @test worst_cleaned_area < 0.12
    @printf("[H3 boundary] cleaned %d / %d sampled cells; area error vs cellAreaRads2: %.2e (intact) %.4f (cleaned)\n",
        n_cleaned, length(samples), worst_intact_area, worst_cleaned_area)
    @test length(H3N.cell_boundary_cartesian(face_crossing; closed_ring=false)) == 7
    @test length(DGG.cell_boundary(S, 1, face_crossing)) == 6
    @test n_cleaned > 0                                   # the pipeline is exercised

    # Native center, not the boundary mean.
    for id in samples[1:20]
        res = H3N.get_resolution(id)
        lon, lat = H3N.cell_center(id)
        center = DGG.cell_center(S, res, id)
        @test center[3] ≈ sin(deg2rad(lat)) atol = 1e-14
        @test sum(abs2, center) ≈ 1.0 atol = 1e-12
        @test atan(center[2], center[1]) ≈ deg2rad(lon) atol = 1e-12
    end
end

@testset "H3 grid convenience constructors" begin
    parent = H3N.lonlat_to_cell(-122.41795063018799, 37.775938728915946, 1)
    ids = DGG.cell_descendants(S, 1, parent, 3)
    lookup = H3Lookup(ids; resolution=3, validate=true)
    grid = DGGSPartialGrid(lookup)
    @test grid.system === S
    @test grid.level == 3
    @test grid.ids === lookup.data                        # by reference, not collected
    @test grid.bucket_size == 0
    @test DGGSPartialGrid(lookup; bucket_size=8).bucket_size == 8
    @test eltype(grid.ids) === DGG.cell_id_type(S)

    chunk = subtree_grid(S, parent; root_level=1, leaf_level=3)
    @test chunk.ids == ids
    @test chunk.root_level == 1 && chunk.root_id == parent
    @test length(chunk.ids) == DGG.subtree_leaf_count(S, 1, parent, 3)

    dense = DGGSGrid(S, 2)
    @test DGG.num_cells(dense.system, dense.level) == 5882
    @test_throws ArgumentError DGGSGrid(S, 16)
end

# ---------------------------------------------------------------------------
# `ConservativeRegridding.Regridder` over the generic grids — the consumer the
# tree layer exists for. Self-regrids, so the answer is known exactly: every
# cell meets only itself, with its own area.
# ---------------------------------------------------------------------------
@testset "H3 Regridder round trips" begin
    dense = CR.Regridder(DGGSGrid(S, 0), DGGSGrid(S, 0); threaded=false, normalize=false)
    @test size(dense.intersections) == (122, 122)
    @test SparseArrays.nnz(dense.intersections) == 122         # diagonal only
    @test isapprox(sum(dense.intersections), sum(dense.dst_areas); rtol=1e-10)

    parent = H3N.lonlat_to_cell(-122.41795063018799, 37.775938728915946, 1)
    lookup = H3Lookup(DGG.cell_descendants(S, 1, parent, 2); resolution=2, validate=true)
    partial = CR.Regridder(DGGSPartialGrid(lookup), DGGSPartialGrid(lookup);
        threaded=false, normalize=false)
    @test size(partial.intersections) == (length(lookup), length(lookup))
    @test isapprox(sum(partial.intersections), sum(partial.dst_areas); rtol=1e-10)

    # The antimeridian seam, which main's H3 suite covered and which the
    # distortion-vertex cleanup above is what makes clippable: `cellToBoundary`
    # hands this cell 7 vertices, the wired ring 6, and only the cleaned convex
    # ring self-intersects to its own area.
    face_crossing = parse(UInt64, "8101bffffffffff"; base=16)
    @test length(H3N.cell_boundary_cartesian(face_crossing; closed_ring=false)) == 7
    @test length(DGG.cell_boundary(S, 1, face_crossing)) == 6
    seam_lookup = H3Lookup([face_crossing]; resolution=1, validate=true)
    seam = CR.Regridder(DGGSPartialGrid(seam_lookup), DGGSPartialGrid(seam_lookup);
        threaded=false, normalize=false)
    @test seam.intersections[1, 1] ≈ only(seam.dst_areas)
end

end # module H3KernelTestSuite
