# The sibling script's question without the sibling's inputs:
#
#   julia --project=<integration-env> test/integration/geomorphometry_synthetic.jl
#
# `geomorphometry_igeo7.jl` needs a Copernicus GeoTIFF, ArchGDAL and a
# conservative regridder, so it can only be run by hand on a machine that has
# all three. Everything it proves about the DGGS side, though, is independent of
# where the elevations came from: what is under test is the cell iterator, the
# neighbour ring, `cellposition` addressing, and the relative-cell direction
# codec the flow routers encode with. So this file synthesises the field
# instead, and needs only an environment holding `Geomorphometry` (the
# `feat/generic` rev that `docs/Project.toml` pins) and `Rasters` alongside this
# package.
#
# It runs the whole battery over TWO cell sets, because they are different code
# underneath and only one of them is easy:
#
#   * one rooted subtree (`PartialGrid`), which compresses to a SINGLE position
#     window. `cellposition` is then one subtraction and `halo` delegates to the
#     system's own `SubtreeHaloIterator`. This is the shape `hydrology.jl` ends
#     up with, and it is the fast path in every direction.
#   * a `MultiOrderCoverage` of a real 1°×1° extent expanded to a working level,
#     which is a few hundred DISJOINT windows. `cellposition` is a binary search
#     over them, membership likewise, and nothing can be recognised as a subtree.
#     This is the shape an arbitrary area of interest produces, and none of the
#     fast paths above apply to it.
#
# A metric that works on the first and not the second is the normal way this
# integration breaks, and testing only the first would not show it.

using Test
import DiscreteGlobalGrids as DGG
import Geomorphometry as GM
using Rasters

const Extents = Rasters.Extents          # a Rasters dep, so no extra import

sys = DGG.IGeo7System()

# The DEM tile `hydrology.jl` works over, and the same two calls it uses to pick
# a cell to work in.
tile = Extents.Extent(X=(10.0, 11.0), Y=(46.0, 47.0))
coverage = DGG.query(sys, DGG.MultiOrderCoverage(tile); level=10)
root = DGG.coarsest_contained(coverage)

"""
A single dome with a small ripple on it, sampled at cell centroids. The dome
gives the flow routers one real drainage tree instead of a plateau of ties; the
ripple keeps the D8 tie-break off the degenerate path where every neighbour is
equidistant in elevation.
"""
function make_dem(cells)
    complete = DGG.levelgrid(sys, DGG.level(cells))
    apex = DGG.cell_centroid(DGG.levelgrid(sys, DGG.level(root)), root)
    elevation = [
        begin
            p = DGG.cell_centroid(complete, c)
            d = acos(clamp(p[1] * apex[1] + p[2] * apex[2] + p[3] * apex[3], -1, 1))
            1000.0 * exp(-(d / 0.02)^2) + 40.0 * sin(120 * p[1]) * cos(97 * p[2])
        end
        for c in cells
    ]
    return Raster(elevation, (DGG.Cells(DGG.CellLookup(cells)),); name=:height)
end

function battery(label, cells)
    dem = make_dem(cells)
    complete = DGG.levelgrid(sys, DGG.level(cells))

    @testset "$label" begin
        @testset "the indexing surface the extension stands on" begin
            indices = eachindex(dem)
            @test indices === cells
            @test indices isa DGG.CellVector
            c = first(indices)
            @test c in indices
            @test dem[c] == parent(dem)[1]

            # `neighbors` on the subset is the system's answer clipped to
            # membership, so it is a SUBSET of the complete level's answer and
            # never a different set.
            clipped = GM.neighbors(dem, c)
            @test !isempty(clipped)
            @test all(in(indices), clipped)
            @test clipped ⊆ DGG.neighbors(complete, c)

            outside = DGG.cellindex(complete, DGG.ncells(complete))
            @test outside ∉ indices
            # String form on purpose: the type form never renders the message,
            # so it cannot catch an unprintable error.
            @test_throws "at index [Z7Cell" dem[outside]
        end

        # Neither the axis nor the ring may put anything on the heap. The axis is
        # compressed windows, so iterating it computes ids rather than reading
        # them; the ring is clipped into the same static container the system
        # answered in. Both are per-cell costs in every loop Geomorphometry runs,
        # and a `Vector` in either place is invisible to every other assertion
        # here — it returns the right cells, just not for free.
        @testset "the hot path allocates nothing" begin
            axisloop(cs) = (n = 0; for c in cs
                n += 1
            end; n)
            ringloop(sub, cs) = (n = 0; for c in cs, _ in DGG.neighbors(sub, c)
                n += 1
            end; n)
            lookup = DGG.CellLookup(cells)
            axisloop(cells)
            ringloop(cells, cells)
            ringloop(lookup, cells)
            @test @allocated(axisloop(cells)) == 0
            @test @allocated(ringloop(cells, cells)) == 0
            @test @allocated(ringloop(lookup, cells)) == 0
        end

        @testset "local metrics" begin
            tpi = GM.topographic_position_index(dem)
            @test size(tpi) == size(dem)
            @test all(isfinite, parent(tpi))

            rough = GM.roughness(dem)
            @test size(rough) == size(dem)
            @test all(>=(0), parent(rough))
        end

        # A MEMBER neighbour: these three take two cells and resolve both through
        # `cellposition`, so on a scattered set the complete level's first
        # neighbour is frequently outside and the call is a `BoundsError` by
        # design. Only the subtree case is forgiving enough to hide that.
        @testset "cell geometry" begin
            neighbor = first(GM.neighbors(dem, first(cells)))
            @test GM.cellarea(dem, first(cells)) > 0
            @test GM.celldistance(dem, first(cells), neighbor) > 0
            @test 0 <= GM.cellbearing(dem, first(cells), neighbor) < 360
            @test GM.cellbearing(dem, first(cells), first(cells)) == 0
        end

        cell_area = GM.cellarea(dem, first(cells))

        @testset "flow routing: $(nameof(typeof(method)))" for method in
                                                               (GM.D8(), GM.DInf(), GM.FD8())
            accumulation, directions = GM.flowaccumulation(dem; method)
            @test size(accumulation) == size(directions) == size(dem)
            @test all(isfinite, parent(accumulation))
            @test all(>=(0.99cell_area), parent(accumulation))

            hand = GM.height_above_nearest_drainage(dem; method, threshold=5cell_area)
            @test size(hand) == size(dem)
            @test all(isfinite, parent(hand))
            @test all(>=(0), parent(hand))

            # The codec closing: every direction the router wrote back decomposes
            # into relative cells that land on real neighbours. A wrong entry in
            # either conversion table shows up here and nowhere in the metrics
            # above, which never read a direction back out.
            bad = 0
            for (i, c) in enumerate(cells)
                iszero(directions[i]) && continue
                for rel in GM.decompose(DGG.RelativeZ7Cell, directions[i], c)
                    target = c + rel
                    target == c || target in DGG.neighbors(complete, c) || (bad += 1)
                end
            end
            @test bad == 0
        end

        @test !isempty(GM.outlets(dem))
    end
    return dem
end

nwindows(cells) = DGG.Fallbacks.nwindows(DGG.Fallbacks.windows(cells))

@testset "Geomorphometry on a synthetic IGeo7 field" begin
    subtree = DGG.CellVector(DGG.PartialGrid(sys, root, DGG.level(root) + 3))
    scattered = DGG.CellVector(coverage; level=10)

    # The two shapes really are different underneath; without this the second
    # case could silently become a rerun of the first.
    @test nwindows(subtree) == 1
    @test nwindows(scattered) > 100
    @test DGG.level(root) == 5

    battery("one rooted subtree ($(length(subtree)) cells, 1 window)", subtree)
    battery("multi-order coverage ($(length(scattered)) cells, $(nwindows(scattered)) windows)",
        scattered)

    # The subset's rim is where a clipped ring makes a local metric answer from
    # fewer samples than the interior. This is not a defect of the walk, it is
    # what `halo` exists to supply — pinned on the small case so that a change in
    # either direction is visible.
    @testset "the rim a halo would complete" begin
        complete = DGG.levelgrid(sys, DGG.level(subtree))
        clipped = count(subtree) do c
            length(DGG.neighbors(subtree, c)) < length(DGG.neighbors(complete, c))
        end
        @test length(subtree) == 343
        @test clipped == 78
        @test length(collect(DGG.halo(subtree))) == 84
    end
end
