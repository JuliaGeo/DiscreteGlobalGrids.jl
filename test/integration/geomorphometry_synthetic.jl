# The sibling script's question without the sibling's inputs:
#
#   julia --project=<integration-env> test/integration/geomorphometry_synthetic.jl
#
# `geomorphometry_igeo7.jl` needs a Copernicus GeoTIFF, ArchGDAL and a
# conservative regridder, so it can only be run by hand on a machine that has
# all three. Everything it proves about the DGGS side, though, is independent of
# where the elevations came from: what is under test is the cell iterator, the
# neighbour ring, `cellposition` addressing, and the direction codecs the flow
# routers encode with. So this file synthesises the field instead, and needs
# only an environment holding `Geomorphometry` (the `clipped-neighbors` rev of
# the fork at https://github.com/asinghvi17/GeoArrayOps.jl.git) and `Rasters`
# alongside this package.
#
# The battery runs over FOUR systems, because the extension serves two
# different backends under one surface:
#
#   * IGeo7 has relative-cell arithmetic, so it is the full-featured backend:
#     LDD directions persisted as `RelativeZ7Cell` codes, DInf/FD8 routing,
#     and height-above-nearest-drainage.
#   * HEALPix, ISEA4R and H3 reach everything through the system-generic
#     `Cells`-lookup path: D8 directions persist the downstream neighbour's
#     RING SLOT (bit `k-1` of a `UInt16` for slot `k`), and DInf/FD8 refuse
#     with an error that says which backend could have answered. (ISEA4R is
#     the ISEA family's aperture-4 grid; there is no separate hexagon system.)
#
# And it runs each system over TWO cell sets, because they are different code
# underneath and only one of them is easy:
#
#   * one rooted subtree (`PartialGrid`), which compresses to a SINGLE position
#     window. `cellposition` is then one subtraction and `halo` delegates to the
#     system's own `SubtreeHaloIterator`. This is the shape `hydrology.jl` ends
#     up with, and it is the fast path in every direction.
#   * a `MultiOrderCoverage` of a real 1°×1° extent expanded to a working level,
#     which is many DISJOINT windows. `cellposition` is a binary search over
#     them, membership likewise, and nothing can be recognised as a subtree.
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

# The DEM tile `hydrology.jl` works over; each system covers it at its own
# working level, chosen so the sets stay in the hundreds-to-thousands. The
# shapes, not the scale, are what changes between systems.
tile = Extents.Extent(X=(10.0, 11.0), Y=(46.0, 47.0))

"""
A single dome with a small ripple on it, sampled at cell centroids. The dome
gives the flow routers one real drainage tree instead of a plateau of ties; the
ripple keeps the D8 tie-break off the degenerate path where every neighbour is
equidistant in elevation.
"""
function make_dem(sys, root, cells)
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

nwindows(cells) = DGG.Fallbacks.nwindows(DGG.Fallbacks.windows(cells))

# --- flow routing, per backend ----------------------------------------------

# The full-featured backend: all three routers run, and every direction they
# persist is a relative-cell code.
function flow_battery(sys::DGG.IGeo7System, dem, cells, complete)
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
end

# The generic backend: D8 only, and its persisted direction is the downstream
# neighbour's slot in the cell's complete-level ring, as one set bit.
# (`height_above_nearest_drainage` is also IGeo7-only upstream, but its
# refusal is a bare `MethodError` on cell subtraction, not a contract worth
# pinning — so it is skipped here, not asserted.)
function flow_battery(sys, dem, cells, complete)
    @testset "flow routing: D8" begin
        accumulation, directions = GM.flowaccumulation(dem; method=GM.D8())
        @test size(accumulation) == size(directions) == size(dem)
        @test all(isfinite, parent(accumulation))
        # Accumulation is seeded with each cell's OWN area before donors are
        # added, so no cell may hold less. Per-cell rather than first-cell,
        # because H3 areas vary a couple of per cent within the tile.
        @test all(1:length(cells)) do i
            accumulation[i] >= 0.99 * GM.cellarea(dem, cells[i])
        end

        # The slot codec closing, decoded here independently from the format's
        # promise (slot `k` is bit `k-1`): a settled cell stores exactly one
        # slot, the slot names a ring member that is IN the subset, and that
        # member — holding its donor's water — accumulates strictly more. An
        # off-by-one in the encoder fails membership at the rim and
        # monotonicity in the interior, where every ring member is a member of
        # the subset and membership alone would not see it.
        nonzero = 0
        bad = 0
        for (i, c) in enumerate(cells)
            bits = Int(directions[i])
            iszero(bits) && continue
            nonzero += 1
            ring = DGG.neighbors(complete, c)
            slot = trailing_zeros(bits) + 1
            if !ispow2(bits) || slot > length(ring)
                bad += 1
                continue
            end
            p = DGG.cellposition(cells, ring[slot])
            (p === nothing || !(accumulation[i] < accumulation[p])) && (bad += 1)
        end
        @test nonzero > 0               # the dome drains; all-zero is vacuous
        @test bad == 0
    end

    # The refusal is part of the surface: it must name the backend that could
    # have answered, or a user on the wrong system gets a wrong-shape error
    # from deep inside a router instead of an instruction.
    @testset "DInf and FD8 refuse, and say why" begin
        for method in (GM.DInf(), GM.FD8())
            @test_throws "only the IGeo7 backend provides" GM.flowaccumulation(dem; method)
        end
    end
end

function battery(sys, root, label, cells)
    dem = make_dem(sys, root, cells)
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
            # so it cannot catch an unprintable error. The fragment is the
            # cell type's own print name, whatever the system calls it.
            @test_throws ("at index [" * String(first(split(repr(outside), '(')))) dem[outside]
        end

        # Neither the axis nor the ring may put anything on the heap. The axis is
        # compressed windows, so iterating it computes ids rather than reading
        # them; the ring is clipped into the same static container the system
        # answered in. Both are per-cell costs in every loop Geomorphometry runs,
        # and a `Vector` in either place is invisible to every other assertion
        # here — it returns the right cells, just not for free. Measured on
        # IGeo7 only: the machinery under measurement is system-generic, so one
        # backend's zero is every backend's zero.
        if sys isa DGG.IGeo7System
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

        flow_battery(sys, dem, cells, complete)

        @test !isempty(GM.outlets(dem))
    end
    return dem
end

@testset "Geomorphometry on synthetic DGGS fields" begin
    @testset "IGeo7System" begin
        sys = DGG.IGeo7System()
        coverage = DGG.query(sys, DGG.MultiOrderCoverage(tile); level=10)
        root = DGG.coarsest_contained(coverage)
        subtree = DGG.CellVector(DGG.PartialGrid(sys, root, DGG.level(root) + 3))
        scattered = DGG.CellVector(coverage; level=10)

        # The two shapes really are different underneath; without this the
        # second case could silently become a rerun of the first.
        @test nwindows(subtree) == 1
        @test nwindows(scattered) > 100
        @test DGG.level(root) == 5

        battery(sys, root, "one rooted subtree ($(length(subtree)) cells, 1 window)",
            subtree)
        battery(sys, root,
            "multi-order coverage ($(length(scattered)) cells, $(nwindows(scattered)) windows)",
            scattered)

        # The subset's rim is where a clipped ring makes a local metric answer
        # from fewer samples than the interior. This is not a defect of the
        # walk, it is what `halo` exists to supply — pinned on the small case
        # so that a change in either direction is visible.
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

    # The generic backends, each over the same tile at its own working level
    # (coverage level, subtree depth below the coverage's coarsest root).
    @testset "$(nameof(typeof(sys)))" for (sys, covlvl, depth) in (
        (DGG.HEALPixSystem(), 11, 4),
        (DGG.ISEA4RSystem(), 11, 4),
        (DGG.H3System(), 7, 3),
    )
        coverage = DGG.query(sys, DGG.MultiOrderCoverage(tile); level=covlvl)
        root = DGG.coarsest_contained(coverage)
        subtree = DGG.CellVector(DGG.PartialGrid(sys, root, DGG.level(root) + depth))
        scattered = DGG.CellVector(coverage; level=covlvl)

        @test nwindows(subtree) == 1
        @test nwindows(scattered) > 1

        battery(sys, root, "one rooted subtree ($(length(subtree)) cells, 1 window)",
            subtree)
        battery(sys, root,
            "multi-order coverage ($(length(scattered)) cells, $(nwindows(scattered)) windows)",
            scattered)
    end
end
