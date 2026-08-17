#   julia --project=<integration-env> test/integration/geomorphometry_synthetic.jl
#
# Run the synthetic Geomorphometry integration battery on IGeo7, HEALPix,
# ISEA4R, and H3. Each system is tested with a rooted subtree and a multi-window
# coverage. IGeo7 exercises every router, HAND, and relative-cell directions;
# the other systems exercise D8 with ring-slot directions and reject
# operations that require relative-cell arithmetic. No system supports the
# stream-mask HAND entry.
#
# The environment needs this package, `Rasters`, and `Geomorphometry` at the
# `clipped-neighbors` rev of https://github.com/asinghvi17/GeoArrayOps.jl.git.

using Test
import DiscreteGlobalGrids as DGG
import Geomorphometry as GM
using Rasters

const Extents = Rasters.Extents          # Available through Rasters.

# Each system covers the same Alpine tile at a suitable working level.
tile = Extents.Extent(X=(10.0, 11.0), Y=(46.0, 47.0))

"""
Return a cell-axis raster sampled from a dome with a small ripple. The field has
a non-flat drainage surface and avoids uniformly tied neighbours.
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

# IGeo7 supports all routers and stores relative-cell directions.
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

        # Every decoded direction must stay at the source or reach a neighbour.
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

    # Even IGeo7 refuses the stream-mask entry: DGGS rasters support only the
    # threshold form exercised above.
    @test_throws "use the threshold form" GM.height_above_nearest_drainage(
        dem, falses(size(dem)))
end

# Other systems support D8 and store one set bit for the downstream ring slot
# (slot `k` as bit `k - 1` of a `UInt16`; the checks below decode it
# independently of the extension). HAND needs IGeo7 relative-cell arithmetic,
# so here it refuses; asserted with DInf/FD8 below.
function flow_battery(sys, dem, cells, complete)
    @testset "flow routing: D8" begin
        accumulation, directions = GM.flowaccumulation(dem; method=GM.D8())
        @test size(accumulation) == size(directions) == size(dem)
        @test all(isfinite, parent(accumulation))
        # Each cell starts with its own area, which may vary across the tile.
        @test all(1:length(cells)) do i
            accumulation[i] >= 0.99 * GM.cellarea(dem, cells[i])
        end

        # A nonzero direction names one in-subset ring slot whose accumulation
        # is greater than the donor's.
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
        @test nonzero > 0               # Require at least one routed cell.
        @test bad == 0
    end

    # Unsupported operations name the backend that could answer; without the
    # guard, HAND would die on cell subtraction with a bare `MethodError`.
    @testset "DInf, FD8 and HAND refuse, and say why" begin
        for method in (GM.DInf(), GM.FD8())
            @test_throws "only the IGeo7 backend provides" GM.flowaccumulation(dem; method)
        end
        @test_throws "needs relative-cell arithmetic" GM.height_above_nearest_drainage(dem)
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

            # Subset neighbours are the complete ring clipped to membership.
            clipped = GM.neighbors(dem, c)
            @test !isempty(clipped)
            @test all(in(indices), clipped)
            @test clipped ⊆ DGG.neighbors(complete, c)

            outside = DGG.cellindex(complete, DGG.ncells(complete))
            @test outside ∉ indices
            # String form on purpose: the type form never renders the message,
            # so it cannot catch an unprintable cell in the error.
            @test_throws ("at index [" * String(first(split(repr(outside), '(')))) dem[outside]
        end

        # The cell-axis and clipped-ring loops allocate nothing. IGeo7 only:
        # the machinery is system-generic, so one backend's zero is every
        # backend's zero.
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

        # Geometry operations use an in-subset neighbour.
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

        # Exercise both one-window and multi-window storage.
        @test nwindows(subtree) == 1
        @test nwindows(scattered) > 100
        @test DGG.level(root) == 5

        battery(sys, root, "one rooted subtree ($(length(subtree)) cells, 1 window)",
            subtree)
        battery(sys, root,
            "multi-order coverage ($(length(scattered)) cells, $(nwindows(scattered)) windows)",
            scattered)

        # The halo contains the neighbours omitted from clipped rim rings.
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

    # Coverage level and subtree depth keep each system's test sets comparable.
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
