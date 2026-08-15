# The sibling script's question without the sibling's inputs:
#
#   julia --project=<integration-env> test/integration/geomorphometry_synthetic.jl
#
# `geomorphometry_igeo7.jl` needs a Copernicus GeoTIFF, ArchGDAL and a
# conservative regridder, so it can only be run by hand on a machine that has
# all three. Everything it proves about the DGGS side, though, is independent of
# where the elevations came from: what is under test is the neighbour iterator,
# `cellposition` addressing, and the relative-cell direction codec the flow
# routers encode with. So this file synthesises the field instead, and needs
# only an environment holding `Geomorphometry` (the `feat/generic` rev that
# `docs/Project.toml` pins) and `Rasters` alongside this package.
#
# The two are complements, not duplicates: the sibling is the only thing that
# checks the regridder feeds Geomorphometry something sane, and this is the only
# thing that runs anywhere.

using Test
import DiscreteGlobalGrids as DGG
import Geomorphometry as GM
using Rasters

sys = DGG.IGeo7System()
root = DGG.cellindex(DGG.levelgrid(sys, 5), 400)
grid = DGG.PartialGrid(sys, root, 8)
cells = DGG.CellVector(grid)
complete = DGG.levelgrid(sys, DGG.level(cells))

# A single dome with a small ripple on it. The dome gives the flow routers one
# real drainage tree instead of a plateau of ties; the ripple keeps the D8
# tie-break off the degenerate path where every neighbour is equidistant in
# elevation.
apex = DGG.cell_centroid(DGG.levelgrid(sys, 5), root)
elevation = [
    begin
        p = DGG.cell_centroid(complete, c)
        d = acos(clamp(p[1] * apex[1] + p[2] * apex[2] + p[3] * apex[3], -1, 1))
        1000.0 * exp(-(d / 0.02)^2) + 40.0 * sin(120 * p[1]) * cos(97 * p[2])
    end
    for c in cells
]

dem = Raster(elevation, (DGG.Cells(DGG.CellLookup(cells)),); name=:height)

@testset "Geomorphometry on a synthetic IGeo7 field" begin
    @testset "the indexing surface the extension stands on" begin
        indices = eachindex(dem)
        @test indices === cells
        @test indices isa DGG.CellVector
        c = first(indices)
        @test c in indices
        @test dem[c] == elevation[1]

        # `neighbors` on the subset is the system's answer clipped to
        # membership, so it is a SUBSET of the complete level's answer and never
        # a different set.
        clipped = GM.neighbors(dem, c)
        @test !isempty(clipped)
        @test all(in(indices), clipped)
        @test clipped ⊆ DGG.neighbors(complete, c)

        # The relative form, which is the half of the interface `directioncode`
        # is for. Taken against the COMPLETE level, so the pentagon and the
        # subset rim are both in scope rather than clipped away.
        for n in DGG.neighbors(complete, c)
            d = n - c
            @test d isa DGG.RelativeZ7Cell
            @test c + d == n
            @test 1 <= DGG.directioncode(d) <= 6
        end
        @test DGG.directioncode(c - c) == 0

        outside = DGG.cellindex(complete, DGG.ncells(complete))
        @test outside ∉ indices
        @test_throws BoundsError dem[outside]
    end

    @testset "local metrics" begin
        tpi = GM.topographic_position_index(dem)
        @test size(tpi) == size(dem)
        @test all(isfinite, parent(tpi))

        rough = GM.roughness(dem)
        @test size(rough) == size(dem)
        @test all(>=(0), parent(rough))
    end

    @testset "cell geometry" begin
        neighbor = first(DGG.neighbors(complete, first(cells)))
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
        for (i, c) in enumerate(cells)
            iszero(directions[i]) && continue
            for rel in GM.decompose(DGG.RelativeZ7Cell, directions[i], c)
                target = c + rel
                @test target == c || target in DGG.neighbors(complete, c)
            end
        end
    end

    # The subset's rim is where a clipped ring makes a local metric answer from
    # fewer samples than the interior. This is not a defect of the walk, it is
    # what `halo` exists to supply — the count is pinned here so that a change
    # in either direction is visible.
    @testset "the rim a halo would complete" begin
        clipped = count(cells) do c
            length(DGG.neighbors(cells, c)) < length(DGG.neighbors(complete, c))
        end
        @test length(cells) == 343
        @test clipped == 78
        @test length(collect(DGG.halo(cells))) == 84
        @test !isempty(GM.outlets(dem))
    end
end
