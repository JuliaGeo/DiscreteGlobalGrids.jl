# test/H3/test_neighbors.jl — edge neighbors for H3: the `grid_disk` /
# `max_grid_disk_size` native wrappers (src/H3/H3Native.jl), the kernel
# wiring (`DGG.cell_neighbors(::H3DGGS, ...)`, src/H3/H3Kernel.jl), and the
# lookup operations built on it (`neighbor_indices`, `stencil`, `zonal` —
# src/core/lookup_ops.jl).
#
# Ground truth for the wiring is libh3 itself (`gridDisk` k=1 is the
# canonical pentagon-safe neighbor enumeration), so beyond set equality with
# the raw disk the sweeps check the properties the *kernel contract* adds and
# libh3 does not promise: ascending order, symmetry, counts (6 per hexagon,
# 5 per pentagon), and a geometric sanity bound tying every claimed neighbor
# to its cell's own scale.

module H3NeighborTests

using Test
using Random
using Statistics
import DimensionalData as DD
import GeoInterface as GI
import GeometryOps as GO
using SmallCollections: SmallVector

using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
using DiscreteGlobalGrids.H3.H3Lookups: H3Lookup, H3Cells

const H3N = DGG.H3.H3Native
const S = H3DGGS()

"All cells of `res`, ascending (root subtrees enumerate in ascending id order)."
level_cells(res) = sort!(reduce(vcat, [DGG.cell_descendants(S, 0, root, res)
                                       for root in DGG.root_ids(S)]))

@testset "H3 neighbors" begin
    @testset "grid_disk native wrappers" begin
        @test H3N.max_grid_disk_size(0) == 1
        @test H3N.max_grid_disk_size(1) == 7
        @test H3N.max_grid_disk_size(3) == 3 * 3 * 4 + 1
        @test_throws ArgumentError H3N.max_grid_disk_size(-1)

        cell = H3N.lonlat_to_cell(-122.41795063018799, 37.775938728915946, 5)
        @test H3N.grid_disk(cell, 0) == [cell]
        disk = H3N.grid_disk(cell, 1)
        @test length(disk) == 7
        @test cell in disk
        @test count(!=(UInt64(0)), disk) == 7        # no pentagon nearby: full disk
        @test all(id -> H3N.get_resolution(id) == 5, filter(!=(UInt64(0)), disk))
        # the pentagon-truncated disk keeps its zero padding
        pentagon = first(H3N.get_pentagons(5))
        pdisk = H3N.grid_disk(pentagon, 1)
        @test count(!=(UInt64(0)), pdisk) == 6       # pentagon + its 5 neighbors
    end

    @testset "kernel wiring against raw gridDisk, res 0:2" begin
        @test DGG.max_neighbors(S) == 6
        for res in 0:2
            cells = level_cells(res)
            @test length(cells) == H3N.num_cells(res)
            nbmap = Dict(id => DGG.cell_neighbors(S, res, id) for id in cells)
            failures = 0
            for id in cells
                nbs = nbmap[id]
                nbs isa SmallVector{6,UInt64} || (failures += 1)
                issorted(nbs) || (failures += 1)
                allunique(nbs) || (failures += 1)
                id in nbs && (failures += 1)
                length(nbs) == (H3N.is_pentagon(id) ? 5 : 6) || (failures += 1)
                # set equality with the disk libh3 answers
                Set(nbs) == Set(c for c in H3N.grid_disk(id, 1)
                                if c != 0 && c != id) || (failures += 1)
                # symmetry, which libh3 promises nowhere on this API
                for nb in nbs
                    id in nbmap[nb] || (failures += 1)
                end
            end
            @test failures == 0
            # 12 pentagons with 5, hexagons with 6: twice the tiling's edge
            # count, E = 3F - 6 (Euler, 3-valent vertices).
            @test sum(length(nbmap[id]) for id in cells) ==
                  2 * (3 * H3N.num_cells(res) - 6)
        end
    end

    @testset "geometric sanity at finer res" begin
        # Around every res-5 pentagon plus a seeded sample: each neighbor's
        # center lies within a small multiple of the cell's own circumradius —
        # a neighbor from a wrong grid distance would sit whole cells away.
        rng = MersenneTwister(11)
        cells = UInt64[]
        append!(cells, H3N.get_pentagons(5))
        for _ in 1:500
            push!(cells, DGG.ordinal_to_cell(S, 5, rand(rng, 1:Int(H3N.num_cells(5)))))
        end
        failures = 0
        for id in cells
            center = DGG.cell_center(S, 5, id)
            radius = maximum(GO.UnitSpherical.spherical_distance(center, p)
                             for p in DGG.cell_boundary(S, 5, id))
            for nb in DGG.cell_neighbors(S, 5, id)
                d = GO.UnitSpherical.spherical_distance(center, DGG.cell_center(S, 5, nb))
                (0 < d < 3 * radius) || (failures += 1)
                id in DGG.cell_neighbors(S, 5, nb) || (failures += 1)
            end
        end
        @test failures == 0
    end

    @testset "neighbor_indices and stencil over an H3Lookup" begin
        res = 4
        origin = H3N.lonlat_to_cell(10.0, 45.0, res)
        ids = sort!(unique(filter(!=(UInt64(0)), H3N.grid_disk(origin, 5))))
        l = H3Lookup(ids; resolution=res, validate=true)
        @test dggs_system(l) === S
        @test dggs_level(l) == res

        nbi = neighbor_indices(l)
        @test eltype(nbi) == SmallVector{6,Int}
        @test length(nbi) == length(ids)
        for (i, id) in enumerate(ids)
            expected = [something(DGG.cell_position(ids, nb), 0)
                        for nb in DGG.cell_neighbors(S, res, id)]
            @test collect(nbi[i]) == expected
        end
        @test any(idx -> all(>(0), idx), nbi)        # disk interior is fully stored
        @test any(idx -> any(==(0), idx), nbi)       # disk rim is a coverage boundary

        A = DD.DimArray(fill(7.0, length(l)), H3Cells(l))
        sm = stencil((c, nbs) -> mean(vcat(c, nbs)), A)
        @test all(==(7.0), parent(sm))
        @test_throws DimensionMismatch stencil((c, nbs) -> c, A; nbidx=nbi[1:end-1])
        B = DD.DimArray([sind(40 * lon) * cosd(40 * lat)
                         for (lon, lat) in (H3N.cell_center(id) for id in ids)],
                        H3Cells(l))
        @test var(parent(stencil((c, nbs) -> mean(vcat(c, nbs)), B; nbidx=nbi))) <
              var(parent(B))
        @test all(==(0.0), parent(stencil((c, nbs) -> isempty(nbs) ? 0.0 : mean(nbs) - c, A)))

        # The globe-complete dimension: every neighbor is stored, and its
        # position resolves through `cell_to_ordinal`, not a binary search.
        globe = H3Lookup(DGGSGlobeIds(H3DGGS(), 0))
        gnbi = neighbor_indices(globe)
        @test length(gnbi) == 122
        @test all(idx -> all(>(0), idx), gnbi)
        Ag = DD.DimArray(fill(1.0, length(globe)), H3Cells(globe))
        @test all(==(1.0), parent(stencil((c, nbs) -> mean(vcat(c, nbs)), Ag; nbidx=gnbi)))
    end

    @testset "zonal over an H3Lookup" begin
        res = 4
        origin = H3N.lonlat_to_cell(10.0, 45.0, res)
        ids = sort!(unique(filter(!=(UInt64(0)), H3N.grid_disk(origin, 6))))
        l = H3Lookup(ids; resolution=res, validate=true)
        centers = [H3N.cell_center(id) for id in ids]
        field(lon, lat) = 20 - 0.5 * lat + 2 * sind(3 * lon)
        A = DD.DimArray([field(lon, lat) for (lon, lat) in centers], H3Cells(l))

        box = GI.Polygon([GI.LinearRing([
            (8.0, 43.5), (12.0, 43.5), (12.0, 46.0), (8.0, 46.0), (8.0, 43.5)])])
        zs = zonal(mean, A; of=[box])
        # spherical oracle: the query engine's own manifold, none of its tree
        prep = GO.prepare(GO.RelateNG(; manifold=GO.Spherical()), box)
        truth = findall(i -> GO.relate_predicate(prep, GO.pred_contains(),
            DGG.cell_center(S, res, ids[i])), eachindex(ids))
        @test !isempty(truth)
        @test zs[1] ≈ mean(parent(A)[truth])
        # `:touches` is a superset of `:center`, and each hit really touches.
        idx_c = DGG._query_positions(l, box, :center)
        @test idx_c == truth
        idx_t = DGG._query_positions(l, box, :touches)
        @test issubset(idx_c, idx_t)
        @test all(i -> GO.relate_predicate(prep, GO.pred_intersects(),
            DGG.cell_polygon_unitsphere(S, res, ids[i])), idx_t)
        far = GI.Polygon([GI.LinearRing([
            (-100.0, -40.0), (-90.0, -40.0), (-90.0, -30.0), (-100.0, -30.0), (-100.0, -40.0)])])
        zs2 = zonal(mean, A; of=[box, far])
        @test zs2[1] == zs[1] && ismissing(zs2[2])
        # extra dimensions survive: f receives the whole cells x Ti selection
        A2 = DD.DimArray([field(lon, lat) + t for (lon, lat) in centers, t in (0.0, 10.0)],
                         (H3Cells(l), DD.Ti([0.0, 10.0])))
        zs3 = zonal(mean, A2; of=[box])
        @test zs3[1] ≈ zs[1] + 5.0
    end
end

end # module H3NeighborTests
