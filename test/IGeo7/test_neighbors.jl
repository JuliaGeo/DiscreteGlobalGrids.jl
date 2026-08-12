# test/IGeo7/test_neighbors.jl — edge neighbors for IGEO7: the native lattice
# implementation (`IGeo7._cell_neighbors`, src/IGeo7/grid.jl), the kernel
# wiring (`DGG.cell_neighbors(::IGEO7DGGS, ...)`, src/IGeo7/IGeo7Kernel.jl),
# and the lookup operations built on it (`neighbor_indices`, `stencil`,
# `zonal` — src/core/lookup_ops.jl).
#
# There is no recorded neighbor oracle in `vectors/`, so ground truth is
# *geometric adjacency* from the module's own oracle-validated boundaries:
# two same-resolution cells share an edge iff their `cell_boundary_cartesian`
# rings share two corners (the corners come from exact lattice arithmetic, so
# shared ones agree to FP roundoff, far below any corner separation). The
# whole-level sweeps then pin the relation completely: every hexagon has
# exactly 6 edges and every pentagon 5, each edge borders exactly two cells,
# so "right count + every claimed neighbor shares an edge + all distinct" is
# the full edge set, and symmetry cross-checks each pair from both sides.
#
# Big loops accumulate a failure counter and assert once (the idiom of
# test_indexing.jl) so the suite's test count stays readable.

module IGeo7NeighborTests

using Test
using Random
using Statistics
import DimensionalData as DD
import GeoInterface as GI
import GeometryOps as GO
using SmallCollections: SmallVector

using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
using DiscreteGlobalGrids: Helpers, ISEA, IGeo7
using DiscreteGlobalGrids.IGeo7.IGeo7Lookups: IGeo7Lookup

const S = IGEO7DGGS()

"All cells of resolution `res` in canonical (ascending id) order."
level_cells(res) = [IGeo7.index_to_cell(i, res) for i in 1:IGeo7.num_cells(res)]

"Number of shared boundary corners between two same-resolution cells."
function shared_corners(a::UInt64, b::UInt64)
    ca = IGeo7.cell_boundary_cartesian(a; closed_ring=false)
    cb = IGeo7.cell_boundary_cartesian(b; closed_ring=false)
    return count(sum(abs2, p .- q) < 1e-18 for p in ca, q in cb)
end

# One cell's neighbor set against every property the contract names. Returns
# a failure count; `0` means: right count (6 hexagon / 5 pentagon), ascending,
# distinct, no self, symmetric, and every claimed neighbor is a true edge
# neighbor (shares two corners).
function neighbor_failures(id::UInt64, res::Int)
    bad = 0
    nbs = DGG.cell_neighbors(S, res, id)
    length(nbs) == (IGeo7.z7_is_pentagon(id) ? 5 : 6) || (bad += 1)
    issorted(nbs) || (bad += 1)
    allunique(nbs) || (bad += 1)
    id in nbs && (bad += 1)
    for nb in nbs
        IGeo7.get_resolution(nb) == res || (bad += 1)
        id in DGG.cell_neighbors(S, res, nb) || (bad += 1)
        shared_corners(id, nb) >= 2 || (bad += 1)
    end
    return bad
end

@testset "IGeo7 neighbors" begin
    @testset "container and wiring" begin
        @test DGG.max_neighbors(S) == 6
        z = IGeo7.lonlat_to_z7(10.0, 45.0, 4)
        nbs = DGG.cell_neighbors(S, 4, z)
        @test nbs isa SmallVector{6,UInt64}
        # ...and the kernel wiring is the native lattice function re-seated.
        @test collect(nbs) == collect(IGeo7._cell_neighbors(z))
        @test IGeo7._cell_neighbors(z) isa Helpers.SmallList{6,UInt64}

        # Validation rides `_geometry_checked`: invalid ids and res-20 ids
        # (prefix arithmetic only, no geometry) keep the native error type.
        @test_throws IGeo7.InvalidZ7Error DGG.cell_neighbors(S, 0, UInt64(12) << 60)
        res20 = z
        for _ in 1:16                               # res 4 -> res 20, one digit at a time
            res20 = IGeo7.z7_child(res20, 1)
        end
        @test_throws IGeo7.InvalidZ7Error DGG.cell_neighbors(S, 20, res20)
    end

    @testset "res 0: the icosahedron itself" begin
        # The twelve res-0 cells are the vertex pentagons, so each one's
        # neighbors must be exactly its vertex's CCW neighbor ring — an
        # independent table (`ISEA.NBRS_CCW`, pinned by the icosahedron suite).
        for base in 0:11
            root = IGeo7.res0_cells()[base+1]
            nbs = DGG.cell_neighbors(S, 0, root)
            @test length(nbs) == 5
            @test sort(IGeo7.z7_base_cell.(collect(nbs))) ==
                  sort(collect(ISEA.NBRS_CCW[base+1]))
        end
    end

    # Exhaustive whole-level sweeps. Res 3 (3432 cells) already contains every
    # relative position a cell can occupy — pentagon rings, base-region fringe
    # across icosahedron edges, cone-cut crossings — at every base.
    @testset "whole-level properties, res 1:3" begin
        for res in 1:3
            cells = level_cells(res)
            failures = sum(neighbor_failures(id, res) for id in cells)
            @test failures == 0
            # 12 pentagons with 5 neighbors, hexagons with 6: the total is
            # twice the tiling's edge count, E = 3F - 6 (Euler, 3-valent
            # vertices) — a global completeness check independent of any
            # single cell's count.
            total = sum(length(DGG.cell_neighbors(S, res, id)) for id in cells)
            @test total == 2 * (3 * IGeo7.num_cells(res) - 6)
        end
    end

    # The fringe zones the exhaustive sweep only reaches at coarse scale:
    # pentagon 2-rings (the cone apex and its collapse), face centers (where
    # three base regions meet) and icosahedron edge midpoints (where two
    # meet), each at res 6, plus a seeded global sample.
    @testset "targeted fringe and random sample, res 6" begin
        res = 6
        targets = Set{UInt64}()
        for base in 0:11
            pent = IGeo7.cell_to_children(IGeo7.res0_cells()[base+1], res)[1]
            @test IGeo7.z7_is_pentagon(pent)
            push!(targets, pent)
            for nb in DGG.cell_neighbors(S, res, pent)
                push!(targets, nb)
                union!(targets, DGG.cell_neighbors(S, res, nb))
            end
        end
        for f in 0:19
            face = ISEA.FACES[f+1]
            push!(targets, IGeo7._xyz_to_z7(face.c, res))
            v = face.verts
            for (i, j) in ((1, 2), (1, 3), (2, 3))
                mid = ISEA.vnormalize(ISEA.vadd(ISEA.VERTICES[v[i]+1], ISEA.VERTICES[v[j]+1]))
                push!(targets, IGeo7._xyz_to_z7(mid, res))
            end
        end
        rng = MersenneTwister(42)
        for _ in 1:2000
            push!(targets, IGeo7.index_to_cell(rand(rng, 1:IGeo7.num_cells(res)), res))
        end
        failures = sum(neighbor_failures(id, res) for id in targets)
        @test failures == 0
    end

    @testset "neighbor_indices and stencil over an IGeo7Lookup" begin
        res = 4
        # A whole res-1 subtree: interior cells have all 6 neighbors stored,
        # the subtree's fractal rim does not.
        root = IGeo7.z7_from_string("053")
        ids = sort(IGeo7.cell_to_children(root, res))
        l = IGeo7Lookup(ids; resolution=res, validate=true)
        @test dggs_system(l) === S
        @test dggs_level(l) == res

        nbi = neighbor_indices(l)
        @test eltype(nbi) == SmallVector{6,Int}
        @test length(nbi) == length(ids)
        # Positions answer the same question `cell_position` does, with `0`
        # standing in for "not stored".
        for (i, id) in enumerate(ids)
            expected = [something(DGG.cell_position(ids, nb), 0)
                        for nb in DGG.cell_neighbors(S, res, id)]
            @test collect(nbi[i]) == expected
        end
        interior = [i for i in eachindex(nbi) if length(nbi[i]) == 6 && all(>(0), nbi[i])]
        @test !isempty(interior)
        @test any(idx -> any(==(0), idx), nbi)      # the rim is really a rim

        A = DD.DimArray(fill(7.0, length(l)), DD.Dim{:cells}(l))
        sm = stencil((c, nbs) -> mean(vcat(c, nbs)), A)
        @test all(==(7.0), parent(sm))
        @test_throws DimensionMismatch stencil((c, nbs) -> c, A; nbidx=nbi[1:end-1])
        # smoothing pulls a bumpy field toward its local mean
        B = DD.DimArray([sind(40 * lon) * cosd(40 * lat)
                         for (lon, lat) in (IGeo7.cell_center(id) for id in ids)],
                        DD.Dim{:cells}(l))
        @test var(parent(stencil((c, nbs) -> mean(vcat(c, nbs)), B; nbidx=nbi))) <
              var(parent(B))
        # Laplacian of a constant is zero
        @test all(==(0.0), parent(stencil((c, nbs) -> isempty(nbs) ? 0.0 : mean(nbs) - c, A)))

        # The globe-complete dimension: every neighbor is stored, positions
        # resolve through `cell_to_ordinal` rather than a binary search.
        globe = IGeo7Lookup(DGGSGlobeIds(IGEO7DGGS(), 1))
        gnbi = neighbor_indices(globe)
        @test all(idx -> all(>(0), idx), gnbi)
        Ag = DD.DimArray(fill(1.0, length(globe)), DD.Dim{:cells}(globe))
        @test all(==(1.0), parent(stencil((c, nbs) -> mean(vcat(c, nbs)), Ag; nbidx=gnbi)))
    end

    @testset "zonal over an IGeo7Lookup" begin
        res = 4
        ids = sort(IGeo7.cell_to_children(IGeo7.z7_from_string("053"), res))
        l = IGeo7Lookup(ids; resolution=res, validate=true)
        centers = [IGeo7.cell_center(id) for id in ids]
        field(lon, lat) = 20 - 0.5 * lat + 2 * sind(3 * lon)
        A = DD.DimArray([field(lon, lat) for (lon, lat) in centers], DD.Dim{:cells}(l))

        lon0, lat0 = centers[length(centers)÷2]
        box = GI.Polygon([GI.LinearRing([
            (lon0 - 4, lat0 - 3), (lon0 + 4, lat0 - 3),
            (lon0 + 4, lat0 + 3), (lon0 - 4, lat0 + 3), (lon0 - 4, lat0 - 3)])])
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
        # a zone fully outside coverage is `missing`, not an error
        far = GI.Polygon([GI.LinearRing([
            (lon0 + 90, -lat0 - 5), (lon0 + 95, -lat0 - 5),
            (lon0 + 95, -lat0), (lon0 + 90, -lat0), (lon0 + 90, -lat0 - 5)])])
        zs2 = zonal(mean, A; of=[box, far])
        @test zs2[1] == zs[1] && ismissing(zs2[2])
    end
end

end # module IGeo7NeighborTests
