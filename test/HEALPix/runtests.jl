module HealpixTestSuite

# Tests for the HealpixLookups submodule of DiscreteGlobalGrids, plus the
# suites included at the bottom: the HEALPix operations-kernel wiring
# (test_healpix_kernel.jl), the closed-form face charts and RING/NESTED index
# maps (test_chart.jl), and the dense face grids as spatial trees
# (test_face_grid.jl).
using Test
using LinearAlgebra: dot
using Random
using Statistics
using ConservativeRegridding
import Healpix
import DimensionalData as DD
using DimensionalData: Lookups
import GeometryOps as GO
import GeoInterface as GI
import Extents
import SparseArrays
using SmallCollections: SmallVector

using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
using DiscreteGlobalGrids.HEALPix.HealpixLookups
using DiscreteGlobalGrids.HEALPix.HealpixLookups: nested_neighbors, HealpixLookups

@testset "nested_neighbors" begin
    res = Healpix.Resolution(4)            # nside=4, 192 pixels, ids 0:191
    npix = 192
    # Brute-force ground truth: same-level pixels are adjacent iff they share
    # >= 1 boundary corner (edge neighbors share 2, diagonal neighbors share 1).
    corners = map(0:npix-1) do p
        ring = Healpix.nest2ring(res, p + 1)
        b = Healpix.boundariesRing(res, ring, 1, Float64)
        [b[i, :] for i in 1:4]
    end
    sharescorner(p, q) = any(dot(a, b) > 1 - 1e-12 for a in corners[p+1], b in corners[q+1])
    sharededge(p, q) = count(dot(a, b) > 1 - 1e-12 for a in corners[p+1], b in corners[q+1]) >= 2
    for p in 0:npix-1
        nbs = filter(>=(0), collect(nested_neighbors(res, p)))
        @test allunique(nbs)
        @test all(q -> sharescorner(p, q), nbs)                       # no false neighbors
        for q in 0:npix-1
            q == p && continue
            if sharededge(p, q)
                @test q in nbs                                        # no missed edge-neighbors
            end
        end
        for q in nbs                                                  # symmetry
            @test p in filter(>=(0), collect(nested_neighbors(res, q)))
        end
    end
    # HEALPix invariant: exactly 24 pixels have 7 neighbors — the base tiling has 8
    # degree-3 vertices (lat ±41.8°), and each of the 3 pixels meeting at one misses its
    # diagonal across it (E/W corners of the 8 polar faces + N/S corners of the 4
    # equatorial faces = 24; matches the -1 rows of NB_FACEARRAY).
    @test count(p -> any(<(0), nested_neighbors(res, p)), 0:npix-1) == 24
    res8 = Healpix.Resolution(8)
    @test count(p -> any(<(0), nested_neighbors(res8, p)), 0:(12*64 - 1)) == 24
end

# The kernel wiring over `nested_neighbors`: same neighbor sets, re-stated in
# the kernel's contract — ascending ids, existing neighbors only, in the
# `SmallVector` container `max_neighbors` sizes (src/HEALPix/HealpixKernel.jl).
@testset "cell_neighbors kernel wiring" begin
    level = 2                                     # nside=4, the brute-forced grid above
    res = Healpix.Resolution(2^level)
    npix = 12 * 4^level
    @test max_neighbors(HEALPixDGGS()) == 8
    sevens = 0
    for p in 0:npix-1
        nbs = cell_neighbors(HEALPixDGGS(), level, p)
        @test nbs isa SmallVector{8,Int64}
        @test issorted(nbs) && allunique(nbs)
        @test !(p in nbs)
        @test sort(collect(nbs)) == sort(filter(>=(0), collect(nested_neighbors(res, p))))
        length(nbs) == 7 && (sevens += 1)
        # symmetry, which the compass tuple never states
        for q in nbs
            @test p in cell_neighbors(HEALPixDGGS(), level, q)
        end
    end
    @test sevens == 24
    @test_throws ArgumentError cell_neighbors(HEALPixDGGS(), level, -1)
    @test_throws ArgumentError cell_neighbors(HEALPixDGGS(), level, npix)
end

@testset "HealpixLookup basics" begin
    level = 3                       # nside=8, 768 cells
    ids = Int64[5, 17, 18, 19, 100, 101, 102, 103, 700]
    l = HealpixLookup(ids; level)
    @test l isa Lookups.Lookup{Int64,1}
    @test length(l) == 9
    @test Lookups.order(l) isa Lookups.ForwardOrdered
    @test l.level == 3 && l.nside == 8
    @test l.metadata["grid_name"] == "healpix"
    @test l.metadata["indexing_scheme"] == "nested"
    @test_throws ArgumentError HealpixLookup([3, 2, 1]; level)        # unsorted
    @test_throws ArgumentError HealpixLookup([1, 1]; level)           # duplicate
    @test_throws ArgumentError HealpixLookup([1, 1000]; level)        # out of range
    @test_throws DomainError HealpixLookup(Int64[]; level=-1)

    source_ids = Int64[5, 17]
    copied = HealpixLookup(source_ids; level)
    source_ids[1] = 17
    @test parent(copied) == Int64[5, 17]

    A = DD.DimArray(Float64.(1:9), Cells(l))
    @test A[Cells(Lookups.At(100))] == 5.0
    @test_throws ArgumentError A[Cells(Lookups.At(6))]                # not stored
    # point selector: center of stored cell 17 maps back to cell 17
    lon, lat = cell_centers(l)[2]
    @test A[Cells(Lookups.Contains((lon, lat)))] == 2.0
    # slicing keeps a working (sorted) lookup
    B = A[2:4]
    @test parent(DD.lookup(B, Cells)) == Int64[17, 18, 19]
    # reverse slicing would break the sorted invariant -> must throw, not corrupt
    @test_throws ArgumentError A[9:-1:1]
    @test_throws ArgumentError DD.rebuild(l; data=Int64[5, 5])
    @test_throws ArgumentError DD.rebuild(l; data=Int64[-1])
    @test_throws ArgumentError DD.rebuild(l; data=Int64[1000])
    # The checks live in the inner constructor, so the positional form runs
    # them too — it used to be the *default* constructor, open beside the
    # keyword one that validated. And `nside` is derived from `level`: there is
    # no constructor that takes one, so it cannot disagree with the ids' grid.
    @test HealpixLookup(Int64[5, 17], 3, Dict{String,Any}()).nside == 8
    @test_throws ArgumentError HealpixLookup(Int64[3, 2], 3, Dict{String,Any}())
    @test_throws ArgumentError HealpixLookup(Int64[1000], 3, Dict{String,Any}())
    @test_throws DomainError HealpixLookup(Int64[], -1, Dict{String,Any}())
    @test_throws MethodError HealpixLookup(Int64[5], 3, 4, Dict{String,Any}())
    # display must not throw (requires show_properties methods)
    @test sprint(show, MIME"text/plain"(), A) isa String
end

# A `HealpixLookup` is regriddable straight out of the box — `treeify` routes it
# through `DGGSPartialGrid` (src/HEALPix/HealpixKernel.jl). The tree layer
# itself is covered by test_healpix_kernel.jl and test/core/test_generic_trees.jl.
@testset "Healpix lookup regridding" begin
    level = 2
    ids = Int64[0, 1, 2, 3, 16, 17, 18, 19]
    l = HealpixLookup(ids; level)

    lookup_regridder = ConservativeRegridding.Regridder(l, l; threaded=false, normalize=false)
    @test size(lookup_regridder.intersections) == (length(ids), length(ids))
    @test SparseArrays.nnz(lookup_regridder.intersections) >= length(ids)
    @test isapprox(sum(lookup_regridder.intersections), sum(lookup_regridder.dst_areas); rtol=1e-10)
end

@testset "bounded high-level query descent" begin
    # A single stored pixel at level 18: the descent walks an 18-level chain
    # of one-child nodes without ever building a subtree outline for it —
    # the depth used to be what forced the planar descent's densification
    # cutoff, and now just has to stay cheap and correct.
    level = 18
    res = Healpix.Resolution(2^level)
    pixel = Int64(
        Healpix.ang2pixNest(res, deg2rad(90 - 20.0), deg2rad(30.0)) - 1,
    )
    lookup = HealpixLookup([pixel]; level)
    lon, lat = only(cell_centers(lookup))
    δ = 1e-5
    box = GI.Polygon([
        GI.LinearRing([
            (lon - δ, lat - δ),
            (lon + δ, lat - δ),
            (lon + δ, lat + δ),
            (lon - δ, lat + δ),
            (lon - δ, lat - δ),
        ]),
    ])
    @test Lookups.selectindices(lookup, Lookups.Contains(box)) == [1]
end

# The queries are spherical now, so the oracle is too: unprepared brute
# force over every stored cell's unit-sphere geometry, through the same
# spherical RelateNG engine but none of the tree, cap or descent machinery.
const SPHERICAL_ALG = GO.RelateNG(; manifold=GO.Spherical())
function sph_truth(l, geom, mode::Symbol)
    prep = GO.prepare(SPHERICAL_ALG, geom)
    mode === :center && return findall(i -> GO.relate_predicate(prep, GO.pred_contains(),
        DGG.cell_center(HEALPixDGGS(), l.level, l.data[i])), eachindex(l.data))
    return findall(i -> GO.relate_predicate(prep, GO.pred_intersects(),
        cell_polygon_unitsphere(HEALPixDGGS(), l.level, l.data[i])), eachindex(l.data))
end

@testset "polygon cover" begin
    level = 5                                    # nside=32, 12288 cells
    # coverage: cells with centers in lon 0..40, lat 30..60
    res = Healpix.Resolution(2^level)
    allcenters = [HealpixLookups._cell_center_lonlat(res, p) for p in 0:(12*4^level - 1)]
    covered = Int64[p for p in 0:(12*4^level - 1) if
        0 <= allcenters[p+1][1] <= 40 && 30 <= allcenters[p+1][2] <= 60]
    l = HealpixLookup(covered; level)
    box = GI.Polygon([GI.LinearRing([(5.0, 40.0), (20.0, 40.0), (20.0, 52.0), (5.0, 52.0), (5.0, 40.0)])])
    idx = Lookups.selectindices(l, Lookups.Contains(box))
    @test sort(idx) == sph_truth(l, box, :center)
    @test !isempty(idx)
    # Touching against its own brute-force ground truth; the subset check
    # alone would pass even if Touching returned every stored cell.
    idxt = Lookups.selectindices(l, Touching(box))
    @test sort(idxt) == sph_truth(l, box, :touches)
    @test issubset(sort(idx), sort(idxt))
    # DD-native Touches with an Extent means the lon/lat-aligned box REGION —
    # parallels top and bottom — which is exactly the densified polygon
    # `_extent_polygon` builds. (NOT the same cells as `Touching(box)`: the
    # plain 4-corner box has great-circle edges that bulge off the parallels.)
    ext = Extents.Extent(X=(5.0, 20.0), Y=(40.0, 52.0))
    @test sort(Lookups.selectindices(l, Lookups.Touches(ext))) ==
          sph_truth(l, HealpixLookups._extent_polygon(ext), :touches)
    @test_throws ArgumentError HealpixLookups._extent_polygon(
        Extents.Extent(X=(-180.0, 180.0), Y=(-10.0, 10.0)))
    # query box entirely outside coverage -> empty
    far = GI.Polygon([GI.LinearRing([(-100.0, -40.0), (-90.0, -40.0), (-90.0, -30.0), (-100.0, -30.0), (-100.0, -40.0)])])
    @test isempty(Lookups.selectindices(l, Lookups.Contains(far)))
    # Fuzz: fixed boxes alone DO NOT catch non-conservative coarse-node
    # pruning (the plan review of the planar descent measured ~4% of random
    # boxes affected before its densification fix). Seeded random boxes,
    # descent vs spherical brute force, must agree exactly.
    rng = Random.MersenneTwister(42)
    for _ in 1:200
        lon0 = rand(rng) * 50 - 5; lat0 = rand(rng) * 30 + 28
        w = rand(rng) * 20 + 0.5; h = rand(rng) * 15 + 0.5
        b = GI.Polygon([GI.LinearRing([(lon0, lat0), (lon0 + w, lat0), (lon0 + w, lat0 + h), (lon0, lat0 + h), (lon0, lat0)])])
        @test sort(Lookups.selectindices(l, Lookups.Contains(b))) == sph_truth(l, b, :center)
        @test sort(Lookups.selectindices(l, Touching(b))) == sph_truth(l, b, :touches)
    end
end

@testset "zonal" begin
    level = 5
    res = Healpix.Resolution(2^level)
    covered = Int64[p for p in 0:(12*4^level - 1) if
        let (lon, lat) = HealpixLookups._cell_center_lonlat(res, p)
            -15 <= lon <= 45 && 30 <= lat <= 65
        end]
    l = HealpixLookup(covered; level)
    field(lon, lat) = 20 - 0.5 * (lat - 30) + 2 * sind(3 * lon)
    A = DD.DimArray([field(lon, lat) for (lon, lat) in cell_centers(l)], Cells(l))
    box = GI.Polygon([GI.LinearRing([(5.0, 40.0), (20.0, 40.0), (20.0, 52.0), (5.0, 52.0), (5.0, 40.0)])])
    zs = zonal(mean, A; of=[box])
    truthidx = sph_truth(l, box, :center)
    @test zs[1] ≈ mean(parent(A)[truthidx])
    # multiple zones incl. one fully outside coverage
    far = GI.Polygon([GI.LinearRing([(-100.0, -40.0), (-90.0, -40.0), (-90.0, -30.0), (-100.0, -30.0), (-100.0, -40.0)])])
    zs2 = zonal(mean, A; of=[box, far])
    @test zs2[1] == zs[1] && ismissing(zs2[2])
    # extra dimensions survive: f receives the whole Cells x Ti selection
    A2 = DD.DimArray([field(lon, lat) + t for (lon, lat) in cell_centers(l), t in (0.0, 10.0)],
                     (Cells(l), DD.Ti([0.0, 10.0])))
    zs3 = zonal(mean, A2; of=[box])
    @test zs3[1] ≈ zs[1] + 5.0
end

# `zonal` is the package-level generic, and so is the spherical tree query
# under it — HEALPix's old planar `_query_indices` quadtree is gone, so
# there is no second descent to cross-check anymore. What replaced the
# old cross-check: the seeded-fuzz spherical-brute-force equivalence in
# "polygon cover" above, and the all-systems (antimeridian/pole included)
# equivalence suite in test/core/test_tree_queries.jl. The accessors the
# generic query reads stay pinned here:
@testset "generic query accessors" begin
    l = HealpixLookup(Int64[3, 5]; level=2)
    @test dggs_system(l) === HEALPixDGGS()
    @test dggs_level(l) == 2
end

@testset "stencil" begin
    level = 5
    res = Healpix.Resolution(2^level)
    covered = Int64[p for p in 0:(12*4^level - 1) if
        let (lon, lat) = HealpixLookups._cell_center_lonlat(res, p)
            -15 <= lon <= 45 && 30 <= lat <= 65
        end]
    l = HealpixLookup(covered; level)
    # constant field: any convex-combination stencil returns the constant on interior cells
    A = DD.DimArray(fill(7.0, length(l)), Cells(l))
    sm = stencil(A) do center, nbs
        mean(vcat(center, nbs))
    end
    @test all(==(7.0), parent(sm))
    # the halo table is the generic one now: `SmallVector` positions in
    # ascending-neighbor-id order, `0` only where coverage ends (nonexistent
    # neighbors are dropped, not sentinelled)
    nbi = HealpixLookups.neighbor_indices(l)
    @test eltype(nbi) == SmallVector{8,Int}
    @test all(idx -> 7 <= length(idx) <= 8, nbi)
    @test_throws DimensionMismatch stencil((c, nbs) -> c + sum(nbs), A; nbidx=nbi[1:end-1])
    interior = [i for i in eachindex(nbi) if all(>(0), nbi[i])]
    @test !isempty(interior)                       # region interior has full 8-neighborhoods
    @test any(idx -> any(==(0), idx), nbi)         # coverage boundary has missing neighbors
    # smoothing pulls a bumpy field toward its local mean: variance must drop
    B = DD.DimArray([sind(40 * lon) * cosd(40 * lat) for (lon, lat) in cell_centers(l)], Cells(l))
    smB = stencil((c, nbs) -> mean(vcat(c, nbs)), B)
    @test var(parent(smB)) < var(parent(B))
    # Laplacian of a constant is zero
    lapA = stencil((c, nbs) -> isempty(nbs) ? 0.0 : mean(nbs) - c, A)
    @test all(==(0.0), parent(lapA))
end

# The complement of this module's partial-coverage premise: the globe-complete
# dimension, an ordinary `HealpixLookup` over computed rather than stored ids
# (docs/design/full_globe_lookups.md §1.2, §1.4, §1.5).
@testset "HEALPix globe lookup" begin
    # Level 29 is 3.5e18 cells — 27 exabytes of `Int64` — so every assertion
    # below doubles as an assertion that nothing on its path materialized: a
    # leak is an `OutOfMemoryError`, not a slow test.
    lookup = HealpixLookup(DGGSGlobeIds(HEALPixDGGS(), 29))
    @test lookup.data isa DGGSGlobeIds
    @test lookup.level == 29                       # defaulted from the ids
    @test lookup.nside == 2^29                     # ...and derived from that
    @test length(lookup) == DiscreteGlobalGrids.num_cells(HEALPixDGGS(), 29)
    @test lookup[1] == 0
    @test lookup[end] == length(lookup) - 1        # ordinal ids: 0:npix-1
    @test lookup.metadata["grid_name"] == "healpix"
    @test lookup.metadata["level"] == 29
    # `show` prints the length and the description, never the ids.
    text = sprint(show, MIME"text/plain"(), lookup)
    @test occursin("level: 29", text)
    @test occursin("ncells: $(length(lookup))", text)
    @test occursin("DGGSGlobeIds", text)
    # Another system's globe is rejected outright rather than falling through
    # to the generic constructor, which would `collect` it to find out.
    @test_throws ArgumentError HealpixLookup(DGGSGlobeIds(H3DGGS(), 3))
    # A declared level disagreeing with the ids' own still fails the O(1)
    # `[0, npix)` endpoint check, as it does for an explicit id vector.
    @test_throws ArgumentError HealpixLookup(DGGSGlobeIds(HEALPixDGGS(), 5); level=3)

    # Laziness survives every rebuild — the invariant most likely to regress
    # silently, since `DD.rebuild` defaults `data=l.data` and routes through the
    # keyword constructor: without its `DGGSGlobeIds` method a `format`, a `set`
    # or a metadata change would quietly materialize the globe, and the failure
    # would be a working answer rather than an error.
    @test DD.rebuild(lookup).data isa DGGSGlobeIds
    relabelled = DD.rebuild(lookup; metadata=Dict{String,Any}("title" => "globe"))
    @test relabelled.data isa DGGSGlobeIds
    @test relabelled.metadata["title"] == "globe"
    dim = Cells(lookup)
    @test DD.val(DD.format(dim, Base.OneTo(length(lookup)))).data isa DGGSGlobeIds
    @test DD.val(DD.set(dim, Cells)).data isa DGGSGlobeIds

    # Degradation needs no new code: a non-scalar index falls through to
    # `AbstractArray`'s default, which materializes a `Vector`, so the rebuild
    # behind it takes the ordinary validating path and lands back on this
    # module's partial-coverage premise. The discriminator is the type of
    # `data`, never "was this a rebuild".
    partial = lookup[1:10]
    @test partial isa HealpixLookup
    @test partial.data isa Vector{Int64}
    @test partial.level == 29
    @test collect(partial) == collect(Int64, 0:9)
    @test Lookups.selectindices(partial, Lookups.At(4)) == 5
    @test_throws ArgumentError Lookups.selectindices(partial, Lookups.At(10))

    # A globe small enough to hold is where the construction can be compared
    # against the explicit one it replaces, dimension and selectors included.
    globe = HealpixLookup(DGGSGlobeIds(HEALPixDGGS(), 1))
    @test collect(globe) == collect(HealpixLookup(collect(Int64, 0:47); level=1))
    array = DD.DimArray(collect(1:length(globe)), Cells(globe))
    @test DD.lookup(array, Cells).data isa DGGSGlobeIds     # survived `format`
    @test DD.lookup(DD.set(array, Cells => Cells), Cells).data isa DGGSGlobeIds
    @test array[Cells(Lookups.At(6))] == 7                  # 0-based ids
    lon, lat = cell_centers(globe)[13]
    @test array[Cells(Lookups.Contains((lon, lat)))] == 13
    sliced = array[Cells(1:5)]
    @test DD.lookup(sliced, Cells).data isa Vector{Int64}
    @test collect(DD.lookup(sliced, Cells)) == collect(Int64, 0:4)
    @test sliced[Cells(Lookups.At(2))] == 3
end

# Operations-kernel wiring of `HEALPixDGGS` (src/HEALPix/HealpixKernel.jl);
# the file wraps itself in a module of its own.
include("test_healpix_kernel.jl")

# Pure closed-form face charts + RING/NESTED index maps (src/HEALPix/chart.jl);
# also its own module.
include("test_chart.jl")

# Dense per-resolution face grids built on the charts, as spatial trees
# (src/HEALPix/face_grid.jl); also its own module.
include("test_face_grid.jl")

println("All testsets finished.")

end # module HealpixTestSuite
