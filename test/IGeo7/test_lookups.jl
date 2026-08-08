# test/test_lookups.jl — DimensionalData integration layer
# (src/IGeo7/IGeo7Lookups.jl, design.md Section 9 task 8).
#
# Ported from the untainted test/runtests.jl of the prototype tree's sibling
# module (an allowed source under CLEANROOM.md), adapted to the clean design:
#
#   * the sibling's `root_map = (1, 6, 2, 7, ...)` is gone — with the Z7 u64 as
#     the cell id (design Section 3), `res0_cells()` ascending *is* bases 0:11,
#     so the root map is the identity;
#   * every ordering assertion is internal-consistency (ascending encoded id ==
#     lookup order == dense-index order), never an oracle-order expectation;
#   * the published Z7 example is ported verbatim — it is a paper/contract
#     fact, not an implementation detail;
#   * exhaustive per-cell loops are aggregated into `@test all(...)` because
#     test_z7.jl / test_indexing.jl already own the per-cell granularity; this
#     file tests the *wiring*.
#
# New here (not in the sibling): the metadata defaults, the mixed-resolution
# `validate=true` rejection, and the `Touching` selector. The layer under test
# sees only the public native API.

using Test

# The suite lives in its own module: `IGeo7Lookups` re-exports the native API
# names that `runtests.jl` has already brought into the suite module from
# `IGeo7`, and only a separate namespace can hold both unambiguously.
module TestLookups

using Test
using DimensionalData
using GeoInterface
using GeometryOps

using DiscreteGlobalGrids: Helpers, ISEA
using DiscreteGlobalGrids.IGeo7.IGeo7Lookups
# Qualified rather than `using`-ed: the package namespace exports `cell_polygon`
# and this module already has `IGeo7Lookups`' one in scope.
import DiscreteGlobalGrids as DGG

const DD = DimensionalData
const GI = GeoInterface
const GO = GeometryOps
"the clean native core the layer delegates to"
const CORE = IGeo7Lookups.IGeo7
const SMALL = Helpers

const NYC = (-73.9857, 40.7484)
const LONDON = (-0.1278, 51.5074)
const TOKYO = (139.6917, 35.6895)

# ---------------------------------------------------------------------------

@testset "IGEO7 lookups: native delegation surface" begin
    @test length(res0_cells()) == 12
    @test num_cells(0) == 12
    @test num_cells(1) == 72
    @test num_cells(2) == 492

    roots = res0_cells()
    @test roots isa SMALL.SmallList{12,UInt64}
    @test all(CORE.is_pentagon, roots)
    @test cell_to_children(roots[1], 0) isa SMALL.SmallList{1,UInt64}
    @test cell_to_children(roots[1], 1) isa SMALL.SmallList{7,UInt64}
    @test cell_to_children(roots[1], 2) isa Vector{UInt64}
    @test length(cell_to_children(roots[1], 1)) == 6
    @test length(cell_to_children(roots[1], 2)) == 41
    @test cell_to_children(roots[1]) == cell_to_children(roots[1], 1)
    @test cell_to_children(roots[1], nothing) == cell_to_children(roots[1], 1)

    child = first(cell_to_children(roots[1], 2))
    @test CORE.get_resolution(child) == 2
    @test cell_to_parent(child, 0) == roots[1]
    @test cell_to_parent(child) == cell_to_parent(child, 1)
    @test child in cell_to_children(cell_to_parent(child, 1), 2)

    cell = lonlat_to_cell(0, 0, 1)
    @test CORE.get_resolution(cell) == 1
    @test CORE.is_valid_cell(cell)
    @test cell_area(cell) > 0
    @test length(cell_boundary(cell; closed_ring=false)) in (5, 6)
    @test length(cell_boundary(cell; closed_ring=true)) ==
          length(cell_boundary(cell; closed_ring=false)) + 1
    @test first(cell_boundary(cell)) == last(cell_boundary(cell))

    # Every delegate is exactly the core function, on ids and on points.
    sample = UInt64[roots[1], cell, lonlat_to_cell(LONDON..., 7),
        lonlat_to_cell(TOKYO..., 12)]
    @test all(id -> cell_area(id) == CORE.cell_area(id), sample)
    @test all(id -> cell_center(id) == CORE.cell_center(id), sample)
    @test all(id -> cell_boundary(id) == CORE.cell_boundary(id), sample)
    @test all(id -> cell_to_index(id) == CORE.cell_to_index(id), sample)
    @test all(id -> cell_to_z7(id) == id && z7_to_cell(id) == id, sample)
    @test all(id -> z7_resolution(id) == CORE.z7_resolution(id), sample)
    @test all(id -> z7_base_cell(id) == CORE.z7_base_cell(id), sample)
    @test all(id -> z7_is_pentagon(id) == CORE.z7_is_pentagon(id), sample)
    @test all(id -> z7_to_string(id) == CORE.z7_to_string(id), sample)
    @test all(id -> z7_to_hex(id) == CORE.z7_to_hex(id), sample)
    @test all(id -> z7_from_hex(z7_to_hex(id)) == id, sample)
    @test all(id -> z7_from_string(z7_to_string(id)) == id, sample)
    @test all(id -> z7_digit(id, 1) == CORE.z7_digit(id, 1), sample)
    @test all(id -> z7_children(z7_parent(id)) == CORE.z7_children(CORE.z7_parent(id)),
        sample[2:end])
    @test all(id -> z7_is_descendant(id, z7_parent(id)), sample[2:end])
    @test z7_parent(sample[3], 2) == CORE.z7_parent(sample[3], 2)
    @test z7_child(roots[1], 1) == CORE.z7_child(roots[1], 1)
    @test lonlat_to_cell(NYC..., 7) == CORE.lonlat_to_cell(NYC..., 7)
    @test lonlat_to_z7(NYC..., 7) == CORE.lonlat_to_z7(NYC..., 7)
    @test lonlat_to_index(NYC..., 7) == CORE.lonlat_to_index(NYC..., 7)
    @test index_to_cell(5, 3) == CORE.index_to_cell(5, 3)

    # id coercion: hex strings and other unsigned widths reach the same cell
    @test cell_area("0x" * z7_to_hex(cell)) == cell_area(cell)
    @test cell_center(z7_to_hex(cell)) == cell_center(cell)
    @test cell_to_index(UInt128(cell)) == cell_to_index(cell)

    # the orientation keyword of the clean core survives the delegation
    rot90z = ISEA.Orientation((0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0))
    @test cell_center(roots[1]; orientation=rot90z) !=
          cell_center(roots[1])
    @test lonlat_to_cell(NYC..., 5; orientation=ISEA.ORIENT_IDENTITY) ==
          lonlat_to_cell(NYC..., 5)
end

@testset "IGEO7 lookups: res-0 root map is the identity" begin
    roots = res0_cells()
    @test issorted(roots)
    @test collect(roots) == [z7_from_string(lpad(string(base), 2, '0')) for base in 0:11]
    @test cell_to_index.(collect(roots)) == 1:12
    @test index_to_cell.(1:12, 0) == collect(roots)

    for (cell, base) in zip(roots, 0:11)          # sibling's root_map, now 0:11
        z7 = cell_to_z7(cell)
        @test z7_base_cell(z7) == base
        @test z7_resolution(z7) == 0
        @test z7_is_pentagon(z7)
        @test z7_to_cell(z7) == cell
        @test z7_to_string(z7) == lpad(string(base), 2, '0')
    end
end

@testset "IGEO7 lookups: published Z7 example" begin
    # Published Z7 example: base cell 08 and resolution digits 0,0,4,3,3.
    published = z7_from_string("0800433")
    @test published == 0x80237fffffffffff
    @test z7_to_hex(published) == "80237fffffffffff"
    @test z7_from_hex("0x80237fffffffffff") == published
    @test z7_to_string(published) == "0800433"
    @test z7_base_cell(published) == 8
    @test z7_resolution(published) == 5
    @test z7_digit(published, 3) == 4
    @test z7_digit(published, 6) == 7
    @test !z7_is_pentagon(published)

    parent = z7_parent(published)
    @test z7_to_string(parent) == "080043"
    @test z7_parent(published, 2) == z7_from_string("0800")
    @test z7_is_descendant(published, parent)
    @test !z7_is_descendant(parent, published)

    north_root = z7_from_string("00")
    @test length(z7_children(north_root)) == 6
    @test z7_to_string.(collect(z7_children(north_root))) ==
          ["000", "001", "003", "004", "005", "006"]
    # these reach the *native* Z7 validators through the delegation layer, so
    # they carry the native error type; the layer's own rejections below stay
    # plain `ArgumentError`s
    @test_throws CORE.InvalidZ7Error z7_child(north_root, 2)

    south_root = z7_from_string("06")
    @test z7_to_string.(collect(z7_children(south_root))) ==
          ["060", "061", "062", "063", "064", "066"]
    @test_throws CORE.InvalidZ7Error z7_child(south_root, 5)
    @test length(z7_children(z7_from_string("001"))) == 7

    @test_throws CORE.InvalidZ7Error z7_from_string("0")
    @test_throws CORE.InvalidZ7Error z7_from_string("12")
    @test_throws CORE.InvalidZ7Error z7_from_string("002")
    @test_throws CORE.InvalidZ7Error z7_from_string("065")
    @test_throws CORE.InvalidZ7Error z7_from_string("0107")
    malformed_prefix = z7_from_string("00") & ~(UInt64(0x07) << 54)
    @test !CORE.is_valid_z7(malformed_prefix)
    @test_throws CORE.InvalidZ7Error z7_resolution(malformed_prefix)

    resolution20 = z7_from_string("01" * repeat("0", 20))
    @test z7_resolution(resolution20) == 20
    @test z7_to_string(resolution20) == "01" * repeat("0", 20)
    @test_throws CORE.InvalidZ7Error z7_to_cell(resolution20)
    @test_throws CORE.InvalidZ7Error z7_children(resolution20)
end

@testset "IGEO7 lookups: full-world indexing through the layer" begin
    for resolution in 0:5
        expected_ids = UInt64[]
        for root in res0_cells()
            append!(expected_ids, cell_to_children(root, resolution))
        end
        # clean-design ordering assertion: base blocks are ascending and each
        # block is emitted ascending, so the concatenation is already sorted
        @test issorted(expected_ids)
        @test length(expected_ids) == num_cells(resolution)
        sort!(expected_ids)

        ids = index_to_cell.(1:length(expected_ids), resolution)
        @test ids == expected_ids
        @test cell_to_index.(ids) == 1:length(ids)
    end

    for resolution in (0, 1, 18, 19)
        ncells = Int(num_cells(resolution))
        for index in unique((1, 2, ncells ÷ 2, ncells - 1, ncells))
            id = index_to_cell(index, resolution)
            @test CORE.get_resolution(id) == resolution
            @test cell_to_index(id) == index
        end
    end

    @test_throws BoundsError index_to_cell(0, 2)
    @test_throws BoundsError index_to_cell(493, 2)
    @test lonlat_to_index(NYC..., 7) == cell_to_index(lonlat_to_cell(NYC..., 7))
end

@testset "IGEO7 lookups: cell/z7 round trip" begin
    # Exhaustive through resolution 3 (aggregated — test_z7.jl owns the
    # per-cell detail), then high-resolution samples.
    for resolution in 0:3
        cells = index_to_cell.(1:Int(num_cells(resolution)), resolution)
        @test all(CORE.is_valid_z7, cells)
        @test all(c -> z7_resolution(c) == resolution, cells)
        @test all(c -> z7_to_cell(cell_to_z7(c)) == c, cells)
        @test all(c -> z7_from_string(z7_to_string(c)) == c, cells)
    end

    for resolution in (4, 7, 18, 19)
        count = Int(num_cells(resolution))
        for index in unique((1, 2, count ÷ 3, count ÷ 2, count - 1, count))
            cell = index_to_cell(index, resolution)
            z7 = cell_to_z7(cell)
            @test z7_to_cell(z7) == cell
            @test z7_resolution(z7) == resolution
        end
    end
end

@testset "IGEO7 lookups: construction and metadata" begin
    ids = collect(cell_to_children(first(res0_cells()), 2))
    lookup = IGeo7Lookup(ids; resolution=2, validate=true)

    @test length(lookup) == 41
    @test size(lookup) == (41,)
    @test lookup[1] == ids[1]
    @test firstindex(lookup) == 1 && lastindex(lookup) == 41
    @test collect(lookup) == ids
    @test DD.parent(lookup) === lookup.data
    @test DD.order(lookup) isa DD.Lookups.ForwardOrdered
    @test eltype(lookup) == UInt64

    md = DD.metadata(lookup)
    @test md["grid_name"] == "igeo7"
    @test md["resolution"] == 2
    @test md["indexing_scheme"] == "z7-u64"
    @test md["external_indexing_schemes"] == ["z7-string"]
    @test md["projection"] == "snyder-isea (standard ISEA placement)"
    # user metadata wins over the defaults, defaults fill the rest
    custom = IGeo7Lookup(ids; resolution=2, metadata=Dict{String,Any}("grid_name" => "x", "note" => 1))
    @test DD.metadata(custom)["grid_name"] == "x"
    @test DD.metadata(custom)["note"] == 1
    @test DD.metadata(custom)["indexing_scheme"] == "z7-u64"

    # ordering contract: lookup order == ascending id == dense index order,
    # and one base's res-2 descendants are a contiguous dense-index block
    @test issorted(lookup.data)
    indices = cell_to_index.(lookup.data)
    @test indices == first(indices):last(indices)

    @test_throws ArgumentError IGeo7Lookup(reverse(ids); resolution=2)
    @test_throws ArgumentError IGeo7Lookup([ids[1], ids[1]]; resolution=2)
    @test_throws ArgumentError IGeo7Lookup(ids; resolution=-1)
    @test_throws ArgumentError IGeo7Lookup(ids; resolution=CORE.MAX_RESOLUTION + 1)
    @test_throws ArgumentError IGeo7Lookup(ids; resolution=3, validate=true)
    # ...and without `validate` too: a Z7 id encodes its own resolution, so the
    # two endpoints are an O(1) test that a whole vector at another resolution
    # cannot pass. This one is what closes `IGeo7Lookup(res2_ids; resolution=3)`
    # constructing a lookup whose selectors then answer "not stored" for cells
    # that are stored.
    @test_throws ArgumentError IGeo7Lookup(ids; resolution=3)

    # mixed resolutions: an id at the *end* is caught by that endpoint check,
    # with or without validate (sorted, unique, and every id individually valid
    # — the resolution is what disagrees)
    mixed = sort(UInt64[ids[1], cell_to_children(ids[2], 3)[1]])
    @test all(CORE.is_valid_cell, mixed)
    @test CORE.get_resolution(last(mixed)) == 3
    @test_throws ArgumentError IGeo7Lookup(mixed; resolution=2, validate=true)
    @test_throws ArgumentError IGeo7Lookup(mixed; resolution=2)
    # ...but an id in the *interior* is invisible to it: both endpoints are
    # res-2, so only the O(n) `validate=true` pass sees the res-3 id between
    # them. That is the documented division of labour, not an oversight.
    interior = sort(UInt64[ids[1], cell_to_children(ids[2], 3)[1], ids[end]])
    @test CORE.get_resolution.(interior) == [2, 3, 2]
    @test_throws ArgumentError IGeo7Lookup(interior; resolution=2, validate=true)
    @test IGeo7Lookup(interior; resolution=2) isa IGeo7Lookup

    # The checks live in the inner constructor, so the positional form runs
    # them too — it used to be the default constructor, open beside the keyword
    # one that validated.
    @test IGeo7Lookup(ids, 2, Dict{String,Any}()).resolution == 2
    @test_throws ArgumentError IGeo7Lookup(reverse(ids), 2, Dict{String,Any}())
    @test_throws ArgumentError IGeo7Lookup(ids, 3, Dict{String,Any}())
    @test_throws ArgumentError IGeo7Lookup(ids, CORE.MAX_RESOLUTION + 1, Dict{String,Any}())

    # hex-string ids go through Helpers.to_uint64_id
    @test IGeo7Lookup(z7_to_hex.(ids[1:3]); resolution=2).data == ids[1:3]

    rebuilt = DD.rebuild(lookup; data=ids[1:5])
    @test rebuilt isa IGeo7Lookup
    @test rebuilt.resolution == 2
    @test rebuilt.data == ids[1:5]
    @test DD.metadata(rebuilt) == md
    @test DD.Lookups.reducelookup(lookup) isa DD.Lookups.NoLookup

    io = IOBuffer()
    show(io, MIME"text/plain"(), lookup)
    text = String(take!(io))
    @test occursin("resolution: 2", text)
    @test occursin("ncells: 41", text)
end

@testset "IGEO7 lookups: geometry helpers" begin
    ids = collect(cell_to_children(first(res0_cells()), 2))
    lookup = IGeo7Lookup(ids; resolution=2, validate=true)

    centers = cell_centers(lookup)
    @test length(centers) == length(lookup)
    @test centers == [cell_center(id) for id in ids]
    # the res-0 pentagon's center child keeps the vertex center
    @test cell_center(cell_to_parent(ids[1], 0)) == cell_center(first(res0_cells()))

    polys = cell_polygons(lookup)
    @test length(polys) == length(lookup)
    @test all(p -> GI.trait(p) isa GI.PolygonTrait, polys)
    @test GI.npoint(first(polys)) in (6, 7)          # closed pentagon/hexagon
    ring = GI.getexterior(first(polys))
    @test GI.getpoint(ring, 1) == GI.getpoint(ring, GI.npoint(ring))
    @test cell_polygon(ids[1]) == first(polys)

    # every cell's own center is inside its own polygon (planar test is safe
    # this far from the antimeridian: base 0 is the north-pole pentagon, so
    # use a mid-latitude base instead)
    mid = collect(cell_to_children(res0_cells()[3], 3))
    @test all(id -> GO.contains(cell_polygon(id), cell_center(id)), mid[1:20])
end

@testset "IGEO7 DimensionalData lookup" begin
    points = [NYC, LONDON, TOKYO]
    point_ids = [lonlat_to_cell(lon, lat, 2) for (lon, lat) in points]
    ids = sort(unique(point_ids))
    @test length(ids) == 3
    lookup = IGeo7Lookup(ids; resolution=2, validate=true)
    dim = Dim{:Geometry}(lookup)
    array = DD.DimArray(collect(1:length(lookup)), (dim,))

    @test array[Geometry=At(ids[1])] == 1
    @test only(DD.dims(array)[1][Contains(points[2])]) == point_ids[2]
    @test only(array[Geometry=Contains(points[3])]) == findfirst(==(point_ids[3]), ids)

    # At: hex string and missing id
    @test array[Geometry=At(z7_to_hex(ids[2]))] == 2
    @test_throws ArgumentError array[Geometry=At(UInt64(0))]

    # Contains: GeoInterface point, tuple, and a point in no listed cell
    @test DD.Lookups.selectindices(lookup, Contains(GI.Point(NYC))) ==
          DD.Lookups.selectindices(lookup, Contains(NYC))
    @test isempty(DD.Lookups.selectindices(lookup, Contains((0.0, 0.0))))

    # StandardIndices pass straight through
    @test DD.Lookups.selectindices(lookup, 2) == 2
    @test array[Geometry=2] == 2
end

@testset "IGEO7 lookups: Touching selector" begin
    ids = collect(cell_to_children(cell_to_parent(lonlat_to_cell(NYC..., 4), 2), 4))
    lookup = IGeo7Lookup(ids; resolution=4, validate=true)
    array = DD.DimArray(collect(1:length(lookup)), (Dim{:Geometry}(lookup),))
    nyc_index = findfirst(==(lonlat_to_cell(NYC..., 4)), ids)

    # a 0.2 deg x 0.2 deg box around Manhattan: much smaller than a res-4 cell
    # (~150 km across), so it can only touch the cell it sits in and any
    # neighbor whose edge crosses it
    box = GI.Polygon([GI.LinearRing([
        (-74.1, 40.6), (-73.9, 40.6), (-73.9, 40.8), (-74.1, 40.8), (-74.1, 40.6),
    ])])
    hits = DD.Lookups.selectindices(lookup, Touching(box))
    @test hits isa Vector{Int}
    @test issorted(hits)
    @test !isempty(hits)
    @test nyc_index in hits
    @test length(hits) <= 4
    @test all(i -> GO.intersects(box, cell_polygon(ids[i])), hits)
    @test collect(array[Geometry=Touching(box)]) == hits

    # a box on the far side of the world touches nothing in this lookup
    far = GI.Polygon([GI.LinearRing([
        (0.0, 0.0), (0.5, 0.0), (0.5, 0.5), (0.0, 0.5), (0.0, 0.0),
    ])])
    @test isempty(DD.Lookups.selectindices(lookup, Touching(far)))

    # Touching is a superset of Contains: a cell whose center is in the
    # geometry certainly intersects it
    big = GI.Polygon([GI.LinearRing([
        (-75.0, 40.0), (-73.0, 40.0), (-73.0, 41.5), (-75.0, 41.5), (-75.0, 40.0),
    ])])
    @test issubset(DD.Lookups.selectindices(lookup, Contains(big)),
        DD.Lookups.selectindices(lookup, Touching(big)))

    # the cell's own polygon always touches itself
    @test DD.Lookups.selectindices(lookup, Touching(cell_polygon(ids[nyc_index]))) ⊇
          [nyc_index]
end

# The globe-complete dimension — an ordinary `IGeo7Lookup` over computed rather
# than stored ids (docs/design/full_globe_lookups.md §1.2, §1.4, §1.5).
@testset "IGEO7 globe lookup" begin
    # Res 19 is 1.1e17 cells; the explicit id vector cannot exist, so every
    # assertion below doubles as an assertion that nothing on its path
    # materialized — a leak is an `OutOfMemoryError`, not a slow test.
    lookup = IGeo7Lookup(DGG.DGGSGlobeIds(DGG.IGEO7DGGS(), 19))
    @test lookup.data isa DGG.DGGSGlobeIds
    @test lookup.resolution == 19               # defaulted from the ids' level
    @test length(lookup) == DGG.num_cells(DGG.IGEO7DGGS(), 19) == num_cells(19)
    @test lookup[1] == DGG.ordinal_to_cell(DGG.IGEO7DGGS(), 19, 1)
    @test lookup.metadata["grid_name"] == "igeo7"
    @test lookup.metadata["indexing_scheme"] == "z7-u64"
    # `validate=true` is a documented no-op here rather than an O(n) decode of
    # every id: each one is `ordinal_to_cell`, a valid res-19 cell by
    # construction.
    @test IGeo7Lookup(DGG.DGGSGlobeIds(DGG.IGEO7DGGS(), 19); validate=true).data isa
          DGG.DGGSGlobeIds
    # `show` prints the length and the description, never the ids.
    text = sprint(show, MIME"text/plain"(), lookup)
    @test occursin("resolution: 19", text)
    @test occursin("ncells: $(length(lookup))", text)
    @test occursin("DGGSGlobeIds", text)
    # Another system's globe is rejected outright rather than falling through
    # to the generic constructor, which would materialize it to find out.
    @test_throws ArgumentError IGeo7Lookup(DGG.DGGSGlobeIds(DGG.H3DGGS(), 3))
    # A declared resolution disagreeing with the ids' level still fails the
    # endpoint check, as it does for an explicit id vector — at O(1) here.
    @test_throws ArgumentError IGeo7Lookup(DGG.DGGSGlobeIds(DGG.IGEO7DGGS(), 19);
        resolution=18)

    # Laziness survives every rebuild — the invariant most likely to regress
    # silently, since `DD.rebuild` defaults `data=l.data` and routes through the
    # keyword constructor: without its `DGGSGlobeIds` method a `format`, a `set`
    # or a metadata change would quietly materialize the globe, and the failure
    # would be a working answer rather than an error.
    @test DD.rebuild(lookup).data isa DGG.DGGSGlobeIds
    relabelled = DD.rebuild(lookup; metadata=Dict{String,Any}("title" => "globe"))
    @test relabelled.data isa DGG.DGGSGlobeIds
    @test relabelled.metadata["title"] == "globe"
    dim = Dim{:Geometry}(lookup)
    @test DD.val(DD.format(dim, Base.OneTo(length(lookup)))).data isa DGG.DGGSGlobeIds
    @test DD.val(DD.set(dim, Dim{:Geometry})).data isa DGG.DGGSGlobeIds

    # Degradation needs no new code: a non-scalar index falls through to
    # `AbstractArray`'s default, which materializes a `Vector`, so the rebuild
    # behind it takes the ordinary validating path to an ordinary partial
    # lookup. The discriminator is the type of `data`, never "was this a
    # rebuild".
    partial = lookup[1:10]
    @test partial isa IGeo7Lookup
    @test partial.data isa Vector{UInt64}
    @test partial.resolution == 19
    @test collect(partial) == [DGG.ordinal_to_cell(DGG.IGEO7DGGS(), 19, i) for i in 1:10]
    @test DD.Lookups.selectindices(partial, At(partial[4])) == 4
    @test_throws ArgumentError DD.Lookups.selectindices(partial, At(lookup[11]))

    # A globe small enough to hold is where the construction can be compared
    # against the explicit one it replaces, dimension and selectors included.
    globe = IGeo7Lookup(DGG.DGGSGlobeIds(DGG.IGEO7DGGS(), 0))
    @test collect(globe) == sort(collect(res0_cells()))
    array = DD.DimArray(collect(1:length(globe)), (Dim{:Geometry}(globe),))
    @test DD.lookup(array, Dim{:Geometry}).data isa DGG.DGGSGlobeIds  # after `format`
    @test DD.lookup(DD.set(array, Dim{:Geometry} => Dim{:Geometry}),
        Dim{:Geometry}).data isa DGG.DGGSGlobeIds
    @test array[Geometry=At(globe[7])] == 7
    @test only(array[Geometry=Contains(NYC)]) ==
          findfirst(==(lonlat_to_cell(NYC..., 0)), collect(globe))
    sliced = array[Geometry=1:5]
    @test DD.lookup(sliced, Dim{:Geometry}).data isa Vector{UInt64}
    @test collect(DD.lookup(sliced, Dim{:Geometry})) == collect(globe)[1:5]
    @test sliced[Geometry=At(globe[3])] == 3
end

end # module TestLookups
