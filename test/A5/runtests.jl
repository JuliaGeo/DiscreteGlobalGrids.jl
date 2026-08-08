module A5TestSuite

using Test
using DimensionalData

using DiscreteGlobalGrids
using DiscreteGlobalGrids.A5.A5Lookups
using DiscreteGlobalGrids: A5

const DD = DimensionalData

@testset "A5 native math" begin
    roots = res0_cells()
    @test length(roots) == 12
    @test roots === res0_cells()
    @test roots isa DiscreteGlobalGrids.Helpers.SmallList{12,UInt64}
    @test issorted(roots)
    @test A5Native.num_cells(-1) == 1
    @test A5Native.num_cells(0) == 12
    @test A5Native.num_cells(1) == 60
    @test A5Native.num_cells(2) == 240
    @test_throws ArgumentError A5Native.num_cells(30)

    cell = lonlat_to_cell(-73.9857, 40.7484, 3)
    @test A5Native.get_resolution(cell) == 3
    @test cell_area(3) > 0
    @test length(cell_boundary(cell; segments=1)) >= 6

    parent = cell_to_parent(cell, 2)
    @test A5Native.get_resolution(parent) == 2
    @test cell in cell_to_children(parent, 3)
    same = cell_to_children(cell, 3)
    @test only(same) == cell
    @test same isa DiscreteGlobalGrids.Helpers.SmallList{1,UInt64}
    @test_throws ArgumentError cell_to_parent(cell, 4)

    supported_parent = A5Native.serialize(
        A5Native.A5Cell(A5Native.ORIGINS[1], 0, 0, 28),
    )
    supported_children = cell_to_children(supported_parent, 30)
    @test length(supported_children) == 16
    @test all(==(30) ∘ A5Native.get_resolution, supported_children)

    unsupported_origin = A5Native.ORIGINS[10]
    unsupported_parent = A5Native.serialize(
        A5Native.A5Cell(unsupported_origin, unsupported_origin.first_quintant, 0, 28),
    )
    @test_throws ArgumentError cell_to_children(unsupported_parent, 30)

    for resolution in 0:8
        sampled = lonlat_to_cell(12.5, -34.25, resolution)
        @test A5Native.get_resolution(sampled) == resolution
        @test cell_to_parent(sampled, resolution) == sampled
        if resolution > 0
            sampled_parent = cell_to_parent(sampled)
            @test sampled in cell_to_children(sampled_parent, resolution)
        end
    end

    tokyo = lonlat_to_cell(139.7623824402441, 35.677369792795794, 4)
    @test tokyo == parse(UInt64, "8708000000000000"; base=16)
    center = cell_center(tokyo)
    @test isapprox(center[1], 141.11036610726177; atol=1e-10)
    @test isapprox(center[2], 34.33639263740687; atol=1e-10)

    antimeridian = parse(UInt64, "eb60000000000000"; base=16)
    ring = cell_boundary(antimeridian; segments=1)
    lons = first.(ring)
    @test maximum(lons) - minimum(lons) < 180
end

@testset "A5 DimensionalData lookup" begin
    @test_throws ArgumentError A5Lookup(UInt64[]; resolution=999)

    points = [
        (-73.9857, 40.7484),
        (-0.1278, 51.5074),
        (139.6917, 35.6895),
    ]
    point_ids = [lonlat_to_cell(lon, lat, 3) for (lon, lat) in points]
    ids = sort(unique(point_ids))
    lookup = A5Lookup(ids; resolution=3, validate=true)
    dim = A5Cells(lookup)
    array = DD.DimArray(collect(1:length(lookup)), (dim,))

    @test array[A5Cells(At(ids[1]))] == 1
    @test only(DD.dims(array)[1][Contains(points[2])]) == point_ids[2]
    @test only(array[A5Cells(Contains(points[3]))]) == findfirst(==(point_ids[3]), ids)

    # The cheap structural checks moved into the inner constructor, so no path
    # skips them — the positional form used to be the open default constructor
    # beside the keyword one that validated. An A5 index carries its resolution
    # in the position of its low marker bit, so checking the sorted vector's two
    # endpoints costs O(1) and rejects a whole vector at the wrong resolution:
    # without it every `Contains` above would answer "not stored" for cells that
    # are stored.
    @test A5Lookup(ids, 3, Dict{String,Any}()).resolution == 3
    @test_throws ArgumentError A5Lookup(ids; resolution=4)
    @test_throws ArgumentError A5Lookup(ids, 4, Dict{String,Any}())
    @test_throws ArgumentError A5Lookup(reverse(ids), 3, Dict{String,Any}())
    @test_throws ArgumentError A5Lookup(ids, A5Native.MAX_RESOLUTION + 1, Dict{String,Any}())
    @test DD.rebuild(lookup; data=ids[1:2]).data == ids[1:2]
    @test_throws ArgumentError DD.rebuild(lookup; data=reverse(ids))

    # `validate=true` is still what checks *every* id — the endpoints cannot
    # see a foreign resolution sitting between them (A5 encodes the resolution
    # in a low marker bit, so mixed resolutions interleave rather than group).
    interior = sort(vcat(ids, cell_to_children(ids[2], 4)[1]))
    @test A5Native.get_resolution.(interior) == [3, 4, 3, 3]
    @test A5Lookup(interior; resolution=3) isa A5Lookup
    @test_throws ArgumentError A5Lookup(interior; resolution=3, validate=true)
    # ...and an endpoint at the wrong resolution is caught either way.
    endpoint = sort(vcat(ids, cell_to_children(ids[1], 4)[1]))
    @test A5Native.get_resolution(first(endpoint)) == 4
    @test_throws ArgumentError A5Lookup(endpoint; resolution=3)
end

# The globe-complete dimension — an ordinary `A5Lookup` over computed rather
# than stored ids (docs/design/full_globe_lookups.md §1.2, §1.4, §1.5).
@testset "A5 globe lookup" begin
    # Res 25 is 1.7e16 cells; the explicit id vector cannot exist, so every
    # assertion below doubles as an assertion that nothing on its path
    # materialized — a leak is an `OutOfMemoryError`, not a slow test.
    lookup = A5Lookup(DGGSGlobeIds(A5DGGS(), 25))
    @test lookup.data isa DGGSGlobeIds
    @test lookup.resolution == 25               # defaulted from the ids' level
    @test length(lookup) == DiscreteGlobalGrids.num_cells(A5DGGS(), 25)
    @test lookup[1] == ordinal_to_cell(A5DGGS(), 25, 1)
    @test lookup.metadata["grid_name"] == "a5"
    # `validate=true` is a documented no-op here rather than an O(n) pass: every
    # id is `ordinal_to_cell`, hence a valid res-25 cell by construction.
    @test A5Lookup(DGGSGlobeIds(A5DGGS(), 25); validate=true).data isa DGGSGlobeIds
    # `show` prints the length and the description, never the ids.
    text = sprint(show, MIME"text/plain"(), lookup)
    @test occursin("resolution: 25", text)
    @test occursin("ncells: $(length(lookup))", text)
    @test occursin("DGGSGlobeIds", text)
    # Another system's globe is rejected outright rather than falling through
    # to the generic constructor, which would materialize it to find out.
    @test_throws ArgumentError A5Lookup(DGGSGlobeIds(H3DGGS(), 3))
    # A declared resolution disagreeing with the ids' level still fails the
    # endpoint check, as it does for an explicit id vector — at O(1) here.
    @test_throws ArgumentError A5Lookup(DGGSGlobeIds(A5DGGS(), 25); resolution=24)

    # Laziness survives every rebuild — the invariant most likely to regress
    # silently, since `DD.rebuild` defaults `data=l.data` and routes through the
    # keyword constructor: without its `DGGSGlobeIds` method a `format`, a `set`
    # or a metadata change would quietly materialize the globe, and the failure
    # would be a working answer rather than an error.
    @test DD.rebuild(lookup).data isa DGGSGlobeIds
    relabelled = DD.rebuild(lookup; metadata=Dict{String,Any}("title" => "globe"))
    @test relabelled.data isa DGGSGlobeIds
    @test relabelled.metadata["title"] == "globe"
    dim = A5Cells(lookup)
    @test DD.val(DD.format(dim, Base.OneTo(length(lookup)))).data isa DGGSGlobeIds
    @test DD.val(DD.set(dim, A5Cells)).data isa DGGSGlobeIds

    # Degradation needs no new code: a non-scalar index falls through to
    # `AbstractArray`'s default, which materializes a `Vector`, so the rebuild
    # behind it takes the ordinary validating path to an ordinary partial
    # lookup. The discriminator is the type of `data`, never "was this a
    # rebuild".
    partial = lookup[1:10]
    @test partial isa A5Lookup
    @test partial.data isa Vector{UInt64}
    @test partial.resolution == 25
    @test collect(partial) == [ordinal_to_cell(A5DGGS(), 25, i) for i in 1:10]
    @test DD.Lookups.selectindices(partial, At(partial[4])) == 4
    @test_throws ArgumentError DD.Lookups.selectindices(partial, At(lookup[11]))

    # A globe small enough to hold is where the construction can be compared
    # against the explicit one it replaces, dimension and selectors included.
    globe = A5Lookup(DGGSGlobeIds(A5DGGS(), 0))
    @test collect(globe) == sort(collect(res0_cells()))
    array = DD.DimArray(collect(1:length(globe)), (A5Cells(globe),))
    @test DD.lookup(array, A5Cells).data isa DGGSGlobeIds   # survived `format`
    @test DD.lookup(DD.set(array, A5Cells => A5Cells), A5Cells).data isa DGGSGlobeIds
    @test array[A5Cells(At(globe[7]))] == 7
    @test only(array[A5Cells(Contains((-73.9857, 40.7484)))]) ==
          findfirst(==(lonlat_to_cell(-73.9857, 40.7484, 0)), collect(globe))
    sliced = array[A5Cells(1:5)]
    @test DD.lookup(sliced, A5Cells).data isa Vector{UInt64}
    @test collect(DD.lookup(sliced, A5Cells)) == collect(globe)[1:5]
    @test sliced[A5Cells(At(globe[3]))] == 3
end

# Operations-kernel wiring of `A5DGGS` (src/A5/A5Kernel.jl); the file wraps
# itself in a module of its own.
include("test_a5_kernel.jl")

end # module A5TestSuite
