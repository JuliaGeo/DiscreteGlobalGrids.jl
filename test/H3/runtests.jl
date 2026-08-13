module H3TestSuite

using Test
using DimensionalData

using DiscreteGlobalGrids
using DiscreteGlobalGrids.H3.H3Lookups

const DD = DimensionalData

@testset "H3 native math" begin
    @test length(res0_cells()) == 122
    @test H3Native.num_cells(0) == 122
    @test H3Native.num_cells(1) == 842
    @test H3Native.num_cells(2) == 5882

    cell = lonlat_to_cell(-122.41795063018799, 37.775938728915946, 9)
    @test cell == parse(UInt64, "8928308280fffff"; base=16)
    @test H3Native.get_resolution(cell) == 9
    @test H3Native.is_valid_cell(cell)
    @test cell_area(cell) > 0

    center = cell_center(cell)
    @test isapprox(center[1], -122.41845932318309; atol=1e-12)
    @test isapprox(center[2], 37.776702349435695; atol=1e-12)

    parent = cell_to_parent(cell, 8)
    @test H3Native.get_resolution(parent) == 8
    children = cell_to_children(parent, 9)
    @test cell in children
    for (position, child) in enumerate(children)
        @test H3Native.cell_to_child_pos(child, 8) == position - 1
        @test H3Native.child_pos_to_cell(position - 1, parent, 9) == child
    end
    @test_throws ArgumentError H3Native.child_pos_to_cell(-1, parent, 9)

    pentagons = H3Native.get_pentagons(2)
    @test length(pentagons) == 12
    @test all(H3Native.is_pentagon, pentagons)
    @test length(cell_boundary(first(pentagons); closed_ring=false)) <= 10

    open_boundary = cell_boundary(cell; closed_ring=false)
    closed_boundary = cell_boundary(cell)
    @test length(closed_boundary) == length(open_boundary) + 1
    @test closed_boundary[1] == closed_boundary[end]
    @test H3Native.get_resolution("0x8928308280fffff") == 9
    @test H3Native.get_resolution("0X8928308280fffff") == 9
end

@testset "H3 DimensionalData lookup" begin
    points = [
        (-122.41795063018799, 37.775938728915946),
        (-0.1278, 51.5074),
        (139.6917, 35.6895),
    ]
    point_ids = [lonlat_to_cell(lon, lat, 4) for (lon, lat) in points]
    ids = sort(unique(point_ids))
    lookup = H3Lookup(ids; resolution=4, validate=true)
    dim = H3Cells(lookup)
    array = DD.DimArray(collect(1:length(lookup)), (dim,))

    @test array[H3Cells(At(ids[1]))] == 1
    @test only(DD.dims(array)[1][Contains(points[2])]) == point_ids[2]
    @test only(array[H3Cells(Contains(points[3]))]) == findfirst(==(point_ids[3]), ids)

    @test_throws ArgumentError H3Lookup(reverse(ids); resolution=4)
    @test_throws ArgumentError H3Lookup([ids[1], ids[1]]; resolution=4)

    # An H3 index encodes its own resolution, so a whole vector at the wrong
    # one is an O(1) rejection — and it has to be, because the failure is
    # silent otherwise: `Contains` hashes the point to a res-5 id, finds no
    # match among the res-4 ids, and reports "not stored" for a stored cell.
    @test_throws ArgumentError H3Lookup(ids; resolution=5)
    # ...and the valid path is untouched: the same ids at their own resolution
    # still resolve a point to the cell that holds it.
    @test only(array[H3Cells(Contains(points[1]))]) == findfirst(==(point_ids[1]), ids)

    # The cheap checks live in the inner constructor, so the positional form
    # (previously the open default constructor) runs them too...
    @test H3Lookup(ids, 4, Dict{String,Any}()).resolution == 4
    @test_throws ArgumentError H3Lookup(ids, 5, Dict{String,Any}())
    @test_throws ArgumentError H3Lookup(reverse(ids), 4, Dict{String,Any}())
    @test_throws ArgumentError H3Lookup(ids, H3Native.MAX_RESOLUTION + 1, Dict{String,Any}())
    # ...and so does `DD.rebuild`, which routes through the keyword form.
    @test DD.rebuild(lookup; data=ids[1:2]).data == ids[1:2]
    @test_throws ArgumentError DD.rebuild(lookup; data=reverse(ids))

    # H3 keeps the resolution in bits 52:55, above the base cell, so sorting
    # groups ids by resolution and a sorted vector's endpoints carry its min
    # and max: here the endpoint check happens to catch *any* mixed-resolution
    # input, with or without `validate`.
    mixed = sort(vcat(ids, cell_to_children(ids[1], 5)[1]))
    @test H3Native.get_resolution(last(mixed)) == 5
    @test_throws ArgumentError H3Lookup(mixed; resolution=4)
    @test_throws ArgumentError H3Lookup(mixed; resolution=4, validate=true)

    # `validate=true` still buys what no endpoint test can: structural validity
    # of every id. This one is tagged res-4, so it passes the endpoint check,
    # but one of its unused digit slots is cleared — not a cell.
    bogus = sort(vcat(ids, ids[1] & ~(UInt64(7) << 30)))
    @test H3Native.get_resolution(first(bogus)) == 4
    @test !H3Native.is_valid_cell(first(bogus))
    @test H3Lookup(bogus; resolution=4) isa H3Lookup
    @test_throws ArgumentError H3Lookup(bogus; resolution=4, validate=true)
end

# The globe-complete dimension — an ordinary `H3Lookup` over computed rather
# than stored ids (docs/design/full_globe_lookups.md §1.2, §1.4, §1.5).
@testset "H3 globe lookup" begin
    # Res 15 is 5.7e14 cells, an explicit id vector of some 4.6 PB. Every
    # assertion below therefore doubles as an assertion that nothing on its
    # path materialized: a leak here is an `OutOfMemoryError`, not a slow test.
    lookup = H3Lookup(DGGSGlobeIds(H3DGGS(), 15))
    @test lookup.data isa DGGSGlobeIds
    @test lookup.resolution == 15               # defaulted from the ids' level
    @test length(lookup) == DiscreteGlobalGrids.num_cells(H3DGGS(), 15)
    @test lookup[1] == ordinal_to_cell(H3DGGS(), 15, 1)
    @test lookup.metadata["grid_name"] == "h3"
    @test lookup.metadata["resolution"] == 15
    # `validate=true` is a documented no-op here, not the O(n) trap it would be:
    # every id is `ordinal_to_cell`, hence a valid res-15 cell by construction.
    @test H3Lookup(DGGSGlobeIds(H3DGGS(), 15); validate=true).data isa DGGSGlobeIds
    # `show` prints the length and the description, never the ids.
    text = sprint(show, MIME"text/plain"(), lookup)
    @test occursin("resolution: 15", text)
    @test occursin("ncells: 569707381193162", text)
    @test occursin("DGGSGlobeIds", text)
    # Another system's globe is rejected outright rather than falling through
    # to the generic constructor, which would materialize it to find out.
    @test_throws ArgumentError H3Lookup(DGGSGlobeIds(HEALPixDGGS(), 3))
    # A declared resolution disagreeing with the ids' level still fails the
    # endpoint check, as it does for an explicit id vector — at O(1) here.
    @test_throws ArgumentError H3Lookup(DGGSGlobeIds(H3DGGS(), 15); resolution=14)

    # Laziness survives every rebuild. This is the invariant most likely to
    # regress silently: `DD.rebuild` defaults `data=l.data` and routes through
    # the keyword constructor, so without that constructor's `DGGSGlobeIds`
    # method a `format`, a `set` or a metadata change would quietly materialize
    # the whole globe — a working answer, not an error.
    @test DD.rebuild(lookup).data isa DGGSGlobeIds
    relabelled = DD.rebuild(lookup; metadata=Dict{String,Any}("title" => "globe"))
    @test relabelled.data isa DGGSGlobeIds
    @test relabelled.metadata["title"] == "globe"
    dim = H3Cells(lookup)
    @test DD.val(DD.format(dim, Base.OneTo(length(lookup)))).data isa DGGSGlobeIds
    @test DD.val(DD.set(dim, H3Cells)).data isa DGGSGlobeIds

    # Degradation, which needs no new code: a non-scalar index falls through to
    # `AbstractArray`'s default and returns a real `Vector`, so `rebuild` sees a
    # `Vector` and takes the ordinary validating path to an ordinary partial
    # lookup. The discriminator is the type of `data`, never "was this a
    # rebuild".
    partial = lookup[1:10]
    @test partial isa H3Lookup
    @test partial.data isa Vector{UInt64}
    @test partial.resolution == 15
    @test collect(partial) == [ordinal_to_cell(H3DGGS(), 15, i) for i in 1:10]
    @test DD.Lookups.selectindices(partial, At(partial[4])) == 4
    @test_throws ArgumentError DD.Lookups.selectindices(partial, At(lookup[11]))

    # A globe small enough to hold is where the construction can be compared
    # against the explicit one it replaces, dimension and selectors included.
    globe = H3Lookup(DGGSGlobeIds(H3DGGS(), 0))
    @test collect(globe) == collect(H3Lookup(sort(res0_cells()); resolution=0))
    array = DD.DimArray(collect(1:length(globe)), (H3Cells(globe),))
    @test DD.lookup(array, H3Cells).data isa DGGSGlobeIds   # survived `format`
    @test DD.lookup(DD.set(array, H3Cells => H3Cells), H3Cells).data isa DGGSGlobeIds
    @test array[H3Cells(At(globe[7]))] == 7
    @test only(array[H3Cells(Contains((-122.4, 37.8)))]) ==
          findfirst(==(lonlat_to_cell(-122.4, 37.8, 0)), collect(globe))
    sliced = array[H3Cells(1:10)]
    @test DD.lookup(sliced, H3Cells).data isa Vector{UInt64}
    @test collect(DD.lookup(sliced, H3Cells)) == collect(globe)[1:10]
    @test sliced[H3Cells(At(globe[3]))] == 3
end

# Operations-kernel wiring of `H3DGGS` (src/H3/H3Kernel.jl); the file wraps
# itself in a module of its own.
include("test_h3_kernel.jl")

# Edge neighbors (`grid_disk` wrappers, `cell_neighbors` wiring) and the
# lookup operations built on them (`neighbor_indices`, `stencil`, `zonal`);
# also its own module.
include("test_neighbors.jl")

# The subtree rim (`DGG.subtree_border` and the digit automaton behind it);
# also its own module. It leans on `cell_neighbors` for ground truth, so it
# runs after the suite that establishes that operation.
include("test_border.jl")

end # module H3TestSuite
