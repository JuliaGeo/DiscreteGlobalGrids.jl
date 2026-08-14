using Test
import DimensionalData as DD

const _CI_SYSTEM = IGEO7DGGS()
const _CI_ROOT = UInt64(0x0c4d9fffffffffff)

function _ci_allocated(f, args...)
    f(args...)
    return @allocated f(args...)
end

@testset "IGEO7 custom indices" begin
    @testset "identity and coordinates" begin
        id = z7_from_string("0800433")
        index = IGEO7Index(id)
        @test index.id == id
        @test IGEO7Index(z7_to_string(id)) == index
        @test UInt64(index) == id
        @test get_resolution(index) == 5
        @test repr(index) == "IGEO7Index(\"0800433\")"
        @test isbitstype(IGEO7Index)
        @test isbitstype(RelativeIGEO7Index)
        @test isbitstype(HexIndex)

        directions = (
            HexIndex(1, -1, 0),
            HexIndex(1, 0, -1),
            HexIndex(0, 1, -1),
            HexIndex(-1, 1, 0),
            HexIndex(-1, 0, 1),
            HexIndex(0, -1, 1),
        )
        @test all(h -> h.i + h.j + h.k == 0, directions)
        @test_throws ArgumentError HexIndex(1, 1, 1)
        for (code, hex) in enumerate(directions)
            relative = RelativeIGEO7Index(hex, 8)
            @test convert(HexIndex, relative) == hex
            @test repr(relative) == "RelativeIGEO7Index($(repr(hex)), 8)"
            @test relative + (-relative) == RelativeIGEO7Index(8)
            @test !iszero(relative)
            @test iszero(zero(relative))
            @test zero(relative).resolution == relative.resolution
            @test directioncode(relative) == code
            @test _ci_allocated(directioncode, relative) == 0
        end
        @test directioncode(RelativeIGEO7Index(8)) == 0
        @test_throws ArgumentError directioncode(RelativeIGEO7Index(2, 0, 8))
        @test_throws ArgumentError directioncode(RelativeIGEO7Index(1, -1, 8))
        @test_throws DimensionMismatch RelativeIGEO7Index(1) + RelativeIGEO7Index(2)
        @test_throws InexactError RelativeIGEO7Index(
            HexIndex(1, typemax(Int64), typemin(Int64)), 8)
    end

    @testset "local lattice arithmetic" begin
        tile = DGGSSubtreeIds(_CI_SYSTEM, 5, _CI_ROOT, 8)
        sample_positions = unique(vcat(1, collect(1:17:length(tile)), length(tile)))
        indices = IGEO7Index[tile[i] for i in sample_positions]

        for (a, b) in zip(indices[2:end], indices[1:end-1])
            displacement = a - b
            @test b + displacement == a
            @test a - displacement == b
            @test (b + displacement) - b == displacement
        end
        a, b = indices[2], indices[1]
        displacement = a - b
        @test _ci_allocated(-, a, b) == 0
        @test _ci_allocated(+, b, displacement) == 0

        zero_index = first(indices)
        @test zero_index - zero_index == RelativeIGEO7Index(8)
        @test zero_index + RelativeIGEO7Index(8) == zero_index
        @test_throws DimensionMismatch trytranslate(
            zero_index, RelativeIGEO7Index(7))
        @test isnothing(trytranslate(
            zero_index, RelativeIGEO7Index(typemax(Int64), typemax(Int64), 8)))

        north_a = IGEO7Index(z7_from_string("0001"))
        north_b = IGEO7Index(z7_from_string("0003"))
        @test_throws ArgumentError north_a - north_b
        @test north_a - north_a == RelativeIGEO7Index(2)

        other_base = IGEO7Index(z7_from_string("0101"))
        @test_throws ArgumentError north_a - other_base
    end

    @testset "neighbors and DimArray indexing" begin
        tile = DGGSSubtreeIds(_CI_SYSTEM, 5, _CI_ROOT, 8)
        lookup = IGeo7.IGeo7Lookups.IGeo7Lookup(tile, 8, Dict{String,Any}())
        values = collect(1:length(tile))
        A = DD.DimArray(values, DD.Dim{:cells}(lookup))

        indices = eachindex(A)
        @test indices isa IGeo7.IGEO7Indices
        @test length(indices) == length(A)
        @test issorted(indices)
        @test all(in(indices), indices)
        @test _ci_allocated(in, indices[10], indices) == 0
        @test all(i -> A[indices[i]] == A[i], eachindex(values))
        @test position_to_cell(A, 10) == indices[10]
        @test cell_to_position(A, indices[10]) == 10
        @test checkbounds(Bool, A, indices[10])

        absent = IGEO7Index(index_to_cell(1, 8))
        @test absent ∉ indices
        @test !checkbounds(Bool, A, absent)
        @test_throws BoundsError A[absent]
        @test_throws BoundsError cellbearing(A, indices[1], absent)
        @test_throws BoundsError celldistance(A, indices[1], absent)

        position = first(subtree_interior_positions(tile))
        index = indices[position]
        global_neighbors = neighbors(_CI_SYSTEM, index)
        stored_neighbors = neighbors(A, index)
        @test length(global_neighbors) == 6
        @test neighbors(index) == global_neighbors
        @test _ci_allocated(neighbors, _CI_SYSTEM, index) == 0
        @test stored_neighbors == global_neighbors
        @test directioncode.(global_neighbors .- Ref(index)) == 1:6
        @test all(n -> n - index in (
            RelativeIGEO7Index(HexIndex(1, -1, 0), 8),
            RelativeIGEO7Index(HexIndex(1, 0, -1), 8),
            RelativeIGEO7Index(HexIndex(0, 1, -1), 8),
            RelativeIGEO7Index(HexIndex(-1, 1, 0), 8),
            RelativeIGEO7Index(HexIndex(-1, 0, 1), 8),
            RelativeIGEO7Index(HexIndex(0, -1, 1), 8),
        ), global_neighbors)
        @test all(n -> A[n] == values[cell_to_position(A, n)], stored_neighbors)
        neighbor = first(stored_neighbors)
        distance = celldistance(A, index, neighbor)
        @test distance > 0
        @test celldistance(A, neighbor, index) ≈ distance
        @test celldistance(A, index, index) == 0.0
        bearing = cellbearing(A, index, neighbor)
        @test 0.0 <= bearing < 360.0
        @test cellbearing(A, index, index) == 0.0

        from_lon, from_lat = cell_center(index.id)
        to_lon, to_lat = cell_center(neighbor.id)
        phi_from, phi_to = deg2rad(from_lat), deg2rad(to_lat)
        delta_lon = deg2rad(to_lon - from_lon)
        expected_bearing = mod(rad2deg(atan(
            sin(delta_lon) * cos(phi_to),
            cos(phi_from) * sin(phi_to) -
            sin(phi_from) * cos(phi_to) * cos(delta_lon),
        )), 360.0)
        @test bearing ≈ expected_bearing

        @test cellarea(A, index) == cell_area(index.id)
        @test cellarea(A, index) > 0
        @test_throws BoundsError cellarea(A, absent)

        border_index = indices[first(subtree_border_positions(tile))]
        expected_stored = filter(n -> checkbounds(Bool, A, n),
            collect(neighbors(_CI_SYSTEM, border_index)))
        @test collect(neighbors(A, border_index)) == expected_stored
        expected_edges = IGEO7Index.(border_descendants(_CI_ROOT, 8))
        @test edges(indices) == expected_edges
        @test edges(A) == expected_edges

        A[indices[1]] = -1
        @test A[1] == -1
        @test sum(copy(A)) == sum(A)
        @test parent(A .+ 1) == parent(A) .+ 1
        @test parent(map(identity, A)) == parent(A)

        plain = DD.DimArray(collect(1:4), DD.Dim{:x}(1:4))
        @test eachindex(plain) == Base.OneTo(4)

        explicit_lookup = IGeo7.IGeo7Lookups.IGeo7Lookup(
            collect(tile); resolution=8, validate=true)
        explicit = DD.DimArray(values, DD.Dim{:cells}(explicit_lookup))
        explicit_indices = eachindex(explicit)
        @test explicit_indices.contiguous
        @test all(in(explicit_indices), explicit_indices)
        @test absent ∉ explicit_indices
        @test edges(explicit) == expected_edges

        sparse_lookup = IGeo7.IGeo7Lookups.IGeo7Lookup(
            collect(tile)[1:2:end]; resolution=8, validate=true)
        sparse = DD.DimArray(values[1:2:end], DD.Dim{:cells}(sparse_lookup))
        sparse_indices = eachindex(sparse)
        @test !sparse_indices.contiguous
        @test first(indices) in sparse_indices
        @test indices[2] ∉ sparse_indices
        sparse_truth = [cell for cell in sparse_indices
                        if any(n -> n ∉ sparse_indices, neighbors(cell))]
        @test edges(sparse) == sparse_truth

        globe_ids = DGGSGlobeIds(_CI_SYSTEM, 2)
        globe_lookup = IGeo7.IGeo7Lookups.IGeo7Lookup(globe_ids)
        globe = DD.DimArray(zeros(length(globe_ids)), DD.Dim{:cells}(globe_lookup))
        @test isempty(edges(globe))
    end
end
