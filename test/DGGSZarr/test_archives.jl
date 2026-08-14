# test_archives.jl — the reader against real DGGS-convention archives
# (src/io/DGGSZarr.jl).
#
# The archives are written by the Python IGEO7/Z7 tooling and live outside this
# repository, so this file is skipped when they are absent (see `ARCHIVE_ROOT`
# in runtests.jl). When they are present it is the only place the reader is
# checked against bytes it did not write itself, which is what makes it worth
# the external dependency.
#
# The reference archives come in (ranges, dense) pairs holding the *same* cells:
# `pori_z7_r10_ranges.zarr` and `pori_z7_r10.zarr` are one dataset stored two
# ways. That pairing is the strongest assertion available here — the range table
# expanded arithmetically must reproduce, bit for bit, the dense `cell_ids`
# array the writer stored independently. If the base-7 number line were even
# slightly wrong, the two would diverge.

using Test

module TestArchives

using Test

using DimensionalData
using DimensionalData.Lookups: At
using GeoInterface

using DiscreteGlobalGrids: Helpers, ISEA, IGeo7
using DiscreteGlobalGrids.DGGSZarr
using DiscreteGlobalGrids.DGGSZarr: DGGSZarrInfo, julia_dimensions
import DiscreteGlobalGrids as DGG
import YAXArrays
import Zarr

const SUITE = parentmodule(@__MODULE__)

if !SUITE.have_archives()
    @info """
    skipping the DGGS-Zarr archive tests: reference archives not found under
    $(SUITE.ARCHIVE_ROOT). Set DGGS_ZARR_TEST_DATA to the directory holding
    pori_z7_r10{,_ranges}.zarr and pori_z7_r12{,_ranges}.zarr to run them.
    """
else

    @testset "convention metadata" begin
        ds = open_dggs_dataset(SUITE.archive_path("pori_z7_r10_ranges.zarr"))
        info = dggs_info(ds)
        @test info.name == "igeo7"
        @test info.level == 10
        @test info.spatial_dimension == "cell_ids"
        @test info.coordinate == "cell_id_ranges"
        @test info.compression == "ranges"
        @test info.ellipsoid["name"] == "wgs84"
        # The grid-specific placement lands in `extra`, which is where the
        # reader has to look for it — it is not a convention-level field.
        @test info.extra["dggs_vert0_lon"] == 11.2
        @test DGGSZarr.igeo7_vert0_lon(info) == 11.2

        # The dense sibling declares the same grid, differently compressed.
        dense = dggs_info(open_dggs_dataset(SUITE.archive_path("pori_z7_r10.zarr")))
        @test dense.compression == "none"
        @test dense.coordinate == "cell_ids"
        @test dense.level == info.level
    end

    @testset "Zarr reverses dimension order" begin
        # The bug this pins: on-disk (R, 2) arrives as (2, R), so a reader that
        # trusts `_ARRAY_DIMENSIONS` as written finds 2 ranges instead of R.
        group = Zarr.zopen(SUITE.archive_path("pori_z7_r10_ranges.zarr"), "r")
        table = group["cell_id_ranges"]
        @test table.attrs["_ARRAY_DIMENSIONS"] == ["ranges", "bounds"]
        @test julia_dimensions(table) == ["bounds", "ranges"]
        @test size(table) == (2, 136)
        # ... and the reader gets R right anyway.
        ids = dggs_cell_ids(open_dggs_dataset(SUITE.archive_path("pori_z7_r10_ranges.zarr")))
        @test IGeo7.z7_nranges(ids) == 136
    end

    @testset "ranges and dense archives decode to identical cells" begin
        for pair in SUITE.ARCHIVE_PAIRS
            ranges_ds = open_dggs_dataset(SUITE.archive_path(pair.ranges))
            dense_ds = open_dggs_dataset(SUITE.archive_path(pair.dense))

            range_ids = dggs_cell_ids(ranges_ds)
            dense_ids = dggs_cell_ids(dense_ds)
            @test range_ids isa IGeo7.Z7RangeIds
            @test dense_ids isa IGeo7.Z7CachedIds
            @test length(range_ids) == length(dense_ids)

            # The assertion this file exists for.
            expanded = collect(range_ids)
            @test expanded == collect(dense_ids)

            @test all(id -> IGeo7.z7_resolution(id) == pair.level, expanded)
            @test Helpers.strictly_increasing(expanded)
            @test all(IGeo7.is_valid_z7, expanded)

            # Positions agree both ways, over every cell.
            @test all(i -> DGG.cell_position(range_ids, expanded[i]) == i, eachindex(expanded))
            @test all(i -> DGG.cell_position(dense_ids, expanded[i]) == i, eachindex(expanded))
        end
    end

    @testset "opening is lazy" begin
        path = SUITE.archive_path("pori_z7_r12_ranges.zarr")
        ds = open_dggs_dataset(path)

        # Data variables are still the archive's Zarr arrays: nothing was read.
        @test ds.elevation.data isa Zarr.ZArray
        @test size(ds.elevation) == (158430,)

        # The index costs O(R). 1073 ranges against 158,430 cells, so the whole
        # coordinate is far smaller than the ids it names would be.
        ids = dggs_cell_ids(ds)
        @test IGeo7.z7_nranges(ids) == 1073
        @test Base.summarysize(ids) < 8 * length(ids)

        # A dense archive defers even its stored coordinate until something asks.
        dense = open_dggs_dataset(SUITE.archive_path("pori_z7_r12.zarr"))
        dense_ids = dggs_cell_ids(dense)
        @test !IGeo7.z7_is_materialized(dense_ids)

        # Printing the dataset's structure is not "something asks": the compact
        # display prints only the two endpoint ids, which the vector serves
        # without a full read. This is the whole point of opening lazily.
        text = sprint(show, MIME"text/plain"(), dense)
        @test occursin("IGeo7Lookup", text)
        @test occursin("ncells: 158430", text)
        @test !IGeo7.z7_is_materialized(dense_ids)

        # Selecting by cell id is.
        DGG.cell_position(dense_ids, dense_ids[1])
        @test IGeo7.z7_is_materialized(dense_ids)
    end

    @testset "selection by cell id" begin
        for pair in SUITE.ARCHIVE_PAIRS
            ranges_ds = open_dggs_dataset(SUITE.archive_path(pair.ranges))
            dense_ds = open_dggs_dataset(SUITE.archive_path(pair.dense))
            ids = collect(dggs_cell_ids(ranges_ds))

            for i in (1, 2, 100, length(ids) ÷ 2, length(ids))
                id = ids[i]
                a = only(ranges_ds.elevation[cell_ids=At(id)])
                b = only(dense_ds.elevation[cell_ids=At(id)])
                # Both archives hold the same DEM, so the two must agree
                # (NaN included — these are float32 elevations with nodata).
                @test (isnan(a) && isnan(b)) || a == b
            end

            # A cell of the right resolution that the archive does not hold.
            absent = IGeo7.index_to_cell(1, pair.level)
            @test absent ∉ ids
            @test_throws ArgumentError ranges_ds.elevation[cell_ids=At(absent)]
        end
    end

    @testset "geometry honours placement and datum" begin
        ds = open_dggs_dataset(SUITE.archive_path("pori_z7_r10_ranges.zarr"))
        ids = collect(dggs_cell_ids(ds))
        centers = cell_centers(ds)
        @test length(centers) == length(ids)

        orientation = IGeo7.vert0_lon_orientation(11.2)
        lon, lat = IGeo7.cell_center(ids[1]; orientation)

        # Longitude is untouched by the datum; latitude is not.
        @test centers[1][1] ≈ lon
        @test centers[1][2] ≈ Helpers.authalic_to_geodeticd(Helpers.WGS84_AUTHALIC, lat)
        @test centers[1][2] != lat
        # The size of the mistake, if the datum were skipped: ~0.115°, ~13 km.
        @test abs(centers[1][2] - lat) > 0.1

        # ... and of skipping the placement: the archive's 11.2 against the
        # package default of 11.25 is a 0.05° shift in longitude.
        default_lon, _ = IGeo7.cell_center(ids[1])
        @test default_lon - lon ≈ ISEA.ISEA_LON0 - 11.2 atol = 1e-9

        # Every center round trips back to its own cell through the same
        # placement and datum — the check that the two conversions compose.
        @test all(eachindex(ids)) do i
            clon, clat = centers[i]
            dggs_cell_at(ds, clon, clat) == ids[i]
        end

        boundaries = cell_boundaries(ds)
        @test length(boundaries) == length(ids)
        # Hexagons, plus the closing repeat; no pentagon in this extent.
        @test all(b -> length(first(GeoInterface.coordinates(b))) == 7, boundaries)
    end

    @testset "sel_latlon and sel_bbox" begin
        for pair in SUITE.ARCHIVE_PAIRS
            ranges_ds = open_dggs_dataset(SUITE.archive_path(pair.ranges))
            dense_ds = open_dggs_dataset(SUITE.archive_path(pair.dense))
            centers = cell_centers(ranges_ds)

            # A point inside a known cell selects that cell, in both forms.
            i = length(centers) ÷ 3
            lon, lat = centers[i]
            a = only(sel_latlon(ranges_ds.elevation, lon, lat))
            b = only(sel_latlon(dense_ds.elevation, lon, lat))
            @test (isnan(a) && isnan(b)) || a == b

            # Somewhere far outside the archive's extent.
            @test_throws ArgumentError sel_latlon(ranges_ds.elevation, 0.0, 0.0)

            # A bbox selects by center, and the two storage forms agree.
            box = (26.6, 26.7, 58.15, 58.2)
            expected = count(c -> box[1] <= c[1] <= box[2] && box[3] <= c[2] <= box[4], centers)
            @test size(sel_bbox(ranges_ds.elevation, box...)) == (expected,)
            @test size(sel_bbox(dense_ds.elevation, box...)) == (expected,)
            @test expected > 0

            # Dataset-wide selection keeps every variable and the properties.
            sub = sel_bbox(ranges_ds, box...)
            @test keys(sub.cubes) == keys(ranges_ds.cubes)
            @test size(first(values(sub.cubes))) == (expected,)
            @test dggs_info(sub).compression == dggs_info(ranges_ds).compression
        end
    end

    @testset "rejections" begin
        # A Zarr group with no dggs block is not a DGGS archive.
        tmp = mktempdir()
        Zarr.zgroup(joinpath(tmp, "plain"))
        @test_throws ArgumentError open_dggs_dataset(joinpath(tmp, "plain"))

        # An unwired grid names itself rather than failing inside the reader.
        info = DGGSZarrInfo(Dict("dggs" => Dict(
            "name" => "healpix", "refinement_level" => 4,
            "spatial_dimension" => "cells", "coordinate" => "cell_ids",
            "compression" => "none")))
        @test info.name == "healpix"
        err = try
            DGGSZarr._decode_coordinate(info, nothing, 1)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("healpix", err.msg)
    end

end # if have_archives

end # module TestArchives
