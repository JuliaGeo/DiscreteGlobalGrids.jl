# Convention detection, decoding, merging and stamping, driven by hand-built
# snapshots that reproduce real stores verbatim (dggs-storage-landscape.md
# sections 4.1-4.5).

module DGGSIOConventionTests

using Test
import DiscreteGlobalGrids as DGG

# The conventions are plain-data metadata logic: no Zarr, no arrays, and no
# lookups. `IO` is the package itself, named here because these tests describe a
# layer rather than a package.
const IO = DGG
const DENSE = IO.ENCODING_REGISTRY["none"]
const RANGES = IO.ENCODING_REGISTRY["ranges"]
const IMPLICIT = IO.ENCODING_REGISTRY["implicit"]

d(pairs...) = Dict{String,Any}(pairs...)

# ---------------------------------------------------------------------------
# Fixtures, verbatim from the landscape report
# ---------------------------------------------------------------------------

# Section 4.1: `pori_z7_r10.zarr/.zattrs`, complete and unedited. Note
# `semimajor_axis` -- the spelling the wild stores use and the schema forbids.
pori_group_attrs(; level=10, compression="none", coordinate="cell_ids") = d(
    "clipper_scale_factor" => 10000000,
    "dggs" => d(
        "compression" => compression,
        "coordinate" => coordinate,
        "dggs_vert0_azimuth" => 0.0,
        "dggs_vert0_lat" => 58.28252559,
        "dggs_vert0_lon" => 11.2,
        "ellipsoid" => d(
            "inverse_flattening" => 298.257223563,
            "name" => "wgs84",
            "semimajor_axis" => 6378137.0),
        "name" => "igeo7",
        "refinement_level" => level,
        "rotation_pattern" => "alternating_cw_odd_ccw_even",
        "spatial_dimension" => "cell_ids"),
    "regridder" => "xdggs_dggrid4py.mapblocks_nearestcentroid",
    "source_crs" => "EPSG:3301",
    "source_path" => "/Users/akmoch/dev/build/igeo7_z7_xarray_paper/data/input/merit_dem_pori_cog.tif",
    "zarr_conventions" => Any[d(
        "description" => "Discrete Global Grid Systems convention for zarr",
        "name" => "dggs",
        "schema_url" => "https://raw.githubusercontent.com/zarr-conventions/dggs/refs/tags/v1/schema.json",
        "spec_url" => "https://github.com/zarr-conventions/dggs/blob/v1/README.md",
        "uuid" => "7b255807-140c-42ca-97f6-7a1cfecdbc38")])

# Section 4.2: the flat coordinate attrs are convention B in the
# xdggs-dggrid4py dialect, matching `IGEO7Info.to_dict()` field for field.
igeo7_coord_attrs(; level=10) = d(
    "_ARRAY_DIMENSIONS" => Any["cell_ids"],
    "grid_name" => "igeo7",
    "level" => level,
    "igeo7_dggs_vert0_lon" => 11.2,
    "igeo7_wgs84_geodetic_conversion" => true)

# Section 4.2: `pori_z7_r10.zarr`, 4 arrays, all 1-D on dimension `cell_ids`.
function pori_dense_snapshot(; level=10, coord_level=level)
    arrays = [
        IO.ArrayEntry(name="cell_ids", attrs=igeo7_coord_attrs(level=coord_level),
            shape=(3101,), eltype=UInt64, dims=["cell_ids"]),
        IO.ArrayEntry(name="elevation", attrs=d("_ARRAY_DIMENSIONS" => Any["cell_ids"]),
            shape=(3101,), eltype=Float32, dims=["cell_ids"]),
        IO.ArrayEntry(name="slope_geodesic",
            attrs=d("_ARRAY_DIMENSIONS" => Any["cell_ids"],
                "long_name" => "slope magnitude (FDA)", "units" => "m/m"),
            shape=(3101,), eltype=Float32, dims=["cell_ids"]),
        IO.ArrayEntry(name="slope_lookup",
            attrs=d("_ARRAY_DIMENSIONS" => Any["cell_ids"],
                "long_name" => "slope magnitude (FDA)", "units" => "m/m"),
            shape=(3101,), eltype=Float32, dims=["cell_ids"])]
    return IO.StoreSnapshot(identifier="gs://geo-assets/igeo7-zarr/pori_z7_r10.zarr",
        attrs=pori_group_attrs(level=level), arrays=arrays)
end

# Section 4.3: the same data with a `(136, 2)` ranges axis. The flat convention
# B attrs ride on that 2-D array, and `spatial_dimension` names a dimension no
# array has.
function pori_ranges_snapshot()
    attrs = pori_group_attrs(compression="ranges", coordinate="cell_id_ranges")
    attrs["coordinates"] = "cell_id_ranges"
    arrays = [
        IO.ArrayEntry(name="cell_id_ranges",
            attrs=merge(igeo7_coord_attrs(),
                d("_ARRAY_DIMENSIONS" => Any["ranges", "bounds"])),
            shape=(136, 2), eltype=UInt64, dims=["ranges", "bounds"]),
        IO.ArrayEntry(name="elevation", attrs=d("_ARRAY_DIMENSIONS" => Any["cell_ids"]),
            shape=(3101,), eltype=Float32, dims=["cell_ids"])]
    return IO.StoreSnapshot(
        identifier="gs://geo-assets/igeo7-zarr/pori_z7_r10_ranges.zarr",
        attrs=attrs, arrays=arrays)
end

# ---------------------------------------------------------------------------
# Convention A -- zarr-conventions/dggs
# ---------------------------------------------------------------------------

@testset "convention A detects on the uuid and decodes the dggs object" begin
    snap = pori_dense_snapshot()
    c = IO.ZarrDGGSConvention()

    det = IO.detect(c, snap)
    @test det isa IO.Detection
    # A UUID declaration outranks any fingerprint.
    @test det.rank === :declared

    desc = IO.decode(c, snap, det)
    @test desc.gridname == "igeo7"
    @test desc.system == DGG.IGeo7System()
    @test desc.idscheme === :z7int
    @test desc.level == 10
    @test desc.encoding === DENSE
    @test desc.coordinate == "cell_ids"
    # `spatial_dimension` names the data dimension, `coordinate` the array.
    @test desc.spatial_dimension == "cell_ids"
    @test desc.variables == ["elevation", "slope_geodesic", "slope_lookup"]
    @test desc.ellipsoid.semi_major_axis == 6378137.0
    @test desc.ellipsoid.inverse_flattening == 298.257223563
    @test desc.orientation.vert0_lon == 11.2
    @test desc.orientation.vert0_lat == 58.28252559
    @test desc.orientation.rotation_pattern == "alternating_cw_odd_ccw_even"
end

@testset "convention A ignores the name field of the declaration" begin
    # The meta-spec: "The `name` MUST NOT be used by tools to identify the
    # Convention". A declaration with the right uuid under a different name is
    # still this convention; one with a foreign uuid is not.
    snap = pori_dense_snapshot()
    snap.attrs["zarr_conventions"][1]["name"] = "not-dggs"
    @test IO.detect(IO.ZarrDGGSConvention(), snap) isa IO.Detection

    other = pori_dense_snapshot()
    other.attrs["zarr_conventions"] = Any[d("name" => "dggs",
        "uuid" => "3e22156d-ea9e-4e01-95fe-e3809a4b41e7")]
    @test IO.detect(IO.ZarrDGGSConvention(), other) === nothing
end

@testset "convention A reads the ranges encoding the flat attrs cannot express" begin
    snap = pori_ranges_snapshot()
    c = IO.ZarrDGGSConvention()
    desc = IO.decode(c, snap, IO.detect(c, snap))
    @test desc.encoding === RANGES
    @test desc.coordinate == "cell_id_ranges"
    @test desc.spatial_dimension == "cell_ids"
    @test desc.variables == ["elevation"]
end

# Section 4.5: `data-taos.ifremer.fr/EGU25_CFOSAT/Sentinel2_test.zarr` -- root
# `.zattrs` is literally `{}`, and the grid is described only by the alias trio
# on the coordinate.
function ifremer_snapshot(; nside=2097152, coord_extra=())
    coord = d("_ARRAY_DIMENSIONS" => Any["cell_ids"], "grid_name" => "healpix",
        "nest" => true, "nside" => nside)
    for (k, v) in coord_extra
        coord[k] = v
    end
    arrays = [
        IO.ArrayEntry(name="cell_ids", attrs=coord, shape=(4194304,),
            eltype=Int64, dims=["cell_ids"]),
        IO.ArrayEntry(name="Sentinel2",
            attrs=d("_ARRAY_DIMENSIONS" => Any["time", "bands", "cell_ids"]),
            shape=(88, 4, 4194304), eltype=Float64,
            dims=["time", "bands", "cell_ids"]),
        IO.ArrayEntry(name="bands", attrs=d("_ARRAY_DIMENSIONS" => Any["bands"]),
            shape=(4,), eltype=String, dims=["bands"]),
        IO.ArrayEntry(name="time", attrs=d("_ARRAY_DIMENSIONS" => Any["time"]),
            shape=(88,), eltype=Int64, dims=["time"])]
    return IO.StoreSnapshot(
        identifier="https://data-taos.ifremer.fr/EGU25_CFOSAT/Sentinel2_test.zarr",
        attrs=d(), arrays=arrays)
end

# Section 2.5: nextGEMS `ngc4008_P1D_7.zarr`. No cell array at all -- the array
# index is the nested index, and `12 * 4^7 == 196608`.
function dkrz_snapshot(; cell_array=false)
    arrays = IO.ArrayEntry[
        IO.ArrayEntry(name="crs",
            attrs=d("_ARRAY_DIMENSIONS" => Any["crs"], "grid_mapping_name" => "healpix",
                "healpix_nside" => 128, "healpix_order" => "nest"),
            shape=(1,), eltype=Int8, dims=["crs"]),
        IO.ArrayEntry(name="tas",
            attrs=d("_ARRAY_DIMENSIONS" => Any["time", "cell"], "grid_mapping" => "crs"),
            shape=(365, 196608), eltype=Float32, dims=["time", "cell"]),
        IO.ArrayEntry(name="time", attrs=d("_ARRAY_DIMENSIONS" => Any["time"]),
            shape=(365,), eltype=Int64, dims=["time"])]
    cell_array && push!(arrays, IO.ArrayEntry(name="cell",
        attrs=d("_ARRAY_DIMENSIONS" => Any["cell"]), shape=(1000,),
        eltype=Int64, dims=["cell"]))
    return IO.StoreSnapshot(
        identifier="https://data.nextgems-h2020.eu/ngc4008_P1D_7.zarr",
        attrs=d(), arrays=arrays)
end

# The Pori store stripped of its group attrs: convention B alone, as on every
# store whose root `.zattrs` is `{}`.
function pori_xdggs_snapshot(; grid_name="igeo7")
    snap = pori_dense_snapshot()
    empty!(snap.attrs)
    IO.getarray(snap, "cell_ids").attrs["grid_name"] = grid_name
    return snap
end

# The ranges store stripped of its group attrs: convention B alone, on the
# `(n, 2)` array it has no vocabulary for.
function pori_ranges_xdggs_snapshot()
    snap = pori_ranges_snapshot()
    empty!(snap.attrs)
    return snap
end

# Convention A on a HEALPix store. `indexing_scheme` is optional there and the
# grid name does not pin the id packing, so A alone may be unable to say it.
function healpix_a_snapshot(; scheme=nothing, coord_extra=())
    dggs = d("name" => "healpix", "refinement_level" => 5,
        "spatial_dimension" => "cells", "coordinate" => "cell_ids",
        "compression" => "none")
    scheme === nothing || (dggs["indexing_scheme"] = scheme)
    coord = d("_ARRAY_DIMENSIONS" => Any["cells"])
    for (k, v) in coord_extra
        coord[k] = v
    end
    arrays = [
        IO.ArrayEntry(name="cell_ids", attrs=coord, shape=(1000,), eltype=Int64,
            dims=["cells"]),
        IO.ArrayEntry(name="sst", attrs=d("_ARRAY_DIMENSIONS" => Any["cells"]),
            shape=(1000,), eltype=Float32, dims=["cells"])]
    return IO.StoreSnapshot(identifier="healpix_a.zarr",
        attrs=d("dggs" => dggs,
            "zarr_conventions" => Any[deepcopy(IO.ZARR_DGGS_DECLARATION)]),
        arrays=arrays)
end

# The same store with every attribute removed, for the write path. Dimension
# names live on the snapshot rather than in attrs, so a blanked snapshot is
# still structurally complete.
function blanked(s::IO.StoreSnapshot)
    out = copy(s)
    empty!(out.attrs)
    for a in out.arrays
        empty!(a.attrs)
    end
    return out
end

# ---------------------------------------------------------------------------
# Convention B -- xdggs flat coordinate attrs
# ---------------------------------------------------------------------------

@testset "convention B decodes the flat coordinate attrs" begin
    snap = pori_xdggs_snapshot()
    c = IO.XdggsConvention()

    det = IO.detect(c, snap)
    @test det isa IO.Detection
    @test det.rank === :fingerprint

    desc = IO.decode(c, snap, det)
    @test desc.gridname == "igeo7"
    @test desc.level == 10
    @test desc.encoding === DENSE
    @test desc.coordinate == "cell_ids"
    @test desc.spatial_dimension == "cell_ids"
    @test desc.variables == ["elevation", "slope_geodesic", "slope_lookup"]
    # The plugin's two extras are grid-defining and must survive decoding.
    @test desc.orientation.vert0_lon == 11.2
    @test desc.geodetic_conversion === true
end

@testset "convention B claims no encoding on a two-dimensional coordinate" begin
    # The ranges stores carry B's flat attrs on an `(n, 2)` array. B has no
    # vocabulary for that layout, so it must leave the encoding and the spatial
    # dimension to convention A rather than assert a dense axis.
    snap = pori_ranges_snapshot()
    c = IO.XdggsConvention()
    desc = IO.decode(c, snap, IO.detect(c, snap))
    @test desc.coordinate == "cell_id_ranges"
    @test desc.encoding === nothing
    @test desc.spatial_dimension === nothing
    @test desc.variables === nothing
end

@testset "the level aliases resolve, and only one may be given" begin
    snap = ifremer_snapshot()
    desc = IO.describe_store(snap)
    @test desc.gridname == "healpix"
    @test desc.system == DGG.HEALPixSystem()
    # nside = 2^21, and reading it as a level would be a silent 2^level error.
    @test desc.level == 21
    @test desc.idscheme === :nested
    @test desc.encoding === DENSE
    @test desc.variables == ["Sentinel2"]

    bad = ifremer_snapshot(nside=1000)
    err = try
        IO.describe_store(bad)
    catch e
        e
    end
    @test err isa IO.DGGSFormatError && err.check === :invalid_nside

    both = ifremer_snapshot(coord_extra=("level" => 21,))
    err2 = try
        IO.describe_store(both)
    catch e
        e
    end
    @test err2 isa IO.DGGSFormatError && err2.check === :duplicate_level_alias
end

@testset "the legacy healpix trio is its own convention and agrees with B" begin
    snap = ifremer_snapshot()
    fired = [IO.conventionname(det.convention)
             for det in IO.detections(snap)]
    @test "xdggs" in fired
    @test "legacy-healpix" in fired
    # Both read the same trio, so the merge must be silent.
    @test IO.describe_store(snap).level == 21
end

# ---------------------------------------------------------------------------
# Merging
# ---------------------------------------------------------------------------

@testset "a dual-stamped store merges field for field" begin
    desc = IO.describe_store(pori_dense_snapshot())
    @test desc.gridname == "igeo7"
    @test desc.level == 10
    @test desc.encoding === DENSE
    # A carries the full orientation, B only the longitude: the two merge.
    @test desc.orientation ==
          IO.GridOrientation(vert0_lon=11.2, vert0_lat=58.28252559,
        vert0_azimuth=0.0, rotation_pattern="alternating_cw_odd_ccw_even")
    @test desc.geodetic_conversion === true
    @test desc.ellipsoid.semi_major_axis == 6378137.0
    @test sort!(collect(keys(desc.provenance))) ==
          ["xdggs", "zarr-conventions/dggs"]

    ranges = IO.describe_store(pori_ranges_snapshot())
    @test ranges.encoding === RANGES
    @test ranges.coordinate == "cell_id_ranges"
    @test ranges.spatial_dimension == "cell_ids"
    @test ranges.variables == ["elevation"]
    @test ranges.geodetic_conversion === true
end

@testset "conventions that contradict each other are refused" begin
    # The store says level 10 in its group attrs and level 12 on its
    # coordinate. Guessing which is right is exactly what we refuse to do.
    snap = pori_dense_snapshot(level=10, coord_level=12)
    err = try
        IO.describe_store(snap)
    catch e
        e
    end
    @test err isa IO.DGGSFormatError
    @test err.check === :level_disagreement
    @test Set([err.declared, err.observed]) == Set([10, 12])
    @test occursin("xdggs", sprint(showerror, err))
    @test occursin("zarr-conventions/dggs", sprint(showerror, err))
end

@testset "a description no convention could complete is refused" begin
    # B's flat attrs on the `(n, 2)` coordinate identify the grid and nothing
    # about the layout. Without A there is no encoding and no spatial
    # dimension, and a description that cannot be opened must say so here
    # rather than fail later as a shape error.
    err = try
        IO.describe_store(pori_ranges_xdggs_snapshot())
    catch e
        e
    end
    @test err isa IO.DGGSFormatError && err.check === :incomplete_description
    msg = sprint(showerror, err)
    @test occursin("xdggs", msg)
    @test occursin("encoding", msg) && occursin("spatial_dimension", msg)
end

@testset "an attribute of the wrong type is refused, not ignored" begin
    # A non-string `spatial_dimension` used to decode as `nothing`, which reads
    # as "the store did not say" -- the one thing it does not mean.
    snap = pori_dense_snapshot()
    snap.attrs["dggs"]["spatial_dimension"] = 3
    err = try
        IO.describe_store(snap)
    catch e
        e
    end
    @test err isa IO.DGGSFormatError && err.check === :invalid_attribute_type
end

@testset "a store with no DGGS metadata is not a DGGS store" begin
    # NICAM `data_healpix_{z}.zarr`: root attrs `{}` and no grid metadata
    # anywhere -- the level exists only in the filename.
    snap = IO.StoreSnapshot(identifier="data_healpix_7.zarr", attrs=d(),
        arrays=[IO.ArrayEntry(name="ta", shape=(196608,), dims=["cell"])])
    err = try
        IO.describe_store(snap)
    catch e
        e
    end
    @test err isa IO.DGGSFormatError && err.check === :no_convention_detected
end

# ---------------------------------------------------------------------------
# The grid name table and its hook
# ---------------------------------------------------------------------------

@testset "an unregistered grid name is refused, with the registry listed" begin
    snap = pori_xdggs_snapshot(grid_name="ISEA7H_Z7")
    err = try
        IO.describe_store(snap)
    catch e
        e
    end
    @test err isa IO.DGGSFormatError && err.check === :unknown_grid_name
    msg = sprint(showerror, err)
    @test occursin("igeo7", msg) && occursin("healpix", msg)
end

struct AliasConvention <: IO.DGGSConvention end
IO.conventionname(::AliasConvention) = "alias"
IO.gridname(::AliasConvention, name, attrs) =
    name == "ISEA7H_Z7" ? "igeo7" : lowercase(String(name))
function IO.detect(c::AliasConvention, s::IO.StoreSnapshot)
    det = IO.detect(IO.XdggsConvention(), s)
    return det === nothing ? nothing : IO.Detection(c, det.rank, det.evidence;
        payload=det.payload)
end
IO.decode(c::AliasConvention, s::IO.StoreSnapshot, det::IO.Detection) =
    IO.decode_flat_attrs(c, s, det)

@testset "a convention may fold its own grid vocabulary in" begin
    snap = pori_xdggs_snapshot(grid_name="ISEA7H_Z7")
    desc = IO.describe_store(snap; conventions=[AliasConvention()])
    @test desc.gridname == "igeo7"
    @test desc.system == DGG.IGeo7System()
    @test desc.level == 10
    # Filled by `applydefaults` from the table, not by the store.
    @test desc.idscheme === :z7int
end

# ---------------------------------------------------------------------------
# The DKRZ dialect
# ---------------------------------------------------------------------------

@testset "DKRZ decodes an implicit axis from the crs variable" begin
    snap = dkrz_snapshot()
    c = IO.DKRZConvention()
    det = IO.detect(c, snap)
    @test det isa IO.Detection

    desc = IO.decode(c, snap, det)
    @test desc.gridname == "healpix"
    # healpix_nside is an nside, where CF's refinement_level is an order.
    @test desc.level == 7
    @test desc.idscheme === :nested       # the legacy spelling is "nest"
    @test desc.coordinate === nothing
    @test desc.encoding === IMPLICIT
    @test desc.spatial_dimension == "cell"
    @test desc.variables == ["tas"]

    # The presence of a cell array is the sparsity signal.
    regional = IO.describe_store(dkrz_snapshot(cell_array=true))
    @test regional.coordinate == "cell"
    @test regional.encoding === DENSE
end

# ---------------------------------------------------------------------------
# Incomplete declarations and their defaults
# ---------------------------------------------------------------------------

@testset "an incomplete convention A stamp is completed, not refused" begin
    # A says which grid and which layout; it need not say how the ids are
    # packed. Refusing there would make a store unreadable that B, or the
    # convention default, describes perfectly well.
    merged = IO.describe_store(healpix_a_snapshot(
        coord_extra=("grid_name" => "healpix", "level" => 5,
            "indexing_scheme" => "ring")))
    @test merged.idscheme === :ring
    @test merged.encoding === DENSE
    @test merged.level == 5
    @test merged.variables == ["sst"]

    # A alone, with nothing to merge: the table default fills it in.
    alone = IO.describe_store(healpix_a_snapshot())
    @test alone.idscheme === :nested
    @test alone.encoding === DENSE
end

@testset "a declared indexing scheme survives the defaults" begin
    # `ring` is not the table default, and reading it as `nested` would
    # misplace every cell of the store.
    @test IO.describe_store(healpix_a_snapshot(scheme="ring")).idscheme === :ring
    @test IO.describe_store(healpix_a_snapshot(scheme="nested")).idscheme === :nested
end

# ---------------------------------------------------------------------------
# Ellipsoid spellings
# ---------------------------------------------------------------------------

@testset "a store that declares no ellipsoid gets the prescribed one" begin
    # The DKRZ crs variable carries no ellipsoid keys at all. Decoding that as
    # an all-`nothing` Ellipsoid would be a declaration, and would shadow the
    # default the convention prescribes.
    desc = IO.describe_store(dkrz_snapshot())
    @test desc.ellipsoid === nothing
    @test IO.ellipsoid(desc) === IO.DEFAULT_ELLIPSOID
end

@testset "both ellipsoid spellings are accepted on read" begin
    wild = pori_dense_snapshot()
    valid = pori_dense_snapshot()
    e = valid.attrs["dggs"]["ellipsoid"]
    e["semi_major_axis"] = pop!(e, "semimajor_axis")
    @test IO.describe_store(wild).ellipsoid == IO.describe_store(valid).ellipsoid
    @test IO.describe_store(valid).ellipsoid.semi_major_axis == 6378137.0
end

# ---------------------------------------------------------------------------
# The write path
# ---------------------------------------------------------------------------

@testset "the default write conventions dual-stamp a schema-valid store" begin
    source = pori_dense_snapshot()
    desc = IO.describe_store(source)

    target = blanked(source)
    for c in IO.DEFAULT_WRITE_CONVENTIONS
        IO.encode!(c, target, desc)
    end

    dggs = target.attrs["dggs"]
    @test dggs["name"] == "igeo7"
    @test dggs["refinement_level"] == 10
    @test dggs["spatial_dimension"] == "cell_ids"
    @test dggs["coordinate"] == "cell_ids"
    @test dggs["compression"] == "none"
    # The schema-valid spelling, against the wild stores' invalid one.
    @test haskey(dggs["ellipsoid"], "semi_major_axis")
    @test !haskey(dggs["ellipsoid"], "semimajor_axis")
    @test any(e -> e["uuid"] == IO.ZARR_DGGS_UUID, target.attrs["zarr_conventions"])

    # xdggs forwards every non-`grid_name` attr as a dataclass keyword, so an
    # extra key on the coordinate is a TypeError in the reader. This is exactly
    # `IGEO7Info.to_dict()`.
    @test Set(keys(IO.getarray(target, "cell_ids").attrs)) ==
          Set(["grid_name", "level", "igeo7_wgs84_geodetic_conversion",
        "igeo7_dggs_vert0_lon"])

    # And the stamp is a fixpoint of detection.
    @test IO.describe_store(target) == desc
end

@testset "stamping is a fixpoint for the computed encodings too" begin
    # The dense case is above. A ranges store's layout is expressible only by
    # convention A, and an implicit one names no coordinate for B to stamp at
    # all: both must still come back as the description they went in as.
    for source in (pori_ranges_snapshot(), dkrz_snapshot())
        desc = IO.describe_store(source)
        target = blanked(source)
        for c in IO.DEFAULT_WRITE_CONVENTIONS
            IO.encode!(c, target, desc)
        end
        @test IO.describe_store(target) == desc
    end
end

@testset "a downstream convention joins the registry" begin
    n = length(IO.CONVENTION_REGISTRY)
    try
        IO.register_convention!(AliasConvention())
        @test last(IO.CONVENTION_REGISTRY) isa AliasConvention
        # Ahead of everything, which is how a store's own dialect outranks the
        # shipped fingerprints.
        IO.register_convention!(AliasConvention(); first=true)
        @test first(IO.CONVENTION_REGISTRY) isa AliasConvention
        @test length(IO.CONVENTION_REGISTRY) == n + 2
    finally
        filter!(c -> !(c isa AliasConvention), IO.CONVENTION_REGISTRY)
    end
    @test length(IO.CONVENTION_REGISTRY) == n
end

@testset "a downstream grid name joins the reference table" begin
    # A grid name pins the id packing, so an unregistered one is refused rather
    # than guessed — and the refusal has to point at the way in. Kills a
    # reference table that only the package itself can extend.
    err = try
        IO.gridreference("isea7h")
    catch e
        e
    end
    @test err isa IO.DGGSFormatError && err.check === :unknown_grid_name
    @test occursin("register_grid!", sprint(showerror, err))

    ref = IO.GridReference("isea7h", DGG.IGeo7System(), :z7int, (:z7int,))
    try
        @test DGG.register_grid!("isea7h", ref) === IO.GRID_REFERENCE
        @test IO.gridreference("isea7h") === ref
    finally
        delete!(IO.GRID_REFERENCE, "isea7h")
    end
    @test !haskey(IO.GRID_REFERENCE, "isea7h")
end

@testset "a read-only convention refuses to stamp" begin
    desc = IO.describe_store(pori_dense_snapshot())
    @test_throws ArgumentError IO.encode!(IO.DKRZConvention(),
        blanked(pori_dense_snapshot()), desc)
end

# ---------------------------------------------------------------------------
# The stubs
# ---------------------------------------------------------------------------

@testset "dggread and dggwrite ask for Zarr, or the extension has them" begin
    # The stub message is only true before the extension loads, and
    # `test/io/runtests.jl` includes this file first so that it usually is. When
    # Zarr came in ahead of the suite the claim to check is the other one: the
    # stub has been taken over rather than left in place.
    ext = Base.get_extension(parentmodule(IO.dggread), :DiscreteGlobalGridsZarrExt)
    fromext(m) = m.module === ext || parentmodule(m.module) === ext
    for f in (IO.dggread, IO.dggwrite)
        if ext === nothing
            err = try
                f("store.zarr")
            catch e
                e
            end
            @test err isa ErrorException && occursin("using Zarr", err.msg)
        else
            @test any(fromext, methods(f))
        end
    end
end

end # module
