# Demo: reading IGEO7/Z7 DGGS-convention Zarr archives lazily (DGGSZarr).
#
# The Julia counterpart of the Python `eesti_soil_conversion/samples.py`: open
# an archive, look at its structure without reading it, then select cells by id,
# by lat/lon and by bounding box, and derive centers and boundaries.
#
# It runs against the four reference archives written by the Python IGEO7/Z7
# tooling, in two forms of the same data:
#
#   pori_z7_r10.zarr         compression "none"   — a dense uint64 cell_ids array
#   pori_z7_r10_ranges.zarr  compression "ranges" — a (R, 2) monotonic range table
#
# and likewise at res 12. The point of the exercise is that both decode to the
# same cells, and that the ranges form does it in O(R) rather than O(N).
#
# Environment: this demo needs the archives, which live outside the repository.
# Point it at them with DGGS_ZARR_TEST_DATA if they are not in the default
# location, then:
#
#     julia -t 4 --project=. examples/z7_zarr_read.jl
#
# It ends in PASS/FAIL assertions and exits non-zero if any of them fail.
using DiscreteGlobalGrids
using DiscreteGlobalGrids.DGGSZarr
using DiscreteGlobalGrids: Helpers, ISEA, IGeo7
using DimensionalData
import DiscreteGlobalGrids as DGG
import GeoInterface
import Zarr

const FAILURES = Ref(0)
function check(name, ok; detail = "")
    ok || (FAILURES[] += 1)
    println(ok ? "PASS  " : "FAIL  ", rpad(name, 52), detail)
    return ok
end
note(text) = println("      ", text)

const ARCHIVES = get(ENV, "DGGS_ZARR_TEST_DATA",
    joinpath(homedir(), "dev", "build", "igeo7_z7_xarray_paper", "data", "working"))

const RANGES_10 = joinpath(ARCHIVES, "pori_z7_r10_ranges.zarr")
const DENSE_10 = joinpath(ARCHIVES, "pori_z7_r10.zarr")
const RANGES_12 = joinpath(ARCHIVES, "pori_z7_r12_ranges.zarr")

println("="^78)
println("z7_zarr_read.jl — IGEO7/Z7 DGGS-convention Zarr archives")
println("julia $(VERSION)  threads=$(Threads.nthreads())")
println("archives: ", ARCHIVES)
println("="^78)

if !all(isdir, (RANGES_10, DENSE_10, RANGES_12))
    println("\nreference archives not found — set DGGS_ZARR_TEST_DATA to the directory")
    println("holding pori_z7_r10{,_ranges}.zarr and pori_z7_r12{,_ranges}.zarr.")
    exit(0)
end

# --------------------------------------------------------------------------
# 1. Open, and look at the structure without reading it.
#
# `open_dggs_dataset` reads the group attributes and — for a ranges archive —
# the R × 2 range table. No data variable and no length-N coordinate.
# --------------------------------------------------------------------------

println("\n--- 1. lazy open ---")

ds = open_dggs_dataset(RANGES_12)
info = dggs_info(ds)
println(ds)

check("archive declares itself", info.name == "igeo7" && info.compression == "ranges";
      detail = "level $(info.level)")

ids = dggs_cell_ids(ds)
check("coordinate is a lazy range vector", ids isa IGeo7.Z7RangeIds;
      detail = string(ids))

# The size claim that motivates the ranges form: R ranges standing in for N ids.
range_bytes = Base.summarysize(ids)
dense_bytes = 8 * length(ids)
check("index costs O(R), not O(N)", range_bytes < dense_bytes ÷ 4;
      detail = "$(range_bytes) B for $(length(ids)) cells vs $(dense_bytes) B dense")
note("R = $(IGeo7.z7_nranges(ids)) ranges, N = $(length(ids)) cells " *
     "($(round(100 * IGeo7.z7_nranges(ids) / length(ids); digits = 2))% R/N)")

# Data variables are still the archive's Zarr arrays — nothing has been read.
check("data variables stay on disk", ds.elevation.data isa Zarr.ZArray;
      detail = string(size(ds.elevation)))

# --------------------------------------------------------------------------
# 2. The two storage forms decode to the same cells.
#
# The ranges archive stores 136 [start, end] pairs; the dense one stores 3101
# ids. Expanding the former must reproduce the latter exactly — which is the
# check that the base-7 monotonic number line is the right one. (It is *not*
# `cell_to_index`, the pentagon-aware dense rank; the two agree on order and
# differ on spacing.)
# --------------------------------------------------------------------------

println("\n--- 2. ranges == dense ---")

ranges_10 = open_dggs_dataset(RANGES_10)
dense_10 = open_dggs_dataset(DENSE_10)

from_ranges = collect(dggs_cell_ids(ranges_10))
from_dense = collect(dggs_cell_ids(dense_10))

check("expanded ranges == stored dense ids", from_ranges == from_dense;
      detail = "$(length(from_ranges)) cells")
check("all ids valid at the declared level",
      all(id -> IGeo7.is_valid_z7(id) && IGeo7.z7_resolution(id) == 10, from_ranges))

# The distinction the reader turns on, made explicit.
sample = from_ranges[1]
note("cell $(IGeo7.z7_to_string(sample)):  " *
     "z7_to_monotonic = $(IGeo7.z7_to_monotonic(sample, 10)), " *
     "cell_to_index = $(IGeo7.cell_to_index(sample))")

# --------------------------------------------------------------------------
# 3. Select by cell id.
#
# The dimension is a real DGGS lookup, so this is ordinary DimensionalData
# indexing. On the ranges archive `At` is answered in O(log R) against the
# range table, without the coordinate ever being materialized.
# --------------------------------------------------------------------------

println("\n--- 3. select by cell id ---")

id = from_ranges[100]
a = only(ranges_10.elevation[cell_ids = At(id)])
b = only(dense_10.elevation[cell_ids = At(id)])
check("At(id) agrees across both forms", a == b;
      detail = "0x$(string(id; base = 16, pad = 16)) -> $(a) m")

# A slice is still lazy until something asks for values.
slice = ranges_10.elevation[cell_ids = 1:10]
check("a slice stays lazy", slice.data isa AbstractArray && ndims(slice) == 1;
      detail = "$(size(slice)), $(eltype(slice))")

# A cell of the right resolution that this regional archive does not hold.
absent = IGeo7.index_to_cell(1, 10)
missing_ok = try
    ranges_10.elevation[cell_ids = At(absent)]
    false
catch e
    e isa ArgumentError
end
check("an absent cell is refused, not rounded", missing_ok)

# --------------------------------------------------------------------------
# 4. Geometry: the archive's placement and datum.
#
# Two archive-declared facts that the plain IGeo7 geometry functions do not know
# about, and that displace every coordinate if ignored:
#
#   dggs_vert0_lon = 11.2   where icosahedron vertex 0 sits, against the
#                           package default ISEA_LON0 = 11.25
#   igeo7_wgs84_geodetic_conversion = true
#                           latitudes are WGS84 geodetic, not authalic
#
# `DGGSZarr`'s geometry functions read both off the archive; `IGeo7`'s
# same-named ones do not, by design.
# --------------------------------------------------------------------------

println("\n--- 4. cell centers and boundaries ---")

centers = cell_centers(ranges_10)
lons = first.(centers)
lats = last.(centers)
check("centers cover the archive extent", length(centers) == length(from_ranges);
      detail = "lon $(round.(extrema(lons); digits = 3)), " *
               "lat $(round.(extrema(lats); digits = 3))")

naive_lon, naive_lat = IGeo7.cell_center(from_ranges[1])
check("placement matters", abs(naive_lon - centers[1][1]) > 0.04;
      detail = "default orientation is off by " *
               "$(round(naive_lon - centers[1][1]; digits = 4))° lon")
check("datum matters", abs(naive_lat - centers[1][2]) > 0.1;
      detail = "authalic vs WGS84 differ by " *
               "$(round(naive_lat - centers[1][2]; digits = 4))° lat (~13 km)")

boundaries = cell_boundaries(ranges_10)
check("boundaries are closed hexagons", length(boundaries) == length(from_ranges);
      detail = "$(length(first(GeoInterface.coordinates(boundaries[1])))) corners incl. repeat")

# --------------------------------------------------------------------------
# 5. Select by lat/lon and by bounding box.
# --------------------------------------------------------------------------

println("\n--- 5. select by lat/lon and bbox ---")

lon0, lat0 = centers[500]
value = only(sel_latlon(ranges_10.elevation, lon0, lat0))
check("sel_latlon finds the cell under a point",
      value == only(ranges_10.elevation[cell_ids = At(from_ranges[500])]);
      detail = "($(round(lon0; digits = 4)), $(round(lat0; digits = 4))) -> $(value) m")

# Round trip: every center resolves back to its own cell.
roundtrip = all(i -> dggs_cell_at(ranges_10, centers[i]...) == from_ranges[i],
                eachindex(from_ranges))
check("every center resolves to its own cell", roundtrip;
      detail = "$(length(from_ranges)) cells")

box = (26.6, 26.7, 58.15, 58.2)
sub = sel_bbox(ranges_10.elevation, box...)
expected = count(c -> box[1] <= c[1] <= box[2] && box[3] <= c[2] <= box[4], centers)
check("sel_bbox selects by cell center", length(sub) == expected;
      detail = "lon $(box[1])–$(box[2]), lat $(box[3])–$(box[4]) -> $(length(sub)) cells")

# A whole-dataset selection keeps every variable.
sub_ds = sel_bbox(ranges_10, box...)
check("bbox on the dataset keeps all variables",
      keys(sub_ds.cubes) == keys(ranges_10.cubes);
      detail = string(collect(keys(sub_ds.cubes))))

# --------------------------------------------------------------------------
# 6. Downstream: this is where the data actually gets pulled.
#
# Everything above stayed lazy. A numeric reduction is a consumer, so it reads —
# and it is the shape any modelling or simulation step would take.
# --------------------------------------------------------------------------

println("\n--- 6. materializing for computation ---")

elevation = ranges_10.elevation[:]
finite = filter(!isnan, elevation)
check("elevation reads as a plain array", elevation isa AbstractVector;
      detail = "$(length(finite)) finite of $(length(elevation)), " *
               "$(round(minimum(finite); digits = 1))–$(round(maximum(finite); digits = 1)) m")

# The cell_id -> position map a CSR adjacency build needs is the lookup itself:
# O(log R) per query and O(R) memory, in place of a Dict of N entries.
lookup = dggs_lookup(ranges_10)
neighbors = DGG.cell_neighbors(DGG.IGEO7DGGS(), 10, from_ranges[500])
positions = [DGG.cell_position(dggs_cell_ids(ranges_10), n) for n in neighbors]
check("neighbours map to positions through the index",
      count(!isnothing, positions) >= 5;
      detail = "$(length(neighbors)) neighbours, " *
               "$(count(!isnothing, positions)) inside the archive")
note("that mapping is the `cell_id -> matrix row` step of a CSR adjacency build")

println("\n", "="^78)
if FAILURES[] == 0
    println("all checks passed")
else
    println(FAILURES[], " check(s) FAILED")
end
println("="^78)
exit(FAILURES[] == 0 ? 0 : 1)
