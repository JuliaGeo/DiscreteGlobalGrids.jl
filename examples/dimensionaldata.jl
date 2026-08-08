# Demo: DimensionalData lookups — from a stored cell axis to a spatial tree.
#
# A vector data cube stores one sparse coverage as a `DimArray` whose dimension
# is a `<X>Lookup` — a sorted vector of cell ids. Everything downstream (a
# `Regridder`, a cap query) has to address the *same* positions, or the answers
# silently transpose. The generic partial grid guarantees that structurally: it
# stores `lookup.data` by reference and never reorders it, so tree leaf i is
# lookup position i is `parent(A)[i]`.
#
# Environment: this demo needs nothing beyond DiscreteGlobalGrids and its own
# dependencies, so it runs in the repository's project environment:
#
#     julia -t 4 --project=. examples/dimensionaldata.jl
#
# It ends in PASS/FAIL assertions and exits non-zero if any of them fail.
#
# Note the imports: the two `<X>Lookups` modules share generic vocabulary
# (`cell_center`, `cell_polygons`, `Touching`, ...), so a script that wants both
# systems at once names what it needs instead of `using` both wholesale.
using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
using DiscreteGlobalGrids.HEALPix.HealpixLookups: HealpixLookup
using DiscreteGlobalGrids.IGeo7.IGeo7Lookups: IGeo7Lookup
const IGeo7 = DGG.IGeo7
import DimensionalData as DD
import Healpix
import ConservativeRegridding as CR
import GeometryOps as GO
import GeometryOps: SpatialTreeInterface as STI
using Random

const FAILURES = Ref(0)
function check(name, ok; detail = "")
    ok || (FAILURES[] += 1)
    println(ok ? "PASS  " : "FAIL  ", rpad(name, 56), detail)
    return ok
end
note(text) = println("      ", text)

const SEED = 20260805
const NCELLS = 3000

# Fixed 10-degree lon/lat destination, as unit-sphere corner points.
const TO_SPHERE = GO.UnitSpherical.UnitSphereFromGeographic()
const DST_LON = collect(range(0, 360; length = 37))
const DST_LAT = collect(range(-90, 90; length = 19))
const DST = [TO_SPHERE((x, y)) for x in DST_LON, y in DST_LAT]
const DST_CELLS = vec(collect(getcell(treeify(DST))))
const INTERSECT = CR.DefaultIntersectionOperator(GO.Spherical())

println("="^78)
println("dimensionaldata.jl — Lookup -> DGGSPartialGrid -> tree, $NCELLS cells per system")
println("julia $(VERSION)  threads=$(Threads.nthreads())  seed=$SEED")
println("="^78)

"`n` distinct ids drawn from a seeded uniform-on-sphere point stream, ascending."
function sampled_ids(encode, T, n; seed = SEED)
    rng = Xoshiro(seed)
    ids = Set{T}()
    while length(ids) < n
        push!(ids, T(encode(360 * rand(rng) - 180, asind(2 * rand(rng) - 1))))
    end
    return sort!(collect(ids))
end

# --------------------------------------------------------------------------
# The generic half of the demo: everything below is system-agnostic. It takes
# the lookup plus the DGGS singleton and level the lookup describes, and never
# names a system-specific tree, grid or polygon helper.
# --------------------------------------------------------------------------

function exercise(label, lookup, system, level)
    # THE CALL SITE — one line from a stored cell axis to a spatial tree:
    tree = treeify(DGGSPartialGrid(lookup))

    grid = DGGSPartialGrid(lookup)
    check("$label: ncells == length(lookup.data)",
          ncells(tree) == length(lookup.data);
          detail = "$(ncells(tree)) cells")
    check("$label: grid.ids === lookup.data (no copy, no reorder)",
          grid.ids === lookup.data)

    sample = (1, 2, 977, NCELLS ÷ 2, NCELLS)
    polygon(id) = cell_polygon_unitsphere(system, level, id)
    check("$label: getcell(tree, i) is lookup position i's cell",
          all(getcell(tree, i) == polygon(lookup.data[i]) for i in sample);
          detail = "sampled i = $(sample)")

    # A DimArray over the same lookup: `parent(A)[i]` is tree leaf i, so a field
    # vector can be handed to `regrid!` with no permutation anywhere.
    array = DD.DimArray(collect(1.0:length(lookup)), (DD.Dim{:cell}(lookup),))
    check("$label: parent(DimArray)[i] pairs with tree leaf i",
          all(parent(array)[i] == Float64(i) &&
              getcell(tree, i) == polygon(lookup.data[i]) for i in sample))

    # Regridder: CR's matrix is (dst cells) x (src cells), so with the DGGS on
    # the source side the lookup axis is the COLUMNS.
    regridder = CR.Regridder(DST, DGGSPartialGrid(lookup))
    A = regridder.intersections
    check("$label: matrix is (dst, lookup positions)",
          size(A) == (length(DST_CELLS), length(lookup.data));
          detail = "$(size(A)), nnz=$(length(A.nzval))")

    # Spot-check: the column with the widest footprint, against the intersection
    # areas computed directly from lookup.data[j]'s polygon — no tree involved.
    j = argmax(diff(A.colptr))
    direct = [INTERSECT(cell, polygon(lookup.data[j])) for cell in DST_CELLS]
    column = Vector(A[:, j])
    coldiff = maximum(abs.(direct .- column) ./ max.(abs.(direct), 1.0))
    check("$label: column $j == direct intersection areas", coldiff <= 1e-12;
          detail = "$(count(!iszero, column)) nonzero dst cells, max rel diff $coldiff")

    # And with the DGGS on the destination side the lookup axis is the ROWS.
    reversed = CR.Regridder(DGGSPartialGrid(lookup), DST)
    R = reversed.intersections
    rowdiff = maximum(abs.(Vector(R[j, :]) .- direct) ./ max.(abs.(direct), 1.0))
    check("$label: reversed regridder row $j == same areas",
          size(R) == (length(lookup.data), length(DST_CELLS)) && rowdiff <= 1e-12;
          detail = "$(size(R)), max rel diff $rowdiff")

    # Value-level alignment: a field that is 1 at lookup position j and 0
    # elsewhere must land exactly on that cell's footprint.
    indicator = zeros(length(lookup.data))
    indicator[j] = 1.0
    destination = zeros(length(regridder.dst_areas))
    CR.regrid!(destination, regridder, indicator)
    check("$label: indicator at $j lights that footprint",
          findall(!iszero, destination) == findall(!iszero, direct);
          detail = "$(count(!iszero, destination)) dst cells lit")

    # A cap query answers in the same index space, so the hits index the
    # DimArray directly. `intersects_cap` is the whole predicate.
    cap = GO.UnitSpherical.SphericalCap(DGG.cell_center(system, level, lookup.data[j]), 0.2)
    hits = STI.query(tree, intersects_cap(cap))
    truth = [i for i in eachindex(lookup.data)
             if intersects_cap(cap, cell_cap(system, level, lookup.data[i]))]
    check("$label: cap query hits index the DimArray",
          j in hits && issubset(hits, truth) && parent(array)[j] == Float64(j);
          detail = "$(length(hits)) hits / $(length(truth)) per-cell candidates")
    return tree
end

const HEALPIX_LEVEL = 6
const HEALPIX_RES = Healpix.Resolution(2^HEALPIX_LEVEL)
healpix_lookup = HealpixLookup(sampled_ids(
        (lon, lat) -> Healpix.ang2pixNest(HEALPIX_RES, deg2rad(90 - lat), deg2rad(mod(lon, 360))) - 1,
        Int64, NCELLS); level = HEALPIX_LEVEL)
exercise("HEALPix 6", healpix_lookup, HEALPixDGGS(), HEALPIX_LEVEL)

println()
const IGEO7_LEVEL = 5
igeo7_lookup = IGeo7Lookup(sampled_ids((lon, lat) -> IGeo7.lonlat_to_cell(lon, lat, IGEO7_LEVEL),
                                       UInt64, NCELLS); resolution = IGEO7_LEVEL)
exercise("IGEO7 5", igeo7_lookup, IGEO7DGGS(), IGEO7_LEVEL)

println()
note("call site, verbatim:  tree = treeify(DGGSPartialGrid(lookup))")
note("the same line for both systems; the Lookup carries level/resolution and eltype")
note("`treeify(lookup)` is the same path in one step, straight off the DD dimension")

# ...and the one-step form really is the same tree.
short = treeify(healpix_lookup)
check("treeify(lookup) == treeify(DGGSPartialGrid(lookup))",
      ncells(short) == length(healpix_lookup.data) &&
      short.grid.ids === healpix_lookup.data &&
      (node_level(short), node_id(short)) == (-1, 0))

println()
println(FAILURES[] == 0 ? "ALL CHECKS PASSED" : "$(FAILURES[]) CHECK(S) FAILED")
exit(FAILURES[] == 0 ? 0 : 1)
