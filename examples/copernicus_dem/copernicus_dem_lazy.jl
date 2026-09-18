# All of GLO-90 as one lazy DiskArrays vector along the level-1 cell order.
#
# Run in the docs environment for ArchGDAL and DiskArrays:
#
#     julia -t auto --project=docs examples/copernicus_dem/copernicus_dem_lazy.jl
#
# One chunk per 1°x1° tile, backed by that tile's COG on AWS and fetched on
# first touch; tiles absent from the bucket's tile list read as `NaN32` with no
# network. The tile list and tiles cache under `tempdir()/copdem90`.

import DiscreteGlobalGrids as DGG
import ConservativeRegridding as CR
import GeometryOps as GO
import GeoInterface as GI
import DimensionalData as DD
import DiskArrays
import ArchGDAL
import Extents
using CopernicusUtils

const CD = DGG.CopernicusDEM

const FAILURES = Ref(0)
function check(name, ok; detail="")
    ok || (FAILURES[] += 1)
    println(ok ? "PASS  " : "FAIL  ", rpad(name, 56), detail)
    return ok
end
note(text) = println("      ", text)

const MANIFOLD = GO.Spherical(; radius=1.0)
const STARTED = time()

println("="^78)
println("copernicus_dem_lazy.jl — all of GLO-90 as one lazy chunked vector")
println("julia $(VERSION)  threads=$(Threads.nthreads())")
println("="^78)

# --- 1. The lazy vector: one chunk per tile. -----------------------------

const DATADIR = joinpath(tempdir(), "copdem90")
const CACHE_TILES = 16

"""
Every GLO-90 pixel as one `TiledDEM` over all 64 800 tiles, so its index order
is the complete level-1 grid's. `loaded` records each listed tile decoded, in
decode order.
"""
function glo90(sys, source)
    loaded = Int[]
    cache = StripedLRUCache{Vector{Float32}}(; slots=CACHE_TILES, stripes=1) do k
        (k - 1) in source.listed && push!(loaded, k - 1)
        loadtile(source, k - 1)
    end
    return TiledDEM(TileIds(sys, 0:64_799), cache), loaded
end

# --- 2. Structure: chunks are exactly the tiles' descendant ranges. ------

sys = DGG.CopernicusDEMSystem(90)
listed = listedtiles(sys, tilelist(DATADIR, 90))
source = CopernicusTiles(sys, listed; cachedir=joinpath(DATADIR, "tiles"))
lazy, loaded = glo90(sys, source)
ec = DiskArrays.eachchunk(lazy)

check("length is ncells(sys, 1)",
    length(lazy) == DGG.ncells(sys, 1) == 68_947_200_000;
    detail="$(length(lazy)) pixels, ~$(round(Int, length(lazy) * 4 / 2^30)) GiB dense")
check("chunks are exactly the tiles' descendant ranges",
    DiskArrays.haschunks(lazy) isa DiskArrays.Chunked &&
    length(ec) == 64_800 && first(ec[1][1]) == 1 && last(ec[end][1]) == length(lazy) &&
    all(ec[k][1] == DGG.descendant_range(sys, DGG.LevelIndex(0, k - 1), 1)
        for k in 1:64_800);
    detail="64800 chunks concatenate to 1:$(length(lazy))")

chunkof(t) = ec[Int(t.index)+1][1]
t49, t50 = CD.tilecell(sys, 49, 6), CD.tilecell(sys, 50, 6)
check("chunk widths step with the band table",
    length(chunkof(t49)) == 1_440_000 && length(chunkof(t50)) == 960_000 &&
    length(chunkof(CD.tilecell(sys, -90, -180))) == 144_000;
    detail="N49 1x, N50 1.5x, S90 10x — 1200/800/120 columns x 1200 rows")

# --- 3. Ocean chunks come from the manifest, not the network. ------------

ocean = CD.tilecell(sys, 0, -30)                   # mid-Atlantic
andes = CD.tilecell(sys, -34, -71)                 # land, exercises the S/W labels
check("the tile list separates land from ocean",
    length(listed) >= 26_000 &&
    Int(t50.index) in source.listed && Int(andes.index) in source.listed &&
    !(Int(ocean.index) in source.listed);
    detail="$(length(listed)) tiles listed, incl. $(tilestem(sys, andes))")
oceanvals = lazy[DGG.descendant_range(sys, ocean, 1)]
check("an ocean chunk is all NaN with zero loads",
    all(isnan, oceanvals) && isempty(loaded);
    detail="$(length(oceanvals)) pixels of $(tilestem(sys, ocean))")

# --- 4. A land chunk equals the COG read directly. -----------------------

directvec(t) = vec(ArchGDAL.read(ds -> ArchGDAL.read(ds, 1),
    tilepath!(source, Int(t.index))))

check("the N50_00_E006_00 chunk equals its COG",
    lazy[DGG.descendant_range(sys, t50, 1)] == directvec(t50) &&
    loaded == [Int(t50.index)];
    detail="960000 pixels, 1 tile loaded")

# --- 5. The cube: Covering selection loads only the tiles it touches. ----

fresh, freshloaded = glo90(sys, source)            # its own load record
dem = DD.DimArray(fresh, DGG.Cells(DGG.CellLookup(DGG.levelgrid(sys, 1))))
window = Extents.Extent(X=(6.4, 6.6), Y=(49.9, 50.1))   # crosses the 50° band edge
sub = dem[DGG.Cells(DGG.Covering(window))]
subids = collect(DD.lookup(sub, DGG.Cells))
touched = unique!([Base.parent(sys, c) for c in subids])
direct = Dict(t => directvec(t) for t in touched)
ranges = Dict(t => DGG.descendant_range(sys, t, 1) for t in touched)

check("the window touches the N49 and N50 tiles", Set(touched) == Set((t49, t50));
    detail="$(length(sub)) cells across two chunk widths")
check("selected values match the direct per-tile reads",
    all(sub[k] == direct[t][Int(subids[k].index) + 2 - first(ranges[t])]
        for k in eachindex(subids) for t in (Base.parent(sys, subids[k]),)))
check("the selection loaded exactly the touched tiles",
    Set(freshloaded) == Set(Int(t.index) for t in touched) &&
    length(freshloaded) == length(touched);
    detail="loaded $(length(freshloaded)) of 64800 chunks")

# --- 6. Conservative regrid of the window onto IGEO7, fed through the cube. ---

igeo7 = DGG.IGeo7System()
pixel_area = DGG.cell_area(DGG.levelgrid(sys, 1), subids[1])
L = argmin(l -> abs(log(4π / DGG.ncells(igeo7, l) / pixel_area)), DGG.levels(igeo7))
kfold = DGG.ncells(igeo7, L + 1) / DGG.ncells(igeo7, L)
ratio = 4π / DGG.ncells(igeo7, L) / pixel_area
check("IGEO7 level $L is the closest level to a GLO-90 pixel",
    1 / sqrt(kfold) <= ratio <= sqrt(kfold);
    detail="area ratio $(round(ratio; digits=4))")

# Densified so the south edge's poleward bow stays inside; padded by one pixel.
function box_polygon(w, e, s, n; pad_lon=0.0, pad_lat=0.0)
    w, e, s, n = w - pad_lon, e + pad_lon, s - pad_lat, n + pad_lat
    steps = max(1, ceil(Int, (e - w) * 64))
    pts = [(w + (e - w) * k / steps, s) for k in 0:steps]
    append!(pts, [(w + (e - w) * k / steps, n) for k in steps:-1:0])
    push!(pts, pts[1])
    return GI.Polygon([GI.LinearRing(pts)])
end

boxes = [CD.cell_box(sys, c) for c in subids]
footprint = box_polygon(
    minimum(b[1] for b in boxes), maximum(b[2] for b in boxes),
    minimum(b[3] for b in boxes), maximum(b[4] for b in boxes);
    pad_lon=maximum(b[2] - b[1] for b in boxes), pad_lat=1 / CD.lat_intervals(sys))
set = DGG.query(igeo7, DGG.MultiOrderCoverage(footprint); level=L)
dst = DGG.PartialGrid(DGG.CellVector(set; level=L))
src = DGG.PartialGrid(DD.lookup(sub, DGG.Cells))
build = @elapsed r = CR.Regridder(MANIFOLD, dst, src)
note("$(DGG.ncells(src)) pixels -> $(DGG.ncells(dst)) IGEO7 level-$L cells, " *
     "Regridder built in $(round(build; digits=2)) s")

col_err = maximum(abs.(vec(sum(r.intersections; dims=1)) .- r.src_areas) ./ r.src_areas)
check("column sums == source cell areas", col_err <= 1e-9;
    detail="max rel err $col_err")

elev, cover = zeros(DGG.ncells(dst)), zeros(DGG.ncells(dst))
CR.regrid!(elev, r, Float64.(collect(sub)))
CR.regrid!(cover, r, ones(DGG.ncells(src)))
inside = findall(>(0.99), cover)
elev[inside] ./= cover[inside]
check("regridded elevation stays inside the lazily read range",
    !isempty(inside) && minimum(sub) - 1 <= minimum(elev[inside]) &&
    maximum(elev[inside]) <= maximum(sub) + 1;
    detail="dst $(round.(extrema(elev[inside]); digits=1)) m over $(length(inside)) " *
           "covered cells, src $(round.(extrema(sub); digits=1)) m")

println()
println("  total wall time $(round(time() - STARTED; digits=1)) s")
println(FAILURES[] == 0 ? "ALL CHECKS PASSED" : "$(FAILURES[]) CHECK(S) FAILED")
exit(FAILURES[] == 0 ? 0 : 1)
