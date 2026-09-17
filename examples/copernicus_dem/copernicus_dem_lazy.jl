# All of GLO-90 as one lazy DiskArrays vector along the level-1 cell order.
#
# Run in the docs environment for ArchGDAL and DiskArrays:
#
#     julia -t auto --project=docs examples/copernicus_dem_lazy.jl
#
# One chunk per 1°x1° tile, backed by that tile's COG on AWS and fetched on
# first touch; tiles absent from the bucket manifest read as `0.0f0` (ocean
# sits at 0 on the EGM2008 geoid) with no network. Tiles cache in `tempdir()`.

import DiscreteGlobalGrids as DGG
import ConservativeRegridding as CR
import GeometryOps as GO
import GeoInterface as GI
import DimensionalData as DD
import DiskArrays
import ArchGDAL
import Downloads
import Extents

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

const BUCKET = "https://copernicus-dem-90m.s3.amazonaws.com"
const CACHE_TILES = 16

"Download `url` to `path` once; atomic rename never caches a partial file."
function fetched(url, path)
    if !isfile(path)
        part = "$path.part-$(getpid())"
        Downloads.download(url, part)
        mv(part, path; force=true)
    end
    return path
end

"The AWS object stem of `tile`, e.g. `Copernicus_DSM_COG_30_N50_00_E006_00_DEM`."
function tilestem(sys, tile)
    lat, lon = CD.tilecorner(sys, tile)
    return string("Copernicus_DSM_COG_30_", lat < 0 ? "S" : "N", lpad(abs(lat), 2, '0'),
        "_00_", lon < 0 ? "W" : "E", lpad(abs(lon), 3, '0'), "_00_DEM")
end

"""
Every GLO-90 pixel as one `Float32` vector in level-1 cell order. Chunk `k` is
tile `k`'s `descendant_range`, so reads are tile-aligned: a touched tile's COG
is downloaded once, decoded, and held in a `CACHE_TILES`-tile LRU. A tile
absent from the bucket manifest reads as `0.0f0` without touching the network.
`loaded` records every tile ordinal actually decoded, in decode order.
"""
struct GLO90Vector{S<:DGG.CopernicusDEMSystem,C} <: DiskArrays.AbstractDiskArray{Float32,1}
    sys::S
    chunks::C
    present::Set{String}                # manifest stems
    cache::Dict{Int,Vector{Float32}}    # decoded tiles by ordinal
    order::Vector{Int}                  # LRU order, oldest first
    loaded::Vector{Int}
end

function GLO90Vector(sys::DGG.CopernicusDEMSystem)
    widths = [length(DGG.descendant_range(sys, t, 1)) for t in DGG.rootcells(sys)]
    chunks = DiskArrays.GridChunks(DiskArrays.IrregularChunks(; chunksizes=widths))
    manifest = fetched("$BUCKET/tileList.txt", joinpath(tempdir(), "copdem90-tileList.txt"))
    present = Set{String}(split(read(manifest, String)))
    return GLO90Vector(sys, chunks, present, Dict{Int,Vector{Float32}}(), Int[], Int[])
end

Base.size(A::GLO90Vector) = (DGG.ncells(A.sys, 1),)
DiskArrays.eachchunk(A::GLO90Vector) = A.chunks
DiskArrays.haschunks(::GLO90Vector) = DiskArrays.Chunked()

"The decoded raster of the tile at `ordinal`, as `vec` of the COG's band 1."
function tilevec!(A::GLO90Vector, ordinal::Int, stem::String)
    if haskey(A.cache, ordinal)
        push!(A.order, splice!(A.order, findfirst(==(ordinal), A.order)))
        return A.cache[ordinal]
    end
    tif = fetched("$BUCKET/$stem/$stem.tif", joinpath(tempdir(), "$stem.tif"))
    v = vec(ArchGDAL.read(ds -> ArchGDAL.read(ds, 1), tif))
    push!(A.loaded, ordinal)
    length(A.cache) >= CACHE_TILES && delete!(A.cache, popfirst!(A.order))
    A.cache[ordinal] = v
    push!(A.order, ordinal)
    return v
end

function DiskArrays.readblock!(A::GLO90Vector, out, r::AbstractUnitRange)
    p = first(r)
    while p <= last(r)
        tile = Base.parent(A.sys, DGG.LevelIndex(1, p - 1))
        window = DGG.descendant_range(A.sys, tile, 1)
        stop = min(last(r), last(window))
        seg = (p - first(r) + 1):(stop - first(r) + 1)
        stem = tilestem(A.sys, tile)
        if stem in A.present
            v = tilevec!(A, Int(tile.index), stem)
            out[seg] .= view(v, (p - first(window) + 1):(stop - first(window) + 1))
        else
            fill!(view(out, seg), 0.0f0)
        end
        p = stop + 1
    end
    return out
end

# --- 2. Structure: chunks are exactly the tiles' descendant ranges. ------

sys = DGG.CopernicusDEMSystem(90)
lazy = GLO90Vector(sys)
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
check("manifest separates land from ocean",
    length(lazy.present) >= 26_000 &&
    tilestem(sys, t50) in lazy.present && tilestem(sys, andes) in lazy.present &&
    !(tilestem(sys, ocean) in lazy.present);
    detail="$(length(lazy.present)) stems listed, incl. $(tilestem(sys, andes))")
oceanvals = lazy[DGG.descendant_range(sys, ocean, 1)]
check("an ocean chunk is all zeros with zero loads",
    all(iszero, oceanvals) && isempty(lazy.loaded);
    detail="$(length(oceanvals)) pixels of $(tilestem(sys, ocean))")

# --- 4. A land chunk equals the COG read directly. -----------------------

directvec(t) = vec(ArchGDAL.read(ds -> ArchGDAL.read(ds, 1),
    fetched("$BUCKET/$(tilestem(sys, t))/$(tilestem(sys, t)).tif",
        joinpath(tempdir(), "$(tilestem(sys, t)).tif"))))

check("the N50_00_E006_00 chunk equals its COG",
    lazy[DGG.descendant_range(sys, t50, 1)] == directvec(t50) &&
    lazy.loaded == [Int(t50.index)];
    detail="960000 pixels, 1 tile loaded")

# --- 5. The cube: Covering selection loads only the tiles it touches. ----

fresh = GLO90Vector(sys)                           # its own load record
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
    Set(fresh.loaded) == Set(Int(t.index) for t in touched) &&
    length(fresh.loaded) == length(touched);
    detail="loaded $(length(fresh.loaded)) of 64800 chunks")

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
