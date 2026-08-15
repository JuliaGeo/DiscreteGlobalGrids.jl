# Demo: one real Copernicus DEM tile as a DGGS chunk, regridded onto IGEO7 and HEALPix.
#
# The claim: a Copernicus DEM tile is an ordinary `AbstractGrid`, so the same three
# lines that move a HEALPix level onto a lon/lat mesh move a 90 m DEM tile onto an
# equal-area DGGS —
#
#     src = DGG.PartialGrid(DGG.CopernicusDEMSystem(90), tile, 1)
#     r   = CR.Regridder(GO.Spherical(; radius = 1.0), dst, src)
#     CR.regrid!(field, r, vec(dem))
#
# — with no `CellBasedGrid` adapter, no corner matrix, and no per-system tree.
#
# Environment: the DOCS project, because reading a COG needs ArchGDAL:
#
#     julia -t auto --project=docs examples/copernicus_dem.jl
#
# It downloads ONE ~3.5 MB GLO-90 tile over plain HTTPS with `Downloads` (not GDAL's
# /vsicurl, whose libcurl does not see the macOS trust store) into `tempdir()`, and
# reuses it if it is already there. It ends in PASS/FAIL assertions and exits non-zero
# if any of them fail.
#
# Cost of this script, ONE run on an M-series laptop with `-t auto` (8 threads):
# 38.4 s wall for the whole tile — 960 000 pixels — onto BOTH destinations at
# MATCHED cell size, of which 11.0 s and 8.3 s are the two `Regridder` builds and
# 4.5 s is the GLO-30 section. Those are machine-local and they move a few percent
# from run to run: `scripts/bench_copdem_cursor.jl`, which is the reproducer for
# this pair and for the source-tree table in `src/systems/CopernicusDEM/cursor.jl`,
# measured the same two builds at 11.5 s and 9.2 s in its own run — against
# 105.7 s and 83.9 s for the generic cursor. `COPDEM_ROWS=n` cuts the source to the
# northernmost `n` raster rows if you want it faster; the default is the whole
# tile, because that is the claim.

import DiscreteGlobalGrids as DGG
import ConservativeRegridding as CR
import GeometryOps as GO
import GeoInterface as GI
import DimensionalData as DD
import ArchGDAL
import Downloads
import Extents

const CD = DGG.CopernicusDEM
const US = GO.UnitSpherical

const FAILURES = Ref(0)
function check(name, ok; detail="")
    ok || (FAILURES[] += 1)
    println(ok ? "PASS  " : "FAIL  ", rpad(name, 56), detail)
    return ok
end
note(text) = println("      ", text)

# THE MANIFOLD, declared once: every grid in this package computes on the UNIT
# sphere, so that is the manifold the regridder must work on. See
# `examples/regridding.jl` for why naming it is not optional.
const MANIFOLD = GO.Spherical(; radius=1.0)
const R_EARTH = 6_371_008.0                # for printing areas in m^2

const STARTED = time()

println("="^78)
println("copernicus_dem.jl — one AWS GLO-90 tile onto IGEO7 and HEALPix")
println("julia $(VERSION)  threads=$(Threads.nthreads())")
println("="^78)

# --------------------------------------------------------------------------
# 1. Fetch the tile, and check the file against the lattice.
#
# The band rule says a GLO-90 tile whose equator-ward edge is at latitude 50 is
# 1.5x reduced: 800 columns of 1200 rows. The geotransform says where its
# pixel-is-point raster sits. Both are predictions this system makes from a
# static table, and the COG on AWS is the oracle.
# --------------------------------------------------------------------------

const TILE_LAT, TILE_LON = 50, 6          # N50_00_E006_00 — band [50,60), 1.5x
const STEM = "Copernicus_DSM_COG_30_N50_00_E006_00_DEM"
const URL = "https://copernicus-dem-90m.s3.amazonaws.com/$STEM/$STEM.tif"
const TIF = joinpath(tempdir(), "$STEM.tif")

if !isfile(TIF)
    println("  downloading $URL")
    # To a scratch name and then `mv`: a download interrupted halfway would
    # otherwise leave a truncated file at the cached path, and every later run
    # would read that instead of re-fetching. `mv` within one directory is atomic.
    part = "$TIF.part-$(getpid())"
    Downloads.download(URL, part)
    mv(part, TIF; force=true)
end
note("tile file: $TIF ($(round(filesize(TIF) / 2^20; digits=2)) MiB)")

sys = DGG.CopernicusDEMSystem(90)
tile = CD.tilecell(sys, TILE_LAT, TILE_LON)
g0 = DGG.levelgrid(sys, 0)
g1 = DGG.levelgrid(sys, 1)
const TILE_W, TILE_E, TILE_S, TILE_N = CD.cell_box(sys, tile)
const NCOLS = Int(CD.ncols_at(sys, TILE_LAT))
const NROWS = Int(CD.lat_intervals(sys))

A, gt, tifsize = ArchGDAL.read(TIF) do ds
    (ArchGDAL.read(ds, 1), ArchGDAL.getgeotransform(ds),
     (ArchGDAL.width(ds), ArchGDAL.height(ds)))
end

check("COG size is the band table's (ncols, nrows)",
    tifsize == (NCOLS, NROWS) == size(A);
    detail="$(tifsize) in a 1.5x band")

# `cell_box` of the TILE is the raster's outer frame — a half pixel outside the
# corner posts — which is exactly GDAL's origin convention for a pixel-is-point
# COG. So the geotransform is not an independent number to look up: it is four
# `cell_box` values and two divisions.
expected_gt = [TILE_W, (TILE_E - TILE_W) / NCOLS, 0.0,
    TILE_N, 0.0, -(TILE_N - TILE_S) / NROWS]
gt_err = maximum(abs.(gt .- expected_gt))
check("geotransform == cell_box(tile) frame", gt_err <= 1e-12;
    detail="max abs err $gt_err")
note("origin ($(gt[1]), $(gt[4]))  dlon $(gt[2])  dlat $(gt[6])")

# --------------------------------------------------------------------------
# 2. The chunk, and the cube axis.
#
# `PartialGrid(sys, tile, 1)` is O(1): sorted subtrees make a tile's pixels a
# lazy window over the level-1 ids, so nothing is materialised. `CellLookup` is
# that window wearing a DimensionalData hat.
# --------------------------------------------------------------------------

src = DGG.PartialGrid(sys, tile, 1)

# THE FLATTENING. Our position `k` for pixel `(j, i)` is `j*ncols + i + 1`, and
# `vec` of an `(ncols, nrows)` matrix gives linear index `j*ncols + i + 1` — so
# `vec(A)` is already in the chunk's position order, with no `permutedims` and
# no copy. The whole demo rests on that line, so it is asserted, not assumed.
values = vec(A)
check("vec(A) is the chunk's position order",
    all(values[DGG.cellposition(src, CD.pixelcell(sys, tile, j, i))] == A[i+1, j+1]
        for (j, i) in ((0, 0), (0, NCOLS - 1), (NROWS - 1, 0),
        (NROWS - 1, NCOLS - 1), (37, 421))))

lk = DGG.CellLookup(src)
dem = DD.DimArray(values, DGG.Cells(lk))

check("chunk holds the whole raster", DGG.ncells(src) == NCOLS * NROWS;
    detail="$(DGG.ncells(src)) cells")
check("positions are chunk-local",
    DGG.cellposition(src, DGG.cellindex(src, 1)) == 1 &&
    DGG.cellposition(src, DGG.cellindex(src, DGG.ncells(src))) == DGG.ncells(src))
check("the axis names level-1 cells of this system",
    DGG.level(lk) == 1 && DGG.system(lk) === sys)
check("PartialGrid(lk) round-trips the chunk",
    DGG.ncells(DGG.PartialGrid(lk)) == DGG.ncells(src))

# Requirement (b) of the brief: select a region out of the cube by covering, and
# land back on a `CellLookup` — the view's axis is still cells, not integers.
window = Extents.Extent(X=(6.10, 6.11), Y=(50.50, 50.51))
sub = dem[DGG.Cells(DGG.Covering(window))]
sublk = DD.lookup(sub, DGG.Cells)
subids = collect(sublk)
check("Covering(extent) selects a CellLookup sub-cube",
    sublk isa DGG.CellLookup && 0 < length(sub) < DGG.ncells(src);
    detail="$(length(sub)) cells over a 0.01 x 0.01 degree window")
check("selected values are those cells' own values",
    all(sub[k] == values[DGG.cellposition(src, subids[k])] for k in eachindex(subids)))
check("every selected cell meets the window",
    all(Extents.intersects(DGG.cell_extent(g1, c), window) for c in subids))

# The axis is windows, not an id vector: `CellVector` stores the chunk's
# position ranges, so the whole 960 000-cell axis is a few hundred bytes.
note("axis memory: Base.summarysize(lk) = $(Base.summarysize(lk)) B, against " *
     "$(8 * DGG.ncells(src)) B for one Int64 id per cell " *
     "($(round(Int, 8 * DGG.ncells(src) / Base.summarysize(lk)))x)")

# --------------------------------------------------------------------------
# 3. Choosing the destination level.
#
# "At comparable cell size" is a question the systems answer themselves: take
# the level whose mean cell area is closest to a pixel's, on a log scale so the
# candidates either side are compared fairly.
# --------------------------------------------------------------------------

pixel_area = DGG.cell_area(g1, DGG.cellindex(src, 1))
matched(dstsys) = argmin(l -> abs(log(4π / DGG.ncells(dstsys, l) / pixel_area)),
    DGG.levels(dstsys))

const DESTINATIONS = (("IGEO7", DGG.IGeo7System()), ("HEALPix", DGG.HEALPixSystem()))
const LEVELS = Dict(name => matched(dstsys) for (name, dstsys) in DESTINATIONS)

println()
println("  destination   level        cells      mean cell m2      pixel m2    ratio")
for (name, dstsys) in DESTINATIONS
    l = LEVELS[name]
    mean_area = 4π / DGG.ncells(dstsys, l)
    ratio = mean_area / pixel_area
    println("  ", rpad(name, 13), lpad(l, 5), lpad(DGG.ncells(dstsys, l), 13),
        lpad(round(mean_area * R_EARTH^2; digits=1), 18),
        lpad(round(pixel_area * R_EARTH^2; digits=1), 14),
        lpad(round(ratio; digits=4), 9))
    # "Within a factor of 2" is not the right bar, and it is not reachable: a
    # `k`-fold system's levels are `k` apart in area, so the closest one can be
    # off by `sqrt(k)` — 2.65 for IGEO7's 7-fold refinement, 2 for HEALPix's
    # 4-fold. What IS assertable is that no neighbouring level is closer, which
    # is what `sqrt(k)` states and what a broken level search would violate.
    k = DGG.ncells(dstsys, l + 1) / DGG.ncells(dstsys, l)
    check("$name level $l is the closest level to a pixel",
        1 / sqrt(k) <= ratio <= sqrt(k);
        detail="ratio $(round(ratio; digits=4)), best possible $(round(sqrt(k); digits=3))x " *
               "for $(round(Int, k))-fold refinement")
end

# --------------------------------------------------------------------------
# 4. The coverings, the regrids, and what they must conserve.
# --------------------------------------------------------------------------

# How many raster rows of the tile the matched-resolution regrids use. The
# default is the whole tile, because "one AWS tile onto a DGGS at comparable
# cell size" is the claim and a fraction of a tile does not make it.
#
# `levels(sys) == 0:1` leaves the generic cursor no interior structure to
# descend — the tile node's children ARE its 960 000 pixels — so the dual tree
# search would cost O(n_src x dst-depth). `src/systems/CopernicusDEM/cursor.jl`
# gives the lattice a real tree instead, by recursively bisecting the pixel
# rectangle, and `treeify` picks it up with nothing named here. The two build
# times this file prints below are the payoff: the same two builds cost 105.7 s
# and 83.9 s on the generic cursor, for an intersection matrix identical to the
# last bit. `scripts/bench_copdem_cursor.jl` is the reproducer for that pair and
# for the whole fanout table in that file's `BlockStrategy` docstring; those
# numbers are one run on the machine named at the top of this file.
const ROWS = let n = parse(Int, get(ENV, "COPDEM_ROWS", string(NROWS)))
    1 <= n <= NROWS || error("COPDEM_ROWS=$n is outside 1:$NROWS — this tile has " *
                             "$NROWS raster rows and the chunk is its northernmost `n`")
    n
end

"The chunk of the tile's northernmost `rows` raster rows, and its values."
function row_band(rows)
    rows == NROWS && return (src, values)
    base = DGG.cellindex(src, 1).index
    ids = [DGG.LevelIndex(1, k) for k in base:(base+rows*NCOLS-1)]
    return (DGG.PartialGrid(sys, 1, ids), values[1:(rows*NCOLS)])
end

"""
The lon/lat box of a contiguous chunk, read off the first and last cell's
`cell_box` — never off a ring, whose poleward edge bows past the box.
"""
function chunk_box(grid)
    w1, e1, s1, n1 = CD.cell_box(sys, DGG.cellindex(grid, 1))
    w2, e2, s2, n2 = CD.cell_box(sys, DGG.cellindex(grid, DGG.ncells(grid)))
    return (min(w1, w2), max(e1, e2), min(s1, s2), max(n1, n2))
end

# A lon/lat box, as a spherical polygon that really contains the box.
#
# Two corrections, and both are needed — an `Extents.Extent` handed straight to
# `query` gets neither, and the regrid then loses whole pixels:
#
#  1. DENSIFY. A parallel is not a great circle, so the arc joining a box's two
#     south corners bows POLEWARD of the parallel by about `(dlam^2/8) sin(phi)
#     cos(phi)`. `query`'s own extent conversion samples at
#     `EXTENT_STEP_DEGREES = 2.0`, i.e. one segment across a 1-degree tile, and
#     at this tile's latitude that bow is 1.87e-5 rad = 0.00107 degrees — 1.3
#     GLO-90 pixel rows. The southernmost raster row then falls OUTSIDE the
#     covering near mid-longitude and regrids to nothing: measured on a 64-row
#     band handed to `query` as a bare `Extents.Extent`, 898 of its 51 200
#     columns came back short of their own area by more than 1e-9 relative, and
#     the worst lost all of it. At 64 segments per degree the same bow is
#     4.58e-9 rad = 2.62e-7 degrees, 3.1e-4 of one pixel row.
#  2. PAD. The published pixel rings bow poleward of their own boxes too, by
#     about 3e-11 rad. A covering of the box exactly would clip those slivers
#     off the outermost rows — 1.2e-5 of a boundary pixel's area, five orders
#     above the interior noise floor. One pixel of pad swallows them.
const DENSIFY_PER_DEGREE = 64

function box_polygon(w, e, s, n; pad_lon=0.0, pad_lat=0.0)
    w, e = w - pad_lon, e + pad_lon
    s, n = max(s - pad_lat, -90.0), min(n + pad_lat, 90.0)
    steps = max(1, ceil(Int, (e - w) * DENSIFY_PER_DEGREE))
    pts = Tuple{Float64,Float64}[]
    for k in 0:steps                        # east along the south edge
        push!(pts, (w + (e - w) * k / steps, s))
    end
    for k in steps:-1:0                     # west along the north edge
        push!(pts, (w + (e - w) * k / steps, n))
    end
    push!(pts, pts[1])
    return GI.Polygon([GI.LinearRing(pts)])
end

"""
How a ring turns, which is what says whether Sutherland-Hodgman can use it as a
clip window: `:ccw_convex` (every turn left, seen from outside the sphere — a
valid window), `:cw_convex` (convex, but wound the other way, which is NOT one),
or `:reflex` (it turns both ways).

ORIENTATION is why this is not the one-line "does it turn right anywhere" test
this repo's regridding suites carry: that test reads a clockwise CONVEX ring as
reflex, so a destination system that emitted its rings clockwise would take the
loose branch below without anything saying which of the two defects it had.
"""
function ring_shape(poly)
    pts = collect(GI.getpoint(GI.getexterior(poly)))
    while length(pts) > 1 && pts[end] == pts[1]
        pop!(pts)
    end
    pts = pts[[i for i in eachindex(pts) if i == 1 || pts[i] != pts[i-1]]]
    n = length(pts)
    n < 3 && return :ccw_convex          # nothing left to turn
    turns = [US.spherical_orient(pts[i], pts[mod1(i + 1, n)], pts[mod1(i + 2, n)])
             for i in 1:n]
    all(>=(0), turns) && return :ccw_convex
    all(<=(0), turns) && return :cw_convex
    return :reflex
end

lat_of(p) = asind(clamp(p[3], -1.0, 1.0))
lon_of(p) = atand(p[2], p[1])

"""
Cover `chunk`'s footprint with `dstsys` cells at level `L`, build the regridder,
and run every law a conservative regrid off this system owes: the matrix's
shape, the column sums against the source areas, the area budget against the
closed-form box areas, and — because a constant field cannot catch a transposed
or flipped raster — the regridded latitude AND longitude against the
destination cells' own centroids.
"""
function regrid_onto(label, dstsys, L, chunk, chunkvalues)
    w, e, s, n = chunk_box(chunk)
    footprint = box_polygon(w, e, s, n;
        pad_lon=(e - w) / NCOLS, pad_lat=(TILE_N - TILE_S) / NROWS)
    set = DGG.query(dstsys, DGG.MultiOrderCoverage(footprint); level=L)
    dst = DGG.PartialGrid(DGG.CellVector(set; level=L))
    build = @elapsed r = CR.Regridder(MANIFOLD, dst, chunk)
    n_src, n_dst = DGG.ncells(chunk), DGG.ncells(dst)
    dstgrid = DGG.levelgrid(dstsys, L)
    println()
    println("  $label: $n_src pixels -> $n_dst $(nameof(typeof(dstsys))) level-$L " *
            "cells, Regridder built in $(round(build; digits=2)) s")

    # WHICH ARM the laws below are checked on is derived, not chosen.
    # Sutherland-Hodgman clips the SOURCE cell against the DESTINATION cell's
    # ring, and only a counter-clockwise CONVEX clip window makes that clip the
    # intersection. IGEO7 cells are plain spherical polygons; HEALPix
    # `cell_boundary` densifies each curved chart edge into eight arcs, so every
    # vertex on a bowed side is reflex and the clipper silently loses area. Same
    # defect as the twelve `@test_broken` cases in
    # `test/systems/crosssystem/regridding_conservation.jl:197-205`, and it
    # belongs to the destination, not to this system — whose own 4-corner quads
    # are convex, which is checked once before the regrids begin.
    #
    # The CONSTANTS are chosen, one per law rather than one reused, each a
    # decade or more above what the whole tile measures on that arm:
    #
    #   law                             exact arm (IGEO7)   clipped arm (HEALPix)
    #   column sums == source areas      2.6e-10 -> 1e-9      6.0e-5 -> 1e-4
    #   no destination over-covered      5.7e-11 -> 1e-9      1.5e-5 -> 1e-4
    #   interior cells covered exactly   5.7e-11 -> 1e-9      5.1e-5 -> 1e-4
    #
    # The three land on the same pair of decades and are not the same number for
    # the same reason: on the exact arm every residue is float rounding over the
    # handful of pieces a pixel is cut into, and on the clipped arm each is a
    # different aggregate of the destination's own bow — which is the defect
    # being reported, not a tolerance being granted.
    shape = ring_shape(DGG.cell_polygon(dstgrid, DGG.cellindex(dst, 1)))
    convex = shape === :ccw_convex
    col_tol = convex ? 1e-9 : 1e-4
    over_tol = convex ? 1e-9 : 1e-4
    deep_tol = convex ? 1e-9 : 1e-4
    nvert = length(DGG.cell_boundary(dstgrid, DGG.cellindex(dst, 1)))
    if convex
        check("$label: destination rings are convex clip windows", true;
            detail="$nvert-vertex rings, counter-clockwise, so the clip is exact " *
                   "and conservation must be too")
    elseif shape === :cw_convex
        note("$label: destination rings are CONVEX but wound CLOCKWISE " *
             "($nvert-vertex rings), which Sutherland-Hodgman cannot use as a clip " *
             "window either — the loose bounds below are taken for that reason, " *
             "which is a different defect from a reflex vertex and is said here " *
             "rather than passed over.")
    else
        note("$label: destination rings have REFLEX vertices ($nvert-vertex rings, " *
             "eight arcs per densified chart edge, so every vertex on a bowed side " *
             "turns the wrong way), so Sutherland-Hodgman loses a little of every " *
             "clipped source cell and conservation below is bounded at $col_tol " *
             "rather than exact. Upstream defect, and it belongs to the " *
             "DESTINATION: see the twelve @test_broken cases in " *
             "test/systems/crosssystem/regridding_conservation.jl:197-205, and the " *
             "same shortfall printed by examples/regridding.jl's last loop.")
    end

    M = r.intersections
    check("$label: matrix is (dst cells, src cells)", size(M) == (n_dst, n_src);
        detail="$(size(M)), nnz=$(length(M.nzval))")

    # THE conservation assertion. The covering OVER-covers the chunk, so the row
    # sums are NOT the destination areas — only the columns are exact. Every
    # source pixel must be fully consumed, and that is the statement that would
    # fail if this system's rings were non-convex, or did not tile, or if the
    # covering left a sliver.
    col_err = maximum(abs.(vec(sum(M; dims=1)) .- r.src_areas) ./ r.src_areas)
    check("$label: column sums == source cell areas", col_err <= col_tol;
        detail="max rel err $col_err (tolerance $col_tol)")
    # The same law summed. Not an independent bound — the per-column errors
    # cancel, so this can only be tighter — but a sign error or a lost block
    # would survive the maximum above and not this.
    check("$label: total intersection area == total source area",
        abs(sum(M) - sum(r.src_areas)) / sum(r.src_areas) <= col_tol;
        detail="$(sum(M)) vs $(sum(r.src_areas)) sr")

    # The area budget against the closed form. `src_areas` measures the published
    # geodesic QUADS; `cell_area` is the exact solid angle of the BOXes those
    # quads approximate. Materialised, not a generator — `cell_area`'s docstring
    # has why: this closed form does not telescope, and Julia's pairwise `sum`
    # over a `Vector` lands some 400x closer to the truth than a generator's
    # sequential accumulation.
    boxes = [DGG.cell_area(g1, DGG.cellindex(chunk, i)) for i in 1:n_src]
    box_sum = sum(boxes)
    ring_gap = abs(sum(r.src_areas) - box_sum) / box_sum
    check("$label: ring areas == box areas to the quad gap", ring_gap <= 1e-9;
        detail="rel gap $ring_gap")
    if n_src == NCOLS * NROWS
        tile_gap = abs(box_sum - DGG.cell_area(g0, tile)) / DGG.cell_area(g0, tile)
        check("$label: the pixel boxes tile the tile's box", tile_gap <= 1e-14;
            detail="rel gap $tile_gap over $(round(box_sum; sigdigits=8)) sr")
    end

    # `regrid!` divides by `dst_areas`, so a regridded field of ones IS each
    # destination cell's covered fraction.
    cover = zeros(n_dst)
    CR.regrid!(cover, r, ones(n_src))
    # Its own bound, and its own physics: this one is about a destination cell
    # being handed MORE than its own area, which the source tessellation's
    # overlaps and the destination's own area formula decide — not about a
    # source pixel being fully consumed.
    check("$label: no destination cell is over-covered", maximum(cover) <= 1 + over_tol;
        detail="max cover - 1 = $(maximum(cover) - 1) (tolerance $over_tol), " *
               "cover in $(round.(extrema(cover); digits=6))")

    # Orientation. Two analytic fields, not one: a row flip and a column flip are
    # different mistakes and a single tilted field would confound them. Both are
    # sampled at SOURCE cell centroids and compared against the DESTINATION
    # cell's own centroid, so nothing in between is trusted. A flip puts the
    # error near half the tile — three orders above the threshold.
    src_lat = [lat_of(DGG.cell_centroid(g1, DGG.cellindex(chunk, i))) for i in 1:n_src]
    src_lon = [lon_of(DGG.cell_centroid(g1, DGG.cellindex(chunk, i))) for i in 1:n_src]
    got_lat, got_lon = zeros(n_dst), zeros(n_dst)
    CR.regrid!(got_lat, r, src_lat)
    CR.regrid!(got_lon, r, src_lon)

    inside = findall(>(0.99), cover)
    width = rad2deg(sqrt(4π / DGG.ncells(dstsys, L)))   # one destination cell, in degrees
    centroids = [DGG.cell_centroid(dstgrid, DGG.cellindex(dst, k)) for k in inside]
    lat_err = maximum(abs(got_lat[k] / cover[k] - lat_of(p))
                      for (k, p) in zip(inside, centroids))
    lon_err = maximum(abs(got_lon[k] / cover[k] - lon_of(p))
                      for (k, p) in zip(inside, centroids))
    check("$label: regridded latitude is the cell's own latitude", lat_err < 5 * width;
        detail="max err $(round(lat_err; sigdigits=4)) deg, 5 cells = $(round(5 * width; sigdigits=4)) deg")
    check("$label: regridded longitude is the cell's own longitude", lon_err < 5 * width;
        detail="max err $(round(lon_err; sigdigits=4)) deg over $(length(inside)) covered cells")

    # A cell well inside the footprint must be covered exactly: a gap in the
    # source tessellation would show here and nowhere else.
    margin = 5 * width
    if n - s > 4 * margin
        deep = [k for (k, p) in zip(inside, centroids)
                if w + margin < lon_of(p) < e - margin &&
                   s + margin < lat_of(p) < n - margin]
        deep_err = isempty(deep) ? NaN : maximum(abs(cover[k] - 1) for k in deep)
        check("$label: interior cells are covered exactly",
            !isempty(deep) && deep_err <= deep_tol;
            detail="$(length(deep)) interior cells, max |cover - 1| = $deep_err " *
                   "(tolerance $deep_tol)")
    else
        note("$label: chunk is thinner than 4 destination cells — no interior " *
             "cells to check coverage on; run with COPDEM_ROWS unset")
    end

    # The elevation the whole exercise is for: the cover-weighted mean of the DEM.
    elevation = zeros(n_dst)
    CR.regrid!(elevation, r, Float64.(chunkvalues))
    elevation[inside] ./= cover[inside]
    note("elevation over covered cells: $(round.(extrema(elevation[inside]); digits=1)) m, " *
         "source range $(round.(extrema(chunkvalues); digits=1)) m")
    check("$label: regridded elevation stays inside the source range",
        minimum(chunkvalues) - 1 <= minimum(elevation[inside]) &&
        maximum(elevation[inside]) <= maximum(chunkvalues) + 1)
    return build
end

# This system's own rings are the convex 4-corner quads it promises — the
# premise the IGEO7 arm above rests on, checked here rather than assumed.
check("source pixel rings are convex quads",
    all(ring_shape(DGG.cell_polygon(g1, CD.pixelcell(sys, tile, j, i))) === :ccw_convex &&
        length(DGG.cell_boundary(g1, CD.pixelcell(sys, tile, j, i))) == 4
        for (j, i) in ((0, 0), (599, 400), (NROWS - 1, NCOLS - 1))))

chunk, chunkvalues = row_band(ROWS)
println()
println("  matched-resolution regrids over $ROWS of $NROWS raster rows " *
        "($(DGG.ncells(chunk)) pixels)")
check("the chunk gets the block cursor", DGG.treeify(chunk) isa CD.BlockCursor;
    detail="$(nameof(typeof(DGG.treeify(chunk)))) from " *
           "src/systems/CopernicusDEM/cursor.jl, chosen by `treeify` with nothing " *
           "named here")
builds = Dict(name => regrid_onto(name, dstsys, LEVELS[name], chunk, chunkvalues)
              for (name, dstsys) in DESTINATIONS)

# --------------------------------------------------------------------------
# 5. The GLO-30 call shape, without a download.
#
# Nothing above is specific to GLO-90. `CopernicusDEMSystem(30)` is the same
# code at N = 3600, its pixels nest 3x3 inside GLO-90's, and a sub-window of one
# is a `PartialGrid` like any other.
# --------------------------------------------------------------------------

sys30 = DGG.CopernicusDEMSystem(30)
tile30 = CD.tilecell(sys30, TILE_LAT, TILE_LON)
pixel = CD.pixelcell(sys, tile, 600, 400)
coarse_box = CD.cell_box(sys, pixel)
fine = CD.refine(sys, sys30, pixel)
fine_boxes = [CD.cell_box(sys30, f) for f in fine]
check("a GLO-90 pixel refines into 9 GLO-30 pixels",
    length(fine) == 9 && all(DGG.level(f) == 1 for f in fine) &&
    all(CD.coarsen(sys30, sys, f) == pixel for f in fine))
check("the 9 fine boxes tile a box of the coarse box's size",
    isapprox(maximum(b[2] for b in fine_boxes) - minimum(b[1] for b in fine_boxes),
        coarse_box[2] - coarse_box[1]; rtol=1e-12) &&
    isapprox(maximum(b[4] for b in fine_boxes) - minimum(b[3] for b in fine_boxes),
        coarse_box[4] - coarse_box[3]; rtol=1e-12))
note("dlon: GLO-90 $(coarse_box[2] - coarse_box[1]) deg, " *
     "GLO-30 $(fine_boxes[1][2] - fine_boxes[1][1]) deg")

# A 256x256 sub-window of the GLO-30 tile: 256 runs of 256 consecutive ids,
# ascending, which is a legal `PartialGrid` like any other.
const WIN = 256
base30 = CD.pixelcell(sys30, tile30, 0, 0).index
nc30 = Int(CD.ncols_at(sys30, TILE_LAT))
src30 = DGG.PartialGrid(sys30, 1,
    [DGG.LevelIndex(1, base30 + Int64(j) * nc30 + i)
     for j in 0:(WIN-1) for i in 0:(WIN-1)])
g1_30 = DGG.levelgrid(sys30, 1)
igeo7 = DGG.IGeo7System()
pixel_area30 = DGG.cell_area(g1_30, DGG.cellindex(src30, 1))
L30 = argmin(l -> abs(log(4π / DGG.ncells(igeo7, l) / pixel_area30)), DGG.levels(igeo7))

nw30 = CD.cell_box(sys30, DGG.cellindex(src30, 1))
se30 = CD.cell_box(sys30, DGG.cellindex(src30, DGG.ncells(src30)))
box30 = box_polygon(min(nw30[1], se30[1]), max(nw30[2], se30[2]),
    min(nw30[3], se30[3]), max(nw30[4], se30[4]);
    pad_lon=1 / nc30, pad_lat=1 / 3600)
set30 = DGG.query(igeo7, DGG.MultiOrderCoverage(box30); level=L30)
dst30 = DGG.PartialGrid(DGG.CellVector(set30; level=L30))
build30 = @elapsed r30 = CR.Regridder(MANIFOLD, dst30, src30)
ramp = [lat_of(DGG.cell_centroid(g1_30, DGG.cellindex(src30, i)))
        for i in 1:DGG.ncells(src30)]
out30, cover30 = zeros(DGG.ncells(dst30)), zeros(DGG.ncells(dst30))
CR.regrid!(out30, r30, ramp)
CR.regrid!(cover30, r30, ones(DGG.ncells(src30)))

println()
println("  GLO-30 $(WIN)x$(WIN) sub-window -> IGEO7 level $L30: $(DGG.ncells(dst30)) " *
        "cells, Regridder built in $(round(build30; digits=2)) s")
col30 = maximum(abs.(vec(sum(r30.intersections; dims=1)) .- r30.src_areas) ./ r30.src_areas)
check("GLO-30: column sums == source cell areas", col30 <= 1e-9;
    detail="max rel err $col30")
inside30 = findall(>(0.99), cover30)
g30 = DGG.levelgrid(igeo7, L30)
ramp_err = maximum(abs(out30[k] / cover30[k] -
                       lat_of(DGG.cell_centroid(g30, DGG.cellindex(dst30, k))))
                   for k in inside30)
check("GLO-30: the latitude ramp lands where it started",
    ramp_err < 5 * rad2deg(sqrt(4π / DGG.ncells(igeo7, L30)));
    detail="max err $(round(ramp_err; sigdigits=4)) deg over $(length(inside30)) cells")

println()
note("call site, verbatim:")
note("  src = DGG.PartialGrid(DGG.CopernicusDEMSystem(90), tile, 1)")
note("  set = DGG.query(dstsys, DGG.MultiOrderCoverage(footprint); level = L)")
note("  dst = DGG.PartialGrid(DGG.CellVector(set; level = L))")
note("  r   = CR.Regridder(GO.Spherical(; radius = 1.0), dst, src)")
note("  CR.regrid!(field, r, vec(dem))")
note("Regridder builds: " *
     join(["$name $(round(builds[name]; digits=1)) s" for (name, _) in DESTINATIONS], ", ") *
     ", GLO-30 window $(round(build30; digits=1)) s")

println()
println("  total wall time $(round(time() - STARTED; digits=1)) s")
println(FAILURES[] == 0 ? "ALL CHECKS PASSED" : "$(FAILURES[]) CHECK(S) FAILED")
exit(FAILURES[] == 0 ? 0 : 1)
