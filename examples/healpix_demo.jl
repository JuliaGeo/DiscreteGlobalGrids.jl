# Demo: a HEALPix vector data cube on PARTIAL coverage — zonal means, stencils.
#
# A cropped cube stores one region's cells, not the globe's. `PartialGrid` is
# that crop, and it is the only object the rest of the demo needs:
#
#     ids  = DGG.query(sys, DGG.Intersects(window); level)
#     grid = DGG.PartialGrid(sys, level, ids)
#
# Indices `1:ncells(grid)` index the data vector; `query` selects zones in
# that same index space through `localindex`; `neighbors` gives the stencil
# ring, with the slots outside the crop simply missing rather than padded.
#
# Environment: needs nothing beyond DiscreteGlobalGrids and its dependencies
# (Healpix.jl is one of them, and is used here as an independent oracle):
#
#     julia -t 4 --project=. examples/healpix_demo.jl
#
# It ends in PASS/FAIL assertions and exits non-zero if any of them fail.

import DiscreteGlobalGrids as DGG
import GeometryOps as GO
import Extents
import Healpix
using Statistics

const FAILURES = Ref(0)
function check(name, ok; detail="")
    ok || (FAILURES[] += 1)
    println(ok ? "PASS  " : "FAIL  ", rpad(name, 52), detail)
    return ok
end
note(text) = println("      ", text)

const SYS = DGG.HEALPixSystem()
const LEVEL = 7                                    # nside = 128, 196608 cells
const WINDOW = Extents.Extent(X=(-12.0, 32.0), Y=(34.0, 62.0))
const LONLAT = GO.UnitSpherical.GeographicFromUnitSphere()

println("="^78)
println("healpix_demo.jl — level-$LEVEL crop over $(WINDOW)")
println("julia $(VERSION)  threads=$(Threads.nthreads())")
println("="^78)

# --------------------------------------------------------------------------
# 1. The crop.
# --------------------------------------------------------------------------

globe = DGG.levelgrid(SYS, LEVEL)
ids = DGG.query(SYS, DGG.Intersects(WINDOW); level=LEVEL)
grid = DGG.PartialGrid(SYS, LEVEL, ids)

check("the crop is a strict subset of the globe",
    0 < DGG.ncells(grid) < DGG.ncells(globe);
    detail="$(DGG.ncells(grid)) of $(DGG.ncells(globe)) cells " *
           "($(round(100 * DGG.ncells(grid) / DGG.ncells(globe); digits=2))%)")
check("indices round trip through the crop",
    all(DGG.localindex(grid, DGG.cellindex(grid, i)) == i
        for i in (1, DGG.ncells(grid) ÷ 2, DGG.ncells(grid))))

# The EOPF claim: the complete level's dense order IS Healpix.jl's nested pixel
# order, so a HealpixMap and a level-grid data vector are the same vector.
res = Healpix.Resolution(2^LEVEL)
theta, phi = Healpix.pix2angNest(res, 1)
lon, lat = LONLAT(DGG.cell_centroid(globe, DGG.cellindex(globe, 1)))
check("globe index 1 is Healpix.jl nested pixel 1",
    rad2deg(phi) ≈ mod(lon, 360) && 90 - rad2deg(theta) ≈ lat;
    detail="($(round(lon; digits=4)), $(round(lat; digits=4)))")

# --------------------------------------------------------------------------
# 2. A field on the crop, and equal-area means.
# --------------------------------------------------------------------------

centers = [LONLAT(DGG.cell_centroid(grid, DGG.cellindex(grid, i)))
           for i in 1:DGG.ncells(grid)]
field(lon, lat) = 20 - 0.5 * (lat - 34) + 2 * sind(3 * lon)
data = [field(c...) for c in centers]

exact = 4pi / DGG.ncells(globe)
areas = [DGG.cell_area(globe, DGG.cellindex(grid, i)) for i in 1:DGG.ncells(grid)]
check("HEALPix cells are exactly equal-area",
    maximum(abs.(areas .- exact)) / exact < 1e-12;
    detail="$(round(exact; sigdigits=6)) sr each, so an unweighted mean IS areal")

# Ask the COMPLETE grid for the area, not the crop: HEALPix's exact
# `4pi/ncells` override is attached to its level grid, and a `PartialGrid` falls
# through to the spherical area of the published boundary ring. That ring is
# already densified — eight great-circle segments per chart edge, 32 vertices —
# but a chord-wise polygon still under-reads a curvilinear diamond, by the
# relative amount printed below.
subset_areas = [DGG.cell_area(grid, DGG.cellindex(grid, i)) for i in 1:DGG.ncells(grid)]
note("cell_area(crop, c) differs from cell_area(globe, c) by " *
     "$(round(maximum(abs.(subset_areas .- areas)) / exact; sigdigits=3)) relative " *
     "— the subset grid does not inherit the system's area override")

# --------------------------------------------------------------------------
# 3. Zonal statistics, as a query in the crop's index space.
#
# `query` returns typed ids; `localindex` turns each into an index into
# `data`. Two predicates bracket the answer: `Within` keeps only cells wholly
# inside the zone, `Intersects` keeps every cell that touches it, and the
# centre-in-zone rule that raster zonal statistics use sits between them.
# --------------------------------------------------------------------------

zones = (Iberia=Extents.Extent(X=(-9.0, 3.0), Y=(36.0, 44.0)),
    France=Extents.Extent(X=(-4.0, 8.0), Y=(43.0, 51.0)),
    Poland=Extents.Extent(X=(14.0, 24.0), Y=(49.0, 55.0)))

zone_indices(g, target) = Int[DGG.localindex(g, c) for c in DGG.query(g, target)]

inside(ext, (lon, lat)) = ext.X[1] <= lon <= ext.X[2] && ext.Y[1] <= lat <= ext.Y[2]

println()
println("  zone      within  centre  touching     mean (within / centre)")
for (name, ext) in pairs(zones)
    within = zone_indices(grid, DGG.Within(ext))
    touching = zone_indices(grid, DGG.Intersects(ext))
    centre = findall(c -> inside(ext, c), centers)
    println("  ", rpad(name, 9), lpad(length(within), 6), lpad(length(centre), 8),
        lpad(length(touching), 8), "     ",
        round(mean(data[within]); digits=3), " / ", round(mean(data[centre]); digits=3))
    check("$name: within ⊆ centre ⊆ touching",
        issubset(within, centre) && issubset(centre, touching))
    check("$name: the two means agree to the cell size",
        abs(mean(data[within]) - mean(data[centre])) < 0.5)
end

# --------------------------------------------------------------------------
# 4. Stencils on partial coverage.
#
# Adjacency is a property of the complete level, so the ring table is built
# from the globe's neighbours and then filtered through `localindex`: a
# neighbour the crop does not own has no index, and is dropped rather than
# padded. Cells that keep all eight are the interior of the crop.
# --------------------------------------------------------------------------

rings = [Int[p for p in (DGG.localindex(grid, nb)
                         for nb in DGG.neighbors(globe, DGG.cellindex(grid, i)))
             if p !== nothing]
         for i in 1:DGG.ncells(grid)]

stencil(f) = [f(data[i], data[rings[i]]) for i in 1:DGG.ncells(grid)]
smoothed = stencil((c, nbs) -> mean(vcat(c, nbs)))
laplacian = stencil((c, nbs) -> isempty(nbs) ? zero(c) : mean(nbs) - c)

full = count(==(8), length.(rings))
println()
println("  field var $(round(var(data); digits=4)) -> smoothed $(round(var(smoothed); digits=4))")
println("  cells with a full 8-neighbourhood: $full / $(length(rings))")
check("smoothing reduces variance", var(smoothed) < var(data))
check("the field is smooth, so its Laplacian is small",
    maximum(abs, laplacian) < 1.0; detail="max |lap| = $(round(maximum(abs, laplacian); digits=5))")
check("only crop-boundary cells lose neighbours", 0 < full < length(rings))

# `neighbors` on the crop answers the same question directly. Subset adjacency
# is the complete level's adjacency clipped to membership, so the call uses
# HEALPix's fast path and drops cells outside the crop. `DGG.adjacency(grid)`
# returns the whole table in one call: its rows are these same clipped rings,
# in the same counter-clockwise order.
sample = 1:50:DGG.ncells(grid)
check("neighbors(crop, c) == the filtered globe ring",
    all(sort(Int[DGG.localindex(grid, nb)
                 for nb in DGG.neighbors(grid, DGG.cellindex(grid, i))]) == sort(rings[i])
        for i in sample); detail="sampled $(length(sample)) cells")

# --------------------------------------------------------------------------
# 5. The same three moves on every system.
# --------------------------------------------------------------------------

println()
println("  system              crop cells   within/touching   full-ring cells")
for sys in (DGG.systems()..., DGG.AuthalicSystem(DGG.HEALPixSystem()))
    base = sys isa DGG.AuthalicSystem ? parent(sys) : sys
    l = base isa Union{DGG.H3System,DGG.IGeo7System} ? 5 : 6
    g = DGG.PartialGrid(sys, l, DGG.query(sys, DGG.Intersects(WINDOW); level=l))
    complete = DGG.levelgrid(sys, l)
    within = zone_indices(g, DGG.Within(zones.France))
    touching = zone_indices(g, DGG.Intersects(zones.France))
    degrees = [count(!isnothing, (DGG.localindex(g, nb)
                                 for nb in DGG.neighbors(complete, DGG.cellindex(g, i))))
               for i in 1:DGG.ncells(g)]
    name = sys isa DGG.AuthalicSystem ?
           "Authalic($(nameof(typeof(base))))" : string(nameof(typeof(sys)))
    println("  ", rpad(name, 20), lpad(DGG.ncells(g), 10), lpad("$(length(within))/$(length(touching))", 18),
        lpad(count(==(maximum(degrees)), degrees), 18))
    check("$name: crop, zonal and stencil all work",
        DGG.ncells(g) > 0 && issubset(within, touching) && maximum(degrees) > 0)
end

println()
note("call site, verbatim:  grid = DGG.PartialGrid(sys, level, DGG.query(sys, DGG.Intersects(window); level))")
note("zonal = query + localindex; stencil = neighbors + localindex")

println()
println(FAILURES[] == 0 ? "ALL CHECKS PASSED" : "$(FAILURES[]) CHECK(S) FAILED")
exit(FAILURES[] == 0 ? 0 : 1)
