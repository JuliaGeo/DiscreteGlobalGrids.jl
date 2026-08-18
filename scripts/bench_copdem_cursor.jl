# Compare source-tree strategies for one GLO-90 tile regridded to IGEO7 and HEALPix.
#
#     julia -t auto --project=docs scripts/bench_copdem_cursor.jl
#
# No network or raster data. `COPDEM_ROWS=n` limits the source to its northernmost
# `n` raster rows. Timings are machine-local.

import DiscreteGlobalGrids as DGG
import ConservativeRegridding as CR
import GeometryOps as GO
import GeoInterface as GI
const CD = DGG.CopernicusDEM

const MANIFOLD = GO.Spherical(; radius=1.0)
const SYS = DGG.CopernicusDEMSystem(90)
const TILE = CD.tilecell(SYS, 50, 6)
const NCOLS = Int(CD.ncols_at(SYS, 50))
const NROWS = Int(CD.lat_intervals(SYS))
const ROWS = clamp(parse(Int, get(ENV, "COPDEM_ROWS", string(NROWS))), 1, NROWS)

rowband(rows) = rows == NROWS ? DGG.subtree(SYS, TILE, 1) :
                DGG.PartialGrid(SYS, 1, [DGG.LevelIndex(1, k) for k in
                                         CD.pixelcell(SYS, TILE, 0, 0).index .+ (0:(rows * NCOLS - 1))])

# Densify the footprint and pad it by one pixel for bowed geodesic edges.
function footprint(chunk)
    w1, e1, s1, n1 = CD.cell_box(SYS, DGG.cellindex(chunk, 1))
    w2, e2, s2, n2 = CD.cell_box(SYS, DGG.cellindex(chunk, DGG.ncells(chunk)))
    w, e, s, n = min(w1, w2), max(e1, e2), min(s1, s2), max(n1, n2)
    w, e = w - (e - w) / NCOLS, e + (e - w) / NCOLS
    s, n = s - 1 / NROWS, n + 1 / NROWS
    steps = max(1, ceil(Int, (e - w) * 64))
    pts = [(w + (e - w) * k / steps, s) for k in 0:steps]
    append!(pts, [(w + (e - w) * k / steps, n) for k in steps:-1:0])
    push!(pts, pts[1])
    return GI.Polygon([GI.LinearRing(pts)])
end

# Destination level with closest mean cell area on a log scale.
matched(dstsys, chunk) = argmin(
    l -> abs(log(4π / DGG.ncells(dstsys, l) /
                 DGG.cell_area(DGG.levelgrid(SYS, 1), DGG.cellindex(chunk, 1)))),
    DGG.levels(dstsys))

const TREES = ("generic `HierarchicalGridCursor`" => g -> DGG.HierarchicalGridCursor(g),
    "`Blocked{3}` — 9 children a node" => g -> CD.BlockCursor(g; strategy=CD.Blocked{3}()),
    "`Blocked{2}` — 4 children a node" => g -> CD.BlockCursor(g; strategy=CD.Blocked{2}()),
    "`Bisected`   — 2 children a node" => g -> CD.BlockCursor(g; strategy=CD.Bisected()))

chunk = rowband(ROWS)
println("julia $(VERSION), $(Threads.nthreads()) threads, $(DGG.ncells(chunk)) source pixels")
times = Dict{String,Vector{Float64}}()
heads = String[]
for (name, dstsys) in (("IGEO7", DGG.IGeo7System()), ("HEALPix", DGG.HEALPixSystem()))
    L = matched(dstsys, chunk)
    push!(heads, "onto $name $L")
    dst = DGG.PartialGrid(DGG.CellVector(
        DGG.query(dstsys, DGG.MultiOrderCoverage(footprint(chunk)); level=L); level=L))
    warm = rowband(1)
    for (_, tree) in TREES     # every cursor type compiled before anything is timed
        CR.Regridder(MANIFOLD, dst, tree(warm))
    end
    println("\n  onto $name level $L: $(DGG.ncells(dst)) cells")
    reference = nothing
    for (label, tree) in TREES
        t = @elapsed r = CR.Regridder(MANIFOLD, dst, tree(chunk))
        agree = reference === nothing ? "(reference)" : string(r.intersections == reference)
        reference === nothing && (reference = r.intersections)
        push!(get!(times, label, Float64[]), t)
        println("    ", rpad(label, 34), lpad(round(t; digits=1), 7), " s   nnz=",
            length(r.intersections.nzval), "  identical=", agree)
    end
end

println("\n| source tree                        | ", rpad(heads[1], 13), " | ", rpad(heads[2], 15), " |")
println("|:-----------------------------------|--------------:|----------------:|")
for (label, _) in TREES
    println("| ", rpad(label, 34), " | ", lpad("$(round(times[label][1]; digits=1)) s", 13),
        " | ", lpad("$(round(times[label][2]; digits=1)) s", 15), " |")
end
