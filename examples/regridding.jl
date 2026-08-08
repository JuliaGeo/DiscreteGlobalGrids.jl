# Demo: conservative regridding from a full-globe DGGS to a lon/lat grid.
#
# The point of the generic grid family: the call site names a system exactly
# once, in the grid value. No per-system tree type, no per-system treeify, no
# per-system polygon helper — `Regridder` picks all of that up through
# `treeify(::DGGSGrid)`.
#
# Environment: this demo needs nothing beyond DiscreteGlobalGrids and its own
# dependencies, so it runs in the repository's project environment:
#
#     julia -t 4 --project=. examples/regridding.jl
#
# It ends in PASS/FAIL assertions and exits non-zero if any of them fail.
using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
import ConservativeRegridding as CR
import GeometryOps as GO

const FAILURES = Ref(0)
function check(name, ok; detail = "")
    ok || (FAILURES[] += 1)
    println(ok ? "PASS  " : "FAIL  ", rpad(name, 52), detail)
    return ok
end
note(text) = println("      ", text)

# Fixed 5-degree lon/lat destination, as a matrix of unit-sphere corner points.
# (The tuple-of-vectors `RegularGrid` path keeps its cells in geographic
# degrees, which the spherical clipper cannot mix with unit-sphere DGGS rings.)
const TO_SPHERE = GO.UnitSpherical.UnitSphereFromGeographic()
const DST_LON = collect(range(0, 360; length = 73))
const DST_LAT = collect(range(-90, 90; length = 37))
const DST = [TO_SPHERE((x, y)) for x in DST_LON, y in DST_LAT]
# Destination cells as polygons, for the tree-free reference intersections.
const DST_CELLS = vec(collect(getcell(treeify(DST))))
const INTERSECT = CR.DefaultIntersectionOperator(GO.Spherical())

println("="^78)
println("regridding.jl — Regridder(lon/lat $(length(DST_LON)-1)x$(length(DST_LAT)-1), DGGS globe)")
println("julia $(VERSION)  threads=$(Threads.nthreads())")
println("="^78)

"""
Structural and conservation checks on one DGGS source grid, tree-free where
possible: the intersection matrix is verified against cell areas and against
directly clipped polygons, never against another tree.
"""
function verify(label, system, level)
    # THE CALL SITE — the singleton is the only system-specific token:
    src = DGGSGrid(system, level)
    regridder = CR.Regridder(DST, src)

    A = regridder.intersections
    n_src = Int(DGG.num_cells(system, level))
    check("$label: matrix is (dst cells, src cells)",
          size(A) == (length(DST_CELLS), n_src);
          detail = "$(size(A)), nnz=$(length(A.nzval))")

    # The 5-degree grid tiles the whole sphere, so every source cell is fully
    # covered: its column of intersection areas sums to its own area. That is
    # the conservation property, read off the matrix.
    column_sums = vec(sum(A; dims = 1))
    col_err = maximum(abs.(column_sums .- regridder.src_areas) ./ regridder.src_areas)
    check("$label: column sums == source cell areas", col_err <= 1e-10;
          detail = "max rel err $col_err")
    row_sums = vec(sum(A; dims = 2))
    row_err = maximum(abs.(row_sums .- regridder.dst_areas) ./ maximum(regridder.dst_areas))
    check("$label: row sums == destination cell areas", row_err <= 1e-10;
          detail = "max rel err $row_err")
    # Source and destination tile the same sphere, so the two area budgets
    # agree. (`Spherical()` measures on the Earth-radius sphere, in m^2.)
    total_err = abs(sum(regridder.src_areas) - sum(regridder.dst_areas)) /
                sum(regridder.dst_areas)
    check("$label: src and dst cover the same sphere", total_err <= 1e-12;
          detail = "rel err $total_err over $(round(sum(regridder.src_areas); sigdigits = 6)) m^2")

    # Spot-check the widest column against intersection areas computed straight
    # from the source cell's polygon — no tree involved on either side.
    j = argmax(diff(A.colptr))
    id = ordinal_to_cell(system, level, j)
    polygon = cell_polygon_unitsphere(system, level, id)
    direct = [INTERSECT(cell, polygon) for cell in DST_CELLS]
    column = Vector(A[:, j])
    coldiff = maximum(abs.(direct .- column) ./ max.(abs.(direct), 1.0))
    check("$label: column $j == direct intersection areas", coldiff <= 1e-12;
          detail = "$(count(!iszero, column)) nonzero dst cells, max rel diff $coldiff")

    # `getcell` on the tree is the same polygon the kernel reports for that
    # ordinal — the leaf index space is the dense ordinal one.
    tree = treeify(src)
    check("$label: tree leaf $j is the ordinal-$j cell",
          getcell(tree, j) == polygon && ncells(tree) == n_src)
    return regridder
end

healpix = verify("HEALPix 4", HEALPixDGGS(), 4)
println()
igeo7 = verify("IGEO7 3", IGEO7DGGS(), 3)

# --------------------------------------------------------------------------
# Regrid a smooth analytic field.
#
# f = 1 + z + 3xy is a degree-2 spherical harmonic combination on the unit
# sphere. Sampling it at cell centers and regridding to the lon/lat grid should
# reproduce it at the destination cells to discretization order, and must
# preserve the area-weighted mean exactly — that is what "conservative" means.
# --------------------------------------------------------------------------

field(p) = 1.0 + p[3] + 3.0 * p[1] * p[2]
const ANALYTIC = vec([field(TO_SPHERE(((DST_LON[i] + DST_LON[i + 1]) / 2,
                                       (DST_LAT[j] + DST_LAT[j + 1]) / 2)))
                      for i in 1:(length(DST_LON) - 1), j in 1:(length(DST_LAT) - 1)])

"Sample `field` at every cell center of `system` at `level`, in ordinal order."
sample(system, level) = [field(DGG.cell_center(system, level, ordinal_to_cell(system, level, i)))
                         for i in 1:Int(DGG.num_cells(system, level))]

function regrid_field(regridder, system, level)
    destination = zeros(length(regridder.dst_areas))
    source = sample(system, level)
    CR.regrid!(destination, regridder, source)
    return source, destination
end

println()
println("  src level    cells    max |regridded - analytic|")
errors = Float64[]
for level in 4:6
    r = level == 4 ? healpix : CR.Regridder(DST, DGGSGrid(HEALPixDGGS(), level))
    source, destination = regrid_field(r, HEALPixDGGS(), level)
    push!(errors, maximum(abs.(destination .- ANALYTIC)))
    println("  $level        $(lpad(12 * 4^level, 9))    $(errors[end])")
    if level == 4
        src_mean = sum(source .* r.src_areas) / sum(r.src_areas)
        dst_mean = sum(destination .* r.dst_areas) / sum(r.dst_areas)
        check("area-weighted mean is preserved",
              abs(src_mean - dst_mean) < 1e-13;
              detail = "src $src_mean vs dst $dst_mean")
    end
end
check("field error shrinks with source resolution", issorted(errors; rev = true);
      detail = "$(round(errors[1] / errors[end]; digits = 1))x from level 4 to 6 " *
               "(field range $(round(maximum(ANALYTIC) - minimum(ANALYTIC); digits = 3)))")

# Two different DGGS, same destination, same field: the regridded results are
# two discretizations of one function and must agree to that order.
_, healpix_dst = regrid_field(healpix, HEALPixDGGS(), 4)
_, igeo7_dst = regrid_field(igeo7, IGEO7DGGS(), 3)
cross = maximum(abs.(healpix_dst .- igeo7_dst))
check("HEALPix 4 and IGEO7 3 agree on the same field", cross < 0.1;
      detail = "max |HEALPix - IGEO7| = $(round(cross; digits = 5)), " *
               "each within $(round(max(errors[1], maximum(abs.(igeo7_dst .- ANALYTIC))); digits = 5)) of analytic")

println()
note("call site, verbatim:  src = DGGSGrid(HEALPixDGGS(), 4); CR.Regridder(DST, src)")
note("swap the singleton for H3DGGS()/IGEO7DGGS() and nothing else changes")

println()
println(FAILURES[] == 0 ? "ALL CHECKS PASSED" : "$(FAILURES[]) CHECK(S) FAILED")
exit(FAILURES[] == 0 ? 0 : 1)
