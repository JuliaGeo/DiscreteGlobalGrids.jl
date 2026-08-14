# Demo: conservative regridding from a full-globe DGGS to a lon/lat grid.
#
# The point of the redesigned grid interface: the call site names a system
# exactly once, in the grid value. No per-system tree type, no per-system
# `treeify`, no per-system polygon helper — `ConservativeRegridding.Regridder`
# picks all of that up through the interface, because `treeify`, `ncells` and
# `getcell` are `ConservativeRegridding.Trees`' own bindings extended for every
# `AbstractGrid`.
#
#     src = DGG.levelgrid(DGG.HEALPixSystem(), 4)
#     regridder = CR.Regridder(MANIFOLD, DST, src)
#
# Swap `HEALPixSystem()` for `IGeo7System()` or `H3System()` and nothing else
# changes. That is the whole claim this file checks.
#
# Environment: needs nothing beyond DiscreteGlobalGrids and its dependencies:
#
#     julia -t 4 --project=. examples/regridding.jl
#
# It ends in PASS/FAIL assertions and exits non-zero if any of them fail.

import DiscreteGlobalGrids as DGG
import ConservativeRegridding as CR
import GeometryOps as GO

const FAILURES = Ref(0)
function check(name, ok; detail="")
    ok || (FAILURES[] += 1)
    println(ok ? "PASS  " : "FAIL  ", rpad(name, 52), detail)
    return ok
end
note(text) = println("      ", text)

# Fixed 5-degree lon/lat destination, as a matrix of unit-sphere corner points.
# (The tuple-of-vectors `RegularGrid` path keeps its cells in geographic
# degrees, which the spherical clipper cannot mix with unit-sphere DGGS rings.)
const TO_SPHERE = GO.UnitSpherical.UnitSphereFromGeographic()
const DST_LON = collect(range(0, 360; length=73))
const DST_LAT = collect(range(-90, 90; length=37))
const DST = [TO_SPHERE((x, y)) for x in DST_LON, y in DST_LAT]
# Destination cells as polygons, for the tree-free reference intersections.
const DST_CELLS = vec(collect(DGG.getcell(DGG.treeify(DST))))

# THE MANIFOLD, declared once. Every grid in this package computes on the UNIT
# sphere — `cell_boundary` returns `UnitSphericalPoint`s and `cell_area` returns
# steradians — so that is the manifold the regridder must work on. A bare
# `Matrix{UnitSphericalPoint}` carries no manifold of its own, and
# `best_manifold` guesses `Spherical()`, whose radius is the WGS84 mean radius;
# mixing the two is a factor of R^2 in every area, so ConservativeRegridding
# refuses rather than silently rescaling. Naming it here is the whole fix.
const MANIFOLD = GO.Spherical(; radius=1.0)
const INTERSECT = CR.DefaultIntersectionOperator(MANIFOLD)

println("="^78)
println("regridding.jl — Regridder(lon/lat $(length(DST_LON)-1)x$(length(DST_LAT)-1), DGGS globe)")
println("julia $(VERSION)  threads=$(Threads.nthreads())")
println("="^78)

"""
Structural and conservation checks on one DGGS source grid, tree-free where
possible: the intersection matrix is verified against cell areas and against
directly clipped polygons, never against another tree.
"""
function verify(label, sys, l)
    # THE CALL SITE — the singleton is the only system-specific token.
    src = DGG.levelgrid(sys, l)
    regridder = CR.Regridder(MANIFOLD, DST, src)

    A = regridder.intersections
    n_src = DGG.ncells(src)
    check("$label: matrix is (dst cells, src cells)",
        size(A) == (length(DST_CELLS), n_src);
        detail="$(size(A)), nnz=$(length(A.nzval))")

    # The 5-degree grid tiles the whole sphere, so every source cell is fully
    # covered: its column of intersection areas sums to its own area. That is
    # the conservation property, read off the matrix.
    column_sums = vec(sum(A; dims=1))
    col_err = maximum(abs.(column_sums .- regridder.src_areas) ./ regridder.src_areas)
    check("$label: column sums == source cell areas", col_err <= 1e-10;
        detail="max rel err $col_err")
    row_sums = vec(sum(A; dims=2))
    row_err = maximum(abs.(row_sums .- regridder.dst_areas) ./ maximum(regridder.dst_areas))
    check("$label: row sums == destination cell areas", row_err <= 1e-10;
        detail="max rel err $row_err")
    # Source and destination tile the same sphere, so the two area budgets
    # agree — and on the unit sphere that budget is 4pi steradians exactly.
    total_err = abs(sum(regridder.src_areas) - sum(regridder.dst_areas)) /
                sum(regridder.dst_areas)
    check("$label: src and dst cover the same sphere", total_err <= 1e-12;
        detail="rel err $total_err over $(round(sum(regridder.src_areas); sigdigits=8)) sr (4pi = $(round(4pi; sigdigits=8)))")

    # Spot-check the widest column against intersection areas computed straight
    # from the source cell's polygon — no tree involved on either side. Note
    # the position/identity split: `j` is a POSITION in the matrix, and
    # `cellindex` is what turns it into the cell's name.
    j = argmax(diff(A.colptr))
    c = DGG.cellindex(src, j)
    polygon = DGG.cell_polygon(src, c)
    direct = [INTERSECT(cell, polygon) for cell in DST_CELLS]
    column = Vector(A[:, j])
    coldiff = maximum(abs.(direct .- column) ./ max.(abs.(direct), 1.0))
    check("$label: column $j == direct intersection areas", coldiff <= 1e-12;
        detail="$(count(!iszero, column)) nonzero dst cells, max rel diff $coldiff")

    # `getcell` on the tree is the same polygon the grid reports for that
    # position — the leaf index space IS the dense position space.
    tree = DGG.treeify(src)
    check("$label: tree leaf $j is the position-$j cell",
        DGG.getcell(tree, j) == polygon && DGG.ncells(tree) == n_src)
    return regridder
end

healpix = verify("HEALPix 4", DGG.HEALPixSystem(), 4)
println()
igeo7 = verify("IGeo7 3", DGG.IGeo7System(), 3)
println()
h3 = verify("H3 2", DGG.H3System(), 2)

# --------------------------------------------------------------------------
# Regrid a smooth analytic field.
#
# f = 1 + z + 3xy is a degree-2 spherical harmonic combination on the unit
# sphere. Sampling it at cell centres and regridding to the lon/lat grid should
# reproduce it at the destination cells to discretisation order, and must
# preserve the area-weighted mean exactly — that is what "conservative" means.
# --------------------------------------------------------------------------

field(p) = 1.0 + p[3] + 3.0 * p[1] * p[2]
const ANALYTIC = vec([field(TO_SPHERE(((DST_LON[i] + DST_LON[i+1]) / 2,
    (DST_LAT[j] + DST_LAT[j+1]) / 2)))
                      for i in 1:(length(DST_LON)-1), j in 1:(length(DST_LAT)-1)])

"Sample `field` at every cell centroid of `grid`, in position order."
sample(grid) = [field(DGG.cell_centroid(grid, DGG.cellindex(grid, i)))
                for i in 1:DGG.ncells(grid)]

function regrid_field(regridder, grid)
    destination = zeros(length(regridder.dst_areas))
    source = sample(grid)
    CR.regrid!(destination, regridder, source)
    return source, destination
end

println()
println("  src level    cells    max |regridded - analytic|")
errors = Float64[]
for l in 4:6
    grid = DGG.levelgrid(DGG.HEALPixSystem(), l)
    r = l == 4 ? healpix : CR.Regridder(MANIFOLD, DST, grid)
    source, destination = regrid_field(r, grid)
    push!(errors, maximum(abs.(destination .- ANALYTIC)))
    println("  $l        $(lpad(DGG.ncells(grid), 9))    $(errors[end])")
    if l == 4
        src_mean = sum(source .* r.src_areas) / sum(r.src_areas)
        dst_mean = sum(destination .* r.dst_areas) / sum(r.dst_areas)
        check("area-weighted mean is preserved",
            abs(src_mean - dst_mean) < 1e-13;
            detail="src $src_mean vs dst $dst_mean")
    end
end
check("field error shrinks with source resolution", issorted(errors; rev=true);
    detail="$(round(errors[1] / errors[end]; digits=1))x from level 4 to 6 " *
           "(field range $(round(maximum(ANALYTIC) - minimum(ANALYTIC); digits=3)))")

# Three different DGGS, same destination, same field: the regridded results are
# three discretisations of one function and must agree to that order.
_, healpix_dst = regrid_field(healpix, DGG.levelgrid(DGG.HEALPixSystem(), 4))
_, igeo7_dst = regrid_field(igeo7, DGG.levelgrid(DGG.IGeo7System(), 3))
_, h3_dst = regrid_field(h3, DGG.levelgrid(DGG.H3System(), 2))
for (name, other) in (("IGeo7 3", igeo7_dst), ("H3 2", h3_dst))
    cross = maximum(abs.(healpix_dst .- other))
    check("HEALPix 4 and $name agree on the same field", cross < 0.2;
        detail="max |HEALPix - $name| = $(round(cross; digits=5))")
end

# --------------------------------------------------------------------------
# The claim, swept: every registered system, plus an ellipsoidal wrap of one.
#
# `AuthalicSystem` re-reads a system's geometry at geodetic latitude. It is a
# system like any other, so it goes through the same call site — which is the
# point of sweeping it here rather than treating it as a special case.
# --------------------------------------------------------------------------

# The wrapper is transparent to everything but geometry, so both the level
# choice and the label read through it.
base(sys) = sys isa DGG.AuthalicSystem ? parent(sys) : sys
demo_level(sys) = base(sys) isa DGG.H3System ? 2 :
                  base(sys) isa DGG.IGeo7System ? 3 : 4
label(sys) = sys isa DGG.AuthalicSystem ?
             "Authalic($(nameof(typeof(parent(sys)))))" : string(nameof(typeof(sys)))

println()
println("  system                     cells    column-sum rel err")
for sys in (DGG.systems()..., DGG.AuthalicSystem(DGG.IGeo7System()))
    l = demo_level(sys)
    grid = DGG.levelgrid(sys, l)
    r = CR.Regridder(MANIFOLD, DST, grid)
    err = maximum(abs.(vec(sum(r.intersections; dims=1)) .- r.src_areas) ./ r.src_areas)
    println("  ", rpad("$(label(sys)) $l", 23), lpad(DGG.ncells(grid), 9), "    ", err)
    check("$(label(sys)) $l: conservative", err <= 1e-10)
end

# --------------------------------------------------------------------------
# The other direction: the lon/lat grid as SOURCE and the DGGS as DESTINATION.
#
# Same two grids, same manifold, arguments swapped — and it is NOT the same
# matrix, because the intersection operator clips the source cell against the
# destination cell and is not symmetric in the two. On a `Spherical` manifold
# that operator is Sutherland-Hodgman, which clips the SUBJECT (the source
# cell) against the half-space of every edge of the CLIP WINDOW (the
# destination cell), and only a CONVEX clip window makes that the intersection.
#
# HEALPix, ISEA4R and A5 cells are curvilinear, so `cell_boundary` densifies
# each chart edge into eight great-circle segments — and a straight chart edge
# is a curve on the sphere, so two of a cell's sides bow inward and every
# densified vertex on them is a reflex vertex. Correct geometry, and it still
# tiles the sphere exactly (which is why the column sums above are exact); but
# as a clip window it is out of contract, and the destination cells lose 0.2%
# to 1.5% of their area — or, where a cell's first ring vertex lands on a
# source cell's boundary, gain a whole extra copy of themselves.
#
# So this loop asserts the law where it holds and PRINTS the shortfall where it
# does not, rather than pretending either. `test/systems/crosssystem/
# regridding_conservation.jl` carries the same law with the failing arms marked
# `@test_broken`, and names the upstream file and lines.
# --------------------------------------------------------------------------

println()
println("  destination direction      cells    row-sum rel err")
for sys in (DGG.systems()..., DGG.AuthalicSystem(DGG.IGeo7System()))
    l = demo_level(sys)
    grid = DGG.levelgrid(sys, l)
    r = CR.Regridder(MANIFOLD, grid, DST)   # DGGS is now the DESTINATION
    err = maximum(abs.(vec(sum(r.intersections; dims=2)) .- r.dst_areas) ./ r.dst_areas)
    ones_back = zeros(DGG.ncells(grid))
    CR.regrid!(ones_back, r, ones(length(r.src_areas)))
    println("  ", rpad("$(label(sys)) $l", 23), lpad(DGG.ncells(grid), 9), "    ", err)
    if err <= 1e-10
        check("$(label(sys)) $l: conservative onto", true)
    else
        note("$(label(sys)) $l: NOT conservative onto — regrid!(ones) spans " *
             "$(round.(extrema(ones_back); digits=3)); non-convex cell rings, " *
             "GeometryOps sutherland_hodgman.jl:306-324")
    end
end

println()
note("call site, verbatim:  src = DGG.levelgrid(DGG.HEALPixSystem(), 4); CR.Regridder(MANIFOLD, DST, src)")
note("swap the singleton for any of systems(), or wrap one in AuthalicSystem")
note("systems() lists them: " * join(string.(nameof.(typeof.(DGG.systems()))), ", "))

println()
println(FAILURES[] == 0 ? "ALL CHECKS PASSED" : "$(FAILURES[]) CHECK(S) FAILED")
exit(FAILURES[] == 0 ? 0 : 1)
