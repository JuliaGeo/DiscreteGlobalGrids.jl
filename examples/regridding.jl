# Demo: conservative regridding from a full-globe DGGS to a lon/lat grid.
#
# The point of the redesigned grid interface: the call site names a system
# exactly once, in the grid value. No per-system tree type, no per-system
# `treeify`, no per-system polygon helper — `DGGSpace` picks all of that up
# through the interface, because `treeify`, `ncells` and `getcell` are
# `ConservativeRegridding.Trees`' own bindings extended for every `AbstractGrid`.
#
#     src = DGG.levelgrid(DGG.HEALPixSystem(), 4)
#     plan = DGG.plan_regrid(values; from = DGG.DGGSpace(src), to = DST)
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
import GlobalRegridding as GR
import ConservativeRegridding as CR
import DimensionalData as DD
import GeometryOps as GO

const FAILURES = Ref(0)
function check(name, ok; detail="")
    ok || (FAILURES[] += 1)
    println(ok ? "PASS  " : "FAIL  ", rpad(name, 52), detail)
    return ok
end
note(text) = println("      ", text)

# Fixed 5-degree lon/lat destination. A `RasterGrid` is a regridding space over
# a dimensional raster's cells and wants the two axes alone: it keeps the cell
# EDGES and synthesizes polygons, caps and areas from them on demand. The axes
# are the cell centres, which midpoint back to exactly these edges.
const DST_LON = collect(range(0, 360; length=73))
const DST_LAT = collect(range(-90, 90; length=37))
midpoints(e) = (e[1:end-1] .+ e[2:end]) ./ 2
const DST = GR.RasterGrid(DD.X(DD.Sampled(midpoints(DST_LON))),
    DD.Y(DD.Sampled(midpoints(DST_LAT))))

"Every cell area of a regridding space, in index order."
areas(space) = [GR.cellarea(space, i) for i in 1:GR.ncells(space)]

# Destination cells as polygons, for the tree-free reference intersections, and
# as areas, for the conservation budget.
const DST_CELLS = [GR.getcell(DST, i) for i in 1:GR.ncells(DST)]
const DST_AREAS = areas(DST)

# The manifold is read off the space rather than declared. Every grid in this
# package computes on the unit sphere — `cell_boundary` returns
# `UnitSphericalPoint`s and `cell_area` returns steradians — and a `RasterGrid`
# is on that same one, which is why the two sides of a plan agree at all;
# `plan_regrid` refuses a pair that disagrees rather than rescaling every area
# by R^2 in silence. The reference intersections must clip where the weights did.
const INTERSECT = CR.DefaultIntersectionOperator(GR.manifold(DST))

println("="^78)
println("regridding.jl — regrid(lon/lat $(length(DST_LON)-1)x$(length(DST_LAT)-1), DGGS globe)")
println("julia $(VERSION)  threads=$(Threads.nthreads())")
println("="^78)

# The field every regrid below carries. f = 1 + z + 3xy is a degree-2 spherical
# harmonic combination on the unit sphere. Sampling it at cell centres and
# regridding to the lon/lat grid should reproduce it at the destination cells to
# discretisation order, and must preserve the area-weighted mean exactly — that
# is what "conservative" means.
const TO_SPHERE = GO.UnitSpherical.UnitSphereFromGeographic()
field(p) = 1.0 + p[3] + 3.0 * p[1] * p[2]
const ANALYTIC = vec([field(TO_SPHERE(((DST_LON[i] + DST_LON[i+1]) / 2,
    (DST_LAT[j] + DST_LAT[j+1]) / 2)))
                      for i in 1:(length(DST_LON)-1), j in 1:(length(DST_LAT)-1)])

"Sample `field` at every cell centroid of `grid`, in index order."
sample(grid) = [field(DGG.cell_centroid(grid, DGG.cellindex(grid, i)))
                for i in 1:DGG.ncells(grid)]

"""
The regridding operator from one DGGS onto `DST`, written once. Planning reads
no source data — weights are geometry — so the array names the shape the plan
will be applied to rather than being measured. `Extensive` leaves the apply a raw
conservative sum: `Weighted` would divide each destination by the area it was
handed, which is the very shortfall the last sweep below is looking for.
"""
planfrom(grid) = DGG.plan_regrid(sample(grid); from=DGG.DGGSpace(grid), to=DST,
    missingpolicy=GR.Extensive())

"""
Structural and conservation checks on one DGGS source grid, tree-free where
possible: the intersection matrix is verified against cell areas and against
directly clipped polygons, never against another tree.
"""
function verify(label, sys, l)
    # The system singleton is the only system-specific token at the call site.
    src = DGG.levelgrid(sys, l)
    plan = planfrom(src)

    # An eager plan holds one whole-domain block, and a whole-domain block's
    # chunk-local indices ARE cell indices — so its weights are the
    # intersection matrix itself, destination cells by source cells.
    A = plan.block.weights
    src_areas = areas(plan.src_space)
    n_src = DGG.ncells(src)
    check("$label: matrix is (dst cells, src cells)",
        size(A) == (length(DST_CELLS), n_src);
        detail="$(size(A)), nnz=$(length(A.nzval))")

    # The 5-degree grid tiles the whole sphere, so every source cell is fully
    # covered: its column of intersection areas sums to its own area. That is
    # the conservation property, read off the matrix.
    column_sums = vec(sum(A; dims=1))
    col_err = maximum(abs.(column_sums .- src_areas) ./ src_areas)
    check("$label: column sums == source cell areas", col_err <= 1e-10;
        detail="max rel err $col_err")
    row_sums = vec(sum(A; dims=2))
    row_err = maximum(abs.(row_sums .- DST_AREAS) ./ maximum(DST_AREAS))
    check("$label: row sums == destination cell areas", row_err <= 1e-10;
        detail="max rel err $row_err")
    # Source and destination tile the same sphere, so the two area budgets
    # agree — and on the unit sphere that budget is 4pi steradians exactly.
    total_err = abs(sum(src_areas) - sum(DST_AREAS)) / sum(DST_AREAS)
    check("$label: src and dst cover the same sphere", total_err <= 1e-12;
        detail="rel err $total_err over $(round(sum(src_areas); sigdigits=8)) sr (4pi = $(round(4pi; sigdigits=8)))")

    # Spot-check the widest column against intersection areas computed straight
    # from the source cell's polygon — no tree involved on either side. Note
    # the index/identity split: `j` is an INDEX in the matrix, and
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
    # index — the leaf index space IS the dense index space, which is why
    # a `DGGSpace` needs no permutation between the two.
    tree = DGG.treeify(src)
    check("$label: tree leaf $j is the index-$j cell",
        DGG.getcell(tree, j) == polygon && DGG.ncells(tree) == n_src)
    return plan
end

healpix = verify("HEALPix 4", DGG.HEALPixSystem(), 4)
println()
igeo7 = verify("IGeo7 3", DGG.IGeo7System(), 3)
println()
h3 = verify("H3 2", DGG.H3System(), 2)

# --------------------------------------------------------------------------
# Apply the plans. A plan takes no keyword arguments — it already carries the
# method, both spaces and the missing policy — and under `Extensive` a
# destination value is the field's integral over that cell, so the cell's own
# area is what turns it back into a mean.
# --------------------------------------------------------------------------

function regrid_field(plan, grid)
    source = sample(grid)
    return source, DGG.regrid(source, plan) ./ DST_AREAS
end

println()
println("  src level    cells    max |regridded - analytic|")
errors = Float64[]
for l in 4:6
    grid = DGG.levelgrid(DGG.HEALPixSystem(), l)
    plan = l == 4 ? healpix : planfrom(grid)
    source, destination = regrid_field(plan, grid)
    push!(errors, maximum(abs.(destination .- ANALYTIC)))
    println("  $l        $(lpad(DGG.ncells(grid), 9))    $(errors[end])")
    if l == 4
        src_areas = areas(plan.src_space)
        src_mean = sum(source .* src_areas) / sum(src_areas)
        dst_mean = sum(destination .* DST_AREAS) / sum(DST_AREAS)
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
    plan = planfrom(grid)
    src_areas = areas(plan.src_space)
    err = maximum(abs.(vec(sum(plan.block.weights; dims=1)) .- src_areas) ./ src_areas)
    println("  ", rpad("$(label(sys)) $l", 23), lpad(DGG.ncells(grid), 9), "    ", err)
    check("$(label(sys)) $l: conservative", err <= 1e-10)
end

# --------------------------------------------------------------------------
# The other direction: the lon/lat grid as SOURCE and the DGGS as DESTINATION.
# A grid is a `to` target as it stands, and the plan hands both resolved spaces
# back, so nothing is built twice to measure them.
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

const COVER = ones(GR.ncells(DST))

println()
println("  destination direction      cells    row-sum rel err")
for sys in (DGG.systems()..., DGG.AuthalicSystem(DGG.IGeo7System()))
    l = demo_level(sys)
    grid = DGG.levelgrid(sys, l)
    plan = DGG.plan_regrid(COVER; from=DST, to=grid, missingpolicy=GR.Extensive())
    dst_areas = areas(plan.dst_space)
    err = maximum(abs.(vec(sum(plan.block.weights; dims=2)) .- dst_areas) ./ dst_areas)
    # A field of ones regridded extensively is the area each destination cell
    # was handed; over the cell's own area that is its coverage, 1 exactly where
    # the law holds. `Weighted` divides by that same area and would answer 1
    # whatever the clipper did.
    ones_back = DGG.regrid(COVER, plan) ./ dst_areas
    println("  ", rpad("$(label(sys)) $l", 23), lpad(DGG.ncells(grid), 9), "    ", err)
    if err <= 1e-10
        check("$(label(sys)) $l: conservative onto", true)
    else
        note("$(label(sys)) $l: NOT conservative onto — regrid(ones) spans " *
             "$(round.(extrema(ones_back); digits=3)); non-convex cell rings, " *
             "GeometryOps sutherland_hodgman.jl:306-324")
    end
end

println()
note("call site, verbatim:  src = DGG.levelgrid(DGG.HEALPixSystem(), 4); DGG.plan_regrid(sample(src); from = DGG.DGGSpace(src), to = DST)")
note("swap the singleton for any of systems(), or wrap one in AuthalicSystem")
note("systems() lists them: " * join(string.(nameof.(typeof.(DGG.systems()))), ", "))

println()
println(FAILURES[] == 0 ? "ALL CHECKS PASSED" : "$(FAILURES[]) CHECK(S) FAILED")
exit(FAILURES[] == 0 ? 0 : 1)
