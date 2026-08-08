# Demo: chunked data — one spatial tree per data chunk (H3).
#
# A datacube stored as H3 res-6 cells grouped by their res-2 ancestor wants a
# tree per chunk: built in O(chunk), indexed 1:chunk_size so the tree's leaf
# indices are positions in the chunk's data array, and never reaching for a cell
# the chunk does not own. `subtree_grid` is that constructor.
#
# Environment: this demo needs nothing beyond DiscreteGlobalGrids and its own
# dependencies, so it runs in the repository's project environment:
#
#     julia -t 4 --project=. examples/chunked_h3.jl
#
# It ends in PASS/FAIL assertions and exits non-zero if any of them fail.
using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
const H3Native = DGG.H3.H3Native
import GeometryOps as GO
import GeometryOps: SpatialTreeInterface as STI

const FAILURES = Ref(0)
function check(name, ok; detail = "")
    ok || (FAILURES[] += 1)
    println(ok ? "PASS  " : "FAIL  ", rpad(name, 56), detail)
    return ok
end
note(text) = println("      ", text)

# Median-of-5 wall time and the allocated bytes of one call.
timed(f) = (sort([(@timed f()).time for _ in 1:5])[3], (@timed f()).bytes)

const CHUNK_LEVEL = 2
const LEAF_LEVEL = 6
# Fixed res-2 hexagon over the western Alps — nothing about the demo depends on
# which cell it is, only that it is the same one on every run.
const CHUNK = H3Native.lonlat_to_cell(10.0, 45.0, CHUNK_LEVEL)

println("="^78)
println("chunked_h3.jl — H3 res-$CHUNK_LEVEL chunk 0x", string(CHUNK; base = 16),
        ", leaves at res $LEAF_LEVEL")
println("julia $(VERSION)  threads=$(Threads.nthreads())")
println("="^78)

# --------------------------------------------------------------------------
# The call site: two lines from "which chunk" to "a spatial tree".
#
# Both levels are keywords, so the call cannot be mis-ordered: `root_id` is the
# only positional argument and the two `Integer` levels name themselves.
# --------------------------------------------------------------------------

grid = subtree_grid(H3DGGS(), CHUNK; root_level = CHUNK_LEVEL, leaf_level = LEAF_LEVEL)
tree = treeify(grid)

# --------------------------------------------------------------------------
# 1. The tree holds exactly the chunk.
# --------------------------------------------------------------------------

expected = subtree_leaf_count(H3DGGS(), CHUNK_LEVEL, CHUNK, LEAF_LEVEL)
check("leaf count == subtree_leaf_count", ncells(tree) == expected;
      detail = "$(ncells(tree)) leaves")
check("grid.ids holds only chunk descendants",
      all(cell_parent(H3DGGS(), LEAF_LEVEL, id, CHUNK_LEVEL) == CHUNK for id in grid.ids))
check("ids ascending and distinct", issorted(grid.ids; lt = (<=)))

# The root cursor stands for the chunk cell itself, not for the whole sphere —
# `node_level` / `node_id` are the public way to ask a node which cell it is.
check("tree roots at the chunk cell, not the sphere",
      (node_level(tree), node_id(tree)) == (CHUNK_LEVEL, CHUNK);
      detail = "node_level=$(node_level(tree)) node_id=0x$(string(node_id(tree); base = 16))")

# The tree's leaf index space is 1:n over the *chunk*, so a per-chunk data
# array indexes straight through it — this is what a global grid cannot give.
leaf_indices = Int[]
function walk_leaves!(out, node)
    if STI.isleaf(node)
        append!(out, first.(STI.child_indices_extents(node)))
    else
        for child in STI.getchild(node)
            walk_leaves!(out, child)
        end
    end
    return out
end
walk_leaves!(leaf_indices, tree)
check("leaf indices are exactly 1:$(expected)", sort(leaf_indices) == collect(1:expected))
check("getcell(tree, i) is the polygon of grid.ids[i]",
      all(getcell(tree, i) == cell_polygon_unitsphere(H3DGGS(), LEAF_LEVEL, grid.ids[i])
          for i in (1, 2, 1200, expected)))

# --------------------------------------------------------------------------
# 2. Cap queries — one small cap inside the chunk, one far outside.
#
# `intersects_cap(cap)` is the query predicate; node extents on a DGGS tree are
# always spherical caps. A bounding-volume query is conservative by
# construction: it answers with the cells whose *caps* survive the descent,
# which is a superset of the cells that truly intersect and a subset of the
# cells a naive per-cell cap test accepts. Both bounds are asserted, and every
# cell in the gap is shown not to touch the cap at all.
# --------------------------------------------------------------------------

const TO_SPHERE = GO.UnitSpherical.UnitSphereFromGeographic()
distance = GO.UnitSpherical.spherical_distance

inside_cap = GO.UnitSpherical.SphericalCap(DGG.cell_center(H3DGGS(), CHUNK_LEVEL, CHUNK), 0.01)
hits = STI.query(tree, intersects_cap(inside_cap))

# Naive brute force: the per-cell caps the tree itself reports at its leaves.
candidates = [i for i in eachindex(grid.ids)
              if intersects_cap(inside_cap, cell_cap(H3DGGS(), LEAF_LEVEL, grid.ids[i]))]
# Sound subset: a cell whose center is inside the cap certainly intersects it.
certain = [i for i in eachindex(grid.ids)
           if distance(inside_cap.point, DGG.cell_center(H3DGGS(), LEAF_LEVEL, grid.ids[i])) <=
              inside_cap.radius]

check("in-chunk query: no false negatives", issubset(certain, hits);
      detail = "$(length(certain)) certain / $(length(hits)) hits")
check("in-chunk query: hits subset the per-cell cap test", issubset(hits, candidates);
      detail = "$(length(hits)) hits / $(length(candidates)) candidates")
dropped = setdiff(candidates, hits)
nearest(i) = minimum(distance(inside_cap.point, v)
                     for v in DGG.cell_boundary(H3DGGS(), LEAF_LEVEL, grid.ids[i]))
check("dropped cells do not touch the query cap",
      all(nearest(i) > inside_cap.radius for i in dropped);
      detail = "$(length(dropped)) dropped; nearest boundary " *
               "$(isempty(dropped) ? "-" : round(minimum(nearest, dropped); digits = 5)) rad " *
               "vs cap radius $(inside_cap.radius)")

outside_cap = GO.UnitSpherical.SphericalCap(TO_SPHERE((-150.0, -40.0)), 0.01)
outside_hits = STI.query(tree, intersects_cap(outside_cap))
check("out-of-chunk query returns nothing", isempty(outside_hits))

# A subtree-rooted chunk of 2,401 leaves is past `SUBTREE_CAP_EXACT_LIMIT`, so
# its root extent is the O(1) `subtree_cap` (the chunk cell's own inflated cap)
# rather than an exact union cap, which at that size gives up and returns the
# whole sphere. So the chunk carries a tight O(1) bound.
root_extent = STI.node_extent(tree)
chunk_cap = subtree_cap(H3DGGS(), CHUNK_LEVEL, CHUNK, LEAF_LEVEL)
check("root extent is the O(1) chunk cap, not the whole sphere",
      root_extent == chunk_cap && root_extent != DGG.full_sphere_extent();
      detail = "radius $(round(root_extent.radius; digits = 4)) rad")
check("one cap test against it settles the far query",
      !intersects_cap(outside_cap, root_extent))

# Worth knowing before wiring a hot query loop: `depth_first_search` (hence
# `STI.query`) starts at the root's *children* and never tests the root's own
# extent, so a far-away query still pays one row of child extents — for this
# chunk, seven `cells_cap` folds over 343 cell boundaries each. Guarding the
# call with the root cap is a caller-side decision the tree makes cheap.
outside_time, outside_bytes = timed(() -> STI.query(tree, intersects_cap(outside_cap)))
guard_time, guard_bytes = timed(() -> intersects_cap(outside_cap, STI.node_extent(tree)))
note("far query, unguarded: $(round(outside_time * 1e3; digits = 2)) ms / " *
     "$(round(outside_bytes / 1024; digits = 0)) KiB " *
     "(STI.query descends from the root's children)")
note("far query, root cap first: $(round(guard_time * 1e6; digits = 2)) us / " *
     "$(guard_bytes) bytes — $(round(Int, outside_time / guard_time))x cheaper")

# --------------------------------------------------------------------------
# 3. O(chunk), not O(globe).
#
# The same construction at one level coarser covers 7x the cells; time and
# allocations follow the chunk, not the 14.1M-cell res-6 globe.
# --------------------------------------------------------------------------

coarse = cell_parent(H3DGGS(), CHUNK_LEVEL, CHUNK, CHUNK_LEVEL - 1)
build_fine = () -> subtree_grid(H3DGGS(), CHUNK;
                                root_level = CHUNK_LEVEL, leaf_level = LEAF_LEVEL)
build_coarse = () -> subtree_grid(H3DGGS(), coarse;
                                  root_level = CHUNK_LEVEL - 1, leaf_level = LEAF_LEVEL)
build_fine(); build_coarse()                      # warm up

fine_time, fine_bytes = timed(build_fine)
coarse_time, coarse_bytes = timed(build_coarse)
n_fine, n_coarse = length(grid.ids), length(build_coarse().ids)
n_globe = DGG.num_cells(H3DGGS(), LEAF_LEVEL)

println()
println("  chunk        cells        build        bytes")
println("  res $CHUNK_LEVEL      $(lpad(n_fine, 9))   $(lpad(round(fine_time * 1e6; digits = 2), 8)) us   $(lpad(fine_bytes, 10))")
println("  res $(CHUNK_LEVEL - 1)      $(lpad(n_coarse, 9))   $(lpad(round(coarse_time * 1e6; digits = 2), 8)) us   $(lpad(coarse_bytes, 10))")
println("  ratio        $(lpad(round(n_coarse / n_fine; digits = 2), 9))   $(lpad(round(coarse_time / fine_time; digits = 2), 11))   $(lpad(round(coarse_bytes / fine_bytes; digits = 2), 10))")
println("  globe        $(lpad(n_globe, 9))   ($(round(n_globe / n_fine; digits = 0))x the res-$CHUNK_LEVEL chunk)")

check("cell count scales 7x with one coarser level", n_coarse == 7 * n_fine)
check("build time scales with the chunk, not the globe",
      2.0 <= coarse_time / fine_time <= 20.0;
      detail = "$(round(coarse_time / fine_time; digits = 2))x for 7x the cells")
check("allocations scale with the chunk",
      2.0 <= coarse_bytes / fine_bytes <= 20.0;
      detail = "$(round(coarse_bytes / fine_bytes; digits = 2))x for 7x the cells")

# --------------------------------------------------------------------------
# 4. The whole-globe alternatives, for contrast.
#
# `DGGSPartialGrid` over every res-6 id is the same tree shape but has to
# materialize 14.1M ids first; `DGGSGrid` is O(1) to build but numbers its
# leaves globally, so a chunk's data array no longer indexes through it.
# --------------------------------------------------------------------------

globe_stats = @timed sort!(reduce(vcat, H3Native.cell_to_children.(H3Native.res0_cells(), LEAF_LEVEL)))
globe_ids = globe_stats.value
globe_grid = DGGSPartialGrid(H3DGGS(), LEAF_LEVEL, globe_ids)
dense_tree = treeify(DGGSGrid(H3DGGS(), LEAF_LEVEL))

println()
note("whole-globe DGGSPartialGrid: $(length(globe_ids)) ids, " *
     "$(round(globe_stats.time; digits = 3)) s / $(round(globe_stats.bytes / 2^20; digits = 1)) MiB to materialize " *
     "($(round(Int, globe_stats.time / fine_time))x the chunk's build time)")
note("whole-globe DGGSGrid: O(1) to build, but ncells = $(ncells(dense_tree)) — " *
     "leaf i is a global ordinal, not chunk position i")

check("chunk leaf 1 and globe leaf 1 are different cells",
      grid.ids[1] != globe_grid.ids[1];
      detail = "chunk id 0x$(string(grid.ids[1]; base = 16)) vs globe id 0x$(string(globe_grid.ids[1]; base = 16))")
check("chunk indices are chunk-local",
      getcell(tree, 1) == cell_polygon_unitsphere(H3DGGS(), LEAF_LEVEL, grid.ids[1]) &&
      getcell(tree, 1) != getcell(dense_tree, 1))

println()
println(FAILURES[] == 0 ? "ALL CHECKS PASSED" : "$(FAILURES[]) CHECK(S) FAILED")
exit(FAILURES[] == 0 ? 0 : 1)
