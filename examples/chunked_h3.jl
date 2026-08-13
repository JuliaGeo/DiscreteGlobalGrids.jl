# Demo: chunked data — one grid, one tree, per data chunk (H3).
#
# A datacube stored as H3 level-6 cells grouped by their level-2 ancestor wants
# a grid per chunk: built in O(1), numbered 1:ncells so a position is an index
# into the chunk's data array, and never naming a cell the chunk does not own.
#
#     grid = DGG.PartialGrid(DGG.H3System(), chunk, LEAF_LEVEL)
#
# is the whole constructor. On a system with sorted subtrees the ids are the
# lazy `descendant_range` window, so nothing is materialised at all.
#
# Environment: needs nothing beyond DiscreteGlobalGrids and its dependencies:
#
#     julia -t 4 --project=. examples/chunked_h3.jl
#
# It ends in PASS/FAIL assertions and exits non-zero if any of them fail.

import DiscreteGlobalGrids as DGG
import GeometryOps as GO
import Extents

const FAILURES = Ref(0)
function check(name, ok; detail="")
    ok || (FAILURES[] += 1)
    println(ok ? "PASS  " : "FAIL  ", rpad(name, 56), detail)
    return ok
end
note(text) = println("      ", text)

# Median-of-5 wall time and the allocated bytes of one call.
timed(f) = (sort([(@timed f()).time for _ in 1:5])[3], (@timed f()).bytes)

const SYS = DGG.H3System()
const CHUNK_LEVEL = 2
const LEAF_LEVEL = 6
# A fixed level-2 hexagon over the western Alps — nothing about the demo depends
# on which cell it is, only that it is the same one on every run.
const CHUNK = DGG.cellat(DGG.levelgrid(SYS, CHUNK_LEVEL), 10.0, 45.0)

println("="^78)
println("chunked_h3.jl — H3 level-$CHUNK_LEVEL chunk $CHUNK, leaves at level $LEAF_LEVEL")
println("julia $(VERSION)  threads=$(Threads.nthreads())")
println("="^78)

# --------------------------------------------------------------------------
# 1. The grid holds exactly the chunk, and numbers it from 1.
# --------------------------------------------------------------------------

grid = DGG.PartialGrid(SYS, CHUNK, LEAF_LEVEL)
tree = DGG.treeify(grid)

expected = length(DGG.descendant_range(SYS, CHUNK, LEAF_LEVEL))
check("ncells == the chunk's descendant range", DGG.ncells(grid) == expected;
    detail="$(DGG.ncells(grid)) cells")
check("every cell is a chunk descendant",
    all(DGG.ancestor(SYS, DGG.cellindex(grid, i), CHUNK_LEVEL) == CHUNK
        for i in 1:DGG.ncells(grid)))
check("positions are chunk-local and round trip",
    all(DGG.cellposition(grid, DGG.cellindex(grid, i)) == i
        for i in (1, 2, 1200, expected)))
check("the tree's leaf i is the grid's position i",
    DGG.ncells(tree) == expected &&
    all(DGG.getcell(tree, i) == DGG.cell_polygon(grid, DGG.cellindex(grid, i))
        for i in (1, 2, 1200, expected)))

# The globe at the same level is O(1) to build too, but its positions are global
# ordinals, so a per-chunk data array no longer indexes through it.
globe = DGG.levelgrid(SYS, LEAF_LEVEL)
check("chunk position 1 is not globe position 1",
    DGG.cellindex(grid, 1) != DGG.cellindex(globe, 1);
    detail="$(DGG.cellindex(grid, 1)) vs $(DGG.cellindex(globe, 1))")
note("globe at level $LEAF_LEVEL: $(DGG.ncells(globe)) cells " *
     "($(round(Int, DGG.ncells(globe) / expected))x the chunk)")

# --------------------------------------------------------------------------
# 2. Queries answer in the chunk's own index space.
#
# `query` returns typed ids; `cellposition` is what turns one into an index into
# the chunk's data array. A query whose target is nowhere near the chunk is
# pruned at the grid's root extent and comes back empty.
# --------------------------------------------------------------------------

near = Extents.Extent(X=(9.9, 10.1), Y=(44.9, 45.1))
far = Extents.Extent(X=(-150.5, -149.5), Y=(-40.5, -39.5))

hits = DGG.query(grid, DGG.Intersects(near))
positions = [DGG.cellposition(grid, c) for c in hits]
check("in-chunk query hits index the chunk", !isempty(hits) &&
                                             all(p -> p isa Int && 1 <= p <= expected, positions);
    detail="$(length(hits)) cells")
check("hits really do meet the target",
    all(DGG.cellindex(grid, p) in hits for p in positions))
check("out-of-chunk query returns nothing",
    isempty(DGG.query(grid, DGG.Intersects(far))))

# `cellat` is the point form of the same question, and answers `nothing` outside
# the chunk's coverage rather than pointing at some other chunk's cell.
check("cellat inside the chunk finds a cell",
    DGG.cellat(grid, 10.0, 45.0) !== nothing)
check("cellat outside the chunk is nothing",
    DGG.cellat(grid, -150.0, -40.0) === nothing)

# --------------------------------------------------------------------------
# 3. The rim, without the interior.
#
# Halo exchange between chunks needs the cells with a neighbour in another
# chunk. `subtree_border` answers that directly, and H3 walks it in O(rim).
# --------------------------------------------------------------------------

rim = DGG.subtree_border(SYS, CHUNK, LEAF_LEVEL)
check("rim is a strict subset of the chunk", 0 < length(rim) < expected;
    detail="$(length(rim)) of $expected cells ($(round(100 * length(rim) / expected; digits=1))%)")
check("every rim cell has a neighbour outside the chunk",
    all(any(DGG.ancestor(SYS, nb, CHUNK_LEVEL) != CHUNK
            for nb in DGG.neighbors(globe, c)) for c in rim))
check("rim and interior partition the chunk",
    length(rim) + length(DGG.subtree_interior(SYS, CHUNK, LEAF_LEVEL)) == expected)

# --------------------------------------------------------------------------
# 4. O(chunk), not O(globe).
#
# The same construction one level coarser covers 7x the cells. On H3 the ids are
# a lazy window over the level grid, so even the 7x chunk allocates nothing much.
# --------------------------------------------------------------------------

coarse = parent(SYS, CHUNK)
build_fine = () -> DGG.PartialGrid(SYS, CHUNK, LEAF_LEVEL)
build_coarse = () -> DGG.PartialGrid(SYS, coarse, LEAF_LEVEL)
build_fine();
build_coarse();                       # warm up

fine_time, fine_bytes = timed(build_fine)
coarse_time, coarse_bytes = timed(build_coarse)
n_coarse = DGG.ncells(build_coarse())

println()
println("  chunk        cells        build        bytes")
println("  level $CHUNK_LEVEL    $(lpad(expected, 9))   $(lpad(round(fine_time * 1e6; digits=2), 8)) us   $(lpad(fine_bytes, 10))")
println("  level $(CHUNK_LEVEL - 1)    $(lpad(n_coarse, 9))   $(lpad(round(coarse_time * 1e6; digits=2), 8)) us   $(lpad(coarse_bytes, 10))")
println("  globe        $(lpad(DGG.ncells(globe), 9))")

check("cell count scales 7x with one coarser level", n_coarse == 7 * expected)
check("build cost does not scale with the chunk", coarse_bytes <= 4 * fine_bytes;
    detail="$(fine_bytes) -> $(coarse_bytes) bytes for 7x the cells")

# --------------------------------------------------------------------------
# 5. The same two lines on every system.
#
# `PartialGrid(sys, cell, leaf)` is generic: sorted-subtree systems get the O(1)
# lazy window, and A5 — which does not claim sorted subtrees — falls back to
# materialising `descendants`. The call site does not change either way.
# --------------------------------------------------------------------------

println()
println("  system            chunk cells   ids materialised")
for sys in (DGG.systems()..., DGG.AuthalicSystem(DGG.H3System()))
    base = sys isa DGG.AuthalicSystem ? parent(sys) : sys
    root_level, leaf_level = 2, base isa Union{DGG.H3System,DGG.IGeo7System} ? 5 : 6
    chunk = DGG.cellat(DGG.levelgrid(sys, root_level), 10.0, 45.0)
    g = DGG.PartialGrid(sys, chunk, leaf_level)
    lazy = !DGG.has_sorted_subtrees(sys)
    name = sys isa DGG.AuthalicSystem ?
           "Authalic($(nameof(typeof(base))))" : string(nameof(typeof(sys)))
    println("  ", rpad(name, 18), lpad(DGG.ncells(g), 11), "   ", lazy ? "yes" : "no")
    check("$name: chunk grid is rooted and local",
        DGG.ncells(g) > 0 && DGG.cellposition(g, DGG.cellindex(g, 1)) == 1)
end

println()
note("call site, verbatim:  grid = DGG.PartialGrid(sys, chunk, LEAF_LEVEL)")
note("`treeify(grid)` is the tree; `cellposition(grid, c)` is the data index")

println()
println(FAILURES[] == 0 ? "ALL CHECKS PASSED" : "$(FAILURES[]) CHECK(S) FAILED")
exit(FAILURES[] == 0 ? 0 : 1)
