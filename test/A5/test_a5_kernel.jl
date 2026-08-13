module A5KernelTestSuite

# Tests for `src/A5/A5Kernel.jl`: the operations-kernel wiring of `A5DGGS`.
#
# A5 is the first real system on the generic cursor's NO-`descendant_range`
# path (`cell_parent` membership filtering into materialized index vectors),
# which until now only mock systems exercised — so the traversal testsets below
# are as much a test of `src/core/generic_cursor.jl` as of this wiring.
#
# It is also the first system whose hierarchy is not a fixed radix (12 roots ->
# 5 quintants -> aperture 4), so every count is checked against an independent
# enumeration rather than against arithmetic, and the mandatory CAP-VALIDATION
# sweep at the end is what justifies the raised `cell_cap_inflation`.
#
# The suite lives in its own module: the generic vocabulary the systems share
# (`cell_center`, `cell_boundary`, ...) must not collide with a sibling's.

using Test
using Printf
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: A5DGGS, DGGSGrid, DGGSPartialGrid, subtree_grid, NotPortedError,
    root_count, radix, max_level, treeify, ncells, getcell, node_level, node_id,
    intersects_cap
using DiscreteGlobalGrids.A5.A5Lookups: A5Lookup, A5Cells
import GeometryOps as GO
import GeoInterface as GI
import GeometryOps: SpatialTreeInterface as STI
import ConservativeRegridding as CR
import ConservativeRegridding: Trees
import DimensionalData as DD
import SparseArrays

const A5N = DGG.A5.A5Native
const S = A5DGGS()
const SD = GO.UnitSpherical.spherical_distance

# Complete levels, built by expanding the hierarchy one level at a time — never
# by the dense-ordinal arithmetic the ordinal testset checks against them.
const ROOTS = DGG.root_ids(S)
const RES1 = reduce(vcat, [DGG.cell_children(S, 0, root) for root in ROOTS])
const RES2 = reduce(vcat, [DGG.cell_children(S, 1, cell) for cell in RES1])
const RES3 = reduce(vcat, [DGG.cell_children(S, 2, cell) for cell in RES2])

ascending(v) = issorted(v; lt=(<=))

# `treeify` is the package's re-export of `Trees.treeify`; the one-argument form
# picks the manifold up from `best_manifold(grid)`.

"Every leaf index the tree yields, in traversal order."
all_leaf_indices(tree) = STI.depth_first_search(Returns(true), tree)

"Leaf indices whose extent meets `cap`, via the public depth-first search."
query(tree, cap) = sort!(STI.depth_first_search(intersects_cap(cap), tree))

"Upper bound on a query: the leaves whose own `cell_cap` meets `cap`."
cap_hits(level, ids, cap) =
    [i for i in eachindex(ids) if GO.UnitSpherical._intersects(cap, DGG.cell_cap(S, level, ids[i]))]

"Lower bound on a query: the leaves with a boundary vertex strictly inside `cap`."
geometry_hits(level, ids, cap) =
    [i for i in eachindex(ids)
     if any(p -> SD(cap.point, p) < cap.radius, DGG.cell_boundary(S, level, ids[i]))]

ring_points(polygon) = collect(GI.getpoint(GI.getexterior(polygon)))

# ---------------------------------------------------------------------------
@testset "A5 kernel facts and traits" begin
    @test DGG.cell_id_type(S) === UInt64
    # All three traits sit at their kernel defaults — this is the conservative
    # corner of the trait space, and the reason the cursor takes the fallback.
    @test DGG.has_ordinal_ids(S) == false
    @test DGG.has_descendant_ranges(S) == false
    @test DGG.has_exact_subtree_cap(S) == false
    @test_throws NotPortedError DGG.descendant_range(S, 0, ROOTS[1], 1)
    @test_throws NotPortedError radix(S)          # no fixed radix to derive from

    @test root_count(S) == 12
    @test max_level(S) == 30                      # the id encoding's deepest level
    @test A5N.MAX_GRID_RESOLUTION == 29           # ... but no uniform grid at 30

    # Raised above the package default; the sweep at the bottom is the evidence.
    @test DGG.cell_cap_inflation(S) == 1.75
    @test DGG.cell_cap_inflation(S) > DGG.CELL_CAP_INFLATION
    @test DGG.cell_cap_inflation(DGG.H3DGGS()) == DGG.CELL_CAP_INFLATION

    # 12 roots, then 60 quintants, then aperture 4: never `root_count * radix^r`.
    @test DGG.num_cells(S, 0) == 12
    @test DGG.num_cells(S, 1) == 60
    @test all(DGG.num_cells(S, r) == 60 * Int64(4)^(r - 1) for r in 1:20)
    @test all(DGG.num_cells(S, r) == Int64(A5N.num_cells(r)) for r in 0:29)
    @test DGG.num_cells(S, 29) == 60 * Int64(4)^28
    # Res 30 exists only for the 42 quintants that fit the encoding, so there is
    # no full-world grid there and the kernel refuses rather than guessing.
    @test_throws ArgumentError DGG.num_cells(S, 30)
    @test_throws ArgumentError DGG.subtree_leaf_count(S, 0, ROOTS[1], 30)
end

# ---------------------------------------------------------------------------
@testset "A5 hierarchy: 12 -> 5 -> 4" begin
    @test length(ROOTS) == 12
    @test ROOTS == collect(A5N.res0_cells())
    @test ascending(ROOTS)
    @test all(root -> A5N.get_resolution(root) == 0, ROOTS)
    @test DGG.root_ids(S) !== DGG.root_ids(S)     # callers get their own copy

    # Level 0 -> 1: five quintants each. The native enumeration walks segments,
    # whose serialized quintant is a rotation of `0:4`, so the kernel has to
    # sort — that is what these two assertions pin.
    for root in ROOTS
        children = DGG.cell_children(S, 0, root)
        @test length(children) == 5
        @test ascending(children)
        @test all(child -> DGG.cell_parent(S, 1, child, 0) == root, children)
        @test all(child -> A5N.get_resolution(child) == 1, children)
    end
    @test length(RES1) == 60
    @test ascending(RES1)

    # Level >= 1: aperture 4, and the native order is already ascending.
    for cell in RES1
        children = DGG.cell_children(S, 1, cell)
        @test length(children) == 4
        @test ascending(children)
        @test all(child -> DGG.cell_parent(S, 2, child, 1) == cell, children)
    end
    @test length(RES2) == 240
    @test ascending(RES2)
    for cell in RES2[1:17:end]
        children = DGG.cell_children(S, 2, cell)
        @test length(children) == 4
        @test ascending(children)
        @test all(child -> DGG.cell_parent(S, 3, child, 2) == cell, children)
    end
    @test length(RES3) == 960
    @test ascending(RES3)

    # Two- and three-level parent hops, back to the roots.
    @test all(cell -> DGG.cell_parent(S, 3, cell, 0) in ROOTS, RES3)
    @test all(cell -> DGG.cell_parent(S, 3, cell, 1) in RES1, RES3)
    @test all(cell -> DGG.cell_parent(S, 2, cell, 2) == cell, RES2)

    # `cell_descendants` against the level-by-level expansion.
    @test DGG.cell_descendants(S, 0, ROOTS[1], 0) == [ROOTS[1]]
    @test reduce(vcat, [DGG.cell_descendants(S, 0, root, 2) for root in ROOTS]) == RES2
    @test reduce(vcat, [DGG.cell_descendants(S, 1, cell, 3) for cell in RES1]) == RES3
    @test all(root -> ascending(DGG.cell_descendants(S, 0, root, 3)), ROOTS)
    @test_throws ArgumentError DGG.cell_descendants(S, 2, RES2[1], 1)
end

# ---------------------------------------------------------------------------
@testset "A5 counts: subtree_leaf_count vs enumeration" begin
    # The two closed forms, against an actual enumeration.
    for root in ROOTS, delta in 0:4
        @test DGG.subtree_leaf_count(S, 0, root, delta) ==
              length(DGG.cell_descendants(S, 0, root, delta))
        @test DGG.subtree_leaf_count(S, 0, root, delta) ==
              (delta == 0 ? 1 : 5 * Int64(4)^(delta - 1))
    end
    for cell in RES1[1:7:end], delta in 0:5
        @test DGG.subtree_leaf_count(S, 1, cell, 1 + delta) ==
              length(DGG.cell_descendants(S, 1, cell, 1 + delta))
        @test DGG.subtree_leaf_count(S, 1, cell, 1 + delta) == Int64(4)^delta
    end
    for cell in RES2[1:31:end], delta in 0:4
        @test DGG.subtree_leaf_count(S, 2, cell, 2 + delta) ==
              length(DGG.cell_descendants(S, 2, cell, 2 + delta))
        @test DGG.subtree_leaf_count(S, 2, cell, 2 + delta) == Int64(4)^delta
    end

    # A complete level decomposes into the roots' subtrees, and into any
    # intermediate level's — anchored to `num_cells`.
    for leaf in 0:6
        @test sum(DGG.subtree_leaf_count(S, 0, root, leaf) for root in ROOTS) ==
              DGG.num_cells(S, leaf)
    end
    for leaf in 1:8
        @test sum(DGG.subtree_leaf_count(S, 1, cell, leaf) for cell in RES1) ==
              DGG.num_cells(S, leaf)
    end
    @test DGG.subtree_leaf_count(S, 0, ROOTS[1], 29) == DGG.num_cells(S, 29) ÷ 12
    @test_throws ArgumentError DGG.subtree_leaf_count(S, 3, RES3[1], 2)
end

# ---------------------------------------------------------------------------
# NEIGHBOR-VALIDATION, the `max_neighbors` counterpart of the CAP-VALIDATION
# battery at the bottom of this file: the bound is a *measurement*, printed, not
# a number read off the cell shape. "A5 cells are pentagons, so 5" is not an
# argument — the res >= 2 answer is a `Set` union of two independently computed
# families (`_find_quintant_neighbor_s`'s within-quintant lattice step and
# `_get_boundary_neighbors`' quintant-seam, face-seam, apex and corner special
# cases), and a union can only ever come out bigger than either part. The sweep
# below is what says it does not.
# ---------------------------------------------------------------------------
@testset "A5 NEIGHBOR-VALIDATION: the degree bound is measured" begin
    bound = DGG.max_neighbors(S)
    @test bound == 5

    # Deep levels are sampled by ordinal, as the cap battery samples them.
    ordinal_sample(level, n) = (total = DGG.num_cells(S, level);
        [DGG.ordinal_to_cell(S, level, o) for o in 1:max(1, (total - 1) ÷ (n - 1)):total])
    RES4 = reduce(vcat, [DGG.cell_children(S, 3, cell) for cell in RES3])

    println("\n  A5 NEIGHBOR-VALIDATION — edge degree (`length(cell_neighbors(...))`) per level")
    groups = (("res 0 (all 12)", 0, ROOTS), ("res 1 (all 60)", 1, RES1),
        ("res 2 (all 240)", 2, RES2), ("res 3 (all 960)", 3, RES3),
        ("res 4 (all 3840)", 4, RES4),
        ("res 9 (400 sample)", 9, ordinal_sample(9, 400)),
        ("res 15 (400 sample)", 15, ordinal_sample(15, 400)),
        ("res 22 (400 sample)", 22, ordinal_sample(22, 400)),
        ("res 29 (400 sample)", 29, ordinal_sample(29, 400)))
    worst = 0
    for (label, level, cells) in groups
        degrees = [length(DGG.cell_neighbors(S, level, id)) for id in cells]
        @printf("  %-20s (%5d cells) degrees observed: %s\n", label, length(cells),
            join(sort(unique(degrees)), ", "))
        worst = max(worst, maximum(degrees))
        # Not merely bounded by 5 — *constant* per regime, which is the stronger
        # statement and the one that makes the bound safe to extrapolate past
        # the sampled levels: 5 at res 0 (a dodecahedron face has 5 face-adjacent
        # faces), 3 at res 1 (a quintant is a triangle: two sibling quintants and
        # one across a face seam), 5 at res >= 2 (the pentagon lattice is
        # edge-to-edge at the apex and along the seams as much as inside a
        # quintant). No cell anywhere below reports 4, 6 or more.
        @test unique(degrees) == [level == 1 ? 3 : 5]
        @test maximum(degrees) <= bound
    end
    @printf("  worst edge degree observed %d against the wired bound %d\n\n", worst, bound)
    # The bound is attained, so it has no slack to trim and none to spare.
    @test worst == bound

    # The container the bound sizes.
    @test DGG.cell_neighbors(S, 0, ROOTS[1]) isa DGG.SmallVector{5,UInt64}
end

# ---------------------------------------------------------------------------
# The kernel contract itself (`cell_neighbors` in `src/core/kernel.jl`):
# ascending, unique, self-excluded, valid cells of the same level, and
# symmetric. Checked exhaustively over complete levels rather than on samples,
# because `edge_only = true` had no caller in the package before this wiring —
# `_spherical_to_cell` takes the native default — so these testsets are the
# whole of that mode's coverage.
# ---------------------------------------------------------------------------
@testset "A5 neighbors: the kernel contract" begin
    RES4 = reduce(vcat, [DGG.cell_children(S, 3, cell) for cell in RES3])
    RES5 = reduce(vcat, [DGG.cell_children(S, 4, cell) for cell in RES4])

    # Complete levels: with every cell of the level in hand, "b in N(a) implies
    # a in N(b)" over all a *is* symmetry, and set membership is exact validity.
    for (level, cells) in ((0, ROOTS), (1, RES1), (2, RES2), (3, RES3), (4, RES4), (5, RES5))
        table = Dict(id => DGG.cell_neighbors(S, level, id) for id in cells)
        stored = Set(cells)
        @test all(id -> ascending(table[id]), cells)          # ascending *and* unique
        @test all(id -> !(id in table[id]), cells)            # never itself
        @test all(id -> all(n -> A5N.get_resolution(n) == level, table[id]), cells)
        @test all(id -> all(n -> n in stored, table[id]), cells)
        @test all(id -> all(n -> id in table[n], table[id]), cells)
    end

    # Deep levels, where no complete level can be held: the same checks on a
    # patch — sampled cells closed under one neighbor step — so that both
    # directions of symmetry are testable for every pair inside the patch.
    ordinal_sample(level, n) = (total = DGG.num_cells(S, level);
        [DGG.ordinal_to_cell(S, level, o) for o in 1:max(1, (total - 1) ÷ (n - 1)):total])
    for level in (9, 15, 22, 29)
        seeds = ordinal_sample(level, 120)
        patch = Set{UInt64}(seeds)
        for id in seeds, n in DGG.cell_neighbors(S, level, id)
            push!(patch, n)
        end
        cells = collect(patch)
        table = Dict(id => DGG.cell_neighbors(S, level, id) for id in cells)
        @test all(id -> ascending(table[id]) && !(id in table[id]), cells)
        @test all(id -> all(n -> A5N.get_resolution(n) == level, table[id]), cells)
        # Validity without a level to enumerate: the ordinal round trip, which
        # only closes for an id that really is a cell of `level`.
        @test all(id -> all(n -> DGG.ordinal_to_cell(S, level,
                DGG.cell_to_ordinal(S, level, n)) == n, table[id]), cells)
        # Symmetry in both directions for every pair with both ends in the patch.
        @test all(id -> all(n -> !(n in patch) || (id in table[n]), table[id]), cells)
    end

    # `edge_only = true` is not the native default, and the difference is not
    # cosmetic: the default is a *vertex* neighborhood, which is what
    # `_spherical_to_cell` wants (more candidates around a point estimate) and
    # is emphatically not what `stencil` or the halo table mean. The wired set
    # is always a strict subset below res 0.
    @test length(A5N._get_global_cell_neighbors(RES1[1]; edge_only=true)) == 3
    @test length(A5N._get_global_cell_neighbors(RES1[1]; edge_only=false)) == 11
    @test all(RES2) do id
        vertexwise = A5N._get_global_cell_neighbors(id; edge_only=false)
        edgewise = DGG.cell_neighbors(S, 2, id)
        length(vertexwise) in 6:8 && Set(edgewise) ⊊ Set(vertexwise)
    end
    # Res 0 is the one level where the keyword changes nothing: a dodecahedron
    # face's 5 face-adjacent faces are also all the faces it shares a vertex
    # with, and `_get_res0_neighbors` ignores `edge_only` outright.
    @test all(root -> A5N._get_global_cell_neighbors(root; edge_only=false) ==
                      collect(DGG.cell_neighbors(S, 0, root)), ROOTS)

    # Guards, the `cell_to_ordinal` discipline: a full-grid level, and an id
    # that is a cell of exactly that level.
    @test_throws ArgumentError DGG.cell_neighbors(S, 30, ROOTS[1])
    @test_throws ArgumentError DGG.cell_neighbors(S, -1, ROOTS[1])
    @test_throws ArgumentError DGG.cell_neighbors(S, 2, ROOTS[1])   # res-0 id at level 2
    @test_throws ArgumentError DGG.cell_neighbors(S, 0, RES2[5])    # res-2 id at level 0
    @test_throws ArgumentError DGG.cell_neighbors(S, 0, UInt64(0))  # the world cell

    # WHY level 30 is refused rather than passed through. Only 42 of the 60
    # quintants fit the res-30 encoding, so a res-30 neighbor landing outside
    # them takes `serialize`'s `S >> 2, MAX_RESOLUTION - 1` fallback and comes
    # back as a *res-29* id. The native answer there mixes resolutions and is
    # not symmetric — this pins that the refusal is a native limit being
    # declined, not a restriction invented by the wiring.
    res30 = A5N.serialize(A5N.A5Cell(A5N.ORIGINS[1], A5N.ORIGINS[1].first_quintant, UInt64(0), 30))
    @test A5N.get_resolution(res30) == 30
    native30 = A5N._get_global_cell_neighbors(res30; edge_only=true)
    @test any(n -> A5N.get_resolution(n) != 30, native30)
end

# ---------------------------------------------------------------------------
# `subtree_border` — the operation `cell_neighbors` unblocks. A5 keeps the
# generic fallback (no index-only shortcut exists for a non-congruent, non-rep-4
# tiling), so what is checked here is the fallback's *pruning*, which is the
# only thing that can make it disagree with the definition.
# ---------------------------------------------------------------------------
@testset "A5 subtree_border vs brute force" begin
    # The definition, straight: a leaf is on the rim iff some edge neighbor of
    # it has a different ancestor at the root's level.
    function brute_border(level, root, leaf)
        out = UInt64[]
        for descendant in DGG.cell_descendants(S, level, root, leaf)
            any(n -> DGG.cell_parent(S, leaf, n, level) != root,
                DGG.cell_neighbors(S, leaf, descendant)) && push!(out, descendant)
        end
        return sort!(out)
    end

    for (level, roots) in ((0, ROOTS), (1, RES1[1:11:end]),
            (2, RES2[1:53:end]), (3, RES3[1:211:end]))
        for delta in 0:5, root in roots
            leaf = level + delta
            border = DGG.subtree_border(S, level, root, leaf)
            @test border == brute_border(level, root, leaf)
            @test ascending(border)
            @test border ⊆ DGG.cell_descendants(S, level, root, leaf)
        end
    end
    # One deeper probe per regime, where the pruning has had room to compound.
    for (level, root, leaf) in ((0, ROOTS[9], 6), (1, RES1[7], 7), (2, RES2[33], 8))
        @test DGG.subtree_border(S, level, root, leaf) == brute_border(level, root, leaf)
    end

    # Depth 0 is the cell itself — and only for an id that is a cell at `level`,
    # which it learns from `cell_descendants` rather than by spelling `[id]`.
    @test DGG.subtree_border(S, 1, RES1[3], 1) == [RES1[3]]
    @test_throws ArgumentError DGG.subtree_border(S, 2, RES2[1], 1)

    # Every quintant of a root is on the rim (each has a cross-face neighbor in
    # another pentagon), and the rim shrinks as a fraction of the subtree by
    # roughly half per level, as a perimeter should against an area.
    @test DGG.subtree_border(S, 0, ROOTS[1], 1) == DGG.cell_children(S, 0, ROOTS[1])
    fractions = [length(DGG.subtree_border(S, 0, ROOTS[1], leaf)) /
                 DGG.subtree_leaf_count(S, 0, ROOTS[1], leaf) for leaf in 2:6]
    @test fractions == sort(fractions; rev=true)
    @test fractions[end] < 0.2

    # THE PRUNING'S PREMISE, measured. The fallback drops the children of any
    # cell that touches nothing outside the subtree, on the premise that a
    # cell's children touch nothing outside the children of that cell and of its
    # edge neighbors. A5 satisfies that at res 0 and from res 2 down — but NOT
    # at res 1, where it fails twice per quintant.
    function closure_violations(level, cells)
        bad = 0
        for cell in cells
            allowed = Set{UInt64}(DGG.cell_children(S, level, cell))
            for neighbor in DGG.cell_neighbors(S, level, cell)
                union!(allowed, DGG.cell_children(S, level, neighbor))
            end
            for child in DGG.cell_children(S, level, cell),
                m in DGG.cell_neighbors(S, level + 1, child)

                m in allowed || (bad += 1)
            end
        end
        return bad
    end
    @test closure_violations(0, ROOTS) == 0
    @test closure_violations(1, RES1) == 120        # 2 per quintant — the gap
    @test closure_violations(2, RES2) == 0
    @test closure_violations(3, RES3[1:7:end]) == 0

    # What the 120 are: res-2 cells at a pentagon apex whose edge neighbors sit
    # in a pentagon that is only a *vertex* neighbor of the parent quintant —
    # the res-1 regime is triangles with 3 edges while the res-2 lattice under
    # it is pentagons with 5, so the child reaches past its parent's edge ring.
    escapes = UInt64[]
    for child in DGG.cell_children(S, 1, RES1[1]), m in DGG.cell_neighbors(S, 2, child)
        parent = DGG.cell_parent(S, 2, m, 1)
        (parent == RES1[1] || parent in DGG.cell_neighbors(S, 1, RES1[1])) && continue
        push!(escapes, parent)
    end
    @test length(escapes) == 2
    @test all(p -> p in A5N._get_global_cell_neighbors(RES1[1]; edge_only=false), escapes)

    # ... and why the gap never bites: a quintant is never *interior* to a
    # subtree, so the fallback never prunes one. Its cross-face neighbor always
    # lies in another pentagon, so for a res-0 root all five quintants stay in
    # the frontier; a res-1 root is its own frontier and expands to res 2
    # directly; a deeper root never traverses res 1 at all. That is the whole
    # reason the brute-force agreement above holds.
    @test all(RES1) do quintant
        root = DGG.cell_parent(S, 1, quintant, 0)
        any(n -> DGG.cell_parent(S, 1, n, 0) != root, DGG.cell_neighbors(S, 1, quintant))
    end
end

# ---------------------------------------------------------------------------
# The lookup layer `cell_neighbors` unblocks: `neighbor_indices` is the halo
# table `stencil` and `zonal` resolve values through.
# ---------------------------------------------------------------------------
@testset "A5 lookup neighbor operations" begin
    ids = DGG.cell_descendants(S, 0, ROOTS[1], 3)
    lookup = A5Lookup(ids; resolution=3, validate=true)
    halo = DGG.neighbor_indices(lookup)
    @test length(halo) == 80
    @test eltype(halo) === DGG.SmallVector{5,Int}
    @test all(v -> length(v) == 5, halo)             # 5 neighbors, stored or not
    # A stored position points at the cell `cell_neighbors` named, in the same
    # ascending order.
    @test all(eachindex(ids)) do i
        [ids[j] for j in halo[i] if j > 0] ==
        [n for n in DGG.cell_neighbors(S, 3, ids[i]) if !isnothing(findfirst(==(n), ids))]
    end
    # The zeros — neighbors outside the stored coverage — are exactly the
    # subtree's rim, which is `subtree_border` seen from the other side.
    @test sort(ids[[i for i in eachindex(halo) if any(==(0), halo[i])]]) ==
          DGG.subtree_border(S, 0, ROOTS[1], 3)

    # `stencil` over that table: sum of the stored neighbors, so an interior
    # cell sees 5 and a rim cell fewer.
    array = DD.DimArray(ones(80), (A5Cells(lookup),))
    counts = parent(DGG.stencil((_, vals) -> length(vals), array; nbidx=halo))
    @test length(counts) == 80
    @test maximum(counts) == 5
    @test minimum(counts) < 5
    @test count(==(5), counts) == 80 - length(DGG.subtree_border(S, 0, ROOTS[1], 3))
end

# ---------------------------------------------------------------------------
@testset "A5 dense ordinals" begin
    # Ascending id order is `(quintant, S)` lexicographic, so the ordinal is a
    # bijection onto `1:num_cells` and monotone in the id.
    for (level, ids) in ((0, ROOTS), (1, RES1), (2, RES2), (3, RES3))
        @test length(ids) == DGG.num_cells(S, level)
        @test [DGG.cell_to_ordinal(S, level, id) for id in ids] == collect(1:length(ids))
        @test [DGG.ordinal_to_cell(S, level, o) for o in 1:length(ids)] == ids
    end

    # Deep levels, without enumerating them: a sweep of ordinals must come back
    # ascending, at the right resolution, and round-trip.
    for level in (5, 9, 15, 29)
        total = DGG.num_cells(S, level)
        # `total` reaches 4.3e18 at level 29, so the sweep steps rather than
        # scaling — `k * total` would overflow `Int64`.
        step = (total - 1) ÷ 17
        ordinals = sort!(unique!(vcat([1, 2, total - 1, total],
            [1 + k * step for k in 1:16])))
        cells = [DGG.ordinal_to_cell(S, level, o) for o in ordinals]
        @test all(cell -> A5N.get_resolution(cell) == level, cells)
        @test ascending(cells)
        @test [DGG.cell_to_ordinal(S, level, cell) for cell in cells] == ordinals
        # Consecutive ordinals really are consecutive cells of the level.
        @test DGG.cell_parent(S, level, DGG.ordinal_to_cell(S, level, 1), 0) == ROOTS[1]
    end

    # The kernel-uniform error for an ordinal that names no cell of its level.
    @test_throws DGG.OrdinalRangeError DGG.ordinal_to_cell(S, 2, 0)
    @test_throws DGG.OrdinalRangeError DGG.ordinal_to_cell(S, 2, DGG.num_cells(S, 2) + 1)
    err = try
        DGG.ordinal_to_cell(S, 2, DGG.num_cells(S, 2) + 1)
    catch e
        e
    end
    @test err.system === :A5
    @test err.level == 2
    @test err.total == DGG.num_cells(S, 2)
    @test occursin("A5", sprint(showerror, err))
    @test_throws ArgumentError DGG.ordinal_to_cell(S, 30, 1)
    @test_throws ArgumentError DGG.cell_to_ordinal(S, 30, ROOTS[1])
    # An A5 id carries its own resolution, so a mislabelled `level` is an error
    # rather than a silently valid ordinal of some other level's grid — the
    # `get_resolution(id) == level` discipline the range systems have.
    @test_throws ArgumentError DGG.cell_to_ordinal(S, 2, ROOTS[1])
    @test_throws ArgumentError DGG.cell_to_ordinal(S, 0, RES2[5])
end

# ---------------------------------------------------------------------------
@testset "A5 geometry" begin
    # Ring shape: a root is a pentagon, a quintant a triangle, everything below
    # a pentagon again — each subdivided by `segments = :auto`.
    for (level, id, corners) in ((0, ROOTS[1], 5), (1, RES1[1], 3),
            (2, RES2[1], 5), (3, RES3[1], 5), (7, DGG.cell_descendants(S, 3, RES3[1], 7)[1], 5))
        open_ring = DGG.cell_boundary(S, level, id)
        closed = DGG.cell_boundary(S, level, id; closed=true)
        @test open_ring isa Vector{GO.UnitSphericalPoint{Float64}}
        @test length(open_ring) % corners == 0
        @test length(closed) == length(open_ring) + 1
        @test closed[1:length(open_ring)] == open_ring
        @test closed[end] == closed[1]
        @test all(p -> abs(sqrt(sum(abs2, p)) - 1) < 1e-12, open_ring)
        @test DGG.cell_boundary(S, level, id) == open_ring   # `closed` never mutates

        # Native center, not the kernel's boundary mean.
        center = DGG.cell_center(S, level, id)
        lon, lat = A5N.cell_to_lonlat(id)
        @test abs(sqrt(sum(abs2, center)) - 1) < 1e-12
        @test center[3] ≈ sin(deg2rad(lat)) atol = 1e-14
        @test atan(center[2], center[1]) ≈ deg2rad(lon) atol = 1e-12

        cap = DGG.cell_cap(S, level, id)
        @test cap.point == center
        @test cap.radius <= Float64(pi)
        @test all(p -> SD(center, p) <= cap.radius, open_ring)
        @test SD(cap.point, center) < cap.radius            # center inside its cap

        polygon = DGG.cell_polygon_unitsphere(S, level, id)
        @test GI.trait(polygon) isa GI.PolygonTrait
        @test ring_points(polygon) == closed
    end

    # A5's projection is equal-area on the ELLIPSOID, and `cell_to_lonlat` /
    # `cell_boundary` report geodetic coordinates — so a whole level's polygons
    # partition the unit sphere (`sum == 4pi`) but their individual areas carry
    # the authalic -> geodetic latitude conversion, ~1% peak to peak. Undoing
    # that one conversion on the wired ring restores equal area to round-off,
    # which is the strongest available check that the boundary is A5's cell and
    # not a projection artifact — and it only passes if the ring really is in
    # the geodetic frame the center and `lonlat_to_cell` use.
    #
    # From level 2 the native ring subdivides edges that are straight in the
    # face plane but not on the sphere, and summing 240 eighty-vertex spherical
    # polygons drifts a fraction of a per mille; that residual is the native
    # geometry's, not the wiring's (the unsubdivided five-chord ring sums
    # exactly but spreads the cell areas twelve times wider, 1.6% against
    # 0.13%), so level 2 gets a tolerance, not an identity.
    level_area(level, ids) = [abs(GO.area(GO.Spherical(; radius=1.0),
        DGG.cell_polygon_unitsphere(S, level, id))) for id in ids]
    # The same ring pulled back to A5's authalic sphere: geodetic latitude
    # forward through `_authalic_forward`, longitude untouched. Closed the way
    # `cell_polygon_unitsphere` closes it (open ring plus its first point), not
    # with the native `closed_ring = true`, which rotates the ring by a vertex.
    function authalic_area(level, id)
        ring = A5N.cell_boundary(UInt64(id); closed_ring=false)
        points = map(ring) do lonlat
            λ = deg2rad(lonlat[1])
            φ = A5N._authalic_forward(deg2rad(lonlat[2]))
            GO.UnitSphericalPoint(cos(φ) * cos(λ), cos(φ) * sin(λ), sin(φ))
        end
        push!(points, points[1])
        return abs(GO.area(GO.Spherical(; radius=1.0), GI.Polygon([GI.LinearRing(points)])))
    end

    for (level, ids) in ((0, ROOTS), (1, RES1))
        areas = level_area(level, ids)
        @test sum(areas) ≈ 4pi rtol = 1e-7               # measured 2e-8 / 5e-8
        @test maximum(areas) - minimum(areas) < 0.015 * maximum(areas)  # 0.82% / 1.07%
        # ... and exactly equal once the latitude conversion is undone.
        authalic = [authalic_area(level, id) for id in ids]
        @test sum(authalic) ≈ 4pi rtol = 1e-12
        @test maximum(authalic) - minimum(authalic) < 1e-12 * maximum(authalic)
    end
    res2_areas = level_area(2, RES2)
    @test sum(res2_areas) ≈ 4pi rtol = 1e-3
    @test (maximum(res2_areas) - minimum(res2_areas)) < 2e-2 * maximum(res2_areas)
    res2_authalic = [authalic_area(2, id) for id in RES2]
    @test sum(res2_authalic) ≈ 4pi rtol = 1e-3
    @test (maximum(res2_authalic) - minimum(res2_authalic)) < 5e-3 * maximum(res2_authalic)
end

# ---------------------------------------------------------------------------
# FRAME REGRESSION. `A5Native` has two spherical frames: the internal one
# `_to_cartesian` produces (no `LONGITUDE_OFFSET` un-rotation, authalic
# latitude) and the geographic one `_to_lonlat` produces. Wiring `cell_boundary`
# to the former while `cell_center` used the latter put the two 93 degrees apart
# in longitude — leaf caps at res 3 were 1.168 rad instead of 0.144, subtree
# containment was meaningless, and every cross-system regrid was silently
# misaligned. The shape, count and ordering testsets above cannot see any of
# that: each of them looks at the boundary alone, never at the boundary and the
# center together. These three assertions do, directly and cheaply.
# ---------------------------------------------------------------------------
@testset "A5 geometry: boundary and center share one frame" begin
    lonlat(p) = (rad2deg(atan(p[2], p[1])), rad2deg(asin(clamp(p[3], -1.0, 1.0))))
    # A point half way from the center to a vertex, back on the sphere.
    function halfway(center, vertex)
        v = ntuple(i -> 0.5 * (center[i] + vertex[i]), 3)
        norm = sqrt(sum(abs2, v))
        return GO.UnitSphericalPoint(v[1] / norm, v[2] / norm, v[3] / norm)
    end

    sample_at(level, n) = level == 0 ? ROOTS :
                          (total = DGG.num_cells(S, level);
                           [DGG.ordinal_to_cell(S, level, o) for o in 1:max(1, (total - 1) ÷ (n - 1)):total])

    previous = Inf
    for level in 0:6
        cells = sample_at(level, 25)
        radii = [maximum(SD(DGG.cell_center(S, level, id), p)
                         for p in DGG.cell_boundary(S, level, id)) for id in cells]
        worst = maximum(radii)
        # (a) The cell's own radius, against an absolute per-resolution budget.
        # A frame mismatch shows up here as radii of order 1 rad at every level
        # (the measured symptom was 1.67-2.06); the real values are 0.65498 at
        # res 0 and then halve, 0.40945 / 0.20490 / 0.10343 / 0.05200 / 0.02599
        # / 0.01301, so the budget below runs ~8% above the measurement.
        @test worst <= (level == 0 ? 0.70 : 0.45 * 0.5^(level - 1))
        @test worst < previous                      # and it shrinks with resolution
        previous = worst

        for id in cells
            center = DGG.cell_center(S, level, id)
            # (b) center -> lonlat -> cell round-trips through the native point
            # lookup, which is the frame `A5Lookups` addresses.
            @test A5N.lonlat_to_cell(lonlat(center)..., level) == id
            # (c) ... and so does the boundary: half way out to every vertex is
            # still inside the same cell. This is the assertion the offset ring
            # fails outright — its vertices are a different hemisphere's.
            @test all(v -> A5N.lonlat_to_cell(lonlat(halfway(center, v))..., level) == id,
                DGG.cell_boundary(S, level, id))
        end
    end

    # Absolute-scale pin for the query sandwich below: a res-3 leaf cap is a
    # small local cap, not the near-hemispheric 1.168 rad the broken frame gave.
    res3_caps = [DGG.cell_cap(S, 3, id).radius for id in RES3[1:37:end]]
    @test maximum(res3_caps) < 0.5
    @test maximum(res3_caps) < 0.30                 # measured 0.181 at inflation 1.75
    @test minimum(res3_caps) > 0.0
end

# ---------------------------------------------------------------------------
@testset "A5 grid construction" begin
    for level in 0:3
        grid = DGGSGrid(S, level)
        @test grid.system === S
        @test DGG.num_cells(grid.system, grid.level) == DGG.num_cells(S, level)
    end
    # `max_level` is the id encoding's bound, and that is all the constructor
    # validates (it cannot call `num_cells` — unwired systems must still
    # construct; see the note at the top of `src/core/grid_types.jl`). Res 30
    # exists only for the 42 quintants that fit 64 bits, so this grid is
    # constructible and every operation that needs a cell count refuses.
    @test DGGSGrid(S, 30).level == 30
    @test_throws ArgumentError DGG.num_cells(S, 30)
    @test_throws ArgumentError treeify(DGGSGrid(S, 30))
    @test_throws ArgumentError DGGSGrid(S, 31)
    @test_throws ArgumentError DGGSGrid(S, -1)

    ids = DGG.cell_descendants(S, 1, RES1[7], 4)
    lookup = A5Lookup(ids; resolution=4, validate=true)
    grid = DGGSPartialGrid(lookup)
    @test grid.system === S
    @test grid.level == 4
    @test grid.ids === lookup.data                       # by reference, not collected
    @test grid.bucket_size == 0
    @test grid.root_level == -1
    @test eltype(grid.ids) === DGG.cell_id_type(S)
    # `kwargs...` forwarding, all three keywords.
    @test DGGSPartialGrid(lookup; bucket_size=8).bucket_size == 8
    rooted = DGGSPartialGrid(lookup; root_level=1, root_id=RES1[7])
    @test rooted.root_level == 1 && rooted.root_id == RES1[7]
    @test_throws ArgumentError DGGSPartialGrid(lookup; root_level=1, root_id=RES1[8])

    chunk = subtree_grid(S, RES1[7]; root_level=1, leaf_level=4)
    @test chunk.ids == ids
    @test chunk.root_level == 1 && chunk.root_id == RES1[7]
    @test length(chunk.ids) == DGG.subtree_leaf_count(S, 1, RES1[7], 4)

    # The generic constructor still enforces its own invariants.
    @test_throws ArgumentError DGGSPartialGrid(S, 3, collect(Int64, 1:4))
    @test_throws ArgumentError DGGSPartialGrid(S, 3, reverse(RES3[1:8]))
end

# ---------------------------------------------------------------------------
# The `cell_parent` fallback path: `treeify` materializes an index vector per
# node instead of the `searchsorted` bounds the `descendant_range` systems get.
# ---------------------------------------------------------------------------
@testset "A5 cursor traversal (fallback path)" begin
    ids = DGG.cell_descendants(S, 0, ROOTS[3], 4)
    sparse = ids[1:13:end]
    for (level, stored, bucket) in ((4, ids, 0), (4, ids, 7), (4, sparse, 0),
            (3, RES3, 0), (2, RES2, 5))
        tree = treeify(DGGSPartialGrid(S, level, stored; bucket_size=bucket))
        @test tree isa DGG.DGGSCursor
        @test tree.selection isa Vector{Int}          # the fallback path, not bounds
        @test tree.level == -1
        @test Trees.ncells(tree) == length(stored)

        # Full leaf coverage, exactly once each.
        @test sort(all_leaf_indices(tree)) == collect(1:length(stored))
        # Leaf order is the stored order: index i is position i of `ids`.
        @test all(i -> ring_points(Trees.getcell(tree, i)) ==
                       ring_points(DGG.cell_polygon_unitsphere(S, level, stored[i])),
            1:min(length(stored), 40))
        @test_throws BoundsError Trees.getcell(tree, length(stored) + 1)
    end

    # An empty grid is a leaf with no entries, not an error.
    empty_tree = treeify(DGGSPartialGrid(S, 3, UInt64[]))
    @test STI.isleaf(empty_tree)
    @test isempty(STI.child_indices_extents(empty_tree))
    @test Trees.ncells(empty_tree) == 0

    # A subtree-rooted chunk answers in chunk-local indices, and its root is the
    # chunk's own cell rather than the synthetic whole-sphere node.
    chunk = subtree_grid(S, RES1[20]; root_level=1, leaf_level=4)
    chunk_tree = treeify(chunk)
    @test chunk_tree.level == 1 && chunk_tree.id == RES1[20]
    @test Trees.ncells(chunk_tree) == 64
    @test sort(all_leaf_indices(chunk_tree)) == collect(1:64)

    # The public node accessors and the re-exported tree vocabulary say the
    # same thing as the fields and the `Trees.`-qualified calls.
    @test (node_level(chunk_tree), node_id(chunk_tree)) == (1, RES1[20])
    @test ncells(chunk_tree) == Trees.ncells(chunk_tree)
    @test getcell(chunk_tree, 1) == Trees.getcell(chunk_tree, 1)
    for child in STI.getchild(chunk_tree)
        @test node_level(child) == 2
        @test DGG.cell_parent(S, 2, node_id(child), 1) == RES1[20]
    end

    # The dense grid takes the same hierarchy with ordinal leaf indices; without
    # `descendant_range` its internal nodes carry no interval at all.
    dense_tree = treeify(DGGSGrid(S, 2))
    @test Trees.ncells(dense_tree) == 240
    @test sort(all_leaf_indices(dense_tree)) == collect(1:240)
    @test all(i -> ring_points(Trees.getcell(dense_tree, i)) ==
                   ring_points(DGG.cell_polygon_unitsphere(S, 2, RES2[i])), 1:240)
    root_child = first(STI.getchild(dense_tree))
    @test root_child.level == 0 && !STI.isleaf(root_child)
    @test root_child.first_index == 0                # no ordinal interval
    @test Trees.ncells(root_child) == DGG.subtree_leaf_count(S, 0, ROOTS[1], 2)
end

# ---------------------------------------------------------------------------
# Query exactness. The kernel's node extents contain their descendants' *cells*
# but not their descendants' inflated leaf *caps*, so the traversal's hits are
# sandwiched: every leaf with a boundary vertex inside the query cap is found,
# and no leaf whose own `cell_cap` misses the query cap ever is.
# ---------------------------------------------------------------------------
@testset "A5 cap queries vs brute force" begin
    ids = DGG.cell_descendants(S, 0, ROOTS[5], 4)
    sparse = ids[1:11:end]
    caps = [DGG.cell_cap(S, 1, RES1[1]), DGG.cell_cap(S, 2, RES2[100]),
        DGG.cell_cap(S, 4, ids[17]),
        GO.UnitSpherical.SphericalCap(DGG.cell_center(S, 0, ROOTS[5]), 0.35),
        GO.UnitSpherical.SphericalCap(GO.UnitSphericalPoint(0.0, 0.0, 1.0), 0.5)]

    nonempty = 0
    for (level, stored) in ((4, ids), (4, sparse), (3, RES3)), bucket in (0, 6)
        tree = treeify(DGGSPartialGrid(S, level, stored; bucket_size=bucket))
        for cap in caps
            hits = query(tree, cap)
            @test hits ⊆ cap_hits(level, stored, cap)
            @test geometry_hits(level, stored, cap) ⊆ hits
            @test allunique(hits)
            nonempty += !isempty(hits)
        end
    end
    @test nonempty > 0

    # Dense grid: leaf indices are ordinals, same sandwich.
    dense_tree = treeify(DGGSGrid(S, 2))
    for cap in caps
        hits = query(dense_tree, cap)
        @test hits ⊆ cap_hits(2, RES2, cap)
        @test geometry_hits(2, RES2, cap) ⊆ hits
    end

    # An off-grid cap: antipodal to the stored subtree, well clear of it.
    center = DGG.cell_center(S, 0, ROOTS[5])
    away = GO.UnitSpherical.SphericalCap(
        GO.UnitSphericalPoint(-center[1], -center[2], -center[3]), 0.05)
    @test isempty(cap_hits(4, ids, away))                  # the ground truth agrees
    for bucket in (0, 6)
        @test isempty(query(treeify(DGGSPartialGrid(S, 4, ids; bucket_size=bucket)), away))
    end
end

# ---------------------------------------------------------------------------
# `ConservativeRegridding.Regridder` over the generic grids — the consumer the
# whole tree layer exists for, and (in the cross-system testset below) the
# second net under the boundary/center frame.
# ---------------------------------------------------------------------------
@testset "A5 Regridder round trips" begin
    # A partial grid against itself: the intersection matrix must add up to the
    # destination's own area budget. Main's A5 suite asserted exactly this over
    # the five quintants of root 1, and it still holds.
    lookup = A5Lookup(DGG.cell_children(S, 0, ROOTS[1]); resolution=1, validate=true)
    partial = CR.Regridder(DGGSPartialGrid(lookup), DGGSPartialGrid(lookup);
        threaded=false, normalize=false)
    @test size(partial.intersections) == (5, 5)
    # 5 diagonal entries plus the round-off slivers the spherical clipper
    # reports where two quintants of the same pentagon share an edge; they
    # cancel to nothing in the total.
    @test SparseArrays.nnz(partial.intersections) >= 5
    @test isapprox(sum(partial.intersections), sum(partial.dst_areas); rtol=1e-10)

    # ... and a subtree grid of small cells, where the native ring is a plain
    # convex five-chord pentagon (`segments = :auto` is one segment per edge
    # from res 6 down) and the identity is exact rather than clipper-limited.
    chunk = subtree_grid(S, RES3[41]; root_level=3, leaf_level=6)
    @test length(chunk.ids) == 64
    deep = CR.Regridder(chunk, chunk; threaded=false, normalize=false)
    @test size(deep.intersections) == (64, 64)
    @test SparseArrays.nnz(deep.intersections) == 64        # diagonal only
    @test isapprox(sum(deep.intersections), sum(deep.dst_areas); rtol=1e-10)

    # The dense grid, at the coarsest level there is. Each cell still regrids
    # onto itself exactly; the off-diagonal entries are GeometryOps' spherical
    # clipper on A5's eighty-vertex res-0 rings, which are not geodesically
    # convex in either of A5's frames, so only the diagonal is asserted (the
    # deep-cell identity above is the exact statement).
    dense = CR.Regridder(DGGSGrid(S, 0), DGGSGrid(S, 0); threaded=false, normalize=false)
    @test size(dense.intersections) == (12, 12)
    @test SparseArrays.nnz(dense.intersections) >= 12
    @test all(i -> isapprox(dense.intersections[i, i], dense.dst_areas[i]; rtol=1e-6), 1:12)
end

# ---------------------------------------------------------------------------
# CROSS-SYSTEM ALIGNMENT — the second net under the frame fix. A5 and HEALPix
# discretize the same longitude-dependent field, both regrid it onto the same
# 15-degree lon/lat destination, and the two answers must agree with each other
# and with the analytic field. A boundary ring that is rotated away from its own
# center (the bug this suite's frame testset pins directly) moves each source
# cell's polygon away from the value sampled at its center, and every number
# below blows up to the field's own scale.
# ---------------------------------------------------------------------------
@testset "A5 cross-system regridding alignment" begin
    to_sphere = GO.UnitSpherical.UnitSphereFromGeographic()
    dst_lon = collect(range(0, 360; length=25))
    dst_lat = collect(range(-90, 90; length=13))
    destination = [to_sphere((x, y)) for x in dst_lon, y in dst_lat]

    # `1 + x` on the unit sphere is `1 + cos(lat)cos(lon)`: a pure longitude
    # signal, range [0, 2], which is what makes it a longitude-offset detector.
    field(p) = 1.0 + p[1]
    analytic = vec([field(to_sphere(((dst_lon[i] + dst_lon[i + 1]) / 2,
                                     (dst_lat[j] + dst_lat[j + 1]) / 2)))
                    for i in 1:(length(dst_lon) - 1), j in 1:(length(dst_lat) - 1)])

    function regrid_from(system, level)
        regridder = CR.Regridder(destination, DGGSGrid(system, level))
        total = Int(DGG.num_cells(system, level))
        source = [field(DGG.cell_center(system, level, DGG.ordinal_to_cell(system, level, i)))
                  for i in 1:total]
        values = zeros(length(regridder.dst_areas))
        CR.regrid!(values, regridder, source)
        return regridder, source, values
    end

    a5, a5_source, a5_values = regrid_from(S, 3)
    healpix, _, healpix_values = regrid_from(DGG.HEALPixDGGS(), 3)

    # Conservation first: the source tiles the same sphere as the destination
    # and every source column is fully consumed. (A5's densified rings cost the
    # spherical clipper about a percent here; the tolerances are that, measured
    # — 0.95% column, 0.037% total.)
    @test size(a5.intersections) == (length(analytic), Int(DGG.num_cells(S, 3)))
    @test maximum(abs.(vec(sum(a5.intersections; dims=1)) .- a5.src_areas) ./ a5.src_areas) < 0.02
    @test abs(sum(a5.src_areas) - sum(a5.dst_areas)) / sum(a5.dst_areas) < 1e-3

    # Mass: the area-weighted mean survives the regrid.
    source_mean = sum(a5_source .* a5.src_areas) / sum(a5.src_areas)
    dst_mean = sum(a5_values .* a5.dst_areas) / sum(a5.dst_areas)
    @test abs(dst_mean - source_mean) < 5e-3

    # Pattern: A5 lands where its longitudes say it does. Measured 0.058 against
    # the analytic field and 0.054 against HEALPix, on a field of range 1.97; a
    # 93-degree boundary offset puts both at the field's own scale.
    @test maximum(abs.(a5_values .- analytic)) < 0.15
    @test maximum(abs.(a5_values .- healpix_values)) < 0.15
    @test maximum(analytic) - minimum(analytic) > 1.9         # the scale it is measured against
end

# ---------------------------------------------------------------------------
# CAP-VALIDATION (`cell_cap` in `src/core/kernel.jl`). "Union ratio" = max
# distance from the wired cap's center to any delta-level descendant vertex,
# divided by the cell's own max center-to-vertex distance — i.e. the inflation
# factor a subtree actually needs. A5's pentagon lattice only approximately
# nests (see the comment on `cell_cap_inflation` in `src/A5/A5Kernel.jl`), so
# the ratio is large and essentially level-independent: 1.45159 worst measured,
# 1.46872 extrapolated, against the shared 1.2 default. That is why
# `A5Kernel.jl` raises `cell_cap_inflation` to 1.75.
# ---------------------------------------------------------------------------
@testset "A5 CAP-VALIDATION: subtree union ratios" begin
    inflation = DGG.cell_cap_inflation(S)

    function union_ratios(cells, level, deltas)
        ratios = zeros(length(deltas))
        of_cap = 0.0
        for id in cells
            cap = DGG.cell_cap(S, level, id)
            raw = maximum(SD(cap.point, p) for p in DGG.cell_boundary(S, level, id))
            for (k, delta) in enumerate(deltas)
                leaf = level + delta
                worst = 0.0
                for descendant in DGG.cell_descendants(S, level, id, leaf)
                    for p in DGG.cell_boundary(S, leaf, descendant)
                        distance = SD(cap.point, p)
                        distance > worst && (worst = distance)
                    end
                end
                ratios[k] = max(ratios[k], worst / raw)
                of_cap = max(of_cap, worst / cap.radius)
            end
        end
        return ratios, of_cap
    end

    # Deep levels are sampled by ordinal so the sweep never enumerates a globe.
    ordinal_sample(level, n) = (total = DGG.num_cells(S, level);
        [DGG.ordinal_to_cell(S, level, o) for o in 1:max(1, (total - 1) ÷ (n - 1)):total])

    groups = (
        ("res 0 (all 12)", ROOTS, 0, 1:5),
        ("res 1 (all 60)", RES1, 1, 1:5),
        ("res 2 (all 240)", RES2, 2, 1:5),
        ("res 3 (96 sample)", RES3[1:10:end], 3, 1:6),
        ("res 5 (200 sample)", ordinal_sample(5, 200), 5, 1:6),
        ("res 8 (60 sample)", ordinal_sample(8, 60), 8, 1:6),
        ("res 2 (deep probe)", RES2[1:20:end], 2, 1:8),
    )

    println("\n  A5 CAP-VALIDATION — union ratio (max descendant-vertex distance / " *
            "the cell's own radius before the $(inflation) inflation)")
    worst_ratio = 0.0
    worst_of_cap = 0.0
    worst_extrapolated = 0.0
    convergence_bad = 0
    for (label, cells, level, deltas) in groups
        ratios, of_cap = union_ratios(cells, level, deltas)
        increments = diff(ratios)
        @printf("  %-22s (%4d cells) deltas %d:%d\n    ratios     %s\n    increments %s\n",
            label, length(cells), first(deltas), last(deltas),
            join((@sprintf("%9.5f", r) for r in ratios), " "),
            join((@sprintf("%9.6f", d) for d in increments), " "))
        # Convergence envelope. A5 refines by aperture 4 below the quintants, so
        # the overhang shrinks by ~1/2 per delta — but it can approach its limit
        # from either side (res 0 and res 1 converge downwards), so the envelope
        # is on the increments' magnitude, not their sign. The first three deltas
        # are exempt: delta 1 spans the 5-way quintant cut, and the lattice's
        # drift only settles into its geometric decay once the subtree is a few
        # levels deep (several groups sit flat at delta 1 -> 2, then jump).
        for k in 4:length(increments)
            abs(increments[k]) <= max(0.55 * abs(increments[k - 1]), 1e-5) ||
                (convergence_bad += 1)
        end
        # Geometric tail beyond the last measured delta at ratio 0.55:
        # sup <= r_last + |inc_last| * 0.55 / (1 - 0.55).
        worst_extrapolated = max(worst_extrapolated,
            ratios[end] + abs(increments[end]) * 1.23)
        worst_ratio = max(worst_ratio, maximum(ratios))
        worst_of_cap = max(worst_of_cap, of_cap)
    end
    @printf("  worst union ratio %.5f | worst fraction of the wired cap radius %.5f | extrapolated supremum %.5f\n\n",
        worst_ratio, worst_of_cap, worst_extrapolated)

    @test convergence_bad == 0
    # A5 needs far more than the shared 1.2 default — that is the finding, and
    # the reason `cell_cap_inflation` is raised. 1.45159 measured / 1.46872
    # extrapolated on this group set; the gates leave the measurement room to
    # drift a few percent without pretending the budget is looser than it is.
    @test worst_ratio > DGG.CELL_CAP_INFLATION * 1.15
    @test worst_ratio <= 1.50
    # ... and the raised budget keeps the ~15% headroom the contract asks for,
    # both at the measured deltas and past them.
    @test worst_extrapolated <= inflation * 0.85
    @test worst_of_cap < 0.85                 # every descendant is inside the cap
end

end # module A5KernelTestSuite
