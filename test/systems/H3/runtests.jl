# ---------------------------------------------------------------------------
# T5 — the H3 system suite.
#
# Two kinds of test, deliberately kept apart:
#
#   * ORACLE tests, which check the wiring against libh3 itself. libh3 is the
#     definition of H3, so anywhere this package could have drifted — the
#     ordinal numbering, the hierarchy, the boundary, the adjacency — is
#     checked against a direct call rather than against a value someone typed
#     in. The sealed constants that remain (the pentagon base cells, a couple
#     of well-known indices) are asserted against libh3 too.
#   * CONFORMANCE tests, the two property suites from
#     `DiscreteGlobalGridsConformanceTesting`, which check the interface
#     contracts every system owes.
# ---------------------------------------------------------------------------

module H3TestSuite

using Test
using Random
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGridsConformanceTesting

const H3 = DGG.H3
const H3N = H3.H3Native
const S = H3.H3System()

# A deterministic sample of cells at `res`: every pentagon, plus an evenly
# spaced sweep of the level's positions. Deterministic by construction rather
# than by seed, so a failure names the same cell on every run.
function sample_cells(res::Int; stride::Int=0)
    grid = DGG.levelgrid(S, res)
    n = DGG.ncells(grid)
    step = stride == 0 ? max(1, n ÷ 40) : stride
    cells = [DGG.cellindex(grid, i) for i in 1:step:n]
    append!(cells, H3.H3Cell.(H3N.get_pentagons(res)))
    return unique!(sort!(cells))
end

# The lon/lat probe grid for `cellat`, avoiding the poles' coordinate
# degeneracy but keeping high latitudes.
const PROBE_LONLAT = [(lon, lat) for lon in -180.0:24.0:170.0, lat in -78.0:12.0:78.0]

@testset "H3" begin

    # =======================================================================
    @testset "native layer (ported oracles)" begin
        @test length(H3N.res0_cells()) == 122
        @test H3N.num_cells(0) == 122
        @test H3N.num_cells(1) == 842
        @test H3N.num_cells(2) == 5882

        # The canonical H3 documentation example.
        cell = H3N.lonlat_to_cell(-122.41795063018799, 37.775938728915946, 9)
        @test cell == parse(UInt64, "8928308280fffff"; base=16)
        @test H3N.get_resolution(cell) == 9

        center = H3N.cell_center(cell)
        @test isapprox(center[1], -122.41845932318309; atol=1e-12)
        @test isapprox(center[2], 37.776702349435695; atol=1e-12)

        # Child position round-trips over a full sibling set.
        parent_cell = H3N.cell_to_parent(cell, 8)
        kids = H3N.cell_to_children(parent_cell, 9)
        for (position, child) in enumerate(kids)
            @test H3N.cell_to_child_pos(child, 8) == position - 1
            @test H3N.child_pos_to_cell(position - 1, parent_cell, 9) == child
        end
        @test_throws ArgumentError H3N.child_pos_to_cell(-1, parent_cell, 9)

        @test length(H3N.get_pentagons(2)) == 12
        @test H3N.get_resolution("0x8928308280fffff") == 9

        # The two wrappers this port added.
        hex = H3N.lonlat_to_cell(10.0, 45.0, 5)
        @test length(H3N.grid_ring_unsafe(hex, 2)) == 12
        @test H3N.grid_ring_unsafe(hex, 0) == [hex]
        # A pentagon defeats the hollow-ring walk, by design.
        @test H3N.grid_ring_unsafe(first(H3N.get_pentagons(5)), 1) === nothing
    end

    # =======================================================================
    @testset "H3Cell encoding" begin
        for res in 0:15, c in sample_cells(res; stride=max(1, DGG.ncells(DGG.levelgrid(S, res)) ÷ 3))
            # `level` reads bits 52-55; libh3 reads the same field.
            @test DGG.level(c) == H3N.get_resolution(DGG.rawid(c))
            @test DGG.level(c) == res
            @test DGG.rawid(c) isa UInt64
            @test H3.H3Cell(DGG.rawid(c)) == c
        end
        a = H3.H3Cell(0x08001fffffffffff)
        b = H3.H3Cell(0x08002fffffffffff)
        @test a < b
        @test isless(a, b)
        @test a == H3.H3Cell(0x08001fffffffffff)
        @test hash(a) == hash(H3.H3Cell(0x08001fffffffffff))
        # Resolution outranks the base cell: the order is resolution-major.
        @test a < H3.H3Cell(0x0810bfffffffffff)
    end

    # =======================================================================
    @testset "sealed tables vs libh3" begin
        # The 122 res-0 indices, built by bit formula rather than by ccall.
        @test collect(H3._H3_ROOT_IDS) == H3N.res0_cells()
        @test [DGG.rawid(c) for c in DGG.rootcells(S)] == H3N.res0_cells()

        # The twelve pentagon base cells.
        @test sort([H3N.get_base_cell(p) for p in H3N.get_pentagons(0)]) ==
              sort(collect(H3.PENTAGON_BASE_CELLS))

        # The prefix sums, computed analytically, against libh3's own counts.
        for res in 0:15
            counts = [H3N.cell_to_children_size(root, res) for root in H3N.res0_cells()]
            @test H3._H3_ROOT_ENDS[res+1] == cumsum(counts)
            @test H3._H3_ROOT_ENDS[res+1][end] == H3N.num_cells(res)
            @test DGG.ncells(DGG.levelgrid(S, res)) == H3N.num_cells(res)
        end
    end

    # =======================================================================
    @testset "dense order: cellindex/cellposition vs libh3" begin
        for res in 0:6
            grid = DGG.levelgrid(S, res)
            n = DGG.ncells(grid)
            step = max(1, n ÷ 200)
            previous = nothing
            for i in 1:step:n
                c = DGG.cellindex(grid, i)
                @test DGG.cellposition(grid, c) == i
                @test H3N.is_valid_cell(DGG.rawid(c))
                @test H3N.get_resolution(DGG.rawid(c)) == res
                # Position order is id order: the base interface requires it.
                previous === nothing || @test previous < c
                previous = c
                # The ordinal is the base-cell prefix plus the child position,
                # which is what makes it hole-free across pentagon gaps.
                b = H3N.get_base_cell(DGG.rawid(c))
                prefix = b == 0 ? 0 : H3._H3_ROOT_ENDS[res+1][b]
                @test i == prefix + H3N.cell_to_child_pos(DGG.rawid(c), 0) + 1
            end
            @test_throws BoundsError DGG.cellindex(grid, 0)
            @test_throws BoundsError DGG.cellindex(grid, n + 1)
        end

        # Exhaustive at res 0 and 1: every cell, every position.
        for res in 0:1
            grid = DGG.levelgrid(S, res)
            cells = [DGG.cellindex(grid, i) for i in 1:DGG.ncells(grid)]
            @test issorted(cells)
            @test allunique(cells)
            @test length(cells) == H3N.num_cells(res)
            @test Set(DGG.rawid.(cells)) ==
                  Set(reduce(vcat, [H3N.cell_to_children(r, res) for r in H3N.res0_cells()]))
        end

        # A cell from another level is simply not in the grid.
        g2 = DGG.levelgrid(S, 2)
        @test DGG.cellposition(g2, DGG.cellindex(DGG.levelgrid(S, 3), 1)) === nothing
        # Neither is a malformed index at the right level. These are the four
        # ways an H3 index can look plausible and name no cell; libh3's own
        # `cellToChildren` would happily enumerate a subtree of any of them,
        # which is why anything that enumerates checks validity first.
        root = DGG.rawid(DGG.cellindex(DGG.levelgrid(S, 3), 1))
        pentagon0 = H3N.get_pentagons(0)[1]
        res1_pent = (pentagon0 & ~(UInt64(0xf) << 52)) | (UInt64(1) << 52)
        malformed = [
            "padding slot cleared" => root & ~(UInt64(7) << 30),
            "digit 7 in an active slot" => root | (UInt64(7) << 36),
            "base cell 125" => (root & ~(UInt64(0x7f) << 45)) | (UInt64(125) << 45),
            "deleted K-axis child of a pentagon" =>
                (res1_pent & ~(UInt64(7) << 42)) | (UInt64(1) << 42),
        ]
        for (name, id) in malformed
            c = H3.H3Cell(id)
            @test !H3N.is_valid_cell(id)
            @test !isvalid(c)
            @test DGG.cellposition(DGG.levelgrid(S, DGG.level(c)), c) === nothing
            @test_throws ArgumentError H3.subtree_border(S, c, DGG.level(c) + 1)
        end
    end

    # =======================================================================
    @testset "hierarchy vs cellToParent/cellToChildren" begin
        for res in 0:6
            for c in sample_cells(res)
                id = DGG.rawid(c)
                kids = DGG.children(S, c)
                @test DGG.rawid.(collect(kids)) == H3N.cell_to_children(id, res + 1)
                @test length(kids) == (H3N.is_pentagon(id) ? 6 : 7)
                @test issorted(kids)
                @test allunique(kids)
                for k in kids
                    @test DGG.level(k) == res + 1
                    @test parent(S, k) == c
                    @test DGG.rawid(parent(S, k)) == H3N.cell_to_parent(DGG.rawid(k), res)
                end
                if res > 0
                    p = parent(S, c)
                    @test DGG.rawid(p) == H3N.cell_to_parent(id, res - 1)
                    @test c in DGG.children(S, p)
                    # `ancestor` agrees with libh3 at every coarser level.
                    for l in 0:(res-1)
                        @test DGG.rawid(DGG.ancestor(S, c, l)) == H3N.cell_to_parent(id, l)
                    end
                end
                @test DGG.ancestor(S, c, res) == c
                # `descendants` is `cellToChildren` across a level gap.
                for d in 1:2
                    res + d > 15 && continue
                    @test DGG.rawid.(DGG.descendants(S, c, res + d)) ==
                          H3N.cell_to_children(id, res + d)
                end
            end
        end
        # Error contract at the two ends of the hierarchy.
        @test_throws ArgumentError parent(S, first(DGG.rootcells(S)))
        @test_throws ArgumentError DGG.children(S, DGG.cellindex(DGG.levelgrid(S, 15), 1))
        @test_throws ArgumentError DGG.levelgrid(S, -1)
        @test_throws ArgumentError DGG.levelgrid(S, 16)
    end

    # =======================================================================
    @testset "descendant_range" begin
        for res in 0:4
            for c in sample_cells(res; stride=max(1, DGG.ncells(DGG.levelgrid(S, res)) ÷ 12))
                id = DGG.rawid(c)
                @test DGG.descendant_range(S, c, res) ==
                      DGG.cellposition(DGG.levelgrid(S, res), c):DGG.cellposition(DGG.levelgrid(S, res), c)
                for target in (res+1):min(res + 3, 15)
                    r = DGG.descendant_range(S, c, target)
                    @test r isa UnitRange{Int}
                    # Size is libh3's own closed-form count, pentagons included.
                    @test length(r) == H3N.cell_to_children_size(id, target)
                    # Two-sided: the range is exactly the descendants' positions.
                    tgrid = DGG.levelgrid(S, target)
                    actual = sort([DGG.cellposition(tgrid, H3.H3Cell(d))
                                   for d in H3N.cell_to_children(id, target)])
                    @test collect(r) == actual
                end
                # Sibling ranges partition the parent's range, in order.
                if res < 15
                    kids = collect(DGG.children(S, c))
                    ranges = [DGG.descendant_range(S, k, res + 1) for k in kids]
                    @test reduce(vcat, collect.(ranges)) ==
                          collect(DGG.descendant_range(S, c, res + 1))
                end
                res > 0 && @test_throws ArgumentError DGG.descendant_range(S, c, res - 1)
            end
        end
        # Deep ranges: no enumeration, and still exactly the libh3 count.
        deep = DGG.cellindex(DGG.levelgrid(S, 3), 7)
        @test length(DGG.descendant_range(S, deep, 15)) ==
              H3N.cell_to_children_size(DGG.rawid(deep), 15)
    end

    # =======================================================================
    @testset "boundary and centroid vs cellToBoundary" begin
        worst_area = 0.0
        for res in 0:5
            grid = DGG.levelgrid(S, res)
            for c in sample_cells(res)
                id = DGG.rawid(c)
                ring = DGG.cell_boundary(grid, c)
                raw = H3N.cell_boundary_cartesian(id; closed_ring=false)
                # The ring is libh3's, vertex for vertex, unchanged.
                @test length(ring) == length(raw)
                @test 5 <= length(ring) <= 10
                for (p, q) in zip(ring, raw)
                    @test p[1] ≈ q[1] atol = 1e-15
                    @test p[2] ≈ q[2] atol = 1e-15
                    @test p[3] ≈ q[3] atol = 1e-15
                end
                # Unit norm, implicitly closed, counter-clockwise from outside.
                for p in ring
                    @test sum(abs2, p) ≈ 1.0 atol = 1e-12
                end
                @test ring[1] != ring[end]
                @test DGG.cell_area(grid, c) > 0

                # The exact ring reproduces libh3's own area, which is the
                # sharpest available statement that no vertex was moved.
                err = abs(DGG.cell_area(grid, c) - H3N.cell_area(id))
                worst_area = max(worst_area, err)

                # The centroid is libh3's cell centre and is strictly interior.
                centroid = DGG.cell_centroid(grid, c)
                @test sum(abs2, centroid) ≈ 1.0 atol = 1e-12
                lon, lat = H3N.cell_center(id)
                @test atand(centroid[2], centroid[1]) ≈ lon atol = 1e-9
                @test asind(clamp(centroid[3], -1, 1)) ≈ lat atol = 1e-9
            end
        end
        @test worst_area < 1e-12
    end

    # =======================================================================
    @testset "cellat vs latLngToCell" begin
        for res in (0, 1, 4, 7, 11)
            grid = DGG.levelgrid(S, res)
            for (lon, lat) in PROBE_LONLAT
                expected = H3.H3Cell(H3N.lonlat_to_cell(lon, lat, res))
                @test DGG.cellat(grid, lon, lat) == expected
                # The unit-sphere primitive agrees with the lon/lat wrapper.
                p = DGG.Fallbacks.unit_point(lon, lat)
                @test DGG.cellat(grid, p) == expected
            end
        end
        # Round trip: the cell at a cell's own centroid is that cell.
        for res in 0:6
            grid = DGG.levelgrid(S, res)
            for c in sample_cells(res)
                @test DGG.cellat(grid, DGG.cell_centroid(grid, c)) == c
            end
        end
    end

    # =======================================================================
    @testset "neighbors vs gridDisk" begin
        for res in 0:2
            grid = DGG.levelgrid(S, res)
            n = DGG.ncells(grid)
            cells = [DGG.cellindex(grid, i) for i in 1:n]
            nbmap = Dict(c => DGG.neighbors(grid, c) for c in cells)
            for c in cells
                id = DGG.rawid(c)
                nbs = nbmap[c]
                @test issorted(nbs)
                @test allunique(nbs)
                @test !(c in nbs)
                @test length(nbs) == (H3N.is_pentagon(id) ? 5 : 6)
                @test eltype(collect(nbs)) === H3.H3Cell
                # The oracle: gridDisk(k=1) minus the origin and the padding.
                @test Set(DGG.rawid.(collect(nbs))) ==
                      Set(x for x in H3N.grid_disk(id, 1) if x != 0 && x != id)
                # Symmetry, both directions.
                for nb in nbs
                    @test c in nbmap[nb]
                end
            end
            # Euler: 12 pentagons with 5 neighbours, hexagons with 6, and the
            # sum is twice the tiling's edge count, E = 3F - 6.
            @test sum(length(nbmap[c]) for c in cells) == 2 * (3 * H3N.num_cells(res) - 6)
        end

        # k = 0 is empty; Vertex() and Edge() coincide on hexagons.
        grid = DGG.levelgrid(S, 5)
        for c in sample_cells(5)
            @test isempty(DGG.neighbors(grid, c, 0))
            @test collect(DGG.neighbors(grid, c, 1; connectivity=DGG.Vertex())) ==
                  collect(DGG.neighbors(grid, c, 1; connectivity=DGG.Edge()))
            # Larger discs are still gridDisk, sorted.
            for k in 2:3
                @test Set(DGG.rawid.(DGG.neighbors(grid, c, k))) ==
                      Set(x for x in H3N.grid_disk(DGG.rawid(c), k) if x != 0 && x != DGG.rawid(c))
            end
        end
    end

    # =======================================================================
    @testset "ring: unsafe walk and pentagon fallback" begin
        grid = DGG.levelgrid(S, 5)
        for c in sample_cells(5)
            @test DGG.ring(grid, c, 0) == [c]
            shells = Set{H3.H3Cell}()
            for k in 1:3
                shell = DGG.ring(grid, c, k)
                @test issorted(shell)
                @test allunique(shell)
                @test isempty(intersect(Set(shell), shells))
                union!(shells, shell)
                # The disc is the union of its shells.
                @test Set(DGG.neighbors(grid, c, k)) == shells
                # And every shell cell is at exactly that libh3 distance.
                cells, dists = H3N.grid_disk_distances(DGG.rawid(c), k)
                @test Set(DGG.rawid.(shell)) ==
                      Set(id for (id, d) in zip(cells, dists) if id != 0 && Int(d) == k)
            end
        end
        # A pentagon is exactly where the fast walk refuses and the fallback
        # has to carry the answer.
        pent = H3.H3Cell(first(H3N.get_pentagons(5)))
        @test H3N.grid_ring_unsafe(DGG.rawid(pent), 1) === nothing
        @test length(DGG.ring(grid, pent, 1)) == 5
        @test Set(DGG.ring(grid, pent, 1)) == Set(DGG.neighbors(grid, pent, 1))
    end

    # =======================================================================
    @testset "subtree border automaton" begin
        # Brute force: enumerate the subtree and keep the cells with a
        # neighbour outside it. This is the definition the automaton replaces.
        function brute_border(c::H3.H3Cell, target::Int)
            id = DGG.rawid(c)
            lvl = DGG.level(c)
            grid = DGG.levelgrid(S, target)
            return H3.H3Cell[d for d in DGG.descendants(S, c, target)
                             if any(nb -> H3N.cell_to_parent(DGG.rawid(nb), lvl) != id,
                                    DGG.neighbors(grid, d))]
        end

        roots = [
            ("hexagon, res 1 (odd)", H3.H3Cell(H3N.lonlat_to_cell(10.0, 45.0, 1))),
            ("hexagon, res 2 (even)", H3.H3Cell(H3N.lonlat_to_cell(10.0, 45.0, 2))),
            ("hexagon, res 3 (odd)", H3.H3Cell(H3N.lonlat_to_cell(10.0, 45.0, 3))),
            ("hexagon, res 4 (even)", H3.H3Cell(H3N.lonlat_to_cell(77.0, 28.0, 4))),
            ("hexagon, res 5, southern", H3.H3Cell(H3N.lonlat_to_cell(-120.0, -30.0, 5))),
            ("hexagon, res 0 base cell", H3.H3Cell(H3N.lonlat_to_cell(-58.0, -15.0, 0))),
            ("pentagon, res 0", H3.H3Cell(H3N.get_pentagons(0)[1])),
            ("pentagon, res 1", H3.H3Cell(H3N.get_pentagons(1)[1])),
            ("pentagon, res 2", H3.H3Cell(H3N.get_pentagons(2)[1])),
            ("pentagon, res 3", H3.H3Cell(H3N.get_pentagons(3)[1])),
        ]

        for (name, root) in roots
            lvl = DGG.level(root)
            @test H3.subtree_border(S, root, lvl) == [root]
            for depth in 1:3
                lvl + depth > 15 && continue
                border = H3.subtree_border(S, root, lvl + depth)
                @test issorted(border)
                @test allunique(border)
                @test all(c -> DGG.level(c) == lvl + depth, border)
                @test all(c -> H3N.is_valid_cell(DGG.rawid(c)), border)
                # The closed-form rim size.
                expected = H3N.is_pentagon(DGG.rawid(root)) ?
                           (5 * (3^depth - 1)) ÷ 2 : 3^(depth + 1) - 3
                @test length(border) == expected
                # And the definition itself.
                @test border == sort(brute_border(root, lvl + depth))
            end
        end

        # All twelve res-0 pentagons, two levels down.
        for p in H3N.get_pentagons(0)
            root = H3.H3Cell(p)
            for depth in 1:2
                border = H3.subtree_border(S, root, depth)
                @test length(border) == (5 * (3^depth - 1)) ÷ 2
                @test border == sort(brute_border(root, depth))
            end
        end

        # A deep rim is O(3^d) where the subtree is O(7^d): 1_594_320 cells
        # against 13_841_287_201, which is why the automaton exists.
        deep = roots[3][2]
        @test length(H3.subtree_border(S, deep, 15)) == 3^13 - 3
        @test H3N.cell_to_children_size(DGG.rawid(deep), 15) == 7^12

        @test_throws ArgumentError H3.subtree_border(S, roots[3][2], 2)
        @test_throws ArgumentError H3.subtree_border(S, roots[1][2], 16)
    end

    # =======================================================================
    @testset "node_extent covers the subtree (measured bound)" begin
        # `cap_inflation` is a claim about H3's aperture-7 overhang. Re-measure
        # it rather than trusting the number: descend along the outermost
        # branch and track how far a descendant's vertices get from the
        # ancestor's cap centre, relative to its radius.
        function overhang(c::H3.H3Cell, depth::Int; beam::Int=60)
            cap = DGG.Fallbacks.cell_cap(DGG.levelgrid(S, DGG.level(c)), c)
            radius = cap.radius
            frontier = [c]
            worst = 1.0
            for _ in 1:depth
                DGG.level(first(frontier)) >= 15 && break
                nxt = H3.H3Cell[]
                for f in frontier
                    append!(nxt, DGG.children(S, f))
                end
                scored = map(nxt) do x
                    g = DGG.levelgrid(S, DGG.level(x))
                    (maximum(p -> DGG.GO.UnitSpherical.spherical_distance(cap.point, p),
                             DGG.cell_boundary(g, x)), x)
                end
                worst = max(worst, maximum(first, scored) / radius)
                sort!(scored; by=first, rev=true)
                frontier = [x for (_, x) in scored[1:min(beam, length(scored))]]
            end
            return worst
        end

        measured = 0.0
        for res in (0, 2), c in sample_cells(res; stride=max(1, DGG.ncells(DGG.levelgrid(S, res)) ÷ 6))
            measured = max(measured, overhang(c, 6))
        end
        # The measured overhang, and the headroom the declared factor keeps.
        @test measured < 1.10
        @test measured < DGG.cap_inflation(S)
        @test DGG.cap_inflation(S) == 1.2

        # The covering law itself, spot-checked deep: a cell's node extent
        # contains every descendant vertex, all the way to max_level.
        for res in (0, 1, 3)
            for c in sample_cells(res; stride=max(1, DGG.ncells(DGG.levelgrid(S, res)) ÷ 4))
                cap = DGG.node_extent(S, c)
                @test cap.radius <= pi / 2      # convex, so vertices imply arcs
                current = c
                while DGG.level(current) < 15
                    kids = DGG.children(S, current)
                    current = kids[1 + (DGG.level(current) % length(kids))]
                    g = DGG.levelgrid(S, DGG.level(current))
                    for p in DGG.cell_boundary(g, current)
                        @test DGG.Fallbacks.cap_contains(cap, p)
                    end
                end
            end
        end
    end

    # =======================================================================
    @testset "traits" begin
        @test DGG.cellindextype(S) === H3.H3Cell
        @test DGG.levels(S) === 0:15
        @test DGG.max_level(S) == 15
        @test DGG.has_sorted_subtrees(S)
        @test DGG.max_neighbors(S) == 6
        @test DGG.max_neighbors(S, DGG.Vertex()) == 6
        @test DGG.max_neighbors(S, DGG.Edge()) == 6
        @test DGG.cellindextypes(S) === (H3.H3Cell,)
        for l in 0:15
            g = DGG.levelgrid(S, l)
            @test DGG.system(g) === S
            @test DGG.level(g) == l
        end
    end

    # =======================================================================
    @testset "conformance" begin
        for l in (0, 1, 3)
            test_grid_interface(DGG.levelgrid(S, l); label="H3Grid(res $l)")
        end
        test_hierarchical_system(S)
    end
end

end # module H3TestSuite
