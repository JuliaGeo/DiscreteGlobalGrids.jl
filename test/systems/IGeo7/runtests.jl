# ---------------------------------------------------------------------------
# IGeo7 (ISEA7H + Z7) — the sealed-oracle suites plus both conformance suites.
#
# The oracle vectors in `test/IGeo7/vectors/` are recorded DGGRID output and are
# the authority for every geometric and combinatorial claim here: the old
# suites' assertions are ported, with the call shapes adapted from the retired
# `(system, level, id)` triple to typed `Z7Cell`s.
#
# Wrapped in a module so the system's names cannot collide with another system's
# in a shared test namespace (T7 includes all three from `test/runtests.jl`).
# ---------------------------------------------------------------------------

module IGeo7SystemTests

using Test
using Random
using Printf

using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
const I = DGG.IGeo7
const ISEA = DGG.ISEA

using DiscreteGlobalGridsConformanceTesting

import GeometryOps as GO
const US = GO.UnitSpherical

const S = I.IGeo7System()

# The sealed vectors still live beside the retired suite. Both locations are
# accepted so that moving them under `test/systems/IGeo7/` in a later task needs
# no edit here.
const VECTORS = let
    here = joinpath(@__DIR__, "vectors")
    legacy = normpath(joinpath(@__DIR__, "..", "..", "IGeo7", "vectors"))
    isdir(here) ? here : legacy
end

const Z7Cell = I.Z7Cell

# ---------------------------------------------------------------------------
# Helpers (parsers ported from the retired suite)
# ---------------------------------------------------------------------------

"Read a headered CSV into a header vector and a vector of row vectors."
function read_csv(path)
    lines = readlines(path)
    header = split(strip(lines[1]), ',')
    rows = [split(strip(l), ',') for l in lines[2:end] if !isempty(strip(l))]
    return header, rows
end

"rows (z7, lon, lat) of dggrid_true_res{r}_centers.txt"
function load_true_centers(r::Int)
    rows = Tuple{UInt64,Float64,Float64}[]
    for line in eachline(joinpath(VECTORS, "dggrid_true_res$(r)_centers.txt"))
        s = strip(line)
        isempty(s) && continue
        p = split(s, ',')
        push!(rows, (I.z7_from_string(String(p[1])), parse(Float64, p[2]),
            parse(Float64, p[3])))
    end
    return rows
end

"angular distance in degrees between a UnitSphericalPoint and an (lon, lat)"
lldist(p, lon, lat) = ISEA.angdist((p[1], p[2], p[3]), ISEA.lonlat_to_xyz(lon, lat))

"Strictly ascending."
ascending(v) = issorted(v; lt=(<=))

"The level-`r` pentagon of base `b` — the all-zero digit string."
pentagon(b::Integer, r::Integer) =
    Z7Cell(I.z7_from_string(lpad(string(b), 2, '0') * repeat("0", r)))

@testset "IGeo7" begin

    @test isdir(VECTORS)

    # =======================================================================
    # 1. Z7 codec and the typed cell wrapper
    # =======================================================================

    @testset "1. Z7Cell: codec round-trips, level, order" begin
        c = Z7Cell("0941054")
        @test DGG.level(c) == 5
        @test DGG.rawid(c) isa UInt64
        @test I.z7_string(c) == "0941054"
        @test Z7Cell(I.z7_string(c)) === c
        @test Z7Cell(I.z7_from_hex(I.z7_hex(c))) === c
        @test I.z7_hex(c; prefix=true) == "0x" * I.z7_hex(c)
        @test length(I.z7_hex(c)) == 16

        # level is total and read from the in-band digits alone
        for r in 0:I.MAX_RESOLUTION
            z = Z7Cell(I.z7_from_string("00" * repeat("0", r)))
            @test DGG.level(z) == r
            @test I.is_pentagon(z)
        end

        # ordering is ascending raw id, and agrees with == / hash
        a, b = Z7Cell("0941054"), Z7Cell("0941055")
        @test a < b
        @test isless(a, b) == isless(DGG.rawid(a), DGG.rawid(b))
        @test a == Z7Cell("0941054")
        @test hash(a) == hash(Z7Cell("0941054"))
        @test isbits(a)

        # a cell index is an identity: no grid needed to read its level
        @test DGG.level(Z7Cell("00")) == 0
        @test DGG.cellindextype(S) === Z7Cell
        @test DGG.cellindextypes(S) === (Z7Cell,)
    end

    # =======================================================================
    # 2. num_cells.csv — the world count at every level
    # =======================================================================

    @testset "2. oracle num_cells.csv vs ncells" begin
        header, rows = read_csv(joinpath(VECTORS, "num_cells.csv"))
        @test header == ["res", "count"]
        @test length(rows) == 20
        nbad = 0
        for row in rows
            r = parse(Int, row[1])
            expected = parse(Int64, row[2])
            DGG.ncells(DGG.levelgrid(S, r)) == expected || (nbad += 1)
        end
        @test nbad == 0
        @test DGG.ncells(DGG.levelgrid(S, 0)) == 12
        @test DGG.ncells(DGG.levelgrid(S, 19)) == 113988951853731432
        @test all(DGG.ncells(DGG.levelgrid(S, r)) == 10 * Int64(7)^r + 2 for r in 0:19)
    end

    # =======================================================================
    # 3. res0_cells.csv — the twelve roots
    # =======================================================================

    @testset "3. oracle res0_cells.csv vs rootcells" begin
        header, rows = read_csv(joinpath(VECTORS, "res0_cells.csv"))
        @test header == ["order", "z7_hex", "z7_string", "center_lon", "center_lat",
            "area_m2", "nverts"]
        @test length(rows) == 12

        roots = collect(DGG.rootcells(S))
        @test length(roots) == 12
        @test ascending(roots)
        @test all(c -> DGG.level(c) == 0, roots)
        # the oracle's row order is its own, so compare as a sorted id set
        @test sort(parse.(UInt64, getindex.(rows, 2); base=16)) == DGG.rawid.(roots)

        # roots ARE positions 1:12 of the level-0 grid
        g0 = DGG.levelgrid(S, 0)
        @test [DGG.cellindex(g0, i) for i in 1:12] == roots

        # every root is a pentagon with five vertices
        #
        # NOTE: only the `z7_hex` column of this file is trustworthy as grid
        # geometry. Its `center_lon`/`center_lat` columns are in a different
        # icosahedron orientation (they put a vertex at the pole and rings at
        # +/-26.565 deg, not the standard ISEA placement's +/-58.28), so the
        # retired suite used this file for the id set alone and took the res-0
        # centres from `dggrid_res0_centers.csv` below. Comparing against
        # columns 4/5 here would fail by ~95 degrees.
        nbad = 0
        for row in rows
            c = Z7Cell(parse(UInt64, row[2]; base=16))
            I.is_pentagon(c) || (nbad += 1)
            parse(Int, row[7]) == 5 || (nbad += 1)
            length(DGG.cell_boundary(g0, c)) == 5 || (nbad += 1)
            I.z7_string(c) == String(row[3]) || (nbad += 1)
        end
        @test nbad == 0
    end

    @testset "3b. oracle dggrid_res0_centers.csv vs cell_centroid" begin
        g0 = DGG.levelgrid(S, 0)
        worst = 0.0
        n = 0
        for line in eachline(joinpath(VECTORS, "dggrid_res0_centers.csv"))
            isempty(strip(line)) && continue
            p = split(line, ',')
            b = parse(Int, p[1])
            c = Z7Cell(I.z7_from_string(lpad(string(b), 2, '0')))
            worst = max(worst, lldist(DGG.cell_centroid(g0, c),
                parse(Float64, p[2]), parse(Float64, p[3])))
            n += 1
        end
        @test n == 12
        @test worst <= 1e-8
        # ...and a LOWER guard, so nobody "fixes" the residual: it belongs to
        # the dump's own print precision (~5.9e-9 deg), not to this package.
        @test worst > 1e-10
    end

    # =======================================================================
    # 4. The DGGRID centre dumps, levels 1-5 — all 196,080 cells
    #
    # Both directions, exactly as the retired suite ran them: the computed
    # centroid must match the published centre to within the dump's own print
    # noise (1e-8 deg of great-circle separation; the file itself carries
    # ~6e-9), and `cellat` must decode the published centre back to the exact
    # id. No sampling.
    # =======================================================================

    @testset "4. oracle DGGRID centre dumps, levels 1-5" begin
        total = 0
        for r in 1:5
            rows = load_true_centers(r)
            grid = DGG.levelgrid(S, r)
            total += length(rows)
            worst = 0.0
            ndecode = 0
            npos = 0
            for (z, lon, lat) in rows
                c = Z7Cell(z)
                worst = max(worst, lldist(DGG.cell_centroid(grid, c), lon, lat))
                DGG.cellat(grid, lon, lat) == c || (ndecode += 1)
                # the published cell is a cell of this grid, at this level
                DGG.cellposition(grid, c) === nothing && (npos += 1)
            end
            @printf("    level %d: %6d cells   centroid max err %.3e deg   decode mismatches %d\n",
                r, length(rows), worst, ndecode)
            @test length(rows) == DGG.ncells(grid)
            @test worst <= 1e-8
            @test ndecode == 0
            @test npos == 0
        end
        @test total == 196080
    end

    # =======================================================================
    # 5. hierarchy.csv — parent chains, children, and the ranges they imply
    # =======================================================================

    @testset "5. oracle hierarchy.csv vs parent/children/ancestor" begin
        header, rows = read_csv(joinpath(VECTORS, "hierarchy.csv"))
        @test header == ["z7_hex", "z7_string", "parent_chain_strings", "children_strings"]
        @test length(rows) == 150

        nbad = 0
        for row in rows
            c = Z7Cell(parse(UInt64, row[1]; base=16))
            text = String(row[2])
            parents = String.(split(row[3], ';'))
            expected_children = String.(split(row[4], ';'))
            r = DGG.level(c)

            I.z7_string(c) == text || (nbad += 1)
            r == length(text) - 2 || (nbad += 1)

            # children: ascending, and exactly the oracle's list
            kids = collect(DGG.children(S, c))
            I.z7_string.(kids) == expected_children || (nbad += 1)
            ascending(kids) || (nbad += 1)
            all(k -> parent(S, k) == c, kids) || (nbad += 1)
            all(k -> DGG.level(k) == r + 1, kids) || (nbad += 1)

            # the parent chain, one level coarser per entry, finest first
            [I.z7_string(DGG.ancestor(S, c, l)) for l in (r-1):-1:0] == parents ||
                (nbad += 1)
            walker = c
            for expected in parents
                walker = parent(S, walker)
                I.z7_string(walker) == expected || (nbad += 1)
            end

            # every ancestor's descendant_range brackets the cell's position
            grid = DGG.levelgrid(S, r)
            pos = DGG.cellposition(grid, c)
            for (k, ptext) in enumerate(parents)
                anc = Z7Cell(I.z7_from_string(ptext))
                rng = DGG.descendant_range(S, anc, r)
                pos in rng || (nbad += 1)
                length(rng) == (I.is_pentagon(anc) ?
                                (5 * 7^(r - (r - k)) + 1) ÷ 6 : 7^(r - (r - k))) || (nbad += 1)
            end
        end
        @test nbad == 0
    end

    # =======================================================================
    # 6. pentagon_chains.csv — six children and the deleted digit
    # =======================================================================

    @testset "6. oracle pentagon_chains.csv vs pentagon children" begin
        header, rows = read_csv(joinpath(VECTORS, "pentagon_chains.csv"))
        @test header == ["base", "depth", "pentagon_z7_string", "children_z7_strings",
            "missing_digits"]
        @test length(rows) == 72          # 12 bases x 6 chain depths

        nbad = 0
        for row in rows
            base = parse(Int, row[1])
            depth = parse(Int, row[2])
            c = Z7Cell(I.z7_from_string(String(row[3])))
            expected = String.(split(row[4], ';'))
            missing_digit = parse(Int, row[5])

            I.is_pentagon(c) || (nbad += 1)
            DGG.level(c) == depth - 1 || (nbad += 1)
            I.z7_base_cell(DGG.rawid(c)) == base || (nbad += 1)
            I.z7_deleted_digit(base) == missing_digit || (nbad += 1)

            # a pentagon has SIX children, not seven, and the missing one is the
            # base's deleted digit
            kids = collect(DGG.children(S, c))
            length(kids) == 6 || (nbad += 1)
            I.z7_string.(kids) == expected || (nbad += 1)
            ascending(kids) || (nbad += 1)
            all(k -> parent(S, k) == c, kids) || (nbad += 1)
            # the deleted digit is rejected outright
            digits_taken = [I.z7_string(k)[end] - '0' for k in kids]
            missing_digit in digits_taken && (nbad += 1)
            sort(digits_taken) == sort([d for d in 0:6 if d != missing_digit]) ||
                (nbad += 1)

            # subtree sizes follow the pentagon recurrence, and the deleted digit
            # leaves NO hole in the position range
            length(DGG.descendants(S, c, DGG.level(c) + 1)) == 6 || (nbad += 1)
            length(DGG.descendants(S, c, DGG.level(c) + 2)) == 41 || (nbad += 1)
            rng = DGG.descendant_range(S, c, DGG.level(c) + 1)
            g = DGG.levelgrid(S, DGG.level(c) + 1)
            extrema(DGG.cellposition(g, k) for k in kids) == (first(rng), last(rng)) ||
                (nbad += 1)
            length(rng) == 6 || (nbad += 1)
        end
        @test nbad == 0

        # the deleted-digit table itself: 2 for bases 0-5, 5 for bases 6-11
        @test I.Z7_DELETED_DIGIT == (2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5)
    end

    # =======================================================================
    # 7. Curve order: positions, ids and subtrees agree
    # =======================================================================

    @testset "7. curve order and position/id consistency" begin
        for r in 0:4
            g = DGG.levelgrid(S, r)
            n = DGG.ncells(g)
            cells = [DGG.cellindex(g, i) for i in 1:n]
            @test ascending(cells)                       # position order IS id order
            @test all(i -> DGG.cellposition(g, cells[i]) == i, 1:n)
            @test allunique(cells)
            @test all(c -> DGG.level(c) == r, cells)
        end

        # a cell at the wrong level is simply not in the grid
        g2 = DGG.levelgrid(S, 2)
        @test DGG.cellposition(g2, Z7Cell("00")) === nothing
        @test DGG.cellposition(g2, Z7Cell("00000")) === nothing
        @test DGG.cellposition(g2, Z7Cell("000")) === nothing
        @test DGG.cellposition(g2, Z7Cell("0000")) == 1

        # bounds
        @test_throws BoundsError DGG.cellindex(g2, 0)
        @test_throws BoundsError DGG.cellindex(g2, DGG.ncells(g2) + 1)

        # descendant_range is two-sided and tight, and siblings partition
        for (lvl, leaf) in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 4))
            parents = [DGG.cellindex(DGG.levelgrid(S, lvl), i)
                       for i in 1:DGG.ncells(DGG.levelgrid(S, lvl))]
            leafgrid = DGG.levelgrid(S, leaf)
            previous_hi = 0
            total = 0
            for p in parents
                desc = DGG.descendants(S, p, leaf)
                rng = DGG.descendant_range(S, p, leaf)
                positions = [DGG.cellposition(leafgrid, d) for d in desc]
                @test rng isa UnitRange{Int}
                @test sort(positions) == collect(rng)      # two-sided, tight
                @test first(rng) > previous_hi             # ordered, disjoint
                previous_hi = last(rng)
                total += length(desc)
            end
            @test total == DGG.ncells(leafgrid)            # the subtrees tile the level
            @test previous_hi == DGG.ncells(leafgrid)
        end

        # self-range is the cell's own single position (the cursor depends on it)
        g3 = DGG.levelgrid(S, 3)
        for i in (1, 17, 500, DGG.ncells(g3))
            c = DGG.cellindex(g3, i)
            @test DGG.descendant_range(S, c, 3) == i:i
        end

        # sibling ranges concatenate, in `children` order, to the parent's
        c = DGG.cellindex(DGG.levelgrid(S, 2), 200)
        kids = DGG.children(S, c)
        @test reduce(vcat, collect.(DGG.descendant_range(S, k, 3) for k in kids)) ==
              collect(DGG.descendant_range(S, c, 3))
    end

    # =======================================================================
    # 8. Error shapes
    # =======================================================================

    @testset "8. error shapes" begin
        @test_throws ArgumentError parent(S, Z7Cell("00"))
        @test_throws ArgumentError DGG.children(S, DGG.cellindex(DGG.levelgrid(S, 19), 1))
        @test_throws ArgumentError DGG.levelgrid(S, -1)
        @test_throws ArgumentError DGG.levelgrid(S, 20)
        @test_throws ArgumentError DGG.ancestor(S, Z7Cell("000"), 2)
        @test_throws ArgumentError DGG.descendants(S, Z7Cell("000"), 0)
        @test_throws ArgumentError DGG.descendant_range(S, Z7Cell("000"), 0)
        @test_throws ArgumentError DGG.neighbors(DGG.levelgrid(S, 1), Z7Cell("000"), -1)
        # a cell from the wrong level handed to a grid operation
        @test_throws ArgumentError DGG.neighbors(DGG.levelgrid(S, 1), Z7Cell("0000"))
        # structurally invalid ids are Z7's own error
        @test_throws I.InvalidZ7Error I.z7_from_string("002")     # deleted digit
        @test_throws I.InvalidZ7Error I.z7_from_string("12")      # base out of range
        @test !I.is_valid_cell(0xffffffffffffffff)
    end

    # =======================================================================
    # 9. Neighbours: counts, determinism, symmetry, winding, pentagon seams
    # =======================================================================

    @testset "9. neighbours" begin
        @test DGG.max_neighbors(S) == 6
        @test DGG.max_neighbors(S, Vertex()) == 6
        @test DGG.max_neighbors(S, Edge()) == 6

        for r in 1:3
            g = DGG.levelgrid(S, r)
            n = DGG.ncells(g)
            npent = 0
            nbad = 0
            for i in 1:n
                c = DGG.cellindex(g, i)
                ns = DGG.neighbors(g, c)
                # hexagons have six neighbours, pentagons five
                expected = I.is_pentagon(c) ? 5 : 6
                length(ns) == expected || (nbad += 1)
                I.is_pentagon(c) && (npent += 1)
                allunique(ns) || (nbad += 1)
                (c in ns) && (nbad += 1)
                all(nb -> DGG.level(nb) == r, ns) || (nbad += 1)
                all(nb -> DGG.cellposition(g, nb) !== nothing, ns) || (nbad += 1)
                # determinism
                collect(DGG.neighbors(g, c)) == collect(ns) || (nbad += 1)
                # symmetry — checked for EVERY cell, so a pentagon seam that
                # forgets a neighbour cannot hide behind sampling
                all(nb -> c in DGG.neighbors(g, nb), ns) || (nbad += 1)
                # Edge() and Vertex() coincide on a hex/pentagon grid
                collect(DGG.neighbors(g, c; connectivity=Edge())) == collect(ns) ||
                    (nbad += 1)
            end
            @test nbad == 0
            @test npent == 12          # exactly twelve pentagons per level
        end

        # the documented order is counter-clockwise seen from outside: the
        # signed volume of consecutive neighbour offsets about the cell centre
        # is positive all the way round
        g = DGG.levelgrid(S, 3)
        nccw = 0
        ntested = 0
        for i in 1:200:DGG.ncells(g)
            c = DGG.cellindex(g, i)
            p = DGG.cell_centroid(g, c)
            ns = collect(DGG.neighbors(g, c))
            cen = [DGG.cell_centroid(g, nb) for nb in ns]
            ok = true
            for j in eachindex(cen)
                k = j == length(cen) ? 1 : j + 1
                a = (cen[j][1] - p[1], cen[j][2] - p[2], cen[j][3] - p[3])
                b = (cen[k][1] - p[1], cen[k][2] - p[2], cen[k][3] - p[3])
                cross = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3],
                    a[1] * b[2] - a[2] * b[1])
                (cross[1] * p[1] + cross[2] * p[2] + cross[3] * p[3]) > 0 || (ok = false)
            end
            ntested += 1
            ok && (nccw += 1)
        end
        @test ntested > 0
        @test nccw == ntested

        # k = 0, and the ring/neighbours relation
        g1 = DGG.levelgrid(S, 2)
        c = DGG.cellindex(g1, 40)
        @test isempty(DGG.neighbors(g1, c, 0))
        @test collect(DGG.ring(g1, c, 0)) == [c]
        @test collect(DGG.ring(g1, c, 1)) == collect(DGG.neighbors(g1, c, 1))
        # The rotational contract, on the sequences and not merely on the sets:
        # the disc IS the rings concatenated outward, so the ring is the disc's
        # tail block element for element. Set equality alone passes happily for
        # a disc that sorts by id, which is what this used to do.
        for k in 1:3
            union_rings = reduce(vcat, [collect(DGG.ring(g1, c, j)) for j in 1:k])
            disc = collect(DGG.neighbors(g1, c, k))
            @test disc == union_rings
            shell = collect(DGG.ring(g1, c, k))
            @test disc[(end-length(shell)+1):end] == shell
            @test allunique(union_rings)
        end

        # pentagon seams specifically: every pentagon at several levels
        for r in 1:4
            g = DGG.levelgrid(S, r)
            for b in 0:11
                c = pentagon(b, r)
                ns = DGG.neighbors(g, c)
                @test length(ns) == 5
                @test all(nb -> c in DGG.neighbors(g, nb), ns)
                @test !I.is_pentagon(first(ns))     # a pentagon's neighbours are hexes
            end
        end
    end

    # =======================================================================
    # 10. Geometry: rings, winding, centroids, areas
    # =======================================================================

    @testset "10. geometry" begin
        for r in 0:3
            g = DGG.levelgrid(S, r)
            n = DGG.ncells(g)
            step = max(1, n ÷ 200)
            nbad = 0
            for i in 1:step:n
                c = DGG.cellindex(g, i)
                ring = DGG.cell_boundary(g, c)
                length(ring) == (I.is_pentagon(c) ? 5 : 6) || (nbad += 1)
                # unit norm
                all(p -> abs(sqrt(p[1]^2 + p[2]^2 + p[3]^2) - 1) < 1e-12, ring) ||
                    (nbad += 1)
                # implicitly closed: the first vertex is NOT repeated
                (first(ring) == last(ring)) && (nbad += 1)
                # counter-clockwise seen from outside => positive spherical area
                GO.area(GO.Spherical(radius=1.0), DGG.cell_polygon(g, c)) > 0 ||
                    (nbad += 1)
                # the centroid is strictly inside, and is where cellat sends it
                DGG.cellat(g, DGG.cell_centroid(g, c)) == c || (nbad += 1)
            end
            @test nbad == 0
        end

        # the published corner rings tile the sphere exactly: adjacent cells
        # share corners, so their great-circle edges coincide and the areas sum
        # to 4pi
        for r in 1:3
            g = DGG.levelgrid(S, r)
            total = sum(DGG.cell_area(g, DGG.cellindex(g, i)) for i in 1:DGG.ncells(g))
            @test total ≈ 4pi rtol = 1e-12
        end

        # the ideal equal-area closed form also sums to 4pi, exactly, and is a
        # DIFFERENT quantity from the ring's geodesic area (see the docstring of
        # `equal_area_steradians`)
        for r in 0:3
            g = DGG.levelgrid(S, r)
            total = sum(I.equal_area_steradians(DGG.cellindex(g, i))
                        for i in 1:DGG.ncells(g))
            @test total ≈ 4pi rtol = 1e-12
        end
        g1 = DGG.levelgrid(S, 1)
        hex = DGG.cellindex(g1, 5)
        @test !I.is_pentagon(hex)
        @test I.equal_area_steradians(hex) ≈ 4pi / 70
        @test I.equal_area_steradians(pentagon(0, 1)) ≈ (5 / 6) * 4pi / 70
        # the two disagree by ~1.6% at level 1 — recorded so a future change
        # that silently conflates them fails here
        @test !isapprox(DGG.cell_area(g1, hex), I.equal_area_steradians(hex); rtol=1e-3)
    end

    # =======================================================================
    # 11. The covering law, exhaustively at coarse levels
    #
    # The conformance suite samples; this checks EVERY descendant of every
    # level-0 and level-1 cell to depth 3 against every ancestor extent, which
    # is what pins `cap_inflation = 1.2` for this system.
    # =======================================================================

    @testset "11. covering law (exhaustive, coarse levels)" begin
        worst_overshoot = -Inf
        worst_ratio = 0.0
        for (lvl, depth) in ((0, 3), (1, 2))
            g = DGG.levelgrid(S, lvl)
            for i in 1:DGG.ncells(g)
                c = DGG.cellindex(g, i)
                cap = DGG.node_extent(S, c)
                @test cap.radius <= pi / 2          # convex, so vertices suffice
                raw = maximum(US.spherical_distance(cap.point, p)
                              for p in DGG.cell_boundary(g, c))
                for d in 1:depth
                    leaf = DGG.levelgrid(S, lvl + d)
                    for dc in DGG.descendants(S, c, lvl + d)
                        for p in DGG.cell_boundary(leaf, dc)
                            over = US.spherical_distance(cap.point, p) - cap.radius
                            worst_overshoot = max(worst_overshoot, over)
                            worst_ratio = max(worst_ratio, US.spherical_distance(cap.point, p) / raw)
                        end
                    end
                end
            end
        end
        @printf("    covering law: worst overshoot %.3e rad, worst union ratio %.5f (budget %.2f)\n",
            worst_overshoot, worst_ratio, DGG.cap_inflation(S))
        @test worst_overshoot <= 0.0                     # nothing escapes
        @test worst_ratio <= 1.10                        # the measured bound
        @test worst_ratio < DGG.cap_inflation(S)         # inside the wired budget
    end

    # =======================================================================
    # 12. Subtree border automaton
    # =======================================================================

    @testset "12. subtree border" begin
        # the rim, cross-checked against a definition that uses adjacency
        # instead of digits: a descendant is on the border iff one of its edge
        # neighbours is not in the subtree
        for (lvl, i, depth) in ((0, 1, 3), (1, 5, 3), (2, 100, 2))
            g = DGG.levelgrid(S, lvl)
            c = DGG.cellindex(g, i)
            for d in 1:depth
                leaf = lvl + d
                leafgrid = DGG.levelgrid(S, leaf)
                rim = I.subtree_border(S, c, leaf)
                @test ascending(rim)
                @test length(rim) == I.subtree_border_count(S, c, leaf)
                inside = Set(DGG.descendants(S, c, leaf))
                @test issubset(Set(rim), inside)
                brute = [x for x in DGG.descendants(S, c, leaf)
                         if any(nb -> !(nb in inside), DGG.neighbors(leafgrid, x))]
                @test rim == brute
            end
        end
        # depth 0 is the cell itself ("015" is at level 1)
        @test DGG.level(Z7Cell("015")) == 1
        @test I.subtree_border(S, Z7Cell("015"), 1) == [Z7Cell("015")]
        # the closed forms: 3^(d+1)-3 for a hexagon, 5*(3^d-1)/2 for a pentagon,
        # with d the depth BELOW the cell (level 1 -> level 6 is d = 5)
        @test I.subtree_border_count(S, Z7Cell("015"), 6) == 3 * 3^5 - 3
        @test I.subtree_border_count(S, pentagon(0, 3), 6) == (5 * (3^3 - 1)) ÷ 2
        # a target coarser than the cell is an ArgumentError
        @test_throws ArgumentError I.subtree_border(S, Z7Cell("015"), 0)
        @test_throws ArgumentError I.subtree_border_count(S, Z7Cell("015"), 0)
    end

    # =======================================================================
    # 13. The generic substrate accepts this system
    #
    # Not covered by the conformance suites (which are deliberately
    # cursor-free), but it is what the whole hierarchy wiring is FOR: the grid
    # must treeify to the package cursor and answer a pruned query exactly.
    # =======================================================================

    @testset "13. tree and query integration" begin
        g = DGG.levelgrid(S, 3)
        tree = DGG.treeify(g)
        @test tree isa HierarchicalGridCursor

        cap = US.SphericalCap(GO.UnitSphericalPoint(0.0, 0.0, 1.0), 0.15)
        hits = DGG.query(g, Intersects(cap))
        @test eltype(hits) === Z7Cell
        @test issorted(hits)
        @test !isempty(hits)
        @test allunique(hits)
        # pruning may only ever over-select relative to a vertex test, never
        # drop a cell: every cell with a vertex in the cap must be in the answer
        brute = [DGG.cellindex(g, i) for i in 1:DGG.ncells(g)
                 if any(p -> US.spherical_distance(cap.point, p) <= cap.radius,
                     DGG.cell_boundary(g, DGG.cellindex(g, i)))]
        @test !isempty(brute)
        @test issubset(Set(brute), Set(hits))

        # the ConservativeRegridding.Trees surface works off dense positions.
        # (`isa Any` was here and asserted nothing — every value is an `Any`.
        # The trait check below is the real assertion.)
        @test GO.GI.geomtrait(DGG.getcell(g, 1)) isa GO.GI.PolygonTrait
        @test DGG.ncells(g) == 3432
    end

    # =======================================================================
    # 14. Conformance suites — both, with default keywords
    # =======================================================================

    @testset "14. conformance" begin
        for l in (0, 1, 3)
            test_grid_interface(DGG.levelgrid(S, l); label="IGeo7Grid(level $l)")
        end
        test_hierarchical_system(S)
    end

end # @testset "IGeo7"

end # module IGeo7SystemTests
