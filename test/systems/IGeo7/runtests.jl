# IGeo7 (ISEA7H + Z7) oracle and conformance tests.
#
# The vectors in `vectors/` contain recorded DGGRID output and independently
# constrain the package's geometry and combinatorics.
#
# Wrapped in a module so the system's names cannot collide with another system's
# in the shared test namespace.

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

# Recorded DGGRID oracle vectors.
const VECTORS = joinpath(@__DIR__, "vectors")

const Z7Cell = I.Z7Cell

@testset "relative Z7 cells" begin
    grid = DGG.PartialGrid(S, Z7Cell("023"), 4)
    c = DGG.cellindex(grid, 10)
    ns = DGG.neighbors(DGG.levelgrid(S, 4), c)
    ds = ns .- Ref(c)
    @test all(d -> d isa DGG.RelativeZ7Cell, ds)
    @test all(d -> d.cell == c, ds)
    @test DGG.directioncode.(ds) == 1:6
    @test all(n -> c + (n - c) == n, ns)
    @test all(n -> n - (n - c) == c, ns)
    @test c - c == DGG.RelativeZ7Cell(c, 0)
    @test c + DGG.RelativeZ7Cell(c, 0) == c
    @test_throws I.RelativeZ7Error first(ns) + (last(ns) - c)
    @test_throws I.RelativeZ7Error Z7Cell("023") - Z7Cell("0230")

    complete = DGG.levelgrid(S, 4)
    first_cell = DGG.cellindex(complete, 1)
    last_cell = DGG.cellindex(complete, DGG.ncells(complete))
    @test first_cell + (last_cell - first_cell) == last_cell
    @test_throws I.RelativeZ7Error first_cell + DGG.RelativeZ7Cell(first_cell, -1)
    # An offset that would overflow a widened target position is still just an
    # out-of-range offset: the check brackets the offset, never the sum.
    @test_throws I.RelativeZ7Error last_cell + DGG.RelativeZ7Cell(last_cell, typemax(Int))
    @test_throws I.RelativeZ7Error DGG.directioncode(last_cell - first_cell)
    # An id that names no cell is Z7's own error, not a displacement error
    @test_throws I.InvalidZ7Error Z7Cell(0xffffffffffffffff) +
                                  DGG.RelativeZ7Cell(Z7Cell(0xffffffffffffffff), 0)

    # Every reason reports itself: the messages are built lazily in `showerror`,
    # so nothing else would notice a field that the formatter reads wrongly.
    for reason in (:level_mismatch, :foreign_origin, :out_of_range, :not_a_neighbor,
        :bogus)
        e = I.RelativeZ7Error(reason, first_cell, last_cell, 3)
        @test occursin("Z7Cell", sprint(showerror, e))
    end
    @test occursin("-1:24010", sprint(showerror,
        I.RelativeZ7Error(:out_of_range, DGG.cellindex(complete, 2), first_cell, -5)))

    # Exhaust every face seam, pentagon, and cone-cut case at small levels.
    for l in 0:2
        level_grid = DGG.levelgrid(S, l)
        for p in 1:DGG.ncells(level_grid)
            origin = DGG.cellindex(level_grid, p)
            for (code, target) in enumerate(DGG.neighbors(level_grid, origin))
                displacement = target - origin
                @test origin + displacement == target
                @test DGG.directioncode(displacement) == code
            end
        end
    end
end

# ---------------------------------------------------------------------------
# Oracle-vector parsers
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

# ---------------------------------------------------------------------------
# `Tally` aggregates large oracle sweeps while retaining the first failing row
# for the test report.
# ---------------------------------------------------------------------------
mutable struct Tally
    n::Int
    first::String
end
Tally() = Tally(0, "")

"Record a failure described by `what` (evaluated only when it fires)."
function bad!(t::Tally, what)
    t.n += 1
    isempty(t.first) && (t.first = string(what))
    return nothing
end

"`cond || bad!(t, what)`, with `what` lazy."
check!(t::Tally, cond::Bool, what) = cond || bad!(t, what)

verdict(t::Tally) = (t.n, t.first)
const CLEAN = (0, "")

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
        t = Tally()
        for row in rows
            r = parse(Int, row[1])
            expected = parse(Int64, row[2])
            check!(t, DGG.ncells(DGG.levelgrid(S, r)) == expected,
                "level $r: ncells $(DGG.ncells(DGG.levelgrid(S, r))) != $expected")
        end
        @test verdict(t) == CLEAN
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
        t = Tally()
        for row in rows
            c = Z7Cell(parse(UInt64, row[2]; base=16))
            check!(t, I.is_pentagon(c), "$(I.z7_string(c)): not a pentagon")
            check!(t, parse(Int, row[7]) == 5, "$(I.z7_string(c)): oracle vertex count != 5")
            check!(t, length(DGG.cell_boundary(g0, c)) == 5,
                "$(I.z7_string(c)): boundary has $(length(DGG.cell_boundary(g0, c))) points, not 5")
            check!(t, I.z7_string(c) == String(row[3]),
                "$(I.z7_string(c)): string form != oracle $(String(row[3]))")
        end
        @test verdict(t) == CLEAN
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

        t = Tally()
        for row in rows
            c = Z7Cell(parse(UInt64, row[1]; base=16))
            text = String(row[2])
            parents = String.(split(row[3], ';'))
            expected_children = String.(split(row[4], ';'))
            r = DGG.level(c)

            check!(t, I.z7_string(c) == text, "$text: string form is $(I.z7_string(c))")
            check!(t, r == length(text) - 2, "$text: level $r != $(length(text) - 2)")

            # children: ascending, and exactly the oracle's list
            kids = collect(DGG.children(S, c))
            check!(t, I.z7_string.(kids) == expected_children,
                "$text: children $(I.z7_string.(kids)) != oracle $expected_children")
            check!(t, ascending(kids), "$text: children not strictly ascending")
            check!(t, all(k -> parent(S, k) == c, kids), "$text: a child's parent is not it")
            check!(t, all(k -> DGG.level(k) == r + 1, kids), "$text: a child is at the wrong level")

            # the parent chain, one level coarser per entry, finest first
            check!(t, [I.z7_string(DGG.ancestor(S, c, l)) for l in (r-1):-1:0] == parents,
                "$text: ancestor chain != oracle $parents")
            walker = c
            for expected in parents
                walker = parent(S, walker)
                check!(t, I.z7_string(walker) == expected,
                    "$text: parent walk reached $(I.z7_string(walker)), expected $expected")
            end

            # every ancestor's descendant_range brackets the cell's position
            grid = DGG.levelgrid(S, r)
            pos = DGG.cellposition(grid, c)
            for (k, ptext) in enumerate(parents)
                anc = Z7Cell(I.z7_from_string(ptext))
                rng = DGG.descendant_range(S, anc, r)
                check!(t, pos in rng, "$text: position $pos outside $ptext's range $rng")
                expected_len = I.is_pentagon(anc) ?
                               (5 * 7^(r - (r - k)) + 1) ÷ 6 : 7^(r - (r - k))
                check!(t, length(rng) == expected_len,
                    "$text: $ptext's range has $(length(rng)) cells, expected $expected_len")
            end
        end
        @test verdict(t) == CLEAN
    end

    # =======================================================================
    # 6. pentagon_chains.csv — six children and the deleted digit
    # =======================================================================

    @testset "6. oracle pentagon_chains.csv vs pentagon children" begin
        header, rows = read_csv(joinpath(VECTORS, "pentagon_chains.csv"))
        @test header == ["base", "depth", "pentagon_z7_string", "children_z7_strings",
            "missing_digits"]
        @test length(rows) == 72          # 12 bases x 6 chain depths

        t = Tally()
        for row in rows
            base = parse(Int, row[1])
            depth = parse(Int, row[2])
            c = Z7Cell(I.z7_from_string(String(row[3])))
            expected = String.(split(row[4], ';'))
            missing_digit = parse(Int, row[5])
            tag = "base $base depth $depth ($(I.z7_string(c)))"

            check!(t, I.is_pentagon(c), "$tag: not a pentagon")
            check!(t, DGG.level(c) == depth - 1, "$tag: level $(DGG.level(c)) != $(depth - 1)")
            check!(t, I.z7_base_cell(DGG.rawid(c)) == base, "$tag: base cell mismatch")
            check!(t, I.z7_deleted_digit(base) == missing_digit,
                "$tag: deleted digit $(I.z7_deleted_digit(base)) != $missing_digit")

            # a pentagon has SIX children, not seven, and the missing one is the
            # base's deleted digit
            kids = collect(DGG.children(S, c))
            check!(t, length(kids) == 6, "$tag: $(length(kids)) children, not 6")
            check!(t, I.z7_string.(kids) == expected, "$tag: children != oracle list")
            check!(t, ascending(kids), "$tag: children not strictly ascending")
            check!(t, all(k -> parent(S, k) == c, kids), "$tag: a child's parent is not it")
            # the deleted digit is rejected outright
            digits_taken = [I.z7_string(k)[end] - '0' for k in kids]
            check!(t, !(missing_digit in digits_taken),
                "$tag: deleted digit $missing_digit appears among the children")
            check!(t, sort(digits_taken) == sort([d for d in 0:6 if d != missing_digit]),
                "$tag: child digits $(sort(digits_taken)) are not 0:6 minus $missing_digit")

            # subtree sizes follow the pentagon recurrence, and the deleted digit
            # leaves NO hole in the position range
            check!(t, length(DGG.descendants(S, c, DGG.level(c) + 1)) == 6,
                "$tag: depth-1 subtree is not 6 cells")
            check!(t, length(DGG.descendants(S, c, DGG.level(c) + 2)) == 41,
                "$tag: depth-2 subtree is not 41 cells")
            rng = DGG.descendant_range(S, c, DGG.level(c) + 1)
            g = DGG.levelgrid(S, DGG.level(c) + 1)
            check!(t, extrema(DGG.cellposition(g, k) for k in kids) == (first(rng), last(rng)),
                "$tag: children do not span descendant_range $rng")
            check!(t, length(rng) == 6, "$tag: descendant_range has $(length(rng)) cells, not 6")
        end
        @test verdict(t) == CLEAN

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
        @test DGG.maxneighbors(S) == 6
        @test DGG.maxneighbors(S, Vertex()) == 6
        @test DGG.maxneighbors(S, Edge()) == 6

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

        # Pin the documented start of the rotational order. Winding is invariant
        # under rotation, but each neighbour position must represent a fixed
        # lattice direction for directional stencil weights.
        #
        # Values produced by the implementation and checked against the
        # documented rule (the six Eisenstein unit steps in lattice order from
        # the dev frame's +1 reference), not typed in from an external oracle.
        let g2 = DGG.levelgrid(S, 2)
            hexagon = DGG.cellindex(g2, 100)
            @test I.z7_string(hexagon) == "0234"
            @test !I.is_pentagon(hexagon)
            @test [I.z7_string(x) for x in DGG.neighbors(g2, hexagon, 1)] ==
                  ["0201", "0203", "0236", "0230", "0235", "0212"]

            # A pentagon yields five, not six: two of the unit directions fold
            # onto the same slot at the cone apex. The start is pinned here too.
            pent = Z7Cell(I.z7_from_string("0400"))
            @test I.is_pentagon(pent)
            @test [I.z7_string(x) for x in DGG.neighbors(g2, pent, 1)] ==
                  ["0405", "0404", "0406", "0403", "0401"]

            # ORACLE PIN on the OUTER rings' start, and on the tolerance that
            # decides it. This pentagon's ring 2 has a cell lying EXACTLY on
            # ring 1's spoke — both sit at compass bearing 198.0 — so which end
            # of the ring it lands on is otherwise a last-bit coin flip.
            # `SPOKE_ATOL` says it starts the ring. Hand-checked: ring 1's
            # bearings are 198, 126, 54, 342, 270 and the five boundary corners
            # are 306, 234, 162, 90, 18, so each neighbour bisects one edge and
            # both sequences decrease by a fifth of a turn — counter-clockwise
            # seen from outside. Ring 2 then steps by 36 degrees.
            @test [I.z7_string(x) for x in DGG.ring(g2, pent, 2)] ==
                  ["0453", "0452", "0441", "0443", "0465",
                      "0461", "0436", "0434", "0412", "0416"]
        end

        # k = 0, and the ring/neighbours relation
        g1 = DGG.levelgrid(S, 2)
        c = DGG.cellindex(g1, 40)
        @test isempty(DGG.neighbors(g1, c, 0))
        @test collect(DGG.ring(g1, c, 0)) == [c]
        @test collect(DGG.ring(g1, c, 1)) == collect(DGG.neighbors(g1, c, 1))
        # The sequence contract requires outward ring concatenation; set
        # equality alone would also accept an id-sorted disc.
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
# 9b. The GBT digit kernel against an independent geometric oracle
    #
    # `_cell_neighbors_ccw` implements the GBT arithmetic documented in
    # `src/systems/IGeo7/gbt.jl`. `_cell_neighbors_ccw_geometric` independently
    # decodes, steps, and re-encodes the lattice. They must return identical ids
    # in identical order for complete levels 0:3, seams, and deep samples.
    # =======================================================================

    @testset "9b. GBT kernel vs geometric oracle" begin
        # The port's own exclusion table and the package's independently fitted
        # deleted-digit table are separate evidence; pin them equal rather than
        # defining one as the other.
        @test I.EXCLUDE_NEIGHBOURS == I.Z7_DELETED_DIGIT

        agrees(z) = collect(I._cell_neighbors_ccw(z)) ==
                    collect(I._cell_neighbors_ccw_geometric(z))

        # every cell of every complete level to 3
        for r in 0:3
            g = DGG.levelgrid(S, r)
            @test count(i -> !agrees(DGG.rawid(DGG.cellindex(g, i))), 1:DGG.ncells(g)) == 0
        end

        # deep samples, where the carry ripples further and the frame rotation
        # has more digits to act on
        rng = Random.MersenneTwister(20260815)
        for r in (4, 6, 8, 10, 12, 15, 19)
            nbad = 0
            for _ in 1:400
                z = (UInt64(rand(rng, 0:11)) << 60) | I.Z7_PAD_MASK
                for _ in 1:r
                    cs = collect(I.z7_children(z))
                    z = cs[rand(rng, eachindex(cs))]
                end
                agrees(z) || (nbad += 1)
            end
            @test nbad == 0
        end

        # the twelve pentagon chains, and the two rings around each link: the
        # exclusion-zone rotation only fires near a chain, so this is where it is
        nbad = 0
        for b in 0:11
            z = (UInt64(b) << 60) | I.Z7_PAD_MASK
            for r in 0:8
                agrees(z) || (nbad += 1)
                for n in I._cell_neighbors_ccw_geometric(z)
                    agrees(n) || (nbad += 1)
                    for m in I._cell_neighbors_ccw_geometric(n)
                        agrees(m) || (nbad += 1)
                    end
                end
                r < 8 && (z = I.z7_child(z, 0))
            end
        end
        @test nbad == 0

        # A cell's six raw steps are distinct. Pentagon handling skips the
        # missing direction rather than computing a duplicate, so a
        # short list would mean a real collision, not a deduplication.
        for r in 1:3
            g = DGG.levelgrid(S, r)
            @test all(1:DGG.ncells(g)) do i
                c = DGG.cellindex(g, i)
                ns = I._cell_neighbors_ccw(DGG.rawid(c))
                allunique(ns) && length(ns) == (I.is_pentagon(c) ? 5 : 6)
            end
        end

        # k = 1 stays allocation-free, which is what makes it usable as the
        # primitive every halo and stencil is built from.
        let g = DGG.levelgrid(S, 8), c = DGG.cellindex(g, 12345)
            DGG.neighbors(g, c, 1)
            @test (@allocated DGG.neighbors(g, c, 1)) == 0
        end

        # ---------------------------------------------------------------
        # k > 1: the whole disc, against a reference walk built on the
        # geometric primitive with the pre-port comparison sort. Same set and
        # same order, both connectivities, k = 0:3.
        # ---------------------------------------------------------------
        # An independent right-handed tangent frame, written out here so the
        # outer-ring order is compared against arithmetic this file owns rather
        # than against the package's own winding helper.
        function ref_frame(centre, toward)
            d = (toward[1] - centre[1], toward[2] - centre[2], toward[3] - centre[3])
            r = d[1] * centre[1] + d[2] * centre[2] + d[3] * centre[3]
            t = (d[1] - r * centre[1], d[2] - r * centre[2], d[3] - r * centre[3])
            n = sqrt(t[1]^2 + t[2]^2 + t[3]^2)
            e1 = (t[1] / n, t[2] / n, t[3] / n)
            return e1, (centre[2] * e1[3] - centre[3] * e1[2],
                centre[3] * e1[1] - centre[1] * e1[3],
                centre[1] * e1[2] - centre[2] * e1[1])
        end
        # A cell exactly on the starting spoke begins the ring: the package's
        # own `SPOKE_ATOL` rule, mirrored here because it is order policy, not
        # neighbour arithmetic, and this oracle is about the arithmetic.
        function ref_azimuth(centre, e1, e2, p)
            d = (p[1] - centre[1], p[2] - centre[2], p[3] - centre[3])
            a = atan(d[1] * e2[1] + d[2] * e2[2] + d[3] * e2[3],
                d[1] * e1[1] + d[2] * e1[2] + d[3] * e1[3])
            a < 0 && (a += 2 * Float64(π))
            return a >= 2 * Float64(π) - DGG.Fallbacks.SPOKE_ATOL ? 0.0 : a
        end

        function reference_shells(g, c, steps)
            shells = Vector{Z7Cell}[]
            ring1(x) = [Z7Cell(z) for z in I._cell_neighbors_ccw_geometric(DGG.rawid(x))]
            first_ring = ring1(c)
            isempty(first_ring) && return shells
            push!(shells, first_ring)
            centre = DGG.cell_centroid(g, c)
            e1, e2 = ref_frame(centre, DGG.cell_centroid(g, first(first_ring)))
            seen = Set{Z7Cell}(first_ring)
            push!(seen, c)
            frontier = first_ring
            for _ in 2:steps
                next = Z7Cell[]
                for x in frontier, y in ring1(x)
                    y in seen && continue
                    push!(seen, y)
                    push!(next, y)
                end
                isempty(next) && break
                sort!(next; by=z -> (ref_azimuth(centre, e1, e2,
                        DGG.cell_centroid(g, z)), z))
                push!(shells, next)
                frontier = next
            end
            return shells
        end

        for r in 1:3
            g = DGG.levelgrid(S, r)
            nbad = 0
            for i in 1:DGG.ncells(g)
                c = DGG.cellindex(g, i)
                shells = reference_shells(g, c, 3)
                for k in 0:3, conn in (Vertex(), Edge())
                    want = k == 0 ? Z7Cell[] :
                           isempty(shells) ? Z7Cell[] :
                           reduce(vcat, shells[1:min(k, length(shells))])
                    collect(DGG.neighbors(g, c, k; connectivity=conn)) == want ||
                        (nbad += 1)
                    wantring = k == 0 ? [c] :
                               k <= length(shells) ? shells[k] : Z7Cell[]
                    collect(DGG.ring(g, c, k; connectivity=conn)) == wantring ||
                        (nbad += 1)
                end
            end
            @test nbad == 0
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
        # the border, cross-checked against a definition that uses adjacency
        # instead of digits: a descendant is on the border iff one of its edge
        # neighbours is not in the subtree
        for (lvl, i, depth) in ((0, 1, 3), (1, 5, 3), (2, 100, 2))
            g = DGG.levelgrid(S, lvl)
            c = DGG.cellindex(g, i)
            for d in 1:depth
                leaf = lvl + d
                leafgrid = DGG.levelgrid(S, leaf)
                border = I.subtree_border(S, c, leaf)
                @test ascending(border)
                @test length(border) == I.subtree_border_count(S, c, leaf)
                inside = Set(DGG.descendants(S, c, leaf))
                @test issubset(Set(border), inside)
                brute = [x for x in DGG.descendants(S, c, leaf)
                         if any(nb -> !(nb in inside), DGG.neighbors(leafgrid, x))]
                @test border == brute
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
        @test tree isa DGG.HierarchicalGridCursor

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
            test_grid_interface(DGG.levelgrid(S, l); label="IGeo7 level $l")
        end
        test_hierarchical_system(S)
    end

end # @testset "IGeo7"

end # module IGeo7SystemTests
