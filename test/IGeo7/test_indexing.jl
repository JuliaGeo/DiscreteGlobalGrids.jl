# test/test_indexing.jl — dense full-world indexing, hierarchy wrappers and
# introspection (design.md Section 9 task 7, Section 4.3).
#
# Authorities:
#   - spec/design.md Section 3 (ordering by encoded id) and Section 4.3
#     (`num_cells`, subtree counts, rank walk)
#   - spec/interface-contract.md ("Native API": hierarchy, dense indexing,
#     introspection, error types)
#   - test/IGeo7/vectors/num_cells.csv, test/IGeo7/vectors/res0_cells.csv
#   - ../IGeo7/test/runtests.jl (untainted) for the exact return types and
#     behavioral shape the wiring layers consume
#
# The canonical full-world order is **ascending encoded Z7 id**; because the
# base cell occupies the top nibble and the padding sentinel 7 sorts after
# every active digit, that is (base, digit-string) lexicographic order. Every
# expectation below is built independently of `src/IGeo7/grid.jl`'s rank
# arithmetic: ids are enumerated by digit-lexicographic DFS over `z7_child`,
# and the resulting order is separately asserted to be ascending.

using Test
using Printf
using Random

const _IX = IGeo7
const _IX_SMALL = Helpers
const _IX_VECTORS = joinpath(@__DIR__, "vectors")
const _IX_RNG = MersenneTwister(20260804)

"""Measure allocations of `f(args...)` after a warm-up call, behind a function
barrier so the measurement is not polluted by boxing of untyped globals."""
function _ix_alloc(f::F, args...) where {F}
    f(args...)
    return @allocated f(args...)
end

"Deleted pentagon digit of `base` — the contract's hemisphere rule, restated."
_ix_deleted(base::Int) = base < 6 ? 2 : 5

"""Descendant count of a *pentagon* prefix at depth `d`: `(5·7^d + 1) / 6`
(6 at d = 1, 41 at d = 2 **[contract]**). A hexagon prefix has `7^d`."""
_ix_pent_count(d::Integer) = (5 * Int64(7)^d + 1) ÷ 6

"res-0 id of base cell `b`"
_ix_root(b::Integer) = (UInt64(b) << 60) | UInt64(0x0fffffffffffffff)

"all-zero-digit (pentagon-chain) id of base `b` at resolution `res`"
_ix_pent(b::Integer, res::Integer) = z7_from_string(lpad(string(b), 2, '0') * repeat("0", res))

"""Every valid id at resolution `res`, enumerated by digit-lexicographic DFS
(base ascending, then digit ascending, deleted pentagon digits skipped).
Independent of the implementation's rank arithmetic."""
function _ix_all_ids(res::Int)
    out = UInt64[]
    sizehint!(out, 10 * 7^res + 2)
    for base in 0:11
        _ix_dfs!(out, _ix_root(base), res, true, _ix_deleted(base))
    end
    return out
end

function _ix_dfs!(out::Vector{UInt64}, z::UInt64, res::Int, pentagon::Bool, deleted::Int)
    if _IX.z7_resolution(z) == res
        push!(out, z)
        return out
    end
    for d in 0:6
        pentagon && d == deleted && continue
        _ix_dfs!(out, _IX.z7_child(z, d), res, pentagon & (d == 0), deleted)
    end
    return out
end

"""Independent (deliberately unoptimised) dense rank of `id`: at every level,
sum the subtree size of each sibling digit that sorts before the taken digit,
looping over the siblings explicitly instead of the implementation's closed
form. Cross-checks the `nhex` shortcut in `src/grid.jl`."""
function _ix_rank(id::UInt64)
    res = get_resolution(id)
    base = z7_base_cell(id)
    deleted = _ix_deleted(base)
    rank = Int64(base) * _ix_pent_count(res)
    pentagon = true
    for k in 1:res
        d = z7_digit(id, k)
        depth = res - k
        for e in 0:(d-1)
            (pentagon && e == deleted) && continue
            rank += (pentagon && e == 0) ? _ix_pent_count(depth) : Int64(7)^depth
        end
        pentagon &= (d == 0)
    end
    return Int(rank) + 1
end

"Read a headered CSV into a header vector and a vector of row vectors."
function _ix_read_csv(path)
    lines = readlines(path)
    header = split(strip(lines[1]), ',')
    rows = [split(strip(l), ',') for l in lines[2:end] if !isempty(strip(l))]
    return header, rows
end

@testset verbose = true "grid.jl — dense indexing, hierarchy and introspection" begin

    # -----------------------------------------------------------------------
    @testset "1. num_cells" begin
        header, rows = _ix_read_csv(joinpath(_IX_VECTORS, "num_cells.csv"))
        @test header == ["res", "count"]
        @test length(rows) == 20
        nbad = 0
        for row in rows
            r = parse(Int, row[1])
            expected = parse(Int64, row[2])
            (num_cells(r) === expected) || (nbad += 1)
        end
        @test nbad == 0                                     # oracle table
        @test num_cells(0) == 12
        @test num_cells(1) == 72
        @test num_cells(2) == 492
        @test num_cells(19) == 113988951853731432
        @test all(num_cells(r) == 10 * Int64(7)^r + 2 for r in 0:19)
        @test all(num_cells(r) == 12 * _ix_pent_count(r) for r in 0:19)
        @test num_cells(3) isa Int64
        @test_throws InvalidZ7Error num_cells(-1)
        @test_throws InvalidZ7Error num_cells(20)
        @test_throws InvalidZ7Error num_cells(MAX_RESOLUTION + 1)
    end

    # -----------------------------------------------------------------------
    @testset "2. res0_cells" begin
        roots = res0_cells()
        @test roots isa _IX_SMALL.SmallList{12,UInt64}
        @test length(roots) == 12
        @test issorted(roots)
        @test allunique(roots)
        @test collect(roots) == [_ix_root(b) for b in 0:11]
        @test all(z7_base_cell(roots[i]) == i - 1 for i in 1:12)
        @test all(get_resolution(z) == 0 for z in roots)
        @test all(is_pentagon, roots)
        @test all(is_valid_cell, roots)
        @test collect(roots) == [z7_from_string(lpad(string(b), 2, '0')) for b in 0:11]
        # the oracle's res-0 set, compared as a set of hex ids (its own row
        # order is the sealed module's internal id order, not ours)
        _, rows = _ix_read_csv(joinpath(_IX_VECTORS, "res0_cells.csv"))
        @test sort(parse.(UInt64, getindex.(rows, 2); base=16)) == collect(roots)
        # index 1..12 at res 0 is exactly the ascending root list
        @test index_to_cell.(1:12, 0) == collect(roots)
        @test cell_to_index.(collect(roots)) == 1:12
    end

    # -----------------------------------------------------------------------
    @testset "3. cell_to_children: counts, types, ordering" begin
        roots = res0_cells()
        root = roots[1]

        # types exactly as the untainted wiring consumes them [contract]
        @test cell_to_children(root, 0) isa _IX_SMALL.SmallList{1,UInt64}
        @test cell_to_children(root, 1) isa _IX_SMALL.SmallList{7,UInt64}
        @test cell_to_children(root, 2) isa Vector{UInt64}
        @test cell_to_children(root) isa _IX_SMALL.SmallList{7,UInt64}

        @test collect(cell_to_children(root, 0)) == [root]
        @test length(cell_to_children(root, 1)) == 6          # pentagon
        @test length(cell_to_children(root, 2)) == 41         # [contract]
        @test length(cell_to_children(root)) == 6
        @test collect(cell_to_children(root)) == collect(cell_to_children(root, 1))
        @test z7_to_string.(collect(cell_to_children(root, 1))) ==
              ["000", "001", "003", "004", "005", "006"]

        # southern base: deleted digit 5
        south = roots[7]
        @test z7_base_cell(south) == 6
        @test z7_to_string.(collect(cell_to_children(south, 1))) ==
              ["060", "061", "062", "063", "064", "066"]

        # a hexagon prefix has the full 7 / 49 descendants
        hexagon = z7_from_string("001")
        @test !is_pentagon(hexagon)
        @test length(cell_to_children(hexagon, 2)) == 7
        @test length(cell_to_children(hexagon, 3)) == 49
        @test length(cell_to_children(hexagon)) == 7
        # a deeper pentagon prefix keeps the pentagon count
        @test length(cell_to_children(z7_from_string("0000"), 6)) == _ix_pent_count(4)

        # every base, every depth 0:4: count == pentagon subtree formula, list
        # ascending, all descendants at the target resolution
        nbadcount = 0
        nbadorder = 0
        nbadres = 0
        nbaddesc = 0
        for b in 0:11, d in 0:4
            r = _ix_root(b)
            kids = cell_to_children(r, d)
            length(kids) == _ix_pent_count(d) || (nbadcount += 1)
            issorted(kids) || (nbadorder += 1)
            all(get_resolution(k) == d for k in kids) || (nbadres += 1)
            all(z7_is_descendant(k, r) for k in kids) || (nbaddesc += 1)
        end
        @test nbadcount == 0
        @test nbadorder == 0
        @test nbadres == 0
        @test nbaddesc == 0

        # hexagon subtrees: 7^d for depth 1:4 over a sample of hexagon prefixes
        nbadhex = 0
        for s in ("001", "0130", "064", "1102", "1162"), d in 1:4
            z = z7_from_string(s)
            length(cell_to_children(z, get_resolution(z) + d)) == 7^d || (nbadhex += 1)
        end
        @test nbadhex == 0

        # child ∈ parent's descendant list, at both one and two levels
        child = first(cell_to_children(root, 2))
        @test get_resolution(child) == 2
        @test cell_to_parent(child, 0) == root
        @test child in cell_to_children(cell_to_parent(child, 1), 2)
        @test child in cell_to_children(cell_to_parent(child, 1))
        nbadround = 0
        for id in cell_to_children(root, 3)
            (id in cell_to_children(cell_to_parent(id, 2), 3)) || (nbadround += 1)
            (id in cell_to_children(cell_to_parent(id, 2))) || (nbadround += 1)
            (id in cell_to_children(cell_to_parent(id, 1), 3)) || (nbadround += 1)
        end
        @test nbadround == 0

        # argument validation
        @test_throws InvalidZ7Error cell_to_children(z7_from_string("001"), 0)
        @test_throws InvalidZ7Error cell_to_children(root, -1)
        @test_throws InvalidZ7Error cell_to_children(root, 20)
        @test_throws InvalidZ7Error cell_to_children(root, MAX_RESOLUTION + 1)
        deep = _ix_pent(1, 19)
        @test get_resolution(deep) == 19
        @test_throws InvalidZ7Error cell_to_children(deep)          # would be res 20
        @test collect(cell_to_children(deep, 19)) == [deep]
        res20 = z7_from_string("01" * repeat("0", 20))
        @test_throws InvalidZ7Error cell_to_children(res20)
        @test_throws InvalidZ7Error cell_to_children(res20, 20)
    end

    # -----------------------------------------------------------------------
    @testset "4. cell_to_parent" begin
        z = z7_from_string("0800433")
        @test cell_to_parent(z, 5) == z
        @test cell_to_parent(z, 4) == z7_from_string("080043")
        @test cell_to_parent(z, 0) == _ix_root(8)
        @test cell_to_parent(z) == z7_from_string("080043")
        @test cell_to_parent(z, 3) == z7_parent(z, 3)
        @test_throws InvalidZ7Error cell_to_parent(z, 6)
        @test_throws InvalidZ7Error cell_to_parent(z, -1)
        @test_throws InvalidZ7Error cell_to_parent(_ix_root(0))       # no parent
        @test_throws InvalidZ7Error cell_to_parent(z7_from_string("01" * repeat("0", 20)), 3)
        # ancestors of every res-3 cell of one base
        nbad = 0
        for id in cell_to_children(_ix_root(3), 3), r in 0:3
            p = cell_to_parent(id, r)
            (get_resolution(p) == r && z7_is_descendant(id, p)) || (nbad += 1)
        end
        @test nbad == 0
    end

    # -----------------------------------------------------------------------
    @testset "5. index <-> cell round trip, exhaustive res 0:4" begin
        for res in 0:4
            expected = _ix_all_ids(res)
            n = Int(num_cells(res))
            @test length(expected) == n
            @test issorted(expected)                       # canonical order
            @test allunique(expected)
            @test index_to_cell.(1:n, res) == expected
            @test cell_to_index.(expected) == 1:n
            # the same set, assembled through the hierarchy API
            assembled = UInt64[]
            for root in res0_cells()
                append!(assembled, cell_to_children(root, res))
            end
            @test assembled == expected
            # each base owns a contiguous block starting at base·p(res) + 1
            @test all(cell_to_index(_ix_pent(b, res)) == b * _ix_pent_count(res) + 1
                      for b in 0:11)
            @test all(z7_base_cell(index_to_cell(b * _ix_pent_count(res) + 1, res)) == b
                      for b in 0:11)
        end
    end

    # -----------------------------------------------------------------------
    @testset "6. spot checks at res 5, 18, 19" begin
        for res in (5, 18, 19)
            n = Int(num_cells(res))
            probes = unique((1, 2, n ÷ 2, n - 1, n))
            for index in probes
                id = index_to_cell(index, res)
                @test get_resolution(id) == res
                @test is_valid_cell(id)
                @test cell_to_index(id) == index
            end
            @test issorted([index_to_cell(i, res) for i in probes])
            @test index_to_cell(1, res) == _ix_pent(0, res)
            @test is_pentagon(index_to_cell(1, res))       # the all-zero chain
            @test z7_base_cell(index_to_cell(n, res)) == 11
            @test cell_to_index(_ix_pent(11, res)) == 11 * _ix_pent_count(res) + 1
        end
        # res 5 exhaustive over base 0's subtree (14,006 cells); the full
        # 168,072-cell sweep is covered by res 0:4 and the spot probes
        base0 = cell_to_children(_ix_root(0), 5)
        @test length(base0) == _ix_pent_count(5)
        @test cell_to_index.(base0) == 1:length(base0)
        @test index_to_cell.(1:length(base0), 5) == base0
    end

    # -----------------------------------------------------------------------
    @testset "7. dense-index bounds" begin
        for res in (0, 1, 2, 5, 19)
            n = Int(num_cells(res))
            @test_throws BoundsError index_to_cell(0, res)
            @test_throws BoundsError index_to_cell(-1, res)
            @test_throws BoundsError index_to_cell(n + 1, res)
        end
        @test_throws BoundsError index_to_cell(0, 2)
        @test_throws BoundsError index_to_cell(493, 2)
        @test_throws InvalidZ7Error index_to_cell(1, 20)
        @test_throws InvalidZ7Error index_to_cell(1, -1)
        # cell_to_index rejects ids without geometry / structure
        @test_throws InvalidZ7Error cell_to_index(z7_from_string("01" * repeat("0", 20)))
        @test_throws InvalidZ7Error cell_to_index((UInt64(12) << 60) | 0x0fffffffffffffff)
    end

    # -----------------------------------------------------------------------
    @testset "8. monotonicity: index order == encoded-id order" begin
        nbad = 0
        npairs = 0
        for res in 0:19, _ in 1:40
            n = num_cells(res)
            a = index_to_cell(rand(_IX_RNG, 1:n), res)
            b = index_to_cell(rand(_IX_RNG, 1:n), res)
            ia = cell_to_index(a)
            ib = cell_to_index(b)
            npairs += 1
            ((a < b) == (ia < ib)) || (nbad += 1)
            ((a == b) == (ia == ib)) || (nbad += 1)
            (1 <= ia <= n) || (nbad += 1)
            index_to_cell(ia, res) == a || (nbad += 1)
            get_resolution(a) == res || (nbad += 1)
            # independent sibling-sum rank (guards the closed-form `nhex`)
            ia == _ix_rank(a) || (nbad += 1)
            ib == _ix_rank(b) || (nbad += 1)
        end
        @test npairs == 20 * 40
        @test nbad == 0

        # the same independent rank over every cell at res 0:3 and over the
        # pentagon chains (whose deleted digit is what `nhex` corrects for)
        nbadrank = 0
        nranked = 0
        for res in 0:3, index in 1:Int(num_cells(res))
            id = index_to_cell(index, res)
            nranked += 1
            cell_to_index(id) == _ix_rank(id) || (nbadrank += 1)
        end
        for b in 0:11, res in 0:19
            id = _ix_pent(b, res)
            nranked += 1
            cell_to_index(id) == _ix_rank(id) || (nbadrank += 1)
            cell_to_index(id) == b * _ix_pent_count(res) + 1 || (nbadrank += 1)
        end
        @test nranked == 12 + 72 + 492 + 3432 + 12 * 20
        @test nbadrank == 0
    end

    # -----------------------------------------------------------------------
    @testset "9. lonlat_to_index" begin
        points = ((-73.9857, 40.7484), (0.0, 0.0), (11.25, 58.28252559),
            (179.9, -12.3), (-180.0, 0.0), (0.0, 90.0), (0.0, -90.0),
            (23.5, -66.0), (-45.0, 45.0), (137.0, 35.0))
        resolutions = (0, 1, 4, 7, 12, 19)
        nbad = 0
        nchecked = 0
        for (lon, lat) in points, res in resolutions
            i = lonlat_to_index(lon, lat, res)
            nchecked += 1
            i == cell_to_index(lonlat_to_cell(lon, lat, res)) || (nbad += 1)
            i == cell_to_index(lonlat_to_z7(lon, lat, res)) || (nbad += 1)
            1 <= i <= num_cells(res) || (nbad += 1)
            index_to_cell(i, res) == lonlat_to_cell(lon, lat, res) || (nbad += 1)
        end
        @test nchecked == length(points) * length(resolutions)
        @test nbad == 0
        @test lonlat_to_index(-73.9857, 40.7484, 7) ==
              cell_to_index(lonlat_to_cell(-73.9857, 40.7484, 7))
        @test lonlat_to_index(0.0, 0.0, 3) isa Int
        @test_throws InvalidZ7Error lonlat_to_index(0.0, 0.0, 20)
        @test_throws InvalidZ7Error lonlat_to_index(0.0, 0.0, -1)
    end

    # -----------------------------------------------------------------------
    @testset "10. introspection: resolution, validity, pentagon, id identity" begin
        @test MAX_RESOLUTION == 19
        z = z7_from_string("0800433")
        @test get_resolution(z) == 5
        @test get_resolution(z) == z7_resolution(z)
        @test !is_pentagon(z)
        @test is_valid_cell(z)
        @test cell_to_z7(z) === z
        @test z7_to_cell(z) === z
        @test z7_to_cell(cell_to_z7(z)) === z

        @test is_pentagon(z7_from_string("00000"))
        @test !is_pentagon(z7_from_string("00001"))
        @test all(is_pentagon(_ix_pent(b, r)) for b in 0:11, r in 0:5)
        @test all(get_resolution(_ix_pent(b, r)) == r for b in 0:11, r in 0:5)

        # validity: structure + geometry resolution gate
        @test !is_valid_cell((UInt64(12) << 60) | 0x0fffffffffffffff)   # base 12
        @test !is_valid_cell(z7_from_string("01" * repeat("0", 20)))    # res 20
        @test is_valid_cell(_ix_pent(1, 19))                            # res 19
        malformed = z7_from_string("00") & ~(UInt64(0x07) << 54)        # gap in the pad
        @test !is_valid_cell(malformed)
        # base 00 with digits 0,2 — the deleted pentagon subsequence "002"
        deleted_chain = (UInt64(2) << 54) | ((UInt64(1) << 54) - UInt64(1))
        @test !is_valid_cell(deleted_chain)
        @test_throws InvalidZ7Error cell_to_index(deleted_chain)

        # geometry-gated wrappers throw on res-20 / malformed ids [contract]
        res20 = z7_from_string("01" * repeat("0", 20))
        @test_throws InvalidZ7Error get_resolution(res20)
        @test_throws InvalidZ7Error is_pentagon(res20)
        @test_throws InvalidZ7Error cell_to_z7(res20)
        @test_throws InvalidZ7Error z7_to_cell(res20)
        @test_throws InvalidZ7Error get_resolution((UInt64(12) << 60) | 0x0fffffffffffffff)
        @test_throws InvalidZ7Error get_resolution(malformed)

        # every res 0:3 cell: identities and agreement with the z7 layer
        nbad = 0
        for res in 0:3, index in 1:Int(num_cells(res))
            id = index_to_cell(index, res)
            is_valid_cell(id) || (nbad += 1)
            get_resolution(id) == res || (nbad += 1)
            cell_to_z7(id) === id || (nbad += 1)
            z7_to_cell(id) === id || (nbad += 1)
            is_pentagon(id) == z7_is_pentagon(id) || (nbad += 1)
            z7_from_string(z7_to_string(id)) === id || (nbad += 1)
        end
        @test nbad == 0
        # exactly 12 pentagons at every resolution
        @test all(count(is_pentagon, index_to_cell.(1:Int(num_cells(r)), r)) == 12
                  for r in 0:3)
    end

    # -----------------------------------------------------------------------
    @testset "11. allocation-free, type-stable scalar paths" begin
        z = z7_from_string("0800433" * repeat("1", 12))       # res 17
        @test get_resolution(z) == 17
        z19 = z7_from_string("0800433" * repeat("1", 14))     # res 19
        @test get_resolution(z19) == 19

        @test (@inferred cell_to_index(z19)) isa Int
        @test (@inferred index_to_cell(7, 3)) isa UInt64
        @test (@inferred lonlat_to_index(12.3, 45.6, 10)) isa Int
        @test (@inferred num_cells(7)) isa Int64
        @test (@inferred get_resolution(z)) isa Int
        @test (@inferred is_valid_cell(z)) isa Bool
        @test (@inferred is_pentagon(z)) isa Bool
        @test (@inferred cell_to_parent(z, 4)) isa UInt64
        @test (@inferred cell_to_z7(z)) isa UInt64
        @test (@inferred z7_to_cell(z)) isa UInt64
        @test (@inferred res0_cells()) isa _IX_SMALL.SmallList{12,UInt64}
        @test (@inferred cell_to_children(z7_from_string("001"))) isa
              _IX_SMALL.SmallList{7,UInt64}

        @test _ix_alloc(cell_to_index, z19) == 0
        @test _ix_alloc(cell_to_index, _ix_root(11)) == 0
        @test _ix_alloc(index_to_cell, 7, 3) == 0
        @test _ix_alloc(index_to_cell, num_cells(19), 19) == 0
        @test _ix_alloc(lonlat_to_index, 12.3, 45.6, 10) == 0
        @test _ix_alloc(lonlat_to_index, -91.0, -37.0, 19) == 0
        @test _ix_alloc(get_resolution, z) == 0
        @test _ix_alloc(num_cells, 12) == 0
        @test _ix_alloc(is_valid_cell, z) == 0
        @test _ix_alloc(is_pentagon, z) == 0
        @test _ix_alloc(cell_to_parent, z, 4) == 0
        @test _ix_alloc(cell_to_z7, z) == 0
        @test _ix_alloc(z7_to_cell, z) == 0
        @test _ix_alloc(res0_cells) == 0

        # timing report (informational; the perf pass is design task 10)
        function bench(f, n)
            t0 = time_ns()
            s = 0
            for _ in 1:n
                s += f()
            end
            return (time_ns() - t0) / n, s
        end
        n19 = Int(num_cells(19))
        t_c2i, _ = bench(() -> cell_to_index(z19), 200_000)
        t_i2c, _ = bench(() -> Int(index_to_cell(n19 ÷ 3, 19) % 7), 200_000)
        @printf("    cell_to_index res 19: %.0f ns/call   index_to_cell res 19: %.0f ns/call\n",
            t_c2i, t_i2c)
    end
end
