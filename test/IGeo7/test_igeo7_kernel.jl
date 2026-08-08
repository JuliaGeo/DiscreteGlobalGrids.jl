# test/IGeo7/test_igeo7_kernel.jl — operations-kernel wiring for `IGEO7DGGS`
# (`src/IGeo7/IGeo7Kernel.jl`), including the MANDATORY descendant_range /
# ordinal checklist and the CAP-VALIDATION union-ratio measurement the
# operation contracts in `src/core/kernel.jl` require.
#
# Ground truth is never the wiring itself:
#   * cell sets are enumerated through `cell_descendants` (digit-lexicographic
#     DFS in the native layer) and separately asserted ascending, so the dense
#     ordinal and the pruning range are each checked against an independent
#     enumeration;
#   * counts are anchored to `IGeo7.num_cells` (`10·7^r + 2`, pinned by
#     `test/IGeo7/vectors/num_cells.csv`);
#   * hierarchy/parent/child relations are cross-checked against the oracle
#     vectors in `test/IGeo7/vectors/`.
#
# Big loops accumulate a failure counter and assert once (the idiom of
# `test/IGeo7/test_indexing.jl`) so the suite's test count stays readable.
#
# The suite lives in its own module: the generic vocabulary the systems share
# (`cell_center`, `num_cells`, ...) must not collide with a sibling's.

module IGeo7KernelTests

using Test
using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
using DiscreteGlobalGrids: IGeo7
import ConservativeRegridding as CR
import GeometryOps as GO
import SparseArrays

const S = IGEO7DGGS()
const SD = GO.UnitSpherical.spherical_distance
const VECTORS = joinpath(@__DIR__, "vectors")

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

"Strictly ascending — the ordering `DGGSPartialGrid` demands of stored ids."
ascending(v) = issorted(v; lt=(<=))

"Read a headered CSV into a header vector and a vector of row vectors."
function read_csv(path)
    lines = readlines(path)
    header = split(strip(lines[1]), ',')
    rows = [split(strip(l), ',') for l in lines[2:end] if !isempty(strip(l))]
    return header, rows
end

"""Every valid id at `res`, ascending, enumerated through the kernel hierarchy
(roots ascend and each root's subtree is a contiguous block of the order), and
independent of the dense-ordinal arithmetic under test."""
function all_ids(res::Int)
    ids = UInt64[]
    sizehint!(ids, Int(IGeo7.num_cells(res)))
    for root in DGG.root_ids(S)
        append!(ids, DGG.cell_descendants(S, 0, root, res))
    end
    return ids
end

"res-`res` pentagon of base `b` — the all-zero digit string."
pentagon(b::Integer, res::Integer) =
    IGeo7.z7_from_string(lpad(string(b), 2, '0') * repeat("0", res))

"""The pentagon neighborhoods at `res`: for every base, the six children of the
res-`res-1` pentagon (the res-`res` pentagon itself plus the five hexagons that
ring it) — the cells whose geometry the missing seventh child distorts."""
function pentagon_neighborhood(res::Integer)
    ids = UInt64[]
    for b in 0:11
        append!(ids, DGG.cell_children(S, res - 1, pentagon(b, res - 1)))
    end
    return ids
end

"""A spread of cells at `res`: the twelve pentagons, both cells straddling each
base-block boundary, and an evenly spaced sweep of the dense order."""
function sample_ids(res::Integer; sweep::Int=24)
    n = IGeo7.num_cells(res)
    ids = UInt64[pentagon(b, res) for b in 0:11]
    block = n ÷ 12
    for b in 1:11
        push!(ids, DGG.ordinal_to_cell(S, res, b * block))
        push!(ids, DGG.ordinal_to_cell(S, res, b * block + 1))
    end
    for k in 1:sweep
        push!(ids, DGG.ordinal_to_cell(S, res, 1 + (k * (n - 1)) ÷ (sweep + 1)))
    end
    return sort!(unique!(ids))
end

@testset verbose = true "IGeo7 kernel wiring (IGEO7DGGS)" begin

    # ----------------------------------------------------------------------
    @testset "1. traits, counts and geometry wiring" begin
        @test DGG.cell_id_type(S) === UInt64
        @test has_ordinal_ids(S) == false          # structural Z7 ids
        @test has_descendant_ranges(S) == true     # enabled by the checklist below
        @test radix(S) == 7
        @test supports_prefix_ranges(S)

        # The root count is a verified trait (12, see testset 2), wired in
        # `src/core/systems/igeo7.jl` — it must agree with the native
        # enumeration `root_ids` hands out.
        @test length(DGG.root_ids(S)) == 12
        @test root_count(S) == 12
        @test max_level(S) == 19        # the Z7 encoding's deepest digit slot

        @test DGG.num_cells(S, 0) == 12
        @test DGG.num_cells(S, 2) == 492
        @test all(DGG.num_cells(S, r) == 10 * Int64(7)^r + 2 for r in 0:19)
        @test all(DGG.num_cells(S, r) === IGeo7.num_cells(r) for r in 0:19)
        @test DGG.num_cells(S, 19) == 113988951853731432

        # geometry: unit-sphere ring, native center, closed-ring kwarg
        root = DGG.root_ids(S)[1]
        hexagon = DGG.cell_descendants(S, 0, root, 1)[2]
        for (level, id, n) in ((0, root, 5), (1, hexagon, 6))
            open_ring = DGG.cell_boundary(S, level, id)
            closed = DGG.cell_boundary(S, level, id; closed=true)
            @test open_ring isa Vector{GO.UnitSphericalPoint{Float64}}
            @test length(open_ring) == n
            @test length(closed) == n + 1
            @test closed[end] == closed[1]
            @test closed[1:n] == open_ring
            @test all(p -> abs(sqrt(sum(abs2, p)) - 1) < 1e-12, open_ring)
            center = DGG.cell_center(S, level, id)
            @test abs(sqrt(sum(abs2, center)) - 1) < 1e-12
            # native center, not the kernel's boundary-mean fallback
            lon, lat = IGeo7.cell_center(id)
            @test all(isapprox.(Tuple(center), IGeo7.lonlat_to_xyz(lon, lat); atol=1e-14))
            cap = cell_cap(S, level, id)
            @test cap.point == center
            @test all(SD(center, p) <= cap.radius for p in open_ring)
            @test DGG.cell_polygon_unitsphere(S, level, id) isa GO.GI.Polygon
        end

        # generic grid types accept the wiring
        grid = DGGSGrid(S, 3)
        @test DGG.num_cells(grid.system, grid.level) == 3432
        sub = subtree_grid(S, root; root_level=0, leaf_level=3)
        @test length(sub.ids) == 286
        @test ascending(sub.ids)

        # the `IGeo7Lookup` convenience constructor (validity trusted)
        ids = collect(DGG.cell_descendants(S, 0, root, 2))
        lookup = IGeo7.IGeo7Lookups.IGeo7Lookup(ids; resolution=2, validate=true)
        pg = DGGSPartialGrid(lookup)
        @test pg.system === S
        @test pg.level == 2
        @test pg.ids === lookup.data
        @test DGGSPartialGrid(lookup; bucket_size=8).bucket_size == 8

        # a lookup treeifies straight through the generic grid
        @test DGG.Trees.treeify(GO.Spherical(), lookup) isa DGG.DGGSCursor
        @test DGG.Trees.ncells(DGG.Trees.treeify(GO.Spherical(), lookup)) == length(lookup)
    end

    # ----------------------------------------------------------------------
    @testset "2. checklist 1: root_ids ascending, bases 0:11" begin
        roots = DGG.root_ids(S)
        @test length(roots) == 12
        @test eltype(roots) === UInt64
        @test ascending(roots)
        @test allunique(roots)
        @test [IGeo7.z7_base_cell(z) for z in roots] == collect(0:11)
        @test all(IGeo7.get_resolution(z) == 0 for z in roots)
        @test all(IGeo7.is_pentagon, roots)
        @test all(IGeo7.is_valid_cell, roots)
        @test [DGG.cell_to_ordinal(S, 0, z) for z in roots] == collect(1:12)
        @test [DGG.ordinal_to_cell(S, 0, i) for i in 1:12] == collect(roots)

        # oracle: the res-0 set as ids (its row order is the sealed module's
        # own, not ours — compare as a sorted set)
        _, rows = read_csv(joinpath(VECTORS, "res0_cells.csv"))
        @test length(rows) == 12
        @test sort(parse.(UInt64, getindex.(rows, 2); base=16)) == collect(roots)
    end

    # ----------------------------------------------------------------------
    @testset "3. checklist 2: children/descendants ascending" begin
        nbad = 0
        nchecked = 0
        function check_ascending(level, id)
            children = DGG.cell_children(S, level, id)
            (ascending(children) && length(children) == (IGeo7.is_pentagon(id) ? 6 : 7)) ||
                (nbad += 1)
            all(c -> DGG.cell_parent(S, level + 1, c, level) == id, children) || (nbad += 1)
            for delta in 1:3
                d = DGG.cell_descendants(S, level, id, level + delta)
                ascending(d) || (nbad += 1)
                allunique(d) || (nbad += 1)
                length(d) == DGG.subtree_leaf_count(S, level, id, level + delta) || (nbad += 1)
                all(c -> IGeo7.get_resolution(c) == level + delta, d) || (nbad += 1)
                all(c -> DGG.cell_parent(S, level + delta, c, level) == id, d) || (nbad += 1)
                nchecked += 1
            end
            return nothing
        end

        for root in DGG.root_ids(S)                       # all 12 roots
            check_ascending(0, root)
        end
        nsampled = 0
        for res in (1, 2, 5)                              # hexagons AND pentagons
            for id in sample_ids(res; sweep=12)
                check_ascending(res, id)
                nsampled += 1
            end
        end
        @test nbad == 0
        @test nsampled >= 90
        @test nchecked == 3 * (12 + nsampled)
    end

    # ----------------------------------------------------------------------
    # Two-sided proof triple (`descendant_range` checklist item 3):
    #   (a) extrema(cell_descendants) == descendant_range   (tight endpoints)
    #   (b) consecutive sorted parents' ranges ordered and disjoint
    #   (c) Σ subtree sizes == num_cells(R) over a complete parent level
    # Jointly: every descendant is in the range, and no non-descendant of the
    # same resolution can be.
    # ----------------------------------------------------------------------
    @testset "4. checklist 3: two-sided range contract" begin
        @testset "complete parent level $level -> res $leaf" for (level, leaf) in
                                                                 ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4))
            parents = all_ids(level)
            @test length(parents) == DGG.num_cells(S, level)
            @test ascending(parents)

            total = Int64(0)
            union_ids = UInt64[]
            sizehint!(union_ids, Int(DGG.num_cells(S, leaf)))
            previous_hi = nothing
            tight = 0
            ordered = 0
            counted = 0
            for p in parents
                desc = DGG.cell_descendants(S, level, p, leaf)
                lo, hi = descendant_range(S, level, p, leaf)
                (minimum(desc), maximum(desc)) == (lo, hi) && (tight += 1)
                (previous_hi === nothing || previous_hi < lo) && (ordered += 1)
                length(desc) == DGG.subtree_leaf_count(S, level, p, leaf) && (counted += 1)
                previous_hi = hi
                total += length(desc)
                append!(union_ids, desc)
            end
            @test tight == length(parents)          # (a)
            @test ordered == length(parents)        # (b)
            @test counted == length(parents)
            @test total == DGG.num_cells(S, leaf)   # (c)
            # (a)+(b)+(c) => the subtrees tile the level exactly:
            @test union_ids == all_ids(leaf)
            @test ascending(union_ids)
        end

        @testset "pentagon neighborhood at res $level" for level in (4, 6)
            parents = pentagon_neighborhood(level)
            @test length(parents) == 72
            nbad = 0
            for leaf in (level + 1):(level + 3)
                for b in 0:11
                    siblings = parents[(6b + 1):(6b + 6)]      # ascending, one base
                    ascending(siblings) || (nbad += 1)
                    previous_hi = nothing
                    subtotal = Int64(0)
                    for p in siblings
                        desc = DGG.cell_descendants(S, level, p, leaf)
                        lo, hi = descendant_range(S, level, p, leaf)
                        (minimum(desc), maximum(desc)) == (lo, hi) || (nbad += 1)    # (a)
                        (previous_hi === nothing || previous_hi < lo) || (nbad += 1) # (b)
                        length(desc) == DGG.subtree_leaf_count(S, level, p, leaf) ||
                            (nbad += 1)
                        previous_hi = hi
                        subtotal += length(desc)
                    end
                    # (c) locally: the six siblings partition their parent
                    parent = pentagon(b, level - 1)
                    subtotal == DGG.subtree_leaf_count(S, level - 1, parent, leaf) ||
                        (nbad += 1)
                end
            end
            @test nbad == 0
        end
    end

    # ----------------------------------------------------------------------
    @testset "5. checklist 4: deep-delta endpoints without enumeration" begin
        nbad = 0
        nchecked = 0
        probes = UInt64[]
        append!(probes, DGG.root_ids(S))
        for res in (1, 2, 5, 10, 15, 18)
            append!(probes, sample_ids(res; sweep=6))
        end
        for id in probes
            res = IGeo7.get_resolution(id)
            for leaf in res:IGeo7.MAX_RESOLUTION
                lo, hi = descendant_range(S, res, id, leaf)
                (IGeo7.is_valid_cell(lo) && IGeo7.is_valid_cell(hi)) || (nbad += 1)
                (IGeo7.get_resolution(lo) == leaf && IGeo7.get_resolution(hi) == leaf) ||
                    (nbad += 1)
                (DGG.cell_parent(S, leaf, lo, res) == id) || (nbad += 1)
                (DGG.cell_parent(S, leaf, hi, res) == id) || (nbad += 1)
                (leaf == res ? (lo == id == hi) : (lo < hi)) || (nbad += 1)
                # the endpoints really are the extremal digit strings
                (IGeo7.z7_to_string(lo) ==
                 IGeo7.z7_to_string(id) * repeat("0", leaf - res)) || (nbad += 1)
                (IGeo7.z7_to_string(hi) ==
                 IGeo7.z7_to_string(id) * repeat("6", leaf - res)) || (nbad += 1)
                nchecked += 1
            end
        end
        @test nbad == 0
        @test nchecked > 500

        # next-sibling separation: consecutive children's ranges never touch,
        # at every leaf resolution down to 19 (no enumeration involved).
        nbad = 0
        for parent in vcat(collect(DGG.root_ids(S)), sample_ids(3; sweep=8))
            level = IGeo7.get_resolution(parent)
            children = DGG.cell_children(S, level, parent)
            for leaf in (level + 1):IGeo7.MAX_RESOLUTION
                previous_hi = nothing
                for c in children
                    lo, hi = descendant_range(S, level + 1, c, leaf)
                    (previous_hi === nothing || lo > previous_hi) || (nbad += 1)
                    previous_hi = hi
                end
                # the children's ranges partition the parent's range
                plo, phi = descendant_range(S, level, parent, leaf)
                (descendant_range(S, level + 1, first(children), leaf)[1] == plo &&
                 descendant_range(S, level + 1, last(children), leaf)[2] == phi) ||
                    (nbad += 1)
            end
        end
        @test nbad == 0

        # the `get_resolution(id) == level` bookkeeping assertion
        root = DGG.root_ids(S)[1]
        @test_throws ArgumentError descendant_range(S, 1, root, 3)
        @test_throws IGeo7.InvalidZ7Error descendant_range(S, 0, root, 20)
        # A negative delta is the kernel's argument error, not Z7 invalidity,
        # so it is an `ArgumentError` here exactly as it is on H3/A5/HEALPix.
        @test_throws ArgumentError descendant_range(S, 1, pentagon(0, 1), 0)
    end

    # ----------------------------------------------------------------------
    @testset "6. checklist 5: ordinal monotonicity" begin
        @testset "exhaustive at res $res" for res in 0:3
            ids = all_ids(res)
            n = Int(DGG.num_cells(S, res))
            @test length(ids) == n
            @test ascending(ids)
            nbad = 0
            for (i, id) in enumerate(ids)
                (DGG.cell_to_ordinal(S, res, id) == i &&
                 DGG.ordinal_to_cell(S, res, i) == id) || (nbad += 1)
            end
            @test nbad == 0
        end

        @testset "sampled at res $res" for res in 5:9
            ids = sample_ids(res; sweep=40)
            append!(ids, DGG.cell_descendants(S, res - 2, pentagon(0, res - 2), res))
            append!(ids, DGG.cell_descendants(S, res - 2, pentagon(7, res - 2), res))
            sort!(unique!(ids))
            ordinals = [DGG.cell_to_ordinal(S, res, id) for id in ids]
            @test ascending(ordinals)                        # monotone in the id
            @test all(1 .<= ordinals .<= DGG.num_cells(S, res))
            @test [DGG.ordinal_to_cell(S, res, o) for o in ordinals] == ids
            # ordinal -> cell -> ordinal on the same sample
            @test [DGG.cell_to_ordinal(S, res, DGG.ordinal_to_cell(S, res, o))
                   for o in ordinals] == ordinals
        end

        # ordinal space is contiguous over a descendant_range (what the dense
        # cursor's O(1) leaf interval relies on)
        nbad = 0
        for level in 0:2, leaf in (level + 1):(level + 3)
            for p in all_ids(level)
                lo, hi = descendant_range(S, level, p, leaf)
                width = DGG.cell_to_ordinal(S, leaf, hi) -
                        DGG.cell_to_ordinal(S, leaf, lo) + 1
                width == DGG.subtree_leaf_count(S, level, p, leaf) || (nbad += 1)
            end
        end
        @test nbad == 0

        # An ordinal outside `1:num_cells` is the kernel's uniform
        # `OrdinalRangeError`, not the `BoundsError` the native
        # `index_to_cell` underneath still raises — that one is pinned as
        # native API shape (`test/IGeo7/test_indexing.jl`), so the wiring
        # range-checks ahead of it rather than relaying it.
        n = Int(DGG.num_cells(S, 2))
        @test_throws OrdinalRangeError DGG.ordinal_to_cell(S, 2, 0)
        @test_throws OrdinalRangeError DGG.ordinal_to_cell(S, 2, n + 1)
        @test_throws BoundsError IGeo7.index_to_cell(n + 1, 2)
        err = try
            DGG.ordinal_to_cell(S, 2, n + 1)
        catch e
            e
        end
        @test err.system === :IGEO7_ISEA7H_Z7
        @test err.level == 2
        @test err.ordinal == n + 1
        @test err.total == n
        @test occursin("IGEO7_ISEA7H_Z7", sprint(showerror, err))
    end

    # ----------------------------------------------------------------------
    @testset "7. subtree_leaf_count is O(1) and exact" begin
        # exact equality against the enumerated children, hexagons and
        # pentagons, deltas 1-4
        nbad = 0
        ncells = 0
        for res in (0, 1, 2, 5)
            for id in (res == 0 ? collect(DGG.root_ids(S)) : sample_ids(res; sweep=4))
                for delta in 1:4
                    length(DGG.cell_descendants(S, res, id, res + delta)) ==
                    DGG.subtree_leaf_count(S, res, id, res + delta) || (nbad += 1)
                end
                # the two closed forms
                for delta in 0:4
                    expected = IGeo7.is_pentagon(id) ?
                               (5 * Int64(7)^delta + 1) ÷ 6 : Int64(7)^delta
                    DGG.subtree_leaf_count(S, res, id, res + delta) == expected ||
                        (nbad += 1)
                end
                ncells += 1
            end
        end
        @test nbad == 0
        @test ncells >= 40

        # anchored to num_cells: a complete level decomposes into 12 pentagon
        # subtrees and hexagon subtrees for the rest
        for r in 1:6
            @test sum(DGG.subtree_leaf_count(S, 0, root, r) for root in DGG.root_ids(S)) ==
                  DGG.num_cells(S, r)
            @test sum(DGG.subtree_leaf_count(S, 1, id, r + 1) for id in all_ids(1)) ==
                  DGG.num_cells(S, r + 1)
            pent = (5 * Int64(7)^r + 1) ÷ 6
            @test 12 * pent == DGG.num_cells(S, r)
            # the res-1 level is 12 pentagons + 60 hexagons
            @test 12 * ((5 * Int64(7)^(r - 1) + 1) ÷ 6) + 60 * Int64(7)^(r - 1) ==
                  DGG.num_cells(S, r)
        end
        @test DGG.subtree_leaf_count(S, 0, DGG.root_ids(S)[1], 19) ==
              113988951853731432 ÷ 12
        @test DGG.subtree_leaf_count(S, 0, DGG.root_ids(S)[1], 0) == 1

        root = DGG.root_ids(S)[1]
        @test_throws IGeo7.InvalidZ7Error DGG.subtree_leaf_count(S, 1, pentagon(0, 1), 0)
        @test_throws IGeo7.InvalidZ7Error DGG.subtree_leaf_count(S, 0, root, 20)

        # Error-type contract (src/core/kernel.jl): a negative delta is a
        # kernel-level argument mistake and is an `ArgumentError` on every
        # system, IGEO7 included — the native `InvalidZ7Error` the Z7 layer
        # would raise there is right about Z7 and wrong about whose mistake it
        # is. A `leaf_level` past the encoding's own maximum is genuine Z7
        # invalidity and keeps the native type.
        @test_throws ArgumentError DGG.cell_descendants(S, 1, pentagon(0, 1), 0)
        @test_throws ArgumentError DGG.cell_descendants(S, 3, root, 2)
        @test_throws IGeo7.InvalidZ7Error DGG.cell_descendants(S, 0, root, 20)
        @test_throws IGeo7.InvalidZ7Error DGG.cell_descendants(S, 0, UInt64(0), 2)
    end

    # ----------------------------------------------------------------------
    @testset "8. oracle vectors (test/IGeo7/vectors)" begin
        header, rows = read_csv(joinpath(VECTORS, "num_cells.csv"))
        @test header == ["res", "count"]
        @test all(DGG.num_cells(S, parse(Int, r[1])) == parse(Int64, r[2]) for r in rows)

        # hierarchy.csv: parent chain, children, and the range/ordinal
        # relations they imply
        header, rows = read_csv(joinpath(VECTORS, "hierarchy.csv"))
        @test header == ["z7_hex", "z7_string", "parent_chain_strings", "children_strings"]
        @test !isempty(rows)
        nbad = 0
        for row in rows
            id = IGeo7.z7_from_hex(String(row[1]))
            res = IGeo7.get_resolution(id)
            parents = String.(split(row[3], ';'))
            children = String.(split(row[4], ';'))

            # cell_children == the oracle's children (ascending)
            wired = DGG.cell_children(S, res, id)
            IGeo7.z7_to_string.(collect(wired)) == children || (nbad += 1)
            ascending(wired) || (nbad += 1)

            # cell_parent along the whole chain, coarser by one per entry
            [IGeo7.z7_to_string(DGG.cell_parent(S, res, id, r)) for r in (res - 1):-1:0] ==
            parents || (nbad += 1)

            # every ancestor's descendant_range brackets the cell, and its
            # ordinal interval brackets the cell's ordinal
            ordinal = DGG.cell_to_ordinal(S, res, id)
            for (k, text) in enumerate(parents)
                ancestor = IGeo7.z7_from_string(text)
                level = res - k
                lo, hi = descendant_range(S, level, ancestor, res)
                (lo <= id <= hi) || (nbad += 1)
                (DGG.cell_to_ordinal(S, res, lo) <= ordinal <=
                 DGG.cell_to_ordinal(S, res, hi)) || (nbad += 1)
                (DGG.cell_to_ordinal(S, res, hi) - DGG.cell_to_ordinal(S, res, lo) + 1 ==
                 DGG.subtree_leaf_count(S, level, ancestor, res)) || (nbad += 1)
            end
        end
        @test nbad == 0

        # pentagon_chains.csv: the six-child pentagons at every chain depth
        header, rows = read_csv(joinpath(VECTORS, "pentagon_chains.csv"))
        @test header == ["base", "depth", "pentagon_z7_string", "children_z7_strings",
            "missing_digits"]
        nbad = 0
        for row in rows
            id = IGeo7.z7_from_string(String(row[3]))
            res = IGeo7.get_resolution(id)
            expected = String.(split(row[4], ';'))
            wired = DGG.cell_children(S, res, id)
            (length(wired) == 6 && IGeo7.z7_to_string.(collect(wired)) == expected) ||
                (nbad += 1)
            DGG.subtree_leaf_count(S, res, id, res + 1) == 6 || (nbad += 1)
            DGG.subtree_leaf_count(S, res, id, res + 2) == 41 || (nbad += 1)
            # the deleted digit leaves no hole in the range: lo/hi still tight
            (minimum(wired), maximum(wired)) == descendant_range(S, res, id, res + 1) ||
                (nbad += 1)
        end
        @test nbad == 0
    end

    # ----------------------------------------------------------------------
    # CAP-VALIDATION (`cell_cap` in `src/core/kernel.jl`). "Union ratio" = max
    # distance from the wired cap's center to any delta-level descendant
    # vertex, divided by the cell's own max center-to-vertex distance — i.e.
    # the inflation factor a subtree actually needs. `CELL_CAP_INFLATION = 1.2`
    # is the wired budget.
    # ----------------------------------------------------------------------
    @testset "9. CAP-VALIDATION: subtree union ratios" begin
        deltas = 1:5

        "cap exactly as wired, plus the *uninflated* max center-to-vertex distance"
        function cap_parts(level, id)
            cap = cell_cap(S, level, id)
            raw = maximum(SD(cap.point, p) for p in DGG.cell_boundary(S, level, id))
            return cap, raw
        end

        function union_ratios(cells, level)
            ratios = zeros(length(deltas))
            of_cap = 0.0
            for id in cells
                cap, raw = cap_parts(level, id)
                for (k, delta) in enumerate(deltas)
                    leaf = level + delta
                    worst = 0.0
                    for c in DGG.cell_descendants(S, level, id, leaf)
                        for p in DGG.cell_boundary(S, leaf, c)
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

        groups = (
            ("res 0 (all 12 cells)", collect(DGG.root_ids(S)), 0),
            ("res 1 (all 72 cells)", all_ids(1), 1),
            ("res 4 pentagon nbhd", pentagon_neighborhood(4), 4),
            ("res 6 pentagon nbhd", pentagon_neighborhood(6), 6),
        )

        println("\n  IGeo7 CAP-VALIDATION — union ratio (max descendant-vertex " *
                "distance / the cell's own radius before the 1.2 inflation)")
        println("  ", rpad("group", 22), join(lpad("d=$d", 11) for d in deltas),
            lpad("in-cap", 11))
        worst_ratio = 0.0
        worst_of_cap = 0.0
        worst_extrapolated = 0.0
        convergence_bad = 0
        for (label, cells, level) in groups
            ratios, of_cap = union_ratios(cells, level)
            println("  ", rpad(label, 22),
                join(lpad(string(round(r; digits=5)), 11) for r in ratios),
                lpad(string(round(of_cap; digits=5)), 11))
            increments = [ratios[k] - ratios[k - 1] for k in 2:length(ratios)]
            println("  ", rpad("  increments", 22), lpad("", 11),
                join(lpad(string(round(i; digits=6)), 11) for i in increments))
            # Convergence is two-step: the Eisenstein lattice's chirality
            # alternates with resolution parity, so a subtree bulges in
            # alternating directions and the increments alternate large/small.
            # Each SAME-PARITY increment must shrink by >= 2x (the 1e-4 floor
            # covers a plateaued maximum, where the previous increment is 0;
            # it is ~2% of the total delta-1..5 growth, far too small to hide
            # divergence).
            for k in 3:length(increments)
                increments[k] <= max(increments[k - 2] / 2, 1e-4) || (convergence_bad += 1)
            end
            # Conservative geometric tail beyond delta 5 at ratio 1/2 per two
            # deltas: sup <= r5 + (inc4 + inc5) * (1/2)/(1 - 1/2).
            worst_extrapolated = max(worst_extrapolated,
                ratios[end] + increments[end - 1] + increments[end])
            worst_ratio = max(worst_ratio, maximum(ratios))
            worst_of_cap = max(worst_of_cap, of_cap)
        end
        println("  worst union ratio ", round(worst_ratio; digits=5),
            " | worst fraction of the wired cap radius ", round(worst_of_cap; digits=5),
            " | extrapolated supremum ", round(worst_extrapolated; digits=5), "\n")

        @test worst_ratio <= 1.10                    # CONTRACT bar (ii)
        @test convergence_bad == 0                   # geometric convergence
        @test worst_extrapolated <= 1.2 * 0.95       # extrapolated < the 1.2 budget
        @test worst_of_cap < 1.0                     # every descendant is inside
    end

    # ----------------------------------------------------------------------
    # `ConservativeRegridding.Regridder` over the generic grids — the consumer
    # the tree layer exists for, and the one case in the package where the two
    # sides are at DIFFERENT resolutions, so the answer is a real partition
    # rather than an identity.
    # ----------------------------------------------------------------------
    @testset "10. Regridder round trips" begin
        dense = CR.Regridder(DGGSGrid(S, 0), DGGSGrid(S, 0); threaded=false, normalize=false)
        @test size(dense.intersections) == (12, 12)
        @test SparseArrays.nnz(dense.intersections) == 12        # diagonal only
        @test isapprox(sum(dense.intersections), sum(dense.dst_areas); rtol=1e-12)

        root = first(DGG.root_ids(S))
        lookup = IGeo7.IGeo7Lookups.IGeo7Lookup(collect(DGG.cell_children(S, 0, root));
            resolution=1, validate=true)
        partial = CR.Regridder(DGGSPartialGrid(lookup), DGGSPartialGrid(lookup);
            threaded=false, normalize=false)
        @test size(partial.intersections) == (length(lookup), length(lookup))
        @test isapprox(sum(partial.intersections), sum(partial.dst_areas); rtol=1e-12)

        # CROSS-LEVEL conservation: res 0 destination, res 1 source. Both tile
        # the same sphere, every res-0 cell is exactly covered by the res-1
        # cells it meets, and every res-1 cell is fully consumed.
        cross = CR.Regridder(DGGSGrid(S, 0), DGGSGrid(S, 1); threaded=false, normalize=false)
        @test size(cross.intersections) == (12, Int(IGeo7.num_cells(1)))
        @test isapprox(sum(cross.dst_areas), 4pi * GO.Spherical().radius^2; rtol=1e-9)
        @test isapprox(sum(cross.src_areas), sum(cross.dst_areas); rtol=1e-9)
        @test isapprox(vec(sum(cross.intersections; dims=2)), cross.dst_areas; rtol=1e-9)
        @test isapprox(vec(sum(cross.intersections; dims=1)), cross.src_areas; rtol=1e-9)
    end
end

end # module IGeo7KernelTests
