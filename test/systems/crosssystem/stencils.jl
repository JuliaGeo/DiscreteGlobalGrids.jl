# Cross-system laws for subset adjacency: clipped rings, and the cached
# `adjacency` table in all three of its row shapes.
# Brute-force oracles derive rows from native one-rings and explicit membership,
# independently of the table implementations.

module StencilTests

using Test
import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import GeometryOps as GO

const FB = DGG.Fallbacks
const EN = DGG.Engine

using DiscreteGlobalGrids: systems, levelgrid, ncells, cellindex, cellposition,
    neighbors, ring, level, levels, maxlevel, descendants,
    descendant_range, has_sorted_subtrees, PartialGrid, CellVector, CellLookup,
    MultiOrderCoverage, member_neighbors, Vertex, Edge,
    Connectivity, cellindextype, query, system,
    subtree, halo, adjacency, AdjacencyTable, halopositions, halocells

include(joinpath(@__DIR__, "..", "..", "helpers.jl"))
using .DGGTestHelpers: syslabel, isquadface, sweepcovers

include(joinpath(@__DIR__, "..", "..", "helpers.jl"))
using .DGGTestHelpers: syslabel, isquadface, sweepcovers

# ---------------------------------------------------------------------------
# Systems, and the depths each is swept at
#
# A root two levels above the leaf gives a subtree big enough to have an
# interior and a border, and small enough that the brute-force oracles below can
# expand everything. The `AuthalicSystem` wrap is over IGEO7, as in the
# neighbouring files.
# ---------------------------------------------------------------------------

const SWEEP = [
    (DGG.IGeo7System(), 1, 3),
    (DGG.H3System(), 1, 3),
    (DGG.HEALPixSystem(), 1, 4),
    (DGG.A5System(), 1, 3),
    (DGG.S2System(), 1, 4),
    (DGG.ISEA4RSystem(), 1, 4),
    (DGG.AuthalicSystem(DGG.IGeo7System()), 1, 3),
]

@testset "the sweep covers every registered system" begin
    sweepcovers(SWEEP)
end

# ---------------------------------------------------------------------------
# The subset shapes
# ---------------------------------------------------------------------------

held(sub, c) = cellposition(sub, c) !== nothing

# A rooted subtree, the shape `subtree(sys, cell, level)` builds.
rooted(sys, base, leaf) =
    subtree(sys, cellindex(levelgrid(sys, base), 3), leaf)

# A scattered set with no root at all: every fifth cell of the whole level.
# Most neighbours are absent, which is the case the clip has to get right.
function scattered(sys, leaf)
    grid = levelgrid(sys, leaf)
    ids = [cellindex(grid, i) for i in 1:5:ncells(grid)]
    return PartialGrid(sys, leaf, ids)
end

# The same subtree with its middle third removed — a hole, and the shape the
# two readings of "distance" part company over.
function holed(sys, base, leaf)
    root = cellindex(levelgrid(sys, base), 3)
    ids = descendants(sys, root, leaf)
    n = length(ids)
    keep = vcat(ids[1:(n÷3)], ids[(2n÷3+1):end])
    return PartialGrid(sys, leaf, keep)
end

shapes(sys, base, leaf) = ("rooted subtree" => rooted(sys, base, leaf),
    "scattered subset" => scattered(sys, leaf),
    "subtree with a hole" => holed(sys, base, leaf))

# A deterministic spread of in-set probes.
function probes(sub, n::Int)
    total = ncells(sub)
    step = max(1, total ÷ n)
    return [cellindex(sub, i) for i in 1:step:total]
end

function induced_ring(sub, complete, c, k::Int; connectivity = Vertex())
    seen = Set([c])
    frontier = [c]
    for _ in 1:k
        next = eltype(frontier)[]
        for x in frontier, y in neighbors(complete, x, 1; connectivity)
            (held(sub, y) && !(y in seen)) || continue
            push!(seen, y)
            push!(next, y)
        end
        frontier = next
    end
    return Set(frontier)
end

# ---------------------------------------------------------------------------
# CLIPPING, and the faces that must agree with it
# ---------------------------------------------------------------------------

@testset "$(syslabel(sys))" for (sys, base, leaf) in SWEEP
    complete = levelgrid(sys, leaf)

    @testset "$label: ring is the complete level's, clipped" for (label, sub) in
                                                                 shapes(sys, base, leaf)
        for c in probes(sub, 6), conn in (Vertex(), Edge()), k in 0:3
            want = [x for x in ring(complete, c, k; connectivity = conn) if held(sub, x)]
            @test ring(sub, c, k; connectivity = conn) == want
        end
        # ...and the disc is still its rings concatenated outward, which the
        # clip preserves because it preserves order.
        for c in probes(sub, 4), conn in (Vertex(), Edge()), k in 0:3
            rings = [ring(sub, c, j; connectivity = conn) for j in 1:k]
            @test neighbors(sub, c, k; connectivity = conn) ==
                  reduce(vcat, rings; init = cellindextype(sys)[])
        end
        cv = CellVector(sub)
        for c in probes(sub, 3), k in 0:2
            @test ring(sub, c, k) == filter(in(sub), ring(complete, c, k))
            @test ring(cv, c, k) == filter(in(cv), ring(complete, c, k))
            @test ring(CellLookup(cv), c, k) ==
                  filter(in(CellLookup(cv)), ring(complete, c, k))
        end
    end

    @testset "$label: the position face is the id face" for (label, sub) in
                                                            shapes(sys, base, leaf)
        for p in 1:max(1, ncells(sub) ÷ 6):ncells(sub), conn in (Vertex(), Edge()), k in 0:2
            c = cellindex(sub, p)
            # Element for element, not as sets: the position form IS the id
            # form read through `cellposition`, and it carries the same
            # counter-clockwise order.
            @test neighbors(sub, p, k; connectivity = conn) ==
                  [cellposition(sub, x)
                   for x in neighbors(sub, c, k; connectivity = conn)]
            @test ring(sub, p, k; connectivity = conn) ==
                  [cellposition(sub, x)
                   for x in ring(sub, c, k; connectivity = conn)]
        end
    end

    @testset "$label: CellVector and CellLookup answer the same" for (label, sub) in
                                                                     shapes(sys, base, leaf)
        cv = CellVector(sub)
        lk = CellLookup(cv)
        for c in probes(sub, 4), conn in (Vertex(), Edge()), k in 0:2
            want = neighbors(sub, c, k; connectivity = conn)
            @test neighbors(cv, c, k; connectivity = conn) == want
            @test neighbors(lk, c, k; connectivity = conn) == want
            @test ring(cv, c, k; connectivity = conn) ==
                  ring(sub, c, k; connectivity = conn)
            @test ring(lk, c, k; connectivity = conn) ==
                  ring(sub, c, k; connectivity = conn)
        end
        for p in 1:max(1, ncells(sub) ÷ 4):ncells(sub)
            @test neighbors(cv, p, 2) == neighbors(sub, p, 2)
            @test neighbors(lk, p, 2) == neighbors(sub, p, 2)
        end
    end

    @testset "$label: the clipped table is the per-cell calls" for (label, sub) in
                                                                   shapes(sys, base, leaf)
        for conn in (Vertex(), Edge())
            table = adjacency(sub; connectivity = conn)
            @test length(table) == ncells(sub)
            @test all(p -> table[p] == neighbors(sub, p, 1; connectivity = conn),
                1:ncells(sub))
        end
        cv = CellVector(sub)
        @test adjacency(cv) == adjacency(sub)
        @test adjacency(CellLookup(cv)) == adjacency(sub)
        # The rows are ROTATIONAL, and that is only visible as an absence of
        # sorting: a `sort!` put back anywhere on this path would pass every
        # other assertion in this file. Some row of a real subset is out of
        # ascending order on every system.
        @test any(!issorted, adjacency(sub))
        @test any(!issorted, adjacency(cv))
        # Threading is a scheduling choice, never an answer.
        @test adjacency(sub; threaded = false) == adjacency(sub)
    end

    @testset "the rooted fast path agrees with the generic route" begin
        sub = rooted(sys, base, leaf)
        for conn in (Vertex(), Edge())
            @test adjacency(sub; connectivity = conn) ==
                  [neighbors(sub, p, 1; connectivity = conn) for p in 1:ncells(sub)]
        end
        @test (EN._whole_subtree_range(sub) !== nothing) == has_sorted_subtrees(sys)
        @test EN._whole_subtree_range(scattered(sys, leaf)) === nothing
        # Rooted but not WHOLE: the root is there, one descendant is not, and
        # "interior" stops meaning "every neighbour is present". The generic
        # route is the only correct one, and the gate has to say so.
        root = cellindex(levelgrid(sys, base), 3)
        ids = descendants(sys, root, leaf)
        part = PartialGrid(sys, leaf, ids[1:(end-1)]; root)
        @test EN._whole_subtree_range(part) === nothing
        @test adjacency(part) == [neighbors(part, p, 1) for p in 1:ncells(part)]
    end

    @testset "a cell outside the subset has no neighbourhood" begin
        sub = rooted(sys, base, leaf)
        cv = CellVector(sub)
        outside = nothing
        for i in 1:ncells(complete)
            c = cellindex(complete, i)
            held(sub, c) || (outside = c; break)
        end
        @test outside !== nothing
        @test_throws ArgumentError neighbors(sub, outside)
        @test_throws ArgumentError ring(sub, outside, 2)
        @test_throws ArgumentError neighbors(cv, outside)
        @test_throws ArgumentError ring(cv, outside, 2)
        @test_throws ArgumentError neighbors(CellLookup(cv), outside)
        # And the position face keeps Base's contract instead.
        @test_throws BoundsError neighbors(sub, ncells(sub) + 1)
        # Directly on `cellindex`, because that is where it is enforced: a
        # subtree's ids are a WINDOW into the level, so an unchecked position
        # past the end reads a real cell of the next subtree and the position
        # form would answer confidently about a cell the subset does not hold.
        @test_throws BoundsError cellindex(sub, 0)
        @test_throws BoundsError cellindex(sub, ncells(sub) + 1)
    end
end

# ---------------------------------------------------------------------------
# Compare `neighborcount` with every ring, including cells below the maximum
# degree.
# ---------------------------------------------------------------------------

@testset "neighborcount is the ring's length, exceptional cells included" begin
    for (sys, base, leaf) in SWEEP
        mismatched = 0
        offdegree = 0
        bound = DGG.maxneighbors(sys, Vertex())
        for l in base:leaf
            grid = levelgrid(sys, l)
            for i in 1:ncells(grid)
                c = cellindex(grid, i)
                for conn in (Vertex(), Edge())
                    DGG.neighborcount(grid, c; connectivity = conn) ==
                        length(neighbors(grid, c, 1; connectivity = conn)) ||
                        (mismatched += 1)
                end
                DGG.neighborcount(grid, c) < bound && (offdegree += 1)
            end
        end
        @test mismatched == 0
        @test offdegree > 0
    end
    # Subsets count clipped rings and reject cells outside the subset.
    sys, base, leaf = SWEEP[1]
    sub = rooted(sys, base, leaf)
    cv = CellVector(sub)
    onborder = first(c for c in (cellindex(sub, i) for i in 1:ncells(sub))
                  if length(neighbors(sub, c)) < DGG.neighborcount(levelgrid(sys, leaf), c))
    @test DGG.neighborcount(sub, onborder) == length(neighbors(sub, onborder))
    @test DGG.neighborcount(cv, onborder) == length(neighbors(cv, onborder))
    @test DGG.neighborcount(CellLookup(cv), onborder) == length(neighbors(cv, onborder))
    outside = cellindex(levelgrid(sys, leaf), ncells(levelgrid(sys, leaf)))
    @test_throws ArgumentError DGG.neighborcount(cv, outside)
end

# ---------------------------------------------------------------------------
# DIVERGENCE
#
# The concrete case: a subtree with its middle third removed. A cell in the
# first third and a cell in the last third can be two SYSTEM steps apart while
# every path between them inside the subset is longer, or absent altogether —
# so `ring(sub, c, 2)` names cells that a breadth-first walk of the subset's own
# adjacency graph puts at distance 3 or at no distance at all. The clipped
# reading keeps them; the induced-subgraph reading drops them. They agree at
# `k == 1` on every subset, because one step cannot leave and come back.
# ---------------------------------------------------------------------------

@testset "clipped distance is not breadth-first distance inside the subset" begin
    divergences = 0
    for (sys, base, leaf) in SWEEP
        complete = levelgrid(sys, leaf)
        for (label, sub) in shapes(sys, base, leaf)
            for c in probes(sub, 8)
                # k == 1: the two readings coincide, always and everywhere.
                @test Set(ring(sub, c, 1)) == induced_ring(sub, complete, c, 1)
                for k in 2:3
                    clipped = Set(ring(sub, c, k))
                    induced = induced_ring(sub, complete, c, k)
                    # One direction is a law rather than a difference: leaving
                    # the subset can only SHORTEN a path, so system distance is
                    # never more than induced distance and a cell the induced
                    # walk puts at `k` is inside the clipped disc of radius `k`.
                    @test induced ⊆ Set(neighbors(sub, c, k))
                    clipped == induced || (divergences += 1)
                end
            end
        end
    end
    # Not "some system somewhere might differ" — the hole makes it happen, and
    # a return to induced-subgraph semantics would make this zero.
    @test divergences > 0
end


# `[chunk; halo]` as ids, and each neighbour looked up as itself. `0` means the
# buffer does not contain it at all, which is the row that could not have been
# completed.
function buffer_rows(sys, l, chunk_ids, halo_ids, conn)
    grid = levelgrid(sys, l)
    slot = Dict{eltype(chunk_ids),Int}()
    for (k, c) in enumerate(chunk_ids)
        slot[c] = k
    end
    for (k, c) in enumerate(halo_ids)
        slot[c] = length(chunk_ids) + k
    end
    return [[get(slot, nb, 0) for nb in neighbors(grid, d, 1; connectivity = conn)]
            for d in chunk_ids]
end

@testset "buffered rows are complete and address [region; halo]" begin
    uncovered = Any[]      # the halo did not cover some neighbour
    disagreed = Any[]      # the table is not the oracle's rows
    layout = Any[]         # row i is not `descendant_range` offset i
    csr = Any[]            # the flat form does not tile into the rows
    tables = 0
    cells = 0
    for (sys, _, _) in SWEEP
        top = first(levels(sys))
        for (rootlevel, depth) in ((top, 0), (top, 1), (top, 2),
                                   (top + 1, 2), (top + 2, 2))
            l = rootlevel + depth
            l <= last(levels(sys)) || continue
            roots = levelgrid(sys, rootlevel)
            grid = levelgrid(sys, l)
            n = ncells(roots)
            # Exhaustive at the shallowest level, where the structure lives;
            # a fixed spread deeper, where the roots are only more of the same.
            step = rootlevel == top ? 1 : max(1, n ÷ 24)
            for i in 1:step:n, conn in (Vertex(), Edge())
                root = cellindex(roots, i)
                pg = subtree(sys, root, l)
                halo_ids = collect(halo(pg; connectivity = conn, cells = true))
                t = adjacency(pg; halo = 1, connectivity = conn)
                ids = [cellindex(pg, p) for p in 1:ncells(pg)]
                want = buffer_rows(sys, l, ids, halo_ids, conn)
                tables += 1
                cells += ncells(pg)

                any(row -> any(iszero, row), want) &&
                    push!(uncovered, (syslabel(sys), root, l, conn))
                (length(t) == ncells(pg) &&
                 all(p -> collect(t[p]) == want[p], 1:ncells(pg))) ||
                    push!(disagreed, (syslabel(sys), root, l, conn))
                (t.offsets[1] == 1 && t.offsets[end] == length(t.indices) + 1 &&
                 all(s -> 1 <= s <= length(t) + length(halopositions(t)), t.indices) &&
                 reduce(vcat, t; init = Int[]) == t.indices) ||
                    push!(csr, (syslabel(sys), root, l, conn))

                # The link between "row i" and "buffer slot i" that nothing else
                # states: the chunk half is read as `descendant_range`, so row i
                # must be that range's i-th cell. Both the table and the oracle
                # above go through `cellindex(pg, ·)`, so they would agree with
                # each other even if the range disagreed with both — and the
                # caller's buffer would then be shifted under a correct table.
                if has_sorted_subtrees(sys)
                    ids == [cellindex(grid, q) for q in descendant_range(sys, root, l)] ||
                        push!(layout, (syslabel(sys), root, l))
                end
            end
        end
    end
    @test isempty(uncovered)
    @test isempty(disagreed)
    @test isempty(layout)
    @test isempty(csr)
    @test tables == 1846
    @test cells == 38828
end

@testset "a hole is addressed, and both paths agree" begin
    holes = Any[]
    hole_in_halo = 0
    for (sys, base, leaf) in SWEEP
        grid = levelgrid(sys, leaf)
        root = cellindex(levelgrid(sys, base), 3)
        ids = descendants(sys, root, leaf)
        n = length(ids)
        punched = Set(ids[(n÷3+1):(2n÷3)])
        sub = holed(sys, base, leaf)
        for conn in (Vertex(), Edge())
            halo_ids = collect(halo(sub; connectivity = conn, cells = true))
            t = adjacency(sub; halo = 1, connectivity = conn)
            want = buffer_rows(sys, leaf, [cellindex(sub, p) for p in 1:ncells(sub)],
                halo_ids, conn)
            all(p -> collect(t[p]) == want[p], 1:ncells(sub)) ||
                push!(holes, (syslabel(sys), conn))
            any(x -> x in punched, halo_ids) && (hole_in_halo += 1)
        end
    end
    @test isempty(holes)
    @test hole_in_halo == 2 * length(SWEEP)

    # The two membership routes are one function. A sorted-subtree system can
    # build the same subtree rooted (contiguous block, integer range) or bare
    # (membership search), and the gate tells them apart — so this is not two
    # spellings of one code path.
    paths = Any[]
    for (sys, base, leaf) in SWEEP
        has_sorted_subtrees(sys) || continue
        grid = levelgrid(sys, leaf)
        root = cellindex(levelgrid(sys, base), 3)
        blocked = subtree(sys, root, leaf)
        member = PartialGrid(sys, leaf, descendants(sys, root, leaf))
        (EN._whole_subtree_range(blocked) !== nothing &&
         EN._whole_subtree_range(member) === nothing) || push!(paths, (syslabel(sys), :gate))
        for conn in (Vertex(), Edge())
            halo_pos = collect(halo(blocked; connectivity = conn))
            adjacency(blocked, halo_pos; connectivity = conn) ==
            adjacency(member, halo_pos; connectivity = conn) ||
                push!(paths, (syslabel(sys), conn))
        end
    end
    @test isempty(paths)

    bridged = Any[]
    for (sys, base, leaf) in SWEEP
        grid = levelgrid(sys, leaf)
        sub = rooted(sys, base, leaf)
        cv = CellVector(sub)
        wrapped = PartialGrid(system(cv), level(cv), cv)
        halo_pos = collect(halo(sub))
        adjacency(sub, halo_pos) == adjacency(wrapped, halo_pos) ||
            push!(bridged, syslabel(sys))
    end
    @test isempty(bridged)
end


@testset "adjacency refuses what would come back short" begin
    sys = DGG.HEALPixSystem()
    leaf = 3
    pg = subtree(sys, cellindex(levelgrid(sys, 1), 3), leaf)
    vertex_halo = collect(halo(pg))
    edge_halo = collect(halo(pg; connectivity = Edge()))
    @test length(vertex_halo) == 19 && length(edge_halo) == 16

    # A row exists only for an in-region position, so a wider receptive field is
    # a wider REGION and the message has to say so.
    @test_throws ArgumentError adjacency(pg; halo = 2)
    @test occursin("grow(region, 1)",
        try
            adjacency(pg; halo = 2)
        catch e
            e.msg
        end)
    @test_throws ArgumentError adjacency(pg; halo = :edge)
    # The lazy halo is not a buffer, and saying so is a method rather than a
    # `MethodError`, because passing it is the mistake the lazy design invites.
    @test_throws ArgumentError adjacency(pg, halo(pg))
    # An unsorted fetch list would misaddress rows rather than fail.
    @test_throws ArgumentError adjacency(pg, reverse(vertex_halo))
    @test all(j -> (try
            adjacency(pg, deleteat!(copy(vertex_halo), j))
            false
        catch e
            e isa ArgumentError
        end), eachindex(vertex_halo))
    # The connectivity has to be the halo's: an `Edge()` margin cannot complete
    # `Vertex()` rows, and it does complete `Edge()` ones.
    @test_throws ArgumentError adjacency(pg, edge_halo; connectivity = Vertex())
    @test adjacency(pg, edge_halo; connectivity = Edge()) isa AdjacencyTable

    # `offsets` is one longer than the table, so an unchecked position past the
    # end reads a real offset and returns a plausible row of somebody else's
    # neighbours.
    t = adjacency(pg, vertex_halo)
    @test_throws BoundsError t[0]
    @test_throws BoundsError t[length(t)+1]
end


# ---------------------------------------------------------------------------
# The three guarantees a caller writes code against
# ---------------------------------------------------------------------------

# THE ANCHOR. Every row is a window onto the ring the COMPLETE level answers.
# The complete-width shapes keep slot indices, so slot k names ring member k in
# any table of any region holding the cell; the clipped shape keeps order and
# drops, which is a different row and deliberately not slot-addressable. A
# re-anchoring per region, or a `sort!`, or a clipped row that padded instead of
# dropping, each breaks exactly one of the three comparisons below.
@testset "every row is a window onto the canonical ring: $(syslabel(sys))" for (sys, base, leaf) in SWEEP
    complete = levelgrid(sys, leaf)
    for (label, sub) in shapes(sys, base, leaf), conn in (Vertex(), Edge())
        n = ncells(sub)
        clipped = adjacency(sub; connectivity = conn)
        marked = adjacency(sub; halo = :mark, connectivity = conn)
        buffered = adjacency(sub; halo = 1, connectivity = conn)
        hpos = halopositions(buffered)
        bad = Int[]
        for p in 1:n
            c = cellindex(sub, p)
            canonical = collect(neighbors(complete, c, 1; connectivity = conn))
            slots = [something(cellposition(sub, nb), 0) for nb in canonical]
            buf = [q == 0 ?
                   n + searchsortedfirst(hpos, cellposition(complete, canonical[j])) : q
                   for (j, q) in enumerate(slots)]
            (collect(marked[p]) == slots &&
             collect(buffered[p]) == buf &&
             collect(clipped[p]) == filter(!iszero, slots)) || push!(bad, p)
        end
        @test isempty(bad)
    end
end

# THE SENTINEL. `:mark` is the slot geometry with NOTHING fetched: `0` is the
# only value outside `1:length(t)`, no halo is walked, and a region that is a
# complete level has no zeros to mark. A `:mark` table that quietly walked a
# halo would still index correctly and would fail here.
@testset "the :mark sentinel names the outside without walking it" begin
    for (sys, base, leaf) in SWEEP
        sub = holed(sys, base, leaf)
        t = adjacency(sub; halo = :mark)
        @test isempty(halopositions(t))
        @test isempty(halocells(t))
        @test all(q -> 0 <= q <= length(t), t.indices)
        @test any(iszero, t.indices)
        # Row width is the complete level's degree, border cells included.
        complete = levelgrid(sys, leaf)
        @test all(p -> length(t[p]) ==
                       length(neighbors(complete, cellindex(sub, p), 1)), 1:length(t))
        # A complete level is all interior, so nothing is marked and the marked
        # table is the clipped one.
        shallow = levelgrid(sys, first(levels(sys)))
        whole = adjacency(shallow; halo = :mark)
        @test !any(iszero, whole.indices)
        @test whole == adjacency(shallow)
    end
end

# THE BUFFER ORDER. Slot `length(t) + j` is the j-th element of `halo(region)`
# under the table's connectivity — the one fact that lets a caller fetch the
# halo separately and index the second half without a lookup table. A build that
# walked the halo in a different order, or numbered it from a `Set`, breaks the
# first comparison; one that renumbered per row breaks the second.
@testset "slot n+j is the j-th halo cell: $(syslabel(sys))" for (sys, base, leaf) in SWEEP
    for (label, sub) in shapes(sys, base, leaf), conn in (Vertex(), Edge())
        t = adjacency(sub; halo = 1, connectivity = conn)
        @test halopositions(t) == collect(halo(sub; connectivity = conn))
        @test halocells(t) == collect(halo(sub; connectivity = conn, cells = true))
        # The caller's own list, passed back in, is the same table.
        @test adjacency(sub, halopositions(t); connectivity = conn) == t
        # And every out-of-region slot resolves to the cell it names.
        n = length(t)
        complete = levelgrid(sys, leaf)
        cells = halocells(t)
        @test all(1:n) do p
            row = collect(neighbors(complete, cellindex(sub, p), 1; connectivity = conn))
            all(zip(t[p], row)) do (q, nb)
                q <= n ? cellindex(sub, q) == nb : cells[q-n] == nb
            end
        end
    end
end

@testset "the buffered table allocates nothing sized by the halo list" begin
    sys = DGG.HEALPixSystem()
    leaf = 5
    grid = levelgrid(sys, leaf)
    root = cellindex(levelgrid(sys, 1), 3)
    pg = subtree(sys, root, leaf)
    block = descendant_range(sys, root, leaf)
    tight = collect(halo(pg))
    wide = [q for q in 1:ncells(grid) if !(first(block) <= q <= last(block))]
    @test length(tight) == 67 && length(wide) == 12032

    a = adjacency(pg, tight; threaded = false)
    b = adjacency(pg, wide; threaded = false)
    @test length(a.indices) == length(b.indices)
    # ...and the two are NOT the same table: the halo half is numbered against
    # the list it was given, so this is two real answers, not one memoised.
    @test a != b

    # The rows themselves cost the same either way — the buffer's second half
    # is addressed by search, never materialised per row. The only thing the
    # public verb allocates in proportion to the list is the copy the table
    # keeps for `halopositions`.
    cv = CellVector(pg)
    n = ncells(pg)
    rows(list) = EN._adjacency_rows(cv, Vertex(), EN.BufferedRows(list, n), false)
    rowbytes(list) = (rows(list); @allocated rows(list))
    @test rowbytes(tight) == rowbytes(wide)
end


const FIXTURE = joinpath(@__DIR__, "..", "..", "fixtures", "california.txt")

# `P` opens a polygon, `R` opens a ring, everything else is `lon lat`. A copy of
# the loader in `multiorder_polygons.jl` rather than an import, by the same
# convention its neighbours keep: these files must be able to disagree.
function load_parts(path)
    parts = Vector{Vector{Vector{Tuple{Float64,Float64}}}}()
    ring = Tuple{Float64,Float64}[]
    for line in eachline(path)
        (isempty(line) || startswith(line, '#')) && continue
        if line == "P"
            push!(parts, Vector{Vector{Tuple{Float64,Float64}}}())
        elseif line == "R"
            ring = Tuple{Float64,Float64}[]
            push!(parts[end], ring)
        else
            lon, lat = split(line)
            push!(ring, (parse(Float64, lon), parse(Float64, lat)))
        end
    end
    return parts
end

const MAINLAND = GI.Polygon([GI.LinearRing(r) for r in load_parts(FIXTURE)[1]])
const DONUT = GI.Polygon([
    GI.LinearRing([(-122.0, 35.0), (-118.0, 35.0), (-118.0, 39.0), (-122.0, 39.0),
        (-122.0, 35.0)]),
    GI.LinearRing([(-121.0, 36.0), (-121.0, 38.0), (-119.0, 38.0), (-119.0, 36.0),
        (-121.0, 36.0)])])

# How many times a sequence of members steps backwards in azimuth about `c`,
# each read on its own level's grid because a coarse member has no centroid at
# the reference level. One for a single counter-clockwise turn, `n - 1` for a
# clockwise one, arbitrary for an id-sorted sequence. Written out here rather
# than taken from the conformance harness, whose `winding_problems` reads one
# grid and so cannot answer about a mixed-level set.
function azimuth_wraps(sys, c, members)
    centroid(x) = DGG.cell_centroid(levelgrid(sys, level(x)), x)
    p = centroid(c)
    ax = abs(p[1]) <= abs(p[2]) ?
         (abs(p[1]) <= abs(p[3]) ? (1.0, 0.0, 0.0) : (0.0, 0.0, 1.0)) :
         (abs(p[2]) <= abs(p[3]) ? (0.0, 1.0, 0.0) : (0.0, 0.0, 1.0))
    s = ax[1] * p[1] + ax[2] * p[2] + ax[3] * p[3]
    t = (ax[1] - s * p[1], ax[2] - s * p[2], ax[3] - s * p[3])
    n = sqrt(t[1]^2 + t[2]^2 + t[3]^2)
    e1 = (t[1] / n, t[2] / n, t[3] / n)
    e2 = (p[2] * e1[3] - p[3] * e1[2], p[3] * e1[1] - p[1] * e1[3],
        p[1] * e1[2] - p[2] * e1[1])
    az = map(members) do m
        q = centroid(m)
        d = (q[1] - p[1], q[2] - p[2], q[3] - p[3])
        r = d[1] * p[1] + d[2] * p[2] + d[3] * p[3]
        u = (d[1] - r * p[1], d[2] - r * p[2], d[3] - r * p[3])
        atan(u[1] * e2[1] + u[2] * e2[2] + u[3] * e2[3],
            u[1] * e1[1] + u[2] * e1[2] + u[3] * e1[3])
    end
    return count(i -> az[mod1(i + 1, length(az))] < az[i] - 1e-9, eachindex(az))
end

function oracle(set, c; connectivity = Vertex(), geometric = false)
    sys = system(set)
    L = set.reference_level
    grid = levelgrid(sys, L)
    tree = geometric ? DGG.treeify(grid) : nothing
    owner = Dict{cellindextype(sys),cellindextype(sys)}()
    for m in set, d in descendants(sys, m, L)
        owner[d] = m
    end
    adjacent(d) = geometric ? FB.adjacent_cells(grid, d, connectivity, tree) :
                  neighbors(grid, d, 1; connectivity)
    out = Set{cellindextype(sys)}()
    for d in descendants(sys, c, L), nb in adjacent(d)
        m = get(owner, nb, nothing)
        (m === nothing || m == c) && continue
        push!(out, m)
    end
    return out
end

# The members to ask about: the coarsest few, which are the ones with a subtree
# and therefore with a border, plus a spread over the rest.
function member_probes(set, n::Int)
    shallowest = minimum(level, set)
    coarse = [i for i in eachindex(set) if level(set[i]) == shallowest]
    step = max(1, length(set) ÷ n)
    return unique(vcat(first(coarse, 3), collect(1:step:length(set))))
end

# The quad-face family is the family whose four children tile their parent
# exactly, so that a member's footprint IS its descendants' union and a
# statement about level-`L` cells is a statement about the drawn polygons.
# IGEO7, H3 and A5 refine non-congruently and are excluded from the GEOMETRIC
# law with that reason — see `member_neighbors`' docstring; they keep every
# other law here.

# The last column is whether `Vertex()` and `Edge()` are different relations at
# all: IGEO7 and H3 have only 3-valent vertices, so on those two they coincide
# and "Edge drops the diagonal contacts" has nothing to drop. Stated rather than
# left to a passing-by-vacuity test.
const MOC_SWEEP = [
    (DGG.IGeo7System(), 6, DONUT, false),
    (DGG.H3System(), 5, DONUT, false),
    (DGG.HEALPixSystem(), 8, MAINLAND, true),
    (DGG.A5System(), 8, DONUT, true),
    (DGG.S2System(), 8, MAINLAND, true),
    (DGG.ISEA4RSystem(), 8, DONUT, true),
    (DGG.AuthalicSystem(DGG.IGeo7System()), 6, DONUT, false),
]

@testset "member_neighbors: $(syslabel(sys))" for (sys, lvl, target, splits) in MOC_SWEEP
    set = query(sys, MultiOrderCoverage(target); level = lvl)

    @testset "the set is genuinely mixed-level" begin
        @test length(set) > 20
        @test minimum(level, set) < maximum(level, set)
    end

    @testset "exact against the expanding oracle" begin
        for i in member_probes(set, 6), conn in (Vertex(), Edge())
            got = member_neighbors(set, set[i]; connectivity = conn)
            @test Set(got) == oracle(set, set[i]; connectivity = conn)
            @test length(got) == length(unique(got))
            @test !(set[i] in got)
        end
    end

    # The order is the package's one order on a mixed-level set: a single
    # counter-clockwise turn about the member, each neighbour read on its own
    # level's grid. `issorted` by `(level, position)` — what this used to be —
    # is exactly the ascending order the contract now refuses.
    @testset "one counter-clockwise turn about the member" begin
        for i in member_probes(set, 6), conn in (Vertex(), Edge())
            c = set[i]
            got = member_neighbors(set, c; connectivity = conn)
            @test got == member_neighbors(set, c; connectivity = conn)
            length(got) < 3 && continue
            @test azimuth_wraps(sys, c, got) == 1
        end
    end

    # `Edge()` is the restriction, so it can only ever drop members — and where
    # a vertex carries more than three cells it does, at every corner contact.
    @testset "Edge is Vertex restricted" begin
        strictly_smaller = 0
        for i in member_probes(set, 6)
            v = member_neighbors(set, set[i])
            e = member_neighbors(set, set[i]; connectivity = Edge())
            @test Set(e) ⊆ Set(v)
            length(e) < length(v) && (strictly_smaller += 1)
        end
        # The exclusion is the systems whose vertices are all 3-valent, where
        # the two connectivities ARE the same relation — not a weaker law.
        splits && @test strictly_smaller > 0
        splits || @test strictly_smaller == 0
    end

    @testset "a non-member has no neighbours here" begin
        outside = cellindex(levelgrid(sys, set.reference_level), 1)
        @test_throws ArgumentError member_neighbors(set, outside)
    end

    if isquadface(sys)
        # The geometric law, and the reason the exclusions above are stated
        # rather than silently taken: `adjacent_cells` reads the boundary RINGS
        # and counts shared vertices, so this is "share a boundary" and "share
        # an edge" in the plain sense, on the systems where a member's footprint
        # is its descendants' union.
        @testset "geometric boundary sharing" begin
            for i in member_probes(set, 4), conn in (Vertex(), Edge())
                @test Set(member_neighbors(set, set[i]; connectivity = conn)) ==
                      oracle(set, set[i]; connectivity = conn, geometric = true)
            end
        end
    end
end


const CAP = GO.UnitSpherical.SphericalCap(
    GO.UnitSpherical.UnitSphereFromGeographic()((10.0, 46.0)), deg2rad(1.0))

coarsest(set) = set[argmin(i -> level(set[i]), eachindex(set))]

function walk_bytes(set)
    c = coarsest(set)
    member_neighbors(set, c)                    # warm up, then measure
    return @allocated member_neighbors(set, c)
end

@testset "the walk is the border, not the subtree" begin
    for (sys, shallow) in ((DGG.HEALPixSystem(), 9), (DGG.S2System(), 9),
        (DGG.ISEA4RSystem(), 9), (DGG.IGeo7System(), 7))
        near = query(sys, MultiOrderCoverage(CAP); level = shallow)
        far = query(sys, MultiOrderCoverage(CAP); level = shallow + 2)
        grow = length(descendant_range(sys, coarsest(far), far.reference_level)) ÷
               length(descendant_range(sys, coarsest(near), near.reference_level))
        @test grow >= 16                        # or the measurement is not one
        @test walk_bytes(far) <= 4 * walk_bytes(near)
    end
end

# The cross-idiom laws compare the subset idioms with each other, which a
# consistent re-anchoring would survive. This compares them with the complete
# level's canonical ring: a clipped ring is that ring with non-members dropped
# in place, rotation pinned.
@testset "clipped rings are pinned-rotation subsequences: $(syslabel(sys))" for (sys, base, leaf) in SWEEP
    pg = subtree(sys, cellindex(levelgrid(sys, base), 3), leaf)
    complete = levelgrid(sys, leaf)
    c = first(c for c in (cellindex(pg, i) for i in 1:ncells(pg))
              if length(neighbors(pg, c, 1)) < length(neighbors(complete, c, 1)))
    for k in 1:2
        @test neighbors(pg, c, k) == filter(in(pg), neighbors(complete, c, k))
    end
    @test adjacency(pg)[cellposition(pg, c)] ==
          [cellposition(pg, x) for x in filter(in(pg), neighbors(complete, c, 1))]
end

end # module StencilTests
