# ---------------------------------------------------------------------------
# T21 — stencils on subsets.
#
# `neighbors` / `ring` on a subset of a level are the COMPLETE level's answer
# clipped to membership. That is a decision, not a discovery: the generic
# geometric path used to walk the subset's own adjacency graph, which agrees at
# `k == 1` and disagrees from `k == 2` wherever the subset has a hole or a
# concave boundary. This file is what makes the chosen reading executable.
#
# The laws, in order:
#
#   * CLIPPING — `ring(sub, c, k) == filter(in(sub), ring(complete, c, k))`,
#     element for element, so the clip preserves the rotational order as well
#     as the membership. Three subset SHAPES (a rooted subtree, a scattered
#     non-rooted id set, a subtree with a hole punched in it), `k` in `0:3`,
#     both connectivities, every system. Asserting the ORDER is what pins that
#     the fast path routes through the system's automaton rather than through
#     the geometric tree walk, which orders by a frame of its own.
#   * DIVERGENCE — the k >= 2 disagreement with breadth-first distance INSIDE
#     the subset, demonstrated on the hole rather than described. Stated as an
#     assertion so that a future return to induced-subgraph semantics fails
#     here instead of silently changing answers.
#   * OUT OF SET — a cell the subset does not hold is an `ArgumentError` on
#     both container types, never the complete level's answer.
#   * FACES — `PartialGrid`, `CellVector` and `CellLookup` answer identically,
#     and the position face is the id face resolved through `cellposition`.
#   * HALO — `halo_table(sub, k)[p] == neighbors(sub, p, k)` for every `p`,
#     including on the rooted-subtree fast path, which takes a different route
#     (T20's rim/interior split) to the same answer.
#   * MEMBER ADJACENCY — `member_neighbors` against a brute-force oracle that
#     expands every member and shares no code with the rim walk: first with the
#     level grid's own adjacency, then — on the systems whose children tile
#     their parent — with the GEOMETRIC shared-vertex/shared-edge test on the
#     boundary rings, which is what "share a boundary" means.
#
# Written against the generic interface only: no system module is imported, and
# every law runs against every system in `systems()` plus an
# `AuthalicSystem`-wrapped one.
# ---------------------------------------------------------------------------

module StencilTests

using Test
import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import GeometryOps as GO

const FB = DGG.Fallbacks

using DiscreteGlobalGrids: systems, levelgrid, ncells, cellindex, cellposition,
    neighbors, ring, halo_table, level, levels, max_level, descendants,
    descendant_range, has_sorted_subtrees, PartialGrid, CellVector, CellLookup,
    MultiOrderCoverage, member_neighbors, AuthalicSystem, Vertex, Edge,
    Connectivity, cellindextype, query, system

# ---------------------------------------------------------------------------
# Systems, and the depths each is swept at
#
# A root two levels above the leaf gives a subtree big enough to have an
# interior and a rim, and small enough that the brute-force oracles below can
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

sysname(sys) = sys isa AuthalicSystem ?
               "Authalic($(nameof(typeof(parent(sys)))))" : string(nameof(typeof(sys)))

@testset "the sweep covers every registered system" begin
    swept = Set(typeof(s) for (s, _, _) in SWEEP)
    for s in systems()
        @test typeof(s) in swept
    end
    @test any(s -> s isa AuthalicSystem, first.(SWEEP))
end

# ---------------------------------------------------------------------------
# The subset shapes
# ---------------------------------------------------------------------------

held(sub, c) = cellposition(sub, c) !== nothing

# A rooted subtree, the shape `PartialGrid(sys, cell, level)` builds.
rooted(sys, base, leaf) =
    PartialGrid(sys, cellindex(levelgrid(sys, base), 3), leaf)

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

# Breadth-first distance INSIDE the subset — the reading this task replaced.
# Kept here as the oracle for DIVERGENCE, and for nothing else.
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

@testset "$(sysname(sys))" for (sys, base, leaf) in SWEEP
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
    end

    @testset "$label: the position face is the id face" for (label, sub) in
                                                            shapes(sys, base, leaf)
        for p in 1:max(1, ncells(sub) ÷ 6):ncells(sub), conn in (Vertex(), Edge()), k in 0:2
            c = cellindex(sub, p)
            @test neighbors(sub, p, k; connectivity = conn) ==
                  sort([cellposition(sub, x)
                        for x in neighbors(sub, c, k; connectivity = conn)])
            @test ring(sub, p, k; connectivity = conn) ==
                  sort([cellposition(sub, x)
                        for x in ring(sub, c, k; connectivity = conn)])
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

    @testset "$label: halo table is the per-cell calls" for (label, sub) in
                                                            shapes(sys, base, leaf)
        for conn in (Vertex(), Edge()), k in (1, 2)
            table = halo_table(sub, k; connectivity = conn)
            @test length(table) == ncells(sub)
            @test all(p -> table[p] == neighbors(sub, p, k; connectivity = conn),
                1:ncells(sub))
        end
        cv = CellVector(sub)
        @test halo_table(cv) == halo_table(sub)
        @test halo_table(CellLookup(cv)) == halo_table(sub)
    end

    # The rooted subtree is the one shape with a fast path of its own — T20's
    # rim/interior split, which never tests membership on an interior cell. It
    # has to reach the same table the generic route does, so both are asked.
    @testset "the rooted fast path agrees with the generic route" begin
        sub = rooted(sys, base, leaf)
        for conn in (Vertex(), Edge())
            @test halo_table(sub, 1; connectivity = conn) ==
                  [neighbors(sub, p, 1; connectivity = conn) for p in 1:ncells(sub)]
        end
        # The split is only claimed for `k == 1`: a cell one step inside the rim
        # still reaches outside at `k == 2`, so the deeper table must NOT take
        # it. Same answer either way, which is what is asserted.
        @test halo_table(sub, 2) == [neighbors(sub, p, 2) for p in 1:ncells(sub)]
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
    end
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

# ---------------------------------------------------------------------------
# MEMBER ADJACENCY on a multi-order set
#
# The fixture is the committed California outline the T15/T18 files use, plus
# their synthetic donut — a box with a square hole, which is the shape that
# produces genuinely mixed levels around a concave boundary.
# ---------------------------------------------------------------------------

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

# The oracle: expand EVERY member to the reference level, own each leaf, then
# walk `c`'s whole subtree — not its rim — and collect the members its leaves'
# neighbours belong to. It materialises what the verb refuses to, and shares no
# arithmetic with the rim walk.
#
# `adjacency` is the level's own relation by default and the GEOMETRIC one when
# asked, so the same oracle states both laws.
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
# and therefore with a rim, plus a spread over the rest.
function member_probes(set, n::Int)
    shallowest = minimum(level, set)
    coarse = [i for i in eachindex(set) if level(set[i]) == shallowest]
    step = max(1, length(set) ÷ n)
    return unique(vcat(first(coarse, 3), collect(1:step:length(set))))
end

# The three systems whose four children tile their parent exactly, so that a
# member's footprint IS its descendants' union and a statement about level-`L`
# cells is a statement about the drawn polygons. IGEO7, H3 and A5 refine
# non-congruently and are excluded from the GEOMETRIC law with that reason —
# see `member_neighbors`' docstring; they keep every other law here.
const CONGRUENT = (DGG.HEALPixSystem, DGG.S2System, DGG.ISEA4RSystem)

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

@testset "member_neighbors: $(sysname(sys))" for (sys, lvl, target, splits) in MOC_SWEEP
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

    @testset "ascending (level, position)" begin
        for i in member_probes(set, 6)
            got = member_neighbors(set, set[i])
            key(m) = (level(m), cellposition(levelgrid(sys, level(m)), m))
            @test issorted(got; by = key)
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

    if any(T -> sys isa T, CONGRUENT)
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

# ---------------------------------------------------------------------------
# The rim walk is a rim walk
#
# The claim `member_neighbors` makes is that a member far above the reference
# level costs its subtree's PERIMETER and not its area. Measured as the growth
# between two coverages of the same cap two levels apart: the coarsest member
# stays where it is, its subtree grows by the aperture squared — sixteenfold on
# the quadtrees, forty-ninefold on IGEO7 — and the walk must not notice. An
# implementation that expanded the subtree would grow with it, so the ratio is
# the assertion and the absolute byte counts are not.
#
# A spherical cap rather than an outline: a short smooth boundary keeps the
# traversal cheap while still leaving a deep level spread between the coarsest
# member and the reference level, which is the whole shape being measured.
# ---------------------------------------------------------------------------

const CAP = GO.UnitSpherical.SphericalCap(
    GO.UnitSpherical.UnitSphereFromGeographic()((10.0, 46.0)), deg2rad(1.0))

coarsest(set) = set[argmin(i -> level(set[i]), eachindex(set))]

function walk_bytes(set)
    c = coarsest(set)
    member_neighbors(set, c)                    # warm up, then measure
    return @allocated member_neighbors(set, c)
end

@testset "the walk is the rim, not the subtree" begin
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

end # module StencilTests
