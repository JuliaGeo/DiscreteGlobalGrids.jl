# The harness is tested against a conforming cube-face quadtree and variants
# that each violate one contract. The base mock implements the required
# interface; an `Extras` variant adds optional methods so guarded checks run.

module TestConformance

using Test
using DiscreteGlobalGrids
using DiscreteGlobalGridsConformanceTesting
import DiscreteGlobalGrids as DGG
import DiscreteGlobalGridsConformanceTesting as Conf
import GeometryOps as GO

const USph = GO.UnitSpherical
# ===========================================================================
# Cube-face quadtree mock
#
# Six root faces of the circumscribed cube, each subdivided 2x2 in face
# coordinates and radially projected onto the sphere. Cell ids are
# `LevelIndex(level, idx)` with `idx` 0-based and face-major, Morton-ordered
# within a face, which makes the hierarchy pure arithmetic: the children of
# `idx` are `4idx .+ (0:3)` and the parent is `idx ÷ 4`. Subtrees are therefore
# contiguous in the dense order, so `has_sorted_subtrees` is `true` and
# `descendant_range` is closed-form.
#
# `Extras` gates the optional methods; `bug` injects exactly one contract
# violation and is `:none` for the correct mock.
# ===========================================================================

struct CubeSystem{Extras} <: AbstractHierarchicalGridSystem
    maxlevel::Int
    bug::Symbol
end

CubeSystem(; maxlevel::Int = 3, bug::Symbol = :none, extras::Bool = false) =
    CubeSystem{extras}(maxlevel, bug)

struct CubeGrid{S<:CubeSystem} <: AbstractGrid
    sys::S
    level::Int
end

# --- geometry --------------------------------------------------------------

# (u, v) -> 3D on face `f`. The bases are ordered so that e_u x e_v is the
# OUTWARD normal, which is what makes the counter-clockwise corner order in
# (u, v) come out counter-clockwise seen from outside the sphere.
_face_point(f::Int, u::Real, v::Real) =
    f == 0 ? (1.0, float(u), float(v)) :      # +X
    f == 1 ? (-1.0, float(v), float(u)) :     # -X
    f == 2 ? (float(v), 1.0, float(u)) :      # +Y
    f == 3 ? (float(u), -1.0, float(v)) :     # -Y
    f == 4 ? (float(u), float(v), 1.0) :      # +Z
             (float(v), float(u), -1.0)       # -Z

function _unit(t)
    n = sqrt(t[1]^2 + t[2]^2 + t[3]^2)
    return GO.UnitSphericalPoint(t[1] / n, t[2] / n, t[3] / n)
end

"Face, and integer cell coordinates within it, of 0-based `idx` at `level`."
function _decode(level::Int, idx::Int)
    f, within = divrem(idx, 4^level)
    x = 0
    y = 0
    for b in 0:(level - 1)
        x |= ((within >> (2b)) & 1) << b
        y |= ((within >> (2b + 1)) & 1) << b
    end
    return f, x, y
end

"Inverse of [`_decode`](@ref): 0-based `idx` of cell `(x, y)` on face `f`."
function _encode(level::Int, f::Int, x::Int, y::Int)
    within = 0
    for b in 0:(level - 1)
        within |= ((x >> b) & 1) << (2b)
        within |= ((y >> b) & 1) << (2b + 1)
    end
    return f * 4^level + within
end

"""
The (u, v) box of a cell: `(f, u0, v0, u1, v1)`.

`w * x` rather than `2x / n` keeps the arithmetic exact for a power-of-two `n`,
which matters because the neighbour relation below identifies shared cube-edge
vertices by exact equality across two different faces' formulas.

With `bug = :deep_drift`, cells at level 2 and deeper move outside their
ancestors while retaining internally consistent cell geometry. With
`bug = :shrunken_cells`, every cell is pulled 1% off its own corners: each cell
is a well-formed ring around its own centroid, and no two cells share a vertex.
"""
function _cell_box(level::Int, idx::Int, bug::Symbol = :none)
    f, x, y = _decode(level, idx)
    w = 2.0 / (1 << level)
    drift = (bug === :deep_drift && level >= 2) ? 2w : 0.0
    pad = bug === :shrunken_cells ? 0.005w : 0.0
    u0 = -1.0 + w * x + drift
    v0 = -1.0 + w * y + drift
    return f, u0 + pad, v0 + pad, u0 + w - pad, v0 + w - pad
end

function _corners(level::Int, idx::Int, bug::Symbol = :none)
    f, u0, v0, u1, v1 = _cell_box(level, idx, bug)
    return [_unit(_face_point(f, u0, v0)), _unit(_face_point(f, u1, v0)),
            _unit(_face_point(f, u1, v1)), _unit(_face_point(f, u0, v1))]
end

function _centroid(level::Int, idx::Int, bug::Symbol = :none)
    f, u0, v0, u1, v1 = _cell_box(level, idx, bug)
    return _unit(_face_point(f, (u0 + u1) / 2, (v0 + v1) / 2))
end

# --- required base grid interface ------------------------------------------

DGG.ncells(g::CubeGrid) = 6 * 4^g.level
DGG.system(g::CubeGrid) = g.sys
DGG.level(g::CubeGrid) = g.level

function DGG.cellindex(g::CubeGrid, i::Int)
    1 <= i <= DGG.ncells(g) || throw(BoundsError(g, i))
    return LevelIndex(g.level, i - 1)
end

function DGG.cellposition(g::CubeGrid, c::LevelIndex)
    level(c) == g.level || return nothing
    idx = rawid(c)
    0 <= idx < DGG.ncells(g) || return nothing
    return g.sys.bug === :bad_position ? Int(idx) + 2 : Int(idx) + 1
end

function DGG.cell_boundary(g::CubeGrid, c::LevelIndex)
    pts = _corners(g.level, Int(rawid(c)), g.sys.bug)
    return g.sys.bug === :cw_boundary ? reverse(pts) : pts
end

function DGG.cell_centroid(g::CubeGrid, c::LevelIndex)
    p = _centroid(g.level, Int(rawid(c)), g.sys.bug)
    # A centroid on the boundary is the subtle one: it is "in" the cell by every
    # closed-region test, and only the strictness clause of the contract rejects it.
    g.sys.bug === :centroid_on_boundary && return first(_corners(g.level, Int(rawid(c))))
    g.sys.bug === :centroid_outside && return GO.UnitSphericalPoint(-p[1], -p[2], -p[3])
    return p
end

# --- required hierarchical system interface --------------------------------

DGG.cellindextype(::CubeSystem) = LevelIndex
DGG.levels(s::CubeSystem) = 0:s.maxlevel
DGG.has_sorted_subtrees(::CubeSystem) = true
DGG.max_neighbors(::CubeSystem, ::Vertex) = 8
DGG.max_neighbors(::CubeSystem, ::Edge) = 4

function DGG.levelgrid(s::CubeSystem, l::Integer)
    l in DGG.levels(s) || throw(ArgumentError("level $l is outside $(DGG.levels(s))"))
    return CubeGrid(s, Int(l))
end

DGG.rootcells(::CubeSystem) = [LevelIndex(0, i) for i in 0:5]

function Base.parent(s::CubeSystem, c::LevelIndex)
    level(c) == 0 && throw(ArgumentError("$c is a root cell and has no parent"))
    return LevelIndex(level(c) - 1, rawid(c) ÷ 4)
end

function DGG.children(s::CubeSystem, c::LevelIndex)
    level(c) == s.maxlevel &&
        throw(ArgumentError("$c is at max_level $(s.maxlevel) and has no children"))
    kids = [LevelIndex(level(c) + 1, 4 * rawid(c) + k) for k in 0:3]
    return s.bug === :unsorted_children ? reverse(kids) : kids
end

# Inflate the cap from the cell centroid to its farthest corner. Each child box
# is a quadrant of its parent, so the parent's spherical quad contains its
# subtree.
function DGG.node_extent(s::CubeSystem, c::LevelIndex)
    l = level(c)
    idx = Int(rawid(c))
    centre = _centroid(l, idx, s.bug)
    tight = maximum(USph.spherical_distance(centre, p) for p in _corners(l, idx, s.bug))
    factor = s.bug === :small_extent ? 0.5 : DGG.cap_inflation(s)
    return USph.SphericalCap(centre, factor * tight)
end

function DGG.descendant_range(s::CubeSystem, c::LevelIndex, l::Integer)
    l < level(c) && throw(ArgumentError("level $l is above $c"))
    d = Int(l) - level(c)
    idx = Int(rawid(c))
    stop = (idx + 1) * 4^d
    return (idx * 4^d + 1):(s.bug === :lying_range ? stop + 1 : stop)
end

# --- neighbours ------------------------------------------------------------
#
# Adjacency is derived from shared boundary vertices. Cube-edge vertices use
# identical arithmetic on both faces and therefore compare exactly. Vertex
# neighbours share at least one corner; edge neighbours share two.

const INCIDENCE = Dict{Int,Dict{NTuple{3,Float64},Vector{Int}}}()

function _incidence(level::Int)
    return get!(INCIDENCE, level) do
        d = Dict{NTuple{3,Float64},Vector{Int}}()
        for idx in 0:(6 * 4^level - 1), p in _corners(level, idx)
            push!(get!(() -> Int[], d, (p[1], p[2], p[3])), idx)
        end
        return d
    end
end

function _adjacent(level::Int, idx::Int, conn::Connectivity)
    d = _incidence(level)
    shared = Dict{Int,Int}()
    for p in _corners(level, idx), other in d[(p[1], p[2], p[3])]
        other == idx && continue
        shared[other] = get(shared, other, 0) + 1
    end
    need = conn isa Edge ? 2 : 1
    return sort!([o for (o, n) in shared if n >= need])
end

function _ball(level::Int, idx::Int, k::Int, conn::Connectivity)
    seen = Set([idx])
    frontier = [idx]
    for _ in 1:k
        next = Int[]
        for i in frontier, j in _adjacent(level, i, conn)
            j in seen && continue
            push!(seen, j)
            push!(next, j)
        end
        frontier = next
    end
    return seen
end

# --- rotational order ------------------------------------------------------
#
# Rings run counter-clockwise as seen from outside the sphere, and
# `neighbors(c, k)` concatenates rings 1:k. Because a cube neighbourhood can
# cross several faces, the mock orders each ring by azimuth in a right-handed
# tangent frame. Every ring starts on the spoke through the first ring-1
# neighbour, with equal azimuths ordered by cell id. This implementation is
# independent of the harness's order checks.

"A right-handed tangent basis at `p`, with `east × north == p` (outward)."
function _frame(p)
    ax = abs(p[1]) <= abs(p[2]) ?
        (abs(p[1]) <= abs(p[3]) ? (1.0, 0.0, 0.0) : (0.0, 0.0, 1.0)) :
        (abs(p[2]) <= abs(p[3]) ? (0.0, 1.0, 0.0) : (0.0, 0.0, 1.0))
    s = ax[1] * p[1] + ax[2] * p[2] + ax[3] * p[3]
    e = (ax[1] - s * p[1], ax[2] - s * p[2], ax[3] - s * p[3])
    n = sqrt(e[1]^2 + e[2]^2 + e[3]^2)
    east = (e[1] / n, e[2] / n, e[3] / n)
    north = (p[2] * east[3] - p[3] * east[2],
             p[3] * east[1] - p[1] * east[3],
             p[1] * east[2] - p[2] * east[1])
    return east, north
end

"Azimuth of `q` about `p`, increasing counter-clockwise seen from outside."
function _azimuth(p, frame, q)
    east, north = frame
    r = q[1] * p[1] + q[2] * p[2] + q[3] * p[3]
    d = (q[1] - r * p[1], q[2] - r * p[2], q[3] - r * p[3])
    return atan(d[1] * north[1] + d[2] * north[2] + d[3] * north[3],
                d[1] * east[1] + d[2] * east[2] + d[3] * east[3])
end

"The ids at adjacency distance exactly `k` from `idx`, counter-clockwise."
function _ring_ids(level::Int, idx::Int, k::Int, conn::Connectivity)
    k == 0 && return [idx]
    shell = collect(setdiff(_ball(level, idx, k, conn), _ball(level, idx, k - 1, conn)))
    length(shell) <= 1 && return shell
    p = _centroid(level, idx)
    frame = _frame(p)
    az(j) = _azimuth(p, frame, _centroid(level, j))
    # Ring 1 starts wherever the measurement does; every outer ring starts on
    # the spoke through the first ring-1 neighbour.
    zero = 0.0
    if k > 1
        inner = _ring_ids(level, idx, 1, conn)
        isempty(inner) || (zero = az(first(inner)))
    end
    return sort!(shell; by = j -> (mod2pi(az(j) - zero), j))
end

function DGG.neighbors(g::CubeGrid, c::LevelIndex, k::Int = 1;
        connectivity::Connectivity = Vertex())
    k >= 0 || throw(ArgumentError("k must be non-negative, got $k"))
    idx = Int(rawid(c))
    k == 0 && return LevelIndex[]
    ids = Int[]
    for j in 1:k
        shell = _ring_ids(g.level, idx, j, connectivity)
        # Sorting each ring by id preserves membership but violates rotational order.
        g.sys.bug === :sorted_neighbors && sort!(shell)
        append!(ids, shell)
    end
    out = [LevelIndex(g.level, i) for i in ids]
    # Cell 0 omits one neighbour that still reports cell 0.
    g.sys.bug === :asymmetric_neighbors && idx == 0 && !isempty(out) && pop!(out)
    return out
end

function DGG.ring(g::CubeGrid, c::LevelIndex, k::Int;
        connectivity::Connectivity = Vertex())
    k >= 0 || throw(ArgumentError("k must be non-negative, got $k"))
    k == 0 && return [c]
    ids = _ring_ids(g.level, Int(rawid(c)), k, connectivity)
    g.sys.bug === :sorted_neighbors && sort!(ids)
    # Rotating the independently computed ring preserves membership and winding
    # but breaks the requirement that it be the disc's tail block.
    g.sys.bug === :rotated_ring && length(ids) > 1 && (ids = circshift(ids, 1))
    return [LevelIndex(g.level, i) for i in ids]
end

# --- optional methods, only on the `Extras` variant ------------------------

"""
Point location by dominant axis: the face is the axis of the largest-magnitude
component, ties broken by taking the first such axis in x, y, z order, which
makes the answer deterministic on a face boundary as the contract requires.
"""
function DGG.cellat(g::CubeGrid{<:CubeSystem{true}}, p::GO.UnitSphericalPoint)
    x, y, z = p[1], p[2], p[3]
    ax, ay, az = abs(x), abs(y), abs(z)
    f, u, v = if ax >= ay && ax >= az
        x > 0 ? (0, y / x, z / x) : (1, -z / x, -y / x)
    elseif ay >= az
        y > 0 ? (2, z / y, x / y) : (3, -x / y, -z / y)
    else
        z > 0 ? (4, x / z, y / z) : (5, -y / z, -x / z)
    end
    n = 1 << g.level
    cx = clamp(floor(Int, (u + 1) / 2 * n), 0, n - 1)
    cy = clamp(floor(Int, (v + 1) / 2 * n), 0, n - 1)
    return LevelIndex(g.level, _encode(g.level, f, cx, cy))
end

function DGG.ancestor(s::CubeSystem{true}, c::LevelIndex, l::Integer)
    l > level(c) && throw(ArgumentError("level $l is below $c"))
    return LevelIndex(l, rawid(c) ÷ 4^(level(c) - Int(l)))
end

function DGG.descendants(s::CubeSystem{true}, c::LevelIndex, l::Integer)
    l < level(c) && throw(ArgumentError("level $l is above $c"))
    d = Int(l) - level(c)
    idx = Int(rawid(c))
    return [LevelIndex(l, i) for i in (idx * 4^d):((idx + 1) * 4^d - 1)]
end

# ===========================================================================
# Failure capture
#
# Captured test sets verify that the public `@testset` entry points report the
# same deliberate failures as the boolean `check_*` predicates.
# ===========================================================================

function count_results(ts)
    # A `DefaultTestSet` keeps failures as objects but collapses passes into a
    # counter, so both have to be read to total a captured run.
    passes = hasproperty(ts, :n_passed) ? ts.n_passed : 0
    fails = 0
    for r in ts.results
        if r isa Test.Pass
            passes += 1
        elseif r isa Test.Fail || r isa Test.Error
            fails += 1
        elseif r isa Test.AbstractTestSet
            p, f = count_results(r)
            passes += p
            fails += f
        end
    end
    return passes, fails
end

"Run `f` under a captured, silenced test set; return `(passes, failures)`."
function capture(f)
    ts = Test.DefaultTestSet("capture")
    Test.push_testset(ts)
    try
        redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                f()
            end
        end
    catch err
        err isa Test.TestSetException || rethrow()
    finally
        Test.pop_testset()
    end
    return count_results(ts)
end

# ===========================================================================
# The suite
# ===========================================================================

const MINIMAL = CubeSystem(; maxlevel = 3)
const FULL = CubeSystem(; maxlevel = 3, extras = true)
broken(bug::Symbol) = CubeSystem(; maxlevel = 3, bug)

@testset "conformance harness" begin

    @testset "correct mock passes test_grid_interface" begin
        for l in 0:3
            test_grid_interface(DGG.levelgrid(MINIMAL, l);
                n_samples = 12, label = "CubeGrid level $l")
        end
        # `Extras` implements `cellat`, enabling the guarded lookup checks.
        test_grid_interface(DGG.levelgrid(FULL, 2);
            n_samples = 12, label = "CubeGrid+cellat level 2")
    end

    @testset "correct mock passes test_hierarchical_system" begin
        test_hierarchical_system(MINIMAL; n_samples = 6, label = "CubeSystem")
        test_hierarchical_system(FULL; n_samples = 6, label = "CubeSystem+extras")
    end

    @testset "harness catches: node_extent too small (covering law)" begin
        @test Conf.check_covering_law(MINIMAL)
        @test !Conf.check_covering_law(broken(:small_extent))
        problems = Conf.covering_law_problems(broken(:small_extent))
        @test !isempty(problems)
        @test all(contains("covering law"), problems)
    end

    # The drifting descendant is internally consistent but falls outside an
    # ancestor's extent, so only the recursive covering check detects it.
    @testset "harness catches: descendant drifting out of its ancestors' extents" begin
        drifted = broken(:deep_drift)
        @test Conf.check_grid_interface(DGG.levelgrid(drifted, 2))   # invisible here
        @test Conf.check_covering_law(MINIMAL)
        @test !Conf.check_covering_law(drifted)                       # caught only here
        problems = Conf.covering_law_problems(drifted)
        @test any(contains("ancestor at level 0, descendant at level 2"), problems)

        # With branch descent disabled, the chain to `max_level` still detects
        # deep drift; disabling the chain removes that coverage.
        @test !Conf.check_covering_law(drifted; descent_depth = 0)
        @test Conf.check_covering_law(drifted; descent_depth = 0, deep_chain = false)
    end

    @testset "harness catches: clockwise boundary ring" begin
        @test Conf.check_grid_interface(DGG.levelgrid(MINIMAL, 2))
        @test !Conf.check_grid_interface(DGG.levelgrid(broken(:cw_boundary), 2))
        problems = Conf.grid_interface_problems(DGG.levelgrid(broken(:cw_boundary), 2))
        @test any(contains("counter-clockwise"), problems)
    end

    # The containment half of the centroid law rests on a predicate that may
    # decline to answer (`spherical_ring_contains` returns `nothing` for
    # near-hemispherical rings), so without these two nothing would notice if it
    # stopped firing altogether.
    @testset "harness catches: centroid on the cell boundary" begin
        @test !Conf.check_grid_interface(DGG.levelgrid(broken(:centroid_on_boundary), 2))
        problems = Conf.grid_interface_problems(DGG.levelgrid(broken(:centroid_on_boundary), 2))
        @test any(contains("strictly interior"), problems)
    end

    @testset "harness catches: centroid outside its cell" begin
        @test !Conf.check_grid_interface(DGG.levelgrid(broken(:centroid_outside), 2))
        problems = Conf.grid_interface_problems(DGG.levelgrid(broken(:centroid_outside), 2))
        @test any(contains("outside its own cell boundary"), problems)
    end

    @testset "harness catches: children out of canonical order" begin
        @test !Conf.check_hierarchical_system(broken(:unsorted_children); n_samples = 3)
        problems = Conf.hierarchy_problems(broken(:unsorted_children), LevelIndex(1, 0))
        @test any(contains("ascending canonical order"), problems)
    end

    @testset "harness catches: cellposition off by one (bijection)" begin
        @test !Conf.check_grid_interface(DGG.levelgrid(broken(:bad_position), 2))
        problems = Conf.grid_interface_problems(DGG.levelgrid(broken(:bad_position), 2))
        @test any(contains("cellposition"), problems)
    end

    @testset "harness catches: descendant_range one too wide" begin
        @test Conf.check_descendant_ranges(MINIMAL)
        @test !Conf.check_descendant_ranges(broken(:lying_range))
        problems = Conf.descendant_range_problems(broken(:lying_range), LevelIndex(0, 0), 1,
            DGG.levelgrid(broken(:lying_range), 1))
        @test any(contains("not descendants"), problems)
    end

    @testset "harness catches: asymmetric neighbours" begin
        @test Conf.check_neighbors(DGG.levelgrid(MINIMAL, 1); n_samples = 24)
        @test !Conf.check_neighbors(DGG.levelgrid(broken(:asymmetric_neighbors), 1);
            n_samples = 24)
        # The break is one-sided: cell 0 forgets a neighbour that still lists it,
        # so the violation is visible from the forgotten cell, not from cell 0.
        dropped = last(DGG.neighbors(DGG.levelgrid(MINIMAL, 1), LevelIndex(1, 0)))
        problems = Conf.neighbor_problems(DGG.levelgrid(broken(:asymmetric_neighbors), 1),
            dropped)
        @test any(contains("not symmetric"), problems)
    end

    # The direct sweep detects an omission only from the omitted cell. The
    # two-hop closure also detects it from the cell that omitted the neighbour.
    @testset "harness catches: asymmetry from the forgetful cell (two-hop closure)" begin
        grid = DGG.levelgrid(broken(:asymmetric_neighbors), 1)
        culprit = LevelIndex(1, 0)   # the cell that drops one of its neighbours

        # Cell 0's remaining neighbours are symmetric, so the direct check passes.
        @test Conf.neighbor_problems(grid, culprit; two_hop = false) == String[]

        # The omitted cell is reachable through a reported neighbour and names cell 0.
        problems = Conf.neighbor_problems(grid, culprit)
        @test any(contains("two-hop closure"), problems)
        @test any(contains("not symmetric"), problems)

        # The reported problem identifies the omitted cell.
        dropped = last(DGG.neighbors(DGG.levelgrid(MINIMAL, 1), culprit))
        @test any(contains(string(dropped)), problems)

        # Both forms accept the conforming mock.
        @test Conf.neighbor_problems(DGG.levelgrid(MINIMAL, 1), culprit) == String[]
        @test Conf.check_neighbors(DGG.levelgrid(MINIMAL, 2); n_samples = 12)
    end

    # Set equality cannot distinguish rotational order from id order.
    @testset "harness catches: neighbours sorted by id instead of wound CCW" begin
        good = DGG.levelgrid(MINIMAL, 2)
        bad = DGG.levelgrid(broken(:sorted_neighbors), 2)
        @test Conf.check_neighbor_order(good; n_samples = 12)
        @test !Conf.check_neighbor_order(bad; n_samples = 12)

        # Sorting by id preserves shell membership and disjointness.
        @test Conf.check_neighbors(bad; n_samples = 12)
        c = DGG.cellindex(bad, 40)
        @test Set(DGG.neighbors(bad, c, 2)) == Set(DGG.neighbors(good, c, 2))

        problems = Conf.neighbor_order_problems(bad, c)
        @test any(contains("single counter-clockwise cycle"), problems)

        # The winding check itself, on the three orders that matter: the correct
        # one passes, the clockwise reversal fails, and a ROTATION passes —
        # where a ring starts is the system's own documented choice, and only its
        # direction is the contract's.
        ccw = collect(DGG.ring(good, c, 1))
        @test Conf.winding_problems(good, c, ccw) == String[]
        @test !isempty(Conf.winding_problems(good, c, reverse(ccw)))
        @test Conf.winding_problems(good, c, circshift(ccw, 3)) == String[]
        @test any(contains("turns clockwise"), Conf.winding_problems(good, c, reverse(ccw)))
        # Ring 1's winding is not negotiable, so opting out of the outer rings
        # does not make an id-sorted neighbourhood conform.
        @test !isempty(Conf.neighbor_order_problems(bad, c; require_rotational_rings = false))
        @test isempty(Conf.neighbor_order_problems(good, c; require_rotational_rings = false))
    end

    # A rotated ring preserves its set and winding but is not the disc's tail block.
    @testset "harness catches: ring that is not the tail block of the disc" begin
        good = DGG.levelgrid(MINIMAL, 2)
        bad = DGG.levelgrid(broken(:rotated_ring), 2)
        c = DGG.cellindex(bad, 40)

        @test Set(DGG.ring(bad, c, 1)) == Set(DGG.ring(good, c, 1))   # same cells
        @test Conf.check_neighbors(bad; n_samples = 12)               # same sets
        @test !Conf.check_neighbor_order(bad; n_samples = 12)

        problems = Conf.neighbor_order_problems(bad, c)
        @test any(contains("tail block"), problems)
        @test any(contains("rings concatenated outward"), problems)
        # ...and it is not caught by mistaking a rotation for a bad winding: a
        # rotated counter-clockwise cycle is still a counter-clockwise cycle.
        @test !any(contains("counter-clockwise cycle"), problems)
    end

    # Default mode asks whether a system's own methods obey; full mode asks
    # whether the methods that answer for a system without them obey. The two
    # can disagree, and this is the disagreement that matters: a tessellation
    # whose cells do not share corner coordinates keeps its own combinatorial
    # adjacency and loses the geometric one entirely — and an empty answer
    # satisfies every relational law there is.
    @testset "full mode catches a fallback that the specialization hides" begin
        gapped = CubeSystem(; maxlevel = 2, bug = :shrunken_cells, extras = true)

        _, specialized = capture() do
            test_hierarchical_system(gapped; n_samples = 4)
        end
        @test specialized == 0

        _, generic = capture() do
            test_generic_fallbacks(gapped; n_samples = 4)
        end
        @test generic > 0

        # The report names the emptiness. Without `require_nonempty` the
        # relational laws pass on the empty answer, which is why full mode
        # turns it on.
        wrapped = DGG.levelgrid(Conf.GenericFallbackSystem(gapped), 1)
        c = DGG.cellindex(wrapped, 1)
        @test isempty(DGG.neighbors(wrapped, c, 1))
        @test Conf.neighbor_problems(wrapped, c) == String[]
        @test any(contains("is empty on a grid"),
            Conf.neighbor_problems(wrapped, c; require_nonempty = true))
    end

    # The wrapper is the mechanism: without it, full mode would re-test the
    # specialized methods and report the fallbacks as certified.
    @testset "the wrapper hides specializations and forwards the contract" begin
        grid = DGG.levelgrid(FULL, 2)
        c = DGG.cellindex(grid, 3)
        p = DGG.cell_centroid(grid, c)
        @test Conf.has_nonfallback_method(DGG.cellat, grid, p)
        @test Conf.has_nonfallback_method(DGG.ancestor, FULL, c, 0)

        wrapped = Conf.GenericFallbackSystem(FULL)
        wgrid = DGG.levelgrid(wrapped, 2)
        @test Conf.dispatches_generically(DGG.cellat, wgrid, p)
        @test Conf.dispatches_generically(DGG.neighbors, wgrid, c, 1)
        @test Conf.dispatches_generically(DGG.ring, wgrid, c, 1)
        @test Conf.dispatches_generically(DGG.ancestor, wrapped, c, 0)
        @test Conf.dispatches_generically(DGG.descendants, wrapped, c, 0)

        # Everything the fallbacks consume is still the system's own answer.
        @test DGG.cell_boundary(wgrid, c) == DGG.cell_boundary(grid, c)
        @test DGG.cellposition(wgrid, c) == DGG.cellposition(grid, c)
        @test collect(DGG.children(wrapped, c)) == collect(DGG.children(FULL, c))
        @test DGG.node_extent(wrapped, c) == DGG.node_extent(FULL, c)
        @test DGG.descendant_range(wrapped, c, 3) == DGG.descendant_range(FULL, c, 3)
    end

    # A skip that does not say what it skipped reads, in a summary, exactly like
    # a known failure.
    @testset "each skip states its reason" begin
        logger = Test.TestLogger()
        Base.CoreLogging.with_logger(logger) do
            capture() do
                test_hierarchical_system(MINIMAL; n_samples = 3)
            end
        end
        msgs = [string(r.message) for r in logger.logs]
        @test any(m -> contains(m, "skipped: ancestor") && contains(m, "does not implement"),
            msgs)
        @test any(m -> contains(m, "skipped: descendants"), msgs)
        # ...and a system that implements them is not told anything.
        empty!(logger.logs)
        Base.CoreLogging.with_logger(logger) do
            capture() do
                test_hierarchical_system(FULL; n_samples = 3)
            end
        end
        @test isempty(logger.logs)
    end

    # The public testset entry point must surface collector failures.
    @testset "the @testset entry points fail on a broken mock" begin
        good_passes, good_fails = capture() do
            test_hierarchical_system(MINIMAL; n_samples = 4)
        end
        @test good_fails == 0
        @test good_passes > 0

        _, bad_fails = capture() do
            test_hierarchical_system(broken(:small_extent); n_samples = 4)
        end
        @test bad_fails > 0

        _, cw_fails = capture() do
            test_grid_interface(DGG.levelgrid(broken(:cw_boundary), 2); n_samples = 4)
        end
        @test cw_fails > 0

        # The order laws too, both halves: an id-sorted neighbourhood and a ring
        # that is not its disc's tail block have to surface through the entry
        # point, or they are laws only the predicates know about.
        _, order_fails = capture() do
            test_hierarchical_system(broken(:sorted_neighbors); n_samples = 4)
        end
        @test order_fails > 0

        _, tail_fails = capture() do
            test_hierarchical_system(broken(:rotated_ring); n_samples = 4)
        end
        @test tail_fails > 0

        # Opting the outer rings out does not opt ring 1 out with them.
        _, ring1_fails = capture() do
            test_hierarchical_system(broken(:sorted_neighbors); n_samples = 4,
                require_rotational_rings = false)
        end
        @test ring1_fails > 0
    end

    # `atol` and `unit_atol` exist because a system with exact (uninflated)
    # subtree caps puts cell corners exactly ON the cap rim, where containment is
    # a floating-point coin toss — HEALPix, in this package. A caller that cannot
    # reach the tolerance cannot use the harness at all, so these assert that the
    # kwargs really reach the checks rather than being decoration.
    @testset "atol and unit_atol reach the checks they name" begin
        # The covering `atol`: an extent that is genuinely too small is forgiven
        # by a slack wide enough to swallow the shortfall, and only by that.
        @test !isempty(Conf.covering_law_problems(broken(:small_extent)))
        @test isempty(Conf.covering_law_problems(broken(:small_extent); atol = 10.0))

        # ...and through the @testset entry point, which is what a caller uses.
        _, tight = capture() do
            test_hierarchical_system(broken(:small_extent); n_samples = 4)
        end
        _, loose = capture() do
            test_hierarchical_system(broken(:small_extent); n_samples = 4, atol = 10.0)
        end
        @test tight > 0
        @test loose == 0

        # `unit_atol` reaches `node_extent_problems`: a tolerance no centre can
        # meet fails the well-formedness set on the *correct* mock.
        _, impossible = capture() do
            test_hierarchical_system(MINIMAL; n_samples = 4, unit_atol = -1.0)
        end
        @test impossible > 0
        @test any(contains("off the unit sphere"),
            Conf.node_extent_problems(DGG.node_extent(MINIMAL, LevelIndex(1, 0));
                unit_atol = -1.0))
    end

    @testset "whole-suite predicates agree with the laws they aggregate" begin
        @test Conf.check_hierarchical_system(MINIMAL; n_samples = 3)
        @test !Conf.check_hierarchical_system(broken(:small_extent); n_samples = 3)
        @test !Conf.check_hierarchical_system(broken(:deep_drift); n_samples = 3)
    end

    @testset "harness is reproducible across runs" begin
        a = Conf.grid_interface_problems(DGG.levelgrid(broken(:bad_position), 2))
        b = Conf.grid_interface_problems(DGG.levelgrid(broken(:bad_position), 2))
        @test a == b
        @test !isempty(a)
        # A caller-supplied seed reaches every law, including the covering walk.
        seeded() = Conf.covering_law_problems(broken(:deep_drift);
            rng = Conf.Random.MersenneTwister(7), n_samples = 3)
        @test seeded() == seeded()
        @test !isempty(seeded())
    end
end

end # module TestConformance
