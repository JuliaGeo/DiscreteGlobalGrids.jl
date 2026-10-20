# Candidate locators: how a pair of spaces answers spatial queries.

"""
    CandidateLocator

How a weight build finds the cells of one space that a cell of the other may
touch. A plan carries one locator and hands it to every block it builds.

  - [`TreeLocator`](@ref) descends the two grids' spatial trees. The default.
  - [`AnalyticLocator`](@ref) places points in closed form and walks cell
    neighbourhoods; it needs `GridInterface.location_is_analytic` on one side.

Every locator answers [`overlappairs`](@ref) and [`locatecell`](@ref).
"""
abstract type CandidateLocator end

"""
    TreeLocator()

The dual-tree `CandidateLocator`: candidate pairs come from a
depth-first descent of both grids' `subtree`s, and a point is located
by descending one tree and testing the cell's chord polygon, the geometry the
conservative clipper measures. Works for every grid.
"""
struct TreeLocator <: CandidateLocator end

"""
    AnalyticLocator()

The closed-form `CandidateLocator`: a point is located with
`GridInterface.cellat`, the system's own answer and the one point methods use,
and candidate pairs come from a breadth-first walk over
`GridInterface.neighbors` seeded under each cell's centroid and corners, pruned
by spherical-cap intersection (`GridInterface.cell_extent`). Builds no tree.

One side must answer `GridInterface.location_is_analytic` with `true`; that
side is walked. When both do, the coarser side (fewer cells) is walked and the
finer side supplies the seeds. Weights built with either locator are identical.
"""
struct AnalyticLocator <: CandidateLocator end

"""
    overlappairs(locator, dst_grid, dst_inds, src_grid, src_inds; threaded)
        -> Vector{Tuple{Int,Int}}

Candidate `(src, dst)` pairs, in the grids' own positions, for cells of
`src_inds` and `dst_inds` whose bounding caps intersect. Same contract as the
dual-tree search of `ConservativeRegridding.intersection_areas`: every pair
that truly overlaps is listed, a listed pair need not overlap, and no pair is
listed twice. `threaded` runs the discovery on several tasks; the result is the
same.
"""
function overlappairs end

function overlappairs(::TreeLocator, dst_space, dst_inds, src_space, src_inds;
    threaded::Bool)
    return ConservativeRegridding.candidate_pairs(GOCore.booltype(threaded),
        _subtree(dst_space, dst_inds; leaf_capacity = GridInterface.ForSweeps()),
        _subtree(src_space, src_inds; leaf_capacity = GridInterface.ForSweeps()))
end

function overlappairs(::AnalyticLocator, dst_space, dst_inds, src_space, src_inds;
    threaded::Bool)
    _walksdestination(dst_space, dst_inds, src_space, src_inds) &&
        return _walkpairs(src_space, src_inds, dst_space, dst_inds, threaded)
    pairs = _walkpairs(dst_space, dst_inds, src_space, src_inds, threaded)
    return map!(reverse, pairs, pairs)
end

# Walk into the side that can place a point and name its neighbours; when both
# can, into the coarser one, so the finer side supplies the seeds.
function _walksdestination(dst_space, dst_inds, src_space, src_inds)
    d = GridInterface.location_is_analytic(dst_space)
    s = GridInterface.location_is_analytic(src_space)
    d || s || throw(ArgumentError(
        "AnalyticLocator needs closed-form `cellat` and `neighbors` on one " *
        "side (`GridInterface.location_is_analytic`), and neither " *
        "$(typeof(dst_space)) nor $(typeof(src_space)) has them"))
    return d && (!s || length(dst_inds) <= length(src_inds))
end

# `(seed, walk)` pairs for every seed cell, in seed order.
function _walkpairs(seedspace, seedinds, walkspace, walkinds, threaded::Bool)
    walkmap = indexmap(walkinds)
    n = length(seedinds)
    ntasks = threaded ? min(n, Threads.nthreads() * 8) : 1
    ntasks <= 1 && return _walkrange!(Tuple{Int,Int}[], Int[],
        seedspace, seedinds, 1:n, walkspace, walkmap)
    tasks = [StableTasks.@spawn _walkrange!(Tuple{Int,Int}[], Int[],
                 seedspace, seedinds, r, walkspace, walkmap)
             for r in _splitrange(n, ntasks)]
    return reduce(vcat, fetch.(tasks))
end

_splitrange(n::Int, k::Int) =
    (round(Int, (i - 1) * n / k) + 1:round(Int, i * n / k) for i in 1:k)

function _walkrange!(out, queue, seedspace, seedinds, positions, walkspace, walkmap)
    for k in positions
        _walk!(out, queue, seedspace, Int(seedinds[k]), walkspace, walkmap)
    end
    return out
end

# One seed cell: breadth-first over the walk side's one-ring from the cells
# under the seed's centroid and corners, expanding a cell only when its cap
# meets the seed's. The queue is never popped, so it doubles as the visited set;
# a walk touches a few dozen cells, and a linear scan of them is cheaper than
# a per-cell stamp array.
function _walk!(out, queue, seedspace, s::Int, walkspace, walkmap)
    cap = _walkcap(seedspace, s)
    empty!(queue)
    _seed!(queue, walkspace, GridInterface.cell_centroid(seedspace, s))
    for p in GridInterface.cell_corners(seedspace, s)
        _seed!(queue, walkspace, p)
    end
    head = 1
    while head <= length(queue)
        w = queue[head]
        head += 1
        Extents.intersects(cap, _walkcap(walkspace, w)) || continue
        localindex(walkmap, w) == 0 || push!(out, (s, w))
        for q in GridInterface.neighbors(walkspace, w)
            Int(q) in queue || push!(queue, Int(q))
        end
    end
    return out
end

# The cap the walk prunes by: the grid's own cell extent, grown by a small
# relative margin so rounding in the cap arithmetic never drops a touching
# cell. A cap wider than a hemisphere is not geodesically convex and becomes
# the whole sphere.
@inline function _walkcap(grid, i::Int)
    cap = convert(Cap, GridInterface.cell_extent(grid, i))
    cap.radius > Float64(pi) / 2 && return _WHOLE_SPHERE
    return Extents.grow(cap, 5e-5)
end

function _seed!(queue, walkspace, p)
    w = GridInterface.cellat(walkspace, p)
    w === nothing && return queue
    Int(w) in queue || push!(queue, Int(w))
    return queue
end

"""
    locatecell(locator, grid, inds, p::UnitSphericalPoint) -> Int

The index within `inds` of the cell of `grid` containing `p`, or `0` when no
cell of `inds` contains it.

The two locators answer from different authorities. An
[`AnalyticLocator`](@ref) returns the system's true cell: the answer
`GridInterface.cellat` gives, and the one point methods use. A
[`TreeLocator`](@ref) returns the cell whose chord polygon
(`GridInterface.cell_geometry`) contains the point, consistent with the
polygons the conservative clipper measures. Where a cell's
edge is a chart curve the polygon replaces by chords, a point in the sliver
between the two lies in the chord polygon of one neighbour and the true cell
of the other, and the two answers differ.
"""
function locatecell end

function locatecell(::AnalyticLocator, space, inds, p::US.UnitSphericalPoint)
    i = GridInterface.cellat(space, p)
    return i === nothing ? 0 : localindex(indexmap(inds), Int(i))
end

# Tree nodes bound their cells loosely, so every leaf the descent reaches is
# confirmed against the cell's own ring. Candidates are taken in index order so
# a point on a shared boundary always names the same cell.
function locatecell(::TreeLocator, space, inds, p::US.UnitSphericalPoint)
    hits = STI.query(_subtree(space, inds; leaf_capacity = GridInterface.ForQueries()),
        node -> _nodecontains(node, p))
    imap = indexmap(inds)
    undecided = 0
    for i in sort!(hits)
        verdict = _cellcontains(GridInterface.cell_geometry(space, i), p)
        verdict === true && return localindex(imap, i)
        verdict === nothing && undecided == 0 && (undecided = i)
    end
    return undecided == 0 ? 0 : localindex(imap, undecided)
end

_nodecontains(cap::SphericalCap, p) = US._contains(cap, p)
_nodecontains(box::Extents.Extent, p) =
    box.X[1] <= p[1] <= box.X[2] && box.Y[1] <= p[2] <= box.Y[2] &&
    box.Z[1] <= p[3] <= box.Z[2]

function _cellcontains(cell, p)
    ring = _ringpoints(cell)
    return US.spherical_ring_contains(ring, length(ring), p)
end

# A cell's exterior ring without its closing repeat.
function _ringpoints(cell)
    ring = GI.getexterior(cell)
    n = GI.npoint(ring) - 1
    return USPoint[USPoint(GI.x(q), GI.y(q), GI.z(q)) for q in Iterators.take(GI.getpoint(ring), n)]
end
