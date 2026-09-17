# Candidate locators: how a pair of spaces answers spatial queries.

"""
    CandidateLocator

How a weight build finds the cells of one space that a cell of the other may
touch. A plan carries one locator and hands it to every block it builds.

  - [`TreeLocator`](@ref) descends the two spaces' spatial trees. The default.
  - [`AnalyticLocator`](@ref) places points in closed form and walks cell
    neighbourhoods; it needs [`hasanalyticlocation`](@ref) on one side.

Every locator answers [`overlappairs`](@ref) and [`locatecell`](@ref).
"""
abstract type CandidateLocator end

"""
    TreeLocator()

The dual-tree `CandidateLocator`: candidate pairs come from a
depth-first descent of both spaces' `subtree`s, and a point is located
by descending one tree and testing the cell's chord polygon, the geometry the
conservative clipper measures. Works for every space.
"""
struct TreeLocator <: CandidateLocator end

"""
    AnalyticLocator()

The closed-form `CandidateLocator`: a point is located with
[`cellat`](@ref), the system's own answer and the one point methods use, and
candidate pairs come from a breadth-first walk over
`cellneighbors` seeded under each cell's centroid and corners, pruned
by spherical-cap intersection (`cellcap`). Builds no tree.

One side must answer `hasanalyticlocation` with `true`; that side is
walked. When both do, the coarser side (fewer cells) is walked and the finer
side supplies the seeds. Weights built with either locator are identical.
"""
struct AnalyticLocator <: CandidateLocator end

"""
    overlappairs(locator, dst_space, dst_inds, src_space, src_inds; threaded)
        -> Vector{Tuple{Int,Int}}

Candidate `(src, dst)` pairs, in the spaces' own local indices, for cells of
`src_inds` and `dst_inds` whose bounding caps intersect. Same contract as
`ConservativeRegridding.get_all_candidate_pairs`: every pair that truly
overlaps is listed, a listed pair need not overlap, and no pair is listed
twice. `threaded` runs the discovery on several tasks; the result is the same.
"""
function overlappairs end

function overlappairs(::TreeLocator, dst_space::RegridSpace, dst_inds,
    src_space::RegridSpace, src_inds; threaded::Bool)
    return ConservativeRegridding.get_all_candidate_pairs(GOCore.booltype(threaded),
        Extents.intersects, subtree(src_space, src_inds), subtree(dst_space, dst_inds))
end

function overlappairs(::AnalyticLocator, dst_space::RegridSpace, dst_inds,
    src_space::RegridSpace, src_inds; threaded::Bool)
    _walksdestination(dst_space, dst_inds, src_space, src_inds) &&
        return _walkpairs(src_space, src_inds, dst_space, dst_inds, threaded)
    pairs = _walkpairs(dst_space, dst_inds, src_space, src_inds, threaded)
    return map!(reverse, pairs, pairs)
end

# Walk into the side that can place a point and name its neighbours; when both
# can, into the coarser one, so the finer side supplies the seeds.
function _walksdestination(dst_space, dst_inds, src_space, src_inds)
    d, s = hasanalyticlocation(dst_space), hasanalyticlocation(src_space)
    d || s || throw(ArgumentError(
        "AnalyticLocator needs `cellat` and `cellneighbors` on one side " *
        "(`hasanalyticlocation`), and neither $(typeof(dst_space)) nor " *
        "$(typeof(src_space)) has them"))
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
    cap = cellcap(seedspace, s)
    empty!(queue)
    _seed!(queue, walkspace, cellcentroid(seedspace, s))
    for p in cellcorners(seedspace, s)
        _seed!(queue, walkspace, p)
    end
    head = 1
    while head <= length(queue)
        w = queue[head]
        head += 1
        Extents.intersects(cap, cellcap(walkspace, w)) || continue
        localindex(walkmap, w) == 0 || push!(out, (s, w))
        for q in cellneighbors(walkspace, w)
            q in queue || push!(queue, q)
        end
    end
    return out
end

function _seed!(queue, walkspace, p)
    w = cellat(walkspace, p)
    w === nothing && return queue
    w in queue || push!(queue, w)
    return queue
end

"""
    locatecell(locator, space, inds, p::UnitSphericalPoint) -> Int

The index within `inds` of the cell of `space` containing `p`, or `0` when no
cell of `inds` contains it.

The two locators answer from different authorities. An
[`AnalyticLocator`](@ref) returns the system's true cell: the answer
[`cellat`](@ref) gives, and the one point methods use. A [`TreeLocator`](@ref)
returns the cell whose chord polygon ([`getcell`](@ref)) contains the point,
consistent with the polygons the conservative clipper measures. Where a cell's
edge is a chart curve the polygon replaces by chords, a point in the sliver
between the two lies in the chord polygon of one neighbour and the true cell
of the other, and the two answers differ.
"""
function locatecell end

function locatecell(::AnalyticLocator, space::RegridSpace, inds, p::US.UnitSphericalPoint)
    i = cellat(space, p)
    return i === nothing ? 0 : localindex(indexmap(inds), i)
end

# Tree nodes bound their cells loosely, so every leaf the descent reaches is
# confirmed against the cell's own ring. Candidates are taken in index order so
# a point on a shared boundary always names the same cell.
function locatecell(::TreeLocator, space::RegridSpace, inds, p::US.UnitSphericalPoint)
    hits = STI.query(subtree(space, inds), node -> _nodecontains(node, p))
    imap = indexmap(inds)
    undecided = 0
    for i in sort!(hits)
        verdict = _cellcontains(getcell(space, i), p)
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
