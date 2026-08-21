# ISEA4T is a true nested partition: twenty Snyder face triangles, each split
# into four equal chart triangles by edge midpoints.  Its canonical face/path
# order is package-defined; the sealed DGGRID SEQNUM dump does not expose a
# face/path crosswalk, so no unsupported SEQNUM compatibility is claimed.

"""`ISEA4TSystem`: nested aperture-4 Snyder ISEA triangular refinement."""
struct ISEA4TSystem <: DGG.AbstractHierarchicalGridSystem
    orientation::Orientation
end
ISEA4TSystem() = ISEA4TSystem(ORIENT_IDENTITY)
Base.show(io::IO, sys::ISEA4TSystem) =
    print(io, "ISEA4TSystem(", sys.orientation.identity ? "" : sys.orientation, ")")

DGG.cellindextype(::ISEA4TSystem) = DGG.LevelIndex
# Float64 Snyder inversion and edge-crossing topology remain reliable through
# level 12. The packed face/path index has spare bits beyond that, but exposing
# numerically invalid geometry would violate the grid interface.
DGG.levels(::ISEA4TSystem) = 0:12
DGG.has_sorted_subtrees(::ISEA4TSystem) = true
DGG.maxneighbors(::ISEA4TSystem, ::DGG.Edge) = 3
DGG.maxneighbors(::ISEA4TSystem, ::DGG.Vertex) = 12
DGG.ncells(::ISEA4TSystem, l::Integer) = Int(20 * POW4[Int(l) + 1])
DGG.rootcells(::ISEA4TSystem) = [DGG.LevelIndex(0, f) for f in 0:19]
DGG.cellindex(::ISEA4TSystem, l::Integer, i::Int) = DGG.LevelIndex(l, i - 1)

function DGG.cellposition(sys::ISEA4TSystem, c::DGG.LevelIndex)
    DGG.level(c) in DGG.levels(sys) || return nothing
    0 <= c.index < 20 * POW4[DGG.level(c) + 1] || return nothing
    return Int(c.index + 1)
end

function Base.parent(::ISEA4TSystem, c::DGG.LevelIndex)
    r = DGG.level(c); r > 0 || throw(ArgumentError("an ISEA4T root has no parent"))
    return DGG.LevelIndex(r - 1, c.index ÷ 4)
end

function DGG.children(sys::ISEA4TSystem, c::DGG.LevelIndex)
    r = DGG.level(c); r < DGG.maxlevel(sys) || throw(ArgumentError("cell is at max level"))
    return [DGG.LevelIndex(r + 1, 4c.index + d) for d in 0:3]
end

function DGG.ancestor(::ISEA4TSystem, c::DGG.LevelIndex, l::Integer)
    target = Int(l); r = DGG.level(c)
    0 <= target <= r || throw(ArgumentError("ancestor level must be in 0:$r"))
    return DGG.LevelIndex(target, c.index >> (2 * (r - target)))
end

function DGG.descendant_range(::ISEA4TSystem, c::DGG.LevelIndex, l::Integer)
    target = Int(l); r = DGG.level(c)
    r <= target <= 12 || throw(ArgumentError("descendant level must be in $r:12"))
    shift = 2 * (target - r)
    return Int((c.index << shift) + 1):Int(((c.index + 1) << shift))
end

@inline function _root_triangle(face::Int)
    c = FACES[face + 1].corner
    # FACES stores vertex ids ascending, not a winding convention.
    return imag(conj(c[2] - c[1]) * (c[3] - c[1])) > 0 ? c : (c[1], c[3], c[2])
end

@inline function _child_triangle(t, d::Int)
    a, b, c = t; ab = (a + b) / 2; bc = (b + c) / 2; ca = (c + a) / 2
    d == 0 && return (a, ab, ca)
    d == 1 && return (ab, b, bc)
    d == 2 && return (ca, bc, c)
    return (ab, bc, ca)
end

function _triangle(c::DGG.LevelIndex)
    r = DGG.level(c); n = POW4[r + 1]
    0 <= c.index < 20n || throw(ArgumentError("ISEA4T index $(c.index) is out of range"))
    face, path = divrem(c.index, n); t = _root_triangle(Int(face))
    for k in r:-1:1
        d = Int((path >> (2 * (k - 1))) & 3)
        t = _child_triangle(t, d)
    end
    return Int(face), t
end

function DGG.cell_centroid(sys::ISEA4TSystem, c::DGG.LevelIndex)
    face, t = _triangle(c); w = (t[1] + t[2] + t[3]) / 3
    p = from_grid(sys.orientation, snyder_inv_xyz(face, w))
    return USPoint(p[1], p[2], p[3])
end

function DGG.cell_boundary(sys::ISEA4TSystem, c::DGG.LevelIndex)
    face, t = _triangle(c)
    segments = 8
    out = Vector{USPoint{Float64}}(undef, 3segments); k = 0
    for j in 1:3
        a, b = t[j], t[mod1(j + 1, 3)]
        for i in 0:(segments - 1)
            w = a + (i / segments) * (b - a)
            p = from_grid(sys.orientation, snyder_inv_xyz(face, w))
            out[k += 1] = USPoint(p[1], p[2], p[3])
        end
    end
    return out
end

function DGG.cell_area(g::LevelGrid{ISEA4TSystem}, c::DGG.LevelIndex)
    DGG.level(c) == g.level || throw(ArgumentError("cell and grid levels differ"))
    DGG.cellposition(g, c) === nothing && throw(ArgumentError("invalid ISEA4T cell"))
    return pi / (5 * 4.0^g.level)
end

@inline function _inside_triangle(w::ComplexF64, t; eps=2e-14)
    a, b, c = t
    s1 = imag(conj(b - a) * (w - a))
    s2 = imag(conj(c - b) * (w - b))
    s3 = imag(conj(a - c) * (w - c))
    return s1 >= -eps && s2 >= -eps && s3 >= -eps
end

function DGG.cellat(g::LevelGrid{ISEA4TSystem}, p::GO.UnitSphericalPoint)
    q = to_grid(g.system.orientation, (Float64(p[1]), Float64(p[2]), Float64(p[3])))
    face, w = snyder_fwd(q); t = _root_triangle(face); path = Int64(0)
    for _ in 1:g.level
        picked = 3
        for d in 0:2
            ct = _child_triangle(t, d)
            if _inside_triangle(w, ct)
                picked = d; break
            end
        end
        t = _child_triangle(t, picked); path = 4path + picked
    end
    return DGG.LevelIndex(g.level, face * POW4[g.level + 1] + path)
end

# Cross each chart edge at its exact inverse-projected midpoint.  The outward
# tangent is obtained from the cell's interior point, so the same construction
# works on icosahedron seams without a face-neighbor table; `cellat` performs
# the global seam ownership step.
function _tri_edge_neighbors(g::LevelGrid{ISEA4TSystem}, c::DGG.LevelIndex)
    DGG.level(c) == g.level || throw(ArgumentError("cell and grid levels differ"))
    DGG.cellposition(g, c) === nothing && throw(ArgumentError("cell is not in the grid"))
    face, t = _triangle(c)
    centre = Tuple(DGG.cell_centroid(g, c))
    out = DGG.LevelIndex[]
    for j in 1:3
        w = (t[j] + t[mod1(j + 1, 3)]) / 2
        m = from_grid(g.system.orientation, snyder_inv_xyz(face, w))
        # `m-centre` has an outward component at the edge midpoint.  A relative
        # 1e-6 step crosses only this edge at every supported level.
        q = vnormalize(vadd(m, vscale(vsub(m, centre), 1e-6)))
        n = DGG.cellat(g, USPoint(q[1], q[2], q[3]))
        n == c && error("ISEA4T edge crossing failed at $c edge $j")
        n in out || push!(out, n)
    end
    return out
end

@inline function _same_vertex(a, b)
    return vdot(a, b) >= 1 - 2e-13
end

"""
    _tri_corners(sys, c) -> NTuple{3,NTuple{3,Float64}}

The cell's three corners on the sphere.

Vertex incidence is a statement about corners, so this is what the vertex star
tests — three inverse projections, against the twenty-four
[`cell_boundary`](@ref) spends densifying the edges between them. The densified
points can only ever confirm a shared edge, which the corners already report,
because a conforming subdivision has no T-junctions.
"""
function _tri_corners(sys::ISEA4TSystem, c::DGG.LevelIndex)
    face, t = _triangle(c)
    return ntuple(j -> from_grid(sys.orientation, snyder_inv_xyz(face, t[j])), 3)
end

function _tri_sort_ccw!(cells, g, subject)
    length(cells) <= 1 && return cells
    o = Tuple(DGG.cell_centroid(g, subject))
    anchor = Tuple(DGG.cell_centroid(g, minimum(cells)))
    t = vsub(anchor, vscale(o, vdot(anchor, o)))
    e1 = vnormalize(t); e2 = vcross(o, e1)
    az(c) = begin
        p = Tuple(DGG.cell_centroid(g, c)); d = vsub(p, o)
        mod(atan(vdot(d, e2), vdot(d, e1)), 2pi)
    end
    zeroaz = az(minimum(cells))
    # Keyed once per cell rather than once per comparison; `az` is a centroid
    # construction, and `by` runs on every comparison the sort makes.
    keys = map(c -> (mod(az(c) - zeroaz, 2pi), c), cells)
    permute!(cells, sortperm(keys))
    return cells
end

function _tri_neighbors1(g::LevelGrid{ISEA4TSystem}, c::DGG.LevelIndex,
        connectivity::DGG.Connectivity)
    edge = _tri_edge_neighbors(g, c)
    connectivity isa DGG.Edge && return _tri_sort_ccw!(edge, g, c)
    # All triangles at an original/refined mesh vertex occur within three
    # edge-adjacency steps.  Filter that bounded disc by exact chart-vertex
    # incidence; this handles valence 5 at icosahedron vertices and valence 6
    # elsewhere without a global mesh table.
    seen = Set{DGG.LevelIndex}((c,)); frontier = DGG.LevelIndex[c]
    candidates = DGG.LevelIndex[]
    for _ in 1:3
        shell = DGG.LevelIndex[]
        for x in frontier, n in _tri_edge_neighbors(g, x)
            if !(n in seen)
                push!(seen, n); push!(shell, n); push!(candidates, n)
            end
        end
        frontier = shell
    end
    ring = _tri_corners(g.system, c)
    # Each candidate's corners are built once. In the flattened generator this
    # replaces, `cell_boundary(g, n)` was the INNER iterable and so was rebuilt
    # for every vertex of the subject's ring — twenty-four densified rings per
    # candidate, for a test that three corners settle.
    out = DGG.LevelIndex[]
    for n in candidates
        corners = _tri_corners(g.system, n)
        any(_same_vertex(a, b) for a in ring, b in corners) && push!(out, n)
    end
    return _tri_sort_ccw!(out, g, c)
end

function _tri_shells(g::LevelGrid{ISEA4TSystem}, c::DGG.LevelIndex, steps::Int,
        connectivity::DGG.Connectivity)
    steps >= 0 || throw(ArgumentError("k must be non-negative"))
    steps == 0 && return Vector{DGG.LevelIndex}[]
    seen = Set{DGG.LevelIndex}((c,)); frontier = DGG.LevelIndex[c]
    shells = Vector{DGG.LevelIndex}[]
    for _ in 1:steps
        shell = DGG.LevelIndex[]
        for x in frontier, n in _tri_neighbors1(g, x, connectivity)
            if !(n in seen)
                push!(seen, n); push!(shell, n)
            end
        end
        _tri_sort_ccw!(shell, g, c)
        push!(shells, shell); frontier = shell
    end
    return shells
end

function DGG.neighbors(g::LevelGrid{ISEA4TSystem}, c::DGG.LevelIndex, k::Integer=1;
        connectivity::DGG.Connectivity=DGG.Vertex())
    shells = _tri_shells(g, c, Int(k), connectivity)
    return isempty(shells) ? DGG.LevelIndex[] : reduce(vcat, shells)
end

function DGG.ring(g::LevelGrid{ISEA4TSystem}, c::DGG.LevelIndex, k::Integer;
        connectivity::DGG.Connectivity=DGG.Vertex())
    steps = Int(k); steps >= 0 || throw(ArgumentError("k must be non-negative"))
    steps == 0 && return DGG.LevelIndex[c]
    return _tri_shells(g, c, steps, connectivity)[steps]
end
