# ISEA3H and ISEA4H are central-place refinements: their chart Voronoi cells
# are not nested.  `parent` below is consequently the published Z3/ZORDER
# prefix parent (a primary indexing tree), not a claim of polygon containment.

"""`ISEA3HSystem`: aperture-3 Snyder ISEA hexagons with canonical Z3 ids."""
struct ISEA3HSystem <: DGG.AbstractHierarchicalGridSystem
    orientation::Orientation
end
ISEA3HSystem() = ISEA3HSystem(ORIENT_IDENTITY)

"""`ISEA4HSystem`: aperture-4 Snyder ISEA hexagons with prefix `LevelIndex` ids."""
struct ISEA4HSystem <: DGG.AbstractHierarchicalGridSystem
    orientation::Orientation
end
ISEA4HSystem() = ISEA4HSystem(ORIENT_IDENTITY)

const POW3 = ntuple(k -> Int64(3)^(k - 1), 31)
const POW4 = ntuple(k -> Int64(4)^(k - 1), 30)

@inline _hexcount(A::Int, level::Int) = 10 * (A == 3 ? POW3[level + 1] : POW4[level + 1]) + 2
@inline _polar(root::Int) = root == 0 || root == 11

DGG.cellindextype(::ISEA3HSystem) = Z3Cell
DGG.levels(::ISEA3HSystem) = 0:30
DGG.has_sorted_subtrees(::ISEA3HSystem) = true
DGG.max_neighbors(::ISEA3HSystem, ::DGG.Connectivity) = 6
# The primary prefix tree is not a spatial containment tree.  In the planar
# chart its limiting center overhang plus the ancestor boundary is at most
# `1 + sqrt(3)/(sqrt(3)-1) = 3.3661` ancestor circumradii.  The 4.5 factor adds
# Snyder distortion/seam margin.  Exhaustive oracle levels 0:5 measured 1.9664;
# 22,977 seeded boundary probes through level 30 measured 1.9704.  Root caps
# consequently exceed a hemisphere (and saturate at the full sphere), which is
# intentional: a convex cap cannot cover these non-spatial prefix subtrees.
DGG.cap_inflation(::ISEA3HSystem) = 4.5
DGG.ncells(::ISEA3HSystem, l::Integer) = Int(_hexcount(3, Int(l)))

DGG.cellindextype(::ISEA4HSystem) = DGG.LevelIndex
DGG.levels(::ISEA4HSystem) = 0:29
DGG.has_sorted_subtrees(::ISEA4HSystem) = true
DGG.max_neighbors(::ISEA4HSystem, ::DGG.Connectivity) = 6
# Planar limiting ratio `1 + sqrt(3) = 2.7321`, with Snyder/seam margin.
# Exhaustive oracle levels 0:4 measured 1.6177; 22,977 seeded probes through
# level 29 measured 1.6983.  See the ISEA3H note above about non-convex roots.
DGG.cap_inflation(::ISEA4HSystem) = 3.5
DGG.ncells(::ISEA4HSystem, l::Integer) = Int(_hexcount(4, Int(l)))

const CentralPlaceHexSystem = Union{ISEA3HSystem,ISEA4HSystem}

function _checked_multiorder(sys::CentralPlaceHexSystem, coverage::DGG.MultiOrderCoverage;
        level::Union{Integer,Nothing}=nothing,
        maxcells::Union{Integer,Nothing}=nothing,
        maxlevel::Union{Integer,Nothing}=nothing)
    maxcells === nothing || throw(ArgumentError(
        "budget MultiOrderCoverage is unavailable for $(nameof(typeof(sys))): " *
        "its canonical prefix parent is non-spatial and its children do not " *
        "cover the parent footprint; use the fixed `level` mode"))
    return DGG.Fallbacks._multi_order_query(sys, coverage.target, level, maxcells, maxlevel)
end

DGG.query(sys::CentralPlaceHexSystem, coverage::DGG.MultiOrderCoverage;
    level::Union{Integer,Nothing}=nothing, maxcells::Union{Integer,Nothing}=nothing,
    maxlevel::Union{Integer,Nothing}=nothing) =
    _checked_multiorder(sys, coverage; level, maxcells, maxlevel)

DGG.MultiOrderCellSet(sys::CentralPlaceHexSystem, coverage::DGG.MultiOrderCoverage;
    level::Union{Integer,Nothing}=nothing, maxcells::Union{Integer,Nothing}=nothing,
    maxlevel::Union{Integer,Nothing}=nothing) =
    _checked_multiorder(sys, coverage; level, maxcells, maxlevel)

is_pentagon(c::Z3Cell) = _valid_z3(c.id) && _z3path(c) == 0

@inline function _i4decode(index::Int64, level::Int)
    width = POW4[level + 1]
    index == 0 && return (0, Int64(0))
    index == 10width + 1 && return (11, Int64(0))
    1 <= index <= 10width || throw(ArgumentError("ISEA4H index $index is out of range"))
    q, p = divrem(index - 1, width)
    return (Int(q) + 1, p)
end

@inline function _i4encode(root::Int, path::Int64, level::Int)
    root == 0 && return Int64(0)
    root == 11 && return 10 * POW4[level + 1] + 1
    return 1 + (root - 1) * POW4[level + 1] + path
end

is_pentagon(c::DGG.LevelIndex) = begin
    root, path = _i4decode(c.index, DGG.level(c)); path == 0
end

function DGG.rootcells(::ISEA3HSystem)
    return [_z3from(root, 0, 0) for root in 0:11]
end
DGG.rootcells(::ISEA4HSystem) = [DGG.LevelIndex(0, root) for root in 0:11]

function DGG.cellindex(::ISEA3HSystem, level::Integer, pos::Int)
    r = Int(level); width = POW3[r + 1]; q = pos - 1
    q == 0 && return _z3from(0, 0, r)
    q == 10width + 1 && return _z3from(11, 0, r)
    block, path = divrem(q - 1, width)
    return _z3from(Int(block) + 1, path, r)
end

function DGG.cellposition(::ISEA3HSystem, c::Z3Cell)
    _valid_z3(c.id) || return nothing
    r = DGG.level(c); root = _z3root(c.id); path = _z3path(c); width = POW3[r + 1]
    root == 0 && return 1
    root == 11 && return Int(10width + 2)
    return Int(1 + (root - 1) * width + path + 1)
end

DGG.cellindex(::ISEA4HSystem, l::Integer, pos::Int) = DGG.LevelIndex(l, pos - 1)
function DGG.cellposition(::ISEA4HSystem, c::DGG.LevelIndex)
    0 <= c.index < _hexcount(4, DGG.level(c)) || return nothing
    return Int(c.index + 1)
end

function Base.parent(::ISEA3HSystem, c::Z3Cell)
    _valid_z3(c.id) || throw(ArgumentError("invalid Z3 cell"))
    r = DGG.level(c); r > 0 || throw(ArgumentError("an ISEA3H root has no parent"))
    return Z3Cell(c.id | (UInt64(3) << _z3shift(r)))
end

function DGG.children(::ISEA3HSystem, c::Z3Cell)
    _valid_z3(c.id) || throw(ArgumentError("invalid Z3 cell"))
    r = DGG.level(c); r < Z3_DIGITS || throw(ArgumentError("cell is at max level"))
    sh = _z3shift(r + 1); cleared = c.id & ~(UInt64(3) << sh)
    ds = _polar(_z3root(c.id)) ? (0:0) : (0:2)
    return [Z3Cell(cleared | (UInt64(d) << sh)) for d in ds]
end

function Base.parent(::ISEA4HSystem, c::DGG.LevelIndex)
    r = DGG.level(c); r > 0 || throw(ArgumentError("an ISEA4H root has no parent"))
    root, path = _i4decode(c.index, r)
    return DGG.LevelIndex(r - 1, _i4encode(root, path ÷ 4, r - 1))
end

function DGG.children(sys::ISEA4HSystem, c::DGG.LevelIndex)
    r = DGG.level(c); r < DGG.max_level(sys) || throw(ArgumentError("cell is at max level"))
    root, path = _i4decode(c.index, r)
    ds = _polar(root) ? (0:0) : (0:3)
    return [DGG.LevelIndex(r + 1, _i4encode(root, 4path + d, r + 1)) for d in ds]
end

function DGG.ancestor(::ISEA3HSystem, c::Z3Cell, l::Integer)
    target = Int(l); r = DGG.level(c)
    0 <= target <= r || throw(ArgumentError("ancestor level must be in 0:$r"))
    z = c.id
    for k in (target + 1):r
        z |= UInt64(3) << _z3shift(k)
    end
    return Z3Cell(z)
end

function DGG.ancestor(::ISEA4HSystem, c::DGG.LevelIndex, l::Integer)
    target = Int(l); r = DGG.level(c)
    0 <= target <= r || throw(ArgumentError("ancestor level must be in 0:$r"))
    root, path = _i4decode(c.index, r)
    return DGG.LevelIndex(target, _i4encode(root, path ÷ Int64(4)^(r - target), target))
end

function DGG.descendant_range(::ISEA3HSystem, c::Z3Cell, l::Integer)
    target = Int(l); r = DGG.level(c); target >= r || throw(ArgumentError("target is above cell"))
    target <= 30 || throw(ArgumentError("target is past max level"))
    root = _z3root(c.id); path = _z3path(c); mul = POW3[target - r + 1]
    lo = if root == 0
        1
    elseif root == 11
        10 * POW3[target + 1] + 2
    else
        2 + (root - 1) * POW3[target + 1] + path * mul
    end
    hi = _polar(root) ? lo : lo + mul - 1
    return Int(lo):Int(hi)
end

function DGG.descendant_range(::ISEA4HSystem, c::DGG.LevelIndex, l::Integer)
    target = Int(l); r = DGG.level(c); target >= r || throw(ArgumentError("target is above cell"))
    target <= 29 || throw(ArgumentError("target is past max level"))
    root, path = _i4decode(c.index, r); mul = POW4[target - r + 1]
    loid = _i4encode(root, path * mul, target)
    hiid = _polar(root) ? loid : loid + mul - 1
    return Int(loid + 1):Int(hiid + 1)
end

@inline function _digits(path::Int64, r::Int, A::Int)
    out = Vector{Int}(undef, r)
    for k in r:-1:1
        path, d = divrem(path, A)
        out[k] = Int(d)
    end
    return out
end

# DGGRID-compatible central-place center construction.  A4 is the ordinary
# class-I triangular lattice.  A3 alternates class II/class I; the even-level
# digit gauge depends on the preceding ternary digit.  These are the
# modified-balanced-ternary directions of Sahr (2008), expressed directly in
# the base vertex's five-face development.
function _hex_dev(A::Int, root::Int, path::Int64, r::Int)
    ds = _digits(path, r, A)
    u = complex(0.0, 0.0)
    south = 6 <= root <= 10
    if A == 4
        dirs = (0.0, 180.0, 60.0, 120.0)
        for k in 1:r
            d = ds[k]; d == 0 && continue
            theta = dirs[d + 1] - (south ? 60.0 : 0.0)
            u += (L_PLANE / 2.0^k) * cis(deg2rad(theta))
        end
    else
        for k in 1:r
            d = ds[k]; d == 0 && continue
            theta = if isodd(k)
                d == 1 ? 90.0 : 150.0
            else
                prev = ds[k - 1]
                prev == 0 ? (d == 1 ? 60.0 : 120.0) :
                    prev == 1 ? (d == 1 ? 180.0 : 0.0) :
                    (d == 1 ? 240.0 : 300.0)
            end
            u += (L_PLANE / sqrt(3.0)^k) * cis(deg2rad(theta - (south ? 60.0 : 0.0)))
        end
    end
    if u != 0
        psi = mod(rad2deg(angle(u)), 360.0)
        psi >= 300.0 - 1e-10 && (u *= cis(-pi / 3))
    end
    return u
end

@inline function _hexparts(c::Z3Cell)
    _valid_z3(c.id) || throw(ArgumentError("invalid Z3 cell"))
    return (_z3root(c.id), _z3path(c), DGG.level(c))
end
@inline function _hexparts(c::DGG.LevelIndex)
    r = DGG.level(c); root, path = _i4decode(c.index, r); return (root, path, r)
end

function _hex_center(sys, c, A::Int)
    root, path, r = _hexparts(c)
    p = path == 0 ? vertex(root) : dev_to_xyz(root, _hex_dev(A, root, path, r))
    p = from_grid(sys.orientation, p)
    return USPoint(p[1], p[2], p[3])
end

DGG.cell_centroid(sys::ISEA3HSystem, c::Z3Cell) = _hex_center(sys, c, 3)
DGG.cell_centroid(sys::ISEA4HSystem, c::DGG.LevelIndex) = _hex_center(sys, c, 4)

function _hex_boundary(sys, c, A::Int)
    root, path, r = _hexparts(c); center = _hex_dev(A, root, path, r)
    pent = path == 0
    n = pent ? 5 : 6
    step = L_PLANE / sqrt(Float64(A))^r
    radius = step / SQRT3
    lattice_rotation = A == 3 && isodd(r) ? 30.0 : 0.0
    corner0 = lattice_rotation + 30.0
    # Starting at 150 degrees agrees with the sealed DGGRID rings for the
    # standard frame.  Rotation of a cyclic ring has no geometric meaning.
    start = corner0
    corners = Vector{ComplexF64}(undef, n)
    centerang = center == 0 ? 0.0 : mod(rad2deg(angle(center)), 360.0)
    for j in 0:(n - 1)
        angledeg = if pent
            if A == 3 && isodd(r)
                (180.0, 240.0, 0.0, 60.0, 120.0)[j + 1]
            else
                (150.0, 210.0, 270.0, 30.0, 90.0)[j + 1]
            end
        else
            start + 60j
        end
        u = center + radius * cis(deg2rad(angledeg))
        psi = dev_angle_deg(u)
        if !pent && psi >= 300.0 - 1e-10
            u *= cis(deg2rad(centerang > 150.0 ? 60.0 : -60.0))
        end
        corners[j + 1] = u
    end
    # Hex cells can straddle the five-face development cut, where the same
    # physical edge has several planar representatives. Canonicalize a paired
    # edge by choosing the lowest icosahedron face containing both endpoints,
    # then interpolate in that one Snyder chart. Both incident cells therefore
    # choose the same curved edge independently. Edges with no common face use
    # the crack-free great-circle fallback.
    segments = 8
    out = Vector{USPoint{Float64}}(undef, n * segments)
    cornergrid = [dev_to_xyz(root, u) for u in corners]
    cornerxyz = [from_grid(sys.orientation, p) for p in cornergrid]
    k = 0
    for j in 1:n
        pa, pb = cornerxyz[j], cornerxyz[mod1(j + 1, n)]
        ga, gb = cornergrid[j], cornergrid[mod1(j + 1, n)]
        besta = maximum(vdot(fc.c, ga) for fc in FACES)
        bestb = maximum(vdot(fc.c, gb) for fc in FACES)
        common = [f for f in 0:19 if besta - vdot(FACES[f+1].c, ga) <= 2e-12 &&
            bestb - vdot(FACES[f+1].c, gb) <= 2e-12]
        face = isempty(common) ? -1 : first(common)
        wa = face < 0 ? 0.0im : snyder_fwd_face(face, ga)
        wb = face < 0 ? 0.0im : snyder_fwd_face(face, gb)
        omega = face < 0 ? acos(clamp(vdot(pa, pb), -1.0, 1.0)) : 0.0
        for i in 0:(segments - 1)
            t = i / segments
            p = if face >= 0
                from_grid(sys.orientation, snyder_inv_xyz(face, wa + t * (wb - wa)))
            elseif omega > eps(Float64)
                so = sin(omega)
                vnormalize(vadd(vscale(pa, sin((1-t)*omega)/so),
                    vscale(pb, sin(t*omega)/so)))
            else
                pa
            end
            out[k += 1] = USPoint(p[1], p[2], p[3])
        end
    end
    return out
end

DGG.cell_boundary(sys::ISEA3HSystem, c::Z3Cell) = _hex_boundary(sys, c, 3)
DGG.cell_boundary(sys::ISEA4HSystem, c::DGG.LevelIndex) = _hex_boundary(sys, c, 4)

function equal_area_steradians(c::Z3Cell)
    _, path, r = _hexparts(c); h = 4pi / (10 * 3.0^r); return path == 0 ? 5h / 6 : h
end
function equal_area_steradians(sys::ISEA4HSystem, c::DGG.LevelIndex)
    _, path, r = _hexparts(c); h = 4pi / (10 * 4.0^r); return path == 0 ? 5h / 6 : h
end

function DGG.cell_area(g::LevelGrid{ISEA3HSystem}, c::Z3Cell)
    DGG.level(c) == g.level || throw(ArgumentError("cell and grid levels differ"))
    DGG.cellposition(g, c) === nothing && throw(ArgumentError("invalid ISEA3H cell"))
    return equal_area_steradians(c)
end
function DGG.cell_area(g::LevelGrid{ISEA4HSystem}, c::DGG.LevelIndex)
    DGG.level(c) == g.level || throw(ArgumentError("cell and grid levels differ"))
    DGG.cellposition(g, c) === nothing && throw(ArgumentError("invalid ISEA4H cell"))
    return equal_area_steradians(g.system, c)
end

# Bounded-candidate inversion of the central-place prefix construction.  At
# every digit it retains the closest small neighborhood in the unfolded chart.
# The neighborhood is much wider than the aperture (24 prefixes per incident
# root), which covers the non-nesting overhang while keeping work O(level) and
# independent of the total global cell count.
@inline function _remaining_radius(A::Int, k::Int, target::Int)
    k == target && return 0.0
    q = inv(sqrt(Float64(A)))
    first = L_PLANE * q^(k + 1)
    return first * (1 - q^(target - k)) / (1 - q)
end

@inline _pathcell(::ISEA3HSystem, root::Int, path::Int64, r::Int) = _z3from(root, path, r)
@inline _pathcell(::ISEA4HSystem, root::Int, path::Int64, r::Int) =
    DGG.LevelIndex(r, _i4encode(root, path, r))

# The five-face development is a 300-degree cone.  A point on its cut has two
# planar representatives related by a 60-degree turn.  Taking the minimum of
# the ordinary and the two adjacent representatives makes face-seam inversion
# deterministic without changing distances in the cone interior.
@inline function _dev_distance(x::ComplexF64, u::ComplexF64)
    d = abs(x - u)
    ax = mod(rad2deg(angle(x)), 360.0)
    au = mod(rad2deg(angle(u)), 360.0)
    # Only the two sides of the missing wedge are identified.  Rotating an
    # interior representative would incorrectly identify different centers.
    au >= 270.0 && ax <= 30.0 && (d = min(d, abs(x - u * cis(pi / 3))))
    ax >= 270.0 && au <= 30.0 && (d = min(d, abs(x - u * cis(-pi / 3))))
    return d
end

function _locate_in_root(A::Int, root::Int, x::ComplexF64, target::Int)
    frontier = Int64[0]
    for k in 1:target
        scored = Tuple{Float64,Int64}[]
        for prefix in frontier
            digits = _polar(root) ? (0:0) : (0:(A - 1))
            for d in digits
                path = prefix * A + d
                dist = _dev_distance(x, _hex_dev(A, root, path, k))
                push!(scored, (dist, path))
            end
        end
        sort!(scored)
        frontier = [scored[i][2] for i in 1:min(24, length(scored))]
    end
    bestpath = frontier[1]
    return bestpath, _dev_distance(x, _hex_dev(A, root, bestpath, target))
end

function _hex_cellat_face(g, q::NTuple{3,Float64}, A::Int)
    face, w = snyder_fwd(q)
    roots = FACES[face + 1].verts
    bestcell = nothing
    bestdist = Inf
    bestpos = typemax(Int)
    for root in roots
        x = face_to_dev(root, face, w)
        path, dist = _locate_in_root(A, root, x, g.level)
        c = _pathcell(g.system, root, path, g.level)
        pos = DGG.cellposition(g.system, c)
        if dist < bestdist - 32eps(L_PLANE) ||
                (abs(dist - bestdist) <= 32eps(L_PLANE) && pos < bestpos)
            bestcell = c; bestdist = dist; bestpos = pos
        end
    end
    return bestcell, bestdist
end

function _hex_cellat(g, p::GO.UnitSphericalPoint, A::Int)
    q = to_grid(g.system.orientation, (Float64(p[1]), Float64(p[2]), Float64(p[3])))
    bestface = maximum(vdot(fc.c, q) for fc in FACES)
    candidates = Tuple{Any,Float64}[]
    for fc in FACES
        bestface - vdot(fc.c, q) <= 2e-12 || continue
        # On an icosahedron edge snyder_fwd deliberately picks one face.  A
        # one-ulp-sized nudge into every tied face exposes the other valid
        # development without changing any Voronoi decision away from the tie.
        qq = vnormalize(vadd(q, vscale(fc.c, 2e-13)))
        push!(candidates, _hex_cellat_face(g, qq, A))
    end
    isempty(candidates) && push!(candidates, _hex_cellat_face(g, q, A))
    bestcell, bestdist = candidates[1]
    bestpos = DGG.cellposition(g.system, bestcell)
    for (c, dist) in candidates[2:end]
        pos = DGG.cellposition(g.system, c)
        if dist < bestdist - 2e-12 || (abs(dist - bestdist) <= 2e-12 && pos < bestpos)
            bestcell = c; bestdist = dist; bestpos = pos
        end
    end
    return bestcell
end
DGG.cellat(g::LevelGrid{ISEA3HSystem}, p::GO.UnitSphericalPoint) = _hex_cellat(g, p, 3)
DGG.cellat(g::LevelGrid{ISEA4HSystem}, p::GO.UnitSphericalPoint) = _hex_cellat(g, p, 4)

function _hex_neighbors1(g, c, A::Int)
    DGG.level(c) == g.level || throw(ArgumentError("cell and grid levels differ"))
    DGG.cellposition(g, c) === nothing && throw(ArgumentError("cell is not in the grid"))
    root, path, r = _hexparts(c)
    center = _hex_dev(A, root, path, r)
    step = L_PLANE / sqrt(Float64(A))^r
    rotation = A == 3 && isodd(r) ? 30.0 : 0.0
    centerang = center == 0 ? 0.0 : mod(rad2deg(angle(center)), 360.0)
    out = typeof(c)[]
    for j in 0:5
        u = center + step * cis(deg2rad(rotation + 60j))
        psi = mod(rad2deg(angle(u)), 360.0)
        if psi >= 300.0 - 1e-10
            u *= cis(deg2rad(centerang > 150.0 ? 60.0 : -60.0))
        end
        xyz = from_grid(g.system.orientation, dev_to_xyz(root, u))
        n = DGG.cellat(g, USPoint(xyz[1], xyz[2], xyz[3]))
        n == c && continue
        n in out || push!(out, n)
    end
    return out
end

function _hex_sort_ccw!(cells, g, subject)
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
    sort!(cells; by=c -> (mod(az(c) - zeroaz, 2pi), c))
    return cells
end

function _hex_shells(g, c, steps::Int, A::Int)
    steps >= 0 || throw(ArgumentError("k must be non-negative"))
    steps == 0 && return Vector{typeof(c)}[]
    seen = Set{typeof(c)}((c,))
    frontier = typeof(c)[c]
    shells = Vector{typeof(c)}[]
    for _ in 1:steps
        shell = typeof(c)[]
        for x in frontier, n in _hex_neighbors1(g, x, A)
            if !(n in seen)
                push!(seen, n); push!(shell, n)
            end
        end
        _hex_sort_ccw!(shell, g, c)
        push!(shells, shell); frontier = shell
    end
    return shells
end

function _hex_neighbors(g, c, steps::Int, A::Int)
    shells = _hex_shells(g, c, steps, A)
    return isempty(shells) ? typeof(c)[] : reduce(vcat, shells)
end

function DGG.neighbors(g::LevelGrid{ISEA3HSystem}, c::Z3Cell, k::Integer=1;
        connectivity::DGG.Connectivity=DGG.Vertex())
    return _hex_neighbors(g, c, Int(k), 3)
end
function DGG.neighbors(g::LevelGrid{ISEA4HSystem}, c::DGG.LevelIndex, k::Integer=1;
        connectivity::DGG.Connectivity=DGG.Vertex())
    return _hex_neighbors(g, c, Int(k), 4)
end
function DGG.ring(g::LevelGrid{ISEA3HSystem}, c::Z3Cell, k::Integer;
        connectivity::DGG.Connectivity=DGG.Vertex())
    steps = Int(k); steps >= 0 || throw(ArgumentError("k must be non-negative"))
    steps == 0 && return Z3Cell[c]
    return _hex_shells(g, c, steps, 3)[steps]
end
function DGG.ring(g::LevelGrid{ISEA4HSystem}, c::DGG.LevelIndex, k::Integer;
        connectivity::DGG.Connectivity=DGG.Vertex())
    steps = Int(k); steps >= 0 || throw(ArgumentError("k must be non-negative"))
    steps == 0 && return DGG.LevelIndex[c]
    return _hex_shells(g, c, steps, 4)[steps]
end
