# Root cut topology is identified from the analytic projection at coarse scale,
# then applied to every depth with exact ternary-lattice arithmetic.  This
# avoids a placement-specific seam table while remaining deterministic for all
# 16 polar-square arrangements.

@inline function _normalize_point(x, y, z)
    n = sqrt(x * x + y * y + z * z)
    return GO.UnitSphericalPoint(x / n, y / n, z / n)
end

function _outside_point(boundary::GO.UnitSphericalPoint,
        inside::GO.UnitSphericalPoint)
    dot = boundary[1] * inside[1] + boundary[2] * inside[2] +
          boundary[3] * inside[3]
    tx = inside[1] - dot * boundary[1]
    ty = inside[2] - dot * boundary[2]
    tz = inside[3] - dot * boundary[3]
    norm = sqrt(tx * tx + ty * ty + tz * tz)
    norm > 16eps(Float64) || return _normalize_point(
        2boundary[1] - inside[1], 2boundary[2] - inside[2],
        2boundary[3] - inside[3])
    tx /= norm; ty /= norm; tz /= norm
    rho = acos(clamp(dot, -1.0, 1.0))
    s, c = sincos(rho)
    return _normalize_point(c * boundary[1] - s * tx,
        c * boundary[2] - s * ty, c * boundary[3] - s * tz)
end

@inline _nside3(level::Integer) = Int64(3)^Int(level)

function _decode_lattice(c::RHEALPixCell)
    l = DGG.level(c)
    root, rest = divrem(c.ordinal, _pow9(l))
    row = col = Int64(0)
    divisor = _pow9(l)
    for _ in 1:l
        divisor ÷= 9
        digit, rest = divrem(rest, divisor)
        row = 3row + digit ÷ 3
        col = 3col + digit % 3
    end
    return Int(root), row, col
end

function _encode_lattice(level::Int, root::Int, row::Int64, col::Int64)
    divisor = _nside3(level)
    ordinal = Int64(root)
    r, c = row, col
    for _ in 1:level
        divisor ÷= 3
        rd, r = divrem(r, divisor)
        cd, c = divrem(c, divisor)
        ordinal = 9ordinal + 3rd + cd
    end
    return RHEALPixCell(level, ordinal)
end

# Derive a root cut's target square, target edge, and orientation at coarse
# scale, where floating point has enormous headroom.  The resulting metadata
# is then applied with integer ternary lattice arithmetic at every depth.
# Direction slots are projected up, left, down, right.
function _seam_metadata(sys::RHEALPixSystem, root::Int, direction::Int)
    ulx, uly = _root_ul(sys, root)
    delta = HALFPI / 32
    function across(t)
        if direction == 1
            bx, by = ulx + t * HALFPI, uly
            ix, iy = bx, by - delta
        elseif direction == 2
            bx, by = ulx, uly - t * HALFPI
            ix, iy = bx + delta, by
        elseif direction == 3
            bx, by = ulx + t * HALFPI, uly - HALFPI
            ix, iy = bx, by + delta
        else
            bx, by = ulx + HALFPI, uly - t * HALFPI
            ix, iy = bx - delta, by
        end
        boundary = _sphere_point(sys, bx, by)
        inside = _sphere_point(sys, ix, iy)
        outside = _outside_point(boundary, inside)
        lon = atan(outside[2], outside[1])
        lat = asin(clamp(outside[3], -1.0, 1.0))
        x, y = _project(sys, lon, lat)
        target = _root_for_plane(x, y)
        tx, ty = _root_ul(sys, target)
        return target, (x - tx) / HALFPI, (ty - y) / HALFPI
    end
    a = across(0.25)
    b = across(0.75)
    a[1] == b[1] || error("rHEALPix seam crosses two roots")
    distances = (abs(a[3]), abs(1 - a[2]), abs(1 - a[3]), abs(a[2]))
    target_edge = argmin(distances) # top, right, bottom, left
    avary = target_edge in (1, 3) ? a[2] : a[3]
    bvary = target_edge in (1, 3) ? b[2] : b[3]
    return a[1], target_edge, bvary < avary
end

"""
    edge_neighbors(grid, cell) -> NTuple{4,RHEALPixCell}

The four edge neighbours in counter-clockwise order from projected up:
`(up, left, down, right)`.  Interior moves and deep seam crossings use exact
integer ternary-lattice arithmetic.  Root seam orientation is derived from the
analytic projection at a fixed coarse scale, making this valid for every polar
square placement without a hand-maintained 16-layout table.
"""
function edge_neighbors(g::LevelGrid, c::RHEALPixCell)
    _checked_index(g, c)
    l = g.level
    root, row, col = _decode_lattice(c)
    side = _nside3(l)
    return ntuple(4) do direction
        # Interior move in up, left, down, right order.
        direction == 1 && row > 0 && return _encode_lattice(l, root, row - 1, col)
        direction == 2 && col > 0 && return _encode_lattice(l, root, row, col - 1)
        direction == 3 && row + 1 < side && return _encode_lattice(l, root, row + 1, col)
        direction == 4 && col + 1 < side && return _encode_lattice(l, root, row, col + 1)

        target, target_edge, reverse = _seam_metadata(g.system, root, direction)
        varying = direction in (1, 3) ? col : row
        reverse && (varying = side - 1 - varying)
        if target_edge == 1
            tr, tc = Int64(0), varying
        elseif target_edge == 2
            tr, tc = varying, side - 1
        elseif target_edge == 3
            tr, tc = side - 1, varying
        else
            tr, tc = varying, Int64(0)
        end
        return _encode_lattice(l, target, tr, tc)
    end
end

@inline function _tangent_frame(point)
    axis = abs(point[3]) < 0.8 ? (0.0, 0.0, 1.0) : (1.0, 0.0, 0.0)
    d = axis[1] * point[1] + axis[2] * point[2] + axis[3] * point[3]
    x = axis[1] - d * point[1]
    y = axis[2] - d * point[2]
    z = axis[3] - d * point[3]
    n = sqrt(x * x + y * y + z * z)
    e1 = (x / n, y / n, z / n)
    e2 = (point[2] * e1[3] - point[3] * e1[2],
          point[3] * e1[1] - point[1] * e1[3],
          point[1] * e1[2] - point[2] * e1[1])
    return e1, e2
end

@inline function _azimuth_about(center, e1, e2, point)
    dot = center[1] * point[1] + center[2] * point[2] + center[3] * point[3]
    x = point[1] - dot * center[1]
    y = point[2] - dot * center[2]
    z = point[3] - dot * center[3]
    return atan(x * e2[1] + y * e2[2] + z * e2[3],
                x * e1[1] + y * e1[2] + z * e1[3])
end

function _sort_ccw!(cells::Vector{RHEALPixCell}, g::LevelGrid,
        c::RHEALPixCell)
    length(cells) <= 1 && return cells
    ulx, uly, width = cell_rectangle(g.system, c)
    center = DGG.cell_centroid(g, c)
    upward = _sphere_point(g.system, ulx + width / 2, uly - width / 4)
    d = center[1] * upward[1] + center[2] * upward[2] + center[3] * upward[3]
    tx = upward[1] - d * center[1]
    ty = upward[2] - d * center[2]
    tz = upward[3] - d * center[3]
    n = sqrt(tx * tx + ty * ty + tz * tz)
    if n <= 16eps(Float64)
        e1, e2 = _tangent_frame(center)
    else
        e1 = (tx / n, ty / n, tz / n)
        e2 = (center[2] * e1[3] - center[3] * e1[2],
              center[3] * e1[1] - center[1] * e1[3],
              center[1] * e1[2] - center[2] * e1[1])
    end
    sort!(cells; by = cell ->
        (mod(_azimuth_about(center, e1, e2, DGG.cell_centroid(g, cell)),
             2 * Float64(pi)), cell))
    return cells
end

function vertex_neighbors(g::LevelGrid, c::RHEALPixCell)
    _checked_index(g, c)
    edges = collect(edge_neighbors(g, c))
    # At level zero each square's four corners collapse pairwise through the
    # projection cuts; the four edge neighbours are the complete vertex star.
    g.level == 0 && return _sort_ccw!(edges, g, c)
    found = copy(edges)
    # At each corner, the diagonal (when the vertex is 4-valent) is the common
    # edge neighbour of the two incident edge neighbours.  At a 3-valent dart
    # vertex the incident pair meet each other and contributes no new cell.
    for i in 1:4
        a = edges[i]
        b = edges[mod1(i + 1, 4)]
        around_a = edge_neighbors(g, a)
        around_b = edge_neighbors(g, b)
        for candidate in around_a
            candidate == c && continue
            candidate in around_b || continue
            candidate in found || push!(found, candidate)
        end
    end
    length(found) <= 8 || error(
        "rHEALPix vertex-star construction found $(length(found)) neighbours for $(suid(c))")
    return _sort_ccw!(found, g, c)
end

function _one_ring(g::LevelGrid, c::RHEALPixCell, ::DGG.Edge)
    result = SmallVector{4,RHEALPixCell}()
    for neighbor in edge_neighbors(g, c)
        result = SmallCollections.push(result, neighbor)
    end
    return result
end

function _one_ring(g::LevelGrid, c::RHEALPixCell, ::DGG.Vertex)
    result = SmallVector{8,RHEALPixCell}()
    for neighbor in vertex_neighbors(g, c)
        result = SmallCollections.push(result, neighbor)
    end
    return result
end

function DGG.neighbors(g::LevelGrid, c::RHEALPixCell, k::Integer=1;
        connectivity::DGG.Connectivity=DGG.Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative"))
    steps == 0 && return SmallVector{8,RHEALPixCell}()
    steps == 1 && return _one_ring(g, c, connectivity)
    shells = _shells(g, c, steps, connectivity)
    return isempty(shells) ? RHEALPixCell[] : reduce(vcat, shells)
end

function DGG.ring(g::LevelGrid, c::RHEALPixCell, k::Integer;
        connectivity::DGG.Connectivity=DGG.Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative"))
    steps == 0 && return RHEALPixCell[c]
    shells = _shells(g, c, steps, connectivity)
    return steps <= length(shells) ? shells[steps] : RHEALPixCell[]
end

function _shells(g::LevelGrid, c::RHEALPixCell, steps::Int,
        connectivity::DGG.Connectivity)
    shells = Vector{RHEALPixCell}[]
    seen = Set{RHEALPixCell}((c,))
    frontier = RHEALPixCell[c]
    for depth in 1:steps
        next = RHEALPixCell[]
        for cell in frontier, candidate in _one_ring(g, cell, connectivity)
            candidate in seen && continue
            push!(seen, candidate)
            push!(next, candidate)
        end
        depth > 1 && _sort_ccw!(next, g, c)
        push!(shells, next)
        isempty(next) && break
        frontier = next
    end
    return shells
end
