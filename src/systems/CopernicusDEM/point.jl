# Point stencils between Copernicus DEM posts.
#
# A pixel is a post, so a point sample interpolates between the posts around it.
# Within a latitude band adjacent sample rows share one longitude lattice and the
# dual cell is a rectangle; at a band edge the two rows' cell edges interleave
# and it is a triangle or a trapezoid. Columns wrap, so a stencil across the
# antimeridian is an ordinary one.
#
# This file dispatches on `DGGSpace`, which `src/regridding.jl` defines after the
# systems are read, so `DiscreteGlobalGrids.jl` reads it into this module there
# rather than `CopernicusDEM.jl` reading it here.

import GlobalRegridding as GR

"""
The collections whose point stencils this lattice builds: a complete level grid,
and any holding of its cells.
"""
const PointGrid = Union{DGG.HierarchicalLevelGrid{<:CopernicusDEMSystem},
    DGG.PartialGrid{<:CopernicusDEMSystem}}

"""
    PointState{N}(ncols, rowbase)

The lattice tables a Copernicus DEM point query reads.

  - `ncols` gives a tile row's columns and `rowbase` the pixels above it, both
    the system's own tables, so a query is arithmetic and no table is built.
  - `N` is the latitude intervals per degree, so row and column arithmetic folds
    to constants.

It is read and never written, so one state serves concurrent queries.
"""
struct PointState{N}
    ncols::Vector{Int64}
    rowbase::Vector{Int64}
end

function _pointstate(sys::CopernicusDEMSystem{N}) where {N}
    t = tables(sys)
    return PointState{N}(t.ncols, t.rowbase)
end

Base.show(io::IO, ::PointState{N}) where {N} =
    print(io, "PointState{", N, "}()")

# ===========================================================================
# Row and column arithmetic on the complete level
# ===========================================================================

# Sample rows are the complete level's pixel rows, numbered north to south.
@inline _samplerows(::PointState{N}) where {N} = 180 * Int64(N)

@inline function _rowcols(st::PointState{N}, J::Int64) where {N}
    r = fld(J, Int64(N))
    return @inbounds st.ncols[Int(r)+1]
end

# The 0-based level-1 id of the cell at sample row `J` and column `K`, `K`
# already reduced modulo the row length.
@inline function _globalid(st::PointState{N}, J::Int64, K::Int64) where {N}
    NI = Int64(N)
    r = fld(J, NI)
    nc = @inbounds st.ncols[Int(r)+1]
    q, i = divrem(K, nc)
    return @inbounds st.rowbase[Int(r)+1] + q * nc * NI + (J - r * NI) * nc + i
end

# The post's latitude, formed exactly as `cell_centroid` forms it, so a sample
# site and the row it names agree bit for bit. The two pole rows are the clamped
# box's midpoint, north of and south of the posts they carry.
@inline function _sitelat(::PointState{N}, J::Int64) where {N}
    NI = Int64(N)
    r = fld(J, NI)
    j = Int(J - r * NI)
    lat_n = 90 - Int(r)
    half = (1 / N) / 2
    north = J == 0 ? 90.0 : (lat_n - j / N) + half
    south = J == 180 * NI - 1 ? -90.0 : (lat_n - (j + 1) / N) + half
    return (north + south) / 2
end

# The index of a node in the space being sampled, or `0` where its collection
# does not hold the node. The node arrives as its tile, its place `off` in that
# tile's raster, its id on the complete level and its tile's pixel count.
#
#   - A complete level is arithmetic: the id.
#   - A holding chunked by tile — the production source, whose chunks are the
#     listed tiles — is one search of the space's ascending tile ids and
#     arithmetic within the tile: a chunk holding every pixel of its tile keeps
#     them in id order, so the node's place in the chunk is its place in the
#     raster. A tile the space has no chunk for holds no pixel of it.
#   - Any other holding — a chunk holding part of its tile, or chunks at another
#     level — is one search of the holding's ascending pixel ids. On a source of
#     tens of billions of pixels named by a lazy id vector that search is the
#     costliest step of a query, which is what the tile route is for.
@inline _localof(::DGG.DGGSpace{<:DGG.HierarchicalLevelGrid{<:CopernicusDEMSystem}},
    ::Int64, ::Int64, id::Int64, ::Int64) = Int(id) + 1

@inline function _localof(space::DGG.DGGSpace{<:DGG.PartialGrid{<:CopernicusDEMSystem}},
    tile::Int64, off::Int64, id::Int64, ntile::Int64)
    if space.chunklevel == 0
        chunks = space.chunkids
        t = DGG.LevelIndex(0, tile)
        k = searchsortedfirst(chunks, t)
        (k <= length(chunks) && @inbounds(chunks[k]) == t) || return 0
        w = @inbounds space.ranges[k]
        length(w) == ntile && return first(w) + Int(off)
    end
    i = DGG.localindex(space.grid, DGG.LevelIndex(1, id))
    return i === nothing ? 0 : Int(i)
end

# The node at sample row `J` and column `K`, `K` not yet reduced, taken to its
# tile, its place in that tile's raster and its id — `_globalid` in parts.
@inline function _nodeindex(st::PointState{N}, space, J::Int64, K::Int64) where {N}
    NI = Int64(N)
    r = fld(J, NI)
    nc = @inbounds st.ncols[Int(r)+1]
    q, i = divrem(mod(K, 360 * nc), nc)
    off = (J - r * NI) * nc + i
    id = @inbounds(st.rowbase[Int(r)+1]) + q * nc * NI + off
    return _localof(space, Int64(tileordinal(r, q)), off, id, nc * NI)
end

# One stencil entry. A node carrying no weight is not required, so it is not
# looked up and its absence is not a rim.
@inline function _emit!(row::GR.WeightRow, st::PointState, space, J::Int64, K::Int64,
    w::Float64)
    w > 0 || return true
    i = _nodeindex(st, space, J, K)
    i == 0 && return false
    GR._addentry!(row, i, w)
    return true
end

# The policy for a point poleward of the outermost post row, `J` that row.
#
# `NearestCell` takes the nearest post of the row. Every post of a row stands at
# one latitude, so the nearest in longitude is the nearest on the sphere too:
# the cosine of the distance grows with the cosine of the longitude difference
# and with nothing else. The column is one rounding, and it wraps at the
# antimeridian like any other.
@inline function _polefill!(row::GR.WeightRow, st::PointState, space,
    ::GR.NearestCell, J::Int64, lon)
    x = Float64(lon) + 180.0
    x -= 360.0 * floor(x / 360.0)
    K = floor(Int64, x * _rowcols(st, J) + 0.5)
    i = _nodeindex(st, space, J, K)
    i == 0 && return GR.WeightsRim
    GR._addentry!(row, i, 1.0)
    return GR.WeightsMapped
end

@inline _polefill!(::GR.WeightRow, ::PointState, _grid, ::Nothing, ::Int64, _lon) =
    GR.WeightsDegenerate

# The four tensor products of a strip's two coordinates, north row first.
@inline function _emitquad!(row::GR.WeightRow, st::PointState, space, JA::Int64,
    K1::Int64, K2::Int64, K3::Int64, K4::Int64, u::Float64, v::Float64)
    ok = _emit!(row, st, space, JA, K1, (1 - u) * (1 - v)) &&
         _emit!(row, st, space, JA, K2, u * (1 - v)) &&
         _emit!(row, st, space, JA + 1, K3, u * v) &&
         _emit!(row, st, space, JA + 1, K4, (1 - u) * v)
    ok || (empty!(row); return GR.WeightsRim)
    return GR.WeightsMapped
end

# A trapezoid between the two sample rows, given its corner longitudes: `x1` and
# `x2` on the north row, `x4` and `x3` below them. Both parallel edges lie on a
# parallel, so the isoparametric map's second coordinate is the latitude
# fraction `v` and the first solves a linear equation along the parallel at `v`.
@inline function _q1strip!(row::GR.WeightRow, st::PointState, space, JA::Int64,
    K1::Int64, K2::Int64, K3::Int64, K4::Int64,
    x1::Float64, x2::Float64, x3::Float64, x4::Float64, v::Float64, x::Float64)
    xl = x1 + v * (x4 - x1)
    xr = x2 + v * (x3 - x2)
    d = xr - xl
    d > 0 || return GR.WeightsDegenerate
    u = GR._snapunit((x - xl) / d)
    (0.0 <= u <= 1.0) || return GR.WeightsOutside
    return _emitquad!(row, st, space, JA, K1, K2, K3, K4, u, v)
end

# Where the segment joining the post of column `KA` to the post of column `KB`
# crosses the parallel at latitude fraction `v`.
@inline _dualx(KA::Int64, KB::Int64, a::Int64, b::Int64, v::Float64) =
    (KA / a) + ((KB / b) - (KA / a)) * v

# The segments in a band-edge strip stand at most one cell of the coarser row
# apart, and the coarsest ratio between adjacent rows is two to one, so the walk
# to the one bounding a point is a handful of steps. Exceeding this is
# arithmetic that did not settle, not a wider strip.
const MAX_STRIP_WALK = 8

# A triangle of three posts, weighted by its barycentric coordinates, in a chart
# whose origin is the query point. The kernel is given the nodes' places in the
# ring and the places it keeps are exchanged for the collection's indices, so a
# node carrying no weight is never looked up.
@inline function _p1strip!(row::GR.WeightRow, st::PointState, space,
    Js::NTuple{3,Int64}, Ks::NTuple{3,Int64}, nodes::NTuple{3,NTuple{2,Float64}})
    status = GR.meanvalueweights!(row, (1, 2, 3), nodes, (0.0, 0.0))
    GR.ismapped(status) || return status
    for k in eachindex(row.indices)
        slot = @inbounds row.indices[k]
        i = _nodeindex(st, space, @inbounds(Js[slot]), @inbounds(Ks[slot]))
        i == 0 && (empty!(row); return GR.WeightsRim)
        @inbounds row.indices[k] = i
    end
    return GR.WeightsMapped
end

# ===========================================================================
# The sampler
# ===========================================================================

"""
    GlobalRegridding.samplerstate(space::DGGSpace{<:PointGrid})

The lattice tables a Copernicus DEM pixel source answers point queries from.

A pixel collection takes the fused row arithmetic below. A collection of tiles
has a lattice of its own and takes whatever the generic space offers.
"""
function GR.samplerstate(space::DGG.DGGSpace{<:PointGrid})
    grid = space.grid
    DGG.level(grid) == 1 || return invoke(GR.samplerstate, Tuple{DGG.DGGSpace}, space)
    return _pointstate(DGG.system(grid))
end

"""
    GlobalRegridding.weightsat!(row, sampler, p)

Weight the Copernicus DEM posts around `p`, at most four of them.

  - Between two rows of one latitude band the posts are a rectangle and the row
    is that rectangle's tensor Q1 weights.
  - At a band edge the rows' cell edges interleave: an edge on one row only
    carries a three-post triangle and its barycentric coordinates, an edge both
    rows share a four-post trapezoid and its Q1 coordinates. Which of the two
    bounds the point is decided in exact integer arithmetic.
  - Columns wrap, so a stencil across the antimeridian is an ordinary one, and a
    1-degree tile seam inside a band changes no stencil.
  - The stencil is built on the complete level and each node exchanged for its
    index in the collection sampled — by tile chunk and raster arithmetic where
    the space is chunked by tile, by a search of the holding's ids otherwise; a
    node the collection does not hold is `WeightsRim` and the row is empty.
  - A point poleward of the outermost post row has no dual cell short of that
    row entire, so the method's `poles` policy answers it: `NearestCell()`, the
    default, takes the nearest post of that row with weight one; `nothing`
    leaves it `WeightsDegenerate` with an empty row.
"""
function GR.weightsat!(row::GR.WeightRow,
    s::GR.Sampler{<:DGG.DGGSpace{<:PointGrid},<:AbstractVector,PointState{N}},
    p) where {N}
    empty!(row)
    st = s.state
    space = s.space
    lon, lat = FROM_SPHERE(p)
    (isfinite(lon) && isfinite(lat)) || return GR.WeightsOutside
    lat = clamp(Float64(lat), -90.0, 90.0)
    nrows = _samplerows(st)

    # The two post rows bracketing the point, and the fraction between them.
    JA = clamp(floor(Int64, (90.0 - lat) * N), Int64(0), nrows - 2)
    latA = _sitelat(st, JA)
    latB = _sitelat(st, JA + 1)
    if lat > latA
        JA == 0 && return _polefill!(row, st, space, s.method.poles, JA, lon)
        JA -= 1
        latB = latA
        latA = _sitelat(st, JA)
    elseif lat < latB
        JA + 1 == nrows - 1 &&
            return _polefill!(row, st, space, s.method.poles, JA + 1, lon)
        JA += 1
        latA = latB
        latB = _sitelat(st, JA + 1)
    end
    v = GR._snapunit((latA - lat) / (latA - latB))
    (0.0 <= v <= 1.0) || return GR.WeightsOutside

    # Longitude east of the antimeridian, in degrees, so a column of a row with
    # `nc` columns per degree stands at `K / nc`.
    x = Float64(lon) + 180.0
    x -= 360.0 * floor(x / 360.0)
    a = _rowcols(st, JA)
    b = _rowcols(st, JA + 1)

    if a == b
        u = x * a
        K = floor(Int64, u)
        u = GR._snapunit(u - K)
        return _emitquad!(row, st, space, JA, K, K + 1, K + 1, K, u, v)
    end

    # A band edge. The strip is divided by the segments joining each row's post
    # to the post facing it, one per stretch of the shared parallel that one cell
    # of each row covers. Walk east or west to the first such segment at or east
    # of the point; the dual cell is the one it and its western neighbour bound.
    KA = floor(Int64, x * a + 0.5)
    KB = floor(Int64, x * b + 0.5)
    settled = false
    for _ in 1:MAX_STRIP_WALK
        if x > _dualx(KA, KB, a, b, v)
            eA = (2 * KA + 1) * b
            eB = (2 * KB + 1) * a
            eA <= eB && (KA += 1)
            eB <= eA && (KB += 1)
            continue
        end
        wA = (2 * KA - 1) * b
        wB = (2 * KB - 1) * a
        qA = wA >= wB ? KA - 1 : KA
        qB = wB >= wA ? KB - 1 : KB
        if x < _dualx(qA, qB, a, b, v)
            KA, KB = qA, qB
            continue
        end
        settled = true
        break
    end
    settled || return GR.WeightsDegenerate

    # The cell's own vertex: an edge of one row alone carries a triangle, an edge
    # both rows share a trapezoid. Cell edges are compared in units of one over
    # twice the product of the two row lengths, so sharing is integer equality.
    wA = (2 * KA - 1) * b
    wB = (2 * KB - 1) * a
    xA = KA / a
    xB = KB / b
    wA == wB && return _q1strip!(row, st, space, JA, KA - 1, KA, KB, KB - 1,
        xA - 1 / a, xA, xB, xB - 1 / b, v, x)
    wA > wB && return _p1strip!(row, st, space,
        (JA, JA, JA + 1), (KA - 1, KA, KB),
        ((xA - 1 / a - x, latA - lat), (xA - x, latA - lat), (xB - x, latB - lat)))
    return _p1strip!(row, st, space,
        (JA, JA + 1, JA + 1), (KA, KB - 1, KB),
        ((xA - x, latA - lat), (xB - 1 / b - x, latB - lat), (xB - x, latB - lat)))
end

"""
    GlobalRegridding.supportradius(::BarycentricPoint, space::DGGSpace{<:PointGrid})

How far a Copernicus DEM point stencil reaches beyond its point, in radians.

A stencil takes the posts around the point: one row away in latitude, and at
most one column of the coarser row away in longitude. A column's east-west
distance is largest at its band's equatorward edge, so the bound is the largest
such diagonal over the six bands, doubled — a bound, not an estimate, and an
overestimate costs only discovery. A collection of tiles is left to the generic
answer.
"""
function GR.supportradius(m::GR.BarycentricPoint, space::DGG.DGGSpace{<:PointGrid})
    DGG.level(space.grid) == 1 || return invoke(GR.supportradius,
        Tuple{GR.BarycentricPoint,DGG.DGGSpace}, m, space)
    n = lat_intervals(DGG.system(space.grid))
    reach = 0.0
    for k in eachindex(BAND_EDGES)
        dlon = BAND_FACTOR2[k] / (2 * n)
        reach = max(reach, hypot(1 / n, dlon * cosd(BAND_EDGES[k])))
    end
    return deg2rad(2 * reach)
end
