module H3Native

using ..Helpers: to_uint64_id
using H3_jll

export MAX_RESOLUTION,
    cell_area,
    cell_boundary,
    cell_boundary_cartesian,
    cell_center,
    cell_to_child_pos,
    cell_to_children,
    cell_to_children_size,
    cell_to_parent,
    child_pos_to_cell,
    get_base_cell,
    get_pentagons,
    get_resolution,
    is_pentagon,
    is_valid_cell,
    lonlat_to_cell,
    num_cells,
    res0_cells

const MAX_RESOLUTION = 15
const H3Index = UInt64
const H3Error = Cint

struct LatLng
    lat::Cdouble
    lng::Cdouble
end

struct CellBoundary
    numVerts::Cint
    verts::NTuple{10,LatLng}
end

function _describe_error(err::Integer)
    msg = ccall((:describeH3Error, H3_jll.libh3), Cstring, (H3Error,), H3Error(err))
    return unsafe_string(msg)
end

function _check(err::Integer, fn::AbstractString)
    err == 0 && return nothing
    throw(ErrorException("$fn failed with H3 error $err: $(_describe_error(err))"))
end

function _check_resolution(resolution::Integer)
    res = Int(resolution)
    0 <= res <= MAX_RESOLUTION ||
        throw(ArgumentError("H3 resolution must be in 0:$MAX_RESOLUTION"))
    return res
end

const _to_id = to_uint64_id

function lonlat_to_cell(lon::Real, lat::Real, resolution::Integer)
    res = _check_resolution(resolution)
    coord = LatLng(deg2rad(Float64(lat)), deg2rad(Float64(lon)))
    out = Ref{H3Index}(0)
    _check(ccall((:latLngToCell, H3_jll.libh3), H3Error,
                 (Ref{LatLng}, Cint, Ref{H3Index}), coord, res, out),
           "latLngToCell")
    return out[]
end

function cell_center(id)
    cell = _to_id(id)
    out = Ref{LatLng}()
    _check(ccall((:cellToLatLng, H3_jll.libh3), H3Error,
                 (H3Index, Ref{LatLng}), cell, out),
           "cellToLatLng")
    return (_wraplon(rad2deg(out[].lng)), rad2deg(out[].lat))
end

function _boundary(id)
    out = Ref{CellBoundary}()
    _check(ccall((:cellToBoundary, H3_jll.libh3), H3Error,
                 (H3Index, Ref{CellBoundary}), _to_id(id), out),
           "cellToBoundary")
    return out[]
end

function _unwrap_lon(lon::Float64, center_lon::Float64)
    return lon - 360 * round((lon - center_lon) / 360)
end

function cell_boundary(id; closed_ring::Bool=true)
    boundary = _boundary(id)
    center_lon, _ = cell_center(id)
    npoints = Int(boundary.numVerts)
    points = Vector{Tuple{Float64,Float64}}(undef, npoints + Int(closed_ring))
    @inbounds for i in 1:npoints
        vertex = boundary.verts[i]
        lon = _unwrap_lon(rad2deg(vertex.lng), center_lon)
        lat = rad2deg(vertex.lat)
        points[i] = (lon, lat)
    end
    if closed_ring
        points[end] = points[1]
    end
    return points
end

function cell_boundary_cartesian(id; closed_ring::Bool=true)
    boundary = _boundary(id)
    npoints = Int(boundary.numVerts)
    points = Vector{NTuple{3,Float64}}(undef, npoints + Int(closed_ring))
    @inbounds for i in 1:npoints
        vertex = boundary.verts[i]
        λ = vertex.lng
        φ = vertex.lat
        cφ = cos(φ)
        points[i] = (cφ * cos(λ), cφ * sin(λ), sin(φ))
    end
    if closed_ring
        points[end] = points[1]
    end
    return points
end

function get_resolution(id)
    return Int(ccall((:getResolution, H3_jll.libh3), Cint, (H3Index,), _to_id(id)))
end

function get_base_cell(id)
    return Int(ccall((:getBaseCellNumber, H3_jll.libh3), Cint, (H3Index,), _to_id(id)))
end

function is_valid_cell(id)
    return ccall((:isValidCell, H3_jll.libh3), Cint, (H3Index,), _to_id(id)) != 0
end

function is_pentagon(id)
    return ccall((:isPentagon, H3_jll.libh3), Cint, (H3Index,), _to_id(id)) != 0
end

function cell_to_parent(id, resolution::Integer)
    res = _check_resolution(resolution)
    out = Ref{H3Index}(0)
    _check(ccall((:cellToParent, H3_jll.libh3), H3Error,
                 (H3Index, Cint, Ref{H3Index}), _to_id(id), res, out),
           "cellToParent")
    return out[]
end

function cell_to_children_size(id, resolution::Integer)
    res = _check_resolution(resolution)
    out = Ref{Int64}(0)
    _check(ccall((:cellToChildrenSize, H3_jll.libh3), H3Error,
                 (H3Index, Cint, Ref{Int64}), _to_id(id), res, out),
           "cellToChildrenSize")
    return out[]
end

function cell_to_children(id, resolution::Union{Nothing,Integer}=nothing)
    cell = _to_id(id)
    child_res = isnothing(resolution) ? get_resolution(cell) + 1 : _check_resolution(resolution)
    n = cell_to_children_size(cell, child_res)
    children = Vector{H3Index}(undef, n)
    _check(ccall((:cellToChildren, H3_jll.libh3), H3Error,
                 (H3Index, Cint, Ptr{H3Index}), cell, child_res, children),
           "cellToChildren")
    return children
end

"""Return the zero-based position of `id` among its ancestor's ordered children."""
function cell_to_child_pos(id, parent_resolution::Integer)
    parent_res = _check_resolution(parent_resolution)
    out = Ref{Int64}(0)
    _check(ccall((:cellToChildPos, H3_jll.libh3), H3Error,
                 (H3Index, Cint, Ref{Int64}), _to_id(id), parent_res, out),
           "cellToChildPos")
    return out[]
end

"""Return the child at a zero-based `position`, in `cell_to_children` order."""
function child_pos_to_cell(position::Integer, id, resolution::Integer)
    position >= 0 || throw(ArgumentError("child position must be non-negative"))
    position <= typemax(Int64) || throw(OverflowError("child position does not fit in Int64"))
    child_res = _check_resolution(resolution)
    out = Ref{H3Index}(0)
    _check(ccall((:childPosToCell, H3_jll.libh3), H3Error,
                 (Int64, H3Index, Cint, Ref{H3Index}),
                 Int64(position), _to_id(id), child_res, out),
           "childPosToCell")
    return out[]
end

function num_cells(resolution::Integer)
    res = _check_resolution(resolution)
    out = Ref{Int64}(0)
    _check(ccall((:getNumCells, H3_jll.libh3), H3Error,
                 (Cint, Ref{Int64}), res, out),
           "getNumCells")
    return out[]
end

function res0_cells()
    n = Int(ccall((:res0CellCount, H3_jll.libh3), Cint, ()))
    cells = Vector{H3Index}(undef, n)
    _check(ccall((:getRes0Cells, H3_jll.libh3), H3Error, (Ptr{H3Index},), cells),
           "getRes0Cells")
    return cells
end

function get_pentagons(resolution::Integer)
    res = _check_resolution(resolution)
    n = Int(ccall((:pentagonCount, H3_jll.libh3), Cint, ()))
    cells = Vector{H3Index}(undef, n)
    _check(ccall((:getPentagons, H3_jll.libh3), H3Error,
                 (Cint, Ptr{H3Index}), res, cells),
           "getPentagons")
    return cells
end

function cell_area(id)
    out = Ref{Cdouble}(0)
    _check(ccall((:cellAreaRads2, H3_jll.libh3), H3Error,
                 (H3Index, Ref{Cdouble}), _to_id(id), out),
           "cellAreaRads2")
    return out[]
end

_wraplon(lon) = lon > 180 ? lon - 360 : lon

end
