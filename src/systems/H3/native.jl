# ---------------------------------------------------------------------------
# The libh3 ccall layer.
#
# Ported wholesale from the pre-redesign `src/H3/H3Native.jl`: these wrappers
# are the oracle, not a reimplementation of one, and they are deliberately
# unchanged. Seven wrappers were *added* on top of that file, none of them
# altering an existing one:
#
#   cell_center_cartesian   the centre without the degree round trip
#   boundary_verts          the boundary into a fixed tuple, no Vector
#   grid_ring_unsafe        the O(k) hollow-ring walk, for `ring`
#   grid_disk_distances     the pentagon-safe fallback `ring` needs
#   grid_disk_1             the k = 1 disk into a stack buffer
#   grid_ring_unsafe_1      the k = 1 hollow ring into a stack buffer
#   cell_to_children_7      one level of children into a stack buffer
#
# The last three exist so that `neighbors(grid, c)` and `children(sys, c)` —
# the two calls a whole-grid sweep makes per cell — do not allocate. libh3
# writes into a caller-supplied array, so a fixed-size `Ref` tuple that never
# escapes is stack-allocated by Julia and the heap is never touched.
#
# Nothing in here knows about the grid interface. It speaks H3's own vocabulary
# — raw `UInt64` indices, degrees, zero-based child positions — and the
# interface wiring lives one directory up in `system.jl` / `geometry.jl`.
# ---------------------------------------------------------------------------

"""
    DiscreteGlobalGrids.H3.H3Native

The raw libh3 ccall layer: a couple of dozen entry points wrapping ~28 `ccall`s
into the C library shipped by
[`H3_jll`](https://github.com/JuliaBinaryWrappers/H3_jll.jl), and the whole of
this package's dependence on it.

It speaks **H3's own vocabulary** — raw `UInt64` indices, degrees, zero-based
child positions — and knows nothing about the grid interface. The interface
wiring lives one directory up, in `system.jl` / `geometry.jl` / `neighbors.jl` /
`border.jl`.

These wrappers are **the oracle, not a reimplementation of one**. They were
ported wholesale from the pre-redesign `src/H3/H3Native.jl` and are deliberately
left unchanged, so that every geometric, hierarchical and adjacency answer this
system gives is libh3's own answer rather than a rewrite that could drift from
it. The only additions are the seven allocation-free wrappers listed in the file
header, none of which alters an existing one.

Callers normally want the interface generics instead — [`cell_boundary`](@ref),
[`cellat`](@ref), [`neighbors`](@ref), [`children`](@ref) and friends on
[`H3System`](@ref) and its [`levelgrid`](@ref) — which is where the canonical
ordering, the dense position numbering and the subtree border walk live. Reach for
`H3Native` only when you specifically want to speak to libh3 directly.
"""
module H3Native

using ...Helpers: to_uint64_id
using H3_jll

export MAX_RESOLUTION,
    cell_area,
    cell_boundary,
    cell_boundary_cartesian,
    cell_center,
    cell_to_child_pos,
    cell_to_children,
    cell_to_children_7,
    cell_to_children_size,
    cell_to_parent,
    child_pos_to_cell,
    get_base_cell,
    get_pentagons,
    get_resolution,
    grid_disk,
    grid_disk_1,
    grid_disk_distances,
    grid_ring_unsafe,
    grid_ring_unsafe_1,
    is_pentagon,
    is_valid_cell,
    lonlat_to_cell,
    max_grid_disk_size,
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

"""
    cell_center_cartesian(id) -> NTuple{3,Float64}

The cell centre straight onto the unit sphere, skipping the degree round trip
`cell_center` makes. This is what [`cell_centroid`](@ref) is built on.
"""
function cell_center_cartesian(id)
    out = Ref{LatLng}()
    _check(ccall((:cellToLatLng, H3_jll.libh3), H3Error,
                 (H3Index, Ref{LatLng}), _to_id(id), out),
           "cellToLatLng")
    λ = out[].lng
    φ = out[].lat
    cφ = cos(φ)
    return (cφ * cos(λ), cφ * sin(λ), sin(φ))
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

"""
    boundary_verts(id) -> (NTuple{10,NTuple{3,Float64}}, Int)

The cell's boundary as unit-sphere `xyz`, in a fixed-size tuple plus the count
of vertices actually used. Allocation-free, which is what
[`cell_boundary`](@ref) wants: an H3 boundary has between 5 and 10 vertices
(the extra ones are the distortion vertices where a cell crosses an
icosahedron face edge), so a `NTuple{10}` is exactly the `CellBoundary` struct
libh3 already filled.
"""
function boundary_verts(id)
    b = _boundary(id)
    n = Int(b.numVerts)
    verts = ntuple(10) do i
        v = b.verts[i]
        cφ = cos(v.lat)
        (cφ * cos(v.lng), cφ * sin(v.lng), sin(v.lat))
    end
    return verts, n
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

"""Number of index slots a grid disk of radius `k` needs (`3k(k+1) + 1`)."""
function max_grid_disk_size(k::Integer)
    k >= 0 || throw(ArgumentError("grid disk radius must be non-negative"))
    out = Ref{Int64}(0)
    _check(ccall((:maxGridDiskSize, H3_jll.libh3), H3Error,
                 (Int64, Ref{Int64}), Int64(k), out),
           "maxGridDiskSize")
    return out[]
end

"""
    grid_disk(id, k) -> Vector{H3Index}

All cells within grid distance `k` of `id`, `id` itself included. The vector
has `max_grid_disk_size(k)` slots and keeps `0` where the disk is truncated
by a pentagon distortion; no cell order is promised. This is `gridDisk`, the
pentagon-safe path — `gridRingUnsafe` fails whenever a pentagon sits inside
the ring.
"""
function grid_disk(id, k::Integer)
    cell = _to_id(id)
    disk = zeros(H3Index, max_grid_disk_size(k))
    _check(ccall((:gridDisk, H3_jll.libh3), H3Error,
                 (H3Index, Cint, Ptr{H3Index}), cell, Cint(k), disk),
           "gridDisk")
    return disk
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

# ---------------------------------------------------------------------------
# Additions for the redesigned `ring` / `neighbors` fast paths.
#
# `gridRingUnsafe` is the O(k) hollow-ring walk, and it is the reason `ring`
# does not have to build and subtract two discs. It refuses — with a non-zero
# error rather than a wrong answer — whenever the walk meets a pentagon, so it
# is wrapped to return `nothing` there and let the caller fall back to the
# pentagon-safe `gridDisk` path.
# ---------------------------------------------------------------------------

"""
    grid_ring_unsafe(id, k) -> Union{Vector{H3Index},Nothing}

The cells at grid distance **exactly** `k` from `id` (`6k` of them, or `id`
itself for `k == 0`), or `nothing` when libh3 reports that the walk met a
pentagon and cannot be completed.

`nothing` is not an error: pentagon distortion genuinely breaks the hollow-ring
walk, and the caller answers those cells from `grid_disk` instead. This is the
whole reason the libh3 function is named "unsafe".
"""
function grid_ring_unsafe(id, k::Integer)
    k >= 0 || throw(ArgumentError("grid ring radius must be non-negative"))
    n = k == 0 ? 1 : 6 * Int(k)
    out = zeros(H3Index, n)
    err = ccall((:gridRingUnsafe, H3_jll.libh3), H3Error,
                (H3Index, Cint, Ptr{H3Index}), _to_id(id), Cint(k), out)
    err == 0 || return nothing
    return out
end

# Seven slots covers both stack-buffer calls below: `maxGridDiskSize(1)` is 7,
# and a cell has at most seven children.
const _ZERO7 = ntuple(_ -> H3Index(0), Val(7))
const _ZERO6 = ntuple(_ -> H3Index(0), Val(6))

"""
    grid_ring_unsafe_1(id) -> Union{NTuple{6,H3Index},Nothing}

The `k = 1` hollow ring — the six neighbours in counter-clockwise order —
**without allocating**, or `nothing` when the walk meets a pentagon.

The stack-buffer form of [`grid_ring_unsafe`](@ref): a `k = 1` ring is always
exactly `6k = 6` slots, so no size call is needed and the `Ref` tuple never
escapes.
"""
function grid_ring_unsafe_1(id)
    buf = Ref(_ZERO6)
    err = GC.@preserve buf ccall((:gridRingUnsafe, H3_jll.libh3), H3Error,
        (H3Index, Cint, Ptr{H3Index}),
        _to_id(id), Cint(1), Base.unsafe_convert(Ptr{H3Index}, buf))
    err == 0 || return nothing
    return buf[]
end

"""
    grid_disk_1(id) -> NTuple{7,H3Index}

The `k = 1` grid disk — `id` and its neighbours — **without allocating**.

`maxGridDiskSize(1)` is always 7, so the size call the general `grid_disk`
makes is skipped and libh3 writes straight into a seven-slot `Ref` tuple. The
`Ref` never escapes, so Julia stack-allocates it and this touches no heap at
all; that is what lets a whole-grid neighbour sweep run garbage-free.

Unused slots — the six neighbours of a pentagon leave one — stay `0`, which is
never a valid H3 index, so the caller filters on it.
"""
function grid_disk_1(id)
    buf = Ref(_ZERO7)
    err = GC.@preserve buf ccall((:gridDisk, H3_jll.libh3), H3Error,
        (H3Index, Cint, Ptr{H3Index}),
        _to_id(id), Cint(1), Base.unsafe_convert(Ptr{H3Index}, buf))
    _check(err, "gridDisk")
    return buf[]
end

"""
    cell_to_children_7(id, child_resolution) -> NTuple{7,H3Index}

One level of children — seven for a hexagon, six for a pentagon — **without
allocating**, by the same stack-buffer route as [`grid_disk_1`](@ref).

`child_resolution` must be exactly one finer than `id`'s, which is the only
case that fits in seven slots. Unused slots stay `0`.
"""
function cell_to_children_7(id, child_resolution::Integer)
    res = _check_resolution(child_resolution)
    buf = Ref(_ZERO7)
    err = GC.@preserve buf ccall((:cellToChildren, H3_jll.libh3), H3Error,
        (H3Index, Cint, Ptr{H3Index}),
        _to_id(id), Cint(res), Base.unsafe_convert(Ptr{H3Index}, buf))
    _check(err, "cellToChildren")
    return buf[]
end

"""
    grid_disk_distances(id, k) -> (Vector{H3Index}, Vector{Cint})

`grid_disk` with each cell's grid distance from `id` alongside it. Slots the
disk does not fill hold `0` and a distance of `-1`; the pentagon-safe path,
like `grid_disk` itself.
"""
function grid_disk_distances(id, k::Integer)
    cell = _to_id(id)
    n = max_grid_disk_size(k)
    cells = zeros(H3Index, n)
    dists = fill(Cint(-1), n)
    _check(ccall((:gridDiskDistances, H3_jll.libh3), H3Error,
                 (H3Index, Cint, Ptr{H3Index}, Ptr{Cint}), cell, Cint(k), cells, dists),
           "gridDiskDistances")
    return cells, dists
end

end # module H3Native
