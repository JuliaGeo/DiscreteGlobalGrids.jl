# Two copies of the audit's OctantSystem toy (audit/dev-experience/toygrid.jl):
# one that DECLARES `maxneighbors`, one that does not. Everything else — the
# geometry, the hierarchy, the id scheme — is shared, so a declared/undeclared
# comparison isolates exactly the neighbour-buffer choice.
#
# OctantSystem: aperture-4, the 8 octants of the sphere refined by spherical
# midpoint subdivision. Canonical id LevelIndex(level, index), zero-based,
# index = octant * 4^level + path.

import DiscreteGlobalGrids as DGG
import GeometryOps as GO
const USP = GO.UnitSpherical.UnitSphericalPoint{Float64}

struct OctantDeclared <: DGG.AbstractHierarchicalGridSystem end
struct OctantBare <: DGG.AbstractHierarchicalGridSystem end

const AnyOctant = Union{OctantDeclared,OctantBare}

# --- identity & levels ------------------------------------------------------

DGG.cellindextype(::AnyOctant) = DGG.LevelIndex
DGG.levels(::AnyOctant) = 0:6
DGG.rootcells(::OctantDeclared) = [DGG.LevelIndex(0, i) for i in 0:7]
DGG.rootcells(::OctantBare) = [DGG.LevelIndex(0, i) for i in 0:7]

# --- hierarchy --------------------------------------------------------------

function Base.parent(sys::AnyOctant, c::DGG.LevelIndex)
    l = DGG.level(c)
    l == first(DGG.levels(sys)) &&
        throw(ArgumentError("root cells have no parent"))
    DGG.LevelIndex(l - 1, DGG.rawid(c) >> 2)
end

function DGG.children(sys::AnyOctant, c::DGG.LevelIndex)
    l = DGG.level(c)
    l == DGG.maxlevel(sys) &&
        throw(ArgumentError("cells at maxlevel have no children"))
    i4 = DGG.rawid(c) << 2
    [DGG.LevelIndex(l + 1, i4 + d) for d in 0:3]
end

# --- the five level-grid primitives -----------------------------------------

DGG.ncells(::AnyOctant, l::Integer) = 8 * 4^Int(l)

DGG.cellindex(::AnyOctant, l::Integer, i::Int) = DGG.LevelIndex(l, i - 1)

function DGG.cellposition(sys::AnyOctant, c::DGG.LevelIndex)
    r = DGG.rawid(c)
    0 <= r < DGG.ncells(sys, DGG.level(c)) || return nothing
    Int(r) + 1
end

_norm3(p) = (n = sqrt(p[1]^2 + p[2]^2 + p[3]^2); (p[1] / n, p[2] / n, p[3] / n))
_mid(a, b) = _norm3((a[1] + b[1], a[2] + b[2], a[3] + b[3]))

function _octant_triangle(r::Integer)
    sx = (r >> 2) & 1 == 0 ? 1.0 : -1.0
    sy = (r >> 1) & 1 == 0 ? 1.0 : -1.0
    sz = r & 1 == 0 ? 1.0 : -1.0
    X = (sx, 0.0, 0.0); Y = (0.0, sy, 0.0); Z = (0.0, 0.0, sz)
    sx * sy * sz > 0 ? (X, Y, Z) : (X, Z, Y)
end

function _triangle(c::DGG.LevelIndex)
    l = DGG.level(c)
    idx = DGG.rawid(c)
    r = idx >> (2l)
    p = idx - (r << (2l))
    (A, B, C) = _octant_triangle(r)
    for k in (l-1):-1:0
        d = (p >> (2k)) & 3
        mab = _mid(A, B); mbc = _mid(B, C); mca = _mid(C, A)
        if d == 0
            B, C = mab, mca
        elseif d == 1
            A, C = mab, mbc
        elseif d == 2
            A, B = mca, mbc
        else
            A, B, C = mab, mbc, mca
        end
    end
    (A, B, C)
end

function DGG.cell_boundary(::AnyOctant, c::DGG.LevelIndex)
    (A, B, C) = _triangle(c)
    [USP(A...), USP(B...), USP(C...)]
end

function DGG.cell_centroid(::AnyOctant, c::DGG.LevelIndex)
    (A, B, C) = _triangle(c)
    USP(_norm3((A[1] + B[1] + C[1], A[2] + B[2] + C[2], A[3] + B[3] + C[3]))...)
end

# --- traits -----------------------------------------------------------------

# THE difference between the two systems: OctantBare declares no bound.
DGG.maxneighbors(::OctantDeclared, ::DGG.Vertex) = 12
DGG.maxneighbors(::OctantDeclared, ::DGG.Edge) = 3

DGG.has_sorted_subtrees(::AnyOctant) = true

function DGG.descendant_range(sys::AnyOctant, c::DGG.LevelIndex, l::Integer)
    l0 = DGG.level(c)
    l < l0 && throw(ArgumentError("descendant level $l is above level($c) = $l0"))
    w = 4^(Int(l) - l0)
    lo = Int(DGG.rawid(c)) * w
    (lo + 1):(lo + w)
end
