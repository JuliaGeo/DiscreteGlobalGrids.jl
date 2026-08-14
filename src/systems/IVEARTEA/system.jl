# Complete package interface for the four rhombic slice-and-dice systems.

const RhombicSystem=Union{IVEA4RSystem,IVEA9RSystem,RTEA4RSystem,RTEA9RSystem}
const RhombicLevelGrid=DGG.HierarchicalLevelGrid{S} where {S<:AbstractSliceDiceRhombicSystem}

@inline _scale(::Union{IVEA4RSystem,RTEA4RSystem})=2
@inline _scale(::Union{IVEA9RSystem,RTEA9RSystem})=3
@inline _aperture(sys)=_scale(sys)^2
@inline _nside(sys,l::Integer)=Int64(_scale(sys))^Int(l)
@inline _ncells(sys,l::Integer)=10*_nside(sys,l)^2
@inline _maxlevel(::Union{IVEA4RSystem,RTEA4RSystem})=25
@inline _maxlevel(::Union{IVEA9RSystem,RTEA9RSystem})=16
@inline _name(::IVEA4RSystem)="IVEA4R"; @inline _name(::IVEA9RSystem)="IVEA9R"
@inline _name(::RTEA4RSystem)="RTEA4R"; @inline _name(::RTEA9RSystem)="RTEA9R"

DGG.cellindextype(::RhombicSystem)=DGG.LevelIndex
DGG.cellindextypes(::RhombicSystem)=(DGG.LevelIndex,)
DGG.levels(sys::RhombicSystem)=0:_maxlevel(sys)
DGG.has_sorted_subtrees(::RhombicSystem)=false
DGG.max_neighbors(::RhombicSystem,::DGG.Vertex)=9
DGG.max_neighbors(::RhombicSystem,::DGG.Edge)=4
DGG.rootcells(::RhombicSystem)=[DGG.LevelIndex(0,r) for r in 0:9]

DGG.ncells(sys::RhombicSystem,l::Integer)=Int(_ncells(sys,l))
DGG.cellindex(::RhombicSystem,l::Integer,i::Int)=DGG.LevelIndex(l,i-1)
function DGG.cellposition(sys::RhombicSystem,c::DGG.LevelIndex)
    l=DGG.level(c); l in DGG.levels(sys) || return nothing
    0<=c.index<_ncells(sys,l) || return nothing
    return Int(c.index+1)
end

function checked(sys::RhombicSystem,c::DGG.LevelIndex)
    l=DGG.level(c); l in DGG.levels(sys) || throw(ArgumentError("$(_name(sys)) level $l is invalid"))
    0<=c.index<_ncells(sys,l) || throw(ArgumentError("$(_name(sys)) id $(c.index) is invalid at level $l"))
    return c.index
end

function Base.parent(sys::RhombicSystem,c::DGG.LevelIndex)
    l=DGG.level(c); l>0 || throw(ArgumentError("$(_name(sys)) root $c has no parent")); checked(sys,c)
    n=_nside(sys,l); ix,iy,r=from_rowmajor(c.index,n); s=_scale(sys); pn=n÷s
    return DGG.LevelIndex(l-1,rowmajor(ix÷s,iy÷s,r,pn))
end

function DGG.children(sys::RhombicSystem,c::DGG.LevelIndex)
    l=DGG.level(c); l<DGG.max_level(sys) || throw(ArgumentError("$c is at max_level")); checked(sys,c)
    n=_nside(sys,l); ix,iy,r=from_rowmajor(c.index,n); s=_scale(sys); cn=n*s
    return [DGG.LevelIndex(l+1,rowmajor(s*ix+dx,s*iy+dy,r,cn)) for dy in 0:s-1 for dx in 0:s-1]
end

function DGG.ancestor(sys::RhombicSystem,c::DGG.LevelIndex,target::Integer)
    t=Int(target); l=DGG.level(c); 0<=t<=l || throw(ArgumentError("ancestor level $t is invalid for level $l")); checked(sys,c)
    t==l && return c
    n=_nside(sys,l); ix,iy,r=from_rowmajor(c.index,n); q=_scale(sys)^(l-t); tn=_nside(sys,t)
    return DGG.LevelIndex(t,rowmajor(ix÷q,iy÷q,r,tn))
end

# Curved inverse-projected chart edges are represented by 32 great-circle
# segments per edge.  This is a polygonal geometry contract (not the analytic
# area calculation); at level zero it reduces midpoint-to-chord deviation from
# about 3.1e-4 rad at 8 segments to below 2e-5 rad.
const BOUNDARY_SEGMENTS=32
function DGG.cell_boundary(sys::RhombicSystem,c::DGG.LevelIndex)
    n=_nside(sys,DGG.level(c)); ix,iy,r=from_rowmajor(checked(sys,c),n)
    return perimeter(sys,ix,iy,r,n,BOUNDARY_SEGMENTS)
end
function DGG.cell_centroid(sys::RhombicSystem,c::DGG.LevelIndex)
    n=_nside(sys,DGG.level(c)); ix,iy,r=from_rowmajor(checked(sys,c),n)
    return cell_center(sys,ix,iy,r,n)
end
function DGG.cell_area(g::RhombicLevelGrid,c::DGG.LevelIndex)
    DGG.level(c)==g.level || throw(ArgumentError("cell level does not match grid")); checked(g.system,c)
    return 4Float64(pi)/_ncells(g.system,g.level)
end

function DGG.node_extent(sys::RhombicSystem,c::DGG.LevelIndex)
    n=_nside(sys,DGG.level(c)); ix,iy,r=from_rowmajor(checked(sys,c),n)
    ctr=cell_center(sys,ix,iy,r,n); pts=perimeter(sys,ix,iy,r,n,8)
    rad=maximum(US.spherical_distance(ctr,p) for p in pts)
    gap=maximum(US.spherical_distance(pts[i],pts[mod1(i+1,length(pts))]) for i in eachindex(pts))
    return SphericalCap(ctr,nextfloat(min(Float64(pi),rad+gap/2)))
end

function DGG.cellat(g::RhombicLevelGrid,p::GO.UnitSphericalPoint)
    x,y,r=point_to_chart(g.system,p); n=_nside(g.system,g.level)
    ix=min(floor(Int,x*n),n-1); iy=min(floor(Int,y*n),n-1)
    return DGG.LevelIndex(g.level,rowmajor(ix,iy,r,n))
end

function one_ring(g::RhombicLevelGrid,c::DGG.LevelIndex,conn::DGG.Connectivity)
    DGG.level(c)==g.level || throw(ArgumentError("cell level does not match grid"))
    n=Int(_nside(g.system,g.level)); ix,iy,r=from_rowmajor(checked(g.system,c),n)
    cells=lattice_neighbors(ix,iy,r,n,conn); sort_ccw!(cells,g.system,(ix,iy,r),n)
    return [DGG.LevelIndex(g.level,rowmajor(x,y,q,n)) for (x,y,q) in cells]
end

function sort_ring!(ids,g::RhombicLevelGrid,subject::DGG.LevelIndex,ref::DGG.LevelIndex)
    isempty(ids) && return ids
    c=_xyz(DGG.cell_centroid(g,subject)); rp=_xyz(DGG.cell_centroid(g,ref))
    rt=vadd(rp,vscale(c,-vdot(rp,c))); vnorm(rt)<1e-12 && return sort!(ids)
    u=vnormalize(rt); w=vcross(c,u)
    az(id)=begin
        p=_xyz(DGG.cell_centroid(g,id)); t=vadd(p,vscale(c,-vdot(p,c)))
        n=vnorm(t); n<1e-12 ? (2pi,DGG.rawid(id)) : (mod(atan(vdot(t,w),vdot(t,u)),2pi),DGG.rawid(id))
    end
    return sort!(ids,by=az)
end

function DGG.neighbors(g::RhombicLevelGrid,c::DGG.LevelIndex,k::Integer=1;connectivity::DGG.Connectivity=DGG.Vertex())
    steps=Int(k); steps>=0 || throw(ArgumentError("k must be non-negative")); steps==0 && return DGG.LevelIndex[]
    seen=Set{DGG.LevelIndex}((c,)); frontier=DGG.LevelIndex[c]; result=DGG.LevelIndex[]
    ref=first(one_ring(g,c,connectivity))
    for j in 1:steps
        nxt=DGG.LevelIndex[]
        for x in frontier, y in one_ring(g,x,connectivity)
            y in seen && continue; push!(seen,y); push!(nxt,y)
        end
        j>1 && sort_ring!(nxt,g,c,ref)
        append!(result,nxt); frontier=nxt
    end
    return result
end

function DGG.ring(g::RhombicLevelGrid,c::DGG.LevelIndex,k::Integer;connectivity::DGG.Connectivity=DGG.Vertex())
    steps=Int(k); steps>=0 || throw(ArgumentError("k must be non-negative")); steps==0 && return [c]
    seen=Set{DGG.LevelIndex}((c,)); frontier=DGG.LevelIndex[c]
    ref=first(one_ring(g,c,connectivity))
    for j in 1:steps
        nxt=DGG.LevelIndex[]
        for x in frontier,y in one_ring(g,x,connectivity)
            y in seen && continue; push!(seen,y); push!(nxt,y)
        end
        j>1 && sort_ring!(nxt,g,c,ref)
        frontier=nxt
    end
    return frontier
end
