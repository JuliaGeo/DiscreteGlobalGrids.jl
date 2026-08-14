# Integer atlas topology.  Edge partners and corner fans are derived from the
# root vertex table, rather than copied from oracle neighbor rows.

const EDGE_S,EDGE_E,EDGE_N,EDGE_W=1,2,3,4
const CORNER_00,CORNER_10,CORNER_11,CORNER_01=1,2,3,4

@inline function edge_pair(root,e)
    a,b,c,d=ROOT_VERTICES[root+1]
    e==EDGE_S && return (a,b); e==EDGE_E && return (b,c)
    e==EDGE_N && return (c,d); return (d,a)
end

const EDGE_PARTNER = ntuple(10) do ri
    r=ri-1
    ntuple(4) do e
        a,b=edge_pair(r,e)
        hits=[(s,f) for s in 0:9 for f in 1:4 if edge_pair(s,f)==(b,a)]
        @assert length(hits)==1
        only(hits)
    end
end

@inline function edge_cell(e,j,n)
    e==EDGE_S && return (j,0); e==EDGE_E && return (n-1,j)
    e==EDGE_N && return (n-1-j,n-1); return (0,n-1-j)
end
@inline function corner_cell(c,n)
    c==CORNER_00 && return (0,0); c==CORNER_10 && return (n-1,0)
    c==CORNER_11 && return (n-1,n-1); return (0,n-1)
end
@inline corner_slot(x,y) = x<0 ? (y<0 ? CORNER_00 : CORNER_01) : (y<0 ? CORNER_10 : CORNER_11)

const OFFSETS=((1,0),(1,1),(0,1),(-1,1),(-1,0),(-1,-1),(0,-1),(1,-1))

function lattice_neighbors(ix,iy,root,n,connectivity)
    slots=connectivity isa DGG.Edge ? (1,3,5,7) : (1,2,3,4,5,6,7,8)
    out=NTuple{3,Int}[]
    add(v)=v in out ? nothing : push!(out,v)
    for k in slots
        dx,dy=OFFSETS[k]; xx=ix+dx; yy=iy+dy
        inx=0<=xx<n; iny=0<=yy<n
        if inx&&iny
            add((xx,yy,root))
        elseif inx||iny
            e,j = xx>=n ? (EDGE_E,yy) : xx<0 ? (EDGE_W,n-1-yy) : yy>=n ? (EDGE_N,n-1-xx) : (EDGE_S,xx)
            r2,e2=EDGE_PARTNER[root+1][e]
            x2,y2=edge_cell(e2,n-1-j,n); add((x2,y2,r2))
        else
            c=corner_slot(xx,yy); v=ROOT_VERTICES[root+1][c]
            for r2 in 0:9, c2 in 1:4
                ROOT_VERTICES[r2+1][c2]==v || continue
                x2,y2=corner_cell(c2,n); (x2,y2,r2)==(ix,iy,root) || add((x2,y2,r2))
            end
        end
    end
    return out
end

@inline _xyz(p)=(Float64(p[1]),Float64(p[2]),Float64(p[3]))

function sort_ccw!(cells,sys,subject,nside)
    ix,iy,r=subject; c=_xyz(cell_center(sys,ix,iy,r,nside))
    # Start ray is the local +x neighbor direction (or its seam continuation).
    refcell=first(lattice_neighbors(ix,iy,r,nside,DGG.Edge()))
    rx,ry,rr=refcell; ref=_xyz(cell_center(sys,rx,ry,rr,nside))
    u=vnormalize(vadd(ref,vscale(c,-vdot(ref,c)))); w=vcross(c,u)
    az(cell)=begin
        x,y,q=cell; p=_xyz(cell_center(sys,x,y,q,nside)); t=vadd(p,vscale(c,-vdot(p,c)))
        mod(atan(vdot(t,w),vdot(t,u)),2pi)
    end
    sort!(cells,by=x->(az(x),rowmajor(x...,nside)))
end
