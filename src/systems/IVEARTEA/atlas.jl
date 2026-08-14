# Ten unit rhombi in the DGGAL 5x6 staircase.  A local `(x,y)` is added to
# ROOT_ORIGINS[root+1].  Each square is split along x == y.

const ROOT_ORIGINS = ((0,0),(0,1),(1,1),(1,2),(2,2),(2,3),(3,3),(3,4),(4,4),(4,5))
const ROOT_LOWER_FACES = (0,10,1,11,2,12,3,13,4,14)
const ROOT_UPPER_FACES = (5,15,6,16,7,17,8,18,9,19)

# `(v00,v10,v11,v01)` at local `(0,0),(1,0),(1,1),(0,1)`.
const ROOT_VERTICES = (
    (1,0,2,6),(6,2,7,11),(2,0,3,7),(7,3,8,11),(3,0,4,8),
    (8,4,9,11),(4,0,5,9),(9,5,10,11),(5,0,1,10),(10,1,6,11),
)

@inline _profile(::Type{<:Union{IVEA4RSystem,IVEA9RSystem}}) = IVEAProfile()
@inline _profile(::Type{<:Union{RTEA4RSystem,RTEA9RSystem}}) = RTEAProfile()
@inline _profile(sys) = _profile(typeof(sys))

function chart_point(profile::ProjectionProfile, x::Real,y::Real,root::Integer)
    ox,oy=ROOT_ORIGINS[Int(root)+1]
    v=inverse_5x6(profile,(ox+Float64(x),oy+Float64(y)))
    return GO.UnitSphericalPoint(v[1],v[2],v[3])
end

chart_point(sys,x,y,root)=chart_point(_profile(sys),x,y,root)

function point_to_chart(profile::ProjectionProfile,p)
    x,y=forward_5x6(profile,p)
    # The ten roots are the staircase squares.  Select by containment, then a
    # deterministic nearest-square fallback for seam roundoff.
    best=Inf; br=0; bx=0.0; by=0.0
    for r in 0:9
        ox,oy=ROOT_ORIGINS[r+1]
        lx=x-ox; ly=y-oy
        d=max(0.0,-lx,lx-1,-ly,ly-1)
        if d<best; best=d; br=r; bx=lx; by=ly end
        d<=2e-13 && break
    end
    return clamp(bx,0.0,1.0),clamp(by,0.0,1.0),br
end

point_to_chart(sys,p)=point_to_chart(_profile(sys),p)

@inline function rowmajor(ix::Integer,iy::Integer,root::Integer,nside::Integer)
    return Int64(root)*Int64(nside)^2 + Int64(iy)*Int64(nside)+Int64(ix)
end

function from_rowmajor(id::Integer,nside::Integer)
    n=Int64(nside); q=Int64(id)
    0<=q<10n*n || throw(ArgumentError("cell id $q is out of range for nside=$n"))
    root,r=divrem(q,n*n); iy,ix=divrem(r,n)
    return Int(ix),Int(iy),Int(root)
end

function cell_corners(sys,ix,iy,root,nside)
    n=Float64(nside); x0=ix/n; x1=(ix+1)/n; y0=iy/n; y1=(iy+1)/n
    # The 5x6 chart has downward-positive orientation relative to conventional
    # lon/lat; this order is CCW on the outward sphere and matches DGGAL roots.
    return (chart_point(sys,x0,y0,root),chart_point(sys,x0,y1,root),
        chart_point(sys,x1,y1,root),chart_point(sys,x1,y0,root))
end

function cell_center(sys,ix,iy,root,nside)
    x=(ix+0.5)/nside; y=(iy+0.5)/nside
    # Deep aperture-9 cells adjacent to an internal fundamental-triangle apex
    # are smaller than the directional condition number of the RTEA forward
    # construction. Return a point one eighth of a cell farther from the apex;
    # `cell_centroid` promises an interior representative, not necessarily the
    # area centroid.
    if sys isa RTEA9RSystem && nside >= 3^15
        d=0.125/Float64(nside)
        if abs(x-0.5) <= 2.5/Float64(nside) && abs(y-0.5) <= 2.5/Float64(nside)
            if x==0.5 && y==0.5
                x+=d; y-=d
            else
                x+=copysign(d,x-0.5); y+=copysign(d,y-0.5)
            end
        end
    end
    return chart_point(sys,x,y,root)
end

function perimeter(sys,ix,iy,root,nside,nseg)
    n=Float64(nside); x0=ix/n; x1=(ix+1)/n; y0=iy/n; y1=(iy+1)/n
    pts=Vector{GO.UnitSphericalPoint{Float64}}(undef,4nseg); k=0
    for i in 0:nseg-1; t=i/nseg; pts[k+=1]=chart_point(sys,x0,y0+t*(y1-y0),root) end
    for i in 0:nseg-1; t=i/nseg; pts[k+=1]=chart_point(sys,x0+t*(x1-x0),y1,root) end
    for i in 0:nseg-1; t=i/nseg; pts[k+=1]=chart_point(sys,x1,y1-t*(y1-y0),root) end
    for i in 0:nseg-1; t=i/nseg; pts[k+=1]=chart_point(sys,x1-t*(x1-x0),y0,root) end
    return pts
end
