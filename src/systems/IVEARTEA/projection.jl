# The generic 120-fundamental-triangle slice-and-dice equal-area kernel.
# References: van Leeuwen & Strebe (2006), Hall et al. (2020), and the
# BSD-3-Clause DGGAL `icoVertexGreatCircle.ec` / `ri5x6.ec` at e16cea7.

abstract type ProjectionProfile end
struct IVEAProfile <: ProjectionProfile end
struct RTEAProfile <: ProjectionProfile end

const Vec3 = NTuple{3,Float64}
const Vec2 = NTuple{2,Float64}
const PHI = (1 + sqrt(5.0)) / 2
const LON0 = 11.2
const LAT_HI = atand(PHI)
const LAT_LO = 90 - LAT_HI
const FUNDAMENTAL_AREA = Float64(pi / 30)

@inline vdot(a::Vec3, b::Vec3) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
@inline vcross(a::Vec3, b::Vec3) = (a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
@inline vadd(a::Vec3, b::Vec3) = (a[1]+b[1], a[2]+b[2], a[3]+b[3])
@inline vscale(a::Vec3, s::Real) = (a[1]*s, a[2]*s, a[3]*s)
@inline vnorm(a::Vec3) = sqrt(vdot(a,a))
@inline vnormalize(a::Vec3) = vscale(a, inv(vnorm(a)))
@inline lonlat(lon, lat) = (cosd(lat)*cosd(lon), cosd(lat)*sind(lon), sind(lat))
@inline angle(a::Vec3, b::Vec3) = atan(vnorm(vcross(a,b)), clamp(vdot(a,b), -1.0, 1.0))

@inline function slerp(a::Vec3, b::Vec3, t::Real)
    d = angle(a,b)
    d < 2eps(Float64) && return a
    sd = sin(d)
    return vnormalize(vadd(vscale(a, sin((1-t)*d)/sd), vscale(b, sin(t*d)/sd)))
end

@inline function spherical_area(a::Vec3, b::Vec3, c::Vec3)
    num = abs(vdot(a, vcross(b,c)))
    den = 1 + vdot(a,b) + vdot(b,c) + vdot(c,a)
    return 2 * atan(num, den)
end

# DGGAL/OGC orientation, on the authalic unit sphere.  These closed-form
# longitudes are the same rotated regular icosahedron constructed by ri5x6.ec.
const ICO_VERTICES = (
    lonlat(LON0, LAT_HI),
    lonlat(LON0-180, LAT_HI),
    lonlat(LON0-90, LAT_LO),
    lonlat(LON0-LAT_LO, 0),
    lonlat(LON0+LAT_LO, 0),
    lonlat(LON0+90, LAT_LO),
    lonlat(LON0-180+LAT_LO, 0),
    lonlat(LON0-90, -LAT_LO),
    lonlat(LON0, -LAT_HI),
    lonlat(LON0+90, -LAT_LO),
    lonlat(LON0+180-LAT_LO, 0),
    lonlat(LON0-180, -LAT_HI),
)

const FACE_VERTICES = (
    (0,1,2),(0,2,3),(0,3,4),(0,4,5),(0,5,1),
    (6,2,1),(7,3,2),(8,4,3),(9,5,4),(10,1,5),
    (2,6,7),(3,7,8),(4,8,9),(5,9,10),(1,10,6),
    (11,7,6),(11,8,7),(11,9,8),(11,10,9),(11,6,10),
)

const FACE_PLANAR = (
    ((1.,0.),(0.,0.),(1.,1.)), ((2.,1.),(1.,1.),(2.,2.)),
    ((3.,2.),(2.,2.),(3.,3.)), ((4.,3.),(3.,3.),(4.,4.)),
    ((5.,4.),(4.,4.),(5.,5.)), ((0.,1.),(1.,1.),(0.,0.)),
    ((1.,2.),(2.,2.),(1.,1.)), ((2.,3.),(3.,3.),(2.,2.)),
    ((3.,4.),(4.,4.),(3.,3.)), ((4.,5.),(5.,5.),(4.,4.)),
    ((1.,1.),(0.,1.),(1.,2.)), ((2.,2.),(1.,2.),(2.,3.)),
    ((3.,3.),(2.,3.),(3.,4.)), ((4.,4.),(3.,4.),(4.,5.)),
    ((5.,5.),(4.,5.),(5.,6.)), ((0.,2.),(1.,2.),(0.,1.)),
    ((1.,3.),(2.,3.),(1.,2.)), ((2.,4.),(3.,4.),(2.,3.)),
    ((3.,5.),(4.,5.),(3.,4.)), ((4.,6.),(5.,6.),(4.,5.)),
)

const FACE_CENTERS = ntuple(20) do i
    a,b,c = FACE_VERTICES[i]
    vnormalize(vadd(vadd(ICO_VERTICES[a+1], ICO_VERTICES[b+1]), ICO_VERTICES[c+1]))
end

@inline midpoint(a::Vec3,b::Vec3) = vnormalize(vadd(a,b))
@inline midpoint(a::Vec2,b::Vec2) = ((a[1]+b[1])/2, (a[2]+b[2])/2)

"Build the six `(sphere vertices, planar vertices)` fundamental triangles of face `f`."
function build_fundamental_triangles(f::Int)
    ia,ib,ic = FACE_VERTICES[f+1]
    va,vb,vc = ICO_VERTICES[ia+1],ICO_VERTICES[ib+1],ICO_VERTICES[ic+1]
    pa,pb,pc = FACE_PLANAR[f+1]
    m = (midpoint(vb,vc), midpoint(vc,va), midpoint(va,vb))
    pm = (midpoint(pb,pc), midpoint(pc,pa), midpoint(pa,pb))
    ctr = FACE_CENTERS[f+1]
    pct = ((pa[1]+pb[1]+pc[1])/3, (pa[2]+pb[2]+pc[2])/3)
    return (
        ((m[1],vc,ctr),(pm[1],pc,pct)), ((m[1],vb,ctr),(pm[1],pb,pct)),
        ((m[2],vc,ctr),(pm[2],pc,pct)), ((m[2],va,ctr),(pm[2],pa,pct)),
        ((m[3],vb,ctr),(pm[3],pb,pct)), ((m[3],va,ctr),(pm[3],pa,pct)),
    )
end

# The twenty faces' triangles, built once. They are a function of the fixed
# icosahedron alone, and rebuilding them per projected point cost six midpoint
# normalisations on every call of the two hottest functions here.
const FUNDAMENTAL_TRIANGLES = ntuple(i -> build_fundamental_triangles(i - 1), 20)

"The six `(sphere vertices, planar vertices)` fundamental triangles of face `f`."
@inline fundamental_triangles(f::Int) = @inbounds FUNDAMENTAL_TRIANGLES[f+1]

@inline _radial_index(::IVEAProfile) = 2 # icosahedron vertex
@inline _radial_index(::RTEAProfile) = 1 # edge midpoint / RT face centre

function oriented_triangle(profile::ProjectionProfile, sv, pv)
    ai = _radial_index(profile)
    rest = ai == 1 ? (2,3) : (1,3)
    # Either order is valid if the corresponding planar vertices follow it.
    # Pick positive spherical winding to keep the inverse's signed triple stable.
    A = sv[ai]; B = sv[rest[1]]; C = sv[rest[2]]
    PA = pv[ai]; PB = pv[rest[1]]; PC = pv[rest[2]]
    if vdot(A,vcross(B,C)) < 0
        B,C = C,B; PB,PC = PC,PB
    end
    return A,B,C,PA,PB,PC
end

@inline function barycentric(p::Vec2, a::Vec2,b::Vec2,c::Vec2)
    det = (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1])
    wb = ((p[1]-a[1])*(c[2]-a[2])-(p[2]-a[2])*(c[1]-a[1]))/det
    wc = ((b[1]-a[1])*(p[2]-a[2])-(b[2]-a[2])*(p[1]-a[1]))/det
    return (1-wb-wc,wb,wc)
end

@inline function cartesian(w, a::Vec2,b::Vec2,c::Vec2)
    return (w[1]*a[1]+w[2]*b[1]+w[3]*c[1], w[1]*a[2]+w[2]*b[2]+w[3]*c[2])
end

function forward_triangle(p::Vec3, A::Vec3,B::Vec3,C::Vec3, PA::Vec2,PB::Vec2,PC::Vec2)
    angle(A,p) < 8eps(Float64) && return PA
    n1 = vcross(A,p); n2 = vcross(B,C)
    dline = vcross(n1,n2)
    vnorm(dline) < 8eps(Float64) && return PA
    D = vnormalize(dline)
    # Great-circle intersections are antipodal; retain the point on minor BC.
    if abs((angle(B,D)+angle(D,C))-angle(B,C)) > 1e-8
        D = vscale(D,-1)
    end
    den = sin(angle(A,D)/2)
    abs(den) < 8eps(Float64) && return PA
    h = clamp(sin(angle(A,p)/2)/den, 0.0, 1.0)
    q = clamp(spherical_area(A,B,D)/FUNDAMENTAL_AREA, 0.0, 1.0)
    return cartesian((1-h,h*(1-q),h*q), PA,PB,PC)
end

"""
    area_parameter(A, B, C, target) -> Float64

The `t` in `[0, 1]` at which `spherical_area(A, B, slerp(B, C, t))` reaches
`target`.

Monotone inversion of the published angular-area coordinate: a safeguarded
evaluation of the analytic construction, independent of any grid fixture, and
avoiding the ill-conditioned closed form at the vertices.

The map is smooth and strictly increasing from `0` to the triangle's own area,
so `[0, 1]` brackets the root from the start. Illinois — regula falsi with the
stale endpoint's value halved — keeps that bracket and converges superlinearly,
where plain bisection spends one area evaluation per bit of `t`. The arc `BC` is
fixed across the whole solve, so its length and the `sin` of it are hoisted out
of the loop; a `slerp` call would recompute both on every step.

This is the inner loop of the whole system: every `cell_centroid`,
`cell_boundary` vertex and `node_extent` is one call of it per point.
"""
function area_parameter(A::Vec3, B::Vec3, C::Vec3, target::Float64)
    d = angle(B,C)
    d < 2eps(Float64) && return 0.0
    sd = sin(d)
    flo = -target
    fhi = spherical_area(A,B,C) - target
    flo >= 0 && return 0.0
    fhi <= 0 && return 1.0
    lo, hi = 0.0, 1.0
    stale = 0
    for _ in 1:64
        hi - lo <= 4eps(hi) && break
        t = (lo*fhi - hi*flo) / (fhi - flo)
        (t > lo && t < hi) || (t = (lo + hi)/2)
        D = vnormalize(vadd(vscale(B, sin((1-t)*d)/sd), vscale(C, sin(t*d)/sd)))
        ft = spherical_area(A,B,D) - target
        ft == 0 && return t
        if ft < 0
            lo, flo = t, ft
            stale == 1 && (fhi /= 2)
            stale = 1
        else
            hi, fhi = t, ft
            stale == -1 && (flo /= 2)
            stale = -1
        end
    end
    return (lo + hi)/2
end

function inverse_triangle(p::Vec2, A::Vec3,B::Vec3,C::Vec3, PA::Vec2,PB::Vec2,PC::Vec2)
    w = barycentric(p,PA,PB,PC)
    w[1] > 1-2e-14 && return A
    w[2] > 1-2e-14 && return B
    w[3] > 1-2e-14 && return C
    h = clamp(w[2]+w[3],0.0,1.0)
    h < 2e-15 && return A
    target = clamp(w[3]/h,0.0,1.0)*FUNDAMENTAL_AREA

    D=slerp(B,C,area_parameter(A,B,C,target))
    ad=angle(A,D)
    ad < 2eps(Float64) && return A
    x=2asin(clamp(h*sin(ad/2),-1.0,1.0))
    return slerp(A,D,x/ad)
end

"Project an authalic unit vector to the DGGAL/OGC 5-by-6 plane."
function forward_5x6(profile::ProjectionProfile, p)
    q=vnormalize((Float64(p[1]),Float64(p[2]),Float64(p[3])))
    f=argmax(ntuple(i->vdot(q,FACE_CENTERS[i]),20))-1
    tris=fundamental_triangles(f)
    k=argmax(ntuple(i->begin sv,_=tris[i]; vdot(q,vnormalize(vadd(vadd(sv[1],sv[2]),sv[3]))) end,6))
    sv,pv=tris[k]
    return forward_triangle(q,oriented_triangle(profile,sv,pv)...)
end

function _face_for_planar(p::Vec2)
    best=-Inf; bf=-1
    for f in 0:19
        a,b,c=FACE_PLANAR[f+1]
        w=barycentric(p,a,b,c)
        m=min(w...)
        if m>best; best=m; bf=f end
    end
    return bf
end

"Inverse-project a point in the 5-by-6 atlas to an authalic unit vector."
function inverse_5x6(profile::ProjectionProfile, p)
    q=(Float64(p[1]),Float64(p[2]))
    f=_face_for_planar(q)
    tris=fundamental_triangles(f)
    best=-Inf; bk=1
    for k in 1:6
        _,pv=tris[k]; m=min(barycentric(q,pv...)...)
        if m>best; best=m; bk=k end
    end
    sv,pv=tris[bk]
    return inverse_triangle(q,oriented_triangle(profile,sv,pv)...)
end
