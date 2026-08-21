module IVEARTEASystemTests

using Test, Random
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
import GeometryOps as GO
using DiscreteGlobalGridsConformanceTesting
const US=GO.UnitSpherical

# Allows this family-owned test to run before the parent branch adds shared
# package wiring.  Once wired, the guard is a no-op.
isdefined(DGG,:IVEARTEA) || Base.include(DGG,joinpath(@__DIR__,"../../../src/systems/IVEARTEA/IVEARTEA.jl"))
const IR=DGG.IVEARTEA

xyz(lon,lat)=(cosd(lat)*cosd(lon),cosd(lat)*sind(lon),sind(lat))
point(v)=GO.UnitSphericalPoint(v...)
sd(a,b)=US.spherical_distance(a,b)

const SYSTEMS=(IR.IVEA4RSystem(),IR.IVEA9RSystem(),IR.RTEA4RSystem(),IR.RTEA9RSystem())

@testset "fundamental slice-and-dice kernel" begin
    rng=MersenneTwister(0x1fea7ea)
    for profile in (IR.IVEAProfile(),IR.RTEAProfile())
        worst=0.0
        for _ in 1:500
            z=2rand(rng)-1; lon=360rand(rng)-180
            p=xyz(lon,asind(z)); q=IR.inverse_5x6(profile,IR.forward_5x6(profile,p))
            worst=max(worst,acos(clamp(sum(p[i]*q[i] for i in 1:3),-1,1)))
        end
        @test worst < 4e-8
        # All special vertices and all 20 face centres are exact/stable cases.
        for p in (IR.ICO_VERTICES...,IR.FACE_CENTERS...)
            q=IR.inverse_5x6(profile,IR.forward_5x6(profile,p))
            @test acos(clamp(IR.vdot(p,q),-1,1)) < 3e-8
        end
    end
end

@testset "DGGAL 0.0.6 level-zero reconnaissance" begin
    # Geodetic oracle latitudes are converted explicitly; the kernel is authalic.
    a=DGG.Helpers.WGS84_AUTHALIC
    expected_centres=[
        (-110.51747441146,54.1219206418414),(-110.51747441146,-18.0755686040723),
        (-47.082525588539,18.0755686040722),(-47.082525588539,-54.1219206418414),
        (11.2,0.0),(69.482525588539,-54.1219206418414),(69.482525588539,18.0755686040723),
        (132.91747441146,-18.0755686040723),(132.91747441146,54.1219206418414),(-168.8,0.0)]
    for sys in (IR.IVEA4RSystem(),IR.IVEA9RSystem())
        g=levelgrid(sys,0)
        for (i,(lon,latg)) in enumerate(expected_centres)
            lat=DGG.Helpers.geodetic_to_authalicd(a,latg)
            @test sd(cell_centroid(g,cellindex(g,i)),point(xyz(lon,lat))) < 3e-13
        end
        @test map(rawid,neighbors(g,cellindex(g,1);connectivity=Edge())) |> Set == Set((1,2,8,9))
    end
    # Root zero corners from cells.jsonl, converted to the spherical kernel.
    expected=[(-168.8,58.3971459074313),(-137.08252558854,0.0),
        (-78.8,31.832359041336),(11.2,58.397145907431)]
    got=IR.cell_corners(IR.IVEA4RSystem(),0,0,0,1)
    for (p,(lon,latg)) in zip(got,expected)
        @test sd(p,point(xyz(lon,DGG.Helpers.geodetic_to_authalicd(a,latg)))) < 5e-13
    end
end

@testset "rhombic system laws" begin
    for sys in SYSTEMS
        s=IR._scale(sys)
        @test ncells(levelgrid(sys,0))==10
        @test ncells(levelgrid(sys,3))==10*(s^2)^3
        @test length(rootcells(sys))==10
        @test !has_sorted_subtrees(sys)
        for l in 0:3
            g=levelgrid(sys,l)
            @test sum(cell_area(g,cellindex(g,i)) for i in 1:ncells(g)) ≈ 4pi atol=2e-12
            for i in unique(round.(Int,range(1,ncells(g),length=min(35,ncells(g)))))
                c=cellindex(g,i)
                @test cellposition(g,c)==i
                @test cellat(g,cell_centroid(g,c))==c
                @test all(sd(p,point(Tuple(p)))==0 for p in cell_boundary(g,c))
                ns=neighbors(g,c;connectivity=Edge())
                @test length(ns)==4
                @test all(c in neighbors(g,n;connectivity=Edge()) for n in ns)
                if l>0
                    p=parent(sys,c); @test c in children(sys,p)
                    @test ancestor(sys,c,0) in rootcells(sys)
                end
            end
        end
    end
end

@testset "profile distinction below root level" begin
    for (a,b) in ((IR.IVEA4RSystem(),IR.RTEA4RSystem()),(IR.IVEA9RSystem(),IR.RTEA9RSystem()))
        ga,gb=levelgrid(a,2),levelgrid(b,2)
        distances=[sd(cell_centroid(ga,cellindex(ga,i)),cell_centroid(gb,cellindex(gb,i))) for i in 1:ncells(ga)]
        @test maximum(distances)>1e-3
        # Root chart centres lie on fundamental-triangle boundaries and remain
        # shared, while refined cell centres expose the radial-profile change.
        @test sd(cell_centroid(levelgrid(a,0),cellindex(levelgrid(a,0),1)),
            cell_centroid(levelgrid(b,0),cellindex(levelgrid(b,0),1)))<1e-12
    end
end

@testset "deep advertised-level centroid lookup" begin
    for sys in SYSTEMS
        l=last(levels(sys)); g=levelgrid(sys,l); n=IR._nside(sys,l)
        probes=unique((0,1,max(0,n÷2-1),n÷2,min(n-1,n÷2+1),n-2,n-1))
        for r in 0:9, ix in probes, iy in probes
            c=DGG.LevelIndex(l,IR.rowmajor(ix,iy,r,n))
            @test cellat(g,cell_centroid(g,c))==c
        end
    end
end

@testset "package interface conformance" begin
    for sys in SYSTEMS
        test_hierarchical_system(sys; levels=0:3, n_levels=4, n_samples=16,
            label=string(nameof(typeof(sys))))
    end
end

end # module
