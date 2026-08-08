module ManifoldAdapterTestSuite

using Test

using DiscreteGlobalGrids
using DiscreteGlobalGrids.Helpers: AuthalicTransform, WGS84_AUTHALIC,
    authalic_radius, geodetic_to_authalicd, authalic_to_geodeticd
import GeometryOps as GO
import GeometryOpsCore as GOCore

@testset "AuthalicTransform from GO manifolds" begin
    @testset "Geodesic carries the ellipsoid through unchanged" begin
        # `Geodesic`'s defaults are WGS84, so this must reproduce the constant
        # the helper ships, field for field.
        t = AuthalicTransform(GO.Geodesic())
        @test t.semimajor_axis === WGS84_AUTHALIC.semimajor_axis
        @test t.eccentricity_squared === WGS84_AUTHALIC.eccentricity_squared
        @test t.authalic_radius === WGS84_AUTHALIC.authalic_radius
        @test t.fwd === WGS84_AUTHALIC.fwd
        @test t.inv === WGS84_AUTHALIC.inv
        @test t === WGS84_AUTHALIC

        # A non-WGS84 ellipsoid must actually be read from the manifold rather
        # than defaulted: GRS80's inverse flattening differs from WGS84's in
        # the 9th significant digit, which is far above the 0.5 ulp the
        # transform itself is accurate to.
        grs80 = GO.Geodesic(; semimajor_axis=6378137.0, inv_flattening=298.257222101)
        g = AuthalicTransform(grs80)
        @test g.semimajor_axis == 6378137.0
        @test g !== WGS84_AUTHALIC
        @test g.eccentricity_squared != WGS84_AUTHALIC.eccentricity_squared
        # ...but the two ellipsoids are close, so the radii must agree to well
        # under a millimetre. This pins that we read `inv_flattening` as an
        # inverse and did not, say, feed it in as `flattening`.
        @test isapprox(g.authalic_radius, WGS84_AUTHALIC.authalic_radius; atol=1e-3)
    end

    @testset "inv_flattening is not confused with flattening" begin
        # The single most damaging way to misread `Geodesic` is to pass
        # `inv_flattening` where `flattening` belongs. That would be rejected
        # outright (298 ∉ [0, 1)), so assert the error rather than a value.
        @test_throws DomainError AuthalicTransform{Float64}(;
            semimajor_axis=6378137.0, flattening=298.257223563)
    end

    @testset "Spherical gives the identity transform" begin
        R = 6371007.180918474
        t = AuthalicTransform(GO.Spherical(; radius=R))
        @test t.eccentricity_squared == 0
        @test t.authalic_radius == R
        @test authalic_radius(t) == R
        # e² = 0 ⇒ every series coefficient vanishes ⇒ both directions are
        # exactly the identity, so a spherical grid needs no branch at the
        # call site.
        @test all(iszero, t.fwd)
        @test all(iszero, t.inv)
        for lat in (-90.0, -37.5, 0.0, 17.25, 45.0, 90.0)
            @test geodetic_to_authalicd(t, lat) === lat
            @test authalic_to_geodeticd(t, lat) === lat
        end
    end

    @testset "element type follows the manifold, and can be overridden" begin
        @test AuthalicTransform(GO.Geodesic()) isa AuthalicTransform{Float64}
        @test AuthalicTransform(GO.Spherical()) isa AuthalicTransform{Float64}
        @test AuthalicTransform(GO.Geodesic(; semimajor_axis=6378137.0f0,
            inv_flattening=298.2572f0)) isa AuthalicTransform{Float32}
        @test AuthalicTransform(GO.Spherical(; radius=6.371f6)) isa AuthalicTransform{Float32}
        @test AuthalicTransform{Float32}(GO.Geodesic()) isa AuthalicTransform{Float32}
        @test AuthalicTransform{Float64}(GO.Spherical(; radius=6.371f6)) isa
              AuthalicTransform{Float64}
    end

    @testset "manifolds without an ellipsoid are rejected with a reason" begin
        # `Planar` has no ellipsoid to read; `AutoManifold` deliberately
        # carries no parameters at all and would have to be guessed at. For a
        # DGGS the manifold is part of the reference-system definition, so
        # guessing is exactly the silent misregistration this adapter exists
        # to prevent.
        @test_throws ArgumentError AuthalicTransform(GO.Planar())
        @test_throws ArgumentError AuthalicTransform(GO.AutoManifold())
        @test_throws ArgumentError AuthalicTransform{Float32}(GO.Planar())
    end
end

@testset "authalic_sphere" begin
    @testset "Geodesic resolves to the equal-area sphere" begin
        s = authalic_sphere(GO.Geodesic())
        @test s isa GOCore.Spherical
        @test s.radius === WGS84_AUTHALIC.authalic_radius

        # This is the whole point of the bridge: the compute manifold must NOT
        # be GO's default `Spherical()`, whose radius is the WGS84 *mean*
        # radius (6371008.8 m). That sphere is not area-preserving; using it
        # would put a 1.6 m scale error into every cell area.
        @test s.radius != GO.Spherical().radius
        @test abs(s.radius - GO.Spherical().radius) > 1.0
    end

    @testset "an AuthalicTransform resolves to its own radius" begin
        @test authalic_sphere(WGS84_AUTHALIC) === authalic_sphere(GO.Geodesic())
        t = AuthalicTransform{Float64}(; semimajor_axis=1.0, eccentricity_squared=0.0)
        @test authalic_sphere(t).radius == 1.0
    end

    @testset "Spherical is returned unchanged" begin
        # Already the compute manifold — resolving must be idempotent, so that
        # a spherical grid and an ellipsoidal one can share one code path.
        s = GO.Spherical(; radius=1.0)
        @test authalic_sphere(s) === s
        @test authalic_sphere(authalic_sphere(GO.Geodesic())) === authalic_sphere(GO.Geodesic())
    end

    @testset "rejects manifolds with no ellipsoid" begin
        @test_throws ArgumentError authalic_sphere(GO.Planar())
        @test_throws ArgumentError authalic_sphere(GO.AutoManifold())
    end
end

@testset "adapter is free at run time" begin
    # The adapter must be constant-folded away: a `const` manifold has to
    # produce a `const` transform with no allocation, else putting it on a
    # grid type would cost something per boundary.
    @test (@allocated AuthalicTransform(GO.Geodesic())) == 0
    @test (@allocated authalic_sphere(GO.Geodesic())) == 0
    @test @inferred(AuthalicTransform(GO.Geodesic())) isa AuthalicTransform{Float64}
    @test @inferred(authalic_sphere(GO.Geodesic())) isa GOCore.Spherical{Float64}
end

@testset "grids carry a manifold" begin
    WGS84 = GO.Geodesic()
    ids = collect(Int64, 0:47)

    @testset "default is the plain sphere, and nothing moves" begin
        # The default must reproduce today's behaviour exactly: every
        # `best_manifold` in this package used to be a bare `GO.Spherical()`.
        g = DGGSGrid(HEALPixDGGS(), 1)
        @test g.manifold === GO.Spherical()
        @test GOCore.best_manifold(g) === GO.Spherical()

        p = DGGSPartialGrid(HEALPixDGGS(), 1, ids)
        @test p.manifold === GO.Spherical()
        @test GOCore.best_manifold(p) === GO.Spherical()
    end

    @testset "a declared ellipsoid is stored, and resolves for compute" begin
        g = DGGSGrid(HEALPixDGGS(), 1; manifold=WGS84)
        # Stored as declared — the ellipsoid is part of the grid's identity...
        @test g.manifold === WGS84
        @test g.manifold isa GOCore.Geodesic
        # ...but `best_manifold` hands the tree/CR layer the authalic sphere,
        # because `Geodesic` is a `MethodError` in `ConservativeRegridding`.
        m = GOCore.best_manifold(g)
        @test m isa GOCore.Spherical
        @test m === authalic_sphere(WGS84)
        @test m.radius === WGS84_AUTHALIC.authalic_radius

        p = DGGSPartialGrid(HEALPixDGGS(), 1, ids; manifold=WGS84)
        @test p.manifold === WGS84
        @test GOCore.best_manifold(p) === authalic_sphere(WGS84)
    end

    @testset "the manifold reaches the type, so grids stay concrete" begin
        g = DGGSGrid(HEALPixDGGS(), 1; manifold=WGS84)
        @test isbits(g)
        @test isconcretetype(typeof(g))
        # Two ellipsoids give two distinct grid types, so a mismatch is at
        # least visible to dispatch rather than hidden in a field.
        @test typeof(g) !== typeof(DGGSGrid(HEALPixDGGS(), 1))
        @test @inferred(GOCore.best_manifold(g)) isa GOCore.Spherical
    end

    @testset "treeify and the cursor inherit it" begin
        for man in (GO.Spherical(), WGS84)
            g = DGGSGrid(HEALPixDGGS(), 1; manifold=man)
            tree = treeify(g)
            @test GOCore.best_manifold(tree) === GOCore.best_manifold(g)
            @test ncells(tree) == 48

            p = DGGSPartialGrid(HEALPixDGGS(), 1, ids; manifold=man)
            ptree = treeify(p)
            @test GOCore.best_manifold(ptree) === GOCore.best_manifold(p)
            @test ncells(ptree) == length(ids)
        end
    end

    @testset "subtree_grid forwards it" begin
        s = subtree_grid(HEALPixDGGS(), Int64(0); root_level=0, leaf_level=2,
            manifold=WGS84)
        @test s.manifold === WGS84
        @test GOCore.best_manifold(s) === authalic_sphere(WGS84)
        @test subtree_grid(HEALPixDGGS(), Int64(0); root_level=0,
            leaf_level=2).manifold === GO.Spherical()
    end

    @testset "the lookup path forwards it" begin
        l = DiscreteGlobalGrids.HEALPix.HealpixLookups.HealpixLookup(ids; level=1)
        @test DGGSPartialGrid(l).manifold === GO.Spherical()
        @test DGGSPartialGrid(l; manifold=WGS84).manifold === WGS84
        # The id vector must still be passed by reference, not copied — the
        # whole Regridder-lines-up-with-data guarantee rests on it.
        @test DGGSPartialGrid(l; manifold=WGS84).ids === l.data
    end

    @testset "a manifold with no ellipsoid is rejected at construction" begin
        # Catching this at construction is the point: a `Planar` grid that
        # built fine and failed later would be a grid whose frame nobody
        # declared.
        @test_throws ArgumentError DGGSGrid(HEALPixDGGS(), 1; manifold=GO.Planar())
        @test_throws ArgumentError DGGSGrid(HEALPixDGGS(), 1; manifold=GO.AutoManifold())
        @test_throws ArgumentError DGGSPartialGrid(HEALPixDGGS(), 1, ids;
            manifold=GO.Planar())
    end
end

end # module ManifoldAdapterTestSuite
