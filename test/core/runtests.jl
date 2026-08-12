module CoreTestSuite

using Test

using DiscreteGlobalGrids
# Only for the `cell_polygon` return-type check below; the rest of this file is
# deliberately dependency-free trait and interface arithmetic.
import GeoInterface as GI

@testset "DGGS registry" begin
    systems = collect(all_systems())
    @test length(systems) == 14 # report 1.8 has ISEA4R and ISEA9R split out
    names = system_name.(systems)
    @test :HEALPix in names
    @test :H3 in names
    @test :S2 in names
    @test root_count(HEALPixDGGS()) == 12
    # Twelve res-0 cells, one per icosahedron vertex, pinned by the sealed
    # oracle vectors (`test/IGeo7/vectors/res0_cells.csv`).
    @test root_count(IGEO7DGGS()) == 12
    @test radix(S2DGGS()) == 4
    @test radix(RHEALPixDGGS()) == 9
    @test root_count(ISEA4RDGGS()) == 10
    @test supports_prefix_ranges(ISEA4RDGGS())
    # ISEA9R carries the same ten roots as its aperture-4 twin, and for a
    # stronger reason than coincidence: OGC 21-038r1 Annex B.2 states them
    # normatively ("The ten root rhombuses are formed by combining two
    # icosahedron triangles at their base") and DGGAL's `countZones(0)` computes
    # `10 * 9^0`. Prefix ranges are true of the PACKAGE ordinal
    # (`diamond * 9^level + base-9 Morton`), not of DGGAL's zone id — see the
    # `ISEA9RDGGS` docstring and `docs/design/isea9r_layout.md`.
    @test root_count(ISEA9RDGGS()) == 10
    @test radix(ISEA9RDGGS()) == 9
    @test supports_prefix_ranges(ISEA9RDGGS())
end

@testset "trait interface" begin
    hp = HEALPixDGGS()
    @test system_name(hp) === :HEALPix
    @test grid_family(hp) === :healpix
    @test base_solid(hp) === :healpix_12_faces
    @test cell_shape(hp) === :curvilinear_quadrilateral
    @test is_equal_area(hp)
    @test aperture(hp) == 4
    @test canonical_index_name(hp) === :nested
    # Nested HEALPix is unbounded in principle; 29 is where `Int64` pixel
    # arithmetic (`12 * 4^level`) stops fitting.
    @test max_level(hp) == 29
    @test root_count(hp) == 12
    @test radix(hp) == 4
    @test supports_prefix_ranges(hp)

    @test max_level(H3DGGS()) == 15
    @test max_level(IGEO7DGGS()) == 19 # the Z7 encoding's deepest digit slot
    @test max_level(RTEADGGS()) === nothing # unbounded is a real answer, not a gap
    @test !is_equal_area(H3DGGS())
    @test canonical_index_name(H3DGGS()) === :h3_index

    # non-integer apertures survive the move off the spec struct
    @test aperture(A5DGGS()) === :implementation_defined
    @test aperture(IVEADGGS()) === :family
    @test aperture(IVEADGGS(:IVEA4R)) == 4

    # the DGGAL families carry a variant on the instance, which drives
    # name / cell shape / aperture / radix
    @test system_name(IVEADGGS()) === :IVEA_family
    @test system_name(IVEADGGS(:IVEA9R)) === :IVEA9R
    @test cell_shape(IVEADGGS(:IVEA9R)) === :rhomb
    @test cell_shape(IVEADGGS()) === :hexagon_or_rhomb
    @test radix(IVEADGGS(:IVEA9R)) == 9
    @test system_name(RTEADGGS()) === :RTEA_family
    @test base_solid(RTEADGGS()) === :rhombic_triacontahedron
    @test root_count(RTEADGGS()) == 30
    @test radix(RTEADGGS(:RTEA3H)) == 3

    # every registered system answers every universally-known fact
    for system in all_systems()
        @test system_name(system) isa Symbol
        @test grid_family(system) isa Symbol
        @test base_solid(system) isa Symbol
        @test cell_shape(system) isa Symbol
        @test is_equal_area(system) isa Bool
        @test aperture(system) isa Union{Int,Symbol,Tuple}
        @test canonical_index_name(system) isa Symbol
        @test max_level(system) isa Union{Int,Nothing}
        @test supports_prefix_ranges(system) isa Bool
    end
end

@testset "unverified facts fall back to NotPortedError" begin
    @test_throws NotPortedError radix(A5DGGS())          # no fixed verified radix
    @test_throws NotPortedError radix(IVEADGGS())        # family pins no aperture
    @test_throws NotPortedError root_count(ISEA3HDGGS())
    @test_throws NotPortedError root_count(ISEA4HDGGS())
    @test_throws NotPortedError root_count(IVEADGGS())

    # `ISEA9RDGGS` used to be the worked example here; its root count is
    # verified now (ten root rhombuses, OGC 21-038r1 Annex B.2 + DGGAL
    # `countZones`), so the example moved to a system whose root layout really
    # is unpinned.
    err = try
        root_count(ISEA3HDGGS())
    catch e
        e
    end
    @test err isa NotPortedError
    @test err.system === :ISEA3H
    @test err.operation === :root_count
    message = sprint(showerror, err)
    @test occursin("root_count", message)
    @test occursin("ISEA3H", message)
    @test occursin("not verified yet", message)
end

@testset "prefix leaf intervals" begin
    hp = HEALPixDGGS()
    @test leaf_count(hp, 0) == 12
    @test leaf_count(hp, 3) == 12 * 4^3
    @test leaf_interval(hp, 0, 2, 3) == 128:191
    @test leaf_interval(hp, 2, 4 * (4 * 2 + 1) + 3, 3) ==
          (4 * (4 * (4 * 2 + 1) + 3)):(4 * (4 * (4 * 2 + 1) + 3) + 3)

    rhp = RHEALPixDGGS()
    @test leaf_count(rhp, 2) == 6 * 9^2
    @test leaf_interval(rhp, 0, 5, 2) == (5 * 81):(6 * 81 - 1)

    isea4r = ISEA4RDGGS()
    @test leaf_count(isea4r, 2) == 10 * 4^2
    @test leaf_interval(isea4r, 0, 9, 2) == (9 * 16):(10 * 16 - 1)

    # Radix 9 over ten roots — the arithmetic
    # `supports_prefix_ranges(ISEA9RDGGS())` unlocked. The interval is the
    # base-9 Morton descendant block `[p * 9^Δ, (p + 1) * 9^Δ)`, exact because
    # the within-diamond lattice nesting is bit-exact
    # (`fl(ix/n) === fl(3ix/3n)`).
    isea9r = ISEA9RDGGS()
    @test leaf_count(isea9r, 0) == 10
    @test leaf_count(isea9r, 2) == 10 * 9^2
    @test leaf_count(isea9r, 5) == 10 * 9^5
    @test leaf_interval(isea9r, 0, 9, 2) == (9 * 81):(10 * 81 - 1)
    @test leaf_interval(isea9r, 1, 5, 2) == (5 * 9):(6 * 9 - 1)
    @test child_ids(isea9r, 0, 3) == collect(27:35)
end

# The same arithmetic, at the edge of `Int64`. `max_level` exists to name the
# wrap point (`12 * 4^30` overflows), so past it there is no cell to count and
# the answer is an error; where a system pins no maximum there is no line to
# draw, so the products are checked instead. Neither ever returns a wrapped
# count — `leaf_count(HEALPixDGGS(), 32)` used to answer a negative number.
@testset "trait arithmetic cannot wrap" begin
    hp = HEALPixDGGS()
    @test max_level(hp) == 29
    @test leaf_count(hp, 29) == 12 * 4^29 > 0
    @test_throws ArgumentError leaf_count(hp, 30)
    @test_throws ArgumentError leaf_count(hp, 32)
    @test_throws ArgumentError leaf_count(hp, -1)

    @test leaf_interval(hp, 0, 0, 29) == 0:(4^29 - 1)
    @test_throws ArgumentError leaf_interval(hp, 0, 0, 30)
    @test_throws ArgumentError leaf_interval(hp, 3, 0, 2)   # unchanged: level > leaf_level

    # Children live one level down, so that is the level that has to be
    # representable: a level-29 HEALPix cell has no level-30 children to name.
    @test child_ids(hp, 28, 3) == collect(12:15)
    @test_throws ArgumentError child_ids(hp, 29, 3)
    @test_throws ArgumentError child_ids(hp, 30, 3)

    # An unbounded system keeps constructing grids at absurd levels (the grid
    # types deliberately allow it) but cannot silently wrap when asked to count.
    rhp = RHEALPixDGGS()
    @test max_level(rhp) === nothing
    @test DGGSGrid(rhp, 40).level == 40
    @test_throws OverflowError leaf_count(rhp, 40)
    @test_throws OverflowError leaf_interval(rhp, 0, 1, 40)
    @test_throws OverflowError child_ids(rhp, 0, typemax(Int))
end

@testset "not ported boundary math is explicit" begin
    # The interface fallback: a system with no wired geometry says so rather
    # than guessing. The example has moved twice as systems got wired — `S2DGGS`
    # first (scaffold ordinals, `src/S2/S2Kernel.jl`), then `ISEA9RDGGS`, whose
    # layout question the OGC 21-038r1 Annex B.2 ten-root statement settled
    # (`src/ISEA9R/Isea9rKernel.jl`). What is left are systems that really are
    # registry-only: `RHEALPixDGGS` carries a written disposition of what a
    # wiring would need, `ISEA4TDGGS` is a plain gap.
    @test_throws NotPortedError cell_polygon(RHEALPixDGGS(), 0, 0)
    @test_throws NotPortedError cell_polygon(ISEA4TDGGS(), 0, 0)
    # ...and the systems whose geometry IS wired answer at this same generic.
    for system in (HEALPixDGGS(), S2DGGS(), ISEA4RDGGS(), ISEA9RDGGS())
        @test cell_polygon(system, 2, 3) isa GI.Polygon
    end
end

# The package namespace promises never to shadow a system submodule's native
# vocabulary, so three kernel generics stay qualified-only. This suite is one of
# the collision oracles: it `using`s the package, and the per-system suites
# `using` the package alongside their submodules.
@testset "kernel export surface" begin
    exported = names(DiscreteGlobalGrids)
    @test :cell_polygon_unitsphere in exported
    for name in (:num_cells, :cell_boundary, :cell_center)
        @test !(name in exported)
        # ...and the reason: a submodule owns a *different* function of that name.
        @test getproperty(DiscreteGlobalGrids, name) !==
              getproperty(DiscreteGlobalGrids.H3.H3Native, name)
    end
    # Re-exported from `ConservativeRegridding.Trees`, as that module's own
    # bindings — so a `using ConservativeRegridding` in the same scope (it
    # exports `ncells`/`getcell` too) cannot make either ambiguous.
    for name in (:treeify, :ncells, :getcell)
        @test name in exported
        @test getproperty(DiscreteGlobalGrids, name) ===
              getproperty(DiscreteGlobalGrids.Trees, name)
    end
    @test :intersects_cap in exported
    @test :node_level in exported && :node_id in exported
    # The lookup layer's two core hooks: the id vector a globe-complete
    # dimension is built from, and the supertype the four `<X>Lookup`s share so
    # generic code can dispatch on "a DGGS lookup" at all.
    @test :DGGSGlobeIds in exported
    @test :AbstractDGGSLookup in exported
    for lookup in (DiscreteGlobalGrids.H3.H3Lookups.H3Lookup,
                   DiscreteGlobalGrids.A5.A5Lookups.A5Lookup,
                   DiscreteGlobalGrids.IGeo7.IGeo7Lookups.IGeo7Lookup,
                   DiscreteGlobalGrids.HEALPix.HealpixLookups.HealpixLookup)
        @test lookup <: AbstractDGGSLookup
    end
end

# Operations kernel (src/core/kernel.jl) and the generic tree family
# (src/core/grid_types.jl, src/core/generic_cursor.jl). Each file wraps itself
# in a module of its own, so its mock systems and the deliberately unexported
# kernel vocabulary (`num_cells`, `cell_boundary`, ...) stay out of this
# suite's namespace.
include("test_kernel_core.jl")
include("test_generic_trees.jl")
include("test_manifolds.jl")
# The lazy globe-complete id vector (src/core/globe_ids.jl), and the tree layer
# and point selectors over the lookups built on it (src/core/lookups.jl); also
# modules of their own, for the same `num_cells` reason.
include("test_globe_ids.jl")
include("test_globe_trees.jl")
include("test_globe_selectors.jl")

end # module CoreTestSuite
