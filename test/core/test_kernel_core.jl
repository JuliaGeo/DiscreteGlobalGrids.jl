# test/core/test_kernel_core.jl — operations-kernel defaults (src/core/kernel.jl)
#
# Three self-contained mocks and one registry singleton cover the ways a
# kernel generic can answer:
#
#   * `OrdinalMock` — `has_ordinal_ids = true`, so the whole hierarchy/ordinal/
#     pruning group must be derivable from `root_count`/`radix` alone. Only
#     geometry (`cell_boundary`) is wired, exactly as the HEALPix wiring will.
#   * `UnboundedMock` — the same id model with `max_level === nothing`, which is
#     the other half of the wrap guards: with no level bound to reject, the
#     radix products themselves are what get checked.
#   * `DegenerateMock` — boundary vertices that average to the origin, the one
#     failure mode of the derived `cell_center`.
#   * `RHEALPixDGGS` — a registered system with no kernel wiring at all: every
#     operation must throw `NotPortedError`, never a silent wrong default.
#     (`S2DGGS` played this part until its geometry group was wired.)
#
# Suite lives in its own module: `cell_boundary`/`cell_center`/`num_cells` are
# deliberately unexported from the package (system submodules own those names),
# so the tests qualify them and never pollute the shared test namespace.
module KernelCoreTests

using Test

using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
import GeometryOps as GO
import GeoInterface as GI

# --------------------------------------------------------------------------
# Mock ordinal system: an equirectangular quadtree.
#
# Level 0 is 4 longitude sectors x 2 latitude bands (8 roots); each cell splits
# into 2x2, with child `4id + k` taking the `(k & 1, k >> 1)` half in
# (lon, lat). Cell ids are dense ordinals by construction, which is exactly the
# nested-HEALPix id model the kernel defaults assume.
# --------------------------------------------------------------------------

struct OrdinalMock <: AbstractDGGS end

const ROOT_COUNT = 8
const RADIX = 4

DGG.system_name(::OrdinalMock) = :OrdinalMock
DGG.root_count(::OrdinalMock) = ROOT_COUNT
DGG.radix(::OrdinalMock) = RADIX
DGG.max_level(::OrdinalMock) = 12
DGG.has_ordinal_ids(::OrdinalMock) = true

"Longitude/latitude box (degrees) of cell `(level, id)`."
function cell_box(level::Integer, id::Integer)
    span = RADIX^Int(level)
    root, within = divrem(Int(id), span)
    x = 0
    y = 0
    for bit in 0:(Int(level) - 1)
        x |= ((within >> (2bit)) & 1) << bit
        y |= ((within >> (2bit + 1)) & 1) << bit
    end
    width = 90.0 / (1 << Int(level))
    lon = -180.0 + 90.0 * (root % 4) + x * width
    lat = -90.0 + 90.0 * (root ÷ 4) + y * width
    return (lon, lat, lon + width, lat + width)
end

function unit_point(lon::Real, lat::Real)
    lambda = deg2rad(lon)
    phi = deg2rad(lat)
    cosphi = cos(phi)
    return GO.UnitSphericalPoint(cosphi * cos(lambda), cosphi * sin(lambda), sin(phi))
end

function DGG.cell_boundary(::OrdinalMock, level::Integer, id; closed::Bool=false)
    lon0, lat0, lon1, lat1 = cell_box(level, id)
    points = [unit_point(lon0, lat0), unit_point(lon1, lat0),
        unit_point(lon1, lat1), unit_point(lon0, lat1)]
    closed && push!(points, points[1])
    return points
end

# No `cell_center` method: the derived one (normalized boundary mean) is under
# test here.

# The same grid with no pinned maximum (`RHEALPixDGGS` et al. in the registry).
# Nothing bounds the exponent of a `radix^delta`, so the guards fall through to
# the checked products and an absurd level is an `OverflowError` rather than a
# wrapped count.
struct UnboundedMock <: AbstractDGGS end

DGG.system_name(::UnboundedMock) = :UnboundedMock
DGG.root_count(::UnboundedMock) = ROOT_COUNT
DGG.radix(::UnboundedMock) = RADIX
DGG.max_level(::UnboundedMock) = nothing
DGG.has_ordinal_ids(::UnboundedMock) = true

struct DegenerateMock <: AbstractDGGS end

DGG.system_name(::DegenerateMock) = :DegenerateMock
DGG.cell_boundary(::DegenerateMock, level::Integer, id; closed::Bool=false) =
    [GO.UnitSphericalPoint(1.0, 0.0, 0.0), GO.UnitSphericalPoint(-1.0, 0.0, 0.0),
        GO.UnitSphericalPoint(0.0, 1.0, 0.0), GO.UnitSphericalPoint(0.0, -1.0, 0.0)]

const MOCK = OrdinalMock()
const UNBOUNDED = UnboundedMock()
# A registered system with NO kernel wiring at all. This used to be `S2DGGS`;
# S2 now answers the geometry group over scaffold ordinals
# (`src/S2/S2Kernel.jl`), so the role moved to `RHEALPixDGGS`, which is
# registry-only and recorded as such in its docstring (its face grid is
# expressible but deferred — piecewise cap charts, unpinned square placement,
# ellipsoidal boundaries). It is the same *kind* of witness: traits wired,
# operations not.
const UNWIRED = RHEALPixDGGS()

@testset "id-model traits" begin
    @test DGG.cell_id_type(MOCK) === Int64
    @test has_ordinal_ids(MOCK)
    # `has_descendant_ranges` defaults to `has_ordinal_ids`, so the pruning
    # path comes along for free with the ordinal id model.
    @test has_descendant_ranges(MOCK)

    @test DGG.cell_id_type(UNWIRED) === Int64
    @test !has_ordinal_ids(UNWIRED)
    @test !has_descendant_ranges(UNWIRED)
end

@testset "derived hierarchy" begin
    @test root_ids(MOCK) == collect(Int64, 0:7)
    @test eltype(root_ids(MOCK)) === Int64

    @test cell_children(MOCK, 0, 3) == collect(Int64, 12:15)
    @test cell_children(MOCK, 2, 5) == collect(Int64, 20:23)
    @test issorted(cell_children(MOCK, 4, 1234))
    @test length(cell_children(MOCK, 4, 1234)) == RADIX

    # `id * radix` names the children of an id that is not a cell just as
    # happily as of one that is, so the ordinal default range-checks the root
    # exactly as `descendant_range` does. The generic cursor's dense descent
    # builds the same block inline as a `UnitRange`, so this costs traversals
    # nothing — only direct callers reach here.
    @test cell_children(MOCK, 0, ROOT_COUNT - 1) ==
        collect(Int64, ((ROOT_COUNT - 1) * RADIX):(ROOT_COUNT * RADIX - 1))
    @test_throws ArgumentError cell_children(MOCK, 0, ROOT_COUNT)
    @test_throws ArgumentError cell_children(MOCK, 0, -1)
    @test_throws ArgumentError cell_children(MOCK, 2, DGG.num_cells(MOCK, 2))

    # Children of every level-1 cell partition that cell's level-2 slice.
    for id in 0:(ROOT_COUNT * RADIX - 1)
        children = cell_children(MOCK, 1, id)
        @test all(child -> cell_parent(MOCK, 2, child, 1) == id, children)
    end

    @test cell_parent(MOCK, 3, 500, 3) == 500
    @test cell_parent(MOCK, 3, 500, 1) == 500 ÷ 16
    @test cell_parent(MOCK, 3, 500, 0) == 500 ÷ 64
    @test_throws ArgumentError cell_parent(MOCK, 2, 5, 3)
    @test_throws ArgumentError cell_parent(MOCK, 2, 5, -1)

    # `parent_level` was checked; the id was not, so an id off the end of its
    # own level divided down to a perfectly ordinary ancestor.
    @test cell_parent(MOCK, 2, DGG.num_cells(MOCK, 2) - 1, 0) == ROOT_COUNT - 1
    @test_throws ArgumentError cell_parent(MOCK, 2, DGG.num_cells(MOCK, 2), 0)
    @test_throws ArgumentError cell_parent(MOCK, 2, -1, 0)

    @test cell_descendants(MOCK, 0, 3, 2) == collect(Int64, 48:63)
    @test cell_descendants(MOCK, 2, 7, 2) == Int64[7]
    @test_throws ArgumentError cell_descendants(MOCK, 3, 5, 2)

    for level in 0:4
        @test DGG.num_cells(MOCK, level) == ROOT_COUNT * RADIX^level
        @test DGG.num_cells(MOCK, level) isa Int64
    end

    @test subtree_leaf_count(MOCK, 1, 5, 4) == RADIX^3
    @test subtree_leaf_count(MOCK, 2, 5, 2) == 1

    # Same id guard again — the count is the same radix arithmetic, and sized
    # the subtree of a nonexistent cell without complaint.
    @test_throws ArgumentError subtree_leaf_count(MOCK, 1, DGG.num_cells(MOCK, 1), 4)
    @test_throws ArgumentError subtree_leaf_count(MOCK, 1, -1, 4)

    # And the same wrap `max_level` names: `Int64(RADIX)^40` is 0, so the count
    # of a 40-level subtree used to be *zero cells*. A pinned maximum rejects
    # the level; where none is pinned the product itself is checked.
    @test subtree_leaf_count(MOCK, 0, 0, 12) == RADIX^12
    @test_throws ArgumentError subtree_leaf_count(MOCK, 0, 0, 40)
    @test Int64(RADIX)^40 == 0
    @test max_level(UNBOUNDED) === nothing
    @test subtree_leaf_count(UNBOUNDED, 0, 0, 10) == RADIX^10
    @test_throws OverflowError subtree_leaf_count(UNBOUNDED, 0, 0, 40)

    # Complete levels: the subtree sizes of one level tile the next.
    for level in 0:2, leaf_level in level:4
        total = sum(subtree_leaf_count(MOCK, level, id, leaf_level)
                    for id in 0:(DGG.num_cells(MOCK, level) - 1))
        @test total == DGG.num_cells(MOCK, leaf_level)
    end
end

@testset "derived ordinals" begin
    # Exhaustive monotonicity: ascending ids map to 1:num_cells in order.
    for level in 0:3
        n = DGG.num_cells(MOCK, level)
        @test [cell_to_ordinal(MOCK, level, id) for id in 0:(n - 1)] == collect(1:n)
        @test [ordinal_to_cell(MOCK, level, ord) for ord in 1:n] == collect(Int64, 0:(n - 1))
    end
    @test ordinal_to_cell(MOCK, 4, cell_to_ordinal(MOCK, 4, 999)) == 999
    @test ordinal_to_cell(MOCK, 4, 1) isa Int64

    # Neither direction was range-checked, and both are the *definition* of the
    # dense leaf numbering the trees index by: `Trees.getcell` resolves a leaf
    # through `ordinal_to_cell`, so an ordinal past the end of a level came
    # back as a cell id past the end of that level rather than as an error.
    # Two integer comparisons each, which is what a hot path can afford.
    total = DGG.num_cells(MOCK, 2)
    @test cell_to_ordinal(MOCK, 2, total - 1) == total
    @test_throws ArgumentError cell_to_ordinal(MOCK, 2, total)
    @test_throws ArgumentError cell_to_ordinal(MOCK, 2, -1)

    # An ordinal is an index into the level, so out of range is an
    # `OrdinalRangeError` here exactly as it is in every wired system (A5, H3,
    # IGEO7). A `BoundsError` against `1:num_cells` is uniform too, and says
    # nothing: no system, no level, no word about ordinals.
    @test ordinal_to_cell(MOCK, 2, total) == total - 1
    @test_throws OrdinalRangeError ordinal_to_cell(MOCK, 2, total + 1)
    @test_throws OrdinalRangeError ordinal_to_cell(MOCK, 2, 0)
    @test_throws OrdinalRangeError ordinal_to_cell(MOCK, 0, 10^9)
    # The audit's own case, on the system the ordinal defaults were written
    # for: level 0 is 12 pixels, and this answered id 999999999.
    @test_throws OrdinalRangeError ordinal_to_cell(HEALPixDGGS(), 0, 10^9)

    # ... and what it says, from that same case. The fields are the facts and
    # the sentence is built in `showerror`, so this is also the check that
    # nothing is formatted on the throwing path.
    err = try
        ordinal_to_cell(HEALPixDGGS(), 0, 10^9)
    catch e
        e
    end
    @test err isa OrdinalRangeError
    @test err.system === :HEALPix
    @test err.level == 0
    @test err.ordinal == 10^9
    @test err.total == 12
    message = sprint(showerror, err)
    @test occursin("HEALPix", message)                  # which system
    @test occursin("level 0", message)                  # which level
    @test occursin("1000000000", message)               # the offending ordinal
    @test occursin("1:12", message)                     # the valid range
    @test occursin("1-based position", message)         # what an ordinal is
end

@testset "derived descendant_range" begin
    @test descendant_range(MOCK, 1, 5, 3) == (Int64(80), Int64(95))
    @test descendant_range(MOCK, 2, 7, 2) == (Int64(7), Int64(7))
    @test_throws ArgumentError descendant_range(MOCK, 3, 5, 2)

    # An ordinal id carries no level of its own, so `0:num_cells - 1` is the
    # only thing there is to check it against — and the radix arithmetic below
    # is happy to return a perfectly well-formed interval of cells that do not
    # exist (`subtree_grid` then builds a grid of them). The structural wirings
    # have had this guard from the start, against the resolution their ids
    # encode; this is its ordinal counterpart.
    @test descendant_range(MOCK, 0, ROOT_COUNT - 1, 2) ==
        (Int64((ROOT_COUNT - 1) * RADIX^2), Int64(ROOT_COUNT * RADIX^2 - 1))
    @test_throws ArgumentError descendant_range(MOCK, 0, ROOT_COUNT, 2)
    @test_throws ArgumentError descendant_range(MOCK, 0, -1, 2)
    @test_throws ArgumentError descendant_range(MOCK, 1, DGG.num_cells(MOCK, 1), 3)
    # `cell_descendants` reaches it through the range branch, and carries the
    # same guard on the expansion branch it falls back to.
    @test_throws ArgumentError cell_descendants(MOCK, 0, ROOT_COUNT, 2)
    @test cell_descendants(MOCK, 0, ROOT_COUNT - 1, 1) ==
        collect(Int64, ((ROOT_COUNT - 1) * RADIX):(ROOT_COUNT * RADIX - 1))

    # Two-sided contract, over complete parent levels: endpoints tight, ranges
    # ordered and disjoint, sizes summing to the level.
    for level in 0:2, leaf_level in level:4
        n = DGG.num_cells(MOCK, level)
        previous_hi = nothing
        total = 0
        for id in 0:(n - 1)
            lo, hi = descendant_range(MOCK, level, id, leaf_level)
            descendants = cell_descendants(MOCK, level, id, leaf_level)
            @test extrema(descendants) == (lo, hi)
            @test hi - lo + 1 == length(descendants)
            previous_hi === nothing || @test lo == previous_hi + 1
            previous_hi = hi
            total += length(descendants)
        end
        @test total == DGG.num_cells(MOCK, leaf_level)
        @test previous_hi == DGG.num_cells(MOCK, leaf_level) - 1
    end

    # The cursor's dense-node interval identity: ordinal monotonicity plus the
    # two-sided range contract make `cell_to_ordinal(lo):cell_to_ordinal(hi)`
    # exactly the node's leaf window.
    for level in 0:2, id in 0:(DGG.num_cells(MOCK, level) - 1)
        leaf_level = level + 2
        lo, hi = descendant_range(MOCK, level, id, leaf_level)
        width = cell_to_ordinal(MOCK, leaf_level, hi) - cell_to_ordinal(MOCK, leaf_level, lo) + 1
        @test width == subtree_leaf_count(MOCK, level, id, leaf_level)
    end
end

@testset "derived geometry" begin
    boundary = DGG.cell_boundary(MOCK, 3, 100)
    @test length(boundary) == 4
    @test all(point -> isapprox(sum(abs2, point), 1.0; atol=1e-12), boundary)

    closed = DGG.cell_boundary(MOCK, 3, 100; closed=true)
    @test length(closed) == 5
    @test closed[1] == closed[end]

    # Derived center: normalized boundary mean, on the sphere and inside the
    # cell's cap.
    center = DGG.cell_center(MOCK, 3, 100)
    @test isapprox(sum(abs2, center), 1.0; atol=1e-12)
    @test_throws ArgumentError DGG.cell_center(DegenerateMock(), 0, 0)

    polygon = DGG.cell_polygon_unitsphere(MOCK, 3, 100)
    @test GI.trait(polygon) isa GI.PolygonTrait
    ring = GI.getexterior(polygon)
    @test collect(GI.getpoint(ring)) == closed

    cap = cell_cap(MOCK, 3, 100)
    @test cap isa GO.UnitSpherical.SphericalCap{Float64}
    @test cap.point == center
    @test all(point -> GO.UnitSpherical._contains(cap, point), boundary)
    # The 1.2 inflation is what makes a leaf cap usable where children overhang.
    tight = maximum(GO.UnitSpherical.spherical_distance(center, p) for p in boundary)
    @test cap.radius > tight
    @test cap.radius <= tight * DGG.CELL_CAP_INFLATION + 1e-6
end

@testset "cells_cap edge cases" begin
    full = DGG.full_sphere_extent()

    @test DGG.cells_cap(MOCK, 3, Int64[]) == full
    @test DGG.cells_cap(MOCK, 3, Int64[7]) == cell_cap(MOCK, 3, 7)

    # Past the exact limit the union cap stops paying for itself.
    @test DGG.cells_cap(MOCK, 6, collect(Int64, 0:DGG.SUBTREE_CAP_EXACT_LIMIT)) == full
    @test length(collect(Int64, 0:DGG.SUBTREE_CAP_EXACT_LIMIT)) == DGG.SUBTREE_CAP_EXACT_LIMIT + 1

    # Exact union: every vertex of every cell in the batch is inside.
    batch = collect(Int64, 96:127)
    batch_cap = DGG.cells_cap(MOCK, 3, batch)
    @test batch_cap != full
    for id in batch, point in DGG.cell_boundary(MOCK, 3, id)
        @test GO.UnitSpherical._contains(batch_cap, point)
    end
    # A view (what the generic partial cursor hands it) behaves the same.
    @test DGG.cells_cap(MOCK, 3, view(batch, 1:length(batch))) == batch_cap
end

@testset "subtree_cap edge cases" begin
    # One answer at every depth — the cell's own cap, O(1) in the node. The
    # subtree size no longer selects between an exact union and the cell cap:
    # the union's O(subtree) boundary calls fell on the small nodes, which is
    # where a traversal spends its time (see the docstring).
    @test DGG.subtree_cap(MOCK, 3, 100, 3) == cell_cap(MOCK, 3, 100)
    @test subtree_leaf_count(MOCK, 0, 2, 5) <= DGG.SUBTREE_CAP_EXACT_LIMIT
    @test DGG.subtree_cap(MOCK, 0, 2, 5) == cell_cap(MOCK, 0, 2)
    @test subtree_leaf_count(MOCK, 0, 2, 6) > DGG.SUBTREE_CAP_EXACT_LIMIT
    @test DGG.subtree_cap(MOCK, 0, 2, 6) == cell_cap(MOCK, 0, 2)

    # ...and it still contains what it claims to: every vertex of every leaf of
    # the subtree, which is the whole contract a traversal rests on.
    cap = DGG.subtree_cap(MOCK, 0, 2, 5)
    @test all(cell_descendants(MOCK, 0, 2, 5)) do id
        all(point -> GO.UnitSpherical._contains(cap, point), DGG.cell_boundary(MOCK, 5, id))
    end

    # The guards did not go with the union cap: a cap is still refused for a
    # subtree that does not exist. `subtree_leaf_count` is what raises them,
    # and it is still asked.
    @test_throws ArgumentError DGG.subtree_cap(MOCK, 1, DGG.num_cells(MOCK, 1), 4)
    @test_throws ArgumentError DGG.subtree_cap(MOCK, 1, -1, 4)
    @test_throws ArgumentError DGG.subtree_cap(MOCK, 0, 0, 40)   # past max_level
    @test_throws OverflowError DGG.subtree_cap(UNBOUNDED, 0, 0, 40)
end

@testset "NotPortedError paths" begin
    # Every kernel generic on an unwired registered system: no silent defaults.
    @test_throws NotPortedError root_ids(UNWIRED)
    @test_throws NotPortedError cell_children(UNWIRED, 0, 0)
    @test_throws NotPortedError cell_parent(UNWIRED, 1, 0, 0)
    @test_throws NotPortedError cell_descendants(UNWIRED, 0, 0, 1)
    @test_throws NotPortedError DGG.num_cells(UNWIRED, 0)
    @test_throws NotPortedError subtree_leaf_count(UNWIRED, 0, 0, 1)
    @test_throws NotPortedError cell_to_ordinal(UNWIRED, 0, 0)
    @test_throws NotPortedError ordinal_to_cell(UNWIRED, 0, 1)
    @test_throws NotPortedError descendant_range(UNWIRED, 0, 0, 1)
    @test_throws NotPortedError DGG.cell_boundary(UNWIRED, 0, 0)
    @test_throws NotPortedError DGG.cell_center(UNWIRED, 0, 0)
    @test_throws NotPortedError DGG.cell_polygon_unitsphere(UNWIRED, 0, 0)
    @test_throws NotPortedError cell_cap(UNWIRED, 0, 0)
    @test_throws NotPortedError DGG.cells_cap(UNWIRED, 0, [0, 1])
    @test_throws NotPortedError DGG.subtree_cap(UNWIRED, 0, 0, 1)

    # The error names the system and the operation, and renders lazily.
    err = try
        root_ids(UNWIRED)
    catch caught
        caught
    end
    @test err isa NotPortedError
    @test err.system === :rHEALPix
    @test err.operation === :root_ids
    message = sprint(showerror, err)
    @test occursin("root_ids", message)
    @test occursin("rHEALPix", message)

    # An unwired *geometry* still short-circuits where the kernel can answer
    # without it: a zero-length batch never touches `cell_boundary`.
    @test DGG.cells_cap(UNWIRED, 0, Int64[]) == DGG.full_sphere_extent()
end

end # module KernelCoreTests
