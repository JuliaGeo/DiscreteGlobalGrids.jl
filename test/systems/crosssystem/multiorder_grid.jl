# ---------------------------------------------------------------------------
# `MultiOrderGrid`: the grid face of a mixed-level container's STORED cells, and
# the regridding route that reads it.
#
#   * counts and ids are the container's — never the reference level's leaves,
#     which is the entire reason the type exists;
#   * geometry is forwarded per cell to its own level's grid, never derived, so
#     a system that overrides area on a level grid keeps that answer;
#   * point location is the covering ancestor, which is what makes a nearest
#     regrid off the stored cells BIT-IDENTICAL to one off the expansion;
#   * the tiling does not conform, so every topology verb throws.
#
# Swept on HEALPix (radix 4, quadrilateral) and IGeo7 (radix 7, hexagonal) —
# the two families the non-congruence argument distinguishes.
# ---------------------------------------------------------------------------

module MultiOrderGridTests

using Test
import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import DimensionalData as DD
import GeometryOps as GO

const EN = DGG.Engine

# A method declaring a sampling a mixed container has no presentation for:
# neither the stored cells (sample sites) nor the expansion (a gap-free cover).
struct OddSampling <: DD.Lookups.Sampling end
struct ThirdSampling <: GR.AbstractRegriddingMethod end
GR.sourcesampling(::ThirdSampling) = OddSampling()

# A container spanning three levels over the whole sphere: every root at level
# `top + 1`, with the first refined two levels deeper and the second one.
function mixed(sys, top::Int)
    roots = collect(DGG.CellVector(DGG.levelgrid(sys, top + 1)))
    kids(c, l) = collect(DGG.CellVector(DGG.subtree(sys, c, l)))
    cells = vcat(kids(roots[1], top + 3), kids(roots[2], top + 2), roots[3:end])
    return DGG.MultiOrderVector(sys, cells; reference_level = top + 3)
end

# Typed wrappers: inference through a closure over an untyped global reports a
# false instability, so every `@inferred` below goes through one of these.
grid_of(mov::DGG.MultiOrderVector) = DGG.MultiOrderGrid(mov)
count_of(g::DGG.MultiOrderGrid) = DGG.ncells(g)
id_at(g::DGG.MultiOrderGrid, i::Int) = DGG.cellindex(g, i)
centroid_of(g::DGG.MultiOrderGrid, c) = DGG.cell_centroid(g, c)
boundary_of(g::DGG.MultiOrderGrid, c) = DGG.cell_boundary(g, c)
area_of(g::DGG.MultiOrderGrid, c) = DGG.cell_area(g, c)
cap_of(g::DGG.MultiOrderGrid, c) = DGG.Fallbacks.cell_cap(g, c)
tree_of(g::DGG.MultiOrderGrid) = DGG.treeify(g)
index_of(g::DGG.MultiOrderGrid, c::DGG.AbstractCellIndex) = DGG.localindex(g, c)
index_at(g::DGG.MultiOrderGrid, p::GO.UnitSphericalPoint) = DGG.localindex(g, p)
space_for(mov::DGG.MultiOrderVector, m::GR.AbstractRegriddingMethod) =
    GR.sourcespacefor(mov, m)

# `congruent` says whether a parent's polygon IS the union of its descendants':
# true for HEALPix's nested quadrilaterals, false for IGeo7's hexagons, where a
# parent and its children cover slightly different ground.
const CASES = ((DGG.HEALPixSystem(), 0, true), (DGG.IGeo7System(), 0, false))

@testset "MultiOrderGrid on $(nameof(typeof(sys)))" for (sys, top, congruent) in CASES

    mov = mixed(sys, top)
    g = DGG.MultiOrderGrid(mov)
    ref = DGG.reference_level(mov)
    leaves = DGG.ncells(DGG.levelgrid(sys, ref))

    @testset "the cells are the stored cells" begin
        # The one thing the type exists for. A grid reporting the leaf count is
        # the mutant every column of a weight matrix would then be sized by.
        @test length(unique(DGG.level, collect(mov))) == 3
        @test DGG.ncells(g) == length(mov) < leaves
        @test all(i -> DGG.cellindex(g, i) == mov[i], 1:DGG.ncells(g))
        @test DGG.system(g) === DGG.system(mov)
        @test DGG.cellset(g) === mov
        # Several levels, so there is no one level — which is what buys the
        # single-chunk fallback in `DGGSpace` with no change there.
        @test DGG.level(g) === nothing
        @test DGG.ncells(DGG.levelgrid(sys, DGG.level(mov[1]))) != DGG.ncells(g)
        @test_throws BoundsError DGG.cellindex(g, DGG.ncells(g) + 1)
        @test occursin("MultiOrderGrid", sprint(show, g))
    end

    @testset "geometry is forwarded to each cell's own level" begin
        # Not derived. HEALPix and ISEA4R override area on their level grid with
        # the exact `4pi/ncells`; a ring-polygon area would be silently wrong on
        # their curvilinear edges, and this is the test that sees it.
        for i in 1:DGG.ncells(g)
            c = DGG.cellindex(g, i)
            own = DGG.levelgrid(sys, DGG.level(c))
            @test DGG.cell_area(g, c) == DGG.cell_area(own, c)
            @test DGG.cell_centroid(g, c) == DGG.cell_centroid(own, c)
            @test collect(DGG.cell_boundary(g, c)) == collect(DGG.cell_boundary(own, c))
            @test DGG.Fallbacks.cell_cap(g, c) == DGG.Fallbacks.cell_cap(own, c)
        end
        # Cells at different levels really do have different areas here, so the
        # forwarding is doing work rather than agreeing by accident.
        areas = [DGG.cell_area(g, DGG.cellindex(g, i)) for i in 1:DGG.ncells(g)]
        @test length(unique(round.(areas; digits = 12))) >= 3
        # `getcell` is the local-index form of the same polygon.
        @test DGG.getcell(g, 1) == DGG.cell_polygon(g, DGG.cellindex(g, 1))

        # Whether the stored cells COVER the sphere is the whole area-method
        # argument, and it is a property of the hierarchy, not of this grid.
        # Expanding to the reference level always tiles; the stored cells tile
        # only where a parent's polygon is the union of its descendants'.
        expanded = sum(DGG.cell_area(DGG.levelgrid(sys, ref), c)
                       for c in DGG.CellVector(mov))
        @test expanded ≈ 4pi rtol = 1e-9
        if congruent
            @test sum(areas) ≈ 4pi rtol = 1e-9
        else
            # The gap Route A exists to avoid: coarse polygons over a
            # non-congruent hierarchy are not a partition, so an area method
            # weighted against them would lose or double-count mass.
            @test !isapprox(sum(areas), 4pi; rtol = 1e-6)
            @test isapprox(sum(areas), 4pi; rtol = 1e-2)
        end
    end

    @testset "the tree is level-agnostic" begin
        # The generic sends any grid with a system to `HierarchicalGridCursor`,
        # whose index windows assume one leaf level. That mutant returns a
        # cursor here and mislabels every leaf.
        t = DGG.treeify(g)
        @test t isa EN.IndexTreeNode
        @test !(t isa EN.HierarchicalGridCursor)
        @test DGG.treeify(t) === t
    end

    @testset "point location is the covering ancestor" begin
        # Exact membership stays exact...
        @test all(i -> DGG.localindex(g, mov[i]) == i, 1:length(mov))
        @test all(i -> mov[i] in g, 1:length(mov))
        # ...and a leaf the container does not store is NOT a member.
        level_ref = DGG.levelgrid(sys, ref)
        deep = DGG.cellindex(level_ref, 1)
        coarse_host = DGG.covering_index(mov, deep)
        @test DGG.localindex(g, deep) === (mov[coarse_host] == deep ?
                                           coarse_host : nothing)

        # But a POINT inside a coarse cell resolves to that coarse cell. A
        # mutant resolving points by exact membership drops every point that
        # lands under a stored cell coarser than the reference level.
        for i in 1:min(200, DGG.ncells(level_ref))
            leaf = DGG.cellindex(level_ref, i)
            p = DGG.cell_centroid(level_ref, leaf)
            @test DGG.localindex(g, p) == DGG.covering_index(mov, leaf)
            @test DGG.cellat(g, p) == mov[DGG.covering_index(mov, leaf)]
        end
        # A container covering the whole sphere leaves no point unmapped.
        @test DGG.localindex(g, DGG.cell_centroid(g, mov[1])) == 1
    end

    @testset "a non-conforming tiling refuses the topology verbs" begin
        # A ring matched on shared vertices is wrong across every T-junction,
        # and on a hexagonal hierarchy parent and child vertices never coincide.
        # A plausible wrong ring is the mutant; a refusal naming the working
        # verb is the answer.
        c = mov[1]
        for call in (() -> DGG.neighbors(g, c), () -> DGG.neighbors(g, c, 2),
                     () -> DGG.neighbors(g, 1), () -> DGG.ring(g, c, 1),
                     () -> DGG.one_ring(g, c), () -> DGG.halo(g),
                     () -> DGG.border(g), () -> DGG.interior(g),
                     () -> DGG.adjacency(g))
            @test_throws ArgumentError call()
            @test_throws "member_neighbors" call()
        end
        # No one level, so no index space carved out of one either.
        @test_throws ArgumentError DGG.globalindex(g, c)
        @test_throws "no global index" DGG.globalindex(g, c)
    end

    @testset "the source space is the stored cells, for point methods only" begin
        for m in (GR.NearestCell(), GR.DirectNearest(), GR.BarycentricPoint())
            space = GR.sourcespacefor(mov, m)
            @test space.grid isa DGG.MultiOrderGrid
            @test DGG.ncells(space) == length(mov)
            # One chunk, over every stored cell, with no `DGGSpace` change.
            @test GR.nchunks(space) == 1
            @test GR.ownedindices(space, 1) == 1:length(mov)
            @test GR.chunkranges(space, 1, (length(mov),)) == (1:length(mov),)
            @test GR.celltree(space) isa EN.IndexTreeNode
            @test GR.cellat(space, DGG.cell_centroid(g, mov[1])) == 1
        end
        # An area method keeps the expansion: a stored cell's descendant leaves
        # are the only gap-free cover of it on a non-congruent hierarchy, and
        # coarse polygons would leave slivers.
        area = GR.sourcespacefor(mov, GR.Conservative())
        @test area.grid isa DGG.PartialGrid
        @test DGG.ncells(area) == leaves

        # A container storing one cell per leaf expands to itself, so it stays
        # on the already-tested subset path whatever the method reads.
        uniform = DGG.MultiOrderVector(sys,
            collect(DGG.CellVector(DGG.levelgrid(sys, top + 1)));
            reference_level = top + 1)
        for m in (GR.NearestCell(), GR.Conservative(), GR.BarycentricPoint())
            u = GR.sourcespacefor(uniform, m)
            @test u.grid isa DGG.PartialGrid
            @test DGG.ncells(u) == length(uniform)
        end

        # A method declaring a third sampling has no presentation to pick.
        @test_throws ArgumentError GR.sourcespacefor(mov, ThirdSampling())
        @test_throws "neither `Points()` nor `Intervals()`" GR.sourcespacefor(mov,
            ThirdSampling())
    end

    @testset "nearest off the stored cells is nearest off the expansion" begin
        # THE equality. Route A locates the destination centroid in the
        # reference-level grid and takes that leaf's replicated value; Route B
        # locates the same point at the same level and takes the covering
        # stored cell's value. Same level, same call, same tie rule on a shared
        # boundary — so it is `isequal`, not `≈`, and every mutant in covering
        # resolution moves at least one number.
        vals = collect(1.0:length(mov))       # distinct per stored cell
        cube = DD.DimArray(vals, DGG.Cells(DGG.MultiOrderLookup(mov)))
        expanded = DGG.expand(cube, ref)

        for dst in (DGG.levelgrid(sys, top + 1), DGG.levelgrid(sys, top + 2))
            for m in (GR.NearestCell(), GR.DirectNearest())
                native = DGG.regrid(cube; to = dst, method = m)
                routeA = DGG.regrid(expanded; to = dst, method = m)
                @test isequal(parent(native), parent(routeA))
                # And the two point methods agree with each other.
                @test isequal(parent(native),
                    parent(DGG.regrid(cube; to = dst,
                        method = m === GR.NearestCell() ? GR.DirectNearest() :
                                 GR.NearestCell())))
                # The plan really did read the stored cells.
                plan = GR.plan_regrid(cube; to = dst, method = m, lazy = false)
                @test DGG.ncells(plan.src_space) == length(mov)
            end
        end

        # Two chunkings of the destination, same answer: the source is one
        # chunk either way, and the destination tiling must not move a value.
        dst = DGG.levelgrid(sys, top + 2)
        base = DGG.regrid(cube; to = dst, method = GR.NearestCell())
        for cells in (8, 4096)
            chunked = DGG.DGGSpace(dst; chunkcells = cells)
            @test isequal(parent(DGG.regrid(cube; to = chunked,
                    method = GR.NearestCell())), parent(base))
        end

        # The weight matrix has one column per STORED cell, not per leaf. This
        # is the shape the perf claim rests on, asserted rather than measured.
        plan = GR.plan_regrid(cube; to = dst, method = GR.NearestCell(),
            lazy = false)
        @test size(plan.block.weights, 2) == length(mov)
        @test size(GR.plan_regrid(expanded; to = dst, method = GR.NearestCell(),
            lazy = false).block.weights, 2) == leaves
    end

    @testset "a hole in the container is unmapped, never nearest" begin
        # Drop one stored cell. Points under it are covered by nothing, and a
        # MOC's holes are meaningful: no stored cell stands in.
        holed = mov[1:length(mov)-1]
        gap = mov[end]
        gh = DGG.MultiOrderGrid(holed)
        p = DGG.cell_centroid(DGG.levelgrid(sys, DGG.level(gap)), gap)
        @test DGG.localindex(gh, p) === nothing

        vals = collect(1.0:length(holed))
        cube = DD.DimArray(vals, DGG.Cells(DGG.MultiOrderLookup(holed)))
        dst = DGG.levelgrid(sys, DGG.level(gap))
        j = DGG.localindex(dst, gap)
        blanked = DGG.regrid(cube; to = dst, method = GR.NearestCell(),
            missingpolicy = GR.Weighted())
        zeroed = DGG.regrid(cube; to = dst, method = GR.NearestCell(),
            missingpolicy = GR.Extensive())
        @test ismissing(parent(blanked)[j]) || isnan(parent(blanked)[j])
        @test parent(zeroed)[j] == 0
        # Everywhere else still answers.
        @test count(x -> !(ismissing(x) || isnan(x)), parent(blanked)) ==
              DGG.ncells(dst) - 1
    end

    @testset "every new entry point infers" begin
        # `juliac` ambition: no dynamic dispatch introduced by the new hooks.
        c = mov[1]
        p = DGG.cell_centroid(g, c)
        @test (@inferred grid_of(mov)) isa DGG.MultiOrderGrid
        @test (@inferred count_of(g)) isa Int
        @test (@inferred id_at(g, 1)) isa DGG.AbstractCellIndex
        @test (@inferred centroid_of(g, c)) isa GO.UnitSphericalPoint
        @test (@inferred area_of(g, c)) isa Float64
        @inferred boundary_of(g, c)
        @inferred cap_of(g, c)
        @test (@inferred tree_of(g)) isa EN.IndexTreeNode
        @test (@inferred Union{Int,Nothing} index_of(g, c)) == 1
        @test (@inferred Union{Int,Nothing} index_at(g, p)) == 1

        # Routing is the one place two space types meet, and the compiler sees
        # exactly those two — never `Any`, so the split is resolved once at plan
        # construction and everything inside a built plan stays concrete.
        for M in (GR.NearestCell, GR.DirectNearest, GR.BarycentricPoint,
                  GR.Conservative)
            rt = Base.infer_return_type(space_for, Tuple{typeof(mov),M})
            @test rt <: DGG.DGGSpace
        end
        # An area method has no choice to make, so its branch is one type.
        @test isconcretetype(Base.infer_return_type(space_for,
            Tuple{typeof(mov),GR.Conservative}))
    end
end

end # module MultiOrderGridTests
