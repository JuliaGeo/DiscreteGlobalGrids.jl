# HEALPix and IGeo7 distinguish congruent from non-congruent stored-cell geometry.

module MultiOrderGridTests

using Test
import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import DimensionalData as DD
import GeometryOps as GO

const EN = DGG.Engine

# A third sampling declaration provides no valid stored or expanded presentation.
struct OddSampling <: DD.Lookups.Sampling end
struct ThirdSampling <: GR.AbstractRegriddingMethod end
GR.sourcesampling(::ThirdSampling) = OddSampling()

# Three stored levels expose accidental reference-level flattening.
function mixed(sys, top::Int)
    roots = collect(DGG.CellVector(DGG.levelgrid(sys, top + 1)))
    kids(c, l) = collect(DGG.CellVector(DGG.subtree(sys, c, l)))
    cells = vcat(kids(roots[1], top + 3), kids(roots[2], top + 2), roots[3:end])
    return DGG.MultiOrderVector(sys, cells; reference_level = top + 3)
end

# Congruence determines whether stored parent polygons tile like their descendants.
const CASES = ((DGG.HEALPixSystem(), 0, true), (DGG.IGeo7System(), 0, false))

@testset "MultiOrderGrid on $(nameof(typeof(sys)))" for (sys, top, congruent) in CASES

    mov = mixed(sys, top)
    g = DGG.MultiOrderGrid(mov)
    ref = DGG.reference_level(mov)
    leaves = DGG.ncells(DGG.levelgrid(sys, ref))

    @testset "the cells are the stored cells" begin
        # Weight-matrix columns are sized by `ncells`: stored cells, not leaves.
        @test length(unique(DGG.level, collect(mov))) == 3
        @test DGG.ncells(g) == length(mov) < leaves
        @test all(i -> DGG.cellindex(g, i) == mov[i], 1:DGG.ncells(g))
        @test DGG.system(g) === DGG.system(mov)
        @test DGG.cellset(g) === mov
        # A `nothing` level selects `DGGSpace`'s single-chunk fallback.
        @test DGG.level(g) === nothing
        @test DGG.ncells(DGG.levelgrid(sys, DGG.level(mov[1]))) != DGG.ncells(g)
        @test_throws BoundsError DGG.cellindex(g, DGG.ncells(g) + 1)
        @test occursin("MultiOrderGrid", sprint(show, g))
    end

    @testset "geometry is forwarded to each cell's own level" begin
        # HEALPix and ISEA4R level grids return the exact `4pi/ncells`; a
        # ring-polygon area is wrong on their curvilinear edges.
        for i in 1:DGG.ncells(g)
            c = DGG.cellindex(g, i)
            own = DGG.levelgrid(sys, DGG.level(c))
            @test DGG.cell_area(g, c) == DGG.cell_area(own, c)
            @test DGG.cell_centroid(g, c) == DGG.cell_centroid(own, c)
            @test collect(DGG.cell_boundary(g, c)) == collect(DGG.cell_boundary(own, c))
            @test DGG.Fallbacks.cell_cap(g, c) == DGG.Fallbacks.cell_cap(own, c)
        end
        # Distinct level areas make forwarding observable.
        areas = [DGG.cell_area(g, DGG.cellindex(g, i)) for i in 1:DGG.ncells(g)]
        @test length(unique(round.(areas; digits = 12))) >= 3
        # `getcell` is the local-index form of the same polygon.
        @test DGG.getcell(g, 1) == DGG.cell_polygon(g, DGG.cellindex(g, 1))

        # The reference-level expansion always tiles the sphere; the stored
        # cells tile only where a parent's polygon is the union of its descendants'.
        expanded = sum(DGG.cell_area(DGG.levelgrid(sys, ref), c)
                       for c in DGG.CellVector(mov))
        @test expanded ≈ 4pi rtol = 1e-9
        if congruent
            @test sum(areas) ≈ 4pi rtol = 1e-9
        else
            # Non-congruent coarse polygons leave gaps or overlaps in area weights.
            @test !isapprox(sum(areas), 4pi; rtol = 1e-6)
            @test isapprox(sum(areas), 4pi; rtol = 1e-2)
        end
    end

    @testset "the tree is level-agnostic" begin
        # `HierarchicalGridCursor`'s index windows assume one leaf level.
        t = DGG.treeify(g)
        @test t isa EN.IndexTreeNode
        @test !(t isa EN.HierarchicalGridCursor)
        @test DGG.treeify(t) === t
    end

    @testset "point location is the covering ancestor" begin
        @test all(i -> DGG.localindex(g, mov[i]) == i, 1:length(mov))
        @test all(i -> mov[i] in g, 1:length(mov))
        # Exact membership excludes unstored covered leaves.
        level_ref = DGG.levelgrid(sys, ref)
        deep = DGG.cellindex(level_ref, 1)
        coarse_host = DGG.covering_index(mov, deep)
        @test DGG.localindex(g, deep) === (mov[coarse_host] == deep ?
                                           coarse_host : nothing)

        # A point inside a coarse cell resolves to that coarse cell.
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
        # A vertex-matched ring is wrong across every T-junction, so the verbs
        # throw and name `member_neighbors`.
        c = mov[1]
        for call in (() -> DGG.neighbors(g, c), () -> DGG.neighbors(g, c, 2),
                     () -> DGG.neighbors(g, 1), () -> DGG.ring(g, c, 1),
                     () -> DGG.one_ring(g, c), () -> DGG.halo(g),
                     () -> DGG.border(g), () -> DGG.interior(g),
                     () -> DGG.adjacency(g))
            @test_throws ArgumentError call()
            @test_throws "member_neighbors" call()
        end
        @test_throws ArgumentError DGG.globalindex(g, c)
        @test_throws "no global index" DGG.globalindex(g, c)
    end

    @testset "the source space is the stored cells, for point methods only" begin
        for m in (GR.NearestCell(), GR.DirectNearest(), GR.BarycentricPoint())
            space = GR.sourcespacefor(mov, m)
            @test space.grid isa DGG.MultiOrderGrid
            @test DGG.ncells(space) == length(mov)
            # One chunk over every stored cell.
            @test GR.nchunks(space) == 1
            @test GR.ownedindices(space, 1) == 1:length(mov)
            @test GR.chunkranges(space, 1, (length(mov),)) == (1:length(mov),)
            @test GR.celltree(space) isa EN.IndexTreeNode
            @test GR.cellat(space, DGG.cell_centroid(g, mov[1])) == 1
        end
        # An area method keeps the expansion: descendant leaves are the only
        # gap-free cover of a stored cell on a non-congruent hierarchy.
        area = GR.sourcespacefor(mov, GR.Conservative())
        @test area.grid isa DGG.PartialGrid
        @test DGG.ncells(area) == leaves

        # A one-cell-per-leaf container expands to itself and takes the subset path.
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
        @test_throws "unsupported source sampling" GR.sourcespacefor(mov,
            ThirdSampling())
    end

    @testset "nearest off the stored cells is nearest off the expansion" begin
        # Both routes locate the destination centroid at the reference level with
        # the same tie rule: one reads the replicated leaf value, the other the
        # covering stored cell's. The results are `isequal`.
        vals = collect(1.0:length(mov))       # distinct per stored cell
        cube = DD.DimArray(vals, DGG.Cells(DGG.MultiOrderLookup(mov)))
        expanded = DGG.expand(cube, ref)

        for dst in (DGG.levelgrid(sys, top + 1), DGG.levelgrid(sys, top + 2))
            for m in (GR.NearestCell(), GR.DirectNearest())
                native = DGG.regrid(cube; to = dst, method = m)
                routeA = DGG.regrid(expanded; to = dst, method = m)
                @test isequal(parent(native), parent(routeA))
                # The two point methods agree.
                @test isequal(parent(native),
                    parent(DGG.regrid(cube; to = dst,
                        method = m === GR.NearestCell() ? GR.DirectNearest() :
                                 GR.NearestCell())))
                plan = GR.plan_regrid(cube; to = dst, method = m, lazy = false)
                @test DGG.ncells(plan.src_space) == length(mov)
            end
        end

        # Destination chunking leaves every value unchanged.
        dst = DGG.levelgrid(sys, top + 2)
        base = DGG.regrid(cube; to = dst, method = GR.NearestCell())
        for cells in (8, 4096)
            chunked = DGG.DGGSpace(dst; chunkcells = cells)
            @test isequal(parent(DGG.regrid(cube; to = chunked,
                    method = GR.NearestCell())), parent(base))
        end

        # Stored-cell columns keep matrix width proportional to compact storage.
        plan = GR.plan_regrid(cube; to = dst, method = GR.NearestCell(),
            lazy = false)
        @test size(plan.block.weights, 2) == length(mov)
        @test size(GR.plan_regrid(expanded; to = dst, method = GR.NearestCell(),
            lazy = false).block.weights, 2) == leaves
    end

    @testset "a hole in the container is unmapped, never nearest" begin
        # Points under a dropped stored cell are covered by nothing.
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
        @test count(x -> !(ismissing(x) || isnan(x)), parent(blanked)) ==
              DGG.ncells(dst) - 1
    end

    @testset "every new entry point infers" begin
        c = mov[1]
        p = DGG.cell_centroid(g, c)
        @test (@inferred DGG.MultiOrderGrid(mov)) isa DGG.MultiOrderGrid
        @test (@inferred DGG.ncells(g)) isa Int
        @test (@inferred DGG.cellindex(g, 1)) isa DGG.AbstractCellIndex
        @test (@inferred DGG.cell_centroid(g, c)) isa GO.UnitSphericalPoint
        @test (@inferred DGG.cell_area(g, c)) isa Float64
        @inferred DGG.cell_boundary(g, c)
        @inferred DGG.Fallbacks.cell_cap(g, c)
        @test (@inferred DGG.treeify(g)) isa EN.IndexTreeNode
        @test (@inferred Union{Int,Nothing} DGG.localindex(g, c)) == 1
        @test (@inferred Union{Int,Nothing} DGG.localindex(g, p)) == 1

        # Routing returns one of two space types, never `Any`.
        for M in (GR.NearestCell, GR.DirectNearest, GR.BarycentricPoint,
                  GR.Conservative)
            rt = Base.infer_return_type(GR.sourcespacefor, Tuple{typeof(mov),M})
            @test rt <: DGG.DGGSpace
        end
        @test isconcretetype(Base.infer_return_type(GR.sourcespacefor,
            Tuple{typeof(mov),GR.Conservative}))
    end
end

end # module MultiOrderGridTests
