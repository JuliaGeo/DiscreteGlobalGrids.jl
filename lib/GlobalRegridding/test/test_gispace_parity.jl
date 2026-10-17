# `GISpace` parity: a GridInterface grid wrapped as a space regrids exactly as
# the space it stands in for, method by method, eager and chunked.
#
# Three pairs. A `ToyLonLatSpace` against a `Structured` twin that shares its
# sample sites; a `RasterGrid` against itself read as a GridInterface grid; and
# the two DGG level grids against their `DGGSpace`s. The weights and the regrid
# outputs must agree; the chunk relations may differ in false positives, so the
# chunked comparison assembles every pair the plan's own relation names.

import DiscreteGlobalGrids as DGG
import GridInterface as GrI
import GridInterface.Structured as GrS
import GeometryOpsCore as GOCore
import DimensionalData as DD
import SparseArrays
using SparseArrays: rowvals, nonzeros, nzrange
import Random

isdefined(@__MODULE__, :ChunkGraphOracles) || include("graphoracles.jl")
using .ChunkGraphOracles: contributing_pairs, graph_pairs, demanded_pairs,
    capjoin_pairs, eager_pairs

# --- the Structured twin of a toy space ------------------------------------
#
# A one-face structured grid whose vertices, chart and sample sites are the
# toy's own formulas, so the two spaces name the same points bit for bit. The
# seams (a periodic longitude, collapsed poles) are a `RegularGrid`'s over the
# same lattice; everything else — rings, extents, the quadtree, tiles, point
# location — is derived by `GridInterface.Structured`.
struct TwinLonLatGrid{R} <: GrS.AbstractStructuredGrid
    toy::ToyLonLatSpace
    reg::R
end

function TwinLonLatGrid(toy::ToyLonLatSpace)
    lons = range(toy.lon0, toy.lon1; length = toy.nlon + 1)
    lats = range(toy.lat0, toy.lat1; length = toy.nlat + 1)
    reg = GrS.RegularGrid(GOCore.Spherical(; radius = 1.0), lons, lats;
        periodic = (toy.lon1 - toy.lon0 >= 360, false),
        collapsed = (toy.lat0 <= -90, toy.lat1 >= 90))
    return TwinLonLatGrid(toy, reg)
end

GOCore.best_manifold(::TwinLonLatGrid) = GOCore.Spherical(; radius = 1.0)
GrS.nfaces(::TwinLonLatGrid) = 1
GrS.facesize(g::TwinLonLatGrid, f) = (g.toy.nlon, g.toy.nlat)
GrS.facevertex(g::TwinLonLatGrid, f, a, b) =
    toy_point(g.toy.lon0 + (a - 1) * dlon(g.toy), g.toy.lat0 + (b - 1) * dlat(g.toy))
GrS.faceposition(g::TwinLonLatGrid, f, i, j) = localindex(g.toy, i, j)
GrS.facecell(g::TwinLonLatGrid, k) = (1, cellsubscript(g.toy, Int(k))...)
GrS.facecross(g::TwinLonLatGrid, f, i, j, s) = GrS.facecross(g.reg, f, i, j, s)
GrS.facefan(g::TwinLonLatGrid, f, a, b) = GrS.facefan(g.reg, f, a, b)
GrI.maxneighbors(g::TwinLonLatGrid, c::GrI.Vertex) = GrI.maxneighbors(g.reg, c)
GrS.hasfaceaxes(::TwinLonLatGrid, f) = true
GrS.faceaxes(g::TwinLonLatGrid, f) =
    ([g.toy.lon0 + (a - 1) * dlon(g.toy) for a in 1:(g.toy.nlon+1)],
        [g.toy.lat0 + (b - 1) * dlat(g.toy) for b in 1:(g.toy.nlat+1)])
GrS.facecoords(g::TwinLonLatGrid, f, p) = GR.chartcoords(g.toy, p)
GrS.faceperiod(::TwinLonLatGrid, f) = (360.0, nothing)
GrI.cell_centroid(g::TwinLonLatGrid, i::Integer) = cellcentroid(g.toy, Int(i))
GrI.chunks(g::TwinLonLatGrid) = GrI.chunks(g, (g.toy.chunklon, g.toy.chunklat))

# --- helpers ------------------------------------------------------------------

# A raster whose chunk grid is out of step with the quadtree's splits.
function parityraster(nx, ny, cx, cy)
    dlon, dlat = 360 / nx, 180 / ny
    array = DD.DimArray(zeros(nx, ny),
        (DD.X(range(-180 + dlon / 2, 180 - dlon / 2; length = nx)),
         DD.Y(range(-90 + dlat / 2, 90 - dlat / 2; length = ny))))
    return RasterGrid(array;
        chunks = ([lo:min(lo + cx - 1, nx) for lo in 1:cx:nx],
            [lo:min(lo + cy - 1, ny) for lo in 1:cy:ny]))
end

# The source field a regrid reads: a plain vector over the source cells, or
# the dimensional raster a `RasterGrid` was built from.
sourcedata(src) = (Random.seed!(7); rand(Int(ncells(src))))
sourcedata(src::ToyLonLatSpace) = (Random.seed!(7); rand(src.nlon, src.nlat))
sourcedata(src::GISpace{<:TwinLonLatGrid}) = sourcedata(src.grid.toy)
sourcedata(src::RasterGrid) = (Random.seed!(7);
    DD.DimArray(rand(GR._nx(src), GR._ny(src)), (src.xdim, src.ydim)))
sourcedata(src::GISpace{<:RasterGrid}) = sourcedata(src.grid)

values_of(x::DD.AbstractDimArray) = vec(collect(parent(x)))
values_of(x) = vec(collect(x))

# The whole operator a chunked plan assembles: every block of every pair its
# relation names, in the spaces' own local indices, with the per-destination
# denominators summed in ascending source-chunk order.
function assembled(plan)
    dst, src = plan.dst_space, plan.src_space
    g = GR.dependencies(plan)
    I, J, V = Int[], Int[], Float64[]
    denom = zeros(Int(ncells(dst)))
    hasdenom = false
    for d in 1:nchunks(dst)
        dinds = ownedindices(dst, d)
        for s in GR.sourcesof(g, d)
            block = GR.buildblock(plan, d, Int(s))
            sinds = ownedindices(src, Int(s))
            W = block.weights
            rows, vals = rowvals(W), nonzeros(W)
            for j in axes(W, 2), t in nzrange(W, j)
                push!(I, Int(dinds[rows[t]]))
                push!(J, Int(sinds[j]))
                push!(V, vals[t])
            end
            if block.denom !== nothing
                hasdenom = true
                for (r, i) in enumerate(dinds)
                    denom[Int(i)] += block.denom[r]
                end
            end
        end
    end
    W = SparseArrays.sparse(I, J, V, Int(ncells(dst)), Int(ncells(src)))
    return W, hasdenom ? denom : nothing
end

# Bit-for-bit: the same sparsity pattern and values, the same denominators, and
# outputs equal under `isequal` (so `NaN` and `missing` match themselves).
sameweights(a, b) = a == b
samedenom(::Nothing, ::Nothing) = true
samedenom(a, b) = a !== nothing && b !== nothing && a == b
samevalues(a, b) = isequal(a, b)
# A lazy read against the eager one on the same side sums its chunk blocks in
# another order, so it agrees to roundoff, with `NaN` in the same places.
closevalues(a, b) = length(a) == length(b) &&
    all(isequal(x, y) || isapprox(x, y; rtol = 1e-12) for (x, y) in zip(a, b))

# The chunk partitions of a native space and its wrapped twin coincide, which
# is what lets their blocks be compared pair by pair.
function samechunks(a, b)
    nchunks(a) == nchunks(b) || return false
    return all(collect(ownedindices(a, c)) == collect(ownedindices(b, c))
               for c in 1:nchunks(a))
end

const PARITY_METHODS = (
    ("Conservative", Conservative(), Weighted(0.5), TreeLocator()),
    ("Conservative, Extensive", Conservative(), Extensive(), TreeLocator()),
    ("NearestCell", NearestCell(), Weighted(0.5), TreeLocator()),
    ("BarycentricPoint", BarycentricPoint(), Weighted(0.5), TreeLocator()),
    ("DirectNearest", DirectNearest(), Weighted(0.5), TreeLocator()),
)

# One pair of (native, wrapped) sides, every method, eager and chunked.
function paritysuite(name, ndst, nsrc, wdst, wsrc; analytic::Bool = false)
    methods = analytic ?
        (PARITY_METHODS..., ("Conservative, AnalyticLocator", Conservative(),
            Weighted(0.5), AnalyticLocator())) : PARITY_METHODS
    @testset "$name" begin
        @test ncells(wdst) == ncells(ndst)
        @test ncells(wsrc) == ncells(nsrc)
        @test samechunks(wdst, ndst)
        @test samechunks(wsrc, nsrc)
        @test manifold(wdst) == manifold(ndst)
        data = sourcedata(nsrc)
        for (label, method, policy, locator) in methods
            @testset "$label" begin
                # Eager: one whole-domain block, and the regrid through it.
                nb = GR.wholeblock(method, ndst, nsrc, locator)
                wb = GR.wholeblock(method, wdst, wsrc, locator)
                @test sameweights(nb.weights, wb.weights)
                @test samedenom(nb.denom, wb.denom)
                nout = regrid(data; to = ndst, from = nsrc, method, missingpolicy = policy,
                    lazy = false, locator)
                wout = regrid(data; to = wdst, from = wsrc, method, missingpolicy = policy,
                    lazy = false, locator)
                @test samevalues(values_of(nout), values_of(wout))

                # Chunked: every block the plan's relation names, and the lazy
                # read through it.
                nplan = ChunkedPlan(method, policy, ndst, nsrc; locator)
                wplan = ChunkedPlan(method, policy, wdst, wsrc; locator)
                nW, nD = assembled(nplan)
                wW, wD = assembled(wplan)
                @test sameweights(nW, wW)
                @test samedenom(nD, wD)
                @test sameweights(nW, nb.weights)
                nlazy = regrid(data; to = ndst, from = nsrc, method, missingpolicy = policy,
                    lazy = true, locator)
                wlazy = regrid(data; to = wdst, from = wsrc, method, missingpolicy = policy,
                    lazy = true, locator)
                @test samevalues(values_of(nlazy), values_of(wlazy))
                @test closevalues(values_of(nlazy), values_of(nout))
            end
        end
    end
end

@testset "GISpace parity" begin

    @testset "the space contract" begin
        toy = ToyLonLatSpace(8, 4; chunks = (4, 2))
        s = GISpace(TwinLonLatGrid(toy))
        @test s isa RegridSpace
        @test ncells(s) == 32
        @test nchunks(s) == 4
        @test all(collect(ownedindices(s, c)) == collect(ownedindices(toy, c))
                  for c in 1:4)
        @test all(GR.chunkat(s, i) == GR.chunkat(toy, i) for i in 1:32)
        @test all(cellcentroid(s, i) == cellcentroid(toy, i) for i in 1:32)
        @test all(cellat(s, cellcentroid(s, i)) == i for i in 1:32)
        @test cellat(s, toy_point(0, 0)) isa Int
        @test hascellchart(s)
        @test GR.chartaxes(s)[1] ≈ GR.chartaxes(toy)[1]
        @test GR.chartperiod(s) == (360.0, nothing)
        @test GR.chartspacing(s) == GR.chartspacing(toy)
        @test GR.chartlocalindex(s, 3, 2) == localindex(toy, 3, 2)
        @test GR.hasdualcells(s)
        @test GR.samplerstate(s) isa GR.ChartState
        @test GR.hasanalyticlocation(s)
        @test Set(GR.cellneighbors(s, 10)) ==
              Set(j for j in 1:32 if j != 10 && !isempty(intersect(
                  Set(round.(Tuple(p); digits = 9) for p in cellcorners(toy, 10)),
                  Set(round.(Tuple(p); digits = 9) for p in cellcorners(toy, j)))))
        caps = GR.chunkextents(s)
        @test length(caps) == 4
        @test all(GR.US._contains(caps[c], p)
                  for c in 1:4 for i in ownedindices(s, c) for p in GR.cellcorners(s, i))
        @test GR.chunkextent(s, 2) == caps[2]
        # The restricted tree of a chunk's own positions is the tile's cursor
        # and the whole space's is its quadtree.
        @test GR.subtree(s, ownedindices(s, 2)) isa GrS.TopDownQuadtreeCursor
        @test GR.subtree(s, 1:32) isa GrS.TopDownQuadtreeCursor
        # Polar cells drop their repeated corner exactly as the toy does.
        @test all(GI.npoint(GI.getexterior(getcell(s, i))) ==
                  GI.npoint(GI.getexterior(getcell(toy, i))) for i in 1:32)

        # A grid with no chunks is one chunk over every cell.
        one = GISpace(GrI.mockchunks(TwinLonLatGrid(toy), GrI.chunks(TwinLonLatGrid(toy), (8, 4))))
        @test nchunks(one) == 1
        @test ownedindices(one, 1) == 1:32

        # A bare GridInterface grid resolves to a `GISpace`; a space stands.
        @test GR._asspace(TwinLonLatGrid(toy), "to") isa GISpace
        @test GR._asspace(s, "to") === s
        @test_throws ArgumentError GR._asspace(1:3, "to")
        @test_throws ArgumentError GISpace(1:3)
        # DGG's own resolution still wins for its grids.
        @test GR._asspace(DGG.levelgrid(DGG.HEALPixSystem(), 1), "to") isa DGG.DGGSpace
    end

    @testset "toy space and its Structured twin" begin
        # Destination centroids avoid the source's parallels and meridians, so
        # the toy's graticule `cellat` and the twin's geodesic one agree.
        dst = ToyLonLatSpace(16, 8; chunks = (4, 4))
        src = ToyLonLatSpace(8, 4; chunks = (4, 2))
        paritysuite("toy 16×8 ← 8×4", dst, src,
            GISpace(TwinLonLatGrid(dst)), GISpace(TwinLonLatGrid(src)))
        dst = ToyLonLatSpace(8, 4; chunks = (4, 2))
        src = ToyLonLatSpace(12, 5; chunks = (3, 2))
        paritysuite("toy 8×4 ← 12×5", dst, src,
            GISpace(TwinLonLatGrid(dst)), GISpace(TwinLonLatGrid(src)))
    end

    @testset "RasterGrid and RasterGrid as a grid" begin
        src = parityraster(36, 18, 7, 5)
        dst = parityraster(24, 12, 5, 4)
        @test GrI.isgrid(src)
        @test GrI.location_is_analytic(src)
        @test length(GrI.chunks(src)) == nchunks(src)
        @test GrI.grid_extent(GrI.chunks(src)[3]) == GR.chunkextent(src, 3)
        @test GR.subtree(GISpace(src), ownedindices(src, 3)) isa GrS.TopDownQuadtreeCursor
        paritysuite("raster 24×12 ← 36×18", dst, src, GISpace(dst), GISpace(src);
            analytic = true)
        toy = ToyLonLatSpace(24, 12; chunks = (2, 2))
        paritysuite("toy 24×12 ← raster 36×18", toy, src, toy, GISpace(src))
    end

    @testset "DGG level grids and their DGGSpaces" begin
        hp = DGG.levelgrid(DGG.HEALPixSystem(), 3)
        i7 = DGG.levelgrid(DGG.IGeo7System(), 2)
        hps, i7s = GISpace(hp), GISpace(i7)
        @test GR.hasdualcells(hps) && GR.hasdualcells(i7s)
        @test GR.samplerstate(hps) isa GR.DualState
        @test !hascellchart(hps)
        # A rooted chunk's tree is the provider's own hierarchy cursor (DGG answers
        # `GridInterface.subtree` for its level grids), the same tree `DGGSpace` built.
        @test GR.subtree(hps, ownedindices(hps, 2)) isa DGG.HierarchicalGridCursor
        @test typeof(GR.subtree(hps, ownedindices(hps, 2))) ==
              typeof(GR.subtree(DGG.DGGSpace(hp; chunklevel = 0), ownedindices(hps, 2)))
        @test GR.supportradius(BarycentricPoint(), hps) ==
              GR.supportradius(BarycentricPoint(), DGG.DGGSpace(hp))
        # The dual cells are the ones the DGG space finds, point for point.
        for (grid, space) in ((hp, hps), (i7, i7s))
            native = GR.sampler(BarycentricPoint(), DGG.DGGSpace(grid))
            wrapped = GR.sampler(BarycentricPoint(), space)
            Random.seed!(11)
            for _ in 1:200
                v = randn(3)
                p = GO.UnitSphericalPoint(v ./ sqrt(sum(abs2, v))...)
                a, b = GR.dualcellat(native, p), GR.dualcellat(wrapped, p)
                @test collect(a.indices) == collect(b.indices)
                @test collect(a.nodes) == collect(b.nodes)
                @test a.kind == b.kind
            end
        end
        paritysuite("HEALPix L3 ← IGeo7 L2", DGG.DGGSpace(hp), DGG.DGGSpace(i7), hps, i7s;
            analytic = true)
        paritysuite("IGeo7 L2 ← HEALPix L3", DGG.DGGSpace(i7), DGG.DGGSpace(hp), i7s, hps;
            analytic = true)
    end

    @testset "the four relation oracles agree for GISpace sources" begin
        cases = (
            ("twin source", ToyLonLatSpace(12, 6; chunks = (3, 2)),
                GISpace(TwinLonLatGrid(ToyLonLatSpace(8, 4; chunks = (2, 1)))), 0.0),
            ("twin source, support", ToyLonLatSpace(12, 6; chunks = (3, 2)),
                GISpace(TwinLonLatGrid(ToyLonLatSpace(8, 4; chunks = (2, 1)))), 0.2),
            ("raster source", ToyLonLatSpace(24, 12; chunks = (2, 2)),
                GISpace(parityraster(36, 18, 7, 5)), 0.0),
            ("HEALPix source", ToyLonLatSpace(24, 12; chunks = (2, 2)),
                GISpace(DGG.levelgrid(DGG.HEALPixSystem(), 2)), 0.0),
        )
        for (name, dst, src, radius) in cases
            @testset "$name" begin
                truth = contributing_pairs(dst, src; radius)
                graph = graph_pairs(GR.chunk_dependency_graph(dst, src; radius))
                @test !isempty(truth)
                @test truth ⊆ graph
                @test graph == demanded_pairs(dst, src; radius)
                @test truth ⊆ capjoin_pairs(dst, src; radius)
                # A `GISpace` keeps the generic packed index, so its relation is
                # the cap join exactly, and the eager weights lie inside it.
                rows = graph_pairs(GR.dependencies(
                    ChunkedPlan(GR.Conservative(), Weighted(0.5), dst, src)))
                eager = eager_pairs(dst, src)
                @test !isempty(eager)
                @test eager ⊆ rows
                @test rows == capjoin_pairs(dst, src)
            end
        end
    end
end
