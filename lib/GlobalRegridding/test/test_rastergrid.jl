# `RasterGrid` cell geometry and orientation. Owned by task T2.

import DimensionalData as DD
import DiskArrays
import ConservativeRegridding as CR
import GeometryOps: SpatialTreeInterface as STI

# A disk-backed parent with the chunking of our choosing, which counts reads:
# a `RasterGrid` must be buildable and fully queryable without ever touching one.
mutable struct CountingChunked{T,N,A<:AbstractArray{T,N}} <: DiskArrays.AbstractDiskArray{T,N}
    data::A
    chunks::DiskArrays.GridChunks{N,NTuple{N,DiskArrays.RegularChunks}}
    reads::Int
end

CountingChunked(a::AbstractArray, cs::Tuple) =
    CountingChunked(a, DiskArrays.GridChunks(a, cs), 0)

Base.size(x::CountingChunked) = size(x.data)
DiskArrays.haschunks(::CountingChunked) = DiskArrays.Chunked()
DiskArrays.eachchunk(x::CountingChunked) = x.chunks
function DiskArrays.readblock!(x::CountingChunked, out, r::AbstractUnitRange...)
    x.reads += 1
    out .= view(x.data, r...)
    return out
end

# A global 8×4 raster of 45° cells, and the same cells addressed however the
# caller's lookups happen to be arranged.
raster_lon() = -157.5:45.0:157.5
raster_lat() = -67.5:45.0:67.5

rg_forward() = RasterGrid(DD.DimArray(zeros(8, 4),
    (DD.X(raster_lon()), DD.Y(raster_lat()))))

cellring(space, i) = collect(GI.getpoint(GI.getexterior(getcell(space, i))))

# The ring as a rotation-invariant sequence: equality of two keys means the same
# corners *and* the same winding, which is what the orientation law needs.
function ringkey(space, i)
    r = cellring(space, i)[1:end-1]
    k = argmin([(p[1], p[2], p[3]) for p in r])
    return circshift(r, 1 - k)
end

# Cells matched across two spaces by centroid, so neither space's numbering is
# assumed.
function same_cells(a, b)
    ncells(a) == ncells(b) || return false
    index = Dict(round.(Tuple(cellcentroid(a, i)), digits = 10) => i for i in 1:ncells(a))
    for j in 1:ncells(b)
        key = round.(Tuple(cellcentroid(b, j)), digits = 10)
        haskey(index, key) || return false
        ringkey(a, index[key]) == ringkey(b, j) || return false
    end
    return true
end

# Every node extent contains the corners of every cell beneath it, and every
# leaf entry's cap contains its own cell. Discovery prunes on these.
function tree_covers(space, node)
    extent = STI.node_extent(node)
    if STI.isleaf(node)
        for (i, cap) in STI.child_indices_extents(node), p in cellring(space, i)
            (GR.US._contains(cap, p) && GR.US._contains(extent, p)) || return false
        end
        return true
    end
    return all(child -> tree_covers(space, child), STI.getchild(node))
end

@testset "RasterGrid" begin

    @testset "cell geometry" begin
        space = rg_forward()
        sphere = manifold(space)
        oriented = GOCore.Spherical(; radius = 1.0, oriented = true)

        @test ncells(space) == 32
        @test hascellchart(space)

        # The cells tile the sphere: edges are shared exactly, so no gap or
        # overlap survives the sum.
        @test sum(GO.area(sphere, getcell(space, i)) for i in 1:ncells(space)) ≈ 4pi

        # Read with orientation honoured, every cell is the region it bounds
        # rather than its complement — counter-clockwise seen from outside.
        @test all(GO.area(oriented, getcell(space, i)) < 2pi for i in 1:ncells(space))

        # `cellat` inverts `cellcentroid`: the chart is not transposed and the
        # edge search lands in the cell the centroid came from.
        @test all(cellat(space, cellcentroid(space, i)) == i for i in 1:ncells(space))

        # Interval bounds and midpointed points describe the same cells.
        intervals = RasterGrid(DD.DimArray(zeros(8, 4), (
            DD.X(DD.Lookups.Sampled(-180.0:45.0:135.0;
                sampling = DD.Lookups.Intervals(DD.Lookups.Start()))),
            DD.Y(DD.Lookups.Sampled(-90.0:45.0:45.0;
                sampling = DD.Lookups.Intervals(DD.Lookups.Start()))))))
        @test intervals.xedges == space.xedges
        @test intervals.yedges == space.yedges

        # A raster need not cover the sphere, and says so.
        patch = RasterGrid(DD.DimArray(zeros(4, 2), (DD.X(5.0:10.0:35.0), DD.Y(5.0:10.0:15.0))))
        @test cellat(patch, GR.LonLatToSphere()(20.0, 10.0)) isa Int
        @test cellat(patch, GR.LonLatToSphere()(20.0, -10.0)) === nothing
        @test cellat(patch, GR.LonLatToSphere()(-20.0, 10.0)) === nothing

        # A raster numbered 0–360 still answers for a negative longitude.
        wrapped = RasterGrid(DD.DimArray(zeros(8, 4), (DD.X(22.5:45.0:337.5), DD.Y(raster_lat()))))
        @test cellat(wrapped, GR.LonLatToSphere()(-170.0, 0.0)) ==
              GR.cellposition(wrapped, 5, 3)
    end

    @testset "lookup order" begin
        # Acceptance law 4 at its root: the geometry is a property of the cells,
        # not of the order the lookups happen to be stored in. A space that read
        # its edges in lookup order without normalizing would wind these
        # backwards; one that assumed X came first would renumber them.
        forward = rg_forward()
        sphere = manifold(forward)
        oriented = GOCore.Spherical(; radius = 1.0, oriented = true)

        variants = (
            "reverse lat" => (DD.X(raster_lon()), DD.Y(reverse(raster_lat()))),
            "reverse lon" => (DD.X(reverse(raster_lon())), DD.Y(raster_lat())),
            "reverse both" => (DD.X(reverse(raster_lon())), DD.Y(reverse(raster_lat()))),
        )
        for (name, ds) in variants
            other = RasterGrid(DD.DimArray(zeros(8, 4), ds))
            @test same_cells(forward, other)
            @test all(GO.area(oriented, getcell(other, i)) < 2pi for i in 1:ncells(other))
            @test sum(GO.area(sphere, getcell(other, i)) for i in 1:ncells(other)) ≈ 4pi
            @test all(cellat(other, cellcentroid(other, i)) == i for i in 1:ncells(other))
        end

        # Y-major storage: the same cells, and the flattening follows the array.
        transposed = RasterGrid(DD.DimArray(zeros(4, 8), (DD.Y(raster_lat()), DD.X(raster_lon()))))
        @test !transposed.xfast
        @test same_cells(forward, transposed)
        @test all(cellat(transposed, cellcentroid(transposed, i)) == i
                  for i in 1:ncells(transposed))
        @test GR.cellposition(transposed, 3, 2) == 2 + (3 - 1) * 4
    end

    @testset "chunks" begin
        parent_x = CountingChunked(zeros(8, 4), (4, 2))
        space = RasterGrid(DD.DimArray(parent_x, (DD.X(raster_lon()), DD.Y(raster_lat()))))

        @test nchunks(space) == 4
        # Chunks partition the cell positions exactly: anything downstream that
        # sums per chunk double-counts or drops cells if this slips.
        covered = reduce(vcat, [collect(cellindices(space, c)) for c in 1:nchunks(space)])
        @test sort(covered) == collect(1:ncells(space))
        @test cellindices(space, 1) == [1, 2, 3, 4, 9, 10, 11, 12]

        # Geometry comes from the lookups alone — a disk-backed raster is never
        # read to build or query a space.
        celltree(space)
        chunktree(space)
        foreach(i -> getcell(space, i), 1:ncells(space))
        @test parent_x.reads == 0

        # A chunk spanning the whole fastest dimension is contiguous in position
        # space, and says so; a partial one cannot be.
        rows = RasterGrid(DD.DimArray(CountingChunked(zeros(8, 4), (8, 2)),
            (DD.X(raster_lon()), DD.Y(raster_lat()))))
        @test cellindices(rows, 2) isa AbstractUnitRange
        @test cellindices(rows, 2) == 17:32
        @test !(cellindices(space, 1) isa AbstractUnitRange)

        # The fastest dimension is the array's, not X's: here full *Y* columns
        # are the contiguous ones, and chunks run along X.
        cols = RasterGrid(DD.DimArray(CountingChunked(zeros(4, 8), (4, 2)),
            (DD.Y(raster_lat()), DD.X(raster_lon()))))
        @test nchunks(cols) == 4
        @test cellindices(cols, 2) isa AbstractUnitRange
        @test cellindices(cols, 2) == 9:16

        # An unchunked parent is one whole-domain chunk.
        @test nchunks(rg_forward()) == 1
        @test cellindices(rg_forward(), 1) == 1:32

        # Every chunk extent covers the geometry of its own cells. Discovery
        # prunes on these, so a cap that does not cover silently drops pairs.
        # Checked on a regional raster, where the caps are tight enough that the
        # law can actually fail.
        region = RasterGrid(DD.DimArray(CountingChunked(zeros(8, 4), (4, 2)),
            (DD.X(-35.0:10.0:35.0), DD.Y(-15.0:10.0:15.0))))
        caps = chunktree(region).caps
        @test all(GR.US._contains(caps[c], p)
                  for c in 1:nchunks(region)
                  for i in cellindices(region, c)
                  for p in cellring(region, i))
        @test maximum(cap.radius for cap in caps) < pi / 2
    end

    @testset "trees" begin
        space = rg_forward()
        regional = RasterGrid(DD.DimArray(zeros(8, 4),
            (DD.X(-35.0:10.0:35.0), DD.Y(-15.0:10.0:15.0))))

        # The covering law, at every level of the derived-extent tree.
        @test tree_covers(space, celltree(space))
        @test tree_covers(regional, celltree(regional))
        @test tree_covers(space, celltree(space, 1))
        @test tree_covers(space, celltree(space, [1, 9, 17, 25]))

        # `ConservativeRegridding` addresses a cell tree as a cell source.
        tree = celltree(space)
        @test Trees.ncells(tree) == ncells(space)
        @test Trees.getcell(tree, 3) == getcell(space, 3)
        @test GOCore.best_manifold(tree) == manifold(space)

        # End to end against the consumer these trees exist for: a coarse raster
        # over a fine one, intersected exactly. Every fine cell is accounted for
        # once and every coarse cell receives its own area, which fails at once
        # if the tree fails to cover, the winding is inverted, or leaf indices
        # and cell positions disagree.
        coarse = RasterGrid(DD.DimArray(zeros(4, 2), (DD.X(-135.0:90.0:135.0), DD.Y(-45.0:90.0:45.0))))
        areas = CR.intersection_areas(manifold(space), GOCore.False(),
            celltree(coarse), celltree(space); progress = false)
        @test size(areas) == (ncells(coarse), ncells(space))
        @test sum(areas) ≈ 4pi
        @test all(≈(4pi / ncells(coarse)), sum(areas; dims = 2))
        @test vec(sum(areas; dims = 1)) ≈
              [GO.area(manifold(space), getcell(space, i)) for i in 1:ncells(space)]

        # And unchanged when the source's latitude lookup runs the other way.
        flipped = RasterGrid(DD.DimArray(zeros(8, 4),
            (DD.X(raster_lon()), DD.Y(reverse(raster_lat())))))
        flipped_areas = CR.intersection_areas(manifold(space), GOCore.False(),
            celltree(coarse), celltree(flipped); progress = false)
        @test sum(flipped_areas) ≈ 4pi
        @test all(≈(4pi / ncells(coarse)), sum(flipped_areas; dims = 2))
    end
end
