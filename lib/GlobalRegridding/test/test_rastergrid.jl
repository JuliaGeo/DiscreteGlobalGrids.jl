# Raster geometry, orientation, and chunking.

import DimensionalData as DD
import DiskArrays
import ConservativeRegridding as CR
import GeometryOps: SpatialTreeInterface as STI

# Disk-backed test array with explicit chunking and read counts.
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

# Equivalent global 8×4 rasters with different lookup layouts.
raster_lon() = -157.5:45.0:157.5
raster_lat() = -67.5:45.0:67.5

rg_forward() = RasterGrid(DD.DimArray(zeros(8, 4),
    (DD.X(raster_lon()), DD.Y(raster_lat()))))

cellring(space, i) = collect(GI.getpoint(GI.getexterior(getcell(space, i))))

# Rotation-invariant ring key that preserves winding.
function ringkey(space, i)
    r = cellring(space, i)[1:end-1]
    k = argmin([(p[1], p[2], p[3]) for p in r])
    return circshift(r, 1 - k)
end

# Match cells across spaces by centroid.
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

# Direct chart constructions used as references for edge-table fast paths.

# The cell-cap fast paths sample the four box corners and nothing between them
# (`_CELL_CAP_SAMPLES` is zero), so a corner reference is the whole reference.
@noinline function reference_corners(space, ix, iy)
    xlo, xhi, ylo, yhi = GR.cellbox(space, ix, iy)
    t = space.transform
    return (t(xlo, ylo), t(xhi, ylo), t(xhi, yhi), t(xlo, yhi))
end

@noinline function reference_cellcap(space, ix, iy)
    corners = reference_corners(space, ix, iy)
    sx = sy = sz = 0.0
    for p in corners
        sx += p[1]
        sy += p[2]
        sz += p[3]
    end
    nrm = sqrt(sx^2 + sy^2 + sz^2)
    nrm <= eps(Float64) && return GR._WHOLE_SPHERE
    centre = USPoint(sx / nrm, sy / nrm, sz / nrm)
    r = 0.0
    for p in corners
        r = max(r, GR.US.spherical_distance(centre, p))
    end
    r = nextfloat(r * 1.0001 + 1e-12)
    r > Float64(pi) / 2 && return GR._WHOLE_SPHERE
    return SphericalCap(centre, r)
end

@noinline function reference_cellring(corners::NTuple{4,USPoint})
    ring = USPoint[]
    for p in corners
        (isempty(ring) || ring[end] != p) && push!(ring, p)
    end
    length(ring) > 1 && ring[end] == ring[1] && pop!(ring)
    push!(ring, ring[1])
    return ring
end

@noinline function reference_cell(space, i)
    ix, iy = GR.cellsubscript(space, i)
    c = reference_corners(space, ix, iy)
    return reference_cellring(space.ccw ? c : (c[4], c[3], c[2], c[1]))
end

# Bitwise comparison includes signed zero and NaN payloads.
bitequal(a::USPoint, b::USPoint) = all(a[k] === b[k] for k in 1:3)
bitequal(a, b) = bitequal(a.point, b.point) && a.radius === b.radius

# Wrapper that forces per-point chart evaluation.
struct OpaqueChart end
(::OpaqueChart)(x, y) = GR.LonLatToSphere()(x, y)

# Sample great-circle edges to test coverage between corners.
function ring_samples(space, i, nseg = 6)
    ring = cellring(space, i)
    out = USPoint[]
    for k in 1:(length(ring)-1)
        p, q = ring[k], ring[k+1]
        for t in range(0, 1; length = nseg + 1)
            v = (1 - t) .* Tuple(p) .+ t .* Tuple(q)
            n = sqrt(sum(v .^ 2))
            n > 0 && push!(out, USPoint(v[1] / n, v[2] / n, v[3] / n))
        end
    end
    return out
end

# Verify caps cover sampled geodesic edges.
function tree_covers_dense(space, node)
    extent = STI.node_extent(node)
    if STI.isleaf(node)
        for (i, cap) in STI.child_indices_extents(node), p in ring_samples(space, i)
            (GR.US._contains(cap, p) && GR.US._contains(extent, p)) || return false
        end
        return true
    end
    return all(child -> tree_covers_dense(space, child), STI.getchild(node))
end

    # Node and leaf extents cover their cells.
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

        # Rings are counter-clockwise from outside.
        @test all(GO.area(oriented, getcell(space, i)) < 2pi for i in 1:ncells(space))

        # Point location inverts centroids.
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
        # Lookup order and array orientation do not change cell geometry.
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
        # Chunks partition cell positions.
        covered = reduce(vcat, [collect(cellindices(space, c)) for c in 1:nchunks(space)])
        @test sort(covered) == collect(1:ncells(space))
        @test cellindices(space, 1) == [1, 2, 3, 4, 9, 10, 11, 12]

        # Geometry queries never read raster values.
        celltree(space)
        chunktree(space)
        foreach(i -> getcell(space, i), 1:ncells(space))
        @test parent_x.reads == 0

        # Full-width chunks are contiguous in position order.
        rows = RasterGrid(DD.DimArray(CountingChunked(zeros(8, 4), (8, 2)),
            (DD.X(raster_lon()), DD.Y(raster_lat()))))
        @test cellindices(rows, 2) isa AbstractUnitRange
        @test cellindices(rows, 2) == 17:32
        @test !(cellindices(space, 1) isa AbstractUnitRange)

        # Y-major arrays make full Y columns contiguous.
        cols = RasterGrid(DD.DimArray(CountingChunked(zeros(4, 8), (4, 2)),
            (DD.Y(raster_lat()), DD.X(raster_lon()))))
        @test nchunks(cols) == 4
        @test cellindices(cols, 2) isa AbstractUnitRange
        @test cellindices(cols, 2) == 9:16

        # An unchunked parent is one whole-domain chunk.
        @test nchunks(rg_forward()) == 1
        @test cellindices(rg_forward(), 1) == 1:32

        # `chunkat` inverts `cellindices` without data reads.
        for s in (space, rows, cols, rg_forward())
            @test all(GR.chunkat(s, i) == c
                      for c in 1:nchunks(s) for i in cellindices(s, c))
        end
        @test parent_x.reads == 0
        @test_throws BoundsError GR.chunkat(space, 0)
        @test_throws BoundsError GR.chunkat(space, ncells(space) + 1)

        # The point form is `cellat` composed with it, and answers nothing
        # exactly where `cellat` does.
        @test GR.chunkat(space, cellcentroid(space, 12)) == GR.chunkat(space, 12)
        patch = RasterGrid(DD.DimArray(zeros(4, 2),
            (DD.X(0.0:10.0:30.0), DD.Y(0.0:10.0:10.0))))
        @test GR.chunkat(patch, GR.LonLatToSphere()(0.0, -80.0)) === nothing

        # Regional chunk extents cover their cell geometry.
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

        # Split weights are the node's own window: children partition the parent.
        @test Trees.split_weight(tree) == ncells(space)
        @test sum(Trees.split_weight, STI.getchild(tree)) == Trees.split_weight(tree)
        @test all(c -> Trees.split_weight(c) < Trees.split_weight(tree), STI.getchild(tree))
        @test Trees.split_weight(celltree(space, [1, 9, 17, 25])) == 4

        # Tree intersections account for every fine cell exactly once.
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

    @testset "edge tables" begin
        # Edge tables match direct chart evaluation bit for bit.
        spaces = (
            "forward" => rg_forward(),
            "reverse lat" => RasterGrid(DD.DimArray(zeros(8, 4),
                (DD.X(raster_lon()), DD.Y(reverse(raster_lat()))))),
            "transposed" => RasterGrid(DD.DimArray(zeros(4, 8),
                (DD.Y(raster_lat()), DD.X(raster_lon())))),
            "0..360" => RasterGrid(DD.DimArray(zeros(8, 4),
                (DD.X(22.5:45.0:337.5), DD.Y(raster_lat())))),
            "fine" => RasterGrid(DD.DimArray(zeros(72, 36),
                (DD.X(-177.5:5.0:177.5), DD.Y(-87.5:5.0:87.5)))),
        )
        for (name, space) in spaces
            @test space.tables isa GR.LonLatEdgeTables
            @test all(1:ncells(space)) do i
                ix, iy = GR.cellsubscript(space, i)
                bitequal(GR._rastercellcap(space, ix, iy), reference_cellcap(space, ix, iy))
            end
            @test all(1:ncells(space)) do i
                r, ref = cellring(space, i), reference_cell(space, i)
                length(r) == length(ref) && all(bitequal(r[k], ref[k]) for k in eachindex(r))
            end
        end

        # Fixed-size ring construction handles repeated corners.
        pts = (USPoint(1.0, 0.0, 0.0), USPoint(0.0, 1.0, 0.0), USPoint(0.0, 0.0, 1.0))
        @test all(Iterators.product(1:3, 1:3, 1:3, 1:3)) do (a, b, c, d)
            corners = (pts[a], pts[b], pts[c], pts[d])
            GR._cellring(corners) == reference_cellring(corners)
        end
        # A cell against the pole loses one corner to the repeat; one away from
        # it keeps all four and closes with a fifth point.
        @test length(cellring(rg_forward(), 1)) == 4
        @test length(cellring(rg_forward(), 9)) == 5

        # Untabulated charts use the equivalent per-point path.
        opaque = RasterGrid(DD.DimArray(zeros(8, 4), (DD.X(raster_lon()), DD.Y(raster_lat())));
            transform = OpaqueChart(), inverse = GR.SphereToLonLat(), xperiod = 360.0,
            ybounds = (-90.0, 90.0))
        @test opaque.tables === nothing
        @test all(1:ncells(opaque)) do i
            ix, iy = GR.cellsubscript(opaque, i)
            bitequal(GR._rastercellcap(opaque, ix, iy), reference_cellcap(opaque, ix, iy))
        end
        @test all(bitequal(GR._rastercellcap(opaque, ix, iy),
                      GR._rastercellcap(rg_forward(), ix, iy)) for ix in 1:8, iy in 1:4)
    end

    @testset "wide chunk extents" begin
        # Full-longitude bands use bounded caps instead of the whole sphere.
        bands = RasterGrid(DD.DimArray(zeros(36, 18),
                (DD.X(-175.0:10.0:175.0), DD.Y(-85.0:10.0:85.0)));
            chunks = ([1:36], [3(k-1)+1:3k for k in 1:6]))
        caps = chunktree(bands).caps
        north = GR.LonLatToSphere()(0.0, 90.0)
        south = GR.LonLatToSphere()(0.0, -90.0)

        # Band caps cover bowed geodesic edges.
        @test all(GR.US._contains(caps[c], p)
                  for c in 1:nchunks(bands)
                  for i in cellindices(bands, c)
                  for p in ring_samples(bands, i))

        # Band caps exclude the opposite pole.
        @test all(cap.radius < Float64(pi) for cap in caps)
        @test !GR.US._contains(caps[1], north)
        @test !GR.US._contains(caps[3], north)
        @test !GR.US._contains(caps[4], south)
        @test !GR.US._contains(caps[6], south)
        # A band away from the equator is bounded by its own pole exactly.
        @test caps[1].radius ≈ deg2rad(30) rtol = 1e-3

        # A band straddling the equator legitimately stays large — no cap holds
        # the whole equator and less than a hemisphere — but it is still a cap.
        straddling = RasterGrid(DD.DimArray(zeros(36, 18),
                (DD.X(-175.0:10.0:175.0), DD.Y(-85.0:10.0:85.0)));
            chunks = ([1:36], [1:6, 7:12, 13:18]))
        mid = chunktree(straddling).caps[2]
        @test mid.radius < Float64(pi)
        @test all(GR.US._contains(mid, p)
                  for i in cellindices(straddling, 2) for p in ring_samples(straddling, i))

        # Pole-to-pole stripes use a mid-meridian cap.
        stripes = RasterGrid(DD.DimArray(zeros(36, 18),
                (DD.X(-175.0:10.0:175.0), DD.Y(-85.0:10.0:85.0)));
            chunks = ([6(k-1)+1:6k for k in 1:6], [1:18]))
        stripecaps = chunktree(stripes).caps
        @test all(cap.radius < Float64(pi) for cap in stripecaps)
        @test all(GR.US._contains(stripecaps[c], p)
                  for c in 1:nchunks(stripes)
                  for i in cellindices(stripes, c)
                  for p in ring_samples(stripes, i))
        # A stripe does not reach the meridian opposite it.
        @test !GR.US._contains(stripecaps[1], GR.LonLatToSphere()(0.0, 0.0))

        # Every node of the tree covers too — the extents the descent prunes on
        # are the same construction — and here against the filled-in geometry.
        @test tree_covers_dense(bands, celltree(bands))
        @test all(tree_covers_dense(bands, celltree(bands, c)) for c in 1:nchunks(bands))
    end
end
