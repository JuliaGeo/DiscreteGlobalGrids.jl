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

# --- the pre-table constructions, kept as the reference the fast paths answer to
#
# `RasterGrid` synthesizes corners, cell caps and rings from per-edge `sin`/`cos`
# tables rather than by evaluating the chart. These are the constructions it
# replaced, transcribed; the "edge tables" testset holds the space to them bit
# for bit, which is the whole claim. `@noinline` so the comparison cannot be
# constant-folded into agreement.

@noinline function reference_boxpoint(t, xlo, xhi, ylo, yhi, m::Int, j::Int)
    e, k = divrem(j, m)
    u = k / m
    return e == 0 ? t(xlo + u * (xhi - xlo), ylo) :
           e == 1 ? t(xhi, ylo + u * (yhi - ylo)) :
           e == 2 ? t(xhi - u * (xhi - xlo), yhi) :
           t(xlo, yhi - u * (yhi - ylo))
end

@noinline function reference_cellcap(space, ix, iy)
    xlo, xhi, ylo, yhi = GR.cellbox(space, ix, iy)
    t = space.transform
    sx = sy = sz = 0.0
    for j in 0:3
        p = reference_boxpoint(t, xlo, xhi, ylo, yhi, 1, j)
        sx += p[1]
        sy += p[2]
        sz += p[3]
    end
    nrm = sqrt(sx^2 + sy^2 + sz^2)
    nrm <= eps(Float64) && return GR._WHOLE_SPHERE
    centre = USPoint(sx / nrm, sy / nrm, sz / nrm)
    r = 0.0
    for j in 0:3
        r = max(r, GR.US.spherical_distance(centre,
            reference_boxpoint(t, xlo, xhi, ylo, yhi, 1, j)))
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
    xlo, xhi, ylo, yhi = GR.cellbox(space, ix, iy)
    t = space.transform
    c = (t(xlo, ylo), t(xhi, ylo), t(xhi, yhi), t(xlo, yhi))
    return reference_cellring(space.ccw ? c : (c[4], c[3], c[2], c[1]))
end

# Equality down to the bit, `-0.0` and `NaN` included — `==` on caps and points
# would let a rounding difference through.
bitequal(a::USPoint, b::USPoint) = all(a[k] === b[k] for k in 1:3)
bitequal(a, b) = bitequal(a.point, b.point) && a.radius === b.radius

# A chart that is `LonLatToSphere` in every respect except that nothing knows it,
# so the space falls back to evaluating the transform per point.
struct OpaqueChart end
(::OpaqueChart)(x, y) = GR.LonLatToSphere()(x, y)

# A cell ring with its great-circle edges filled in — the geometry an extent has
# to cover, as opposed to the corners, which bow outside the graticule box.
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

# The same law against the filled-in geometry rather than the corners: a cap that
# holds a cell's corners but not the arcs between them does not cover.
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

        # `chunkat` inverts `cellindices` — every cell of every chunk is placed
        # back in the chunk it came from, on all three chunkings, so a lattice
        # arithmetic slip cannot hide behind a symmetric case. It answers
        # without building a `cellindices` vector, so the disk-backed parent is
        # still untouched.
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

    @testset "edge tables" begin
        # Corners, cell caps and rings come from `cos`/`sin` tabulated once per
        # cell edge rather than from evaluating the chart per point. The table
        # holds the chart's own `Float64`s, so this is not an approximation and
        # the assertion is equality of bits, not of value: a plan built through
        # the tables is the plan built through the chart, entry for entry.
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

        # Polar cells repeat a corner and cells elsewhere do not, so the ring
        # length varies; every way four corners can repeat gives the ring the
        # growing construction gave.
        pts = (USPoint(1.0, 0.0, 0.0), USPoint(0.0, 1.0, 0.0), USPoint(0.0, 0.0, 1.0))
        @test all(Iterators.product(1:3, 1:3, 1:3, 1:3)) do (a, b, c, d)
            corners = (pts[a], pts[b], pts[c], pts[d])
            GR._cellring(corners) == reference_cellring(corners)
        end
        # A cell against the pole loses one corner to the repeat; one away from
        # it keeps all four and closes with a fifth point.
        @test length(cellring(rg_forward(), 1)) == 4
        @test length(cellring(rg_forward(), 9)) == 5

        # A chart nothing can tabulate keeps the per-point path, and answers the
        # same. This is the seam a projected raster arrives through.
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
        # The shape a global raster is normally chunked in: full-longitude bands.
        # A cap through such a box's own boundary has to reach around the sphere,
        # so the sampled construction gives up and reports the whole sphere —
        # which makes every band appear to touch every source chunk. The band cap
        # is what stops that, and it has to cover.
        bands = RasterGrid(DD.DimArray(zeros(36, 18),
                (DD.X(-175.0:10.0:175.0), DD.Y(-85.0:10.0:85.0)));
            chunks = ([1:36], [3(k-1)+1:3k for k in 1:6]))
        caps = chunktree(bands).caps
        north = GR.LonLatToSphere()(0.0, 90.0)
        south = GR.LonLatToSphere()(0.0, -90.0)

        # (a) covering, against the cells' *geometry* and not just their corners:
        # an east-west cell edge is a great-circle arc that bows poleward of the
        # parallel it joins, and the band cap's latitude bound accounts for it.
        @test all(GR.US._contains(caps[c], p)
                  for c in 1:nchunks(bands)
                  for i in cellindices(bands, c)
                  for p in ring_samples(bands, i))

        # (b) meaningfully smaller than the whole sphere: no band's cap reaches
        # the pole it faces away from, which a whole-sphere cap does.
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

        # The other shape that defeats a boundary-sampled cap: a chunk running
        # pole to pole. No cap about a pole helps there, and the one on the
        # chunk's own mid-meridian is what keeps it off the whole sphere.
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
