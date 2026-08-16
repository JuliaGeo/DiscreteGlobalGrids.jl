# `Conservative()` weight construction. Owned by task T3.

import ConservativeRegridding as CR

"""
    DensifiedCellSpace(lon, lat, n)

A one-cell space whose single cell is a graticule box with its parallel edges
sampled at `n` segments — a cell whose ring is **not convex**, since a
great-circle chord between two points of a parallel bows poleward of it.

Everything a conservative weight build touches is here: `ncells`, `getcell`,
`manifold`, and a `celltree`. It is the shape a densified DGGS cell has, and it
is the shape the spherical clip's convexity precondition fails on.
"""
struct DensifiedCellSpace <: RegridSpace
    lon::Tuple{Float64,Float64}
    lat::Tuple{Float64,Float64}
    n::Int
end

ncells(::DensifiedCellSpace) = 1
nchunks(::DensifiedCellSpace) = 1
cellindices(::DensifiedCellSpace, ::Int) = 1:1
manifold(::DensifiedCellSpace) = GOCore.Spherical(; radius = 1.0)
celltree(space::DensifiedCellSpace) = GR.CellCapTree(space, 1:1)

function getcell(space::DensifiedCellSpace, i::Int)
    i == 1 || throw(BoundsError(space, i))
    lon0, lon1 = space.lon
    lat0, lat1 = space.lat
    ring = USPoint[]
    for k in 0:space.n          # south edge, west to east
        push!(ring, toy_point(lon0 + (lon1 - lon0) * k / space.n, lat0))
    end
    for k in 0:space.n          # north edge, east to west
        push!(ring, toy_point(lon1 - (lon1 - lon0) * k / space.n, lat1))
    end
    push!(ring, ring[1])
    return GI.Polygon([GI.LinearRing(ring)])
end

# --- helpers ---------------------------------------------------------------

conservative_block(dst, dst_inds, src, src_inds) =
    WeightBlock(
        build_weights!(WeightCOO(length(dst_inds)), Conservative(),
            dst, dst_inds, src, src_inds),
        length(dst_inds), length(src_inds))

cellareas(space, inds) = [GO.area(manifold(space), getcell(space, i)) for i in inds]

@testset "Conservative weights" begin
    @test Conservative() isa AbstractRegriddingMethod

    @testset "identity" begin
        space = ToyLonLatSpace(4, 2)
        inds = cellindices(space, 1)
        block = conservative_block(space, inds, space, inds)
        areas = cellareas(space, inds)
        W = Matrix(block.weights)

        # A space regridded onto itself pairs each cell with itself and nothing
        # else: the diagonal is the cell's area and the row holds nothing more,
        # so neighbours sharing an edge contribute no sliver.
        @test size(W) == (8, 8)
        @test LinearAlgebra.diag(W) ≈ areas rtol = 1e-8
        @test vec(sum(W; dims = 2)) ≈ areas rtol = 1e-8

        # The denominator is accumulated coverage, which here is full coverage:
        # a builder that reported cell areas from the wrong side, or reported
        # the diagonal only once per block, fails here and nowhere else.
        @test block.denom ≈ areas rtol = 1e-8
        @test sum(W) ≈ 4pi rtol = 1e-8
    end

    @testset "refinement" begin
        coarse = ToyLonLatSpace(4, 2)
        fine = ToyLonLatSpace(8, 4)
        dst_inds = cellindices(fine, 1)
        src_inds = cellindices(coarse, 1)
        block = conservative_block(fine, dst_inds, coarse, src_inds)
        W = Matrix(block.weights)
        areas = cellareas(fine, dst_inds)

        # Every fine cell nests in exactly one coarse cell, so every
        # destination row is a single entry carrying the whole fine cell.
        @test size(W) == (32, 8)
        @test all(count(>(1e-12), view(W, j, :)) == 1 for j in 1:32)
        @test vec(sum(W; dims = 2)) ≈ areas rtol = 1e-8

        # Nothing is dropped and nothing is counted twice.
        @test sum(W) ≈ 4pi rtol = 1e-8
        @test block.denom ≈ areas rtol = 1e-8
    end

    @testset "chunked source partitions the block" begin
        fine = ToyLonLatSpace(8, 4)
        whole = ToyLonLatSpace(4, 2)
        chunked = ToyLonLatSpace(4, 2; chunks = (2, 2))
        dst_inds = cellindices(fine, 1)

        reference = conservative_block(fine, dst_inds, whole, cellindices(whole, 1))

        # Two source chunks, neither contiguous in position space — so a
        # builder that emitted cell positions instead of positions within
        # `src_inds` writes outside the block or into the wrong column.
        @test nchunks(chunked) == 2
        parts = [cellindices(chunked, c) for c in 1:nchunks(chunked)]
        @test !(parts[1] isa AbstractUnitRange)
        @test sort(vcat(parts...)) == collect(1:8)

        assembled = zeros(32, 8)
        denom = zeros(32)
        for inds in parts
            block = conservative_block(fine, dst_inds, chunked, inds)
            @test size(block) == (32, length(inds))
            assembled[:, inds] .+= Matrix(block.weights)
            denom .+= block.denom
        end

        @test assembled ≈ Matrix(reference.weights) rtol = 1e-8
        @test denom ≈ reference.denom rtol = 1e-8
    end

    @testset "cells across the antimeridian" begin
        fine = ToyLonLatSpace(8, 4)
        # A band whose cells straddle 180°: its second column is the same
        # geometry as the destination's first. Nothing here wraps longitude —
        # the cells are unit-sphere points — and this is what says so.
        band = ToyLonLatSpace(2, 2; lon = (135.0, 225.0), lat = (-45.0, 45.0))
        dst_inds = cellindices(fine, 1)
        src_inds = cellindices(band, 1)
        block = conservative_block(fine, dst_inds, band, src_inds)
        W = Matrix(block.weights)
        areas = cellareas(band, src_inds)

        # The band's cells coincide with four destination cells exactly.
        @test count(>(1e-12), W) == 4
        for k in 1:4
            j = cellat(fine, cellcentroid(band, k))
            @test W[j, k] ≈ areas[k] rtol = 1e-8
        end
        @test sum(W) ≈ sum(areas) rtol = 1e-8
    end

    @testset "non-convex cells" begin
        dense = DensifiedCellSpace((0.0, 90.0), (0.0, 45.0), 12)
        tiling = ToyLonLatSpace(16, 8)
        exact = GO.area(manifold(dense), getcell(dense, 1))
        tiling_inds = cellindices(tiling, 1)

        # A cell's intersections with a full tiling of the sphere sum to its
        # own area. As the *source* the non-convex cell is the clipped subject
        # and this holds.
        as_source = conservative_block(tiling, tiling_inds, dense, 1:1)
        @test sum(as_source.weights) ≈ exact rtol = 1e-6

        # As the *destination* it is the clip ring, and
        # `ConvexConvexSutherlandHodgman` requires a convex one. Broken until
        # the GeometryOps Sutherland-Hodgman non-convex-clip fix lands
        # (`GeometryOps/src/methods/clipping/sutherland_hodgman.jl`); this is
        # the same defect that makes conservative regridding *onto* densified
        # DGGS cells wrong.
        as_destination = conservative_block(dense, 1:1, tiling, tiling_inds)
        @test_broken isapprox(sum(as_destination.weights), exact; rtol = 1e-6)
    end

    @testset "the block is the assembled matrix, entry for entry" begin
        # `build_weights!` reads `ConservativeRegridding`'s assembled block
        # straight into the `WeightCOO`, in the order the entries are stored,
        # rather than through a `findnz` triple. That is a memory economy only,
        # and this is what says so: the block a whole-domain build produces is
        # the matrix the descent assembled, with no tolerance.
        coarse = ToyLonLatSpace(4, 2)
        fine = ToyLonLatSpace(8, 4)
        dst_inds, src_inds = cellindices(coarse, 1), cellindices(fine, 1)
        block = conservative_block(coarse, dst_inds, fine, src_inds)
        areas = CR.intersection_areas(manifold(coarse), GOCore.False(),
            GR.subtree(coarse, dst_inds), GR.subtree(fine, src_inds); progress = false)
        @test block.weights == areas
        # And the denominator is the same sum over the same entries.
        @test block.denom == vec(sum(areas; dims = 2))
    end

    @testset "banded chunk extents" begin
        # A toy space chunked into full-longitude rows — the shape a global grid
        # is normally stored in. A cap through the chunk's cell corners loses
        # convexity past π/2 and degenerates to the whole sphere, at which point
        # every band appears to touch every other; the band cap does not.
        banded = ToyLonLatSpace(16, 8; chunks = (16, 2))
        caps = chunktree(banded).caps
        @test all(cap.radius < Float64(pi) for cap in caps)

        # Covering, against the arcs and not only the corners.
        function samples(space, i, nseg = 6)
            ring = collect(GI.getpoint(GI.getexterior(getcell(space, i))))
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
        @test all(US._contains(caps[c], p)
                  for c in 1:nchunks(banded)
                  for i in cellindices(banded, c)
                  for p in samples(banded, i))

        # The southern band does not reach the north pole, which is exactly what
        # a whole-sphere extent did.
        @test !US._contains(caps[1], toy_point(0.0, 90.0))
        @test !US._contains(caps[4], toy_point(0.0, -90.0))

        # A regional space's chunks never degenerate, and keep the cap they had.
        region = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (4, 2))
        @test all(c -> c.radius < Float64(pi) / 2, chunktree(region).caps)
    end

    @testset "disjoint chunks and mismatched manifolds" begin
        north = ToyLonLatSpace(2, 1; lat = (60.0, 90.0))
        south = ToyLonLatSpace(2, 1; lat = (-90.0, -60.0))
        block = conservative_block(north, 1:2, south, 1:2)

        # Chunks that do not meet still produce a block, and it still carries a
        # denominator — zero coverage is an answer, not "no denominator, so
        # finalize as a raw sum".
        @test GR.SparseArrays.nnz(block.weights) == 0
        @test block.denom == zeros(2)
    end
end
