# `NearestCell()` and `BilinearPoint()` weight construction. Owned by task T4.
#
# `toyspaces.jl` is frozen, so the chart accessor the bilinear stencil is
# written against is implemented for `ToyLonLatSpace` here. These five methods
# are the entire contract a charted space owes `BilinearPoint`; a real raster
# space implements exactly this list.

# Cell centres, in the space's native degrees. `ix` runs west to east, `iy`
# south to north, matching `cellposition`.
GR.chartaxes(space::ToyLonLatSpace) =
    ([space.lon0 + (ix - 0.5) * dlon(space) for ix in 1:space.nlon],
        [space.lat0 + (iy - 0.5) * dlat(space) for iy in 1:space.nlat])

# Onto the branch the axes are written on. A periodic axis is reduced again by
# the locator, so any representative would do there; a partial span needs this.
function GR.chartcoords(space::ToyLonLatSpace, p)
    lon, lat = toy_lonlat(p)
    return (space.lon0 + mod(lon - space.lon0, 360.0), lat)
end

GR.chartposition(space::ToyLonLatSpace, ix::Int, iy::Int) = cellposition(space, ix, iy)

GR.chartperiod(space::ToyLonLatSpace) =
    (space.lon1 - space.lon0 >= 360 ? 360.0 : nothing, nothing)

# Along a parallel the centres are Δλ·cos φ apart, so Δλ bounds them.
GR.chartspacing(space::ToyLonLatSpace) = (deg2rad(dlon(space)), deg2rad(dlat(space)))

# --- helpers ---------------------------------------------------------------

function t4_build(method, dst, dst_inds, src, src_inds)
    coo = WeightCOO(length(dst_inds))
    build_weights!(coo, method, dst, dst_inds, src, src_inds)
    return coo
end

"""
    t4_entries(coo, dst_inds, src_inds) -> Dict{(dst position, src position), weight}

A block's chunk-local entries lifted back to cell positions, which is the only
form in which two blocks over different source chunks are comparable.
"""
function t4_entries(coo, dst_inds, src_inds)
    entries = Dict{Tuple{Int,Int},Float64}()
    for k in eachindex(coo.vals)
        key = (Int(dst_inds[coo.rows[k]]), Int(src_inds[coo.cols[k]]))
        entries[key] = get(entries, key, 0.0) + coo.vals[k]
    end
    return entries
end

t4_rowsum(entries, dst_position) =
    sum(w for ((j, _), w) in entries if j == dst_position; init = 0.0)

# The coarse cell a twice-refined fine cell sits inside.
function t4_parent(coarse, fine, i)
    ix, iy = cellsubscript(fine, i)
    return cellposition(coarse, cld(ix, 2), cld(iy, 2))
end

struct T4NoChartSpace <: RegridSpace end

struct T4BareChartSpace <: RegridSpace end
GR.hascellchart(::T4BareChartSpace) = true

@testset "Interpolation weights" begin

    @testset "NearestCell" begin
        # A space onto itself is the identity: each centroid lands in its own
        # cell, weight 1, and nowhere else.
        space = ToyLonLatSpace(8, 4)
        inds = cellindices(space, 1)
        block = WeightBlock(t4_build(NearestCell(), space, inds, space, inds),
            length(inds), length(inds))
        @test Matrix(block.weights) == Matrix(LinearAlgebra.I, 32, 32)

        # No denominator: a point sample is a value, not a coverage, and the
        # executor must finalize this block as the raw weighted value.
        @test block.denom === nothing

        # Coarse source, twice-refined destination: every fine cell takes the
        # coarse cell its centroid sits in, and takes nothing else.
        coarse = ToyLonLatSpace(4, 2)
        fine = ToyLonLatSpace(8, 4)
        dst_inds, src_inds = cellindices(fine, 1), cellindices(coarse, 1)
        entries = t4_entries(t4_build(NearestCell(), fine, dst_inds, coarse, src_inds),
            dst_inds, src_inds)
        @test length(entries) == ncells(fine)
        @test all(entries[(i, t4_parent(coarse, fine, i))] == 1.0
                  for i in 1:ncells(fine))

        # A centroid outside the source's coverage emits nothing at all — the
        # missing policy decides, not the weights.
        north = ToyLonLatSpace(4, 1; lat = (0.0, 90.0))
        global_dst = ToyLonLatSpace(4, 2)
        dst_inds = cellindices(global_dst, 1)
        entries = t4_entries(
            t4_build(NearestCell(), global_dst, dst_inds, north, cellindices(north, 1)),
            dst_inds, cellindices(north, 1))
        @test sort(unique(first.(keys(entries)))) == collect(5:8)
    end

    @testset "BilinearPoint stencils" begin
        # Source centres: lon -135, -45, 45, 135; lat -45, 45.
        src = ToyLonLatSpace(4, 2)
        src_inds = cellindices(src, 1)

        # A destination centroid a quarter of the way east and three quarters
        # north between four centres. Hand values, asymmetric in both axes, so a
        # transposed lattice or a flipped fraction fails here.
        dst = ToyLonLatSpace(1, 1; lon = (-27.5, -17.5), lat = (17.5, 27.5))
        entries = t4_entries(t4_build(BilinearPoint(), dst, [1], src, src_inds),
            [1], src_inds)
        @test length(entries) == 4
        @test entries[(1, cellposition(src, 2, 1))] ≈ 0.1875
        @test entries[(1, cellposition(src, 3, 1))] ≈ 0.0625
        @test entries[(1, cellposition(src, 2, 2))] ≈ 0.5625
        @test entries[(1, cellposition(src, 3, 2))] ≈ 0.1875
        @test t4_rowsum(entries, 1) ≈ 1.0

        # The seam: a centroid at lon 180 lies between the last centre and the
        # first, not off the edge of the lattice. Without the period this
        # collapses onto column 4 alone.
        seam = ToyLonLatSpace(1, 1; lon = (175.0, 185.0), lat = (-5.0, 5.0))
        entries = t4_entries(t4_build(BilinearPoint(), seam, [1], src, src_inds),
            [1], src_inds)
        @test length(entries) == 4
        @test all(entries[(1, cellposition(src, ix, iy))] ≈ 0.25
                  for ix in (1, 4), iy in (1, 2))

        # Poleward of the outermost latitude centre the stencil degrades to
        # nearest in latitude and stays linear in longitude — it never
        # extrapolates, and the weights still sum to 1.
        polar = ToyLonLatSpace(1, 1; lon = (-5.0, 5.0), lat = (75.0, 85.0))
        entries = t4_entries(t4_build(BilinearPoint(), polar, [1], src, src_inds),
            [1], src_inds)
        @test length(entries) == 2
        @test entries[(1, cellposition(src, 2, 2))] ≈ 0.5
        @test entries[(1, cellposition(src, 3, 2))] ≈ 0.5

        # Outside both axes of a non-periodic chart: one point, weight 1.
        patch = ToyLonLatSpace(4, 2; lon = (-40.0, 40.0), lat = (-20.0, 20.0))
        corner = ToyLonLatSpace(1, 1; lon = (38.0, 40.0), lat = (18.0, 20.0))
        entries = t4_entries(
            t4_build(BilinearPoint(), corner, [1], patch, cellindices(patch, 1)),
            [1], cellindices(patch, 1))
        @test entries == Dict((1, cellposition(patch, 4, 2)) => 1.0)

        # A field linear in lon and lat is reproduced exactly at every
        # destination centroid strictly inside the source's centre lattice.
        fsrc = ToyLonLatSpace(36, 18)
        fdst = ToyLonLatSpace(17, 8; lon = (-170.0, 170.0), lat = (-80.0, 80.0))
        fdst_inds, fsrc_inds = cellindices(fdst, 1), cellindices(fsrc, 1)
        function linear(p)
            lon, lat = toy_lonlat(p)
            return 2.0 + 0.01 * lon + 0.03 * lat
        end
        field = [linear(cellcentroid(fsrc, i)) for i in fsrc_inds]
        weights = WeightBlock(
            t4_build(BilinearPoint(), fdst, fdst_inds, fsrc, fsrc_inds),
            length(fdst_inds), length(fsrc_inds)).weights
        @test weights * field ≈ [linear(cellcentroid(fdst, i)) for i in fdst_inds]
        @test all(≈(1.0), sum(weights; dims = 2))

        # The gate, and the chart contract behind it.
        @test_throws "hascellchart" build_weights!(
            WeightCOO(1), BilinearPoint(), src, [1], T4NoChartSpace(), [1])
        @test_throws "chartaxes" build_weights!(
            WeightCOO(1), BilinearPoint(), src, [1], T4BareChartSpace(), [1])
    end

    @testset "stencils partition across source chunks" begin
        # Source split east/west between longitude centres -45 and 45; every
        # destination centroid sits between them, so every stencil straddles the
        # split.
        src = ToyLonLatSpace(4, 2; chunks = (2, 2))
        dst = ToyLonLatSpace(3, 2; lon = (-30.0, 30.0), lat = (-30.0, 30.0))
        dst_inds = cellindices(dst, 1)
        @test nchunks(src) == 2

        whole_inds = cellindices(ToyLonLatSpace(4, 2), 1)
        whole = t4_entries(t4_build(BilinearPoint(), dst, dst_inds, src, whole_inds),
            dst_inds, whole_inds)

        blocks = [t4_entries(
                      t4_build(BilinearPoint(), dst, dst_inds, src, cellindices(src, c)),
                      dst_inds, cellindices(src, c)) for c in 1:nchunks(src)]

        # Each block emits only its own chunk's source cells, and both blocks
        # have something to say — otherwise this proves nothing.
        @test all(!isempty(b) for b in blocks)
        @test all(key[2] in cellindices(src, c)
                  for (c, b) in enumerate(blocks) for key in keys(b))

        # The union of the blocks is exactly the unchunked stencil: nothing
        # dropped at the boundary, nothing counted twice.
        merged = merge(+, blocks...)
        @test keys(merged) == keys(whole)
        @test all(merged[key] ≈ whole[key] for key in keys(whole))

        # And the invariant that matters to the answer: a straddling stencil
        # still sums to 1 once the executor has accumulated both source chunks.
        @test all(t4_rowsum(merged, i) ≈ 1.0 for i in dst_inds)
    end

    @testset "support radius" begin
        # Nearest-cell weights vanish outside the geometric overlap; a bilinear
        # stencil reaches into the neighbouring chunk, and says so.
        space = ToyLonLatSpace(8, 4)
        @test support_radius(NearestCell(), space) == 0.0
        @test support_radius(BilinearPoint(), space) > 0
        @test support_radius(BilinearPoint(), space) ≈ deg2rad(45.0)

        # The larger axis bounds both, so no stencil point can be missed.
        oblong = ToyLonLatSpace(36, 6)
        @test support_radius(BilinearPoint(), oblong) >= deg2rad(dlon(oblong))
        @test support_radius(BilinearPoint(), oblong) >= deg2rad(dlat(oblong))
    end
end
