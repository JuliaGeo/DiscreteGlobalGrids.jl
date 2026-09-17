# Acceptance: a tiled global DEM onto an area-matched DGGS level, including the
# destination chunk that holds the south pole.
#
# That chunk is the fan-in case: a cap barely a degree across whose cells are
# spanned by a whole latitude row of source tiles. It is the case that has to
# stream, the case whose weights are worth spilling, and the case where a
# missed chunk pair is invisible in the answer's shape and obvious in its
# values.
#
# The source is a disk array that synthesizes each block it is asked for and
# records the request, so every claim here is about what was read and what was
# held — not about a timing. The whole destination is 8.2 million cells and
# never materializes: one chunk of it is read, and nothing else is computed.

module RegridAcceptanceTests

using Test
import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import DimensionalData as DD
import DiskArrays
import GeometryOps as GO

# ===========================================================================
# The source: a tiled global DEM
# ===========================================================================

const NX, NY = 3600, 1800                       # 0.1°, a global DEM mosaic
const CHUNK = 300                               # 12 × 6 tiles
const CHUNKBYTES = CHUNK * CHUNK * sizeof(Float32)

# Analytic terrain, everything below sea level read as nodata. The ocean
# patches are irregular and reach into the polar band, so the coverage
# normalization is exercised where the cells are smallest.
elevation(lon, lat) = Float32(1800 * sinpi(lon / 60) * cospi(lat / 50) +
                              900 * cospi(lon / 23 + lat / 17) + 1200)
isocean(lon, lat) = elevation(lon, lat) < 0

mutable struct DEMTiles{T} <: DiskArrays.AbstractDiskArray{T,2}
    lon::Vector{Float64}
    lat::Vector{Float64}
    chunks::DiskArrays.GridChunks{2,NTuple{2,DiskArrays.RegularChunks}}
    reads::Vector{NTuple{2,UnitRange{Int}}}
    nodata::T
end

Base.size(d::DEMTiles) = (length(d.lon), length(d.lat))
DiskArrays.haschunks(::DEMTiles) = DiskArrays.Chunked()
DiskArrays.eachchunk(d::DEMTiles) = d.chunks

function DiskArrays.readblock!(d::DEMTiles{T}, out, r::AbstractUnitRange...) where {T}
    push!(d.reads, map(UnitRange{Int}, r))
    @inbounds for (jj, j) in enumerate(r[2]), (ii, i) in enumerate(r[1])
        x, y = d.lon[i], d.lat[j]
        out[ii, jj] = isocean(x, y) ? d.nodata : elevation(x, y)
    end
    return out
end

_axis(D, centres, step) = D(DD.Sampled(centres; span = DD.Regular(step),
    sampling = DD.Intervals(DD.Center()),
    order = step > 0 ? DD.ForwardOrdered() : DD.ReverseOrdered()))

function demtiles(nodata)
    dx, dy = 360 / NX, 180 / NY
    lon = collect((-180 + dx / 2):dx:180)
    lat = collect((90 - dy / 2):-dy:-90)
    tiles = DEMTiles(lon, lat, DiskArrays.GridChunks((NX, NY), (CHUNK, CHUNK)),
        NTuple{2,UnitRange{Int}}[], nodata)
    return DD.DimArray(tiles, (_axis(DD.X, lon, dx), _axis(DD.Y, lat, -dy)))
end

reads(raster) = parent(raster).reads
resetreads!(raster) = (empty!(reads(raster)); raster)

# ===========================================================================
# The destination: the chunk that holds the south pole
# ===========================================================================

const SYS = DGG.IGeo7System()
const DEM = demtiles(NaN32)

# Neither of these reads a source block: the level comes from cell areas and the
# chunking from the hierarchy.
const PLAN = DGG.plan_regrid(DEM; to = SYS, lazy = true)
const DST = PLAN.dst_space
const POLE = GO.UnitSphericalPoint(0.0, 0.0, -1.0)
const POLARCHUNK = GR.chunkat(DST, POLE)
const POLARCELLS = GR.ownedindices(DST, POLARCHUNK)

# The tile row those cells lie in — the fan-in the chunk pairs have to find.
const POLARTILE = (NY - CHUNK + 1):NY

# One destination chunk, through the bare-system target the page-level API uses.
function polarchunk(raster; kwargs...)
    A = parent(DGG.regrid(raster; to = SYS, lazy = true, kwargs...))
    return A[POLARCELLS], GR.residency(A)
end

# Blanked destinations are NaN on both sides; the rest agree to fp noise. A
# missed chunk pair moves a value by a whole source tile's worth of terrain, so
# the tolerance is not what is under test.
agree(a, b) = length(a) == length(b) &&
    all(isnan(x) ? isnan(y) : isapprox(x, y; rtol = 1e-5) for (x, y) in zip(a, b))

# The same cells computed eagerly off a source cropped to the latitude band they
# can reach — whole source rows, so the retained cells' weights are unchanged
# and the two answers are comparable exactly.
function eagerreference()
    g = DST.grid
    tolonlat = GO.UnitSpherical.GeographicFromUnitSphere()
    top = maximum(tolonlat(p)[2] for i in POLARCELLS
                  for p in DGG.cell_boundary(g, DGG.cellindex(g, i)))
    rows = findall(<(top + 1.0), collect(DD.lookup(DEM, DD.Y)))
    sub = DGG.PartialGrid(DGG.CellVector(g)[POLARCELLS])
    return parent(DGG.regrid(DEM[:, rows]; to = sub, lazy = false))
end

@testset "a tiled DEM onto an area-matched DGGS" begin

    @testset "planning reads nothing" begin
        # Weights are geometry-only and a chunked plan builds them on first
        # touch, so planning a regrid of a grid this size costs no IO at all.
        @test isempty(reads(DEM))
        @test PLAN isa GR.ChunkedPlan
        # A 0.1° source area-matches IGeo7 level 7.
        @test DGG.level(DST.grid) == 7
        A = parent(DGG.regrid(DEM; to = SYS, lazy = true))
        @test size(A) == (DGG.ncells(DST),)
        @test isempty(reads(DEM))
        @test GR.residency(A).loads == 0

        # The chunk under test is named by `chunkat`, which is one binary search
        # over the chunk windows. Against the scan it replaces: a wrong window
        # would move the whole file onto a chunk that is not the pole's, where
        # every read assertion below would still pass.
        @test POLARCHUNK == findfirst(c -> GR.cellat(DST, POLE) in GR.ownedindices(DST, c),
            1:GR.nchunks(DST))
        @test GR.cellat(DST, POLE) in POLARCELLS
    end

    held, heldstats = (resetreads!(DEM); polarchunk(DEM; budget = 2^30))
    heldreads = copy(reads(DEM))
    streamed, streamstats = (resetreads!(DEM); polarchunk(DEM; budget = 2^20))

    @testset "the south-pole chunk streams" begin
        # Every read is one whole source tile of the southernmost tile row, and
        # all twelve of them are read: over-connection would pull in tiles from
        # the row above, and a missed dilation or a cap that does not cover the
        # pole would drop one.
        @test all(r[2] == POLARTILE for r in heldreads)
        @test length(unique(first, heldreads)) == NX ÷ CHUNK
        @test heldstats.loads == NX ÷ CHUNK

        # The budget changes what is resident and nothing else.
        @test isequal(streamed, held)
        @test streamstats.loads == heldstats.loads
        @test streamstats.peakbytes <= 2 * CHUNKBYTES
        @test heldstats.peakbytes >= 8 * CHUNKBYTES

        # The blanking is real: the polar band has nodata in it, and the cells
        # that survive it are the majority.
        @test 0 < count(isnan, held) < length(held) ÷ 2
    end

    @testset "the streamed chunk is the eager answer" begin
        @test agree(held, eagerreference())
    end

    @testset "weights spill to disk" begin
        # A one-byte memory bound evicts every block as the next is built, so
        # each of the destination chunk's pairs — one per source tile across the
        # southernmost row — leaves a file behind for the second read to reload.
        dir = mktempdir()
        storage = GR.Spilled(dir; maxbytes = 1)
        A = parent(DGG.regrid(DEM; to = SYS, lazy = true, storage, budget = 2^20))
        @test isequal(A[POLARCELLS], held)
        @test length(GR.spilledfiles(storage)) == NX ÷ CHUNK
        @test isequal(A[POLARCELLS], held)
    end

    @testset "a sentinel is nodata" begin
        sentinel = demtiles(-32767.0f0)
        declared, _ = polarchunk(sentinel; budget = 2^20, missingval = -32767.0f0)
        # The sentinel is nodata on both sides: the same cells survive, and the
        # blanked ones come back holding it rather than NaN.
        unsentinel(x) = isequal(x, -32767.0f0) ? oftype(x, NaN) : x
        @test isequal(map(unsentinel, declared), held)
        @test count(isequal(-32767.0f0), declared) == count(isnan, held)
        # And it is `missingval` that says so: undeclared, the sentinel is
        # terrain, and the cells over the ocean patches come back wrong.
        undeclared, _ = polarchunk(sentinel; budget = 2^20)
        @test !isequal(undeclared, held)
    end
end

end # module RegridAcceptanceTests
