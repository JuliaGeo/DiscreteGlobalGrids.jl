# A fabricated DEM with an analytic answer: the field, the land mask that gives
# it coastlines, one tile's worth of it, and the per-cell oracle predicates.

"""
    synthetic_elevation(lon, lat) -> Float64
    synthetic_elevation(p::UnitSphericalPoint) -> Float64

Stand-in elevation, degrees in and metres out:
`1000 sin(3λ) cos(2φ) + 500 cos(7λ) sin(5φ) + 100`.

The field varies on a scale four orders of magnitude larger than an IGeo7
level-12 cell, so an area-weighted cell mean equals the field at the cell centre
to well within a metre. That makes it an oracle, including over a part-ocean
cell: a renormalised mean of a near-constant field is still that constant.
"""
synthetic_elevation(lon, lat) = let λ = deg2rad(lon), φ = deg2rad(lat)
    1000 * sin(3λ) * cos(2φ) + 500 * cos(7λ) * sin(5φ) + 100
end

synthetic_elevation(p::GO.UnitSphericalPoint) =
    synthetic_elevation(US.GeographicFromUnitSphere()(p)...)

"""
    LandMask(ncols, nrows, dlon, dlat, bits)

A global land/ocean bitmap on the plate-carrée lattice: `bits[col, row]` is true
over land, column 1 starting at 180°W and row 1 at 90°N.

The mask decides which synthetic posts are nodata, which puts partly-fed cells
at every coastline. Real tiles carry their own nodata and never consult it.
"""
struct LandMask
    ncols::Int
    nrows::Int
    dlon::Float64
    dlat::Float64
    bits::BitMatrix
end

"No mask: every post of a listed tile is valid."
const NOMASK = nothing

"Is (`lon`, `lat`) over land? True everywhere under [`NOMASK`](@ref)."
@inline function island(m::LandMask, lon::Real, lat::Real)
    c = clamp(floor(Int, (lon + 180) / m.dlon) + 1, 1, m.ncols)
    r = clamp(floor(Int, (90 - lat) / m.dlat) + 1, 1, m.nrows)
    return @inbounds m.bits[c, r]
end

island(::Nothing, lon, lat) = true

# Every ring of every polygon, holes included: the even-odd fill reads an
# interior ring as the hole it is.
function landrings(path::AbstractString)
    out = Vector{Vector{Tuple{Float64,Float64}}}()
    ArchGDAL.read(path) do ds
        for feat in ArchGDAL.getlayer(ds, 0)
            _collectrings!(out, ArchGDAL.getgeom(feat))
        end
    end
    return out
end

function _collectrings!(out, g)
    t = GI.geomtrait(g)
    if t isa GI.PolygonTrait
        for r in GI.getring(g)
            push!(out, [(Float64(GI.x(p)), Float64(GI.y(p))) for p in GI.getpoint(r)])
        end
    elseif t isa GI.MultiPolygonTrait
        for p in GI.getgeom(g)
            _collectrings!(out, p)
        end
    end
    return out
end

"""
    rasterize_land(rings, arcsec) -> LandMask

Even-odd scanline fill of `rings` onto the global `arcsec` lattice. A row is
sampled at its centre latitude and every pixel a span touches in longitude is
filled, end pixels included.

`Rasters.rasterize(boundary = :center)` gives a different mask: this fill is a
one-pixel dilation of every coastline in longitude (0.017 % of pixels at 15
arcsec). A store written against this mask resumes only against this mask.

The rings must stay inside [-180, 180]; the sweep has no antimeridian handling.
"""
function rasterize_land(rings, arcsec::Int)
    nrows = (180 * 3600) ÷ arcsec
    ncols = (360 * 3600) ÷ arcsec
    dlat = 180 / nrows
    dlon = 360 / ncols
    X1 = Float64[]; Y1 = Float64[]; X2 = Float64[]; Y2 = Float64[]
    for r in rings, k in 1:(length(r) - 1)
        (xa, ya), (xb, yb) = r[k], r[k + 1]
        ya == yb && continue
        if ya < yb
            push!(X1, xa); push!(Y1, ya); push!(X2, xb); push!(Y2, yb)
        else
            push!(X1, xb); push!(Y1, yb); push!(X2, xa); push!(Y2, ya)
        end
    end
    ne = length(X1)
    rowof(y) = clamp(floor(Int, (90 - y) / dlat) + 1, 1, nrows)
    starts = [rowof(Y2[e]) for e in 1:ne]     # `Y2` is the northern endpoint
    order = sortperm(starts)
    bits = falses(ncols, nrows)
    active = Int[]
    xs = Float64[]
    p = 1
    for j in 1:nrows
        ylat = 90 - (j - 0.5) * dlat
        while p <= ne && starts[order[p]] <= j
            push!(active, order[p]); p += 1
        end
        isempty(active) && continue
        empty!(xs)
        i = 1
        while i <= length(active)
            e = active[i]
            if Y1[e] > ylat                    # wholly north of the scanline now
                active[i] = active[end]; pop!(active); continue
            end
            Y1[e] <= ylat < Y2[e] &&
                push!(xs, X1[e] + (ylat - Y1[e]) * (X2[e] - X1[e]) / (Y2[e] - Y1[e]))
            i += 1
        end
        length(xs) < 2 && continue
        sort!(xs)
        for k in 1:2:(length(xs) - 1)
            c1 = clamp(floor(Int, (xs[k] + 180) / dlon) + 1, 1, ncols)
            c2 = clamp(ceil(Int, (xs[k + 1] + 180) / dlon), 1, ncols)
            c1 <= c2 && (@inbounds bits[c1:c2, j] .= true)
        end
    end
    return LandMask(ncols, nrows, dlon, dlat, bits)
end

"""
    landmask(shapefile, arcsec) -> LandMask or NOMASK

The land polygons of `shapefile` (Natural Earth `ne_10m_land.shp`) rasterised at
`arcsec`. `arcsec = 0` returns [`NOMASK`](@ref) and reads nothing.
"""
function landmask(shapefile::AbstractString, arcsec::Integer)
    arcsec <= 0 && return NOMASK
    isfile(shapefile) || throw(ArgumentError(
        "no land shapefile at $shapefile; pass arcsec = 0 to run unmasked"))
    return rasterize_land(landrings(shapefile), Int(arcsec))
end

"""
    synthetic_tile(sys, tile, mask) -> Vector{Float32}

[`synthetic_elevation`](@ref) at every post of `tile`, ocean posts `NaN32`, in
[`loadtile`](@ref) order. Posts are pixel-is-point: column `i` sits at
`lon_w + i/ncols`, row `j` at `lat_s + 1 - j/nrows`.
"""
function synthetic_tile(sys, tile::DGG.LevelIndex, mask)
    lat_s, lon_w = CD.tilecorner(sys, tile)
    nc = Int(CD.ncols_at(sys, lat_s))
    nrows = Int(CD.lat_intervals(sys))
    out = Vector{Float32}(undef, nc * nrows)
    lons = [lon_w + i / nc for i in 0:(nc - 1)]
    s3 = [sin(3 * deg2rad(l)) for l in lons]
    c7 = [cos(7 * deg2rad(l)) for l in lons]
    @inbounds for j in 0:(nrows - 1)
        lat = (lat_s + 1) - j / nrows
        φ = deg2rad(lat)
        c2, s5 = cos(2φ), sin(5φ)
        base = j * nc
        for i in 1:nc
            out[base + i] = island(mask, lons[i], lat) ?
                            1000 * s3[i] * c2 + 500 * c7[i] * s5 + 100 : NaN32
        end
    end
    return out
end

"""
    SyntheticTiles(sys, mask = NOMASK)

A tile source that fabricates every tile from [`synthetic_elevation`](@ref),
with `mask` deciding which posts are ocean.
"""
struct SyntheticTiles{S<:DGG.CopernicusDEMSystem,M}
    sys::S
    mask::M
end

SyntheticTiles(sys) = SyntheticTiles(sys, NOMASK)

loadtile(source::SyntheticTiles, ordinal::Int) =
    synthetic_tile(source.sys, DGG.LevelIndex(0, ordinal), source.mask)

"""
    SourceMask(mask, g0, listed)

Whether a point has a valid synthetic post under it. A point has none when its
tile is unlisted or when the land mask calls it ocean; the two disagree often
enough (lake islands, mostly-sea tiles, region filters) that an oracle needs
both. `g0` is the level-0 tile grid.
"""
struct SourceMask{M,G}
    mask::M
    g0::G
    listed::Set{Int}
end

@inline function hassource(s::SourceMask, p::GO.UnitSphericalPoint)
    lon, lat = US.GeographicFromUnitSphere()(p)
    island(s.mask, lon, lat) || return false
    t = DGG.cellat(s.g0, p)
    return t !== nothing && Int(t.index) in s.listed
end

"""
    cellsource(grid, c, sm; samples = 8) -> :full, :none or :mixed

Whether destination cell `c` is fed everywhere, nowhere, or partly, sampled at
its centroid and at `samples` points along each boundary edge. `:mixed` is a
coastal cell.
"""
function cellsource(grid, c, sm::SourceMask; samples::Integer = 8)
    has = hassource(sm, DGG.cell_centroid(grid, c))
    ring = DGG.cell_boundary(grid, c)
    for j in eachindex(ring)
        a, b = ring[j], ring[mod1(j + 1, length(ring))]
        for s in 0:(samples - 1)
            hassource(sm, US.slerp(a, b, s / samples)) == has || return :mixed
        end
    end
    return has ? :full : :none
end
