# Everything about the run that is FABRICATED, so it can be reviewed apart from
# the real-data path in `copdem_production.jl`:
#
#   SYNTHETIC       the analytic elevation field that stands in for a DEM tile.
#   LandMask        a global land/ocean bitmap.  It exists ONLY to decide which
#                   fabricated pixels are nodata; a real Copernicus GeoTIFF
#                   carries its own nodata and never consults it.
#   synthetic_tile  one tile's worth of SYNTHETIC, ocean posts set to NaN32.
#   verify          the oracle: read written chunks back and hold them against
#                   SYNTHETIC.
#
# Included by `copdem_production.jl`, which supplies `say`, `check` and the
# imports.

# ===========================================================================
# The field
# ===========================================================================

"""
    SYNTHETIC(lon, lat) -> Float64

Stand-in elevation: degrees in, metres out,
`1000 sin(3λ) cos(2φ) + 500 cos(7λ) sin(5φ) + 100`.

It varies on a scale four orders of magnitude larger than a level-12 cell, so a
conservative (area-weighted) cell mean is the field's value at the cell centre to
well within a metre. That is what makes it an oracle — and it stays one over a
part-ocean cell, since a weighted mean of a near-constant field is that constant
however the weights renormalise.
"""
SYNTHETIC(lon, lat) = let λ = deg2rad(lon), φ = deg2rad(lat)
    1000 * sin(3λ) * cos(2φ) + 500 * cos(7λ) * sin(5φ) + 100
end

SYNTHETIC(p::GO.UnitSphericalPoint) = SYNTHETIC(US.GeographicFromUnitSphere()(p)...)

# ===========================================================================
# The land mask
# ===========================================================================

"""
    LandMask(ncols, nrows, dlon, dlat, bits)

A global land/ocean bitmap on the plate-carrée lattice: `bits[col, row]` is true
over land, column 1 starting at 180°W and row 1 at 90°N. A lookup is two
divisions and a bit read.

Synthetic-only. It decides which fabricated pixels are nodata, which is what
puts partly-fed cells at every coastline — the case that makes the conservative
weights renormalise over a cell's land fraction.
"""
struct LandMask
    ncols::Int
    nrows::Int
    dlon::Float64
    dlat::Float64
    bits::BitMatrix
end

"No mask: every pixel of a listed tile is valid."
const NOMASK = nothing

"Is (`lon`, `lat`) over land? True everywhere when there is no mask."
@inline function island(m::LandMask, lon::Real, lat::Real)
    c = clamp(floor(Int, (lon + 180) / m.dlon) + 1, 1, m.ncols)
    r = clamp(floor(Int, (90 - lat) / m.dlat) + 1, 1, m.nrows)
    return @inbounds m.bits[c, r]
end

island(::Nothing, lon, lat) = true

"""
    landrings(path) -> Vector{Vector{Tuple{Float64,Float64}}}

Every ring of every polygon in a land shapefile, exterior and interior alike.
Holes need no special case: the fill below is even-odd, which reads an interior
ring as the hole it is.
"""
function landrings(path::AbstractString)
    out = Vector{Vector{Tuple{Float64,Float64}}}()
    ArchGDAL.read(path) do ds
        layer = ArchGDAL.getlayer(ds, 0)
        for feat in layer
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

Even-odd scanline fill of `rings` onto the global `arcsec` lattice.

NOT `Rasters.rasterize`, and deliberately so. This fill samples a row at its
centre latitude but fills every pixel a span TOUCHES in longitude, end pixels
included — `:center` in y, `:touches` in x. Measured on the production 15-arcsec
lattice the two masks differ in 646 280 of 3.73e9 pixels (0.017 %): all but one
are pixels this fill calls land and `boundary = :center` calls ocean, they run
1.01 pixels long on average, and 98.8 % of them touch a pixel `:center` also
calls land. It is a one-pixel dilation of every coastline in longitude, nothing
else. The store already on disk was written against THIS mask and resuming it
has to reproduce the same coastal values, so the scanline stays for the current
store; `Rasters.rasterize` is the right choice for the next one.

Edges are bucketed by the first scanline they reach and swept with an active
list, so cost goes as (edge, row) incidences rather than rings x rows. Natural
Earth is clipped to [-180, 180], so no edge wraps the antimeridian and the sweep
needs no seam handling.
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

Build the mask and report it. `arcsec = 0` disables masking entirely.
"""
function landmask(shapefile::AbstractString, arcsec::Integer)
    if arcsec <= 0
        say("land mask: DISABLED (maskarcsec=0), every pixel of a listed tile is valid")
        return NOMASK
    end
    t0 = time()
    rings = landrings(shapefile)
    m = rasterize_land(rings, Int(arcsec))
    say(@sprintf("land mask: %d x %d at %d arcsec from %d rings / %d vertices, %.2f%% land by lattice cell, %.0f MiB, %s",
        m.ncols, m.nrows, arcsec, length(rings), sum(length, rings),
        100 * count(m.bits) / length(m.bits), sizeof(m.bits.chunks) / 2^20,
        secs(time() - t0)))
    return m
end

# ===========================================================================
# Fabricating one tile
# ===========================================================================

"""
    synthetic_tile(sys, tile, mask) -> (Vector{Float32}, nland)

[`SYNTHETIC`](@ref) at every post of `tile`, ocean posts set to `NaN32`, in the
tile's own index order: raster row `j` north first, column `i` west to east,
`i` fastest. Also returns how many posts came out land.

Posts are pixel-is-point: column `i` sits at `lon_w + i/ncols`, row `j` at
`lat_s + 1 - j/nrows`.
"""
function synthetic_tile(sys, tile::DGG.LevelIndex, mask)
    lat_s, lon_w = CD.tilecorner(sys, tile)
    nc = Int(CD.ncols_at(sys, lat_s))
    nrows = Int(CD.lat_intervals(sys))
    out = Vector{Float32}(undef, nc * nrows)
    lons = [lon_w + i / nc for i in 0:(nc - 1)]
    s3 = [sin(3 * deg2rad(l)) for l in lons]
    c7 = [cos(7 * deg2rad(l)) for l in lons]
    nland = 0
    @inbounds for j in 0:(nrows - 1)
        lat = (lat_s + 1) - j / nrows
        φ = deg2rad(lat)
        c2, s5 = cos(2φ), sin(5φ)
        base = j * nc
        for i in 1:nc
            if island(mask, lons[i], lat)
                out[base + i] = 1000 * s3[i] * c2 + 500 * c7[i] * s5 + 100
                nland += 1
            else
                out[base + i] = NaN32
            end
        end
    end
    return out, nland
end

# ===========================================================================
# The oracle
# ===========================================================================

"How many points to sample along each edge of a destination cell's boundary."
const BOUNDARY_SAMPLES = 8

"""
    SourceMask(mask, g0, listed)

Does a point have a valid source pixel under it? Two independent ways to have
none, and the oracle needs both: the point's tile may not be on the Copernicus
list at all, or the land mask may call the post ocean.

Either test alone is wrong. Natural Earth calls points land that Copernicus
never listed a tile for — lake islands, coastlines the two datasets disagree
about, tiles a region filter cut — and Copernicus lists tiles that are mostly
sea. A destination cell over either is `NaN`, correctly.
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
    cellsource(g12, c, sm) -> :full, :none or :mixed

Is destination cell `c` fed everywhere, nowhere, or only partly? Sampled at its
centroid and at [`BOUNDARY_SAMPLES`](@ref) points along each boundary edge.

`:mixed` — a coastal cell — is the case this run exists to exercise, and the one
the analytic oracle is only approximately right about, so those cells are checked
against the field's RANGE across the cell rather than its centre value: a
renormalised mean is a mean over some subset of the cell, and the field is not
exactly constant on one.
"""
function cellsource(g12, c, sm::SourceMask)
    has = hassource(sm, DGG.cell_centroid(g12, c))
    ring = DGG.cell_boundary(g12, c)
    for j in eachindex(ring)
        a, b = ring[j], ring[mod1(j + 1, length(ring))]
        for s in 0:(BOUNDARY_SAMPLES - 1)
            hassource(sm, US.slerp(a, b, s / BOUNDARY_SAMPLES)) == has || return :mixed
        end
    end
    return has ? :full : :none
end

"""
    verify(config, sys7, layout, chunks, sm, real, dem, written)

Read written chunks back out of the store and hold them against
[`SYNTHETIC`](@ref). Five claims, in the order they can fail:

  * every value is finite or `NaN`, nothing else;
  * a cell with no valid source pixel under it is `NaN`, and a fully-fed cell
    matches the field at its centre to within the field's curvature across a
    cell;
  * a coastal cell lies inside the field's range across that cell — which is
    what says the weights renormalised over the fed fraction rather than
    diluting toward zero;
  * a real GeoTIFF tile's own cells are finite, the S90 pole tiles included;
  * a chunk nobody wrote reads back all `NaN`.
"""
function verify(config, sys7, layout, chunks, sm, real, dem, written)
    level, ancestor, storepath = config.level, config.ancestor, config.store
    g12 = DGG.levelgrid(sys7, level)
    g0 = DGG.levelgrid(dem.builder.sys, 0)
    n = min(config.checkchunks, length(chunks))
    picked = chunks[round.(Int, range(1, length(chunks); length = n))]
    say("verify: reading back $n of $(length(chunks)) chunks")
    anymixed = false
    for ch in picked
        a = DGG.columncell(layout, ch)
        stack = DGG.dggread(storepath; ancestors = [a])
        vals = collect(stack[:elevation])
        cells = collect(DD.lookup(stack[:elevation], DGG.Cells))
        h = DGG.columnlength(layout, ch)
        check("chunk $ch: read back $h values",
            length(vals) == h && length(cells) == h;
            detail = "$(length(vals)) values, $(length(cells)) cells")

        nfin = count(isfinite, vals)
        nnan = count(isnan, vals)
        check("chunk $ch: every value is finite or NaN", nfin + nnan == length(vals);
            detail = "$nfin finite, $nnan NaN")

        # A stride, not every cell: the boundary walk is 8 slerps per edge and a
        # chunk is 823 543 cells. The stride is prime to 7 so it does not land on
        # one subtree's worth of siblings.
        stride = max(1, length(vals) ÷ 400)
        nfull = nnone = nmixed = 0
        errs = Float64[]
        bracketfail = 0
        nanfail = Int[]
        finfail = Int[]
        for k in 1:stride:length(vals)
            c = cells[k]
            # Only cells fed exclusively by synthetic tiles have an oracle.
            Int(DGG.cellat(g0, DGG.cell_centroid(g12, c)).index) in keys(real) && continue
            kind = cellsource(g12, c, sm)
            if kind === :none
                nnone += 1
                isnan(vals[k]) || push!(nanfail, k)
            elseif kind === :full
                nfull += 1
                isfinite(vals[k]) || (push!(finfail, k); continue)
                push!(errs, abs(vals[k] - SYNTHETIC(DGG.cell_centroid(g12, c))))
            else
                nmixed += 1
                isfinite(vals[k]) || continue
                samples = SYNTHETIC.(DGG.cell_boundary(g12, c))
                push!(samples, SYNTHETIC(DGG.cell_centroid(g12, c)))
                lo, hi = extrema(samples)
                tol = 1e-3 + 0.05 * (hi - lo)
                (lo - tol <= vals[k] <= hi + tol) || (bracketfail += 1)
            end
        end
        nmixed > 0 && (anymixed = true)
        say("chunk $ch sampled: $nfull fully fed / $nnone unfed / $nmixed coastal")
        check("chunk $ch: unfed cells are NaN", isempty(nanfail);
            detail = "$(length(nanfail)) of $nnone unfed cells are not NaN")
        check("chunk $ch: fully-fed cells are finite", isempty(finfail);
            detail = "$(length(finfail)) of $nfull fully-fed cells are NaN")
        if !isempty(errs)
            mx = maximum(errs)
            check("chunk $ch: fully-fed cells match the analytic field", mx < 1.0;
                detail = @sprintf("max %.3e m, RMS %.3e m over %d cells", mx,
                    sqrt(Statistics.mean(abs2, errs)), length(errs)))
        end
        check("chunk $ch: coastal cells bracket the field", bracketfail == 0;
            detail = "$bracketfail of $nmixed coastal cells outside their own field range")
        if nnan > 0 && nfin > 0
            say("chunk $ch holds BOTH land and nodata: $nfin finite, $nnan NaN " *
                @sprintf("(%.1f%%)", 100 * nnan / length(vals)))
        end
    end
    # A global run must meet the coast; a region need not — an all-Antarctica box
    # is 100 % valid and has nothing to renormalise.
    if config.region === nothing
        check("some sampled cell straddles the coast", anymixed;
            detail = anymixed ? "the renormalising path ran" : "no coastal cell was sampled")
    else
        say("coastal cells sampled in this region: " *
            (anymixed ? "yes, the renormalising path ran" :
             "none — a region with no coastline has nothing to renormalise"))
    end

    # Real tiles, the pole in particular: a GLO tile keeps its full degree of
    # longitude but that degree pinches to nothing at 90°, so the S90 row is 360
    # slivers meeting at a point — the one degenerate source geometry.
    for (o, _) in sort!(collect(real))
        t = DGG.LevelIndex(0, o)
        lat, lon = CD.tilecorner(dem.builder.sys, t)
        c = DGG.cellat(g12, US.UnitSphereFromGeographic()((lon + 0.5, lat + 0.5)))
        c === nothing && continue
        ch = DGG.columnindex(layout, DGG.ancestor(sys7, c, ancestor))
        ch in written || continue
        stack = DGG.dggread(storepath; ancestors = [DGG.columncell(layout, ch)])
        vals = collect(stack[:elevation])
        k = findfirst(==(c), collect(DD.lookup(stack[:elevation], DGG.Cells)))
        check("real tile $(tilestem(dem.builder.sys, t)) is finite at its centre",
            k !== nothing && isfinite(vals[k]);
            detail = k === nothing ? "its centre cell is not in chunk $ch" :
                     @sprintf("%.2f m in chunk %d", vals[k], ch))
    end

    # A chunk nobody wrote: the first level-`ancestor` chunk outside the
    # covering, which by construction meets no listed tile.
    inset = Set(chunks)
    empty = findfirst(i -> !(i in inset), 1:DGG.ncells(sys7, ancestor))
    if empty !== nothing
        a = DGG.columncell(layout, empty)
        vals = collect(DGG.dggread(storepath; ancestors = [a])[:elevation])
        check("unwritten ocean chunk $empty reads back NaN", all(isnan, vals);
            detail = "$(count(isnan, vals)) of $(length(vals)) NaN")
    end
    return nothing
end
