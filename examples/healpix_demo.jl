# Demo: EOPF-convention HEALPix vector-data-cube — zonal stats + stencils on partial coverage.
#
# Environment: this demo additionally needs NaturalEarth and Rasters, which are
# intentionally NOT DiscreteGlobalGrids dependencies (they are only used here, to
# fetch reference country polygons and to cross-check against a regular lon/lat
# grid). Run it in an environment that has DiscreteGlobalGrids plus those two.
# NaturalEarth downloads its data at runtime, so the first run needs network access.
import Healpix, NaturalEarth, GeoInterface as GI, GeometryOps as GO
import DimensionalData as DD
using DimensionalData: Lookups
using Statistics
using DiscreteGlobalGrids
using DiscreteGlobalGrids.HEALPix.HealpixLookups

level = 7
res = Healpix.Resolution(2^level)
npix = 12 * 4^level

# --- build partial coverage: cells over a European window (the "country crop" case) ---
inwindow(lon, lat) = -12 <= lon <= 32 && 34 <= lat <= 62
covered = Int64[p for p in 0:npix-1 if inwindow(HealpixLookups._cell_center_lonlat(res, p)...)]
l = HealpixLookup(covered; level)
println("coverage: $(length(covered)) of $npix cells ($(round(100length(covered)/npix; digits=2))%)")

# --- synthetic temperature-like field on cell centers ---
field(lon, lat) = 20 - 0.5 * (lat - 34) + 2 * sind(3 * lon)
A = DD.DimArray([field(lon, lat) for (lon, lat) in cell_centers(l)], Cells(l); name=:tsynth)

# --- zonal means over real country polygons ---
countries = NaturalEarth.naturalearth("admin_0_countries", 110)
sel = ["Germany", "France", "Spain", "Italy", "Poland", "Austria", "Switzerland"]
rows = findall(in(sel), countries.NAME)
geoms = countries.geometry[rows]
zs = zonal(mean, A; of=geoms)
println("\nzonal means (HEALPix, equal-area => unweighted):")
for (n, z) in zip(countries.NAME[rows], zs)
    println("  ", rpad(n, 14), round(z; digits=3))
end

# --- independent check #1: brute force over stored centers ---
centers = cell_centers(l)
for (n, g, z) in zip(countries.NAME[rows], geoms, zs)
    bf = mean(parent(A)[findall(c -> GO.contains(g, c), centers)])
    @assert z ≈ bf "zonal mismatch for $n"
end
println("\nbrute-force agreement: OK")

# --- independent check #2: Rasters.zonal on a regular lon/lat grid of the same field ---
import Rasters
using Rasters: X, Y
lons, lats = -12:0.25:32, 34:0.25:62
rast = Rasters.Raster([field(lon, lat) for lon in lons, lat in lats],
                      (X(lons), Y(lats)); name=:tsynth, missingval=NaN)
rz = Rasters.zonal(mean, rast; of=geoms, boundary=:center, progress=false)
println("\nHEALPix zonal vs Rasters.zonal on a 0.25° lon/lat grid:")
for (n, a, b) in zip(countries.NAME[rows], zs, rz)
    println("  ", rpad(n, 14), " healpix: ", round(a; digits=3), "  rasters: ", round(b; digits=3))
    @assert abs(a - b) < 0.2 "HEALPix and lon/lat zonal means disagree for $n"
    # (different discretizations of the same smooth field; equal-area vs lon/lat-weighted)
end

# --- selector sugar from the lookup itself ---
germany = countries.geometry[findfirst(==("Germany"), countries.NAME)]
@assert mean(A[Cells(Lookups.Contains(germany))]) ≈ zs[findfirst(==("Germany"), countries.NAME[rows])]
println("\nA[Cells(Contains(germany))] == zonal(mean) for Germany: OK")

# --- stencil ops: smoothing + Laplacian ---
nbi = HealpixLookups.neighbor_indices(l)
smooth = stencil((c, nbs) -> mean(vcat(c, nbs)), A; nbidx=nbi)
lap = stencil((c, nbs) -> isempty(nbs) ? zero(c) : mean(nbs) - c, A; nbidx=nbi)
println("\nstencil ops:")
println("  field var: ", round(var(parent(A)); digits=4), " -> smoothed: ", round(var(parent(smooth)); digits=4))
println("  max |laplacian| (smooth field, should be small): ", round(maximum(abs, parent(lap)); digits=5))
nfull = count(t -> all(>(0), t), nbi)
println("  cells with full 8-neighborhood: $nfull / $(length(nbi)) (rest are coverage-boundary cells)")
