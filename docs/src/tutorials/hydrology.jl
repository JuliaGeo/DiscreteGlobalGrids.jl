# # Hydrology: a DEM on an IGEO7 grid
#
# Terrain analysis wants an equal-area grid: slope, flow accumulation and
# catchment area are all areal quantities, and on a lon/lat raster every one of
# them is a function of latitude. This page moves a Copernicus 30 m DEM tile
# over the Alps onto IGEO7 — hexagons, equal-area by construction — and does
# the first step of a flow-routing model on it. The worked example reads the
# native 30 m tile and works at a level whose cells are a couple of pixels
# across — as fine as a tile's worth of them will go on a standard CI runner.
#
# Three calls carry the page: `MultiOrderCoverage` names the cells the tile
# touches, `regrid` fills them from the raster, and `mapneighbors` routes water
# out of every cell.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

import DiscreteGlobalGrids as DGG
import Geomorphometry as GM
import GeoInterface as GI, GeometryOps as GO
using Rasters, RasterDataSources
import ArchGDAL
import Extents
using Statistics
using GLMakie, GeoMakie
using DiscreteGlobalGridsVisualization: dggsurface, dggsurface!, dggpoly, dggpoly!
GLMakie.activate!(inline = true)

# ## Acquiring data

# Let's first get a DEM tile and regrid it to IGeo7.
# We'll use a Copernicus DEM 30m-tile over the Alps for this example.

centre = GI.extent((10.5, 46.5))
path = only(skipmissing(RasterDataSources.getraster(CopernicusDEM; extent = centre)))
Sys.isapple() && Rasters.checkmem!(false) # needed for Apple systems
dem = Raster(path; lazy = false)
## dem = aggregate(mean, dem, 8; progress = false)
# This is what the raster looks like:
plot(dem; axis = (; aspect = DataAspect()))
# Now, let's regrid it to IGeo7.  First, the system:
sys = DGG.IGeo7System()
# We can get the grid of IGeo7 cells that would cover the DEM tile,
# using a `MultiOrderCoverage` query.
# You can specify an integer here, or find the closest matching level automatically:
leaf_level = DGG.levelfor(sys, dem)
# Once we know the leaf level we want, we can query the system at that level with the
# [`MultiOrderCoverage`](@ref) query, to get a cell set that we can use.
region = @time DGG.query(
    sys, 
    DGG.MultiOrderCoverage(Rasters.extent(dem)); 
    level = leaf_level
)
# Here's what this looks like:
f, a, p = plot(dem; axis = (; aspect = DataAspect()))
poly!(a, region; color = :transparent, strokewidth = 2, strokecolor = (:black, 0.5))
f
# This is a nice way to compress the set of cells that would be covered in memory.
# Note that `region` says it has ~44,000 cells.  But when you look at the number of
# cells at level 13,
DGG.CellLookup(region) |> length
# That's a lot of cells!  This optimization helps to decrease memory pressure,
# especially on datasets that don't fit in memory in the first place.

# When we construct a DimArray with this, it will interpret `region` as a one
# level cell axis.  Indices will run from `1:length(CellLookup(region))`,
# linearly, so indexing is done as in a regular vector.

# ## Regridding the DEM
# DiscreteGlobalGrids provides a `regrid` function that will take in a raster
# and some sort of grid - a full level grid, or a region, or a partial grid -
# and return a Raster whose axis is this new grid.
# By default, it returns a DimArray, which we can promote to a Raster (those
# have better handling for missing / NODATA values).

igeo7_dem = @time DGG.regrid(dem; to = region)
igeo7_dem = @time DGG.regrid(dem; to = region, method = DGG.BarycentricPoint())
elevation = Raster(igeo7_dem; missingval = oftype(first(igeo7_dem), NaN))
# Let's now plot this too, using the specialized [`dggsurface`](@ref) recipe for efficiency:
f, a, p = dggpoly(lookup(elevation, DGG.Cells); color = vec(elevation), axis = (; aspect = DataAspect()))
f

# There are also some nice overloads to make this really feel like a surface plot.
# To enhance realism, we'll transform it to "real" coordinates at least.
f, a, p = dggpoly(
    elevation .* 2; # just for effect, since this will be a static plot
    color = vec(elevation), 
    axis = (; type = Axis3, aspect = :data, clip = false)
);
p.transformation.transform_func[] = GeoMakie.create_transform("+proj=ortho +lon_0=10.5 +lat_0=46.5 +datum=WGS84", "+proj=longlat +datum=WGS84")
f
# We can compute the elevation of the cells pretty easily, as well.
extrema(skipmissing(elevation))

# ## Flow direction
#
# Each cell sends its water to the lowest of its neighbours — the first step of
# every flow-routing model. `mapneighbors` walks every cell with its clipped
# one-ring as handles that index the cube; a neighbour whose
# elevation is missing (the coverage overhangs the tile) is skipped. A cell
# with no lower neighbour is a pit.

function downhill(elev, cell, nbrs)
    dest = 0
    zmin = typemax(eltype(elev))
    for n in nbrs
        z = elev[n]
        isnan(z) && continue
        if z < zmin
            zmin = z
            dest = DGG.localindex(n)
        end
    end
    zc = elev[cell]
    (dest == 0 || zmin >= zc) && return (0, NaN32)
    return (dest, zc - zmin)
end

flow, drop = @time DGG.mapneighbors(elevation) do cell, nbrs
    (isnan(elevation[cell]) ? (0, NaN32) : downhill(elevation, cell, nbrs))::Tuple{Int,Float32}
end

# (; n_pits = count(i -> !isnan(elevation[i]) && flow[i] == 0, eachindex(flow)),
#    max_drop = maximum(filter(!isnan, drop)))

# Elevation, and the drop to the downhill neighbour — the drop map picks out
# valley floors as the flat regions and headwalls as the steep ones.
# `CellVector(grid)[shown]` is the covered cells as a vector, which `poly!`
# draws directly, in the same order `igeo7_dem[shown]` follows.


fig = Figure(size = (900, 430))
ax1 = GeoAxis(fig[1, 1]; dest = "+proj=longlat +datum=WGS84", title = "elevation (m)")
p1 = dggsurface!(ax1, lookup(elevation, DGG.Cells); color = vec(elevation), colormap = :terrain)
Colorbar(fig[2, 1], p1; vertical = false)
ax2 = GeoAxis(fig[1, 2]; dest = "+proj=longlat +datum=WGS84", title = "drop to downhill neighbour (m)")
p2 = dggsurface!(ax2, lookup(elevation, DGG.Cells); color = drop, colormap = :magma,
    nan_color = :gray80)
Colorbar(fig[2, 2], p2; vertical = false)
fig

# ## Terrain analysis with Geomorphometry
#
# The cube's axis is already a `CellLookup`: canonical IGEO7 cell identities
# over indexed vector storage. Geomorphometry reads the DGGS neighbourhood
# and physical cell geometry off it directly. TPI is a single local call; flow
# accumulation uses D8 here, with each result expressed as upstream area in
# square metres.

tpi = @time GM.topographic_position_index(elevation)
#
accumulation, directions = @time GM.flowaccumulation(elevation; method = GM.D8())
# accumulation, directions = @time GM.flowaccumulation(elevation; method = GM.DInf())

cell_area = @time GM.cellarea(elevation, first(eachindex(elevation)))
log_cells = log10.(accumulation ./ cell_area)

fig = Figure(size = (900, 430))
ax1 = GeoAxis(fig[1, 1]; dest = "+proj=longlat +datum=WGS84",
    title = "topographic position index (m)")
p1 = dggsurface!(ax1, lookup(tpi, DGG.Cells); color = vec(tpi), colorrange = (-25, 25),
    colormap = :delta)
Colorbar(fig[2, 1], p1; vertical = false)
ax2 = GeoAxis(fig[1, 2]; dest = "+proj=longlat +datum=WGS84",
    title = "D8 flow accumulation")
p2 = dggsurface!(ax2, lookup(log_cells, DGG.Cells); color = vec(log_cells),
    colormap = :devon)
Colorbar(fig[2, 2], p2; vertical = false,
    label = "log₁₀(upstream cell equivalents)")
fig
save("geomorphometry_igeo7.png", fig)

# `MultiOrderCoverage`, `subtree`, `regrid` and `mapneighbors` are interface
# methods, so the regridding and routing portions can use another system.
# `RelativeZ7Cell` and the lazy cell-index iterator provide the IGEO7-specific
# Geomorphometry integration.
