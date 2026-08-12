# # Zonal statistics
#
# Which cells of a discrete global grid lie inside a polygon, and what is the mean of a
# field over them? This tutorial answers that for Texas on a whole-globe HEALPix grid,
# twice: first with the built-in `zonal` and `Touching` machinery, then with a custom
# descent of the grid's spatial tree using GeometryOps' spherical predicates — written
# entirely against `SpatialTreeInterface`, so it runs on any grid the package can treeify.

using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids.HEALPix.HealpixLookups
import DimensionalData as DD
import NaturalEarth
import GeometryOps as GO
import GeometryOps.SpatialTreeInterface as STI
using Statistics
using CairoMakie, GeoMakie
CairoMakie.activate!()
nothing

# ## Texas and a field on the globe
#
# Texas comes from Natural Earth's admin-1 dataset. The grid is every HEALPix cell at
# level 7 (`nside = 128`, 196608 cells), which `DGGSGlobeIds` provides lazily — the id
# vector is computed, not stored — and the data is a smooth synthetic field sampled at
# the cell centers.

fc = NaturalEarth.naturalearth("admin_1_states_provinces", 50)
texas = fc.geometry[findfirst(==("Texas"), fc.name)]

level = 7
lookup = HealpixLookup(DGGSGlobeIds(HEALPixDGGS(), level))
field(lon, lat) = 20 - 0.5 * (lat - 25) + 2 * sind(3 * lon)
A = DD.DimArray([field(lon, lat) for (lon, lat) in cell_centers(lookup)], Cells(lookup); name = :tsynth)
nothing

# ## The built-in way
#
# `zonal` reduces over the cells whose center falls inside the geometry, and the
# `Touching` selector subsets to every cell whose polygon intersects it. HEALPix cells
# are equal-area, so the unweighted mean *is* the areal mean — no latitude weights.

tex_mean = only(zonal(mean, A; of = texas))
touching = A[Cells(Touching(texas))]
(; tex_mean, n_touching = length(touching))

# ## The generic descent
#
# What `zonal` does internally generalizes to every DGGS in this package, because the id
# hierarchy *is* a spatial tree — `treeify` hands the grid to you as one, speaking
# `SpatialTreeInterface`. First prepare Texas once, on the spherical manifold:
# this builds the point locators and edge index (~100 ms), and by default also validates
# the rings for spherical use, so do it outside any loop.

prep = GO.prepare(GO.RelateNG(; manifold = GO.Spherical()), texas)
nothing

# Spherical predicates treat polygon edges as great-circle arcs, so a coarse cell's
# four-vertex polygon is already an honest region on the sphere — no densification
# hacks. They also consume the cells' native unit-sphere geometry directly, so the
# whole walk stays on the sphere; Texas is the only lon/lat object, and `prepare`
# ingested it once. The walk prunes any subtree disjoint from Texas, accepts a covered
# subtree wholesale — `STI.depth_first_search` with an always-true predicate enumerates
# its leaf indices, which for a dense whole-globe tree are the cell ordinals — and only
# cells reached at the leaves, the ones straddling the outline, are classified
# individually, by their center. The tree's root is a synthetic whole-sphere node
# (`node_level` of `-1`), which the walk recurses through without a polygon test.

function descend!(interior, boundary, prep, sys, leaf_level, node)
    if node_level(node) >= 0  # the root is synthetic: no cell, no polygon test
        poly = cell_polygon_unitsphere(sys, node_level(node), node_id(node))
        GO.relate_predicate(prep, GO.pred_disjoint(), poly) && return nothing
        if GO.relate_predicate(prep, GO.pred_covers(), poly)
            return STI.depth_first_search(i -> push!(interior, i), Returns(true), node)
        end
    end
    if STI.isleaf(node)
        for (i, _) in STI.child_indices_extents(node)
            center = DGG.cell_center(sys, leaf_level, ordinal_to_cell(sys, leaf_level, i))
            GO.relate_predicate(prep, GO.pred_covers(), center) && push!(boundary, i)
        end
    else
        for child in STI.getchild(node)
            descend!(interior, boundary, prep, sys, leaf_level, child)
        end
    end
    return nothing
end

tree = treeify(DGGSGrid(HEALPixDGGS(), level))
interior, boundary = Int[], Int[]
descend!(interior, boundary, prep, HEALPixDGGS(), level, tree)
(; n_interior = length(interior), n_boundary = length(boundary))

# The descent selects exactly the cells `zonal` averaged over, so the means agree — and
# of 196608 cells on the globe, only the handful of leaf cells along the outline were
# ever handed to an exact predicate; the interior arrived wholesale as accepted
# subtrees, and everything else was pruned at coarse levels.

ords = sort!(vcat(interior, boundary))
descent_mean = mean(A[Cells(ords)])
(; descent_mean, tex_mean, agree = descent_mean ≈ tex_mean)

# Boundary cells are included whenever their center is inside, so they straddle the
# outline; the interior arrived as accepted subtrees. The walk only spoke
# `SpatialTreeInterface` — `isleaf`, `getchild`, `child_indices_extents`,
# `depth_first_search` — plus `node_level`/`node_id` to ask the kernel for each node's
# polygon, so it runs unchanged on anything the package can `treeify`. Its accept and
# prune steps do assume a parent cell geographically contains its descendants — true
# for a congruent grid like HEALPix; for overhanging hierarchies like H3, test the
# conservative `subtree_cap` instead of the parent polygon. And to classify many
# geometries at once against the same tree, `STI.dual_depth_first_search` walks a tree
# of geometries against the grid tree so both sides prune.
#
# The descent never left the sphere; only the plot wants lon/lat.

cellpoly(o) = GO.transform(GO.GeographicFromUnitSphere(),
    cell_polygon_unitsphere(HEALPixDGGS(), level, ordinal_to_cell(HEALPixDGGS(), level, o)))

fig = Figure(size = (700, 620))
ax = GeoAxis(fig[1, 1]; dest = "+proj=longlat +datum=WGS84",
    limits = ((-107.5, -92.5), (25.0, 37.5)), title = "HEALPix level-7 cells selected for Texas")
poly!(ax, cellpoly.(interior); color = :steelblue, strokecolor = :white, strokewidth = 0.5)
poly!(ax, cellpoly.(boundary); color = :orange, strokecolor = :white, strokewidth = 0.5)
poly!(ax, texas; color = :transparent, strokecolor = :black, strokewidth = 2)
axislegend(ax, [PolyElement(color = :steelblue), PolyElement(color = :orange)],
    ["whole subtrees accepted", "boundary cells (center inside)"]; position = :lt)
fig
