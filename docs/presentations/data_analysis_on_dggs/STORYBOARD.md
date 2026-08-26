# What is a DGGS?

Discrete global grids are mosaics of cells that cover the entire globe.

The simplest example is a longitude-latitude grid, that's the most familiar grid that exists.

But the interesting DGGSes are all different - Uber's H3, made of hexagons, which they use for movement analysis, comes to mind.

# What are the properties of DGGSes?

They may be hierarchical, meaning that each cell is "composed" of smaller cells (think of a quadtree).

Some of them try to have each cell be equal area on the ellipsoid, others try to have directional fixes, others try to be well positioned for simulation.

DGGSes can be generally better representations of data on the earth, especially when dealing with spread out locations.

# How do people use them?

- Binning point data into cells
- ML style analysis in stencil operators
- Sparse or dense data cubes for traditional zonal etc analysis

# State of DGGS in Julia