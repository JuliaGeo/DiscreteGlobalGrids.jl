Global grid systems (DGGS), or non-planar grids, are the new hotness - but how can you get data onto them and analyze? In this talk, we'll show the high-level way to do this, which packages you'll need, as well as how to handle issues like geometries on the boundaries of faces, regridding error, and more.

In all likelihood you are familiar with DGGS already - in simulation, tripolar and cubed-sphere grids are common, as are HEALPIX and other formulations.

Julia now has an up-and-coming ecosystem for global grid systems (commonly called DGGS). The question of getting data to and from such a grid is mostly solved. But once data is on that grid, how do you interpret and analyse it?

One of the most interesting applications here, especially for earth scientists, might be analysing data on the grid it's simulated on - thus removing regridding error, and allowing easier debugging at the simulation level.

By "data analysis", we mean both "traditional" zonal statistics and similar methods as well as more interesting things like applying stencil operations, movement analysis and more.

