---
theme: none
title: Dagger for Global Grids
titleTemplate: '%s'
author: Anshul Singhvi
info: false

aspectRatio: 16/9
canvasWidth: 960
colorSchema: light
transition: fade

lineNumbers: false
selectable: true
drawings:
  persist: false

fonts:
  sans: 'Inter'
  serif: 'Rale Grotesk'
  mono: 'JuliaMono'
  local: 'Inter,Rale Grotesk,JuliaMono'
  weights: '400,500,700'
  italic: false

layout: title
tripod:
  - ANSHUL SINGHVI
  - DAGGER × JULIAGEO
  - MEETING PRIMER
---

# Dagger for<br>Global Grids.

::sub::

A ten-minute introduction to the work hidden inside global grids.

---
layout: figure
title: A DGGS gives every part of the globe an address
figure: /figures/10-more-dggs.html
label: DISCRETE GLOBAL GRID SYSTEMS
preload: false
---

DIFFERENT CELLS · DIFFERENT TOPOLOGIES · THE SAME GLOBAL IDEA

---
layout: figure
title: Each cell sits inside a hierarchy
figure: /figures/10b-igeo7-hierarchy.html
label: IGEO7 · THREE LEVELS
preload: false
---

COARSE ADDRESSES CONTAIN FINER ADDRESSES

---
layout: figure
title: The hierarchy becomes a storage layout
figure: /figures/11a-storage-subtree.html
label: LEVEL 5 → LEVEL 12
preload: false
---

ONE PARENT CELL · ONE CONTIGUOUS DESCENDANT INTERVAL

---
layout: figure
title: Address order becomes a space-filling curve
figure: /figures/11c-z7-curve.html
label: Z7 ADDRESS ORDER
preload: false
---

SUCCESSIVE CENTROIDS · ONE CONTIGUOUS WALK

---
layout: figure
title: Storage chunks are polygons, not raster blocks
figure: /figures/11b-storage-globe.html?v=2
label: GLOBAL CHUNKING
preload: false
---

LEVEL-5 CHUNKS · LEVEL-12 CELL BOUNDARIES

---
layout: figure
title: First, build and apply one regridder
figure: /figures/12-regrid-single.html?v=2
label: ONE-MAP REGRIDDING
preload: false
---

SOURCE PIXELS → SPARSE WEIGHTS → DESTINATION CELLS

---
layout: figure
title: Then reuse it across a long data cube
figure: /figures/13-regrid-timeseries.html?v=3
label: ITS_LIVE-SCALE CUBE
preload: false
---

120 M PIXELS · MANY VARIABLES · MANY TIMESTAMPS

---
layout: figure
title: We already know the chunk dependencies
figure: /figures/14-regrid-chunkgraph.html?v=4
label: REGRIDDING DEPENDENCY GRAPH
preload: false
---

DESTINATION CHUNK → CANDIDATE SOURCE CHUNKS

---
layout: video
title: A stencil advances in storage order
video: /video/15-stencil-order.mp4
label: ONE-RING MEAN
preload: false
---

READ HALO · COMPUTE · WRITE · ADVANCE

---
layout: figure
title: The grid can describe every halo read
figure: /figures/16-halo-mapping.html?v=4
label: CHUNK AND HALO MAPPING
preload: false
---

OWNED CELLS · HALO CELLS · SUPPLYING CHUNKS

---
layout: video
title: Both workloads have the same runtime shape
video: /video/17-common-shape.mp4
label: COMMON TASK SHAPE
preload: false
---

DOMAIN PLANS THE WORK · DAGGER RUNS IT

---
layout: figure
title: What should cross the boundary?
figure: /figures/18-boundary-options.html
label: INTEGRATION CONTRACTS
preload: false
---

CHUNK GRAPH · CELL MAPPING · TASK PLAN

---
layout: end
links:
  - { label: GRID, value: JuliaGeo/DiscreteGlobalGrids.jl }
  - { label: RUNTIME, value: JuliaParallel/Dagger.jl }
---

# Start with one<br>thin vertical slice.

::sub::

One grid · one operation · one multiprocess execution path.
