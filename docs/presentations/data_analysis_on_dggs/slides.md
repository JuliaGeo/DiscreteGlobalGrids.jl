---
# ============================================================
# Data analysis on DGGS — JuliaCon 2026
#
# Headmatter is the JuliaGeo template's, unchanged. Canvas is
# 960x540, which is also the canvas every figure in figures/ is
# drawn on, so an exported figure drops in at exactly 1:1.
# ============================================================
theme: none
title: Data analysis on discrete global grids
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
  - JULIACON 2026
  - 2026-07-29
---

# Data analysis on<br>discrete global grids.

::sub::

Getting data onto a global grid is mostly solved. This talk is about
what happens next — zonal statistics, stencils and movement analysis
on a sphere that no longer pretends to be a rectangle.

---
layout: section
---

# What is a DGGS?

---
title: A mosaic of cells
---

- A **discrete global grid** is a mosaic of cells that covers the
  entire globe.
- The simplest one is the standard **longitude–latitude** grid, which lots of global datasets exist on.
  - UTM and Web Mercator grids are also common.
- But those are unsuitable for many applications, cell areas vary, and the grid has a singularity at the poles.
- DGGS are the extension of this - hexagons, triangles, pentagons, quadrilaterals, etc.

---
layout: figure
title: Latitude-longitude grid
figure: /figures/01-longlat-grid.html
label: LATITUDE-LONGITUDE GRID
preload: false
---

---
layout: figure
title: H3 grid
figure: /figures/02-h3-grid.html
label: H3 GRID
preload: false
---

---
layout: figure
title: Projection distortion
figure: /figures/03-projection-distortion.html
label: PROJECTION DISTORTION
preload: false
---

---
layout: figure
title: More DGGS'es
figure: /figures/10-more-dggs.html
label: MORE DGGSES
preload: false
---

---
layout: section
---

# What are the properties<br>of a DGGS?

---
title: Hierarchy, area, direction
---

- **Hierarchical** — each cell is composed of smaller cells, the way a
  quadtree is
- **Equal area** — every cell measuring the same on the ellipsoid
- **Direction**, or a layout a simulation can step across

<hr class="dy-rule dy-rule--hair">

<p class="dy-note">
No grid gives you all of them; each one picks. What they share is
being a generally better representation of data on the earth —
especially when the locations are spread out.
</p>

---
layout: section
---

# Every cell the same size.<br>Or not.

::note::

One shared colour scale across all five grids: purple is smaller than
the global mean, green is larger, white is the mean itself.

---
layout: figure
title: Cell area — latitude-longitude
figure: /figures/04-cell-area-longlat.html
label: CELL AREA · LATITUDE-LONGITUDE
preload: false
---

---
layout: figure
title: Cell area — H3
figure: /figures/04-cell-area-h3.html
label: CELL AREA · H3
preload: false
---

---
layout: figure
title: Cell area — S2
figure: /figures/04-cell-area-s2.html
label: CELL AREA · S2
preload: false
---

---
layout: figure
title: Cell area — IGEO7
figure: /figures/04-cell-area-igeo7.html
label: CELL AREA · IGEO7
preload: false
---

---
layout: figure
title: Cell area — ISEA4R
figure: /figures/04-cell-area-isea4r.html
label: CELL AREA · ISEA4R
preload: false
---

---
layout: section
---

# How do people use them?

---
layout: cards
title: Three shapes of analysis
cols: 3
---

<DyCard label="01" title="Binning">
Point data binned into cells: counts, sums, densities, etc.
</DyCard>

<DyCard label="02" title="Stencils" feature>
Convolution-style analysis over a cell and its neighbours, reading pixel neighbourhoods.
</DyCard>

<DyCard label="03" title="Cubes">
Sparse or dense data cubes, for more traditional analysis.
</DyCard>

---
layout: section
---

# State of DGGS in Julia.

---