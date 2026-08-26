<!-- Sample Venn backbone for the domain-split layout.

     Three ink-stroked lobes, every label outside on a horizontal or
     vertical muted leader, one green-tint band in the triple
     intersection, and a single accent dot as the apex. The centre
     needs no label of its own — the slide's compound title above
     the diagram is the centre label.

     The labels are placeholders. Redraw the geometry for your own
     content; keep the classes.

     Geometry: r = 70, centres on an equilateral triangle of
     circumradius 42 about (240, 152). -->
<template>
  <svg class="dy-diagram" viewBox="0 0 480 300" role="img"
       aria-label="Three-lobe domain diagram">
    <!-- Lobes -->
    <circle class="dy-dg-lobe" cx="240" cy="110" r="70" />
    <circle class="dy-dg-lobe" cx="204" cy="173" r="70" />
    <circle class="dy-dg-lobe" cx="276" cy="173" r="70" />

    <!-- Structural band in the triple intersection, then the apex.
         The centroid sits 42 from each centre, so any radius up to
         70 − 42 = 28 stays inside all three lobes; 22 leaves air. -->
    <circle class="dy-dg-band" cx="240" cy="152" r="22" />
    <circle class="dy-dg-apex" cx="240" cy="152" r="5.5" />

    <!-- Leaders: horizontal or vertical only, 4px clear of the label. -->
    <path class="dy-dg-leader" d="M240 26 V40" />
    <path class="dy-dg-leader" d="M118 173 H134" />
    <path class="dy-dg-leader" d="M346 173 H362" />

    <!-- Lobe labels. The Dyad original hung a domain marker glyph off
         two of the three; that rule was about AI framing and does not
         carry over, so the labels stand alone.

         `text-anchor` does the horizontal placement — middle for the
         top label, end for the left one, start for the right — so a
         label can change length without anything needing recomputing. -->
    <text class="dy-dg-label dy-dg-label--outside" x="240" y="20"
          text-anchor="middle">VECTOR</text>

    <text class="dy-dg-label dy-dg-label--outside" x="112" y="177"
          text-anchor="end">RASTER</text>

    <text class="dy-dg-label dy-dg-label--outside" x="368" y="177">CRS</text>
  </svg>
</template>
