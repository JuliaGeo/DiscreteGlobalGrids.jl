<!-- Domain marker B: 8-ray asterisk, stroked, ink.
     Renders at 1.2x the surrounding type. A block that carries a
     domain gets exactly one marker — this one or <Gen />.

     INHERITED from the Dyad template, where it meant "Scientific
     AI". Generic here — see the note in Gen.vue. The stroke
     treatment is the point, so a Unicode glyph will not do. -->
<template>
  <svg
    class="dy-sci"
    viewBox="0 0 32 32"
    fill="none"
    stroke="currentColor"
    stroke-width="2"
    stroke-linecap="round"
    width="1.2em"
    height="1.2em"
    aria-hidden="true"
  >
    <line x1="16" y1="3" x2="16" y2="29" />
    <line x1="3" y1="16" x2="29" y2="16" />
    <line x1="7" y1="7" x2="25" y2="25" />
    <line x1="25" y1="7" x2="7" y2="25" />
  </svg>
</template>
