<!-- Domain marker A: 4-point sparkle, filled, accent green.
     Renders at 1.2x the surrounding type. A block that carries a
     domain gets exactly one marker — this one or <Sci />.

     INHERITED, and worth a look before you use it. In the Dyad
     template this meant "Generative AI", and a filled sparkle still
     reads that way to most audiences. JuliaGeo has no two-AI split,
     so it is now just "the filled marker". If your deck's two
     domains are geospatial, swap this path for something that says
     so — a vertex, a tile, a graticule — and leave everything else
     alone; the size, colour and alignment all come from CSS.

     Filled as a path rather than typed as U+2726 so it is genuinely
     filled instead of whatever the system symbol font decides. -->
<template>
  <svg
    class="dy-gen"
    viewBox="0 0 32 32"
    fill="currentColor"
    width="1.2em"
    height="1.2em"
    aria-hidden="true"
  >
    <path
      d="M16 2
         C16 2, 17.5 10, 20 13
         C23 16.5, 30 16, 30 16
         C30 16, 23 15.5, 20 19
         C17.5 22, 16 30, 16 30
         C16 30, 14.5 22, 12 19
         C9 15.5, 2 16, 2 16
         C2 16, 9 16.5, 12 13
         C14.5 10, 16 2, 16 2
         Z"
    />
  </svg>
</template>
