<!-- Sample structural diagram: a closed-loop control schematic.

     It exists to show the diagram rules in one place — orthogonal
     ink wires, mono block labels inside, outside labels on 1px
     muted leaders, closed-triangle arrow markers, an accent dashed
     feedback path, and exactly one apex (the accent tap dot).

     Redraw the geometry for your own content; keep the classes. -->
<template>
  <svg class="dy-diagram" viewBox="0 0 420 268" role="img"
       aria-label="Closed-loop control diagram">
    <defs>
      <marker id="dy-arrow" viewBox="0 0 6 4" refX="6" refY="2"
              markerWidth="6" markerHeight="4" orient="auto-start-reverse">
        <path class="dy-dg-marker" d="M0 0 L6 2 L0 4 Z" />
      </marker>
      <!-- The feedback path is accent, so its head is too. -->
      <marker id="dy-arrow-accent" viewBox="0 0 6 4" refX="6" refY="2"
              markerWidth="6" markerHeight="4" orient="auto-start-reverse">
        <path class="dy-dg-marker-accent" d="M0 0 L6 2 L0 4 Z" />
      </marker>
    </defs>

    <!-- Forward path. Orthogonal only; every segment is horizontal. -->
    <path class="dy-dg-wire" d="M16 128 H46" marker-end="url(#dy-arrow)" />
    <path class="dy-dg-wire" d="M82 128 H126" marker-end="url(#dy-arrow)" />
    <path class="dy-dg-wire" d="M214 128 H250" marker-end="url(#dy-arrow)" />
    <path class="dy-dg-wire" d="M338 128 H404" marker-end="url(#dy-arrow)" />

    <!-- Summing junction: arrows in from left and bottom, out right.
         Each sign sits next to the input it belongs to. -->
    <circle cx="64" cy="128" r="18" class="dy-dg-block" />
    <text class="dy-dg-sign" x="55" y="126" text-anchor="middle">+</text>
    <text class="dy-dg-sign" x="70" y="143" text-anchor="middle">−</text>

    <!-- Blocks: rx=1, 1.5px ink, paper fill, mono label centred inside. -->
    <rect class="dy-dg-block" x="126" y="106" width="88" height="44" rx="1" />
    <text class="dy-dg-label" x="170" y="132" text-anchor="middle">CONTROL</text>

    <rect class="dy-dg-block" x="250" y="106" width="88" height="44" rx="1" />
    <text class="dy-dg-label" x="294" y="132" text-anchor="middle">PLANT</text>
    <text class="dy-dg-sub" x="294" y="166" text-anchor="middle">G(s)</text>

    <!-- Apex: the one accent mark on the diagram. -->
    <circle class="dy-dg-apex" cx="370" cy="128" r="5.5" />

    <!-- Feedback: accent, dashed, 0.45 opacity, orthogonal. One path,
         so the dash pattern runs continuously round both corners and
         the head lands on the junction's bottom edge. -->
    <path
      class="dy-dg-feedback"
      d="M370 128 V214 H64 V150"
      marker-end="url(#dy-arrow-accent)"
    />

    <!-- Outside labels on 1px muted leaders. Never diagonal,
         minimum 4px between the label and the leader endpoint. -->
    <path class="dy-dg-leader" d="M370 122 V96" />
    <text class="dy-dg-label--outside dy-dg-label" x="370" y="88"
          text-anchor="middle">MEASURE</text>

    <text class="dy-dg-label--outside dy-dg-label" x="16" y="116">REF</text>
    <text class="dy-dg-label--outside dy-dg-label" x="404" y="116"
          text-anchor="end">OUT</text>
    <text class="dy-dg-label--outside dy-dg-label" x="217" y="234"
          text-anchor="middle">SENSOR</text>
  </svg>
</template>
