<!-- Layout 3 — Headline plus hero image.
     Eyebrow and accent H2 at the top, one full-bleed grayscale
     image below, a three-column mono-kicker triad under the image.

     Photography is grayscale and hairlined; no colour overlays and
     no duotones. Drop <Placeholder /> into ::image:: until the real
     asset exists.

     Frontmatter:
       layout: hero-image
       eyebrow: CUSTOMERS
     Slots:  default → the H2.  ::image:: → <img class="dy-photo">
             or <Placeholder />.  ::triad:: → three columns. -->
<script setup lang="ts">
defineProps<{ eyebrow?: string }>()
</script>

<template>
  <div class="slidev-layout dy-slide dy-hero">
    <Eyebrow v-if="eyebrow" class="dy-slide__eyebrow">{{ eyebrow }}</Eyebrow>

    <div class="dy-hero__h2">
      <slot />
    </div>

    <div class="dy-hero__figure">
      <slot name="image">
        <Placeholder />
      </slot>
    </div>

    <div v-if="$slots.triad" class="dy-hero__triad">
      <slot name="triad" />
    </div>
  </div>
</template>
