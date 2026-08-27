<!-- Addition — the closing slide.
     The title layout inverted: mark at the top, the closing line in
     accent, a row of contact labels on a hairline baseline.
     CTAs take a trailing › — never an arrow, never an exclamation.

     Frontmatter:
       layout: end
       links:
         - { label: DOCS,   value: juliageo.org }
         - { label: GITHUB, value: JuliaGeo/GeometryOps.jl }
     Slots:  default → the closing line.  ::sub:: → one line under it. -->
<script setup lang="ts">
defineProps<{
  links?: { label: string, value: string }[]
}>()
</script>

<template>
  <div class="slidev-layout dy-slide dy-end">
    <div class="dy-title__mark" style="margin-bottom: var(--dy-s-8)">
      <img class="dy-title__logo" src="/logos/juliageo-logo.svg" alt="">
      <span class="dy-title__wordmark">JuliaGeo</span>
    </div>

    <div class="dy-end__h1">
      <slot />
    </div>

    <div v-if="$slots.sub" class="dy-title__sub">
      <slot name="sub" />
    </div>

    <div v-if="links?.length" class="dy-end__links">
      <div v-for="link in links" :key="link.label">
        <div class="dy-end__link-label">{{ link.label }}</div>
        <div class="dy-end__link-value">{{ link.value }}</div>
      </div>
    </div>
  </div>
</template>
