<!-- One feature row: domain marker, title, a line of support.
     `domain` picks the marker — "gen" for the filled one, "sci" for
     the stroked one. Omit it for a plain accent dot.

     The two markers are inherited from the Dyad template, where they
     named Generative and Scientific AI. JuliaGeo has no equivalent
     split, so they are now generic: use them when a deck genuinely
     runs on two domains (vector/raster, algorithm/interface) and
     wants a consistent mark for each. See components/Gen.vue. -->
<script setup lang="ts">
withDefaults(defineProps<{
  title: string
  domain?: 'gen' | 'sci' | 'none'
}>(), { domain: 'none' })
</script>

<template>
  <div class="dy-feature">
    <div class="dy-feature__mark">
      <Gen v-if="domain === 'gen'" />
      <Sci v-else-if="domain === 'sci'" />
      <span v-else class="dy-dot dy-dot--accent" />
    </div>
    <div>
      <div class="dy-feature__title">{{ title }}</div>
      <div v-if="$slots.default" class="dy-feature__body">
        <slot />
      </div>
    </div>
  </div>
</template>
