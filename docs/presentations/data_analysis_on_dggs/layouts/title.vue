<!-- Layout 1 — Title.
     Wordmark top-left, huge accent H1 bottom-left, mono tripod on a
     single baseline below it. No page number.

     Frontmatter:
       layout: title
       tripod: [SPEAKER, VENUE, DATE]   # three short mono labels
       lockup: true                     # logo mark beside the wordmark
     Slots:  default → the H1.  ::sub:: → one supporting line. -->
<script setup lang="ts">
withDefaults(defineProps<{
  tripod?: string[]
  lockup?: boolean
}>(), { lockup: true })
</script>

<template>
  <div class="slidev-layout dy-slide dy-title">
    <div class="dy-title__mark">
      <img
        v-if="lockup"
        class="dy-title__logo"
        src="/logos/juliageo-logo.svg"
        alt=""
      >
      <span class="dy-title__wordmark">JuliaGeo</span>
    </div>

    <div class="dy-title__foot">
      <div class="dy-title__h1">
        <slot />
      </div>

      <div v-if="$slots.sub" class="dy-title__sub">
        <slot name="sub" />
      </div>

      <div v-if="tripod?.length" class="dy-tripod">
        <div v-for="item in tripod" :key="item" class="dy-tripod__item">
          {{ item }}
        </div>
      </div>
    </div>
  </div>
</template>
