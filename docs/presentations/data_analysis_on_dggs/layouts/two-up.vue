<!-- Layout 4 — Two-up diagram.
     Left column: eyebrow, accent H2, feature rows carrying the
     domain markers. Right column: exactly one structural diagram.
     Columns split the content width 50/50 with one col-gap between.

     Frontmatter:
       layout: two-up
       eyebrow: TECHNOLOGY
       wideLeft: true      # 1.15 / 1 split for a narrow diagram
     Slots:  default → the H2 and feature rows.
             ::right:: → the diagram. -->
<script setup lang="ts">
withDefaults(defineProps<{
  eyebrow?: string
  wideLeft?: boolean
}>(), { wideLeft: false })
</script>

<template>
  <div
    class="slidev-layout dy-slide dy-twoup"
    :class="{ 'dy-twoup--wide-left': wideLeft }"
  >
    <div class="dy-twoup__cols">
      <div class="dy-twoup__left">
        <Eyebrow v-if="eyebrow" class="dy-slide__eyebrow">{{ eyebrow }}</Eyebrow>
        <slot />
      </div>

      <div class="dy-twoup__right">
        <slot name="right" />
      </div>
    </div>
  </div>
</template>
