<!-- Layout 2 — Section divider.
     Eyebrow plus a massive accent H2. Whitespace does the work:
     no images, no diagrams, no supporting bullets.

     Frontmatter:
       layout: section
       index: "01"        # optional running number
       eyebrow: MISSION   # single mono uppercase word
       quiet: true        # sit on paper-off instead of paper
     Slots:  default → the H2.  ::note:: → at most one line. -->
<script setup lang="ts">
withDefaults(defineProps<{
  index?: string
  eyebrow?: string
  quiet?: boolean
}>(), { quiet: true })
</script>

<template>
  <div
    class="slidev-layout dy-slide dy-section"
    :class="{ 'dy-slide--quiet': quiet }"
  >
    <div v-if="index || eyebrow" class="dy-section__index">
      <template v-if="index">{{ index }}</template>
      <template v-if="index && eyebrow"> · </template>
      <template v-if="eyebrow">{{ eyebrow }}</template>
    </div>

    <div class="dy-section__h2">
      <slot />
    </div>

    <div v-if="$slots.note" class="dy-section__note">
      <slot name="note" />
    </div>
  </div>
</template>
