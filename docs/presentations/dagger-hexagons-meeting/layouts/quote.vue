<!-- Layout 6 — Full-bleed quote.
     Mono eyebrow at the top, large ink quote centred, mono
     attribution at the bottom. The quote is ink, not accent — accent
     is for headlines, never for body-weight copy.

     Frontmatter:
       layout: quote
       eyebrow: CUSTOMERS
       attribution: NAME · ROLE, COMPANY
       quiet: true      # paper-off full-bleed (default)
     Slots:  default → the quote. -->
<script setup lang="ts">
withDefaults(defineProps<{
  eyebrow?: string
  attribution?: string
  quiet?: boolean
}>(), { quiet: true })
</script>

<template>
  <div
    class="slidev-layout dy-slide dy-quote"
    :class="{ 'dy-slide--quiet': quiet }"
  >
    <Eyebrow v-if="eyebrow">{{ eyebrow }}</Eyebrow>

    <div class="dy-quote__body">
      <div class="dy-quote__text">
        <slot />
      </div>
    </div>

    <div v-if="attribution" class="dy-quote__attr">
      {{ attribution }}
    </div>
  </div>
</template>
