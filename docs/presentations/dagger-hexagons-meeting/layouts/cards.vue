<!-- Addition — a row of hairline cards under the standard header.
     Cards are the least decorative container available: hairline
     border, 2px radius, no coloured left bar. Mark at most one as
     `feature` — the accent-tint fill is the only emphasis on offer.

     Frontmatter:
       layout: cards
       eyebrow: TECHNOLOGY
       title: What ships in the box
       cols: 3            # 2, 3 or 4
     Slots:  default → <DyCard> elements. -->
<script setup lang="ts">
import { computed } from 'vue'
import { useSlideContext } from '@slidev/client'

const props = withDefaults(defineProps<{
  eyebrow?: string
  heading?: string
  cols?: number
}>(), { cols: 3 })

// `title` is reserved by Slidev — see layouts/default.vue.
const { $frontmatter } = useSlideContext()
const title = computed(() => props.heading ?? $frontmatter.title)
</script>

<template>
  <div class="slidev-layout dy-slide dy-cards">
    <div v-if="eyebrow || title" class="dy-slide__head">
      <Eyebrow v-if="eyebrow" class="dy-slide__eyebrow">{{ eyebrow }}</Eyebrow>
      <div v-if="title" class="dy-slide__title">{{ title }}</div>
      <hr class="dy-rule dy-rule--hair dy-slide__rule">
    </div>

    <div class="dy-cards__grid" :class="`dy-cards__grid--${cols}`">
      <slot />
    </div>
  </div>
</template>
