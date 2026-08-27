<!-- Addition — a code slide.
     Header, then one code block owning the rest of the canvas, with
     an optional note column at the right. Code sits on paper-off
     with a hairline edge, same as any other quiet surface.

     Frontmatter:
       layout: code
       eyebrow: DYAD
       title: A component in 12 lines
     Slots:  default → the fenced code block.
             ::note:: → optional right-hand commentary. -->
<script setup lang="ts">
import { computed } from 'vue'
import { useSlideContext } from '@slidev/client'

const props = defineProps<{
  eyebrow?: string
  heading?: string
}>()

// `title` is reserved by Slidev — see layouts/default.vue.
const { $frontmatter } = useSlideContext()
const title = computed(() => props.heading ?? $frontmatter.title)
</script>

<template>
  <div class="slidev-layout dy-slide dy-code">
    <div v-if="eyebrow || title" class="dy-slide__head">
      <Eyebrow v-if="eyebrow" class="dy-slide__eyebrow">{{ eyebrow }}</Eyebrow>
      <div v-if="title" class="dy-slide__title">{{ title }}</div>
      <hr class="dy-rule dy-rule--hair dy-slide__rule">
    </div>

    <div
      class="dy-code__body"
      :class="{ 'dy-code__body--with-note': !!$slots.note }"
    >
      <div class="dy-code__main">
        <slot />
      </div>
      <div v-if="$slots.note" class="dy-code__note">
        <slot name="note" />
      </div>
    </div>
  </div>
</template>
