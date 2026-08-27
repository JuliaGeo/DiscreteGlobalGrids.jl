<script setup lang="ts">
import { computed } from 'vue'
import { useSlideContext } from '@slidev/client'

const props = withDefaults(defineProps<{
  video: string
  label?: string
  eyebrow?: string
  heading?: string
  rule?: boolean
}>(), { label: 'VIDEO', rule: true })

const { $frontmatter } = useSlideContext()
const title = computed(() => props.heading ?? $frontmatter.title)
</script>

<template>
  <div class="slidev-layout dy-slide dy-figureslide" :class="{ 'dy-figureslide--bare': !title }">
    <div v-if="eyebrow || title" class="dy-slide__head">
      <Eyebrow v-if="eyebrow" class="dy-slide__eyebrow">{{ eyebrow }}</Eyebrow>
      <div v-if="title" class="dy-slide__title">{{ title }}</div>
      <hr v-if="rule" class="dy-rule dy-rule--hair dy-slide__rule">
    </div>

    <div class="dy-figureslide__body">
      <VideoFrame :src="video" :label="label" />
    </div>

    <div v-if="$slots.default" class="dy-figureslide__note">
      <slot />
    </div>
  </div>
</template>
