<script setup lang="ts">
import { onMounted, ref } from 'vue'

withDefaults(defineProps<{
  src: string
  label?: string
}>(), { label: 'VIDEO' })

const video = ref<HTMLVideoElement>()
const loaded = ref(false)

onMounted(() => {
  video.value?.play().catch(() => {})
})
</script>

<template>
  <div class="dy-figure">
    <video
      ref="video"
      class="dy-video__media"
      :src="src"
      :aria-label="label"
      autoplay
      loop
      muted
      playsinline
      preload="metadata"
      @canplay="loaded = true"
    />
    <div class="dy-figure__cover" :class="{ 'dy-figure__cover--loaded': loaded }">
      <span class="dy-figure__cover-label">{{ label }}</span>
    </div>
  </div>
</template>
