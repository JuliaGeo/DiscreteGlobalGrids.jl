<!-- A colour chip and its token name.

     Takes both the light and dark hex so the printed value always
     matches the chip you are looking at — a swatch that renders the
     dark palette under a light-mode hex is worse than no swatch. Pass
     only `light` for a token that has one value in both modes.

     Values come from the design system's colors_and_type.css; do not
     type a hex that is not in it (or derived from it in tokens.css). -->
<script setup lang="ts">
import { computed } from 'vue'
import { useDarkMode } from '@slidev/client'

const props = defineProps<{
  token: string
  light: string
  dark?: string
}>()

const { isDark } = useDarkMode()
const hex = computed(() => (isDark.value ? props.dark ?? props.light : props.light))
</script>

<template>
  <div class="dy-swatch-row">
    <div class="dy-swatch" :style="{ background: hex }" />
    <div class="dy-swatch-meta">{{ token }} · {{ hex }}</div>
  </div>
</template>
