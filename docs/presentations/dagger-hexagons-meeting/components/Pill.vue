<!-- Mono pill. The one place a 999px radius is allowed. -->
<script setup lang="ts">
withDefaults(defineProps<{
  variant?: 'default' | 'accent' | 'solid'
}>(), { variant: 'default' })
</script>

<template>
  <span
    class="dy-pill"
    :class="{
      'dy-pill--accent': variant === 'accent',
      'dy-pill--solid': variant === 'solid',
    }"
  >
    <slot />
  </span>
</template>
