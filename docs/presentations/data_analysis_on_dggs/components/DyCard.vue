<!-- Hairline card. `feature` fills it with the softest accent tint;
     that fill is the only thing that ever marks a card as special —
     no coloured left border, no heavy shadow. -->
<script setup lang="ts">
withDefaults(defineProps<{
  label?: string
  title?: string
  feature?: boolean
}>(), { feature: false })
</script>

<template>
  <div class="dy-card" :class="{ 'dy-card--feature': feature }">
    <div v-if="label" class="dy-card__label">{{ label }}</div>
    <div v-if="title" class="dy-card__title">{{ title }}</div>
    <div class="dy-card__body">
      <slot />
    </div>
  </div>
</template>
