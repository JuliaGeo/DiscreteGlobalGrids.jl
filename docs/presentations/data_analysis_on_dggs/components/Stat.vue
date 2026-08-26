<!-- A number and what it counts. Specific numbers beat vague claims. -->
<script setup lang="ts">
defineProps<{ value: string, label: string }>()
</script>

<template>
  <div class="dy-stat">
    <div class="dy-stat__value">{{ value }}</div>
    <div class="dy-stat__label">{{ label }}</div>
  </div>
</template>
