<!-- Stand-in for a photo that does not exist yet: hairline
     rectangle, one centred mono caption. Swap it for an <img>
     with class="dy-photo" when the real asset arrives —
     customer photography is grayscale, never duotoned. -->
<script setup lang="ts">
withDefaults(defineProps<{ caption?: string }>(), {
  caption: 'IMAGE · 1920×1080',
})
</script>

<template>
  <div class="dy-placeholder">
    <div class="dy-placeholder__caption">{{ caption }}</div>
  </div>
</template>
