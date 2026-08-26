<!-- Page number, drawn on top of every slide.
     Bottom-right of the safe area, JuliaMono Medium 10.5pt, fg-3.
     No prefix, no slash, no total — just `07`.

     Suppressed on the title and end slides, and anywhere the
     frontmatter says `hidePage: true`. -->
<script setup lang="ts">
import { computed } from 'vue'
import { useSlideContext } from '@slidev/client'

const { $page, $frontmatter, $renderContext } = useSlideContext()

const HIDDEN_LAYOUTS = ['title', 'end']

const show = computed(() =>
  $renderContext.value === 'slide'
  && !$frontmatter.hidePage
  && !HIDDEN_LAYOUTS.includes($frontmatter.layout),
)
</script>

<template>
  <div v-if="show" class="dy-pagenum">
    {{ String($page).padStart(2, '0') }}
  </div>
</template>
