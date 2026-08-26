<!-- Addition — the plain content slide.
     Not one of the six canonical layouts; a deck needs somewhere to
     put prose, a list or a table. It reuses the same frame and type
     scale rather than inventing a second system.

     Frontmatter:
       eyebrow: TECHNOLOGY
       title: The slide title      # or just write `# Title` in body
       rule: false                 # drop the hairline under the head
       quiet: true                 # paper-off background
     Slots:  default → the body. -->
<script setup lang="ts">
import { computed } from 'vue'
import { useSlideContext } from '@slidev/client'

const props = withDefaults(defineProps<{
  eyebrow?: string
  heading?: string
  rule?: boolean
  quiet?: boolean
}>(), { rule: true, quiet: false })

// Slidev reserves `title` for its own navigation, so it never reaches
// a layout as a prop — read it off the frontmatter instead. That lets
// one `title:` serve as both the slide headline and the TOC entry.
// `heading:` overrides it when the two should differ.
const { $frontmatter } = useSlideContext()
const title = computed(() => props.heading ?? $frontmatter.title)
</script>

<template>
  <div
    class="slidev-layout dy-slide"
    :class="{ 'dy-slide--quiet': quiet }"
  >
    <div v-if="eyebrow || title" class="dy-slide__head">
      <Eyebrow v-if="eyebrow" class="dy-slide__eyebrow">{{ eyebrow }}</Eyebrow>
      <div v-if="title" class="dy-slide__title">{{ title }}</div>
      <hr v-if="rule" class="dy-rule dy-rule--hair dy-slide__rule">
    </div>

    <div class="dy-slide__body">
      <slot />
    </div>
  </div>
</template>
