<!-- Addition — one exported figure under a slide title.

     The pages in figures/html are the figure and nothing else: no
     title, no subtitle, no rule. Those belong to the slide, so this
     layout is the `default` header — title, hairline — over a frame
     that takes the rest of the canvas. A figure slide carries a title
     and no subtitle; if it needs a sentence, it is not a figure slide.

     Omit `title:` and the header goes with it and the figure runs
     full bleed.

     Pair it with `preload: false`. Slidev keeps every visited slide
     mounted, and these pages carry a WebGL context and several
     megabytes of geometry each; `preload: false` is what makes a
     figure mount on entry and unmount on exit, so only the figure on
     screen is alive.

     Frontmatter:
       layout: figure
       title: H3 grid                   # headline and TOC entry
       figure: /figures/02-h3-grid.html # public/figures → figures/html
       label: H3 GRID                   # mono caption during load
       preload: false
       heading: ...                     # headline, when it differs
       eyebrow: GRIDS
       rule: false                      # drop the hairline
     Slots:  default → an optional mono note, bottom-left.

     The path key is `figure:`, not `src:` — Slidev reserves `src:` in
     frontmatter for pulling in another markdown file, and strips it
     before the layout ever sees it. Same reason the component is
     <FigureFrame> and not <Figure>: layouts are auto-registered as
     components too, so a <Figure> here would resolve to this file and
     recurse until the stack gives out. -->
<script setup lang="ts">
import { computed } from 'vue'
import { useSlideContext } from '@slidev/client'

const props = withDefaults(defineProps<{
  figure: string
  label?: string
  eyebrow?: string
  heading?: string
  rule?: boolean
}>(), { label: 'FIGURE', rule: true })

// Slidev reserves `title` for its own navigation, so it never reaches
// a layout as a prop — read it off the frontmatter, exactly as the
// default layout does, and one `title:` serves as both the headline
// and the TOC entry.
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
      <FigureFrame :src="figure" :label="label" />
    </div>

    <div v-if="$slots.default" class="dy-figureslide__note">
      <slot />
    </div>
  </div>
</template>
