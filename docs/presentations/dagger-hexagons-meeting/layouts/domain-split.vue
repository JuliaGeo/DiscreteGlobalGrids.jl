<!-- Layout 5 — Domain split.
     Centred compound title, one three-lobe diagram, label stacks at
     gutter-left and gutter-right. The title is the centre label for
     the diagram, so the intersection itself carries no text.

     Frontmatter:
       layout: domain-split
       eyebrow: AI
     Slots:  default → the compound title.
             ::left:: / ::right:: → the gutter stacks.
             ::figure:: → the diagram (defaults to <VennDiagram />). -->
<script setup lang="ts">
defineProps<{ eyebrow?: string }>()
</script>

<template>
  <div class="slidev-layout dy-slide dy-split">
    <Eyebrow v-if="eyebrow" class="dy-slide__eyebrow">{{ eyebrow }}</Eyebrow>

    <div class="dy-split__title">
      <slot />
    </div>

    <div class="dy-split__stage">
      <div class="dy-split__aside">
        <slot name="left" />
      </div>

      <div class="dy-split__figure">
        <slot name="figure">
          <VennDiagram />
        </slot>
      </div>

      <div class="dy-split__aside dy-split__aside--right">
        <slot name="right" />
      </div>
    </div>
  </div>
</template>
