<!-- An exported WGLMakie/Bonito figure, embedded as a frame.

     The pages in figures/html are not images — they are standalone
     documents holding a live WebGL scene, drawn out of the same tokens
     as this deck. They carry no title of their own; the slide does
     that.

     The files reach the browser through public/figures, a symlink to
     figures/html, so `src` is always an absolute served path:
     /figures/02-h3-grid.html.

     Named FigureFrame, not Figure: Slidev registers layouts as
     components alongside these, so <Figure> inside layouts/figure.vue
     would resolve to the layout and recurse.

     Props:
       src    served path to the exported page
       label  mono caption shown over the frame until it draws -->
<script setup lang="ts">
import { onBeforeUnmount, onMounted, ref, computed } from 'vue'

withDefaults(defineProps<{
  src: string
  label?: string
}>(), { label: 'FIGURE' })

const box = ref<HTMLElement>()
const frame = ref<HTMLIFrameElement>()
const loaded = ref(false)

// The exported page pins itself to a fixed pixel size — Bonito's
// `resize_to = :parent` measures .dggs-root, which the export styles
// at the Figure's own size. So the frame is laid out at that size and
// scaled to fit whatever the slide gives it, rather than being sized
// to the slide and silently cropped. 960x540 is the fallback for a
// page that reports nothing, and matches the deck canvas.
const natural = ref({ w: 960, h: 540 })
const scale = ref(1)

// Never above 1: the canvas holds real pixels at its own size, so
// scaling up is blur. A figure exported at the body's exact size
// therefore lands 1:1.
const frameStyle = computed(() => ({
  width: `${natural.value.w}px`,
  height: `${natural.value.h}px`,
  transform: `scale(${scale.value})`,
}))

function fit() {
  const el = box.value
  if (!el)
    return
  // clientWidth, not getBoundingClientRect: Slidev scales the whole
  // slide, and the rect would carry that scale into the arithmetic.
  scale.value = Math.min(
    1,
    el.clientWidth / natural.value.w,
    el.clientHeight / natural.value.h,
  )
}

function measureNatural() {
  try {
    const root = frame.value?.contentDocument?.querySelector('.dggs-root') as HTMLElement | null
    if (root?.offsetWidth && root?.offsetHeight)
      natural.value = { w: root.offsetWidth, h: root.offsetHeight }
  }
  catch {
    // cross-origin: keep the fallback
  }
  fit()
}

// `load` fires when the document is parsed, which is a beat before
// WGLMakie has a canvas up: uncovering there shows an empty frame.
// Same origin, so the frame's own document can be polled for one, and
// a ceiling means a figure that never draws still uncovers rather than
// sitting under a label forever.
let timer: number | undefined

function onLoad() {
  measureNatural()
  const until = performance.now() + 20000

  const poll = () => {
    let ready = performance.now() > until
    try {
      ready ||= !!frame.value?.contentDocument?.querySelector('canvas')
    }
    catch {
      ready = true
    }
    if (ready) {
      measureNatural()
      loaded.value = true
    }
    else {
      timer = setTimeout(poll, 100) as unknown as number
    }
  }

  poll()
}

let observer: ResizeObserver | undefined

onMounted(() => {
  fit()
  observer = new ResizeObserver(fit)
  if (box.value)
    observer.observe(box.value)
})

onBeforeUnmount(() => {
  clearTimeout(timer)
  observer?.disconnect()
})
</script>

<template>
  <div ref="box" class="dy-figure">
    <iframe
      ref="frame"
      class="dy-figure__frame"
      :style="frameStyle"
      :src="src"
      :title="label"
      @load="onLoad"
    />
    <div class="dy-figure__cover" :class="{ 'dy-figure__cover--loaded': loaded }">
      <span class="dy-figure__cover-label">{{ label }}</span>
    </div>
  </div>
</template>
