import { defineShikiSetup } from '@slidev/types'

/**
 * Syntax highlighting in the JuliaGeo palette.
 *
 * Stock themes bring six or seven hues to a system that allows
 * three families. This one keeps to ink, paper and accent: ink for
 * the code you read, muted for the parts you skim, accent for the
 * parts that carry meaning. Nothing is italic and only definitions
 * go bold.
 *
 * Both themes run on a dark slab. That is not a dark-mode default —
 * the JuliaGeo docs site puts code on `hsl(220, 20%, 9%)` in both
 * modes, and on a projector it is also the better surface: a light
 * block glares and loses its edge against the paper around it.
 * The `light` theme below is therefore light-MODE, not light-BG.
 *
 * Values come from the design system's colors_and_type.css. Do not
 * add a hue that is not in that file.
 */

interface Tokens {
  bg: string
  fg: string
  quiet: string
  muted: string
  accent: string
  accentInk: string
}

function juliageoTheme(name: string, t: Tokens) {
  return {
    name,
    type: 'dark' as const,
    colors: {
      'editor.background': t.bg,
      'editor.foreground': t.fg,
    },
    settings: [
      { scope: ['source', 'text'], settings: { foreground: t.fg } },

      // Skimmed: comments and punctuation recede.
      { scope: ['comment', 'punctuation.definition.comment'], settings: { foreground: t.quiet } },
      { scope: ['punctuation', 'meta.brace', 'meta.delimiter'], settings: { foreground: t.quiet } },

      // Structural: keywords and operators carry the accent.
      {
        scope: ['keyword', 'storage', 'storage.type', 'storage.modifier', 'keyword.control'],
        settings: { foreground: t.accentInk },
      },
      { scope: ['keyword.operator'], settings: { foreground: t.muted } },

      // Data: literals in solid accent, strings quiet.
      { scope: ['constant.numeric', 'constant.language', 'constant.character'], settings: { foreground: t.accent } },
      { scope: ['string', 'string.quoted', 'punctuation.definition.string'], settings: { foreground: t.muted } },

      // Definitions: the one place weight changes.
      { scope: ['entity.name.function', 'support.function', 'meta.function-call'], settings: { foreground: t.fg, fontStyle: 'bold' } },
      { scope: ['entity.name.type', 'entity.name.class', 'support.type', 'support.class'], settings: { foreground: t.accentInk } },
      { scope: ['variable', 'variable.parameter', 'variable.other'], settings: { foreground: t.fg } },
      { scope: ['entity.name.tag'], settings: { foreground: t.accentInk } },
      { scope: ['entity.other.attribute-name'], settings: { foreground: t.muted } },

      { scope: ['invalid'], settings: { foreground: t.fg } },
    ],
  }
}

export default defineShikiSetup(() => ({
  themes: {
    // Light MODE — still a dark slab, per the docs site.
    light: juliageoTheme('juliageo-light', {
      bg: '#12151c',        // hsl(220, 20%, 9%) — --dark-bg
      fg: '#e9ecef',        // gray-200
      quiet: '#7d848f',     // hsl(220, 8%, 53%)
      muted: '#adb5bd',     // gray-500
      accent: '#7dda71',    // green-300 — literals
      accentInk: '#91dd33', // the source's dark-mode lime, for keywords
    }) as any,
    // Dark mode — the slab steps up to --dark-bg-soft so it still
    // separates from a page that is already near-black.
    dark: juliageoTheme('juliageo-dark', {
      bg: '#1c1f26',        // hsl(220, 16%, 13%) — --dark-bg-soft
      fg: '#e9ecef',        // gray-200
      quiet: '#6b7280',     // hsl(220, 9%, 45%)
      muted: '#9ca3ad',     // hsl(220, 9%, 64%)
      accent: '#7dda71',    // green-300
      accentInk: '#91dd33',
    }) as any,
  },
}))
