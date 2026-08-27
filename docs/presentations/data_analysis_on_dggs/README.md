# Data analysis on discrete global grids

JuliaCon 2026. A [Slidev](https://sli.dev) deck wired to the JuliaGeo
design system — every colour, radius, shadow and rule comes from that
system's `colors_and_type.css` and `README.md`; nothing here is a
chosen value.

```bash
npm install
npm run dev          # localhost:3030
npm run build        # static site → dist/
npm run export       # PDF → exports/
npm run export:png   # one PNG per slide → exports/png/
```

The figures are Julia, and they are built separately. Instantiate once after
checking out the repository so the deck resolves both local packages:

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project figures/01-what-is_dggs-series.jl   # all four, → figures/html/
```

## Where this came from

Copied from `../../templates/juliageo-template`, which was in turn
forked from `../talk-template`, the Dyad deck. The frame, layouts,
components, canvas maths and type are unchanged; the palette is not.
Token *names* still carry a `--dy-` prefix so the two directories stay
diffable — port a fix from one to the other and you get a clean diff.
Rename them if that stops mattering:

```bash
sed -i '' 's/--dy-/--jg-/g; s/\bdy-/jg-/g' \
  slides.md slide-top.vue layouts/*.vue components/*.vue styles/*.css
```

## The canvas

960×540 — the half-res 16:9 fallback, chosen for one reason: 960px
across 13.333in is exactly 72px/in, so **1pt = 1px**. Every point size
in a presentation spec drops into `styles/tokens.css` unconverted, and
pt is what slide type is specified in because it survives the jump to
projector, PDF and print handout.

To move to the full-res 1920×1080 canvas, set `canvasWidth: 1920` in
the headmatter and double every size in the slide-canvas block at the
bottom of `styles/tokens.css`. Nothing else changes.

Safe area is 32px — 6% of the 540px short side. Nothing crosses it. The
content column sits just outside it at 40px horizontal / 36px top; the
safe line is a floor, not a target. The bottom inset is 48px because the
page number sits on the safe line and copy has to clear it. All four are
`--dy-pad-*` in `styles/tokens.css`.

## Colour

Three families — ink, paper, Julia green — plus two hues that each
appear in exactly one role.

| Role | Light | Dark | Source |
| --- | --- | --- | --- |
| Ink | `#212529` | `#e9ecef` | gray-900 / gray-200 |
| Paper | `#ffffff` | `hsl(220,20%,9%)` | gray-0 / `--dark-bg` |
| Paper-off | `#f8f9fa` | `hsl(220,16%,13%)` | gray-50 / `--dark-bg-soft` |
| Accent (solid) | `#389826` | `#4dc43d` | Julia green / green-400 |
| Accent (letterforms) | `#2c7a1e` | `#91dd33` | green-600 / dark link colour |
| Alt — 2-domain only | `#9558b2` | `#b07fc9` | Julia purple |
| Inline code | `#cb3c33` | `#e8837b` | Julia red |

**Two accent variants, and they invert between modes.** Solid shapes
take `--dy-accent`, thin letterforms take `--dy-accent-ink`. In light
mode the letterform variant is the deeper of the two; in dark mode it
is the brighter. That flip is what keeps single-line type readable
whichever way paper goes.

**Julia purple is the only permitted fourth hue**, and only for a
genuine two-domain split. Blue reads as water on a geospatial slide and
red is already spoken for by inline code.

**Green tints never sit behind a slide.** They appear only inside
diagrams, as structural bands, and as the fill of one feature card.

### Three documented deviations

Everything else is quoted from the source. These three are not, and
each has a reason:

1. **`--dy-hairline` is gray-300 `#dee2e6`, not the source's
   `#e9ecef`.** A projector eats low-contrast detail that a monitor at
   arm's length resolves fine; at `#e9ecef` the table rules and card
   borders vanish off the back wall. Set it to `--dy-hairline-soft`
   for an exact match with the docs site.
2. **`--dy-fg-3` maps to `--dy-muted-2` (gray-600), not `--dy-muted-3`
   (gray-500).** Same reason: `fg-3` carries page numbers, captions,
   card labels and diagram leader lines, all at 10–11px. The raw ramp
   is still 1:1 with the source — only the semantic mapping steps off
   it.
3. **`--dy-alt-ink`, the dark `--dy-alt-*`, and the dark
   `--dy-code-inline-ink` are derived.** The source ships purple-500
   and purple-100 only, and Julia red drops to 2.6:1 on a dark page.
   Each derivation is written out in `styles/tokens.css`.

## Typography

Three faces, split by job rather than by level.

| Role | Face | Token | Where it comes from |
| --- | --- | --- | --- |
| Headlines, all headings, card and feature titles | Rale Grotesk 400/500/700 | `--dy-font-display` | Bundled OTF |
| Body, lists, tables, notes, captions | Inter 100–900 variable | `--dy-font-sans` | Bundled TTF |
| Labels, code, page numbers | JuliaMono 400/500/700 | `--dy-font-mono` | Bundled woff2 |

**Inter is the design system's specified sans**, and it is the better
reading face at 17–19px on a wall: bigger x-height, wider apertures,
drawn for screens. It ships as a single variable file covering the
whole 100–900 axis, so weights interpolate rather than snapping to a
bundled cut, and one request replaces three.

**Rale Grotesk is not in the design system at all.** It is kept for
display type because the deck's headline sizes, tracking and hand-set
line breaks were drawn against it, and because a grotesk headline over
an Inter body is a sharper pairing than Inter over Inter.

The split is enforced by token, not by convention: anything that reads
as a heading takes `--dy-font-display` explicitly, down to `h3` and the
card titles. If you add a heading, say so — inheriting drops it onto
the body face.

**JuliaMono** is the Julia community's typeface and it carries the
glyphs Julia code actually contains — `∂`, `∇`, `⊗`, `α`, subscripts.

Nothing is fetched at runtime. All three faces are bundled in
`public/fonts/` and marked `local` in the headmatter, so the deck
renders with no network at all. Both Inter and JuliaMono are SIL OFL
1.1; their licence files ship alongside them as the licence requires.

### No font-metric constants

Alignment is left to the browser. List markers are typed characters in
the normal inline flow using a hanging indent, so a dash lands where
the body face draws its own dash and a numeral lands on the same
baseline as the words beside it — at any size, in any face, at any
nesting depth. Marker cells in feature rows are one line tall and
centre their contents.

Earlier revisions positioned these absolutely and placed them with
measured offsets (`0.282em` above the baseline, and so on). That works
until the type changes, and then every constant is silently wrong.
There are none left; changing a face requires changing the face.

## Layouts

Six canonical layouts. **Pick one per slide. Do not combine them.**

| `layout:` | Shape | Frontmatter |
| --- | --- | --- |
| `title` | Lockup, accent display H1, mono tripod | `tripod: [A, B, C]`, `lockup` |
| `section` | Eyebrow and a massive accent H2, nothing else | `index`, `eyebrow`, `quiet` |
| `hero-image` | Headline, grayscale hero, kicker triad | `eyebrow` |
| `two-up` | Text column plus one structural diagram | `eyebrow`, `wideLeft` |
| `domain-split` | Compound title over a three-lobe Venn | `eyebrow` |
| `quote` | Mono eyebrow, ink quote, mono attribution | `eyebrow`, `attribution`, `quiet` |

Four additions a talk needs that the canonical six do not cover. They
reuse the same frame and type scale rather than starting a second
system:

| `layout:` | Shape | Frontmatter |
| --- | --- | --- |
| `default` | Header, hairline, prose or a table | `eyebrow`, `title`, `rule`, `quiet` |
| `cards` | A row of hairline cards | `eyebrow`, `title`, `cols` |
| `code` | One code block, optional note column | `eyebrow`, `title` |
| `figure` | Header, then one exported figure | `title`, `figure:`, `label`, `preload: false` |
| `end` | Closing frame and contact labels | `links: [{label, value}]` |

### Nothing sits above the headline

An eyebrow (mono, uppercase) can sit over the headline on several
layouts, and section dividers can carry a running number. **This deck
uses neither** — no slide passes `eyebrow:` or `index:`, so nothing
renders above the headline and it sits flush at the top of the frame.

Both props survive on every layout. Add `eyebrow: WHATEVER` or
`index: '01'` to a slide's frontmatter and it comes back, correctly
spaced, with no other change.

Named slots go in with Slidev's `::name::` syntax:

```markdown
---
layout: two-up
eyebrow: ALGORITHMS
---

# The headline

<Feature title="Vector" domain="gen">One line of support.</Feature>

::right::

<BlockDiagram />
```

### `title:` is special

Slidev reserves `title` for its own navigation, so it never reaches a
layout as a prop. The layouts read it off the frontmatter instead, which
means one `title:` serves as both the slide headline and the TOC entry.
Pass `heading:` when the two should differ.

## Components

| Component | What it is |
| --- | --- |
| `<Gen />` `<Sci />` | The two domain markers — filled green, stroked ink |
| `<Eyebrow>` `<Kicker>` | Mono uppercase signage, tracked 0.16em / 0.12em |
| `<Feature title domain>` | Marker, title, a line of support |
| `<DyCard label title feature>` | Hairline card; `feature` fills it with the softest tint |
| `<Pill variant>` | `default`, `accent`, `solid` |
| `<Stat value label>` | A number and what it counts |
| `<Placeholder caption>` | Hairline stand-in for a figure |
| `<FigureFrame src label>` | An exported `figures/html` page, scaled to fit its box |
| `<Swatch token light dark>` | Colour chip that prints the hex for the active mode |
| `<BlockDiagram />` `<VennDiagram />` | Sample diagrams — redraw the geometry, keep the classes |

Utility classes usable straight from markdown: `dy-eyebrow`, `dy-label`,
`dy-kicker`, `dy-caption` (mono, for figures), `dy-note` (sans, for
running commentary), `dy-rule` (`--hair`, `--thick`, `--dashed`,
`--accent`), `dy-grid` + `dy-grid--2|3|4`, `dy-compact`, `dy-cta`.

### The two domain markers are inherited

`<Gen />` and `<Sci />` come from the Dyad template, where they named
Generative and Scientific AI. JuliaGeo has no equivalent split, so they
are now just "the filled marker" and "the stroked marker" — useful when
a deck genuinely runs on two domains (vector/raster,
algorithms/interface).

A filled 4-point sparkle still reads as "AI" to most audiences. If
that is wrong for your talk, replace the `<path>` in
`components/Gen.vue` with something geospatial — a vertex, a tile, a
graticule. Size, colour and alignment all come from CSS, so nothing
else has to change.

## Figures

The figures are not images. Each one is a Julia script in `figures/`
that exports a standalone WGLMakie/Bonito page into `figures/html/` —
a live WebGL scene, drawn out of the same tokens as the deck
(`figures/00-dggs-theme.jl` mirrors the palette and loads the three
bundled faces).

Package grids are passed directly from `levelgrid(...)` to
`DiscreteGlobalGridsVisualization.dggpoly!`; the recipe reads the current
`cell_boundary` interface and emits one mesh for the complete cell set. Only
the non-DGGS Oceananigans comparison grids and hand-built Tissot geometry use
Makie's generic `poly!` path.

**The figure is the figure; the slide is the slide.** An exported page
carries no title, no subtitle and no rule — the `figure` layout draws
the standard header and hands the rest of the canvas over. A figure
slide takes a title and nothing else; if it needs a sentence of
argument, it is not a figure slide.

They reach the browser through `public/figures`, a symlink to
`figures/html`. Vite serves it in dev and dereferences it on build, so
`dist/` carries real copies and the deck still renders with no network
at all. `figures/html/` is gitignored alongside `dist/` — it is ~100 MB
of generated pages, so rebuild it rather than commit it.

```markdown
---
layout: figure
title: H3 grid                        # the headline, and the TOC entry
figure: /figures/02-h3-grid.html      # served path, i.e. public/figures
label: H3 GRID                        # mono caption held during load
preload: false
---
```

Drop the `title:` and the header goes with it: the figure runs full
bleed, edge to edge.

### Export at the size of the hole

The frame lays the exported page out at its own pixel size and scales
it down to fit, so nothing is ever cropped — but a page larger than
the space it lands in is a page rendered smaller than it could be.

The layout gives the figure everything it can: the rule sits tight
under the title, the body starts at the rule, the figure bleeds past
the content inset to both canvas edges, and the bottom stops on the
safe line rather than the copy inset. That leaves **960×415** under a
title, and the whole 960×540 canvas on a bare slide.

Export at those numbers — `Figure(size = (960, 415))` — and the scale
factor is 1 and the plot fills the hole. A 960×540 page under a title
scales to 0.769 and gives back a third of its area as white margin,
because 16:9 does not fit a 2.3:1 hole. No layout tweak reaches that
last third; only the export size does.

Three sharp edges worth knowing before you add a figure slide:

1. **`preload: false` is not optional.** Slidev keeps every visited
   slide mounted, and each of these pages is 10–15 MB carrying its own
   WebGL context. `preload: false` is what makes a figure mount on
   entry and unmount on exit, so exactly one is ever alive. The cost is
   a beat of load on entry, which the mono `label:` covers until the
   scene has actually drawn.
2. **The path key is `figure:`, not `src:`.** Slidev reserves `src:` in
   frontmatter for pulling in another markdown file and strips it
   before a layout sees it.
3. **The component is `<FigureFrame>`, not `<Figure>`.** Layouts are
   registered as components alongside `components/`, so a `<Figure>`
   inside `layouts/figure.vue` resolves to the layout itself and
   recurses until the stack gives out.

## Files

```
slides.md              the deck
slide-top.vue          page number, drawn over every slide
layouts/               one file per layout
components/            reusable pieces
figures/*.jl           the Julia figure scripts
figures/html/          their exported pages — built, not checked in
public/figures         symlink → figures/html
styles/
  tokens.css           the design system mirrored — the only place values live
  fonts.css            Rale Grotesk + JuliaMono @font-face
  base.css             element defaults
  components.css       reusable pieces
  diagrams.css         diagram primitives
  layouts.css          per-layout frames
setup/shiki.ts         syntax highlighting in the JuliaGeo palette
public/fonts           Inter variable TTF + Rale Grotesk OTFs
                       + JuliaMono woff2, with OFL licences
public/logos           JuliaGeo org logo, GeometryOps package logo
```

## Rules worth not breaking

1. **Rale Grotesk for headings, Inter for body, JuliaMono for
   signage.** Never a system font, never Roboto. All three are bundled
   in `public/fonts/` and marked `local` in the headmatter, so none is
   ever fetched and the deck renders offline.
2. **Mono labels are uppercase and tracked.** 0.12em for labels, 0.16em
   for eyebrows. Signage, not decoration; no weight-mixing in one label.
3. **One accent: Julia green.** Headlines take `accent-ink` (`#2c7a1e`
   light / `#91dd33` dark), solid shapes take `accent` (`#389826` /
   `#4dc43d`). The two invert between modes so single-line type always
   has contrast. Body copy is never accent.
4. **Small radii, not zero radii.** 0 full-bleed, 3px badges and inline
   code, 8px cards and code blocks, 12px modals. Oversized rounding
   reads as dated.
5. **Green tints never sit behind a slide.** They appear only inside
   diagrams, as structural bands.
6. **One apex per diagram.** Exactly one accent dot, thick accent arrow
   or filled accent node claims the eye. If two things compete, the
   diagram is wrong.
7. **Fades and 1–2px translates only.** No bounce, no spring, no
   overshoot — 150ms fast, 220ms base, plain `ease`.
8. **No emoji.** CTAs end in ` ›`, never `→`.

## Code blocks are dark in both modes

The design system puts code on `hsl(220, 20%, 9%)` with light text
whatever the page is doing — it is the most recognisable thing about
the JuliaGeo/VitePress look. It also happens to be the better choice
on a projector, where a light block glares and loses its edge against
the paper around it.

On dark mode the slab steps up to `--dark-bg-soft` and takes a
hairline, so it still separates from a page that is already near-black.
Both are `--dy-code-*` in `styles/tokens.css`; `setup/shiki.ts` carries
the matching token colours.

Inline code is the exception, and the one place a third hue appears: a
gray-100 chip with Julia red type, exactly as the docs site sets it.

## Colour scheme

**This deck is pinned to light mode** — `colorSchema: light` in the
headmatter. That hides Slidev's dark toggle, so the deck renders
identically whatever the presenting machine is set to.

The dark half of the palette is still fully specified in
`styles/tokens.css` and every layout and component handles it. Set
`colorSchema: auto` (follow the system, toggle available) or `dark` (pin
dark) and it all comes back; nothing else needs to change.

## Where the design system lives

`../../.claude/skills/juliageo-design-system/`

Vendored into the repo rather than referenced from wherever it was
downloaded, so the pointer cannot go stale and both templates cite
their design system the same way.

- `colors_and_type.css` — the token file. If `styles/tokens.css`
  disagrees with it, fix `styles/tokens.css`.
- `README.md` — voice, radii, shadows, borders, iconography, imagery.
- `preview/` — live specimens for colour, type, spacing, components.
- `assets/` — the org logo, the GeometryOps logo, the favicon.
- `ui_kits/docs/` — the docs-site prototype the palette was drawn from.
