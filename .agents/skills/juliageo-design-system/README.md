# JuliaGeo Design System

A clean, technical design system for the JuliaGeo open-source geospatial ecosystem. Built primarily for documentation websites, infographics, and presentation slides. Spring green themed, highly readable, developer-first.

---

## Sources

- **Logo SVG**: `assets/juliageo-logo.svg` — JuliaGeo organization logo (three circles + geographic map outline)
- **Codebase**: [JuliaGeo/GeometryOps.jl](https://github.com/JuliaGeo/GeometryOps.jl) — primary package in the JuliaGeo ecosystem
- **Docs CSS**: `docs/src/.vitepress/theme/style.css` within the repo — VitePress theme overrides with brand colors
- **Docs assets**: `docs/src/assets/logo.png` — GeometryOps package logo

---

## Product Context

**JuliaGeo** is the GitHub organization that hosts Julia's geospatial software ecosystem. It includes packages for:
- **GeometryOps.jl** — blazing-fast 2D geometry operations in pure Julia (GIS use cases: intersection, union, simplification, projection, polygonization)
- **GeoInterface.jl** — a common interface standard that all geo packages implement
- Various extension packages (LibGEOS, Proj, FlexiJoins)

The primary surfaces this design system targets:
1. **Documentation website** — VitePress-based, with tutorials, API reference, and literate source code
2. **Infographics & diagrams** — geometry illustrations, benchmark charts
3. **Presentations** — technical talks for the Julia/geospatial community

---

## CONTENT FUNDAMENTALS

### Voice & Tone
- **Direct and precise** — technical writing, no fluff. "Blazing fast geometry operations in pure Julia."
- **Welcoming to contributors** — "We welcome contributions, either as pull requests or discussion on issues!"
- **Honest about limitations** — warnings are upfront: "This package is still under heavy development! Use with care."
- **Scientifically grounded** — references to OGC methods, CRS, GeoInterface; assumes a technically literate reader

### Casing
- Product names: Title Case — `GeometryOps.jl`, `GeoInterface.jl`, `Julia`
- Function/method names: `snake_case` in code, backtick-wrapped in prose
- Section headers in docs: Sentence case — "How to navigate the docs"
- UI labels: Sentence case

### Person / Voice
- Uses **"we"** for the project team ("We are focusing primarily on 2/2.5D geometries")
- Addresses users in **second person** implicitly through imperative docs ("See the API page")
- No emoji in prose; emoji appear only in CI badges and icon slots on the homepage

### Writing Style
- Short paragraphs, bulleted lists for method catalogs
- Code blocks are preferred for examples over inline references
- Warnings use `> [!WARNING]` callout syntax
- Links are descriptive: "See its integrations page" not "click here"

---

## VISUAL FOUNDATIONS

### Color
- **Primary brand**: Julia Green `#389826` — used for buttons, links, accents, borders
- **Light green**: `#dcf5d7` — background tints, logo circle fill
- **Julia Blue** `#4063D8`, **Julia Purple** `#9558B2`, **Julia Red** `#CB3C33` — secondary accents; rarely used for UI, but appear in the logo and hero imagery
- **Link color (light mode)**: `#CB3C33` (Julia Red) — a deliberate departure from the green primary
- **Link color (dark mode)**: `#91dd33` — bright lime green
- **Warning**: `#b45309` (dark amber) — high-contrast against light backgrounds; old `#dccc50` was too pale to read

### Typography
- **Body / UI sans**: `Inter` (user specified) — clean, highly readable, variable weight
- **Docs sans**: historically `Barlow` + `Space Grotesk` in the VitePress site
- **Monospace**: `Space Mono` — used for all code, terminals, filenames
- Type scale is conservative: body 16px, relaxed leading (1.65), tight headings
- Headings: semibold/bold Inter, tracking tight (`-0.025em`)
- Labels: uppercase, widest tracking, small size — used for category metadata

### Backgrounds
- Light mode: pure white `#ffffff` or near-white `#f8f9fa`
- Dark mode: `hsl(220, 20%, 9%)` — cool near-black (blue-tinted)
- No texture, no patterns, no gradient backgrounds in UI
- Hero section only: gradient glow blur behind logo image (decorative, blurred 40–72px)

### Borders & Radius
- Border color: `#e9ecef` light, `hsl(220,12%,23%)` dark
- Radii: small (3px) for code/badges, medium (6px) for inputs/buttons, large (10px) for cards
- No pill-shaped buttons; prefer `6px` radius on CTAs

### Shadows
- Subtle: `0 1px 4px rgba(0,0,0,0.09)` for cards and inputs
- No heavy drop shadows; the aesthetic is flat-to-barely-raised

### Animation
- Minimal — transitions on hover states only
- Duration: 150ms (fast) for hover color changes, 220ms (base) for reveals
- Easing: `ease` — no bounces, no springs in UI elements
- No page-entry animations in docs

### Cards
- Background: white with `1px solid #e9ecef` border + `shadow-sm`
- Radius: `10px`
- Padding: `24px`
- Feature cards on homepage: icon + title + detail text; no colored borders

### Hover / Active States
- Links: color shift only (no underline → underline on hover)
- Buttons: slightly darker background shade on hover; same border
- No opacity tricks; no scale transforms on UI elements

### Iconography
See ICONOGRAPHY section below.

### Imagery
- Geometry plots generated by **CairoMakie** / **Makie.jl** in SVG format (inline)
- Benchmark charts as PNG/SVG images
- No photography; no illustrations beyond diagrams
- Logo mark is the three-circle Julia logo adapted with JuliaGeo geographic map overlay

### Code Blocks
- Dark background (`hsl(220,20%,9%)`), light text
- `Space Mono` at 14px, relaxed line height
- Inline code: light gray background, red text (`#CB3C33`)

### Corner Radii Summary
| Element     | Radius |
|-------------|--------|
| Badges/tags | 3px    |
| Buttons     | 5px    |
| Inputs      | 5px    |
| Cards       | 8px    |
| Modals      | 12px   |

Intentionally tight — oversized rounding feels dated; small radii read as technical and current.

---

## ICONOGRAPHY

### Approach
JuliaGeo does not ship a custom icon system. Icons are used sparingly — primarily on the documentation homepage feature cards and in navigation.

### Usage in Docs
- Feature icons on the homepage are **inline `<img>` tags** pointing to external CDN-hosted PNGs/SVGs:
  - Julia dots logo (JuliaLang GitHub CDN)
  - Literate.jl logo (fredrikekre.github.io)
  - JuliaGeo SVG logo (juliageo.github.io CDN)
- No icon font is bundled
- No emoji in navigation or headings

### Recommendation for New Designs
Use **Lucide Icons** (CDN: `https://unpkg.com/lucide@latest`) — stroke-style, 1.5px weight, matches the clean technical aesthetic. Use at 16px or 20px in UI, 24px for feature callouts.

### Brand Assets (in `assets/`)
| File | Description |
|------|-------------|
| `assets/juliageo-logo.svg` | JuliaGeo org logo — three circles + geo map |
| `assets/geometryops-logo.png` | GeometryOps.jl package logo |
| `assets/favicon.ico` | Browser favicon |

---

## FILE INDEX

```
README.md                  ← This file
SKILL.md                   ← Agent skill descriptor
colors_and_type.css        ← All CSS design tokens + semantic styles
assets/
  juliageo-logo.svg        ← JuliaGeo organization logo
  geometryops-logo.png     ← GeometryOps package logo
  favicon.ico              ← Favicon
preview/
  colors-brand.html        ← Brand + Julia color swatches
  colors-semantic.html     ← Semantic color tokens
  colors-neutral.html      ← Gray scale
  type-scale.html          ← Heading + body type specimens
  type-mono.html           ← Monospace / code specimens
  spacing-tokens.html      ← Spacing scale
  radius-shadow.html       ← Border radius + shadow system
  components-buttons.html  ← Button variants
  components-badges.html   ← Badges + labels
  components-code.html     ← Code block + inline code
  components-cards.html    ← Card patterns
  components-inputs.html   ← Form inputs
  brand-logo.html          ← Logo usage
ui_kits/
  docs/
    README.md              ← Docs UI kit notes
    index.html             ← Docs site prototype
    Sidebar.jsx            ← Navigation sidebar
    Header.jsx             ← Top nav + search
    ContentPage.jsx        ← Doc page with prose + code
    HomePage.jsx           ← Hero + feature cards
```
