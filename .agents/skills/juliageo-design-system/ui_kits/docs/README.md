# JuliaGeo Docs UI Kit

A high-fidelity recreation of the GeometryOps.jl documentation site (VitePress-based).

## Screens
- `index.html` — Full interactive prototype (home → doc page → API reference)

## Components
- `Header.jsx` — Top nav bar with logo, search, dark mode toggle, GitHub link
- `Sidebar.jsx` — Collapsible navigation sidebar with sections
- `ContentPage.jsx` — Documentation page with prose, code blocks, callouts
- `HomePage.jsx` — Hero section + feature cards

## Design Notes
- Font: Inter (body/UI) + Space Mono (code)
- Brand green: #389826
- VitePress-style layout: fixed sidebar left, content center, TOC right
- Code blocks: dark bg `hsl(220,20%,9%)` with Space Mono 13px
- Inline code: light gray bg, red (#CB3C33) text
- Max content width: ~860px
