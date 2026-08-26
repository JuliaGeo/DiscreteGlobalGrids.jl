# Dagger for Global Grids

A ten-minute Slidev meeting primer titled "Dagger for Global Grids." It uses
the same JuliaGeo visual system as `../data_analysis_on_dggs` and retains that
deck's Julia figure environment for reuse.

```bash
npm install
npm run dev
npm run build
```

The existing WGLMakie exports live in `figures/html/` and are served through
`public/figures`. To rebuild the retained Julia figures, instantiate the copied
environment and run the relevant series script:

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project figures/01-what-is_dggs-series.jl
```

`STORYBOARD.md` is the narrative source. `slides.md` follows it one beat per
visual, with interactive exports under `figures/html/` and rendered camera
sequences under `figures/video/`.
