// HomePage.jsx — JuliaGeo docs hero + feature cards
// Load with <script type="text/babel" src="HomePage.jsx"></script>

const FEATURES = [
  {
    icon: <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="#389826" strokeWidth="1.75"><polygon points="12 2 22 8.5 22 15.5 12 22 2 15.5 2 8.5"/><line x1="12" y1="2" x2="12" y2="22"/><line x1="2" y1="8.5" x2="22" y2="8.5"/></svg>,
    title: 'Pure Julia code',
    detail: 'Fast, understandable, extensible functions written entirely in Julia — no C bindings required.',
    link: '/introduction',
  },
  {
    icon: <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="#389826" strokeWidth="1.75"><path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"/><polyline points="14 2 14 8 20 8"/><line x1="16" y1="13" x2="8" y2="13"/><line x1="16" y1="17" x2="8" y2="17"/></svg>,
    title: 'Literate programming',
    detail: 'Documented source code with embedded examples via Literate.jl — code and explanation together.',
    link: '/source/methods/clipping/cut',
  },
  {
    icon: <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="#389826" strokeWidth="1.75"><circle cx="12" cy="12" r="10"/><line x1="2" y1="12" x2="22" y2="12"/><path d="M12 2a15.3 15.3 0 0 1 4 10 15.3 15.3 0 0 1-4 10 15.3 15.3 0 0 1-4-10 15.3 15.3 0 0 1 4-10z"/></svg>,
    title: 'GeoInterface integration',
    detail: 'Use any GeoInterface.jl-compatible geometry — ArchGDAL, Shapefile, GeoJSON, and more.',
    link: 'https://juliageo.org/GeoInterface.jl/stable',
  },
  {
    icon: <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="#389826" strokeWidth="1.75"><polyline points="22 12 18 12 15 21 9 3 6 12 2 12"/></svg>,
    title: 'Blazing fast',
    detail: 'Outperforms Python and R equivalents on most geometry benchmarks. Pure Julia means full JIT optimization.',
    link: '#benchmarks',
  },
];

const METHODS = [
  'equals', 'extent', 'distance', 'crosses', 'contains', 'intersects',
  'intersection', 'difference', 'union', 'simplify', 'centroid',
  'signed_area', 'segmentize', 'polygonize', 'barycentric_coordinates',
];

const HomePage = ({ dark, onNav }) => {
  const s = homeStyles;
  return (
    <div style={s.page}>
      {/* Hero */}
      <section style={s.hero(dark)}>
        <div style={s.heroInner}>
          <img src="../../assets/geometryops-logo.png" alt="GeometryOps" style={s.heroLogo} />
          <div>
            <h1 style={s.heroTitle(dark)}>GeometryOps<span style={{color:'#389826'}}>.jl</span></h1>
            <p style={s.heroTagline(dark)}>Blazing fast geometry operations in pure Julia</p>
            <div style={s.heroActions}>
              <button style={s.btnPrimary} onClick={() => onNav('intro')}>Introduction</button>
              <button style={s.btnSecondary(dark)} onClick={() => onNav('api')}>API Reference</button>
              <a href="https://github.com/JuliaGeo/GeometryOps.jl" target="_blank" rel="noreferrer" style={s.btnGhost(dark)}>View on GitHub</a>
            </div>
          </div>
        </div>
      </section>

      {/* Content */}
      <div style={s.content}>
        {/* Features */}
        <section style={s.section}>
          <div style={s.featureGrid}>
            {FEATURES.map(f => (
              <div key={f.title} style={s.featureCard(dark)}>
                <div style={s.featureIcon(dark)}>{f.icon}</div>
                <div style={s.featureTitle(dark)}>{f.title}</div>
                <div style={s.featureDetail(dark)}>{f.detail}</div>
              </div>
            ))}
          </div>
        </section>

        {/* What is it */}
        <section style={s.section}>
          <h2 style={s.h2(dark)}>What is GeometryOps.jl?</h2>
          <p style={s.p(dark)}>
            GeometryOps.jl is a package for geometric calculations on (primarily 2D) geometries. The driving idea is to unify all the disparate packages for geometric calculations in Julia, and make them <a href="#" style={s.link}>GeoInterface.jl</a>-compatible.
          </p>
          <p style={s.p(dark)}>
            Most use cases are driven by GIS and similar Earth data workflows. Methods are always general to any coordinate space.
          </p>
        </section>

        {/* Methods grid */}
        <section style={s.section}>
          <h2 style={s.h2(dark)}>Available methods</h2>
          <div style={s.methodGrid}>
            {METHODS.map(m => (
              <span key={m} style={s.methodTag(dark)}>{m}</span>
            ))}
          </div>
        </section>
      </div>
    </div>
  );
};

const homeStyles = {
  page: { fontFamily: 'Inter, sans-serif' },
  hero: (dark) => ({
    padding: '72px 48px 56px',
    background: dark ? 'hsl(220,20%,9%)' : '#ffffff',
    borderBottom: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`,
  }),
  heroInner: { maxWidth: 860, display: 'flex', gap: 48, alignItems: 'center' },
  heroLogo: { width: 120, height: 'auto', flexShrink: 0 },
  heroTitle: (dark) => ({ fontSize: 48, fontWeight: 700, color: dark ? '#f8f9fa' : '#212529', letterSpacing: '-0.03em', margin: '0 0 10px', lineHeight: 1.05 }),
  heroTagline: (dark) => ({ fontSize: 18, color: dark ? 'hsl(220,8%,56%)' : '#6c757d', margin: '0 0 28px', lineHeight: 1.4 }),
  heroActions: { display: 'flex', gap: 10, flexWrap: 'wrap' },
  btnPrimary: { height: 38, padding: '0 20px', background: '#389826', color: '#fff', border: '1.5px solid #389826', borderRadius: 6, fontSize: 14, fontWeight: 500, cursor: 'pointer', fontFamily: 'Inter, sans-serif' },
  btnSecondary: (dark) => ({ height: 38, padding: '0 20px', background: 'transparent', color: '#389826', border: '1.5px solid #389826', borderRadius: 6, fontSize: 14, fontWeight: 500, cursor: 'pointer', fontFamily: 'Inter, sans-serif' }),
  btnGhost: (dark) => ({ height: 38, padding: '0 20px', background: 'transparent', color: dark ? '#adb5bd' : '#495057', border: `1.5px solid ${dark ? 'hsl(220,12%,23%)' : '#dee2e6'}`, borderRadius: 6, fontSize: 14, fontWeight: 400, cursor: 'pointer', fontFamily: 'Inter, sans-serif', textDecoration: 'none', display: 'inline-flex', alignItems: 'center' }),
  content: { maxWidth: 860, padding: '48px 48px 80px' },
  section: { marginBottom: 52 },
  featureGrid: { display: 'grid', gridTemplateColumns: 'repeat(2, 1fr)', gap: 16 },
  featureCard: (dark) => ({ background: dark ? 'hsl(220,16%,13%)' : '#ffffff', border: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`, borderRadius: 10, padding: '20px 24px', boxShadow: '0 1px 4px rgba(0,0,0,0.06)' }),
  featureIcon: (dark) => ({ width: 44, height: 44, background: dark ? 'hsl(130,20%,12%)' : '#f0faea', borderRadius: 8, display: 'flex', alignItems: 'center', justifyContent: 'center', marginBottom: 12 }),
  featureTitle: (dark) => ({ fontSize: 15, fontWeight: 600, color: dark ? '#f8f9fa' : '#212529', marginBottom: 6 }),
  featureDetail: (dark) => ({ fontSize: 13, color: dark ? 'hsl(220,8%,56%)' : '#6c757d', lineHeight: 1.6 }),
  h2: (dark) => ({ fontSize: 28, fontWeight: 700, color: dark ? '#f8f9fa' : '#212529', letterSpacing: '-0.025em', marginBottom: 16 }),
  p:  (dark) => ({ fontSize: 15, color: dark ? '#adb5bd' : '#495057', lineHeight: 1.7, marginBottom: 14 }),
  link: { color: '#2c7a1e', textDecoration: 'none' },
  methodGrid: { display: 'flex', flexWrap: 'wrap', gap: 8 },
  methodTag: (dark) => ({ display: 'inline-flex', padding: '4px 12px', borderRadius: 6, fontFamily: 'Space Mono, monospace', fontSize: 12, background: dark ? 'hsl(220,14%,17%)' : '#f1f3f5', color: dark ? '#adb5bd' : '#495057', border: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}` }),
};

Object.assign(window, { HomePage });
