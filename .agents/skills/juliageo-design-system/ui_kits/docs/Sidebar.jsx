// Sidebar.jsx — JuliaGeo Docs navigation sidebar
// Load with <script type="text/babel" src="Sidebar.jsx"></script>

const NAV_SECTIONS = [
  {
    title: null,
    items: [{ label: 'Introduction', id: 'intro' }, { label: 'API Reference', id: 'api' }],
  },
  {
    title: 'Tutorials',
    items: [
      { label: 'Creating Geometry', id: 'tut-creating' },
      { label: 'Spatial Joins', id: 'tut-joins' },
    ],
  },
  {
    title: 'Explanations',
    items: [
      { label: 'Paradigms', id: 'exp-paradigms' },
      { label: 'Manifolds', id: 'exp-manifolds' },
      { label: 'Performance', id: 'exp-perf' },
      { label: 'Peculiarities', id: 'exp-peculiar' },
    ],
  },
  {
    title: 'GIS Terminology',
    items: [
      { label: 'CRS', id: 'gis-crs' },
      { label: 'Winding Order', id: 'gis-winding' },
    ],
  },
  {
    title: 'Source Code',
    items: [
      { label: 'Methods / Clipping', id: 'src-clipping' },
      { label: 'Methods / Simplify', id: 'src-simplify' },
      { label: 'Methods / Distance', id: 'src-distance' },
      { label: 'GeometryOpsCore', id: 'src-core' },
    ],
  },
];

const Sidebar = ({ darkMode, currentPage, onNav }) => {
  const s = sidebarStyles;
  return (
    <aside style={s.sidebar(darkMode)}>
      <nav style={s.nav}>
        {NAV_SECTIONS.map((section, i) => (
          <div key={i} style={s.section}>
            {section.title && <div style={s.sectionTitle(darkMode)}>{section.title}</div>}
            {section.items.map(item => (
              <button key={item.id}
                style={s.item(darkMode, currentPage === item.id)}
                onClick={() => onNav(item.id)}>
                {item.label}
              </button>
            ))}
          </div>
        ))}
      </nav>
    </aside>
  );
};

const sidebarStyles = {
  sidebar: (dark) => ({
    position: 'fixed', top: 56, left: 0, bottom: 0, width: 260,
    background: dark ? 'hsl(220,20%,9%)' : '#ffffff',
    borderRight: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`,
    overflowY: 'auto', padding: '20px 0 40px',
  }),
  nav: { display: 'flex', flexDirection: 'column' },
  section: { marginBottom: 4 },
  sectionTitle: (dark) => ({
    fontSize: 11, fontWeight: 700, textTransform: 'uppercase', letterSpacing: '0.1em',
    color: dark ? 'hsl(220,8%,40%)' : '#adb5bd',
    padding: '14px 20px 4px', fontFamily: 'Inter, sans-serif',
  }),
  item: (dark, active) => ({
    display: 'block', width: '100%', textAlign: 'left',
    padding: '5px 20px 5px 24px', border: 'none', cursor: 'pointer',
    fontSize: 13, fontFamily: 'Inter, sans-serif',
    fontWeight: active ? 500 : 400,
    color: active ? '#389826' : (dark ? '#adb5bd' : '#495057'),
    background: active ? (dark ? 'hsl(220,14%,17%)' : '#f0faea') : 'none',
    borderLeft: active ? '2px solid #389826' : '2px solid transparent',
    transition: 'all 120ms ease',
  }),
};

Object.assign(window, { Sidebar, NAV_SECTIONS });
