// Header.jsx — JuliaGeo Docs top navigation
// Load with <script type="text/babel" src="Header.jsx"></script>

const Header = ({ darkMode, onToggleDark, currentPage, onNav }) => {
  const s = headerStyles;
  return (
    <header style={s.header(darkMode)}>
      <div style={s.inner}>
        {/* Logo */}
        <button style={s.logo} onClick={() => onNav('home')}>
          <img src="../../assets/juliageo-logo.svg" alt="JuliaGeo" style={{height: 28, width: 'auto'}} />
          <span style={s.logoName(darkMode)}>GeometryOps<span style={s.logoJl(darkMode)}>.jl</span></span>
        </button>

        {/* Search */}
        <div style={s.searchWrap(darkMode)}>
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke={darkMode ? '#6c757d' : '#adb5bd'} strokeWidth="2">
            <circle cx="11" cy="11" r="8"/><line x1="21" y1="21" x2="16.65" y2="16.65"/>
          </svg>
          <span style={s.searchText(darkMode)}>Search docs…</span>
          <span style={s.searchKbd(darkMode)}>⌘K</span>
        </div>

        {/* Nav links */}
        <nav style={s.nav}>
          {['Introduction', 'API', 'Tutorials', 'Source'].map(label => (
            <button key={label}
              style={s.navLink(darkMode, currentPage === label.toLowerCase())}
              onClick={() => onNav(label.toLowerCase())}>
              {label}
            </button>
          ))}
          <div style={s.divider(darkMode)} />
          {/* GitHub */}
          <a href="https://github.com/JuliaGeo/GeometryOps.jl" target="_blank" rel="noreferrer" style={s.iconBtn(darkMode)} title="GitHub">
            <svg width="18" height="18" viewBox="0 0 24 24" fill={darkMode ? '#adb5bd' : '#6c757d'}>
              <path d="M12 2C6.477 2 2 6.484 2 12.017c0 4.425 2.865 8.18 6.839 9.504.5.092.682-.217.682-.483 0-.237-.008-.868-.013-1.703-2.782.605-3.369-1.343-3.369-1.343-.454-1.158-1.11-1.466-1.11-1.466-.908-.62.069-.608.069-.608 1.003.07 1.531 1.032 1.531 1.032.892 1.53 2.341 1.088 2.91.832.092-.647.35-1.088.636-1.338-2.22-.253-4.555-1.113-4.555-4.951 0-1.093.39-1.988 1.029-2.688-.103-.253-.446-1.272.098-2.65 0 0 .84-.27 2.75 1.026A9.564 9.564 0 0112 6.844c.85.004 1.705.115 2.504.337 1.909-1.296 2.747-1.027 2.747-1.027.546 1.379.202 2.398.1 2.651.64.7 1.028 1.595 1.028 2.688 0 3.848-2.339 4.695-4.566 4.943.359.309.678.92.678 1.855 0 1.338-.012 2.419-.012 2.747 0 .268.18.58.688.482A10.019 10.019 0 0022 12.017C22 6.484 17.522 2 12 2z"/>
            </svg>
          </a>
          {/* Dark mode toggle */}
          <button style={s.iconBtn(darkMode)} onClick={onToggleDark} title="Toggle dark mode">
            {darkMode
              ? <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="#adb5bd" strokeWidth="2"><circle cx="12" cy="12" r="5"/><line x1="12" y1="1" x2="12" y2="3"/><line x1="12" y1="21" x2="12" y2="23"/><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/><line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/><line x1="1" y1="12" x2="3" y2="12"/><line x1="21" y1="12" x2="23" y2="12"/><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/><line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/></svg>
              : <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="#6c757d" strokeWidth="2"><path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/></svg>
            }
          </button>
        </nav>
      </div>
    </header>
  );
};

const headerStyles = {
  header: (dark) => ({
    position: 'fixed', top: 0, left: 0, right: 0, zIndex: 100,
    height: 56,
    background: dark ? 'hsl(220,20%,9%)' : '#ffffff',
    borderBottom: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`,
    display: 'flex', alignItems: 'center',
  }),
  inner: {
    width: '100%', maxWidth: 1280, margin: '0 auto',
    padding: '0 24px', display: 'flex', alignItems: 'center', gap: 16,
  },
  logo: {
    display: 'flex', alignItems: 'center', gap: 10,
    background: 'none', border: 'none', cursor: 'pointer', padding: 0,
    textDecoration: 'none',
  },
  logoName: (dark) => ({
    fontSize: 16, fontWeight: 700, color: dark ? '#f8f9fa' : '#212529',
    fontFamily: 'Inter, sans-serif', letterSpacing: '-0.02em', whiteSpace: 'nowrap',
  }),
  logoJl: (dark) => ({ color: '#389826' }),
  searchWrap: (dark) => ({
    flex: 1, maxWidth: 280, height: 34,
    background: dark ? 'hsl(220,16%,13%)' : '#f8f9fa',
    border: `1.5px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`,
    borderRadius: 6, display: 'flex', alignItems: 'center', gap: 8,
    padding: '0 10px', cursor: 'text',
  }),
  searchText: (dark) => ({
    flex: 1, fontSize: 13, color: dark ? 'hsl(220,8%,56%)' : '#adb5bd',
    fontFamily: 'Inter, sans-serif',
  }),
  searchKbd: (dark) => ({
    fontSize: 11, color: dark ? 'hsl(220,8%,56%)' : '#ced4da',
    background: dark ? 'hsl(220,14%,17%)' : '#f1f3f5',
    border: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#dee2e6'}`,
    borderRadius: 4, padding: '1px 5px', fontFamily: 'Space Mono, monospace',
  }),
  nav: { display: 'flex', alignItems: 'center', gap: 4, marginLeft: 'auto' },
  navLink: (dark, active) => ({
    background: 'none', border: 'none', cursor: 'pointer',
    fontSize: 13, fontWeight: active ? 600 : 400,
    color: active ? '#389826' : (dark ? '#adb5bd' : '#495057'),
    padding: '4px 10px', borderRadius: 5,
    fontFamily: 'Inter, sans-serif',
    transition: 'color 150ms, background 150ms',
  }),
  divider: (dark) => ({
    width: 1, height: 20, background: dark ? 'hsl(220,12%,23%)' : '#e9ecef', margin: '0 4px',
  }),
  iconBtn: (dark) => ({
    display: 'flex', alignItems: 'center', justifyContent: 'center',
    width: 32, height: 32, borderRadius: 6, background: 'none', border: 'none',
    cursor: 'pointer', textDecoration: 'none',
    transition: 'background 150ms',
  }),
};

Object.assign(window, { Header });
