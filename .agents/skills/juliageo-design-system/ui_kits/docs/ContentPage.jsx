// ContentPage.jsx — Doc page with prose, code blocks, callouts
// Load with <script type="text/babel" src="ContentPage.jsx"></script>

const InlineCode = ({ children, dark }) => (
  <code style={{
    fontFamily: 'Space Mono, monospace', fontSize: '0.875em',
    background: dark ? 'hsl(220,14%,17%)' : '#f1f3f5',
    color: dark ? '#f472b6' : '#CB3C33',
    padding: '0.15em 0.4em', borderRadius: 3,
    border: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`,
  }}>{children}</code>
);

const CodeBlock = ({ dark, children, lang = 'julia' }) => (
  <pre style={{
    background: 'hsl(220,20%,9%)', color: '#d4d4d4',
    fontFamily: 'Space Mono, monospace', fontSize: 13, lineHeight: 1.65,
    padding: '20px 24px', borderRadius: 8, overflowX: 'auto',
    margin: '16px 0', border: '1px solid hsl(220,12%,23%)',
  }}>
    <code dangerouslySetInnerHTML={{__html: children}} />
  </pre>
);

const Callout = ({ type = 'tip', title, children, dark }) => {
  const colors = {
    tip:     { border: '#389826', bg: dark ? 'hsl(130,30%,10%)' : '#f0faea', title: '#2c7a1e' },
    warning: { border: '#b45309', bg: dark ? 'hsl(35,40%,12%)'  : '#fef3c7', title: '#78350f' },
    danger:  { border: '#CB3C33', bg: dark ? 'hsl(5,30%,10%)'   : '#fef2f2', title: '#9a2e27' },
  };
  const c = colors[type] || colors.tip;
  return (
    <div style={{
      borderLeft: `3px solid ${c.border}`, background: c.bg,
      padding: '12px 16px', borderRadius: '0 6px 6px 0', margin: '16px 0',
    }}>
      <div style={{ fontSize: 13, fontWeight: 600, color: c.title, marginBottom: 4, fontFamily: 'Inter, sans-serif' }}>{title || type.toUpperCase()}</div>
      <div style={{ fontSize: 13, color: dark ? '#adb5bd' : '#495057', lineHeight: 1.6, fontFamily: 'Inter, sans-serif' }}>{children}</div>
    </div>
  );
};

const Badge = ({ children, color = 'green' }) => {
  const colors = {
    green:  { bg: '#dcf5d7', fg: '#2c7a1e' },
    blue:   { bg: '#d4dcf9', fg: '#2d47a8' },
    purple: { bg: '#eadaf2', fg: '#6b3a85' },
  };
  const c = colors[color] || colors.green;
  return (
    <span style={{
      background: c.bg, color: c.fg, fontSize: 11, fontWeight: 600,
      padding: '2px 8px', borderRadius: 3, fontFamily: 'Inter, sans-serif',
      display: 'inline-block', marginRight: 6,
    }}>{children}</span>
  );
};

const ContentPage = ({ dark }) => {
  const s = contentStyles;
  return (
    <article style={s.article}>
      {/* Breadcrumb */}
      <div style={s.breadcrumb(dark)}>Source Code <span style={{color:'#adb5bd'}}>/</span> Methods <span style={{color:'#adb5bd'}}>/</span> Clipping</div>

      {/* Title */}
      <h1 style={s.h1(dark)}>Polygon Clipping</h1>
      <div style={{marginBottom: 20}}>
        <Badge>Source Code</Badge><Badge color="blue">Pure Julia</Badge>
      </div>

      <p style={s.p(dark)}>
        GeometryOps provides polygon clipping via three operations: <InlineCode dark={dark}>intersection</InlineCode>, <InlineCode dark={dark}>difference</InlineCode>, and <InlineCode dark={dark}>union</InlineCode>. All methods accept any <InlineCode dark={dark}>GeoInterface.jl</InlineCode>-compatible geometry type and return the same.
      </p>

      <Callout type="warning" title="Warning" dark={dark}>
        This package is still under heavy development! Results may change between minor versions.
      </Callout>

      <h2 style={s.h2(dark)}>Usage</h2>
      <p style={s.p(dark)}>The simplest way to compute an intersection is with the <InlineCode dark={dark}>GO.intersection</InlineCode> function. Specify a <InlineCode dark={dark}>target</InlineCode> trait to control what geometry type is returned.</p>

      <CodeBlock dark={dark}>{`<span style="color:#91dd33">using</span> GeometryOps <span style="color:#91dd33">as</span> GO
<span style="color:#91dd33">using</span> GeoInterface <span style="color:#91dd33">as</span> GI

<span style="color:#6c8a6c; font-style:italic"># Two overlapping polygons</span>
poly1 = GI.Polygon([[(<span style="color:#f9a8d4">0</span>, <span style="color:#f9a8d4">0</span>), (<span style="color:#f9a8d4">2</span>, <span style="color:#f9a8d4">0</span>), (<span style="color:#f9a8d4">2</span>, <span style="color:#f9a8d4">2</span>), (<span style="color:#f9a8d4">0</span>, <span style="color:#f9a8d4">2</span>), (<span style="color:#f9a8d4">0</span>, <span style="color:#f9a8d4">0</span>)]])
poly2 = GI.Polygon([[(<span style="color:#f9a8d4">1</span>, <span style="color:#f9a8d4">0</span>), (<span style="color:#f9a8d4">3</span>, <span style="color:#f9a8d4">0</span>), (<span style="color:#f9a8d4">3</span>, <span style="color:#f9a8d4">2</span>), (<span style="color:#f9a8d4">1</span>, <span style="color:#f9a8d4">2</span>), (<span style="color:#f9a8d4">1</span>, <span style="color:#f9a8d4">0</span>)]])

result = GO.<span style="color:#7dd3fc">intersection</span>(poly1, poly2; target = GI.PolygonTrait())
GO.<span style="color:#7dd3fc">area</span>(result)  <span style="color:#6c8a6c; font-style:italic"># → 2.0</span>`}</CodeBlock>

      <Callout type="tip" title="Tip" dark={dark}>
        Use <InlineCode dark={dark}>GO.apply</InlineCode> with a target trait to apply operations over large nested geometry collections efficiently.
      </Callout>

      <h2 style={s.h2(dark)}>Algorithm</h2>
      <p style={s.p(dark)}>
        Clipping is implemented using the Sutherland–Hodgman algorithm for convex polygons, and the Greiner–Hormann algorithm for the general case. The implementation is in pure Julia with no external C dependencies.
      </p>

      <h3 style={s.h3(dark)}>Available functions</h3>
      <div style={s.apiTable(dark)}>
        {[
          ['intersection(a, b; target)', 'Returns the intersection of two geometries'],
          ['difference(a, b; target)',   'Returns a minus b'],
          ['union(a, b; target)',        'Returns the union of two geometries'],
        ].map(([fn, desc]) => (
          <div key={fn} style={s.apiRow(dark)}>
            <div style={s.apiFn(dark)}><InlineCode dark={dark}>{fn}</InlineCode></div>
            <div style={s.apiDesc(dark)}>{desc}</div>
          </div>
        ))}
      </div>
    </article>
  );
};

const contentStyles = {
  article: { maxWidth: 860, padding: '40px 48px 80px', fontFamily: 'Inter, sans-serif' },
  breadcrumb: (dark) => ({ fontSize: 13, color: dark ? 'hsl(220,8%,56%)' : '#adb5bd', marginBottom: 16, display: 'flex', gap: 8, alignItems: 'center' }),
  h1: (dark) => ({ fontSize: 36, fontWeight: 700, color: dark ? '#f8f9fa' : '#212529', letterSpacing: '-0.025em', marginBottom: 12, lineHeight: 1.15 }),
  h2: (dark) => ({ fontSize: 24, fontWeight: 600, color: dark ? '#f8f9fa' : '#212529', letterSpacing: '-0.02em', marginTop: 36, marginBottom: 12, lineHeight: 1.25, borderTop: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`, paddingTop: 28 }),
  h3: (dark) => ({ fontSize: 18, fontWeight: 600, color: dark ? '#e9ecef' : '#343a40', marginTop: 24, marginBottom: 10 }),
  p:  (dark) => ({ fontSize: 15, color: dark ? '#adb5bd' : '#495057', lineHeight: 1.7, marginBottom: 16 }),
  apiTable: (dark) => ({ border: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`, borderRadius: 8, overflow: 'hidden' }),
  apiRow: (dark) => ({ display: 'flex', gap: 16, padding: '10px 16px', borderBottom: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`, alignItems: 'baseline', background: 'transparent' }),
  apiFn:   (dark) => ({ minWidth: 260, flexShrink: 0 }),
  apiDesc: (dark) => ({ fontSize: 13, color: dark ? '#6c757d' : '#6c757d', lineHeight: 1.5 }),
};

Object.assign(window, { ContentPage, InlineCode, CodeBlock, Callout, Badge });
