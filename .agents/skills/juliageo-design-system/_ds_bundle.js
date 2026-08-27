/* @ds-bundle: {"format":3,"namespace":"JuliaGeoDesignSystem_019dcb","components":[],"sourceHashes":{"animations.jsx":"2b1e0ed2c732","slides/deck-stage.js":"ad1c016a6256","ui_kits/docs/ContentPage.jsx":"0070e70fc2e8","ui_kits/docs/Header.jsx":"372d0d0c325a","ui_kits/docs/HomePage.jsx":"37d01ce49d16","ui_kits/docs/Sidebar.jsx":"0b9a94ebcc7b"},"inlinedExternals":[],"unexposedExports":[]} */

(() => {

const __ds_ns = (window.JuliaGeoDesignSystem_019dcb = window.JuliaGeoDesignSystem_019dcb || {});

const __ds_scope = {};

(__ds_ns.__errors = __ds_ns.__errors || []);

// animations.jsx
try { (() => {
// animations.jsx
// Reusable animation starter: Stage, Timeline, Sprite, easing helpers.
// Usage (in an HTML file that loads React + Babel):
//
//   <Stage width={1280} height={720} duration={10} background="#f6f4ef">
//     <MyScene />
//   </Stage>
//
// Inside <Stage>, any child can call useTime() to read the current
// playhead (seconds). Or wrap content in <Sprite start={1} end={4}>...</Sprite>
// to only render during that window -- children receive a `localTime` and
// `progress` via the useSprite() hook.
//
// ─────────────────────────────────────────────────────────────────────────────

// ── Easing functions (hand-rolled, Popmotion-style) ─────────────────────────
// All easings take t ∈ [0,1] and return eased t ∈ [0,1] (may overshoot for back/elastic).
const Easing = {
  linear: t => t,
  // Quad
  easeInQuad: t => t * t,
  easeOutQuad: t => t * (2 - t),
  easeInOutQuad: t => t < 0.5 ? 2 * t * t : -1 + (4 - 2 * t) * t,
  // Cubic
  easeInCubic: t => t * t * t,
  easeOutCubic: t => --t * t * t + 1,
  easeInOutCubic: t => t < 0.5 ? 4 * t * t * t : (t - 1) * (2 * t - 2) * (2 * t - 2) + 1,
  // Quart
  easeInQuart: t => t * t * t * t,
  easeOutQuart: t => 1 - --t * t * t * t,
  easeInOutQuart: t => t < 0.5 ? 8 * t * t * t * t : 1 - 8 * --t * t * t * t,
  // Expo
  easeInExpo: t => t === 0 ? 0 : Math.pow(2, 10 * (t - 1)),
  easeOutExpo: t => t === 1 ? 1 : 1 - Math.pow(2, -10 * t),
  easeInOutExpo: t => {
    if (t === 0) return 0;
    if (t === 1) return 1;
    if (t < 0.5) return 0.5 * Math.pow(2, 20 * t - 10);
    return 1 - 0.5 * Math.pow(2, -20 * t + 10);
  },
  // Sine
  easeInSine: t => 1 - Math.cos(t * Math.PI / 2),
  easeOutSine: t => Math.sin(t * Math.PI / 2),
  easeInOutSine: t => -(Math.cos(Math.PI * t) - 1) / 2,
  // Back (overshoot)
  easeOutBack: t => {
    const c1 = 1.70158,
      c3 = c1 + 1;
    return 1 + c3 * Math.pow(t - 1, 3) + c1 * Math.pow(t - 1, 2);
  },
  easeInBack: t => {
    const c1 = 1.70158,
      c3 = c1 + 1;
    return c3 * t * t * t - c1 * t * t;
  },
  easeInOutBack: t => {
    const c1 = 1.70158,
      c2 = c1 * 1.525;
    return t < 0.5 ? Math.pow(2 * t, 2) * ((c2 + 1) * 2 * t - c2) / 2 : (Math.pow(2 * t - 2, 2) * ((c2 + 1) * (t * 2 - 2) + c2) + 2) / 2;
  },
  // Elastic
  easeOutElastic: t => {
    const c4 = 2 * Math.PI / 3;
    if (t === 0) return 0;
    if (t === 1) return 1;
    return Math.pow(2, -10 * t) * Math.sin((t * 10 - 0.75) * c4) + 1;
  }
};

// ── Core interpolation helpers ──────────────────────────────────────────────

// Clamp a value to [min, max]
const clamp = (v, min, max) => Math.max(min, Math.min(max, v));

// interpolate([0, 0.5, 1], [0, 100, 50], ease?) -> fn(t)
// Popmotion-style: linearly maps t across input keyframes to output values,
// with optional easing per segment (single fn or array of fns).
function interpolate(input, output, ease = Easing.linear) {
  return t => {
    if (t <= input[0]) return output[0];
    if (t >= input[input.length - 1]) return output[output.length - 1];
    for (let i = 0; i < input.length - 1; i++) {
      if (t >= input[i] && t <= input[i + 1]) {
        const span = input[i + 1] - input[i];
        const local = span === 0 ? 0 : (t - input[i]) / span;
        const easeFn = Array.isArray(ease) ? ease[i] || Easing.linear : ease;
        const eased = easeFn(local);
        return output[i] + (output[i + 1] - output[i]) * eased;
      }
    }
    return output[output.length - 1];
  };
}

// animate({from, to, start, end, ease})(t) — simpler single-segment tween.
// Returns `from` before `start`, `to` after `end`.
function animate({
  from = 0,
  to = 1,
  start = 0,
  end = 1,
  ease = Easing.easeInOutCubic
}) {
  return t => {
    if (t <= start) return from;
    if (t >= end) return to;
    const local = (t - start) / (end - start);
    return from + (to - from) * ease(local);
  };
}

// ── Timeline context ────────────────────────────────────────────────────────

const TimelineContext = React.createContext({
  time: 0,
  duration: 10,
  playing: false
});
const useTime = () => React.useContext(TimelineContext).time;
const useTimeline = () => React.useContext(TimelineContext);

// ── Sprite ──────────────────────────────────────────────────────────────────
// Renders children only when the playhead is inside [start, end]. Provides
// a sub-context with `localTime` (seconds since start) and `progress` (0..1).
//
//   <Sprite start={2} end={5}>
//     {({ localTime, progress }) => <Thing x={progress * 100} />}
//   </Sprite>
//
// Or as a plain wrapper — children can call useSprite() themselves.

const SpriteContext = React.createContext({
  localTime: 0,
  progress: 0,
  duration: 0
});
const useSprite = () => React.useContext(SpriteContext);
function Sprite({
  start = 0,
  end = Infinity,
  children,
  keepMounted = false
}) {
  const {
    time
  } = useTimeline();
  const visible = time >= start && time <= end;
  if (!visible && !keepMounted) return null;
  const duration = end - start;
  const localTime = Math.max(0, time - start);
  const progress = duration > 0 && isFinite(duration) ? clamp(localTime / duration, 0, 1) : 0;
  const value = {
    localTime,
    progress,
    duration,
    visible
  };
  return /*#__PURE__*/React.createElement(SpriteContext.Provider, {
    value: value
  }, typeof children === 'function' ? children(value) : children);
}

// ── Sample sprite components ────────────────────────────────────────────────

// TextSprite: fades/slides text in on entry, holds, then fades out on exit.
// Props: text, x, y, size, color, font, entryDur, exitDur, align
function TextSprite({
  text,
  x = 0,
  y = 0,
  size = 48,
  color = '#111',
  font = 'Inter, system-ui, sans-serif',
  weight = 600,
  entryDur = 0.45,
  exitDur = 0.35,
  entryEase = Easing.easeOutBack,
  exitEase = Easing.easeInCubic,
  align = 'left',
  letterSpacing = '-0.01em'
}) {
  const {
    localTime,
    duration
  } = useSprite();
  const exitStart = Math.max(0, duration - exitDur);
  let opacity = 1;
  let ty = 0;
  if (localTime < entryDur) {
    const t = entryEase(clamp(localTime / entryDur, 0, 1));
    opacity = t;
    ty = (1 - t) * 16;
  } else if (localTime > exitStart) {
    const t = exitEase(clamp((localTime - exitStart) / exitDur, 0, 1));
    opacity = 1 - t;
    ty = -t * 8;
  }
  const translateX = align === 'center' ? '-50%' : align === 'right' ? '-100%' : '0';
  return /*#__PURE__*/React.createElement("div", {
    style: {
      position: 'absolute',
      left: x,
      top: y,
      transform: `translate(${translateX}, ${ty}px)`,
      opacity,
      fontFamily: font,
      fontSize: size,
      fontWeight: weight,
      color,
      letterSpacing,
      whiteSpace: 'pre',
      lineHeight: 1.1,
      willChange: 'transform, opacity'
    }
  }, text);
}

// ImageSprite: scales + fades in; optional Ken Burns drift during hold.
function ImageSprite({
  src,
  x = 0,
  y = 0,
  width = 400,
  height = 300,
  entryDur = 0.6,
  exitDur = 0.4,
  kenBurns = false,
  kenBurnsScale = 1.08,
  radius = 12,
  fit = 'cover',
  placeholder = null // {label: string} for striped placeholder
}) {
  const {
    localTime,
    duration
  } = useSprite();
  const exitStart = Math.max(0, duration - exitDur);
  let opacity = 1;
  let scale = 1;
  if (localTime < entryDur) {
    const t = Easing.easeOutCubic(clamp(localTime / entryDur, 0, 1));
    opacity = t;
    scale = 0.96 + 0.04 * t;
  } else if (localTime > exitStart) {
    const t = Easing.easeInCubic(clamp((localTime - exitStart) / exitDur, 0, 1));
    opacity = 1 - t;
    scale = (kenBurns ? kenBurnsScale : 1) + 0.02 * t;
  } else if (kenBurns) {
    const holdSpan = exitStart - entryDur;
    const holdT = holdSpan > 0 ? (localTime - entryDur) / holdSpan : 0;
    scale = 1 + (kenBurnsScale - 1) * holdT;
  }
  const content = placeholder ? /*#__PURE__*/React.createElement("div", {
    style: {
      width: '100%',
      height: '100%',
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      background: 'repeating-linear-gradient(135deg, #e9e6df 0 10px, #dcd8cf 10px 20px)',
      color: '#6b6458',
      fontFamily: 'JetBrains Mono, ui-monospace, monospace',
      fontSize: 13,
      letterSpacing: '0.04em',
      textTransform: 'uppercase'
    }
  }, placeholder.label || 'image') : /*#__PURE__*/React.createElement("img", {
    src: src,
    alt: "",
    style: {
      width: '100%',
      height: '100%',
      objectFit: fit,
      display: 'block'
    }
  });
  return /*#__PURE__*/React.createElement("div", {
    style: {
      position: 'absolute',
      left: x,
      top: y,
      width,
      height,
      opacity,
      transform: `scale(${scale})`,
      transformOrigin: 'center',
      borderRadius: radius,
      overflow: 'hidden',
      willChange: 'transform, opacity'
    }
  }, content);
}

// RectSprite: simple rectangle that animates position/size/color via props.
// Useful demo primitive — takes a `render` fn for per-frame customization.
function RectSprite({
  x = 0,
  y = 0,
  width = 100,
  height = 100,
  color = '#111',
  radius = 8,
  entryDur = 0.4,
  exitDur = 0.3,
  render // optional: (ctx) => style overrides
}) {
  const spriteCtx = useSprite();
  const {
    localTime,
    duration
  } = spriteCtx;
  const exitStart = Math.max(0, duration - exitDur);
  let opacity = 1;
  let scale = 1;
  if (localTime < entryDur) {
    const t = Easing.easeOutBack(clamp(localTime / entryDur, 0, 1));
    opacity = clamp(localTime / entryDur, 0, 1);
    scale = 0.4 + 0.6 * t;
  } else if (localTime > exitStart) {
    const t = Easing.easeInQuad(clamp((localTime - exitStart) / exitDur, 0, 1));
    opacity = 1 - t;
    scale = 1 - 0.15 * t;
  }
  const overrides = render ? render(spriteCtx) : {};
  return /*#__PURE__*/React.createElement("div", {
    style: {
      position: 'absolute',
      left: x,
      top: y,
      width,
      height,
      background: color,
      borderRadius: radius,
      opacity,
      transform: `scale(${scale})`,
      transformOrigin: 'center',
      willChange: 'transform, opacity',
      ...overrides
    }
  });
}
function Stage({
  width = 1280,
  height = 720,
  duration = 10,
  background = '#f6f4ef',
  fps = 60,
  loop = true,
  autoplay = true,
  persistKey = 'animstage',
  children
}) {
  const [time, setTime] = React.useState(() => {
    try {
      const v = parseFloat(localStorage.getItem(persistKey + ':t') || '0');
      return isFinite(v) ? clamp(v, 0, duration) : 0;
    } catch {
      return 0;
    }
  });
  const [playing, setPlaying] = React.useState(autoplay);
  const [hoverTime, setHoverTime] = React.useState(null);
  const [scale, setScale] = React.useState(1);
  const stageRef = React.useRef(null);
  const canvasRef = React.useRef(null);
  const rafRef = React.useRef(null);
  const lastTsRef = React.useRef(null);

  // Persist playhead
  React.useEffect(() => {
    try {
      localStorage.setItem(persistKey + ':t', String(time));
    } catch {}
  }, [time, persistKey]);

  // Auto-scale to fit viewport
  React.useEffect(() => {
    if (!stageRef.current) return;
    const el = stageRef.current;
    const measure = () => {
      const barH = 44; // playback bar height
      const s = Math.min(el.clientWidth / width, (el.clientHeight - barH) / height);
      setScale(Math.max(0.05, s));
    };
    measure();
    const ro = new ResizeObserver(measure);
    ro.observe(el);
    window.addEventListener('resize', measure);
    return () => {
      ro.disconnect();
      window.removeEventListener('resize', measure);
    };
  }, [width, height]);

  // Animation loop
  React.useEffect(() => {
    if (!playing) {
      lastTsRef.current = null;
      return;
    }
    const step = ts => {
      if (lastTsRef.current == null) lastTsRef.current = ts;
      const dt = (ts - lastTsRef.current) / 1000;
      lastTsRef.current = ts;
      setTime(t => {
        let next = t + dt;
        if (next >= duration) {
          if (loop) next = next % duration;else {
            next = duration;
            setPlaying(false);
          }
        }
        return next;
      });
      rafRef.current = requestAnimationFrame(step);
    };
    rafRef.current = requestAnimationFrame(step);
    return () => {
      if (rafRef.current) cancelAnimationFrame(rafRef.current);
      lastTsRef.current = null;
    };
  }, [playing, duration, loop]);

  // Keyboard: space = play/pause, ← → = seek
  React.useEffect(() => {
    const onKey = e => {
      if (e.target && (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA')) return;
      if (e.code === 'Space') {
        e.preventDefault();
        setPlaying(p => !p);
      } else if (e.code === 'ArrowLeft') {
        setTime(t => clamp(t - (e.shiftKey ? 1 : 0.1), 0, duration));
      } else if (e.code === 'ArrowRight') {
        setTime(t => clamp(t + (e.shiftKey ? 1 : 0.1), 0, duration));
      } else if (e.key === '0' || e.code === 'Home') {
        setTime(0);
      }
    };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, [duration]);
  const displayTime = hoverTime != null ? hoverTime : time;
  const ctxValue = React.useMemo(() => ({
    time: displayTime,
    duration,
    playing,
    setTime,
    setPlaying
  }), [displayTime, duration, playing]);
  return /*#__PURE__*/React.createElement("div", {
    ref: stageRef,
    style: {
      position: 'absolute',
      inset: 0,
      display: 'flex',
      flexDirection: 'column',
      alignItems: 'center',
      background: '#0a0a0a',
      fontFamily: 'Inter, system-ui, sans-serif'
    }
  }, /*#__PURE__*/React.createElement("div", {
    style: {
      flex: 1,
      width: '100%',
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      overflow: 'hidden',
      minHeight: 0
    }
  }, /*#__PURE__*/React.createElement("div", {
    ref: canvasRef,
    style: {
      width,
      height,
      background,
      position: 'relative',
      transform: `scale(${scale})`,
      transformOrigin: 'center',
      flexShrink: 0,
      boxShadow: '0 20px 60px rgba(0,0,0,0.4)',
      overflow: 'hidden'
    }
  }, /*#__PURE__*/React.createElement(TimelineContext.Provider, {
    value: ctxValue
  }, children))), /*#__PURE__*/React.createElement(PlaybackBar, {
    time: displayTime,
    actualTime: time,
    duration: duration,
    playing: playing,
    onPlayPause: () => setPlaying(p => !p),
    onReset: () => {
      setTime(0);
    },
    onSeek: t => setTime(t),
    onHover: t => setHoverTime(t)
  }));
}

// ── Playback bar ────────────────────────────────────────────────────────────
// Play/pause, return-to-begin, scrub track, time display.
// Uses fixed-width time fields so layout doesn't thrash.

function PlaybackBar({
  time,
  duration,
  playing,
  onPlayPause,
  onReset,
  onSeek,
  onHover
}) {
  const trackRef = React.useRef(null);
  const [dragging, setDragging] = React.useState(false);
  const timeFromEvent = React.useCallback(e => {
    const rect = trackRef.current.getBoundingClientRect();
    const x = clamp((e.clientX - rect.left) / rect.width, 0, 1);
    return x * duration;
  }, [duration]);
  const onTrackMove = e => {
    if (!trackRef.current) return;
    const t = timeFromEvent(e);
    if (dragging) {
      onSeek(t);
    } else {
      onHover(t);
    }
  };
  const onTrackLeave = () => {
    if (!dragging) onHover(null);
  };
  const onTrackDown = e => {
    setDragging(true);
    const t = timeFromEvent(e);
    onSeek(t);
    onHover(null);
  };
  React.useEffect(() => {
    if (!dragging) return;
    const onUp = () => setDragging(false);
    const onMove = e => {
      if (!trackRef.current) return;
      const t = timeFromEvent(e);
      onSeek(t);
    };
    window.addEventListener('mouseup', onUp);
    window.addEventListener('mousemove', onMove);
    return () => {
      window.removeEventListener('mouseup', onUp);
      window.removeEventListener('mousemove', onMove);
    };
  }, [dragging, timeFromEvent, onSeek]);
  const pct = duration > 0 ? time / duration * 100 : 0;
  const fmt = t => {
    const total = Math.max(0, t);
    const m = Math.floor(total / 60);
    const s = Math.floor(total % 60);
    const cs = Math.floor(total * 100 % 100);
    return `${String(m).padStart(1, '0')}:${String(s).padStart(2, '0')}.${String(cs).padStart(2, '0')}`;
  };
  const mono = 'JetBrains Mono, ui-monospace, SFMono-Regular, monospace';
  return /*#__PURE__*/React.createElement("div", {
    style: {
      display: 'flex',
      alignItems: 'center',
      gap: 12,
      padding: '8px 16px',
      background: 'rgba(20,20,20,0.92)',
      borderTop: '1px solid rgba(255,255,255,0.08)',
      width: '100%',
      maxWidth: 680,
      alignSelf: 'center',
      borderRadius: 8,
      color: '#f6f4ef',
      fontFamily: 'Inter, system-ui, sans-serif',
      userSelect: 'none',
      flexShrink: 0
    }
  }, /*#__PURE__*/React.createElement(IconButton, {
    onClick: onReset,
    title: "Return to start (0)"
  }, /*#__PURE__*/React.createElement("svg", {
    width: "14",
    height: "14",
    viewBox: "0 0 14 14",
    fill: "none"
  }, /*#__PURE__*/React.createElement("path", {
    d: "M3 2v10M12 2L5 7l7 5V2z",
    stroke: "currentColor",
    strokeWidth: "1.5",
    strokeLinejoin: "round",
    strokeLinecap: "round"
  }))), /*#__PURE__*/React.createElement(IconButton, {
    onClick: onPlayPause,
    title: "Play/pause (space)"
  }, playing ? /*#__PURE__*/React.createElement("svg", {
    width: "14",
    height: "14",
    viewBox: "0 0 14 14",
    fill: "none"
  }, /*#__PURE__*/React.createElement("rect", {
    x: "3",
    y: "2",
    width: "3",
    height: "10",
    fill: "currentColor"
  }), /*#__PURE__*/React.createElement("rect", {
    x: "8",
    y: "2",
    width: "3",
    height: "10",
    fill: "currentColor"
  })) : /*#__PURE__*/React.createElement("svg", {
    width: "14",
    height: "14",
    viewBox: "0 0 14 14",
    fill: "none"
  }, /*#__PURE__*/React.createElement("path", {
    d: "M3 2l9 5-9 5V2z",
    fill: "currentColor"
  }))), /*#__PURE__*/React.createElement("div", {
    style: {
      fontFamily: mono,
      fontSize: 12,
      fontVariantNumeric: 'tabular-nums',
      width: 64,
      textAlign: 'right',
      color: '#f6f4ef'
    }
  }, fmt(time)), /*#__PURE__*/React.createElement("div", {
    ref: trackRef,
    onMouseMove: onTrackMove,
    onMouseLeave: onTrackLeave,
    onMouseDown: onTrackDown,
    style: {
      flex: 1,
      height: 22,
      position: 'relative',
      cursor: 'pointer',
      display: 'flex',
      alignItems: 'center'
    }
  }, /*#__PURE__*/React.createElement("div", {
    style: {
      position: 'absolute',
      left: 0,
      right: 0,
      height: 4,
      background: 'rgba(255,255,255,0.12)',
      borderRadius: 2
    }
  }), /*#__PURE__*/React.createElement("div", {
    style: {
      position: 'absolute',
      left: 0,
      width: `${pct}%`,
      height: 4,
      background: 'oklch(72% 0.12 250)',
      borderRadius: 2
    }
  }), /*#__PURE__*/React.createElement("div", {
    style: {
      position: 'absolute',
      left: `${pct}%`,
      top: '50%',
      width: 12,
      height: 12,
      marginLeft: -6,
      marginTop: -6,
      background: '#fff',
      borderRadius: 6,
      boxShadow: '0 2px 4px rgba(0,0,0,0.4)'
    }
  })), /*#__PURE__*/React.createElement("div", {
    style: {
      fontFamily: mono,
      fontSize: 12,
      fontVariantNumeric: 'tabular-nums',
      width: 64,
      textAlign: 'left',
      color: 'rgba(246,244,239,0.55)'
    }
  }, fmt(duration)));
}
function IconButton({
  children,
  onClick,
  title
}) {
  const [hover, setHover] = React.useState(false);
  return /*#__PURE__*/React.createElement("button", {
    onClick: onClick,
    title: title,
    onMouseEnter: () => setHover(true),
    onMouseLeave: () => setHover(false),
    style: {
      width: 28,
      height: 28,
      display: 'flex',
      alignItems: 'center',
      justifyContent: 'center',
      background: hover ? 'rgba(255,255,255,0.12)' : 'rgba(255,255,255,0.04)',
      border: '1px solid rgba(255,255,255,0.1)',
      borderRadius: 6,
      color: '#f6f4ef',
      cursor: 'pointer',
      padding: 0,
      transition: 'background 120ms'
    }
  }, children);
}
Object.assign(window, {
  Easing,
  interpolate,
  animate,
  clamp,
  TimelineContext,
  useTime,
  useTimeline,
  Sprite,
  SpriteContext,
  useSprite,
  TextSprite,
  ImageSprite,
  RectSprite,
  Stage,
  PlaybackBar
});
})(); } catch (e) { __ds_ns.__errors.push({ path: "animations.jsx", error: String((e && e.message) || e) }); }

// slides/deck-stage.js
try { (() => {
/**
 * <deck-stage> — reusable web component for HTML decks.
 *
 * Handles:
 *  (a) speaker notes — reads <script type="application/json" id="speaker-notes">
 *      and posts {slideIndexChanged: N} to the parent window on nav.
 *  (b) keyboard navigation — ←/→, PgUp/PgDn, Space, Home/End, number keys.
 *  (c) press R to reset to slide 0 (with a tasteful keyboard hint).
 *  (d) bottom-center overlay showing slide count + hints, fades out on idle.
 *  (e) auto-scaling — inner canvas is a fixed design size (default 1920×1080)
 *      scaled with `transform: scale()` to fit the viewport, letterboxed.
 *      Set the `noscale` attribute to render at authored size (1:1) — the
 *      PPTX exporter sets this so its DOM capture sees unscaled geometry.
 *  (f) print — `@media print` lays every slide out as its own page at the
 *      design size, so the browser's Print → Save as PDF produces a clean
 *      one-page-per-slide PDF with no extra setup.
 *
 * Slides are HIDDEN, not unmounted. Non-active slides stay in the DOM with
 * `visibility: hidden` + `opacity: 0`, so their state (videos, iframes,
 * form inputs, React trees) is preserved across navigation.
 *
 * Lifecycle event — the component dispatches a `slidechange` CustomEvent on
 * itself whenever the active slide changes (including the initial mount).
 * The event bubbles and composes out of shadow DOM, so you can listen on
 * the <deck-stage> element or on document:
 *
 *   document.querySelector('deck-stage').addEventListener('slidechange', (e) => {
 *     e.detail.index         // new 0-based index
 *     e.detail.previousIndex // previous index, or -1 on init
 *     e.detail.total         // total slide count
 *     e.detail.slide         // the new active slide element
 *     e.detail.previousSlide // the prior slide element, or null on init
 *     e.detail.reason        // 'init' | 'keyboard' | 'click' | 'tap' | 'api'
 *   });
 *
 * Persistence: none at the deck level. The host app keeps the current slide
 * in its own URL (?slide=) and re-delivers it via location.hash on load, so a
 * bare load with no hash always starts at slide 1.
 *
 * Usage:
 *   <deck-stage width="1920" height="1080">
 *     <section data-label="Title">...</section>
 *     <section data-label="Agenda">...</section>
 *   </deck-stage>
 *
 * Slides are the direct element children of <deck-stage>. Each slide is
 * automatically tagged with:
 *   - data-screen-label="NN Label"   (1-indexed, for comment flow)
 *   - data-om-validate="no_overflowing_text,no_overlapping_text,slide_sized_text"
 */

(() => {
  const DESIGN_W_DEFAULT = 1920;
  const DESIGN_H_DEFAULT = 1080;
  const OVERLAY_HIDE_MS = 1800;
  const VALIDATE_ATTR = 'no_overflowing_text,no_overlapping_text,slide_sized_text';
  const pad2 = n => String(n).padStart(2, '0');
  const stylesheet = `
    :host {
      position: fixed;
      inset: 0;
      display: block;
      background: #000;
      color: #fff;
      font-family: -apple-system, BlinkMacSystemFont, "Helvetica Neue", Helvetica, Arial, sans-serif;
      overflow: hidden;
    }

    .stage {
      position: absolute;
      inset: 0;
      display: flex;
      align-items: center;
      justify-content: center;
    }

    .canvas {
      position: relative;
      transform-origin: center center;
      flex-shrink: 0;
      background: #fff;
      will-change: transform;
    }

    /* Slides live in light DOM (via <slot>) so authored CSS still applies.
       We absolutely position each slotted child to stack them. */
    ::slotted(*) {
      position: absolute !important;
      inset: 0 !important;
      width: 100% !important;
      height: 100% !important;
      box-sizing: border-box !important;
      overflow: hidden;
      opacity: 0;
      pointer-events: none;
      visibility: hidden;
    }
    ::slotted([data-deck-active]) {
      opacity: 1;
      pointer-events: auto;
      visibility: visible;
    }

    /* Tap zones for mobile — back/forward thirds like Stories.
       Transparent, no visible UI, don't block the overlay. */
    .tapzones {
      position: fixed;
      inset: 0;
      display: flex;
      z-index: 2147482000;
      pointer-events: none;
    }
    .tapzone {
      flex: 1;
      pointer-events: auto;
      -webkit-tap-highlight-color: transparent;
    }
    /* Only activate tap zones on coarse pointers (touch devices). */
    @media (hover: hover) and (pointer: fine) {
      .tapzones { display: none; }
    }

    .overlay {
      position: fixed;
      left: 50%;
      bottom: 22px;
      transform: translate(-50%, 6px) scale(0.92);
      filter: blur(6px);
      display: flex;
      align-items: center;
      gap: 4px;
      padding: 4px;
      background: #000;
      color: #fff;
      border-radius: 999px;
      font-size: 12px;
      font-feature-settings: "tnum" 1;
      letter-spacing: 0.01em;
      opacity: 0;
      pointer-events: none;
      transition: opacity 260ms ease, transform 260ms cubic-bezier(.2,.8,.2,1), filter 260ms ease;
      transform-origin: center bottom;
      z-index: 2147483000;
      user-select: none;
    }
    .overlay[data-visible] {
      opacity: 1;
      pointer-events: auto;
      transform: translate(-50%, 0) scale(1);
      filter: blur(0);
    }

    .btn {
      appearance: none;
      -webkit-appearance: none;
      background: transparent;
      border: 0;
      margin: 0;
      padding: 0;
      color: inherit;
      font: inherit;
      cursor: default;
      display: inline-flex;
      align-items: center;
      justify-content: center;
      height: 28px;
      min-width: 28px;
      border-radius: 999px;
      color: rgba(255,255,255,0.72);
      transition: background 140ms ease, color 140ms ease;
      -webkit-tap-highlight-color: transparent;
    }
    .btn:hover { background: rgba(255,255,255,0.12); color: #fff; }
    .btn:active { background: rgba(255,255,255,0.18); }
    .btn:focus { outline: none; }
    .btn:focus-visible { outline: none; }
    .btn::-moz-focus-inner { border: 0; }
    .btn svg { width: 14px; height: 14px; display: block; }
    .btn.reset {
      font-size: 11px;
      font-weight: 500;
      letter-spacing: 0.02em;
      padding: 0 10px 0 12px;
      gap: 6px;
      color: rgba(255,255,255,0.72);
    }
    .btn.reset .kbd {
      display: inline-flex;
      align-items: center;
      justify-content: center;
      min-width: 16px;
      height: 16px;
      padding: 0 4px;
      font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace;
      font-size: 10px;
      line-height: 1;
      color: rgba(255,255,255,0.88);
      background: rgba(255,255,255,0.12);
      border-radius: 4px;
    }

    .count {
      font-variant-numeric: tabular-nums;
      color: #fff;
      font-weight: 500;
      padding: 0 8px;
      min-width: 42px;
      text-align: center;
      font-size: 12px;
    }
    .count .sep { color: rgba(255,255,255,0.45); margin: 0 3px; font-weight: 400; }
    .count .total { color: rgba(255,255,255,0.55); }

    .divider {
      width: 1px;
      height: 14px;
      background: rgba(255,255,255,0.18);
      margin: 0 2px;
    }

    /* ── Print: one page per slide, no chrome ────────────────────────────
       The screen layout stacks every slide at inset:0 inside a scaled
       canvas; for print we want them in document flow at the authored
       design size so the browser paginates one slide per sheet. The
       @page size is set from the width/height attributes via the inline
       <style id="deck-stage-print-page"> that connectedCallback injects
       into <head> (the @page at-rule has no effect inside shadow DOM). */
    @media print {
      :host {
        position: static;
        inset: auto;
        background: none;
        overflow: visible;
        color: inherit;
      }
      .stage { position: static; display: block; }
      .canvas {
        transform: none !important;
        width: auto !important;
        height: auto !important;
        background: none;
        will-change: auto;
      }
      ::slotted(*) {
        position: relative !important;
        inset: auto !important;
        width: var(--deck-design-w) !important;
        height: var(--deck-design-h) !important;
        box-sizing: border-box !important;
        opacity: 1 !important;
        visibility: visible !important;
        pointer-events: auto;
        break-after: page;
        page-break-after: always;
        break-inside: avoid;
        overflow: hidden;
      }
      ::slotted(*:last-child) {
        break-after: auto;
        page-break-after: auto;
      }
      .overlay, .tapzones { display: none !important; }
    }
  `;
  class DeckStage extends HTMLElement {
    static get observedAttributes() {
      return ['width', 'height', 'noscale'];
    }
    constructor() {
      super();
      this._root = this.attachShadow({
        mode: 'open'
      });
      this._index = 0;
      this._slides = [];
      this._notes = [];
      this._hideTimer = null;
      this._mouseIdleTimer = null;
      this._onKey = this._onKey.bind(this);
      this._onResize = this._onResize.bind(this);
      this._onSlotChange = this._onSlotChange.bind(this);
      this._onMouseMove = this._onMouseMove.bind(this);
      this._onTapBack = this._onTapBack.bind(this);
      this._onTapForward = this._onTapForward.bind(this);
    }
    get designWidth() {
      return parseInt(this.getAttribute('width'), 10) || DESIGN_W_DEFAULT;
    }
    get designHeight() {
      return parseInt(this.getAttribute('height'), 10) || DESIGN_H_DEFAULT;
    }
    connectedCallback() {
      this._render();
      this._loadNotes();
      this._syncPrintPageRule();
      window.addEventListener('keydown', this._onKey);
      window.addEventListener('resize', this._onResize);
      window.addEventListener('mousemove', this._onMouseMove, {
        passive: true
      });
      // Initial collection + layout happens via slotchange, which fires on mount.
    }
    disconnectedCallback() {
      window.removeEventListener('keydown', this._onKey);
      window.removeEventListener('resize', this._onResize);
      window.removeEventListener('mousemove', this._onMouseMove);
      if (this._hideTimer) clearTimeout(this._hideTimer);
      if (this._mouseIdleTimer) clearTimeout(this._mouseIdleTimer);
    }
    attributeChangedCallback() {
      if (this._canvas) {
        this._canvas.style.width = this.designWidth + 'px';
        this._canvas.style.height = this.designHeight + 'px';
        this._canvas.style.setProperty('--deck-design-w', this.designWidth + 'px');
        this._canvas.style.setProperty('--deck-design-h', this.designHeight + 'px');
        this._fit();
        this._syncPrintPageRule();
      }
    }
    _render() {
      const style = document.createElement('style');
      style.textContent = stylesheet;
      const stage = document.createElement('div');
      stage.className = 'stage';
      const canvas = document.createElement('div');
      canvas.className = 'canvas';
      canvas.style.width = this.designWidth + 'px';
      canvas.style.height = this.designHeight + 'px';
      canvas.style.setProperty('--deck-design-w', this.designWidth + 'px');
      canvas.style.setProperty('--deck-design-h', this.designHeight + 'px');
      const slot = document.createElement('slot');
      slot.addEventListener('slotchange', this._onSlotChange);
      canvas.appendChild(slot);
      stage.appendChild(canvas);

      // Tap zones (mobile): left third = back, right third = forward.
      const tapzones = document.createElement('div');
      tapzones.className = 'tapzones export-hidden';
      tapzones.setAttribute('aria-hidden', 'true');
      tapzones.setAttribute('data-noncommentable', '');
      const tzBack = document.createElement('div');
      tzBack.className = 'tapzone tapzone--back';
      const tzMid = document.createElement('div');
      tzMid.className = 'tapzone tapzone--mid';
      tzMid.style.pointerEvents = 'none';
      const tzFwd = document.createElement('div');
      tzFwd.className = 'tapzone tapzone--fwd';
      tzBack.addEventListener('click', this._onTapBack);
      tzFwd.addEventListener('click', this._onTapForward);
      tapzones.append(tzBack, tzMid, tzFwd);

      // Overlay: compact, solid black, with clickable controls.
      const overlay = document.createElement('div');
      overlay.className = 'overlay export-hidden';
      overlay.setAttribute('role', 'toolbar');
      overlay.setAttribute('aria-label', 'Deck controls');
      overlay.setAttribute('data-noncommentable', '');
      overlay.innerHTML = `
        <button class="btn prev" type="button" aria-label="Previous slide" title="Previous (←)">
          <svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M10 3L5 8l5 5"/></svg>
        </button>
        <span class="count" aria-live="polite"><span class="current">1</span><span class="sep">/</span><span class="total">1</span></span>
        <button class="btn next" type="button" aria-label="Next slide" title="Next (→)">
          <svg viewBox="0 0 16 16" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true"><path d="M6 3l5 5-5 5"/></svg>
        </button>
        <span class="divider"></span>
        <button class="btn reset" type="button" aria-label="Reset to first slide" title="Reset (R)">Reset<span class="kbd">R</span></button>
      `;
      overlay.querySelector('.prev').addEventListener('click', () => this._go(this._index - 1, 'click'));
      overlay.querySelector('.next').addEventListener('click', () => this._go(this._index + 1, 'click'));
      overlay.querySelector('.reset').addEventListener('click', () => this._go(0, 'click'));
      this._root.append(style, stage, tapzones, overlay);
      this._canvas = canvas;
      this._slot = slot;
      this._overlay = overlay;
      this._countEl = overlay.querySelector('.current');
      this._totalEl = overlay.querySelector('.total');
    }

    /** @page must live in the document stylesheet — it's a no-op inside
     *  shadow DOM. Inject/update a single <head> style tag so the print
     *  sheet matches the design size and Save-as-PDF yields one slide per
     *  page with no margins. */
    _syncPrintPageRule() {
      const id = 'deck-stage-print-page';
      let tag = document.getElementById(id);
      if (!tag) {
        tag = document.createElement('style');
        tag.id = id;
        document.head.appendChild(tag);
      }
      tag.textContent = '@page { size: ' + this.designWidth + 'px ' + this.designHeight + 'px; margin: 0; } ' + '@media print { html, body { margin: 0 !important; padding: 0 !important; background: none !important; overflow: visible !important; height: auto !important; } ' + '* { -webkit-print-color-adjust: exact; print-color-adjust: exact; } }';
    }
    _onSlotChange() {
      this._collectSlides();
      this._restoreIndex();
      this._applyIndex({
        showOverlay: false,
        broadcast: true,
        reason: 'init'
      });
      this._fit();
    }
    _collectSlides() {
      const assigned = this._slot.assignedElements({
        flatten: true
      });
      this._slides = assigned.filter(el => {
        // Skip template/style/script nodes even if someone slots them.
        const tag = el.tagName;
        return tag !== 'TEMPLATE' && tag !== 'SCRIPT' && tag !== 'STYLE';
      });
      this._slides.forEach((slide, i) => {
        const n = i + 1;
        // Determine a label for comment flow: prefer explicit data-label,
        // then an existing data-screen-label, then first heading, else "Slide".
        let label = slide.getAttribute('data-label');
        if (!label) {
          const existing = slide.getAttribute('data-screen-label');
          if (existing) {
            // Strip any leading number the author may have included.
            label = existing.replace(/^\s*\d+\s*/, '').trim() || existing;
          }
        }
        if (!label) {
          const h = slide.querySelector('h1, h2, h3, [data-title]');
          if (h) label = (h.textContent || '').trim().slice(0, 40);
        }
        if (!label) label = 'Slide';
        slide.setAttribute('data-screen-label', `${pad2(n)} ${label}`);

        // Validation attribute for comment flow / auto-checks.
        if (!slide.hasAttribute('data-om-validate')) {
          slide.setAttribute('data-om-validate', VALIDATE_ATTR);
        }
        slide.setAttribute('data-deck-slide', String(i));
      });
      if (this._totalEl) this._totalEl.textContent = String(this._slides.length || 1);
      if (this._index >= this._slides.length) this._index = Math.max(0, this._slides.length - 1);
    }
    _loadNotes() {
      const tag = document.getElementById('speaker-notes');
      if (!tag) {
        this._notes = [];
        return;
      }
      try {
        const parsed = JSON.parse(tag.textContent || '[]');
        if (Array.isArray(parsed)) this._notes = parsed;
      } catch (e) {
        console.warn('[deck-stage] Failed to parse #speaker-notes JSON:', e);
        this._notes = [];
      }
    }
    _restoreIndex() {
      // The host's ?slide= param is delivered as a #<int> hash (1-indexed) on
      // the iframe src. No hash → slide 1; the deck itself keeps no position
      // state across loads.
      const h = (location.hash || '').match(/^#(\d+)$/);
      if (h) {
        const n = parseInt(h[1], 10) - 1;
        if (n >= 0 && n < this._slides.length) this._index = n;
      }
    }
    _applyIndex({
      showOverlay = true,
      broadcast = true,
      reason = 'init'
    } = {}) {
      if (!this._slides.length) return;
      const prev = this._prevIndex == null ? -1 : this._prevIndex;
      const curr = this._index;
      // Keep the iframe's own hash in sync so an in-iframe location.reload()
      // (reload banner path in viewer-handle.ts) lands on the current slide,
      // not the stale deep-link hash from initial load.
      try {
        history.replaceState(null, '', '#' + (curr + 1));
      } catch (e) {}
      this._slides.forEach((s, i) => {
        if (i === curr) s.setAttribute('data-deck-active', '');else s.removeAttribute('data-deck-active');
      });
      if (this._countEl) this._countEl.textContent = String(curr + 1);
      if (broadcast) {
        // (1) Legacy: host-window postMessage for speaker-notes renderers.
        try {
          window.postMessage({
            slideIndexChanged: curr
          }, '*');
        } catch (e) {}

        // (2) In-page CustomEvent on the <deck-stage> element itself.
        //     Bubbles and composes out of shadow DOM so slide code can listen:
        //       document.querySelector('deck-stage').addEventListener('slidechange', e => {
        //         e.detail.index, e.detail.previousIndex, e.detail.total, e.detail.slide, e.detail.reason
        //       });
        const detail = {
          index: curr,
          previousIndex: prev,
          total: this._slides.length,
          slide: this._slides[curr] || null,
          previousSlide: prev >= 0 ? this._slides[prev] || null : null,
          reason: reason // 'init' | 'keyboard' | 'click' | 'tap' | 'api'
        };
        this.dispatchEvent(new CustomEvent('slidechange', {
          detail,
          bubbles: true,
          composed: true
        }));
      }
      this._prevIndex = curr;
      if (showOverlay) this._flashOverlay();
    }
    _flashOverlay() {
      if (!this._overlay) return;
      this._overlay.setAttribute('data-visible', '');
      if (this._hideTimer) clearTimeout(this._hideTimer);
      this._hideTimer = setTimeout(() => {
        this._overlay.removeAttribute('data-visible');
      }, OVERLAY_HIDE_MS);
    }
    _fit() {
      if (!this._canvas) return;
      // PPTX export sets noscale so the DOM capture sees authored-size
      // geometry — the scaled canvas is in shadow DOM, so the exporter's
      // resetTransformSelector can't reach .canvas.style.transform directly.
      if (this.hasAttribute('noscale')) {
        this._canvas.style.transform = 'none';
        return;
      }
      const vw = window.innerWidth;
      const vh = window.innerHeight;
      const s = Math.min(vw / this.designWidth, vh / this.designHeight);
      this._canvas.style.transform = `scale(${s})`;
    }
    _onResize() {
      this._fit();
    }
    _onMouseMove() {
      // Keep overlay visible while mouse moves; hide after idle.
      this._flashOverlay();
    }
    _onTapBack(e) {
      e.preventDefault();
      this._go(this._index - 1, 'tap');
    }
    _onTapForward(e) {
      e.preventDefault();
      this._go(this._index + 1, 'tap');
    }
    _onKey(e) {
      // Ignore when the user is typing.
      const t = e.target;
      if (t && (t.isContentEditable || /^(INPUT|TEXTAREA|SELECT)$/.test(t.tagName))) return;
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      const key = e.key;
      let handled = true;
      if (key === 'ArrowRight' || key === 'PageDown' || key === ' ' || key === 'Spacebar') {
        this._go(this._index + 1, 'keyboard');
      } else if (key === 'ArrowLeft' || key === 'PageUp') {
        this._go(this._index - 1, 'keyboard');
      } else if (key === 'Home') {
        this._go(0, 'keyboard');
      } else if (key === 'End') {
        this._go(this._slides.length - 1, 'keyboard');
      } else if (key === 'r' || key === 'R') {
        this._go(0, 'keyboard');
      } else if (/^[0-9]$/.test(key)) {
        // 1..9 jump to that slide; 0 jumps to 10.
        const n = key === '0' ? 9 : parseInt(key, 10) - 1;
        if (n < this._slides.length) this._go(n, 'keyboard');
      } else {
        handled = false;
      }
      if (handled) {
        e.preventDefault();
        this._flashOverlay();
      }
    }
    _go(i, reason = 'api') {
      if (!this._slides.length) return;
      const clamped = Math.max(0, Math.min(this._slides.length - 1, i));
      if (clamped === this._index) {
        this._flashOverlay();
        return;
      }
      this._index = clamped;
      this._applyIndex({
        showOverlay: true,
        broadcast: true,
        reason
      });
    }

    // Public API ------------------------------------------------------------

    /** Current slide index (0-based). */
    get index() {
      return this._index;
    }
    /** Total slide count. */
    get length() {
      return this._slides.length;
    }
    /** Programmatically navigate. */
    goTo(i) {
      this._go(i, 'api');
    }
    next() {
      this._go(this._index + 1, 'api');
    }
    prev() {
      this._go(this._index - 1, 'api');
    }
    reset() {
      this._go(0, 'api');
    }
  }
  if (!customElements.get('deck-stage')) {
    customElements.define('deck-stage', DeckStage);
  }
})();
})(); } catch (e) { __ds_ns.__errors.push({ path: "slides/deck-stage.js", error: String((e && e.message) || e) }); }

// ui_kits/docs/ContentPage.jsx
try { (() => {
// ContentPage.jsx — Doc page with prose, code blocks, callouts
// Load with <script type="text/babel" src="ContentPage.jsx"></script>

const InlineCode = ({
  children,
  dark
}) => /*#__PURE__*/React.createElement("code", {
  style: {
    fontFamily: 'Space Mono, monospace',
    fontSize: '0.875em',
    background: dark ? 'hsl(220,14%,17%)' : '#f1f3f5',
    color: dark ? '#f472b6' : '#CB3C33',
    padding: '0.15em 0.4em',
    borderRadius: 3,
    border: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`
  }
}, children);
const CodeBlock = ({
  dark,
  children,
  lang = 'julia'
}) => /*#__PURE__*/React.createElement("pre", {
  style: {
    background: 'hsl(220,20%,9%)',
    color: '#d4d4d4',
    fontFamily: 'Space Mono, monospace',
    fontSize: 13,
    lineHeight: 1.65,
    padding: '20px 24px',
    borderRadius: 8,
    overflowX: 'auto',
    margin: '16px 0',
    border: '1px solid hsl(220,12%,23%)'
  }
}, /*#__PURE__*/React.createElement("code", {
  dangerouslySetInnerHTML: {
    __html: children
  }
}));
const Callout = ({
  type = 'tip',
  title,
  children,
  dark
}) => {
  const colors = {
    tip: {
      border: '#389826',
      bg: dark ? 'hsl(130,30%,10%)' : '#f0faea',
      title: '#2c7a1e'
    },
    warning: {
      border: '#b45309',
      bg: dark ? 'hsl(35,40%,12%)' : '#fef3c7',
      title: '#78350f'
    },
    danger: {
      border: '#CB3C33',
      bg: dark ? 'hsl(5,30%,10%)' : '#fef2f2',
      title: '#9a2e27'
    }
  };
  const c = colors[type] || colors.tip;
  return /*#__PURE__*/React.createElement("div", {
    style: {
      borderLeft: `3px solid ${c.border}`,
      background: c.bg,
      padding: '12px 16px',
      borderRadius: '0 6px 6px 0',
      margin: '16px 0'
    }
  }, /*#__PURE__*/React.createElement("div", {
    style: {
      fontSize: 13,
      fontWeight: 600,
      color: c.title,
      marginBottom: 4,
      fontFamily: 'Inter, sans-serif'
    }
  }, title || type.toUpperCase()), /*#__PURE__*/React.createElement("div", {
    style: {
      fontSize: 13,
      color: dark ? '#adb5bd' : '#495057',
      lineHeight: 1.6,
      fontFamily: 'Inter, sans-serif'
    }
  }, children));
};
const Badge = ({
  children,
  color = 'green'
}) => {
  const colors = {
    green: {
      bg: '#dcf5d7',
      fg: '#2c7a1e'
    },
    blue: {
      bg: '#d4dcf9',
      fg: '#2d47a8'
    },
    purple: {
      bg: '#eadaf2',
      fg: '#6b3a85'
    }
  };
  const c = colors[color] || colors.green;
  return /*#__PURE__*/React.createElement("span", {
    style: {
      background: c.bg,
      color: c.fg,
      fontSize: 11,
      fontWeight: 600,
      padding: '2px 8px',
      borderRadius: 3,
      fontFamily: 'Inter, sans-serif',
      display: 'inline-block',
      marginRight: 6
    }
  }, children);
};
const ContentPage = ({
  dark
}) => {
  const s = contentStyles;
  return /*#__PURE__*/React.createElement("article", {
    style: s.article
  }, /*#__PURE__*/React.createElement("div", {
    style: s.breadcrumb(dark)
  }, "Source Code ", /*#__PURE__*/React.createElement("span", {
    style: {
      color: '#adb5bd'
    }
  }, "/"), " Methods ", /*#__PURE__*/React.createElement("span", {
    style: {
      color: '#adb5bd'
    }
  }, "/"), " Clipping"), /*#__PURE__*/React.createElement("h1", {
    style: s.h1(dark)
  }, "Polygon Clipping"), /*#__PURE__*/React.createElement("div", {
    style: {
      marginBottom: 20
    }
  }, /*#__PURE__*/React.createElement(Badge, null, "Source Code"), /*#__PURE__*/React.createElement(Badge, {
    color: "blue"
  }, "Pure Julia")), /*#__PURE__*/React.createElement("p", {
    style: s.p(dark)
  }, "GeometryOps provides polygon clipping via three operations: ", /*#__PURE__*/React.createElement(InlineCode, {
    dark: dark
  }, "intersection"), ", ", /*#__PURE__*/React.createElement(InlineCode, {
    dark: dark
  }, "difference"), ", and ", /*#__PURE__*/React.createElement(InlineCode, {
    dark: dark
  }, "union"), ". All methods accept any ", /*#__PURE__*/React.createElement(InlineCode, {
    dark: dark
  }, "GeoInterface.jl"), "-compatible geometry type and return the same."), /*#__PURE__*/React.createElement(Callout, {
    type: "warning",
    title: "Warning",
    dark: dark
  }, "This package is still under heavy development! Results may change between minor versions."), /*#__PURE__*/React.createElement("h2", {
    style: s.h2(dark)
  }, "Usage"), /*#__PURE__*/React.createElement("p", {
    style: s.p(dark)
  }, "The simplest way to compute an intersection is with the ", /*#__PURE__*/React.createElement(InlineCode, {
    dark: dark
  }, "GO.intersection"), " function. Specify a ", /*#__PURE__*/React.createElement(InlineCode, {
    dark: dark
  }, "target"), " trait to control what geometry type is returned."), /*#__PURE__*/React.createElement(CodeBlock, {
    dark: dark
  }, `<span style="color:#91dd33">using</span> GeometryOps <span style="color:#91dd33">as</span> GO
<span style="color:#91dd33">using</span> GeoInterface <span style="color:#91dd33">as</span> GI

<span style="color:#6c8a6c; font-style:italic"># Two overlapping polygons</span>
poly1 = GI.Polygon([[(<span style="color:#f9a8d4">0</span>, <span style="color:#f9a8d4">0</span>), (<span style="color:#f9a8d4">2</span>, <span style="color:#f9a8d4">0</span>), (<span style="color:#f9a8d4">2</span>, <span style="color:#f9a8d4">2</span>), (<span style="color:#f9a8d4">0</span>, <span style="color:#f9a8d4">2</span>), (<span style="color:#f9a8d4">0</span>, <span style="color:#f9a8d4">0</span>)]])
poly2 = GI.Polygon([[(<span style="color:#f9a8d4">1</span>, <span style="color:#f9a8d4">0</span>), (<span style="color:#f9a8d4">3</span>, <span style="color:#f9a8d4">0</span>), (<span style="color:#f9a8d4">3</span>, <span style="color:#f9a8d4">2</span>), (<span style="color:#f9a8d4">1</span>, <span style="color:#f9a8d4">2</span>), (<span style="color:#f9a8d4">1</span>, <span style="color:#f9a8d4">0</span>)]])

result = GO.<span style="color:#7dd3fc">intersection</span>(poly1, poly2; target = GI.PolygonTrait())
GO.<span style="color:#7dd3fc">area</span>(result)  <span style="color:#6c8a6c; font-style:italic"># → 2.0</span>`), /*#__PURE__*/React.createElement(Callout, {
    type: "tip",
    title: "Tip",
    dark: dark
  }, "Use ", /*#__PURE__*/React.createElement(InlineCode, {
    dark: dark
  }, "GO.apply"), " with a target trait to apply operations over large nested geometry collections efficiently."), /*#__PURE__*/React.createElement("h2", {
    style: s.h2(dark)
  }, "Algorithm"), /*#__PURE__*/React.createElement("p", {
    style: s.p(dark)
  }, "Clipping is implemented using the Sutherland\u2013Hodgman algorithm for convex polygons, and the Greiner\u2013Hormann algorithm for the general case. The implementation is in pure Julia with no external C dependencies."), /*#__PURE__*/React.createElement("h3", {
    style: s.h3(dark)
  }, "Available functions"), /*#__PURE__*/React.createElement("div", {
    style: s.apiTable(dark)
  }, [['intersection(a, b; target)', 'Returns the intersection of two geometries'], ['difference(a, b; target)', 'Returns a minus b'], ['union(a, b; target)', 'Returns the union of two geometries']].map(([fn, desc]) => /*#__PURE__*/React.createElement("div", {
    key: fn,
    style: s.apiRow(dark)
  }, /*#__PURE__*/React.createElement("div", {
    style: s.apiFn(dark)
  }, /*#__PURE__*/React.createElement(InlineCode, {
    dark: dark
  }, fn)), /*#__PURE__*/React.createElement("div", {
    style: s.apiDesc(dark)
  }, desc)))));
};
const contentStyles = {
  article: {
    maxWidth: 860,
    padding: '40px 48px 80px',
    fontFamily: 'Inter, sans-serif'
  },
  breadcrumb: dark => ({
    fontSize: 13,
    color: dark ? 'hsl(220,8%,56%)' : '#adb5bd',
    marginBottom: 16,
    display: 'flex',
    gap: 8,
    alignItems: 'center'
  }),
  h1: dark => ({
    fontSize: 36,
    fontWeight: 700,
    color: dark ? '#f8f9fa' : '#212529',
    letterSpacing: '-0.025em',
    marginBottom: 12,
    lineHeight: 1.15
  }),
  h2: dark => ({
    fontSize: 24,
    fontWeight: 600,
    color: dark ? '#f8f9fa' : '#212529',
    letterSpacing: '-0.02em',
    marginTop: 36,
    marginBottom: 12,
    lineHeight: 1.25,
    borderTop: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`,
    paddingTop: 28
  }),
  h3: dark => ({
    fontSize: 18,
    fontWeight: 600,
    color: dark ? '#e9ecef' : '#343a40',
    marginTop: 24,
    marginBottom: 10
  }),
  p: dark => ({
    fontSize: 15,
    color: dark ? '#adb5bd' : '#495057',
    lineHeight: 1.7,
    marginBottom: 16
  }),
  apiTable: dark => ({
    border: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`,
    borderRadius: 8,
    overflow: 'hidden'
  }),
  apiRow: dark => ({
    display: 'flex',
    gap: 16,
    padding: '10px 16px',
    borderBottom: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`,
    alignItems: 'baseline',
    background: 'transparent'
  }),
  apiFn: dark => ({
    minWidth: 260,
    flexShrink: 0
  }),
  apiDesc: dark => ({
    fontSize: 13,
    color: dark ? '#6c757d' : '#6c757d',
    lineHeight: 1.5
  })
};
Object.assign(window, {
  ContentPage,
  InlineCode,
  CodeBlock,
  Callout,
  Badge
});
})(); } catch (e) { __ds_ns.__errors.push({ path: "ui_kits/docs/ContentPage.jsx", error: String((e && e.message) || e) }); }

// ui_kits/docs/Header.jsx
try { (() => {
// Header.jsx — JuliaGeo Docs top navigation
// Load with <script type="text/babel" src="Header.jsx"></script>

const Header = ({
  darkMode,
  onToggleDark,
  currentPage,
  onNav
}) => {
  const s = headerStyles;
  return /*#__PURE__*/React.createElement("header", {
    style: s.header(darkMode)
  }, /*#__PURE__*/React.createElement("div", {
    style: s.inner
  }, /*#__PURE__*/React.createElement("button", {
    style: s.logo,
    onClick: () => onNav('home')
  }, /*#__PURE__*/React.createElement("img", {
    src: "../../assets/juliageo-logo.svg",
    alt: "JuliaGeo",
    style: {
      height: 28,
      width: 'auto'
    }
  }), /*#__PURE__*/React.createElement("span", {
    style: s.logoName(darkMode)
  }, "GeometryOps", /*#__PURE__*/React.createElement("span", {
    style: s.logoJl(darkMode)
  }, ".jl"))), /*#__PURE__*/React.createElement("div", {
    style: s.searchWrap(darkMode)
  }, /*#__PURE__*/React.createElement("svg", {
    width: "14",
    height: "14",
    viewBox: "0 0 24 24",
    fill: "none",
    stroke: darkMode ? '#6c757d' : '#adb5bd',
    strokeWidth: "2"
  }, /*#__PURE__*/React.createElement("circle", {
    cx: "11",
    cy: "11",
    r: "8"
  }), /*#__PURE__*/React.createElement("line", {
    x1: "21",
    y1: "21",
    x2: "16.65",
    y2: "16.65"
  })), /*#__PURE__*/React.createElement("span", {
    style: s.searchText(darkMode)
  }, "Search docs\u2026"), /*#__PURE__*/React.createElement("span", {
    style: s.searchKbd(darkMode)
  }, "\u2318K")), /*#__PURE__*/React.createElement("nav", {
    style: s.nav
  }, ['Introduction', 'API', 'Tutorials', 'Source'].map(label => /*#__PURE__*/React.createElement("button", {
    key: label,
    style: s.navLink(darkMode, currentPage === label.toLowerCase()),
    onClick: () => onNav(label.toLowerCase())
  }, label)), /*#__PURE__*/React.createElement("div", {
    style: s.divider(darkMode)
  }), /*#__PURE__*/React.createElement("a", {
    href: "https://github.com/JuliaGeo/GeometryOps.jl",
    target: "_blank",
    rel: "noreferrer",
    style: s.iconBtn(darkMode),
    title: "GitHub"
  }, /*#__PURE__*/React.createElement("svg", {
    width: "18",
    height: "18",
    viewBox: "0 0 24 24",
    fill: darkMode ? '#adb5bd' : '#6c757d'
  }, /*#__PURE__*/React.createElement("path", {
    d: "M12 2C6.477 2 2 6.484 2 12.017c0 4.425 2.865 8.18 6.839 9.504.5.092.682-.217.682-.483 0-.237-.008-.868-.013-1.703-2.782.605-3.369-1.343-3.369-1.343-.454-1.158-1.11-1.466-1.11-1.466-.908-.62.069-.608.069-.608 1.003.07 1.531 1.032 1.531 1.032.892 1.53 2.341 1.088 2.91.832.092-.647.35-1.088.636-1.338-2.22-.253-4.555-1.113-4.555-4.951 0-1.093.39-1.988 1.029-2.688-.103-.253-.446-1.272.098-2.65 0 0 .84-.27 2.75 1.026A9.564 9.564 0 0112 6.844c.85.004 1.705.115 2.504.337 1.909-1.296 2.747-1.027 2.747-1.027.546 1.379.202 2.398.1 2.651.64.7 1.028 1.595 1.028 2.688 0 3.848-2.339 4.695-4.566 4.943.359.309.678.92.678 1.855 0 1.338-.012 2.419-.012 2.747 0 .268.18.58.688.482A10.019 10.019 0 0022 12.017C22 6.484 17.522 2 12 2z"
  }))), /*#__PURE__*/React.createElement("button", {
    style: s.iconBtn(darkMode),
    onClick: onToggleDark,
    title: "Toggle dark mode"
  }, darkMode ? /*#__PURE__*/React.createElement("svg", {
    width: "18",
    height: "18",
    viewBox: "0 0 24 24",
    fill: "none",
    stroke: "#adb5bd",
    strokeWidth: "2"
  }, /*#__PURE__*/React.createElement("circle", {
    cx: "12",
    cy: "12",
    r: "5"
  }), /*#__PURE__*/React.createElement("line", {
    x1: "12",
    y1: "1",
    x2: "12",
    y2: "3"
  }), /*#__PURE__*/React.createElement("line", {
    x1: "12",
    y1: "21",
    x2: "12",
    y2: "23"
  }), /*#__PURE__*/React.createElement("line", {
    x1: "4.22",
    y1: "4.22",
    x2: "5.64",
    y2: "5.64"
  }), /*#__PURE__*/React.createElement("line", {
    x1: "18.36",
    y1: "18.36",
    x2: "19.78",
    y2: "19.78"
  }), /*#__PURE__*/React.createElement("line", {
    x1: "1",
    y1: "12",
    x2: "3",
    y2: "12"
  }), /*#__PURE__*/React.createElement("line", {
    x1: "21",
    y1: "12",
    x2: "23",
    y2: "12"
  }), /*#__PURE__*/React.createElement("line", {
    x1: "4.22",
    y1: "19.78",
    x2: "5.64",
    y2: "18.36"
  }), /*#__PURE__*/React.createElement("line", {
    x1: "18.36",
    y1: "5.64",
    x2: "19.78",
    y2: "4.22"
  })) : /*#__PURE__*/React.createElement("svg", {
    width: "18",
    height: "18",
    viewBox: "0 0 24 24",
    fill: "none",
    stroke: "#6c757d",
    strokeWidth: "2"
  }, /*#__PURE__*/React.createElement("path", {
    d: "M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"
  }))))));
};
const headerStyles = {
  header: dark => ({
    position: 'fixed',
    top: 0,
    left: 0,
    right: 0,
    zIndex: 100,
    height: 56,
    background: dark ? 'hsl(220,20%,9%)' : '#ffffff',
    borderBottom: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`,
    display: 'flex',
    alignItems: 'center'
  }),
  inner: {
    width: '100%',
    maxWidth: 1280,
    margin: '0 auto',
    padding: '0 24px',
    display: 'flex',
    alignItems: 'center',
    gap: 16
  },
  logo: {
    display: 'flex',
    alignItems: 'center',
    gap: 10,
    background: 'none',
    border: 'none',
    cursor: 'pointer',
    padding: 0,
    textDecoration: 'none'
  },
  logoName: dark => ({
    fontSize: 16,
    fontWeight: 700,
    color: dark ? '#f8f9fa' : '#212529',
    fontFamily: 'Inter, sans-serif',
    letterSpacing: '-0.02em',
    whiteSpace: 'nowrap'
  }),
  logoJl: dark => ({
    color: '#389826'
  }),
  searchWrap: dark => ({
    flex: 1,
    maxWidth: 280,
    height: 34,
    background: dark ? 'hsl(220,16%,13%)' : '#f8f9fa',
    border: `1.5px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`,
    borderRadius: 6,
    display: 'flex',
    alignItems: 'center',
    gap: 8,
    padding: '0 10px',
    cursor: 'text'
  }),
  searchText: dark => ({
    flex: 1,
    fontSize: 13,
    color: dark ? 'hsl(220,8%,56%)' : '#adb5bd',
    fontFamily: 'Inter, sans-serif'
  }),
  searchKbd: dark => ({
    fontSize: 11,
    color: dark ? 'hsl(220,8%,56%)' : '#ced4da',
    background: dark ? 'hsl(220,14%,17%)' : '#f1f3f5',
    border: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#dee2e6'}`,
    borderRadius: 4,
    padding: '1px 5px',
    fontFamily: 'Space Mono, monospace'
  }),
  nav: {
    display: 'flex',
    alignItems: 'center',
    gap: 4,
    marginLeft: 'auto'
  },
  navLink: (dark, active) => ({
    background: 'none',
    border: 'none',
    cursor: 'pointer',
    fontSize: 13,
    fontWeight: active ? 600 : 400,
    color: active ? '#389826' : dark ? '#adb5bd' : '#495057',
    padding: '4px 10px',
    borderRadius: 5,
    fontFamily: 'Inter, sans-serif',
    transition: 'color 150ms, background 150ms'
  }),
  divider: dark => ({
    width: 1,
    height: 20,
    background: dark ? 'hsl(220,12%,23%)' : '#e9ecef',
    margin: '0 4px'
  }),
  iconBtn: dark => ({
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    width: 32,
    height: 32,
    borderRadius: 6,
    background: 'none',
    border: 'none',
    cursor: 'pointer',
    textDecoration: 'none',
    transition: 'background 150ms'
  })
};
Object.assign(window, {
  Header
});
})(); } catch (e) { __ds_ns.__errors.push({ path: "ui_kits/docs/Header.jsx", error: String((e && e.message) || e) }); }

// ui_kits/docs/HomePage.jsx
try { (() => {
// HomePage.jsx — JuliaGeo docs hero + feature cards
// Load with <script type="text/babel" src="HomePage.jsx"></script>

const FEATURES = [{
  icon: /*#__PURE__*/React.createElement("svg", {
    width: "24",
    height: "24",
    viewBox: "0 0 24 24",
    fill: "none",
    stroke: "#389826",
    strokeWidth: "1.75"
  }, /*#__PURE__*/React.createElement("polygon", {
    points: "12 2 22 8.5 22 15.5 12 22 2 15.5 2 8.5"
  }), /*#__PURE__*/React.createElement("line", {
    x1: "12",
    y1: "2",
    x2: "12",
    y2: "22"
  }), /*#__PURE__*/React.createElement("line", {
    x1: "2",
    y1: "8.5",
    x2: "22",
    y2: "8.5"
  })),
  title: 'Pure Julia code',
  detail: 'Fast, understandable, extensible functions written entirely in Julia — no C bindings required.',
  link: '/introduction'
}, {
  icon: /*#__PURE__*/React.createElement("svg", {
    width: "24",
    height: "24",
    viewBox: "0 0 24 24",
    fill: "none",
    stroke: "#389826",
    strokeWidth: "1.75"
  }, /*#__PURE__*/React.createElement("path", {
    d: "M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"
  }), /*#__PURE__*/React.createElement("polyline", {
    points: "14 2 14 8 20 8"
  }), /*#__PURE__*/React.createElement("line", {
    x1: "16",
    y1: "13",
    x2: "8",
    y2: "13"
  }), /*#__PURE__*/React.createElement("line", {
    x1: "16",
    y1: "17",
    x2: "8",
    y2: "17"
  })),
  title: 'Literate programming',
  detail: 'Documented source code with embedded examples via Literate.jl — code and explanation together.',
  link: '/source/methods/clipping/cut'
}, {
  icon: /*#__PURE__*/React.createElement("svg", {
    width: "24",
    height: "24",
    viewBox: "0 0 24 24",
    fill: "none",
    stroke: "#389826",
    strokeWidth: "1.75"
  }, /*#__PURE__*/React.createElement("circle", {
    cx: "12",
    cy: "12",
    r: "10"
  }), /*#__PURE__*/React.createElement("line", {
    x1: "2",
    y1: "12",
    x2: "22",
    y2: "12"
  }), /*#__PURE__*/React.createElement("path", {
    d: "M12 2a15.3 15.3 0 0 1 4 10 15.3 15.3 0 0 1-4 10 15.3 15.3 0 0 1-4-10 15.3 15.3 0 0 1 4-10z"
  })),
  title: 'GeoInterface integration',
  detail: 'Use any GeoInterface.jl-compatible geometry — ArchGDAL, Shapefile, GeoJSON, and more.',
  link: 'https://juliageo.org/GeoInterface.jl/stable'
}, {
  icon: /*#__PURE__*/React.createElement("svg", {
    width: "24",
    height: "24",
    viewBox: "0 0 24 24",
    fill: "none",
    stroke: "#389826",
    strokeWidth: "1.75"
  }, /*#__PURE__*/React.createElement("polyline", {
    points: "22 12 18 12 15 21 9 3 6 12 2 12"
  })),
  title: 'Blazing fast',
  detail: 'Outperforms Python and R equivalents on most geometry benchmarks. Pure Julia means full JIT optimization.',
  link: '#benchmarks'
}];
const METHODS = ['equals', 'extent', 'distance', 'crosses', 'contains', 'intersects', 'intersection', 'difference', 'union', 'simplify', 'centroid', 'signed_area', 'segmentize', 'polygonize', 'barycentric_coordinates'];
const HomePage = ({
  dark,
  onNav
}) => {
  const s = homeStyles;
  return /*#__PURE__*/React.createElement("div", {
    style: s.page
  }, /*#__PURE__*/React.createElement("section", {
    style: s.hero(dark)
  }, /*#__PURE__*/React.createElement("div", {
    style: s.heroInner
  }, /*#__PURE__*/React.createElement("img", {
    src: "../../assets/geometryops-logo.png",
    alt: "GeometryOps",
    style: s.heroLogo
  }), /*#__PURE__*/React.createElement("div", null, /*#__PURE__*/React.createElement("h1", {
    style: s.heroTitle(dark)
  }, "GeometryOps", /*#__PURE__*/React.createElement("span", {
    style: {
      color: '#389826'
    }
  }, ".jl")), /*#__PURE__*/React.createElement("p", {
    style: s.heroTagline(dark)
  }, "Blazing fast geometry operations in pure Julia"), /*#__PURE__*/React.createElement("div", {
    style: s.heroActions
  }, /*#__PURE__*/React.createElement("button", {
    style: s.btnPrimary,
    onClick: () => onNav('intro')
  }, "Introduction"), /*#__PURE__*/React.createElement("button", {
    style: s.btnSecondary(dark),
    onClick: () => onNav('api')
  }, "API Reference"), /*#__PURE__*/React.createElement("a", {
    href: "https://github.com/JuliaGeo/GeometryOps.jl",
    target: "_blank",
    rel: "noreferrer",
    style: s.btnGhost(dark)
  }, "View on GitHub"))))), /*#__PURE__*/React.createElement("div", {
    style: s.content
  }, /*#__PURE__*/React.createElement("section", {
    style: s.section
  }, /*#__PURE__*/React.createElement("div", {
    style: s.featureGrid
  }, FEATURES.map(f => /*#__PURE__*/React.createElement("div", {
    key: f.title,
    style: s.featureCard(dark)
  }, /*#__PURE__*/React.createElement("div", {
    style: s.featureIcon(dark)
  }, f.icon), /*#__PURE__*/React.createElement("div", {
    style: s.featureTitle(dark)
  }, f.title), /*#__PURE__*/React.createElement("div", {
    style: s.featureDetail(dark)
  }, f.detail))))), /*#__PURE__*/React.createElement("section", {
    style: s.section
  }, /*#__PURE__*/React.createElement("h2", {
    style: s.h2(dark)
  }, "What is GeometryOps.jl?"), /*#__PURE__*/React.createElement("p", {
    style: s.p(dark)
  }, "GeometryOps.jl is a package for geometric calculations on (primarily 2D) geometries. The driving idea is to unify all the disparate packages for geometric calculations in Julia, and make them ", /*#__PURE__*/React.createElement("a", {
    href: "#",
    style: s.link
  }, "GeoInterface.jl"), "-compatible."), /*#__PURE__*/React.createElement("p", {
    style: s.p(dark)
  }, "Most use cases are driven by GIS and similar Earth data workflows. Methods are always general to any coordinate space.")), /*#__PURE__*/React.createElement("section", {
    style: s.section
  }, /*#__PURE__*/React.createElement("h2", {
    style: s.h2(dark)
  }, "Available methods"), /*#__PURE__*/React.createElement("div", {
    style: s.methodGrid
  }, METHODS.map(m => /*#__PURE__*/React.createElement("span", {
    key: m,
    style: s.methodTag(dark)
  }, m))))));
};
const homeStyles = {
  page: {
    fontFamily: 'Inter, sans-serif'
  },
  hero: dark => ({
    padding: '72px 48px 56px',
    background: dark ? 'hsl(220,20%,9%)' : '#ffffff',
    borderBottom: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`
  }),
  heroInner: {
    maxWidth: 860,
    display: 'flex',
    gap: 48,
    alignItems: 'center'
  },
  heroLogo: {
    width: 120,
    height: 'auto',
    flexShrink: 0
  },
  heroTitle: dark => ({
    fontSize: 48,
    fontWeight: 700,
    color: dark ? '#f8f9fa' : '#212529',
    letterSpacing: '-0.03em',
    margin: '0 0 10px',
    lineHeight: 1.05
  }),
  heroTagline: dark => ({
    fontSize: 18,
    color: dark ? 'hsl(220,8%,56%)' : '#6c757d',
    margin: '0 0 28px',
    lineHeight: 1.4
  }),
  heroActions: {
    display: 'flex',
    gap: 10,
    flexWrap: 'wrap'
  },
  btnPrimary: {
    height: 38,
    padding: '0 20px',
    background: '#389826',
    color: '#fff',
    border: '1.5px solid #389826',
    borderRadius: 6,
    fontSize: 14,
    fontWeight: 500,
    cursor: 'pointer',
    fontFamily: 'Inter, sans-serif'
  },
  btnSecondary: dark => ({
    height: 38,
    padding: '0 20px',
    background: 'transparent',
    color: '#389826',
    border: '1.5px solid #389826',
    borderRadius: 6,
    fontSize: 14,
    fontWeight: 500,
    cursor: 'pointer',
    fontFamily: 'Inter, sans-serif'
  }),
  btnGhost: dark => ({
    height: 38,
    padding: '0 20px',
    background: 'transparent',
    color: dark ? '#adb5bd' : '#495057',
    border: `1.5px solid ${dark ? 'hsl(220,12%,23%)' : '#dee2e6'}`,
    borderRadius: 6,
    fontSize: 14,
    fontWeight: 400,
    cursor: 'pointer',
    fontFamily: 'Inter, sans-serif',
    textDecoration: 'none',
    display: 'inline-flex',
    alignItems: 'center'
  }),
  content: {
    maxWidth: 860,
    padding: '48px 48px 80px'
  },
  section: {
    marginBottom: 52
  },
  featureGrid: {
    display: 'grid',
    gridTemplateColumns: 'repeat(2, 1fr)',
    gap: 16
  },
  featureCard: dark => ({
    background: dark ? 'hsl(220,16%,13%)' : '#ffffff',
    border: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`,
    borderRadius: 10,
    padding: '20px 24px',
    boxShadow: '0 1px 4px rgba(0,0,0,0.06)'
  }),
  featureIcon: dark => ({
    width: 44,
    height: 44,
    background: dark ? 'hsl(130,20%,12%)' : '#f0faea',
    borderRadius: 8,
    display: 'flex',
    alignItems: 'center',
    justifyContent: 'center',
    marginBottom: 12
  }),
  featureTitle: dark => ({
    fontSize: 15,
    fontWeight: 600,
    color: dark ? '#f8f9fa' : '#212529',
    marginBottom: 6
  }),
  featureDetail: dark => ({
    fontSize: 13,
    color: dark ? 'hsl(220,8%,56%)' : '#6c757d',
    lineHeight: 1.6
  }),
  h2: dark => ({
    fontSize: 28,
    fontWeight: 700,
    color: dark ? '#f8f9fa' : '#212529',
    letterSpacing: '-0.025em',
    marginBottom: 16
  }),
  p: dark => ({
    fontSize: 15,
    color: dark ? '#adb5bd' : '#495057',
    lineHeight: 1.7,
    marginBottom: 14
  }),
  link: {
    color: '#2c7a1e',
    textDecoration: 'none'
  },
  methodGrid: {
    display: 'flex',
    flexWrap: 'wrap',
    gap: 8
  },
  methodTag: dark => ({
    display: 'inline-flex',
    padding: '4px 12px',
    borderRadius: 6,
    fontFamily: 'Space Mono, monospace',
    fontSize: 12,
    background: dark ? 'hsl(220,14%,17%)' : '#f1f3f5',
    color: dark ? '#adb5bd' : '#495057',
    border: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`
  })
};
Object.assign(window, {
  HomePage
});
})(); } catch (e) { __ds_ns.__errors.push({ path: "ui_kits/docs/HomePage.jsx", error: String((e && e.message) || e) }); }

// ui_kits/docs/Sidebar.jsx
try { (() => {
// Sidebar.jsx — JuliaGeo Docs navigation sidebar
// Load with <script type="text/babel" src="Sidebar.jsx"></script>

const NAV_SECTIONS = [{
  title: null,
  items: [{
    label: 'Introduction',
    id: 'intro'
  }, {
    label: 'API Reference',
    id: 'api'
  }]
}, {
  title: 'Tutorials',
  items: [{
    label: 'Creating Geometry',
    id: 'tut-creating'
  }, {
    label: 'Spatial Joins',
    id: 'tut-joins'
  }]
}, {
  title: 'Explanations',
  items: [{
    label: 'Paradigms',
    id: 'exp-paradigms'
  }, {
    label: 'Manifolds',
    id: 'exp-manifolds'
  }, {
    label: 'Performance',
    id: 'exp-perf'
  }, {
    label: 'Peculiarities',
    id: 'exp-peculiar'
  }]
}, {
  title: 'GIS Terminology',
  items: [{
    label: 'CRS',
    id: 'gis-crs'
  }, {
    label: 'Winding Order',
    id: 'gis-winding'
  }]
}, {
  title: 'Source Code',
  items: [{
    label: 'Methods / Clipping',
    id: 'src-clipping'
  }, {
    label: 'Methods / Simplify',
    id: 'src-simplify'
  }, {
    label: 'Methods / Distance',
    id: 'src-distance'
  }, {
    label: 'GeometryOpsCore',
    id: 'src-core'
  }]
}];
const Sidebar = ({
  darkMode,
  currentPage,
  onNav
}) => {
  const s = sidebarStyles;
  return /*#__PURE__*/React.createElement("aside", {
    style: s.sidebar(darkMode)
  }, /*#__PURE__*/React.createElement("nav", {
    style: s.nav
  }, NAV_SECTIONS.map((section, i) => /*#__PURE__*/React.createElement("div", {
    key: i,
    style: s.section
  }, section.title && /*#__PURE__*/React.createElement("div", {
    style: s.sectionTitle(darkMode)
  }, section.title), section.items.map(item => /*#__PURE__*/React.createElement("button", {
    key: item.id,
    style: s.item(darkMode, currentPage === item.id),
    onClick: () => onNav(item.id)
  }, item.label))))));
};
const sidebarStyles = {
  sidebar: dark => ({
    position: 'fixed',
    top: 56,
    left: 0,
    bottom: 0,
    width: 260,
    background: dark ? 'hsl(220,20%,9%)' : '#ffffff',
    borderRight: `1px solid ${dark ? 'hsl(220,12%,23%)' : '#e9ecef'}`,
    overflowY: 'auto',
    padding: '20px 0 40px'
  }),
  nav: {
    display: 'flex',
    flexDirection: 'column'
  },
  section: {
    marginBottom: 4
  },
  sectionTitle: dark => ({
    fontSize: 11,
    fontWeight: 700,
    textTransform: 'uppercase',
    letterSpacing: '0.1em',
    color: dark ? 'hsl(220,8%,40%)' : '#adb5bd',
    padding: '14px 20px 4px',
    fontFamily: 'Inter, sans-serif'
  }),
  item: (dark, active) => ({
    display: 'block',
    width: '100%',
    textAlign: 'left',
    padding: '5px 20px 5px 24px',
    border: 'none',
    cursor: 'pointer',
    fontSize: 13,
    fontFamily: 'Inter, sans-serif',
    fontWeight: active ? 500 : 400,
    color: active ? '#389826' : dark ? '#adb5bd' : '#495057',
    background: active ? dark ? 'hsl(220,14%,17%)' : '#f0faea' : 'none',
    borderLeft: active ? '2px solid #389826' : '2px solid transparent',
    transition: 'all 120ms ease'
  })
};
Object.assign(window, {
  Sidebar,
  NAV_SECTIONS
});
})(); } catch (e) { __ds_ns.__errors.push({ path: "ui_kits/docs/Sidebar.jsx", error: String((e && e.message) || e) }); }

})();
