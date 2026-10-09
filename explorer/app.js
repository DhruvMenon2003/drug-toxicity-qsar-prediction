// Drug Hansch Space: a scroll-driven 3D scatter of the curated drugs in RDKit descriptor space, with Hansch-type fits.
// Data come from explorer/build_data.py (data/drugs.json, data/hansch.json, data/molecules.json); nothing is computed
// from the network at view time except the 2D depictions, which SmilesDrawer draws from the curated SMILES.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';

const $ = (s, el = document) => el.querySelector(s);
const S = 10;                                   // half-width of the plot cube, scene units
const reduceMotion = matchMedia('(prefers-reduced-motion: reduce)');
const narrow = matchMedia('(max-width: 900px)');
const darkQuery = matchMedia('(prefers-color-scheme: dark)');
const fin = Number.isFinite;
const neg = s => String(s).replace(/^-/, '−');
const fx = (v, k = 2) => (v == null || !fin(v) ? '–' : neg(Number(v).toFixed(k)));
const cap = s => (s ? s.charAt(0).toUpperCase() + s.slice(1).toLowerCase() : '');
const esc = s => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
const median = a => { const b = [...a].sort((x, y) => x - y), m = b.length >> 1; return b.length % 2 ? b[m] : (b[m - 1] + b[m]) / 2; };
const quant = (a, q) => { const b = [...a].sort((x, y) => x - y), i = (b.length - 1) * q, lo = Math.floor(i); return b[lo] + (b[Math.ceil(i)] - b[lo]) * (i - lo); };
const pearson = (xs, ys) => {
  const n = xs.length, mx = xs.reduce((a, b) => a + b) / n, my = ys.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { const dx = xs[i] - mx, dy = ys[i] - my; sxy += dx * dy; sxx += dx * dx; syy += dy * dy; }
  return sxy / Math.sqrt(sxx * syy);
};

// ---------------------------------------------------------------- axes the viewer can choose
const AX = {
  pc1: { label: 'PC1', get: d => d.pc?.[0], dec: 1 },
  pc2: { label: 'PC2', get: d => d.pc?.[1], dec: 1 },
  pc3: { label: 'PC3', get: d => d.pc?.[2], dec: 1 },
  logp: { label: 'log P', unit: '', get: d => d.d.logp, dec: 1, long: 'log P (Crippen, RDKit)' },
  xlogp: { label: 'XLogP3', get: d => d.d.xlogp, dec: 1, long: 'XLogP3 (PubChem)' },
  ts1: { label: 'Structure map 1', get: d => d.ts?.[0], dec: 0, arbitrary: true, long: 'Structure map, t-SNE 1 (Morgan count fingerprints)' },
  ts2: { label: 'Structure map 2', get: d => d.ts?.[1], dec: 0, arbitrary: true, long: 'Structure map, t-SNE 2 (Morgan count fingerprints)' },
  mr: { label: 'MR', unit: 'cm³/mol', get: d => d.d.mr, dec: 0, long: 'Molar refractivity (cm³/mol)' },
  tpsa: { label: 'TPSA', unit: 'Å²', get: d => d.d.tpsa, dec: 0, long: 'Topological polar surface area (Å²)' },
  act: { label: 'pChEMBL', get: d => d.act?.p, dec: 1, long: 'Activity at the mechanism target (pChEMBL)' },
  mw: { label: 'MW', unit: 'Da', get: d => d.d.mw, dec: 0, long: 'Molecular weight (Da)' },
  hbd: { label: 'HBD', get: d => d.d.hbd, dec: 0, long: 'Hydrogen-bond donors (NH + OH)' },
  hba: { label: 'HBA', get: d => d.d.hba, dec: 0, long: 'Hydrogen-bond acceptors (N + O)' },
  rotb: { label: 'Rotatable bonds', get: d => d.d.rotb, dec: 0 },
  arom: { label: 'Aromatic rings', get: d => d.d.arom, dec: 0 },
  fsp3: { label: 'Fsp3', get: d => d.d.fsp3, dec: 2, long: 'Fraction of sp3 carbons' },
  cx: { label: 'Complexity', get: d => d.d.cx, dec: 0, long: 'PubChem complexity' },
  none: { label: 'Flat (no axis)', pseudo: true },
  jitter: { label: 'Random spread (no meaning)', pseudo: true, get: d => d._jit },
};
const axLong = k => AX[k].long || AX[k].label;

// ---------------------------------------------------------------- state
const app = {
  drugs: [], byCid: new Map(), meta: null, hansch: null, groups: new Map(), mols: null, molsPromise: null,
  view: { x: 'pc1', y: 'pc2', z: 'pc3', color: 'atc', atc: [null, null, null], l1: '', l2: '', l3: '', route: '', small: false, group: null, veber: false },
  scales: null, visible: [], selected: -1, hover: -1, tokens: {}, step: null, seriesId: null,
};

// ---------------------------------------------------------------- colours from the CSS tokens (both themes)
function readTokens() {
  const cs = getComputedStyle(document.documentElement), t = {};
  for (const k of ['panel', 'ink', 'ink-2', 'muted', 'line', 'grid', 'accent', 'mark', 'point', 's1', 's2', 's3', 'q1', 'q2', 'q3', 'q4', 'q5', 'surface-alpha'])
    t[k] = cs.getPropertyValue('--' + k).trim();
  app.tokens = t;
}
const isDark = () => {
  const th = document.documentElement.getAttribute('data-theme');
  return th ? th === 'dark' : darkQuery.matches;
};

// ---------------------------------------------------------------- 3D stage
const stage = { renderer: null, labels: null, scene: null, camera: null, controls: null, points: null, axes: null, extras: null, ring: null, selLabel: null };
const anim = { from: null, to: null, t: 1, tween: null };
const _m = new THREE.Matrix4(), _c = new THREE.Color();

function initStage() {
  const canvas = $('#space'), host = $('.stage-inner');
  const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  renderer.setClearColor(0x000000, 0);
  const labels = new CSS2DRenderer();
  labels.domElement.style.position = 'absolute';
  labels.domElement.style.inset = '0';
  $('#labels').appendChild(labels.domElement);
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(40, 1, 0.1, 500);
  camera.position.set(24, 14, 26);
  scene.add(new THREE.HemisphereLight(0xffffff, 0x8a9a96, 2.1));
  const sun = new THREE.DirectionalLight(0xffffff, 1.8);
  sun.position.set(12, 24, 16);
  scene.add(sun);
  const controls = new OrbitControls(camera, canvas);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;
  controls.minDistance = 8;
  controls.maxDistance = 90;
  controls.autoRotateSpeed = 0.6;
  // Plain wheel scrolls the page (the tour is scroll-driven); Ctrl/⌘ + wheel or a trackpad pinch zooms.
  host.addEventListener('wheel', e => { if (!e.ctrlKey && !e.metaKey) e.stopImmediatePropagation(); }, { capture: true });

  const n = app.drugs.length;
  const geo = new THREE.SphereGeometry(0.22, 18, 12);
  const mat = new THREE.MeshStandardMaterial({ roughness: 0.5, metalness: 0.0 });
  const points = new THREE.InstancedMesh(geo, mat, n);
  points.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
  points.frustumCulled = false;                 // instances move between layouts; keep culling and picking bounds fixed to the plot cube
  points.boundingSphere = new THREE.Sphere(new THREE.Vector3(), S * 2);
  for (let i = 0; i < n; i++) points.setColorAt(i, _c.set('#888888'));
  scene.add(points);

  const ring = new THREE.Mesh(new THREE.RingGeometry(0.36, 0.48, 40), new THREE.MeshBasicMaterial({ side: THREE.DoubleSide, depthTest: false, transparent: true }));
  ring.renderOrder = 10;
  ring.visible = false;
  scene.add(ring);
  const selDiv = document.createElement('div');
  selDiv.className = 'sel-label';
  const selLabel = new CSS2DObject(selDiv);
  selLabel.center.set(0, 1.4);
  selLabel.visible = false;
  scene.add(selLabel);

  Object.assign(stage, { renderer, labels, scene, camera, controls, points, ring, selLabel });
  anim.from = new Float32Array(n * 4);
  anim.to = new Float32Array(n * 4);
  anim.cur = new Float32Array(n * 4);

  new ResizeObserver(resize).observe(host);
  resize();
  bindPointer(canvas, host);
  renderer.setAnimationLoop(frame);
}

function resize() {
  const host = $('.stage-inner'), w = host.clientWidth, h = host.clientHeight;
  if (!w || !h) return;
  stage.renderer.setSize(w, h, false);
  stage.labels.setSize(w, h);
  stage.camera.aspect = w / h;
  stage.camera.updateProjectionMatrix();
}

function frame() {
  const { points, controls, camera, ring, selLabel } = stage;
  controls.update();
  const t = anim.t, n = app.drugs.length, F = anim.from, T = anim.to, C = anim.cur;
  for (let i = 0; i < n; i++) {
    const k = i * 4;
    for (let j = 0; j < 4; j++) C[k + j] = F[k + j] + (T[k + j] - F[k + j]) * t;
    let s = C[k + 3];
    if (i === app.hover && s > 0) s *= 1.7;
    else if (i === app.selected && s > 0) s *= 1.35;
    _m.makeScale(s, s, s).setPosition(C[k], C[k + 1], C[k + 2]);
    points.setMatrixAt(i, _m);
  }
  points.instanceMatrix.needsUpdate = true;
  if (app.selected >= 0 && C[app.selected * 4 + 3] > 0.5) {
    const k = app.selected * 4;
    ring.position.set(C[k], C[k + 1], C[k + 2]);
    ring.quaternion.copy(camera.quaternion);
    ring.visible = selLabel.visible = true;
    selLabel.position.copy(ring.position);
  } else ring.visible = selLabel.visible = false;
  stage.renderer.render(stage.scene, camera);
  stage.labels.render(stage.scene, camera);
  keepTitlesInside();
}

// Axis titles hang off the axis ends; slide any that would run past the stage edge back inside it, and below the legend
// or the buttons if they would sit on them. The CSS translate property composes with the transform CSS2DRenderer
// writes every frame, so the two do not fight.
function keepTitlesInside() {
  const host = stage.labels.domElement, box = host.getBoundingClientRect();
  const keepOut = ['.hud', '.stage-tools'].map(q => $(q).getBoundingClientRect()).filter(r => r.width && r.height);
  for (const el of host.querySelectorAll('.ax-title')) {
    const r = el.getBoundingClientRect();
    if (!r.width) continue;
    const [wx, wy] = (el.dataset.shift || '0 0').split(' ').map(Number);
    const left = r.left - wx, right = r.right - wx, top = r.top - wy, bottom = r.bottom - wy;
    const x = Math.round(Math.max(box.left + 8 - left, Math.min(0, box.right - 8 - right)));
    let y = 0;
    for (const k of keepOut) if (left + x < k.right && right + x > k.left && top + y < k.bottom && bottom + y > k.top) y = Math.round(k.bottom + 6 - top);
    const shift = `${x} ${y}`;
    if (shift !== (el.dataset.shift || '0 0')) { el.dataset.shift = shift; el.style.translate = x || y ? `${x}px ${y}px` : ''; }
  }
}

// ---------------------------------------------------------------- scales and layout
function niceStep(range, count) {
  const raw = range / Math.max(1, count - 1), mag = 10 ** Math.floor(Math.log10(raw)), f = raw / mag;
  return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * mag;
}
function makeScale(key, list) {
  if (key === 'none') return { pseudo: true, map: () => -S };
  if (key === 'jitter') return { pseudo: true, map: v => v * S * 0.3 };
  const vals = list.map(AX[key].get).filter(fin);
  let lo = vals.length ? Math.min(...vals) : 0, hi = vals.length ? Math.max(...vals) : 1;
  if (lo === hi) { lo -= 1; hi += 1; }
  const step = niceStep(hi - lo, 5);
  lo = Math.floor(lo / step + 1e-9) * step;
  hi = Math.ceil(hi / step - 1e-9) * step;
  const ticks = [];
  for (let v = lo; v <= hi + step * 1e-6; v += step) ticks.push(+v.toFixed(10));
  return { lo, hi, step, ticks, map: v => ((v - lo) / (hi - lo)) * 2 * S - S };
}

function passes(d, axes = true) {
  const v = app.view;
  if (v.l1 && !d.atc.some(a => a.code.startsWith(v.l1))) return false;
  if (v.l2 && !d.atc.some(a => a.code.startsWith(v.l2))) return false;
  if (v.l3 && !d.atc.some(a => a.code.startsWith(v.l3))) return false;
  if (v.route && !d.routes.includes(v.route)) return false;
  if (v.small && !(d.d.mw <= app.hansch.max_mw)) return false;
  if (v.group && !app.groups.get(v.group)?.cids.has(d.cid)) return false;
  if (axes) for (const k of [v.x, v.y, v.z]) if (!AX[k].pseudo && !fin(AX[k].get(d))) return false;
  return true;
}

function layout({ instant = false } = {}) {
  const v = app.view, list = app.drugs.filter(d => passes(d));
  app.visible = list;
  const vis = new Set(list);
  const sx = makeScale(v.x, list), sy = makeScale(v.y, list), sz = makeScale(v.z, list);
  app.scales = { x: sx, y: sy, z: sz };
  // start from wherever the points are now, so a change in mid-flight does not jump
  anim.from.set(anim.cur);
  app.drugs.forEach((d, i) => {
    const k = i * 4;
    if (vis.has(d)) {
      anim.to[k] = sx.map(AX[v.x].get?.(d) ?? 0);
      anim.to[k + 1] = sy.map(AX[v.y].get?.(d) ?? 0);
      anim.to[k + 2] = sz.map(AX[v.z].get?.(d) ?? 0);
      anim.to[k + 3] = 1;
      if (anim.from[k + 3] === 0) { anim.from[k] = anim.to[k]; anim.from[k + 1] = anim.to[k + 1]; anim.from[k + 2] = anim.to[k + 2]; }
    } else {
      anim.to[k] = anim.from[k]; anim.to[k + 1] = anim.from[k + 1]; anim.to[k + 2] = anim.from[k + 2];
      anim.to[k + 3] = 0;
    }
  });
  anim.tween?.kill();
  if (instant || reduceMotion.matches || !window.gsap) { anim.t = 1; anim.from.set(anim.to); anim.cur.set(anim.to); }
  else { anim.t = 0; anim.tween = gsap.to(anim, { t: 1, duration: 1.1, ease: 'power2.inOut' }); }
  buildAxes();
  buildExtras();
  recolor();
  renderAxisKey();
  renderCount();
  if ($('#table-details').open) renderTable();
}

// ---------------------------------------------------------------- axes, grid and tick labels
function disposeGroup(g) {
  if (!g) return;
  g.parent?.remove(g);
  // CSS2D labels only drop their DOM element when they themselves are removed, not when an ancestor group is
  g.traverse(o => { o.geometry?.dispose(); o.material?.dispose?.(); if (o.isCSS2DObject) o.element.remove(); });
}
function label(text, cls, pos, center = [0.5, 0.5]) {
  const el = document.createElement('div');
  el.className = cls;
  el.textContent = text;
  const o = new CSS2DObject(el);
  o.position.copy(pos);
  o.center.set(...center);
  return o;
}
const tickText = (key, v) => neg(Math.abs(v) >= 1000 ? v.toFixed(0) : +v.toPrecision(4) + '');

function buildAxes() {
  disposeGroup(stage.axes);
  const g = new THREE.Group(), { x: sx, y: sy, z: sz } = app.scales, v = app.view, T = app.tokens;
  const grid = [], frame = [];
  const seg = (arr, a, b) => arr.push(a[0], a[1], a[2], b[0], b[1], b[2]);
  // floor (y = −S): x ticks run along z, z ticks run along x
  if (!sx.pseudo) for (const t of sx.ticks) seg(grid, [sx.map(t), -S, -S], [sx.map(t), -S, S]);
  if (!sz.pseudo) for (const t of sz.ticks) seg(grid, [-S, -S, sz.map(t)], [S, -S, sz.map(t)]);
  if (!sy.pseudo) {
    for (const t of sy.ticks) { seg(grid, [-S, sy.map(t), -S], [S, sy.map(t), -S]); seg(grid, [-S, sy.map(t), -S], [-S, sy.map(t), S]); }
    if (!sx.pseudo) for (const t of sx.ticks) seg(grid, [sx.map(t), -S, -S], [sx.map(t), S, -S]);
    if (!sz.pseudo) for (const t of sz.ticks) seg(grid, [-S, -S, sz.map(t)], [-S, S, sz.map(t)]);
  }
  seg(frame, [-S, -S, S], [S, -S, S]);
  seg(frame, [S, -S, -S], [S, -S, S]);
  if (!sy.pseudo) seg(frame, [-S, -S, S], [-S, S, S]);
  const mk = (arr, color, opacity) => {
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.Float32BufferAttribute(arr, 3));
    return new THREE.LineSegments(geo, new THREE.LineBasicMaterial({ color, transparent: true, opacity }));
  };
  g.add(mk(grid, T.line, 0.75), mk(frame, T['ink-2'], 0.9));
  const P = (x, y, z) => new THREE.Vector3(x, y, z);
  // t-SNE units mean nothing, so those axes get no numbers and a short name (the axis key carries the long one);
  // phone-width stages get the short names throughout
  const narrow = $('.stage-inner').clientWidth < 600, title = k => narrow || AX[k].arbitrary ? AX[k].label + (AX[k].unit ? ` (${AX[k].unit})` : '') : axLong(k), ticks = (k, sc) => AX[k].arbitrary ? [] : sc.ticks;
  if (!sx.pseudo) {
    for (const t of ticks(v.x, sx)) g.add(label(tickText(v.x, t), 'ax-label', P(sx.map(t), -S, S + 0.9)));
    g.add(label(title(v.x), 'ax-title', P(0, -S, S + 2.4)));
  }
  if (!sz.pseudo) {
    for (const t of ticks(v.z, sz)) g.add(label(tickText(v.z, t), 'ax-label', P(S + 0.9, -S, sz.map(t)), [0, 0.5]));
    g.add(label(title(v.z), 'ax-title', P(S + 2.2, -S, 0), [0, 0.5]));
  }
  if (!sy.pseudo) {
    for (const t of ticks(v.y, sy)) g.add(label(tickText(v.y, t), 'ax-label', P(-S - 0.6, sy.map(t), S), [1, 0.5]));
    g.add(label(title(v.y), 'ax-title', P(-S + 0.3, S, S), [0, 0]));   // beside the top of the axis, clear of the legend and the canvas edge
  }
  stage.scene.add(g);
  stage.axes = g;
}

// ---------------------------------------------------------------- Veber plane, fitted surface, residual sticks
function buildExtras() {
  disposeGroup(stage.extras);
  const g = new THREE.Group(), v = app.view, T = app.tokens, sc = app.scales;
  const fade = (mat, to) => {
    mat.transparent = true;
    if (reduceMotion.matches || !window.gsap) { mat.opacity = to; return; }
    mat.opacity = 0;
    gsap.to(mat, { opacity: to, duration: 0.8, delay: 0.6 });
  };
  // Veber et al. 2002: oral bioavailability falls above 140 Å² of polar surface
  const axisOf = ['x', 'y', 'z'].find(a => v[a] === 'tpsa');
  if (v.veber && axisOf && 140 >= sc[axisOf].lo && 140 <= sc[axisOf].hi) {
    const at = sc[axisOf].map(140);
    const plane = new THREE.Mesh(new THREE.PlaneGeometry(2 * S, 2 * S), new THREE.MeshBasicMaterial({ color: T.mark, side: THREE.DoubleSide, depthWrite: false }));
    if (axisOf === 'y') { plane.rotation.x = -Math.PI / 2; plane.position.y = at; }
    else if (axisOf === 'x') { plane.rotation.y = Math.PI / 2; plane.position.x = at; }
    else plane.position.z = at;
    fade(plane.material, 0.13);
    const edge = new THREE.LineSegments(new THREE.EdgesGeometry(plane.geometry), new THREE.LineBasicMaterial({ color: T.mark }));
    edge.rotation.copy(plane.rotation); edge.position.copy(plane.position);
    fade(edge.material, 0.9);
    const lp = axisOf === 'y' ? new THREE.Vector3(S, at, -S) : axisOf === 'x' ? new THREE.Vector3(at, S, -S) : new THREE.Vector3(S, S, at);
    g.add(plane, edge, label('TPSA 140 Å² (Veber)', 'ax-plane', lp, [1, 1.2]));
  }
  const grp = v.group && app.groups.get(v.group);
  if (grp && v.x === 'logp' && v.z === 'mr' && v.y === 'act') {
    const h = grp.fit, members = app.visible.filter(d => grp.cids.has(d.cid));
    const lps = members.map(d => d.d.logp), mrs = members.map(d => d.d.mr);
    const tpsaMean = members.reduce((a, d) => a + (d.d.tpsa ?? 0), 0) / Math.max(1, members.length);
    const f = (lp, mr, tp) => h.coef.const + (h.coef.logP ?? 0) * lp + (h.coef['logP²'] ?? 0) * lp * lp + (h.coef.MR ?? 0) * mr + (h.coef.TPSA ?? 0) * tp;
    const [l0, l1] = [Math.min(...lps), Math.max(...lps)], [m0, m1] = [Math.min(...mrs), Math.max(...mrs)];
    const N = 28, geo = new THREE.PlaneGeometry(1, 1, N, N), pos = geo.attributes.position;
    for (let i = 0; i < pos.count; i++) {
      const u = pos.getX(i) + 0.5, w = pos.getY(i) + 0.5, lp = l0 + (l1 - l0) * u, mr = m0 + (m1 - m0) * w;
      pos.setXYZ(i, sc.x.map(lp), THREE.MathUtils.clamp(sc.y.map(f(lp, mr, tpsaMean)), -S, S), sc.z.map(mr));
    }
    geo.computeVertexNormals();
    const surf = new THREE.Mesh(geo, new THREE.MeshBasicMaterial({ color: T.accent, side: THREE.DoubleSide, depthWrite: false }));
    fade(surf.material, parseFloat(T['surface-alpha']) || 0.3);
    const wire = new THREE.Mesh(geo, new THREE.MeshBasicMaterial({ color: T.accent, wireframe: true }));
    fade(wire.material, 0.35);
    // residual sticks: observed pChEMBL to the model's fitted value for that drug
    const fitBy = new Map(h.points.map(p => [p.cid, p.fit])), sticks = [];
    for (const d of members) {
      const x = sc.x.map(d.d.logp), z = sc.z.map(d.d.mr);
      sticks.push(x, sc.y.map(d.act.p), z, x, THREE.MathUtils.clamp(sc.y.map(fitBy.get(d.cid)), -S, S), z);
    }
    const sg = new THREE.BufferGeometry();
    sg.setAttribute('position', new THREE.Float32BufferAttribute(sticks, 3));
    const stickLines = new THREE.LineSegments(sg, new THREE.LineBasicMaterial({ color: T.ink }));
    fade(stickLines.material, 0.55);
    g.add(surf, wire, stickLines);
  }
  stage.scene.add(g);
  stage.extras = g;
}

// ---------------------------------------------------------------- colour modes (at most three categorical hues)
const ACT_EDGES = [6, 7, 8, 9];
function l1sOf(d) { return [...new Set(d.atc.map(a => a.code[0]))]; }
function colorOf(d) {
  const v = app.view, T = app.tokens;
  switch (v.color) {
    case 'atc': {
      const mine = l1sOf(d), slot = v.atc.findIndex(c => c && mine.includes(c));
      return slot >= 0 ? T['s' + (slot + 1)] : T.point;
    }
    case 'oral': return !d.routes.length ? T.point : d.routes.includes('ORAL') ? T.s1 : T.s2;
    case 'act': {
      if (!d.act) return T.point;
      return T['q' + (1 + ACT_EDGES.filter(e => d.act.p >= e).length)];
    }
    case 'ro5': return d.d.ro5 == null ? T.point : T[['q2', 'q4', 'q5'][Math.min(2, d.d.ro5)]];
    default: return T.s1;
  }
}
function recolor() {
  app.drugs.forEach((d, i) => stage.points.setColorAt(i, _c.set(colorOf(d) || '#888')));
  stage.points.instanceColor.needsUpdate = true;
  stage.ring.material.color.set(app.tokens.mark);
  renderLegend();
}

function renderLegend() {
  const v = app.view, T = app.tokens, sw = (c, t) => `<span><span class="sw" style="background:${c}"></span>${esc(t)}</span>`;
  let h = '';
  if (v.color === 'atc') {
    v.atc.forEach((c, i) => { if (c) h += sw(T['s' + (i + 1)], `${c} ${cap(app.meta.atc_names[c])}`); });
    h += sw(T.point, 'Other ATC groups');
  } else if (v.color === 'oral') h = sw(T.s1, 'Oral route') + sw(T.s2, 'No oral route') + sw(T.point, 'No route listed');
  else if (v.color === 'act') h = `<span>pChEMBL &lt;6<span class="ramp"></span>≥9</span>` + sw(T.point, 'No activity value');
  else if (v.color === 'ro5') h = sw(T.q2, '0 rule-of-five violations') + sw(T.q4, '1') + sw(T.q5, '2 or more');
  else h = `<span>Each sphere is one drug</span>`;
  $('#legend').innerHTML = h;
  // chips mirror the slots
  document.querySelectorAll('#atc-pick .chip').forEach(ch => {
    const slot = v.atc.indexOf(ch.dataset.code);
    ch.setAttribute('aria-pressed', slot >= 0);
    ch.querySelector('.sw').style.background = slot >= 0 ? T['s' + (slot + 1)] : '';
    ch.disabled = slot < 0 && !v.atc.includes(null);
  });
  $('#atc-pick-wrap').hidden = v.color !== 'atc';
}

function renderAxisKey() {
  const v = app.view;
  $('#axis-key').innerHTML = ['x', 'y', 'z'].filter(a => v[a] !== 'none')
    .map(a => `<span><b>${a.toUpperCase()}</b> ${esc(axLong(v[a]))}</span>`).join('');
}
function renderCount() {
  const v = app.view, n = app.visible.length, total = app.drugs.length;
  const grp = v.group && app.groups.get(v.group);
  const big = d => v.small && !(d.d.mw <= app.hansch.max_mw);
  const missing = app.drugs.filter(d => !(v.group && !grp?.cids.has(d.cid)) && !big(d) && [v.x, v.y, v.z].some(k => !AX[k].pseudo && !fin(AX[k].get(d)))).length;
  let nBig = 0;                                  // drugs the size filter alone removes from this view
  if (v.small && !v.group) { v.small = false; nBig = app.drugs.filter(d => !(d.d.mw <= app.hansch.max_mw) && passes(d, false)).length; v.small = true; }
  let s = `Showing ${n} of ${total} drugs`;
  if (grp) s += ` in the fit for ${grp.fit.id === 'ALL' ? 'all drugs with an activity value' : grp.fit.label}`;
  if (nBig) s += `; ${nBig} drugs over ${app.hansch.max_mw} Da are left out`;
  if (missing) s += `; ${missing} have no value on one of the axes`;
  $('#count').textContent = s + '.';
}

// ---------------------------------------------------------------- hover, click, tooltip
const ray = new THREE.Raycaster(), ndc = new THREE.Vector2();
function pick(ev, canvas) {
  const r = canvas.getBoundingClientRect();
  ndc.set(((ev.clientX - r.left) / r.width) * 2 - 1, -((ev.clientY - r.top) / r.height) * 2 + 1);
  ray.setFromCamera(ndc, stage.camera);
  const hits = ray.intersectObject(stage.points);
  for (const h of hits) if (anim.cur[h.instanceId * 4 + 3] > 0.5) return h.instanceId;
  return -1;
}
function bindPointer(canvas, host) {
  const tip = $('#tooltip');
  let down = null, pending = null;
  canvas.addEventListener('pointermove', ev => {
    if (pending) return;
    pending = requestAnimationFrame(() => {
      pending = null;
      const i = ev.buttons ? -1 : pick(ev, canvas);
      app.hover = i;
      canvas.style.cursor = i >= 0 ? 'pointer' : '';
      if (i < 0) { tip.hidden = true; return; }
      const d = app.drugs[i], v = app.view, hr = host.getBoundingClientRect();
      const rows = ['x', 'y', 'z'].filter(a => !AX[v[a]].pseudo)
        .map(a => `<div class="t-row">${esc(AX[v[a]].label)} ${fx(AX[v[a]].get(d), AX[v[a]].dec)}${AX[v[a]].unit ? ' ' + AX[v[a]].unit : ''}</div>`).join('');
      tip.innerHTML = `<div class="t-name">${esc(cap(d.name))}</div>${rows}<div class="t-row">ATC ${esc(d.atc.map(a => a.code).slice(0, 3).join(', '))}${d.atc.length > 3 ? '…' : ''}</div>`;
      tip.hidden = false;
      const x = ev.clientX - hr.left, y = ev.clientY - hr.top, w = tip.offsetWidth, h = tip.offsetHeight;
      tip.style.left = Math.min(x + 14, hr.width - w - 8) + 'px';
      tip.style.top = (y + h + 18 > hr.height ? y - h - 12 : y + 16) + 'px';
    });
  });
  canvas.addEventListener('pointerleave', () => { app.hover = -1; tip.hidden = true; });
  canvas.addEventListener('pointerdown', ev => { down = [ev.clientX, ev.clientY]; });
  canvas.addEventListener('pointerup', ev => {
    if (!down || Math.hypot(ev.clientX - down[0], ev.clientY - down[1]) > 5) return;
    const i = pick(ev, canvas);
    if (i >= 0) selectDrug(app.drugs[i].cid);
  });
}

// ---------------------------------------------------------------- camera moves
function camTo(pos, { instant = false } = {}) {
  const c = stage.camera, k = Math.max(1, 1.2 / Math.max(0.3, c.aspect));
  const p = new THREE.Vector3(...pos).multiplyScalar(k);
  if (instant || reduceMotion.matches || !window.gsap) { c.position.copy(p); stage.controls.target.set(0, 0, 0); return; }
  gsap.killTweensOf(c.position);
  gsap.killTweensOf(stage.controls.target);
  gsap.to(c.position, { x: p.x, y: p.y, z: p.z, duration: 1.4, ease: 'power2.inOut' });
  gsap.to(stage.controls.target, { x: 0, y: 0, z: 0, duration: 1.4, ease: 'power2.inOut' });
}

// ---------------------------------------------------------------- the guided tour
const STEPS = {
  space: { view: { x: 'pc1', y: 'pc2', z: 'pc3', color: 'atc' }, cam: [24, 14, 27] },
  structure: { view: { x: 'ts1', y: 'none', z: 'ts2', color: 'atc' }, cam: [0, 37, 3] },
  logp: { view: { x: 'logp', y: 'none', z: 'jitter', color: 'single', small: true }, cam: [0, 22, 26] },
  mr: { view: { x: 'logp', y: 'none', z: 'mr', color: 'single', small: true }, cam: [0, 37, 3] },
  tpsa: { view: { x: 'logp', y: 'tpsa', z: 'mr', color: 'oral', veber: true, small: true }, cam: [27, 12, 24] },
  activity: { view: { x: 'logp', y: 'act', z: 'mr', color: 'act', small: true }, cam: [27, 11, 24] },
  series: { view: { x: 'logp', y: 'act', z: 'mr', color: 'single', group: 'SERIES', small: true }, cam: [25, 12, 26] },
  all: { view: { x: 'logp', y: 'act', z: 'mr', color: 'single', group: 'ALL', small: true }, cam: [27, 12, 23] },
};
function goStep(name) {
  if (app.step === name) return;
  app.step = name;
  document.querySelectorAll('.step').forEach(s => s.classList.toggle('is-active', s.dataset.step === name));
  const st = STEPS[name], gv = st.view.group === 'SERIES' ? app.seriesId : st.view.group ?? null;
  setView({ l1: '', l2: '', l3: '', route: '', small: false, veber: false, ...st.view, group: gv });
  camTo(st.cam);
}
function initTour() {
  if (!window.ScrollTrigger) { goStep('space'); return; }
  gsap.registerPlugin(ScrollTrigger);
  document.querySelectorAll('.step').forEach(el => {
    ScrollTrigger.create({
      trigger: el, start: () => (narrow.matches ? 'top 66%' : 'top 62%'), end: () => (narrow.matches ? 'bottom 66%' : 'bottom 62%'),
      onToggle: self => { if (self.isActive) goStep(el.dataset.step); },
    });
  });
  if (!app.step) goStep('space');
}

// ---------------------------------------------------------------- view changes from the tools
function setView(patch, opts) {
  Object.assign(app.view, patch);
  syncControls();
  layout(opts);
}
function syncControls() {
  const v = app.view;
  for (const a of ['x', 'y', 'z']) $('#ax-' + a).value = v[a];
  $('#color-by').value = v.color;
  fillAtcSelects();
  $('#f-route').value = v.route;
  $('#f-small').checked = v.small;
  $('#h-group').value = v.group && app.groups.has(v.group) ? v.group : $('#h-group').value;
}

function initControls() {
  for (const a of ['x', 'y', 'z']) {
    const sel = $('#ax-' + a);
    sel.innerHTML = Object.entries(AX).filter(([k]) => a === 'x' ? !AX[k].pseudo : true)
      .map(([k, o]) => `<option value="${k}">${esc(o.long || o.label)}</option>`).join('');
    sel.addEventListener('change', () => setView({ [a]: sel.value, group: null, veber: false }));
  }
  $('#color-by').addEventListener('change', e => setView({ color: e.target.value }));
  // ATC chips: L1 groups by drug count; colours follow the slot a group was given, never its rank
  const count = new Map();
  for (const d of app.drugs) for (const c of l1sOf(d)) count.set(c, (count.get(c) || 0) + 1);
  const l1 = [...count.keys()].sort();
  $('#atc-pick').innerHTML = l1.map(c => `<button type="button" class="chip" data-code="${c}" aria-pressed="false" title="${esc(cap(app.meta.atc_names[c]))}"><span class="sw"></span>${c} ${esc(cap(app.meta.atc_names[c]))} <span class="mono">${count.get(c)}</span></button>`).join('');
  const top3 = [...count.entries()].sort((a, b) => b[1] - a[1]).slice(0, 3).map(e => e[0]);
  app.view.atc = top3;
  $('#atc-pick').addEventListener('click', e => {
    const ch = e.target.closest('.chip');
    if (!ch) return;
    const atc = [...app.view.atc], slot = atc.indexOf(ch.dataset.code);
    if (slot >= 0) atc[slot] = null;
    else { const free = atc.indexOf(null); if (free < 0) return; atc[free] = ch.dataset.code; }
    app.view.atc = atc;
    recolor();
  });
  for (const id of ['f-l1', 'f-l2', 'f-l3']) {
    $('#' + id).addEventListener('change', e => {
      const lvl = id.slice(-2), patch = { [lvl]: e.target.value };
      if (lvl === 'l1') Object.assign(patch, { l2: '', l3: '' });
      if (lvl === 'l2') patch.l3 = '';
      setView(patch);
    });
  }
  const routes = [...new Set(app.drugs.flatMap(d => d.routes))].sort();
  $('#f-route').innerHTML = '<option value="">All routes</option>' + routes.map(r => `<option value="${esc(r)}">${esc(cap(r))}</option>`).join('');
  $('#f-route').addEventListener('change', e => setView({ route: e.target.value }));
  $('#f-small').addEventListener('change', e => setView({ small: e.target.checked, group: null }));
  // search
  $('#drug-names').innerHTML = [...app.drugs].sort((a, b) => a.name.localeCompare(b.name)).map(d => `<option value="${esc(cap(d.name))}">`).join('');
  const find = () => {
    const q = $('#search').value.trim().toLowerCase();
    if (!q) return;
    const d = app.drugs.find(x => x.name.toLowerCase() === q) || app.drugs.find(x => x.name.toLowerCase().startsWith(q));
    if (d) {
      if (!passes(d)) setView({ l1: '', l2: '', l3: '', route: '', small: false, group: null });
      selectDrug(d.cid);
    }
  };
  $('#search').addEventListener('change', find);
  $('#search').addEventListener('keydown', e => { if (e.key === 'Enter') { e.preventDefault(); find(); } });
  $('#btn-reset').addEventListener('click', () => camTo(STEPS[app.step]?.cam || STEPS.space.cam));
  $('#btn-spin').addEventListener('click', e => {
    stage.controls.autoRotate = !stage.controls.autoRotate;
    e.currentTarget.setAttribute('aria-pressed', stage.controls.autoRotate);
  });
  $('#table-details').addEventListener('toggle', () => { if ($('#table-details').open) renderTable(); });
  $('#show-h').addEventListener('change', () => { if (app.selected >= 0) showMolecule(app.drugs[app.selected]); });
}

function fillAtcSelects() {
  const v = app.view, names = app.meta.atc_names, codes = new Set();
  for (const d of app.drugs) for (const a of d.atc) { codes.add(a.code[0]); codes.add(a.code.slice(0, 3)); codes.add(a.code.slice(0, 4)); }
  const opt = (c, sel) => `<option value="${c}"${c === sel ? ' selected' : ''}>${c} ${esc(cap(names[c] || ''))}</option>`;
  const all = [...codes].sort();
  $('#f-l1').innerHTML = '<option value="">All groups</option>' + all.filter(c => c.length === 1).map(c => opt(c, v.l1)).join('');
  $('#f-l2').innerHTML = '<option value="">All</option>' + all.filter(c => c.length === 3 && (!v.l1 || c.startsWith(v.l1))).map(c => opt(c, v.l2)).join('');
  $('#f-l3').innerHTML = '<option value="">All</option>' + all.filter(c => c.length === 4 && (!v.l2 || c.startsWith(v.l2)) && (!v.l1 || c.startsWith(v.l1))).map(c => opt(c, v.l3)).join('');
}

// ---------------------------------------------------------------- data table (the accessible view of the scatter)
function renderTable() {
  const v = app.view, ax = ['x', 'y', 'z'].filter(a => !AX[v[a]].pseudo);
  $('#table thead').innerHTML = `<tr><th scope="col">Drug</th><th scope="col">ATC</th>${ax.map(a => `<th scope="col">${esc(AX[v[a]].label)}</th>`).join('')}<th scope="col">pChEMBL</th></tr>`;
  const rows = [...app.visible].sort((a, b) => a.name.localeCompare(b.name));
  $('#table tbody').innerHTML = rows.map(d => `<tr data-cid="${d.cid}" tabindex="0" aria-selected="${app.drugs[app.selected]?.cid === d.cid}"><td>${esc(cap(d.name))}</td><td>${esc(d.atc[0]?.code || '')}</td>${ax.map(a => `<td>${fx(AX[v[a]].get(d), AX[v[a]].dec)}</td>`).join('')}<td>${d.act ? fx(d.act.p, 1) : '–'}</td></tr>`).join('');
}
$('#table').addEventListener('click', e => { const tr = e.target.closest('tr[data-cid]'); if (tr) selectDrug(+tr.dataset.cid); });
$('#table').addEventListener('keydown', e => { const tr = e.target.closest('tr[data-cid]'); if (tr && (e.key === 'Enter' || e.key === ' ')) { e.preventDefault(); selectDrug(+tr.dataset.cid); } });

// ---------------------------------------------------------------- drug detail
function selectDrug(cid) {
  const i = app.drugs.findIndex(d => d.cid === cid);
  if (i < 0) return;
  app.selected = i;
  const d = app.drugs[i];
  stage.selLabel.element.textContent = cap(d.name);
  document.querySelectorAll('#table tbody tr').forEach(tr => tr.setAttribute('aria-selected', +tr.dataset.cid === cid));
  $('#detail-empty').hidden = true;
  $('#detail-body').hidden = false;
  $('#d-name').textContent = d.name;
  const pc = `https://pubchem.ncbi.nlm.nih.gov/compound/${d.pcid}`;
  $('#d-ids').innerHTML = `<a href="${pc}" target="_blank" rel="noopener">PubChem CID ${d.pcid}</a>` +
    (d.pcid !== d.cid ? ` (dataset CID ${d.cid})` : '') +
    (d.chembl ? ` · <a href="https://www.ebi.ac.uk/chembl/compound_report_card/${esc(d.chembl)}/" target="_blank" rel="noopener">${esc(d.chembl)}</a>` : '') +
    `<br>${esc(d.ik || '')}${d.d.formula ? ' · ' + esc(d.d.formula) : ''}`;
  const D = d.d, kv = (k, v, small = '') => `<div><dt>${k}</dt><dd>${v}${small ? ` <small>${small}</small>` : ''}</dd></div>`;
  $('#d-desc').innerHTML = [
    kv('log P', fx(D.logp, 2), fin(D.xlogp) ? `Crippen; XLogP3 ${fx(D.xlogp, 2)}` : 'Crippen'),
    kv('MR', fx(D.mr, 1), 'cm³/mol'), kv('TPSA', fx(D.tpsa, 1), 'Å²'), kv('MW', fx(D.mw, 1), 'Da'),
    kv('H-bond donors', D.hbd ?? '–'), kv('H-bond acceptors', D.hba ?? '–'), kv('Rotatable bonds', D.rotb ?? '–'),
    kv('Aromatic rings', D.arom ?? '–'), kv('Fsp3', fx(D.fsp3, 2)), kv('Complexity', fx(D.cx, 0), 'PubChem'),
    kv('Rule-of-five violations', D.ro5 ?? '–'), kv('Formal charge', D.charge ?? '–'),
  ].join('');
  const a = d.act;
  $('#d-act').innerHTML = a
    ? `<div class="act-box"><span class="big">${fx(a.p, 2)}</span> pChEMBL <small>(about ${fmtConc(a.p)})</small><br>` +
      `at <a href="https://www.ebi.ac.uk/chembl/target_report_card/${esc(a.target)}/" target="_blank" rel="noopener">${esc(a.target)}</a> ${esc(a.tname || '')}` +
      `<br><small>Median of ${esc(a.types || 'IC50/Ki/EC50/Kd')} values from ${a.ndoc ?? '≥2'} ChEMBL documents.</small></div>`
    : `<div class="act-box">No activity value: fewer than two ChEMBL documents agreed on a pChEMBL at the mechanism target.</div>`;
  const names = app.meta.atc_names;
  $('#d-atc').innerHTML = d.atc.map(e => {
    const path = [e.code[0], e.code.slice(0, 3), e.code.slice(0, 4), e.code.slice(0, 5)].filter(c => names[c]).map(c => cap(names[c])).join(' › ');
    return `<li><span class="atc-code">${esc(e.code)}</span> <span class="atc-path">${esc(path)}</span><br>${e.routes.length ? e.routes.map(r => `<span class="route">${esc(cap(r))}</span>`).join('') : '<span class="route">No route listed</span>'}</li>`;
  }).join('');
  $('#d-flags').textContent = d.flags.length ? 'Curation notes: ' + d.flags.join('; ') + '.' : '';
  $('#d-src').textContent = d.conf ? d.conf : 'No 3D conformer';
  drawDepiction(d);
  showMolecule(d);
}
function fmtConc(p) {
  const nM = 10 ** (9 - p);
  return nM >= 1000 ? `${+(nM / 1000).toPrecision(2)} µM` : nM >= 1 ? `${+nM.toPrecision(2)} nM` : `${+(nM * 1000).toPrecision(2)} pM`;
}

function drawDepiction(d) {
  const host = $('#d-2d');
  host.innerHTML = '';
  if (!window.SmilesDrawer || !d.smiles) return;
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('role', 'img');
  svg.setAttribute('aria-label', `2D structure of ${d.name}`);
  host.appendChild(svg);
  const drawer = new SmilesDrawer.SvgDrawer({ width: 360, height: 240, bondThickness: 1.1, padding: 12, compactDrawing: false });
  SmilesDrawer.parse(d.smiles, tree => {
    try { drawer.draw(tree, svg, isDark() ? 'dark' : 'light'); svg.removeAttribute('style'); }
    catch { host.textContent = '2D depiction not available for this structure.'; }
  }, () => { host.textContent = '2D depiction not available for this structure.'; });
}

// ---------------------------------------------------------------- ball-and-stick conformer viewer
const CPK = { 1: '#e8e8e8', 3: '#cc80ff', 5: '#ffb5b5', 6: '#909090', 7: '#3050f8', 8: '#ff0d0d', 9: '#90e050', 11: '#ab5cf2', 12: '#8aff00',
  13: '#bfa6a6', 14: '#f0c8a0', 15: '#ff8000', 16: '#e6c629', 17: '#1ff01f', 19: '#8f40d4', 20: '#3dff00', 26: '#e06633', 27: '#f090a0',
  29: '#c88033', 30: '#7d80b0', 31: '#c28f8f', 33: '#bd80e3', 34: '#ffa100', 35: '#a62929', 43: '#3b9e9e', 47: '#c0c0c0', 51: '#9e63b5',
  53: '#940094', 57: '#70d4ff', 64: '#45ffc7', 78: '#d0d0e0', 79: '#ffd123', 80: '#b8b8d0', 83: '#9e4fb5' };
const molv = { renderer: null, scene: null, camera: null, controls: null, group: null };
function initMolViewer() {
  const canvas = $('#mol');
  molv.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
  molv.renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  molv.scene = new THREE.Scene();
  molv.scene.add(new THREE.HemisphereLight(0xffffff, 0x666f6d, 2.2));
  const sun = new THREE.DirectionalLight(0xffffff, 1.6);
  sun.position.set(5, 8, 10);
  molv.camera = new THREE.PerspectiveCamera(35, 1, 0.1, 500);
  molv.camera.add(sun);
  molv.scene.add(molv.camera);
  molv.controls = new OrbitControls(molv.camera, canvas);
  molv.controls.enableDamping = true;
  molv.controls.enableZoom = false;              // the page scrolls over this canvas
  molv.controls.enablePan = false;
  molv.controls.autoRotate = !reduceMotion.matches;
  molv.controls.autoRotateSpeed = 1.2;
  const fit = () => {
    const w = canvas.clientWidth, h = canvas.clientHeight;
    if (!w || !h) return;
    molv.renderer.setSize(w, h, false);
    molv.camera.aspect = w / h;
    molv.camera.updateProjectionMatrix();
  };
  new ResizeObserver(fit).observe(canvas);
  molv.renderer.setAnimationLoop(() => {
    if (!molv.group || $('#detail-body').hidden) return;
    molv.controls.update();
    molv.renderer.render(molv.scene, molv.camera);
  });
}
function loadMolecules() {
  app.molsPromise ??= fetch('data/molecules.json').then(r => r.json()).then(m => (app.mols = m));
  return app.molsPromise;
}
async function showMolecule(d) {
  if (!molv.renderer) initMolViewer();
  $('#d-src').textContent = 'Loading conformer…';
  let m;
  try { m = (await loadMolecules())[String(d.cid)]; }
  catch { $('#d-src').textContent = 'Could not load the conformer file.'; return; }
  if (app.drugs[app.selected] !== d) return;
  disposeGroup(molv.group);
  molv.group = null;
  if (!m) { $('#d-src').textContent = 'No 3D conformer (PubChem has none and RDKit could not embed it).'; molv.renderer.clear(); return; }
  $('#d-src').textContent = m.src;
  const showH = $('#show-h').checked, n = m.z.length, keep = [];
  for (let i = 0; i < n; i++) if (showH || m.z[i] !== 1) keep.push(i);
  const P = i => new THREE.Vector3(m.xyz[3 * i], m.xyz[3 * i + 1], m.xyz[3 * i + 2]);
  const g = new THREE.Group();
  const rad = z => (z === 1 ? 0.17 : z <= 10 ? 0.3 : z <= 18 ? 0.36 : 0.42);
  const atoms = new THREE.InstancedMesh(new THREE.SphereGeometry(1, 20, 14), new THREE.MeshStandardMaterial({ roughness: 0.45 }), keep.length);
  keep.forEach((i, j) => {
    const r = rad(m.z[i]);
    atoms.setMatrixAt(j, _m.makeScale(r, r, r).setPosition(P(i)));
    atoms.setColorAt(j, _c.set(CPK[m.z[i]] || '#ff1493'));
  });
  g.add(atoms);
  // bonds: each half takes its atom's colour; double, triple and aromatic bonds as parallel sticks
  const halves = [], up = new THREE.Vector3(0, 1, 0), q = new THREE.Quaternion();
  for (let k = 0; k < m.b.length; k += 3) {
    const [i, j, o] = [m.b[k], m.b[k + 1], m.b[k + 2]];
    if (!showH && (m.z[i] === 1 || m.z[j] === 1)) continue;
    const a = P(i), b = P(j), dir = b.clone().sub(a);
    let side = new THREE.Vector3().crossVectors(dir, Math.abs(dir.y) < 0.9 * dir.length() ? up : new THREE.Vector3(1, 0, 0)).normalize();
    const sticks = o === 2 ? [[-0.07, 0.055], [0.07, 0.055]] : o === 3 ? [[-0.1, 0.045], [0, 0.045], [0.1, 0.045]] : o === 4 ? [[0, 0.07], [0.12, 0.035]] : [[0, 0.08]];
    for (const [off, r] of sticks) {
      const sh = side.clone().multiplyScalar(off), a2 = a.clone().add(sh), b2 = b.clone().add(sh), mid = a2.clone().add(b2).multiplyScalar(0.5);
      halves.push([a2, mid, r, m.z[i]], [mid, b2, r, m.z[j]]);
    }
  }
  const bonds = new THREE.InstancedMesh(new THREE.CylinderGeometry(1, 1, 1, 10, 1), new THREE.MeshStandardMaterial({ roughness: 0.5 }), Math.max(1, halves.length));
  halves.forEach(([a, b, r, z], k) => {
    const dir = b.clone().sub(a), len = dir.length();
    q.setFromUnitVectors(up, dir.normalize());
    _m.compose(a.clone().add(b).multiplyScalar(0.5), q, new THREE.Vector3(r, len, r));
    bonds.setMatrixAt(k, _m);
    bonds.setColorAt(k, _c.set(CPK[z] || '#ff1493'));
  });
  bonds.count = halves.length;
  g.add(bonds);
  molv.scene.add(g);
  molv.group = g;
  let R = 1;
  for (const i of keep) R = Math.max(R, P(i).length());
  const dist = (R + 0.6) / Math.sin(THREE.MathUtils.degToRad(molv.camera.fov / 2)) * (molv.camera.aspect < 1 ? 1 / molv.camera.aspect : 1);
  molv.camera.position.set(0, 0, dist);
  molv.controls.target.set(0, 0, 0);
  molv.controls.update();
}

// ---------------------------------------------------------------- Hansch panel
const TERM = { logP: 'log P', 'logP²': '(log P)²', MR: 'MR', TPSA: 'TPSA' };
const sig = c => (Math.abs(c) >= 10 ? c.toFixed(1) : Math.abs(c) >= 0.1 ? c.toFixed(3) : c.toPrecision(2));
function equation(g, withSe = true) {
  let s = `pChEMBL = ${neg(sig(g.coef.const))}`;
  if (withSe && g.se.const != null) s += ` (±${sig(g.se.const)})`;
  for (const t of g.terms) {
    const c = g.coef[t];
    s += ` ${c < 0 ? '−' : '+'} ${sig(Math.abs(c))}${withSe && g.se[t] != null ? ` (±${sig(g.se[t])})` : ''} ${TERM[t]}`;
  }
  return s;
}
function verdict(g) {
  if (g.q2 >= 0.5 && g.yrand_beats === 0)
    return ['Predictive within this set', `Leave-one-out q² is ${fx(g.q2)} and none of 100 shuffled activity vectors reached r² = ${fx(g.r2)}.`];
  if (g.q2 >= 0.2)
    return ['Weak', `q² is ${fx(g.q2)}: better than the mean, but below the usual 0.5 bar${g.yrand_beats ? `, and ${g.yrand_beats} of 100 shuffled fits matched its r²` : ''}.`];
  return ['Not predictive', g.q2 < 0 ? `q² is ${fx(g.q2)}: predicting every drug with the group mean does better.` : `q² is ${fx(g.q2)}, close to zero.`];
}
function statsTable(g) {
  return `<table class="stats"><tbody>
    <tr><th scope="row">Drugs</th><td>${g.n}</td><th scope="row">Model</th><td>${esc(g.terms.map(t => TERM[t]).join(' + '))}</td></tr>
    <tr><th scope="row">r²</th><td>${fx(g.r2, 3)}</td><th scope="row">q² (LOO)</th><td>${fx(g.q2, 3)}</td></tr>
    <tr><th scope="row">s</th><td>${fx(g.s, 3)}</td><th scope="row">F</th><td>${fx(g.F, 1)}</td></tr>
    <tr><th scope="row">Shuffled r²</th><td colspan="3">mean ${fx(g.yrand_r2_mean, 3)}, max ${fx(g.yrand_r2_max, 3)} (100 runs)</td></tr>
  </tbody></table>`;
}
function initHansch() {
  const sel = $('#h-group');
  sel.innerHTML = app.hansch.groups.map(g => `<option value="${esc(g.id)}">${g.id === 'ALL' ? 'All drugs with an activity value' : esc(cap(g.label))} (n = ${g.n}, q² = ${fx(g.q2)})</option>`).join('');
  sel.value = app.seriesId || app.hansch.groups[0].id;
  sel.addEventListener('change', renderHansch);
  $('#h-show').addEventListener('click', () => {
    setView({ x: 'logp', y: 'act', z: 'mr', color: 'single', group: sel.value, l1: '', l2: '', l3: '', route: '', small: true, veber: false });
    camTo(STEPS.series.cam);
    if (narrow.matches) $('#stage').scrollIntoView({ behavior: reduceMotion.matches ? 'auto' : 'smooth' });
  });
  renderHansch();
}
function renderHansch() {
  const g = app.groups.get($('#h-group').value).fit, [word, why] = verdict(g);
  const head = g.id === 'ALL' ? '' : `<p class="lede">First mechanism target <a href="https://www.ebi.ac.uk/chembl/target_report_card/${esc(g.id)}/" target="_blank" rel="noopener">${esc(g.id)}</a>, ${esc(g.label)}.</p>`;
  $('#h-out').innerHTML = `${head}<p class="fig eq-block">${esc(equation(g))}</p>${statsTable(g)}<p class="verdict"><b>${word}.</b> ${esc(why)}</p>`;
  $('#h-cands').innerHTML = `<table class="stats"><thead><tr><th>Model</th><th>r²</th><th>q²</th></tr></thead><tbody>${g.candidates.map(c => `<tr><td>${esc(c.terms.map(t => TERM[t]).join(' + '))}</td><td>${fx(c.r2, 3)}</td><td>${fx(c.q2, 3)}</td></tr>`).join('')}</tbody></table>`;
  drawObsPred(g);
}
function drawObsPred(g) {
  const svg = $('#h-plot'), W = 320, H = 260, m = { l: 40, r: 12, t: 12, b: 34 };
  const vals = g.points.flatMap(p => [p.obs, p.loo]).filter(fin);
  let lo = Math.floor(Math.min(...vals)), hi = Math.ceil(Math.max(...vals));
  if (lo === hi) hi = lo + 1;
  const X = v => m.l + ((v - lo) / (hi - lo)) * (W - m.l - m.r), Y = v => H - m.b - ((v - lo) / (hi - lo)) * (H - m.t - m.b);
  const T = app.tokens;
  let h = '';
  for (let v = lo; v <= hi; v++) {
    h += `<line x1="${X(v)}" x2="${X(v)}" y1="${m.t}" y2="${H - m.b}" stroke="${T.grid}"/><line x1="${m.l}" x2="${W - m.r}" y1="${Y(v)}" y2="${Y(v)}" stroke="${T.grid}"/>`;
    h += `<text x="${X(v)}" y="${H - m.b + 14}" text-anchor="middle">${v}</text><text x="${m.l - 6}" y="${Y(v) + 3}" text-anchor="end">${v}</text>`;
  }
  h += `<line x1="${X(lo)}" y1="${Y(lo)}" x2="${X(hi)}" y2="${Y(hi)}" stroke="${T['ink-2']}" stroke-dasharray="4 4"/>`;
  h += `<text x="${(m.l + W - m.r) / 2}" y="${H - 4}" text-anchor="middle">Predicted pChEMBL (leave-one-out)</text>`;
  h += `<text transform="translate(11 ${(m.t + H - m.b) / 2}) rotate(-90)" text-anchor="middle">Observed pChEMBL</text>`;
  for (const p of g.points) {
    if (!fin(p.loo)) continue;
    const d = app.byCid.get(p.cid);
    h += `<circle cx="${X(p.loo)}" cy="${Y(p.obs)}" r="5" fill="${T.s1}" stroke="${T.panel || 'transparent'}" stroke-width="1.5" data-cid="${p.cid}" style="cursor:pointer"><title>${esc(cap(d?.name))}: observed ${fx(p.obs)}, predicted ${fx(p.loo)}</title></circle>`;
  }
  svg.innerHTML = h;
}
$('#h-plot').addEventListener('click', e => { const c = e.target.closest('circle[data-cid]'); if (c) selectDrug(+c.dataset.cid); });

// ---------------------------------------------------------------- tour text from the data (so the numbers always match)
function fillTourText() {
  const D = app.drugs, M = app.meta;
  const pcaNames = { logp: 'log P', mr: 'MR', tpsa: 'TPSA', mw: 'weight', hbd: 'donors', hba: 'acceptors', rotb: 'rotatable bonds', arom: 'aromatic rings', fsp3: 'Fsp3' };
  const ex = M.pca.explained.map(v => (v * 100).toFixed(0) + '%');
  const top = k => Object.entries(M.pca.loadings).sort((a, b) => Math.abs(b[1][k]) - Math.abs(a[1][k])).slice(0, 2).map(e => pcaNames[e[0]] || e[0]).join(' and ');
  $('#t-pca').textContent = `PC1, PC2 and PC3 explain ${ex[0]}, ${ex[1]} and ${ex[2]} of the variance across ${D.filter(d => d.pc).length} drugs. PC1 is driven mostly by ${top(0)}, PC2 by ${top(1)}, PC3 by ${top(2)}. Colours mark the three largest ATC groups.`;
  $('#t-ts').textContent = structureMapText(D);
  const lp = D.filter(d => fin(d.d.logp) && d.d.mw <= app.hansch.max_mw), lpv = lp.map(d => d.d.logp), nBig = D.length - lp.length;
  const lo = lp.reduce((a, b) => (b.d.logp < a.d.logp ? b : a)), hi = lp.reduce((a, b) => (b.d.logp > a.d.logp ? b : a));
  const over5 = lpv.filter(v => v > 5).length;
  $('#t-logp').textContent = `Among the ${lp.length} drugs up to ${app.hansch.max_mw} Da, log P runs from ${fx(lo.d.logp, 1)} (${cap(lo.name)}) to ${fx(hi.d.logp, 1)} (${cap(hi.name)}); the median is ${fx(median(lpv), 1)}. ${over5} (${(100 * over5 / lp.length).toFixed(0)}%) are above 5, Lipinski's limit. The ${nBig} larger drugs, mostly peptides, are left out from here on: Crippen log P falls as low as ${fx(Math.min(...D.map(d => d.d.logp)), 0)} for them, and Hansch analysis is about small molecules. Depth is random spread so the points do not overlap.`;
  const both = D.filter(d => fin(d.d.logp) && fin(d.d.mr));
  const small = both.filter(d => d.d.mw <= 500), r = a => fx(pearson(a.map(d => d.d.logp), a.map(d => d.d.mr)), 2);
  $('#t-mr').textContent = `Among the ${small.length} drugs up to 500 Da, Pearson r between log P and MR is ${r(small)}: larger small molecules tend to be greasier. Across all ${both.length} it is ${r(both)}, because peptides and other large polar drugs have a high MR and a very low log P.`;
  const tp = lp.filter(d => fin(d.d.tpsa) && d.routes.length), oral = tp.filter(d => d.routes.includes('ORAL')), non = tp.filter(d => !d.routes.includes('ORAL'));
  const share = a => (100 * a.filter(d => d.d.tpsa <= 140).length / Math.max(1, a.length)).toFixed(0);
  $('#t-tpsa').textContent = `${share(oral)}% of the ${oral.length} drugs with an oral route sit at or below 140 Å², against ${share(non)}% of the ${non.length} drugs given only by other routes.`;
  const ac = D.filter(d => d.act).map(d => d.act.p);
  $('#t-act').textContent = `${ac.length} of ${D.length} drugs have a value. Median pChEMBL ${fx(median(ac), 1)} (interquartile range ${fx(quant(ac, 0.25), 1)} to ${fx(quant(ac, 0.75), 1)}), about ${fmtConc(median(ac))}.`;
  const sg = app.groups.get(app.seriesId)?.fit;
  if (sg) {
    const [word, why] = verdict(sg);
    $('#t-series-h').textContent = `${cap(sg.label)}: ${sg.n} drugs`;
    $('#t-series').textContent = `These ${sg.n} drugs share ${sg.id} (${sg.label}) as their first mechanism target. Of the ${sg.candidates.length} models allowed at this size (one descriptor per five drugs), ${sg.terms.map(t => TERM[t]).join(' + ')} gave the best leave-one-out q².`;
    $('#t-series-eq').innerHTML = `${esc(equation(sg))}<br>r² ${fx(sg.r2)} · q² ${fx(sg.q2)} · s ${fx(sg.s)} · n ${sg.n}<br><b>${word}.</b> ${esc(why)}`;
    const sticks = document.querySelector('[data-step="series"] p:last-of-type');
    if (sg.terms.includes('TPSA')) sticks.textContent = 'The surface is the fitted equation with TPSA held at the group mean. Each thin line joins a drug to its own fitted value, so its length is the residual.';
  }
  const all = app.groups.get('ALL').fit;
  $('#t-all').textContent = `Pool all ${all.n} drugs with an activity value, whatever their target, and the best model is ${all.terms.map(t => TERM[t]).join(' + ')} with r² = ${fx(all.r2)} and q² = ${fx(all.q2)}.`;
  $('#t-all-eq').innerHTML = `${esc(equation(all))}<br><b>${verdict(all)[0]}.</b> ${esc(verdict(all)[1])}`;
}

// How often a drug's nearest neighbour shares its ATC anatomical group: in the fingerprint map, in descriptor (PCA) space,
// and by chance. Fingerprints see scaffolds; bulk descriptors do not, which is why Hansch analysis needs a congeneric series.
function structureMapText(D) {
  const pts = D.filter(d => d.ts && d.pc), L1 = pts.map(d => new Set(d.atc.map(e => e.code[0])));
  const share = (i, j) => [...L1[i]].some(c => L1[j].has(c));
  const nnShare = key => {
    let hit = 0;
    for (let i = 0; i < pts.length; i++) {
      let best = -1, bd = Infinity;
      for (let j = 0; j < pts.length; j++) {
        if (j === i || pts[j].smiles === pts[i].smiles) continue;
        let s = 0;
        for (let k = 0; k < key(pts[i]).length; k++) { const t = key(pts[i])[k] - key(pts[j])[k]; s += t * t; }
        if (s < bd) { bd = s; best = j; }
      }
      if (best >= 0 && share(i, best)) hit++;
    }
    return hit / pts.length;
  };
  let chance = 0;
  for (let i = 0; i < pts.length; i++) { let c = 0; for (let j = 0; j < pts.length; j++) if (j !== i && share(i, j)) c++; chance += c / (pts.length - 1); }
  chance /= pts.length;
  const pct = v => (100 * v).toFixed(0) + '%';
  return `In this map ${pct(nnShare(d => d.ts))} of ${pts.length} drugs have a nearest neighbour from the same ATC anatomical group, against ${pct(nnShare(d => d.pc))} in the descriptor space of step 1 and ${pct(chance)} by chance.`;
}

function fillMethod() {
  const D = app.drugs, M = app.meta, ms = D.filter(d => d.conf);
  const pcs = ms.filter(d => d.conf.startsWith('PubChem')).length;
  $('#m-count').textContent = D.length;
  $('#m-checks').textContent = `${D.length} drugs, all described by RDKit; ${D.filter(d => fin(d.d.xlogp)).length} with a PubChem XLogP3; ${pcs} PubChem 3D conformers and ${ms.length - pcs} RDKit conformers; ` +
    `${D.filter(d => d.act).length} with an activity value; ${app.hansch.groups.length - 1} target groups with at least ${app.hansch.min_group} drugs. Data built ${M.generated_utc.slice(0, 10)}.`;
}

// ---------------------------------------------------------------- theme changes
function onTheme() {
  readTokens();
  recolor();
  buildAxes();
  buildExtras();
  if ($('#h-group').value) renderHansch();
  if (app.selected >= 0) drawDepiction(app.drugs[app.selected]);
}

// ---------------------------------------------------------------- boot
async function boot() {
  const status = $('#stage-status');
  try {
    const [dj, hj] = await Promise.all(['data/drugs.json', 'data/hansch.json'].map(u => fetch(u).then(r => { if (!r.ok) throw new Error(`${u}: ${r.status}`); return r.json(); })));
    app.meta = dj.meta;
    app.drugs = dj.drugs;
    app.hansch = hj;
    for (const d of app.drugs) { d._jit = (Math.imul(d.cid, 2654435761) >>> 0) / 2 ** 32 * 2 - 1; app.byCid.set(d.cid, d); }
    for (const g of hj.groups) app.groups.set(g.id, { fit: g, cids: new Set(g.points.map(p => p.cid)) });
    // the tour's one-target example: glucocorticoid receptor if it is there and fits at all, otherwise the best q² with n >= 8
    const cands = hj.groups.filter(g => g.id !== 'ALL');
    const gr = cands.find(g => g.id === 'CHEMBL2034' && g.q2 >= 0.3);
    app.seriesId = (gr || [...cands].filter(g => g.n >= 8).sort((a, b) => b.q2 - a.q2)[0] || cands[0])?.id ?? 'ALL';
    readTokens();
    const head = $('.masthead');
    new ResizeObserver(() => document.documentElement.style.setProperty('--stage-top', `${head.offsetHeight}px`)).observe(head);
    initStage();
    initControls();
    initHansch();
    fillTourText();
    fillMethod();
    syncControls();
    initTour();
    status.hidden = true;
    darkQuery.addEventListener('change', onTheme);
    new MutationObserver(onTheme).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });
    reduceMotion.addEventListener('change', () => { molv.controls && (molv.controls.autoRotate = !reduceMotion.matches); });
    narrow.addEventListener('change', () => window.ScrollTrigger?.refresh());
  } catch (err) {
    console.error(err);
    status.textContent = `Could not load the data (${err.message}).`;
  }
}
boot();
