// The binding movie: a receptor built from its UniProt topology, a ligand built from an RDKit
// conformer, and a four-beat sequence that carries the ligand from the extracellular bulk into the
// pocket and holds it there against the residues UniProt annotates as essential for binding.
//
// Everything with a number behind it is data: helix lengths are residue counts, residue markers sit
// at the residue's own fractional depth in its helix, the affinity is a ChEMBL median, the last beat
// is the ChEMBL action type. The helices' positions around the bundle are a schematic class-A
// arrangement, and the page says so in the caveat panel.

import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

const RISE = 1.5;          // A per residue along an alpha helix
const RES_PER_TURN = 3.6;
const HELIX_R = 2.3;       // A, Calpha radius of an alpha helix
const MEMBRANE = 30;       // A, hydrophobic span drawn for the bilayer
const SPREAD = 1.55;       // helices pushed apart from the bundle axis so all seven stay readable

// Class-A seven-transmembrane bundle seen from outside the cell: angle around the bundle, distance
// from its axis, and tilt. Schematic, not coordinates from a structure.
const BUNDLE = {
  1: { deg: 152, r: 11.0, tilt: 6 },
  2: { deg: 104, r: 9.2, tilt: -8 },
  3: { deg: 34, r: 6.6, tilt: 12 },
  4: { deg: 327, r: 11.6, tilt: 4 },
  5: { deg: 264, r: 9.6, tilt: -10 },
  6: { deg: 214, r: 7.6, tilt: 14 },
  7: { deg: 178, r: 5.6, tilt: -6 },
};

const CPK = { C: 0x44524e, N: 0x2a6fd6, O: 0xd14b3a, S: 0xd8b13a, F: 0x4fbf8b, Cl: 0x4aa05a,
  Br: 0xa35a2a, I: 0x8e4fb0, P: 0xe0893b, H: 0xd9e2df };
const RADIUS = { H: 0.30, C: 0.52, N: 0.50, O: 0.48, S: 0.62, F: 0.44, Cl: 0.60, Br: 0.66, I: 0.76, P: 0.64 };

const BEATS = [
  { at: 0.00, name: "approach" },
  { at: 0.34, name: "pocket mouth" },
  { at: 0.60, name: "contacts" },
  { at: 0.84, name: "verdict" },
];

const el = (id) => document.getElementById(id);
const dark = () => matchMedia("(prefers-color-scheme: dark)").matches;

let scene, camera, renderer, controls, clock;
let receptorGroup, ligandGroup, contactGroup, labelGroup, tm6Group;
let data, group, ligand, t = 0, playing = true, orbiting = true;

init();

async function init() {
  data = window.RECEPTORS || await fetch("data/receptors.json").then((r) => r.json());
  if (!data || !data.groups || !data.groups.length) {
    el("hud").textContent = "No receptor data: run curation/h2_features.py then explorer/binding/build_binding.py.";
    return;
  }
  setupScene();
  renderGroups();
  pickGroup(0);
  bindControls();
  clock = new THREE.Clock();
  animate();
}

/* ------------------------------------------------------------------ scene */
function setupScene() {
  const canvas = el("scene");
  renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  scene = new THREE.Scene();
  camera = new THREE.PerspectiveCamera(42, 1, 0.5, 600);
  camera.position.set(52, 27, 60);
  controls = new OrbitControls(camera, canvas);
  controls.enableDamping = true;
  controls.target.set(0, 2, 0);

  scene.add(new THREE.HemisphereLight(0xffffff, dark() ? 0x0b1316 : 0xc8d4d0, 0.85));
  const key = new THREE.DirectionalLight(0xffffff, 1.1);
  key.position.set(30, 44, 28);
  scene.add(key);
  const rim = new THREE.DirectionalLight(0x8fc4ff, 0.5);
  rim.position.set(-34, -12, -26);
  scene.add(rim);

  scene.add(membrane());
  receptorGroup = new THREE.Group();
  tm6Group = new THREE.Group();
  ligandGroup = new THREE.Group();
  contactGroup = new THREE.Group();
  labelGroup = new THREE.Group();
  scene.add(receptorGroup, tm6Group, ligandGroup, contactGroup, labelGroup);
  resize();
  addEventListener("resize", resize);
}

function resize() {
  const canvas = el("scene");
  const w = canvas.clientWidth, h = canvas.clientHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}

function membrane() {
  const g = new THREE.Group();
  const tails = new THREE.Mesh(
    new THREE.CylinderGeometry(27, 27, MEMBRANE, 72, 1, true),
    new THREE.MeshBasicMaterial({ color: dark() ? 0x1b2a2e : 0xdfe8e5, transparent: true,
      opacity: 0.32, side: THREE.DoubleSide, depthWrite: false }));
  g.add(tails);
  for (const y of [MEMBRANE / 2, -MEMBRANE / 2]) {        // phosphate head planes
    const ring = new THREE.Mesh(
      new THREE.RingGeometry(4, 27, 72),
      new THREE.MeshBasicMaterial({ color: dark() ? 0x2a3f45 : 0xc3d2ce, transparent: true,
        opacity: 0.5, side: THREE.DoubleSide, depthWrite: false }));
    ring.rotation.x = Math.PI / 2;
    ring.position.y = y;
    g.add(ring);
  }
  g.add(label("extracellular", new THREE.Vector3(21, MEMBRANE / 2 + 8, 0), 0.5));
  g.add(label("cytoplasm", new THREE.Vector3(21, -MEMBRANE / 2 - 8, 0), 0.5));
  return g;
}

function label(text, pos, scale = 1, color) {
  const pad = 8, font = 30;
  const c = document.createElement("canvas");
  const ctx = c.getContext("2d");
  ctx.font = `500 ${font}px "IBM Plex Mono", monospace`;
  c.width = ctx.measureText(text).width + pad * 2;
  c.height = font + pad * 2;
  const ctx2 = c.getContext("2d");
  ctx2.font = `500 ${font}px "IBM Plex Mono", monospace`;
  ctx2.fillStyle = color || (dark() ? "#dce8e5" : "#1d2c28");
  ctx2.textBaseline = "middle";
  ctx2.fillText(text, pad, c.height / 2);
  const sprite = new THREE.Sprite(new THREE.SpriteMaterial({
    map: new THREE.CanvasTexture(c), transparent: true, depthTest: false }));
  sprite.scale.set((c.width / c.height) * 3.4 * scale, 3.4 * scale, 1);
  sprite.position.copy(pos);
  return sprite;
}

/* -------------------------------------------------------------- receptor */
function helixCurve(nRes, origin, tilt) {
  const pts = [];
  const steps = Math.max(8, nRes * 2);
  for (let i = 0; i <= steps; i++) {
    const res = (i / steps) * nRes;
    const a = (res / RES_PER_TURN) * Math.PI * 2;
    const y = origin.y - res * RISE;                      // run from the extracellular side downwards
    const lean = Math.tan((tilt * Math.PI) / 180) * (res * RISE);
    pts.push(new THREE.Vector3(origin.x + Math.cos(a) * HELIX_R + lean * origin.ux,
                               y,
                               origin.z + Math.sin(a) * HELIX_R + lean * origin.uz));
  }
  return new THREE.CatmullRomCurve3(pts);
}

function helixOrigin(n, nRes) {
  const b = BUNDLE[n] || { deg: (n * 51) % 360, r: 9, tilt: 0 };
  const rad = (b.deg * Math.PI) / 180;
  const span = nRes * RISE;
  return { x: Math.cos(rad) * b.r * SPREAD, z: Math.sin(rad) * b.r * SPREAD,
           y: Math.min(MEMBRANE / 2 + 2, span / 2 + 1),
           ux: Math.cos(rad), uz: Math.sin(rad), tilt: b.tilt, r: b.r, deg: b.deg };
}

function buildReceptor(rec) {
  receptorGroup.clear(); tm6Group.clear(); labelGroup.clear();
  const helices = (rec && rec.helices) || [];
  const marks = [];
  helices.forEach((h, i) => {
    const n = Number(h.name) || i + 1;
    const nRes = h.end - h.start + 1;
    const o = helixOrigin(n, nRes);
    const curve = helixCurve(nRes, o, o.tilt);
    const mesh = new THREE.Mesh(
      new THREE.TubeGeometry(curve, Math.max(24, nRes * 3), 0.8, 10, false),
      new THREE.MeshStandardMaterial({
        color: n === 6 ? 0x2f7fc7 : (dark() ? 0x5d6f6b : 0x8fa09b),
        roughness: 0.55, metalness: 0.05,
        transparent: true, opacity: n === 4 || n === 1 || n === 2 ? 0.45 : 0.9 }));
    (n === 6 ? tm6Group : receptorGroup).add(mesh);
    const top = curve.getPoint(0).clone()
      .add(new THREE.Vector3(o.ux * 4.5, 2.4 + (n % 3) * 2.2, o.uz * 4.5));
    const tag = label(`TM${n} ${h.start}-${h.end}`, top, 0.5);
    labelGroup.add(tag);
    marks.push({ n, h, curve, o });
  });

  // pocket residues: placed at the residue's own depth inside its helix, pushed toward the bundle axis
  const sites = [];
  for (const s of (rec && rec.sites) || []) {
    const host = marks.find((m) => s.pos >= m.h.start && s.pos <= m.h.end);
    if (!host) continue;
    const f = (s.pos - host.h.start) / Math.max(1, host.h.end - host.h.start);
    const p = host.curve.getPoint(Math.min(0.98, f)).clone();
    const inward = new THREE.Vector3(-host.o.ux, 0, -host.o.uz).multiplyScalar(2.1);
    p.add(inward);
    const dot = new THREE.Mesh(new THREE.SphereGeometry(1.15, 24, 18),
      new THREE.MeshStandardMaterial({ color: 0xc94f24, roughness: 0.4 }));
    dot.position.copy(p);
    receptorGroup.add(dot);
    labelGroup.add(label(`${s.label} TM${host.n}`, p.clone().add(new THREE.Vector3(0, 2.2, 0)), 0.58, "#c94f24"));
    sites.push({ ...s, point: p, tm: host.n });
  }
  // disulfide between the extracellular ends, where UniProt reports one
  for (const [a, b] of (rec && rec.disulfides) || []) {
    const pa = marks.find((m) => Math.abs(m.h.start - a) < 12), pb = marks.find((m) => Math.abs(m.h.start - b) < 14);
    if (!pa || !pb) continue;
    const g = new THREE.BufferGeometry().setFromPoints([pa.curve.getPoint(0), pb.curve.getPoint(0)]);
    receptorGroup.add(new THREE.Line(g, new THREE.LineDashedMaterial({ color: 0xd8b13a, dashSize: 1.1, gapSize: 0.8 })
      .setValues({})) );
    receptorGroup.children.at(-1).computeLineDistances();
  }
  return { sites, pocket: pocketCentre(sites) };
}

function pocketCentre(sites) {
  if (!sites.length) return new THREE.Vector3(0, 6, 0);
  const c = new THREE.Vector3();
  sites.forEach((s) => c.add(s.point));
  return c.multiplyScalar(1 / sites.length);
}

/* ---------------------------------------------------------------- ligand */
function buildLigand(lig) {
  ligandGroup.clear();
  const atoms = lig.atoms.filter((a) => a.el !== "H");
  const index = new Map();
  lig.atoms.forEach((a, i) => index.set(i, a));
  const centre = new THREE.Vector3();
  atoms.forEach((a) => centre.add(new THREE.Vector3(...a.p)));
  centre.multiplyScalar(1 / Math.max(1, atoms.length));

  const pos = (a) => new THREE.Vector3(...a.p).sub(centre);
  lig.atoms.forEach((a, i) => {
    if (a.el === "H") return;
    const s = new THREE.Mesh(new THREE.SphereGeometry((RADIUS[a.el] || 0.55) * 1.35, 26, 20),
      new THREE.MeshStandardMaterial({ color: CPK[a.el] ?? 0x9aa5a2, roughness: 0.35, metalness: 0.08 }));
    s.position.copy(pos(a));
    s.userData.el = a.el;
    ligandGroup.add(s);
  });
  for (const [i, j] of lig.bonds) {
    const ai = index.get(i), aj = index.get(j);
    if (!ai || !aj || ai.el === "H" || aj.el === "H") continue;
    const a = pos(ai), b = pos(aj);
    const dir = new THREE.Vector3().subVectors(b, a);
    const tube = new THREE.Mesh(new THREE.CylinderGeometry(0.26, 0.26, dir.length(), 10),
      new THREE.MeshStandardMaterial({ color: dark() ? 0x9fb0ac : 0x6d7b77, roughness: 0.5 }));
    tube.position.copy(a).add(dir.clone().multiplyScalar(0.5));
    tube.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.clone().normalize());
    ligandGroup.add(tube);
  }
  return ligandGroup;
}

/* ------------------------------------------------------------- sequencing */
let staged = null;

function stage() {
  const rec = buildReceptor(group.receptor);
  buildLigand(ligand);
  staged = { ...rec };
  contactGroup.clear();
  // contact lines: nearest heavy atom of the ligand to each annotated residue, drawn at the pocket pose
  staged.contacts = staged.sites.map((s) => ({ site: s, line: null }));
}

function poseAt(u) {                                     // u in 0..1 over the whole sequence
  const pocket = staged.pocket;
  const entry = pocket.clone().setY(MEMBRANE / 2 + 13);
  const mouth = pocket.clone().setY(pocket.y + 7);
  if (u < BEATS[1].at) {                                 // approach through the extracellular bulk
    const k = u / BEATS[1].at;
    const swirl = new THREE.Vector3(Math.cos(k * 4.2) * 16 * (1 - k), 0, Math.sin(k * 4.2) * 16 * (1 - k));
    return { p: entry.clone().add(swirl).lerp(mouth, k * 0.55), spin: k * 5 };
  }
  if (u < BEATS[2].at) {                                 // slide down between TM5, TM6 and TM7
    const k = (u - BEATS[1].at) / (BEATS[2].at - BEATS[1].at);
    return { p: mouth.clone().lerp(pocket, k), spin: 2.7 + k * 1.4 };
  }
  const k = (u - BEATS[2].at) / (1 - BEATS[2].at);
  const settle = Math.sin(k * Math.PI * 2) * 0.22;       // small thermal jiggle, held in the pocket
  return { p: pocket.clone().add(new THREE.Vector3(settle, settle * 0.6, -settle)), spin: 4.1 };
}

function drawContacts(show) {
  contactGroup.clear();
  if (!show) return;
  const ligandAtoms = ligandGroup.children.filter((c) => c.userData.el);
  staged.sites.forEach((s, i) => {
    let best = null;
    for (const a of ligandAtoms) {
      const w = a.getWorldPosition(new THREE.Vector3());
      const d = w.distanceTo(s.point);
      if (!best || d < best.d) best = { d, w };
    }
    if (!best) return;
    const g = new THREE.BufferGeometry().setFromPoints([best.w, s.point]);
    const line = new THREE.Line(g, new THREE.LineDashedMaterial({
      color: 0xc94f24, dashSize: 0.8, gapSize: 0.55, transparent: true, opacity: 0.95 }));
    line.computeLineDistances();
    contactGroup.add(line);
    contactGroup.add(label(`${best.d.toFixed(1)} A to ${s.label}`,
      best.w.clone().lerp(s.point, 0.28 + 0.22 * i).add(new THREE.Vector3(0, 1.1 + i * 1.9, 0)),
      0.44, "#c94f24"));
  });
}

function beatName(u) {
  let name = BEATS[0].name;
  for (const b of BEATS) if (u >= b.at) name = b.name;
  return name;
}

function hud(u) {
  const rec = group.receptor || {};
  const action = (ligand.action || "").toLowerCase();
  const holds = /antagonist|inverse/.test(action);
  const lines = [
    `<b>${ligand.name}</b> &rarr; ${group.target || "target"} <span>(${rec.accession || "?"})</span>`,
    `beat: ${beatName(u)}`,
  ];
  if (u >= BEATS[2].at) {
    lines.push(`pChEMBL ${ligand.pchembl ?? "not established"}${ligand.pchembl ? ` &asymp; ${fmtKd(ligand.pchembl)}` : ""}`);
    lines.push(`contacts drawn to ${staged.sites.length} residue(s) UniProt calls essential for binding`);
  }
  if (u >= BEATS[3].at) {
    lines.push(holds
      ? `ChEMBL action: ${ligand.action}. TM6 (blue) stays where it is: the pocket is occupied, the receptor is not switched on.`
      : `ChEMBL action: ${ligand.action || "not stated"}. TM6 (blue) swings out: the hallmark of activation.`);
  }
  el("hud").innerHTML = lines.join("<br>");
}

function fmtKd(p) {
  const kd = Math.pow(10, -p);                            // mol/L
  if (kd < 1e-9) return `${(kd * 1e12).toFixed(0)} pM`;
  if (kd < 1e-6) return `${(kd * 1e9).toFixed(0)} nM`;
  return `${(kd * 1e6).toFixed(1)} uM`;
}

/* ----------------------------------------------------------------- render */
function animate() {
  requestAnimationFrame(animate);
  const dt = clock.getDelta();
  if (playing) {
    t = (t + dt / 16) % 1;                                // one pass about every sixteen seconds
    el("scrub").value = Math.round(t * 1000);
  }
  if (!staged) return;
  const { p, spin } = poseAt(t);
  ligandGroup.position.copy(p);
  ligandGroup.rotation.set(spin * 0.3, spin, spin * 0.17);
  const contactsOn = t >= BEATS[2].at;
  if (contactsOn !== contactGroup.userData.on) {
    contactGroup.userData.on = contactsOn;
    drawContacts(contactsOn);
  } else if (contactsOn && playing) {
    drawContacts(true);                                   // keep the distances honest while it jiggles
  }
  const holds = /antagonist|inverse/.test((ligand.action || "").toLowerCase());
  const open = t >= BEATS[3].at && !holds ? (t - BEATS[3].at) / (1 - BEATS[3].at) : 0;
  tm6Group.rotation.z = -open * 0.16;
  tm6Group.position.x = open * 1.6;
  if (orbiting) {                                         // slow orbit, closing in once the ligand is in the pocket
    controls.target.lerp(t >= BEATS[2].at ? staged.pocket : new THREE.Vector3(0, 2, 0), 0.04);
    camera.position.applyAxisAngle(new THREE.Vector3(0, 1, 0), dt * 0.12);
    const want = t >= BEATS[2].at ? 36 : 62;
    const flat = new THREE.Vector3(camera.position.x, 0, camera.position.z);
    camera.position.setLength(camera.position.length() + (want - camera.position.length()) * 0.03);
    void flat;
  }
  controls.update();
  el("beat").textContent = beatName(t);
  hud(t);
  renderer.render(scene, camera);
}

/* ------------------------------------------------------------------- rail */
function renderGroups() {
  const box = el("groups");
  box.innerHTML = "";
  data.groups.forEach((g, i) => {
    const b = document.createElement("button");
    b.className = "chip";
    b.textContent = `${g.target || g.target_chembl_id} (${g.ligands.length})`;
    b.setAttribute("aria-pressed", "false");
    b.onclick = () => pickGroup(i);
    box.appendChild(b);
  });
}

function pickGroup(i) {
  group = data.groups[i];
  [...el("groups").children].forEach((c, k) => c.setAttribute("aria-pressed", String(k === i)));
  const box = el("ligands");
  box.innerHTML = "";
  group.ligands.forEach((l, k) => {
    const b = document.createElement("button");
    b.className = "chip";
    b.textContent = `${l.name}${l.pchembl ? ` · pChEMBL ${l.pchembl}` : ""}`;
    b.onclick = () => pickLigand(k);
    box.appendChild(b);
  });
  pickLigand(0);
}

function pickLigand(k) {
  ligand = group.ligands[k];
  [...el("ligands").children].forEach((c, i) => c.setAttribute("aria-pressed", String(i === k)));
  stage();
  t = 0;
  contactGroup.userData.on = false;
  facts();
}

function facts() {
  const rec = group.receptor || {};
  const rows = [
    ["drug", `${ligand.name} · ATC ${ligand.atc} · CID ${ligand.cid}`],
    ["receptor", `${rec.name || group.target} · ${rec.length ? rec.length + " residues" : ""} · ${(rec.helices || []).length} TM helices`],
    ["pocket residues", (rec.sites || []).map((s) => s.label).join(", ") || "none annotated"],
    ["affinity", ligand.pchembl ? `pChEMBL ${ligand.pchembl} (&asymp; ${fmtKd(ligand.pchembl)}), ChEMBL median` : "not established by two documents"],
    ["action", ligand.action || "not stated"],
    ["mechanism", ligand.mechanism || "-"],
    ["properties", `MW ${ligand.mw} · cLogP ${ligand.clogp} · QED ${ligand.qed}`],
    ["structures", (rec.pdb || []).map((p) => p.id).join(", ") || "none for this accession"],
    ["sources", `<a href="${rec.uniprot_url || "#"}">UniProt</a> · <a href="${group.chembl_url}">ChEMBL target</a> · <a href="https://pubchem.ncbi.nlm.nih.gov/compound/${ligand.cid}">PubChem</a>`],
  ];
  el("factlist").innerHTML = rows.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join("");
  const sites = (rec.sites || []).map((s) => `${s.label}: ${s.note}`).join("; ");
  el("caveat").innerHTML = `Helix <em>lengths</em> are the UniProt residue ranges and each orange marker sits at
    its residue's own depth inside its helix. Where the helices sit <em>around</em> the bundle is a schematic
    class-A arrangement, not coordinates from a solved structure, so the contact distances on screen are
    indicative, not measured: treat them as a picture of which residues matter, not as docking output.
    UniProt's own words for them &mdash; ${sites || "none"}.
    ${(rec.pdb || []).length ? `For measured geometry see PDB ${(rec.pdb || []).map((p) => p.id).join(", ")}.` : ""}`;
}

function bindControls() {
  el("play").onclick = () => {
    playing = !playing;
    el("play").textContent = playing ? "Pause" : "Play";
  };
  el("scrub").oninput = (e) => {
    playing = false;
    el("play").textContent = "Play";
    t = Number(e.target.value) / 1000;
    contactGroup.userData.on = null;
  };
  el("spin").onclick = () => {
    orbiting = !orbiting;
    el("spin").textContent = orbiting ? "Orbit on" : "Orbit off";
    el("spin").setAttribute("aria-pressed", String(orbiting));
  };
}
