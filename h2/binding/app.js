/* The binding movie: a docked H2 antagonist travelling into the orthosteric pocket of PDB 7UL3.
 *
 * Everything with coordinates comes from h2/data/binding.json, written by h2/dock.py: the receptor's CA trace,
 * the pocket residues as they sit in the structure, the docked pose, the crystallographic famotidine and the
 * hydrogen bonds measured on the docked pose. The only invented part is the path the ligand takes to get there,
 * which runs along the receptor's own extracellular axis and ends exactly on the docked pose.
 */

import * as THREE from "three";

const DATA = "../data/binding.json";
const FEATURES = "../data/h2_features.json";

const ELEMENT = {
  C: { color: 0x9fb4d6, radius: 0.30 },
  N: { color: 0x5b8dff, radius: 0.30 },
  O: { color: 0xff6b6b, radius: 0.30 },
  S: { color: 0xffd166, radius: 0.36 },
  H: { color: 0xdfe7f5, radius: 0.18 },
};
const RESIDUE_COLOR = { ASP: 0xff8585, GLU: 0xff8585, LYS: 0x7aa5ff, ARG: 0x7aa5ff, TYR: 0x8fe3c3, PHE: 0x8fe3c3,
  TRP: 0x8fe3c3, CYS: 0xffd166, THR: 0xc3cfe6, VAL: 0xc3cfe6, LEU: 0xc3cfe6 };

const canvas = document.getElementById("scene");
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x0a0d14);
scene.fog = new THREE.Fog(0x0a0d14, 60, 150);

const camera = new THREE.PerspectiveCamera(42, 1, 0.1, 500);
scene.add(new THREE.AmbientLight(0xffffff, 0.55));
const key = new THREE.DirectionalLight(0xffffff, 1.5);
key.position.set(20, 30, 25);
scene.add(key);
const rim = new THREE.DirectionalLight(0x6fd3ff, 0.7);
rim.position.set(-25, -10, -20);
scene.add(rim);

const world = new THREE.Group();          // everything, re-centred on the pocket
scene.add(world);
const traceGroup = new THREE.Group();
const pocketGroup = new THREE.Group();
const ligandGroup = new THREE.Group();
const ghostGroup = new THREE.Group();
const bondGroup = new THREE.Group();      // hydrogen bonds, drawn each frame
world.add(traceGroup, pocketGroup, ligandGroup, ghostGroup, bondGroup);

const state = { data: null, features: {}, ligand: null, frame: 89, playing: true, spin: true, centre: new THREE.Vector3() };

const sphere = new THREE.SphereGeometry(1, 20, 16);
const cylinder = new THREE.CylinderGeometry(1, 1, 1, 12, 1, true);

function atomMesh(xyz, element, scale = 1) {
  const spec = ELEMENT[element] || ELEMENT.C;
  const mesh = new THREE.Mesh(sphere, new THREE.MeshStandardMaterial({
    color: spec.color, roughness: 0.42, metalness: 0.08, emissive: spec.color, emissiveIntensity: 0,
  }));
  mesh.scale.setScalar(spec.radius * scale);
  mesh.position.set(...xyz);
  return mesh;
}

function stick(a, b, color = 0x8ea4c8, radius = 0.1, opacity = 1) {
  const start = new THREE.Vector3(...a);
  const end = new THREE.Vector3(...b);
  const mesh = new THREE.Mesh(cylinder, new THREE.MeshStandardMaterial({
    color, roughness: 0.5, transparent: opacity < 1, opacity,
  }));
  place(mesh, start, end, radius);
  return mesh;
}

function place(mesh, start, end, radius) {
  const dir = new THREE.Vector3().subVectors(end, start);
  mesh.scale.set(radius, dir.length(), radius);
  mesh.position.copy(start).addScaledVector(dir, 0.5);
  mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.clone().normalize());
}

/* ---------------------------------------------------------------- receptor */

function buildReceptor(receptor) {
  const points = receptor.trace.map((r) => new THREE.Vector3(...r.xyz));
  const curve = new THREE.CatmullRomCurve3(points);
  const tube = new THREE.Mesh(
    new THREE.TubeGeometry(curve, points.length * 4, 0.30, 8, false),
    new THREE.MeshStandardMaterial({ color: 0x2f4466, roughness: 0.8, transparent: true, opacity: 0.5 }),
  );
  traceGroup.add(tube);

  const residues = new Map();
  for (const atom of receptor.pocket) {
    if (!residues.has(atom.resseq)) residues.set(atom.resseq, []);
    residues.get(atom.resseq).push(atom);
  }
  for (const [, atoms] of residues) {
    const color = RESIDUE_COLOR[atoms[0].resname] || 0xb9c6de;
    for (const atom of atoms) {
      if (atom.name === "C" || atom.name === "N" || atom.name === "O") continue;   // backbone stays in the tube
      const mesh = atomMesh(atom.xyz, atom.element, 0.62);
      mesh.material.color.setHex(color);
      pocketGroup.add(mesh);
    }
    for (let i = 0; i < atoms.length; i += 1) {
      for (let j = i + 1; j < atoms.length; j += 1) {
        const d = new THREE.Vector3(...atoms[i].xyz).distanceTo(new THREE.Vector3(...atoms[j].xyz));
        if (d < 1.95) pocketGroup.add(stick(atoms[i].xyz, atoms[j].xyz, color, 0.075));
      }
    }
  }
}

function buildGhost(crystal) {
  for (let i = 0; i < crystal.xyz.length; i += 1) {
    const mesh = atomMesh(crystal.xyz[i], crystal.elements[i], 0.55);
    mesh.material.transparent = true;
    mesh.material.opacity = 0.35;
    ghostGroup.add(mesh);
  }
  for (const [i, j] of crystal.bonds) {
    ghostGroup.add(stick(crystal.xyz[i], crystal.xyz[j], 0x7f8ba3, 0.07, 0.3));
  }
}

/* ---------------------------------------------------------------- ligand */

let ligandAtoms = [];
let ligandBonds = [];

function buildLigand(ligand) {
  ligandGroup.clear();
  ligandAtoms = ligand.elements.map((el, i) => {
    const mesh = atomMesh(ligand.pose[i], el, 1.35);
    mesh.material.emissiveIntensity = 0.28;
    ligandGroup.add(mesh);
    return mesh;
  });
  ligandBonds = ligand.bonds.map(([i, j]) => {
    const mesh = stick(ligand.pose[i], ligand.pose[j], 0xd8e3f7, 0.11);
    ligandGroup.add(mesh);
    return { mesh, i, j };
  });
}

function setFrame(frameIndex) {
  const ligand = state.data.ligands[state.ligand];
  const frame = ligand.path[frameIndex];
  frame.xyz.forEach((xyz, i) => ligandAtoms[i].position.set(...xyz));
  for (const bond of ligandBonds) {
    place(bond.mesh, new THREE.Vector3(...frame.xyz[bond.i]), new THREE.Vector3(...frame.xyz[bond.j]), 0.11);
  }
  drawHydrogenBonds(ligand, frame);
  document.getElementById("scrub").value = String(frameIndex);
  const t = frame.t;
  document.getElementById("caption").textContent = t < 0.98
    ? `Approach: ${(100 * t).toFixed(0)}% of the way in from the extracellular vestibule`
    : `Docked pose: ${ligand.score} kcal/mol, ${ligand.hbonds.length} hydrogen bonds, ${ligand.contacts.length} residues within 4 Å`;
}

function drawHydrogenBonds(ligand, frame) {
  bondGroup.clear();
  if (frame.t < 0.8) return;
  const strength = (frame.t - 0.8) / 0.2;
  for (const hbond of ligand.hbonds) {
    const from = new THREE.Vector3(...frame.xyz[hbond.ligand_atom]);
    const to = new THREE.Vector3(...hbond.xyz);
    const steps = 7;
    for (let s = 0; s < steps; s += 1) {
      if (s % 2) continue;
      const a = from.clone().lerp(to, s / steps);
      const b = from.clone().lerp(to, (s + 1) / steps);
      const dash = stick([a.x, a.y, a.z], [b.x, b.y, b.z], 0x6fd3ff, 0.045, 0.25 + 0.65 * strength);
      bondGroup.add(dash);
    }
  }
}

/* ---------------------------------------------------------------- panel */

function showStats() {
  const ligand = state.data.ligands[state.ligand];
  const rows = state.features[state.ligand] || [];
  const stats = {
    "Vina score": `${ligand.score} kcal/mol`,
    "Pose rank": `${ligand.selected_rank} of 20`,
    "Hydrogen bonds": ligand.hbonds.length,
    "Residues within 4 Å": ligand.contacts.length,
  };
  if (rows[0]) {
    Object.assign(stats, {
      QED: rows[0].qed,
      "pChEMBL at H2": rows[0].pchembl ?? "—",
      "Rat oral LD50": rows[0].ld50 ? `${rows[0].ld50} mg/kg` : "—",
      Routes: rows.map((r) => r.route).join(", "),
    });
  }
  document.getElementById("stats").innerHTML = Object.entries(stats)
    .map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join("");

  const bonded = new Set(ligand.hbonds.map((h) => h.residue));
  document.getElementById("contacts").innerHTML = "<b>Pocket contacts (4 Å)</b>" + ligand.contacts
    .map((c) => `<span class="chip${bonded.has(c) ? " bonded" : ""}">${c}</span>`).join("");
  document.getElementById("hbonds").innerHTML = ligand.hbonds.length
    ? "<b>Hydrogen bonds</b>" + ligand.hbonds
      .map((h) => `<span class="chip bonded">${h.residue} ${h.atom} · ${h.distance} Å</span>`).join("")
    : "<b>Hydrogen bonds</b>none within 3.5 Å in this pose";

  const famotidine = state.data.ligands.famotidine;
  document.getElementById("validation").innerHTML =
    `Redocking check: famotidine's selected pose sits ${famotidine.redock_rmsd} Å from its experimental pose in
     7UL3 (closest of the twenty poses: ${famotidine.redock_rmsd_best} Å). ${state.data.method.docking}.`;
}

/* ---------------------------------------------------------------- loop */

function resize() {
  const { clientWidth: w, clientHeight: h } = canvas.parentElement;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}

let angle = 0.6;
function animate() {
  requestAnimationFrame(animate);
  if (state.playing) {
    state.frame = Math.min(state.frame + 1, 89);
    setFrame(state.frame);
    if (state.frame === 89) state.playing = false;
  }
  if (state.spin) angle += 0.0022;
  const radius = 48;
  camera.position.set(Math.cos(angle) * radius, 13, Math.sin(angle) * radius);
  camera.lookAt(0, 0, 0);
  renderer.render(scene, camera);
}

function selectLigand(name) {
  state.ligand = name;
  buildLigand(state.data.ligands[name]);
  state.frame = 0;
  state.playing = true;
  showStats();
}

async function main() {
  const [data, features] = await Promise.all([
    fetch(DATA).then((r) => r.json()),
    fetch(FEATURES).then((r) => r.json()).catch(() => ({})),
  ]);
  state.data = data;
  state.features = features;
  buildReceptor(data.receptor);
  buildGhost(data.receptor.crystal_ligand);
  state.centre.set(...data.receptor.centre);
  world.position.copy(state.centre).multiplyScalar(-1);

  const select = document.getElementById("ligand");
  select.innerHTML = Object.keys(data.ligands)
    .map((n) => `<option value="${n}">${n[0].toUpperCase() + n.slice(1)}</option>`).join("");
  select.onchange = () => selectLigand(select.value);
  document.getElementById("play").onclick = () => { state.frame = 0; state.playing = true; };
  document.getElementById("ghost").onchange = (e) => { ghostGroup.visible = e.target.checked; };
  document.getElementById("spin").onchange = (e) => { state.spin = e.target.checked; };
  document.getElementById("pocketOnly").onchange = (e) => { traceGroup.visible = !e.target.checked; };
  document.getElementById("scrub").oninput = (e) => {
    state.playing = false;
    state.frame = Number(e.target.value);
    setFrame(state.frame);
  };
  document.getElementById("loading").remove();

  selectLigand(Object.keys(data.ligands)[0]);
  resize();
  addEventListener("resize", resize);
  animate();
}

main();
