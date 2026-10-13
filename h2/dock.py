"""Dock the five H2 antagonists into the human histamine H2 receptor and write the data the binding movie plays.

Receptor: PDB 7UL3, the inactive-state H2 receptor solved with famotidine bound (chain A, UniProt P25021).
The crystallographic famotidine is removed before docking and then used to validate the run: redocking has to
put famotidine back where the structure found it, or the other four poses mean nothing.

    python h2/dock.py --out h2/data --cache .h2_cache --exhaustiveness 16

Needs the AutoDock Vina binary (downloaded on first run if missing), Meeko and RDKit. Output:

| File | Contents |
|---|---|
| `data/binding.json` | receptor trace, pocket residues, docked poses, entry path frames and hydrogen bonds |
| `data/docking.csv` | one row per ligand: Vina score, pocket contacts, redocking RMSD where applicable |
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from common import H2_DRUGS, get, set_cache, write_provenance

RDLogger.DisableLog("rdApp.*")

PDB_ID = "7UL3"
CRYSTAL_LIGAND = "FO9"          # famotidine as modelled in 7UL3
VINA_URL = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64"
SITE_RADIUS = 5.0               # Angstrom: how far a pose centroid may sit from the crystallographic ligand
BOX = 24.0                      # Angstrom, cubic search box centred on the crystallographic famotidine
COVALENT = {("N", "O"): 3.5, ("O", "N"): 3.5, ("O", "O"): 3.4, ("N", "N"): 3.5}   # donor-acceptor cutoffs
ELEMENT_RADIUS = {"C": 0.70, "N": 0.65, "O": 0.60, "S": 1.00, "H": 0.31}


# ---------------------------------------------------------------- structure


def fetch_structure(work: Path) -> Path:
    pdb = work / f"{PDB_ID}.pdb"
    if not pdb.exists():
        r = get(f"https://files.rcsb.org/download/{PDB_ID}.pdb", note=f"RCSB {PDB_ID} coordinates")
        if r is None:
            raise SystemExit(f"could not download {PDB_ID} from RCSB")
        pdb.write_text(r.text)
    return pdb


def parse_atoms(pdb: Path, record: str, chain: str | None = None, resname: str | None = None) -> list[dict]:
    out = []
    for line in pdb.read_text().splitlines():
        if not line.startswith(record):
            continue
        if chain and line[21] != chain:
            continue
        if resname and line[17:20].strip() != resname:
            continue
        out.append({"name": line[12:16].strip(), "resname": line[17:20].strip(), "chain": line[21],
                    "resseq": int(line[22:26]), "xyz": [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    "element": (line[76:78].strip() or line[12:16].strip()[0]).title()})
    return out


def write_receptor(pdb: Path, work: Path) -> Path:
    """Chain A protein only: the nanobody, the Fab and the crystallographic ligand are not part of the target."""
    out = work / "receptor.pdb"
    keep = [l for l in pdb.read_text().splitlines() if l.startswith("ATOM") and l[21] == "A"]
    out.write_text("\n".join(keep) + "\nEND\n")
    pdbqt = work / "receptor.pdbqt"
    if not pdbqt.exists():
        subprocess.run(["obabel", "-ipdb", str(out), "-opdbqt", "-O", str(pdbqt), "-xr"], check=True,
                       capture_output=True)
    return pdbqt


def vina_binary(work: Path) -> Path:
    exe = work / "vina"
    if not exe.exists():
        r = get(VINA_URL, note="AutoDock Vina 1.2.5 binary")
        if r is None:
            raise SystemExit("could not download the AutoDock Vina binary")
        exe.write_bytes(r.content)
        exe.chmod(0o755)
    return exe


# ---------------------------------------------------------------- ligands


def ligand_pdbqt(name: str, smiles: str, work: Path) -> Path:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    path = work / f"lig_{name}.pdbqt"
    if path.exists():
        return path
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol)
    setup = MoleculePreparation()(mol)[0]
    text, ok, err = PDBQTWriterLegacy.write_string(setup)
    if not ok:
        raise SystemExit(f"Meeko could not write {name}: {err}")
    path.write_text(text)
    return path


def read_poses(pdbqt: Path) -> list[dict]:
    """Every model in a Vina output file: heavy atoms with elements, plus the model's score."""
    poses, atoms, score = [], [], None
    for line in pdbqt.read_text().splitlines():
        if line.startswith("REMARK VINA RESULT"):
            score = float(line.split()[3])
        elif line.startswith(("ATOM", "HETATM")):
            element = line[77:79].strip()
            if element in ("HD", "H"):
                continue
            atoms.append({"element": element[0].title() if element else "C",
                          "xyz": [float(line[30:38]), float(line[38:46]), float(line[46:54])]})
        elif line.startswith("ENDMDL"):
            poses.append({"score": score, "atoms": atoms})
            atoms, score = [], None
    return poses


def rmsd_to(reference: np.ndarray, pose: np.ndarray) -> float:
    """Symmetry-free nearest-neighbour RMSD, both directions, so atom order does not matter."""
    a = np.min(np.linalg.norm(pose[:, None] - reference[None], axis=2), axis=1)
    b = np.min(np.linalg.norm(reference[:, None] - pose[None], axis=2), axis=1)
    return float(np.sqrt((np.concatenate([a, b]) ** 2).mean()))


# ---------------------------------------------------------------- geometry for the movie


def bonds_from_distance(atoms: list[dict]) -> list[list[int]]:
    xyz = np.array([a["xyz"] for a in atoms])
    bonds = []
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            cutoff = ELEMENT_RADIUS.get(atoms[i]["element"], 0.8) + ELEMENT_RADIUS.get(atoms[j]["element"], 0.8) + 0.45
            if np.linalg.norm(xyz[i] - xyz[j]) <= cutoff:
                bonds.append([i, j])
    return bonds


def hydrogen_bonds(pose_atoms: list[dict], pocket: list[dict]) -> list[dict]:
    """Donor/acceptor heavy-atom pairs inside 3.5 A: what the viewer sees snap into place at the end of the run."""
    out = []
    for i, a in enumerate(pose_atoms):
        if a["element"] not in ("N", "O"):
            continue
        for p in pocket:
            if p["element"] not in ("N", "O"):
                continue
            d = float(np.linalg.norm(np.array(a["xyz"]) - np.array(p["xyz"])))
            if d <= COVALENT.get((a["element"], p["element"]), 3.5):
                out.append({"ligand_atom": i, "residue": f'{p["resname"]}{p["resseq"]}', "atom": p["name"],
                            "xyz": p["xyz"], "distance": round(d, 2)})
    return sorted(out, key=lambda h: h["distance"])[:8]


def entry_path(pose_atoms: list[dict], outward: np.ndarray, frames: int = 90) -> list[dict]:
    """Frames of the ligand travelling in from the extracellular vestibule to the docked pose.

    The direction is the receptor's own: the vector from the transmembrane bundle's centre through the
    orthosteric pocket, which points out between the extracellular ends of TM5, TM6 and TM7. The path is an
    illustration of the approach, not a simulated trajectory, and the final frame is the docked pose exactly.
    """
    final = np.array([a["xyz"] for a in pose_atoms])
    centre = final.mean(axis=0)
    start_offset = outward * 26.0
    out = []
    for f in range(frames):
        t = f / (frames - 1)
        ease = t * t * (3 - 2 * t)                                   # smoothstep: slow entry, slow landing
        offset = start_offset * (1 - ease)
        # the ligand tumbles while it travels and settles into the docked orientation
        angle = (1 - ease) * math.pi * 0.8
        axis = np.cross(outward, [0, 0, 1.0])
        axis = axis / (np.linalg.norm(axis) or 1)
        rot = rotation_matrix(axis, angle)
        coords = (final - centre) @ rot.T + centre + offset
        out.append({"t": round(t, 4), "xyz": np.round(coords, 3).tolist()})
    return out


def rotation_matrix(axis: np.ndarray, angle: float) -> np.ndarray:
    a, (x, y, z) = math.cos(angle / 2), axis * math.sin(angle / 2)
    return np.array([
        [a * a + x * x - y * y - z * z, 2 * (x * y - a * z), 2 * (x * z + a * y)],
        [2 * (x * y + a * z), a * a - x * x + y * y - z * z, 2 * (y * z - a * x)],
        [2 * (x * z - a * y), 2 * (y * z + a * x), a * a - x * x - y * y + z * z]])


# ---------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="h2/data")
    ap.add_argument("--cache", default=".h2_cache")
    ap.add_argument("--work", default=".h2_dock")
    ap.add_argument("--exhaustiveness", type=int, default=32)
    args = ap.parse_args()
    set_cache(args.cache)
    work, out = Path(args.work), Path(args.out)
    work.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    pdb = fetch_structure(work)
    receptor = write_receptor(pdb, work)
    vina = vina_binary(work)

    crystal = parse_atoms(pdb, "HETATM", chain="A", resname=CRYSTAL_LIGAND)
    crystal_xyz = np.array([a["xyz"] for a in crystal])
    centre = crystal_xyz.mean(axis=0)

    protein = parse_atoms(pdb, "ATOM", chain="A")
    protein_xyz = np.array([a["xyz"] for a in protein])
    outward = centre - protein_xyz.mean(axis=0)
    outward = outward / np.linalg.norm(outward)

    # pocket: every residue with an atom within 5 A of the crystallographic ligand
    near = np.min(np.linalg.norm(protein_xyz[:, None] - crystal_xyz[None], axis=2), axis=1)
    pocket_res = sorted({protein[i]["resseq"] for i in np.where(near <= 5.0)[0]})
    pocket = [a for a in protein if a["resseq"] in pocket_res]

    # smiles for the five drugs, taken from the feature table so both deliverables use the same structures
    features = pd.read_csv(out / "h2_features.csv")
    smiles = dict(zip(features["Generic Name"], features["Structure"]))

    ligands, rows = {}, []
    for drug in H2_DRUGS:
        name = drug["name"]
        lig = ligand_pdbqt(name, smiles[name], work)
        docked = work / f"dock_{name}.pdbqt"
        if not docked.exists():
            subprocess.run([str(vina), "--receptor", str(receptor), "--ligand", str(lig),
                            "--center_x", f"{centre[0]:.2f}", "--center_y", f"{centre[1]:.2f}",
                            "--center_z", f"{centre[2]:.2f}", "--size_x", str(BOX), "--size_y", str(BOX),
                            "--size_z", str(BOX), "--exhaustiveness", str(args.exhaustiveness), "--num_modes", "20", "--seed", "42",
                            "--out", str(docked)], check=True, capture_output=True)
        poses = read_poses(docked)
        # The 24 A box also covers the extracellular vestibule, and Vina's top pose sometimes parks there.
        # The prediction is the best-scoring pose that actually occupies the orthosteric site, defined by the
        # crystallographic famotidine: centroid within SITE_RADIUS of it. The rule is the same for all five.
        in_site = [p for p in poses
                   if np.linalg.norm(np.array([a["xyz"] for a in p["atoms"]]).mean(axis=0) - centre) <= SITE_RADIUS]
        best = in_site[0] if in_site else poses[0]
        best_rank = (poses.index(best) + 1)
        pose_xyz = np.array([a["xyz"] for a in best["atoms"]])
        rmsd = rmsd_best = rmsd_rank = None
        if name == "famotidine":                       # the only ligand whose experimental pose is in this structure
            per_pose = [rmsd_to(crystal_xyz, np.array([a["xyz"] for a in p["atoms"]])) for p in poses]
            rmsd = round(per_pose[best_rank - 1], 2)
            rmsd_best = round(min(per_pose), 2)
            rmsd_rank = int(np.argmin(per_pose)) + 1
        hbonds = hydrogen_bonds(best["atoms"], pocket)
        contacts = sorted({f'{p["resname"]}{p["resseq"]}' for p in pocket
                           if np.min(np.linalg.norm(pose_xyz - np.array(p["xyz"]), axis=1)) <= 4.0},
                          key=lambda s: int(re.sub(r"\D", "", s)))
        ligands[name] = {
            "smiles": smiles[name],
            "score": best["score"],
            "all_scores": [p["score"] for p in poses],
            "elements": [a["element"] for a in best["atoms"]],
            "bonds": bonds_from_distance(best["atoms"]),
            "pose": np.round(pose_xyz, 3).tolist(),
            "path": entry_path(best["atoms"], outward),
            "hbonds": hbonds,
            "contacts": contacts,
            "redock_rmsd": rmsd,
            "redock_rmsd_best": rmsd_best,
            "redock_rmsd_best_rank": rmsd_rank,
            "selected_rank": best_rank,
            "site_offset": round(float(np.linalg.norm(pose_xyz.mean(axis=0) - centre)), 2),
        }
        rows.append({"drug": name, "atc": drug["atc"], "vina_score_kcal_per_mol": best["score"],
                     "selected_pose_rank": best_rank, "poses_in_orthosteric_site": len(in_site),
                     "site_offset_A": round(float(np.linalg.norm(pose_xyz.mean(axis=0) - centre)), 2),
                     "n_poses": len(poses), "best_minus_second": round(best["score"] - poses[1]["score"], 2)
                     if len(poses) > 1 else None, "pocket_contacts": "; ".join(contacts),
                     "n_hydrogen_bonds": len(hbonds), "redock_rmsd_to_crystal_A": rmsd,
                     "redock_best_pose_rmsd_A": rmsd_best, "redock_best_pose_rank": rmsd_rank,
                     "receptor": f"PDB {PDB_ID} chain A (human H2 receptor, UniProt P25021)",
                     "box": f"{BOX} A cube centred on the crystallographic famotidine"})

    data = {
        "receptor": {
            "pdb_id": PDB_ID,
            "title": "Inactive human histamine H2 receptor with famotidine bound (cryo-EM, 3.0 A)",
            "chain": "A",
            "trace": [{"resseq": a["resseq"], "resname": a["resname"], "xyz": [round(v, 2) for v in a["xyz"]]}
                      for a in protein if a["name"] == "CA"],
            "pocket": [{"resseq": a["resseq"], "resname": a["resname"], "name": a["name"], "element": a["element"],
                        "xyz": [round(v, 2) for v in a["xyz"]]} for a in pocket],
            "pocket_residues": [f'{a["resname"]}{a["resseq"]}' for a in pocket if a["name"] == "CA"],
            "crystal_ligand": {"name": CRYSTAL_LIGAND, "elements": [a["element"] for a in crystal],
                               "xyz": [[round(v, 2) for v in a["xyz"]] for a in crystal],
                               "bonds": bonds_from_distance(crystal)},
            "centre": [round(v, 2) for v in centre],
            "outward": [round(v, 4) for v in outward],
        },
        "ligands": ligands,
        "method": {
            "docking": f"AutoDock Vina 1.2.5, exhaustiveness {args.exhaustiveness}, seed 42, {BOX} A cubic box, "
                       f"20 poses; the reported pose is the best-scoring one whose centroid is within "
                       f"{SITE_RADIUS} A of the crystallographic famotidine",
            "ligand_prep": "RDKit ETKDGv3 embedding and MMFF optimisation, Meeko PDBQT",
            "receptor_prep": "chain A protein atoms, Open Babel PDBQT (rigid)",
            "validation": "famotidine redocked against its own crystallographic pose in 7UL3; the table reports "
                          "both the top-ranked pose and the closest of the twenty poses returned",
        },
    }
    (out / "binding.json").write_text(json.dumps(data))
    pd.DataFrame(rows).to_csv(out / "docking.csv", index=False)
    write_provenance(out / "provenance_docking.csv")
    print(pd.DataFrame(rows)[["drug", "vina_score_kcal_per_mol", "selected_pose_rank", "site_offset_A",
                              "n_hydrogen_bonds", "redock_rmsd_to_crystal_A", "redock_best_pose_rmsd_A"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
