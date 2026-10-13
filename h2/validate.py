"""Cross-check the H2 run and write `h2/data/validation.md`.

The point is that nothing here trusts the pipeline that produced the numbers: each check recomputes a value
from a second source, or from the structure itself, and the report says plainly when a check fails.

    python h2/validate.py --out h2
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED

RDLogger.DisableLog("rdApp.*")


def check(name: str, ok: bool, detail: str) -> dict:
    return {"check": name, "result": "pass" if ok else "LOOK", "detail": detail}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="h2")
    args = ap.parse_args()
    root = Path(args.out)
    data = root / "data"
    features = pd.read_csv(data / "h2_features.csv")
    docking = pd.read_csv(data / "docking.csv")
    binding = json.loads((data / "binding.json").read_text())
    rows: list[dict] = []

    # 1. every SMILES parses, and the stored InChIKey is the one RDKit derives from it
    bad = [r["Generic Name"] for _, r in features.iterrows()
           if Chem.MolFromSmiles(r["Structure"]) is None
           or Chem.MolToInchiKey(Chem.MolFromSmiles(r["Structure"])) != r["InChIKey"]]
    rows.append(check("SMILES parse and InChIKey agree", not bad, f"{len(features)} rows, mismatches: {bad or 'none'}"))

    # 2. no two different drugs share a structure (the reason the curated SMILES column exists)
    pairs = features.drop_duplicates("Generic Name").groupby("InChIKey")["Generic Name"].nunique()
    rows.append(check("one structure per drug", (pairs <= 1).all(),
                      f'{len(pairs)} distinct structures for {features["Generic Name"].nunique()} drugs'))

    # 3. RDKit weight against the PubChem weight for the same CID
    deltas = []
    for _, r in features.drop_duplicates("Generic Name").iterrows():
        mol = Chem.MolFromSmiles(r["Structure"])
        deltas.append(abs(Descriptors.MolWt(mol) - r["Molecular Weight"]))
    rows.append(check("RDKit molecular weight reproducible", max(deltas) < 0.05,
                      f"largest difference {max(deltas):.3f} Da"))

    # 4. QED recomputed from the SMILES matches the stored value, and sits near ChEMBL's own QED
    qed_delta, chembl_delta = [], []
    for _, r in features.drop_duplicates("Generic Name").iterrows():
        mol = Chem.MolFromSmiles(r["Structure"])
        qed_delta.append(abs(QED.qed(mol) - r["Quantitative Estimate of Drug Likeness"]))
        if pd.notna(r.get("ChEMBL QED")):
            chembl_delta.append((r["Generic Name"], round(abs(QED.qed(mol) - float(r["ChEMBL QED"])), 3)))
    rows.append(check("QED reproducible from the SMILES", max(qed_delta) < 0.002,
                      f"largest difference {max(qed_delta):.4f}"))
    worst = max(chembl_delta, key=lambda x: x[1]) if chembl_delta else ("none", 0)
    rows.append(check("QED close to ChEMBL's own value", worst[1] <= 0.1,
                      f"largest difference {worst[1]} ({worst[0]}); ChEMBL rounds to two decimals and uses its "
                      f"own parent structure"))

    # 5. Crippen log P against PubChem XLogP3 (two different estimators, so a loose bound)
    gaps = [(r["Generic Name"], round(abs(Crippen.MolLogP(Chem.MolFromSmiles(r["Structure"]))
                                          - r["Lipophilicity(XLogP3, PubChem)"]), 2))
            for _, r in features.drop_duplicates("Generic Name").iterrows()
            if pd.notna(r["Lipophilicity(XLogP3, PubChem)"])]
    worst_logp = max(gaps, key=lambda x: x[1]) if gaps else ("none", 0)
    rows.append(check("Crippen log P and PubChem XLogP3 agree within 1.5 log units", worst_logp[1] <= 1.5,
                      f"largest gap {worst_logp[1]} ({worst_logp[0]}) over {len(gaps)} drugs"))

    # 6. redocking: the famotidine pose has to come back where the structure put it
    fam = binding["ligands"]["famotidine"]
    rows.append(check("redocking recovers the crystallographic famotidine pose (< 2 A)", fam["redock_rmsd"] < 2.0,
                      f'selected pose {fam["redock_rmsd"]} A (rank {fam["selected_rank"]} of 20); closest pose '
                      f'{fam["redock_rmsd_best"]} A at rank {fam["redock_rmsd_best_rank"]}'))

    # 7. every docked pose sits in the orthosteric site, and touches the conserved Asp98 of aminergic receptors
    offsets = {n: l["site_offset"] for n, l in binding["ligands"].items()}
    rows.append(check("all poses inside the orthosteric site (centroid < 5 A from the crystal ligand)",
                      max(offsets.values()) <= 5.0, json.dumps(offsets)))
    asp = {n: ("ASP98" in l["contacts"]) for n, l in binding["ligands"].items()}
    rows.append(check("poses contact Asp98, the conserved aminergic anchor", all(asp.values()),
                      f'contacting: {[n for n, v in asp.items() if v]}; not: {[n for n, v in asp.items() if not v]}'))

    # 8. routes: each drug keeps at most two non-intravenous routes, and each row is one ATC code x route pair
    per_drug = features[~features["Mode of Administration"].str.startswith("INTRAVENOUS", na=False)]
    counts = per_drug.groupby(["Generic Name", "Anatomical Therapeutic Chemical Code"])["Mode of Administration"].nunique()
    rows.append(check("at most two non-intravenous routes per ATC code (rules 3 and 4)", (counts <= 2).all(),
                      counts.to_dict().__str__()))
    rows.append(check("one row per ATC code and route (rule 5)",
                      not features.duplicated(["Generic Name", "Anatomical Therapeutic Chemical Code",
                                               "Mode of Administration"]).any(),
                      f"{len(features)} rows for {features['Generic Name'].nunique()} drugs"))

    # 9. bioavailability: intravenous rows are 1.0 by definition, oral rows below it
    iv = features[features["Mode of Administration"].str.startswith("INTRAVENOUS", na=False)]
    oral = features[features["Mode of Administration"] == "ORAL"]
    rows.append(check("intravenous bioavailability is 1 and oral is below it",
                      (iv["Bioavailability(Fraction Absorbed)"] == 1.0).all()
                      and (oral["Bioavailability(Fraction Absorbed)"].dropna() < 1.0).all(),
                      f'IV rows {list(iv["Bioavailability(Fraction Absorbed)"])}, '
                      f'oral rows {list(oral["Bioavailability(Fraction Absorbed)"].dropna())}'))

    # 10. every filled cell carries a source, every empty one carries a reason
    unsourced = []
    for value_col, source_col in [("In Vitro Activity (pChEMBL at H2)", "In Vitro Activity_source"),
                                  ("In Vivo Toxicity (rat oral LD50, mg/kg)", "In Vivo Toxicity_source"),
                                  ("Bioavailability(Fraction Absorbed)", "Bioavailability_source"),
                                  ("Half Life (h)", "Half Life_source"), ("Dosage (mg)", "Dosage_source"),
                                  ("Mode of Administration", "Mode of Administration_source")]:
        for _, r in features.iterrows():
            if pd.isna(r[source_col]) or not str(r[source_col]).strip():
                unsourced.append((r["Generic Name"], value_col))
            elif pd.isna(r[value_col]) and "left empty" not in str(r[source_col]) and "no route" not in str(r[source_col]):
                unsourced.append((r["Generic Name"], value_col + " (empty without a reason)"))
    rows.append(check("every value has a source and every gap has a reason", not unsourced, str(unsourced or "none")))

    # 11. therapeutic index: the stored value matches the formula the spreadsheet carries
    import math
    bad_ti = []
    for _, r in features.iterrows():
        tox, act = r["In Vivo Toxicity (rat oral LD50, mg/kg)"], r["In Vivo Activity (1/effective dose, kg/mg)"]
        if pd.notna(tox) and pd.notna(act):
            expected = round(math.log10(tox) - math.log10(1 / act), 3)
            if abs(expected - r["Therapeutic Index_value"]) > 1e-6:
                bad_ti.append(r["Generic Name"])
    rows.append(check("therapeutic index equals log(toxicity) - log(1/activity)", not bad_ti,
                      f"checked {features['Therapeutic Index_value'].notna().sum()} rows, off: {bad_ti or 'none'}"))

    report = pd.DataFrame(rows)
    lines = ["# Validation of the H2 antagonist run", "",
             f"Rows: {len(features)} (one per drug x ATC code x route) for {features['Generic Name'].nunique()} drugs.",
             f"Receptor: {binding['receptor']['title']} ({binding['receptor']['pdb_id']}).", "",
             report.to_markdown(index=False), "",
             "## Docking", "", docking.to_markdown(index=False), "",
             "## What is measured and what is illustrated", "",
             "- Pose coordinates, Vina scores, pocket contacts and hydrogen-bond distances are computed from the",
             "  structure and the docking run; the binding page shows them unchanged.",
             "- The path the ligand takes into the pocket is an illustration along the receptor's extracellular",
             "  axis. It is not a molecular dynamics trajectory and no energy is computed along it.",
             "- FAERS counts are spontaneous reports: a reporting signal, never an incidence rate.",
             "- The rat oral LD50 is an acute single-dose animal value; the therapeutic index built from it and",
             "  the human label dose is a ratio for comparing these five drugs with each other, not a clinical",
             "  safety margin.",
             ]
    (data / "validation.md").write_text("\n".join(lines) + "\n")
    print(report.to_string(index=False))
    if (report.result == "LOOK").any():
        print("\nchecks needing a look:", list(report[report.result == "LOOK"]["check"]))


if __name__ == "__main__":
    main()
