"""Static version of the QED notebook: the figures and the structure grid, without running marimo.

This is the part of deliverable 2 you can check at a glance. It writes, for the five A02BA antagonists against
the 759 curated drugs:

| File | What it shows |
|---|---|
| `figures/qed_distributions.png` | one histogram per QED term with each H2 drug marked on it |
| `figures/qed_terms_heatmap.png` | the eight desirability terms per drug, scaled over the curated range |
| `figures/h2_structures.png` | the five structures drawn from the curated SMILES, with QED under each |
| `figures/h2_grid.html` | the same five as an interactive mols2grid page |
| `data/qed_table.csv` | every term, QED, and each drug's percentile among the curated drugs |

    python h2/qed_report.py --out h2
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                    # noqa: E402
import numpy as np                                                 # noqa: E402
import pandas as pd                                                # noqa: E402
from rdkit import Chem, RDLogger                                   # noqa: E402
from rdkit.Chem import Crippen, Descriptors, Draw, QED             # noqa: E402

RDLogger.DisableLog("rdApp.*")
TERMS = ["MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"]
TERM_LABEL = {"MW": "molecular weight", "ALOGP": "ALOGP", "HBA": "H-bond acceptors", "HBD": "H-bond donors",
              "PSA": "polar surface area", "ROTB": "rotatable bonds", "AROM": "aromatic rings",
              "ALERTS": "structural alerts"}


def describe(smiles: str) -> dict | None:
    mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if mol is None:
        return None
    p = QED.properties(mol)
    return {"MW": p.MW, "ALOGP": p.ALOGP, "HBA": p.HBA, "HBD": p.HBD, "PSA": p.PSA, "ROTB": p.ROTB,
            "AROM": p.AROM, "ALERTS": p.ALERTS, "QED": QED.qed(mol), "MolMR": Crippen.MolMR(mol),
            "FractionCSP3": Descriptors.FractionCSP3(mol), "InChIKey": Chem.MolToInchiKey(mol)}


def frame(records: list[tuple[str, str]]) -> pd.DataFrame:
    rows = [{"name": name, "SMILES": smiles, **d} for name, smiles in records if (d := describe(smiles))]
    return pd.DataFrame(rows).drop_duplicates("InChIKey").reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="h2")
    ap.add_argument("--curated", default="explorer/data/drugs.json")
    args = ap.parse_args()
    root = Path(args.out)
    figures, data = root / "figures", root / "data"
    figures.mkdir(parents=True, exist_ok=True)

    h2 = pd.read_csv(data / "h2_features.csv").drop_duplicates("Generic Name")
    h2_df = frame(list(zip(h2["Generic Name"], h2["Structure"])))
    curated = json.loads(Path(args.curated).read_text())["drugs"]
    all_df = frame([(d["name"], d["smiles"]) for d in curated])
    all_df["is_h2"] = all_df.InChIKey.isin(set(h2_df.InChIKey))

    # ---- one histogram per term, with the five drugs marked
    columns = ["QED", *TERMS]
    fig, axes = plt.subplots(3, 3, figsize=(15, 10.5), facecolor="white")
    colours = plt.cm.turbo(np.linspace(0.12, 0.92, len(h2_df)))
    for ax, column in zip(axes.ravel(), columns):
        # peptides stretch MW, PSA and HBA by an order of magnitude; the axis is clipped to the 1st-99th
        # percentile of the curated set so the small-molecule region stays readable, and the title says so
        lo, hi = np.percentile(all_df[column], [1, 99])
        span = max(hi - lo, 1e-6)
        lo, hi = lo - 0.05 * span, hi + 0.05 * span
        lo = min(lo, float(h2_df[column].min()) - 0.02 * span)      # never clip a marked drug off the axis
        hi = max(hi, float(h2_df[column].max()) + 0.02 * span)
        clipped = int(((all_df[column] < lo) | (all_df[column] > hi)).sum())
        ax.hist(all_df[column].clip(lo, hi), bins=40, range=(lo, hi), color="#9fb6d8", edgecolor="white",
                linewidth=0.4)
        ax.set_xlim(lo, hi)
        for colour, (_, row) in zip(colours, h2_df.sort_values(column).iterrows()):
            ax.axvline(row[column], color=colour, linewidth=1.8, label=row["name"])
        title = f'{column}  -  {TERM_LABEL.get(column, "quantitative estimate of drug likeness")}'
        ax.set_title(title + (f"  ({clipped} drugs outside this range)" if clipped else ""), fontsize=10.5)
        ax.set_ylabel("curated drugs", fontsize=9)
        ax.tick_params(labelsize=8.5)
    axes.ravel()[0].legend(fontsize=8, frameon=False, loc="upper right")
    fig.suptitle(f"QED and its eight terms: five A02BA H2 antagonists against {len(all_df)} curated drugs",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(figures / "qed_distributions.png", dpi=150)
    plt.close(fig)

    # ---- the eight terms per drug, scaled over the curated range
    scaled = pd.DataFrame(
        [[0 if all_df[t].max() == all_df[t].min() else
          (row[t] - all_df[t].min()) / (all_df[t].max() - all_df[t].min()) for t in TERMS]
         for _, row in h2_df.iterrows()], index=h2_df["name"], columns=TERMS)
    fig, ax = plt.subplots(figsize=(9, 3.4), facecolor="white")
    im = ax.imshow(scaled.values, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(TERMS)), TERMS, fontsize=9)
    ax.set_yticks(range(len(h2_df)), h2_df["name"], fontsize=9)
    for i in range(len(h2_df)):
        for j, term in enumerate(TERMS):
            ax.text(j, i, f"{h2_df.iloc[i][term]:.0f}" if term in ("HBA", "HBD", "ROTB", "AROM", "ALERTS")
                    else f"{h2_df.iloc[i][term]:.1f}", ha="center", va="center", fontsize=8,
                    color="white" if scaled.values[i, j] < 0.6 else "black")
    fig.colorbar(im, ax=ax, label="position in the curated range")
    ax.set_title("QED terms, scaled over the curated set (cell text is the raw value)", fontsize=11)
    fig.tight_layout()
    fig.savefig(figures / "qed_terms_heatmap.png", dpi=150)
    plt.close(fig)

    # ---- the structures, so a wrong SMILES is visible rather than hidden in a number
    mols = [Chem.MolFromSmiles(s) for s in h2_df.SMILES]
    legends = [f'{r["name"]}  QED {r.QED:.2f}  MW {r.MW:.0f}' for _, r in h2_df.iterrows()]
    Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(300, 250), legends=legends).save(
        figures / "h2_structures.png")

    try:
        import mols2grid
        grid = mols2grid.MolGrid(h2_df.assign(QED=h2_df.QED.round(3), MW=h2_df.MW.round(1)), smiles_col="SMILES")
        (figures / "h2_grid.html").write_text(
            grid.render(subset=["img", "name", "QED"], tooltip=["MW", "ALOGP", "PSA", "HBD", "HBA", "ROTB"],
                        n_items_per_page=10))
    except Exception as exc:                                        # mols2grid is optional for the static run
        print(f"mols2grid page skipped: {exc}")

    # ---- the numbers behind the figures
    out = h2_df.copy()
    for column in ["QED", *TERMS, "MolMR", "FractionCSP3"]:
        out[f"{column}_percentile"] = [round(100 * (all_df[column] < v).mean(), 1) for v in out[column]]
    out.to_csv(data / "qed_table.csv", index=False)
    print(out[["name", "QED", "QED_percentile", "MW", "ALOGP", "PSA", "ALERTS"]].to_string(index=False))


if __name__ == "__main__":
    main()
