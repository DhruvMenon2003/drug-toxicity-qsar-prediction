# Every physicochemical property on its own, and QED taken apart into the eight properties it is made of.
#
#   python3 curation/qed_properties.py --group-features h2_features/features.csv
#
# Output (folder qed_properties/):
#   property_<name>.png     one figure per property: the whole curated set's distribution, the group marked
#   qed_desirability.png    QED's eight sub-properties per group drug, each on the 0-1 desirability scale
#   qed_vs_property.png     QED against each property across the whole set, the group marked
#   grid_qed.html           mols2grid deck of the group sorted by QED, with every property in the tooltip
#   properties.csv          the numbers behind all of it
#
# RDKit's QED is a weighted average of eight desirability functions (Bickerton et al. 2012), so the
# honest way to show "the distribution of each property individually" is to show both: the raw property
# and the desirability QED assigns it. Both are plotted.

import argparse, json, os
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PROPS = [("MW", "Molecular weight (Da)"), ("ALOGP", "log P (Crippen)"), ("HBA", "H-bond acceptors"),
         ("HBD", "H-bond donors"), ("PSA", "Topological polar surface area (A^2)"),
         ("ROTB", "Rotatable bonds"), ("AROM", "Aromatic rings"), ("ALERTS", "Structural alerts")]


def row(smiles, name, source):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    p = QED.properties(m)                                   # the eight properties QED is built from
    return {"name": name, "smiles": smiles, "source": source, "QED": round(QED.qed(m), 4),
            "MW": round(p.MW, 2), "ALOGP": round(p.ALOGP, 3), "HBA": p.HBA, "HBD": p.HBD,
            "PSA": round(p.PSA, 2), "ROTB": p.ROTB, "AROM": p.AROM, "ALERTS": p.ALERTS,
            "MR": round(Crippen.MolMR(m), 2), "heavy_atoms": m.GetNumHeavyAtoms(),
            "fsp3": round(rdMolDescriptors.CalcFractionCSP3(m), 3),
            "formal_charge": Chem.GetFormalCharge(m)}


def desirability(df):
    """each property's QED desirability, 0 to 1, from RDKit's own ADS parameters."""
    out = []
    for _, r in df.iterrows():
        m = Chem.MolFromSmiles(r.smiles)
        props = QED.properties(m)
        d = {"name": r["name"]}
        for i, (key, _) in enumerate(PROPS):
            d[key] = round(QED.ads(getattr(props, key), QED.adsParameters[key]), 3)
        out.append(d)
    return pd.DataFrame(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group-features", default="h2_features/features.csv")
    ap.add_argument("--all-drugs", default=os.path.join(HERE, "..", "explorer", "data", "drugs.json"),
                    help="the curated whole-set file the Hansch explorer uses")
    ap.add_argument("--out", default="qed_properties")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    rows = []
    allj = json.load(open(args.all_drugs))
    for d in allj["drugs"]:
        r = row(d["smiles"], d["name"], "curated set")
        if r:
            rows.append(r)
    group = pd.read_csv(args.group_features).drop_duplicates("Generic Name")
    for _, g in group.iterrows():
        if isinstance(g["Structure"], str):
            r = row(g["Structure"], g["Generic Name"], "group")
            if r:
                r["atc"] = g["Anatomical Therapeutic Chemical Code"]
                rows.append(r)
    df = pd.DataFrame(rows).drop_duplicates(["name", "source"])
    df.to_csv(f"{args.out}/properties.csv", index=False)
    whole = df[df.source == "curated set"]
    grp = df[df.source == "group"]
    print(f"{len(whole)} drugs in the curated set, {len(grp)} in the group")

    # one figure per property: where the group sits in the whole set
    for key, label in [("QED", "Quantitative estimate of drug likeness")] + PROPS + \
                      [("MR", "Molar refractivity"), ("fsp3", "Fraction of sp3 carbons")]:
        fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.2, 5.4), height_ratios=[2.2, 1],
                                     gridspec_kw={"hspace": 0.32})
        vals = pd.to_numeric(whole[key], errors="coerce").dropna()
        ax.hist(vals, bins=40, color="#9ecae1", edgecolor="#5d8bb3", linewidth=0.4)
        for _, g in grp.iterrows():
            ax.axvline(g[key], color="#c1443a", lw=1.2)
            ax.annotate(g["name"], (g[key], ax.get_ylim()[1]), rotation=90, va="top", ha="right",
                        fontsize=7, color="#c1443a")
        q = vals.quantile([0.25, 0.5, 0.75])
        ax.set_title(f"{label}\ncurated set n={len(vals)}, median {q[0.5]:.2f} (IQR {q[0.25]:.2f}-{q[0.75]:.2f})",
                     fontsize=10)
        ax.set_ylabel("drugs")
        bx.scatter(pd.to_numeric(grp[key], errors="coerce"), range(len(grp)), color="#c1443a", s=26)
        bx.set_yticks(range(len(grp)))
        bx.set_yticklabels(grp["name"], fontsize=8)
        bx.set_xlim(ax.get_xlim())
        bx.grid(axis="x", color="#e3e9e7")
        bx.set_xlabel(label, fontsize=9)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False); bx.spines[s].set_visible(False)
        fig.savefig(f"{args.out}/property_{key}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # QED taken apart: the eight desirability values per group drug
    des = desirability(grp)
    fig, axes = plt.subplots(2, 4, figsize=(13, 6), sharex=True)
    for ax, (key, label) in zip(axes.ravel(), PROPS):
        ax.barh(des["name"], des[key], color="#2f7fc7")
        ax.set_xlim(0, 1)
        ax.set_title(f"{label}\ndesirability", fontsize=9)
        ax.tick_params(labelsize=8)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.suptitle("QED is the weighted average of these eight desirability functions (Bickerton et al. 2012); "
                 "1.0 means the property is in the range approved drugs live in", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{args.out}/qed_desirability.png", dpi=150)
    plt.close(fig)
    des.merge(grp[["name", "QED"]], on="name").to_csv(f"{args.out}/qed_desirability.csv", index=False)

    # QED against each property, so a low QED can be traced to the property that caused it
    fig, axes = plt.subplots(2, 4, figsize=(13.5, 6.4))
    for ax, (key, label) in zip(axes.ravel(), PROPS):
        ax.scatter(pd.to_numeric(whole[key], errors="coerce"), whole["QED"], s=5, alpha=0.3, color="#9aaaa5")
        ax.scatter(pd.to_numeric(grp[key], errors="coerce"), grp["QED"], s=30, color="#c1443a")
        for _, g in grp.iterrows():
            ax.annotate(g["name"][:10], (g[key], g["QED"]), fontsize=6.5, color="#c1443a",
                        xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel(label, fontsize=8.5); ax.set_ylabel("QED", fontsize=8.5)
        ax.tick_params(labelsize=7.5)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.suptitle("QED against each property it is made of: grey the curated set, red the group", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{args.out}/qed_vs_property.png", dpi=150)
    plt.close(fig)

    try:
        import mols2grid
        deck = grp.sort_values("QED", ascending=False)
        html = mols2grid.MolGrid(deck.rename(columns={"smiles": "SMILES"}), smiles_col="SMILES",
                                 size=(240, 190)).display(
            subset=["img", "name", "QED"],
            tooltip=[k for k, _ in PROPS] + ["MR", "fsp3", "heavy_atoms"]).data
        open(f"{args.out}/grid_qed.html", "w").write(html)
    except Exception as e:
        print(f"  mols2grid deck skipped: {type(e).__name__}: {e}")

    print(f"figures and tables -> {args.out}/")


if __name__ == "__main__":
    main()
