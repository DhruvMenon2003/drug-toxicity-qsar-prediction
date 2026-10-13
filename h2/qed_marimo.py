# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "marimo-chem-utils",
#     "mols2grid",
#     "rdkit",
#     "pandas",
#     "numpy",
#     "altair",
# ]
# ///
"""One physicochemical property at a time: where the H2 antagonists sit among the curated drugs.

QED is one number made of eight desirability functions (molecular weight, ALOGP, hydrogen-bond donors and
acceptors, polar surface area, rotatable bonds, aromatic rings and structural alerts). Looking only at the
number hides which term is pulling it down, so every term gets its own distribution here, with the five A02BA
drugs marked on it, and the molecules behind any slice shown as structures through mols2grid.

Run it with `marimo edit h2/qed_marimo.py`, or `marimo run h2/qed_marimo.py` for the app view.
Data: `h2/data/h2_features.csv` (this thread's five drugs) and `explorer/data/drugs.json` (the curated set),
read from disk when present and from GitHub otherwise.
"""

import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        # QED, one term at a time

        The five H2-receptor antagonists (ATC **A02BA**) against 759 curated drugs. Pick a property to see its
        distribution, where each H2 drug falls in it, and the structures behind any part of the curve.
        """
    )
    return


@app.cell
def _():
    import json
    import urllib.request
    from pathlib import Path

    import altair as alt
    import marimo as mo
    import marimo_chem_utils as mcu
    import mols2grid
    import numpy as np
    import pandas as pd
    from rdkit import Chem, RDLogger
    from rdkit.Chem import QED, Crippen, Descriptors

    RDLogger.DisableLog("rdApp.*")
    return (Chem, Crippen, Descriptors, QED, alt, json, mcu, mo, mols2grid, np, pd, Path, urllib)


@app.cell
def _(Path, json, mo, pd, urllib):
    _root = mo.notebook_dir().parent
    _raw = "https://raw.githubusercontent.com/DhruvMenon2003/drug-toxicity-qsar-prediction/main/"

    def _load(rel, parser):
        local = _root / rel
        if local.exists():
            return parser(local.read_text())
        with urllib.request.urlopen(_raw + rel) as fh:
            return parser(fh.read().decode())

    curated = pd.DataFrame(_load("explorer/data/drugs.json", json.loads)["drugs"])
    h2 = _load("h2/data/h2_features.csv", lambda t: pd.read_csv(__import__("io").StringIO(t)))
    return curated, h2


@app.cell
def _(Chem, Crippen, Descriptors, QED, curated, h2, pd):
    TERMS = ["MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"]

    def describe(smiles: str) -> dict | None:
        """The eight QED terms as RDKit computes them, plus QED itself and the two other weightings."""
        mol = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
        if mol is None:
            return None
        p = QED.properties(mol)
        return {
            "MW": p.MW, "ALOGP": p.ALOGP, "HBA": p.HBA, "HBD": p.HBD, "PSA": p.PSA,
            "ROTB": p.ROTB, "AROM": p.AROM, "ALERTS": p.ALERTS,
            "QED": QED.qed(mol), "QED_max": QED.weights_max(mol), "QED_none": QED.weights_none(mol),
            "MolMR": Crippen.MolMR(mol), "FractionCSP3": Descriptors.FractionCSP3(mol),
            "InChIKey": Chem.MolToInchiKey(mol),
        }

    def table(frame: pd.DataFrame, smiles_col: str, name_col: str) -> pd.DataFrame:
        rows = []
        for _, r in frame.iterrows():
            d = describe(r[smiles_col])
            if d:
                rows.append({"name": r[name_col], "SMILES": r[smiles_col], **d})
        return pd.DataFrame(rows).drop_duplicates("InChIKey").reset_index(drop=True)

    all_drugs = table(curated, "smiles", "name")
    h2_drugs = table(h2.drop_duplicates("Generic Name"), "Structure", "Generic Name")
    all_drugs["is_h2"] = all_drugs.InChIKey.isin(set(h2_drugs.InChIKey))
    return TERMS, all_drugs, h2_drugs, table


@app.cell(hide_code=True)
def _(TERMS, mo):
    prop = mo.ui.dropdown(["QED", *TERMS, "QED_max", "QED_none", "MolMR", "FractionCSP3"], value="QED",
                          label="Property")
    bins = mo.ui.slider(10, 80, value=40, label="Bins")
    mo.hstack([prop, bins], justify="start", gap=2)
    return bins, prop


@app.cell(hide_code=True)
def _(all_drugs, alt, bins, h2_drugs, mo, prop):
    _column = prop.value
    _hist = (
        alt.Chart(all_drugs)
        .mark_bar(color="#3f5f8f")
        .encode(
            x=alt.X(f"{_column}:Q", bin=alt.Bin(maxbins=bins.value), title=_column),
            y=alt.Y("count()", title="curated drugs"),
            tooltip=["count()"],
        )
        .properties(height=240, title=f"{_column} across {len(all_drugs)} curated drugs")
    )
    _marks = (
        alt.Chart(h2_drugs)
        .mark_rule(color="#6fd3ff", strokeWidth=2)
        .encode(x=f"{_column}:Q", tooltip=["name", f"{_column}:Q"])
    )
    _labels = (
        alt.Chart(h2_drugs)
        .mark_text(align="left", dx=4, dy=-6, color="#6fd3ff", fontSize=11)
        .encode(x=f"{_column}:Q", y=alt.value(10), text="name")
    )
    mo.ui.altair_chart(_hist + _marks + _labels)
    return


@app.cell(hide_code=True)
def _(all_drugs, h2_drugs, mo, prop):
    _column = prop.value
    _rank = [
        {
            "drug": r["name"],
            _column: round(r[_column], 3),
            "percentile among curated drugs": round(100 * (all_drugs[_column] < r[_column]).mean(), 1),
        }
        for _, r in h2_drugs.sort_values(_column).iterrows()
    ]
    mo.ui.table(_rank, selection=None)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### The eight terms behind one QED

        Each row is one drug, each column one desirability term, scaled so 0 is the worst and 1 the best value
        seen in the curated set. A low QED with one dark cell is a different molecule from a low QED with eight.
        """
    )
    return


@app.cell(hide_code=True)
def _(TERMS, all_drugs, alt, h2_drugs, mo, pd):
    _rows = []
    for _, r in h2_drugs.iterrows():
        for term in TERMS:
            lo, hi = all_drugs[term].min(), all_drugs[term].max()
            _rows.append({"drug": r["name"], "term": term,
                          "scaled": 0 if hi == lo else float((r[term] - lo) / (hi - lo)),
                          "value": round(float(r[term]), 2)})
    _heat = (
        alt.Chart(pd.DataFrame(_rows))
        .mark_rect()
        .encode(x=alt.X("term:N", sort=TERMS), y=alt.Y("drug:N"),
                color=alt.Color("scaled:Q", scale=alt.Scale(scheme="viridis"), title="position in range"),
                tooltip=["drug", "term", "value"])
        .properties(height=160)
    )
    mo.ui.altair_chart(_heat)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### The molecules themselves (mols2grid)""")
    return


@app.cell(hide_code=True)
def _(h2_drugs, mo, mols2grid):
    _grid = mols2grid.display(
        h2_drugs.assign(QED=h2_drugs.QED.round(3), MW=h2_drugs.MW.round(1), ALOGP=h2_drugs.ALOGP.round(2)),
        smiles_col="SMILES", subset=["img", "name", "QED"], tooltip=["MW", "ALOGP", "PSA", "HBD", "HBA", "ROTB"],
        size=(210, 160), n_cols=5,
    )
    mo.iframe(_grid.data) if hasattr(_grid, "data") else _grid
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### Where they sit in chemical space

        Morgan count fingerprints and a t-SNE projection from `marimo_chem_utils`, coloured by QED. The five
        H2 antagonists are drawn larger: structural neighbours should land next to each other, which is the
        quickest check that the fingerprints and the structures agree.
        """
    )
    return


@app.cell(hide_code=True)
def _(all_drugs, mcu, mo):
    _df = mcu.add_fingerprint_column(all_drugs.copy(), fp_type="counts_fp")
    _df = mcu.add_image_column(_df, smiles_column="SMILES")
    _df = mcu.add_tsne_columns(_df)
    _chart = mcu.interactive_chart(_df, x_col="TSNE_x", y_col="TSNE_y", color_col="QED", image_col="image")
    mo.ui.altair_chart(_chart)
    return
