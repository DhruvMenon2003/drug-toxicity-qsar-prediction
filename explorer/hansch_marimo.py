# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "marimo-chem-utils",
#     "pandas",
#     "altair",
# ]
# ///
"""Hansch Space in marimo: the structure map and the Hansch fits from explorer/data, with structures on hover and
molecule grids for any selection. Run with `marimo edit explorer/hansch_marimo.py` (or `marimo run` for app view).
The data files are written by explorer/build_data.py; when they are not next to the notebook it reads them from GitHub."""

import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Hansch Space in marimo
    The curated DrugBank drugs, placed by structure (Morgan count fingerprints and t-SNE from `marimo_chem_utils`) and
    fitted with Hansch-type regressions of ChEMBL potency on RDKit Crippen log P, log P², molar refractivity and TPSA.
    The same data drive the 3D page in `explorer/index.html`. The fits are exploratory: drugs that share a target are not
    a congeneric series.
    """)
    return


@app.cell
def _():
    import json
    import urllib.request

    import altair as alt
    import marimo as mo
    import marimo_chem_utils as mcu
    import pandas as pd
    return alt, json, mcu, mo, pd, urllib


@app.cell
def _(json, mo, urllib):
    _local = mo.notebook_dir() / "data"
    _remote = "https://raw.githubusercontent.com/DhruvMenon2003/drug-toxicity-qsar-prediction/main/explorer/data/"

    def load(name):
        if (_local / name).exists():
            return json.loads((_local / name).read_text())
        with urllib.request.urlopen(_remote + name) as r:
            return json.load(r)

    drugs_json, hansch_json = load("drugs.json"), load("hansch.json")
    return drugs_json, hansch_json


@app.cell
def _(drugs_json, mcu, pd):
    atc1 = {k: v for k, v in drugs_json["meta"]["atc_names"].items() if len(k) == 1}
    df = pd.DataFrame([{
        "CID": d["cid"], "Name": d["name"].capitalize(), "SMILES": d["smiles"],
        "ATC group": "; ".join(sorted({atc1.get(e["code"][0], e["code"][0]).capitalize() for e in d["atc"]})),
        "Crippen logP": d["d"].get("logp"), "XLogP3": d["d"].get("xlogp"), "MR": d["d"].get("mr"), "TPSA": d["d"].get("tpsa"),
        "pChEMBL": (d["act"] or {}).get("p"), "Target": (d["act"] or {}).get("tname"),
        "TSNE_x": (d.get("ts") or [None, None])[0], "TSNE_y": (d.get("ts") or [None, None])[1],
    } for d in drugs_json["drugs"] if d.get("ts")])
    mcu.add_image_column(df)
    return (df,)


@app.cell(hide_code=True)
def _(df, mo):
    mo.md(f"""
    ### Structure map
    {df["pChEMBL"].notna().sum()} of {len(df)} drugs have a ChEMBL potency at their mechanism target; they are coloured
    by pChEMBL. Hover for the structure, drag a box to draw the selected drugs below.
    """)
    return


@app.cell
def _(df, mcu):
    tsne_chart = mcu.interactive_chart(df[df["pChEMBL"].notna()], "TSNE_x", "TSNE_y", color_col="pChEMBL",
                                       x_title="t-SNE 1", y_title="t-SNE 2")
    return (tsne_chart,)


@app.cell
def _(mcu, mo, tsne_chart):
    _sel = tsne_chart.value
    _grid = (mcu.draw_molecule_grid(_sel.assign(Label=_sel["Name"] + " " + _sel["pChEMBL"].round(1).astype(str)),
                                    legend_column="Label", num_cols=4, max_to_show=24)
             if len(_sel) else mo.md("Drag a box on the map to see those structures."))
    mo.hstack([tsne_chart, _grid], widths="equal")
    return


@app.cell(hide_code=True)
def _(hansch_json, mo):
    _opts = {f"{g['label']} ({g['id']}, n = {g['n']})": g["id"] for g in hansch_json["groups"]}
    group = mo.ui.dropdown(options=_opts, value=next(iter(_opts)), label="Target group")
    mo.vstack([mo.md("### Hansch fits\nOne regression per ChEMBL mechanism target shared by at least "
                     f"{hansch_json['min_group']} drugs, chosen by leave-one-out q²."), group])
    return (group,)


@app.cell
def _(alt, df, group, hansch_json, mo, pd):
    fit = next(g for g in hansch_json["groups"] if g["id"] == group.value)
    _terms = " ".join(f"{'+' if c >= 0 else '−'} {abs(c):.3f}·{t}" for t, c in fit["coef"].items() if t != "const")
    _verdict = ("predictive within this set" if fit["q2"] >= 0.5 and fit["yrand_beats"] == 0   # same rule as the 3D page
                else "weak" if fit["q2"] >= 0.2 else "not predictive")
    pts = pd.DataFrame(fit["points"]).merge(df[["CID", "Name", "SMILES", "image"]], left_on="cid", right_on="CID")
    _diag = alt.Chart(pd.DataFrame({"v": [pts[["obs", "loo"]].min().min(), pts[["obs", "loo"]].max().max()]})).mark_line(
        strokeDash=[4, 4], color="gray").encode(x="v:Q", y="v:Q")
    _dots = alt.Chart(pts).mark_circle(size=70).encode(
        x=alt.X("loo:Q", title="Leave-one-out predicted pChEMBL", scale=alt.Scale(zero=False)),
        y=alt.Y("obs:Q", title="Observed pChEMBL", scale=alt.Scale(zero=False)),
        tooltip=["Name", "image", alt.Tooltip("obs:Q", format=".2f"), alt.Tooltip("loo:Q", format=".2f")])
    mo.vstack([
        mo.md(f"**pChEMBL = {fit['coef']['const']:.3f} {_terms}**  \n"
              f"n {fit['n']} · r² {fit['r2']} · q² {fit['q2']} · s {fit['s']} · shuffled-activity r² mean {fit['yrand_r2_mean']}, "
              f"max {fit['yrand_r2_max']} · **{_verdict}**"),
        mo.ui.altair_chart(_diag + _dots),
        mo.ui.table(pd.DataFrame(fit["candidates"]).assign(terms=lambda t: t.terms.map(" + ".join)), selection=None,
                    label="Every model tried"),
    ])
    return (pts,)


@app.cell
def _(mcu, pts):
    mcu.draw_molecule_grid(pts.sort_values("obs", ascending=False).assign(Label=lambda t: t.Name + " " + t.obs.round(1).astype(str)),
                           legend_column="Label", num_cols=5, max_to_show=30)
    return


if __name__ == "__main__":
    app.run()
