# /// script
# requires-python = ">=3.11"
# dependencies = ["marimo", "pandas", "rdkit", "mols2grid", "altair", "marimo-chem-utils"]
# ///
# Interactive version of curation/qed_properties.py: pick a property, see its distribution across the
# whole curated set, see where the chosen ATC group sits in it, and look at the molecules behind any
# slice of that distribution.
#
#   marimo edit explorer/qed_properties_marimo.py
#
# marimo-chem-utils supplies the structure images, count fingerprints and the t-SNE map;
# mols2grid supplies the card deck; RDKit's QED module supplies both the score and the eight
# desirability values it is the average of.

import marimo

__generated_with = "0.9.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import json, os
    import marimo as mo
    import pandas as pd
    from rdkit import Chem, RDLogger
    from rdkit.Chem import QED, Crippen, rdMolDescriptors

    RDLogger.DisableLog("rdApp.*")
    return Chem, Crippen, QED, json, mo, os, pd, rdMolDescriptors


@app.cell
def _(mo):
    mo.md(
        """
        # Physicochemical properties, one at a time

        QED is a single number standing in for eight properties at once. Picking a property below shows
        what that one property does on its own across every curated drug, and where the drugs of one ATC
        chemical group fall inside it.
        """
    )
    return


@app.cell
def _(Chem, Crippen, QED, json, os, pd, rdMolDescriptors):
    HERE = os.path.dirname(os.path.abspath(__file__))
    ALL = os.path.join(HERE, "data", "drugs.json")
    GROUP = os.environ.get("GROUP_FEATURES", "h2_features/features.csv")

    def describe(smiles, name, source, atc=""):
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        p = QED.properties(m)
        d = {"name": name, "SMILES": smiles, "source": source, "atc": atc,
             "QED": round(QED.qed(m), 4), "MW": round(p.MW, 2), "ALOGP": round(p.ALOGP, 3),
             "HBA": p.HBA, "HBD": p.HBD, "PSA": round(p.PSA, 2), "ROTB": p.ROTB,
             "AROM": p.AROM, "ALERTS": p.ALERTS, "MR": round(Crippen.MolMR(m), 2),
             "fsp3": round(rdMolDescriptors.CalcFractionCSP3(m), 3)}
        for key in ["MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"]:
            d[f"des_{key}"] = round(QED.ads(getattr(p, key), QED.adsParameters[key]), 3)
        return d

    rows = [describe(d["smiles"], d["name"], "curated set") for d in json.load(open(ALL))["drugs"]]
    if os.path.exists(GROUP):
        g = pd.read_csv(GROUP).drop_duplicates("Generic Name")
        rows += [describe(r["Structure"], r["Generic Name"], "group",
                          r["Anatomical Therapeutic Chemical Code"])
                 for _, r in g.iterrows() if isinstance(r["Structure"], str)]
    df = pd.DataFrame([r for r in rows if r]).drop_duplicates(["name", "source"])
    return GROUP, df, describe


@app.cell
def _(df, mo):
    prop = mo.ui.dropdown(
        options=["QED", "MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS", "MR", "fsp3"],
        value="QED", label="Property")
    rng = mo.ui.range_slider(start=0.0, stop=1.0, step=0.01, value=[0.0, 1.0],
                             label="Percentile slice of the curated set")
    mo.hstack([prop, rng, mo.md(f"**{len(df)}** molecules loaded")])
    return prop, rng


@app.cell
def _(df, mo, prop, rng):
    import altair as alt

    whole = df[df.source == "curated set"]
    grp = df[df.source == "group"]
    lo, hi = whole[prop.value].quantile(rng.value[0]), whole[prop.value].quantile(rng.value[1])
    slice_ = whole[(whole[prop.value] >= lo) & (whole[prop.value] <= hi)]

    hist = alt.Chart(whole).mark_bar(color="#9ecae1").encode(
        alt.X(f"{prop.value}:Q", bin=alt.Bin(maxbins=40), title=prop.value),
        alt.Y("count()", title="drugs"))
    marks = alt.Chart(grp).mark_rule(color="#c1443a", size=2).encode(
        x=f"{prop.value}:Q", tooltip=["name", prop.value, "QED", "atc"])
    chart = mo.ui.altair_chart((hist + marks).properties(height=260, width="container"))
    mo.vstack([chart,
               mo.md(f"Slice holds **{len(slice_)}** drugs, {prop.value} from "
                     f"**{lo:.2f}** to **{hi:.2f}**; the red rules are the group.")])
    return chart, grp, slice_, whole


@app.cell
def _(grp, mo, slice_):
    import mols2grid

    deck = slice_.sort_values("QED", ascending=False).head(60)
    mo.vstack([
        mo.md("### The group, sorted by QED"),
        mo.Html(mols2grid.MolGrid(grp.sort_values("QED", ascending=False), smiles_col="SMILES",
                                  size=(230, 180)).display(
            subset=["img", "name", "QED"],
            tooltip=["atc", "MW", "ALOGP", "HBA", "HBD", "PSA", "ROTB", "AROM", "ALERTS"]).data),
        mo.md("### The slice of the curated set behind the bars above (first 60)"),
        mo.Html(mols2grid.MolGrid(deck, smiles_col="SMILES", size=(210, 165)).display(
            subset=["img", "name", "QED"], tooltip=["MW", "ALOGP", "PSA"]).data),
    ])
    return deck, mols2grid


@app.cell
def _(grp, mo):
    import altair as alt2

    des = grp.melt(id_vars=["name"], value_vars=[c for c in grp.columns if c.startswith("des_")],
                   var_name="property", value_name="desirability")
    des["property"] = des["property"].str.replace("des_", "", regex=False)
    mo.vstack([
        mo.md("### QED taken apart\nEach bar is one of the eight desirability functions QED averages; "
              "a short bar is the property that is holding the score down."),
        mo.ui.altair_chart(alt2.Chart(des).mark_bar(color="#2f7fc7").encode(
            x=alt2.X("desirability:Q", scale=alt2.Scale(domain=[0, 1])),
            y=alt2.Y("name:N", title=""),
            row=alt2.Row("property:N", title="")).properties(height=70, width=420)),
    ])
    return


@app.cell
def _(grp, mo, whole):
    # structure map: marimo-chem-utils count fingerprints + t-SNE, the group marked on the whole set
    try:
        from marimo_chem_utils import add_fingerprint_column, add_image_column, add_tsne_columns, interactive_chart
        import pandas as pd2

        both = pd2.concat([whole.sample(min(400, len(whole)), random_state=0), grp])
        both = add_fingerprint_column(both, fp_type="counts_fp", smiles_column="SMILES")
        both = add_image_column(both, smiles_column="SMILES")
        both = add_tsne_columns(both)
        out = mo.ui.altair_chart(interactive_chart(both, x_col="TSNE_x", y_col="TSNE_y",
                                                   color_col="QED", image_col="image"))
    except Exception as e:
        out = mo.md(f"marimo-chem-utils map skipped: `{type(e).__name__}: {e}`")
    mo.vstack([mo.md("### Where the group sits in structure space"), out])
    return


if __name__ == "__main__":
    app.run()
