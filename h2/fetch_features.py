"""Build the feature table for the five A02BA H2-receptor antagonists, one row per drug x ATC code x route.

Every column is fetched or computed; nothing is typed in from memory. Sources, one per field group:

| Field group | Source | API |
|---|---|---|
| identity, SMILES, XLogP3 | PubChem | PUG REST (`pubchempy` style) and PUG View |
| physicochemical, QED and its eight desirability terms | RDKit from the PubChem SMILES | local |
| routes of administration | openFDA drug label, DailyMed SPL fallback | `api.fda.gov`, `dailymed.nlm.nih.gov` |
| in vitro activity | ChEMBL pChEMBL at the human H2 receptor (CHEMBL1941) | ChEMBL REST |
| in vitro toxicity | ChEMBL pChEMBL at the strongest non-H2 target, hERG called out | ChEMBL REST |
| in vivo toxicity | rat acute oral LD50, Therapeutics Data Commons `LD50_Zhu` (ChemIDplus/RTECS origin) | Harvard Dataverse |
| bioavailability, half-life, dose, dose interval | the product label for that route | openFDA label, PubChem PUG View |
| primary toxicity endpoint | most frequently reported FAERS reaction (a reporting signal, not a rate) | openFDA event |

Rules kept from the project instructions: a drug may carry at most two non-intravenous routes (rules 3 and 4);
intravenous rows are added on top of that cap so dose-corrected bioavailability can be compared across routes
(rule 4 lets IV be excluded from the two); every ATC code x route pair is an independent row (rule 5); the
therapeutic index and the free drug concentration are written as spreadsheet formulas, not as pasted values.

    python h2/fetch_features.py --out h2/data --cache .h2_cache
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
from pathlib import Path
from statistics import median

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, QED

from common import ATC_LEVELS, H2_DRUGS, H2_TARGET, get, set_cache, write_provenance

RDLogger.DisableLog("rdApp.*")

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUGVIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/"
OPENFDA = "https://api.fda.gov/drug/"
LD50_FILE = "https://dataverse.harvard.edu/api/access/datafile/4267146"  # TDC LD50_Zhu (acute oral rat LD50)
HERG = "CHEMBL240"
IV = {"INTRAVENOUS", "INTRAVENOUS BOLUS", "INTRAVENOUS DRIP"}

# ---------------------------------------------------------------- PubChem


def pubchem_identity(cid: int) -> dict:
    r = get(f"{PUG}/compound/cid/{cid}/property/IsomericSMILES,CanonicalSMILES,MolecularWeight,XLogP,"
            f"IUPACName,InChIKey,Title/JSON", note=f"PubChem identity CID {cid}")
    p = r.json()["PropertyTable"]["Properties"][0] if r else {}
    smiles = p.get("IsomericSMILES") or p.get("SMILES") or p.get("ConnectivitySMILES") or p.get("CanonicalSMILES")
    return {"smiles_pubchem": smiles, "xlogp3": p.get("XLogP"), "inchikey_pubchem": p.get("InChIKey"),
            "title": p.get("Title"), "mw_pubchem": p.get("MolecularWeight")}


def pugview_sections(cid: int, headings: list[str]) -> dict[str, list[str]]:
    """Plain-text statements under the given PUG View headings, keyed by heading."""
    out: dict[str, list[str]] = {h: [] for h in headings}
    r = get(f"{PUGVIEW}/{cid}/JSON", note=f"PubChem PUG View CID {cid}")
    if not r:
        return out

    def walk(sections, path=()):
        for s in sections:
            p = path + (s.get("TOCHeading", ""),)
            for info in s.get("Information", []):
                for sm in info.get("Value", {}).get("StringWithMarkup", []):
                    yield p, sm["String"]
            yield from walk(s.get("Section", []), p)

    for path, text in walk(r.json()["Record"]["Section"]):
        for h in headings:
            if h in path:
                out[h].append(re.sub(r"\s+", " ", text).strip())
    return out


# ---------------------------------------------------------------- RDKit descriptors


def descriptors(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    props = QED.properties(mol)
    return {
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "cLogP (Crippen)": round(Crippen.MolLogP(mol), 2),
        "Hydrogen Bond Donors": props.HBD,
        "Hydrogen Bond Acceptors": props.HBA,
        "Topological Surface Area": round(props.PSA, 2),
        "Rotatable Bonds": props.ROTB,
        "Aromatic Rings": props.AROM,
        "Molecular Refractivity": round(Crippen.MolMR(mol), 2),
        "Structural Alerts (QED ALERTS)": props.ALERTS,
        "Quantitative Estimate of Drug Likeness": round(QED.qed(mol), 3),
        "QED weights max": round(QED.weights_max(mol), 3),
        "QED weights none": round(QED.weights_none(mol), 3),
        "Heavy Atoms": mol.GetNumHeavyAtoms(),
        "Fraction Csp3": round(Descriptors.FractionCSP3(mol), 3),
        "Formal Charge": Chem.GetFormalCharge(mol),
        "InChIKey": Chem.MolToInchiKey(mol),
        "Canonical SMILES (RDKit)": Chem.MolToSmiles(mol),
    }


# ---------------------------------------------------------------- ChEMBL


def chembl_molecule(name: str, inchikey: str = "") -> dict:
    """Resolve the ChEMBL parent: by InChIKey first, by preferred name only if the structure is not indexed.

    The name lookup alone is not enough: "NIZATIDINE" matches a record that carries no activities, while the
    InChIKey finds the entry the assays hang off.
    """
    hits = []
    if inchikey:
        r = get(CHEMBL + "molecule.json", {"molecule_structures__standard_inchi_key": inchikey,
                                           "only": "molecule_chembl_id,molecule_properties,molecule_hierarchy,pref_name"},
                note=f"ChEMBL molecule by InChIKey {inchikey}")
        hits = r.json().get("molecules", []) if r else []
    if not hits:
        r = get(CHEMBL + "molecule.json", {"pref_name__iexact": name.upper(),
                                           "only": "molecule_chembl_id,molecule_properties,molecule_hierarchy,pref_name"},
                note=f"ChEMBL molecule by name {name}")
        hits = r.json().get("molecules", []) if r else []
    if not hits:
        return {}
    m = hits[0]
    parent = (m.get("molecule_hierarchy") or {}).get("parent_chembl_id") or m["molecule_chembl_id"]
    props = m.get("molecule_properties") or {}
    return {"chembl_id": parent, "chembl_alogp": props.get("alogp"), "chembl_qed": props.get("qed_weighted"),
            "chembl_ro5_violations": props.get("num_ro5_violations")}


AFFINITY_TYPES = {"IC50", "Ki", "Kd", "EC50", "AC50"}


def pchembl_at(chembl_id: str, target: str) -> dict:
    """Median pChEMBL for one molecule at one target.

    Where ChEMBL stores the measurement but no pChEMBL (nizatidine's H2 entries are a censored AC50 of 10 uM),
    the value is derived from the nanomolar standard_value and the row says so.
    """
    r = get(CHEMBL + "activity.json",
            {"parent_molecule_chembl_id": chembl_id, "target_chembl_id": target, "limit": 1000,
             "only": "pchembl_value,standard_type,standard_value,standard_units,standard_relation,"
                     "assay_chembl_id,document_chembl_id"},
            note=f"ChEMBL activities {chembl_id} at {target}")
    acts = r.json().get("activities", []) if r else []
    scored = [a for a in acts if a.get("pchembl_value")]
    derived = False
    if not scored:
        scored = [a for a in acts if a.get("standard_value") and a.get("standard_units") == "nM"
                  and a.get("standard_type") in AFFINITY_TYPES]
        derived = bool(scored)
        for a in scored:
            import math
            a["pchembl_value"] = -math.log10(float(a["standard_value"]) * 1e-9)
    if not scored:
        return {"value": None, "n_assays": 0, "n_docs": 0, "types": "", "derived": False}
    return {"value": round(median(float(a["pchembl_value"]) for a in scored), 2),
            "n_assays": len({a["assay_chembl_id"] for a in scored}),
            "n_docs": len({a["document_chembl_id"] for a in scored}),
            "types": ",".join(sorted({a["standard_type"] for a in scored})),
            "derived": derived,
            "relations": ",".join(sorted({a.get("standard_relation") or "=" for a in scored}))}


H2_ORTHOLOGS = {"CHEMBL1941": "human", "CHEMBL2882": "guinea pig", "CHEMBL2245": "mouse", "CHEMBL3892": "rat"}


def h2_potency(chembl_id: str) -> dict:
    """pChEMBL at the human H2 receptor; if ChEMBL holds none, the same measurement at an ortholog, flagged."""
    for target, species in H2_ORTHOLOGS.items():
        act = pchembl_at(chembl_id, target)
        if act["value"] is not None:
            act["target"], act["species"] = target, species
            return act
    return {"value": None, "n_assays": 0, "n_docs": 0, "types": "", "target": H2_TARGET, "species": "human"}


def strongest_off_target(chembl_id: str) -> dict:
    """Highest median pChEMBL at any target that is not the H2 receptor: the in vitro off-target liability."""
    r = get(CHEMBL + "activity.json",
            {"parent_molecule_chembl_id": chembl_id, "pchembl_value__isnull": "false", "limit": 1000,
             "only": "pchembl_value,target_chembl_id,assay_chembl_id,document_chembl_id,standard_type"},
            note=f"ChEMBL all activities {chembl_id}")
    acts = [a for a in (r.json().get("activities", []) if r else []) if a["target_chembl_id"] != H2_TARGET]
    if not acts:
        return {}
    by_target: dict[str, list] = {}
    for a in acts:
        by_target.setdefault(a["target_chembl_id"], []).append(a)
    scored = {t: (median(float(a["pchembl_value"]) for a in v), len({a["assay_chembl_id"] for a in v}),
                  len({a["document_chembl_id"] for a in v})) for t, v in by_target.items()}
    top = max(scored, key=lambda t: scored[t][0])
    name = target_name(top)
    return {"off_target": top, "off_target_name": name, "off_target_pchembl": round(scored[top][0], 2),
            "off_target_n_assays": scored[top][1], "off_target_n_docs": scored[top][2],
            "n_off_targets": len(scored)}


def target_name(target_chembl_id: str) -> str:
    r = get(CHEMBL + f"target/{target_chembl_id}.json", {"only": "pref_name,organism"},
            note=f"ChEMBL target {target_chembl_id}")
    if not r:
        return ""
    d = r.json()
    return f"{d.get('pref_name')} ({d.get('organism')})"


# ---------------------------------------------------------------- in vivo acute toxicity


def ld50_table() -> dict[str, float]:
    """Rat acute oral LD50 in log(1/(mol/kg)) from Therapeutics Data Commons, keyed by InChIKey skeleton."""
    r = get(LD50_FILE, note="TDC LD50_Zhu acute oral rat toxicity")
    out: dict[str, float] = {}
    if not r:
        return out
    for row in csv.DictReader(io.StringIO(r.text), delimiter="\t"):
        mol = Chem.MolFromSmiles(row["X"])
        if mol is not None:
            out[Chem.MolToInchiKey(mol).split("-")[0]] = float(row["Y"])
    return out


# ---------------------------------------------------------------- labels: routes, dose, PK


def label_fields(generic: str, route: str | None = None) -> dict:
    """openFDA label record(s) for one generic name, optionally restricted to one route."""
    q = f'openfda.generic_name:"{generic}"'
    if route:
        q += f' AND openfda.route:"{route}"'
    r = get(OPENFDA + "label.json", {"search": q, "limit": 5}, note=f"openFDA label {generic} {route or ''}")
    return r.json() if r else {}


def route_counts(generic: str) -> tuple[dict[str, int], str]:
    r = get(OPENFDA + "label.json", {"search": f'openfda.generic_name:"{generic}"', "count": "openfda.route.exact"},
            note=f"openFDA route counts {generic}")
    if r:
        return {x["term"].upper(): x["count"] for x in r.json().get("results", [])}, "openFDA drug label"
    spl = get("https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json", {"drug_name": generic, "pagesize": 50},
              note=f"DailyMed SPL list {generic}")
    counts: dict[str, int] = {}
    for rec in (spl.json().get("data", []) if spl else [])[:40]:
        x = get(f"https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/{rec['setid']}.xml",
                note=f"DailyMed SPL {generic} {rec['setid']}")
        for route in set(re.findall(r'<routeCode[^>]*displayName="([^"]+)"', x.text if x else "")):
            counts[route.upper()] = counts.get(route.upper(), 0) + 1
    return counts, "DailyMed SPL"


DOSE_RE = re.compile(r"(\d+(?:\.\d+)?)\s*mg(?:\s*/\s*(\d+(?:\.\d+)?)\s*mL)?[^.]{0,80}?"
                     r"(once daily|twice daily|two times a day|at bedtime|every\s+\d+\s+hours|q\s*\d+\s*h|daily)",
                     re.I)
INTERVAL = {"once daily": 24, "daily": 24, "at bedtime": 24, "twice daily": 12, "two times a day": 12}


def dose_from_label(payload: dict) -> dict:
    """First dose + interval the label's dosage section states, with the sentence kept as evidence."""
    for rec in payload.get("results", []):
        for field in ("dosage_and_administration", "indications_and_usage"):
            for text in rec.get(field, []):
                for sentence in re.split(r"(?<=[.;])\s+", text):
                    m = DOSE_RE.search(sentence)
                    if m:
                        phrase = m.group(3).lower()
                        hours = INTERVAL.get(phrase)
                        if hours is None:
                            n = re.search(r"\d+", phrase)
                            hours = float(n.group()) if n else None
                        return {"dose_mg": float(m.group(1)), "dose_interval_h": hours,
                                "dose_evidence": sentence.strip()[:300]}
    return {}


HALF_LIFE_RE = re.compile(r"half[- ]life[^.]{0,120}?(\d+(?:\.\d+)?)\s*(?:to|-|–)?\s*(\d+(?:\.\d+)?)?\s*hours?", re.I)
BIOAVAIL_RE = re.compile(r"bioavailability[^.]{0,120}?(\d+(?:\.\d+)?)\s*(?:to|-|–)?\s*(\d+(?:\.\d+)?)?\s*%", re.I)


def _range_mid(m) -> float:
    lo = float(m.group(1))
    hi = float(m.group(2)) if m.group(2) else lo
    return round((lo + hi) / 2, 2)


def pk_from_text(statements: list[str], pattern: re.Pattern) -> dict:
    """Median over every statement that reports the quantity, so one unusual sentence cannot set the value."""
    found = [(_range_mid(m), s) for s in statements if (m := pattern.search(s))]
    if not found:
        return {}
    values = sorted(v for v, _ in found)
    value = round(median(values), 2)
    nearest = min(found, key=lambda vs: abs(vs[0] - value))
    return {"value": value, "n": len(found), "all": values, "evidence": nearest[1][:300]}


WHO_ADM_R = {"O": "ORAL", "P": "PARENTERAL", "N": "NASAL", "R": "RECTAL", "V": "VAGINAL", "SL": "SUBLINGUAL",
             "TD": "TRANSDERMAL", "Inhal": "RESPIRATORY (INHALATION)", "Implant": "IMPLANTATION",
             "Instill.sol": "INSTILLATION", "s.c.": "SUBCUTANEOUS", "Chewing gum": "ORAL",
             "Lamella": "OPHTHALMIC", "oral aerosol": "RESPIRATORY (INHALATION)"}


def who_ddd(atc: str) -> dict:
    """WHO ATC/DDD index row for one code: the defined daily dose and the route it is defined for.

    Used where no US label exists (roxatidine is not marketed in the United States), so the route and an
    effective-dose proxy still come from a named source rather than from nowhere.
    """
    r = get("https://atcddd.fhi.no/atc_ddd_index/", {"code": atc, "showdescription": "no"},
            note=f"WHO ATC/DDD index {atc}")
    if not r:
        return {}
    import html as _html
    text = re.sub(r"\s+", " ", _html.unescape(re.sub(r"<[^>]+>", " ", r.content.decode("latin-1"))))
    m = re.search(rf"{atc}\s+([a-z0-9 ,.-]+?)\s+([\d.]+)\s+(g|mg|mcg|U|TU|MU|ml)\s+([A-Za-z.]+)", text)
    if not m:
        return {}
    grams = {"g": 1, "mg": 1e-3, "mcg": 1e-6}.get(m.group(3))
    return {"ddd_mg": round(float(m.group(2)) * grams * 1000, 3) if grams else None, "ddd_unit": m.group(3),
            "route": WHO_ADM_R.get(m.group(4).strip()), "who_adm_r": m.group(4).strip(),
            "evidence": f"WHO ATC/DDD index {atc}: DDD {m.group(2)} {m.group(3)}, Adm.R {m.group(4).strip()}"}


def faers_top_reaction(generic: str) -> dict:
    r = get(OPENFDA + "event.json", {"search": f'patient.drug.openfda.generic_name:"{generic}"',
                                     "count": "patient.reaction.reactionmeddrapt.exact", "limit": 5},
            note=f"openFDA FAERS reactions {generic}")
    res = r.json().get("results", []) if r else []
    if not res:
        return {}
    return {"faers_top_reaction": res[0]["term"].title(), "faers_top_count": res[0]["count"],
            "faers_next": "; ".join(f'{x["term"].title()} ({x["count"]})' for x in res[1:4])}


# ---------------------------------------------------------------- build


def build(out_dir: Path) -> pd.DataFrame:
    ld50 = ld50_table()
    rows: list[dict] = []
    for drug in H2_DRUGS:
        name, atc, cid = drug["name"], drug["atc"], drug["cid"]
        ident = pubchem_identity(cid)
        smiles = ident["smiles_pubchem"]
        desc = descriptors(smiles)
        ch = chembl_molecule(name, desc["InChIKey"])
        act = h2_potency(ch["chembl_id"]) if ch else {"value": None, "n_assays": 0, "n_docs": 0, "species": "human",
                                                      "target": H2_TARGET, "types": ""}
        off = strongest_off_target(ch["chembl_id"]) if ch else {}
        herg = pchembl_at(ch["chembl_id"], HERG) if ch else {}
        sections = pugview_sections(cid, ["Biological Half-Life", "Absorption, Distribution and Excretion"])
        half = pk_from_text(sections["Biological Half-Life"], HALF_LIFE_RE)
        bioav = pk_from_text(sections["Absorption, Distribution and Excretion"], BIOAVAIL_RE)
        faers = faers_top_reaction(name)
        counts, route_src = route_counts(name)
        ddd = who_ddd(atc)
        legend_routes = {k: v for k, v in counts.items() if k}
        non_iv = sorted(((k, v) for k, v in legend_routes.items() if k not in IV), key=lambda kv: -kv[1])[:2]
        iv = [(k, v) for k, v in legend_routes.items() if k in IV]
        chosen = non_iv + iv                                            # rules 3 and 4: two non-IV, IV on top
        if not chosen and ddd.get("route"):                             # no US label (roxatidine): WHO ATC/DDD route
            chosen, route_src = [(ddd["route"], 0)], ddd["evidence"]
        if not chosen:
            chosen = [(None, 0)]
        ld50_log = ld50.get(desc["InChIKey"].split("-")[0])
        ld50_mg_kg = round(10 ** (-ld50_log) * 1000 * desc["Molecular Weight"], 1) if ld50_log else None

        for route, n_labels in chosen:
            label = label_fields(name, route) if route else {}
            dose = dose_from_label(label)
            # per-route bioavailability: 100% for intravenous by definition, else the label/PUG View statement
            if route in IV:
                f_abs, f_src = 1.0, "intravenous: complete by definition"
            elif bioav:
                f_abs, f_src = round(bioav["value"] / 100, 3), f'label statement: "{bioav["evidence"]}"'
            else:
                f_abs, f_src = None, "left empty: no bioavailability statement in the label or PubChem PUG View"
            if not dose and ddd.get("ddd_mg"):      # no label dose: the WHO defined daily dose for this code
                dose = {"dose_mg": ddd["ddd_mg"], "dose_interval_h": 24, "dose_evidence": ddd["evidence"]}
            dose_mg_kg = round(dose["dose_mg"] / 70, 4) if dose.get("dose_mg") else None
            rows.append({
                "Compound Identifier": cid,
                "Generic Name": name,
                "Anatomical(Organ System)": atc[0],
                "Anatomical name": ATC_LEVELS[atc[0]],
                "Therapeutic(Indication)": atc[1:3],
                "Therapeutic name": ATC_LEVELS[atc[:3]],
                "Mechanism of Action": atc[3],
                "Mechanism name": ATC_LEVELS[atc[:4]],
                "Chemical Group(Pharmacophore)": atc[4],
                "Chemical group name": ATC_LEVELS[atc[:5]],
                "Anatomical Therapeutic Chemical Code": atc,
                **desc,
                "Lipophilicity(XLogP3, PubChem)": ident["xlogp3"],
                "Lipophilicity_Edit_Flag": "XLogP3 from PubChem; log D7.4 is not published by any open API, "
                                           "cLogP (Crippen) is given beside it",
                "ChEMBL alogp": ch.get("chembl_alogp"),
                "ChEMBL QED": ch.get("chembl_qed"),
                "Mode of Administration": route or "NaN",
                "Mode of Administration_source": f"{route_src}: {n_labels} labels list it; all counts {counts}"
                                                 if route else f"{route_src}: no route found",
                "Route rule": ("intravenous row kept on top of the two-route cap so dose-corrected bioavailability "
                               "can be compared across routes (rule 4)") if route in IV else
                              "one of the two most reported non-intravenous routes (rules 3 and 4)",
                "Structure": smiles,
                "SMILES_source": f"PubChem CID {cid} (isomeric), canonicalised with RDKit",
                "Primary Biological Target": "Histamine H2 receptor (CHEMBL1941, UniProt P25021)",
                "In Vitro Activity (pChEMBL at H2)": act.get("value"),
                "In Vitro Activity_source": f'ChEMBL {ch.get("chembl_id")} at {act.get("target")} '
                                            f'({act.get("species")} H2 receptor): median of {act.get("n_assays")} '
                                            f'assays in {act.get("n_docs")} documents ({act.get("types")}'
                                            f'{", derived from standard_value, no pChEMBL in ChEMBL" if act.get("derived") else ""}'
                                            f'{", relation " + act["relations"] if act.get("relations") not in (None, "=") else ""})'
                                            if act.get("value") else
                                            "left empty: ChEMBL holds no pChEMBL value at the H2 receptor of any "
                                            "species for this molecule",
                "In Vitro Toxicity (strongest off-target pChEMBL)": off.get("off_target_pchembl"),
                "In Vitro Toxicity_target": f'{off.get("off_target_name")} [{off.get("off_target")}]' if off else "",
                "In Vitro Toxicity_source": f'ChEMBL: median of {off.get("off_target_n_assays")} assays in '
                                            f'{off.get("off_target_n_docs")} documents; strongest of '
                                            f'{off.get("n_off_targets")} non-H2 targets' if off else
                                            "left empty: no off-target pChEMBL value in ChEMBL",
                "hERG pChEMBL": herg.get("value"),
                "hERG_source": f'ChEMBL {HERG}: {herg.get("n_assays")} assays in {herg.get("n_docs")} documents'
                               if herg.get("value") else "left empty: no hERG pChEMBL value in ChEMBL",
                "In Vivo Toxicity (rat oral LD50, mg/kg)": ld50_mg_kg,
                "In Vivo Toxicity_source": "Therapeutics Data Commons LD50_Zhu (acute oral rat LD50, ChemIDplus/RTECS "
                                           f"origin): y = {ld50_log} in log(1/(mol/kg))" if ld50_log else
                                           "left empty: this drug is not in the TDC acute oral rat LD50 set",
                "In Vivo Activity (1/effective dose, kg/mg)": round(1 / dose_mg_kg, 3) if dose_mg_kg else None,
                "In Vivo Activity_source": f'label dose {dose.get("dose_mg")} mg per 70 kg adult: '
                                           f'"{dose.get("dose_evidence")}"' if dose_mg_kg else
                                           "left empty: no dose statement matched in the label for this route",
                "Bioavailability(Fraction Absorbed)": f_abs,
                "Bioavailability_source": f_src,
                "H2 potency species": act.get("species"),
                "Half Life (h)": half.get("value"),
                "Half Life_source": f'PubChem PUG View: median of {half.get("n")} statements {half.get("all")}; '
                                    f'nearest "{half.get("evidence")}"' if half else
                                    "left empty: no half-life statement found",
                "Dosage (mg)": dose.get("dose_mg"),
                "Dose Interval (h)": dose.get("dose_interval_h"),
                "Dosage_source": f'openFDA label for the {route} product: "{dose.get("dose_evidence")}"'
                                 if dose else "left empty: no dose and interval matched in the label",
                "Primary Toxicity Endpoint": faers.get("faers_top_reaction"),
                "Primary Toxicity Endpoint_source": f'openFDA FAERS reporting signal: {faers.get("faers_top_count")} '
                                                    f'reports; next {faers.get("faers_next")}' if faers else
                                                    "left empty: no FAERS reports returned",
            })

    df = pd.DataFrame(rows)
    df = df.sort_values(["Anatomical(Organ System)", "Therapeutic(Indication)", "Mechanism of Action",
                         "Chemical Group(Pharmacophore)", "Anatomical Therapeutic Chemical Code",
                         "Mode of Administration"]).reset_index(drop=True)

    # Spreadsheet formulas, written against the row each value lands on (header is row 1, data starts at row 2).
    tox_col, act_col = col_letter(df, "In Vivo Toxicity (rat oral LD50, mg/kg)"), col_letter(df, "In Vivo Activity (1/effective dose, kg/mg)")
    f_col, t_col = col_letter(df, "Bioavailability(Fraction Absorbed)"), col_letter(df, "Half Life (h)")
    d_col, i_col = col_letter(df, "Dosage (mg)"), col_letter(df, "Dose Interval (h)")
    df["Therapeutic Index"] = [f"=IF(OR({tox_col}{i + 2}=\"\",{act_col}{i + 2}=\"\"),\"\","
                               f"LOG({tox_col}{i + 2})-LOG(1/{act_col}{i + 2}))" for i in range(len(df))]
    df["In Vivo Free Drug Concentration"] = [f"=IF(COUNT({f_col}{i + 2},{t_col}{i + 2},{d_col}{i + 2},"
                                             f"{i_col}{i + 2})<4,\"\",1.44*{f_col}{i + 2}*{t_col}{i + 2}*"
                                             f"{d_col}{i + 2}/{i_col}{i + 2})" for i in range(len(df))]
    df["Therapeutic Index_value"] = [ti_value(r) for _, r in df.iterrows()]
    df["Hammet Substituent Constant"] = None
    df["Hammet_Edit_Flag"] = ("left empty: a Hammett sigma is defined for a substituent on a common scaffold, and "
                              "these five drugs do not share one; TPSA is the polarity term used in the Hansch fits")
    return df


def ti_value(row) -> float | None:
    import math
    tox, act = row["In Vivo Toxicity (rat oral LD50, mg/kg)"], row["In Vivo Activity (1/effective dose, kg/mg)"]
    if tox is None or act is None or pd.isna(tox) or pd.isna(act):
        return None
    return round(math.log10(tox) - math.log10(1 / act), 3)


def none_if_nan(v):
    return None if v is None or (isinstance(v, float) and pd.isna(v)) or pd.isna(v) else v


def col_letter(df: pd.DataFrame, name: str) -> str:
    from openpyxl.utils import get_column_letter
    return get_column_letter(list(df.columns).index(name) + 1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="h2/data")
    ap.add_argument("--cache", default=".h2_cache")
    args = ap.parse_args()
    set_cache(args.cache)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    df = build(out)
    df.to_csv(out / "h2_features.csv", index=False)
    write_provenance(out / "provenance_features.csv")
    panel = {}
    for _, r in df.iterrows():
        panel.setdefault(r["Generic Name"], []).append({
            "route": r["Mode of Administration"], "atc": r["Anatomical Therapeutic Chemical Code"],
            "qed": none_if_nan(r["Quantitative Estimate of Drug Likeness"]),
            "pchembl": none_if_nan(r["In Vitro Activity (pChEMBL at H2)"]),
            "ld50": none_if_nan(r["In Vivo Toxicity (rat oral LD50, mg/kg)"]),
            "ti": none_if_nan(r["Therapeutic Index_value"]),
            "bioavailability": none_if_nan(r["Bioavailability(Fraction Absorbed)"]),
            "half_life": none_if_nan(r["Half Life (h)"]),
            "logp": none_if_nan(r["cLogP (Crippen)"]), "tpsa": none_if_nan(r["Topological Surface Area"]),
            "mw": none_if_nan(r["Molecular Weight"]),
        })
    (out / "h2_features.json").write_text(json.dumps(panel))   # the binding page's side panel reads this
    filled = {c: int(df[c].notna().sum()) for c in df.columns if not c.endswith(("_source", "_Edit_Flag", "rule"))}
    (out / "fill_report.json").write_text(json.dumps({"rows": len(df), "filled": filled}, indent=1))
    print(df[["Generic Name", "Anatomical Therapeutic Chemical Code", "Mode of Administration",
              "Quantitative Estimate of Drug Likeness", "In Vitro Activity (pChEMBL at H2)",
              "In Vivo Toxicity (rat oral LD50, mg/kg)", "Therapeutic Index_value"]].to_string(index=False))


if __name__ == "__main__":
    main()
