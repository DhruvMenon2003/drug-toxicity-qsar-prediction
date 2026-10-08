"""Merge spider outputs into the curated TSV: one row per drug x mode of administration.

Inputs (any may be missing; columns are then left blank and the gap is logged in Evidence_Conflicts):
  output/pubchem.jsonl  output/drugbank.jsonl  output/chembl.jsonl  output/openfda.jsonl  output/whocc.jsonl
  inputs/admetsar3_results.csv   (optional: batch export from admetsar3 web server, keyed by SMILES or name)

Orthogonal-assay rules (from the curation brief):
  * two assays per value, SAME units, NEVER mixing in vitro with in vivo
  * "orthogonal" = different ChEMBL assay AND different source document (independent experiments);
    for Activity we also prefer different readouts (e.g. Ki binding vs IC50 functional)
  * in vivo toxicity is matched to the row's route (oral LD50 for the oral row, i.v. LD50 for the i.v. row)

Usage: python merge.py   -> output/ATC_A_curated.tsv
"""
import json
import math
import os
import re
from collections import defaultdict
from urllib.parse import urlparse

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "output")

COLUMNS = [
    "Compound_Identifier", "Generic_Name", "ATC_Code", "DrugBank_ID", "ChEMBL_ID", "Mode_of_Administration",
    "SMILES_PubChem", "Molecular_Formula", "MW_PubChem", "MW_ChEMBL_FreeBase", "Mono_IsotopicMass",
    "IUPACName_PubChem", "XLogP_PubChem", "ALogP_ChEMBL", "LogP_DrugBank_Experimental",
    "Melting_Point_C", "Solubility_mgL", "Solubility_Raw_DrugBank",
    "DrugBank_Kingdom", "DrugBank_Superclass", "DrugBank_Class", "DrugBank_Subclass", "DrugBank_Direct_Parent",
    "DrugBank_Molecular_Framework", "Pharmacophore_DrugBank_Substituents",
    "Biological_Target", "Mechanism_of_Action_ChEMBL",
    "First_Approval_Year", "Withdrawn_Flag", "Black_Box_Warning", "Orphan_Drug", "Ro5_Violations", "Prodrug",
    "Toxicity_Context", "Toxicity_Assay1_Value", "Toxicity_Assay1_Type", "Toxicity_Assay1_Source",
    "Toxicity_Assay2_Value", "Toxicity_Assay2_Type", "Toxicity_Assay2_Source", "Toxicity_Units",
    "Toxicity_DrugBank_Text",
    "Activity_Context", "Activity_Assay1_Value", "Activity_Assay1_Type", "Activity_Assay1_Source",
    "Activity_Assay2_Value", "Activity_Assay2_Type", "Activity_Assay2_Source", "Activity_Units",
    "Therapeutic_Index_Computed", "Therapeutic_Index_Basis",
    "Highest_Freq_ADR_FAERS", "FAERS_Reports_For_Top_ADR", "FAERS_Total_Route_Reports",
    "WHO_DDD", "WHO_DDD_Unit",
    "PK_Absorption", "PK_Volume_of_Distribution", "PK_Protein_Binding", "PK_Half_Life", "PK_Clearance",
    "PK_Metabolism", "PK_Route_of_Elimination", "PD_Pharmacodynamics",
    "ADMET_Physchem_RDKit", "ADMET_ADMETsar3",
    "Evidence_Conflicts", "Filtered_Domains_Log",
]

ROUTE_NORMAL = {
    "oral": "Oral", "intravenous": "Intravenous", "subcutaneous": "Subcutaneous",
    "intramuscular": "Intramuscular", "rectal": "Rectal", "topical": "Topical", "transdermal": "Transdermal",
    "sublingual": "Sublingual", "buccal": "Buccal", "ophthalmic": "Ophthalmic", "nasal": "Nasal",
    "intranasal": "Nasal", "parenteral": "Parenteral", "inhalation": "Inhalation", "respiratory": "Inhalation",
    "intraperitoneal": "Intraperitoneal", "dental": "Dental", "vaginal": "Vaginal",
}
ROUTE_PATTERNS = [
    ("Oral", r"\b(oral(ly)?|p\.?\s?o\.?|per os|gavage)\b"),
    ("Intravenous", r"\b(intravenous(ly)?|i\.?\s?v\.?)\b"),
    ("Intraperitoneal", r"\b(intraperitoneal(ly)?|i\.?\s?p\.?)\b"),
    ("Subcutaneous", r"\b(subcutaneous(ly)?|s\.?\s?c\.?)\b"),
    ("Intramuscular", r"\b(intramuscular(ly)?|i\.?\s?m\.?)\b"),
    ("Rectal", r"\brectal(ly)?\b"),
    ("Topical", r"\b(topical(ly)?|dermal)\b"),
]
# MedDRA PTs that describe product use, efficacy or outcome rather than an adverse reaction; skipped
# when picking the highest-frequency ADR (they are still in output/openfda.jsonl).
NON_ADR_TERMS = {
    "DRUG INEFFECTIVE", "OFF LABEL USE", "PRODUCT USE IN UNAPPROVED INDICATION", "PRODUCT USE ISSUE",
    "INTENTIONAL PRODUCT USE ISSUE", "INTENTIONAL PRODUCT MISUSE", "DRUG INTERACTION", "NO ADVERSE EVENT",
    "INAPPROPRIATE SCHEDULE OF PRODUCT ADMINISTRATION", "WRONG TECHNIQUE IN PRODUCT USAGE PROCESS",
    "PRODUCT DOSE OMISSION ISSUE", "PRODUCT DOSE OMISSION", "INCORRECT DOSE ADMINISTERED", "DRUG DOSE OMISSION",
    "TREATMENT FAILURE", "THERAPY NON-RESPONDER", "CONDITION AGGRAVATED", "DISEASE PROGRESSION", "DEATH",
    "MALIGNANT NEOPLASM PROGRESSION", "TOXICITY TO VARIOUS AGENTS", "OVERDOSE", "EXPOSURE DURING PREGNANCY",
    "MATERNAL EXPOSURE DURING PREGNANCY", "THERAPEUTIC RESPONSE UNEXPECTED", "PRODUCT QUALITY ISSUE",
    "PRODUCT SUBSTITUTION ISSUE", "DRUG HYPERSENSITIVITY", "ILLNESS", "HOSPITALISATION",
}
# Indication terms by ATC prefix: confounding by indication, not ADRs; skipped and noted.
INDICATION_TERMS = {
    "A02": {"GASTROOESOPHAGEAL REFLUX DISEASE", "DYSPEPSIA", "GASTRIC ULCER", "DUODENAL ULCER"},
    "A03": {"IRRITABLE BOWEL SYNDROME", "ABDOMINAL PAIN"}, "A04": {"NAUSEA", "VOMITING"},
    "A06": {"CONSTIPATION"}, "A07E": {"COLITIS ULCERATIVE", "CROHN'S DISEASE", "RHEUMATOID ARTHRITIS"},
    "A08": {"WEIGHT INCREASED", "OBESITY"},
    "A10": {"BLOOD GLUCOSE INCREASED", "DIABETES MELLITUS", "TYPE 2 DIABETES MELLITUS",
            "DIABETES MELLITUS INADEQUATE CONTROL", "GLYCOSYLATED HAEMOGLOBIN INCREASED"},
}
# Known stimulated-reporting (litigation) signals in FAERS; flagged, not removed.
LITIGATION_SIGNALS = {("A02B", "CHRONIC KIDNEY DISEASE"), ("A02B", "ACUTE KIDNEY INJURY"),
                      ("A02B", "END STAGE RENAL DISEASE"), ("A10BK", "DIABETIC KETOACIDOSIS"),
                      ("A10BK", "FOURNIER'S GANGRENE"), ("A10BG", "CARDIAC FAILURE CONGESTIVE"),
                      ("A10BG", "BLADDER CANCER")}
IN_VIVO_TYPES = {"LD50", "TD50", "MTD", "NOAEL", "LOAEL", "ED50", "LC50"}
IN_VITRO_TOX_TYPES = {"CC50", "TC50"}


def read_jsonl(name):
    p = os.path.join(OUT_DIR, name)
    if not os.path.exists(p):
        return []
    with open(p, encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def by_key(rows, key="generic_name"):
    d = defaultdict(list)
    for r in rows:
        d[str(r.get(key, "")).lower()].append(r)
    return d


def norm_route(text):
    t = (text or "").strip().lower()
    for k, v in ROUTE_NORMAL.items():
        if k in t:
            return v
    return t.title() if t else None


def route_from_description(desc):
    d = (desc or "").lower()
    for name, pat in ROUTE_PATTERNS:
        if re.search(pat, d):
            return name
    return None


def first_number(s):
    m = re.search(r"-?\d+(?:\.\d+)?", s or "")
    return float(m.group()) if m else None


def solubility_mg_l(raw):
    """Convert DrugBank water-solubility strings to mg/L when the unit is explicit; else None."""
    if not raw:
        return None
    s = raw.replace(",", "").lower()
    m = re.search(r"(\d+(?:\.\d+)?(?:e[-+]?\d+)?)\s*(mg/ml|mg/l|g/l|ug/ml|µg/ml|mcg/ml|g/100\s*ml)", s)
    if not m:
        return None
    v, u = float(m.group(1)), m.group(2)
    factor = {"mg/ml": 1000, "mg/l": 1, "g/l": 1000, "ug/ml": 1, "µg/ml": 1, "mcg/ml": 1}.get(u)
    if factor is None and u.startswith("g/100"):
        factor = 10000
    return v * factor


def exp_prop(db, key):
    for row in db.get("experimental_properties") or []:
        if key in (row.get("Property") or "").lower():
            return row.get("Value"), row.get("Source")
    return None, None


def classify_tox(rec):
    st = (rec.get("standard_type") or "").upper()
    if st in IN_VITRO_TOX_TYPES or rec.get("assay_cell_type"):
        return "in vitro"
    if st in IN_VIVO_TYPES:
        return "in vivo"
    return "in vitro" if rec.get("assay_type") in ("B", "A") else "in vivo"


def src(rec):
    return f"ChEMBL {rec.get('assay_chembl_id')} / {rec.get('document_chembl_id')} ({rec.get('document_year') or 'n.d.'})"


def pick_orthogonal(records, prefer_distinct_type=False):
    """Return up to 2 records from different assays AND different documents (and, if possible, readouts)."""
    recs = [r for r in records if r.get("standard_value") not in (None, "") and not r.get("data_validity_comment")
            and not r.get("potential_duplicate")]
    recs.sort(key=lambda r: (-(r.get("document_year") or 0), r.get("assay_chembl_id") or ""))
    if not recs:
        return []
    first = recs[0]
    def ok(r, strict):
        diff = r["assay_chembl_id"] != first["assay_chembl_id"] and r["document_chembl_id"] != first["document_chembl_id"]
        return diff and (not strict or r.get("standard_type") != first.get("standard_type"))
    for strict in ([True, False] if prefer_distinct_type else [False]):
        for r in recs[1:]:
            if ok(r, strict):
                return [first, r]
    return [first]


def fmt_val(r, pchembl=False):
    if pchembl:
        return f"{float(r['pchembl_value']):.2f}"
    rel = r.get("standard_relation") or "="
    rel = "" if rel == "=" else rel
    return f"{rel}{r['standard_value']}"


def toxicity_for_route(chembl, route):
    tox = [r for r in (chembl or {}).get("toxicity_records", []) if (r.get("standard_type") or "").upper() != "ED50"]
    groups = defaultdict(list)
    for r in tox:
        ctx = classify_tox(r)
        r_route = route_from_description(r.get("assay_description")) if ctx == "in vivo" else "n/a"
        if ctx == "in vivo" and r_route != route:
            continue
        key = (ctx, r.get("standard_type"), r.get("standard_units"), r.get("assay_organism") if ctx == "in vivo" else "")
        groups[key].append(r)
    if not groups:
        return None
    # prefer in vivo route-matched, then largest group
    key = max(groups, key=lambda k: (k[0] == "in vivo", len(groups[k])))
    picks = pick_orthogonal(groups[key])
    return key, picks


def therapeutic_index(chembl, route):
    recs = (chembl or {}).get("toxicity_records", [])
    ld = defaultdict(list); ed = defaultdict(list)
    for r in recs:
        if route_from_description(r.get("assay_description")) != route or r.get("standard_value") in (None, ""):
            continue
        k = (r.get("assay_organism"), r.get("standard_units"))
        st = (r.get("standard_type") or "").upper()
        if st in ("LD50", "TD50"):
            ld[k].append(float(r["standard_value"]))
        elif st == "ED50":
            ed[k].append(float(r["standard_value"]))
    for k in ld:
        if k in ed:
            gm = lambda xs: math.exp(sum(math.log(x) for x in xs if x > 0) / len(xs))
            return round(gm(ld[k]) / gm(ed[k]), 2), f"LD50/ED50 geometric means, {k[0]}, {route}, {k[1]} (n={len(ld[k])}/{len(ed[k])})"
    return None, None


def main():
    drugs = pd.read_csv(os.path.join(HERE, "inputs", "atc_a_small_molecules.csv"), dtype=str)
    pub = by_key(read_jsonl("pubchem.jsonl"))
    # DrugBank, best source first: live scrape > official XML release > browser capture from 2026-10-07
    dbk = by_key(read_jsonl("drugbank.jsonl"), "drugbank_id")
    for fallback in ["drugbank_xml.jsonl"] + sorted(f for f in os.listdir(OUT_DIR) if f.startswith("drugbank_seed")):
        for r in read_jsonl(fallback):
            k = str(r.get("drugbank_id", "")).lower()
            if not any(x.get("status") == "ok" for x in dbk.get(k, [])):
                r.setdefault("capture", fallback)
                dbk[k].append(r)
    chb = by_key(read_jsonl("chembl.jsonl"))
    fda = by_key(read_jsonl("openfda.jsonl"))
    who = by_key(read_jsonl("whocc.jsonl"))
    admetsar = None
    ap = os.path.join(HERE, "inputs", "admetsar3_results.csv")
    if os.path.exists(ap):
        admetsar = pd.read_csv(ap)

    rows = []
    for _, d in drugs.iterrows():
        name = d.generic_name.lower()
        p = next((x for x in pub.get(name, []) if x.get("status") == "ok"), {})
        db = next((x for x in dbk.get(str(d.drugbank_id).lower(), []) if x.get("status") == "ok"), {})
        ch = next((x for x in chb.get(name, []) if x.get("status") == "ok"), {})
        fd = next(iter(fda.get(name, [])), {})
        wh = [w for w in who.get(name, [])]
        conflicts, domains = [], set()
        for src_rec in [p, db, ch, fd, *wh]:
            domains.update(urlparse(u).netloc for u in src_rec.get("source_urls", []) if u)
        for label, rec in (("PubChem", p), ("DrugBank", db), ("ChEMBL", ch)):
            if not rec:
                conflicts.append(f"{label}: no data retrieved")

        # routes: DrugBank dosage forms > WHO DDD Adm.R > ChEMBL flags
        routes = sorted({norm_route(r.get("Route")) for r in db.get("dosage_forms", []) if r.get("Route")} - {None})
        if not routes:
            routes = sorted({norm_route(x["adm_route"]) for w in wh for x in w.get("ddd_rows", [])} - {None})
        if not routes and ch:
            routes = [r for r, f in (("Oral", ch.get("oral")), ("Parenteral", ch.get("parenteral")),
                                     ("Topical", ch.get("topical"))) if f]
        if not routes and fd.get("routes"):
            top_r = norm_route(fd["routes"][0]["route"])
            if top_r and top_r != "Unknown":
                routes = [top_r]
                conflicts.append(f"Route: taken from most-reported FAERS route ({top_r}); confirm on label")
        if not routes:
            routes = ["Not reported"]
            conflicts.append("Route: none found in DrugBank/WHO/ChEMBL")

        mp_raw, _ = exp_prop(db, "melting point")
        sol_raw, _ = exp_prop(db, "water solubility")
        logp_raw, _ = exp_prop(db, "logp")
        sol = solubility_mg_l(sol_raw)
        if sol_raw and sol is None:
            conflicts.append(f"Solubility qualitative only, no mg/L value ('{sol_raw}')")

        subs = list(db.get("substituents") or [])
        if db.get("substituents_truncated_more"):
            subs.append(f"(+ {db['substituents_truncated_more']} more; full list needs live page or XML)")
        if db.get("capture") and db.get("capture") != "drugbank.jsonl":
            conflicts.append(f"DrugBank fields from {db['capture']} (live scrape unavailable)")
        mechs = ch.get("mechanisms") or []
        rd = p.get("rdkit") or {}
        adm_txt = ""
        if admetsar is not None:
            hit = admetsar[(admetsar.astype(str).apply(lambda c: c.str.lower()) == name).any(axis=1)]
            if len(hit):
                adm_txt = "; ".join(f"{k}={v}" for k, v in hit.iloc[0].items())

        # activity: in vitro pChEMBL on mechanism targets (route independent)
        act = pick_orthogonal([r for r in ch.get("activity_records", []) if r.get("pchembl_value")],
                              prefer_distinct_type=True)
        if len(act) < 2:
            conflicts.append(f"Activity: only {len(act)} independent assay(s) in ChEMBL; literature needed")

        for route in routes:
            row = dict.fromkeys(COLUMNS, "")
            row.update({
                "Compound_Identifier": d.pubchem_cid, "Generic_Name": d.generic_name, "ATC_Code": d.atc_code,
                "DrugBank_ID": d.drugbank_id, "ChEMBL_ID": ch.get("chembl_id", ""), "Mode_of_Administration": route,
                "SMILES_PubChem": p.get("smiles", ""), "Molecular_Formula": p.get("molecular_formula", ""),
                "MW_PubChem": p.get("molecular_weight", ""), "MW_ChEMBL_FreeBase": ch.get("mw_freebase", ""),
                "Mono_IsotopicMass": p.get("monoisotopic_mass", ""), "IUPACName_PubChem": p.get("iupac_name", ""),
                "XLogP_PubChem": p.get("xlogp", ""), "ALogP_ChEMBL": ch.get("alogp", ""),
                "LogP_DrugBank_Experimental": logp_raw or "",
                "Melting_Point_C": first_number(mp_raw) if mp_raw else "",
                "Solubility_mgL": sol if sol is not None else "", "Solubility_Raw_DrugBank": sol_raw or "",
                "DrugBank_Kingdom": db.get("kingdom") or "", "DrugBank_Superclass": db.get("super_class") or "",
                "DrugBank_Class": db.get("class") or "", "DrugBank_Subclass": db.get("sub_class") or "",
                "DrugBank_Direct_Parent": db.get("direct_parent") or "",
                "DrugBank_Molecular_Framework": db.get("molecular_framework") or "",
                "Pharmacophore_DrugBank_Substituents": " / ".join(subs),
                "Biological_Target": "; ".join(f"{t.get('Target')} ({t.get('Actions')})" for t in db.get("moa_targets", [])),
                "Mechanism_of_Action_ChEMBL": "; ".join(m["mechanism_of_action"] or "" for m in mechs),
                "First_Approval_Year": ch.get("first_approval") or "", "Withdrawn_Flag": ch.get("withdrawn_flag", ""),
                "Black_Box_Warning": ch.get("black_box_warning", ""), "Orphan_Drug": ch.get("orphan", ""),
                "Ro5_Violations": ch.get("num_ro5_violations", ""), "Prodrug": ch.get("prodrug", ""),
                "Toxicity_DrugBank_Text": db.get("toxicity") or "",
                "PK_Absorption": db.get("absorption") or "", "PK_Volume_of_Distribution": db.get("volume_of_distribution") or "",
                "PK_Protein_Binding": db.get("protein_binding") or "", "PK_Half_Life": db.get("half_life") or "",
                "PK_Clearance": db.get("clearance") or "", "PK_Metabolism": db.get("metabolism") or "",
                "PK_Route_of_Elimination": db.get("route_of_elimination") or "",
                "PD_Pharmacodynamics": db.get("pharmacodynamics") or "",
                "ADMET_Physchem_RDKit": "; ".join(f"{k}={v}" for k, v in rd.items()),
                "ADMET_ADMETsar3": adm_txt,
            })
            row_conf = list(conflicts)

            t = toxicity_for_route(ch, route)
            if t:
                (ctx, stype, units, org), picks = t
                row["Toxicity_Context"] = f"{ctx}; {stype}" + (f"; {org}; {route}" if ctx == "in vivo" else "")
                row["Toxicity_Units"] = units or ""
                for i, r in enumerate(picks, 1):
                    row[f"Toxicity_Assay{i}_Value"] = fmt_val(r)
                    row[f"Toxicity_Assay{i}_Type"] = f"{r.get('standard_type')} - {(r.get('assay_description') or '')[:120]}"
                    row[f"Toxicity_Assay{i}_Source"] = src(r)
                if len(picks) < 2:
                    row_conf.append("Toxicity: only 1 independent assay for this route/context; literature needed")
            else:
                row_conf.append(f"Toxicity: no ChEMBL toxicity assay for route '{route}'")

            if act:
                row["Activity_Context"] = "in vitro; pChEMBL on mechanism target(s)"
                row["Activity_Units"] = "pChEMBL (-log10 M)"
                for i, r in enumerate(act, 1):
                    row[f"Activity_Assay{i}_Value"] = fmt_val(r, pchembl=True)
                    row[f"Activity_Assay{i}_Type"] = f"{r.get('standard_type')} vs {r.get('target_pref_name')} ({r.get('target_organism')})"
                    row[f"Activity_Assay{i}_Source"] = src(r)

            ti, basis = therapeutic_index(ch, route)
            row["Therapeutic_Index_Computed"] = ti if ti is not None else ""
            row["Therapeutic_Index_Basis"] = basis or ("No same-species/same-route LD50+ED50 pair in ChEMBL" if route != "Not reported" else "")

            rx = fd.get("reactions", {})
            route_key = route.lower()
            if route == "Parenteral" and not rx.get(route_key):
                route_key = next((k for k in ("intravenous", "subcutaneous", "intramuscular") if rx.get(k)), route_key)
            ind = set().union(*[v for a, v in INDICATION_TERMS.items()
                                for code in str(d.atc_code).split("|") if code.startswith(a)] or [set()])
            raw = rx.get(route_key) or rx.get("all") or []
            skipped = [r["term"].title() for r in raw[:5] if r["term"].upper() in ind]
            rlist = [r for r in raw if r["term"].upper() not in NON_ADR_TERMS | ind]
            if skipped:
                row_conf.append("FAERS: indication term(s) skipped as confounding by indication: " + ", ".join(skipped))
            if rlist:
                top = rlist[0]
                row["Highest_Freq_ADR_FAERS"] = top["term"].title()
                row["FAERS_Reports_For_Top_ADR"] = top["reports"]
                if not rx.get(route_key):
                    row_conf.append("FAERS: route-specific counts unavailable; all-route counts used")
                elif route_key != route.lower():
                    row_conf.append(f"FAERS: '{route}' row uses {route_key} reports")
                if any(str(d.atc_code).startswith(a) and top["term"].upper() == t for a, t in LITIGATION_SIGNALS):
                    row_conf.append(f"FAERS: '{top['term'].title()}' is a known litigation-stimulated reporting signal")
                row_conf.append("FAERS: report counts (not incidence); drug matched as suspect or concomitant")
            rt = [r for r in fd.get("routes", []) if norm_route(r["route"]) == route]
            row["FAERS_Total_Route_Reports"] = rt[0]["reports"] if rt else ""

            ddd = [x for w in wh for x in w.get("ddd_rows", []) if norm_route(x["adm_route"]) == route
                   or (route in ("Intravenous", "Subcutaneous", "Intramuscular") and x["adm_route"] == "parenteral")]
            if ddd:
                row["WHO_DDD"], row["WHO_DDD_Unit"] = ddd[0]["ddd"], ddd[0]["unit"]

            row["Evidence_Conflicts"] = " | ".join(row_conf)
            row["Filtered_Domains_Log"] = "; ".join(sorted(domains))
            rows.append(row)

    out = pd.DataFrame(rows, columns=COLUMNS)
    path = os.path.join(OUT_DIR, "ATC_A_curated.tsv")
    out.to_csv(path, sep="\t", index=False)
    print(f"wrote {path}: {len(out)} rows x {len(out.columns)} cols")
    fill = (out.replace("", pd.NA).notna().mean() * 100).round(0)
    print(fill.to_string())


if __name__ == "__main__":
    main()
