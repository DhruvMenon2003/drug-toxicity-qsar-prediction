# Curate the DrugBank dataset in one run (Google Colab or local Python).
#
# In Colab:
#   !pip -q install rdkit
#   !wget -q -O colab_curate_all.py https://raw.githubusercontent.com/DhruvMenon2003/drug-toxicity-qsar-prediction/main/curation/colab_curate_all.py
#   %run colab_curate_all.py
#
# It asks you to upload "Drug Bank Dataset.tsv.md" if that file is not already in the working folder.
# Everything else (mode-of-administration table, FDA route legend) is read from this GitHub repo.
#
# Output (folder curated/):
#   drugbank_curated.csv      one row per drug x ATC code x route, sorted by ATC level 1 > 2 > 3 > code
#   structures_<ATC2>.png     2D structures of each ATC level-2 group, to check by eye
#   same_smiles.png           drugs that share a structure (should be empty)
#   Activity                  median pChEMBL at the drug's mechanism target (ChEMBL), only when >= 2 documents agree
#   provenance.csv            which API answered each lookup, with URL and time

import io, re, time, json, os, datetime as dt
from collections import Counter
import requests, pandas as pd

REPO = "https://raw.githubusercontent.com/DhruvMenon2003/drug-toxicity-qsar-prediction/main/curation/inputs/"
DRUGBANK_FILE = "Drug Bank Dataset.tsv.md"
OUT = "curated"
DROP_EXACT_DUPLICATES = True    # Star, 2026-10-09: the source file repeats its table; keep one copy of each identical row
os.makedirs(OUT, exist_ok=True)
try:   # optional openFDA key from Colab Secrets (key icon in the left bar, name FDA_API_KEY); never paste it into code
    from google.colab import userdata
    os.environ.setdefault("OPENFDA_API_KEY", userdata.get("FDA_API_KEY") or "")
except Exception:
    pass
os.environ.setdefault("OPENFDA_API_KEY", os.environ.get("FDA_API_KEY", ""))
S = requests.Session(); S.headers["User-Agent"] = "drug-curation (research; contact via GitHub)"
PROV = []

def get(url, **kw):
    for i in range(4):
        try:
            r = S.get(url, timeout=90, **kw)
            if r.status_code in (429, 500, 502, 503, 504) and i < 3:   # busy server: wait and retry
                time.sleep(2 ** (i + 1)); continue
            PROV.append({"url": re.sub(r"api_key=[^&]+", "api_key=hidden", r.url), "status": r.status_code, "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")})
            return r
        except requests.RequestException:
            time.sleep(2 ** i)
    return None

def read_md_table(text):
    rows = [l.strip()[1:-1].strip().split(" | ") for l in text.splitlines() if l.startswith("|") and not l.startswith("| ---")]
    df = pd.DataFrame(rows[1:], columns=[c.strip() for c in rows[0]])
    return df[df.iloc[:, 0] != df.columns[0]]          # drop repeated header lines

# ---------- inputs ----------
if not os.path.exists(DRUGBANK_FILE):
    from google.colab import files
    files.upload()
db = read_md_table(open(DRUGBANK_FILE, encoding="utf-8").read())
mode = read_md_table(get(REPO + "mode_of_administration.md").text)
legend = pd.read_csv(io.StringIO(get(REPO + "fda_route_of_administration.tsv").text), sep="\t", dtype=str, keep_default_na=False)
who2fda = pd.read_csv(io.StringIO(get(REPO + "who_to_fda_route.csv").text), dtype=str).set_index("who_admr").fda_route
L = legend.set_index("NAME")
print(f"DrugBank rows {len(db)} (exact repeats {db.duplicated().sum()}), mode rows {len(mode)}, legend routes {len(legend)}")
if DROP_EXACT_DUPLICATES:
    db = db.drop_duplicates()

CID, NAME, ATCCOL = "Compound Identifier", "Generic Name", "Anatomical Therapeutic Chemical Code"

# ---------- deliverable 2: one row per ATC code, salt-free names ----------
db = db.assign(ATC=db[ATCCOL].str.split("|")).explode("ATC")
# Only these name parts are counter-ions. Everything else after a "." is part of the drug name
# (acids, esters/prodrugs, insulin analogues, combinations, radiopharmaceuticals) and is left as it is.
COUNTER_IONS = {"bromide", "isethionate", "acetate", "sodium", "calcium"}
KEEP_WHOLE = {"sodium.bicarbonate", "sodium.phosphate"}          # here the sodium salt is the drug itself
def salt_free(name):
    if name in KEEP_WHOLE or "." not in name:
        return name, ""
    parts = [p for p in name.split(".") if p not in COUNTER_IONS]
    new = ".".join(parts)
    return (new, f"salt removed: {name} -> {new}") if new != name else (name, "")
db[[NAME, "Name_Edit_Flag"]] = db[NAME].apply(lambda n: pd.Series(salt_free(n)))

# ---------- deliverable 1: ATC levels 1-3 with WHO names ----------
atc_names = {}
for l3 in sorted(db.ATC.str[:4].unique()):
    r = get(f"https://atcddd.fhi.no/atc_ddd_index/?code={l3}&showdescription=no")
    for code, label in re.findall(r'<a href="\./\?code=([A-Z0-9]+)&showdescription=(?:no|yes)">([^<]+)</a>', r.text if r else ""):
        atc_names.setdefault(code, label.strip())
for n, k in [(1, 1), (2, 3), (3, 4)]:
    db[f"ATC_L{n}"] = db.ATC.str[:k]
    db[f"ATC_L{n}_name"] = db[f"ATC_L{n}"].map(atc_names)

# ---------- deliverable 4 (fetched first, used again in deliverable 3): isomeric canonical SMILES (PubChem by CID, canonicalised with RDKit) ----------
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.MolStandardize import rdMolStandardize
parent, neutral = rdMolStandardize.LargestFragmentChooser(), rdMolStandardize.Uncharger()
smiles, inchikey, title = {}, {}, {}
cids = sorted(db[CID].unique(), key=int)
for i in range(0, len(cids), 100):
    r = get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/" + ",".join(cids[i:i + 100]) + "/property/SMILES,Title/JSON")
    for p in (r.json()["PropertyTable"]["Properties"] if r is not None and r.ok else []):
        title[str(p["CID"])] = p.get("Title", "")
        m = Chem.MolFromSmiles(p.get("SMILES", ""))
        if m:                                   # strip counter-ions, neutralise, canonicalise (stereo kept)
            m = neutral.uncharge(parent.choose(m))
            smiles[str(p["CID"])], inchikey[str(p["CID"])] = Chem.MolToSmiles(m), Chem.MolToInchiKey(m)
# ---------- deliverable 3: Mode of Administration from the FDA route legend ----------
IV = {"INTRAVENOUS", "INTRAVENOUS BOLUS", "INTRAVENOUS DRIP"}
# Labels are not tagged with ATC codes, so for a drug with several ATC codes a label route only counts
# for the ATC group whose body site it fits (rule 4). Groups not listed accept every route in the legend.
ATC_ROUTES = {
    "S01": {"OPHTHALMIC", "INTRAOCULAR", "CONJUNCTIVAL", "INTRAVITREAL", "INTRACAMERAL", "SUBCONJUNCTIVAL", "RETROBULBAR"},
    "S02": {"AURICULAR (OTIC)"}, "S03": {"OPHTHALMIC", "AURICULAR (OTIC)"},
    "R01": {"NASAL"}, "R03": {"RESPIRATORY (INHALATION)", "ORAL", "NASAL"},
    "D": {"TOPICAL", "CUTANEOUS", "TRANSDERMAL", "INTRADERMAL", "INTRALESIONAL", "SUBCUTANEOUS"},
    "G01": {"VAGINAL"}, "A01": {"DENTAL", "ORAL", "BUCCAL", "OROPHARYNGEAL", "SUBLINGUAL", "TOPICAL"},
}
def allowed_for(atc, n_atc):
    if n_atc == 1:
        return None
    return ATC_ROUTES.get(atc[:3]) or ATC_ROUTES.get(atc[:1])

def label_route_counts(drug, us_name=""):
    """Route counts over all US labels: openFDA first, DailyMed SPL XML if openFDA has none.
    WHO names (aciclovir, indometacin) are retried under the PubChem title, which is usually the US name."""
    for q in dict.fromkeys([drug.replace(".", " "), us_name.lower()]):
        if q:
            c, src = _route_counts(q)
            if c:
                return c, src + ("" if q == drug.replace(".", " ") else f" (searched as {q})")
    return Counter(), "no US label found"

def _route_counts(q):
    params = {"search": f'openfda.generic_name:"{q}"', "count": "openfda.route.exact"}
    if os.environ.get("OPENFDA_API_KEY"):                       # optional; without a key openFDA allows 1,000 calls a day
        params["api_key"] = os.environ["OPENFDA_API_KEY"]
    r = get("https://api.fda.gov/drug/label.json", params=params)
    if r is not None and r.ok:
        return Counter({x["term"]: x["count"] for x in r.json()["results"]}), "openFDA drug label"
    r = get("https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json", params={"drug_name": q, "pagesize": 50})
    c = Counter()
    for spl in (r.json().get("data", []) if r is not None and r.ok else []):
        x = get(f"https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/{spl['setid']}.xml")
        if x is not None and x.ok:
            c.update(set(re.findall(r'<routeCode[^>]*displayName="([^"]+)"', x.text)))
    return c, "DailyMed SPL"

n_atc = db.groupby(CID).ATC.nunique()
cache, rows = {}, []
for (cid, atc), grp in mode.groupby(["CID", "ATC"], sort=False):
    known = [a for a in grp.AdmR if a != "NaN"]
    drug = grp.GenericName.iloc[0]
    if known:   # WHO codes already in the table -> FDA legend names
        for a in dict.fromkeys(known):
            rows.append({"CID": cid, "ATC": atc, "route": who2fda.get(a), "AdmR_Edit_Flag": f"WHO code {a} mapped to FDA legend", "AdmR_evidence": "source table"})
        continue
    if drug not in cache:
        cache[drug] = label_route_counts(drug, title.get(cid, ""))
    counts, src = cache[drug]
    ok = allowed_for(atc, n_atc.get(cid, 1))
    counts = Counter({k.upper(): v for k, v in counts.items() if k.upper() in L.index and (ok is None or k.upper() in ok)})
    top = Counter({k: v for k, v in counts.items() if k not in IV}).most_common(2)          # rules 3 and 4
    flag = f"filled from {src}"
    if not top and counts:                                                               # IV-only drug (decided 2026-10-09)
        top, flag = [("INTRAVENOUS", sum(v for k, v in counts.items() if k in IV))], flag + "; iv_only"
    if not top:
        rows.append({"CID": cid, "ATC": atc, "route": None, "AdmR_Edit_Flag": f"left empty: {src}" + (" has no route fitting this ATC group" if src != "no US label found" else ""), "AdmR_evidence": ""})
    for route, n in top:                                                                 # rule 5: one row per route
        rows.append({"CID": cid, "ATC": atc, "route": route, "AdmR_Edit_Flag": flag, "AdmR_evidence": f"{n} labels list {route}; all counts {dict(counts)}"})
adm = pd.DataFrame(rows).drop_duplicates(["CID", "ATC", "route"])
adm["AdmR_fda_code"] = adm.route.map(L["FDA CODE"])
adm["AdmR_nci_concept_id"] = adm.route.map(L["NCI CONCEPT ID"])

db = db.merge(adm, left_on=[CID, "ATC"], right_on=["CID", "ATC"], how="left").drop(columns="CID")
db["Mode of Administration"] = db.route.fillna("NaN")

db["2D Chemical Structure"] = db[CID].map(smiles).fillna("NaN")
db["SMILES_source"] = db[CID].map(lambda c: "PubChem CID " + c + " (isomeric; salt stripped, neutralised, RDKit canonical)" if smiles.get(c) else "skipped: no small-molecule structure in PubChem (biologic, polymer or mixture)")
db["InChIKey"] = db[CID].map(inchikey)
same = db.dropna(subset=["InChIKey"]).groupby("InChIKey")[NAME].nunique()
db["SMILES_shared_with_other_drug"] = db.InChIKey.isin(same[same > 1].index)   # flagged only, no rows removed

# ---------- Activity: ChEMBL potency (pChEMBL = -log10 molar IC50/Ki/EC50/Kd) at the drug's mechanism target ----------
# Mechanism and activities are read for the parent compound and all its salt forms.
# Filled only when at least two independent assays (different ChEMBL documents) report a value; otherwise left
# empty and the single finding is noted as unconfirmed. Value = median pChEMBL over all those measurements.
from concurrent.futures import ThreadPoolExecutor
from statistics import median
CH = "https://www.ebi.ac.uk/chembl/api/data/"
keys = sorted({k for k in inchikey.values() if k})
chembl_id = {}
for i in range(0, len(keys), 40):
    r = get(CH + "molecule.json", params={"molecule_structures__standard_inchi_key__in": ",".join(keys[i:i + 40]), "limit": 1000,
                                          "only": "molecule_chembl_id,molecule_structures,molecule_hierarchy"})
    for m in (r.json()["molecules"] if r is not None and r.ok else []):
        par = (m.get("molecule_hierarchy") or {}).get("parent_chembl_id") or m["molecule_chembl_id"]
        chembl_id[m["molecule_structures"]["standard_inchi_key"]] = par
names = db.drop_duplicates(CID).set_index(CID)[NAME]
for c, k in inchikey.items():          # structure not matched exactly (tautomer, stereo): try the drug name
    if k not in chembl_id and c in names:
        r = get(CH + "molecule.json", params={"pref_name__iexact": names[c].replace(".", " "), "only": "molecule_chembl_id,molecule_hierarchy"})
        hit = r.json()["molecules"] if r is not None and r.ok else []
        if hit:
            chembl_id[k] = (hit[0].get("molecule_hierarchy") or {}).get("parent_chembl_id") or hit[0]["molecule_chembl_id"]

def potency(mol):
    r = get(CH + "mechanism.json", params={"parent_molecule_chembl_id": mol, "only": "target_chembl_id,mechanism_of_action", "limit": 100})
    mech = [x for x in (r.json()["mechanisms"] if r is not None and r.ok else []) if x["target_chembl_id"]]
    if not mech:
        return {"Activity_Edit_Flag": f"left empty: no mechanism target in ChEMBL ({mol})"}
    tg = sorted({x["target_chembl_id"] for x in mech})
    r = get(CH + "activity.json", params={"parent_molecule_chembl_id": mol, "target_chembl_id__in": ",".join(tg), "pchembl_value__isnull": "false",
                                          "limit": 1000, "only": "pchembl_value,standard_type,assay_chembl_id,document_chembl_id,target_chembl_id"})
    acts = r.json()["activities"] if r is not None and r.ok else []
    docs = {a["document_chembl_id"] for a in acts}
    out = {"Activity_target": "; ".join(f'{x["target_chembl_id"]} ({x["mechanism_of_action"]})' for x in mech),
           "Activity_n_measurements": len(acts), "Activity_n_documents": len(docs),
           "Activity_types": ", ".join(sorted({a["standard_type"] for a in acts}))}
    if len(docs) >= 2 and len({a["assay_chembl_id"] for a in acts}) >= 2:
        out.update(Activity=round(median(float(a["pchembl_value"]) for a in acts), 2), Activity_Edit_Flag=f"median pChEMBL, ChEMBL {mol}")
    else:
        out["Activity_Edit_Flag"] = (f"left empty: unconfirmed, one report only (pChEMBL {acts[0]['pchembl_value']})" if acts
                                     else "left empty: no pChEMBL value at the mechanism target")
    return out

mols = sorted(set(chembl_id.values()))
with ThreadPoolExecutor(8) as ex:
    pot = dict(zip(mols, ex.map(potency, mols)))
act = pd.DataFrame(columns=[CID, "Activity"]) if not inchikey else pd.DataFrame([{CID: c, "Activity_chembl_id": chembl_id.get(k), **(pot.get(chembl_id.get(k)) or
                     {"Activity_Edit_Flag": "left empty: structure not found in ChEMBL"})} for c, k in inchikey.items()])
db = db.drop(columns="Activity").merge(act, on=CID, how="left")
db["Activity"] = db["Activity"].astype(object).where(db["Activity"].notna(), "NaN")
db["Activity_Edit_Flag"] = db["Activity_Edit_Flag"].fillna("left empty: no small-molecule structure")

# ---------- group by ATC level 1 > 2 > 3 and save ----------
db = db.sort_values(["ATC_L1", "ATC_L2", "ATC_L3", "ATC", NAME, "Mode of Administration"]).reset_index(drop=True)
cols = [CID, NAME, "Name_Edit_Flag", "ATC", "ATC_L1", "ATC_L1_name", "ATC_L2", "ATC_L2_name", "ATC_L3", "ATC_L3_name",
        "2D Chemical Structure", "InChIKey", "SMILES_source", "SMILES_shared_with_other_drug", "Biological Target",
        "Mode of Administration", "AdmR_fda_code", "AdmR_nci_concept_id", "AdmR_Edit_Flag", "AdmR_evidence",
        "Activity", "Activity_target", "Activity_types", "Activity_n_measurements", "Activity_n_documents",
        "Activity_chembl_id", "Activity_Edit_Flag", "Toxicity", "Therapeutic Index", "Highest Frequency Adverse Drug Reaction"]
db[cols].to_csv(f"{OUT}/drugbank_curated.csv", index=False)
pd.DataFrame(PROV).to_csv(f"{OUT}/provenance.csv", index=False)

# ---------- validation you can see ----------
uniq = db.drop_duplicates([CID, "ATC"])
for l2, g in uniq.groupby("ATC_L2"):
    g = g[g["2D Chemical Structure"] != "NaN"].drop_duplicates(CID).head(120)
    if len(g):
        Draw.MolsToGridImage([Chem.MolFromSmiles(s) for s in g["2D Chemical Structure"]], molsPerRow=6, subImgSize=(220, 180),
                             legends=[f"{n} {a}" for n, a in zip(g[NAME], g.ATC)]).save(f"{OUT}/structures_{l2}.png")
dup = uniq[uniq.SMILES_shared_with_other_drug].drop_duplicates(CID)
if len(dup):
    Draw.MolsToGridImage([Chem.MolFromSmiles(s) for s in dup["2D Chemical Structure"]], molsPerRow=6,
                         legends=list(dup[NAME])).save(f"{OUT}/same_smiles.png")

print(f"\nrows {len(db)} | ATC codes {db.ATC.nunique()} | L1/L2/L3 groups {db.ATC_L1.nunique()}/{db.ATC_L2.nunique()}/{db.ATC_L3.nunique()}"
      f" | missing ATC names {db[['ATC_L1_name', 'ATC_L2_name', 'ATC_L3_name']].isna().any(axis=1).sum()}")
print(f"names changed {(db.Name_Edit_Flag != '').sum()} rows: {sorted(set(db.loc[db.Name_Edit_Flag != '', 'Name_Edit_Flag']))}")
print("Activity filled:", (db.Activity != "NaN").sum(), "rows;", db.loc[db.Activity == "NaN", "Activity_Edit_Flag"].str.split(r" \(|,|:").str[1].str.strip().value_counts().to_dict())
print("Mode of Administration:", db["Mode of Administration"].value_counts().head(12).to_dict())
print("AdmR still empty:", (db["Mode of Administration"] == "NaN").sum(), "rows | SMILES missing:",
      (db["2D Chemical Structure"] == "NaN").sum(), "rows | drugs sharing a structure:", dup[NAME].tolist())
try:
    from IPython.display import Image, display
    display(db[cols].head(15))
    display(Image(f"{OUT}/structures_{uniq.ATC_L2.iloc[0]}.png"))
except Exception:
    pass
