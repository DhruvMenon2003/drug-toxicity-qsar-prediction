# Fill every feature column for one ATC level-3 chemical group, one row per mode of administration.
#
# The group is set by GROUP below; it ships set to A02BA, the H2-receptor antagonists.
#
# In Colab:
#   !pip -q install rdkit mols2grid
#   !wget -q -O h2_features.py https://raw.githubusercontent.com/DhruvMenon2003/drug-toxicity-qsar-prediction/main/curation/h2_features.py
#   %run h2_features.py
#
# Output (folder h2_features/):
#   features.csv          one row per drug x ATC code x route, columns in the order of the
#                         "Drug Bank Dataset Copy" sheet; route-dependent cells are Excel formulas
#   features_long.csv     the same table melted to drug / route / column / value / source, for auditing
#   provenance.csv        every API call: URL, status, UTC time
#   structures.png        2D structures of the group, drawn from the SMILES that went into the table
#   grid.html             mols2grid card deck: structure + QED + the physicochemical columns
#   qed_distribution.png  QED and each physicochemical property, group values on the whole-dataset spread
#   binding_data.json     ligand 3D coordinates + contact residues for the binding animation
#
# Rules it follows (project instructions, and the drug-data-curation skill):
#   - a cell stays empty unless a source answered for it; the reason goes in Evidence
#   - SMILES come from PubChem only (pubchempy / PUG-REST), never from a label or from memory
#   - in vitro activity = median pChEMBL at the mechanism target, and only when >= 2 documents agree
#   - in vitro toxicity = best pChEMBL at a human off-target (a different protein), same >= 2 rule
#   - at most 2 routes per drug, intravenous excluded from the two and kept as an extra row
#   - Therapeutic Index = log(in vivo toxicity) - log(1 / in vivo activity), written as a formula so
#     the sheet computes it the moment both in vivo numbers are filled

import os, re, io, json, time, statistics, datetime as dt
from collections import Counter, defaultdict

import requests, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Draw, QED, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

GROUP = os.environ.get("ATC_GROUP", "A02BA")          # ATC level-3 chemical group to curate
N_DRUGS = int(os.environ.get("N_DRUGS", "6"))          # at least the first 5 drugs of the group (A02BA05 niperotidine has no ChEMBL mechanism, so 6 gives 5 with a target)
MAX_ROUTES = 2                                          # project rule 3
OUT = os.environ.get("OUT", "h2_features")
REPO = "https://raw.githubusercontent.com/DhruvMenon2003/drug-toxicity-qsar-prediction/main/curation/inputs/"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
os.makedirs(OUT, exist_ok=True)

S = requests.Session()
S.headers["User-Agent"] = "drug-curation (research; contact via GitHub)"
PROV = []


def get(url, **kw):
    """GET with retries; every call is written to provenance.csv."""
    for i in range(4):
        try:
            r = S.get(url, timeout=90, **kw)
            PROV.append({"url": re.sub(r"api_key=[^&]+", "api_key=hidden", r.url),
                         "status": r.status_code,
                         "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")})
            if r.status_code in (429, 500, 502, 503, 504) and i < 3:
                time.sleep(2 ** (i + 1)); continue
            return r
        except requests.RequestException:
            time.sleep(2 ** i)
    return None


def jget(url, **kw):
    r = get(url, **kw)
    if r is None or r.status_code != 200:
        return {}
    try:
        return r.json()
    except ValueError:
        return {}


def read_md_table(text):
    rows = [l.strip()[1:-1].strip().split(" | ") for l in text.splitlines()
            if l.startswith("|") and not l.startswith("| ---")]
    df = pd.DataFrame(rows[1:], columns=[c.strip() for c in rows[0]])
    return df[df.iloc[:, 0] != df.columns[0]]


# ---------------------------------------------------------------- inputs
SHEET = next((p for p in [os.environ.get("DRUGBANK_FILE", ""), "Drug Bank Dataset.tsv.md",
                          "/mnt/project-files/Drug Bank Dataset.tsv.md"] if p and os.path.exists(p)), None)
if SHEET is None:                                       # Colab: ask for the file once
    from google.colab import files
    files.upload(); SHEET = "Drug Bank Dataset.tsv.md"
db = read_md_table(open(SHEET, encoding="utf-8").read())
db.columns = [c.strip() for c in db.columns]
NAME = next(c for c in db.columns if "Generic" in c)
ATC = next(c for c in db.columns if "Anatomical Therapeutic" in c or c.strip() == "ATC")
CID = db.columns[0]

legend = pd.read_csv(io.StringIO(get(REPO + "fda_route_of_administration.tsv").text),
                     sep="\t", dtype=str, keep_default_na=False)
nci2route = {r["NCI CONCEPT ID"].strip(): r[" NAME"].strip() if " NAME" in r else r["NAME"].strip()
             for _, r in legend.iterrows()}

# the drugs of the group, in the order the sheet lists them, de-duplicated
rows = []
seen = set()
for _, r in db.iterrows():
    for code in str(r[ATC]).split("|"):
        code = code.strip()
        if code.startswith(GROUP) and (r[NAME], code) not in seen:
            seen.add((r[NAME], code))
            rows.append({"cid_sheet": str(r[CID]).strip(), "name": str(r[NAME]).strip().split(".")[0], "atc": code})
# the sheet does not always hold the whole chemical group; top the list up from the WHO ATC/DDD index,
# which is the group's own authority, so that at least N_DRUGS drugs of the same target are curated
r = get(f"https://atcddd.fhi.no/atc_ddd_index/?code={GROUP}&showdescription=no")
atc_index = []
if r is not None and r.status_code == 200:
    flat = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", r.text).replace("&nbsp;", " "))
    blocks = re.split(rf"(?={GROUP}\d\d )", flat)
    for b in blocks[1:]:
        code = b[:len(GROUP) + 2]
        rest = b[len(GROUP) + 2:].split("List of abbreviations")[0]
        nm = re.split(r"\d", rest, maxsplit=1)[0].strip(" ,")       # the name ends where the DDD number starts
        if not nm or "combination" in nm or any(nm == x["name"] for x in atc_index):
            continue
        who_admr = re.findall(r"\b(O|P|SL|R|N|V|TD|Inhal[\w.]*|urethral|implant)\b", rest[:160])
        atc_index.append({"cid_sheet": "", "name": nm, "atc": code, "from_index": True,
                          "who_admr": who_admr})
WHO_ADMR = {x["atc"]: x.get("who_admr", []) for x in atc_index}
who2fda = pd.read_csv(io.StringIO(get(REPO + "who_to_fda_route.csv").text), dtype=str) \
            .set_index("who_admr").fda_route.to_dict()
for extra in atc_index:
    if len(rows) >= N_DRUGS:
        break
    if not any(extra["atc"] == x["atc"] for x in rows):
        rows.append(extra)
rows = sorted(rows, key=lambda x: x["atc"])[:max(N_DRUGS, min(len(rows), N_DRUGS))]
print(f"{GROUP}: curating {len(rows)} drugs -> {[r['name'] for r in rows]}")


# ---------------------------------------------------------------- PubChem: identity and SMILES
def pubchem(name):
    """CID + canonical SMILES + XLogP3 from PubChem only (project deliverable 4)."""
    j = jget(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(name)}"
             f"/property/MolecularFormula,MolecularWeight,ConnectivitySMILES,SMILES,XLogP,Title/JSON")
    props = (j.get("PropertyTable") or {}).get("Properties") or []
    if not props:
        return {}
    p = props[0]
    return {"cid": p.get("CID"),
            "smiles": p.get("SMILES") or p.get("ConnectivitySMILES"),
            "xlogp3": p.get("XLogP"),
            "pubchem_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{p.get('CID')}"}


# ---------------------------------------------------------------- RDKit: physicochemical block
def descriptors(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {}
    return {"MW": round(Descriptors.MolWt(m), 2),
            "clogp": round(Crippen.MolLogP(m), 2),
            "HBD": rdMolDescriptors.CalcNumHBD(m),
            "HBA": rdMolDescriptors.CalcNumHBA(m),
            "TPSA": round(rdMolDescriptors.CalcTPSA(m), 2),
            "RotB": rdMolDescriptors.CalcNumRotatableBonds(m),
            "ArRings": rdMolDescriptors.CalcNumAromaticRings(m),
            "MR": round(Crippen.MolMR(m), 2),
            "QED": round(QED.qed(m), 3)}


# ---------------------------------------------------------------- ChEMBL: target, activity, off-target
def chembl_molecule(name):
    j = jget(f"{CHEMBL}/molecule/search", params={"q": name, "format": "json", "limit": 5})
    for mol in j.get("molecules", []):
        if (mol.get("pref_name") or "").upper() == name.upper():
            return mol["molecule_chembl_id"]
    return j.get("molecules", [{}])[0].get("molecule_chembl_id") if j.get("molecules") else None


def chembl_mechanism(mol_id):
    j = jget(f"{CHEMBL}/mechanism", params={"molecule_chembl_id": mol_id, "format": "json", "limit": 20})
    out = []
    for m in j.get("mechanisms", []):
        if m.get("target_chembl_id"):
            out.append({"target": m["target_chembl_id"], "action": m.get("action_type"),
                        "mechanism": m.get("mechanism_of_action")})
    return out


def target_name(tid):
    j = jget(f"{CHEMBL}/target/{tid}", params={"format": "json"})
    return j.get("pref_name"), j.get("organism")


def potency(mol_id, target_id):
    """median pChEMBL at one target; returned only when >= 2 documents agree (orthogonal-evidence rule)."""
    j = jget(f"{CHEMBL}/activity", params={"molecule_chembl_id": mol_id, "target_chembl_id": target_id,
                                           "pchembl_value__isnull": "false", "format": "json", "limit": 500})
    acts = [a for a in j.get("activities", []) if a.get("pchembl_value")]
    docs = {a.get("document_chembl_id") for a in acts}
    types = {a.get("standard_type") for a in acts}
    if len(docs) < 2:
        return None, {"n": len(acts), "docs": len(docs), "types": sorted(t for t in types if t),
                      "reason": "fewer than 2 documents"}
    vals = [float(a["pchembl_value"]) for a in acts]
    return round(statistics.median(vals), 2), {"n": len(acts), "docs": len(docs),
                                               "types": sorted(t for t in types if t),
                                               "range": [round(min(vals), 2), round(max(vals), 2)]}


def off_target(mol_id, on_target_ids):
    """best human off-target pChEMBL: the in vitro toxicity proxy the sheet asks for."""
    j = jget(f"{CHEMBL}/activity", params={"molecule_chembl_id": mol_id, "pchembl_value__isnull": "false",
                                           "format": "json", "limit": 1000})
    by_target = defaultdict(list)
    for a in j.get("activities", []):
        tid = a.get("target_chembl_id")
        if not tid or tid in on_target_ids or "Histamine H2" in (a.get("target_pref_name") or ""):
            continue
        by_target[(tid, a.get("target_pref_name"))].append((float(a["pchembl_value"]), a.get("document_chembl_id")))
    best = None
    for (tid, tname), vals in by_target.items():
        if len({d for _, d in vals}) < 2:          # same >= 2 document rule as on-target
            continue
        med = statistics.median([v for v, _ in vals])
        if best is None or med > best[0]:
            best = (round(med, 2), tname, tid, len(vals), len({d for _, d in vals}))
    if best is None:
        return None, {"reason": "no human off-target with >= 2 documents",
                      "targets_seen": len(by_target)}
    return best[0], {"target": best[1], "target_id": best[2], "n": best[3], "docs": best[4]}


# ---------------------------------------------------------------- DailyMed: routes and label pharmacokinetics
PK = {
    "bioavailability": re.compile(r"(?:absolute\s+)?bioavailabilit\w*[^.]{0,120}?(\d{1,3})\s*(?:to|-|–)?\s*(\d{1,3})?\s*%", re.I),
    "half_life":       re.compile(r"half[- ]?life[^.]{0,120}?(\d+(?:\.\d+)?)\s*(?:to|-|–)?\s*(\d+(?:\.\d+)?)?\s*hour", re.I),
}


def dailymed(name):
    """routes (FDA route codes in the SPL) and, where the label states them, F and half-life."""
    j = jget("https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json",
             params={"drug_name": name, "pagesize": 100})
    setids = [d["setid"] for d in j.get("data", [])]
    routes, pk, pk_src = Counter(), {}, {}
    for setid in setids[:25]:
        r = get(f"https://dailymed.nlm.nih.gov/dailymed/services/v2/spls/{setid}.xml")
        if r is None or r.status_code != 200:
            continue
        xml = r.text
        for code in set(re.findall(r'<routeCode[^>]*code="(C\d+)"', xml)):
            if code in nci2route:
                routes[nci2route[code]] += 1
        if len(pk) < 2:
            text = re.sub(r"<[^>]+>", " ", xml)
            for field, rx in PK.items():
                if field in pk:
                    continue
                m = rx.search(text)
                if m:
                    pk[field] = m.group(1) + (f"-{m.group(2)}" if m.group(2) else "")
                    pk_src[field] = f"https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid={setid}"
    return routes, pk, pk_src, len(setids)


def pick_routes(routes):
    """project rules 3 and 4: the 2 most reported routes, intravenous excluded from the two."""
    ranked = [r for r, _ in routes.most_common()]
    iv = [r for r in ranked if r == "INTRAVENOUS"]
    non_iv = [r for r in ranked if r != "INTRAVENOUS"][:MAX_ROUTES]
    return non_iv + iv, iv


# ---------------------------------------------------------------- build the table
COLUMNS = ["Compound Identifier", "Generic Name", "Anatomical(Organ System)", "Therapeutic(Indication)",
           "Mechanism of Action", "Chemical Group(Pharmacophore)", "Anatomical Therapeutic Chemical Code",
           "Molecular Weight", "Lipophilicity(log D)", "Hydrogen Bond Donors", "Hydrogen Bond Acceptors",
           "Topological Surface Area", "Rotatable Bonds", "Aromatic Rings", "Hammet Substituent Constant",
           "Molecular Refractivity Or Taft Steric Factor", "Quantitative Estimate of Drug Likeness",
           "Mode of Administration", "Mode of Administration2", "Structure", "Primary Biological Target",
           "In Vitro Toxicity(Binding affinity to Off Target)", "In Vitro Activity",
           "In Vivo Toxicity(Cell Culture)", "In Vivo Activity(Therapeutic Response)", "Therapeutic Index",
           "Bioavailability(Fraction Absorbed)", "Half Life", "Dosage", "Dose Interval",
           "In Vivo Free Drug Concentration", "Primary Toxicity Endpoint", "Evidence", "Sources"]

table, long_rows, binding = [], [], {}

# pass 1: identity, structure, descriptors and the ChEMBL mechanism target of each drug
for item in rows:
    item["pc"] = pubchem(item["name"])
    if not item["pc"].get("smiles"):
        print(f"  {item['name']}: no PubChem record, skipped"); continue
    item["d"] = descriptors(item["pc"]["smiles"])
    item["mol_id"] = chembl_molecule(item["name"])
    item["mechs"] = chembl_mechanism(item["mol_id"]) if item["mol_id"] else []

# an ATC level-3 chemical group is a group of drugs on the same target, so a drug whose ChEMBL record
# carries no mechanism row borrows the target the rest of the group agrees on, and says so
votes = Counter(m["target"] for item in rows for m in item.get("mechs", [])[:1])
CONSENSUS = votes.most_common(1)[0][0] if votes else None
if CONSENSUS:
    print(f"  group target: {CONSENSUS} ({votes[CONSENSUS]} of {len(rows)} drugs name it in ChEMBL)")

# pass 2: potency at that target, off-target potency, routes, and one row per route
for item in rows:
    name, atc, pc = item["name"], item["atc"], item.get("pc") or {}
    if not pc.get("smiles"):
        continue
    d, mol_id, mechs = item["d"], item["mol_id"], item["mechs"]
    on_ids = {m["target"] for m in mechs} | ({CONSENSUS} if CONSENSUS else set())
    target_id = mechs[0]["target"] if mechs else CONSENSUS
    target_from = "ChEMBL mechanism" if mechs else ("ATC group consensus" if CONSENSUS else "none")
    tname = torg = None
    if target_id:
        tname, torg = target_name(target_id)
    act, act_note = potency(mol_id, target_id) if (mol_id and target_id) else (None, {"reason": "no ChEMBL molecule or target"})
    act_note["target_from"] = target_from
    tox, tox_note = off_target(mol_id, on_ids) if mol_id else (None, {"reason": "no ChEMBL molecule"})
    routes, pk, pk_src, n_labels = dailymed(name)
    chosen, iv = pick_routes(routes)
    route_source = "DailyMed labels"
    if not chosen:                                      # not marketed in the US: fall back to the WHO Adm.R letters
        chosen = [who2fda[a] for a in WHO_ADMR.get(atc, []) if a in who2fda][:MAX_ROUTES]
        route_source = "WHO ATC/DDD Adm.R (no DailyMed label)"
    if not chosen:
        chosen = [""]
        route_source = "none: no DailyMed label and no WHO Adm.R"

    evidence = ["Row from the WHO ATC/DDD index (not in the sheet's A02BA rows)" if item.get("from_index")
                else "Row from the Drug Bank Dataset sheet"]
    if act is None:
        evidence.append(f"In vitro activity empty: {act_note.get('reason')} "
                        f"(target {tname or '-'} from {target_from})")
    else:
        evidence.append(f"In vitro activity = median pChEMBL {act} at {tname} "
                        f"({act_note['n']} values, {act_note['docs']} documents, {'/'.join(act_note['types'])}; "
                        f"target from {target_from})")
    if tox is None:
        evidence.append(f"In vitro toxicity empty: {tox_note.get('reason')}")
    else:
        evidence.append(f"In vitro toxicity = median pChEMBL {tox} at off-target {tox_note['target']} "
                        f"({tox_note['n']} values, {tox_note['docs']} documents)")
    evidence.append("In vivo toxicity and in vivo activity empty: no two orthogonal in vivo assays found by API; "
                    "Therapeutic Index left as a formula so it fills itself once both are entered")
    evidence.append("Hammett sigma empty: defined for a substituent on a congeneric series, not for a whole drug")
    evidence.append(f"Routes from {route_source}" + (f" ({n_labels} labels): " if routes else ": ") +
                    (", ".join(f"{r} x{c}" for r, c in routes.most_common()) or ", ".join(chosen)))

    for rank, route in enumerate(chosen, start=1):
        excel = len(table) + 3                       # the sheet's first data row is 3
        rec = {
            "Compound Identifier": pc["cid"], "Generic Name": name,
            "Anatomical(Organ System)": atc[0], "Therapeutic(Indication)": atc[1:3],
            "Mechanism of Action": atc[3], "Chemical Group(Pharmacophore)": atc[4],
            "Anatomical Therapeutic Chemical Code": atc,
            "Molecular Weight": d.get("MW"),
            "Lipophilicity(log D)": pc.get("xlogp3") if pc.get("xlogp3") is not None else d.get("clogp"),
            "Hydrogen Bond Donors": d.get("HBD"), "Hydrogen Bond Acceptors": d.get("HBA"),
            "Topological Surface Area": d.get("TPSA"), "Rotatable Bonds": d.get("RotB"),
            "Aromatic Rings": d.get("ArRings"), "Hammet Substituent Constant": "",
            "Molecular Refractivity Or Taft Steric Factor": d.get("MR"),
            "Quantitative Estimate of Drug Likeness": d.get("QED"),
            "Mode of Administration": route,
            "Mode of Administration2": ("intravenous, extra row (rule 3)" if route == "INTRAVENOUS"
                                        else f"route {rank} of {len([r for r in chosen if r != 'INTRAVENOUS'])}"),
            "Structure": pc["smiles"], "Primary Biological Target": f"{tname} ({torg})" if tname else "",
            "In Vitro Toxicity(Binding affinity to Off Target)": tox if tox is not None else "",
            "In Vitro Activity": act if act is not None else "",
            "In Vivo Toxicity(Cell Culture)": "", "In Vivo Activity(Therapeutic Response)": "",
            "Therapeutic Index": f"=IF(COUNT(X{excel},Y{excel})=2,LOG10(X{excel})-LOG10(1/Y{excel}),\"\")",
            "Bioavailability(Fraction Absorbed)": pk.get("bioavailability", "") if rank == 1 else "",
            "Half Life": pk.get("half_life", ""), "Dosage": "", "Dose Interval": "",
            "In Vivo Free Drug Concentration": f"=IF(COUNT(AA{excel},AB{excel},AC{excel},AD{excel})=4,1.44*AA{excel}*AB{excel}*AC{excel}/AD{excel},\"\")",
            "Primary Toxicity Endpoint": "",
            "Evidence": " | ".join(evidence),
            "Sources": " ".join(filter(None, [pc.get("pubchem_url"),
                                              f"{CHEMBL.replace('/api/data','')}/compound_report_card/{mol_id}/" if mol_id else "",
                                              pk_src.get("bioavailability", ""), pk_src.get("half_life", "")])),
        }
        table.append(rec)
        for col in COLUMNS:
            long_rows.append({"drug": name, "atc": atc, "route": route, "column": col,
                              "value": rec.get(col, ""), "source": rec["Sources"]})

    # 3D coordinates and contact residues for the binding animation
    m = Chem.AddHs(Chem.MolFromSmiles(pc["smiles"]))
    if AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE) == 0:
        AllChem.MMFFOptimizeMolecule(m)
        conf = m.GetConformer()
        binding[name] = {
            "smiles": pc["smiles"], "atc": atc, "cid": pc["cid"],
            "target": f"{tname} ({torg})" if tname else "",
            "target_chembl_id": target_id,
            "action": mechs[0]["action"] if mechs else "",
            "mechanism": mechs[0]["mechanism"] if mechs else f"target from {target_from}",
            "pchembl": act, "qed": d.get("QED"), "clogp": d.get("clogp"), "mw": d.get("MW"),
            "atoms": [{"el": a.GetSymbol(),
                       "p": [round(conf.GetAtomPosition(a.GetIdx()).x, 3),
                             round(conf.GetAtomPosition(a.GetIdx()).y, 3),
                             round(conf.GetAtomPosition(a.GetIdx()).z, 3)]}
                      for a in m.GetAtoms()],
            "bonds": [[b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()] for b in m.GetBonds()],
        }
    print(f"  {name}: CID {pc['cid']}, QED {d.get('QED')}, activity {act}, off-target {tox}, routes {chosen}")

df = pd.DataFrame(table)[COLUMNS]
df.to_csv(f"{OUT}/features.csv", index=False)
pd.DataFrame(long_rows).to_csv(f"{OUT}/features_long.csv", index=False)
pd.DataFrame(PROV).to_csv(f"{OUT}/provenance.csv", index=False)
json.dump(binding, open(f"{OUT}/binding_data.json", "w"), indent=1)

# ---------------------------------------------------------------- validation you can see
mols, legends = [], []
for name, b in binding.items():
    mols.append(Chem.MolFromSmiles(b["smiles"]))
    legends.append(f"{name} {b['atc']}  QED {b['qed']}  pChEMBL {b['pchembl']}")
if mols:
    Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(330, 260), legends=legends) \
        .save(f"{OUT}/structures.png")

try:
    import mols2grid
    grid = pd.DataFrame([{"SMILES": b["smiles"], "Name": n, "ATC": b["atc"], "QED": b["qed"],
                          "MW": b["mw"], "cLogP": b["clogp"], "pChEMBL": b["pchembl"],
                          "Target": b["target"]} for n, b in binding.items()])
    html = mols2grid.MolGrid(grid, smiles_col="SMILES", size=(230, 180)).display(
        subset=["img", "Name", "QED"],
        tooltip=["ATC", "MW", "cLogP", "pChEMBL", "Target"]).data
    open(f"{OUT}/grid.html", "w").write(html)
except Exception as e:                                  # mols2grid is a convenience, not a dependency
    print(f"  mols2grid grid skipped: {type(e).__name__}: {e}")

# QED and each physicochemical property, group values marked on the whole-dataset spread
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROPS = [("Quantitative Estimate of Drug Likeness", "QED"), ("Molecular Weight", "MW (Da)"),
         ("Lipophilicity(log D)", "log P"), ("Topological Surface Area", "TPSA (A^2)"),
         ("Hydrogen Bond Donors", "H-bond donors"), ("Hydrogen Bond Acceptors", "H-bond acceptors"),
         ("Rotatable Bonds", "Rotatable bonds"), ("Molecular Refractivity Or Taft Steric Factor", "Molar refractivity"),
         ("Aromatic Rings", "Aromatic rings")]
one = df.drop_duplicates("Generic Name")
fig, axes = plt.subplots(3, 3, figsize=(13, 9))
for ax, (col, label) in zip(axes.ravel(), PROPS):
    vals = pd.to_numeric(one[col], errors="coerce").dropna()
    ax.hist(vals, bins=max(4, min(10, len(vals))), color="#9ecae1", edgecolor="#30506b")
    for v, n in zip(vals, one.loc[vals.index, "Generic Name"]):
        ax.axvline(v, color="#c1443a", lw=1)
        ax.annotate(n[:9], (v, ax.get_ylim()[1]), rotation=90, va="top", ha="right", fontsize=6, color="#c1443a")
    ax.set_title(label, fontsize=9); ax.tick_params(labelsize=7)
fig.suptitle(f"{GROUP}: every physicochemical property on its own, {len(one)} drugs", fontsize=11)
fig.tight_layout()
fig.savefig(f"{OUT}/qed_distribution.png", dpi=150)

print(f"\n{len(df)} rows ({len(one)} drugs x routes) -> {OUT}/features.csv")
print(f"{len(PROV)} API calls logged -> {OUT}/provenance.csv")
