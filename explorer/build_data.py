#!/usr/bin/env python3
"""Build the data files for the Hansch Space 3D explorer (explorer/index.html).

Reads the curated table written by curation/colab_curate_all.py and adds PubChem data for every drug.

Colab:
    !pip -q install rdkit
    !wget -q -O build_data.py https://raw.githubusercontent.com/DhruvMenon2003/drug-toxicity-qsar-prediction/main/explorer/build_data.py
    %run build_data.py --curated curated/drugbank_curated_dedup.csv --out explorer_data

Local:
    python explorer/build_data.py --curated path/to/drugbank_curated_dedup.csv --out explorer/data

Sources, in this order (APIs only, no scraping of drug pages):
    PubChem PUG REST   parent CID by InChIKey, descriptors (XLogP3, TPSA, MW, H-bond counts, ...), 3D conformers
    WHO ATC index      ATC level 1-3 names (the curated table's level-3 names were wrong before this fix)
    RDKit              Crippen logP and molar refractivity, and a 3D conformer when PubChem has none

Outputs (in --out):
    drugs.json       one record per drug: ATC codes and names, routes, descriptors, PCA scores, ChEMBL activity
    molecules.json   one 3D conformer per drug: atomic numbers, coordinates in angstrom, bonds
    hansch.json      Hansch-type regressions of activity on logP, MR and TPSA per ChEMBL target, and for all drugs
    provenance.csv   every HTTP request with URL, status, size and UTC time (cached answers are marked)
    validation.md    cross-checks to read before trusting the page
"""
import argparse, csv, datetime as dt, hashlib, io, json, math, os, re, threading, time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
WHO = "https://atcddd.fhi.no/atc_ddd_index/?code={}&showdescription=no"
PROPS = ("XLogP,TPSA,MolecularWeight,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,"
         "HeavyAtomCount,Complexity,Charge,MolecularFormula,Title,InChIKey")
CID, NAME = "Compound Identifier", "Generic Name"
MIN_GROUP = 5            # smallest ChEMBL-target group that gets a regression
SEED = 20261009          # fixed seed: y-randomisation and RDKit embedding give the same answer on every run


# ---------------------------------------------------------------- HTTP: cached, rate-limited, logged
class Http:
    """GET with a disk cache, <= 4.5 requests/s, retries on busy answers, and a provenance log.
    PubChem reports this client's request budget in X-Throttling-Control; when it turns Red or Black, every thread pauses."""

    def __init__(self, cache_dir):
        self.cache_dir, self.log, self.lock, self.last, self.pause_until, self.down = cache_dir, [], threading.Lock(), 0.0, 0.0, {}
        os.makedirs(cache_dir, exist_ok=True)
        self.s = requests.Session()
        self.s.headers["User-Agent"] = "drug-toxicity-qsar-prediction explorer (research; github.com/DhruvMenon2003)"

    def _wait(self):
        with self.lock:
            gap = max(0.22 - (time.monotonic() - self.last), self.pause_until - time.monotonic())
            if gap > 0:
                time.sleep(gap)
            self.last = time.monotonic()

    def _throttle(self, r):
        """Back off from PubChem's own load report and Retry-After before the server starts refusing."""
        mine = re.findall(r"Request (?:Count|Time) status: (\w+)", r.headers.get("X-Throttling-Control", ""))
        wait = float(r.headers.get("Retry-After", 0) or 0) if r.status_code in (429, 503) else 0
        if "Black" in mine:                      # our own request budget; "Service status" is the server's overall load
            wait = max(wait, 30)
        elif "Red" in mine:
            wait = max(wait, 5)
        if wait:
            with self.lock:
                self.pause_until = max(self.pause_until, time.monotonic() + wait)

    def get(self, url):
        """Return (status, text). 404 is a normal answer ("not found") and is cached too."""
        path = os.path.join(self.cache_dir, hashlib.sha1(url.encode()).hexdigest() + ".json")
        if os.path.exists(path):
            c = json.load(open(path, encoding="utf-8"))
            self.log.append({"url": url, "status": c["status"], "bytes": len(c["text"]), "utc": c["utc"], "cached": True})
            return c["status"], c["text"]
        host = url.split("/")[2]
        if self.down.get(host, 0) >= 3:      # three requests in a row ran out of retries: stop hammering a server that is down
            self.log.append({"url": url, "status": None, "bytes": 0, "utc": "", "cached": False})
            return None, "skipped: host unavailable during this run"
        status, text = None, ""
        for attempt in range(6):
            self._wait()
            try:
                r = self.s.get(url, timeout=120)
                status, text = r.status_code, r.text
                self._throttle(r)
            except requests.RequestException as e:
                status, text = None, str(e)
            busy = status in (429, 500, 502, 503, 504) or "ServerBusy" in text[:300]
            if not busy and status is not None:
                break
            time.sleep(min(60, 2 ** (attempt + 1)))
        utc = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
        with self.lock:
            self.down[host] = 0 if status in (200, 404) else self.down.get(host, 0) + 1
        self.log.append({"url": url, "status": status, "bytes": len(text), "utc": utc, "cached": False})
        if status in (200, 404):
            json.dump({"status": status, "text": text, "utc": utc}, open(path, "w", encoding="utf-8"))
        return status, text


def chunks(xs, n):
    return [xs[i:i + n] for i in range(0, len(xs), n)]


# ---------------------------------------------------------------- ChEMBL target text in the curated table
def first_target(text):
    """'CHEMBL245 (Muscarinic M3 antagonist); CHEMBL216 (...)' -> ('CHEMBL245', 'Muscarinic M3 antagonist').
    Mechanism texts can contain their own brackets, e.g. 'Serotonin 2a (5-HT2a) receptor antagonist'."""
    if not isinstance(text, str):
        return None, None
    m = re.search(r"(CHEMBL\d+) \(", text)
    if not m:
        return None, None
    depth, i = 1, m.end()
    while i < len(text) and depth:
        depth += {"(": 1, ")": -1}.get(text[i], 0)
        i += 1
    return m.group(1), text[m.end():i - 1 if depth == 0 else i].strip()


# ---------------------------------------------------------------- RDKit
def rdkit_descriptors(mol):
    return {"clogp": Crippen.MolLogP(mol), "mr": Crippen.MolMR(mol), "tpsa_rd": rdMolDescriptors.CalcTPSA(mol),
            "mw_rd": Descriptors.MolWt(mol), "hbd_rd": Lipinski.NumHDonors(mol), "hba_rd": Lipinski.NumHAcceptors(mol),
            "rotb_rd": rdMolDescriptors.CalcNumRotatableBonds(mol), "arom": rdMolDescriptors.CalcNumAromaticRings(mol),
            "fsp3": rdMolDescriptors.CalcFractionCSP3(mol), "heavy_rd": mol.GetNumHeavyAtoms()}


def embed_3d(smiles):
    """ETKDGv3 conformer + MMFF (UFF if MMFF has no parameters). None if embedding fails or the molecule is huge."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumHeavyAtoms() > 150:
        return None
    mol = Chem.AddHs(mol)
    p = AllChem.ETKDGv3()
    p.randomSeed = SEED
    if AllChem.EmbedMolecule(mol, p) != 0:
        p.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, p) != 0:
            return None
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
    except Exception:
        pass
    return mol


def compact(mol, source):
    """Molecule -> {'z': atomic numbers, 'xyz': centred coordinates (2 dp), 'b': [i, j, order, ...], 'src'}.
    Bond order 4 means aromatic (RDKit conformers); PubChem SDF bonds are Kekule 1/2/3."""
    pos = mol.GetConformer().GetPositions()
    pos = pos - pos.mean(axis=0)
    bonds = []
    for b in mol.GetBonds():
        o = b.GetBondTypeAsDouble()
        bonds += [b.GetBeginAtomIdx(), b.GetEndAtomIdx(), 4 if o == 1.5 else int(o)]
    return {"z": [a.GetAtomicNum() for a in mol.GetAtoms()], "xyz": [round(float(v), 2) for v in pos.ravel()],
            "b": bonds, "src": source}


# ---------------------------------------------------------------- regression helpers
def ols(X, y):
    """Least squares with intercept. Returns coef (intercept first), standard errors, fitted values, s."""
    A = np.column_stack([np.ones(len(y)), X])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    fit = A @ coef
    dof = len(y) - A.shape[1]
    s2 = float(((y - fit) ** 2).sum() / dof) if dof > 0 else float("nan")
    try:
        se = np.sqrt(np.diag(s2 * np.linalg.inv(A.T @ A)))
    except np.linalg.LinAlgError:
        se = np.full(A.shape[1], np.nan)
    return coef, se, fit, math.sqrt(s2) if dof > 0 else float("nan")


def r2(y, fit):
    ss = ((y - y.mean()) ** 2).sum()
    return float(1 - ((y - fit) ** 2).sum() / ss) if ss > 0 else float("nan")


def loo(X, y):
    pred = np.empty(len(y))
    for i in range(len(y)):
        keep = np.arange(len(y)) != i
        coef, *_ = ols(X[keep], y[keep])
        pred[i] = coef[0] + X[i] @ coef[1:]
    return pred


MODELS = [("logP",), ("MR",), ("TPSA",), ("logP", "logP²"), ("logP", "MR"), ("logP", "TPSA"),
          ("logP", "logP²", "MR"), ("logP", "MR", "TPSA"), ("logP", "logP²", "MR", "TPSA")]


def design(rows, terms):
    col = {"logP": [r["logp"] for r in rows], "logP²": [r["logp"] ** 2 for r in rows],
           "MR": [r["mr"] for r in rows], "TPSA": [r["tpsa"] for r in rows]}
    return np.column_stack([col[t] for t in terms]).astype(float)


def hansch(rows, rng):
    """Pick the model with the best leave-one-out q2 among those with n >= 5 drugs per descriptor (Topliss & Costello).
    Ties go to the smaller model. Then y-randomisation: refit on 100 shuffled activity vectors."""
    y = np.array([r["p"] for r in rows], float)
    tried = []
    for terms in MODELS:
        if len(y) < 5 * len(terms):
            continue
        X = design(rows, terms)
        if np.linalg.matrix_rank(np.column_stack([np.ones(len(y)), X])) < X.shape[1] + 1:
            continue
        coef, se, fit, s = ols(X, y)
        q = r2(y, loo(X, y))
        tried.append({"terms": terms, "coef": coef, "se": se, "fit": fit, "s": s, "r2": r2(y, fit), "q2": q})
    if not tried:
        return None
    best = max(tried, key=lambda m: (round(m["q2"], 3), -len(m["terms"])))
    X = design(rows, best["terms"])
    scr = [r2(ys, ols(X, ys)[2]) for ys in (rng.permutation(y) for _ in range(100))]
    p_, n = len(best["terms"]), len(y)
    F = (best["r2"] / p_) / ((1 - best["r2"]) / (n - p_ - 1)) if best["r2"] < 1 and n > p_ + 1 else float("nan")
    pred = loo(X, y)
    rnd = lambda v, k=3: None if v is None or not np.isfinite(v) else round(float(v), k)
    return {"n": n, "terms": list(best["terms"]),
            "coef": {t: rnd(c, 4) for t, c in zip(["const", *best["terms"]], best["coef"])},
            "se": {t: rnd(c, 4) for t, c in zip(["const", *best["terms"]], best["se"])},
            "r2": rnd(best["r2"]), "q2": rnd(best["q2"]), "s": rnd(best["s"]), "F": rnd(F, 2),
            "yrand_r2_mean": rnd(np.mean(scr)), "yrand_r2_max": rnd(np.max(scr)),
            "yrand_beats": int(sum(v >= best["r2"] for v in scr)),
            "candidates": [{"terms": list(m["terms"]), "r2": rnd(m["r2"]), "q2": rnd(m["q2"])} for m in tried],
            "points": [{"cid": r["cid"], "obs": r["p"], "fit": rnd(f), "loo": rnd(l)} for r, f, l in zip(rows, best["fit"], pred)]}


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--curated", default="curated/drugbank_curated_dedup.csv")
    ap.add_argument("--out", default="explorer/data")
    ap.add_argument("--cache", default=".pubchem_cache")
    ap.add_argument("--fix-curated-l3", action="store_true", help="write the corrected ATC_L3_name back into --curated")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    http = Http(a.cache)
    started = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    df = pd.read_csv(a.curated, dtype={CID: str, "ATC": str})
    print(f"curated rows {len(df)}, drugs {df[CID].nunique()}, ATC codes {df.ATC.nunique()}")

    # ---- ATC names (WHO). The page lists the level-1/2/3 names above the level-4 table.
    names = {}
    l3s = sorted(df.ATC.str[:4].unique())
    for l3 in l3s:
        st, text = http.get(WHO.format(l3))
        for code, label in re.findall(r'<a href="\./\?code=([A-Z0-9]+)&showdescription=(?:no|yes)">([^<]+)</a>', text or ""):
            if "text from Guidelines" not in label:
                names.setdefault(code, label.strip())
    for n, k in [(1, 1), (2, 3)]:              # keep the curated names if the WHO page did not answer
        for code, nm in df[[f"ATC_L{n}", f"ATC_L{n}_name"]].drop_duplicates().itertuples(index=False):
            names.setdefault(code, nm)
    missing_l3 = [c for c in l3s if c not in names]
    print(f"ATC level-3 names: {len(l3s) - len(missing_l3)}/{len(l3s)}")
    if a.fix_curated_l3:
        df["ATC_L3_name"] = df.ATC.str[:4].map(names)
        df.to_csv(a.curated, index=False)
        print(f"wrote corrected ATC_L3_name into {a.curated}")

    # ---- one record per drug (dataset CID)
    drugs = {}
    for cid, g in df.groupby(CID, sort=False):
        r0 = g.iloc[0]
        atc = []
        for code, gg in g.groupby("ATC", sort=False):
            atc.append({"code": code, "routes": sorted(gg["Mode of Administration"].dropna().unique().tolist())})
        tid, tname = first_target(r0.get("Activity_target"))
        act = None
        if pd.notna(r0.get("Activity")) and str(r0.get("Activity")) not in ("", "NaN"):
            act = {"p": float(r0.Activity), "target": tid, "tname": tname,
                   "types": r0.get("Activity_types") if pd.notna(r0.get("Activity_types")) else "",
                   "ndoc": int(r0.Activity_n_documents) if pd.notna(r0.get("Activity_n_documents")) else None,
                   "chembl": r0.get("Activity_chembl_id") if pd.notna(r0.get("Activity_chembl_id")) else None}
        flags = [f for f in [r0.get("Name_Edit_Flag")] if isinstance(f, str) and f]
        if bool(r0.get("SMILES_shared_with_other_drug")):
            flags.append("structure shared with another drug in the dataset")
        drugs[cid] = {"cid": int(cid), "name": r0[NAME], "smiles": r0["2D Chemical Structure"], "ik": r0.InChIKey,
                      "atc": atc, "routes": sorted({x for e in atc for x in e["routes"]}), "act": act,
                      "chembl": r0.get("Activity_chembl_id") if pd.notna(r0.get("Activity_chembl_id")) else None,
                      "flags": flags}

    # ---- PubChem descriptors, 100 CIDs per request (dataset CIDs first)
    props = {}
    def fetch_props(cids):
        for ch in chunks(sorted(cids), 100):
            st, text = http.get(f"{PUG}/compound/cid/{','.join(map(str, ch))}/property/{PROPS}/JSON")
            if st == 200:
                for p in json.loads(text)["PropertyTable"]["Properties"]:
                    props[p["CID"]] = p
    fetch_props({d["cid"] for d in drugs.values()})

    # ---- parent CID: the curated SMILES is the neutral parent (salts stripped). Where the dataset CID already has that
    # InChIKey it is the parent; otherwise look the InChIKey up (one request each, so only for the salt/mixture records).
    keys = sorted({d["ik"] for d in drugs.values() if isinstance(d["ik"], str) and props.get(d["cid"], {}).get("InChIKey") != d["ik"]})
    def lookup(k):
        st, text = http.get(f"{PUG}/compound/inchikey/{k}/cids/JSON")
        try:
            return k, json.loads(text)["IdentifierList"]["CID"][0] if st == 200 else None
        except (ValueError, KeyError, IndexError):
            return k, None
    with ThreadPoolExecutor(2) as ex:
        parent = dict(ex.map(lookup, keys))
    for d in drugs.values():
        if props.get(d["cid"], {}).get("InChIKey") == d["ik"]:
            d["pcid"], d["pcid_src"] = d["cid"], "dataset CID (same InChIKey as the curated parent)"
        elif parent.get(d["ik"]):
            d["pcid"], d["pcid_src"] = parent[d["ik"]], "PubChem record for the curated parent's InChIKey"
        else:
            d["pcid"], d["pcid_src"] = d["cid"], "dataset CID (curated parent's InChIKey not in PubChem)"
    pcids = sorted({d["pcid"] for d in drugs.values()})
    print(f"parent CIDs {len(pcids)} ({sum(d['pcid'] != d['cid'] for d in drugs.values())} differ from the dataset CID)")

    fetch_props(set(pcids) - set(props))
    print(f"PubChem property records {sum(c in props for c in pcids)}/{len(pcids)}")

    # ---- PubChem 3D conformers (first conformer of each record), 50 CIDs per request; RDKit fills the gaps
    mols = {}
    for ch in chunks(pcids, 50):
        st, text = http.get(f"{PUG}/compound/cid/{','.join(map(str, ch))}/SDF?record_type=3d")
        if st != 200:
            continue
        for m in Chem.ForwardSDMolSupplier(io.BytesIO(text.encode()), removeHs=False, sanitize=False):
            if m is not None and m.HasProp("PUBCHEM_COMPOUND_CID"):
                mols[int(m.GetProp("PUBCHEM_COMPOUND_CID"))] = compact(m, "PubChem 3D conformer")
    print(f"PubChem 3D conformers {len(mols)}/{len(pcids)}")
    molecules, n_rd = {}, 0
    for d in drugs.values():
        m = mols.get(d["pcid"])
        if m is None:
            em = embed_3d(d["smiles"])
            if em is not None:
                m = compact(em, "RDKit ETKDGv3 + force field (PubChem has no 3D record)")
                n_rd += 1
        if m is not None:
            molecules[str(d["cid"])] = m
        d["conf"] = m["src"] if m else None
    print(f"3D conformers: PubChem {sum(1 for m in molecules.values() if m['src'].startswith('PubChem'))}, RDKit {n_rd}, none {len(drugs) - len(molecules)}")

    # ---- descriptors (PubChem first, RDKit for Crippen logP/MR and as a cross-check)
    for d in drugs.values():
        mol = Chem.MolFromSmiles(d["smiles"]) if isinstance(d["smiles"], str) else None
        rd = rdkit_descriptors(mol) if mol else {}
        p = props.get(d["pcid"], {})
        xl = p.get("XLogP")
        desc = {"logp": xl if xl is not None else rd.get("clogp"), "logp_src": "XLogP3 (PubChem)" if xl is not None else "Crippen (RDKit)",
                "xlogp": xl, "clogp": rd.get("clogp"), "mr": rd.get("mr"),
                "tpsa": p.get("TPSA", rd.get("tpsa_rd")), "mw": float(p["MolecularWeight"]) if "MolecularWeight" in p else rd.get("mw_rd"),
                "hbd": p.get("HBondDonorCount", rd.get("hbd_rd")), "hba": p.get("HBondAcceptorCount", rd.get("hba_rd")),
                "rotb": p.get("RotatableBondCount", rd.get("rotb_rd")), "heavy": p.get("HeavyAtomCount", rd.get("heavy_rd")),
                "cx": p.get("Complexity"), "charge": p.get("Charge"), "arom": rd.get("arom"), "fsp3": rd.get("fsp3"),
                "formula": p.get("MolecularFormula"), "_mw_rd": rd.get("mw_rd"), "_tpsa_rd": rd.get("tpsa_rd"), "_heavy_rd": rd.get("heavy_rd")}
        if desc["logp"] is not None:
            desc["ro5"] = int((desc["mw"] or 0) > 500) + int(desc["logp"] > 5) + int((desc["hbd"] or 0) > 5) + int((desc["hba"] or 0) > 10)
        d["title"] = p.get("Title")
        d["desc_src"] = "PubChem" if p else "RDKit (no PubChem record in this build)"
        d["d"] = desc

    # ---- PCA on standardised descriptors (drugs with all of them)
    feats = ["logp", "mr", "tpsa", "mw", "hbd", "hba", "rotb", "arom", "fsp3"]
    ok = [d for d in drugs.values() if all(d["d"].get(f) is not None for f in feats)]
    X = np.array([[d["d"][f] for f in feats] for d in ok], float)
    mu, sd = X.mean(0), X.std(0)
    U, S, Vt = np.linalg.svd((X - mu) / sd, full_matrices=False)
    scores = U[:, :3] * S[:3]
    for d, sc in zip(ok, scores):
        d["pc"] = [round(float(v), 3) for v in sc]
    explained = (S ** 2 / (S ** 2).sum())[:3]

    # ---- Hansch regressions: one structure per group (shared SMILES counted once)
    rng = np.random.default_rng(SEED)
    act_rows, seen = [], set()
    for d in drugs.values():
        if d["act"] and d["d"].get("logp") is not None and d["d"].get("mr") is not None and d["ik"] not in seen:
            seen.add(d["ik"])
            act_rows.append({"cid": d["cid"], "p": d["act"]["p"], "target": d["act"]["target"], "tname": d["act"]["tname"],
                             "logp": d["d"]["logp"], "mr": d["d"]["mr"], "tpsa": d["d"]["tpsa"]})
    groups = []
    for tid, g in pd.DataFrame(act_rows).groupby("target"):
        if len(g) >= MIN_GROUP:
            res = hansch(g.to_dict("records"), rng)
            if res:
                groups.append({"id": tid, "label": g.tname.iloc[0], **res})
    groups.sort(key=lambda x: -x["n"])
    overall = hansch(act_rows, rng)
    overall.update(id="ALL", label="All drugs with an activity value (unrelated targets)")
    print(f"Hansch fits: {len(groups)} target groups + all drugs (n={overall['n']})")

    # ---- write
    clean = lambda v: None if isinstance(v, float) and not math.isfinite(v) else (round(v, 3) if isinstance(v, float) else v)
    out_drugs = []
    for d in drugs.values():
        dd = {k: clean(v) for k, v in d["d"].items() if not k.startswith("_")}
        out_drugs.append({**{k: v for k, v in d.items() if k not in ("d",)}, "d": dd})
    meta = {"generated_utc": started, "n_drugs": len(out_drugs), "atc_names": names,
            "pca": {"features": feats, "explained": [round(float(v), 4) for v in explained],
                    "loadings": {f: [round(float(v), 3) for v in Vt[:3, i]] for i, f in enumerate(feats)}},
            "sources": {"pubchem": PUG, "who_atc": "https://atcddd.fhi.no/atc_ddd_index/", "curated": os.path.basename(a.curated)}}
    json.dump({"meta": meta, "drugs": out_drugs}, open(os.path.join(a.out, "drugs.json"), "w"), separators=(",", ":"))
    json.dump(molecules, open(os.path.join(a.out, "molecules.json"), "w"), separators=(",", ":"))
    json.dump({"min_group": MIN_GROUP, "seed": SEED, "groups": [overall, *groups]}, open(os.path.join(a.out, "hansch.json"), "w"), separators=(",", ":"))
    with open(os.path.join(a.out, "provenance.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, ["url", "status", "bytes", "utc", "cached"])
        w.writeheader()
        w.writerows(http.log)
    write_validation(a, drugs, molecules, groups, overall, names, l3s, missing_l3, http.log, explained, feats, started)
    print("done:", ", ".join(sorted(os.listdir(a.out))))


def write_validation(a, drugs, molecules, groups, overall, names, l3s, missing_l3, log, explained, feats, started):
    D = list(drugs.values())
    pair = lambda k1, k2: np.array([(d["d"][k1], d["d"][k2]) for d in D if d["d"].get(k1) is not None and d["d"].get(k2) is not None], float)
    xl = pair("xlogp", "clogp")
    mw = pair("mw", "_mw_rd")
    tp = pair("tpsa", "_tpsa_rd")
    worst_logp = sorted((d for d in D if d["d"].get("xlogp") is not None and d["d"].get("clogp") is not None),
                        key=lambda d: -abs(d["d"]["xlogp"] - d["d"]["clogp"]))[:8]
    mw_bad = [d for d in D if d["d"].get("mw") is not None and d["d"].get("_mw_rd") is not None and abs(d["d"]["mw"] - d["d"]["_mw_rd"]) > 0.5]
    heavy_bad = []
    for d in D:
        m = molecules.get(str(d["cid"]))
        if m and d["d"].get("_heavy_rd") is not None and sum(z > 1 for z in m["z"]) != d["d"]["_heavy_rd"]:
            heavy_bad.append(d)
    hosts = pd.DataFrame(log).assign(host=lambda t: t.url.str.extract(r"https?://([^/]+)")[0]).groupby(["host", "status", "cached"]).size()
    L = [f"# Hansch Space data checks", "", f"Built {started} from `{os.path.basename(a.curated)}`.", "",
         "## Counts", "",
         f"- Drugs (dataset CIDs): {len(D)}",
         f"- Parent CID differs from the dataset CID (salt or mixture record in the source): {sum(d['pcid'] != d['cid'] for d in D)}",
         f"- PubChem descriptors found: {sum(d['desc_src'] == 'PubChem' for d in D)} (the rest use RDKit values from the curated SMILES)",
         f"- PubChem requests that failed after retries or were skipped while PubChem was down: "
         f"{sum(1 for r in log if 'pubchem' in r['url'] and r['status'] not in (200, 404))} (rerun the build to fill them from PubChem)",
         f"- 3D conformers from PubChem: {sum(1 for m in molecules.values() if m['src'].startswith('PubChem'))}; "
         f"from RDKit: {sum(1 for m in molecules.values() if m['src'].startswith('RDKit'))}; none: {len(D) - len(molecules)}",
         f"- Drugs with a ChEMBL activity value: {sum(d['act'] is not None for d in D)}",
         f"- ATC level-3 names found: {len(l3s) - len(missing_l3)}/{len(l3s)}" + (f" (missing: {', '.join(missing_l3)})" if missing_l3 else ""),
         "", "## Cross-checks (PubChem against RDKit on the curated SMILES)", "",
         f"- Molecular weight: {len(mw)} pairs, max difference {np.abs(mw[:, 0] - mw[:, 1]).max():.2f} Da; "
         f"{len(mw_bad)} drugs differ by more than 0.5 Da (would mean the parent CID is the wrong record).",
         f"- TPSA: {len(tp)} pairs, Pearson r {np.corrcoef(tp.T)[0, 1]:.3f}, mean absolute difference {np.abs(tp[:, 0] - tp[:, 1]).mean():.2f} A^2 "
         "(PubChem counts S and P polar surface; RDKit's default does not).",
         f"- logP: XLogP3 against Crippen, {len(xl)} pairs, Pearson r {np.corrcoef(xl.T)[0, 1]:.3f}, "
         f"mean absolute difference {np.abs(xl[:, 0] - xl[:, 1]).mean():.2f}. Two different logP models; large gaps are listed below.",
         f"- 3D conformer heavy-atom count equals the curated SMILES: {len(molecules) - len(heavy_bad)}/{len(molecules)}"
         + (f"; mismatches: {', '.join(d['name'] for d in heavy_bad[:15])}" if heavy_bad else ""), ""]
    if mw_bad:
        L += ["Molecular weight mismatches: " + ", ".join(f"{d['name']} ({d['d']['mw']:.1f} vs {d['d']['_mw_rd']:.1f})" for d in mw_bad[:20]), ""]
    L += ["Largest XLogP3 vs Crippen gaps:", "", "| Drug | XLogP3 | Crippen |", "|---|---|---|"]
    L += [f"| {d['name']} | {d['d']['xlogp']:.2f} | {d['d']['clogp']:.2f} |" for d in worst_logp]
    L += ["", f"PCA on standardised {', '.join(feats)}: PC1-3 explain " + ", ".join(f"{v:.1%}" for v in explained) + ".", "",
          "## Hansch-type fits", "",
          "Activity is the median ChEMBL pChEMBL at the drug's mechanism target(s). Drugs are grouped by their first listed target. "
          "These groups are not congeneric series, so the fits are exploratory. A fit is credible only when q2 is well above 0 and "
          "the shuffled-activity r2 (y-randomisation, 100 runs) stays far below the real r2.", "",
          "| Group | n | Model | r2 | q2 (LOO) | s | shuffled r2 mean / max |", "|---|---|---|---|---|---|---|"]
    for g in [overall, *groups]:
        L.append(f"| {g['label']} ({g['id']}) | {g['n']} | {' + '.join(g['terms'])} | {g['r2']} | {g['q2']} | {g['s']} | {g['yrand_r2_mean']} / {g['yrand_r2_max']} |")
    L += ["", "## Requests", "", "| Host | Status | From cache | Count |", "|---|---|---|---|"]
    L += [f"| {h} | {s} | {c} | {n} |" for (h, s, c), n in hosts.items()]
    open(os.path.join(a.out, "validation.md"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
