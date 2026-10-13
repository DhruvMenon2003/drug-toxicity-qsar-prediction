# Build the data the ligand-binding animation plays: one receptor per target group, with its real
# topology, plus the drugs of that group with real 3D ligand coordinates.
#
#   python3 explorer/binding/build_binding.py --features h2_features/binding_data.json
#
# Nothing here is drawn from memory:
#   ligand 3D coordinates   RDKit ETKDG + MMFF, from the PubChem SMILES in the feature table
#   receptor topology       UniProt (transmembrane helix ranges, pocket residues, disulfide, sequence)
#   affinity                ChEMBL pChEMBL carried over from the feature table
#   structures to cite      RCSB PDB entries for the same UniProt accession
#
# The helix *positions* on screen are a schematic class-A arrangement, not coordinates from a structure;
# the page says so, and every residue label it prints is the one UniProt gives for that receptor.

import argparse, json, os, re, sys, time, datetime as dt
import requests

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
S = requests.Session()
S.headers["User-Agent"] = "drug-curation (research; contact via GitHub)"
PROV = []

AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
       "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
       "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"}


def jget(url, **kw):
    for i in range(4):
        try:
            r = S.get(url, timeout=90, **kw)
            PROV.append({"url": r.url, "status": r.status_code,
                         "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")})
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 500, 502, 503, 504) and i < 3:
                time.sleep(2 ** (i + 1)); continue
            return {}
        except requests.RequestException:
            time.sleep(2 ** i)
    return {}


def uniprot_accession(target_chembl_id):
    t = jget(f"{CHEMBL}/target/{target_chembl_id}", params={"format": "json"})
    for c in t.get("target_components", []):
        if c.get("accession"):
            return c["accession"], t.get("pref_name"), t.get("target_type")
    return None, t.get("pref_name"), t.get("target_type")


def receptor(accession):
    """transmembrane helices, pocket residues, disulfide bonds, straight from UniProt."""
    u = jget(f"https://rest.uniprot.org/uniprotkb/{accession}.json")
    if not u:
        return None
    seq = u.get("sequence", {}).get("value", "")
    helices, sites, disulfides, topo = [], [], [], []
    for f in u.get("features", []):
        loc = f.get("location", {})
        a, b = loc.get("start", {}).get("value"), loc.get("end", {}).get("value")
        if f["type"] == "Transmembrane":
            name = re.search(r"Name=(\w+)", f.get("description", "") or "")
            helices.append({"start": a, "end": b, "name": name.group(1) if name else str(len(helices) + 1)})
        elif f["type"] in ("Binding site", "Site"):
            res = seq[a - 1] if seq and a and a <= len(seq) else ""
            sites.append({"pos": a, "residue": AA3.get(res, res),
                          "label": f"{AA3.get(res, res)}{a}", "note": f.get("description", "")})
        elif f["type"] == "Disulfide bond":
            disulfides.append([a, b])
        elif f["type"] == "Topological domain":
            topo.append({"start": a, "end": b, "side": f.get("description", "")})
    pdbs = []
    for x in u.get("uniProtKBCrossReferences", []):
        if x.get("database") == "PDB":
            method = next((p["value"] for p in x.get("properties", []) if p["key"] == "Method"), "")
            res = next((p["value"] for p in x.get("properties", []) if p["key"] == "Resolution"), "")
            pdbs.append({"id": x["id"], "method": method, "resolution": res})
    return {"accession": accession,
            "name": u.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""),
            "gene": (u.get("genes") or [{}])[0].get("geneName", {}).get("value", ""),
            "length": u.get("sequence", {}).get("length"),
            "family": next((c["texts"][0]["value"] for c in u.get("comments", [])
                            if c.get("commentType") == "SIMILARITY" and c.get("texts")), ""),
            "helices": helices, "sites": sites, "disulfides": disulfides, "topology": topo,
            "pdb": pdbs[:8],
            "uniprot_url": f"https://www.uniprot.org/uniprotkb/{accession}"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", default="h2_features/binding_data.json",
                    help="binding_data.json written by curation/h2_features.py")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "data"))
    args = ap.parse_args()
    if not os.path.exists(args.features):
        sys.exit(f"{args.features} not found: run curation/h2_features.py first")
    ligands = json.load(open(args.features))
    os.makedirs(args.out, exist_ok=True)

    groups = {}
    for name, lig in ligands.items():
        tid = lig.get("target_chembl_id")
        if not tid:
            print(f"  {name}: no mechanism target, left out of the animation"); continue
        if tid not in groups:
            acc, pref, ttype = uniprot_accession(tid)
            rec = receptor(acc) if acc else None
            groups[tid] = {"target_chembl_id": tid, "target": pref, "target_type": ttype,
                           "chembl_url": f"https://www.ebi.ac.uk/chembl/target_report_card/{tid}/",
                           "receptor": rec, "ligands": []}
            print(f"  {pref}: UniProt {acc}, "
                  f"{len(rec['helices']) if rec else 0} TM helices, "
                  f"{len(rec['sites']) if rec else 0} annotated pocket residues")
        groups[tid]["ligands"].append({"name": name, **{k: v for k, v in lig.items()
                                                        if k != "target_chembl_id"}})

    for g in groups.values():
        g["ligands"].sort(key=lambda l: (l.get("pchembl") is None, -(l.get("pchembl") or 0)))
    out = {"built_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
           "groups": list(groups.values())}
    json.dump(out, open(os.path.join(args.out, "receptors.json"), "w"), indent=1)
    import csv
    with open(os.path.join(args.out, "provenance.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, ["url", "status", "utc"]); w.writeheader(); w.writerows(PROV)
    print(f"{len(groups)} target group(s), {sum(len(g['ligands']) for g in groups.values())} ligands "
          f"-> {args.out}/receptors.json")


if __name__ == "__main__":
    main()
