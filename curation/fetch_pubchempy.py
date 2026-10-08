"""SMILES + PubChem properties strictly via pubchempy (project rule: SMILES only from pubchempy),
plus RDKit-computed physicochemical descriptors from that SMILES.

Usage: python fetch_pubchempy.py            -> output/pubchem.jsonl
"""
import json
import os
import time

import pubchempy as pcp

from drugcuration.inputs import load_drugs

try:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
except ImportError:  # descriptors are optional
    Chem = None

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output", "pubchem.jsonl")


def smiles_of(c):
    # PubChem renamed IsomericSMILES -> SMILES in 2025; support both pubchempy generations.
    for attr in ("isomeric_smiles", "smiles", "canonical_smiles", "connectivity_smiles"):
        v = getattr(c, attr, None)
        if v:
            return v, attr
    return None, None


def rdkit_panel(smi):
    if not (Chem and smi):
        return {}
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return {"rdkit_error": "unparseable SMILES"}
    return {
        "MW": round(Descriptors.MolWt(m), 3), "ExactMW": round(Descriptors.ExactMolWt(m), 4),
        "nHA": Lipinski.NumHAcceptors(m), "nHD": Lipinski.NumHDonors(m),
        "nRot": Lipinski.NumRotatableBonds(m), "nRing": rdMolDescriptors.CalcNumRings(m),
        "nAromRing": rdMolDescriptors.CalcNumAromaticRings(m), "nHet": rdMolDescriptors.CalcNumHeteroatoms(m),
        "nHeavy": m.GetNumHeavyAtoms(), "fChar": Chem.GetFormalCharge(m),
        "TPSA": round(rdMolDescriptors.CalcTPSA(m), 2), "logP_Crippen": round(Crippen.MolLogP(m), 3),
        "MR": round(Crippen.MolMR(m), 3), "Fsp3": round(rdMolDescriptors.CalcFractionCSP3(m), 3),
        "nStereo": len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)),
    }


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as out:
        for d in load_drugs():
            rec = {"generic_name": d["generic_name"], "pubchem_cid": d["pubchem_cid"], "source": "pubchempy"}
            try:
                c = pcp.Compound.from_cid(int(d["pubchem_cid"]))
                smi, attr = smiles_of(c)
                rec.update({
                    "smiles": smi, "smiles_field": attr, "molecular_formula": c.molecular_formula,
                    "molecular_weight": c.molecular_weight, "monoisotopic_mass": c.monoisotopic_mass,
                    "xlogp": c.xlogp, "iupac_name": c.iupac_name, "tpsa": c.tpsa,
                    "inchikey": c.inchikey, "status": "ok",
                    "source_urls": [f"https://pubchem.ncbi.nlm.nih.gov/compound/{d['pubchem_cid']}"],
                })
                rec["rdkit"] = rdkit_panel(smi)
            except Exception as e:  # keep going; record the failure
                rec.update({"status": "error", "error": repr(e)[:300]})
            out.write(json.dumps(rec) + "\n")
            print(rec["generic_name"], rec["status"])
            time.sleep(0.25)  # PubChem limit: <= 5 requests/second


if __name__ == "__main__":
    main()
