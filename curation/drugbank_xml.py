"""Alternative to scraping: read DrugBank's official full-database XML release.

DrugBank offers the complete database as XML (CC BY-NC 4.0) to anyone with a free academic account:
https://go.drugbank.com/releases/latest  ->  "Full database" (drugbank_all_full_database.xml.zip).
This gives the same fields the drugbank spider scrapes (ClassyFire taxonomy, all substituents,
experimental properties, toxicity, PK/PD text, dosage routes, targets) without loading any web page,
and it is the route DrugBank's terms intend for research datasets.

Usage:  python drugbank_xml.py path\\to\\drugbank_all_full_database.xml.zip   -> output/drugbank_xml.jsonl
merge.py uses this file automatically for any drug the live scrape did not return.
"""
import json
import os
import sys
import zipfile
import xml.etree.ElementTree as ET

from drugcuration.inputs import load_drugs

NS = "{http://www.drugbank.ca}"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output", "drugbank_xml.jsonl")
TEXT = {"toxicity": "toxicity", "pharmacodynamics": "pharmacodynamics", "mechanism-of-action": "mechanism_of_action",
        "absorption": "absorption", "volume-of-distribution": "volume_of_distribution",
        "protein-binding": "protein_binding", "metabolism": "metabolism",
        "route-of-elimination": "route_of_elimination", "half-life": "half_life", "clearance": "clearance"}


def t(el, tag):
    x = el.find(NS + tag)
    return (x.text or "").strip() if x is not None and x.text else None


def open_xml(path):
    if path.lower().endswith(".zip"):
        z = zipfile.ZipFile(path)
        name = next(n for n in z.namelist() if n.lower().endswith(".xml"))
        return z.open(name)
    return open(path, "rb")


def drug_record(d):
    ids = [x.text for x in d.findall(NS + "drugbank-id")]
    primary = next((x.text for x in d.findall(NS + "drugbank-id") if x.get("primary") == "true"), ids[0])
    rec = {"drugbank_id": primary, "generic_name": t(d, "name"), "status": "ok", "capture": "drugbank_xml.jsonl",
           "source_urls": [f"https://go.drugbank.com/drugs/{primary}"], "groups":
           [g.text for g in d.findall(f"{NS}groups/{NS}group")]}
    c = d.find(NS + "classification")
    if c is not None:
        rec.update({"kingdom": t(c, "kingdom"), "super_class": t(c, "superclass"), "class": t(c, "class"),
                    "sub_class": t(c, "subclass"), "direct_parent": t(c, "direct-parent"),
                    "molecular_framework": t(c, "molecular-framework"),
                    "alternative_parents": [x.text for x in c.findall(NS + "alternative-parent")],
                    "substituents": [x.text for x in c.findall(NS + "substituent")]})
    for tag, key in TEXT.items():
        rec[key] = t(d, tag)
    rec["experimental_properties"] = [
        {"Property": t(p, "kind"), "Value": t(p, "value"), "Source": t(p, "source")}
        for p in d.findall(f"{NS}experimental-properties/{NS}property")]
    rec["dosage_forms"] = [{"Form": t(x, "form"), "Route": t(x, "route"), "Strength": t(x, "strength")}
                           for x in d.findall(f"{NS}dosages/{NS}dosage")]
    rec["moa_targets"] = [{"Target": t(x, "name"), "Actions": ", ".join(a.text for a in x.findall(f"{NS}actions/{NS}action"))}
                          for x in d.findall(f"{NS}targets/{NS}target")]
    rec["atc_codes_text"] = " | ".join(a.get("code") for a in d.findall(f"{NS}atc-codes/{NS}atc-code"))
    return rec, set(ids)


def main(path):
    wanted = {d["drugbank_id"] for d in load_drugs() if d.get("drugbank_id")}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    found = 0
    with open_xml(path) as fh, open(OUT, "w", encoding="utf-8") as out:
        depth = 0
        for event, el in ET.iterparse(fh, events=("start", "end")):
            if el.tag == NS + "drug":
                depth += 1 if event == "start" else -1
                # only top-level <drug> elements (nested ones appear inside drug-interactions)
                if event == "end" and depth == 0:
                    rec, ids = drug_record(el)
                    if ids & wanted:
                        out.write(json.dumps(rec, ensure_ascii=False) + "\n")
                        found += 1
                    el.clear()
    print(f"wrote {OUT}: {found}/{len(wanted)} drugs")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
