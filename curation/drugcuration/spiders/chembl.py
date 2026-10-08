"""ChEMBL REST API -> molecule metadata, mechanism, and bioactivities used for Toxicity / Activity assays.

Per drug:
  1. resolve ChEMBL ID by PubChem CID-independent name search (pref_name / synonyms)
  2. molecule record (max_phase, first_approval, black_box_warning, withdrawn_flag, ro5, prodrug...)
  3. mechanism(s) -> target ChEMBL IDs
  4. activities on mechanism targets with pChEMBL values (Activity candidates; in vitro)
  5. toxicity-type activities: LD50/TD50/LC50/CC50/TC50 (Toxicity candidates; in vivo vs in vitro tagged)

Usage: scrapy crawl chembl -O output/chembl.jsonl
"""
import json
from urllib.parse import urlencode

import scrapy

from drugcuration.inputs import load_drugs

API = "https://www.ebi.ac.uk/chembl/api/data"
TOX_TYPES = "LD50,TD50,LC50,CC50,TC50,MTD,NOAEL,LOAEL,ED50"
ACT_FIELDS = ["activity_id", "assay_chembl_id", "assay_type", "assay_description", "assay_organism",
              "assay_cell_type", "bao_label", "standard_type", "standard_relation", "standard_value",
              "standard_units", "pchembl_value", "target_chembl_id", "target_pref_name",
              "target_organism", "document_chembl_id", "document_year", "src_id",
              "data_validity_comment", "potential_duplicate"]


def url(path, **params):
    params.setdefault("format", "json")
    return f"{API}/{path}?{urlencode(params)}"


class ChemblSpider(scrapy.Spider):
    name = "chembl"
    allowed_domains = ["www.ebi.ac.uk"]
    custom_settings = {"DOWNLOAD_DELAY": 0.5, "ROBOTSTXT_OBEY": False}

    def __init__(self, input=None, *a, **kw):
        super().__init__(*a, **kw)
        self.drugs = load_drugs(input)

    async def start(self):
        for r in self._requests():
            yield r

    def start_requests(self):
        yield from self._requests()

    def _requests(self):
        for d in self.drugs:
            name = d["generic_name"]
            yield scrapy.Request(url("molecule", pref_name__iexact=name),
                                 cb_kwargs={"drug": d, "tried_syn": False}, callback=self.parse_lookup)

    def parse_lookup(self, response, drug, tried_syn):
        mols = json.loads(response.text).get("molecules", [])
        if not mols and not tried_syn:
            yield scrapy.Request(url("molecule", molecule_synonyms__molecule_synonym__iexact=drug["generic_name"]),
                                 cb_kwargs={"drug": drug, "tried_syn": True}, callback=self.parse_lookup)
            return
        if not mols:
            yield {"generic_name": drug["generic_name"], "status": "not_found", "source_urls": [response.url]}
            return
        # prefer the parent (non-salt) form
        mols.sort(key=lambda m: (m.get("molecule_hierarchy") or {}).get("parent_chembl_id") != m["molecule_chembl_id"])
        m = mols[0]
        props = m.get("molecule_properties") or {}
        item = {
            "generic_name": drug["generic_name"], "drugbank_id": drug.get("drugbank_id"),
            "chembl_id": m["molecule_chembl_id"], "pref_name": m.get("pref_name"),
            "max_phase": m.get("max_phase"), "first_approval": m.get("first_approval"),
            "black_box_warning": m.get("black_box_warning"), "withdrawn_flag": m.get("withdrawn_flag"),
            "orphan": m.get("orphan"), "prodrug": m.get("prodrug"), "oral": m.get("oral"),
            "parenteral": m.get("parenteral"), "topical": m.get("topical"),
            "alogp": props.get("alogp"), "mw_freebase": props.get("mw_freebase"),
            "num_ro5_violations": props.get("num_ro5_violations"), "psa": props.get("psa"),
            "hba": props.get("hba"), "hbd": props.get("hbd"), "rtb": props.get("rtb"),
            "aromatic_rings": props.get("aromatic_rings"), "qed_weighted": props.get("qed_weighted"),
            "source_urls": [response.url],
        }
        yield scrapy.Request(url("mechanism", molecule_chembl_id=item["chembl_id"]),
                             cb_kwargs={"item": item}, callback=self.parse_mech)

    def parse_mech(self, response, item):
        mechs = json.loads(response.text).get("mechanisms", [])
        item["mechanisms"] = [{"mechanism_of_action": x.get("mechanism_of_action"),
                               "action_type": x.get("action_type"),
                               "target_chembl_id": x.get("target_chembl_id")} for x in mechs]
        item["source_urls"].append(response.url)
        targets = sorted({x["target_chembl_id"] for x in item["mechanisms"] if x["target_chembl_id"]})
        params = dict(molecule_chembl_id=item["chembl_id"], pchembl_value__isnull="false", limit=1000)
        if targets:
            params["target_chembl_id__in"] = ",".join(targets)
        yield scrapy.Request(url("activity", **params), cb_kwargs={"item": item}, callback=self.parse_act)

    def parse_act(self, response, item):
        acts = json.loads(response.text).get("activities", [])
        item["activity_records"] = [{k: a.get(k) for k in ACT_FIELDS} for a in acts]
        item["source_urls"].append(response.url)
        yield scrapy.Request(url("activity", molecule_chembl_id=item["chembl_id"],
                                 standard_type__in=TOX_TYPES, limit=1000),
                             cb_kwargs={"item": item}, callback=self.parse_tox)

    def parse_tox(self, response, item):
        acts = json.loads(response.text).get("activities", [])
        item["toxicity_records"] = [{k: a.get(k) for k in ACT_FIELDS} for a in acts]
        item["source_urls"].append(response.url)
        item["status"] = "ok"
        yield item
