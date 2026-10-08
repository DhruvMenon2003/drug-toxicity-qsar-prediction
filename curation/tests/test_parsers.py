"""Offline tests: DrugBank parser on a fixture built from markup verified on go.drugbank.com (2026-10-07),
and helper functions in merge.py. Run: python -m pytest -q tests  (or python tests/test_parsers.py)"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scrapy.http import HtmlResponse, Request

import merge
from drugcuration.spiders.drugbank import DrugBankSpider

FIXTURE = """<html><head><title>Tegaserod: Uses | DrugBank</title></head><body><dl>
<dt id="drugbank-accession-number">DrugBank Accession Number</dt><dd>DB01079</dd>
<dt id="toxicity">Toxicity</dt><dd class="col-xl-10"><p>Single oral doses of 120 mg were administered. Oral LD50 rat &gt;2000 mg/kg.</p></dd>
<dt id="half-life">Half-life</dt><dd><p>11 +/- 5 hours</p></dd>
<dt id="kingdom">Kingdom</dt><dd class="col-md-8"><a class="classyfire-taxnode" href="x">Organic compounds</a></dd>
<dt id="super-class">Super Class</dt><dd><a>Organoheterocyclic compounds</a></dd>
<dt id="class">Class</dt><dd><a>Indoles and derivatives</a></dd>
<dt id="sub-class">Sub Class</dt><dd><a>Hydroxyindoles</a></dd>
<dt id="direct-parent">Direct Parent</dt><dd><a>Hydroxyindoles</a></dd>
<dt id="substituents">Substituents</dt><dd><span class="separated-list-container"><span class="separated-list">
<span class="separated-list-item">Alkyl aryl ether</span><span class="list-separator"> / </span>
<span class="separated-list-item">Anisole</span><span class="list-separator"> / </span>
<span class="separated-list-item">Azacycle</span></span></span></dd>
<dt id="molecular-framework">Molecular Framework</dt><dd>Aromatic heteropolycyclic compounds</dd>
<dt id="dosage-forms">Dosage Forms</dt><dd><table id="dosages"><thead><tr><th>Form</th><th>Route</th><th>Strength</th></tr></thead>
<tbody><tr><td>Tablet</td><td>Oral</td><td>6 mg</td></tr></tbody></table></dd>
<dt id="experimental-properties">Experimental Properties</dt><dd><table class="table" id="experimental-properties">
<thead><tr><th>Property</th><th>Value</th><th>Source</th></tr></thead><tbody>
<tr><td>melting point (°C)</td><td>155 °C</td><td><span>Not Available</span></td></tr>
<tr><td>logP</td><td>2.6</td><td><span>Not Available</span></td></tr></tbody></table></dd>
</dl></body></html>"""


def test_drugbank_parse():
    sp = DrugBankSpider(ids="DB01079")
    req = Request("https://go.drugbank.com/drugs/DB01079")
    resp = HtmlResponse(url=req.url, body=FIXTURE.encode(), encoding="utf-8", request=req)
    item = next(iter(sp.parse(resp, drug={"drugbank_id": "DB01079", "generic_name": "tegaserod"})))
    assert item["status"] == "ok"
    assert item["kingdom"] == "Organic compounds"
    assert item["super_class"] == "Organoheterocyclic compounds"
    assert item["sub_class"] == "Hydroxyindoles"
    assert item["substituents"] == ["Alkyl aryl ether", "Anisole", "Azacycle"]
    assert item["dosage_forms"] == [{"Form": "Tablet", "Route": "Oral", "Strength": "6 mg"}]
    assert merge.exp_prop(item, "melting point")[0] == "155 °C"
    assert "LD50" in item["toxicity"]
    return item


def test_helpers():
    assert merge.solubility_mg_l("3.56 mg/mL") == 3560
    assert merge.solubility_mg_l(">20 mg/mL") == 20000
    assert merge.solubility_mg_l("Practically insoluble") is None
    assert merge.route_from_description("Acute toxicity in rat assessed as LD50 after oral administration") == "Oral"
    assert merge.route_from_description("LD50 in mouse administered i.v.") == "Intravenous"
    assert merge.route_from_description("Cytotoxicity against HepG2 cells") is None
    assert merge.norm_route("ORAL") == "Oral"
    a = {"standard_value": "1", "assay_chembl_id": "A1", "document_chembl_id": "D1", "standard_type": "Ki", "document_year": 2010}
    b = {"standard_value": "2", "assay_chembl_id": "A2", "document_chembl_id": "D1", "standard_type": "IC50", "document_year": 2009}
    c = {"standard_value": "3", "assay_chembl_id": "A3", "document_chembl_id": "D2", "standard_type": "Ki", "document_year": 2008}
    d = {"standard_value": "4", "assay_chembl_id": "A4", "document_chembl_id": "D3", "standard_type": "IC50", "document_year": 2005}
    picks = merge.pick_orthogonal([a, b, c, d], prefer_distinct_type=True)
    assert [p["assay_chembl_id"] for p in picks] == ["A1", "A4"]  # different assay, document AND readout
    picks = merge.pick_orthogonal([a, b, c])
    assert [p["assay_chembl_id"] for p in picks] == ["A1", "A3"]  # same-document assay b rejected


def test_openfda_search_term():
    from drugcuration.spiders.openfda_faers import OpenFdaSpider
    q = OpenFdaSpider.search_term("mesalazine")
    assert 'openfda.generic_name:"MESALAMINE"' in q and 'medicinalproduct:"MESALAZINE"' in q
    assert q.startswith("(") and q.endswith(")")  # safe to AND with a route filter


def test_playwright_meta():
    meta = DrugBankSpider._pw_meta()
    assert meta.get("playwright") is True
    assert meta["playwright_page_methods"][0].args[0] == "dt#drugbank-accession-number"


if __name__ == "__main__":
    test_drugbank_parse(); test_helpers(); test_openfda_search_term(); test_playwright_meta(); print("all tests passed")
