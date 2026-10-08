"""DrugBank public drug pages -> taxonomy, pharmacophore (substituents), experimental properties,
toxicity text, PK text, dosage forms/routes, predicted ADMET.

Selectors verified against go.drugbank.com/drugs/DB01079 on 2026-10-07:
  <dt id="kingdom"> ... <dd>, <dt id="substituents"> -> span.separated-list-item,
  table#experimental-properties, table#dosages (Form, Route, Strength), table#drug-predicted-admet.

Rendering: scrapy-playwright headless Chromium (see settings.py) because DrugBank serves a JS challenge.

Usage:  scrapy crawl drugbank -O output/drugbank.jsonl
        scrapy crawl drugbank -a ids=DB01079,DB00331 -O output/test.jsonl
"""
import re

import scrapy

try:
    from scrapy_playwright.page import PageMethod
except ImportError:  # scrapy-playwright not installed -> plain HTTP (will usually be challenged)
    PageMethod = None

from drugcuration.inputs import load_drugs

TEXT_FIELDS = [
    "kingdom", "super-class", "class", "sub-class", "direct-parent", "molecular-framework",
    "toxicity", "pharmacodynamics", "mechanism-of-action", "absorption", "volume-of-distribution",
    "protein-binding", "metabolism", "route-of-elimination", "half-life", "clearance",
    "groups", "contraindications-blackbox-warnings", "adverse-effects", "state",
]
LIST_FIELDS = ["substituents", "alternative-parents"]


def clean(s):
    return re.sub(r"\s+", " ", s or "").strip()


def table_rows(sel):
    if not sel:
        return []
    heads = [clean(" ".join(h.css("::text").getall())) for h in sel.css("thead th")]
    rows = []
    for tr in sel.css("tbody tr"):
        cells = [clean(" ".join(td.css("::text").getall())) for td in tr.css("td")]
        if cells:
            rows.append(dict(zip(heads, cells)) if heads else {"cells": cells})
    return rows


class DrugBankSpider(scrapy.Spider):
    name = "drugbank"
    allowed_domains = ["go.drugbank.com"]
    custom_settings = {"DOWNLOAD_DELAY": 4.0, "CONCURRENT_REQUESTS_PER_DOMAIN": 1,
                       # the first navigation can return 403/503 while the JS challenge runs;
                       # let parse() judge by page content instead of dropping the response
                       "HTTPERROR_ALLOWED_CODES": [403, 503]}

    def __init__(self, ids=None, input=None, *a, **kw):
        super().__init__(*a, **kw)
        if ids:
            self.targets = [{"drugbank_id": i.strip(), "generic_name": ""} for i in ids.split(",")]
        else:
            self.targets = [d for d in load_drugs(input) if d.get("drugbank_id")]

    async def start(self):  # Scrapy >= 2.13
        for r in self._requests():
            yield r

    def start_requests(self):  # older Scrapy
        yield from self._requests()

    def _requests(self):
        for d in self.targets:
            yield scrapy.Request(
                f"https://go.drugbank.com/drugs/{d['drugbank_id']}",
                cb_kwargs={"drug": d},
                meta=self._pw_meta(),
                errback=self.on_error,
            )

    @staticmethod
    def _pw_meta():
        """Render with Chromium and wait until the challenge has cleared and the taxonomy block exists."""
        if PageMethod is None:
            return {}
        return {
            "playwright": True,
            "playwright_page_methods": [
                # the challenge page reloads itself; wait for real drug-page markup (up to 45 s)
                PageMethod("wait_for_selector", "dt#drugbank-accession-number", timeout=45_000),
                PageMethod("wait_for_load_state", "domcontentloaded"),
            ],
        }

    def on_error(self, failure):
        req = failure.request
        yield {
            "drugbank_id": req.cb_kwargs["drug"]["drugbank_id"],
            "status": "error",
            "error": repr(failure.value)[:300],
            "source_urls": [req.url],
        }

    def parse(self, response, drug):
        body = response.text
        has_drug = bool(response.css("dt#drugbank-accession-number, dt#kingdom"))
        if not has_drug and ("Just a moment" in body[:5000] or "Securing Connection" in body[:3000]
                             or response.status in (403, 503)):
            yield {"drugbank_id": drug["drugbank_id"], "status": f"blocked_{response.status}",
                   "source_urls": [response.url]}
            return
        if "/login" in response.url or not response.css("dt#kingdom, dt#drugbank-accession-number"):
            yield {"drugbank_id": drug["drugbank_id"], "status": "no_content_or_login",
                   "source_urls": [response.url]}
            return

        def dd(id_):
            return response.xpath(f'//dt[@id="{id_}"]/following-sibling::dd[1]')

        item = {
            "drugbank_id": drug["drugbank_id"],
            "generic_name": drug.get("generic_name", ""),
            "status": "ok",
            "page_title": clean(response.css("title::text").get()),
        }
        for f in TEXT_FIELDS:
            node = dd(f)
            item[f.replace("-", "_")] = clean(" ".join(node.css("::text").getall())) if node else None
        for f in LIST_FIELDS:
            node = dd(f)
            vals = [clean(x) for x in node.css("span.separated-list-item ::text").getall()] if node else []
            if node and not vals:  # fallback: plain text split on " / "
                vals = [v for v in (clean(t) for t in " ".join(node.css("::text").getall()).split(" / ")) if v]
            item[f.replace("-", "_")] = vals

        item["experimental_properties"] = table_rows(response.css("table#experimental-properties"))
        item["dosage_forms"] = table_rows(response.css("table#dosages"))
        item["approved_products_routes"] = sorted({r.get("Route", "") for r in
                                                   table_rows(response.css("table#approved-products"))} - {""})
        item["predicted_admet"] = table_rows(response.css("table#drug-predicted-admet"))
        item["predicted_properties"] = table_rows(response.css("table#drug-moldb-properties"))
        item["moa_targets"] = table_rows(response.css("table#drug-moa-target-table"))
        item["atc_codes_text"] = clean(" ".join(dd("atc-codes").css("::text").getall())) if dd("atc-codes") else None
        item["source_urls"] = [response.url]
        yield item
