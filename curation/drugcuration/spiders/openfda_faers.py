"""FDA Adverse Event Reporting System (FAERS) via the official openFDA API (same data as the FAERS
public dashboard, which is a Qlik JS app that Scrapy cannot render).

For each drug: total reports, top routes, top 10 MedDRA reaction terms overall and per top-3 route.
Counts are REPORT counts, not incidence frequencies.

Usage: scrapy crawl openfda -O output/openfda.jsonl     (optional env OPENFDA_API_KEY raises rate limits)
"""
import json
import os
from urllib.parse import urlencode

import scrapy

from drugcuration.inputs import load_drugs

BASE = "https://api.fda.gov/drug/event.json"
# FAERS route codes (E2B): https://open.fda.gov/apis/drug/event/searchable-fields/
ROUTE_CODES = {"048": "oral", "042": "intravenous", "058": "subcutaneous", "030": "intramuscular",
               "061": "topical", "067": "transdermal", "054": "rectal", "065": "unknown",
               "026": "intra-arterial", "047": "ophthalmic", "032": "intranasal"}


def q(params):
    key = os.environ.get("OPENFDA_API_KEY")
    if key:
        params = {**params, "api_key": key}
    return f"{BASE}?{urlencode(params)}"


# INN -> US adopted / label names used in FAERS
US_NAMES = {
    "dicycloverine": ["dicyclomine"], "mesalazine": ["mesalamine"], "mercaptamine": ["cysteamine"],
    "ursodeoxycholic acid": ["ursodiol"], "levocarnitine": ["levocarnitine", "carnitine"],
    "glycerol phenylbutyrate": ["glycerol phenylbutyrate", "ravicti"], "sodium phenylbutyrate": ["buphenyl"],
    "amphotericin b": ["amphotericin b liposome", "amphotericin b lipid complex"],
    "methylnaltrexone bromide": ["methylnaltrexone"], "hyoscyamine": ["hyoscyamine sulfate"],
    "sapropterin": ["sapropterin dihydrochloride"], "tegaserod": ["tegaserod maleate", "zelnorm"],
    "dolasetron": ["dolasetron mesylate", "anzemet"], "lorcaserin": ["lorcaserin hydrochloride", "belviq"],
    "rosiglitazone": ["rosiglitazone maleate", "avandia"], "nabilone": ["cesamet"],
}


class OpenFdaSpider(scrapy.Spider):
    name = "openfda"
    allowed_domains = ["api.fda.gov"]
    custom_settings = {"DOWNLOAD_DELAY": 0.4, "ROBOTSTXT_OBEY": False,  # API endpoint, no robots.txt
                       "HTTPERROR_ALLOWED_CODES": [404]}

    def __init__(self, input=None, *a, **kw):
        super().__init__(*a, **kw)
        self.drugs = load_drugs(input)

    async def start(self):
        for r in self._requests():
            yield r

    def start_requests(self):
        yield from self._requests()

    @staticmethod
    def search_term(name):
        """Match the drug by its openFDA-harmonised generic name OR the reporter's free-text product name,
        including US adopted names (INN 'mesalazine' = USAN 'mesalamine'). The free-text field is needed for
        drugs with no current US label (e.g. gliclazide, vildagliptin, rosiglitazone), which carry no
        openfda annotation. The suspect-drug restriction (drugcharacterization:1) cannot be tied to the same
        array element in openFDA, so concomitant mentions are included; merge.py records that caveat."""
        names = [name] + US_NAMES.get(name.lower(), [])
        parts = []
        for n in names:
            n = n.upper()
            parts += [f'patient.drug.openfda.generic_name:"{n}"', f'patient.drug.medicinalproduct:"{n}"']
        return "(" + " OR ".join(parts) + ")"

    def _requests(self):
        for d in self.drugs:
            s = self.search_term(d["generic_name"])
            yield scrapy.Request(q({"search": s, "count": "patient.drug.drugadministrationroute.exact"}),
                                 cb_kwargs={"drug": d, "s": s}, callback=self.parse_routes)

    def parse_routes(self, response, drug, s):
        data = json.loads(response.text) if response.status == 200 else {}
        routes = data.get("results", [])
        item = {"generic_name": drug["generic_name"], "drugbank_id": drug.get("drugbank_id"),
                "search": s, "routes": [{"code": r["term"], "route": ROUTE_CODES.get(r["term"], r["term"]),
                                          "reports": r["count"]} for r in routes[:6]],
                "source_urls": [response.url], "reactions": {}}
        if not routes:
            item["status"] = "no_reports"
            yield item
            return
        yield scrapy.Request(q({"search": s, "count": "patient.reaction.reactionmeddrapt.exact", "limit": 25}),
                             cb_kwargs={"item": item, "pending": [r["term"] for r in routes[:3]], "label": "all"},
                             callback=self.parse_reactions)

    def parse_reactions(self, response, item, pending, label):
        data = json.loads(response.text) if response.status == 200 else {}
        item["reactions"][label] = [{"term": r["term"], "reports": r["count"]} for r in data.get("results", [])]
        item["source_urls"].append(response.url)
        if pending:
            code, rest = pending[0], pending[1:]
            s = f'{item["search"]} AND patient.drug.drugadministrationroute:"{code}"'
            yield scrapy.Request(q({"search": s, "count": "patient.reaction.reactionmeddrapt.exact", "limit": 25}),
                                 cb_kwargs={"item": item, "pending": rest,
                                            "label": ROUTE_CODES.get(code, code)},
                                 callback=self.parse_reactions, dont_filter=True)
        else:
            item["status"] = "ok"
            yield item
