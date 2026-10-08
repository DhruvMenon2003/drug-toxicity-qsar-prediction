"""WHO Collaborating Centre ATC/DDD index -> DDD, unit, administration route per ATC code.
Academic use with citation of WHOCC; commercial reuse of the index is not permitted.

Usage: scrapy crawl whocc -O output/whocc.jsonl
"""
import re

import scrapy

from drugcuration.inputs import load_drugs

ADM_R = {"O": "oral", "P": "parenteral", "R": "rectal", "SL": "sublingual", "TD": "transdermal",
         "N": "nasal", "V": "vaginal", "Inhal": "inhalation", "Chewing gum": "chewing gum",
         "Inhal.powder": "inhalation powder", "Inhal.solution": "inhalation solution",
         "Inhal.aerosol": "inhalation aerosol", "implant": "implant", "lamella": "lamella",
         "ointment": "ointment", "oral aerosol": "oral aerosol", "s.c. implant": "s.c. implant",
         "urethral": "urethral", "intravesical": "intravesical", "instill.sol.": "instillation"}


def clean(s):
    return re.sub(r"\s+", " ", s or "").strip()


class WhoccSpider(scrapy.Spider):
    name = "whocc"
    allowed_domains = ["atcddd.fhi.no", "www.whocc.no"]
    custom_settings = {"DOWNLOAD_DELAY": 3.0}

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
            for code in str(d["atc_code"]).split("|"):
                code = code.strip()
                if code.startswith("A"):
                    yield scrapy.Request(f"https://atcddd.fhi.no/atc_ddd_index/?code={code}&showdescription=no",
                                         cb_kwargs={"drug": d, "code": code})

    def parse(self, response, drug, code):
        rows = []
        for tr in response.css("table tr"):
            cells = [clean(" ".join(td.css("::text").getall())) for td in tr.css("td")]
            if len(cells) >= 5 and cells[0].upper() != "ATC CODE":
                rows.append(cells)
        # rows: ATC code | Name | DDD | U | Adm.R | Note ; continuation rows have empty first cells
        ddds, cur = [], None
        for c in rows:
            if c[0]:
                cur = c[0]
            if cur == code and len(c) >= 5:
                ddds.append({"atc_code": code, "name": c[1] or drug["generic_name"], "ddd": c[2], "unit": c[3],
                             "adm_route_code": c[4], "adm_route": ADM_R.get(c[4], c[4]),
                             "note": c[5] if len(c) > 5 else ""})
        yield {"generic_name": drug["generic_name"], "atc_code": code, "ddd_rows": ddds,
               "status": "ok" if ddds else "no_ddd", "source_urls": [response.url]}
