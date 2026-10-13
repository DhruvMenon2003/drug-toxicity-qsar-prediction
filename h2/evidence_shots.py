"""Screenshots of the pages this run drew on, plus three frames of the binding movie.

The feature table is built from APIs rather than scraped pages (they are faster and kinder to the sites), so
the audit trail is `data/provenance_*.csv`. This script adds the visual half of that trail: the same records
as a human sees them in a browser, captured with Playwright, so the numbers in the table can be eyeballed
against the source without re-running anything.

    python h2/evidence_shots.py --out h2/figures/evidence          # source pages
    python h2/evidence_shots.py --out h2/figures/evidence --serve  # also the binding page, via a local server
"""

from __future__ import annotations

import argparse
import contextlib
import functools
import http.server
import socketserver
import threading
from pathlib import Path

CHROMIUM_CANDIDATES = [
    "/opt/pw-browsers/chromium-1194/chrome-linux/chrome",
    "/opt/google/chrome/chrome",
    "/usr/bin/chromium",
]

PAGES = [
    ("pubchem_famotidine", "https://pubchem.ncbi.nlm.nih.gov/compound/5702160", "PubChem record: identity and SMILES"),
    ("dailymed_famotidine", "https://dailymed.nlm.nih.gov/dailymed/search.cfm?labeltype=all&query=famotidine",
     "DailyMed labels: the routes and doses behind the route rows"),
    ("who_atc_a02ba", "https://atcddd.fhi.no/atc_ddd_index/?code=A02BA&showdescription=no",
     "WHO ATC/DDD index: the A02BA group, its codes and defined daily doses"),
    ("chembl_h2_target", "https://www.ebi.ac.uk/chembl/explore/target/CHEMBL1941",
     "ChEMBL target CHEMBL1941: the H2 receptor the activity column is measured at"),
    ("rcsb_7ul3", "https://www.rcsb.org/structure/7UL3", "PDB 7UL3: the receptor structure used for docking"),
]


def chromium() -> str | None:
    return next((p for p in CHROMIUM_CANDIDATES if Path(p).exists()), None)


@contextlib.contextmanager
def static_server(root: Path, port: int = 8777):
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(root))
    with socketserver.TCPServer(("", port), handler) as httpd:
        thread = threading.Thread(target=httpd.serve_forever, daemon=True)
        thread.start()
        try:
            yield f"http://localhost:{port}"
        finally:
            httpd.shutdown()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="h2/figures/evidence")
    ap.add_argument("--serve", action="store_true", help="also capture the local binding page")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    from playwright.sync_api import sync_playwright

    executable = chromium()
    index = []
    with sync_playwright() as pw:
        browser = pw.chromium.launch(executable_path=executable) if executable else pw.chromium.launch()
        page = browser.new_page(viewport={"width": 1440, "height": 950})
        for name, url, caption in PAGES:
            try:
                page.goto(url, wait_until="domcontentloaded", timeout=90_000)
                page.wait_for_timeout(5_000)
                page.screenshot(path=out / f"{name}.png")
                index.append((name, url, caption, "captured"))
            except Exception as exc:
                index.append((name, url, caption, f"not captured: {type(exc).__name__}"))

        if args.serve:
            with static_server(Path(__file__).resolve().parent) as base:
                page.goto(f"{base}/binding/", wait_until="networkidle", timeout=60_000)
                page.wait_for_timeout(3_000)
                for label, frame in [("approach", 25), ("landing", 70), ("bound", 89)]:
                    page.evaluate(f"const s=document.getElementById('scrub'); s.value={frame};"
                                  " s.dispatchEvent(new Event('input'))")
                    page.wait_for_timeout(900)
                    page.screenshot(path=out.parent / f"binding_{label}.png")
                    index.append((f"binding_{label}", "h2/binding/index.html",
                                  f"binding movie, frame {frame} of 89", "captured"))
        browser.close()

    lines = ["# Evidence screenshots", "",
             "Captured with Playwright (Chromium). The feature values themselves come from the APIs logged in",
             "`data/provenance_features.csv` and `data/provenance_docking.csv`; these pages are the same records",
             "as a person would read them.", "",
             "| File | Source | What it shows | Status |", "|---|---|---|---|"]
    lines += [f"| `{n}.png` | {u} | {c} | {s} |" for n, u, c, s in index]
    (out / "README.md").write_text("\n".join(lines) + "\n")
    for row in index:
        print(row[0], row[3])


if __name__ == "__main__":
    main()
