"""One-click runner for Spyder (press F5), Jupyter, or a terminal.

Why subprocesses?  Scrapy runs on a Twisted reactor that can be started only ONCE per Python process,
and Spyder's IPython console already runs its own asyncio event loop. Calling CrawlerProcess inside
the console works the first time and then fails with `ReactorNotRestartable`. Running each spider
as `python -m scrapy crawl ...` in a child process avoids both problems, and every step can be
re-run as often as you like.

Edit the CONFIG block, then press F5 in Spyder (or run `python run_pipeline.py`).
"""
import os
import subprocess
import sys
import time
from pathlib import Path

# ----------------------------- CONFIG ---------------------------------------
CONTACT_EMAIL = "you@university.edu"   # sent in the User-Agent so site operators can reach you
STEPS = ["pubchem", "drugbank", "chembl", "openfda", "whocc", "merge"]   # remove steps to skip them
DRUGBANK_IDS = ""          # e.g. "DB01079,DB00331" to test a few drugs; "" = all drugs in inputs/
PLAYWRIGHT_HEADLESS = True  # False opens a visible Chromium window for DrugBank (use if challenged)
DRUGBANK_XML = ""          # path to drugbank_all_full_database.xml.zip if you have an academic account
OPENFDA_API_KEY = ""       # optional; free key at https://open.fda.gov/apis/authentication/
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"


def run(label, args):
    env = dict(os.environ, CURATION_CONTACT=CONTACT_EMAIL,
               PLAYWRIGHT_HEADLESS="1" if PLAYWRIGHT_HEADLESS else "0")
    if OPENFDA_API_KEY:
        env["OPENFDA_API_KEY"] = OPENFDA_API_KEY
    print(f"\n=== {label}: {' '.join(args)}", flush=True)
    t0 = time.time()
    proc = subprocess.run([sys.executable, *args], cwd=HERE, env=env)
    print(f"=== {label} finished in {time.time() - t0:.0f} s (exit code {proc.returncode})", flush=True)
    if proc.returncode != 0:
        print(f"!!! {label} failed. Fix the error above and re-run; finished pages are cached in httpcache/.")
    return proc.returncode == 0


def crawl(name, outfile, *extra):
    return run(name, ["-m", "scrapy", "crawl", name, "-O", str(OUT / outfile), *extra])


def main():
    OUT.mkdir(exist_ok=True)
    ok = True
    for step in STEPS:
        if step == "pubchem":
            ok &= run("PubChem (pubchempy)", ["fetch_pubchempy.py"])
        elif step == "drugbank":
            if DRUGBANK_XML:
                ok &= run("DrugBank XML", ["drugbank_xml.py", DRUGBANK_XML])
            extra = ["-a", f"ids={DRUGBANK_IDS}"] if DRUGBANK_IDS else []
            ok &= crawl("drugbank", "drugbank.jsonl", *extra)
        elif step in ("chembl", "openfda", "whocc"):
            ok &= crawl(step, f"{step}.jsonl")
        elif step == "merge":
            ok &= run("Merge", ["merge.py"])
    tsv = OUT / "ATC_A_curated.tsv"
    print("\nAll steps OK." if ok else "\nSome steps failed (see above).",
          f"Curated table: {tsv}" if tsv.exists() else "")
    return tsv


if __name__ == "__main__":
    curated_path = main()
    # In Spyder this leaves `curated` in the Variable Explorer for inspection.
    try:
        import pandas as pd
        curated = pd.read_csv(curated_path, sep="\t")
        print(curated.shape)
    except Exception as e:  # file missing if merge was skipped
        print("Could not load the curated TSV:", e)
