# curation/ — Scrapy + scrapy-playwright pipeline for the ATC-A small-molecule dataset

This pipeline fills the incomplete `Drug_Bank_Dataset.tsv.xlsx` for ATC group A (60 small molecules) and writes
`output/ATC_A_curated.tsv`. It produces one row per drug × mode of administration, with 69 columns, a provenance log and an evidence-gap log.
The layout follows the [scrapy/quotesbot](https://github.com/scrapy/quotesbot) template: `scrapy.cfg`, `drugcuration/settings.py`
and `drugcuration/spiders/`.

| Step | Tool | Source | Feeds columns |
|---|---|---|---|
| `fetch_pubchempy.py` | pubchempy + RDKit | pubchem.ncbi.nlm.nih.gov | SMILES (pubchempy only), formula, MW, monoisotopic mass, IUPAC, XLogP, RDKit physchem panel |
| `scrapy crawl drugbank` | Scrapy + **scrapy-playwright** (Chromium) | go.drugbank.com/drugs/DBxxxxx | ClassyFire taxonomy, substituents (pharmacophore), MP, solubility, logP, toxicity text, PK/PD, dosage routes, targets |
| `drugbank_xml.py` (optional) | ElementTree | DrugBank full XML release (free academic account) | same fields as the drugbank spider, without loading web pages |
| `scrapy crawl chembl` | Scrapy (REST JSON) | www.ebi.ac.uk/chembl/api | ChEMBL ID, approval / withdrawn / BBW / orphan / Ro5 / prodrug flags, MoA, activity + toxicity assays, ED50 for TI |
| `scrapy crawl openfda` | Scrapy (REST JSON) | api.fda.gov/drug/event (FAERS) | highest-frequency ADR per route (report counts) |
| `scrapy crawl whocc` | Scrapy | atcddd.fhi.no | WHO DDD, unit, administration route |
| `merge.py` | pandas | — | `output/ATC_A_curated.tsv` |
| `run_pipeline.py` | — | — | runs all of the above from **Spyder (F5)**, Jupyter or a terminal |

## Step by step on Windows with Spyder

The commands below go in the **Anaconda Prompt** (Start menu → "Anaconda Prompt"). Install
[Anaconda](https://www.anaconda.com/download) or Miniconda first if you don't have it.

**1. Get the code** into your internship folder:
```
cd /d "D:\Downloads\Drug Toxicity Internship"
git clone https://github.com/DhruvMenon2003/drug-toxicity-qsar-prediction.git
cd drug-toxicity-qsar-prediction\curation
```
(No git? On GitHub click **Code → Download ZIP**, unzip into the same folder, and `cd` into `curation`.)

**2. Create an environment and install everything.** This is a one-time step of about 5 minutes.
```
conda create -n drugcur python=3.11 -y
conda activate drugcur
conda install -c conda-forge rdkit spyder -y
pip install -r requirements.txt
playwright install chromium
```
`playwright install chromium` downloads the browser that scrapy-playwright drives (about 150 MB).
Spyder is installed *inside* `drugcur` so that its console sees Scrapy, RDKit and Playwright. This avoids
the "spyder-kernels version mismatch" error you get when an outside Spyder points at a different environment.

**3. Check the install** (you should see 4 passed):
```
python -m pytest -q tests
```

**4. Start Spyder from the same prompt:**
```
spyder
```

**5. Open `run_pipeline.py`** with *File → Open…* (it is in the `curation` folder). In the CONFIG block at the top, set:
```python
CONTACT_EMAIL = "your.name@university.edu"
```

**6. Run a quick DrugBank test first.** Temporarily set:
```python
STEPS = ["drugbank"]
DRUGBANK_IDS = "DB01079"
PLAYWRIGHT_HEADLESS = False
```
Press **F5**. A Chromium window opens on the tegaserod page. If DrugBank shows a "Verify you are human"
box, tick it in that window; you have 45 seconds. Your choice is kept in `.pw-profile/` for later runs. When it finishes,
`output\drugbank.jsonl` should contain `"status": "ok"`.
* A status of `"error"` (Timeout) or `"blocked_403"` means the site check was not passed. Run it again with the window visible.

**7. Run the full pipeline.** Set it back to:
```python
STEPS = ["pubchem", "drugbank", "chembl", "openfda", "whocc", "merge"]
DRUGBANK_IDS = ""
PLAYWRIGHT_HEADLESS = True   # or False if step 6 needed a click
```
Press **F5**. It takes about 15–20 minutes for 60 drugs, most of it the polite 4-second delay between DrugBank pages.
Progress prints in the IPython console. Every page is cached in `httpcache\`, so a re-run after an
error only fetches what is missing.

**8. Open the result.** The result is `output\ATC_A_curated.tsv`. After the run, `curated` (a pandas DataFrame) appears
in Spyder's **Variable Explorer**; double-click it to browse. Read the `Evidence_Conflicts` column to see why
a cell is empty.

**9. Add ADMETsar 3.0** (optional; it has no API). Open the admetsar3 web server and paste the SMILES column from
`output\pubchem.jsonl` (or the TSV) into its batch mode. Download the result as CSV and save it as
`inputs\admetsar3_results.csv`. Then set `STEPS = ["merge"]` and press F5 again.

**10. Use the official DrugBank XML instead of scraping** (optional, and the most complete route). With a free academic
account, download *drugbank_all_full_database.xml.zip* from go.drugbank.com/releases/latest. Then set
`DRUGBANK_XML = r"D:\Downloads\drugbank_all_full_database.xml.zip"` and run with `STEPS = ["drugbank", "merge"]`.

### Why the runner uses subprocesses
Scrapy's Twisted reactor can start only once per Python process. Spyder's IPython console also runs its own
event loop. Calling `CrawlerProcess` directly in the console therefore works once and then fails with
`ReactorNotRestartable`. `run_pipeline.py` runs every spider with `python -m scrapy crawl …` in a child process, so you
can press F5 as often as you like. On Windows, scrapy-playwright runs Playwright in its own `ProactorEventLoop`
thread, so no event-loop changes are needed.

### Terminal alternative
Run `run_all.bat` from the Anaconda Prompt, or run single steps directly, e.g. `scrapy crawl drugbank -a ids=DB01079 -O output\test.jsonl`.

## How scrapy-playwright is configured (`drugcuration/settings.py`)
* `DOWNLOAD_HANDLERS` routes http/https through `ScrapyPlaywrightDownloadHandler`. Only requests with
  `meta={"playwright": True}` (DrugBank) open Chromium; ChEMBL, openFDA and WHOCC stay on plain Scrapy HTTP.
* `TWISTED_REACTOR = AsyncioSelectorReactor` is required by scrapy-playwright.
* `PageMethod("wait_for_selector", "dt#drugbank-accession-number", timeout=45000)` waits until the real drug page has
  replaced DrugBank's JavaScript "Securing Connection" page before the HTML is handed to `parse()`.
* `PLAYWRIGHT_ABORT_REQUEST` skips images, fonts and trackers, which makes pages faster and puts less load on the site.
* `PLAYWRIGHT_CONTEXTS` uses a persistent profile (`.pw-profile/`) so cookies persist between runs.
* `PLAYWRIGHT_HEADLESS=0` (or `PLAYWRIGHT_HEADLESS = False` in the runner) shows the browser window. `USE_PLAYWRIGHT=0` disables the browser entirely.

## Curation rules implemented in merge.py
* **Orthogonal assays.** Each value needs two records from different ChEMBL assays *and* different source documents, in the same units.
  In vitro and in vivo are never mixed. For Activity, different readouts (e.g. Ki vs IC50) are preferred.
* **Route matching.** In vivo toxicity is matched to the row's route, which is parsed from the assay description.
* **Therapeutic index.** It is computed only from a same-species, same-route, same-units LD50 (or TD50) and ED50 pair, using geometric means.
* **Solubility.** It is converted to mg/L only when the unit is explicit. The raw DrugBank string is always kept, and qualitative terms are flagged.
* **FAERS ADR.** Product-use and outcome terms (Drug ineffective, Off label use, Death, …) and the drug's own indication
  terms (e.g. Colitis ulcerative for mesalazine) are skipped. Known litigation-stimulated signals (CKD with PPIs/H2
  blockers, DKA with SGLT2 inhibitors) are flagged. Drugs with no US label are matched by reporter product name.
* **DrugBank fallback order.** The live scrape is used first, then `output/drugbank_xml.jsonl`, then `output/drugbank_seed_2026-10-07.jsonl`.
  The seed is a capture of the public pages taken in the previous session; it lists only the first 10 substituents.
* **Logs.** Gaps and caveats go in `Evidence_Conflicts`. Every domain contacted goes in `Filtered_Domains_Log`.

## Known limits
* **FAERS.** Counts are spontaneous-report counts, not incidence rates. "Suspect drug" cannot be tied to the same drug entry, so
  reports where the drug was only a concomitant drug are included.
* **ChEMBL.** The API returned HTTP 500 (an EBI outage) on 2026-10-08. Re-run `chembl` + `merge` to fill the ChEMBL columns.
* **Literature.** Full-text curation for the 2nd orthogonal assay, where ChEMBL has only one, is not automated. Those rows are flagged.
* **DrugBank terms.** The website data is for internal, non-commercial research use, and the terms restrict building a
  redistributed database from it. For a published TSV, use the DrugBank academic XML release (CC BY-NC 4.0) and cite DrugBank.

## Reproducible single-target Hansch QSAR

Open [`hansch_drd2.ipynb`](hansch_drd2.ipynb) in a fresh Google Colab runtime and run all cells.
It uses one human target record (DRD2 / CHEMBL217 / UniProt P14416), binding assays (`B`),
and `Ki` measurements expressed only as pChEMBL. No API key or GPU is needed.

The first cell installs and prints exact dependency versions. By default, the notebook
downloads the immutable, SHA-256-verified [`hansch_drd2_cache.zip`](hansch_drd2_cache.zip),
then performs the analysis without live API requests. A local copy of the archive is
also supported. Set `QSAR_OUTPUT_DIR` to change the output folder; the default is
`hansch_drd2/`. `QSAR_LIVE_REFRESH=1` opts into a new live ChEMBL snapshot in a separate
cache, which may yield different results and must be audited as a new dataset.

The frozen dataset contains 250 standardized compounds, each supported by at least two
documents with replicate spread at most one pChEMBL unit. The notebook performs nested
descriptor selection for the scaffold holdout, five-fold scaffold CV, LOO and all 100
Y-scrambles. The requested 80/20 holdout achieves 195/55 (78/22) using the fixed whole-scaffold
selection rule. The final cell repeats the complete pipeline from cache and asserts
byte-identical `model.json` output. The model is **not predictive**, as documented in the
validation results; the full-data fit is an exploratory description rather than a validated
prediction tool.

Committed outputs are under [`results/hansch_drd2/`](results/hansch_drd2/): raw and curated
CSV files, candidate and retained descriptors, exclusion reasons, model and prediction
tables, validation metadata, provenance, five diagnostic figures, and the verdict.
The compact `qsar_drd2.json` is copied without changes into both `data/qsar_drd2.json`
and `site/data/qsar_drd2.json` in `drug-hansch-space`, where the existing home explorer
loads it in its QSAR overlay. Training-model predictions, surface and applicability domain
are distinguished from the full-data descriptive coefficient table.

## One chemical group at a time: `h2_features.py` and `qed_properties.py`

Two scripts fill every column of the "Drug Bank Dataset Copy" layout for a single ATC level-3 chemical
group (drugs that share a target), one row per mode of administration, and then show the
physicochemical columns one property at a time.

```
pip install rdkit mols2grid pandas requests matplotlib
python3 curation/h2_features.py                 # ATC_GROUP=A02BA by default (H2-receptor antagonists)
python3 curation/qed_properties.py --group-features h2_features/features.csv
```

In Colab, under the **updated workflow** heading:

```python
!pip -q install rdkit mols2grid
!wget -q -O h2_features.py https://raw.githubusercontent.com/DhruvMenon2003/drug-toxicity-qsar-prediction/main/curation/h2_features.py
!wget -q -O qed_properties.py https://raw.githubusercontent.com/DhruvMenon2003/drug-toxicity-qsar-prediction/main/curation/qed_properties.py
%env ATC_GROUP=A02BA
%run h2_features.py
%run qed_properties.py --group-features h2_features/features.csv --all-drugs <path to explorer/data/drugs.json>
```

| Script | Writes | What it is for |
|---|---|---|
| `h2_features.py` | `h2_features/features.csv`, `features_long.csv`, `provenance.csv`, `structures.png`, `grid.html`, `qed_distribution.png`, `binding_data.json` | the filled table; every API call logged; a 2D structure sheet to check the group by eye |
| `qed_properties.py` | `qed_properties/property_<name>.png` (one per property), `qed_desirability.png`, `qed_vs_property.png`, `grid_qed.html`, `properties.csv` | each property on its own against the whole curated set, and QED taken apart into the eight desirability functions it averages |
| `explorer/qed_properties_marimo.py` | — | the same thing interactively (marimo + mols2grid + marimo-chem-utils), with a property picker and a percentile slice you can look at as structures |
| `explorer/binding/build_binding.py` | `explorer/binding/data/receptors.json` | the receptor topology and ligand conformers the [Ligand Binding Theatre](../explorer/binding/README.md) animates |

Rules the table follows, with the reason written into its own `Evidence` column:

- a cell stays **empty** unless a source answered for it, and the reason is logged; nothing is filled from memory
- SMILES come from PubChem only
- in vitro activity is the **median pChEMBL at the mechanism target**, and only when at least two documents agree
- in vitro toxicity is the best human **off-target** pChEMBL under the same two-document rule
- a drug whose ChEMBL record has no mechanism row borrows the target the rest of the ATC group agrees on, and the row says so
- at most two routes per drug, intravenous excluded from the two and kept as an extra row (project rules 3 and 4)
- `Therapeutic Index` and `In Vivo Free Drug Concentration` are written as **Excel formulas**, so they compute themselves the moment the in vivo numbers are entered: `=IF(COUNT(X3,Y3)=2,LOG10(X3)-LOG10(1/Y3),"")` and `=IF(COUNT(AA3,AB3,AC3,AD3)=4,1.44*AA3*AB3*AC3/AD3,"")`
