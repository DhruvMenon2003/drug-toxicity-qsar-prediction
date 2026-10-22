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

Open [`hansch_drd2.ipynb`](hansch_drd2.ipynb) in a fresh Google Colab runtime and run every
cell in order. The notebook now explains the complete path from collection to interpretation,
with the exact attrition counts, descriptor-removal decisions, validation rationale and a
worked record trail. Its scope is one human target record (DRD2 / CHEMBL217 / UniProt P14416),
binding assays (`B`), and `Ki` expressed only as pChEMBL (pKi). No API key or GPU is needed.

| Question | Read in the notebook | Inspect in the committed results |
|---|---|---|
| What was collected and excluded? | A1 endpoint, confidence, species and pKi explanation | `selection.json`, `raw_activities.csv`, source provenance |
| What is one compound and one response? | A2 normalization, two-document rule and spread cutoff | `dataset.csv`, `dropped.csv` |
| Why these descriptors and which survived? | A3/A4 descriptor meanings and removed-pair table | `candidate_features.csv`, `features.csv`, `feature_selection.json` |
| Does the association predict unseen compounds? | Validation design and interpretation | `model.json`, `predictions.csv`, five diagnostic figures |
| Can I replay or trace a point to its source? | A6 and the final worked CHEMBL8946 record trail | `reproducibility.json`, cached responses, `source_provenance.csv` |

The frozen ChEMBL 37 snapshot starts with 32,602 target activity rows. Selecting Ki and
the exact eligibility rules leaves 10,132 rows. After 22 missing-structure rows are
removed, 7,980 standardized compound groups are evaluated: 7,558 fail only the document
rule, 79 fail only the spread rule, and 93 fail both. **250 compounds remain**; none is
removed by the 1000 Da cutoff. Those counting units differ: activity rows are not compounds.
Confidence ≥8 is applied exactly as requested. It includes homologous or unspecified-source
assays assigned to the human target, which the notebook reports rather than conceals.

The model begins with 11 whole-molecule descriptors. The full-data and scaffold-training
fits each remove MolMR in favor of MolWt, and MolLogP² in favor of MolLogP, under the
pairwise-correlation rule. Both retain nine descriptors. Classical Hansch substituent π
is distinguished from this notebook's whole-molecule MolLogP proxy; MolMR is explained as
a size/polarizability proxy and QED as a composite, not an independent potency measurement.
Neutralization provides a consistent representation, not a physiological protonation state.

Descriptor selection is repeated within the scaffold holdout, five-fold scaffold CV,
every LOO fold and every one of 100 y-scrambles. The requested 80/20 holdout achieves
195/55 compounds (78/22) under the fixed whole-scaffold rule. The frozen model is
**not predictive**: scaffold-test r² is −0.089821, five-fold scaffold q² is −0.284990 and
LOO q² is 0.071419. The real in-sample r² (0.151520) exceeds the largest scrambled fit
(0.126145), but passing this null check does not repair poor prediction. Coefficient
p-values are exploratory after selection; associations and conditional surface changes
do not establish mechanisms or useful chemical interventions.

The first cell installs and prints exact dependency versions. By default the notebook
downloads the immutable, SHA-256-verified [`hansch_drd2_cache.zip`](hansch_drd2_cache.zip),
then uses cached responses for the entire analysis. A local archive is also supported
via `QSAR_SNAPSHOT_ARCHIVE`. `QSAR_OUTPUT_DIR` changes the output folder (default
`hansch_drd2/`). Leave `QSAR_LIVE_REFRESH` unset to reproduce the frozen result;
`QSAR_LIVE_REFRESH=1` deliberately creates a new live snapshot in a separate default cache.

The final code cell repeats the complete pipeline from cache and requires identical
`model.json` bytes. The committed verified model SHA-256 is
`622ba5d1d93e3ff4eee49b665e0ceeb8645785713ce5660eb8765f2b513aef23`.
The full execution and replay were checked under the pinned Python 3.12 runtime;
the hosted Colab browser runtime was not exercised. Original source requests are preserved
separately and imported once into provenance; later cache reads append audit rows.

Committed evidence is under [`results/hansch_drd2/`](results/hansch_drd2/). To trace one
point, follow InChIKey `AADCDMQTJNYOSS-LBPRGKRZSA-N` from the exported predictions into
`dataset.csv`, then molecule CHEMBL8946 into `raw_activities.csv`: its three documents
support pChEMBL values 10.07, 9.59 and 9.93, a median of 9.93 and a spread of 0.48.
The notebook continues the trail through descriptors, test prediction, domain flag,
cached request metadata and publication IDs. Rejected cases have their own trail in
`dropped.csv`.

The unchanged compact `qsar_drd2.json` is copied into `data/qsar_drd2.json` and
`site/data/qsar_drd2.json` in `drug-hansch-space`. The existing home explorer reads it
in the QSAR overlay. Its primary predictions, surface and applicability domain use the
training-only model; the full-data coefficient table is separately descriptive. Surface
grids hold other descriptors fixed and need not represent chemically possible molecules.

Official methods references are linked at their relevant explanations and collected at
the notebook's end: ChEMBL data definitions, the original Hansch/Fujita and Wildman/Crippen
papers, the original QED paper, RDKit descriptor/standardization documentation, scikit-learn
leakage guidance and statsmodels OLS documentation. The executed code and data remain the
source of all numerical claims above.

### Teaching reference reviewed

The user supplied [Santiago Bolivar’s interactive QSAR notebook](https://github.com/Sbolivar16/QSAR_Tools/blob/cf2514d10038252627d397d2164741077ff5d06f/QSAR_Universal_Interactive.ipynb).
Its stage explanations, curation report and export guide are useful presentation references.
We inspected its source without executing or importing it. Its classification workflow differs
from this fixed DRD2 pKi regression: units are not filtered before naming activity values nM,
and preprocessing is fitted before the holdout split/CV. We retain explicit unit eligibility,
fold-local selection and measured-residual Williams diagnostics. No code or model values changed.

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
