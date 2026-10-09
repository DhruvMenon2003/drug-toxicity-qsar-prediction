# Drug Hansch Space

A browser page that places every drug in the curated DrugBank table (`curation/colab_curate_all.py`) in a 3D descriptor space built from PubChem, then walks through a Hansch-style QSAR analysis on it. It uses plain Three.js and GSAP from a CDN, so there is no build step.

## What it shows

- **Tour.** Seven scroll steps: the PCA chemical space, then log P, molar refractivity and TPSA one axis at a time (with the Veber 140 Å² line), then ChEMBL potency, a Hansch fit for one target group, and the same fit across all targets.
- **Explore.** Any three descriptors as axes, colour by ATC group (up to three at a time), oral route, potency or rule-of-five violations, and filters for ATC levels 1 to 3 and route.
- **Drug panel.** The PubChem 3D conformer as ball and stick (RDKit ETKDGv3 where PubChem has none), a 2D depiction of the curated SMILES, the descriptors, the ChEMBL activity value and every ATC code with its routes.
- **Hansch fits.** For each ChEMBL mechanism target shared by at least five drugs: the model with the best leave-one-out q² among log P, log P², MR and TPSA terms (one descriptor per five drugs), its statistics, a y-randomisation check, an observed against predicted plot and the fitted surface in 3D.

The fits are exploratory. Drugs grouped by target are not a congeneric series, the activity values are medians across assays, and TPSA is a whole-molecule polarity proxy rather than a Hammett σ. The page labels each fit as predictive, weak or not predictive from q² and the shuffled-activity runs.

## Rebuild the data

```bash
pip install rdkit pandas numpy requests
python explorer/build_data.py --curated path/to/drugbank_curated_dedup.csv --out explorer/data --cache .pubchem_cache
```

The script reads only the curated CSV and public APIs (PubChem PUG REST and the WHO ATC index). It caches every response under `--cache`, so a second run is offline and gives the same output. It writes:

| File | Contents |
|---|---|
| `data/drugs.json` | One record per drug: identifiers, ATC codes and routes, descriptors, PCA scores, ChEMBL activity |
| `data/molecules.json` | 3D coordinates and bonds for the ball-and-stick viewer, loaded when a drug is first opened |
| `data/hansch.json` | Every target-group fit and the all-drugs fit, with candidates and leave-one-out predictions |
| `data/validation.md` | Counts and cross-checks: PubChem against RDKit weight, TPSA and log P; conformer atom counts; fit table |
| `data/provenance.csv` | Every request the build made, with status and time |

`--fix-curated-l3` also writes the corrected ATC level-3 names back into the curated CSV (the earlier curation run stored the WHO page's "Show text from Guidelines" link text instead of the group name).

## View it locally

```bash
cd explorer && python -m http.server 8000
```

Then open http://localhost:8000. The page fetches the JSON files next to it, so opening `index.html` straight from disk will not work.
