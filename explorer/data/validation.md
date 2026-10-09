# Hansch Space data checks

Built 2026-10-09T18:26:26+00:00 from `drugbank_curated_dedup.csv`.

## Counts

- Drugs (dataset CIDs): 759
- Parent CID differs from the dataset CID (salt or mixture record in the source): 123
- RDKit descriptors (used for the PCA and every fit): 759
- Curated InChIKey reproduced from the curated SMILES (marimo-chem-utils smi2inchi_key): 759/759
- t-SNE structure-map coordinates: 759
- PubChem records found (XLogP3, complexity, cross-checks): 759; with XLogP3: 752
- PubChem requests that failed after retries or were skipped while PubChem was down: 0 (rerun the build to fill them from PubChem)
- 3D conformers from PubChem: 649; from RDKit: 89; none: 21
- Drugs with a ChEMBL activity value: 308
- ATC level-3 names found: 162/162

## Cross-checks (PubChem against RDKit on the curated SMILES)

- Molecular weight: 759 pairs, max difference 1.09 Da; 3 drugs differ by more than 0.5 Da (would mean the parent CID is the wrong record).
- TPSA: 759 pairs, Pearson r 0.998, mean absolute difference 7.52 A^2 (PubChem counts S and P polar surface; RDKit's default, used here, does not).
- logP: XLogP3 against Crippen, 752 pairs, Pearson r 0.903, mean absolute difference 0.92. Two different logP models; the fits use Crippen for every drug, and the largest gaps are listed below.
- 3D conformer heavy-atom count equals the curated SMILES: 738/738

Molecular weight mismatches (PubChem vs RDKit): desirudin (6963.0 vs 6963.5), lepirudin (6985.0 vs 6985.6), technetium.(99mTc).exametazime (96.9 vs 98.0)

Largest XLogP3 vs Crippen gaps:

| Drug | XLogP3 | Crippen |
|---|---|---|
| mipomersen | -36.50 | 4.99 |
| insulin.detemir | -3.50 | -13.38 |
| insulin.glargine | -14.10 | -21.12 |
| liraglutide | -3.40 | -10.16 |
| pegvisomant | -13.10 | -19.37 |
| gadoteridol | -8.20 | -2.16 |
| gadobenic.acid | -7.10 | -1.10 |
| gadobutrol | -9.80 | -3.82 |

PCA on standardised logp, mr, tpsa, mw, hbd, hba, rotb, arom, fsp3 (mr, tpsa, mw, hbd, hba, rotb, arom as log(1 + x)): PC1-3 explain 59.2%, 18.3%, 10.3%.

## Hansch-type fits

Activity is the median ChEMBL pChEMBL at the drug's mechanism target(s). Drugs are grouped by their first listed target. Only drugs up to 1000 Da are fitted (6 larger drugs with an activity value, mostly peptides, are left out: Crippen logP reaches -43 for them and would set every slope). These groups are not congeneric series, so the fits are exploratory. A fit is credible only when q2 is well above 0 and the shuffled-activity r2 (y-randomisation, 100 runs) stays far below the real r2.

| Group | n | Model | r2 | q2 (LOO) | s | shuffled r2 mean / max |
|---|---|---|---|---|---|---|
| All drugs with an activity value (unrelated targets) (ALL) | 299 | logP + MR + TPSA | 0.153 | 0.126 | 1.123 | 0.01 / 0.031 |
| Serotonin transporter inhibitor (CHEMBL228) | 14 | logP | 0.371 | 0.216 | 0.836 | 0.073 / 0.563 |
| Beta-1 adrenergic receptor antagonist (CHEMBL213) | 11 | TPSA | 0.214 | -0.1 | 1.074 | 0.093 / 0.575 |
| Muscarinic acetylcholine receptor M3 antagonist (CHEMBL245) | 10 | TPSA | 0.176 | -0.271 | 1.181 | 0.094 / 0.717 |
| Human immunodeficiency virus type 1 reverse transcriptase inhibitor (CHEMBL247) | 9 | logP | 0.154 | -0.431 | 1.074 | 0.101 / 0.622 |
| Glucocorticoid receptor agonist (CHEMBL2034) | 8 | logP | 0.358 | -0.143 | 0.593 | 0.134 / 0.602 |
| Serotonin 2a (5-HT2a) receptor antagonist (CHEMBL224) | 8 | MR | 0.102 | -0.275 | 0.799 | 0.141 / 0.656 |
| Mu opioid receptor antagonist (CHEMBL233) | 8 | logP | 0.673 | 0.424 | 0.32 | 0.161 / 0.843 |
| Human immunodeficiency virus type 1 protease inhibitor (CHEMBL243) | 8 | TPSA | 0.601 | 0.196 | 0.446 | 0.162 / 0.632 |
| Beta-2 adrenergic receptor agonist (CHEMBL210) | 7 | TPSA | 0.819 | 0.689 | 0.884 | 0.18 / 0.814 |
| Type-1 angiotensin II receptor antagonist (CHEMBL227) | 6 | logP | 0.584 | 0.084 | 0.26 | 0.192 / 0.677 |
| HMG-CoA reductase inhibitor (CHEMBL402) | 6 | TPSA | 0.082 | -0.566 | 0.635 | 0.209 / 0.776 |
| Farnesyl diphosphate synthase inhibitor (CHEMBL1782) | 5 | TPSA | 0.575 | -0.985 | 0.487 | 0.268 / 0.998 |
| Angiotensin-converting enzyme inhibitor (CHEMBL1808) | 5 | MR | 0.03 | -1.648 | 0.912 | 0.239 / 0.824 |
| Tyrosine-protein kinase ABL inhibitor (CHEMBL1862) | 5 | logP | 0.165 | -1.557 | 0.891 | 0.226 / 0.935 |
| Androgen Receptor agonist (CHEMBL1871) | 5 | logP | 0.058 | -1.685 | 1.076 | 0.26 / 0.783 |
| Serotonin 1b (5-HT1b) receptor agonist (CHEMBL1898) | 5 | logP | 0.316 | -0.56 | 0.387 | 0.213 / 0.908 |
| Dopamine D2 receptor agonist (CHEMBL217) | 5 | logP | 0.478 | -0.552 | 0.537 | 0.25 / 0.958 |
| Histamine H1 receptor antagonist (CHEMBL231) | 5 | TPSA | 0.312 | -0.229 | 0.968 | 0.236 / 0.815 |
| Dopamine transporter inhibitor (CHEMBL238) | 5 | MR | 0.668 | 0.301 | 0.545 | 0.223 / 0.988 |

## Requests

| Host | Status | From cache | Count |
|---|---|---|---|
| atcddd.fhi.no | 200 | True | 162 |
| pubchem.ncbi.nlm.nih.gov | 200 | True | 149 |
