# Internship Project: Predicting Drug Toxicity Endpoints using QSAR 

Multi-stage ML pipeline for toxicity prediction with mathematical modeling and parallel computing optimization.

## High Level Workflow
<img width="2591" height="2666" alt="High Level Workflow" src="https://github.com/user-attachments/assets/c2034162-5a96-431b-bed8-893a49675e94" />



## Curated drug dataset (ATC group A, small molecules)
`curation/` contains a Scrapy + scrapy-playwright pipeline that fills the drug-level dataset (PubChem, DrugBank, ChEMBL, FAERS, WHO ATC/DDD) with orthogonal-assay and provenance rules. See [curation/README.md](curation/README.md) for the step-by-step Spyder instructions and section 4 of `Exploratory Data Analysis.ipynb` for the analysis.
