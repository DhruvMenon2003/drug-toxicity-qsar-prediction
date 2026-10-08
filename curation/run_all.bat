@echo off
REM Runs the full ATC-A curation pipeline on Windows. Run from this folder in an Anaconda Prompt.
REM First time only:   conda create -n drugcur python=3.11 -y  &&  conda activate drugcur  &&  pip install -r requirements.txt  &&  playwright install chromium
REM Set your contact address once (sent in the User-Agent so site operators can reach you):
REM   setx CURATION_CONTACT "you@university.edu"

if not exist output mkdir output

echo [1/6] PubChem via pubchempy (SMILES + properties + RDKit descriptors)
python fetch_pubchempy.py || goto :err

echo [2/6] DrugBank public pages via scrapy-playwright Chromium (4 s between requests, ~6 min for 60 drugs)
scrapy crawl drugbank -O output\drugbank.jsonl || goto :err

echo [3/6] ChEMBL REST (metadata, mechanism, bioactivity, toxicity)
scrapy crawl chembl -O output\chembl.jsonl || goto :err

echo [4/6] openFDA FAERS (adverse events by route)
scrapy crawl openfda -O output\openfda.jsonl || goto :err

echo [5/6] WHO ATC/DDD index
scrapy crawl whocc -O output\whocc.jsonl || goto :err

echo [6/6] Merge into output\ATC_A_curated.tsv
python merge.py || goto :err

echo Done. Check output\ATC_A_curated.tsv and the Evidence_Conflicts column.
goto :eof

:err
echo A step failed - see the messages above. Re-running is safe: finished pages are cached in httpcache\.
exit /b 1
