"""Shared helpers for the H2-antagonist feature run: cached HTTP with a provenance log.

Every request any script in `h2/` makes goes through `get()`, so `h2/data/provenance.csv` is a complete
record of what was fetched, when, and with what status. Responses are cached under `--cache` so a second
run is offline and reproduces the same table.
"""

from __future__ import annotations

import csv
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

CACHE = Path(".h2_cache")
PROV: list[dict] = []
SESSION = requests.Session()
SESSION.headers["User-Agent"] = "drug-toxicity-qsar-prediction/h2 (research; contact via GitHub)"


def set_cache(path: str | Path) -> None:
    global CACHE
    CACHE = Path(path)
    CACHE.mkdir(parents=True, exist_ok=True)


def get(url: str, params: dict | None = None, *, note: str = "", tries: int = 3, timeout: int = 60):
    """GET with on-disk cache, retries and a provenance row. Returns a requests.Response or None."""
    key = hashlib.sha1(f"{url}?{sorted((params or {}).items())}".encode()).hexdigest()[:16]
    blob = CACHE / f"{key}.json"
    if blob.exists():
        d = json.loads(blob.read_text())
        PROV.append({"when": d["when"], "url": url, "params": json.dumps(params or {}), "status": d["status"],
                     "note": note, "from_cache": True})
        r = requests.Response()
        r.status_code = d["status"]
        r._content = d["text"].encode()
        r.url = url
        return r if d["status"] else None
    last = None
    for attempt in range(tries):
        try:
            last = SESSION.get(url, params=params, timeout=timeout)
            if last.status_code < 500:
                break
        except Exception as exc:  # network flake
            last, _ = None, exc
        time.sleep(2 ** attempt)
    when = datetime.now(timezone.utc).isoformat(timespec="seconds")
    status = last.status_code if last is not None else 0
    PROV.append({"when": when, "url": url, "params": json.dumps(params or {}), "status": status,
                 "note": note, "from_cache": False})
    if last is not None:
        CACHE.mkdir(parents=True, exist_ok=True)
        blob.write_text(json.dumps({"when": when, "status": status, "text": last.text}))
    return last if last is not None and last.ok else None


def write_provenance(path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["when", "url", "params", "status", "note", "from_cache"])
        w.writeheader()
        w.writerows(PROV)


# The five A02BA (H2-receptor antagonist) drugs in the source table's chemical group, in ATC order.
H2_DRUGS = [
    {"name": "cimetidine", "atc": "A02BA01", "cid": 2756},
    {"name": "ranitidine", "atc": "A02BA02", "cid": 3001055},
    {"name": "famotidine", "atc": "A02BA03", "cid": 5702160},
    {"name": "nizatidine", "atc": "A02BA04", "cid": 4513},
    {"name": "roxatidine", "atc": "A02BA06", "cid": 91276},
]

ATC_LEVELS = {
    "A": "Alimentary tract and metabolism",
    "A02": "Drugs for acid related disorders",
    "A02B": "Drugs for peptic ulcer and gastro-oesophageal reflux disease (GORD)",
    "A02BA": "H2-receptor antagonists",
}

H2_TARGET = "CHEMBL1941"  # Histamine H2 receptor, Homo sapiens
