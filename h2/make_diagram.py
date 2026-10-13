"""Draw the H2 run as an Excalidraw file (and an SVG), so the pipeline is readable before the code is.

    python h2/make_diagram.py --out h2/figures/h2_workflow

Writes `<out>.excalidraw` (open at excalidraw.com) and `<out>.svg`.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

W, H, GAP_X, GAP_Y = 250, 92, 86, 44
PALETTE = {"source": ("#e7f5ff", "#1971c2"), "step": ("#fff9db", "#e8590c"), "out": ("#ebfbee", "#2f9e44"),
           "check": ("#f3f0ff", "#6741d9")}

COLUMNS: list[tuple[str, str, list[tuple[str, str]]]] = [
    ("Sources", "source", [
        ("PubChem", "CID, isomeric SMILES, XLogP3,\nPK statements (PUG REST + PUG View)"),
        ("ChEMBL", "pChEMBL at the H2 receptor,\nstrongest off-target, molecule properties"),
        ("openFDA + DailyMed", "label routes, dose and interval,\nFAERS reporting signal"),
        ("WHO ATC/DDD", "A02BA codes, defined daily dose,\nroute where no US label exists"),
        ("TDC LD50_Zhu", "acute oral rat LD50\n(ChemIDplus/RTECS origin)"),
        ("RCSB PDB 7UL3", "inactive human H2 receptor\nwith famotidine bound"),
    ]),
    ("Pipeline", "step", [
        ("fetch_features.py", "one row per drug x ATC code x route;\nrules 3, 4 and 5 applied"),
        ("RDKit descriptors", "MW, cLogP, MR, TPSA, QED\nand its eight terms"),
        ("dock.py", "Meeko + AutoDock Vina into the\northosteric site; 20 poses each"),
        ("qed_report.py", "per-property distributions,\nmols2grid structure grid"),
        ("validate.py", "13 cross-checks against a\nsecond source or the structure"),
    ]),
    ("Deliverables", "out", [
        ("h2_features.csv", "the feature table, every cell\nwith its source"),
        ("binding/ movie", "ligand entering the pocket,\nhydrogen bonds snapping in"),
        ("qed_marimo.py", "one property at a time,\nstructures on hover"),
        ("workbook sheets", "appended to the source workbook,\nformulas live"),
    ]),
    ("Checks you can see", "check", [
        ("redocking 1.15 A", "famotidine lands back on its\ncrystallographic pose"),
        ("figures/", "structures, distributions,\nbinding frames"),
        ("provenance_*.csv", "every request, status and time"),
        ("evidence/", "screenshots of the source pages"),
    ]),
]


def rounded(seed: int, x: float, y: float, w: float, h: float, bg: str, stroke: str, index: str) -> dict:
    return {"id": f"r{seed}", "type": "rectangle", "x": x, "y": y, "width": w, "height": h, "angle": 0,
            "strokeColor": stroke, "backgroundColor": bg, "fillStyle": "solid", "strokeWidth": 1.5,
            "strokeStyle": "solid", "roughness": 1, "opacity": 100, "groupIds": [], "roundness": {"type": 3},
            "seed": seed, "version": 1, "versionNonce": seed, "isDeleted": False, "boundElements": [],
            "updated": 1, "link": None, "locked": False, "index": index}


def text(seed: int, x: float, y: float, w: float, content: str, size: int, colour: str, index: str) -> dict:
    lines = content.split("\n")
    return {"id": f"t{seed}", "type": "text", "x": x, "y": y, "width": w, "height": size * 1.25 * len(lines),
            "angle": 0, "strokeColor": colour, "backgroundColor": "transparent", "fillStyle": "solid",
            "strokeWidth": 1, "strokeStyle": "solid", "roughness": 1, "opacity": 100, "groupIds": [],
            "roundness": None, "seed": seed, "version": 1, "versionNonce": seed, "isDeleted": False,
            "boundElements": [], "updated": 1, "link": None, "locked": False, "text": content,
            "fontSize": size, "fontFamily": 2, "textAlign": "left", "verticalAlign": "top",
            "containerId": None, "originalText": content, "lineHeight": 1.25, "index": index}


def arrow(seed: int, x: float, y: float, dx: float, dy: float, index: str) -> dict:
    return {"id": f"a{seed}", "type": "arrow", "x": x, "y": y, "width": abs(dx), "height": abs(dy), "angle": 0,
            "strokeColor": "#868e96", "backgroundColor": "transparent", "fillStyle": "solid", "strokeWidth": 1.2,
            "strokeStyle": "solid", "roughness": 1, "opacity": 100, "groupIds": [], "roundness": {"type": 2},
            "seed": seed, "version": 1, "versionNonce": seed, "isDeleted": False, "boundElements": [],
            "updated": 1, "link": None, "locked": False, "points": [[0, 0], [dx, dy]], "lastCommittedPoint": None,
            "startBinding": None, "endBinding": None, "startArrowhead": None, "endArrowhead": "arrow",
            "index": index}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="h2/figures/h2_workflow")
    args = ap.parse_args()
    random.seed(7)

    elements: list[dict] = []
    order = iter(f"a{i:04d}" for i in range(2000))
    column_x = []
    for c, (heading, kind, boxes) in enumerate(COLUMNS):
        x = c * (W + GAP_X)
        column_x.append(x)
        elements.append(text(random.randint(1, 10**9), x, -56, W, heading, 22, "#212529", next(order)))
        bg, stroke = PALETTE[kind]
        for b, (title, body) in enumerate(boxes):
            y = b * (H + GAP_Y)
            seed = random.randint(1, 10**9)
            elements.append(rounded(seed, x, y, W, H, bg, stroke, next(order)))
            elements.append(text(seed + 1, x + 14, y + 12, W - 28, title, 16, "#212529", next(order)))
            elements.append(text(seed + 2, x + 14, y + 36, W - 28, body, 12, "#495057", next(order)))

    for c in range(len(COLUMNS) - 1):
        elements.append(arrow(random.randint(1, 10**9), column_x[c] + W + 10, 150,
                              GAP_X - 20, 0, next(order)))

    scene = {"type": "excalidraw", "version": 2, "source": "drug-toxicity-qsar-prediction/h2/make_diagram.py",
             "elements": elements, "appState": {"gridSize": None, "viewBackgroundColor": "#ffffff"}, "files": {}}
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.with_suffix(".excalidraw").write_text(json.dumps(scene, indent=1))

    # a plain SVG of the same layout, for anyone who will not open Excalidraw
    width = len(COLUMNS) * (W + GAP_X)
    height = max(len(c[2]) for c in COLUMNS) * (H + GAP_Y) + 90
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
           f'viewBox="-20 -80 {width} {height}" font-family="Inter, system-ui, sans-serif">',
           '<rect x="-20" y="-80" width="100%" height="100%" fill="#ffffff"/>']
    for c, (heading, kind, boxes) in enumerate(COLUMNS):
        x = c * (W + GAP_X)
        bg, stroke = PALETTE[kind]
        svg.append(f'<text x="{x}" y="-40" font-size="20" font-weight="600" fill="#212529">{heading}</text>')
        for b, (title, body) in enumerate(boxes):
            y = b * (H + GAP_Y)
            svg.append(f'<rect x="{x}" y="{y}" width="{W}" height="{H}" rx="10" fill="{bg}" stroke="{stroke}"/>')
            svg.append(f'<text x="{x + 14}" y="{y + 26}" font-size="15" font-weight="600" fill="#212529">{title}</text>')
            for k, line in enumerate(body.split("\n")):
                svg.append(f'<text x="{x + 14}" y="{y + 48 + k * 15}" font-size="11.5" fill="#495057">{line}</text>')
        if c < len(COLUMNS) - 1:
            svg.append(f'<path d="M{x + W + 12} 150 L{x + W + GAP_X - 14} 150" stroke="#868e96" stroke-width="1.5" '
                       f'marker-end="url(#head)"/>')
    svg.insert(1, '<defs><marker id="head" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto">'
                  '<path d="M0,0 L6,3 L0,6 Z" fill="#868e96"/></marker></defs>')
    svg.append("</svg>")
    out.with_suffix(".svg").write_text("\n".join(svg))
    print(f"wrote {out.with_suffix('.excalidraw')} and {out.with_suffix('.svg')}")


if __name__ == "__main__":
    main()
