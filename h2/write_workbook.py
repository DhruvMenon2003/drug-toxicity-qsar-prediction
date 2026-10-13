"""Add this run's sheets to a copy of the source workbook, leaving every existing sheet as it was.

The source workbook keeps its sheets, its values and its formulas; four sheets are appended:

| Sheet | Contents |
|---|---|
| `H2 per route` | the feature table, one row per drug x ATC code x route, with live formulas for the therapeutic index and the free drug concentration |
| `H2 QED terms` | QED, its eight terms and each drug's percentile among the curated drugs |
| `H2 docking` | the Vina run: score, selected pose, contacts, redocking check |
| `H2 validation` | the cross-checks from `validate.py` |

    python h2/write_workbook.py --source "Copy of Drug Bank Dataset.xlsx" --out out.xlsx
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

HEADER_FILL = PatternFill("solid", fgColor="1F3864")
SOURCE_FILL = PatternFill("solid", fgColor="EEF2F9")


def write_sheet(book, title: str, frame: pd.DataFrame, formulas: dict[str, str] | None = None) -> None:
    if title in book.sheetnames:
        del book[title]
    sheet = book.create_sheet(title)
    for j, column in enumerate(frame.columns, start=1):
        cell = sheet.cell(row=1, column=j, value=column)
        cell.font = Font(bold=True, color="FFFFFF", size=10)
        cell.fill = HEADER_FILL
        cell.alignment = Alignment(wrap_text=True, vertical="top")
        sheet.column_dimensions[get_column_letter(j)].width = 42 if column.endswith(("_source", "_Edit_Flag")) else 17
    letters = {column: get_column_letter(j) for j, column in enumerate(frame.columns, start=1)}
    for i, (_, row) in enumerate(frame.iterrows(), start=2):
        for column, letter in letters.items():
            value = row[column]
            cell = sheet.cell(row=i, column=letters[column] and frame.columns.get_loc(column) + 1)
            if pd.isna(value):
                value = None
            cell.value = value
            if column.endswith(("_source", "_Edit_Flag", "rule")):
                cell.fill = SOURCE_FILL
                cell.font = Font(size=9, color="44506B")
                cell.alignment = Alignment(wrap_text=True, vertical="top")
        for column, template in (formulas or {}).items():
            sheet.cell(row=i, column=frame.columns.get_loc(column) + 1).value = template.format(
                row=i, **{re.sub(r"\W", "_", k): v for k, v in letters.items()})
    sheet.freeze_panes = "C2"
    sheet.auto_filter.ref = f"A1:{get_column_letter(len(frame.columns))}{len(frame) + 1}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, help="the workbook to extend (it is not modified in place)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--data", default="h2/data")
    args = ap.parse_args()
    data = Path(args.data)

    features = pd.read_csv(data / "h2_features.csv")
    # the two spreadsheet formulas are rewritten against this sheet's own columns
    features = features.drop(columns=[c for c in ("Therapeutic Index", "In Vivo Free Drug Concentration")
                                      if c in features.columns])
    features["Therapeutic Index"] = ""
    features["In Vivo Free Drug Concentration"] = ""
    formulas = {
        "Therapeutic Index":
            '=IF(OR({In_Vivo_Toxicity__rat_oral_LD50__mg_kg_}{row}="",'
            '{In_Vivo_Activity__1_effective_dose__kg_mg_}{row}=""),"",'
            'LOG({In_Vivo_Toxicity__rat_oral_LD50__mg_kg_}{row})-'
            'LOG(1/{In_Vivo_Activity__1_effective_dose__kg_mg_}{row}))',
        "In Vivo Free Drug Concentration":
            '=IF(COUNT({Bioavailability_Fraction_Absorbed_}{row},{Half_Life__h_}{row},'
            '{Dosage__mg_}{row},{Dose_Interval__h_}{row})<4,"",'
            '1.44*{Bioavailability_Fraction_Absorbed_}{row}*{Half_Life__h_}{row}*'
            '{Dosage__mg_}{row}/{Dose_Interval__h_}{row})',
    }

    book = load_workbook(args.source)
    write_sheet(book, "H2 per route", features, formulas)
    write_sheet(book, "H2 QED terms", pd.read_csv(data / "qed_table.csv"))
    write_sheet(book, "H2 docking", pd.read_csv(data / "docking.csv"))

    # the first markdown table in validation.md is the checks; the second one is the docking table, already a sheet
    block: list[str] = []
    for line in (data / "validation.md").read_text().splitlines():
        if line.startswith("|"):
            block.append(line)
        elif block:
            break
    rows = [re.split(r"\s*\|\s*", line.strip("| ")) for line in block]
    table = [r for r in rows if r and not set("".join(r)) <= set("-: ")]
    if table:
        write_sheet(book, "H2 validation", pd.DataFrame(table[1:], columns=table[0]))

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    book.save(args.out)
    print(f"wrote {args.out} with sheets: {book.sheetnames}")


if __name__ == "__main__":
    main()
