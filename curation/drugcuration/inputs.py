import csv
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_INPUT = os.path.join(HERE, "inputs", "atc_a_small_molecules.csv")


def load_drugs(path=None):
    with open(path or DEFAULT_INPUT, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))
