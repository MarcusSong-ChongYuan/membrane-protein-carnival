import pandas as pd

import run_pagtn_sample_experiment as experiment


def structural_class_nan_safe(row: pd.Series) -> str:
    def integer(name: str) -> int:
        value = row.get(name, 0)
        return 0 if pd.isna(value) else int(value)

    formula_value = row.get("molecular_formula", "")
    formula = "" if pd.isna(formula_value) else str(formula_value)
    rings = integer("ring_count")
    max_ring = integer("max_ring_size")
    heavy = integer("heavy_atom_count")
    amide = integer("amide_bond_count")
    if "C" not in formula:
        return "inorganic/no carbon"
    if heavy >= 25 and amide >= 4:
        return "peptide-like"
    if max_ring >= 12:
        return "macrocycle"
    if rings >= 4:
        return "polycyclic (>=4 rings)"
    if rings >= 2:
        return "polycyclic (2-3 rings)"
    if rings == 1:
        return "monocyclic"
    return "acyclic organic"


experiment.structural_class = structural_class_nan_safe

if __name__ == "__main__":
    experiment.main()
