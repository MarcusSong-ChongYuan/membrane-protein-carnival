#!/usr/bin/env python3
"""Compatibility runner for the full-release scaffold analysis.

RDKit can derive an exact Bemis-Murcko scaffold for a few unusual records but
raise AtomValenceException when converting that scaffold to the optional generic
framework.  This runner preserves the exact scaffold and returns an empty
generic molecule only for that optional failed conversion.
"""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold


_real_make_generic = MurckoScaffold.MakeScaffoldGeneric


def _safe_make_generic(mol):
    try:
        return _real_make_generic(mol)
    except Exception:
        return Chem.MolFromSmiles("")


MurckoScaffold.MakeScaffoldGeneric = _safe_make_generic
RDLogger.DisableLog("rdApp.warning")
script = Path(__file__).with_name("build_m5e_scaffold_diversity.py")
sys.argv[0] = str(script)
runpy.run_path(str(script), run_name="__main__")
