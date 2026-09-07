#!/usr/bin/env python3
"""Run the V6.3 classifier with NumPy scalar JSON compatibility."""

from __future__ import annotations

import json
import runpy
import sys
from pathlib import Path

import numpy as np


_default = json.JSONEncoder.default


def _numpy_default(self, value):
    if isinstance(value, np.generic):
        return value.item()
    return _default(self, value)


json.JSONEncoder.default = _numpy_default
script = Path(__file__).with_name("build_protein_cross_classification_v63.py")
sys.argv[0] = str(script)
runpy.run_path(str(script), run_name="__main__")
