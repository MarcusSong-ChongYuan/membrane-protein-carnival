"""Generate only the remaining V8--V10 and optional Butina outputs.

This intentionally imports the established plotting implementation without
rerunning or overwriting V1--V7.
"""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "02_generate_visuals.py"
spec = importlib.util.spec_from_file_location("mempro_compound_visuals", SOURCE)
module = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(module)

analysis = module.pd.read_csv(module.ANALYSIS, sep="\t", low_memory=False)
module.plot_v8(analysis)
module.plot_v9(analysis)
module.plot_v10(analysis)
module.butina_benchmark(analysis)
print({"status": "PASS", "analysis_compounds": len(analysis), "generated": ["V8", "V9", "V10", "Butina"]})
