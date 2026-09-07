#!/usr/bin/env python3
"""Start BRENDA identity resolution after the PubChem P3 job finishes."""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RUN = ROOT / "runs" / "incremental_v61_20260727"
P3 = RUN / "qa" / "PUBCHEM_IDENTITY_FETCH_P3_V61_PROGRESS.json"
BRENDA = ROOT / "scripts" / "resolve_brenda_names_v61.py"


if __name__ == "__main__":
    while True:
        if P3.exists():
            try:
                state = json.loads(P3.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError):
                state = {}
            if state.get("status") == "complete":
                break
            if state.get("status") == "incomplete":
                raise RuntimeError(f"PubChem P3 did not complete: {state}")
        time.sleep(30)
    subprocess.run([sys.executable, str(BRENDA)], check=True)
