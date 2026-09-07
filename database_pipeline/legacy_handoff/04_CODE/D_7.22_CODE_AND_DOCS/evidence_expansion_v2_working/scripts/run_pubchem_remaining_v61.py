#!/usr/bin/env python3
"""Wait for P1, then fetch P2 and P3 with independent checkpoints."""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RUN = ROOT / "runs" / "incremental_v61_20260727"
SCRIPT = ROOT / "scripts" / "fetch_pubchem_compound_properties_v61.py"
MAP_SCRIPT = ROOT / "scripts" / "map_pubchem_properties_v61.py"
BRENDA_SCRIPT = ROOT / "scripts" / "resolve_brenda_names_v61.py"
P1_PROGRESS = RUN / "qa" / "PUBCHEM_IDENTITY_FETCH_V61_PROGRESS.json"


def wait_for_p1() -> None:
    while True:
        if P1_PROGRESS.exists():
            try:
                state = json.loads(P1_PROGRESS.read_text(encoding="utf-8"))
            except (json.JSONDecodeError, OSError):
                state = {}
            if state.get("status") == "complete":
                return
            if state.get("status") == "incomplete":
                raise RuntimeError(f"P1 did not complete: {state}")
        time.sleep(30)


def run(priority: str, tag: str) -> None:
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--priority",
            priority,
            "--job-tag",
            tag,
            "--batch-size",
            "100",
            "--requests-per-second",
            "3",
        ],
        check=True,
    )
    subprocess.run(
        [
            sys.executable,
            str(MAP_SCRIPT),
            "--job-tag",
            tag,
            "--fetch-progress-name",
            f"PUBCHEM_IDENTITY_FETCH_{tag.upper()}_V61_PROGRESS.json",
            "--properties-name",
            f"pubchem_properties_{tag}_v61.tsv.gz",
            "--wait-hours",
            "1",
        ],
        check=True,
    )


if __name__ == "__main__":
    wait_for_p1()
    run("P2_ACTIVE_OTHER", "p2")
    run("P3_NONACTIVE_OR_UNSPECIFIED", "p3")
    subprocess.run([sys.executable, str(BRENDA_SCRIPT)], check=True)
