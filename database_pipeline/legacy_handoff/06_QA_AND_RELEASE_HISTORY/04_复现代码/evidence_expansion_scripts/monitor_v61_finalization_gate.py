#!/usr/bin/env python3
"""Track when all identity inputs are ready for the formal V6.1 rebuild."""

from __future__ import annotations

import datetime as dt
import json
import time
from pathlib import Path


ROOT = Path(r"D:\7.22\evidence_expansion_v2_working")
RUN = ROOT / "runs" / "incremental_v61_20260727"
QA = RUN / "qa"
OUT = QA / "V61_FINALIZATION_GATE.json"
REQUIRED = {
    "PubChem_P1_map": "PUBCHEM_IDENTITY_MAP_P1_V61_PROGRESS.json",
    "PubChem_P2_map": "PUBCHEM_IDENTITY_MAP_P2_V61_PROGRESS.json",
    "PubChem_P3_map": "PUBCHEM_IDENTITY_MAP_P3_V61_PROGRESS.json",
    "BRENDA_identity": "BRENDA_IDENTITY_V61_PROGRESS.json",
    "PDBe_CCD": "PDBE_CCD_IDENTITY_V61_PROGRESS.json",
    "PDBbind_identity": "PDBBIND_IDENTITY_V61_PROGRESS.json",
}


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def write(payload: dict) -> None:
    temporary = OUT.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    temporary.replace(OUT)


def read_state(filename: str) -> dict:
    path = QA / filename
    if not path.exists():
        return {"status": "not_started", "path": str(path)}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        return {"status": "unreadable", "path": str(path), "error": str(exc)}
    return {
        "status": payload.get("status", payload.get("state", "unknown")),
        "updated_utc": payload.get(
            "completed_utc", payload.get("updated_utc", payload.get("updated_at_utc"))
        ),
        "path": str(path),
    }


if __name__ == "__main__":
    while True:
        states = {name: read_state(path) for name, path in REQUIRED.items()}
        ready = all(value["status"] == "complete" for value in states.values())
        write(
            {
                "status": (
                    "ready_for_final_refresh" if ready else "waiting_for_upstream"
                ),
                "updated_utc": now(),
                "required_inputs": states,
                "next_action": (
                    "run full identity refresh, evidence rebuild, validation and freeze"
                    if ready
                    else "continue checkpointed upstream jobs"
                ),
            }
        )
        if ready:
            break
        time.sleep(60)
