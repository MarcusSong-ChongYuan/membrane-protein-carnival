#!/usr/bin/env python3
"""Wait for evidence construction, then build, validate and freeze V6."""

from __future__ import annotations

import datetime as dt
import json
import pathlib
import subprocess
import sys
import time
import traceback


ROOT = pathlib.Path(__file__).resolve().parents[1]
RUN = ROOT / "runs" / "incremental_20260727"
QA = RUN / "qa"
LOG = RUN / "logs" / "finalize_v6_pipeline.log"
STATUS = QA / "PIPELINE_V6_STATUS.json"
PYTHON = pathlib.Path(sys.executable)


def write_status(state: str, detail: str = "") -> None:
    payload = {
        "state": state,
        "detail": detail,
        "updated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    temporary = STATUS.with_suffix(".json.tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    temporary.replace(STATUS)


def run_script(name: str) -> None:
    script = ROOT / "scripts" / name
    write_status("running", name)
    with LOG.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(
            f"\n[{dt.datetime.now(dt.timezone.utc).isoformat()}] "
            f"START {name}\n"
        )
        handle.flush()
        result = subprocess.run(
            [str(PYTHON), str(script)],
            cwd=str(ROOT),
            stdout=handle,
            stderr=handle,
            text=True,
            check=False,
        )
        handle.write(
            f"[{dt.datetime.now(dt.timezone.utc).isoformat()}] "
            f"END {name} exit={result.returncode}\n"
        )
    if result.returncode:
        raise RuntimeError(f"{name} failed with exit code {result.returncode}")


def main() -> int:
    try:
        write_status("waiting", "INCREMENTAL_EVIDENCE_BUILD_V6_QA.json")
        deadline = time.time() + 12 * 60 * 60
        evidence_qa = QA / "INCREMENTAL_EVIDENCE_BUILD_V6_QA.json"
        while time.time() < deadline:
            if evidence_qa.exists():
                payload = json.loads(evidence_qa.read_text(encoding="utf-8"))
                if payload.get("status") == "passed":
                    break
                if payload.get("status") == "failed":
                    raise RuntimeError("evidence build QA failed")
            time.sleep(30)
        else:
            raise TimeoutError("timed out waiting for evidence build QA")

        run_script("build_incremental_master_release_v6.py")
        run_script("validate_and_freeze_release_v6.py")
        write_status(
            "completed",
            str(ROOT / "releases" / "release_mempro_v6_0_20260727"),
        )
        return 0
    except Exception:
        detail = traceback.format_exc()
        write_status("failed", detail)
        with LOG.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(detail + "\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
