#!/usr/bin/env python3
"""Download and validate the selected HPA v25.1 physiological matrices."""

from __future__ import annotations

import concurrent.futures
import datetime as dt
import hashlib
import json
import os
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path


RUN = Path(
    r"D:\7.22\evidence_expansion_v2_working\runs\v62_completion_20260729"
)
RAW = RUN / "raw" / "hpa_v25_1"
QA = RUN / "qa"
PROGRESS = QA / "HPA_V25_1_DOWNLOAD_PROGRESS.json"
BASE_URL = "https://www.proteinatlas.org/download/tsv"
FILES = {
    "rna_tissue_consensus.tsv.zip": 5_293_680,
    "rna_tissue_hpa.tsv.zip": 6_653_963,
    "rna_single_cell_type.tsv.zip": 16_346_850,
    "rna_single_cell_cluster.tsv.zip": 196_793_490,
    "normal_ihc_data.tsv.zip": 5_732_831,
    "ms_tissue.tsv.zip": 2_215_153,
    "subcellular_location.tsv.zip": 252_331,
}


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    os.replace(tmp, path)


def validate_zip(path: Path) -> tuple[bool, list[str]]:
    try:
        with zipfile.ZipFile(path) as archive:
            bad = archive.testzip()
            names = archive.namelist()
        return bad is None and bool(names), names
    except (OSError, zipfile.BadZipFile):
        return False, []


def download_one(name: str) -> dict:
    destination = RAW / name
    expected = FILES[name]
    if destination.exists():
        valid, members = validate_zip(destination)
        if valid:
            return {
                "status": "complete_cached",
                "bytes": destination.stat().st_size,
                "expected_bytes": expected,
                "sha256": sha256(destination),
                "members": members,
                "path": str(destination),
            }
        destination.unlink()

    part = destination.with_suffix(destination.suffix + ".part")
    url = f"{BASE_URL}/{name}"
    last_error = ""
    for attempt in range(1, 7):
        try:
            offset = part.stat().st_size if part.exists() else 0
            headers = {
                "User-Agent": "MemPro-V6.2-HPA-import/1.0",
                "Accept": "application/zip",
            }
            if offset:
                headers["Range"] = f"bytes={offset}-"
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=120) as response:
                append = offset > 0 and response.status == 206
                if offset and not append:
                    offset = 0
                mode = "ab" if append else "wb"
                with part.open(mode) as output:
                    while True:
                        block = response.read(4 * 1024 * 1024)
                        if not block:
                            break
                        output.write(block)
            os.replace(part, destination)
            valid, members = validate_zip(destination)
            if not valid:
                raise zipfile.BadZipFile("downloaded archive failed CRC validation")
            return {
                "status": "complete",
                "bytes": destination.stat().st_size,
                "expected_bytes": expected,
                "sha256": sha256(destination),
                "members": members,
                "path": str(destination),
                "url": url,
            }
        except Exception as exc:  # network retry is intentional
            last_error = f"{type(exc).__name__}: {exc}"
            time.sleep(min(30, 2**attempt))
    return {
        "status": "failed",
        "expected_bytes": expected,
        "error": last_error,
        "path": str(destination),
        "url": url,
    }


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    started = now()
    state = {
        "status": "running",
        "started_utc": started,
        "updated_utc": now(),
        "hpa_version": "25.1",
        "ensembl_version": "109",
        "files": {},
    }
    write_json(PROGRESS, state)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        futures = {
            executor.submit(download_one, name): name for name in FILES
        }
        for future in concurrent.futures.as_completed(futures):
            name = futures[future]
            state["files"][name] = future.result()
            state["updated_utc"] = now()
            write_json(PROGRESS, state)

    failed = [
        name
        for name, result in state["files"].items()
        if result["status"] == "failed"
    ]
    state["status"] = "complete" if not failed else "incomplete"
    state["failed_files"] = failed
    state["completed_utc"] = now()
    state["total_downloaded_bytes"] = sum(
        result.get("bytes", 0) for result in state["files"].values()
    )
    write_json(PROGRESS, state)
    print(json.dumps(state, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
