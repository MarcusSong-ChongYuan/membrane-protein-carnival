#!/usr/bin/env python3
"""Create a relative-path SHA-256 manifest for the complete module package."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    files = []
    for path in sorted(p for p in args.root.rglob("*") if p.is_file() and p != args.output):
        relative = path.relative_to(args.root).as_posix()
        if "/__pycache__/" in f"/{relative}/" or relative.endswith(".pyc"):
            continue
        files.append({
            "relative_path": relative,
            "bytes": path.stat().st_size,
            "sha256": digest(path),
            "package_section": relative.split("/", 1)[0],
        })
    manifest = {
        "package": "MemPro V6.3 candidate identity layer",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "file_count": len(files),
        "total_bytes": sum(item["bytes"] for item in files),
        "final_validation": "qa/IDENTITY_LAYER_V0_3_VALIDATION.json",
        "files": files,
    }
    args.output.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"file_count": manifest["file_count"], "total_bytes": manifest["total_bytes"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
