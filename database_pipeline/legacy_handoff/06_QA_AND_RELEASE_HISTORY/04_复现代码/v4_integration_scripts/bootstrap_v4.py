from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "v4_working"
BASELINE = ROOT / "upload" / "normalized_tables" / "releases" / "partner_handoff_v3_1" / "v3_1"
DIRECTORIES = [
    "raw",
    "intermediate",
    "releases",
    "reports",
    "logs",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    for relative in DIRECTORIES:
        (WORK / relative).mkdir(parents=True, exist_ok=True)

    files = []
    checksum_lines = []
    for path in sorted(BASELINE.iterdir()):
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT).as_posix()
        checksum = sha256(path)
        files.append(
            {
                "path": relative,
                "bytes": path.stat().st_size,
                "sha256": checksum,
            }
        )
        checksum_lines.append(f"{checksum}  {relative}")

    manifest = {
        "manifest_version": "v4-baseline-manifest-v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "baseline_release": "partner_handoff_v3_1/v3_1",
        "baseline_is_immutable": True,
        "file_count": len(files),
        "files": files,
    }
    (WORK / "baseline_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    (WORK / "checksums.sha256").write_text(
        "\n".join(checksum_lines) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": "passed", "file_count": len(files)}, indent=2))


if __name__ == "__main__":
    main()
