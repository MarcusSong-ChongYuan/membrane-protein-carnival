from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "membrane_master_v5_working"
V4 = ROOT / "v4_working" / "releases" / "release_v4"
DIRECTORIES = [
    "raw",
    "normalized",
    "crosswalks",
    "review",
    "reports",
    "release",
    "logs",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    for name in DIRECTORIES:
        (WORK / name).mkdir(parents=True, exist_ok=True)
    files = []
    for path in sorted(V4.iterdir()):
        if path.is_file() and not path.name.startswith("~$"):
            files.append(
                {
                    "path": path.relative_to(ROOT).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    manifest = {
        "manifest_version": "membrane-master-v5-baseline-v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "baseline_release": "v4_working/releases/release_v4",
        "baseline_is_immutable": True,
        "file_count": len(files),
        "files": files,
    }
    (WORK / "baseline_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    (WORK / "baseline_checksums.sha256").write_text(
        "\n".join(f"{item['sha256']}  {item['path']}" for item in files) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": "passed", "baseline_file_count": len(files)}, indent=2))


if __name__ == "__main__":
    main()
