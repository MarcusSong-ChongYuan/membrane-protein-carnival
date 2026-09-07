from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "v4_working"
RELEASE = WORK / "releases" / "release_v4"
WORKBOOK = WORK / "outputs" / "v4_excel" / "membrane_protein_database_v4.xlsx"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    destination = RELEASE / WORKBOOK.name
    shutil.copy2(WORKBOOK, destination)
    manifest_path = RELEASE / "RELEASE_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["status"] = "validated"
    manifest["finalized_at_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["files"] = [
        {
            "filename": path.name,
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for path in sorted(RELEASE.iterdir())
        if path.is_file()
        and path.name not in {"RELEASE_MANIFEST.json", "checksums.sha256"}
    ]
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    checksum_path = RELEASE / "checksums.sha256"
    checksum_path.write_text(
        "\n".join(
            f"{sha256(path)}  {path.name}"
            for path in sorted(RELEASE.iterdir())
            if path.is_file() and path.name != "checksums.sha256"
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": "validated",
                "release": str(RELEASE),
                "file_count": len([p for p in RELEASE.iterdir() if p.is_file()]),
                "workbook_bytes": destination.stat().st_size,
                "workbook_sha256": sha256(destination),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
