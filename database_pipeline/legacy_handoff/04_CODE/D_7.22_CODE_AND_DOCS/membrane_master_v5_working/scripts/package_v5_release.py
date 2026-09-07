from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(r"C:\Users\Administrator\7.22\membrane_master_v5_working")
PACKAGE = ROOT / "releases" / "release_v5"
WORKBOOK = Path(
    r"C:\Users\Administrator\outputs\membrane-v5-20260724"
    r"\human_membrane_protein_master_v5.xlsx"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_file(source: Path, relative_destination: Path) -> None:
    destination = PACKAGE / relative_destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def main() -> None:
    if PACKAGE.exists():
        raise FileExistsError(
            f"Release directory already exists and was not overwritten: {PACKAGE}"
        )

    PACKAGE.mkdir(parents=True)

    for source in sorted((ROOT / "release").glob("*")):
        if source.is_file():
            copy_file(source, Path(source.name))

    copy_file(WORKBOOK, Path(WORKBOOK.name))

    for folder in ("crosswalks", "review", "normalized", "reports"):
        for source in sorted((ROOT / folder).glob("*")):
            if source.is_file():
                copy_file(source, Path(folder) / source.name)

    for source in sorted((ROOT / "scripts").glob("*")):
        if source.is_file() and source.suffix.lower() in {".py", ".ps1"}:
            copy_file(source, Path("reproducibility") / "scripts" / source.name)

    copy_file(
        ROOT / "spreadsheet_build" / "build_v5_workbook.mjs",
        Path("reproducibility") / "scripts" / "build_v5_workbook.mjs",
    )
    copy_file(
        ROOT / "configs" / "release_tiers.yaml",
        Path("reproducibility") / "configs" / "release_tiers.yaml",
    )
    copy_file(
        ROOT / "raw" / "source_download_manifest.json",
        Path("reproducibility") / "source_download_manifest.json",
    )
    copy_file(
        ROOT / "baseline_manifest.json",
        Path("reproducibility") / "baseline_manifest.json",
    )
    copy_file(
        ROOT / "baseline_checksums.sha256",
        Path("reproducibility") / "baseline_checksums.sha256",
    )

    package_readme = """# Human Membrane Protein Master v5 release package

Release date: 2026-07-24

## Entry points

- `human_membrane_protein_master_v5.xlsx`: human-readable review workbook.
- `human_membrane_protein_master_v5.tsv`: full machine-readable master table.
- `human_integral_membrane_view_v5.tsv`: Tier A.
- `human_core_membrane_view_v5.tsv`: Tiers A+B.
- `human_extended_membrane_view_v5.tsv`: Tiers A+B+C.
- `review/human_membrane_review_queue_v5.tsv`: Tier R records.

## Supporting folders

- `crosswalks/`: TrEMBL canonicalization, legacy-record, topology,
  source-evidence, and classification audit tables.
- `normalized/`: normalized source-specific evidence tables used by the build.
- `review/`: manual-review and source-conflict queues.
- `reports/`: machine-readable build statistics and validation results.
- `reproducibility/`: scripts, tier policy, source-download manifest, and
  baseline hashes.

The raw third-party downloads are not redistributed in this package. Their
official URLs, release identifiers, roles, and license notes are recorded in
`SOURCE_REGISTRY_v5.tsv` and the reproducibility manifest.
"""
    (PACKAGE / "README_PACKAGE.md").write_text(package_readme, encoding="utf-8")

    files = [path for path in PACKAGE.rglob("*") if path.is_file()]
    records = []
    for path in sorted(files):
        relative = path.relative_to(PACKAGE).as_posix()
        records.append(
            {
                "path": relative,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )

    manifest = {
        "release": "v5",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "file_count_excluding_manifests": len(records),
        "files": records,
    }
    manifest_path = PACKAGE / "RELEASE_MANIFEST_v5.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    checksum_lines = [
        f"{record['sha256']}  {record['path']}" for record in records
    ]
    checksum_lines.append(f"{sha256(manifest_path)}  RELEASE_MANIFEST_v5.json")
    (PACKAGE / "SHA256SUMS_v5.txt").write_text(
        "\n".join(checksum_lines) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "package": str(PACKAGE),
                "file_count": len(records) + 1,
                "bytes": sum(path.stat().st_size for path in PACKAGE.rglob("*") if path.is_file()),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
