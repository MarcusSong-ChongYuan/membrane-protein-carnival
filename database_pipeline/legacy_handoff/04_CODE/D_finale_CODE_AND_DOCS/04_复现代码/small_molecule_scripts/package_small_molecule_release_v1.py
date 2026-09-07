from __future__ import annotations

import csv
import hashlib
import json
import shutil
import zipfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(r"D:\7.22")
RELEASE = (
    ROOT
    / "small_molecule_v1_working"
    / "releases"
    / "release_small_molecule_v1_0"
)
OUTPUTS = ROOT / "outputs" / "small_molecule_v1_0_20260727"
SEMIFINAL = ROOT / "semifinal"
WORKBOOK_NAME = "human_membrane_small_molecule_master_v1_0.xlsx"
ZIP_NAME = "human_membrane_small_molecule_release_v1_0.zip"
MANIFEST_NAME = "RELEASE_MANIFEST_small_molecule_v1_0.tsv"
CHECKSUM_NAME = "SHA256SUMS_small_molecule_v1_0.txt"
INFO_NAME = "RELEASE_INFO_small_molecule_v1_0.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def role_for(name: str) -> str:
    roles = {
        WORKBOOK_NAME: "curated Excel review workbook",
        "small_molecule_master_v1_0.tsv": "canonical parent-compound master",
        "compound_form_hierarchy_v1_0.tsv": "exact form and parent hierarchy",
        "compound_source_xref_v1_0.tsv": "source-record mapping and provenance",
        "compound_synonyms_v1_0.tsv": "compound names and synonyms",
        "compound_unresolved_candidates_v1_0.tsv": "unresolved identity audit",
        "compound_exclusion_audit_v1_0.tsv": "excluded non-small-molecule audit",
        "compound_merge_audit_v1_0.tsv": "identity and merge conflict audit",
        "compound_connectivity_families_v1_0.tsv": "connectivity-only related-structure index",
        "binding_evidence_master_v4_2.tsv": "full binding evidence with canonical compound IDs",
        "binding_site_instances_v4_2.tsv": "binding sites with canonical compound IDs",
        "protein_compound_summary_v4_2.tsv": "canonical protein-compound pairs",
        "compound_protein_summary_v1_0.tsv": "compound-level membrane-protein summary",
        "protein_binding_summary_v4_2.tsv": "protein-level canonical compound summary",
        "human_membrane_audit_master_v5_5.tsv": "protein master with canonical compound counts",
        "SMALL_MOLECULE_V1_VALIDATION_REPORT.json": "machine-readable validation report",
    }
    if name.endswith(".md"):
        return "release documentation"
    if name.startswith("SOURCE_REGISTRY"):
        return "source versions and URLs"
    if "excel_review" in name:
        return "compact source view used in the Excel review workbook"
    return roles.get(name, "release payload")


def main() -> None:
    SEMIFINAL.mkdir(parents=True, exist_ok=True)
    workbook_source = OUTPUTS / WORKBOOK_NAME
    if not workbook_source.exists():
        raise FileNotFoundError(workbook_source)
    shutil.copy2(workbook_source, RELEASE / WORKBOOK_NAME)
    shutil.copy2(workbook_source, SEMIFINAL / WORKBOOK_NAME)

    excluded_names = {MANIFEST_NAME, CHECKSUM_NAME, INFO_NAME}
    payload = sorted(
        path
        for path in RELEASE.iterdir()
        if path.is_file() and path.name not in excluded_names
    )
    entries = []
    for path in payload:
        entries.append(
            {
                "file_name": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "role": role_for(path.name),
            }
        )

    manifest = RELEASE / MANIFEST_NAME
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["file_name", "bytes", "sha256", "role"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(entries)

    checksum = RELEASE / CHECKSUM_NAME
    with checksum.open("w", encoding="utf-8", newline="\n") as handle:
        for entry in entries:
            handle.write(f"{entry['sha256']}  {entry['file_name']}\n")

    validation = json.loads(
        (RELEASE / "SMALL_MOLECULE_V1_VALIDATION_REPORT.json").read_text(
            encoding="utf-8"
        )
    )
    info = RELEASE / INFO_NAME
    info.write_text(
        json.dumps(
            {
                "release": "human-membrane-small-molecule-v1.0",
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "payload_file_count": len(entries),
                "validation_status": validation.get("status"),
                "canonical_parent_compounds": validation.get("parent_master_rows"),
                "canonical_forms": validation.get("form_rows"),
                "binding_evidence_rows": validation.get(
                    "binding_evidence_rows_v4_2"
                ),
                "protein_compound_pairs": validation.get(
                    "protein_compound_pair_rows_v4_2"
                ),
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    zip_path = SEMIFINAL / ZIP_NAME
    with zipfile.ZipFile(
        zip_path,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=6,
        allowZip64=True,
    ) as archive:
        for path in payload:
            archive.write(
                path,
                arcname=f"human_membrane_small_molecule_v1_0/{path.name}",
            )
        for path in (manifest, checksum, info):
            archive.write(
                path,
                arcname=f"human_membrane_small_molecule_v1_0/{path.name}",
            )

    print(
        json.dumps(
            {
                "release_dir": str(RELEASE),
                "workbook": str(SEMIFINAL / WORKBOOK_NAME),
                "zip_path": str(zip_path),
                "zip_bytes": zip_path.stat().st_size,
                "payload_file_count": len(entries),
                "validation_status": validation.get("status"),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
