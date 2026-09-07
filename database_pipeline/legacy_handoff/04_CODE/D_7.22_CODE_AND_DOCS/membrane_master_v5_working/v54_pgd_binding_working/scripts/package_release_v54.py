from __future__ import annotations

import csv
import hashlib
import json
import shutil
import zipfile
from datetime import datetime, timezone
from pathlib import Path


RELEASE_DIR = Path(
    r"D:\7.22\membrane_master_v5_working\releases\release_v5_4_pgd_binding"
)
OUTPUT_DIR = Path(r"D:\7.22\outputs\membrane_v54_pgd_binding_20260724")
SEMIFINAL_DIR = Path(r"D:\7.22\semifinal")

PROTEIN_WORKBOOK = "human_membrane_protein_master_v5_4_with_disease_binding.xlsx"
BINDING_WORKBOOK = "human_membrane_small_molecule_binding_evidence_v4_1.xlsx"
ZIP_NAME = "human_membrane_release_v5_4_pgd_binding.zip"

RELEASE_PAYLOAD = [
    PROTEIN_WORKBOOK,
    BINDING_WORKBOOK,
    "human_membrane_audit_master_v5_4.tsv",
    "protein_gene_disease_relations_v5_4.tsv",
    "protein_primary_disease_v5_4.tsv",
    "protein_disease_summary_v5_4.tsv",
    "binding_evidence_master_v4_1.tsv",
    "binding_evidence_excel_review_v4_1.tsv",
    "binding_site_instances_v4_1.tsv",
    "protein_compound_summary_v4_1.tsv",
    "protein_binding_summary_v4_1.tsv",
    "small_molecule_index_v4_1.tsv",
    "v3_pgd_relations_not_in_v5_3_audit.tsv",
    "v3_binding_evidence_not_in_v5_3_audit.tsv",
    "README_v5_4.md",
    "METHODS_v5_4.md",
    "INCLUSION_POLICY_v5_4.md",
    "DATA_DICTIONARY_v5_4.md",
    "SOURCE_REGISTRY_v5_4.tsv",
    "V54_VALIDATION_REPORT.json",
]

ROLES = {
    PROTEIN_WORKBOOK: "review workbook: protein master with disease and binding summaries",
    BINDING_WORKBOOK: "review workbook: representative binding evidence and sites",
    "human_membrane_audit_master_v5_4.tsv": "canonical protein master",
    "protein_gene_disease_relations_v5_4.tsv": "canonical protein-gene-disease relation table",
    "protein_primary_disease_v5_4.tsv": "one primary disease row per supported protein",
    "protein_disease_summary_v5_4.tsv": "protein-level disease summary for master merge",
    "binding_evidence_master_v4_1.tsv": "canonical full binding evidence table",
    "binding_evidence_excel_review_v4_1.tsv": "capped Excel review view",
    "binding_site_instances_v4_1.tsv": "compound-specific and curated binding-site instances",
    "protein_compound_summary_v4_1.tsv": "protein-compound pair summary",
    "protein_binding_summary_v4_1.tsv": "protein-level binding summary for master merge",
    "small_molecule_index_v4_1.tsv": "compound identifier and provenance index",
    "v3_pgd_relations_not_in_v5_3_audit.tsv": "legacy disease rows outside the v5.3 protein universe",
    "v3_binding_evidence_not_in_v5_3_audit.tsv": "legacy binding rows outside the v5.3 protein universe",
    "README_v5_4.md": "release overview",
    "METHODS_v5_4.md": "methods",
    "INCLUSION_POLICY_v5_4.md": "evidence inclusion policy",
    "DATA_DICTIONARY_v5_4.md": "data dictionary",
    "SOURCE_REGISTRY_v5_4.tsv": "source versions and URLs",
    "V54_VALIDATION_REPORT.json": "machine-readable release validation",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    SEMIFINAL_DIR.mkdir(parents=True, exist_ok=True)

    for workbook_name in (PROTEIN_WORKBOOK, BINDING_WORKBOOK):
        source = OUTPUT_DIR / workbook_name
        if not source.exists():
            raise FileNotFoundError(source)
        shutil.copy2(source, RELEASE_DIR / workbook_name)
        shutil.copy2(source, SEMIFINAL_DIR / workbook_name)

    entries = []
    for name in RELEASE_PAYLOAD:
        path = RELEASE_DIR / name
        if not path.exists():
            raise FileNotFoundError(path)
        entries.append(
            {
                "file_name": name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
                "role": ROLES[name],
            }
        )

    generated_at = datetime.now(timezone.utc).isoformat()
    manifest_path = RELEASE_DIR / "RELEASE_MANIFEST_v5_4.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["file_name", "bytes", "sha256", "role"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(entries)

    checksums_path = RELEASE_DIR / "SHA256SUMS_v5_4.txt"
    with checksums_path.open("w", encoding="utf-8", newline="\n") as handle:
        for entry in entries:
            handle.write(f"{entry['sha256']}  {entry['file_name']}\n")

    release_info_path = RELEASE_DIR / "RELEASE_INFO_v5_4.json"
    validation = json.loads(
        (RELEASE_DIR / "V54_VALIDATION_REPORT.json").read_text(encoding="utf-8")
    )
    release_info_path.write_text(
        json.dumps(
            {
                "release": "v5.4-pgd-binding",
                "generated_at_utc": generated_at,
                "payload_file_count": len(entries),
                "validation_status": validation.get("status"),
                "validation": validation,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    zip_path = SEMIFINAL_DIR / ZIP_NAME
    with zipfile.ZipFile(
        zip_path,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=6,
        allowZip64=True,
    ) as archive:
        for entry in entries:
            source = RELEASE_DIR / entry["file_name"]
            archive.write(
                source,
                arcname=f"human_membrane_release_v5_4/{entry['file_name']}",
            )
        for extra in (manifest_path, checksums_path, release_info_path):
            archive.write(
                extra,
                arcname=f"human_membrane_release_v5_4/{extra.name}",
            )

    result = {
        "release_dir": str(RELEASE_DIR),
        "semifinal_dir": str(SEMIFINAL_DIR),
        "zip_path": str(zip_path),
        "zip_bytes": zip_path.stat().st_size,
        "payload_file_count": len(entries),
        "validation_status": validation.get("status"),
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
