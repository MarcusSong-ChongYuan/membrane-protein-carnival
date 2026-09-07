from __future__ import annotations

import csv
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path


ROOT = Path(r"C:\Users\Administrator\7.22\membrane_master_v5_working")
V52 = ROOT / "releases" / "release_v5_2_evidence_tiered"
RELEASE = ROOT / "releases" / "release_v5_3_sequences"
WORKING = ROOT / "v53_sequence_working"
AUDIT_V52 = V52 / "human_membrane_audit_master_v5_2.tsv"
LOOKUP = WORKING / "uniprot_sequence_lookup_v53.tsv"


SEQUENCE_FIELDS = [
    "canonical_sequence",
    "sequence_length_v53",
    "sequence_length_match_v53",
    "sequence_version_v53",
    "sequence_sha256_v53",
    "sequence_status_v53",
    "sequence_resolved_accession_v53",
    "isoform_scope_v53",
    "sequence_source_v53",
    "sequence_source_url_v53",
    "uniprot_release_v53",
    "uniprot_release_date_v53",
    "sequence_retrieval_date_v53",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def wrap_sequence(sequence: str, width: int = 60) -> str:
    return "\n".join(sequence[start : start + width] for start in range(0, len(sequence), width))


def write_fasta(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for row in rows:
            prefix = "sp" if row["reviewed"].lower() == "reviewed" else "tr"
            header = (
                f">{prefix}|{row['target_uniprot_id']}|{row['uniprot_entry_name']} "
                f"HMP={row['membrane_protein_id']} "
                f"class={row['membrane_class_v52']} "
                f"evidence={row['evidence_level_v52']} "
                f"isoform=canonical SV={row['sequence_version_v53']}"
            )
            handle.write(header + "\n")
            handle.write(wrap_sequence(row["canonical_sequence"]) + "\n")


def main() -> None:
    RELEASE.mkdir(parents=True, exist_ok=True)
    (RELEASE / "reports").mkdir(exist_ok=True)
    (RELEASE / "reproducibility" / "scripts").mkdir(parents=True, exist_ok=True)

    source_rows = read_tsv(AUDIT_V52)
    sequence_rows = {
        row["requested_accession"]: row for row in read_tsv(LOOKUP)
    }
    fields = list(source_rows[0]) + SEQUENCE_FIELDS
    merged: list[dict[str, str]] = []
    for row in source_rows:
        accession = row["target_uniprot_id"]
        sequence = sequence_rows[accession]
        current_length = sequence["sequence_length_api"]
        item = dict(row)
        item.update(
            {
                "canonical_sequence": sequence["canonical_sequence"],
                "sequence_length_v53": current_length,
                "sequence_length_match_v53": "1"
                if row["sequence_length"] == current_length
                else "0",
                "sequence_version_v53": sequence["sequence_version"],
                "sequence_sha256_v53": sequence["sequence_sha256"],
                "sequence_status_v53": sequence["sequence_status"],
                "sequence_resolved_accession_v53": sequence["resolved_accession"],
                "isoform_scope_v53": sequence["isoform_scope"],
                "sequence_source_v53": sequence["sequence_source"],
                "sequence_source_url_v53": sequence["sequence_source_url"],
                "uniprot_release_v53": sequence["uniprot_release"],
                "uniprot_release_date_v53": sequence["uniprot_release_date"],
                "sequence_retrieval_date_v53": sequence[
                    "sequence_retrieval_date"
                ],
            }
        )
        merged.append(item)

    subsets = {
        "human_membrane_audit_master_v5_3.tsv": merged,
        "human_membrane_experimental_core_E1_v5_3.tsv": [
            row for row in merged if row["evidence_level_v52"] == "E1"
        ],
        "human_membrane_strong_supported_E2_v5_3.tsv": [
            row for row in merged if row["evidence_level_v52"] == "E2"
        ],
        "human_membrane_default_E1_E2_v5_3.tsv": [
            row for row in merged if row["evidence_level_v52"] in {"E1", "E2"}
        ],
        "human_membrane_prediction_candidates_E3_v5_3.tsv": [
            row for row in merged if row["evidence_level_v52"] == "E3"
        ],
        "human_membrane_excluded_E0_v5_3.tsv": [
            row for row in merged if row["evidence_level_v52"] == "E0"
        ],
        "human_membrane_class_A_E1_E2_v5_3.tsv": [
            row
            for row in merged
            if row["evidence_level_v52"] in {"E1", "E2"}
            and row["membrane_class_v52"] == "A"
        ],
        "human_membrane_class_B_E1_E2_v5_3.tsv": [
            row
            for row in merged
            if row["evidence_level_v52"] in {"E1", "E2"}
            and row["membrane_class_v52"] == "B"
        ],
        "human_membrane_class_C_E1_E2_v5_3.tsv": [
            row
            for row in merged
            if row["evidence_level_v52"] in {"E1", "E2"}
            and row["membrane_class_v52"] == "C"
        ],
    }
    for name, rows in subsets.items():
        write_tsv(RELEASE / name, rows, fields)

    default_rows = subsets["human_membrane_default_E1_E2_v5_3.tsv"]
    write_fasta(RELEASE / "human_membrane_all_audit_v5_3.fasta", merged)
    write_fasta(RELEASE / "human_membrane_default_E1_E2_v5_3.fasta", default_rows)

    long_rows = [
        {
            "target_uniprot_id": row["target_uniprot_id"],
            "uniprot_entry_name": row["uniprot_entry_name"],
            "sequence_length_v53": row["sequence_length_v53"],
            "excel_storage_action": "split_across_sequence_part_1_and_part_2",
            "tsv_fasta_status": "full_sequence_preserved",
        }
        for row in merged
        if len(row["canonical_sequence"]) > 32_767
    ]
    write_tsv(
        RELEASE / "sequence_excel_storage_exceptions_v5_3.tsv",
        long_rows,
        [
            "target_uniprot_id",
            "uniprot_entry_name",
            "sequence_length_v53",
            "excel_storage_action",
            "tsv_fasta_status",
        ],
    )

    shutil.copy2(V52 / "EVIDENCE_POLICY_v5_2.md", RELEASE / "EVIDENCE_POLICY_v5_2.md")
    shutil.copy2(
        ROOT / "scripts" / "fetch_uniprot_sequences_v53.py",
        RELEASE / "reproducibility" / "scripts" / "fetch_uniprot_sequences_v53.py",
    )
    shutil.copy2(
        ROOT / "scripts" / "build_v53_sequence_release.py",
        RELEASE / "reproducibility" / "scripts" / "build_v53_sequence_release.py",
    )
    shutil.copy2(
        ROOT / "scripts" / "build_v52_evidence_tiers.py",
        RELEASE / "reproducibility" / "scripts" / "build_v52_evidence_tiers.py",
    )

    readme = """# Human Membrane Protein Master v5.3 — sequence enhanced

This release preserves the v5.2 ABC × E0–E3 biological/evidence model and adds
the current UniProtKB canonical amino-acid sequence to every audited record.

## Coverage

- All audited records: 10,997
- E1 experimental core: 3,451
- E2 strongly supported: 4,453
- Default website release (E1+E2): 7,904
- E3 prediction candidates: 2,652
- E0 excluded audit records: 441
- Canonical sequences retrieved: 10,997 (100%)
- Missing sequences: 0
- Sequence-length mismatches against the previous table: 0

## Sequence files

- `human_membrane_all_audit_v5_3.fasta` contains all 10,997 canonical sequences.
- `human_membrane_default_E1_E2_v5_3.fasta` contains the 7,904 default website records.
- Every TSV includes one full `canonical_sequence` field per record.
- Excel stores sequences in `Canonical sequence part 1` and `part 2`. Only Q8WZ42
  requires part 2 because its 34,350-aa sequence exceeds Excel's per-cell limit.

The v5.2 release is not modified.
"""
    (RELEASE / "README_v5_3.md").write_text(readme, encoding="utf-8")

    policy = """# Sequence policy v5.3

## Scope

The canonical sequence for each `target_uniprot_id` was retrieved from the
official UniProtKB REST API on 2026-07-24. Isoform-specific sequences are not
included in this release.

## Provenance

- UniProtKB release: 2026_02
- UniProtKB release date: 10-June-2026
- Retrieval endpoint: https://rest.uniprot.org/uniprotkb/search
- Returned fields: accession, entry name, reviewed status, sequence, length,
  and sequence version

## Integrity

`sequence_sha256_v53` is a SHA-256 digest calculated locally from the unwrapped
uppercase amino-acid sequence. `sequence_length_match_v53` compares the
retrieved sequence length with the length already stored in v5.2.

## Excel limitation

Excel permits at most 32,767 characters in one cell. The workbook therefore
splits sequences after 32,000 residues. Concatenating the two sequence-part
columns reconstructs the exact canonical sequence. TSV and FASTA outputs
always contain the unsplit full sequence.
"""
    (RELEASE / "SEQUENCE_POLICY_v5_3.md").write_text(policy, encoding="utf-8")

    dictionary = """# Sequence fields added in v5.3

| Field | Meaning |
|---|---|
| canonical_sequence | Complete unwrapped UniProtKB canonical sequence |
| sequence_length_v53 | Length returned by UniProtKB |
| sequence_length_match_v53 | 1 when the v5.3 length matches the prior table |
| sequence_version_v53 | UniProtKB sequence version |
| sequence_sha256_v53 | SHA-256 digest of the canonical sequence |
| sequence_status_v53 | Retrieval/mapping result |
| sequence_resolved_accession_v53 | Current accession returned by UniProtKB |
| isoform_scope_v53 | Sequence scope; canonical in this release |
| sequence_source_v53 | Source service |
| sequence_source_url_v53 | Direct UniProtKB FASTA endpoint |
| uniprot_release_v53 | UniProtKB data release |
| uniprot_release_date_v53 | UniProtKB release date |
| sequence_retrieval_date_v53 | Date the sequence was retrieved |
"""
    (RELEASE / "DATA_DICTIONARY_v5_3.md").write_text(dictionary, encoding="utf-8")

    registry = """source\tversion_or_date\trole\turl
UniProtKB REST\t2026_02; released 10-June-2026; retrieved 2026-07-24\tCanonical amino-acid sequence, length, sequence version, accession resolution\thttps://rest.uniprot.org/uniprotkb/search
Human Membrane Protein Master v5.2\tfrozen v5.2 evidence-tiered release\tProtein identity, ABC biological class, E0-E3 evidence level, membrane evidence and external cross-references\tlocal release_v5_2_evidence_tiered
"""
    (RELEASE / "SOURCE_REGISTRY_v5_3.tsv").write_text(registry, encoding="utf-8")

    counts = Counter(row["evidence_level_v52"] for row in merged)
    report = {
        "release": "v5.3-sequence-enhanced",
        "audit_rows": len(merged),
        "unique_accessions": len({row["target_uniprot_id"] for row in merged}),
        "sequence_rows": sum(bool(row["canonical_sequence"]) for row in merged),
        "missing_sequences": sum(not row["canonical_sequence"] for row in merged),
        "length_matches": sum(row["sequence_length_match_v53"] == "1" for row in merged),
        "length_mismatches": sum(row["sequence_length_match_v53"] == "0" for row in merged),
        "total_amino_acids": sum(int(row["sequence_length_v53"]) for row in merged),
        "longest_sequence_length": max(
            int(row["sequence_length_v53"]) for row in merged
        ),
        "excel_split_records": len(long_rows),
        "evidence_level_counts": dict(sorted(counts.items())),
        "default_E1_E2_rows": len(default_rows),
        "uniprot_release": "2026_02",
        "uniprot_release_date": "10-June-2026",
        "retrieval_date": "2026-07-24",
    }
    (RELEASE / "reports" / "V53_SEQUENCE_STATS.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
