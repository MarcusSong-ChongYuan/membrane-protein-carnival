from __future__ import annotations

import hashlib
import json
from pathlib import Path


BASE = Path(r"C:\github-repos\upload\normalized_tables")


def read_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines()


def write_lines(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def replace_file_summary_line(lines: list[str], file_name: str, rows: str, cols: str, size: str, meaning_tail: str) -> list[str]:
    prefix = f"| `{file_name}` |"
    out = []
    for line in lines:
        if line.startswith(prefix):
            out.append(f"| `{file_name}` | {rows} | {cols} | {size} | {meaning_tail} |")
        else:
            out.append(line)
    return out


def remove_evidence_levels_lines(lines: list[str]) -> list[str]:
    return [line for line in lines if "`evidence_levels`" not in line]


def replace_exact(lines: list[str], old: str, new: str) -> list[str]:
    return [new if line == old else line for line in lines]


def replace_contains(lines: list[str], needle: str, new: str) -> list[str]:
    return [new if needle in line else line for line in lines]


def renumber_markdown_rows(lines: list[str], mapping: dict[str, int]) -> list[str]:
    out = []
    for line in lines:
        replaced = False
        for col, num in mapping.items():
            marker = f"| `{col}` |"
            if marker in line:
                parts = line.split("|")
                if len(parts) > 2:
                    parts[1] = f" {num} "
                    out.append("|".join(parts))
                    replaced = True
                    break
        if not replaced:
            out.append(line)
    return out


def update_data_dictionary() -> None:
    path = BASE / "DATA_DICTIONARY.md"
    lines = read_lines(path)
    lines = replace_file_summary_line(
        lines,
        "protein_database.tsv",
        "1,773",
        "25",
        "2.8 MB",
        "One row per UniProt protein. Drug counts and source summaries are recomputed from the binary table.",
    )
    lines = replace_file_summary_line(
        lines,
        "small_molecule_database.tsv",
        "87,552",
        "32",
        "42.7 MB",
        "One row per unique compound/drug ID. Identity, activity, PubChem properties, safety/toxicity, clinical fields.",
    )
    lines = remove_evidence_levels_lines(lines)
    lines = replace_exact(lines, "### Relationship Summary (14-18)", "### Relationship Summary (14-17)")
    lines = replace_exact(lines, "### UniProt Enrichment (19-26)", "### UniProt Enrichment (18-25)")
    lines = replace_contains(lines, "### Source Evidence Summary", "### Source Evidence Summary (6-9)")
    lines = replace_contains(lines, "### Activity Summary", "### Activity Summary (10-12)")
    lines = replace_contains(lines, "### PubChem Chemical Properties", "### PubChem Chemical Properties (13-25)")
    lines = replace_contains(lines, "### Safety/Toxicity", "### Safety/Toxicity (26-32)")
    lines = renumber_markdown_rows(
        lines,
        {
            "has_matched_evidence": 17,
            "reviewed": 18,
            "go_cellular_component": 19,
            "go_molecular_function": 20,
            "pdb_structures": 21,
            "disease_association": 22,
            "transmembrane_count": 23,
            "subcellular_location": 24,
            "protein_function": 25,
            "clinical_statuses": 9,
            "activity_types": 10,
            "activity_value_uM_min": 11,
            "activity_value_uM_median": 12,
            "molecular_formula": 13,
            "molecular_weight": 14,
            "canonical_smiles": 15,
            "isomeric_smiles": 16,
            "inchikey": 17,
            "iupac_name": 18,
            "xlogp": 19,
            "tpsa": 20,
            "hbond_donor_count": 21,
            "hbond_acceptor_count": 22,
            "rotatable_bond_count": 23,
            "complexity": 24,
            "formal_charge": 25,
            "ghs_hazard_classes": 26,
            "ghs_signal_words": 27,
            "toxicity_summary": 28,
            "livertox": 29,
            "drug_classes": 30,
            "pharmacodynamics": 31,
            "drug_indication": 32,
        },
    )
    lines = replace_contains(
        lines,
        "| 17 | `has_matched_evidence` |",
        "| 17 | `has_matched_evidence` | `1` if any matched PDB/STITCH/structure-context evidence exists. | `0` or `1` |",
    )
    write_lines(path, lines)


def update_column_guide() -> None:
    path = BASE / "SHEET2_SHEET3_COLUMN_GUIDE.md"
    lines = read_lines(path)
    lines = replace_exact(lines, "Current shape: 1,773 rows x 26 columns. One row is one UniProt protein.", "Current shape: 1,773 rows x 25 columns. One row is one UniProt protein.")
    lines = replace_exact(lines, "Current shape: 87,552 rows x 33 columns. One row is one compound/drug identifier.", "Current shape: 87,552 rows x 32 columns. One row is one compound/drug identifier.")
    lines = remove_evidence_levels_lines(lines)
    lines = renumber_markdown_rows(
        lines,
        {
            "has_matched_evidence": 17,
            "reviewed": 18,
            "go_cellular_component": 19,
            "go_molecular_function": 20,
            "pdb_structures": 21,
            "disease_association": 22,
            "transmembrane_count": 23,
            "subcellular_location": 24,
            "protein_function": 25,
            "clinical_statuses": 9,
            "activity_types": 10,
            "activity_value_uM_min": 11,
            "activity_value_uM_median": 12,
            "molecular_formula": 13,
            "molecular_weight": 14,
            "canonical_smiles": 15,
            "isomeric_smiles": 16,
            "inchikey": 17,
            "iupac_name": 18,
            "xlogp": 19,
            "tpsa": 20,
            "hbond_donor_count": 21,
            "hbond_acceptor_count": 22,
            "rotatable_bond_count": 23,
            "complexity": 24,
            "formal_charge": 25,
            "ghs_hazard_classes": 26,
            "ghs_signal_words": 27,
            "toxicity_summary": 28,
            "livertox": 29,
            "drug_classes": 30,
            "pharmacodynamics": 31,
            "drug_indication": 32,
        },
    )
    write_lines(path, lines)


def update_process_log() -> None:
    path = BASE / "PROCESS_LOG.md"
    text = path.read_text(encoding="utf-8")
    text = text.replace("| `upload/normalized_tables/drug_protein_binary_relationships.tsv` | 95,234 | 29 | 50.0 MB |", "| `upload/normalized_tables/drug_protein_binary_relationships.tsv` | 101,171 | 29 | 53.0 MB |")
    text = text.replace("| `upload/normalized_tables/protein_database.tsv` | 1,773 | 26 | 5.4 MB |", "| `upload/normalized_tables/protein_database.tsv` | 1,773 | 25 | 2.8 MB |")
    text = text.replace("| `upload/normalized_tables/small_molecule_database.tsv` | 87,552 | 33 | 44.2 MB |", "| `upload/normalized_tables/small_molecule_database.tsv` | 87,552 | 32 | 42.7 MB |")
    text = text.replace("- Recomputed protein-level `unique_drug_count`, `unique_drug_cid_sample`, `source_databases`, `evidence_levels`, and `has_matched_evidence` from the merged binary table.", "- Recomputed protein-level `unique_drug_count`, `unique_drug_cid_sample`, `source_databases`, and `has_matched_evidence` from the merged binary table.")
    text = text.replace("- Removed `cid_mapping_method` from `small_molecule_database.tsv`; mapping method is internal provenance and remains documented in process/audit artifacts rather than the final result table.", "- Removed `cid_mapping_method` from `small_molecule_database.tsv`; mapping method is internal provenance and remains documented in process/audit artifacts rather than the final result table.\n- Removed `evidence_levels` from `protein_database.tsv` and `small_molecule_database.tsv`; it duplicated broad source/category information and did not provide useful evidence grading for these entity summary tables.")
    text = text.replace("| `protein_database.tsv` | 1,773 | 26 |", "| `protein_database.tsv` | 1,773 | 25 |")
    text = text.replace("| `small_molecule_database.tsv` | 87,552 | 33 |", "| `small_molecule_database.tsv` | 87,552 | 32 |")
    path.write_text(text, encoding="utf-8")


def update_full_work_log() -> None:
    path = BASE / "FULL_WORK_LOG.md"
    text = path.read_text(encoding="utf-8")
    text = text.replace("drug_protein_binary_relationships.tsv   95,234 rows", "drug_protein_binary_relationships.tsv   101,171 rows")
    text = text.replace("protein_database.tsv                     1,773 rows 脳 26 cols", "protein_database.tsv                     1,773 rows 脳 25 cols")
    text = text.replace("small_molecule_database.tsv             87,552 rows 脳 33 cols", "small_molecule_database.tsv             87,552 rows 脳 32 cols")
    text = text.replace("| Binary relationships | 95,234 unique (drug, protein) pairs |", "| Binary relationships | 101,171 relationship/context rows |")
    text = text.replace("- Removed `cid_mapping_method` from `small_molecule_database.tsv`; it is internal CID-mapping provenance rather than a final result field.", "- Removed `cid_mapping_method` from `small_molecule_database.tsv`; it is internal CID-mapping provenance rather than a final result field.\n- Removed `evidence_levels` from `protein_database.tsv` and `small_molecule_database.tsv`; it was too coarse for the current entity-summary tables and duplicated source/context information.")
    text = text.replace("| `protein_database.tsv` | 1,773 | 26 |", "| `protein_database.tsv` | 1,773 | 25 |")
    text = text.replace("| `small_molecule_database.tsv` | 87,552 | 33 |", "| `small_molecule_database.tsv` | 87,552 | 32 |")
    path.write_text(text, encoding="utf-8")


def update_manifest() -> None:
    path = BASE / "MANIFEST.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    data.setdefault("summary", {})
    data["summary"].update(
        {
            "binary_relationship_rows": 101171,
            "protein_rows": 1773,
            "small_molecule_rows": 87552,
            "pocket_instance_rows": 6815,
        }
    )
    tracked = [
        "drug_protein_binary_relationships.tsv",
        "protein_database.tsv",
        "protein_source_summaries.tsv",
        "normalized_tables_v3.xlsx",
        "small_molecule_database.tsv",
        "pocket_instances.tsv",
        "DATA_DICTIONARY.md",
        "PROCESS_LOG.md",
        "FULL_WORK_LOG.md",
        "SHEET2_SHEET3_COLUMN_GUIDE.md",
    ]
    rows = []
    for name in tracked:
        p = BASE / name
        rows.append(
            {
                "path": f"upload/normalized_tables/{name}",
                "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                "bytes": p.stat().st_size,
            }
        )
    existing = {item["path"]: item for item in data.get("files", [])}
    for item in rows:
        existing[item["path"]] = item
    data["files"] = list(existing.values())
    data["outputs"] = [item for item in rows if item["path"].endswith((".tsv", ".xlsx")) and Path(item["path"]).name in {
        "drug_protein_binary_relationships.tsv",
        "protein_database.tsv",
        "small_molecule_database.tsv",
        "pocket_instances.tsv",
        "normalized_tables_v3.xlsx",
    }]
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    update_data_dictionary()
    update_column_guide()
    update_process_log()
    update_full_work_log()
    update_manifest()


if __name__ == "__main__":
    main()
