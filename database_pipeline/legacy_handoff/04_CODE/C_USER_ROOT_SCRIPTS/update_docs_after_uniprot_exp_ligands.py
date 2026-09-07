from __future__ import annotations

import hashlib
import json
from pathlib import Path


BASE = Path(r"C:\github-repos\upload\normalized_tables")


def hash_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def update_text_files() -> None:
    replacements = {
        "101,171 | 31": "101,802 | 31",
        "87,552 | 35": "87,611 | 35",
        "101,171 rows": "101,802 rows",
        "87,552 rows": "87,611 rows",
        "101171": "101802",
        "87552": "87611",
    }
    note = (
        "- Split explicit UniProt experimental binding-site ligands into binary rows with "
        "`source_database = UniProt` and `compound_source = uniprot_binding_site_annotation`; "
        "generic/non-specific terms were excluded to `uniprot_exp_ligand_excluded_terms.tsv`.\n"
    )
    for name in ["DATA_DICTIONARY.md", "SHEET2_SHEET3_COLUMN_GUIDE.md", "PROCESS_LOG.md", "FULL_WORK_LOG.md"]:
        p = BASE / name
        if not p.exists():
            continue
        text = p.read_text(encoding="utf-8")
        for old, new in replacements.items():
            text = text.replace(old, new)
        if name in {"PROCESS_LOG.md", "FULL_WORK_LOG.md"} and "uniprot_binding_site_annotation" not in text:
            marker = "## Known Limitations" if "## Known Limitations" in text else "## Rebuild Command"
            if marker in text:
                text = text.replace(marker, note + "\n" + marker)
            else:
                text += "\n" + note
        if name == "DATA_DICTIONARY.md":
            text = text.replace(
                "`curated_target_relation`, `bioassay_active_relation`, or `structure_ligand_context`",
                "`curated_target_relation`, `bioassay_active_relation`, `structure_ligand_context`, or `uniprot_binding_site_annotation`",
            )
            text = text.replace(
                "`PubChem CID`, original type for unmapped biologics, or `unmapped_structure_ligand`",
                "`PubChem CID`, original type for unmapped biologics, `unmapped_structure_ligand`, or `unmapped_uniprot_ligand`",
            )
            text = text.replace(
                "`structure_affinity_ligand`, `approved_or_clinical_target_ligand`, `cofactor_or_functional_ligand`, `membrane_lipid_or_sterol`, `curated_target_ligand`, `structure_bound_ligand`, `bioassay_active_ligand`, `ion_or_metal`, `buffer_salt_solvent`, and `unknown_structure_ligand`",
                "`structure_affinity_ligand`, `approved_or_clinical_target_ligand`, `cofactor_or_functional_ligand`, `membrane_lipid_or_sterol`, `uniprot_experimental_site_ligand`, `curated_target_ligand`, `structure_bound_ligand`, `bioassay_active_ligand`, `ion_or_metal`, `buffer_salt_solvent`, and `unknown_structure_ligand`",
            )
        p.write_text(text, encoding="utf-8")


def update_manifest() -> None:
    path = BASE / "MANIFEST.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    data.setdefault("summary", {}).update(
        {
            "binary_relationship_rows": 101802,
            "small_molecule_rows": 87611,
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
        "data_dictionary_review.tsv",
        "data_dictionary_review.xlsx",
        "data_dictionary_review_zh.tsv",
        "data_dictionary_review_zh.xlsx",
        "uniprot_exp_ligand_excluded_terms.tsv",
        "uniprot_exp_ligand_merge_report.json",
        "drug_name_normalization_report.tsv",
        "protein_structure_ligand_context_review.tsv",
        "protein_structure_ligand_context_review.xlsx",
        "protein_structure_ligand_context_review_summary.json",
        "protein_structure_ligand_context_binary_format_review.tsv",
        "protein_structure_ligand_context_binary_format_review.xlsx",
        "structure_context_merge_report.json",
        "pubchem_structure_context_properties_cache.json",
        "stitch_cid_only_normalization_report.json",
        "structure_ligand_cid_supplement_report.json",
        "structure_ligand_remaining_unmapped_nonempty.tsv",
        "uniprot_exp_ligand_pubchem_cache.json",
    ]
    items = []
    for name in tracked:
        p = BASE / name
        if p.exists():
            items.append({"path": f"upload/normalized_tables/{name}", "sha256": hash_file(p), "bytes": p.stat().st_size})
    data["files"] = items
    data["outputs"] = [item for item in items if Path(item["path"]).name in {
        "drug_protein_binary_relationships.tsv",
        "protein_database.tsv",
        "small_molecule_database.tsv",
        "pocket_instances.tsv",
        "normalized_tables_v3.xlsx",
    }]
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    update_text_files()
    update_manifest()


if __name__ == "__main__":
    main()
