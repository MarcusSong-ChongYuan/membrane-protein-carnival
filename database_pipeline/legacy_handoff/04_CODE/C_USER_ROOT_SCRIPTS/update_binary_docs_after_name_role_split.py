from __future__ import annotations

import hashlib
import json
from pathlib import Path


BASE = Path(r"C:\github-repos\upload\normalized_tables")

BINARY_SECTION = """## `drug_protein_binary_relationships.tsv`

Column order: protein identity -> compound identity/role -> relationship/evidence -> binding-site details.

### Protein Identity (1-4)
| # | Column | Meaning |
|---|---|---|
| 1 | `target_uniprot_id` | UniProt accession; primary protein identifier. |
| 2 | `approved_symbol` | HGNC-approved gene symbol. |
| 3 | `ncbi_gene_id` | NCBI Gene identifier mapped for the protein. |
| 4 | `protein_name` | Human-readable protein name. |

### Compound Identity and Role (5-10)
| # | Column | Meaning |
|---|---|---|
| 5 | `drug_id` | Compound identifier. Usually PubChem CID; unmapped biologics retain ChEMBL/DrugCentral ID; unmapped structure ligands use stable IDs such as `UNMAPPED_PDB_LIGAND:ITC`. |
| 6 | `compound_id_type` | `PubChem CID`, original type for unmapped biologics, or `unmapped_structure_ligand`. |
| 7 | `drug_name` | Preferred compound display name. |
| 8 | `compound_name_category` | Type of display name: `common_or_drug_name`, `research_code`, `systematic_or_chemical_name`, `pubchem_synonym_name`, `pdb_ligand_code_or_name`, or `unknown_name_type`. |
| 9 | `compound_biological_role` | Broad role copied from the small-molecule entity table, such as `therapeutic_or_clinical_compound`, `bioactive_research_ligand`, `structure_affinity_ligand`, `endogenous_ligand_or_cofactor`, `membrane_lipid_or_sterol`, `biologic_or_large_molecule`, or `unknown_structure_ligand`. |
| 10 | `compound_source` | Relationship evidence origin: `curated_target_relation`, `bioassay_active_relation`, or `structure_ligand_context`. |

### Relationship Evidence (11-15)
| # | Column | Meaning |
|---|---|---|
| 11 | `source_database` | Originating database(s) for this compound-protein row. |
| 12 | `assay_or_mechanism` | Assay type, mechanism-of-action text, or structure-context descriptor. |
| 13 | `activity_type` | `IC50`, `Ki`, `Kd`, `EC50`, etc. |
| 14 | `activity_value_uM` | Activity value in uM when available. |
| 15 | `clinical_or_approval_status` | Row-level captured status from source relationship data: `approved_or_listed`, `clinical_trial`, combined value, or empty. |

### Matched Binding Evidence and Details (16-30)
Boolean flags use `1` = present and `0` = absent.
| # | Column | Meaning |
|---|---|---|
| 16 | `has_pdb_biolip_matched` | BioLiP co-crystal/context evidence for this row. |
| 17 | `pdb_biolip_matched_sites` | BioLiP site detail, usually `LIGAND(CODE)@PDB:CHAIN=RESIDUES`. |
| 18 | `has_pdb_scpdb_matched` | sc-PDB co-crystal/context evidence for this row. |
| 19 | `pdb_scpdb_matched_sites` | sc-PDB site detail. |
| 20 | `has_pdbbind_matched` | PDBbind affinity/structure context for this row. |
| 21 | `pdbbind_matched_sites` | PDBbind site/affinity detail. |
| 22 | `has_stitch_matched` | STITCH evidence for this compound-protein pair. |
| 23 | `stitch_matched_compounds` | `CIDxxx:confidence`; ChEMBL prefixes from the original STITCH mapping are removed. |
| 24 | `has_exp_binding_site` | UniProt literature binding-site annotation exists for this protein. |
| 25 | `exp_binding_sites` | Protein-level residue binding-site annotations. |
| 26 | `exp_ligands` | Known ligands from UniProt experimental annotation. |
| 27 | `exp_evidence_quality` | UniProt/ECO evidence quality tiers. |
| 28 | `exp_pubmed` | PubMed IDs for binding-site evidence. |
| 29 | `exp_pdb` | Cross-referenced PDB entries when captured. |
| 30 | `target_group_id` | Group ID preserving source rows that originally listed multiple UniProt targets. |
"""


def replace_between(text: str, start: str, end: str, replacement: str) -> str:
    a = text.index(start)
    b = text.index(end, a)
    return text[:a] + replacement.rstrip() + "\n\n" + text[b:]


def update_manifest() -> None:
    path = BASE / "MANIFEST.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    tracked = [Path(item["path"]).name for item in data.get("files", [])]
    items = []
    for name in tracked:
        p = BASE / name
        if p.exists():
            items.append(
                {
                    "path": f"upload/normalized_tables/{name}",
                    "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                    "bytes": p.stat().st_size,
                }
            )
    data["files"] = items
    data["outputs"] = [
        item
        for item in items
        if Path(item["path"]).name
        in {
            "drug_protein_binary_relationships.tsv",
            "protein_database.tsv",
            "small_molecule_database.tsv",
            "pocket_instances.tsv",
            "normalized_tables_v3.xlsx",
        }
    ]
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    data_dict = BASE / "DATA_DICTIONARY.md"
    text = data_dict.read_text(encoding="utf-8")
    text = text.replace(
        "| `drug_protein_binary_relationships.tsv` | 101,171 | 29 | 50 MB | One row per unique drug-protein relationship, including curated structure ligand context rows. |",
        "| `drug_protein_binary_relationships.tsv` | 101,171 | 30 | 53.0 MB | One row per compound-protein relationship/context row, including curated structure ligand context rows. |",
    )
    text = replace_between(text, "## `drug_protein_binary_relationships.tsv`", "---\n\n## `protein_database.tsv`", BINARY_SECTION)
    data_dict.write_text(text, encoding="utf-8")

    for name in ["PROCESS_LOG.md", "FULL_WORK_LOG.md"]:
        p = BASE / name
        text = p.read_text(encoding="utf-8")
        text = text.replace("| `upload/normalized_tables/drug_protein_binary_relationships.tsv` | 101,171 | 29 | 53.0 MB |", "| `upload/normalized_tables/drug_protein_binary_relationships.tsv` | 101,171 | 30 | 53.0 MB |")
        text = text.replace("| `drug_protein_binary_relationships.tsv` | 101,171 | 29 |", "| `drug_protein_binary_relationships.tsv` | 101,171 | 30 |")
        text = text.replace("drug_protein_binary_relationships.tsv   101,171 rows 脳 29 cols", "drug_protein_binary_relationships.tsv   101,171 rows 脳 30 cols")
        text = text.replace(
            "- Split `drug_name_type` into `compound_name_category` and `compound_biological_role`.",
            "- Split `drug_name_type` into `compound_name_category` and `compound_biological_role` in both the binary and small-molecule tables.",
        )
        text = text.replace(
            "split `drug_name_type` into `compound_name_category` plus `compound_biological_role`.",
            "split `drug_name_type` into `compound_name_category` plus `compound_biological_role` in both the binary and small-molecule tables.",
        )
        p.write_text(text, encoding="utf-8")

    update_manifest()


if __name__ == "__main__":
    main()
