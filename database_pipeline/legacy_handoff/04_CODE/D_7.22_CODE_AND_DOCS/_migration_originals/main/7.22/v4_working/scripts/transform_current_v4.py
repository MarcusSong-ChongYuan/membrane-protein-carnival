from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "v4_working"
INTERMEDIATE = WORK / "intermediate"
REPORTS = WORK / "reports"
V31 = (
    ROOT
    / "upload"
    / "normalized_tables"
    / "releases"
    / "partner_handoff_v3_1"
    / "v3_1"
)
COMPOUNDS = V31 / "small_molecule_database_v3_1.tsv"
RELATIONSHIPS = V31 / "drug_protein_binary_relationships_v3_1.tsv"
MEMBRANE_MASTER = INTERMEDIATE / "human_membrane_protein_master.tsv"
GENERATED_AT = datetime.now(timezone.utc).isoformat()


def stable_id(prefix: str, *values: str) -> str:
    payload = "\x1f".join(values).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(payload).hexdigest()[:18].upper()}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if fields is None:
        if not rows:
            raise ValueError(f"No rows and no fields for {path}")
        fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def classify_compound(row: dict[str, str]) -> tuple[str, str, str]:
    role = row["compound_biological_role"]
    id_type = row["compound_id_type"]
    if role == "biologic_or_large_molecule":
        return "excluded", "biologic_or_large_molecule", "v3.1 biological-role classification"
    if role == "ion_or_metal":
        return "separate_not_core", "inorganic_ion_or_metal", "v3.1 biological-role classification"
    if id_type in {"unmapped_structure_ligand", "unmapped_uniprot_ligand"}:
        if row["canonical_smiles"] or row["inchikey"] or row["molecular_formula"]:
            return "accepted_core", "defined_structure_ligand", "defined chemical structure present"
        return "review", "unresolved_ligand_identity", "unmapped ligand without normalized chemical structure"
    mapping = {
        "therapeutic_or_clinical_compound": "therapeutic_or_clinical_compound",
        "bioactive_research_ligand": "research_ligand",
        "structure_affinity_ligand": "defined_structure_ligand",
        "endogenous_ligand_or_cofactor": "endogenous_metabolite_or_cofactor",
        "membrane_lipid_or_sterol": "lipid_or_sterol",
        "unknown_structure_ligand": "defined_structure_ligand",
    }
    return "accepted_core", mapping.get(role, "defined_chemical_entity"), "normalized v3.1 compound record"


def evidence_rank(row: dict[str, str]) -> tuple[str, int, str]:
    evidence_type = row["evidence_type"]
    activity_type = row["activity_type"].strip().lower()
    activity_value = row["activity_value_uM"].strip()
    source_count = len([x for x in row["source_database"].split(";") if x.strip()])
    if evidence_type == "assay_activity" and activity_value and activity_type in {"ki", "kd"}:
        if source_count >= 2:
            return "SM5", 5, "quantitative Ki/Kd with multiple source labels"
        return "SM4", 4, "quantitative Ki/Kd assay"
    if evidence_type == "curated_target_assertion":
        return "SM3", 3, "curated pharmacological target assertion"
    if evidence_type == "assay_activity":
        return "SM2", 2, "functional or bioactivity assay"
    if evidence_type in {"structure_ligand_observation", "uniprot_binding_site_annotation"}:
        return "SM1", 1, "structure context or binding-site annotation"
    return "SM0", 0, "unclassified evidence type"


def parse_pmids(value: str) -> list[str]:
    return sorted(set(re.findall(r"\b\d{6,9}\b", value or "")))


def main() -> None:
    INTERMEDIATE.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    compound_rows = read_tsv(COMPOUNDS)
    relation_rows = read_tsv(RELATIONSHIPS)
    membrane_rows = read_tsv(MEMBRANE_MASTER)
    membrane_ids = {row["target_uniprot_id"] for row in membrane_rows}

    compounds = []
    compound_status: dict[str, str] = {}
    for row in compound_rows:
        status, scope_class, reason = classify_compound(row)
        compound_status[row["drug_id"]] = status
        compounds.append(
            {
                **row,
                "small_molecule_scope_status": status,
                "small_molecule_scope_class": scope_class,
                "small_molecule_scope_reason": reason,
                "is_core_small_molecule": "1" if status == "accepted_core" else "0",
                "v4_scope_policy_version": "membrane-small-molecule-scope-v1",
            }
        )
    write_tsv(INTERMEDIATE / "small_molecule_master.tsv", compounds)
    write_tsv(
        INTERMEDIATE / "excluded_non_small_molecules.tsv",
        [row for row in compounds if row["small_molecule_scope_status"] == "excluded"],
        list(compounds[0]),
    )
    write_tsv(
        INTERMEDIATE / "compound_review_queue.tsv",
        [
            row
            for row in compounds
            if row["small_molecule_scope_status"] in {"review", "separate_not_core"}
        ],
        list(compounds[0]),
    )

    evidence_export = []
    assertions = []
    measurements = []
    structures = []
    binding_sites = []
    assays_by_id: dict[str, dict] = {}
    publications_by_id: dict[str, dict] = {}
    sources = set()
    relationship_qc = Counter()
    evidence_levels = Counter()
    mapped_relation_count = 0
    core_relation_count = 0

    for row in relation_rows:
        level, level_rank, reason = evidence_rank(row)
        protein_in_reference = row["target_uniprot_id"] in membrane_ids
        is_core_compound = compound_status.get(row["drug_id"]) == "accepted_core"
        if protein_in_reference:
            mapped_relation_count += 1
        if protein_in_reference and is_core_compound:
            core_relation_count += 1
        if not protein_in_reference:
            relationship_qc["protein_not_in_reviewed_membrane_reference"] += 1
        if row["drug_id"] not in compound_status:
            relationship_qc["compound_not_in_v3_1_master"] += 1
        if not is_core_compound:
            relationship_qc["noncore_or_review_compound"] += 1
        if not row["source_database"]:
            relationship_qc["missing_source_database"] += 1
        evidence_levels[level] += 1

        assay_id = ""
        if row["evidence_type"] == "assay_activity":
            assay_id = stable_id(
                "ASSAY",
                row["source_database"],
                row["assay_or_mechanism"],
                row["activity_type"],
                row["target_uniprot_id"],
            )
            assays_by_id.setdefault(
                assay_id,
                {
                    "assay_id": assay_id,
                    "source_database": row["source_database"],
                    "source_assay_id": "",
                    "target_uniprot_id": row["target_uniprot_id"],
                    "assay_description": row["assay_or_mechanism"],
                    "assay_type_original": row["activity_type"],
                    "assay_category": "bioactivity_assay",
                    "target_type": "single_protein_assumed_from_v3_1",
                    "target_organism_id": "9606",
                    "assay_organism_id": "",
                    "source_record_qc_status": "missing_source_assay_id",
                },
            )
            measurements.append(
                {
                    "activity_measurement_id": stable_id(
                        "ACT",
                        row["relationship_context_id"],
                        row["activity_type"],
                        row["activity_value_uM"],
                    ),
                    "relationship_context_id": row["relationship_context_id"],
                    "assay_id": assay_id,
                    "target_uniprot_id": row["target_uniprot_id"],
                    "drug_id": row["drug_id"],
                    "activity_type_original": row["activity_type"],
                    "activity_relation_original": row["activity_relation"],
                    "activity_value_original": row["activity_value_uM"],
                    "activity_unit_original": row["activity_value_unit"],
                    "standard_value_uM": row["activity_value_uM"],
                    "standardization_method": "retained_from_v3_1",
                    "direct_binding_candidate": (
                        "1"
                        if row["activity_type"].strip().lower() in {"ki", "kd"}
                        and row["activity_value_uM"].strip()
                        else "0"
                    ),
                    "measurement_qc_status": (
                        "missing_original_source_value_or_unit"
                        if not row["activity_value_uM"] or not row["activity_value_unit"]
                        else "v3_1_standardized_value_available"
                    ),
                }
            )
        elif row["evidence_type"] == "curated_target_assertion":
            assertions.append(
                {
                    "assertion_id": stable_id("ASSERT", row["relationship_context_id"]),
                    "relationship_context_id": row["relationship_context_id"],
                    "target_uniprot_id": row["target_uniprot_id"],
                    "drug_id": row["drug_id"],
                    "assertion_type": "curated_pharmacological_target",
                    "mechanism_or_description": row["assay_or_mechanism"],
                    "clinical_or_approval_status": row["clinical_or_approval_status"],
                    "source_database": row["source_database"],
                    "source_record_id": "",
                    "assertion_qc_status": "missing_source_record_id",
                }
            )
        elif row["evidence_type"] == "structure_ligand_observation":
            structures.append(
                {
                    "structure_observation_id": stable_id("STRUCT", row["relationship_context_id"]),
                    "relationship_context_id": row["relationship_context_id"],
                    "target_uniprot_id": row["target_uniprot_id"],
                    "drug_id": row["drug_id"],
                    "source_database": row["source_database"],
                    "pdb_biolip_matched_sites": row["pdb_biolip_matched_sites"],
                    "pdb_scpdb_matched_sites": row["pdb_scpdb_matched_sites"],
                    "pdbbind_matched_sites": row["pdbbind_matched_sites"],
                    "evidence_interpretation": "structure_context_not_automatically_pharmacological",
                    "observation_qc_status": row["record_qc_status"],
                }
            )
        elif row["evidence_type"] == "uniprot_binding_site_annotation":
            binding_sites.append(
                {
                    "binding_site_annotation_id": stable_id("SITE", row["relationship_context_id"]),
                    "relationship_context_id": row["relationship_context_id"],
                    "target_uniprot_id": row["target_uniprot_id"],
                    "drug_id": row["drug_id"],
                    "binding_sites": row["exp_binding_sites"],
                    "ligands": row["exp_ligands"],
                    "evidence_quality": row["exp_evidence_quality"],
                    "pubmed": row["exp_pubmed"],
                    "pdb_ids": row["exp_pdb"],
                    "source_database": row["source_database"],
                    "annotation_qc_status": row["record_qc_status"],
                }
            )

        for source in (x.strip() for x in row["source_database"].split(";")):
            if source:
                sources.add(source)
        for pmid in parse_pmids(row["exp_pubmed"]):
            publications_by_id.setdefault(
                pmid,
                {
                    "publication_id": f"PMID:{pmid}",
                    "pmid": pmid,
                    "doi": "",
                    "patent_id": "",
                    "citation": "",
                    "source": "v3.1 exp_pubmed",
                },
            )
        evidence_export.append(
            {
                **row,
                "protein_in_reviewed_membrane_reference": "1" if protein_in_reference else "0",
                "small_molecule_scope_status": compound_status.get(row["drug_id"], "missing"),
                "qualifies_for_core_small_molecule_flag": (
                    "1" if protein_in_reference and is_core_compound else "0"
                ),
                "candidate_sm_level": level,
                "candidate_sm_level_rank": level_rank,
                "candidate_sm_level_reason": reason,
                "assay_id": assay_id,
                "v4_evidence_policy_version": "membrane-small-molecule-scope-v1",
            }
        )

    write_tsv(INTERMEDIATE / "compound_protein_evidence_export.tsv", evidence_export)
    write_tsv(
        INTERMEDIATE / "compound_protein_assertions.tsv",
        assertions,
        list(assertions[0]),
    )
    write_tsv(
        INTERMEDIATE / "activity_measurements.tsv",
        measurements,
        list(measurements[0]),
    )
    write_tsv(
        INTERMEDIATE / "structure_ligand_observations.tsv",
        structures,
        list(structures[0]),
    )
    write_tsv(
        INTERMEDIATE / "binding_site_annotations.tsv",
        binding_sites,
        list(binding_sites[0]),
    )
    write_tsv(
        INTERMEDIATE / "assays.tsv",
        sorted(assays_by_id.values(), key=lambda row: row["assay_id"]),
        list(next(iter(assays_by_id.values()))),
    )
    write_tsv(
        INTERMEDIATE / "publications.tsv",
        sorted(publications_by_id.values(), key=lambda row: row["publication_id"]),
        ["publication_id", "pmid", "doi", "patent_id", "citation", "source"],
    )
    source_rows = [
        {
            "source_id": stable_id("SOURCE", source),
            "source_name": source,
            "source_version": "",
            "retrieval_date": "",
            "source_role": "retained_from_v3_1",
            "license_status": "requires_source_registry_review",
        }
        for source in sorted(sources)
    ]
    write_tsv(
        INTERMEDIATE / "sources.tsv",
        source_rows,
        list(source_rows[0]),
    )

    compound_scope_counts = Counter(
        row["small_molecule_scope_status"] for row in compounds
    )
    report = {
        "status": "passed",
        "generated_at_utc": GENERATED_AT,
        "compound_rows": len(compounds),
        "compound_scope_counts": dict(sorted(compound_scope_counts.items())),
        "relationship_rows": len(evidence_export),
        "relationships_mapping_to_reviewed_membrane_reference": mapped_relation_count,
        "relationships_qualifying_for_core_small_molecule_flag": core_relation_count,
        "evidence_type_counts": dict(Counter(row["evidence_type"] for row in relation_rows)),
        "candidate_sm_level_counts": dict(sorted(evidence_levels.items())),
        "assay_rows": len(assays_by_id),
        "activity_measurement_rows": len(measurements),
        "curated_assertion_rows": len(assertions),
        "structure_observation_rows": len(structures),
        "binding_site_annotation_rows": len(binding_sites),
        "publication_rows": len(publications_by_id),
        "source_rows": len(source_rows),
        "relationship_qc_counts": dict(sorted(relationship_qc.items())),
        "known_limitations": [
            "Most v3.1 rows do not retain a source record identifier.",
            "Many activity rows retain standardized uM values but not original source values and units.",
            "Assay organism, construct, mutation, and target-complex context are largely unavailable in v3.1.",
            "SM4 candidates are conservatively restricted to quantitative Ki or Kd rows.",
        ],
    }
    (REPORTS / "CURRENT_DATA_TRANSFORMATION_REPORT.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
