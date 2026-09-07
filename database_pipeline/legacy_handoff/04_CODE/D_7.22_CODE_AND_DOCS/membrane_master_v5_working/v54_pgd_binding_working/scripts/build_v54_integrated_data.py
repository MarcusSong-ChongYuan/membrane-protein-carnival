#!/usr/bin/env python3
"""Build v5.4 protein/disease and v4.1 binding-evidence data releases."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import shutil
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import median


ROOT = Path(r"D:\7.22")
WORK = ROOT / "membrane_master_v5_working" / "v54_pgd_binding_working"
INTERMEDIATE = WORK / "intermediate"
REPORTS = WORK / "reports"
RELEASE = (
    ROOT
    / "membrane_master_v5_working"
    / "releases"
    / "release_v5_4_pgd_binding"
)

V53_MASTER = (
    ROOT
    / "membrane_master_v5_working"
    / "releases"
    / "release_v5_3_sequences"
    / "human_membrane_audit_master_v5_3.tsv"
)
V3_PGD_GZ = (
    ROOT
    / "upload"
    / "normalized_tables"
    / "releases"
    / "protein_gene_disease_v1"
    / "protein_gene_disease_rule_b_scored_v1.tsv.gz"
)
UNIPROT_DISEASE = WORK / "raw" / "uniprot_disease_annotations_v54.tsv"

V4_RELEASE = ROOT / "v4_working" / "releases" / "release_v4"
V4_COMPOUNDS = V4_RELEASE / "small_molecule_master_v4.tsv"
V4_ACTIVITY = V4_RELEASE / "activity_measurements_v4.tsv"
V4_ASSERTIONS = V4_RELEASE / "compound_protein_assertions_v4.tsv"
V4_STRUCTURES = V4_RELEASE / "structure_ligand_observations_v4.tsv"
V4_SITES = V4_RELEASE / "binding_site_annotations_v4.tsv"
V4_GTOPDB = V4_RELEASE / "external_gtopdb_interactions_v4.tsv"
V4_BDB_INDEX = V4_RELEASE / "external_bindingdb_target_index_v4.tsv"
V4_CHEMBL_INDEX = V4_RELEASE / "external_chembl_target_index_v4.tsv"
CHEMBL_INDEX_SUPPLEMENT = (
    INTERMEDIATE / "external_chembl_target_index_v54_supplement.tsv"
)
V3_POCKETS = ROOT / "upload" / "normalized_tables" / "pocket_instances.tsv"

BDB_EVIDENCE = INTERMEDIATE / "bindingdb_article_evidence_v54.tsv"
BDB_COMPOUNDS = INTERMEDIATE / "bindingdb_compounds_v54.tsv"
CHEMBL_ACTIVITY = INTERMEDIATE / "chembl_37_binding_activities_v54.tsv"
CHEMBL_ACTIVITY_SUPPLEMENT = (
    INTERMEDIATE / "chembl_37_binding_activities_v54_supplement.tsv"
)

PGD_OUT = RELEASE / "protein_gene_disease_relations_v5_4.tsv"
PGD_EXCLUDED_OUT = RELEASE / "v3_pgd_relations_not_in_v5_3_audit.tsv"
PGD_PRIMARY_OUT = RELEASE / "protein_primary_disease_v5_4.tsv"
PGD_SUMMARY_OUT = RELEASE / "protein_disease_summary_v5_4.tsv"
BINDING_OUT = RELEASE / "binding_evidence_master_v4_1.tsv"
BINDING_WORKBOOK_VIEW_OUT = RELEASE / "binding_evidence_workbook_view_v4_1.tsv"
BINDING_EXCEL_REVIEW_OUT = RELEASE / "binding_evidence_excel_review_v4_1.tsv"
BINDING_EXCLUDED_OUT = RELEASE / "v3_binding_evidence_not_in_v5_3_audit.tsv"
SITE_OUT = RELEASE / "binding_site_instances_v4_1.tsv"
PAIR_OUT = RELEASE / "protein_compound_summary_v4_1.tsv"
BINDING_SUMMARY_OUT = RELEASE / "protein_binding_summary_v4_1.tsv"
COMPOUND_OUT = RELEASE / "small_molecule_index_v4_1.tsv"
MASTER_OUT = RELEASE / "human_membrane_audit_master_v5_4.tsv"
VALIDATION_OUT = RELEASE / "V54_VALIDATION_REPORT.json"

DISEASE_BLOCK_RE = re.compile(
    r"(?:^|;\s*)DISEASE:\s*(?P<name>.+?)"
    r"(?:\s+\((?P<abbr>[^()]*)\))?"
    r"(?:\s+\[(?P<id>MIM:\d+)\])?:\s*"
    r"(?P<body>.*?)(?=(?:;\s*DISEASE:)|$)"
)
PUBMED_RE = re.compile(r"PubMed:(\d+)")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def iter_tsv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def stable_id(prefix: str, *parts: str) -> str:
    digest = hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:20].upper()
    return f"{prefix}_{digest}"


def unique_tokens(value: str) -> list[str]:
    return sorted(
        {
            token.strip()
            for token in re.split(r"[;,|]", value or "")
            if token.strip()
        }
    )


def is_true(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def number_or_blank(value: str):
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return ""


def master_universe():
    rows = read_tsv(V53_MASTER)
    by_accession = {row["target_uniprot_id"]: row for row in rows}
    return rows, by_accession


def build_disease_relations(master_rows, master_by_accession):
    relations: list[dict] = []
    v3_keys: set[tuple[str, str, str]] = set()
    v3_not_in_master: list[dict] = []
    with gzip.open(V3_PGD_GZ, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        v3_header = list(reader.fieldnames or [])
        for row in reader:
            accession = row["target_uniprot_id"].strip()
            if accession not in master_by_accession:
                v3_not_in_master.append(
                    {
                        **row,
                        "v54_exclusion_reason": (
                            "target accession is not present in the audited "
                            "v5.3 membrane-protein universe"
                        ),
                    }
                )
                continue
            master = master_by_accession[accession]
            relation_id = stable_id(
                "PGD",
                "OpenTargets_v3",
                accession,
                row["ensembl_gene_id"],
                row["disease_id"],
            )
            evidence_level = row["evidence_level"]
            relations.append(
                {
                    "disease_relation_id": relation_id,
                    "target_uniprot_id": accession,
                    "approved_symbol": master.get("approved_symbol", ""),
                    "hgnc_ids": master.get("hgnc_ids", ""),
                    "ncbi_gene_ids": master.get("ncbi_gene_ids", ""),
                    "ensembl_gene_id": row["ensembl_gene_id"],
                    "disease_id": row["disease_id"],
                    "disease_name": row["disease_name"],
                    "disease_ontology": (
                        "MONDO" if row["disease_id"].startswith("MONDO") else ""
                    ),
                    "association_scope": row["association_scope"],
                    "relation_context": row["relation_context"],
                    "disease_evidence_level": evidence_level,
                    "disease_evidence_rank": row["evidence_level_rank"],
                    "highest_rule_passed": row["highest_rule_passed"],
                    "ot_overall_score": row["ot_overall_score"],
                    "primary_disease_flag_v3": row["primary_disease_flag"],
                    "human_genetic_flag": row["human_genetic_flag"],
                    "clinical_genetic_flag": row["clinical_genetic_flag"],
                    "somatic_mutation_flag": row["somatic_mutation_flag"],
                    "functional_flag": row["functional_flag"],
                    "animal_model_flag": row["animal_model_flag"],
                    "expression_flag": row["expression_flag"],
                    "literature_flag": row["literature_flag"],
                    "uniprot_disease_flag": row["uniprot_disease_flag"],
                    "evidence_count": row["evidence_count"],
                    "independent_source_family_count": row[
                        "independent_source_family_count"
                    ],
                    "independent_source_families": row[
                        "independent_source_families"
                    ],
                    "source_database": "Open Targets",
                    "source_version": row["opentargets_release"],
                    "source_record_id": (
                        f"{row['ensembl_gene_id']}|{row['disease_id']}"
                    ),
                    "source_url": (
                        "https://platform.opentargets.org/evidence/"
                        f"{row['ensembl_gene_id']}/{row['disease_id']}"
                    ),
                    "pubmed_ids": "",
                    "record_qc_status": "retained_from_v3_rule_b_release",
                    "default_relation_inclusion": row[
                        "default_release_inclusion"
                    ],
                    "protein_release_tier": master.get("evidence_level_v52", ""),
                    "protein_website_default": master.get("website_default_v52", ""),
                    "relation_source_layer": "v3_rule_b",
                }
            )
            v3_keys.add((accession, row["ensembl_gene_id"], row["disease_id"]))

    write_tsv(
        PGD_EXCLUDED_OUT,
        v3_not_in_master,
        v3_header + ["v54_exclusion_reason"],
    )

    uniprot_by_accession = {
        row["target_uniprot_id"]: row for row in iter_tsv(UNIPROT_DISEASE)
    }
    for accession, source in uniprot_by_accession.items():
        if accession not in master_by_accession:
            continue
        comment = source.get("uniprot_disease_comment", "")
        if not comment:
            continue
        master = master_by_accession[accession]
        ensembl_ids = unique_tokens(master.get("ensembl_gene_ids", ""))
        gene_id = ensembl_ids[0] if ensembl_ids else ""
        for index, match in enumerate(DISEASE_BLOCK_RE.finditer(comment), start=1):
            disease_name = match.group("name").strip()
            if disease_name.lower().startswith("note=") or len(disease_name) > 200:
                # UniProt sometimes emits a disease-section Note as a separate
                # "DISEASE:" block; it is supporting prose, not a disease entity.
                continue
            disease_id = (match.group("id") or "").replace(":", "_")
            if not disease_id:
                disease_id = stable_id("UNIPROT_DISEASE", disease_name)
            body = match.group("body").strip()
            pubmed_ids = sorted(set(PUBMED_RE.findall(body)), key=int)
            causal = "caused by variants affecting" in body.lower()
            evidence_level = "high" if causal and pubmed_ids else "medium"
            evidence_rank = 3 if evidence_level == "high" else 2
            relations.append(
                {
                    "disease_relation_id": stable_id(
                        "PGD", "UniProt", accession, gene_id, disease_id, str(index)
                    ),
                    "target_uniprot_id": accession,
                    "approved_symbol": master.get("approved_symbol", ""),
                    "hgnc_ids": master.get("hgnc_ids", ""),
                    "ncbi_gene_ids": master.get("ncbi_gene_ids", ""),
                    "ensembl_gene_id": gene_id,
                    "disease_id": disease_id,
                    "disease_name": disease_name,
                    "disease_ontology": (
                        "OMIM" if disease_id.startswith("MIM_") else "UniProt"
                    ),
                    "association_scope": "direct_disease_annotation",
                    "relation_context": (
                        "causal_variant_annotation"
                        if causal
                        else "curated_disease_involvement"
                    ),
                    "disease_evidence_level": evidence_level,
                    "disease_evidence_rank": evidence_rank,
                    "highest_rule_passed": "UniProt_curated",
                    "ot_overall_score": "",
                    "primary_disease_flag_v3": "False",
                    "human_genetic_flag": "True" if causal else "False",
                    "clinical_genetic_flag": "True" if causal else "False",
                    "somatic_mutation_flag": "False",
                    "functional_flag": "False",
                    "animal_model_flag": "False",
                    "expression_flag": "False",
                    "literature_flag": "True" if pubmed_ids else "False",
                    "uniprot_disease_flag": "True",
                    "evidence_count": len(pubmed_ids) if pubmed_ids else 1,
                    "independent_source_family_count": 1,
                    "independent_source_families": "UniProtKB",
                    "source_database": "UniProtKB",
                    "source_version": master.get("uniprot_release_v53", "2026_02"),
                    "source_record_id": f"{accession}|{disease_id}",
                    "source_url": source.get("uniprot_disease_source_url", ""),
                    "pubmed_ids": ";".join(pubmed_ids),
                    "record_qc_status": "current_uniprot_curated_disease_comment",
                    "default_relation_inclusion": "True",
                    "protein_release_tier": master.get("evidence_level_v52", ""),
                    "protein_website_default": master.get("website_default_v52", ""),
                    "relation_source_layer": "uniprot_2026_02_supplement",
                }
            )

    relations.sort(
        key=lambda row: (
            row["target_uniprot_id"],
            -int(row["disease_evidence_rank"] or 0),
            row["disease_id"],
            row["source_database"],
        )
    )
    pgd_header = [
        "disease_relation_id",
        "target_uniprot_id",
        "approved_symbol",
        "hgnc_ids",
        "ncbi_gene_ids",
        "ensembl_gene_id",
        "disease_id",
        "disease_name",
        "disease_ontology",
        "association_scope",
        "relation_context",
        "disease_evidence_level",
        "disease_evidence_rank",
        "highest_rule_passed",
        "ot_overall_score",
        "primary_disease_flag_v3",
        "human_genetic_flag",
        "clinical_genetic_flag",
        "somatic_mutation_flag",
        "functional_flag",
        "animal_model_flag",
        "expression_flag",
        "literature_flag",
        "uniprot_disease_flag",
        "evidence_count",
        "independent_source_family_count",
        "independent_source_families",
        "source_database",
        "source_version",
        "source_record_id",
        "source_url",
        "pubmed_ids",
        "record_qc_status",
        "default_relation_inclusion",
        "protein_release_tier",
        "protein_website_default",
        "relation_source_layer",
    ]
    write_tsv(PGD_OUT, relations, pgd_header)

    by_protein: dict[str, list[dict]] = defaultdict(list)
    for relation in relations:
        if is_true(relation["default_relation_inclusion"]):
            by_protein[relation["target_uniprot_id"]].append(relation)

    primary_rows: list[dict] = []
    summary_rows: list[dict] = []
    for master in master_rows:
        accession = master["target_uniprot_id"]
        protein_relations = by_protein.get(accession, [])
        v3_primary = [
            row for row in protein_relations if is_true(row["primary_disease_flag_v3"])
        ]
        primary = (
            sorted(
                v3_primary,
                key=lambda row: (
                    -int(row["disease_evidence_rank"] or 0),
                    -(number_or_blank(row["ot_overall_score"]) or 0),
                    row["disease_id"],
                ),
            )[0]
            if v3_primary
            else (protein_relations[0] if protein_relations else None)
        )
        source_names = sorted({row["source_database"] for row in protein_relations})
        best_rank = max(
            [int(row["disease_evidence_rank"] or 0) for row in protein_relations],
            default=0,
        )
        rank_label = {0: "", 1: "low", 2: "medium", 3: "high", 4: "very_high"}
        summary = {
            "target_uniprot_id": accession,
            "approved_symbol": master.get("approved_symbol", ""),
            "primary_gene_id": (
                primary["ensembl_gene_id"]
                if primary
                else (
                    unique_tokens(master.get("ensembl_gene_ids", ""))[0]
                    if unique_tokens(master.get("ensembl_gene_ids", ""))
                    else ""
                )
            ),
            "has_supported_disease": 1 if protein_relations else 0,
            "supported_disease_count": len(
                {row["disease_id"] for row in protein_relations}
            ),
            "disease_relation_count": len(protein_relations),
            "primary_disease_id": primary["disease_id"] if primary else "",
            "primary_disease_name": primary["disease_name"] if primary else "",
            "best_disease_evidence_level": rank_label[best_rank],
            "best_disease_evidence_rank": best_rank,
            "disease_relation_ids": ";".join(
                row["disease_relation_id"] for row in protein_relations
            ),
            "disease_source_count": len(source_names),
            "disease_sources": ";".join(source_names),
        }
        summary_rows.append(summary)
        if primary:
            primary_rows.append(
                {
                    **summary,
                    "primary_relation_source": primary["source_database"],
                    "primary_relation_context": primary["relation_context"],
                    "primary_relation_id": primary["disease_relation_id"],
                }
            )

    summary_header = list(summary_rows[0])
    write_tsv(PGD_SUMMARY_OUT, summary_rows, summary_header)
    primary_header = list(primary_rows[0]) if primary_rows else summary_header
    write_tsv(PGD_PRIMARY_OUT, primary_rows, primary_header)
    return relations, summary_rows, uniprot_by_accession


def build_compound_index():
    compounds: dict[str, dict] = {}
    for row in iter_tsv(V4_COMPOUNDS):
        key = row["drug_id"]
        compounds[key] = {
            "compound_id": key,
            "compound_id_type": row.get("compound_id_type", ""),
            "compound_name": row.get("drug_name", ""),
            "canonical_smiles": row.get("canonical_smiles", ""),
            "inchikey": row.get("inchikey", ""),
            "pubchem_cids": key if row.get("compound_id_type") == "PubChem CID" else "",
            "chembl_ids": "",
            "small_molecule_scope_status": row.get(
                "small_molecule_scope_status", ""
            ),
            "small_molecule_scope_class": row.get("small_molecule_scope_class", ""),
            "is_core_small_molecule": row.get("is_core_small_molecule", "0"),
            "source_databases": row.get("source_databases", ""),
            "source_versions": "v3.1;GtoPdb 2026.2",
            "record_qc_status": row.get("record_qc_status", ""),
        }
    if BDB_COMPOUNDS.exists():
        for row in iter_tsv(BDB_COMPOUNDS):
            key = row["compound_id"]
            if key not in compounds:
                compounds[key] = {
                    "compound_id": key,
                    "compound_id_type": row["compound_id_type"],
                    "compound_name": row["compound_name"],
                    "canonical_smiles": row["canonical_smiles"],
                    "inchikey": row["inchikey"],
                    "pubchem_cids": row["pubchem_cids"],
                    "chembl_ids": row["chembl_ids"],
                    "small_molecule_scope_status": row[
                        "small_molecule_scope_status"
                    ],
                    "small_molecule_scope_class": row[
                        "small_molecule_scope_class"
                    ],
                    "is_core_small_molecule": row["is_core_small_molecule"],
                    "source_databases": "BindingDB",
                    "source_versions": "2026-07",
                    "record_qc_status": "supplemental_bindingdb_compound",
                }
    for chembl_path in (CHEMBL_ACTIVITY, CHEMBL_ACTIVITY_SUPPLEMENT):
        if chembl_path.exists():
            for row in iter_tsv(chembl_path):
                molecule = row["chembl_molecule_id"]
                key = f"CHEMBL:{molecule}"
                if key not in compounds:
                    compounds[key] = {
                        "compound_id": key,
                        "compound_id_type": "ChEMBL ID",
                        "compound_name": molecule,
                        "canonical_smiles": "",
                        "inchikey": "",
                        "pubchem_cids": "",
                        "chembl_ids": molecule,
                        "small_molecule_scope_status": "accepted_core_provisional",
                        "small_molecule_scope_class": "defined_chembl_entity",
                        "is_core_small_molecule": 1,
                        "source_databases": "ChEMBL",
                        "source_versions": "ChEMBL 37",
                        "record_qc_status": "properties_pending_compound_enrichment",
                    }
    rows = [compounds[key] for key in sorted(compounds)]
    write_tsv(COMPOUND_OUT, rows, list(rows[0]))
    return compounds


def evidence_template(
    evidence_id,
    accession,
    compound_id,
    compound,
    tier,
    evidence_type,
    source_database,
):
    return {
        "evidence_id": evidence_id,
        "target_uniprot_id": accession,
        "approved_symbol": "",
        "protein_release_tier": "",
        "protein_website_default": "",
        "compound_id": compound_id,
        "compound_id_type": compound.get("compound_id_type", ""),
        "compound_name": compound.get("compound_name", ""),
        "small_molecule_scope_status": compound.get(
            "small_molecule_scope_status", ""
        ),
        "is_core_small_molecule": compound.get("is_core_small_molecule", 0),
        "evidence_tier": tier,
        "evidence_type": evidence_type,
        "evidence_directness": "",
        "target_assignment_status": "single_protein_direct",
        "relationship_context_id": "",
        "activity_type": "",
        "activity_relation": "",
        "activity_value": "",
        "activity_unit": "",
        "standard_value_nM": "",
        "assay_or_mechanism": "",
        "pdb_ids": "",
        "binding_site_residues": "",
        "ligand_het_id": "",
        "pubmed_ids": "",
        "doi": "",
        "source_database": source_database,
        "source_version": "",
        "source_record_id": "",
        "source_url": "",
        "record_qc_status": "",
        "default_release_inclusion": 0,
        "originating_dataset": "",
    }


def build_binding_evidence(master_rows, master_by_accession, compounds):
    evidence: list[dict] = []
    binding_excluded: list[dict] = []

    for row in iter_tsv(V4_ACTIVITY):
        accession = row["target_uniprot_id"]
        if accession not in master_by_accession:
            binding_excluded.append(
                {
                    "legacy_record_id": row["activity_measurement_id"],
                    "legacy_record_type": "activity_measurement",
                    "target_uniprot_id": accession,
                    "drug_id": row["drug_id"],
                    "relationship_context_id": row["relationship_context_id"],
                    "v54_exclusion_reason": "target not present in v5.3 protein universe",
                }
            )
            continue
        compound_id = row["drug_id"]
        compound = compounds.get(compound_id, {})
        tier = "BE2" if is_true(row["direct_binding_candidate"]) else "BE3"
        record = evidence_template(
            stable_id("BE", "v3_activity", row["activity_measurement_id"]),
            accession,
            compound_id,
            compound,
            tier,
            (
                "direct_quantitative_binding"
                if tier == "BE2"
                else "functional_or_bioactivity_measurement"
            ),
            "v3.1_mixed_sources",
        )
        value_um = number_or_blank(row["standard_value_uM"])
        record.update(
            {
                "evidence_directness": "direct" if tier == "BE2" else "indirect_or_functional",
                "relationship_context_id": row["relationship_context_id"],
                "activity_type": row["activity_type_original"],
                "activity_relation": row["activity_relation_original"],
                "activity_value": row["activity_value_original"],
                "activity_unit": row["activity_unit_original"],
                "standard_value_nM": value_um * 1000 if value_um != "" else "",
                "source_version": "v3.1 legacy snapshot",
                "source_record_id": row["activity_measurement_id"],
                "record_qc_status": row["measurement_qc_status"],
                "default_release_inclusion": int(
                    is_true(compound.get("is_core_small_molecule", 0))
                ),
                "originating_dataset": "normalized_tables_v3_1",
            }
        )
        evidence.append(record)

    for row in iter_tsv(V4_ASSERTIONS):
        accession = row["target_uniprot_id"]
        if accession not in master_by_accession:
            binding_excluded.append(
                {
                    "legacy_record_id": row["assertion_id"],
                    "legacy_record_type": "curated_assertion",
                    "target_uniprot_id": accession,
                    "drug_id": row["drug_id"],
                    "relationship_context_id": row["relationship_context_id"],
                    "v54_exclusion_reason": "target not present in v5.3 protein universe",
                }
            )
            continue
        compound_id = row["drug_id"]
        compound = compounds.get(compound_id, {})
        record = evidence_template(
            stable_id("BE", "v3_assertion", row["assertion_id"]),
            accession,
            compound_id,
            compound,
            "BE3",
            "curated_pharmacological_target_assertion",
            row["source_database"],
        )
        record.update(
            {
                "evidence_directness": "curated_target_relation",
                "relationship_context_id": row["relationship_context_id"],
                "assay_or_mechanism": row["mechanism_or_description"],
                "source_version": "v3.1 legacy snapshot",
                "source_record_id": row["source_record_id"] or row["assertion_id"],
                "record_qc_status": row["assertion_qc_status"],
                "default_release_inclusion": int(
                    is_true(compound.get("is_core_small_molecule", 0))
                ),
                "originating_dataset": "normalized_tables_v3_1",
            }
        )
        evidence.append(record)

    site_instances: list[dict] = []
    for row in iter_tsv(V4_STRUCTURES):
        accession = row["target_uniprot_id"]
        if accession not in master_by_accession:
            binding_excluded.append(
                {
                    "legacy_record_id": row["structure_observation_id"],
                    "legacy_record_type": "structure_observation",
                    "target_uniprot_id": accession,
                    "drug_id": row["drug_id"],
                    "relationship_context_id": row["relationship_context_id"],
                    "v54_exclusion_reason": "target not present in v5.3 protein universe",
                }
            )
            continue
        compound_id = row["drug_id"]
        compound = compounds.get(compound_id, {})
        all_sites = ";".join(
            value
            for value in [
                row["pdb_biolip_matched_sites"],
                row["pdb_scpdb_matched_sites"],
                row["pdbbind_matched_sites"],
            ]
            if value
        )
        pdb_ids = sorted(
            set(re.findall(r"@([0-9A-Za-z]{4})", all_sites))
            | set(re.findall(r"\b([0-9][A-Za-z0-9]{3})\b", all_sites))
        )
        record = evidence_template(
            stable_id("BE", "v3_structure", row["structure_observation_id"]),
            accession,
            compound_id,
            compound,
            "BE1",
            "compound_specific_structure_context",
            row["source_database"],
        )
        record.update(
            {
                "evidence_directness": "direct_structural_observation",
                "relationship_context_id": row["relationship_context_id"],
                "pdb_ids": ";".join(pdb_ids),
                "binding_site_residues": all_sites,
                "source_version": "v3.1 legacy snapshot",
                "source_record_id": row["structure_observation_id"],
                "record_qc_status": row["observation_qc_status"],
                "default_release_inclusion": int(
                    is_true(compound.get("is_core_small_molecule", 0))
                ),
                "originating_dataset": "normalized_tables_v3_1",
            }
        )
        evidence.append(record)
        site_instances.append(
            {
                "binding_site_instance_id": stable_id(
                    "BSI", row["structure_observation_id"]
                ),
                "evidence_id": record["evidence_id"],
                "target_uniprot_id": accession,
                "compound_id": compound_id,
                "site_type": "experimental_structure_match",
                "pdb_ids": ";".join(pdb_ids),
                "residue_or_site_description": all_sites,
                "site_compound_specificity": "compound_specific",
                "source_database": row["source_database"],
                "record_qc_status": row["observation_qc_status"],
            }
        )

    for row in iter_tsv(V4_SITES):
        accession = row["target_uniprot_id"]
        if accession not in master_by_accession:
            binding_excluded.append(
                {
                    "legacy_record_id": row["binding_site_annotation_id"],
                    "legacy_record_type": "binding_site_annotation",
                    "target_uniprot_id": accession,
                    "drug_id": row["drug_id"],
                    "relationship_context_id": row["relationship_context_id"],
                    "v54_exclusion_reason": "target not present in v5.3 protein universe",
                }
            )
            continue
        compound_id = row["drug_id"]
        compound = compounds.get(compound_id, {})
        exact_compound = not compound_id.startswith("UNMAPPED_")
        tier = "BE1" if exact_compound and row["pdb_ids"] else "BE3"
        record = evidence_template(
            stable_id("BE", "v3_site", row["binding_site_annotation_id"]),
            accession,
            compound_id,
            compound,
            tier,
            (
                "compound_specific_curated_site"
                if tier == "BE1"
                else "protein_level_curated_site_annotation"
            ),
            row["source_database"],
        )
        record.update(
            {
                "evidence_directness": (
                    "direct_site_annotation" if tier == "BE1" else "site_context_only"
                ),
                "relationship_context_id": row["relationship_context_id"],
                "pdb_ids": row["pdb_ids"],
                "binding_site_residues": row["binding_sites"],
                "pubmed_ids": row["pubmed"],
                "source_version": "v3.1 legacy snapshot",
                "source_record_id": row["binding_site_annotation_id"],
                "record_qc_status": row["annotation_qc_status"],
                "default_release_inclusion": int(
                    tier == "BE1"
                    and exact_compound
                    and is_true(compound.get("is_core_small_molecule", 0))
                ),
                "originating_dataset": "normalized_tables_v3_1",
            }
        )
        evidence.append(record)
        site_instances.append(
            {
                "binding_site_instance_id": stable_id(
                    "BSI", row["binding_site_annotation_id"]
                ),
                "evidence_id": record["evidence_id"],
                "target_uniprot_id": accession,
                "compound_id": compound_id,
                "site_type": "uniprot_curated_binding_site",
                "pdb_ids": row["pdb_ids"],
                "residue_or_site_description": row["binding_sites"],
                "site_compound_specificity": (
                    "compound_specific"
                    if exact_compound
                    else "protein_level_or_unmapped_ligand"
                ),
                "source_database": row["source_database"],
                "record_qc_status": row["annotation_qc_status"],
            }
        )

    for row in iter_tsv(V4_GTOPDB):
        accession = row["target_uniprot_id"]
        if accession not in master_by_accession:
            continue
        compound_id = (
            row["pubchem_cid"]
            if row["pubchem_cid"]
            else f"GTOPDB:{row['ligand_id']}"
        )
        compound = compounds.get(compound_id, {})
        tier = "BE2" if row["candidate_sm_level"] == "SM4" else "BE3"
        record = evidence_template(
            stable_id("BE", "gtopdb", row["external_evidence_id"]),
            accession,
            compound_id,
            compound,
            tier,
            (
                "direct_quantitative_binding"
                if tier == "BE2"
                else "curated_pharmacological_interaction"
            ),
            row["source_database"],
        )
        relation = row["original_affinity_relation"]
        value_nm = row["original_affinity_median_nm"]
        if not value_nm:
            value_nm = row["original_affinity_low_nm"] or row["original_affinity_high_nm"]
        record.update(
            {
                "evidence_directness": "direct" if tier == "BE2" else "curated_or_functional",
                "activity_type": row["original_affinity_units"],
                "activity_relation": relation,
                "activity_value": value_nm,
                "activity_unit": "nM" if value_nm else "",
                "standard_value_nM": value_nm,
                "assay_or_mechanism": row["assay_description"]
                or row["interaction_type"],
                "pubmed_ids": row["pubmed_id"],
                "source_version": row["source_version"],
                "source_record_id": row["external_evidence_id"],
                "source_url": row["source_url"],
                "record_qc_status": "ok",
                "default_release_inclusion": int(
                    is_true(compound.get("is_core_small_molecule", 0))
                ),
                "originating_dataset": "GtoPdb_2026_2",
            }
        )
        evidence.append(record)

    if BDB_EVIDENCE.exists():
        for row in iter_tsv(BDB_EVIDENCE):
            accession = row["target_uniprot_id"]
            if accession not in master_by_accession:
                continue
            compound_id = row["compound_id"]
            compound = compounds.get(compound_id, {})
            record = evidence_template(
                row["evidence_id"],
                accession,
                compound_id,
                compound,
                row["evidence_tier"],
                row["evidence_type"],
                row["source_database"],
            )
            value_nm = row["activity_value"] if row["activity_unit"] == "nM" else ""
            record.update(
                {
                    "evidence_directness": (
                        "direct"
                        if row["evidence_tier"] in {"BE1", "BE2"}
                        else "binding_assay_activity"
                    ),
                    "target_assignment_status": row["target_assignment_status"],
                    "activity_type": row["activity_type"],
                    "activity_relation": row["activity_relation"],
                    "activity_value": row["activity_value"],
                    "activity_unit": row["activity_unit"],
                    "standard_value_nM": value_nm,
                    "assay_or_mechanism": row["target_context"],
                    "pdb_ids": row["pdb_ids"],
                    "ligand_het_id": row["ligand_het_id"],
                    "pubmed_ids": row["pubmed_id"],
                    "doi": row["doi"],
                    "source_version": row["source_version"],
                    "source_record_id": row["source_record_id"],
                    "source_url": row["source_url"],
                    "record_qc_status": row["record_qc_status"],
                    "default_release_inclusion": row["default_release_inclusion"],
                    "originating_dataset": "BindingDB_curated_articles_202607",
                }
            )
            evidence.append(record)
            if row["evidence_tier"] == "BE1":
                site_instances.append(
                    {
                        "binding_site_instance_id": stable_id(
                            "BSI", row["evidence_id"]
                        ),
                        "evidence_id": row["evidence_id"],
                        "target_uniprot_id": accession,
                        "compound_id": compound_id,
                        "site_type": "bindingdb_ligand_target_complex",
                        "pdb_ids": row["pdb_ids"],
                        "residue_or_site_description": "",
                        "site_compound_specificity": "compound_specific",
                        "source_database": "BindingDB",
                        "record_qc_status": row["record_qc_status"],
                    }
                )

    for chembl_path in (CHEMBL_ACTIVITY, CHEMBL_ACTIVITY_SUPPLEMENT):
        if chembl_path.exists():
            for row in iter_tsv(chembl_path):
                accession = row["target_uniprot_id"]
                if accession not in master_by_accession:
                    continue
                compound_id = f"CHEMBL:{row['chembl_molecule_id']}"
                compound = compounds.get(compound_id, {})
                tier = "BE2" if row["activity_type"] in {"Ki", "Kd"} else "BE3"
                record = evidence_template(
                    stable_id(
                        "BE", "chembl", str(row["chembl_activity_id"]), accession
                    ),
                    accession,
                    compound_id,
                    compound,
                    tier,
                    (
                        "direct_quantitative_binding"
                        if tier == "BE2"
                        else "binding_assay_activity"
                    ),
                    "ChEMBL",
                )
                record.update(
                    {
                        "evidence_directness": (
                            "direct" if tier == "BE2" else "binding_assay_activity"
                        ),
                        "target_assignment_status": row[
                            "target_assignment_status"
                        ],
                        "relationship_context_id": row["chembl_assay_id"],
                        "activity_type": row["activity_type"],
                        "activity_relation": row["activity_relation"],
                        "activity_value": row["activity_value_nM"],
                        "activity_unit": row["activity_unit"],
                        "standard_value_nM": row["activity_value_nM"],
                        "assay_or_mechanism": row["bao_label"],
                        "source_version": row["source_version"],
                        "source_record_id": row["chembl_activity_id"],
                        "source_url": row["source_url"],
                        "record_qc_status": (
                            "ok"
                            if not row["data_validity_comment"]
                            else f"review:{row['data_validity_comment']}"
                        ),
                        "default_release_inclusion": int(
                            row["target_assignment_status"]
                            == "single_uniprot_mapping"
                            and not row["data_validity_comment"]
                            and str(row["potential_duplicate"]).lower()
                            not in {"1", "true"}
                        ),
                        "originating_dataset": "ChEMBL_37_binding_assays",
                    }
                )
                evidence.append(record)

    # Add protein fields and remove exact duplicate source records.
    deduplicated: dict[str, dict] = {}
    for row in evidence:
        master = master_by_accession[row["target_uniprot_id"]]
        row["approved_symbol"] = master.get("approved_symbol", "")
        row["protein_release_tier"] = master.get("evidence_level_v52", "")
        row["protein_website_default"] = master.get("website_default_v52", "")
        if not is_true(row["protein_website_default"]):
            row["default_release_inclusion"] = 0
        deduplicated[row["evidence_id"]] = row
    evidence = sorted(
        deduplicated.values(),
        key=lambda row: (
            row["target_uniprot_id"],
            {"BE1": 1, "BE2": 2, "BE3": 3}.get(row["evidence_tier"], 9),
            row["compound_id"],
            row["source_database"],
            row["evidence_id"],
        ),
    )
    evidence_header = list(evidence[0])
    write_tsv(BINDING_OUT, evidence, evidence_header)
    workbook_view = []
    for row in evidence:
        if row["originating_dataset"] != "ChEMBL_37_binding_assays":
            workbook_view.append(row)
            continue
        value_nm = number_or_blank(row["standard_value_nM"])
        if (
            is_true(row["default_release_inclusion"])
            and value_nm != ""
            and value_nm <= 10000
            and row["activity_relation"] not in {">", ">="}
        ):
            workbook_view.append(row)
    write_tsv(BINDING_WORKBOOK_VIEW_OUT, workbook_view, evidence_header)
    review_caps = {"BE1": 5, "BE2": 7, "BE3": 3}
    source_priority = {
        "BioLiP": 1,
        "sc-PDB": 2,
        "PDBbind": 3,
        "BindingDB": 4,
        "IUPHAR/BPS Guide to PHARMACOLOGY": 5,
        "ChEMBL": 6,
        "UniProt": 7,
        "v3.1_mixed_sources": 8,
    }
    review_groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in evidence:
        if is_true(row["default_release_inclusion"]):
            review_groups[
                (row["target_uniprot_id"], row["evidence_tier"])
            ].append(row)
    excel_review: list[dict] = []
    for key in sorted(review_groups):
        tier = key[1]
        candidates = sorted(
            review_groups[key],
            key=lambda row: (
                source_priority.get(row["source_database"], 50),
                number_or_blank(row["standard_value_nM"])
                if number_or_blank(row["standard_value_nM"]) != ""
                else float("inf"),
                row["compound_id"],
                row["evidence_id"],
            ),
        )
        excel_review.extend(candidates[: review_caps[tier]])
    excel_review.sort(
        key=lambda row: (
            row["target_uniprot_id"],
            {"BE1": 1, "BE2": 2, "BE3": 3}[row["evidence_tier"]],
            row["source_database"],
            row["compound_id"],
        )
    )
    write_tsv(BINDING_EXCEL_REVIEW_OUT, excel_review, evidence_header)
    write_tsv(
        BINDING_EXCLUDED_OUT,
        binding_excluded,
        [
            "legacy_record_id",
            "legacy_record_type",
            "target_uniprot_id",
            "drug_id",
            "relationship_context_id",
            "v54_exclusion_reason",
        ],
    )
    site_instances.sort(
        key=lambda row: (
            row["target_uniprot_id"],
            row["compound_id"],
            row["binding_site_instance_id"],
        )
    )
    write_tsv(SITE_OUT, site_instances, list(site_instances[0]))

    pocket_counts: Counter[str] = Counter()
    if V3_POCKETS.exists():
        for row in iter_tsv(V3_POCKETS):
            if row["target_uniprot_id"] in master_by_accession:
                pocket_counts[row["target_uniprot_id"]] += 1

    bdb_index = {row["target_uniprot_id"] for row in iter_tsv(V4_BDB_INDEX)}
    chembl_index = {row["target_uniprot_id"] for row in iter_tsv(V4_CHEMBL_INDEX)}
    if CHEMBL_INDEX_SUPPLEMENT.exists():
        chembl_index.update(
            row["target_uniprot_id"] for row in iter_tsv(CHEMBL_INDEX_SUPPLEMENT)
        )
    by_protein: dict[str, list[dict]] = defaultdict(list)
    by_pair: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in evidence:
        if is_true(row["default_release_inclusion"]):
            by_protein[row["target_uniprot_id"]].append(row)
            by_pair[(row["target_uniprot_id"], row["compound_id"])].append(row)

    pair_rows: list[dict] = []
    tier_rank = {"BE1": 1, "BE2": 2, "BE3": 3}
    for (accession, compound_id), rows in sorted(by_pair.items()):
        best = min(rows, key=lambda row: tier_rank[row["evidence_tier"]])
        sources = sorted({row["source_database"] for row in rows})
        numeric_nm = [
            float(row["standard_value_nM"])
            for row in rows
            if number_or_blank(row["standard_value_nM"]) != ""
        ]
        pair_rows.append(
            {
                "target_uniprot_id": accession,
                "approved_symbol": master_by_accession[accession].get(
                    "approved_symbol", ""
                ),
                "compound_id": compound_id,
                "compound_name": compounds.get(compound_id, {}).get(
                    "compound_name", ""
                ),
                "best_binding_evidence_level": best["evidence_tier"],
                "has_BE1_structural_binding": int(
                    any(row["evidence_tier"] == "BE1" for row in rows)
                ),
                "has_BE2_direct_binding": int(
                    any(row["evidence_tier"] == "BE2" for row in rows)
                ),
                "has_BE3_pharmacology": int(
                    any(row["evidence_tier"] == "BE3" for row in rows)
                ),
                "binding_evidence_count": len(rows),
                "independent_source_count": len(sources),
                "independent_sources": ";".join(sources),
                "best_standard_value_nM": min(numeric_nm) if numeric_nm else "",
                "median_standard_value_nM": median(numeric_nm) if numeric_nm else "",
                "evidence_ids": ";".join(row["evidence_id"] for row in rows),
            }
        )
    write_tsv(PAIR_OUT, pair_rows, list(pair_rows[0]))

    protein_summary: list[dict] = []
    for master in master_rows:
        accession = master["target_uniprot_id"]
        rows = by_protein.get(accession, [])
        sources = sorted({row["source_database"] for row in rows})
        best = (
            min(rows, key=lambda row: tier_rank[row["evidence_tier"]])[
                "evidence_tier"
            ]
            if rows
            else ""
        )
        protein_summary.append(
            {
                "target_uniprot_id": accession,
                "approved_symbol": master.get("approved_symbol", ""),
                "has_BE1_structural_binding": int(
                    any(row["evidence_tier"] == "BE1" for row in rows)
                ),
                "has_BE2_direct_binding": int(
                    any(row["evidence_tier"] == "BE2" for row in rows)
                ),
                "has_BE3_pharmacology": int(
                    any(row["evidence_tier"] == "BE3" for row in rows)
                ),
                "has_predicted_pocket": int(pocket_counts[accession] > 0),
                "predicted_pocket_count": pocket_counts[accession],
                "experimental_small_molecule_count": len(
                    {row["compound_id"] for row in rows}
                ),
                "best_binding_evidence_level": best,
                "binding_evidence_count": len(rows),
                "binding_source_count": len(sources),
                "binding_sources": ";".join(sources),
                "bindingdb_target_index_hit": int(accession in bdb_index),
                "chembl_target_index_hit": int(accession in chembl_index),
            }
        )
    write_tsv(BINDING_SUMMARY_OUT, protein_summary, list(protein_summary[0]))
    return evidence, site_instances, pair_rows, protein_summary


def merge_master(
    master_rows,
    disease_summary,
    binding_summary,
    uniprot_by_accession,
):
    disease_by_accession = {
        row["target_uniprot_id"]: row for row in disease_summary
    }
    binding_by_accession = {
        row["target_uniprot_id"]: row for row in binding_summary
    }
    new_fields = [
        "primary_gene_id",
        "has_supported_disease",
        "supported_disease_count",
        "disease_relation_count",
        "primary_disease_id",
        "primary_disease_name",
        "best_disease_evidence_level",
        "best_disease_evidence_rank",
        "disease_relation_ids",
        "disease_source_count",
        "disease_sources",
        "uniprot_mim_ids_v54",
        "uniprot_orphanet_ids_v54",
        "uniprot_disgenet_ids_v54",
        "uniprot_malacards_ids_v54",
        "has_BE1_structural_binding",
        "has_BE2_direct_binding",
        "has_BE3_pharmacology",
        "has_predicted_pocket",
        "predicted_pocket_count",
        "experimental_small_molecule_count",
        "best_binding_evidence_level",
        "binding_evidence_count",
        "binding_source_count",
        "binding_sources",
        "bindingdb_target_index_hit",
        "chembl_target_index_hit",
        "v54_release_date",
        "v54_disease_policy_version",
        "v54_binding_policy_version",
    ]
    merged: list[dict] = []
    for row in master_rows:
        accession = row["target_uniprot_id"]
        disease = disease_by_accession[accession]
        binding = binding_by_accession[accession]
        uniprot = uniprot_by_accession.get(accession, {})
        out = dict(row)
        for field in [
            "primary_gene_id",
            "has_supported_disease",
            "supported_disease_count",
            "disease_relation_count",
            "primary_disease_id",
            "primary_disease_name",
            "best_disease_evidence_level",
            "best_disease_evidence_rank",
            "disease_relation_ids",
            "disease_source_count",
            "disease_sources",
        ]:
            out[field] = disease[field]
        out["uniprot_mim_ids_v54"] = uniprot.get("uniprot_mim_ids", "")
        out["uniprot_orphanet_ids_v54"] = uniprot.get(
            "uniprot_orphanet_ids", ""
        )
        out["uniprot_disgenet_ids_v54"] = uniprot.get(
            "uniprot_disgenet_ids", ""
        )
        out["uniprot_malacards_ids_v54"] = uniprot.get(
            "uniprot_malacards_ids", ""
        )
        for field in [
            "has_BE1_structural_binding",
            "has_BE2_direct_binding",
            "has_BE3_pharmacology",
            "has_predicted_pocket",
            "predicted_pocket_count",
            "experimental_small_molecule_count",
            "best_binding_evidence_level",
            "binding_evidence_count",
            "binding_source_count",
            "binding_sources",
            "bindingdb_target_index_hit",
            "chembl_target_index_hit",
        ]:
            out[field] = binding[field]
        out["v54_release_date"] = "2026-07-24"
        out["v54_disease_policy_version"] = "pgd_v54_ot_rule_b_plus_uniprot_v1"
        out["v54_binding_policy_version"] = "binding_BE123_v1"
        merged.append(out)
    fieldnames = list(master_rows[0]) + new_fields
    write_tsv(MASTER_OUT, merged, fieldnames)
    return merged


def validate(
    master_rows,
    relations,
    disease_summary,
    evidence,
    site_instances,
    pair_rows,
    binding_summary,
    compounds,
):
    master_ids = {row["target_uniprot_id"] for row in master_rows}
    checks = {
        "master_row_count": len(master_rows),
        "master_key_unique": len(master_ids) == len(master_rows),
        "disease_relation_rows": len(relations),
        "disease_relation_key_unique": len(
            {row["disease_relation_id"] for row in relations}
        )
        == len(relations),
        "disease_foreign_key_closed": all(
            row["target_uniprot_id"] in master_ids for row in relations
        ),
        "disease_summary_rows": len(disease_summary),
        "disease_summary_key_closed": {
            row["target_uniprot_id"] for row in disease_summary
        }
        == master_ids,
        "binding_evidence_rows": len(evidence),
        "binding_evidence_key_unique": len(
            {row["evidence_id"] for row in evidence}
        )
        == len(evidence),
        "binding_foreign_key_closed": all(
            row["target_uniprot_id"] in master_ids for row in evidence
        ),
        "binding_site_rows": len(site_instances),
        "protein_compound_summary_rows": len(pair_rows),
        "binding_summary_rows": len(binding_summary),
        "binding_summary_key_closed": {
            row["target_uniprot_id"] for row in binding_summary
        }
        == master_ids,
        "compound_index_rows": len(compounds),
        "binding_tier_counts": dict(Counter(row["evidence_tier"] for row in evidence)),
        "binding_source_counts": dict(
            Counter(row["source_database"] for row in evidence)
        ),
        "disease_source_counts": dict(
            Counter(row["source_database"] for row in relations)
        ),
        "proteins_with_supported_disease": sum(
            int(row["has_supported_disease"]) for row in disease_summary
        ),
        "proteins_with_binding_evidence": sum(
            bool(row["best_binding_evidence_level"]) for row in binding_summary
        ),
    }
    boolean_checks = [
        value for key, value in checks.items() if isinstance(value, bool)
    ]
    checks["status"] = "passed" if all(boolean_checks) else "failed"
    checks["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    VALIDATION_OUT.write_text(json.dumps(checks, indent=2), encoding="utf-8")
    REPORTS.mkdir(parents=True, exist_ok=True)
    shutil.copy2(VALIDATION_OUT, REPORTS / VALIDATION_OUT.name)
    return checks


def main() -> None:
    RELEASE.mkdir(parents=True, exist_ok=True)
    if not UNIPROT_DISEASE.exists():
        raise RuntimeError("UniProt disease retrieval output is missing")
    if not BDB_EVIDENCE.exists():
        raise RuntimeError("BindingDB normalized evidence output is missing")
    if not CHEMBL_ACTIVITY.exists():
        raise RuntimeError("ChEMBL normalized activity output is missing")
    if not CHEMBL_ACTIVITY_SUPPLEMENT.exists():
        raise RuntimeError("Supplemental ChEMBL activity output is missing")

    master_rows, master_by_accession = master_universe()
    relations, disease_summary, uniprot_by_accession = build_disease_relations(
        master_rows, master_by_accession
    )
    compounds = build_compound_index()
    evidence, site_instances, pair_rows, binding_summary = build_binding_evidence(
        master_rows, master_by_accession, compounds
    )
    merge_master(
        master_rows, disease_summary, binding_summary, uniprot_by_accession
    )
    report = validate(
        master_rows,
        relations,
        disease_summary,
        evidence,
        site_instances,
        pair_rows,
        binding_summary,
        compounds,
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
