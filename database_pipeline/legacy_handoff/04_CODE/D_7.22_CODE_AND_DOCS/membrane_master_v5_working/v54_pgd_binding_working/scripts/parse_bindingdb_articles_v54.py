#!/usr/bin/env python3
"""Normalize BindingDB-curated article measurements for membrane proteins."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(r"D:\7.22")
MASTER = (
    ROOT
    / "membrane_master_v5_working"
    / "releases"
    / "release_v5_3_sequences"
    / "human_membrane_audit_master_v5_3.tsv"
)
SOURCE = (
    ROOT
    / "membrane_master_v5_working"
    / "v54_pgd_binding_working"
    / "raw"
    / "bindingdb_articles_202607"
    / "BindingDB_BindingDB_Articles.tsv"
)
WORK = ROOT / "membrane_master_v5_working" / "v54_pgd_binding_working"
EVIDENCE_OUT = WORK / "intermediate" / "bindingdb_article_evidence_v54.tsv"
COMPOUND_OUT = WORK / "intermediate" / "bindingdb_compounds_v54.tsv"
REPORT = WORK / "reports" / "BINDINGDB_ARTICLE_INTEGRATION_REPORT.json"

MEASUREMENT_FIELDS = {
    "Ki (nM)": ("Ki", "BE2"),
    "Kd (nM)": ("Kd", "BE2"),
    "IC50 (nM)": ("IC50", "BE3"),
    "EC50 (nM)": ("EC50", "BE3"),
}
VALUE_RE = re.compile(
    r"^\s*(?P<relation>[<>=~]{0,2})\s*(?P<value>[0-9]+(?:\.[0-9]+)?(?:[Ee][+-]?[0-9]+)?)"
)


def read_master() -> set[str]:
    with MASTER.open("r", encoding="utf-8-sig", newline="") as handle:
        return {
            row["target_uniprot_id"].strip()
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("target_uniprot_id", "").strip()
        }


def split_ids(value: str) -> list[str]:
    return [
        token.strip()
        for token in re.split(r"[,;\s]+", value or "")
        if token.strip()
    ]


def parse_value(raw: str) -> tuple[str, str] | None:
    match = VALUE_RE.search(raw or "")
    if not match:
        return None
    return match.group("relation") or "=", match.group("value")


def compound_identity(row: dict[str, str]) -> tuple[str, str, str]:
    pubchem_ids = split_ids(row.get("PubChem CID", ""))
    chembl_ids = split_ids(row.get("ChEMBL ID of Ligand", ""))
    inchikey = (row.get("Ligand InChI Key", "") or "").strip()
    monomer = (row.get("BindingDB MonomerID", "") or "").strip()
    if pubchem_ids:
        # v3.1 and GtoPdb use the bare PubChem CID as the canonical key.
        return pubchem_ids[0], "PubChem CID", pubchem_ids[0]
    if chembl_ids:
        return f"CHEMBL:{chembl_ids[0]}", "ChEMBL ID", chembl_ids[0]
    if inchikey:
        return f"INCHIKEY:{inchikey}", "InChIKey", inchikey
    return f"BINDINGDB:{monomer}", "BindingDB MonomerID", monomer


def compound_scope(row: dict[str, str]) -> tuple[str, str, str, int]:
    smiles = (row.get("Ligand SMILES", "") or "").strip()
    name = (row.get("BindingDB Ligand Name", "") or "").lower()
    if not smiles:
        return "review", "missing_structure", "no ligand SMILES", 0
    if len(smiles) > 1000:
        return "review", "large_or_polymeric_entity", "SMILES length > 1000", 0
    excluded_tokens = ("antibody", "immunoglobulin", "protein ", "oligonucleotide")
    if any(token in name for token in excluded_tokens):
        return "excluded", "biologic_or_polymer", "name indicates non-small-molecule", 0
    if smiles.count("C(=O)N") + smiles.count("C(N)=O") >= 15:
        return "review", "peptide_like", "high amide-bond count", 0
    return "accepted_core", "defined_chemical_entity", "defined structure", 1


def matching_accessions(row: dict[str, str], universe: set[str]) -> set[str]:
    matches: set[str] = set()
    for chain in range(1, 51):
        for prefix in ("UniProt (SwissProt)", "UniProt (TrEMBL)"):
            value = row.get(f"{prefix} Primary ID of Target Chain {chain}", "")
            for accession in split_ids(value):
                if accession in universe:
                    matches.add(accession)
    return matches


def evidence_id(parts: list[str]) -> str:
    digest = hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()[:20].upper()
    return f"BDBE_{digest}"


def main() -> None:
    csv.field_size_limit(min(sys.maxsize, 2_147_483_647))
    EVIDENCE_OUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    universe = read_master()

    evidence_header = [
        "evidence_id",
        "target_uniprot_id",
        "compound_id",
        "compound_id_type",
        "compound_source_id",
        "compound_name",
        "evidence_tier",
        "evidence_type",
        "activity_type",
        "activity_relation",
        "activity_value",
        "activity_unit",
        "pdb_ids",
        "ligand_het_id",
        "target_context",
        "target_assignment_status",
        "small_molecule_scope_status",
        "is_core_small_molecule",
        "pubmed_id",
        "doi",
        "source_record_id",
        "source_database",
        "source_version",
        "source_url",
        "record_qc_status",
        "default_release_inclusion",
    ]
    compounds: dict[str, dict[str, str | int]] = {}
    source_rows = 0
    matched_source_rows = 0
    evidence_rows = 0
    unique_proteins: set[str] = set()
    tier_counts: Counter[str] = Counter()
    assignment_counts: Counter[str] = Counter()

    with SOURCE.open("r", encoding="utf-8-sig", newline="") as source_handle, EVIDENCE_OUT.open(
        "w", encoding="utf-8", newline=""
    ) as evidence_handle:
        reader = csv.DictReader(source_handle, delimiter="\t")
        writer = csv.DictWriter(
            evidence_handle,
            fieldnames=evidence_header,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            source_rows += 1
            accessions = matching_accessions(row, universe)
            if not accessions:
                continue
            matched_source_rows += 1
            compound_id, compound_id_type, compound_source_id = compound_identity(row)
            scope_status, scope_class, scope_reason, is_core = compound_scope(row)
            compounds.setdefault(
                compound_id,
                {
                    "compound_id": compound_id,
                    "compound_id_type": compound_id_type,
                    "compound_source_id": compound_source_id,
                    "compound_name": (row.get("BindingDB Ligand Name", "") or "").strip(),
                    "canonical_smiles": (row.get("Ligand SMILES", "") or "").strip(),
                    "inchi": (row.get("Ligand InChI", "") or "").strip(),
                    "inchikey": (row.get("Ligand InChI Key", "") or "").strip(),
                    "pubchem_cids": (row.get("PubChem CID", "") or "").strip(),
                    "chembl_ids": (row.get("ChEMBL ID of Ligand", "") or "").strip(),
                    "chebi_ids": (row.get("ChEBI ID of Ligand", "") or "").strip(),
                    "drugbank_ids": (row.get("DrugBank ID of Ligand", "") or "").strip(),
                    "small_molecule_scope_status": scope_status,
                    "small_molecule_scope_class": scope_class,
                    "small_molecule_scope_reason": scope_reason,
                    "is_core_small_molecule": is_core,
                    "source_database": "BindingDB",
                    "source_version": "2026-07",
                },
            )

            try:
                chain_count = int(
                    (row.get("Number of Protein Chains in Target (>1 implies a multichain complex)", "") or "0").strip()
                    or "0"
                )
            except ValueError:
                chain_count = 0
            assignment = (
                "single_protein_direct"
                if chain_count == 1 and len(accessions) == 1
                else "complex_context_ambiguous"
            )
            assignment_counts[assignment] += 1
            target_context = (
                (row.get("Target Name", "") or "").strip()
                + " | chains="
                + str(chain_count)
            )
            source_record_id = (row.get("BindingDB Reactant_set_id", "") or "").strip()
            pubmed_id = (row.get("PMID", "") or "").strip()
            doi = (row.get("Article DOI", "") or "").strip()
            pdb_ids = (row.get("PDB ID(s) for Ligand-Target Complex", "") or "").strip()
            ligand_het = (row.get("Ligand HET ID in PDB", "") or "").strip()
            source_url = (
                row.get("Link to Ligand-Target Pair in BindingDB", "") or ""
            ).strip()
            default_ok = int(assignment == "single_protein_direct" and is_core == 1)

            measurements: list[tuple[str, str, str, str]] = []
            for field, (activity_type, tier) in MEASUREMENT_FIELDS.items():
                parsed = parse_value(row.get(field, ""))
                if parsed:
                    relation, value = parsed
                    measurements.append((activity_type, tier, relation, value))

            for accession in sorted(accessions):
                unique_proteins.add(accession)
                if pdb_ids and ligand_het:
                    record = {
                        "evidence_id": evidence_id(
                            [source_record_id, accession, compound_id, "STRUCTURE"]
                        ),
                        "target_uniprot_id": accession,
                        "compound_id": compound_id,
                        "compound_id_type": compound_id_type,
                        "compound_source_id": compound_source_id,
                        "compound_name": (row.get("BindingDB Ligand Name", "") or "").strip(),
                        "evidence_tier": "BE1",
                        "evidence_type": "compound_specific_structure",
                        "activity_type": "",
                        "activity_relation": "",
                        "activity_value": "",
                        "activity_unit": "",
                        "pdb_ids": pdb_ids,
                        "ligand_het_id": ligand_het,
                        "target_context": target_context,
                        "target_assignment_status": assignment,
                        "small_molecule_scope_status": scope_status,
                        "is_core_small_molecule": is_core,
                        "pubmed_id": pubmed_id,
                        "doi": doi,
                        "source_record_id": source_record_id,
                        "source_database": "BindingDB",
                        "source_version": "2026-07",
                        "source_url": source_url,
                        "record_qc_status": "ok",
                        "default_release_inclusion": default_ok,
                    }
                    writer.writerow(record)
                    evidence_rows += 1
                    tier_counts["BE1"] += 1

                for activity_type, tier, relation, value in measurements:
                    record = {
                        "evidence_id": evidence_id(
                            [
                                source_record_id,
                                accession,
                                compound_id,
                                activity_type,
                                relation,
                                value,
                            ]
                        ),
                        "target_uniprot_id": accession,
                        "compound_id": compound_id,
                        "compound_id_type": compound_id_type,
                        "compound_source_id": compound_source_id,
                        "compound_name": (row.get("BindingDB Ligand Name", "") or "").strip(),
                        "evidence_tier": tier,
                        "evidence_type": (
                            "direct_quantitative_binding"
                            if tier == "BE2"
                            else "binding_assay_activity"
                        ),
                        "activity_type": activity_type,
                        "activity_relation": relation,
                        "activity_value": value,
                        "activity_unit": "nM",
                        "pdb_ids": pdb_ids,
                        "ligand_het_id": ligand_het,
                        "target_context": target_context,
                        "target_assignment_status": assignment,
                        "small_molecule_scope_status": scope_status,
                        "is_core_small_molecule": is_core,
                        "pubmed_id": pubmed_id,
                        "doi": doi,
                        "source_record_id": source_record_id,
                        "source_database": "BindingDB",
                        "source_version": "2026-07",
                        "source_url": source_url,
                        "record_qc_status": "ok",
                        "default_release_inclusion": default_ok,
                    }
                    writer.writerow(record)
                    evidence_rows += 1
                    tier_counts[tier] += 1

    compound_header = [
        "compound_id",
        "compound_id_type",
        "compound_source_id",
        "compound_name",
        "canonical_smiles",
        "inchi",
        "inchikey",
        "pubchem_cids",
        "chembl_ids",
        "chebi_ids",
        "drugbank_ids",
        "small_molecule_scope_status",
        "small_molecule_scope_class",
        "small_molecule_scope_reason",
        "is_core_small_molecule",
        "source_database",
        "source_version",
    ]
    with COMPOUND_OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=compound_header, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for compound_id in sorted(compounds):
            writer.writerow(compounds[compound_id])

    report = {
        "status": "passed",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": "BindingDB curated from articles",
        "source_version": "2026-07",
        "source_rows": source_rows,
        "matched_source_rows": matched_source_rows,
        "evidence_rows": evidence_rows,
        "unique_membrane_proteins": len(unique_proteins),
        "unique_compounds": len(compounds),
        "evidence_tier_counts": dict(sorted(tier_counts.items())),
        "target_assignment_counts": dict(sorted(assignment_counts.items())),
        "default_inclusion_policy": (
            "single-protein target assignment and accepted-core small molecule"
        ),
        "evidence_output": str(EVIDENCE_OUT),
        "compound_output": str(COMPOUND_OUT),
    }
    REPORT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
