from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "v4_working"
INTERMEDIATE = WORK / "intermediate"
REPORTS = WORK / "reports"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if fields is None:
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


def compound_key(drug_id: str, pubchem_cid: str = "", fallback: str = "") -> str:
    if pubchem_cid:
        return f"CID:{pubchem_cid}"
    if drug_id and drug_id.isdigit():
        return f"CID:{drug_id}"
    return drug_id or fallback


def main() -> None:
    master = read_tsv(INTERMEDIATE / "human_membrane_protein_master.tsv")
    current = read_tsv(INTERMEDIATE / "compound_protein_evidence_export.tsv")
    gtopdb = read_tsv(INTERMEDIATE / "external_gtopdb_interactions.tsv")
    bindingdb = read_tsv(INTERMEDIATE / "external_bindingdb_target_index.tsv")
    chembl = read_tsv(INTERMEDIATE / "external_chembl_target_index.tsv")
    compounds = read_tsv(INTERMEDIATE / "small_molecule_master.tsv")
    gtopdb_additions = read_tsv(INTERMEDIATE / "external_gtopdb_compound_additions.tsv")

    bindingdb_ids = {row["target_uniprot_id"] for row in bindingdb}
    chembl_ids = {row["target_uniprot_id"] for row in chembl}
    evidence_by_protein: dict[str, list[dict]] = defaultdict(list)
    direct_pair_sources: dict[tuple[str, str], set[str]] = defaultdict(set)
    compound_sets: dict[str, set[str]] = defaultdict(set)
    source_sets: dict[str, set[str]] = defaultdict(set)
    current_proteins = set()
    gtopdb_proteins = set()

    for row in current:
        if row["qualifies_for_core_small_molecule_flag"] != "1":
            continue
        accession = row["target_uniprot_id"]
        current_proteins.add(accession)
        level = row["candidate_sm_level"]
        key = compound_key(row["drug_id"])
        sources = [x.strip() for x in row["source_database"].split(";") if x.strip()]
        compound_sets[accession].add(key)
        source_sets[accession].update(sources)
        evidence_by_protein[accession].append(
            {
                "level": level,
                "rank": int(row["candidate_sm_level_rank"]),
                "source": row["source_database"],
                "compound_key": key,
                "evidence_id": row["relationship_context_id"],
            }
        )
        if level == "SM4":
            direct_pair_sources[(accession, key)].update(sources)

    for row in gtopdb:
        accession = row["target_uniprot_id"]
        gtopdb_proteins.add(accession)
        key = compound_key(
            row.get("drug_id", ""),
            row["pubchem_cid"],
            f"GTOPDB:{row['ligand_id']}",
        )
        level = row["candidate_sm_level"]
        rank = int(level.replace("SM", ""))
        compound_sets[accession].add(key)
        source_sets[accession].add("IUPHAR/BPS Guide to PHARMACOLOGY")
        evidence_by_protein[accession].append(
            {
                "level": level,
                "rank": rank,
                "source": "IUPHAR/BPS Guide to PHARMACOLOGY",
                "compound_key": key,
                "evidence_id": row["external_evidence_id"],
            }
        )
        if level == "SM4":
            direct_pair_sources[(accession, key)].add(
                "IUPHAR/BPS Guide to PHARMACOLOGY"
            )

    sm5_proteins = {
        accession
        for (accession, _), sources in direct_pair_sources.items()
        if len(sources) >= 2
    }
    status_rows = []
    for row in master:
        accession = row["target_uniprot_id"]
        detailed = evidence_by_protein.get(accession, [])
        has_detailed = bool(detailed)
        has_bindingdb = accession in bindingdb_ids
        has_chembl = accession in chembl_ids
        has_index = has_bindingdb or has_chembl
        if accession in sm5_proteins:
            level = "SM5"
            rank = 5
            level_reason = "same protein-compound direct evidence supported by multiple independent sources"
        elif detailed:
            best = max(detailed, key=lambda item: item["rank"])
            level = best["level"]
            rank = best["rank"]
            level_reason = "highest qualifying detailed evidence"
        else:
            level = "SM0"
            rank = 0
            level_reason = "no qualifying detailed small-molecule evidence in integrated release"

        if has_detailed:
            evidence_status = "confirmed_detailed_evidence"
        elif has_index:
            evidence_status = "external_index_candidate_pending_measurement_import"
        else:
            evidence_status = "current_sources_no_evidence_found"
        level_counts = Counter(item["level"] for item in detailed)
        status_rows.append(
            {
                "membrane_protein_id": row["membrane_protein_id"],
                "target_uniprot_id": accession,
                "approved_symbol": row["approved_symbol"],
                "protein_name": row["protein_name"],
                "reviewed": row["reviewed"],
                "membrane_scope": row["membrane_scope"],
                "membrane_topology": row["membrane_topology"],
                "transmembrane_count": row["transmembrane_count"],
                "functional_primary_class": row["functional_primary_class"],
                "functional_subclass": row["functional_subclass"],
                "has_small_molecule_evidence": "1" if has_detailed else "0",
                "has_small_molecule_source_index_hit": "1" if has_index else "0",
                "small_molecule_evidence_status": evidence_status,
                "small_molecule_evidence_level": level,
                "small_molecule_evidence_level_rank": rank,
                "small_molecule_evidence_level_reason": level_reason,
                "has_direct_quantitative_binding": (
                    "1" if any(item["rank"] >= 4 for item in detailed) else "0"
                ),
                "has_multi_source_validated_direct_evidence": (
                    "1" if accession in sm5_proteins else "0"
                ),
                "has_functional_activity": (
                    "1" if level_counts["SM2"] else "0"
                ),
                "has_curated_target_relation": (
                    "1" if level_counts["SM3"] else "0"
                ),
                "has_structure_or_site_evidence": (
                    "1" if level_counts["SM1"] else "0"
                ),
                "present_in_current_v3_1_detailed_evidence": (
                    "1" if accession in current_proteins else "0"
                ),
                "present_in_gtopdb_detailed_evidence": (
                    "1" if accession in gtopdb_proteins else "0"
                ),
                "present_in_bindingdb_target_index": "1" if has_bindingdb else "0",
                "present_in_chembl_target_index": "1" if has_chembl else "0",
                "small_molecule_count": len(compound_sets[accession]),
                "detailed_evidence_record_count": len(detailed),
                "independent_source_count": len(source_sets[accession]),
                "independent_sources": ";".join(sorted(source_sets[accession])),
                "sm1_record_count": level_counts["SM1"],
                "sm2_record_count": level_counts["SM2"],
                "sm3_record_count": level_counts["SM3"],
                "sm4_record_count": level_counts["SM4"],
                "status_policy_version": "membrane-small-molecule-scope-v1",
            }
        )
    write_tsv(
        INTERMEDIATE / "human_membrane_protein_small_molecule_status.tsv",
        status_rows,
    )

    current_keys = {
        (
            row["drug_id"]
            if row["compound_id_type"] == "PubChem CID"
            else row["inchikey"] or row["drug_id"]
        )
        for row in compounds
    }
    consolidated_compounds = list(compounds)
    additions_added = 0
    for row in gtopdb_additions:
        key = row["pubchem_cid"] or row["inchikey"] or row["drug_id"]
        if key in current_keys:
            continue
        consolidated_compounds.append(
            {
                "drug_id": row["drug_id"],
                "compound_id_type": row["compound_id_type"],
                "drug_name": row["drug_name"],
                "compound_name_category": "gtopdb_preferred_name",
                "compound_biological_role": "bioactive_research_ligand",
                "ligand_interpretation_classes": "curated_target_ligand",
                "best_ligand_interpretation_class": "curated_target_ligand",
                "compound_source": "curated_target_relation",
                "source_databases": row["source_database"],
                "source_row_count": "",
                "unique_target_count": "",
                "compound_evidence_status": "gtopdb_curated_interaction",
                "activity_types": "",
                "activity_value_uM_min": "",
                "activity_value_uM_median": "",
                "molecular_formula": "",
                "molecular_weight": "",
                "canonical_smiles": row["smiles"],
                "isomeric_smiles": "",
                "inchikey": row["inchikey"],
                "iupac_name": "",
                "xlogp": "",
                "tpsa": "",
                "hbond_donor_count": "",
                "hbond_acceptor_count": "",
                "rotatable_bond_count": "",
                "complexity": "",
                "formal_charge": "",
                "ghs_hazard_classes": "",
                "ghs_signal_words": "",
                "toxicity_summary": "",
                "livertox": "",
                "drug_classes": "",
                "pharmacodynamics": "",
                "drug_indication": "",
                "record_qc_status": "external_addition_requires_property_enrichment",
                "record_qc_note": f"GtoPdb {row['source_version']}",
                "small_molecule_scope_status": "accepted_core",
                "small_molecule_scope_class": "defined_chemical_entity",
                "small_molecule_scope_reason": "IUPHAR/BPS curated small-molecule interaction",
                "is_core_small_molecule": "1",
                "v4_scope_policy_version": "membrane-small-molecule-scope-v1",
            }
        )
        current_keys.add(key)
        additions_added += 1
    write_tsv(
        INTERMEDIATE / "small_molecule_master_consolidated.tsv",
        consolidated_compounds,
        list(compounds[0]),
    )

    evidence_status_counts = Counter(
        row["small_molecule_evidence_status"] for row in status_rows
    )
    level_counts = Counter(row["small_molecule_evidence_level"] for row in status_rows)
    report = {
        "status": "passed",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "membrane_protein_rows": len(status_rows),
        "confirmed_detailed_evidence_proteins": sum(
            row["has_small_molecule_evidence"] == "1" for row in status_rows
        ),
        "external_index_only_candidate_proteins": sum(
            row["small_molecule_evidence_status"]
            == "external_index_candidate_pending_measurement_import"
            for row in status_rows
        ),
        "no_evidence_or_index_hit_proteins": sum(
            row["small_molecule_evidence_status"] == "current_sources_no_evidence_found"
            for row in status_rows
        ),
        "evidence_status_counts": dict(sorted(evidence_status_counts.items())),
        "protein_sm_level_counts": dict(sorted(level_counts.items())),
        "sm5_protein_count": len(sm5_proteins),
        "current_v3_1_detailed_proteins": len(current_proteins),
        "gtopdb_detailed_proteins": len(gtopdb_proteins),
        "gtopdb_new_detailed_proteins": len(gtopdb_proteins - current_proteins),
        "bindingdb_index_proteins": len(bindingdb_ids),
        "chembl_index_proteins": len(chembl_ids),
        "consolidated_compound_rows": len(consolidated_compounds),
        "gtopdb_compound_additions": additions_added,
    }
    (REPORTS / "PROTEIN_SMALL_MOLECULE_STATUS_REPORT.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
