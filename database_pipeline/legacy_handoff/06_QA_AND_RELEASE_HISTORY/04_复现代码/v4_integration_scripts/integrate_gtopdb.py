from __future__ import annotations

import csv
import hashlib
import json
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "v4_working"
RAW = WORK / "raw" / "gtopdb_2026_2"
INTERMEDIATE = WORK / "intermediate"
REPORTS = WORK / "reports"
BASE_URL = "https://www.guidetopharmacology.org/DATA"
FILES = [
    "interactions.tsv",
    "ligands.tsv",
    "targets_and_families.tsv",
    "GtP_to_UniProt_mapping.tsv",
]
VERSION = "2026.2"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(name: str) -> Path:
    RAW.mkdir(parents=True, exist_ok=True)
    path = RAW / name
    if not path.exists():
        request = urllib.request.Request(
            f"{BASE_URL}/{name}",
            headers={"User-Agent": "mempro1-v4-build/1.0"},
        )
        with urllib.request.urlopen(request, timeout=180) as response:
            path.write_bytes(response.read())
    return path


def read_comment_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        lines = (line for line in handle if not line.startswith('"#'))
        return list(csv.DictReader(lines, delimiter="\t"))


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


def ligand_scope(ligand_type: str) -> tuple[str, str]:
    low = ligand_type.lower()
    if any(term in low for term in ("antibody", "protein", "peptide", "nucleic acid")):
        return "excluded", "biologic_or_polymer"
    if any(term in low for term in ("inorganic", "ion", "metal")):
        return "separate_not_core", "inorganic_ion_or_metal"
    return "accepted_core", "defined_chemical_entity"


def main() -> None:
    paths = {name: download(name) for name in FILES}
    interactions = read_comment_tsv(paths["interactions.tsv"])
    ligands = read_comment_tsv(paths["ligands.tsv"])
    targets = read_comment_tsv(paths["targets_and_families.tsv"])
    membrane = read_tsv(INTERMEDIATE / "human_membrane_protein_master.tsv")
    current_compounds = read_tsv(INTERMEDIATE / "small_molecule_master.tsv")
    membrane_ids = {row["target_uniprot_id"] for row in membrane}
    ligand_by_id = {row["Ligand ID"]: row for row in ligands}
    current_pubchem = {
        row["drug_id"]
        for row in current_compounds
        if row["compound_id_type"] == "PubChem CID"
    }
    current_inchikey = {row["inchikey"] for row in current_compounds if row["inchikey"]}

    accepted = []
    excluded = []
    additions: dict[str, dict] = {}
    target_classification = []
    for row in targets:
        accession = row["Human SwissProt"].strip()
        if accession and accession in membrane_ids:
            target_classification.append(
                {
                    "target_uniprot_id": accession,
                    "gtopdb_target_id": row["Target id"],
                    "gtopdb_target_name": row["Target name"],
                    "gtopdb_type": row["Type"],
                    "gtopdb_family_id": row["Family id"],
                    "gtopdb_family_name": row["Family name"],
                    "hgnc_id": row["HGNC id"],
                    "hgnc_symbol": row["HGNC symbol"],
                    "source_version": VERSION,
                }
            )

    for row in interactions:
        if row["Target Species"] != "Human":
            continue
        accession = row["Target UniProt ID"].strip()
        if accession not in membrane_ids:
            continue
        ligand = ligand_by_id.get(row["Ligand ID"], {})
        scope_status, scope_class = ligand_scope(row["Ligand Type"])
        pubchem_cid = ligand.get("PubChem CID", "")
        inchikey = ligand.get("InChIKey", "")
        original_type = row["Original Affinity Units"].strip().lower()
        has_quantitative = any(
            row[field].strip()
            for field in (
                "Original Affinity Low nm",
                "Original Affinity Median nm",
                "Original Affinity High nm",
            )
        )
        if scope_status == "accepted_core" and original_type in {"ki", "kd"} and has_quantitative:
            candidate_level = "SM4"
        else:
            candidate_level = "SM3"
        normalized = {
            "external_evidence_id": f"GTOPDB_{row['Target ID']}_{row['Ligand ID']}",
            "target_uniprot_id": accession,
            "target_gene_symbol": row["Target Gene Symbol"],
            "gtopdb_target_id": row["Target ID"],
            "ligand_id": row["Ligand ID"],
            "ligand_name": row["Ligand"],
            "ligand_type": row["Ligand Type"],
            "pubchem_cid": pubchem_cid,
            "inchikey": inchikey,
            "smiles": ligand.get("SMILES", ""),
            "small_molecule_scope_status": scope_status,
            "small_molecule_scope_class": scope_class,
            "approved": row["Approved"],
            "interaction_type": row["Type"],
            "action": row["Action"],
            "selectivity": row["Selectivity"],
            "endogenous": row["Endogenous"],
            "primary_target": row["Primary Target"],
            "affinity_units": row["Affinity Units"],
            "affinity_high": row["Affinity High"],
            "affinity_median": row["Affinity Median"],
            "affinity_low": row["Affinity Low"],
            "original_affinity_units": row["Original Affinity Units"],
            "original_affinity_low_nm": row["Original Affinity Low nm"],
            "original_affinity_median_nm": row["Original Affinity Median nm"],
            "original_affinity_high_nm": row["Original Affinity High nm"],
            "original_affinity_relation": row["Original Affinity Relation"],
            "assay_description": row["Assay Description"],
            "pubmed_id": row["PubMed ID"],
            "source_url": row["Webpage URLs"],
            "patent_numbers": row["Patent Numbers"],
            "candidate_sm_level": candidate_level,
            "source_database": "IUPHAR/BPS Guide to PHARMACOLOGY",
            "source_version": VERSION,
        }
        if scope_status == "accepted_core":
            accepted.append(normalized)
            compound_key = pubchem_cid or inchikey or f"GTOPDB:{row['Ligand ID']}"
            if pubchem_cid not in current_pubchem and inchikey not in current_inchikey:
                additions.setdefault(
                    compound_key,
                    {
                        "external_compound_id": f"GTOPDB:{row['Ligand ID']}",
                        "drug_id": pubchem_cid or f"GTOPDB:{row['Ligand ID']}",
                        "compound_id_type": "PubChem CID" if pubchem_cid else "GtoPdb ligand ID",
                        "drug_name": row["Ligand"],
                        "ligand_type": row["Ligand Type"],
                        "pubchem_cid": pubchem_cid,
                        "inchikey": inchikey,
                        "smiles": ligand.get("SMILES", ""),
                        "approved": row["Approved"],
                        "small_molecule_scope_status": "accepted_core",
                        "source_database": "IUPHAR/BPS Guide to PHARMACOLOGY",
                        "source_version": VERSION,
                    },
                )
        else:
            excluded.append(normalized)

    write_tsv(
        INTERMEDIATE / "external_gtopdb_interactions.tsv",
        accepted,
        list(accepted[0]),
    )
    write_tsv(
        INTERMEDIATE / "external_gtopdb_noncore_interactions.tsv",
        excluded,
        list(excluded[0]) if excluded else list(accepted[0]),
    )
    write_tsv(
        INTERMEDIATE / "external_gtopdb_compound_additions.tsv",
        sorted(additions.values(), key=lambda row: row["drug_id"]),
        list(next(iter(additions.values()))) if additions else ["drug_id"],
    )
    write_tsv(
        INTERMEDIATE / "gtopdb_target_classification.tsv",
        target_classification,
        list(target_classification[0]),
    )

    report = {
        "status": "passed",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": "IUPHAR/BPS Guide to PHARMACOLOGY",
        "source_version": VERSION,
        "download_files": {
            name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for name, path in paths.items()
        },
        "all_interaction_rows": len(interactions),
        "human_membrane_core_small_molecule_interactions": len(accepted),
        "human_membrane_noncore_ligand_interactions": len(excluded),
        "unique_membrane_targets": len({row["target_uniprot_id"] for row in accepted}),
        "new_compound_candidates": len(additions),
        "target_classification_rows": len(target_classification),
        "candidate_sm_level_counts": dict(
            Counter(row["candidate_sm_level"] for row in accepted)
        ),
        "ligand_type_counts": dict(Counter(row["ligand_type"] for row in accepted)),
    }
    (REPORTS / "GTOPDB_INTEGRATION_REPORT.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
