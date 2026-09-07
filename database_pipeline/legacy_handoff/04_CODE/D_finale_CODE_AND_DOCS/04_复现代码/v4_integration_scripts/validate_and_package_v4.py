from __future__ import annotations

import csv
import gzip
import hashlib
import json
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "v4_working"
INTERMEDIATE = WORK / "intermediate"
REPORTS = WORK / "reports"
RELEASE = WORK / "releases" / "release_v4"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy(source: str, target: str) -> Path:
    destination = RELEASE / target
    shutil.copy2(INTERMEDIATE / source, destination)
    return destination


def main() -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    RELEASE.mkdir(parents=True, exist_ok=True)
    master = read_tsv(INTERMEDIATE / "human_membrane_protein_master.tsv")
    status = read_tsv(INTERMEDIATE / "human_membrane_protein_small_molecule_status.tsv")
    compounds = read_tsv(INTERMEDIATE / "small_molecule_master_consolidated.tsv")
    evidence = read_tsv(INTERMEDIATE / "compound_protein_evidence_export.tsv")
    gtopdb = read_tsv(INTERMEDIATE / "external_gtopdb_interactions.tsv")
    assertions = read_tsv(INTERMEDIATE / "compound_protein_assertions.tsv")
    activities = read_tsv(INTERMEDIATE / "activity_measurements.tsv")
    structures = read_tsv(INTERMEDIATE / "structure_ligand_observations.tsv")
    sites = read_tsv(INTERMEDIATE / "binding_site_annotations.tsv")

    master_ids = [row["target_uniprot_id"] for row in master]
    status_ids = [row["target_uniprot_id"] for row in status]
    compound_ids = [row["drug_id"] for row in compounds]
    relation_ids = [row["relationship_context_id"] for row in evidence]
    require(len(master_ids) == len(set(master_ids)), "Duplicate protein in membrane master")
    require(len(status_ids) == len(set(status_ids)), "Duplicate protein in status table")
    require(set(master_ids) == set(status_ids), "Status and master protein sets differ")
    require(len(compound_ids) == len(set(compound_ids)), "Duplicate drug_id in compound master")
    require(len(relation_ids) == len(set(relation_ids)), "Duplicate relationship_context_id")
    require(all(row["organism_id"] == "9606" for row in master), "Non-human master row")
    require(
        all(row["small_molecule_evidence_level"] in {f"SM{i}" for i in range(6)} for row in status),
        "Invalid SM level",
    )
    require(
        all(
            row["has_small_molecule_evidence"] == "1"
            or row["small_molecule_evidence_level"] == "SM0"
            for row in status
        ),
        "Non-SM0 row lacks detailed-evidence flag",
    )
    require(
        all(row["target_uniprot_id"] in set(master_ids) for row in gtopdb),
        "GtoPdb target outside membrane master",
    )
    require(
        len(assertions) + len(activities) + len(structures) + len(sites) == len(evidence),
        "Evidence child tables do not partition v3.1 relationships",
    )

    baseline = json.loads((WORK / "baseline_manifest.json").read_text(encoding="utf-8"))
    for item in baseline["files"]:
        path = ROOT / item["path"]
        require(path.exists(), f"Missing baseline file: {item['path']}")
        require(sha256(path) == item["sha256"], f"Baseline changed: {item['path']}")

    validation = {
        "status": "passed",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "master_protein_rows": len(master),
        "status_rows": len(status),
        "compound_rows": len(compounds),
        "v3_1_evidence_rows": len(evidence),
        "gtopdb_external_evidence_rows": len(gtopdb),
        "activity_rows": len(activities),
        "assertion_rows": len(assertions),
        "structure_rows": len(structures),
        "binding_site_rows": len(sites),
        "protein_key_unique": True,
        "status_key_unique": True,
        "compound_key_unique": True,
        "relationship_key_unique": True,
        "status_foreign_key_closed": True,
        "gtopdb_target_foreign_key_closed": True,
        "evidence_partition_complete": True,
        "baseline_checksums_unchanged": True,
    }
    validation_path = REPORTS / "V4_VALIDATION_REPORT.json"
    validation_path.write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")

    level_counts = Counter(row["small_molecule_evidence_level"] for row in status)
    evidence_status_counts = Counter(row["small_molecule_evidence_status"] for row in status)
    scope_counts = Counter(row["membrane_scope"] for row in master)
    coverage = [
        "# v4 Completeness Report",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Reference universe",
        "",
        f"- Reviewed human membrane or membrane-associated proteins: {len(master):,}",
        *[f"- {key}: {value:,}" for key, value in sorted(scope_counts.items())],
        "",
        "## Small-molecule status",
        "",
        *[f"- {key}: {value:,}" for key, value in sorted(evidence_status_counts.items())],
        "",
        "## Final detailed-evidence levels",
        "",
        *[f"- {key}: {value:,}" for key, value in sorted(level_counts.items())],
        "",
        "## Interpretation",
        "",
        "Index-only candidates are not promoted to a final detailed-evidence level.",
        "SM0 means no qualifying detailed evidence in this release, not proof of no interaction.",
        "Peripheral membrane-associated proteins are reported separately from integral proteins.",
        "",
        "## Remaining completeness work",
        "",
        "- Import BindingDB measurement-level records for index-only candidates.",
        "- Import ChEMBL 37 activity-level records for index-only candidates.",
        "- Resolve 82 unreviewed transmembrane candidates without a reviewed gene-symbol match.",
        "- Review v3.1 proteins outside the reviewed UniProt membrane reference.",
        "- Add independent Human Protein Atlas, OPM/PDBTM, GPCRdb, and TCDB cross-checks.",
    ]
    completeness_path = REPORTS / "COMPLETENESS_REPORT.md"
    completeness_path.write_text("\n".join(coverage) + "\n", encoding="utf-8")

    packaged = [
        copy("human_membrane_protein_master.tsv", "human_membrane_protein_master_v4.tsv"),
        copy("human_membrane_protein_small_molecule_status.tsv", "human_membrane_protein_small_molecule_status_v4.tsv"),
        copy("small_molecule_master_consolidated.tsv", "small_molecule_master_v4.tsv"),
        copy("activity_measurements.tsv", "activity_measurements_v4.tsv"),
        copy("compound_protein_assertions.tsv", "compound_protein_assertions_v4.tsv"),
        copy("structure_ligand_observations.tsv", "structure_ligand_observations_v4.tsv"),
        copy("binding_site_annotations.tsv", "binding_site_annotations_v4.tsv"),
        copy("assays.tsv", "assays_v4.tsv"),
        copy("publications.tsv", "publications_v4.tsv"),
        copy("sources.tsv", "sources_v4.tsv"),
        copy("protein_classification.tsv", "protein_classification_v4.tsv"),
        copy("external_gtopdb_interactions.tsv", "external_gtopdb_interactions_v4.tsv"),
        copy("external_bindingdb_target_index.tsv", "external_bindingdb_target_index_v4.tsv"),
        copy("external_chembl_target_index.tsv", "external_chembl_target_index_v4.tsv"),
        copy("unreviewed_transmembrane_candidate_review.tsv", "unreviewed_transmembrane_candidate_review_v4.tsv"),
        copy("v3_1_proteins_not_in_reference.tsv", "v3_1_proteins_not_in_reference_v4.tsv"),
    ]
    evidence_gz = RELEASE / "compound_protein_evidence_v4.tsv.gz"
    with (INTERMEDIATE / "compound_protein_evidence_export.tsv").open("rb") as source:
        with gzip.GzipFile(filename=str(evidence_gz), mode="wb", compresslevel=9, mtime=0) as target:
            shutil.copyfileobj(source, target)
    packaged.append(evidence_gz)
    for source, target in [
        (WORK / "INCLUSION_POLICY.md", RELEASE / "INCLUSION_POLICY.md"),
        (WORK / "METHODS.md", RELEASE / "METHODS.md"),
        (WORK / "DATA_DICTIONARY.md", RELEASE / "DATA_DICTIONARY.md"),
        (completeness_path, RELEASE / "COMPLETENESS_REPORT.md"),
        (validation_path, RELEASE / "V4_VALIDATION_REPORT.json"),
    ]:
        shutil.copy2(source, target)
        packaged.append(target)

    manifest = {
        "release": "v4",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "validated_pre_excel",
        "files": [
            {
                "filename": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in sorted(packaged)
        ],
        "summary": {
            "membrane_proteins": len(master),
            "confirmed_detailed_evidence_proteins": sum(
                row["has_small_molecule_evidence"] == "1" for row in status
            ),
            "index_only_candidate_proteins": sum(
                row["small_molecule_evidence_status"]
                == "external_index_candidate_pending_measurement_import"
                for row in status
            ),
            "small_molecules": len(compounds),
            "v3_1_evidence_rows": len(evidence),
            "gtopdb_external_evidence_rows": len(gtopdb),
        },
    }
    manifest_path = RELEASE / "RELEASE_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    checksum_path = RELEASE / "checksums.sha256"
    checksum_path.write_text(
        "\n".join(
            f"{sha256(path)}  {path.name}"
            for path in sorted(RELEASE.iterdir())
            if path.is_file() and path.name != "checksums.sha256"
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"validation": validation, "manifest_summary": manifest["summary"]}, indent=2))


if __name__ == "__main__":
    main()
