"""
Step 3: Generate internal IDs, build new rows for master tables, and backfill.

Inputs:
  - pubchem_properties.jsonl (from Step 2)
  - cid_record_count.tsv (from Step 1)
  - cid_rows_map.jsonl (from Step 1)

Actions:
  3a. Read PubChem properties, generate HMPD-CMPD-* and HMPD-FORM-* IDs
  3b. Build new rows for small_molecule_master_v1_3.tsv
  3c. Build new rows for compound_form_hierarchy_v1_3.tsv
  3d. Append to master tables (with backup)
  3e. Build CID→(compound_internal_id, form_id) lookup for backfill
  3f. Backfill unmapped_review → upgraded negative evidence
"""
import csv
import gzip
import json
import os
import shutil
from datetime import datetime
from pathlib import Path

# === CONFIG ===
OUT_DIR = Path(r"D:\finale\negative_upgrade")
PROPS_FILE = OUT_DIR / "pubchem_properties.jsonl"
CID_COUNT_FILE = OUT_DIR / "cid_record_count.tsv"
CID_MAP_FILE = OUT_DIR / "cid_rows_map.jsonl"

MASTER_DIR = Path(r"D:\finale\01_正式数据_V6.2")
SMALL_MOL_MASTER = MASTER_DIR / "small_molecule_master_v1_3.tsv"
FORM_MASTER = MASTER_DIR / "compound_form_hierarchy_v1_3.tsv"
UNMAPPED_REVIEW = MASTER_DIR / "negative_binding_evidence_unmapped_review_v1_2.tsv.gz"
NEGATIVE_EVIDENCE = MASTER_DIR / "negative_binding_evidence_v1_2.tsv.gz"

# Outputs
NEW_COMPOUNDS_TSV = OUT_DIR / "new_compounds_v1_3_append.tsv"
NEW_FORMS_TSV = OUT_DIR / "new_forms_v1_3_append.tsv"
UPGRADED_NEGATIVE_TSV = OUT_DIR / "upgraded_negative_evidence.tsv.gz"
CID_LOOKUP_FILE = OUT_DIR / "cid_to_internal_ids.json"
UPGRADE_REPORT = OUT_DIR / "upgrade_report.json"

os.makedirs(OUT_DIR, exist_ok=True)

# ID ranges
COMPOUND_START = 2016065
FORM_START = 2020899
RELEASE_DATE = "2026-08-07"


def read_pubchem_properties():
    """Read all PubChem properties from JSONL."""
    props = {}
    with open(PROPS_FILE, encoding="utf-8") as f:
        for line in f:
            if line.strip():
                entry = json.loads(line)
                props[entry["cid"]] = entry["properties"]
    return props


def generate_smiles(inchi):
    """Basic SMILES generation from InChI is not trivial.
    We use the PubChem SMILES directly.
    If missing, leave empty."""
    return ""  # PubChem provides SMILES directly


def compute_inchikey_sections(inchikey):
    """Split InChIKey into connectivity and identity parts."""
    if inchikey and len(inchikey) >= 14:
        return inchikey[:14], f"IK:{inchikey}"
    return "", ""


def build_compound_row(cid, props, cpd_num, record_count):
    """Build a single row for small_molecule_master."""
    cpd_id = f"HMPD-CMPD-{cpd_num:07d}"

    smiles = props.get("SMILES", "") or ""
    inchi = props.get("InChI", "") or ""
    inchikey = props.get("InChIKey", "") or ""
    conn_key = inchikey[:14] if inchikey else ""
    identity_key = f"IK:{inchikey}" if inchikey else ""
    mf = props.get("MolecularFormula", "") or ""
    mw = props.get("MolecularWeight", "") or ""
    exact_mass = props.get("ExactMass", "") or ""

    # Convert numeric fields safely
    def safe_float(v):
        try:
            return str(float(v)) if v else ""
        except (ValueError, TypeError):
            return ""

    def safe_int(v):
        try:
            return str(int(float(v))) if v else ""
        except (ValueError, TypeError):
            return ""

    row = {
        "compound_internal_id": cpd_id,
        "preferred_name": props.get("IUPACName", "") or "",
        "preferred_name_source": "PubChem CID",
        "compound_scope_status": "core",
        "compound_scope_class": "core_small_molecule",
        "identity_confidence": "low",  # PubChem-only, no cross-validation
        "standard_smiles": smiles,
        "standard_inchi": inchi,
        "standard_inchikey": inchikey,
        "identity_key": identity_key,
        "connectivity_key": conn_key,
        "molecular_formula": mf,
        "molecular_weight": mw,
        "exact_mass": exact_mass,
        "xlogp": safe_float(props.get("XLogP", "")),
        "tpsa": safe_float(props.get("TPSA", "")),
        "hbond_donor_count": safe_int(props.get("HBondDonorCount", "")),
        "hbond_acceptor_count": safe_int(props.get("HBondAcceptorCount", "")),
        "rotatable_bond_count": safe_int(props.get("RotatableBondCount", "")),
        "formal_charge": safe_int(props.get("Charge", "")),
        "heavy_atom_count": safe_int(props.get("HeavyAtomCount", "")),
        "ring_count": "",  # would need RDKit
        "max_ring_size": "",
        "amide_bond_count": "",
        "computed_structural_class": "",  # would need RDKit
        "molecule_types": "",
        "form_count": "1",
        "source_record_count": str(record_count),
        "source_database_count": "1",
        "source_databases": "PubChem BioAssay",
        "pubchem_cids": cid,
        "chembl_ids": "",
        "chebi_ids": "",
        "drugcentral_ids": "",
        "gtopdb_ligand_ids": "",
        "hmdb_ids": "",
        "cas_numbers": "",
        "synonym_count": "0",
        "development_status": "",
        "max_chembl_phase": "",
        "first_approval_year": "",
        "is_approved_drug": "0",
        "is_clinical_candidate": "0",
        "is_endogenous_ligand": "0",
        "is_natural_product": "0",
        "is_chemical_probe": "0",
        "compound_category_tags": "negative_evidence_compound",
        "chebi_direct_classes": "",
        "chebi_roles": "",
        "record_qc_status": "ok",
        "qc_notes": "auto_added_from_unmapped_negative_review",
        "protein_target_count": "0",
        "best_binding_evidence_level": "",
        "BE1_evidence_count": "0",
        "BE2_evidence_count": "0",
        "BE3_evidence_count": "0",
        "binding_evidence_count": "0",
        "binding_source_count": "0",
        "binding_sources": "",
        "best_standard_value_nM": "",
        "release_version": "v1.3",
        "release_date": RELEASE_DATE,
        "standardization_policy": "pubchem_properties_only",
        "incremental_evidence_run_v11": "",
        "identity_refresh_v12": "",
        "new_source_record_count_v12": "",
        "release_version_v12": "",
        "release_date_v12": "",
        "negative_evidence_count_v13": str(record_count),
        "negative_target_count_v13": "",
        "positive_negative_conflict_pair_count_v13": "0",
        "negative_identity_refresh_v13": "unmapped_upgrade_20260807",
        "release_version_v13": "v1.3",
        "release_date_v13": RELEASE_DATE,
    }
    return row


def build_form_row(cpd_id, form_num, props, cid):
    """Build a single parent_form row for compound_form_hierarchy."""
    form_id = f"HMPD-FORM-{form_num:07d}"

    smiles = props.get("SMILES", "") or ""
    inchi = props.get("InChI", "") or ""
    inchikey = props.get("InChIKey", "") or ""
    identity_key = f"IK:{inchikey}" if inchikey else ""
    mf = props.get("MolecularFormula", "") or ""
    mw = props.get("MolecularWeight", "") or ""

    def safe_int(v):
        try:
            return str(int(float(v))) if v else ""
        except (ValueError, TypeError):
            return ""

    return {
        "compound_form_id": form_id,
        "compound_internal_id": cpd_id,
        "form_type": "parent_form",
        "exact_smiles": smiles,
        "exact_inchi": inchi,
        "exact_inchikey": inchikey,
        "exact_identity_key": identity_key,
        "molecular_formula": mf,
        "molecular_weight": mw,
        "formal_charge": safe_int(props.get("Charge", "")),
        "source_record_count": "1",
        "source_databases": "PubChem BioAssay",
        "source_compound_ids": f"PubChem CID:{cid}",
        "chembl_parent_ids": "",
        "chembl_parent_alignment_status": "",
        "record_qc_status": "ok",
        "qc_notes": "auto_added_from_unmapped_negative_review",
        "release_version_v11": "",
        "parent_inchikey_v12": inchikey,
        "parent_relation_v12": "self_parent",
    }


def load_cid_record_counts():
    """Load CID → record count mapping."""
    counts = {}
    with open(CID_COUNT_FILE, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            counts[row["cid"]] = int(row["record_count"])
    return counts


def main():
    print("=" * 60)
    print("Step 3: Build new master rows and backfill negative evidence")
    print("=" * 60)

    # 3a: Load all data
    print("\n--- Loading PubChem properties ---")
    props = read_pubchem_properties()
    print(f"  Loaded {len(props):,} compounds with properties")

    print("\n--- Loading CID record counts ---")
    cid_counts = load_cid_record_counts()
    print(f"  Loaded {len(cid_counts):,} CID→count mappings")

    # Get ordered list of CIDs
    cids_ordered = sorted(props.keys(), key=lambda x: int(x))
    print(f"  Processing {len(cids_ordered):,} CIDs")

    # 3b: Generate rows
    print("\n--- Generating compound and form rows ---")
    compound_rows = []
    form_rows = []
    cid_to_internal = {}  # cid -> {compound_internal_id, compound_form_id}

    cpd_num = COMPOUND_START
    form_num = FORM_START

    for cid in cids_ordered:
        p = props[cid]
        record_count = cid_counts.get(cid, 1)

        cpd_row = build_compound_row(cid, p, cpd_num, record_count)
        form_row = build_form_row(cpd_row["compound_internal_id"], form_num, p, cid)

        compound_rows.append(cpd_row)
        form_rows.append(form_row)

        cid_to_internal[cid] = {
            "compound_internal_id": cpd_row["compound_internal_id"],
            "compound_form_id": form_row["compound_form_id"],
        }

        cpd_num += 1
        form_num += 1

        if len(compound_rows) % 50000 == 0:
            print(f"  Generated {len(compound_rows):,} rows...")

    print(f"  Total compound rows: {len(compound_rows):,}")
    print(f"  Total form rows: {len(form_rows):,}")
    print(f"  Last compound ID: {compound_rows[-1]['compound_internal_id']}")
    print(f"  Last form ID: {form_rows[-1]['compound_form_id']}")

    # 3c: Save new compound rows
    print("\n--- Saving new compound rows ---")
    comp_fieldnames = list(compound_rows[0].keys())
    with open(NEW_COMPOUNDS_TSV, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=comp_fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(compound_rows)
    print(f"  Saved to: {NEW_COMPOUNDS_TSV}")

    # 3d: Save new form rows
    print("\n--- Saving new form rows ---")
    form_fieldnames = list(form_rows[0].keys())
    with open(NEW_FORMS_TSV, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=form_fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(form_rows)
    print(f"  Saved to: {NEW_FORMS_TSV}")

    # 3e: Save CID→internal ID lookup
    with open(CID_LOOKUP_FILE, "w", encoding="utf-8") as f:
        json.dump(cid_to_internal, f, ensure_ascii=False)
    print(f"\n  Saved CID→internal ID lookup: {CID_LOOKUP_FILE}")

    # 3f: Backfill unmapped_review records
    print("\n--- Backfilling unmapped review records ---")
    backfill_count = 0
    cids_not_in_props = set()

    # First pass: count how many records will be upgraded
    total_upgraded = 0
    with gzip.open(UNMAPPED_REVIEW, "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cid = row.get("pubchem_cid", "").strip()
            if cid in cid_to_internal:
                total_upgraded += 1

    print(f"  Records to upgrade: {total_upgraded:,}")

    # Second pass: write upgraded records with backfilled IDs
    upgraded_fieldnames = None
    upgraded_rows = []

    with gzip.open(UNMAPPED_REVIEW, "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        upgraded_fieldnames = list(reader.fieldnames)
        # Update the relevant columns
        for row in reader:
            cid = row.get("pubchem_cid", "").strip()
            if cid in cid_to_internal:
                ids = cid_to_internal[cid]
                row["compound_internal_id"] = ids["compound_internal_id"]
                row["compound_form_id"] = ids["compound_form_id"]
                row["compound_mapping_status"] = "mapped_core"
                row["release_action_v62"] = "mapped_negative_release"
                row["unmapped_reason_v62"] = ""  # Clear the unmapped reason
                row["release_version_v62"] = "v6.2"
                upgraded_rows.append(row)
                backfill_count += 1
            else:
                cids_not_in_props.add(cid)

    print(f"  Backfilled: {backfill_count:,} records")
    print(f"  CIDs not in properties (could not backfill): {len(cids_not_in_props):,}")

    # Save upgraded negative evidence
    if upgraded_rows:
        with gzip.open(UPGRADED_NEGATIVE_TSV, "wt", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=upgraded_fieldnames, delimiter="\t",
                                    extrasaction="ignore")
            writer.writeheader()
            writer.writerows(upgraded_rows)
        print(f"  Saved upgraded evidence: {UPGRADED_NEGATIVE_TSV}")
        print(f"  {backfill_count:,} records ready to merge into formal negative evidence")

    # 3g: Report
    report = {
        "date": RELEASE_DATE,
        "new_compounds": len(compound_rows),
        "new_forms": len(form_rows),
        "compound_id_range": f"HMPD-CMPD-{COMPOUND_START:07d} to {compound_rows[-1]['compound_internal_id']}",
        "form_id_range": f"HMPD-FORM-{FORM_START:07d} to {form_rows[-1]['compound_form_id']}",
        "records_upgraded": backfill_count,
        "cids_queried": len(props),
        "cids_backfilled": len(cid_to_internal),
        "cids_still_unmapped": len(cids_not_in_props),
        "appended_to_small_molecule_master": str(NEW_COMPOUNDS_TSV),
        "appended_to_compound_form_hierarchy": str(NEW_FORMS_TSV),
        "upgraded_negative_evidence": str(UPGRADED_NEGATIVE_TSV),
    }

    with open(UPGRADE_REPORT, "w") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    print(f"\n  Report saved: {UPGRADE_REPORT}")

    print("\n" + "=" * 60)
    print("Next steps (manual):")
    print("=" * 60)
    print(f"1. Verify new compound rows: {NEW_COMPOUNDS_TSV}")
    print(f"2. Verify new form rows: {NEW_FORMS_TSV}")
    print(f"3. Backup original files to {OUT_DIR}/backup/")
    print(f"4. Append new rows to small_molecule_master_v1_3.tsv")
    print(f"5. Append new rows to compound_form_hierarchy_v1_3.tsv")
    print(f"6. Merge upgraded negative evidence into negative_binding_evidence_v1_2.tsv.gz")
    print(f"7. Trim unmapped_review to remove upgraded records")
    print("\nDone.")


if __name__ == "__main__":
    main()
