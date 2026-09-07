from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


BASE = Path(r"D:\finale\20_MemPro_V7.1.1_final_20260814")
SRC = BASE / "01_release_tables"
XWALK = Path(r"D:\7.22\v71_publication_candidate_working\work\v3_1_activity_source_crosswalk_v71.tsv.gz")
OUT = Path(r"D:\finale\21_MemPro_V7.2_candidate_20260817")
TABLES = OUT / "01_release_tables"
META = OUT / "02_metadata"
QA = OUT / "03_QA"
for d in (TABLES, META, QA):
    d.mkdir(parents=True, exist_ok=True)

VERSION = "MemPro_V7.2-candidate"
NOW = datetime.now(timezone.utc).isoformat()


def write(df: pd.DataFrame, name: str) -> Path:
    p = TABLES / name
    df.to_csv(p, sep="\t", index=False, compression="gzip" if name.endswith(".gz") else None)
    return p


def v72_name(name: str) -> str:
    return name.replace("_v71", "_v72")


def source_split(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(";") if x.strip() and x.strip() != "nan"]


# 1. Frozen baseline and rules.
rules = pd.DataFrame([
    ("protein_release", "E1_or_E2", "CONFIRMED_CORE", "default website and manuscript core"),
    ("protein_release", "E3", "CANDIDATE_E3", "searchable and downloadable; excluded from default core statistics"),
    ("form_specificity", "source_does_not_identify_isoform", "CANONICAL_PLUS_ISOFORM_UNSPECIFIED", "retain; never delete"),
    ("complex", "all_complex_records", "UNIFIED_COMPLEX_MODULE", "no formal/beta/hidden release split; retain evidence-state fields"),
    ("site", "exact_cocrystal_and_complete_coordinate", "S1", "query ligand, structure, chain and coordinates resolved"),
    ("site", "coordinate_complete_ligand_unresolved", "S2", "coordinate-usable, not an exact query-ligand cocrystal claim"),
    ("site", "partial_or_source_only", "S3", "source site retained with explicit coordinate limitation"),
    ("interaction", "no_reported_residue", "NO_REPORTED_RESIDUE", "valid interaction evidence; residue not required"),
    ("negative", "nonconflicting_negative", "ARCHIVE_ONLY", "not part of positive interaction core"),
], columns=["module", "condition", "v72_status", "interpretation"])
rules.to_csv(META / "V72_RELEASE_RULES.tsv", sep="\t", index=False)

baseline = []
for p in sorted(SRC.iterdir()):
    if p.is_file():
        baseline.append((p.name, p.stat().st_size, hashlib.sha256(p.read_bytes()).hexdigest()))
pd.DataFrame(baseline, columns=["file", "bytes", "sha256"]).to_csv(META / "V711_BASELINE_MANIFEST.tsv", sep="\t", index=False)

# 2. Protein confirmed/candidate release roles.
cls = pd.read_csv(SRC / "canonical_protein_membrane_class_v71.tsv.gz", sep="\t", low_memory=False)
cls["release_role_v72"] = cls.evidence_level.map(lambda x: "CONFIRMED_CORE" if x in {"E1", "E2"} else "CANDIDATE_E3")
cls["default_web_inclusion_v72"] = cls.evidence_level.isin(["E1", "E2"]).astype(int)
cls["searchable_inclusion_v72"] = 1
cls["release_version_v72"] = VERSION
write(cls, "canonical_protein_membrane_class_v72.tsv.gz")
prot_role = cls.set_index("canonical_uniprot_accession")["release_role_v72"]
prot_default = cls.set_index("canonical_uniprot_accession")["default_web_inclusion_v72"]

pm = pd.read_csv(SRC / "protein_master_v71.tsv.gz", sep="\t", low_memory=False)
pm["release_role_v72"] = pm.canonical_uniprot_accession.map(prot_role)
pm["default_web_inclusion_v72"] = pm.canonical_uniprot_accession.map(prot_default).fillna(0).astype(int)
pm["searchable_inclusion_v72"] = pm.canonical_uniprot_accession.isin(prot_role.index).astype(int)
pm["release_version_v72"] = VERSION
write(pm, "protein_master_v72.tsv.gz")

# 3. Recover legacy v3.1 assay provenance and rebuild positive evidence.
xw = pd.read_csv(XWALK, sep="\t", dtype=str, low_memory=False)
xw = xw.drop_duplicates("activity_measurement_id").set_index("activity_measurement_id")
ev = pd.read_csv(SRC / "positive_interaction_evidence_v71.tsv.gz", sep="\t", low_memory=False)
legacy = ev.source_database.eq("v3.1_mixed_sources")
mapped = ev.loc[legacy, "source_record_id"].map(xw.source_database)
assay_id = ev.loc[legacy, "source_record_id"].map(xw.assay_id)
source_assay = ev.loc[legacy, "source_record_id"].map(xw.source_assay_id)
ev["legacy_source_database_v72"] = ev.source_database.where(legacy)
ev["legacy_source_record_id_v72"] = ev.source_record_id.where(legacy)
ev.loc[legacy, "source_database"] = mapped.fillna("LEGACY_SOURCE_UNRESOLVED")
resolved_record = source_assay.where(source_assay.notna() & source_assay.ne(""), assay_id)
ev.loc[legacy, "source_record_id"] = resolved_record.values
ev["source_provenance_resolution_v72"] = "SOURCE_ALREADY_SPECIFIC"
ev.loc[legacy & mapped.notna(), "source_provenance_resolution_v72"] = "RECOVERED_HISTORICAL_ASSAY_REGISTRY"
ev.loc[legacy & mapped.isna(), "source_provenance_resolution_v72"] = "LEGACY_SOURCE_UNRESOLVED"
ev["protein_release_role_v72"] = ev.target_uniprot_id.map(prot_role)
ev["default_web_inclusion_v72"] = ev.target_uniprot_id.map(prot_default).fillna(0).astype(int)
ev["searchable_inclusion_v72"] = ev.target_uniprot_id.isin(prot_role.index).astype(int)
ev["protein_form_specificity_v72"] = "ISOFORM_UNSPECIFIED"
ev["target_assignment_policy_v72"] = "DIRECT_TO_CANONICAL_RETAINED"
has_residue = ev.binding_site_residues.fillna("").astype(str).str.strip().ne("")
ev["site_availability_v72"] = has_residue.map({True: "SOURCE_REPORTED_RESIDUE", False: "NO_REPORTED_RESIDUE"})
ev["release_version_v72"] = VERSION
write(ev, "positive_interaction_evidence_v72.tsv.gz")

# 4. Canonical/form assignment policy: retain every evidence record.
assign = pd.read_csv(SRC / "evidence_target_entity_assignment_v71.tsv.gz", sep="\t", low_memory=False)
assign["legacy_source_database_v72"] = assign.source_database.where(assign.source_database.eq("v3.1_mixed_sources"))
assign["legacy_source_record_id_v72"] = assign.source_record_id.where(assign.source_database.eq("v3.1_mixed_sources"))
assign["source_database"] = assign.evidence_id.map(ev.set_index("evidence_id").source_database).fillna(assign.source_database)
assign["source_record_id"] = assign.evidence_id.map(ev.set_index("evidence_id").source_record_id).fillna(assign.source_record_id)
unspec = assign.protein_form_specificity_v71.eq("ISOFORM_UNSPECIFIED") | assign.protein_isoform_entity_id.isna()
assign["target_entity_type_v72"] = assign.target_entity_type_v71
assign["target_entity_id_v72"] = assign.target_entity_id_v71
assign["protein_form_specificity_v72"] = assign.protein_form_specificity_v71
assign.loc[unspec, "target_entity_type_v72"] = "canonical_protein"
assign.loc[unspec, "target_entity_id_v72"] = assign.loc[unspec, "canonical_uniprot_accession"]
assign.loc[unspec, "protein_form_specificity_v72"] = "ISOFORM_UNSPECIFIED"
assign["assignment_policy_v72"] = "EXPLICIT_FORM_WHEN_SOURCE_SPECIFIC_ELSE_CANONICAL"
assign["record_retained_v72"] = 1
assign["release_version_v72"] = VERSION
write(assign, "evidence_target_entity_assignment_v72.tsv.gz")

iso = pd.read_csv(SRC / "isoform_membrane_class_v71.tsv.gz", sep="\t", low_memory=False)
iso["canonical_interaction_fallback_v72"] = "DIRECT_TO_CANONICAL_WHEN_ISOFORM_UNSPECIFIED"
iso["record_retained_v72"] = 1
iso["release_version_v72"] = VERSION
write(iso, "isoform_membrane_class_v72.tsv.gz")

pf = pd.read_csv(SRC / "processed_form_membrane_class_v71.tsv.gz", sep="\t", low_memory=False)
pf["canonical_interaction_fallback_v72"] = "DIRECT_TO_CANONICAL_WHEN_PROCESSED_FORM_UNSPECIFIED"
pf["record_retained_v72"] = 1
pf["release_version_v72"] = VERSION
write(pf, "processed_form_membrane_class_v72.tsv.gz")

# 5. Site semantic tiers.
site = pd.read_csv(SRC / "interaction_site_v71.tsv.gz", sep="\t", low_memory=False)
complete = site.coordinate_mapping_status.eq("COORDINATE_COMPLETE") & site.chain_assignment_status.isin(["UNIQUE_CHAIN", "MULTIPLE_EQUIVALENT_CHAINS"])
exact = site.structure_ligand_relationship.eq("EXACT_COCRYSTAL_LIGAND")
site["site_tier_v72"] = "S3_SOURCE_OR_PARTIAL"
site.loc[complete, "site_tier_v72"] = "S2_COORDINATE_COMPLETE_LIGAND_UNRESOLVED"
site.loc[complete & exact, "site_tier_v72"] = "S1_EXACT_COCRYSTAL_COORDINATE_COMPLETE"
site["query_ligand_identity_confirmed_v72"] = exact.astype(int)
site["residue_assertion_retained_v72"] = 1
site["site_default_web_inclusion_v72"] = 1
site["release_version_v72"] = VERSION
write(site, "interaction_site_v72.tsv.gz")

# 6. Unified complex module: no release-layer split, retain uncertainty fields.
complex_files = [
    "protein_complex_v71.tsv.gz", "complex_component_v71.tsv.gz", "complex_assembly_v71.tsv.gz",
    "complex_stoichiometry_v71.tsv.gz", "complex_state_v71.tsv.gz", "complex_binding_evidence_v71.tsv.gz"
]
for name in complex_files:
    d = pd.read_csv(SRC / name, sep="\t", low_memory=False)
    d["complex_module_status_v72"] = "UNIFIED_COMPLEX_MODULE"
    d["complex_module_inclusion_v72"] = 1
    d["record_retained_v72"] = 1
    d["release_version_v72"] = VERSION
    write(d, v72_name(name))

# 7. Evidence lineage source fields follow the corrected evidence registry.
lin = pd.read_csv(SRC / "evidence_lineage_v71.tsv.gz", sep="\t", low_memory=False)
srcmap = ev.set_index("evidence_id")[["source_database", "source_record_id", "source_provenance_resolution_v72"]]
lin["legacy_source_database_v72"] = lin.source_database.where(lin.evidence_id.isin(ev.loc[legacy, "evidence_id"]))
lin["legacy_source_record_id_v72"] = lin.source_record_id.where(lin.evidence_id.isin(ev.loc[legacy, "evidence_id"]))
lin["source_database"] = lin.evidence_id.map(srcmap.source_database).fillna(lin.source_database)
lin["source_record_id"] = lin.evidence_id.map(srcmap.source_record_id).fillna(lin.source_record_id)
lin["source_provenance_resolution_v72"] = lin.evidence_id.map(srcmap.source_provenance_resolution_v72)
lin["distinct_database_is_not_independent_experiment_v72"] = 1
lin["release_version_v72"] = VERSION
write(lin, "evidence_lineage_v72.tsv.gz")

# 8. Pair source lists and release roles rebuilt from corrected evidence.
pair = pd.read_csv(SRC / "protein_compound_pair_v71.tsv.gz", sep="\t", low_memory=False)
px = ev[["target_uniprot_id", "compound_internal_id", "source_database"]].copy()
px["source_database"] = px.source_database.fillna("UNSPECIFIED").map(source_split)
px = px.explode("source_database").drop_duplicates()
pg = px.groupby(["target_uniprot_id", "compound_internal_id"])["source_database"].agg(lambda z: ";".join(sorted(set(z)))).rename("source_databases_v72")
pair = pair.merge(pg, on=["target_uniprot_id", "compound_internal_id"], how="left")
pair["distinct_database_count_v72"] = pair.source_databases_v72.fillna("").map(lambda z: len(source_split(z)))
pair["legacy_source_databases_v72"] = pair.source_databases
pair["source_databases"] = pair.source_databases_v72
pair["distinct_database_count"] = pair["distinct_database_count_v72"]
pair["protein_release_role_v72"] = pair.target_uniprot_id.map(prot_role)
pair["default_web_inclusion_v72"] = pair.target_uniprot_id.map(prot_default).fillna(0).astype(int)
pair["searchable_inclusion_v72"] = pair.target_uniprot_id.isin(prot_role.index).astype(int)
pair["release_version_v72"] = VERSION
write(pair, "protein_compound_pair_v72.tsv.gz")

# 9. Preserve cross-classification unresolved roles explicitly.
cross = pd.read_csv(SRC / "protein_cross_classification_summary_v71.tsv.gz", sep="\t", low_memory=False)
cross["membrane_role_resolution_v72"] = cross.membrane_role_primary_v63.eq("family_defined_membrane_role_unresolved").map({True: "FAMILY_KNOWN_ROLE_UNRESOLVED", False: "ROLE_RESOLVED"})
cross["release_version_v72"] = VERSION
write(cross, "protein_cross_classification_summary_v72.tsv.gz")

# 10. HPA status semantics retained; no unsafe split of missing vs not measured.
expr_src = SRC / "expression_measurement_v71.tsv.gz"
expr_out = TABLES / "expression_measurement_v72.tsv.gz"
first = True
for chunk in pd.read_csv(expr_src, sep="\t", chunksize=200_000, low_memory=False):
    lim = chunk.missingness_limitation_v71.fillna("") if "missingness_limitation_v71" in chunk else pd.Series("", index=chunk.index)
    chunk["hpa_denominator_policy_v72"] = "MAPPED_AND_MEASURED_ONLY"
    chunk["hpa_missingness_status_v72"] = lim.eq("source_combines_missing_and_not_measured").map({True: "SOURCE_UNRESOLVED_MISSING_OR_NOT_MEASURED", False: "NO_COMBINED_MISSINGNESS_LIMITATION"})
    chunk["release_version_v72"] = VERSION
    chunk.to_csv(expr_out, sep="\t", index=False, compression="gzip", mode="wt" if first else "at", header=first)
    first = False

# Copy unchanged relation/entity tables into the V7.2 package with explicit routing metadata.
modified_v71 = {
    "canonical_protein_membrane_class_v71.tsv.gz", "protein_master_v71.tsv.gz", "positive_interaction_evidence_v71.tsv.gz",
    "evidence_target_entity_assignment_v71.tsv.gz", "isoform_membrane_class_v71.tsv.gz", "processed_form_membrane_class_v71.tsv.gz",
    "interaction_site_v71.tsv.gz", "evidence_lineage_v71.tsv.gz", "protein_compound_pair_v71.tsv.gz",
    "protein_cross_classification_summary_v71.tsv.gz", "expression_measurement_v71.tsv.gz", *complex_files
}
routing = []
for p in sorted(SRC.iterdir()):
    if not p.is_file() or p.name in modified_v71:
        continue
    target = TABLES / v72_name(p.name)
    if target.exists():
        target.chmod(0o666)
    shutil.copy2(p, target)
    target.chmod(0o666)
    routing.append((p.name, target.name, "INHERITED_UNCHANGED_FROM_V7.1.1"))
for p in sorted(TABLES.iterdir()):
    if p.is_file() and not any(r[1] == p.name for r in routing):
        routing.append((p.name.replace("_v72", "_v71"), p.name, "REBUILT_FOR_V7.2"))
pd.DataFrame(routing, columns=["v711_source_table", "v72_release_table", "routing_status"]).to_csv(META / "V72_TABLE_ROUTING.tsv", sep="\t", index=False)

# 11. Automated QA.
canonical = set(cls.canonical_uniprot_accession)
compound = pd.read_csv(SRC / "compound_master_v71.tsv.gz", sep="\t", usecols=["compound_internal_id"])
compound_set = set(compound.compound_internal_id)
qa = {
    "version": VERSION, "generated_at_utc": NOW,
    "baseline_files": len(baseline),
    "canonical_proteins": len(cls),
    "confirmed_core_E1_E2": int(cls.default_web_inclusion_v72.sum()),
    "candidate_E3": int((cls.release_role_v72 == "CANDIDATE_E3").sum()),
    "E3_in_default_core": int(((cls.evidence_level == "E3") & (cls.default_web_inclusion_v72 == 1)).sum()),
    "positive_evidence": len(ev),
    "mixed_source_remaining": int(ev.source_database.eq("v3.1_mixed_sources").sum()),
    "legacy_source_rows_recovered": int((ev.source_provenance_resolution_v72 == "RECOVERED_HISTORICAL_ASSAY_REGISTRY").sum()),
    "legacy_source_rows_unresolved": int((ev.source_provenance_resolution_v72 == "LEGACY_SOURCE_UNRESOLVED").sum()),
    "evidence_unknown_protein_fk": int((~ev.target_uniprot_id.isin(canonical)).sum()),
    "evidence_unknown_compound_fk": int((~ev.compound_internal_id.isin(compound_set)).sum()),
    "pair_rows": len(pair),
    "pair_mixed_source_remaining": int(pair.source_databases_v72.fillna("").str.contains("v3.1_mixed_sources", regex=False).sum()),
    "pair_formal_source_field_mixed_remaining": int(pair.source_databases.fillna("").str.contains("v3.1_mixed_sources", regex=False).sum()),
    "assignment_mixed_source_remaining": int(assign.source_database.eq("v3.1_mixed_sources").sum()),
    "isoform_rows_retained": len(iso),
    "isoform_unresolved_retained": int(iso.primary_membrane_class.eq("UNRESOLVED").sum()),
    "processed_form_rows_retained": len(pf),
    "processed_form_unresolved_retained": int(pf.primary_membrane_class.eq("UNRESOLVED").sum()),
    "site_rows": len(site),
    "site_tiers": site.site_tier_v72.value_counts().to_dict(),
    "complex_module_all_included": True,
    "complex_binding_rows": int(len(pd.read_csv(TABLES / "complex_binding_evidence_v72.tsv.gz", sep="\t", usecols=["complex_module_inclusion_v72"]))),
    "blocking_errors": [],
}
if qa["E3_in_default_core"] != 0: qa["blocking_errors"].append("E3_PRESENT_IN_DEFAULT_CORE")
if any(qa[k] != 0 for k in ["mixed_source_remaining", "pair_mixed_source_remaining", "pair_formal_source_field_mixed_remaining", "assignment_mixed_source_remaining"]): qa["blocking_errors"].append("MIXED_SOURCE_REMAINS")
if qa["legacy_source_rows_unresolved"] != 0: qa["blocking_errors"].append("LEGACY_SOURCE_UNRESOLVED")
if qa["evidence_unknown_protein_fk"] != 0 or qa["evidence_unknown_compound_fk"] != 0: qa["blocking_errors"].append("EVIDENCE_FOREIGN_KEY_FAILURE")
qa["status"] = "PASS" if not qa["blocking_errors"] else "FAIL"
(QA / "V72_VALIDATION_REPORT.json").write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding="utf-8")

# 12. Manifest and checksums.
manifest = []
for p in sorted(OUT.rglob("*")):
    if p.is_file() and p.name not in {"V72_MANIFEST.tsv", "SHA256SUMS.tsv"}:
        manifest.append((str(p.relative_to(OUT)), p.stat().st_size, hashlib.sha256(p.read_bytes()).hexdigest()))
mf = pd.DataFrame(manifest, columns=["file", "bytes", "sha256"])
mf.to_csv(META / "V72_MANIFEST.tsv", sep="\t", index=False)
mf[["sha256", "file"]].to_csv(META / "SHA256SUMS.tsv", sep="\t", index=False, header=False)
(OUT / "RELEASE_INFO.json").write_text(json.dumps({
    "release": VERSION, "status": "candidate", "generated_at_utc": NOW,
    "baseline": str(BASE), "validation": qa["status"],
    "policy": "E1/E2 confirmed core; E3 searchable candidate; unified complex module; canonical fallback for unspecified isoform; no record deletion"
}, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(qa, ensure_ascii=False, indent=2))
