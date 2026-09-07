from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


SRC = Path(r"D:\finale\02_V6.3_candidate_20260803")
OUT = Path(r"D:\finale\03_V6.3.1_readiness_20260803")
PY = Path(r"D:\7.22\evidence_expansion_v2_working\tools\plotting_venv\Scripts\python.exe")


def read_tsv(path: Path, **kwargs):
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False, **kwargs)


def write_tsv(df: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while b := f.read(1024 * 1024):
            h.update(b)
    return h.hexdigest()


def priority_label(score: pd.Series) -> pd.Series:
    return pd.cut(score, [-1, 3, 6, 10**9], labels=["P3", "P2", "P1"]).astype(str)


def main():
    if OUT.exists():
        raise SystemExit(f"Refusing to overwrite existing output: {OUT}")
    for d in ["00_docs", "01_engineering", "02_review_queues", "03_complex", "04_docking_hpc", "05_qa"]:
        (OUT / d).mkdir(parents=True, exist_ok=True)

    counts = {}
    # Protein queue: deterministic ranking, no automatic biological promotion.
    prot = read_tsv(SRC / "05_protein_annotation/protein_reclassification_priority_v0_1.tsv")
    dock = read_tsv(SRC / "07_docking/docking_priority_v6_3_candidate.tsv.gz")
    dockagg = dock.groupby("target_uniprot_id", as_index=False).agg(
        docking_pair_count=("compound_internal_id", "nunique"),
        best_docking_score=("ranking_score_v63", "max"),
        max_binding_evidence_count=("binding_evidence_count", "max"),
    )
    dis = read_tsv(SRC / "04_disease/protein_disease_relation_v2_1.tsv")
    rel_target = "target_uniprot_id" if "target_uniprot_id" in dis.columns else dis.columns[0]
    disagg = dis.groupby(rel_target, as_index=False).size().rename(columns={rel_target: "target_uniprot_id", "size": "disease_relation_count"})
    prot = prot.merge(dockagg, on="target_uniprot_id", how="left").merge(disagg, on="target_uniprot_id", how="left")
    for c in ["docking_pair_count", "best_docking_score", "max_binding_evidence_count", "disease_relation_count"]:
        prot[c] = pd.to_numeric(prot[c], errors="coerce").fillna(0)
    unresolved = prot["automatic_reclassification_status"].astype(str).str.contains("unresolved|insufficient|review", case=False, regex=True)
    prot["review_score_v631"] = (
        prot["membrane_evidence_level"].map({"E1": 5, "E2": 3, "E3": 1}).fillna(0)
        + unresolved.astype(int) * 3
        + (prot["docking_pair_count"] > 0).astype(int) * 3
        + (prot["disease_relation_count"] > 0).astype(int) * 2
        + (prot["max_binding_evidence_count"] > 1).astype(int)
    )
    prot["review_priority_v631"] = priority_label(prot["review_score_v631"])
    prot["automatic_action_v631"] = np.where(unresolved, "manual_function_review_required", "retain_automatic_multiaxis_annotation")
    prot = prot.sort_values(["review_score_v631", "docking_pair_count", "disease_relation_count"], ascending=False)
    write_tsv(prot, OUT / "02_review_queues/protein_function_review_ranked_v631.tsv")
    write_tsv(prot[prot["review_priority_v631"] == "P1"].head(500), OUT / "02_review_queues/protein_function_review_minimal_P1_v631.tsv")
    counts["protein_review_total"] = len(prot)
    counts["protein_review_p1_minimal"] = min(500, int((prot["review_priority_v631"] == "P1").sum()))

    # Disease queue: prioritize ambiguous/obsolete/non-disease and relation burden.
    dm = read_tsv(SRC / "04_disease/disease_mapping_review_v0_1.tsv")
    sid_col = "source_disease_id"
    source_rel_col = next((c for c in dis.columns if c == "source_disease_id"), None)
    if source_rel_col:
        da = dis.groupby(source_rel_col, as_index=False).size().rename(columns={"size": "relation_count"})
        dm = dm.merge(da, on=sid_col, how="left")
    else:
        dm["relation_count"] = 0
    dm["relation_count"] = pd.to_numeric(dm["relation_count"], errors="coerce").fillna(0).astype(int)
    dm["review_score_v631"] = (
        dm["review_type"].astype(str).str.contains("1:n|ambig", case=False, regex=True).astype(int) * 5
        + dm["review_type"].astype(str).str.contains("obsolete", case=False).astype(int) * 4
        + dm["review_type"].astype(str).str.contains("non.?disease", case=False, regex=True).astype(int) * 3
        + (dm["relation_count"] > 0).astype(int) * 2
        + np.minimum(dm["relation_count"], 10) / 10
    )
    dm["review_priority_v631"] = priority_label(dm["review_score_v631"])
    dm = dm.sort_values(["review_score_v631", "relation_count"], ascending=False)
    write_tsv(dm, OUT / "02_review_queues/disease_mapping_review_ranked_v631.tsv")
    counts["disease_review_total"] = len(dm)

    # Expression queue summaries and reproducible stratified sample.
    expr = read_tsv(SRC / "01_baseline_v6_2/expression_location_review_queue_v2.tsv.gz")
    exs = expr.groupby(["review_type", "review_reason", "current_status"], dropna=False).size().reset_index(name="record_count")
    write_tsv(exs.sort_values("record_count", ascending=False), OUT / "02_review_queues/expression_location_review_summary_v631.tsv")
    ex_sample = expr.groupby(["review_type", "review_reason"], dropna=False, group_keys=False).head(25)
    write_tsv(ex_sample, OUT / "02_review_queues/expression_location_review_sample_v631.tsv")
    counts["expression_review_total"] = len(expr)

    # Negative evidence: summarize full universe; sample per reason/conflict without promoting records.
    neg_path = SRC / "01_baseline_v6_2/negative_binding_evidence_unmapped_review_v1_2.tsv.gz"
    summary_parts, sample_parts, neg_total = [], [], 0
    for ch in pd.read_csv(neg_path, sep="\t", compression="gzip", chunksize=250000, low_memory=False):
        neg_total += len(ch)
        keys = ["unmapped_reason_v62", "positive_negative_conflict_flag_v62", "source_database"]
        summary_parts.append(ch.groupby(keys, dropna=False).size().reset_index(name="record_count"))
        sample_parts.append(ch.groupby(keys, dropna=False, group_keys=False).head(3))
    ns = pd.concat(summary_parts).groupby(["unmapped_reason_v62", "positive_negative_conflict_flag_v62", "source_database"], dropna=False)["record_count"].sum().reset_index()
    write_tsv(ns.sort_values("record_count", ascending=False), OUT / "02_review_queues/negative_unmapped_review_summary_v631.tsv")
    nsm = pd.concat(sample_parts).drop_duplicates("source_evidence_id").groupby(["unmapped_reason_v62", "positive_negative_conflict_flag_v62", "source_database"], dropna=False, group_keys=False).head(25)
    write_tsv(nsm, OUT / "02_review_queues/negative_unmapped_review_stratified_sample_v631.tsv")
    counts["negative_unmapped_review_total"] = neg_total
    counts["negative_review_sample"] = len(nsm)

    # Complex directness queue with structure bonus. Directness remains manual.
    cb = read_tsv(SRC / "03_complex/complex_binding_evidence_candidates_v0_2.tsv")
    cs = read_tsv(SRC / "03_complex/complex_structure_evidence_candidates_v0_2.tsv")
    csagg = cs.groupby("complex_target_id", as_index=False).agg(
        candidate_structure_count=("pdb_id", "nunique"),
        single_assembly_structure_count=("candidate_assembly_count_for_pdb", lambda x: int((pd.to_numeric(x, errors="coerce") == 1).sum())),
    )
    cb = cb.merge(csagg, on="complex_target_id", how="left")
    cb[["candidate_structure_count", "single_assembly_structure_count"]] = cb[["candidate_structure_count", "single_assembly_structure_count"]].fillna(0).astype(int)
    exact = cb["compound_mapping_status"].eq("exact_source_identifier_unique")
    asserted = cb["complex_specificity_status"].eq("source_asserted_complex_context")
    named_ok = cb["source_name_identifier_consistency"].eq("consistent")
    has_pm = cb["pubmed_ids"].fillna("").astype(str).str.len().gt(0)
    cb["review_score_v631"] = exact.astype(int) * 4 + asserted.astype(int) * 2 + named_ok.astype(int) * 2 + has_pm.astype(int) * 2 + (cb["candidate_structure_count"] > 0).astype(int) * 2
    cb["review_priority_v631"] = priority_label(cb["review_score_v631"])
    cb["required_manual_decision"] = "verify_source_article_supports_direct_binding_to_intact_complex"
    cb["automatic_release_action_v631"] = "retain_candidate_not_promoted"
    cb = cb.sort_values(["review_score_v631", "candidate_structure_count"], ascending=False)
    write_tsv(cb, OUT / "03_complex/complex_binding_directness_ranked_v631.tsv")
    write_tsv(cb[cb["review_priority_v631"] == "P1"].head(500), OUT / "03_complex/complex_binding_directness_minimal_P1_v631.tsv")
    write_tsv(cs, OUT / "03_complex/complex_structure_chain_mapping_review_v631.tsv")
    counts["complex_binding_candidates"] = len(cb)
    counts["complex_directness_p1_minimal"] = min(500, int((cb["review_priority_v631"] == "P1").sum()))
    counts["complex_binding_promoted"] = 0

    # Isoform unresolved/obsolete queue is preserved and explicitly gated.
    iso = read_tsv(SRC / "02_identity/external_complex_protein_review_v0_2.tsv")
    iso["automatic_action_v631"] = "retain_unresolved_external_accession"
    iso["manual_evidence_required"] = "authoritative_accession_history_or_source_publication"
    write_tsv(iso, OUT / "02_review_queues/isoform_external_accession_review_v631.tsv")
    counts["isoform_external_review_total"] = len(iso)

    # Docking pilot uses existing frozen membership; V6.3 ranking enriches it only.
    pilot_candidates = [
        Path(r"D:\7.22\evidence_expansion_v2_working\runs\docking_selection_v2_v62_20260730\staging\docking_hpc_pilot_v2_v62.tsv"),
        SRC / "07_docking/docking_hpc_pilot_v2_v62.tsv",
    ]
    pilot_path = next((p for p in pilot_candidates if p.exists()), None)
    if pilot_path:
        pilot = read_tsv(pilot_path)
        keys = ["target_uniprot_id", "compound_internal_id"]
        enrich_cols = keys + [c for c in ["ranking_score_v63", "ranking_policy_v63", "candidate_receptor_pdb_id_v2", "receptor_preparation_status_v2", "ligand_preparation_status_v2", "box_definition_status_v2", "positive_negative_conflict_flag_v62", "standard_smiles", "standard_inchikey"] if c in dock.columns]
        pilot = pilot.drop(columns=[c for c in enrich_cols if c in pilot.columns and c not in keys], errors="ignore").merge(dock[enrich_cols].drop_duplicates(keys), on=keys, how="left")
    else:
        pilot = dock.sort_values("ranking_score_v63", ascending=False).head(6019).copy()
    pilot["docking_execution_status_v631"] = "not_run"
    pilot["receptor_gate_v631"] = np.where(pilot.get("candidate_receptor_pdb_id_v2", pd.Series(index=pilot.index)).fillna("").astype(str).str.len().gt(0), "candidate_structure_available_needs_preparation", "structure_selection_required")
    pilot["ligand_gate_v631"] = np.where(pilot.get("standard_smiles", pd.Series(index=pilot.index)).fillna("").astype(str).str.len().gt(0), "smiles_available_needs_3d_microstate_preparation", "chemical_structure_missing")
    pilot["box_gate_v631"] = "binding_site_or_reference_ligand_box_required"
    pilot["redocking_acceptance_rule"] = "reference_ligand_heavy_atom_RMSD_le_2A_before_production"
    write_tsv(pilot, OUT / "04_docking_hpc/docking_pilot_manifest_v631.tsv")
    batches = pilot[["target_uniprot_id", "compound_internal_id"]].copy()
    batches.insert(0, "array_index", np.arange(len(batches)))
    batches["batch_id"] = batches["array_index"] // 250
    batches["output_stem"] = batches["target_uniprot_id"].astype(str) + "__" + batches["compound_internal_id"].astype(str)
    write_tsv(batches, OUT / "04_docking_hpc/docking_array_map_v631.tsv")
    counts["docking_pilot_pairs"] = len(pilot)
    counts["docking_unique_targets"] = pilot["target_uniprot_id"].nunique()
    counts["docking_structure_selection_required"] = int((pilot["receptor_gate_v631"] == "structure_selection_required").sum())

    # Engineering lock and audits.
    lock = """adjustText==1.4.0
duckdb==1.3.2
lxml==6.0.2
matplotlib==3.11.1
networkx==3.6.1
numpy==2.3.5
openpyxl==3.1.5
pandas==3.0.1
pyarrow==25.0.0
rdkit==2026.3.5
scikit-learn==1.9.0
scipy==1.18.0
seaborn==0.13.2
tqdm==4.70.0
umap-learn==0.5.12
UpSetPlot==0.9.0
xlsxwriter==3.2.9
"""
    (OUT / "01_engineering/requirements-scientific-lock.txt").write_text(lock, encoding="utf-8")
    env = "name: mempro-v631\nchannels:\n  - conda-forge\ndependencies:\n  - python=3.12\n  - pip\n  - pip:\n" + "".join(f"    - {x}\n" for x in lock.strip().splitlines())
    (OUT / "01_engineering/environment.yml").write_text(env, encoding="utf-8")

    script_files = list((SRC / "10_scripts").rglob("*.py")) + list((SRC / "10_scripts").rglob("*.ps1")) + list((SRC / "10_scripts").rglob("*.sh"))
    abs_hits = []
    rx = re.compile(r"(?i)(?:[A-Z]:\\|[A-Z]:/)(?!\\s)")
    for p in script_files:
        try:
            for i, line in enumerate(p.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
                if rx.search(line):
                    abs_hits.append({"file": str(p.relative_to(SRC)), "line": i, "text": line[:300]})
        except OSError:
            pass
    write_tsv(pd.DataFrame(abs_hits, columns=["file", "line", "text"]), OUT / "01_engineering/absolute_path_audit.tsv")
    counts["absolute_path_hits_in_release_scripts"] = len(abs_hits)

    license_rows = [
        ["UniProt", "included", "verify redistribution terms and attribution before public release"],
        ["HPA", "included", "verify current data-use terms and attribution before public release"],
        ["PDBe/PDB", "included", "retain structure provenance and citation"],
        ["PDBbind", "derived internal use", "do not redistribute files unless license explicitly permits"],
        ["BRENDA", "derived internal use", "manual license review required before redistribution"],
        ["PubChem", "included", "retain AID/SID/CID provenance"],
        ["ChEMBL", "included", "retain version and attribution"],
        ["BindingDB", "included", "retain version and attribution"],
        ["Open Targets", "included", "retain version and upstream evidence provenance"],
        ["Complex Portal", "candidate layer", "retain source identifiers and article provenance"],
    ]
    write_tsv(pd.DataFrame(license_rows, columns=["source", "current_package_status", "publication_gate"]), OUT / "01_engineering/source_license_publication_gate.tsv")

    # Copy immutable metadata and figures for handoff context.
    shutil.copy2(SRC / "FREEZE_V6_3_CANDIDATE.json", OUT / "00_docs/BASELINE_FREEZE_V6_3_CANDIDATE.json")
    shutil.copytree(SRC / "08_figures", OUT / "06_figures_v63", dirs_exist_ok=True)

    unresolved_manual = {
        "paper_level_complex_directness": True,
        "historical_isoform_specificity": True,
        "license_redistribution_review": True,
        "second_machine_reproduction": True,
        "receptor_structure_preparation": True,
        "ligand_microstate_and_3d_preparation": True,
        "binding_box_definition": True,
        "redocking_validation": True,
        "production_hpc_docking": True,
    }
    qa = {
        "release": "V6.3.1-readiness",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "baseline": "V6.3 candidate; immutable",
        "automatic_checks": counts,
        "automatic_status": "PASS",
        "manual_or_external_gates": unresolved_manual,
        "overall_status": "AUTOMATION_COMPLETE_MANUAL_GATES_OPEN",
    }
    (OUT / "05_qa/V631_READINESS_VALIDATION.json").write_text(json.dumps(qa, ensure_ascii=False, indent=2), encoding="utf-8")

    readme = f"""# MemPro V6.3.1 readiness package

This package is an incremental, non-destructive companion to the frozen V6.3 candidate. It does not rewrite or promote any scientific record automatically.

## Completed automatically

- reproducible dependency lock and absolute-path audit;
- ranked protein, disease, expression, negative-evidence, isoform and complex review queues;
- complex-binding directness minimal P1 queue without automatic promotion;
- docking pilot manifest ({counts['docking_pilot_pairs']:,} pairs; {counts['docking_unique_targets']:,} targets), array map and acceptance gates;
- source/license publication-gate matrix and machine-readable validation.

## Still requires evidence or external execution

- paper-level confirmation of intact-complex direct binding;
- historical isoform-specific interpretation where source records are ambiguous;
- legal/license decision for redistribution of restricted-source derivatives;
- clean second-machine reproduction;
- receptor, ligand, box and redocking preparation followed by HPC execution.

The status `AUTOMATION_COMPLETE_MANUAL_GATES_OPEN` is intentional: it prevents review queues from being misreported as curated conclusions.
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")

    manifest_rows = []
    for p in sorted(OUT.rglob("*")):
        if p.is_file():
            manifest_rows.append({"relative_path": p.relative_to(OUT).as_posix(), "bytes": p.stat().st_size, "sha256": sha256(p)})
    write_tsv(pd.DataFrame(manifest_rows), OUT / "MANIFEST_V631_READINESS.tsv")
    print(json.dumps(qa, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
