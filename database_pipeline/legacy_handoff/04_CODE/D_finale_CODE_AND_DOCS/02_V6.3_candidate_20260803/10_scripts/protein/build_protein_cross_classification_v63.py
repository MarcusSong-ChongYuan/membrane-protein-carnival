#!/usr/bin/env python3
"""Build a non-destructive V6.3 five-axis protein annotation candidate.

The frozen V6.2 protein table is never modified.  This module joins a frozen
UniProt reviewed-human snapshot to GO, InterPro/Pfam and Reactome resources,
emits long-form source annotations and five orthogonal classification axes,
and writes a V6.3 candidate view with provenance-bearing added columns.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def args():
    p = argparse.ArgumentParser()
    p.add_argument("--protein-master", type=Path, required=True)
    p.add_argument("--uniprot", type=Path, required=True)
    p.add_argument("--uniprot-headers", type=Path, required=True)
    p.add_argument("--go-obo", type=Path, required=True)
    p.add_argument("--interpro-list", type=Path, required=True)
    p.add_argument("--pfam-clans", type=Path, required=True)
    p.add_argument("--reactome-pathways", type=Path, required=True)
    p.add_argument("--reactome-relations", type=Path, required=True)
    p.add_argument("--reactome-version", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--qa-dir", type=Path, required=True)
    return p.parse_args()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def clean(value) -> str:
    return "" if pd.isna(value) else str(value).strip()


def split_ids(value: str, prefix: str) -> list[str]:
    return sorted(set(re.findall(rf"{re.escape(prefix)}\d+", clean(value))))


def parse_go(path: Path):
    terms, parents, obsolete = {}, defaultdict(set), set()
    block = {}
    def commit():
        tid = block.get("id", "")
        if tid.startswith("GO:"):
            terms[tid] = (block.get("name", ""), block.get("namespace", ""))
            parents[tid].update(block.get("is_a", []))
            if block.get("is_obsolete") == "true":
                obsolete.add(tid)
    with path.open(encoding="utf-8") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line == "[Term]":
                commit(); block = {"is_a": []}
            elif not line:
                commit(); block = {}
            elif block is not None and ": " in line:
                key, value = line.split(": ", 1)
                if key == "is_a":
                    block.setdefault("is_a", []).append(value.split(" ! ", 1)[0])
                elif key in {"id", "name", "namespace", "is_obsolete"}:
                    block[key] = value
        commit()
    return terms, parents, obsolete


def ancestor_resolver(parents):
    cache = {}
    def resolve(term):
        if term in cache:
            return cache[term]
        seen, stack = set(), [term]
        while stack:
            current = stack.pop()
            for parent in parents.get(current, ()):
                if parent not in seen:
                    seen.add(parent); stack.append(parent)
        cache[term] = seen
        return seen
    return resolve


MF_ANCHORS = [
    ("ion_channel_activity", "GO:0005216"),
    ("transmembrane_transporter_activity", "GO:0022857"),
    ("signaling_receptor_activity", "GO:0038023"),
    ("catalytic_activity", "GO:0003824"),
    ("molecular_function_regulator", "GO:0098772"),
    ("structural_molecule_activity", "GO:0005198"),
    ("binding_activity", "GO:0005488"),
]

BP_ANCHORS = [
    ("cell_adhesion", "GO:0007155"),
    ("immune_system_process", "GO:0002376"),
    ("signal_transduction", "GO:0007165"),
    ("transmembrane_transport", "GO:0055085"),
    ("transport_and_localization", "GO:0006810"),
    ("developmental_process", "GO:0032502"),
    ("cell_death", "GO:0008219"),
    ("cell_cycle", "GO:0007049"),
    ("organelle_organization", "GO:0006996"),
    ("protein_metabolic_process", "GO:0019538"),
    ("lipid_metabolic_process", "GO:0006629"),
    ("carbohydrate_metabolic_process", "GO:0005975"),
    ("nucleic_acid_metabolic_process", "GO:0090304"),
    ("response_to_stimulus", "GO:0050896"),
    ("metabolic_process", "GO:0008152"),
]


def controlled_go_labels(ids, terms, ancestors, anchors, namespace):
    result = defaultdict(list)
    direct_in_namespace = []
    for go_id in ids:
        if go_id not in terms or terms[go_id][1] != namespace:
            continue
        direct_in_namespace.append(go_id)
        lineage = ancestors(go_id) | {go_id}
        for label, anchor in anchors:
            if anchor in lineage:
                result[label].append(go_id)
    if not result and direct_in_namespace:
        result[f"other_{namespace}_annotated"] = direct_in_namespace
    if not result:
        result["unclassified"] = []
    return result


def parse_semicolon(value: str) -> list[str]:
    return [item.strip() for item in clean(value).split(";") if item.strip()]


def normalize_family_text(value: str) -> list[str]:
    values = []
    for item in parse_semicolon(value):
        item = re.sub(r"^Belongs to the\s+", "", item, flags=re.I).strip(" .")
        item = re.sub(r"\s+family\.?$", " family", item, flags=re.I)
        if item and item.lower() not in {"protein family", "uncharacterized protein family"}:
            values.append(item)
    return list(dict.fromkeys(values))


def determine_membrane_roles(row, mf_labels, family_labels):
    roles, support = [], defaultdict(list)
    def add(label, reason):
        if label not in roles:
            roles.append(label)
        support[label].append(reason)
    legacy = clean(row.get("functional_primary_class_v5", "")).lower()
    sub = clean(row.get("functional_subclass_v5", "")).lower()
    gpcr = "|".join([clean(row.get("gpcrdb_class_v5", "")), clean(row.get("gpcrdb_family_v5", "")), clean(row.get("gtopdb_type_v5", ""))]).lower()
    tcdb = clean(row.get("tcdb_tcids_v5", ""))
    family_text = " | ".join(family_labels).lower()
    if gpcr and any(token in gpcr for token in ("gpcr", "class a", "class b", "class c", "adhesion receptor")):
        add("receptor", "GPCRdb/IUPHAR annotation")
    if tcdb:
        add("transporter", "TCDB identifier")
    legacy_map = {
        "receptor": "receptor", "ion_channel": "ion_channel", "transporter": "transporter",
        "enzyme": "membrane_associated_enzyme", "adhesion_recognition": "adhesion_recognition",
        "immune_or_cell_recognition": "immune_or_cell_recognition", "membrane_trafficking": "membrane_trafficking",
        "cytoskeleton_membrane_linker": "membrane_scaffold_or_linker", "signaling_regulator": "signaling_regulator",
    }
    if legacy in legacy_map:
        add(legacy_map[legacy], f"legacy V5 functional class={legacy}")
    for mf in mf_labels:
        if mf == "ion_channel_activity": add("ion_channel", "GO ion-channel activity")
        elif mf == "transmembrane_transporter_activity": add("transporter", "GO transmembrane-transporter activity")
        elif mf == "signaling_receptor_activity": add("receptor", "GO signaling-receptor activity")
        elif mf == "catalytic_activity": add("membrane_associated_enzyme", "GO catalytic activity")
        elif mf == "molecular_function_regulator": add("signaling_regulator", "GO molecular-function regulator")
        elif mf == "structural_molecule_activity": add("membrane_scaffold_or_linker", "GO structural-molecule activity")
    family_rules = [
        ("tetraspan", "membrane_organizer"), ("claudin", "junction_or_adhesion"),
        ("caveolin", "membrane_organizer"), ("snare", "membrane_trafficking"),
        ("syntaxin", "membrane_trafficking"), ("integrin", "adhesion_recognition"),
        ("cadherin", "adhesion_recognition"), ("connexin", "ion_channel"),
        ("aquaporin", "transporter"), ("solute carrier", "transporter"),
    ]
    for token, role in family_rules:
        if token in family_text:
            add(role, f"family rule:{token}")
    if not roles and family_labels:
        add("family_defined_membrane_role_unresolved", "family annotation present; role not over-inferred")
    if not roles:
        add("unclassified", "no role rule satisfied")
    priority = ["ion_channel", "transporter", "receptor", "membrane_associated_enzyme", "adhesion_recognition", "junction_or_adhesion", "immune_or_cell_recognition", "membrane_trafficking", "membrane_organizer", "membrane_scaffold_or_linker", "signaling_regulator", "family_defined_membrane_role_unresolved", "unclassified"]
    roles.sort(key=lambda x: priority.index(x) if x in priority else 999)
    return roles, support


def main():
    a = args(); a.output_dir.mkdir(parents=True, exist_ok=True); a.qa_dir.mkdir(parents=True, exist_ok=True)
    protein = pd.read_csv(a.protein_master, sep="\t", dtype=str, keep_default_na=False)
    uni = pd.read_csv(a.uniprot, sep="\t", compression="gzip", dtype=str, keep_default_na=False)
    uni = uni.rename(columns={
        "Entry": "target_uniprot_id", "Protein families": "uniprot_protein_families",
        "Gene Ontology (GO)": "go_all", "Gene Ontology (biological process)": "go_bp",
        "Gene Ontology (molecular function)": "go_mf", "Gene Ontology (cellular component)": "go_cc",
        "InterPro": "interpro_ids", "Pfam": "pfam_ids", "Reactome": "reactome_ids",
        "Function [CC]": "uniprot_function_comment", "Pathway": "uniprot_pathway_comment",
        "Sequence similarities": "uniprot_similarity_comment",
    })
    uni = uni.drop_duplicates("target_uniprot_id", keep="first").set_index("target_uniprot_id")

    go_terms, go_parents, go_obsolete = parse_go(a.go_obo)
    ancestors = ancestor_resolver(go_parents)
    interpro = pd.read_csv(a.interpro_list, sep="\t", dtype=str, keep_default_na=False)
    interpro_name = dict(zip(interpro["ENTRY_AC"], interpro["ENTRY_NAME"]))
    interpro_type = dict(zip(interpro["ENTRY_AC"], interpro["ENTRY_TYPE"]))
    pfam = pd.read_csv(a.pfam_clans, sep="\t", compression="gzip", dtype=str, header=None, names=["pfam_id", "clan_id", "clan_name", "pfam_name", "pfam_description"], keep_default_na=False)
    pfam_name = dict(zip(pfam.pfam_id, pfam.pfam_description))
    pathways = pd.read_csv(a.reactome_pathways, sep="\t", dtype=str, header=None, names=["reactome_id", "reactome_name", "species"], keep_default_na=False)
    human_path = pathways[pathways.species.eq("Homo sapiens")]
    reactome_name = dict(zip(human_path.reactome_id, human_path.reactome_name))
    relations = pd.read_csv(a.reactome_relations, sep="\t", dtype=str, header=None, names=["parent", "child"], keep_default_na=False)
    parents_r = defaultdict(set)
    for r in relations.itertuples(index=False):
        if r.parent.startswith("R-HSA-") and r.child.startswith("R-HSA-"):
            parents_r[r.child].add(r.parent)
    r_anc = ancestor_resolver(parents_r)
    human_ids = set(reactome_name)
    top_level = {rid for rid in human_ids if not (parents_r.get(rid, set()) & human_ids)}

    bridge_rows, summary_rows, source_rows, review_rows = [], [], [], []
    unmapped = []
    for row in protein.to_dict("records"):
        acc = row["target_uniprot_id"]
        if acc not in uni.index:
            unmapped.append({"target_uniprot_id": acc, "reason": "not_in_frozen_reviewed_human_uniprot_snapshot"})
            u = {}
        else:
            u = uni.loc[acc].to_dict()
        go_ids = split_ids(u.get("go_all", ""), "GO:")
        mf = controlled_go_labels(go_ids, go_terms, ancestors, MF_ANCHORS, "molecular_function")
        bp = controlled_go_labels(go_ids, go_terms, ancestors, BP_ANCHORS, "biological_process")
        ipr_ids = split_ids(u.get("interpro_ids", ""), "IPR")
        pf_ids = split_ids(u.get("pfam_ids", ""), "PF")
        reactome_ids = sorted(set(re.findall(r"R-HSA-\d+", clean(u.get("reactome_ids", "")))))

        family_labels = normalize_family_text(u.get("uniprot_protein_families", ""))
        family_support = {label: ["UniProt protein-family annotation"] for label in family_labels}
        if not family_labels:
            for ipr in ipr_ids:
                if interpro_type.get(ipr) in {"Family", "Homologous_superfamily"} and interpro_name.get(ipr):
                    label = interpro_name[ipr]
                    if label not in family_labels:
                        family_labels.append(label); family_support[label] = [ipr]
        if not family_labels:
            for pf in pf_ids[:5]:
                if pfam_name.get(pf):
                    family_labels.append(pfam_name[pf]); family_support[pfam_name[pf]] = [pf]
        if not family_labels:
            family_labels = ["unclassified"]; family_support["unclassified"] = []

        specialist = []
        for label, value in [
            ("GPCRdb", clean(row.get("gpcrdb_family_v5", "")) or clean(row.get("gpcrdb_class_v5", ""))),
            ("IUPHAR_GtoPdb", clean(row.get("gtopdb_family_name_v5", "")) or clean(row.get("gtopdb_type_v5", ""))),
            ("TCDB", clean(row.get("tcdb_tcids_v5", ""))),
        ]:
            if value:
                specialist.append((f"{label}: {value}", label))
        ec_ids = sorted(set(re.findall(r"\b\d+\.\d+\.\d+\.[\d-]+\b", clean(u.get("Protein names", "")) + " " + clean(row.get("protein_name", "")))))
        for ec in ec_ids:
            specialist.append((f"EC: {ec}", "UniProt EC"))
        if not specialist:
            specialist = [("no_specialist_classification", "none")]

        roles, role_support = determine_membrane_roles(row, list(mf), family_labels)
        reactome_tops = set()
        for rid in reactome_ids:
            lineage = r_anc(rid) | {rid}
            roots = lineage & top_level
            reactome_tops.update(roots or ({rid} if rid in reactome_name else set()))
        reactome_top_names = sorted({reactome_name[rid] for rid in reactome_tops if rid in reactome_name})

        axes = {
            "structural_family": [(label, ";".join(family_support.get(label, []))) for label in family_labels],
            "molecular_function": [(label, ";".join(ids)) for label, ids in mf.items()],
            "biological_process": [(label, ";".join(ids)) for label, ids in bp.items()],
            "membrane_role": [(label, ";".join(role_support.get(label, []))) for label in roles],
            "specialist_classification": [(label, source) for label, source in specialist],
        }
        for axis, labels in axes.items():
            for index, (label, support_value) in enumerate(labels):
                bridge_rows.append({
                    "target_uniprot_id": acc, "classification_axis": axis, "classification_label": label,
                    "primary_flag": int(index == 0), "supporting_source_or_ids": support_value,
                    "classification_rule_version": "protein_cross_classification_v0.1.0",
                    "source_release": "UniProt 2026_02; GO 2026-06-15; Reactome 97; InterPro/Pfam frozen 2026-08-03",
                })
        unresolved_axes = sum(1 for axis, labels in axes.items() if labels[0][0] in {"unclassified", "no_specialist_classification"})
        display_status = "five_axis_annotated" if unresolved_axes == 0 else ("partially_annotated" if unresolved_axes < 4 else "functionally_unresolved")
        db_sources = sorted(set(parse_semicolon(row.get("independent_membrane_sources_v5", ""))))
        summary = {
            "target_uniprot_id": acc,
            "structural_family_primary_v63": family_labels[0], "structural_family_all_v63": "; ".join(family_labels),
            "molecular_function_primary_v63": next(iter(mf)), "molecular_function_all_v63": "; ".join(mf),
            "biological_process_primary_v63": next(iter(bp)), "biological_process_all_v63": "; ".join(bp),
            "membrane_role_primary_v63": roles[0], "membrane_role_all_v63": "; ".join(roles),
            "specialist_classification_primary_v63": specialist[0][0], "specialist_classification_all_v63": "; ".join(x[0] for x in specialist),
            "reactome_top_level_v63": "; ".join(reactome_top_names),
            "display_classification_status_v63": display_status,
            "distinct_contributing_membrane_database_count_v63": len(db_sources),
            "contributing_membrane_databases_v63": "; ".join(db_sources),
            "membrane_database_count_semantics_v63": "distinct contributing database labels; not independent experiments",
            "protein_annotation_module_version_v63": "v0.1.0",
        }
        summary_rows.append(summary)
        source_rows.append({
            "target_uniprot_id": acc, "uniprot_mapping_status": "mapped" if u else "unmapped",
            "go_ids": ";".join(go_ids), "interpro_ids": ";".join(ipr_ids), "pfam_ids": ";".join(pf_ids),
            "reactome_ids": ";".join(reactome_ids), "reactome_top_level_ids": ";".join(sorted(reactome_tops)),
            "uniprot_protein_families": clean(u.get("uniprot_protein_families", "")),
            "uniprot_function_comment": clean(u.get("uniprot_function_comment", "")),
            "uniprot_pathway_comment": clean(u.get("uniprot_pathway_comment", "")),
            "annotation_snapshot": "UniProt 2026_02",
        })
        legacy_status = clean(row.get("classification_status_v5", ""))
        e = clean(row.get("evidence_level_v52", ""))
        if legacy_status in {"family_level_classified_function_unresolved", "unclassified"} and e in {"E1", "E2"}:
            changed = roles[0] not in {"family_defined_membrane_role_unresolved", "unclassified"}
            review_rows.append({
                "target_uniprot_id": acc, "approved_symbol": clean(row.get("approved_symbol", "")),
                "membrane_evidence_level": e, "legacy_classification_status": legacy_status,
                "legacy_primary_class": clean(row.get("functional_primary_class_v5", "")),
                "proposed_membrane_role_primary_v63": roles[0],
                "proposed_structural_family_primary_v63": family_labels[0],
                "proposed_molecular_function_primary_v63": next(iter(mf)),
                "proposed_biological_process_primary_v63": next(iter(bp)),
                "automatic_reclassification_status": "specific_role_proposed" if changed else "multiaxis_context_added_role_unresolved",
                "manual_review_priority": "high" if clean(row.get("canonical_best_binding_evidence_level_v61", "")) in {"BE1", "BE2"} or clean(row.get("best_disease_evidence_level", "")) in {"DE1", "DE2"} else "standard",
            })

    source_df = pd.DataFrame(source_rows)
    bridge_df = pd.DataFrame(bridge_rows)
    summary_df = pd.DataFrame(summary_rows)
    review_df = pd.DataFrame(review_rows)
    source_df.to_csv(a.output_dir / "protein_function_source_annotation_v0_1.tsv.gz", sep="\t", index=False, compression="gzip")
    bridge_df.to_csv(a.output_dir / "protein_cross_classification_v0_1.tsv.gz", sep="\t", index=False, compression="gzip")
    summary_df.to_csv(a.output_dir / "protein_cross_classification_summary_v0_1.tsv", sep="\t", index=False)
    review_df.to_csv(a.output_dir / "protein_reclassification_priority_v0_1.tsv", sep="\t", index=False)
    pd.DataFrame(unmapped, columns=["target_uniprot_id", "reason"]).to_csv(a.output_dir / "protein_annotation_unmapped_review_v0_1.tsv", sep="\t", index=False)
    candidate = protein.merge(summary_df, on="target_uniprot_id", how="left", validate="one_to_one")
    candidate.to_csv(a.output_dir / "human_membrane_protein_master_v6_3_candidate.tsv.gz", sep="\t", index=False, compression="gzip")

    headers = a.uniprot_headers.read_text(encoding="utf-8", errors="replace")
    uniprot_release = re.search(r"X-UniProt-Release:\s*([^\r\n]+)", headers, re.I)
    registry = [
        {"source": "UniProtKB", "version": uniprot_release.group(1).strip() if uniprot_release else "2026_02", "local_file": str(a.uniprot), "sha256": sha256(a.uniprot), "role": "canonical reviewed-human annotation and external cross-references"},
        {"source": "Gene Ontology", "version": "2026-06-15", "local_file": str(a.go_obo), "sha256": sha256(a.go_obo), "role": "ontology-derived molecular-function and biological-process axes"},
        {"source": "InterPro", "version": "current_release snapshot 2026-08-03", "local_file": str(a.interpro_list), "sha256": sha256(a.interpro_list), "role": "family/superfamily names"},
        {"source": "Pfam", "version": "current_release snapshot 2026-08-03", "local_file": str(a.pfam_clans), "sha256": sha256(a.pfam_clans), "role": "fallback structural-family names"},
        {"source": "Reactome", "version": clean(a.reactome_version.read_text()), "local_file": str(a.reactome_pathways), "sha256": sha256(a.reactome_pathways), "role": "pathway identifiers and hierarchy"},
    ]
    pd.DataFrame(registry).to_csv(a.output_dir / "protein_annotation_source_registry_v0_1.tsv", sep="\t", index=False)

    checks = {
        "protein_count_preserved": len(candidate) == len(protein) == 10997,
        "protein_primary_key_unique": candidate.target_uniprot_id.is_unique,
        "all_proteins_have_five_axes": bridge_df.groupby("target_uniprot_id").classification_axis.nunique().eq(5).all(),
        "all_bridge_foreign_keys_valid": set(bridge_df.target_uniprot_id) <= set(protein.target_uniprot_id),
        "candidate_columns_populated": summary_df.notna().all().all(),
        "legacy_columns_preserved": all(col in candidate.columns for col in protein.columns),
        "review_scope_only_e1_e2": set(review_df.membrane_evidence_level) <= {"E1", "E2"},
    }
    counts = {
        "proteins": len(protein), "uniprot_mapped": len(protein) - len(unmapped), "uniprot_unmapped": len(unmapped),
        "classification_bridge_rows": len(bridge_df), "priority_reclassification_rows": len(review_df),
        "specific_role_proposed": int((review_df.automatic_reclassification_status == "specific_role_proposed").sum()),
        "role_unresolved_after_multiaxis": int((summary_df.membrane_role_primary_v63 == "unclassified").sum()),
        "family_defined_role_unresolved": int((summary_df.membrane_role_primary_v63 == "family_defined_membrane_role_unresolved").sum()),
        "display_status": summary_df.display_classification_status_v63.value_counts().to_dict(),
        "membrane_role_primary": summary_df.membrane_role_primary_v63.value_counts().to_dict(),
    }
    validation = {"module": "protein_cross_classification", "version": "0.1.0", "generated_at_utc": datetime.now(timezone.utc).isoformat(), "status": "PASS" if all(checks.values()) else "FAIL", "checks": checks, "counts": counts, "input_sha256": {"protein_master_v6_2": sha256(a.protein_master), "uniprot": sha256(a.uniprot), "go": sha256(a.go_obo), "interpro": sha256(a.interpro_list), "pfam": sha256(a.pfam_clans), "reactome_pathways": sha256(a.reactome_pathways), "reactome_relations": sha256(a.reactome_relations)}}
    (a.qa_dir / "PROTEIN_CROSS_CLASSIFICATION_V0_1_VALIDATION.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(validation, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
