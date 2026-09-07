from __future__ import annotations

import csv
import json
import shutil
from collections import Counter
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
V5 = ROOT / "releases" / "release_v5"
RELEASE = ROOT / "releases" / "release_v5_1"
MASTER_V5 = V5 / "human_membrane_protein_master_v5.tsv"
REVIEW = ROOT / "v51_review_working" / "analysis" / "r1109_final_review_v51.tsv"
OVERRIDES = ROOT / "configs" / "v51_review_overrides.tsv"
CHAIN = (
    ROOT
    / "v51_review_working"
    / "normalized"
    / "pdbtm_target_chain_topology_v51.tsv"
)
TODAY = date(2026, 7, 24).isoformat()


NEW_FIELDS = [
    "candidate_membrane_class_v51",
    "final_membrane_class_v51",
    "evidence_status_v51",
    "release_disposition_v51",
    "decision_code_v51",
    "decision_basis_v51",
    "cautions_v51",
    "review_date_v51",
    "review_method_v51",
    "v51_decision_layer",
    "manual_review_flag_v51",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def main() -> None:
    if RELEASE.exists():
        shutil.rmtree(RELEASE)
    for subdir in ("review", "normalized", "reports", "reproducibility/scripts", "reproducibility/configs"):
        (RELEASE / subdir).mkdir(parents=True, exist_ok=True)

    master = read_tsv(MASTER_V5)
    review_rows = read_tsv(REVIEW)
    review_by_accession = {row["target_uniprot_id"]: row for row in review_rows}
    if len(review_by_accession) != 1109:
        raise RuntimeError("Expected exactly 1,109 unique v5 R review rows")

    master_fields = list(master[0])
    output_fields = master_fields + [field for field in NEW_FIELDS if field not in master_fields]
    audit_rows: list[dict[str, str]] = []
    reviewed_r = 0
    for source in master:
        row = dict(source)
        accession = row["target_uniprot_id"]
        if row["release_tier_v5"] != "R":
            cls = row["release_tier_v5"]
            row.update(
                {
                    "candidate_membrane_class_v51": cls,
                    "final_membrane_class_v51": cls,
                    "evidence_status_v51": "confirmed" if cls in {"A", "B"} else "probable",
                    "release_disposition_v51": "included",
                    "decision_code_v51": "inherited_v5_class",
                    "decision_basis_v51": "Retained from the audited v5 A/B/C layer; v5.1 specifically resolves the former Tier R workflow queue.",
                    "cautions_v51": "",
                    "review_date_v51": TODAY,
                    "review_method_v51": "v5_inheritance_plus_v51_R_resolution",
                    "v51_decision_layer": "inherited_v5",
                    "manual_review_flag_v51": "0",
                }
            )
        else:
            review = review_by_accession.get(accession)
            if not review:
                raise RuntimeError(f"Missing v5.1 review for {accession}")
            reviewed_r += 1
            row.update(
                {
                    "candidate_membrane_class_v51": review[
                        "candidate_membrane_class_v51"
                    ],
                    "final_membrane_class_v51": review["final_membrane_class_v51"],
                    "evidence_status_v51": review["evidence_status_v51"],
                    "release_disposition_v51": review["inclusion_status_v51"],
                    "decision_code_v51": review["decision_code_v51"],
                    "decision_basis_v51": review["decision_basis_v51"],
                    "cautions_v51": review["cautions_v51"],
                    "review_date_v51": review["review_date_v51"],
                    "review_method_v51": review["review_method_v51"],
                    "v51_decision_layer": review["v51_decision_layer"],
                    "manual_review_flag_v51": review["manual_review_flag_v51"],
                }
            )
        audit_rows.append(row)
    if reviewed_r != 1109:
        raise RuntimeError(f"Resolved {reviewed_r}, expected 1,109")

    included = [row for row in audit_rows if row["release_disposition_v51"] == "included"]
    candidates = [row for row in audit_rows if row["release_disposition_v51"] == "candidate"]
    excluded = [row for row in audit_rows if row["release_disposition_v51"] == "excluded"]
    class_views = {
        cls: [row for row in included if row["final_membrane_class_v51"] == cls]
        for cls in ("A", "B", "C")
    }

    write_tsv(RELEASE / "human_membrane_protein_audit_master_v5_1.tsv", audit_rows, output_fields)
    write_tsv(RELEASE / "human_membrane_protein_master_v5_1.tsv", included, output_fields)
    write_tsv(RELEASE / "human_integral_membrane_view_v5_1.tsv", class_views["A"], output_fields)
    write_tsv(RELEASE / "human_lipid_anchored_membrane_view_v5_1.tsv", class_views["B"], output_fields)
    write_tsv(RELEASE / "human_peripheral_membrane_view_v5_1.tsv", class_views["C"], output_fields)
    write_tsv(RELEASE / "human_membrane_candidate_view_v5_1.tsv", candidates, output_fields)
    write_tsv(RELEASE / "human_membrane_excluded_audit_v5_1.tsv", excluded, output_fields)

    shutil.copy2(REVIEW, RELEASE / "review" / REVIEW.name)
    shutil.copy2(OVERRIDES, RELEASE / "review" / OVERRIDES.name)
    shutil.copy2(CHAIN, RELEASE / "normalized" / CHAIN.name)
    for script_name in (
        "download_v51_review_sources.py",
        "download_v51_aux_sources.py",
        "download_v51_pdbe_mappings.py",
        "build_v51_preliminary_review.py",
        "finalize_v51_review.py",
        "build_v51_release.py",
    ):
        shutil.copy2(
            ROOT / "scripts" / script_name,
            RELEASE / "reproducibility" / "scripts" / script_name,
        )
    shutil.copy2(
        OVERRIDES,
        RELEASE / "reproducibility" / "configs" / OVERRIDES.name,
    )

    stats = {
        "release_date": TODAY,
        "v5_audit_master_rows": len(audit_rows),
        "published_membrane_master_rows": len(included),
        "published_class_counts": {key: len(value) for key, value in class_views.items()},
        "former_tier_R_review_rows": len(review_rows),
        "former_tier_R_outcomes": {
            "included": sum(row["inclusion_status_v51"] == "included" for row in review_rows),
            "candidate": sum(row["inclusion_status_v51"] == "candidate" for row in review_rows),
            "excluded": sum(row["inclusion_status_v51"] == "excluded" for row in review_rows),
        },
        "former_tier_R_final_classes": dict(
            Counter(row["final_membrane_class_v51"] for row in review_rows)
        ),
        "former_tier_R_evidence_status": dict(
            Counter(row["evidence_status_v51"] for row in review_rows)
        ),
        "former_tier_R_review_reasons": dict(
            Counter(row["v5_review_reason"] for row in review_rows)
        ),
        "pdbtm_target_accessions_with_TM_chain": sum(
            row["pdbtm_target_chain_tm_flag"] == "1" for row in review_rows
        ),
        "manual_override_rows": sum(
            row["manual_review_flag_v51"] == "1" for row in review_rows
        ),
    }
    (RELEASE / "reports" / "V51_BUILD_STATS.json").write_text(
        json.dumps(stats, indent=2, ensure_ascii=False, sort_keys=True),
        encoding="utf-8",
    )

    source_rows = [
        {
            "source": "UniProtKB REST",
            "version_or_date": "retrieved 2026-07-24",
            "role": "Current reviewed sequence, feature, localization, GO cross-reference, domain, structure and literature fields for all 1,109 former R accessions",
            "url": "https://rest.uniprot.org/",
        },
        {
            "source": "Human Protein Atlas",
            "version_or_date": "v5 normalized source snapshot; methods checked 2026-07-24",
            "role": "Predicted membrane class and immunofluorescence subcellular location with Enhanced/Supported/Approved/Uncertain reliability",
            "url": "https://www.proteinatlas.org/",
        },
        {
            "source": "UniTmp high-throughput topology",
            "version_or_date": "v5 normalized source snapshot",
            "role": "Legacy high-throughput TM topology candidates",
            "url": "https://www.unitmp.org/",
        },
        {
            "source": "TMbed human proteome predictions",
            "version_or_date": "2022 precomputed predictions; exact sequence matched",
            "role": "Independent sequence-based TM/SP check; used only when the stored sequence exactly matched current UniProt",
            "url": "https://github.com/BernhoferM/TMbed",
        },
        {
            "source": "GOA Human and Gene Ontology",
            "version_or_date": "GOA 2026-07-08; GO basic 2026-06-26",
            "role": "Evidence-coded cellular-component localization",
            "url": "https://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/",
        },
        {
            "source": "InterPro",
            "version_or_date": "entry list 2026-06-10",
            "role": "Membrane-binding-capable domain context; never sufficient alone",
            "url": "https://www.ebi.ac.uk/interpro/",
        },
        {
            "source": "OPM",
            "version_or_date": "v5 normalized source snapshot",
            "role": "Structure-based peripheral membrane placement",
            "url": "https://opm.phar.umich.edu/",
        },
        {
            "source": "PDBTM plus PDBe/SIFTS",
            "version_or_date": "PDBTM snapshot plus mappings retrieved 2026-07-24",
            "role": "Chain-level correction: distinguish the reviewed target chain from a different membrane chain in the same PDB entry",
            "url": "https://www.ebi.ac.uk/pdbe/",
        },
    ]
    write_tsv(
        RELEASE / "SOURCE_REGISTRY_v5_1.tsv",
        source_rows,
        ["source", "version_or_date", "role", "url"],
    )

    (RELEASE / "README_v5_1.md").write_text(
        f"""# Human Membrane Protein Master v5.1

Release date: {TODAY}

v5.1 resolves every one of the 1,109 records that v5 placed in workflow Tier R. The biological classes and evidence state are now separate:

- **A**: integral membrane protein with one or more membrane-spanning/intramembrane segments.
- **B**: directly membrane-inserted without a transmembrane span, chiefly lipid-anchored proteins.
- **C**: peripheral or stable membrane-associated protein without a membrane-spanning segment.
- **Evidence status**: confirmed, probable, uncertain, or excluded.

The publication-ready master contains {len(included):,} included proteins: A {len(class_views['A']):,}, B {len(class_views['B']):,}, and C {len(class_views['C']):,}. The audit master retains all {len(audit_rows):,} v5 union records. Of the former R records, {stats['former_tier_R_outcomes']['included']:,} are included, {stats['former_tier_R_outcomes']['candidate']:,} remain explicit evidence-limited candidates, and {stats['former_tier_R_outcomes']['excluded']:,} are excluded from the publication-ready membrane set.

The candidate set is a completed review outcome, not an unfinished queue. It records proteins for which a candidate biological class can be proposed but current evidence does not justify a final A/B/C assignment.

## Principal files

- `human_membrane_protein_master_v5_1.tsv`: publication-ready included A/B/C set.
- `human_membrane_protein_audit_master_v5_1.tsv`: all v5 union records with v5.1 disposition.
- `human_integral_membrane_view_v5_1.tsv`, `human_lipid_anchored_membrane_view_v5_1.tsv`, `human_peripheral_membrane_view_v5_1.tsv`: class-specific views.
- `human_membrane_candidate_view_v5_1.tsv`: evidence-limited candidates.
- `human_membrane_excluded_audit_v5_1.tsv`: excluded records and reasons.
- `review/r1109_final_review_v51.tsv`: full evidence matrix for the former R queue.
""",
        encoding="utf-8",
    )

    (RELEASE / "INCLUSION_POLICY_v5_1.md").write_text(
        """# Inclusion policy v5.1

## Separation of biology and evidence

`candidate_membrane_class_v51` may be A, B, or C even when `final_membrane_class_v51` is `unknown`. Only records with `release_disposition_v51=included` and a final class A/B/C enter the publication-ready master.

## Evidence thresholds

- **Confirmed**: curated/experimental cellular membrane evidence, or a target chain itself demonstrated as transmembrane by structure-level mapping.
- **Probable**: concordant independent evidence, such as external topology plus a full current-sequence TM prediction, experimentally supported membrane-organelle localization for a non-integral protein, an experimental lipid anchor plus membrane localization, or OPM placement at a defined cellular membrane.
- **Uncertain**: a plausible candidate class with incomplete or ambiguous evidence, including HPA Approved localization alone, isolated hydropathy, or an N-terminal signal-anchor/signal-peptide ambiguity.
- **Excluded**: prediction not reproducible on the current sequence, HPA localization rated Uncertain with no corroboration, secreted/lumen proteins without stable cellular membrane association, or a PDB/PDBTM association belonging to another chain.

InterPro membrane-binding domains are contextual evidence and never sufficient alone. HPA `Approved` does not mean higher confidence than `Supported`; it denotes a different validation situation and is retained conservatively as uncertain when uncorroborated.
""",
        encoding="utf-8",
    )

    (RELEASE / "METHODS_v5_1.md").write_text(
        f"""# Methods v5.1

## Scope

The v5 Tier R queue contained 1,109 reviewed human accessions assembled from HPA location/prediction, UniTmp high-throughput topology, OPM, and PDBTM. v5.1 audited all 1,109 without changing the immutable v5 release.

## Current evidence refresh

Current UniProt records were retrieved on {TODAY}. TMbed 2022 human-proteome predictions were accepted only when their stored amino-acid sequence exactly matched the current UniProt sequence. GOA Human 2026-07-08 and GO basic 2026-06-26 supplied evidence-coded cellular-component annotations. InterPro 2026-06-10 supplied domain context.

## Chain-level structural correction

All PDB identifiers attached to former Tier R records were mapped to UniProt accessions through PDBe/SIFTS and compared with PDBTM chain topology. No reviewed target accession in Tier R mapped to a PDBTM chain that itself contained a TM segment. Therefore, PDBTM membership alone was not used to promote any R record; in many cases the membrane chain was a receptor or partner in the same complex.

## Decision process

An initial deterministic evidence matrix integrated current UniProt features, TMbed, hydropathy, GOA, HPA, UniTmp, OPM, InterPro, and target-chain mapping. Refined rules rejected short TM calls, non-reproduced external predictions, and signal-peptide confusion. Nineteen targeted overrides document high-impact edge cases such as retroviral products, lipid anchors, and false signal-peptide assignments. Every record retains the rule, evidence basis, cautions, source links, and review date.
""",
        encoding="utf-8",
    )

    (RELEASE / "DATA_DICTIONARY_v5_1.md").write_text(
        """# Data dictionary additions in v5.1

| Field | Meaning |
|---|---|
| `candidate_membrane_class_v51` | Best biological hypothesis (A/B/C), including evidence-limited candidates. |
| `final_membrane_class_v51` | Release-level biological class A/B/C, or `unknown` if evidence is insufficient/excluded. |
| `evidence_status_v51` | `confirmed`, `probable`, `uncertain`, or `excluded`. |
| `release_disposition_v51` | `included`, `candidate`, or `excluded`; controls publication-ready membership. |
| `decision_code_v51` | Machine-readable reason for the final disposition. |
| `decision_basis_v51` | Human-readable evidence summary. |
| `cautions_v51` | Conflicts and limitations that should accompany the record. |
| `v51_decision_layer` | Inherited v5, refined deterministic rule, or targeted manual override. |
| `manual_review_flag_v51` | `1` for targeted edge-case override, otherwise `0`. |

Legacy v5 fields are preserved for provenance. `release_tier_v5=R` is historical workflow state only and must not be interpreted as a biological class in v5.1.
""",
        encoding="utf-8",
    )

    (RELEASE / "REVIEW_REPORT_v5_1.md").write_text(
        f"""# Former Tier R review report

All 1,109 former Tier R records received a completed disposition.

| Outcome | Count |
|---|---:|
| Included A | {sum(r['final_membrane_class_v51'] == 'A' for r in review_rows)} |
| Included B | {sum(r['final_membrane_class_v51'] == 'B' for r in review_rows)} |
| Included C | {sum(r['final_membrane_class_v51'] == 'C' for r in review_rows)} |
| Evidence-limited candidate | {sum(r['inclusion_status_v51'] == 'candidate' for r in review_rows)} |
| Excluded | {sum(r['inclusion_status_v51'] == 'excluded' for r in review_rows)} |

The most consequential quality correction was chain-level PDBTM validation: zero former R target accessions were themselves a transmembrane PDBTM chain. The structural false-positive route has therefore been closed.
""",
        encoding="utf-8",
    )

    print(json.dumps(stats, indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
