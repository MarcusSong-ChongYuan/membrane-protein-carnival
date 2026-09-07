from __future__ import annotations

import csv
import json
import re
import shutil
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_RELEASE = ROOT / "releases" / "release_v5_1"
SOURCE_MASTER = SOURCE_RELEASE / "human_membrane_protein_audit_master_v5_1.tsv"
RELEASE = ROOT / "releases" / "release_v5_2_evidence_tiered"
RELEASE_DATE = "2026-07-24"

NEW_FIELDS = [
    "membrane_class_v52",
    "evidence_level_v52",
    "evidence_label_zh_v52",
    "evidence_label_en_v52",
    "direct_experimental_flag_v52",
    "release_scope_v52",
    "website_default_v52",
    "evidence_rule_v52",
    "evidence_basis_v52",
]

LABELS = {
    "E1": ("实验证实", "experimentally confirmed", "core", "1", "1"),
    "E2": ("强证据支持", "strongly supported", "supported", "0", "1"),
    "E3": ("预测候选", "prediction candidate", "candidate", "0", "0"),
    "E0": ("排除", "excluded", "excluded", "0", "0"),
}

MEMBRANE_CONTEXT = (
    r"cell membrane|plasma membrane|(?:endoplasmic|sarcoplasmic) reticulum(?: membrane)?|"
    r"golgi apparatus(?: membrane)?|endosome(?: membrane)?|lysosome(?: membrane)?|"
    r"nucle(?:us|ar) envelope|mitochondri(?:on|al) (?:inner |outer )?membrane|"
    r"peroxisom(?:e|al) membrane|lipid droplet|autophagosome membrane|vesicle membrane|"
    r"membrane raft|postsynaptic membrane|presynaptic membrane"
)


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


def has_experimental_membrane_location(text: str) -> bool:
    return bool(
        re.search(
            rf"(?is)(?:{MEMBRANE_CONTEXT})[^.;]{{0,220}}ECO:0000269",
            text or "",
        )
    )


def has_experimental_lipid_anchor(text: str) -> bool:
    return bool(
        re.search(r"(?is)Lipid-anchor[^.;]{0,160}ECO:0000269", text or "")
        or re.search(r"(?is)ECO:0000269[^.;]{0,160}Lipid-anchor", text or "")
    )


def assign_inherited(row: dict[str, str]) -> tuple[str, str, str]:
    cls = row["release_tier_v5"]
    basis = row.get("membrane_evidence_basis", "")
    subcell = row.get("subcellular_location", "")
    transmem = row.get("transmembrane_features", "")
    intramem = row.get("intramembrane_features", "")
    hpa_reliability = row.get("hpa_reliability_if_v5", "")
    hpa_plasma = row.get("hpa_plasma_membrane_location_v5", "") == "1"
    external_count = int(row.get("independent_membrane_source_count_v5", "") or 0)
    external_sources = row.get("independent_membrane_sources_v5", "")
    opm = row.get("opm_present_v5", "") == "1"
    membranome = row.get("membranome_present_v5", "") == "1"

    direct_topology = "ECO:0000269" in transmem or "ECO:0000269" in intramem
    direct_membrane_location = has_experimental_membrane_location(subcell)
    direct_lipid_anchor = has_experimental_lipid_anchor(subcell)

    if cls == "A" and direct_topology:
        return (
            "E1",
            "E1_experimental_topology_feature",
            "UniProt TM/intramembrane feature carries ECO:0000269 experimental evidence.",
        )
    if cls == "B" and direct_lipid_anchor:
        return (
            "E1",
            "E1_experimental_lipid_anchor",
            "UniProt reports an experimentally supported lipid-anchor annotation.",
        )
    if direct_membrane_location:
        return (
            "E1",
            "E1_experimental_membrane_location",
            "UniProt reports membrane or membrane-organelle localization with ECO:0000269 experimental evidence.",
        )

    curated_topology = basis.startswith("uniprot_") and cls in {"A", "B"}
    reliable_hpa = hpa_plasma and hpa_reliability in {"Enhanced", "Supported"}
    multi_source = external_count >= 2
    structural_or_family_support = opm or membranome
    no_evidence_codes = "ECO:" not in subcell

    if curated_topology:
        return (
            "E2",
            "E2_curated_UniProt_topology_or_anchor",
            "UniProt provides a curated TM/intramembrane/lipid-anchor classification, but direct experimental evidence is not explicit in the exported evidence text.",
        )
    if reliable_hpa:
        return (
            "E2",
            "E2_reliable_HPA_membrane_location",
            f"HPA plasma-membrane localization has {hpa_reliability} reliability.",
        )
    if structural_or_family_support:
        source = ";".join(
            value
            for value, present in (("OPM", opm), ("Membranome", membranome))
            if present
        )
        return (
            "E2",
            "E2_structure_or_membrane_family_support",
            f"Strong membrane support from {source}; direct experimental evidence for the reviewed protein is not explicit.",
        )
    if multi_source:
        return (
            "E2",
            "E2_multiple_independent_sources",
            f"At least two independent membrane sources agree: {external_sources}.",
        )
    if cls == "C" and basis == "uniprot_membrane_location_annotation" and no_evidence_codes:
        return (
            "E2",
            "E2_curated_UniProt_membrane_location",
            "UniProt provides a curated membrane-location statement without an explicit computational/by-similarity evidence code.",
        )

    return (
        "E3",
        "E3_single_or_inferred_support",
        "Membrane assignment is supported only by a single source, prediction, by-similarity statement, or limited localization evidence.",
    )


def assign_tier(row: dict[str, str]) -> tuple[str, str, str]:
    if row["release_tier_v5"] == "R":
        status = row["evidence_status_v51"]
        code = row["decision_code_v51"]
        basis = row["decision_basis_v51"]
        if status == "probable" and code in {
            "B_probable_experimental_myristoylation_plus_HPA",
            "C_probable_experimental_organelle_membrane_context",
        }:
            return (
                "E1",
                f"E1_former_R_{code}",
                basis
                + " Reclassified to E1 because the underlying UniProt evidence carries direct experimental support.",
            )
        mapping = {
            "confirmed": ("E1", f"E1_former_R_{code}", basis),
            "probable": ("E2", f"E2_former_R_{code}", basis),
            "uncertain": ("E3", f"E3_former_R_{code}", basis),
            "excluded": ("E0", f"E0_former_R_{code}", basis),
        }
        return mapping[status]
    return assign_inherited(row)


def main() -> None:
    source_rows = read_tsv(SOURCE_MASTER)
    if len(source_rows) != 10997:
        raise RuntimeError("Expected 10,997 rows in the v5.1 audit master")
    if RELEASE.exists():
        shutil.rmtree(RELEASE)
    for name in ("reports", "reproducibility/scripts", "review"):
        (RELEASE / name).mkdir(parents=True, exist_ok=True)

    fields = list(source_rows[0]) + [
        field for field in NEW_FIELDS if field not in source_rows[0]
    ]
    rows: list[dict[str, str]] = []
    for source in source_rows:
        row = dict(source)
        level, rule, evidence_basis = assign_tier(row)
        zh, en, scope, direct_flag, website_default = LABELS[level]
        if level == "E0":
            membrane_class = "unknown"
        elif row["final_membrane_class_v51"] in {"A", "B", "C"}:
            membrane_class = row["final_membrane_class_v51"]
        else:
            membrane_class = row["candidate_membrane_class_v51"]
        row.update(
            {
                "membrane_class_v52": membrane_class,
                "evidence_level_v52": level,
                "evidence_label_zh_v52": zh,
                "evidence_label_en_v52": en,
                "direct_experimental_flag_v52": direct_flag,
                "release_scope_v52": scope,
                "website_default_v52": website_default,
                "evidence_rule_v52": rule,
                "evidence_basis_v52": evidence_basis,
            }
        )
        rows.append(row)

    e1 = [row for row in rows if row["evidence_level_v52"] == "E1"]
    e2 = [row for row in rows if row["evidence_level_v52"] == "E2"]
    e3 = [row for row in rows if row["evidence_level_v52"] == "E3"]
    e0 = [row for row in rows if row["evidence_level_v52"] == "E0"]
    default = e1 + e2
    class_views = {
        cls: [row for row in default if row["membrane_class_v52"] == cls]
        for cls in ("A", "B", "C")
    }

    write_tsv(RELEASE / "human_membrane_audit_master_v5_2.tsv", rows, fields)
    write_tsv(RELEASE / "human_membrane_experimental_core_E1_v5_2.tsv", e1, fields)
    write_tsv(RELEASE / "human_membrane_strong_supported_E2_v5_2.tsv", e2, fields)
    write_tsv(RELEASE / "human_membrane_default_E1_E2_v5_2.tsv", default, fields)
    write_tsv(RELEASE / "human_membrane_prediction_candidates_E3_v5_2.tsv", e3, fields)
    write_tsv(RELEASE / "human_membrane_excluded_E0_v5_2.tsv", e0, fields)
    for cls, view in class_views.items():
        write_tsv(RELEASE / f"human_membrane_class_{cls}_E1_E2_v5_2.tsv", view, fields)

    stats = {
        "release": "v5.2-evidence-tiered",
        "release_date": RELEASE_DATE,
        "audit_rows": len(rows),
        "evidence_level_counts": dict(
            Counter(row["evidence_level_v52"] for row in rows)
        ),
        "evidence_level_by_class": {
            level: dict(
                Counter(
                    row["membrane_class_v52"]
                    for row in rows
                    if row["evidence_level_v52"] == level
                )
            )
            for level in ("E1", "E2", "E3", "E0")
        },
        "default_E1_E2_rows": len(default),
        "default_E1_E2_class_counts": {
            cls: len(view) for cls, view in class_views.items()
        },
        "former_R_evidence_level_counts": dict(
            Counter(
                row["evidence_level_v52"]
                for row in rows
                if row["release_tier_v5"] == "R"
            )
        ),
        "rule_counts": dict(Counter(row["evidence_rule_v52"] for row in rows)),
    }
    (RELEASE / "reports" / "V52_EVIDENCE_TIER_STATS.json").write_text(
        json.dumps(stats, indent=2, ensure_ascii=False, sort_keys=True),
        encoding="utf-8",
    )

    (RELEASE / "README_v5_2.md").write_text(
        f"""# Human Membrane Protein Master v5.2 — evidence-tiered

This release separates biological membrane class from evidence strength.

## Biological class

- **A** — integral membrane protein.
- **B** — directly membrane-inserted without a TM span, chiefly lipid-anchored.
- **C** — peripheral or stable membrane-associated protein.

## Evidence level

- **E1 / 实验证实** — direct experimental topology, lipid-anchor, or membrane-location evidence.
- **E2 / 强证据支持** — curated UniProt topology/location, reliable HPA localization, structure/family support, or multiple independent sources, without explicit direct evidence in the available record.
- **E3 / 预测候选** — prediction, single-source, by-similarity, or limited localization evidence.
- **E0 / 排除** — not retained as a membrane protein.

The strict experimental core contains {len(e1):,} records. The recommended website default contains E1+E2 ({len(default):,} records). E3 contains {len(e3):,} separately searchable candidates, and E0 contains {len(e0):,} internal exclusion records.
""",
        encoding="utf-8",
    )
    (RELEASE / "EVIDENCE_POLICY_v5_2.md").write_text(
        """# Evidence policy v5.2

The fields `membrane_class_v52` and `evidence_level_v52` answer different questions. A/B/C describes how a protein associates with membrane; E1/E2/E3/E0 describes confidence.

## E1

E1 requires an explicit ECO:0000269 experimental qualifier attached to a UniProt TM/intramembrane feature, lipid-anchor statement, or membrane/membrane-organelle localization. Former-R records previously marked `confirmed` also map to E1.

## E2

E2 includes curated UniProt TM/intramembrane/lipid-anchor assignments without explicit direct evidence in the exported text; HPA Enhanced/Supported plasma-membrane localization; OPM/Membranome support; at least two independent external sources; and curated UniProt membrane-location statements without computational/by-similarity evidence codes. Former-R `probable` records map to E2.

## E3

E3 contains single-source predictions, by-similarity statements, HPA Approved-only localization, and unresolved signal-anchor/signal-peptide cases. It is not included in default statistics.

## E0

E0 contains secreted/signal-peptide false positives, non-reproduced predictions, target-chain structural mismatches, and other excluded records. It is retained only for audit and re-import prevention.
""",
        encoding="utf-8",
    )
    (RELEASE / "DATA_DICTIONARY_v5_2.md").write_text(
        """# v5.2 added fields

| Field | Meaning |
|---|---|
| `membrane_class_v52` | A/B/C biological class; `unknown` for E0. |
| `evidence_level_v52` | E1, E2, E3, or E0. |
| `evidence_label_zh_v52` | Chinese evidence label. |
| `evidence_label_en_v52` | English evidence label. |
| `direct_experimental_flag_v52` | `1` only for E1. |
| `release_scope_v52` | core, supported, candidate, or excluded. |
| `website_default_v52` | `1` for E1+E2 default website display. |
| `evidence_rule_v52` | Machine-readable rule used for tiering. |
| `evidence_basis_v52` | Human-readable rationale. |
""",
        encoding="utf-8",
    )
    (RELEASE / "SOURCE_REGISTRY_v5_2.tsv").write_text(
        (SOURCE_RELEASE / "SOURCE_REGISTRY_v5_1.tsv").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    shutil.copy2(
        ROOT / "scripts" / "build_v52_evidence_tiers.py",
        RELEASE / "reproducibility" / "scripts" / "build_v52_evidence_tiers.py",
    )
    shutil.copy2(
        SOURCE_RELEASE / "review" / "r1109_final_review_v51.tsv",
        RELEASE / "review" / "r1109_final_review_v51.tsv",
    )
    print(json.dumps(stats, indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
