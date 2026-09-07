from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "v51_review_working"
PRELIM = WORK / "analysis" / "r1109_preliminary_review_v51.tsv"
OVERRIDES = ROOT / "configs" / "v51_review_overrides.tsv"
OUTPUT = WORK / "analysis" / "r1109_final_review_v51.tsv"
STATS = WORK / "reports" / "V51_FINAL_REVIEW_STATS.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def segment_lengths(value: str) -> list[int]:
    lengths: list[int] = []
    for start, end in re.findall(r"(\d+)-(\d+)", value or ""):
        lengths.append(int(end) - int(start) + 1)
    return lengths


def set_decision(
    row: dict[str, str],
    candidate: str,
    final: str,
    status: str,
    inclusion: str,
    code: str,
    basis: str,
    layer: str,
) -> None:
    row["candidate_membrane_class_v51"] = candidate
    row["final_membrane_class_v51"] = final
    row["evidence_status_v51"] = status
    row["inclusion_status_v51"] = inclusion
    row["decision_code_v51"] = code
    row["decision_basis_v51"] = basis
    row["v51_decision_layer"] = layer


def refine(row: dict[str, str]) -> dict[str, str]:
    reason = row["v5_review_reason"]
    status = row["evidence_status_v51"]
    subcell = row["uniprot_subcellular_location"]
    max_tmbed = max(segment_lengths(row["tmbed_tm_segments"]) or [0])
    max_hydro = max(segment_lengths(row["hydropathy_non_signal_segments"]) or [0])
    tm_count = int(row["tmbed_tm_segment_count"] or 0)
    htp_count = int(row["htp_num_tm"] or 0)
    exact_tmbed = row["tmbed_sequence_exact"] == "1"
    hpa_rel = row["hpa_reliability"]
    secreted = bool(re.search(r"(?i)\bsecreted\b|extracellular space", subcell))
    context_pattern = (
        r"endoplasmic reticulum|golgi apparatus|endosome|lipid droplet|"
        r"nucle(?:us|ar) envelope|cell membrane|plasma membrane|"
        r"cytoplasmic vesicle|clathrin-coated vesicle"
    )
    membrane_context = bool(re.search(rf"(?i){context_pattern}", subcell))
    # Keep the evidence qualifier tied to the membrane-context clause.  A
    # protein may have experimental nuclear localization and only a
    # by-similarity Golgi annotation elsewhere in the same UniProt comment.
    experimental_context = bool(
        re.search(rf"(?i)(?:{context_pattern})[^.;]{{0,180}}ECO:0000269", subcell)
    )
    lumen_context = bool(re.search(r"(?i)\blumen\b", subcell))

    row["v51_decision_layer"] = "preliminary_rule_retained"
    row["manual_review_flag_v51"] = "0"

    # An experimentally supported non-lumenal membrane-organelle localization
    # can rescue a false integral-protein prediction as peripheral class C.
    if (
        not secreted
        and not lumen_context
        and experimental_context
        and membrane_context
        and status in {"uncertain", "excluded"}
    ):
        set_decision(
            row,
            "C",
            "C",
            "probable",
            "included",
            "C_probable_experimental_organelle_membrane_context",
            "Current UniProt reports experimentally supported localization to a membrane-bounded organelle or membrane surface; no integral topology is reproduced, so the protein is assigned to peripheral class C.",
            "v51_refined_rule",
        )

    if reason == "new_hpa_location_only":
        if secreted:
            return row
        if hpa_rel == "Uncertain" and row["final_membrane_class_v51"] == "unknown":
            set_decision(
                row,
                "C",
                "unknown",
                "excluded",
                "excluded",
                "excluded_HPA_uncertain_location_only",
                "The only membrane evidence is an HPA localization explicitly rated Uncertain.",
                "v51_refined_rule",
            )
        return row

    if reason in {
        "new_hpa_predicted_membrane_only",
        "new_htp_entry_without_experimental_or_membranome_support",
    }:
        if row["final_membrane_class_v51"] == "C":
            return row

        reproduced_full_tm = exact_tmbed and max_tmbed >= 15
        strong_sequence = reproduced_full_tm or max_hydro >= 19
        multi_segment = (tm_count >= 2 and max_tmbed >= 15) or htp_count >= 2

        if row["final_membrane_class_v51"] == "A":
            if max_tmbed < 15 and max_hydro < 19:
                set_decision(
                    row,
                    "A",
                    "unknown",
                    "uncertain",
                    "candidate",
                    "A_candidate_short_external_TM_call",
                    "The external predictors agree, but the reproduced TMbed segment is shorter than 15 residues and no full hydrophobic window is present.",
                    "v51_refined_rule",
                )
            elif (
                max_tmbed >= 15
                and row["tmbed_tm_segments"].split(";")[0].startswith(
                    tuple(str(i) + "-" for i in range(1, 26))
                )
                and tm_count == 1
                and not membrane_context
                and not multi_segment
            ):
                set_decision(
                    row,
                    "A",
                    "unknown",
                    "uncertain",
                    "candidate",
                    "A_candidate_N_terminal_signal_anchor_ambiguity",
                    "A single N-terminal hydrophobic segment is reproduced, but no independent localization evidence distinguishes a signal anchor from a cleaved signal peptide.",
                    "v51_refined_rule",
                )
            return row

        if strong_sequence:
            set_decision(
                row,
                "A",
                "unknown",
                "uncertain",
                "candidate",
                "A_candidate_external_prediction_plus_sequence_support",
                "An external membrane prediction has hydrophobic sequence support, but current curated topology and independent experimental membrane evidence are absent.",
                "v51_refined_rule",
            )
        else:
            set_decision(
                row,
                "A",
                "unknown",
                "excluded",
                "excluded",
                "excluded_external_TM_prediction_not_reproduced",
                "The legacy HPA/UniTmp integral-membrane prediction is not reproduced by current-sequence TMbed or a full hydrophobic window, and current UniProt has no curated TM feature.",
                "v51_refined_rule",
            )
        return row

    if reason == "new_opm_membrane_association_without_integral_segments":
        if row["evidence_status_v51"] == "uncertain":
            set_decision(
                row,
                "C",
                "unknown",
                "excluded",
                "excluded",
                "excluded_OPM_undefined_or_noncellular_context",
                "The OPM entry does not place the target at a defined non-secreted cellular membrane, and the target chain is non-transmembrane.",
                "v51_refined_rule",
            )
        return row

    return row


def main() -> None:
    rows = [refine(dict(row)) for row in read_tsv(PRELIM)]
    overrides = {row["target_uniprot_id"]: row for row in read_tsv(OVERRIDES)}
    seen: set[str] = set()
    for row in rows:
        accession = row["target_uniprot_id"]
        override = overrides.get(accession)
        if not override:
            continue
        seen.add(accession)
        for field in (
            "candidate_membrane_class_v51",
            "final_membrane_class_v51",
            "evidence_status_v51",
            "inclusion_status_v51",
            "decision_code_v51",
        ):
            row[field] = override[field]
        row["decision_basis_v51"] = override["override_rationale"]
        row["v51_decision_layer"] = "targeted_manual_override"
        row["manual_review_flag_v51"] = "1"

    missing = sorted(set(overrides) - seen)
    if missing:
        raise RuntimeError(f"Override accessions absent from review: {missing}")

    fieldnames = list(rows[0])
    for field in ("v51_decision_layer", "manual_review_flag_v51"):
        if field not in fieldnames:
            fieldnames.append(field)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

    stats = {
        "review_rows": len(rows),
        "candidate_class_counts": Counter(
            row["candidate_membrane_class_v51"] for row in rows
        ),
        "final_class_counts": Counter(row["final_membrane_class_v51"] for row in rows),
        "evidence_status_counts": Counter(
            row["evidence_status_v51"] for row in rows
        ),
        "inclusion_status_counts": Counter(
            row["inclusion_status_v51"] for row in rows
        ),
        "decision_code_counts": Counter(row["decision_code_v51"] for row in rows),
        "review_reason_outcomes": Counter(
            f'{row["v5_review_reason"]}|{row["evidence_status_v51"]}|{row["inclusion_status_v51"]}'
            for row in rows
        ),
        "manual_override_rows": sum(
            row["manual_review_flag_v51"] == "1" for row in rows
        ),
        "pdbtm_target_chain_tm_rows": sum(
            row["pdbtm_target_chain_tm_flag"] == "1" for row in rows
        ),
    }
    STATS.parent.mkdir(parents=True, exist_ok=True)
    STATS.write_text(
        json.dumps(stats, indent=2, ensure_ascii=False, sort_keys=True),
        encoding="utf-8",
    )
    print(f"Wrote {OUTPUT} ({len(rows)} rows)")
    print(json.dumps(stats, indent=2, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
