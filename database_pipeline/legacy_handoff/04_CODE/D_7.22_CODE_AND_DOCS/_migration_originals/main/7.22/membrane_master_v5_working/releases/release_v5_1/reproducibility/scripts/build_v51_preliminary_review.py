from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


ROOT = Path(r"C:\Users\Administrator\7.22\membrane_master_v5_working")
V51 = ROOT / "v51_review_working"
RAW = V51 / "raw"
NORMALIZED = V51 / "normalized"
ANALYSIS = V51 / "analysis"
REPORTS = V51 / "reports"

QUEUE = ROOT / "review" / "human_membrane_review_queue_v5.tsv"
UNIPROT = RAW / "uniprot_r_review_2026-07-24.tsv"
TM_BED = RAW / "tmbed_human_2022_predictions.txt"
GO_OBO = RAW / "go_basic_2026-06-26.obo"
GOA = RAW / "goa_human_2026-07-08.gaf.gz"
INTERPRO = RAW / "interpro_entry_list_2026-06-10.tsv"
HPA = ROOT / "normalized" / "hpa_membrane_evidence_v5.tsv"
HTP = ROOT / "normalized" / "unitmp_htp_human_topology_v5.tsv"
OPM = ROOT / "normalized" / "opm_human_structures_v5.tsv"
PDBTM_CHAINS = NORMALIZED / "pdbtm_target_chain_topology_v51.tsv"

EXPERIMENTAL_GO_CODES = {
    "EXP",
    "IDA",
    "IPI",
    "IMP",
    "IGI",
    "IEP",
    "HTP",
    "HDA",
    "HMP",
    "HGI",
    "HEP",
}

HYDROPATHY = {
    "I": 4.5,
    "V": 4.2,
    "L": 3.8,
    "F": 2.8,
    "C": 2.5,
    "M": 1.9,
    "A": 1.8,
    "G": -0.4,
    "T": -0.7,
    "S": -0.8,
    "W": -0.9,
    "Y": -1.3,
    "P": -1.6,
    "H": -3.2,
    "E": -3.5,
    "Q": -3.5,
    "D": -3.5,
    "N": -3.5,
    "K": -3.9,
    "R": -4.5,
}

DIRECT_MEMBRANE_TEXT = re.compile(
    r"(?i)("
    r"(binds?|binding|associates?|association|recruited|recruitment|anchors?|"
    r"anchoring|targets?|targeting|inserts?|insertion).{0,90}"
    r"(cell membrane|plasma membrane|organelle membrane|lipid bilayer|"
    r"phospholipid|phosphoinositide)"
    r"|"
    r"(cell membrane|plasma membrane|organelle membrane|lipid bilayer|"
    r"phospholipid|phosphoinositide).{0,90}"
    r"(binds?|binding|associates?|association|recruited|recruitment|anchors?|"
    r"anchoring|targets?|targeting|inserts?|insertion)"
    r")"
)

MEMBRANE_BINDING_DOMAIN_TEXT = re.compile(
    r"(?i)("
    r"pleckstrin homology domain|"
    r"\bC2 domain\b|"
    r"\bPX domain\b|"
    r"\bFYVE domain\b|"
    r"\bF-BAR domain\b|"
    r"\bI-BAR domain\b|"
    r"\bBAR domain\b|"
    r"\bFERM domain\b|"
    r"\bENTH domain\b|"
    r"\bGRAM domain\b|"
    r"annexin|"
    r"CRAL-TRIO|"
    r"SEC14-like"
    r")"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_obo(path: Path) -> dict[str, dict[str, object]]:
    terms: dict[str, dict[str, object]] = {}
    current: dict[str, object] | None = None
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "[Term]":
                current = {
                    "id": "",
                    "name": "",
                    "namespace": "",
                    "parents": [],
                    "obsolete": False,
                }
                continue
            if line.startswith("[") and line != "[Term]":
                current = None
                continue
            if current is None:
                continue
            if line.startswith("id: "):
                current["id"] = line[4:]
            elif line.startswith("name: "):
                current["name"] = line[6:]
            elif line.startswith("namespace: "):
                current["namespace"] = line[11:]
            elif line.startswith("is_a: "):
                current["parents"].append(line[6:].split()[0])
            elif line == "is_obsolete: true":
                current["obsolete"] = True
            elif line == "" and current["id"]:
                terms[str(current["id"])] = current
                current = None
    if current and current["id"]:
        terms[str(current["id"])] = current
    return terms


def descendants(
    root_id: str,
    terms: dict[str, dict[str, object]],
) -> set[str]:
    children: dict[str, set[str]] = defaultdict(set)
    for term_id, term in terms.items():
        for parent in term["parents"]:
            children[str(parent)].add(term_id)
    result: set[str] = set()
    stack = [root_id]
    while stack:
        term_id = stack.pop()
        if term_id in result:
            continue
        result.add(term_id)
        stack.extend(children.get(term_id, set()))
    return result


def parse_goa(
    path: Path,
    accessions: set[str],
    terms: dict[str, dict[str, object]],
) -> dict[str, list[dict[str, str]]]:
    annotations: dict[str, list[dict[str, str]]] = defaultdict(list)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 15:
                continue
            accession = fields[1].strip().upper()
            if accession not in accessions or fields[8] != "C":
                continue
            qualifiers = fields[3].split("|")
            if "NOT" in qualifiers:
                continue
            go_id = fields[4]
            if go_id not in terms:
                continue
            annotations[accession].append(
                {
                    "go_id": go_id,
                    "go_name": str(terms[go_id]["name"]),
                    "evidence_code": fields[6],
                    "reference": fields[5],
                    "qualifier": fields[3],
                    "annotation_date": fields[13],
                    "assigned_by": fields[14],
                }
            )
    return annotations


def contiguous_segments(prediction: str, pattern: str) -> list[tuple[int, int, str]]:
    segments = []
    for match in re.finditer(pattern, prediction):
        segments.append((match.start() + 1, match.end(), match.group(0)[0]))
    return segments


def parse_tmbed(path: Path, relevant_accessions: set[str]) -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}
    with path.open("r", encoding="utf-8") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip()
            prediction = handle.readline().strip()
            match = re.match(r">sp\|([^|]+)\|", header)
            if not match:
                continue
            accession = match.group(1).upper()
            if accession not in relevant_accessions:
                continue
            tm_segments = contiguous_segments(prediction, r"[Hh]+")
            beta_segments = contiguous_segments(prediction, r"[Bb]+")
            signal_segments = contiguous_segments(prediction, r"S+")
            result[accession] = {
                "sequence": sequence,
                "prediction": prediction,
                "tm_segments": tm_segments,
                "beta_segments": beta_segments,
                "signal_segments": signal_segments,
            }
    return result


def parse_signal_end(value: str) -> int:
    match = re.search(r"\bSIGNAL\s+\d+\.\.(\d+)", value or "")
    return int(match.group(1)) if match else 0


def hydrophobic_segments(sequence: str, window: int = 19) -> list[tuple[int, int, float]]:
    if len(sequence) < window:
        return []
    qualifying = []
    for start in range(len(sequence) - window + 1):
        peptide = sequence[start : start + window]
        score = sum(HYDROPATHY.get(residue, 0.0) for residue in peptide) / window
        charged = sum(residue in "DEKR" for residue in peptide)
        prolines = peptide.count("P")
        if score >= 1.8 and charged <= 2 and prolines <= 2:
            qualifying.append((start + 1, start + window, score))
    if not qualifying:
        return []
    merged: list[list[float]] = []
    for start, end, score in qualifying:
        if not merged or start > merged[-1][1] + 2:
            merged.append([start, end, score])
        else:
            merged[-1][1] = max(merged[-1][1], end)
            merged[-1][2] = max(merged[-1][2], score)
    return [(int(start), int(end), round(score, 3)) for start, end, score in merged]


def segment_text(segments: list[tuple[int, int, object]]) -> str:
    return ";".join(f"{start}-{end}" for start, end, _ in segments)


def pubmed_ids(text: str) -> list[str]:
    return sorted(set(re.findall(r"PubMed:(\d+)", text or "")), key=int)


def classify_record(
    base: dict[str, str],
    uniprot: dict[str, str],
    hpa: dict[str, str],
    htp: dict[str, str],
    opm_rows: list[dict[str, str]],
    chain_rows: list[dict[str, str]],
    go_annotations: list[dict[str, str]],
    membrane_go_terms: set[str],
    interpro_names: dict[str, str],
    tmbed: dict[str, object] | None,
) -> dict[str, object]:
    accession = base["target_uniprot_id"].upper()
    sequence = uniprot.get("Sequence", "")
    signal_feature = uniprot.get("Signal peptide", "") or base.get(
        "signal_peptide_features", ""
    )
    signal_end = parse_signal_end(signal_feature)
    subcellular = uniprot.get("Subcellular location [CC]", "")
    keywords = uniprot.get("Keywords", "")
    function = uniprot.get("Function [CC]", "")
    domain_comment = uniprot.get("Domain [CC]", "")
    ptm_comment = uniprot.get("Post-translational modification", "")
    go_function = uniprot.get("Gene Ontology (molecular function)", "")
    lipidation = uniprot.get("Lipidation", "")

    secreted = bool(
        re.search(r"(?i)\bSecreted\b|extracellular space|extracellular region", subcellular)
        or re.search(r"(?i)(^|;)Secreted(;|$)", keywords)
    )
    lumen_only = bool(
        re.search(r"(?i)\blumen\b", subcellular)
        and not re.search(r"(?i)cytoplasm|cytosol|cell membrane|plasma membrane", subcellular)
    )

    hydro_segments = hydrophobic_segments(sequence)
    hydro_after_signal = [
        segment
        for segment in hydro_segments
        if segment[0] > signal_end + 5 or signal_end == 0
    ]

    tmbed_available = bool(tmbed)
    tmbed_sequence_exact = bool(tmbed and tmbed["sequence"] == sequence)
    tm_segments = list(tmbed["tm_segments"]) if tmbed_sequence_exact else []
    beta_segments = list(tmbed["beta_segments"]) if tmbed_sequence_exact else []
    sp_segments = list(tmbed["signal_segments"]) if tmbed_sequence_exact else []
    tmbed_tm = bool(tm_segments or beta_segments)
    tmbed_sp_only = bool(sp_segments and not tmbed_tm)

    htp_present = bool(htp)
    htp_evidence = htp.get("htp_evidence", "")
    htp_tm_count = htp.get("htp_num_tm", "")
    try:
        htp_first_tm_start = int(htp.get("htp_first_tm_start", "") or 0)
    except ValueError:
        htp_first_tm_start = 0
    htp_overlaps_signal = bool(
        signal_end
        and htp_first_tm_start
        and htp_first_tm_start <= signal_end + 5
    )

    hpa_main = hpa.get("hpa_subcellular_main_location", "")
    hpa_additional = hpa.get("hpa_subcellular_additional_location", "")
    hpa_reliability = hpa.get("hpa_reliability_if", "") or base.get(
        "hpa_reliability_if_v5", ""
    )
    hpa_main_plasma = "Plasma membrane" in hpa_main
    hpa_additional_plasma = "Plasma membrane" in hpa_additional
    hpa_predicted = base.get("hpa_predicted_membrane_v5", "") == "1"

    membrane_go = [
        annotation
        for annotation in go_annotations
        if annotation["go_id"] in membrane_go_terms
    ]
    experimental_membrane_go = [
        annotation
        for annotation in membrane_go
        if annotation["evidence_code"] in EXPERIMENTAL_GO_CODES
    ]

    interpro_ids = [
        value.strip()
        for value in uniprot.get("InterPro", "").split(";")
        if value.strip()
    ]
    matching_domains = sorted(
        {
            f"{interpro_id}:{interpro_names[interpro_id]}"
            for interpro_id in interpro_ids
            if interpro_id in interpro_names
            and MEMBRANE_BINDING_DOMAIN_TEXT.search(interpro_names[interpro_id])
        }
    )
    direct_text_blob = " ".join(
        [function, domain_comment, ptm_comment, go_function]
    )
    direct_membrane_text = bool(DIRECT_MEMBRANE_TEXT.search(direct_text_blob))

    anchor_type = ""
    if re.search(r"(?i)GPI-anchor", lipidation):
        anchor_type = "GPI_anchor"
    elif re.search(r"(?i)S-geranylgeranyl|S-farnesyl", lipidation):
        anchor_type = "prenyl_anchor"
    elif re.search(r"(?i)N-myristoyl", lipidation):
        anchor_type = "myristoyl_anchor"
    elif re.search(r"(?i)S-palmitoyl", lipidation):
        anchor_type = "palmitoyl_anchor"

    opm_membranes = sorted(
        {row.get("opm_membrane", "") for row in opm_rows if row.get("opm_membrane", "")}
    )
    opm_nonsecreted_context = any(
        value not in {"Secreted", "Granule", "Undefined"}
        for value in opm_membranes
    )
    target_chain_tm = any(row.get("pdbtm_target_has_tm") == "1" for row in chain_rows)
    target_chain_non_tm = any(
        row.get("pdbtm_target_all_non_tm") == "1" for row in chain_rows
    )

    evidence: list[str] = []
    warnings: list[str] = []
    candidate_class = ""
    final_class = "unknown"
    evidence_status = "uncertain"
    inclusion_status = "candidate"
    decision_code = ""

    # No R entry should be promoted by a PDBTM record unless the reviewed
    # accession maps to the transmembrane chain itself.
    if target_chain_tm:
        candidate_class = "A"
        final_class = "A"
        evidence_status = "confirmed"
        inclusion_status = "included"
        decision_code = "A_confirmed_target_pdbtm_chain"
        evidence.append("target accession maps to a PDBTM chain with TM segment(s)")

    if not decision_code and tmbed_tm and htp_present and not secreted:
        candidate_class = "A"
        final_class = "A"
        evidence_status = "probable"
        inclusion_status = "included"
        decision_code = "A_probable_TUniTmp_TTMBed_consensus"
        evidence.extend(
            [
                "UniTmp HTP predicts transmembrane topology",
                "TMbed predicts a transmembrane segment on an identical sequence",
            ]
        )
        if hydro_after_signal:
            evidence.append("independent hydropathy window supports a non-signal TM segment")

    if (
        not decision_code
        and tmbed_tm
        and hpa_predicted
        and not secreted
    ):
        candidate_class = "A"
        final_class = "A"
        evidence_status = "probable"
        inclusion_status = "included"
        decision_code = "A_probable_HPA_TTMBed_consensus"
        evidence.extend(
            [
                "HPA predicted-membrane class",
                "TMbed predicts a transmembrane segment on an identical sequence",
            ]
        )

    if (
        not decision_code
        and anchor_type
        and (hpa_main_plasma or hpa_additional_plasma or membrane_go or opm_rows)
        and not secreted
    ):
        candidate_class = "B"
        final_class = "B"
        evidence_status = "probable"
        inclusion_status = "included"
        decision_code = "B_probable_lipid_anchor_plus_membrane_evidence"
        evidence.extend(
            [
                f"UniProt lipidation feature: {anchor_type}",
                "independent membrane localization/structure evidence",
            ]
        )

    if not decision_code and experimental_membrane_go:
        candidate_class = "C"
        final_class = "C"
        evidence_status = "confirmed"
        inclusion_status = "included"
        decision_code = "C_confirmed_experimental_GOA_membrane_location"
        evidence.append("experimental GO cellular-component membrane annotation")

    if (
        not decision_code
        and opm_nonsecreted_context
        and not secreted
        and not lumen_only
    ):
        candidate_class = "C"
        final_class = "C"
        evidence_status = "probable"
        inclusion_status = "included"
        decision_code = "C_probable_OPM_peripheral_structure"
        evidence.append(
            "OPM places the non-integral target structure at a non-secreted cellular membrane"
        )
        if target_chain_non_tm:
            evidence.append("PDBTM confirms the target chain itself is non-transmembrane")

    if (
        not decision_code
        and (hpa_main_plasma or hpa_additional_plasma)
        and hpa_reliability in {"Enhanced", "Supported"}
        and not secreted
    ):
        candidate_class = "C"
        final_class = "C"
        evidence_status = "probable"
        inclusion_status = "included"
        decision_code = "C_probable_HPA_validated_plasma_membrane_location"
        evidence.append(
            f"HPA plasma-membrane localization with {hpa_reliability} reliability"
        )
        if hpa_additional_plasma and not hpa_main_plasma:
            warnings.append("plasma membrane is an additional rather than main HPA location")

    if (
        not decision_code
        and direct_membrane_text
        and (hpa_main_plasma or hpa_additional_plasma or opm_rows or membrane_go)
        and not secreted
    ):
        candidate_class = "C"
        final_class = "C"
        evidence_status = "probable"
        inclusion_status = "included"
        decision_code = "C_probable_direct_membrane_binding_text"
        evidence.extend(
            [
                "UniProt function/domain text describes direct membrane or phospholipid interaction",
                "independent membrane source is present",
            ]
        )

    signal_confusion = bool(
        signal_end
        and not tmbed_tm
        and (
            tmbed_sp_only
            or htp_overlaps_signal
            or (
                htp_present
                and hydro_segments
                and not hydro_after_signal
            )
        )
    )
    if (
        not decision_code
        and secreted
        and not experimental_membrane_go
        and not (
            direct_membrane_text
            and (hpa_main_plasma or hpa_additional_plasma)
            and hpa_reliability in {"Enhanced", "Supported"}
        )
    ):
        candidate_class = (
            "A"
            if base["membrane_decision_v5"]
            in {
                "new_hpa_predicted_membrane_only",
                "new_htp_entry_without_experimental_or_membranome_support",
            }
            else "C"
        )
        final_class = "unknown"
        evidence_status = "excluded"
        inclusion_status = "excluded"
        decision_code = (
            "excluded_signal_peptide_or_secreted_soluble"
            if signal_confusion or signal_end
            else "excluded_secreted_without_stable_membrane_association"
        )
        evidence.append("UniProt describes the protein as secreted/extracellular")
        if signal_confusion:
            warnings.append("external TM prediction is consistent with a cleaved signal peptide")

    if (
        not decision_code
        and base["membrane_decision_v5"] == "new_hpa_location_only"
    ):
        candidate_class = "C"
        final_class = "unknown"
        evidence_status = "uncertain"
        inclusion_status = "candidate"
        if hpa_main_plasma:
            decision_code = "C_candidate_HPA_main_location_only"
            evidence.append("HPA reports plasma membrane as a main location")
        else:
            decision_code = "C_candidate_HPA_additional_location_only"
            evidence.append("HPA reports plasma membrane only as an additional location")
        evidence.append(f"HPA reliability: {hpa_reliability or 'not available'}")

    if (
        not decision_code
        and base["membrane_decision_v5"] == "new_hpa_predicted_membrane_only"
    ):
        candidate_class = "A"
        final_class = "unknown"
        evidence_status = "uncertain"
        inclusion_status = "candidate"
        decision_code = "A_candidate_HPA_prediction_without_sequence_support"
        evidence.append("HPA predicted-membrane class")
        warnings.append("no concordant current TMbed prediction or curated TM feature")

    if (
        not decision_code
        and base["membrane_decision_v5"]
        == "new_htp_entry_without_experimental_or_membranome_support"
    ):
        candidate_class = "A"
        if signal_confusion:
            final_class = "unknown"
            evidence_status = "excluded"
            inclusion_status = "excluded"
            decision_code = "excluded_HTP_signal_peptide_confusion"
            evidence.append("UniTmp HTP TM call overlaps a UniProt signal peptide")
        elif hydro_after_signal and not secreted:
            final_class = "unknown"
            evidence_status = "uncertain"
            inclusion_status = "candidate"
            decision_code = "A_candidate_UniTmp_plus_hydropathy"
            evidence.extend(
                [
                    f"UniTmp HTP evidence class: {htp_evidence or 'unknown'}",
                    "hydropathy scan finds a non-signal hydrophobic segment",
                ]
            )
            warnings.append("no concordant current curated or TMbed TM annotation")
        else:
            final_class = "unknown"
            evidence_status = "uncertain"
            inclusion_status = "candidate"
            decision_code = "A_candidate_UniTmp_only"
            evidence.append(f"UniTmp HTP evidence class: {htp_evidence or 'unknown'}")
            warnings.append("no independent sequence/topology confirmation")

    if (
        not decision_code
        and base["membrane_decision_v5"]
        == "new_opm_membrane_association_without_integral_segments"
    ):
        candidate_class = "C"
        final_class = "unknown"
        if secreted or lumen_only:
            evidence_status = "excluded"
            inclusion_status = "excluded"
            decision_code = "excluded_OPM_secreted_or_lumen_protein"
            evidence.append("OPM includes a non-integral structure")
            warnings.append(
                "UniProt describes a secreted/lumen protein rather than a stable cellular membrane protein"
            )
        else:
            evidence_status = "uncertain"
            inclusion_status = "candidate"
            decision_code = "C_candidate_OPM_nonintegral_structure"
            evidence.append("OPM includes a non-integral target structure")

    if not decision_code:
        candidate_class = "C"
        final_class = "unknown"
        evidence_status = "uncertain"
        inclusion_status = "candidate"
        decision_code = "unresolved_external_membrane_candidate"
        warnings.append("no rule reached a defensible membrane-class conclusion")

    if target_chain_non_tm:
        warnings.append("PDBTM target chain is explicitly non-transmembrane")
    if chain_rows and not target_chain_tm and base.get("pdbtm_present_v5") == "1":
        warnings.append(
            "PDBTM membrane status belongs to another chain or the target chain is non-TM"
        )
    if tmbed_available and not tmbed_sequence_exact:
        warnings.append("TMbed 2022 prediction sequence differs from current UniProt sequence")
    if not tmbed_available:
        warnings.append("no TMbed 2022 prediction matched the accession")
    if matching_domains:
        evidence.append("membrane-binding-capable domain(s) detected")

    relevant_text = " ".join(
        [subcellular, function, domain_comment, ptm_comment]
    )
    literature = pubmed_ids(relevant_text)
    go_summary = ";".join(
        f"{item['go_id']}:{item['go_name']}|{item['evidence_code']}|{item['reference']}"
        for item in membrane_go
    )

    priority = 5
    if int(base.get("independent_membrane_source_count_v5", "0") or 0) >= 2:
        priority = 1
    elif base["membrane_decision_v5"].startswith("new_htp"):
        priority = 2
    elif base["membrane_decision_v5"].startswith("new_opm"):
        priority = 3
    elif base["membrane_decision_v5"].startswith("new_hpa_predicted"):
        priority = 4

    return {
        "target_uniprot_id": accession,
        "approved_symbol": base.get("approved_symbol", ""),
        "protein_name": base.get("protein_name", ""),
        "v5_review_reason": base.get("membrane_decision_v5", ""),
        "v5_independent_source_count": base.get(
            "independent_membrane_source_count_v5", ""
        ),
        "v5_independent_sources": base.get("independent_membrane_sources_v5", ""),
        "review_priority": priority,
        "candidate_membrane_class_v51": candidate_class,
        "final_membrane_class_v51": final_class,
        "evidence_status_v51": evidence_status,
        "inclusion_status_v51": inclusion_status,
        "decision_code_v51": decision_code,
        "decision_basis_v51": "; ".join(dict.fromkeys(evidence)),
        "cautions_v51": "; ".join(dict.fromkeys(warnings)),
        "hpa_reliability": hpa_reliability,
        "hpa_main_location": hpa_main,
        "hpa_additional_location": hpa_additional,
        "htp_evidence": htp_evidence,
        "htp_num_tm": htp_tm_count,
        "htp_reliability": htp.get("htp_reliability", ""),
        "htp_first_tm_start": htp_first_tm_start or "",
        "uniprot_signal_feature": signal_feature,
        "uniprot_lipidation": lipidation,
        "uniprot_subcellular_location": subcellular,
        "uniprot_go_cellular_component": uniprot.get(
            "Gene Ontology (cellular component)", ""
        ),
        "goa_membrane_annotations": go_summary,
        "goa_experimental_membrane_count": len(experimental_membrane_go),
        "tmbed_prediction_available": "1" if tmbed_available else "0",
        "tmbed_sequence_exact": "1" if tmbed_sequence_exact else "0",
        "tmbed_tm_segment_count": len(tm_segments) + len(beta_segments),
        "tmbed_tm_segments": segment_text(tm_segments + beta_segments),
        "tmbed_signal_segments": segment_text(sp_segments),
        "hydropathy_segment_count": len(hydro_segments),
        "hydropathy_segments": segment_text(hydro_segments),
        "hydropathy_non_signal_segment_count": len(hydro_after_signal),
        "hydropathy_non_signal_segments": segment_text(hydro_after_signal),
        "membrane_binding_domains": ";".join(matching_domains),
        "direct_membrane_binding_text_flag": "1" if direct_membrane_text else "0",
        "opm_membrane_contexts": ";".join(opm_membranes),
        "opm_pdb_ids": base.get("opm_pdb_ids_v5", ""),
        "pdbtm_pdb_ids": base.get("pdbtm_pdb_ids_v5", ""),
        "pdbtm_target_chain_tm_flag": "1" if target_chain_tm else "0",
        "pdbtm_target_chain_non_tm_flag": "1" if target_chain_non_tm else "0",
        "supporting_pubmed_ids": ";".join(literature),
        "uniprot_url": f"https://www.uniprot.org/uniprotkb/{accession}/entry",
        "hpa_url": (
            f"https://www.proteinatlas.org/{hpa.get('hpa_ensembl', '')}"
            if hpa.get("hpa_ensembl", "")
            else ""
        ),
        "review_date_v51": date.today().isoformat(),
        "review_method_v51": "versioned_rule_based_evidence_integration",
    }


def main() -> None:
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)

    queue_rows = load_tsv(QUEUE)
    if len(queue_rows) != 1109:
        raise ValueError("Frozen R queue row count is not 1,109")
    accessions = {row["target_uniprot_id"].upper() for row in queue_rows}

    uniprot_index = {
        row["Entry"].upper(): row for row in load_tsv(UNIPROT)
    }
    hpa_index = {
        row["target_uniprot_id"].upper(): row for row in load_tsv(HPA)
    }
    htp_index = {
        row["uniprot_entry_name"].upper(): row for row in load_tsv(HTP)
    }
    opm_index: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in load_tsv(OPM):
        opm_index[row["target_uniprot_id"].upper()].append(row)
    chain_index: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in load_tsv(PDBTM_CHAINS):
        chain_index[row["target_uniprot_id"].upper()].append(row)

    terms = parse_obo(GO_OBO)
    membrane_go_terms = descendants("GO:0016020", terms)
    go_annotations = parse_goa(GOA, accessions, terms)
    interpro_names = {
        row["ENTRY_AC"]: row["ENTRY_NAME"] for row in load_tsv(INTERPRO)
    }
    tmbed_index = parse_tmbed(TM_BED, accessions)

    output_rows = []
    for base in queue_rows:
        accession = base["target_uniprot_id"].upper()
        if accession not in uniprot_index:
            raise KeyError(f"Missing current UniProt record: {accession}")
        entry_name = uniprot_index[accession].get("Entry Name", "").upper()
        output_rows.append(
            classify_record(
                base=base,
                uniprot=uniprot_index[accession],
                hpa=hpa_index.get(accession, {}),
                htp=htp_index.get(entry_name, {}),
                opm_rows=opm_index.get(accession, []),
                chain_rows=chain_index.get(accession, []),
                go_annotations=go_annotations.get(accession, []),
                membrane_go_terms=membrane_go_terms,
                interpro_names=interpro_names,
                tmbed=tmbed_index.get(accession),
            )
        )

    columns = list(output_rows[0])
    output_path = ANALYSIS / "r1109_preliminary_review_v51.tsv"
    write_tsv(output_path, output_rows, columns)

    stats = {
        "review_date": date.today().isoformat(),
        "baseline_rows": len(queue_rows),
        "unique_accessions": len(accessions),
        "candidate_class_counts": Counter(
            row["candidate_membrane_class_v51"] for row in output_rows
        ),
        "final_class_counts": Counter(
            row["final_membrane_class_v51"] for row in output_rows
        ),
        "evidence_status_counts": Counter(
            row["evidence_status_v51"] for row in output_rows
        ),
        "inclusion_status_counts": Counter(
            row["inclusion_status_v51"] for row in output_rows
        ),
        "decision_code_counts": Counter(
            row["decision_code_v51"] for row in output_rows
        ),
        "priority_counts": Counter(str(row["review_priority"]) for row in output_rows),
        "tmbed_accession_coverage": len(tmbed_index),
        "tmbed_exact_sequence_count": sum(
            row["tmbed_sequence_exact"] == "1" for row in output_rows
        ),
        "pdbtm_target_chain_tm_hits": sum(
            row["pdbtm_target_chain_tm_flag"] == "1" for row in output_rows
        ),
        "pdbtm_target_chain_non_tm_hits": sum(
            row["pdbtm_target_chain_non_tm_flag"] == "1" for row in output_rows
        ),
        "input_hashes": {
            str(path.relative_to(ROOT)): sha256(path)
            for path in [
                QUEUE,
                UNIPROT,
                TM_BED,
                GO_OBO,
                GOA,
                INTERPRO,
                HPA,
                HTP,
                OPM,
                PDBTM_CHAINS,
            ]
        },
        "output_file": str(output_path),
        "output_sha256": sha256(output_path),
    }
    stats_path = REPORTS / "V51_PRELIMINARY_REVIEW_STATS.json"
    stats_path.write_text(
        json.dumps(stats, indent=2, ensure_ascii=False, default=dict) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(stats, indent=2, ensure_ascii=False, default=dict))


if __name__ == "__main__":
    main()
