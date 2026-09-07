#!/usr/bin/env python3
"""Build the publication-oriented v5 human membrane-protein master layer.

The script uses only Python's standard library. It keeps source-specific
evidence separate, resolves identifiers deterministically, and never promotes
location-only evidence to integral-membrane status.
"""

from __future__ import annotations

import csv
import html
import json
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


WORK = Path(__file__).resolve().parents[1]
REPO = WORK.parent
RAW = WORK / "raw"
NORMALIZED = WORK / "normalized"
CROSSWALKS = WORK / "crosswalks"
REVIEW = WORK / "review"
REPORTS = WORK / "reports"
RELEASE = WORK / "release"
V4_RELEASE = REPO / "v4_working" / "releases" / "release_v4"
V4_RAW = REPO / "v4_working" / "raw"
TODAY = date.today().isoformat()

for directory in (NORMALIZED, CROSSWALKS, REVIEW, REPORTS, RELEASE):
    directory.mkdir(parents=True, exist_ok=True)


V4_FIELDS = [
    "membrane_protein_id",
    "target_uniprot_id",
    "uniprot_entry_name",
    "reviewed",
    "protein_name",
    "approved_symbol",
    "gene_names",
    "hgnc_ids",
    "ncbi_gene_ids",
    "ensembl_gene_ids",
    "organism_id",
    "sequence_length",
    "membrane_scope",
    "membrane_topology",
    "transmembrane_count",
    "intramembrane_count",
    "membrane_inclusion_status",
    "membrane_evidence_basis",
    "transmembrane_features",
    "intramembrane_features",
    "signal_peptide_features",
    "subcellular_location",
    "go_cellular_component",
    "uniprot_keywords",
    "protein_families",
    "pdb_ids",
    "alphafolddb_ids",
    "functional_primary_class",
    "functional_subclass",
    "functional_class_rule",
    "present_in_v3_1",
    "v3_1_unique_drug_count",
    "source",
    "source_query",
    "retrieval_date",
]

V5_FIELDS = [
    "release_tier_v5",
    "release_tier_label_v5",
    "record_status_v5",
    "membrane_decision_v5",
    "membrane_confidence_v5",
    "membrane_scope_v5",
    "membrane_topology_v5",
    "transmembrane_count_v5",
    "independent_membrane_source_count_v5",
    "independent_membrane_sources_v5",
    "hpa_predicted_membrane_v5",
    "hpa_plasma_membrane_location_v5",
    "hpa_reliability_if_v5",
    "hpa_subcellular_location_v5",
    "htp_present_v5",
    "htp_evidence_v5",
    "htp_num_tm_v5",
    "htp_reliability_v5",
    "htp_n_terminal_side_v5",
    "htp_c_terminal_side_v5",
    "membranome_present_v5",
    "membranome_family_v5",
    "membranome_segment_v5",
    "opm_present_v5",
    "opm_integral_structure_v5",
    "opm_pdb_ids_v5",
    "pdbtm_present_v5",
    "pdbtm_pdb_ids_v5",
    "single_pass_resolution_v5",
    "single_pass_resolution_basis_v5",
    "functional_primary_class_v5",
    "functional_subclass_v5",
    "classification_source_v5",
    "classification_status_v5",
    "gpcrdb_class_v5",
    "gpcrdb_family_v5",
    "gpcrdb_ligand_type_v5",
    "gpcrdb_subfamily_v5",
    "gtopdb_type_v5",
    "gtopdb_family_id_v5",
    "gtopdb_family_name_v5",
    "tcdb_tcids_v5",
    "cross_source_conflict_v5",
    "review_reason_v5",
]


def clean(value) -> str:
    if value is None:
        return ""
    return str(value).strip()


def strip_html(value: str) -> str:
    value = html.unescape(clean(value))
    return re.sub(r"<[^>]+>", "", value).strip()


def split_tokens(value: str) -> list[str]:
    return [x.strip() for x in re.split(r"[;,|]", clean(value)) if x.strip()]


def ensembl_genes(value: str) -> set[str]:
    return set(re.findall(r"ENSG\d+", clean(value).upper()))


def numeric_ids(value: str) -> set[str]:
    return set(re.findall(r"\d+", clean(value)))


def read_tsv(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def normalize_pdb(value: str) -> str:
    value = clean(value).lower()
    value = value.replace('="', "").replace('"', "")
    match = re.search(r"\b[0-9][a-z0-9]{3}\b", value)
    return match.group(0) if match else ""


def parse_hpa(reviewed_accessions: set[str]) -> tuple[dict[str, dict], list[dict]]:
    index: dict[str, dict] = {}
    normalized: list[dict] = []
    archive = RAW / "hpa_proteinatlas_v25_1.tsv.zip"
    with zipfile.ZipFile(archive) as zf:
        member = zf.namelist()[0]
        with zf.open(member) as binary:
            text = (line.decode("utf-8-sig") for line in binary)
            reader = csv.DictReader(text, delimiter="\t")
            for row in reader:
                protein_class = clean(row.get("Protein class"))
                subcell = clean(row.get("Subcellular location"))
                main_loc = clean(row.get("Subcellular main location"))
                add_loc = clean(row.get("Subcellular additional location"))
                predicted = "predicted membrane proteins" in protein_class.lower()
                plasma = "plasma membrane" in " ".join(
                    [subcell, main_loc, add_loc]
                ).lower()
                accessions = [
                    token
                    for token in split_tokens(row.get("Uniprot", ""))
                    if token in reviewed_accessions
                ]
                for accession in accessions:
                    evidence = {
                        "target_uniprot_id": accession,
                        "hpa_gene": clean(row.get("Gene")),
                        "hpa_ensembl": clean(row.get("Ensembl")),
                        "hpa_predicted_membrane": "1" if predicted else "0",
                        "hpa_plasma_membrane_location": "1" if plasma else "0",
                        "hpa_reliability_if": clean(row.get("Reliability (IF)")),
                        "hpa_subcellular_location": subcell,
                        "hpa_subcellular_main_location": main_loc,
                        "hpa_subcellular_additional_location": add_loc,
                        "hpa_protein_class": protein_class,
                        "source_version": "HPA 25.1",
                    }
                    index[accession] = evidence
                    if predicted or plasma:
                        normalized.append(evidence)
    return index, normalized


def parse_htp() -> tuple[dict[str, dict], list[dict]]:
    index: dict[str, dict] = {}
    rows: list[dict] = []
    path = RAW / "unitmp_htp_all_d2_2.xml"
    for _, element in ET.iterparse(path, events=("end",)):
        if element.tag.split("}")[-1] != "HtpItem":
            continue
        entry_name = clean(element.attrib.get("id")).upper()
        topology = element.find("./Topology")
        num_tm = clean(topology.attrib.get("numTm")) if topology is not None else ""
        reliability = (
            clean(topology.attrib.get("reliability")) if topology is not None else ""
        )
        regions = []
        if topology is not None:
            for region in topology.findall("./Region"):
                regions.append(
                    {
                        "from": clean(region.attrib.get("from")),
                        "to": clean(region.attrib.get("to")),
                        "loc": clean(region.attrib.get("loc")),
                    }
                )
        non_tm = [r for r in regions if r["loc"] in {"I", "O"}]
        tm_regions = [r for r in regions if r["loc"] == "M"]
        n_side = non_tm[0]["loc"] if non_tm else ""
        c_side = non_tm[-1]["loc"] if non_tm else ""
        tm_start = tm_regions[0]["from"] if tm_regions else ""
        evidence = {
            "uniprot_entry_name": entry_name,
            "htp_protein_name": clean(element.findtext("./Name")),
            "htp_transmembrane": clean(element.attrib.get("transmembrane")),
            "htp_evidence": clean(element.attrib.get("evidence")),
            "htp_num_tm": num_tm,
            "htp_reliability": reliability,
            "htp_n_terminal_side": n_side,
            "htp_c_terminal_side": c_side,
            "htp_first_tm_start": tm_start,
            "source_version": "UniTmp HTP d.2.2",
        }
        index[entry_name] = evidence
        rows.append(evidence)
        element.clear()
    return index, rows


def parse_membranome(
    accession_set: set[str], entry_to_accession: dict[str, str]
) -> tuple[dict[str, list[dict]], list[dict]]:
    index: dict[str, list[dict]] = defaultdict(list)
    rows: list[dict] = []
    path = RAW / "membranome_proteins_2026-07-24.csv"
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            if "homo sapiens" not in clean(row.get("species_name_cache")).lower():
                continue
            accession = clean(row.get("uniprotcode")).upper()
            entry_name = clean(row.get("uniprot_id")).upper()
            if accession not in accession_set:
                accession = entry_to_accession.get(entry_name, "")
            out = {
                "target_uniprot_id": accession,
                "uniprot_entry_name": entry_name,
                "membranome_name": clean(row.get("name")),
                "membranome_gene": clean(row.get("genename")),
                "membranome_family": clean(row.get("protein_family_name_cache")),
                "membranome_segment": clean(row.get("segment")),
                "membranome_membrane": clean(row.get("membrane_name_cache")),
                "membranome_opm_links_count": clean(row.get("opm_links_count")),
                "membranome_pdb_links_count": clean(row.get("pdb_links_count")),
                "source_version": "BioMembHub Membranome current 2026-07-24",
            }
            rows.append(out)
            if accession:
                index[accession].append(out)
    return index, rows


def parse_opm(
    accession_set: set[str], entry_to_accession: dict[str, str]
) -> tuple[dict[str, list[dict]], list[dict]]:
    index: dict[str, list[dict]] = defaultdict(list)
    rows: list[dict] = []
    path = RAW / "opm_primary_structures_2026-07-24.csv"
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            if "homo sapiens" not in clean(row.get("species_name_cache")).lower():
                continue
            code = clean(row.get("uniprotcode")).upper()
            accession = code if code in accession_set else entry_to_accession.get(code, "")
            pdb_id = normalize_pdb(row.get("pdbid", ""))
            topology_subunit = clean(row.get("topology_subunit"))
            segments = clean(row.get("subunit_segments"))
            integral = bool(topology_subunit and (segments or clean(row.get("thickness"))))
            out = {
                "target_uniprot_id": accession,
                "uniprot_entry_name_or_accession": code,
                "opm_pdb_id": pdb_id,
                "opm_name": clean(row.get("name")),
                "opm_family": clean(row.get("family_name_cache")),
                "opm_membrane": clean(row.get("membrane_name_cache")),
                "opm_topology_subunit": topology_subunit,
                "opm_subunit_segments": segments,
                "opm_thickness": clean(row.get("thickness")),
                "opm_integral_structure": "1" if integral else "0",
                "source_version": "OPM current 2026-07-24",
            }
            rows.append(out)
            if accession:
                index[accession].append(out)
    return index, rows


def parse_pdbtm() -> tuple[set[str], list[dict]]:
    pdb_ids: set[str] = set()
    rows: list[dict] = []
    path = RAW / "unitmp_pdbtm_all_2026-07-24.xml"
    for _, element in ET.iterparse(path, events=("end",)):
        if element.tag.split("}")[-1] != "pdbtm":
            continue
        pdb_id = clean(element.attrib.get("ID")).lower()
        chains = element.findall("./CHAIN")
        tm_counts = []
        chain_types = []
        for chain in chains:
            tm_counts.append(clean(chain.attrib.get("NUM_TM")))
            chain_types.append(clean(chain.attrib.get("TYPE")))
        row = {
            "pdb_id": pdb_id,
            "pdbtm_tmp": clean(element.attrib.get("TMP")),
            "pdbtm_chain_count": str(len(chains)),
            "pdbtm_chain_tm_counts": ";".join(tm_counts),
            "pdbtm_chain_types": ";".join(chain_types),
            "source_version": "UniTmp PDBTM current 2026-07-24",
        }
        rows.append(row)
        if pdb_id:
            pdb_ids.add(pdb_id)
        element.clear()
    return pdb_ids, rows


def parse_tcdb(accession_set: set[str]) -> tuple[dict[str, list[dict]], list[dict]]:
    index: dict[str, list[dict]] = defaultdict(list)
    rows: list[dict] = []
    path = RAW / "tcdb_human_2026-07-24.csv"
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            accession = clean(row.get("Accession")).upper()
            out = {
                "target_uniprot_id": accession,
                "tcdb_tcid": clean(row.get("TCID")),
                "tcdb_name": clean(row.get("Name")),
                "tcdb_symbol": clean(row.get("Symbol")),
                "tcdb_aliases": clean(row.get("Aliases")),
                "source_version": "TCDB human.csv retrieved 2026-07-24",
            }
            rows.append(out)
            if accession in accession_set:
                index[accession].append(out)
    return index, rows


def parse_gpcrdb(accession_set: set[str]) -> tuple[dict[str, dict], list[dict]]:
    payload = json.loads(
        (RAW / "gpcrdb_receptorlist_2026-07-24.json").read_text(encoding="utf-8")
    )
    index: dict[str, dict] = {}
    rows: list[dict] = []
    for item in payload:
        if clean(item.get("species")) != "Homo sapiens":
            continue
        accession = clean(item.get("accession")).upper()
        if accession not in accession_set:
            continue
        row = {
            "target_uniprot_id": accession,
            "gpcrdb_entry_name": clean(item.get("entry_name")).upper(),
            "gpcrdb_name": strip_html(item.get("name", "")),
            "gpcrdb_class": strip_html(item.get("receptor_class", "")),
            "gpcrdb_family": strip_html(item.get("receptor_family", "")),
            "gpcrdb_ligand_type": strip_html(item.get("ligand_type", "")),
            "gpcrdb_subfamily": strip_html(item.get("subfamily", "")),
            "source_version": "GPCRdb current 2026-07-24",
        }
        index[accession] = row
        rows.append(row)
    return index, rows


def parse_gtopdb(accession_set: set[str]) -> tuple[dict[str, list[dict]], list[dict]]:
    index: dict[str, list[dict]] = defaultdict(list)
    rows: list[dict] = []
    path = RAW / "gtopdb_targets_and_families_2026_2.tsv"
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        if not first.startswith('"#'):
            handle.seek(0)
        for row in csv.DictReader(handle, delimiter="\t"):
            accessions = split_tokens(row.get("Human SwissProt", ""))
            for accession in accessions:
                accession = accession.upper()
                if accession not in accession_set:
                    continue
                out = {
                    "target_uniprot_id": accession,
                    "gtopdb_type": clean(row.get("Type")),
                    "gtopdb_family_id": clean(row.get("Family id")),
                    "gtopdb_family_name": strip_html(row.get("Family name", "")),
                    "gtopdb_target_id": clean(row.get("Target id")),
                    "gtopdb_target_name": strip_html(row.get("Target name", "")),
                    "gtopdb_hgnc_symbol": clean(row.get("HGNC symbol")),
                    "source_version": "GtoPdb 2026.2",
                }
                index[accession].append(out)
                rows.append(out)
    return index, rows


def parse_feature_count(value: str, feature_name: str) -> int:
    return len(re.findall(rf"\b{re.escape(feature_name)}\b", clean(value)))


def reviewed_record(raw: dict) -> dict:
    transmembrane = clean(raw.get("Transmembrane"))
    intramembrane = clean(raw.get("Intramembrane"))
    return {
        "membrane_protein_id": "HMP_" + clean(raw.get("Entry")),
        "target_uniprot_id": clean(raw.get("Entry")),
        "uniprot_entry_name": clean(raw.get("Entry Name")),
        "reviewed": "reviewed",
        "protein_name": clean(raw.get("Protein names")),
        "approved_symbol": clean(raw.get("Gene Names (primary)")),
        "gene_names": clean(raw.get("Gene Names")),
        "hgnc_ids": clean(raw.get("HGNC")),
        "ncbi_gene_ids": clean(raw.get("GeneID")),
        "ensembl_gene_ids": clean(raw.get("Ensembl")),
        "organism_id": clean(raw.get("Organism (ID)")),
        "sequence_length": clean(raw.get("Length")),
        "membrane_scope": "uncertain_membrane_scope",
        "membrane_topology": "unresolved",
        "transmembrane_count": str(parse_feature_count(transmembrane, "TRANSMEM")),
        "intramembrane_count": str(parse_feature_count(intramembrane, "INTRAMEM")),
        "membrane_inclusion_status": "review_external_database_entry",
        "membrane_evidence_basis": "",
        "transmembrane_features": transmembrane,
        "intramembrane_features": intramembrane,
        "signal_peptide_features": clean(raw.get("Signal peptide")),
        "subcellular_location": clean(raw.get("Subcellular location [CC]")),
        "go_cellular_component": clean(raw.get("Gene Ontology (cellular component)")),
        "uniprot_keywords": clean(raw.get("Keywords")),
        "protein_families": clean(raw.get("Protein families")),
        "pdb_ids": clean(raw.get("PDB")),
        "alphafolddb_ids": clean(raw.get("AlphaFoldDB")),
        "functional_primary_class": "unclassified",
        "functional_subclass": "",
        "functional_class_rule": "",
        "present_in_v3_1": "0",
        "v3_1_unique_drug_count": "",
        "source": "UniProt reviewed reference proteome + external membrane source",
        "source_query": "UP000005640; organism_id:9606; reviewed:true",
        "retrieval_date": TODAY,
    }


def source_flags(
    accession: str,
    record: dict,
    hpa_index: dict[str, dict],
    htp_index: dict[str, dict],
    membranome_index: dict[str, list[dict]],
    opm_index: dict[str, list[dict]],
    pdbtm_ids: set[str],
) -> dict:
    hpa = hpa_index.get(accession, {})
    htp = htp_index.get(clean(record.get("uniprot_entry_name")).upper(), {})
    membranome = membranome_index.get(accession, [])
    opm = opm_index.get(accession, [])
    record_pdbs = {normalize_pdb(x) for x in split_tokens(record.get("pdb_ids", ""))}
    record_pdbs.discard("")
    pdbtm_hits = sorted(record_pdbs & pdbtm_ids)
    opm_pdbs = sorted({x["opm_pdb_id"] for x in opm if x["opm_pdb_id"]})
    opm_integral = any(x["opm_integral_structure"] == "1" for x in opm)

    sources = []
    if hpa.get("hpa_predicted_membrane") == "1":
        sources.append("HPA_predicted_membrane")
    if hpa.get("hpa_plasma_membrane_location") == "1":
        sources.append("HPA_plasma_membrane_location")
    if htp:
        sources.append("UniTmp_HTP")
    if membranome:
        sources.append("Membranome")
    if opm:
        sources.append("OPM")
    if pdbtm_hits:
        sources.append("PDBTM")

    return {
        "hpa": hpa,
        "htp": htp,
        "membranome": membranome,
        "opm": opm,
        "pdbtm_hits": pdbtm_hits,
        "opm_pdbs": opm_pdbs,
        "opm_integral": opm_integral,
        "sources": sources,
    }


def decide_tier(record: dict, flags: dict, existed_v4: bool) -> tuple[str, str, str, str]:
    scope = clean(record.get("membrane_scope"))
    htp_evidence = clean(flags["htp"].get("htp_evidence"))
    membranome = bool(flags["membranome"])
    opm_integral = bool(flags["opm_integral"])
    strong_integral = membranome or opm_integral or htp_evidence in {"3D", "Experiment"}

    if existed_v4:
        if scope == "integral_membrane":
            return "A", "core_integral_membrane", "retained_v4_integral", "high"
        if scope in {"integral_monotopic", "lipid_anchored"}:
            if strong_integral:
                return "A", "core_integral_membrane", "external_integral_evidence_promoted", "high"
            return "B", "core_extended_membrane", "retained_v4_monotopic_or_lipid_anchored", "high"
        if scope == "peripheral_membrane_associated":
            if strong_integral:
                return "A", "core_integral_membrane", "external_integral_evidence_promoted", "high"
            return "C", "peripheral_membrane_associated", "retained_v4_peripheral", "moderate"

    if strong_integral:
        return "A", "core_integral_membrane", "new_external_integral_evidence", "high"
    if flags["htp"]:
        return "R", "review", "new_htp_entry_without_experimental_or_membranome_support", "review"
    if flags["hpa"].get("hpa_predicted_membrane") == "1":
        return "R", "review", "new_hpa_predicted_membrane_only", "review"
    if flags["hpa"].get("hpa_plasma_membrane_location") == "1":
        return "R", "review", "new_hpa_location_only", "review"
    if flags["opm"]:
        return "R", "review", "new_opm_membrane_association_without_integral_segments", "review"
    return "R", "review", "external_database_entry_requires_review", "review"


def resolve_single_pass(record: dict, htp: dict) -> tuple[str, str]:
    if clean(record.get("membrane_topology")) not in {
        "single_pass_unresolved",
        "single_pass_unresolved_type",
    }:
        return clean(record.get("membrane_topology")), "not_applicable_or_already_curated"
    if clean(htp.get("htp_num_tm")) != "1":
        return "single_pass_unresolved", "no_matching_single_pass_htp_topology"
    n_side = clean(htp.get("htp_n_terminal_side"))
    c_side = clean(htp.get("htp_c_terminal_side"))
    try:
        tm_start = int(clean(htp.get("htp_first_tm_start")) or "0")
    except ValueError:
        tm_start = 0
    has_signal = bool(clean(record.get("signal_peptide_features")))
    if n_side == "I" and c_side == "O":
        return "single_pass_type_ii", "UniTmp HTP N-in/C-out orientation"
    if n_side == "O" and c_side == "I":
        if has_signal:
            return "single_pass_type_i", "UniTmp HTP N-out/C-in plus UniProt signal peptide"
        if tm_start and tm_start <= 60:
            return "single_pass_type_iii", "UniTmp HTP N-out/C-in and N-terminal signal anchor"
        return (
            "single_pass_type_i_or_iii",
            "UniTmp HTP N-out/C-in; type I versus III not distinguishable",
        )
    return "single_pass_unresolved", "HTP terminal orientation unavailable_or_noncanonical"


def classify(
    accession: str,
    record: dict,
    gpcrdb_index: dict[str, dict],
    gtopdb_index: dict[str, list[dict]],
    tcdb_index: dict[str, list[dict]],
) -> dict:
    gpcr = gpcrdb_index.get(accession, {})
    gtop = gtopdb_index.get(accession, [])
    tcdb = tcdb_index.get(accession, [])
    old_primary = clean(record.get("functional_primary_class")) or "unclassified"
    old_subclass = clean(record.get("functional_subclass"))
    unresolved_old = old_primary == "unclassified" or old_subclass == "unclassified"
    primary = "unclassified" if unresolved_old else old_primary
    subclass = "" if unresolved_old else old_subclass
    source = clean(record.get("functional_class_rule")) or "v4"
    status = "retained_v4"

    if gpcr:
        primary = "receptor"
        subclass = "GPCR | " + " | ".join(
            x
            for x in [
                gpcr.get("gpcrdb_class"),
                gpcr.get("gpcrdb_ligand_type"),
                gpcr.get("gpcrdb_family"),
                gpcr.get("gpcrdb_subfamily"),
            ]
            if x
        )
        source = "GPCRdb"
        status = "specialist_database_aligned"
    elif gtop:
        target_type = clean(gtop[0].get("gtopdb_type")).lower()
        type_map = {
            "gpcr": "receptor",
            "catalytic_receptor": "receptor",
            "ion_channel": "ion_channel",
            "vgic": "ion_channel",
            "lgic": "ion_channel",
            "other_ic": "ion_channel",
            "transporter": "transporter",
            "enzyme": "enzyme",
            "nuclear_receptor": "receptor",
            "nhr": "receptor",
            "other_protein": "other_protein_target",
        }
        mapped = type_map.get(target_type, "")
        if mapped and (
            unresolved_old or mapped in {"receptor", "ion_channel", "transporter"}
        ):
            primary = mapped
            subclass = "GtoPdb | " + " | ".join(
                x
                for x in [
                    gtop[0].get("gtopdb_type"),
                    gtop[0].get("gtopdb_family_name"),
                    gtop[0].get("gtopdb_target_name"),
                ]
                if x
            )
            source = "IUPHAR/BPS Guide to PHARMACOLOGY 2026.2"
            status = "specialist_database_aligned"
    elif tcdb:
        primary = "transporter"
        subclass = "TCDB | " + ";".join(
            sorted({clean(x.get("tcdb_tcid")) for x in tcdb if clean(x.get("tcdb_tcid"))})
        )
        source = "TCDB"
        status = "specialist_database_aligned"
    elif unresolved_old:
        family = clean(record.get("protein_families"))
        keywords = clean(record.get("uniprot_keywords"))
        text = " ".join(
            [
                clean(record.get("protein_name")),
                family,
                keywords,
                clean(record.get("approved_symbol")),
            ]
        ).lower()
        if re.search(r"\b(ion channel|channel subunit|aquaporin|connexin|pannexin|porin)\b", text):
            primary = "ion_channel"
            subclass = "name_or_family_inference"
        elif re.search(
            r"\b(transporter|solute carrier|permease|antiporter|symporter|exchanger|abc transporter)\b",
            text,
        ):
            primary = "transporter"
            subclass = "name_or_family_inference"
        elif "receptor" in text:
            primary = "receptor"
            subclass = "name_or_family_inference"
        elif re.search(r"\b(cadherin|integrin|cell adhesion|adhesion molecule)\b", text):
            primary = "cell_adhesion"
            subclass = "name_or_family_inference"
        elif re.search(r"\b(snare|vesicle-associated|vesicular trafficking)\b", text):
            primary = "membrane_trafficking"
            subclass = "name_or_family_inference"
        elif re.search(
            r"\b(sorting nexin|synaptotagmin|syntaxin|adaptor complex|rab family|arf family)\b",
            text,
        ):
            primary = "membrane_trafficking"
            subclass = "name_or_family_inference"
        elif re.search(
            r"\b(myosin|kinesin|intermediate filament|annexin|cytoskeleton)\b",
            text,
        ):
            primary = "cytoskeleton_membrane_linker"
            subclass = "name_or_family_inference"
        elif re.search(
            r"\b(mhc class|major histocompatibility|immunoglobulin|tumor necrosis factor|immune)\b",
            text,
        ):
            primary = "immune_or_cell_recognition"
            subclass = "name_or_family_inference"
        elif re.search(
            r"\b(g protein gamma|small gtpase|signal transduction|apoptosis)\b",
            text,
        ):
            primary = "signaling_regulator"
            subclass = "name_or_family_inference"
        elif family:
            primary = "other_family_defined"
            subclass = "UniProt family | " + family
            source = "UniProt family annotation"
            status = "family_level_classified_function_unresolved"
        if primary != "unclassified":
            if status != "family_level_classified_function_unresolved":
                source = "conservative UniProt name/family rule"
                status = "rule_classified"
        else:
            source = "unresolved"
            status = "unclassified"

    return {
        "primary": primary,
        "subclass": subclass,
        "source": source,
        "status": status,
        "gpcr": gpcr,
        "gtop": gtop[0] if gtop else {},
        "tcdb_tcids": ";".join(
            sorted({clean(x.get("tcdb_tcid")) for x in tcdb if clean(x.get("tcdb_tcid"))})
        ),
    }


def build_indices(reviewed_raw: list[dict]) -> dict[str, dict[str, set[str]]]:
    indices = {
        "symbol": defaultdict(set),
        "hgnc": defaultdict(set),
        "geneid": defaultdict(set),
        "ensembl": defaultdict(set),
    }
    for row in reviewed_raw:
        accession = clean(row.get("Entry"))
        symbol = clean(row.get("Gene Names (primary)")).upper()
        if symbol:
            indices["symbol"][symbol].add(accession)
        for value in numeric_ids(row.get("HGNC", "")):
            indices["hgnc"][value].add(accession)
        for value in numeric_ids(row.get("GeneID", "")):
            indices["geneid"][value].add(accession)
        for value in ensembl_genes(row.get("Ensembl", "")):
            indices["ensembl"][value].add(accession)
    return indices


def canonical_matches(candidate: dict, indices) -> tuple[list[str], list[str]]:
    method_sets: list[tuple[str, set[str]]] = []
    symbol = clean(candidate.get("Gene Names (primary)")).upper()
    if symbol and symbol in indices["symbol"]:
        method_sets.append(("exact_primary_gene_symbol", set(indices["symbol"][symbol])))
    for value in numeric_ids(candidate.get("HGNC", "")):
        if value in indices["hgnc"]:
            method_sets.append(("HGNC", set(indices["hgnc"][value])))
    for value in numeric_ids(candidate.get("GeneID", "")):
        if value in indices["geneid"]:
            method_sets.append(("NCBI_GeneID", set(indices["geneid"][value])))
    for value in ensembl_genes(candidate.get("Ensembl", "")):
        if value in indices["ensembl"]:
            method_sets.append(("Ensembl_gene", set(indices["ensembl"][value])))
    if not method_sets:
        return [], []
    scores: Counter[str] = Counter()
    methods_by_accession: dict[str, set[str]] = defaultdict(set)
    for method, accessions in method_sets:
        for accession in accessions:
            scores[accession] += 1
            methods_by_accession[accession].add(method)
    best_score = max(scores.values())
    best = sorted([accession for accession, score in scores.items() if score == best_score])
    methods = sorted(set().union(*(methods_by_accession[x] for x in best)))
    return best, methods


def main() -> None:
    v4_rows = read_tsv(V4_RELEASE / "human_membrane_protein_master_v4.tsv")
    reviewed_raw = read_tsv(
        V4_RAW / "uniprot_human_reviewed_reference_proteome_identity_v2_2026-07-24.tsv"
    )
    unreviewed_raw = read_tsv(
        V4_RAW / "uniprot_human_unreviewed_transmembrane_candidates_identity_v2_2026-07-24.tsv"
    )
    v4_unreviewed_review = read_tsv(
        V4_RELEASE / "unreviewed_transmembrane_candidate_review_v4.tsv"
    )
    v4_unreviewed_status = {
        clean(row.get("target_uniprot_id")): clean(
            row.get("canonical_resolution_status")
        )
        for row in v4_unreviewed_review
    }
    legacy_rows = read_tsv(V4_RELEASE / "v3_1_proteins_not_in_reference_v4.tsv")

    reviewed_by_accession = {clean(r["Entry"]): r for r in reviewed_raw}
    accession_set = set(reviewed_by_accession)
    entry_to_accession = {
        clean(r["Entry Name"]).upper(): clean(r["Entry"]) for r in reviewed_raw
    }
    v4_by_accession = {clean(r["target_uniprot_id"]): dict(r) for r in v4_rows}

    print("Parsing HPA...")
    hpa_index, hpa_rows = parse_hpa(accession_set)
    print("Parsing UniTmp HTP...")
    htp_index, htp_rows = parse_htp()
    print("Parsing Membranome...")
    membranome_index, membranome_rows = parse_membranome(
        accession_set, entry_to_accession
    )
    print("Parsing OPM...")
    opm_index, opm_rows = parse_opm(accession_set, entry_to_accession)
    print("Parsing PDBTM...")
    pdbtm_ids, pdbtm_rows = parse_pdbtm()
    print("Parsing TCDB, GPCRdb and GtoPdb...")
    tcdb_index, tcdb_rows = parse_tcdb(accession_set)
    gpcrdb_index, gpcrdb_rows = parse_gpcrdb(accession_set)
    gtopdb_index, gtopdb_rows = parse_gtopdb(accession_set)

    write_tsv(NORMALIZED / "hpa_membrane_evidence_v5.tsv", hpa_rows)
    write_tsv(NORMALIZED / "unitmp_htp_human_topology_v5.tsv", htp_rows)
    write_tsv(NORMALIZED / "membranome_human_proteins_v5.tsv", membranome_rows)
    write_tsv(NORMALIZED / "opm_human_structures_v5.tsv", opm_rows)
    write_tsv(NORMALIZED / "unitmp_pdbtm_structures_v5.tsv", pdbtm_rows)
    write_tsv(NORMALIZED / "tcdb_human_classification_v5.tsv", tcdb_rows)
    write_tsv(NORMALIZED / "gpcrdb_classification_v5.tsv", gpcrdb_rows)
    write_tsv(NORMALIZED / "gtopdb_classification_v5.tsv", gtopdb_rows)

    # Build reviewed union: v4 plus every reviewed entry found in an independent
    # membrane-specific source. HPA location-only evidence enters Tier R.
    selected_accessions = set(v4_by_accession)
    for accession, raw in reviewed_by_accession.items():
        probe = v4_by_accession.get(accession) or reviewed_record(raw)
        flags = source_flags(
            accession,
            probe,
            hpa_index,
            htp_index,
            membranome_index,
            opm_index,
            pdbtm_ids,
        )
        if (
            flags["htp"]
            or flags["membranome"]
            or flags["opm"]
            or flags["hpa"].get("hpa_predicted_membrane") == "1"
            or flags["hpa"].get("hpa_plasma_membrane_location") == "1"
        ):
            selected_accessions.add(accession)

    master_rows: list[dict] = []
    evidence_matrix: list[dict] = []
    conflicts: list[dict] = []
    topology_rows: list[dict] = []
    classification_rows: list[dict] = []

    for accession in sorted(selected_accessions):
        existed_v4 = accession in v4_by_accession
        record = (
            dict(v4_by_accession[accession])
            if existed_v4
            else reviewed_record(reviewed_by_accession[accession])
        )
        flags = source_flags(
            accession,
            record,
            hpa_index,
            htp_index,
            membranome_index,
            opm_index,
            pdbtm_ids,
        )
        tier, tier_label, decision, confidence = decide_tier(record, flags, existed_v4)
        htp = flags["htp"]
        membrane_sources = list(flags["sources"])
        if existed_v4:
            membrane_sources.insert(0, "UniProt_v4")

        topology, topology_basis = resolve_single_pass(record, htp)
        if clean(record.get("membrane_topology")) in {
            "single_pass_unresolved",
            "single_pass_unresolved_type",
        }:
            topology_rows.append(
                {
                    "target_uniprot_id": accession,
                    "approved_symbol": clean(record.get("approved_symbol")),
                    "original_topology": clean(record.get("membrane_topology")),
                    "resolved_topology_v5": topology,
                    "resolution_basis": topology_basis,
                    "htp_evidence": clean(htp.get("htp_evidence")),
                    "htp_num_tm": clean(htp.get("htp_num_tm")),
                    "htp_n_terminal_side": clean(htp.get("htp_n_terminal_side")),
                    "htp_c_terminal_side": clean(htp.get("htp_c_terminal_side")),
                    "uniprot_signal_peptide": clean(record.get("signal_peptide_features")),
                }
            )
            if topology in {
                "single_pass_type_i",
                "single_pass_type_ii",
                "single_pass_type_iii",
            }:
                record["membrane_topology"] = topology

        if not existed_v4:
            if tier == "A":
                record["membrane_scope"] = "integral_membrane"
                record["membrane_inclusion_status"] = "accepted_external_core"
                record["membrane_evidence_basis"] = "+".join(membrane_sources)
                try:
                    htp_tm = int(clean(htp.get("htp_num_tm")) or "0")
                except ValueError:
                    htp_tm = 0
                if htp_tm:
                    record["transmembrane_count"] = str(htp_tm)
                    record["membrane_topology"] = (
                        topology if htp_tm == 1 else "multi_pass"
                    )
            else:
                record["membrane_scope"] = "uncertain_membrane_scope"
                record["membrane_inclusion_status"] = "review_external_database_entry"
                record["membrane_evidence_basis"] = "+".join(membrane_sources)

        classification = classify(
            accession, record, gpcrdb_index, gtopdb_index, tcdb_index
        )
        gpcr = classification["gpcr"]
        gtop = classification["gtop"]

        conflict_reasons = []
        original_scope = clean(v4_by_accession.get(accession, {}).get("membrane_scope"))
        if existed_v4 and original_scope == "peripheral_membrane_associated" and tier == "A":
            conflict_reasons.append("v4_peripheral_vs_external_integral")
        try:
            v4_tm = int(clean(v4_by_accession.get(accession, {}).get("transmembrane_count")) or "0")
            htp_tm = int(clean(htp.get("htp_num_tm")) or "0")
            if existed_v4 and v4_tm and htp_tm and abs(v4_tm - htp_tm) > 1:
                conflict_reasons.append("uniprot_vs_htp_tm_count_difference_gt_1")
        except ValueError:
            pass
        if tier == "R":
            conflict_reasons.append(decision)

        record.update(
            {
                "release_tier_v5": tier,
                "release_tier_label_v5": tier_label,
                "record_status_v5": "retained_from_v4" if existed_v4 else "new_reviewed_external_union",
                "membrane_decision_v5": decision,
                "membrane_confidence_v5": confidence,
                "membrane_scope_v5": (
                    "integral_membrane"
                    if tier == "A"
                    else (
                        clean(record.get("membrane_scope"))
                        if tier in {"B", "C"}
                        else "uncertain_membrane_scope"
                    )
                ),
                "membrane_topology_v5": clean(record.get("membrane_topology")),
                "transmembrane_count_v5": clean(record.get("transmembrane_count")),
                "independent_membrane_source_count_v5": str(
                    len(set(flags["sources"]))
                ),
                "independent_membrane_sources_v5": ";".join(
                    sorted(set(flags["sources"]))
                ),
                "hpa_predicted_membrane_v5": clean(
                    flags["hpa"].get("hpa_predicted_membrane")
                ),
                "hpa_plasma_membrane_location_v5": clean(
                    flags["hpa"].get("hpa_plasma_membrane_location")
                ),
                "hpa_reliability_if_v5": clean(
                    flags["hpa"].get("hpa_reliability_if")
                ),
                "hpa_subcellular_location_v5": clean(
                    flags["hpa"].get("hpa_subcellular_location")
                ),
                "htp_present_v5": "1" if htp else "0",
                "htp_evidence_v5": clean(htp.get("htp_evidence")),
                "htp_num_tm_v5": clean(htp.get("htp_num_tm")),
                "htp_reliability_v5": clean(htp.get("htp_reliability")),
                "htp_n_terminal_side_v5": clean(htp.get("htp_n_terminal_side")),
                "htp_c_terminal_side_v5": clean(htp.get("htp_c_terminal_side")),
                "membranome_present_v5": "1" if flags["membranome"] else "0",
                "membranome_family_v5": ";".join(
                    sorted(
                        {
                            clean(x.get("membranome_family"))
                            for x in flags["membranome"]
                            if clean(x.get("membranome_family"))
                        }
                    )
                ),
                "membranome_segment_v5": ";".join(
                    sorted(
                        {
                            clean(x.get("membranome_segment"))
                            for x in flags["membranome"]
                            if clean(x.get("membranome_segment"))
                        }
                    )
                ),
                "opm_present_v5": "1" if flags["opm"] else "0",
                "opm_integral_structure_v5": "1" if flags["opm_integral"] else "0",
                "opm_pdb_ids_v5": ";".join(flags["opm_pdbs"]),
                "pdbtm_present_v5": "1" if flags["pdbtm_hits"] else "0",
                "pdbtm_pdb_ids_v5": ";".join(flags["pdbtm_hits"]),
                "single_pass_resolution_v5": topology,
                "single_pass_resolution_basis_v5": topology_basis,
                "functional_primary_class_v5": classification["primary"],
                "functional_subclass_v5": classification["subclass"],
                "classification_source_v5": classification["source"],
                "classification_status_v5": classification["status"],
                "gpcrdb_class_v5": clean(gpcr.get("gpcrdb_class")),
                "gpcrdb_family_v5": clean(gpcr.get("gpcrdb_family")),
                "gpcrdb_ligand_type_v5": clean(gpcr.get("gpcrdb_ligand_type")),
                "gpcrdb_subfamily_v5": clean(gpcr.get("gpcrdb_subfamily")),
                "gtopdb_type_v5": clean(gtop.get("gtopdb_type")),
                "gtopdb_family_id_v5": clean(gtop.get("gtopdb_family_id")),
                "gtopdb_family_name_v5": clean(gtop.get("gtopdb_family_name")),
                "tcdb_tcids_v5": classification["tcdb_tcids"],
                "cross_source_conflict_v5": "1" if conflict_reasons else "0",
                "review_reason_v5": ";".join(conflict_reasons),
            }
        )
        master_rows.append(record)

        matrix_row = {
            "target_uniprot_id": accession,
            "uniprot_entry_name": clean(record.get("uniprot_entry_name")),
            "approved_symbol": clean(record.get("approved_symbol")),
            "release_tier_v5": tier,
            "present_in_v4": "1" if existed_v4 else "0",
            "uniprot_membrane_scope_v4": original_scope,
            "hpa_predicted_membrane": record["hpa_predicted_membrane_v5"],
            "hpa_plasma_membrane_location": record["hpa_plasma_membrane_location_v5"],
            "htp_present": record["htp_present_v5"],
            "htp_evidence": record["htp_evidence_v5"],
            "htp_num_tm": record["htp_num_tm_v5"],
            "membranome_present": record["membranome_present_v5"],
            "opm_present": record["opm_present_v5"],
            "opm_integral_structure": record["opm_integral_structure_v5"],
            "pdbtm_present": record["pdbtm_present_v5"],
            "independent_source_count": record["independent_membrane_source_count_v5"],
            "decision": decision,
            "conflict": record["cross_source_conflict_v5"],
            "review_reason": record["review_reason_v5"],
        }
        evidence_matrix.append(matrix_row)
        if conflict_reasons:
            conflicts.append(matrix_row)
        classification_rows.append(
            {
                "target_uniprot_id": accession,
                "approved_symbol": clean(record.get("approved_symbol")),
                "v4_primary_class": clean(record.get("functional_primary_class")),
                "v5_primary_class": classification["primary"],
                "v5_subclass": classification["subclass"],
                "classification_source": classification["source"],
                "classification_status": classification["status"],
                "gpcrdb_class": record["gpcrdb_class_v5"],
                "gpcrdb_family": record["gpcrdb_family_v5"],
                "gtopdb_type": record["gtopdb_type_v5"],
                "gtopdb_family": record["gtopdb_family_name_v5"],
                "tcdb_tcids": record["tcdb_tcids_v5"],
            }
        )

    # Canonical resolution for every TrEMBL transmembrane candidate.
    canonical_indices = build_indices(reviewed_raw)
    unreviewed_resolution: list[dict] = []
    reviewed_entry_to_accession = entry_to_accession
    for candidate in unreviewed_raw:
        matches, methods = canonical_matches(candidate, canonical_indices)
        candidate_accession = clean(candidate.get("Entry"))
        candidate_entry = clean(candidate.get("Entry Name")).upper()
        htp = htp_index.get(candidate_entry, {})
        membranome_hit = candidate_accession in membranome_index
        external = []
        if htp:
            external.append("UniTmp_HTP")
        if membranome_hit:
            external.append("Membranome")
        if len(matches) == 1:
            status = "absorbed_to_reviewed_canonical"
            canonical = matches[0]
        elif len(matches) > 1:
            status = "ambiguous_reviewed_canonical_match"
            canonical = ";".join(matches)
        elif external:
            status = "retain_review_unreviewed_external_support"
            canonical = ""
        else:
            status = "retain_review_unreviewed_unconfirmed"
            canonical = ""
        unreviewed_resolution.append(
            {
                "candidate_uniprot_id": candidate_accession,
                "candidate_entry_name": candidate_entry,
                "candidate_gene_symbol": clean(candidate.get("Gene Names (primary)")),
                "candidate_hgnc": clean(candidate.get("HGNC")),
                "candidate_geneid": clean(candidate.get("GeneID")),
                "candidate_ensembl": clean(candidate.get("Ensembl")),
                "v4_gene_symbol_resolution_status": v4_unreviewed_status.get(
                    candidate_accession, ""
                ),
                "canonical_reviewed_uniprot_id": canonical,
                "canonical_match_methods": ";".join(methods),
                "canonical_resolution_status": status,
                "external_membrane_support": ";".join(external),
                "htp_evidence": clean(htp.get("htp_evidence")),
                "htp_num_tm": clean(htp.get("htp_num_tm")),
            }
        )

    # Resolve all 377 legacy rows using current accession first and approved
    # gene symbol second; distinguish identifier resolution from membrane scope.
    legacy_resolution: list[dict] = []
    symbol_index = canonical_indices["symbol"]
    final_by_accession = {r["target_uniprot_id"]: r for r in master_rows}
    for legacy in legacy_rows:
        old_id = clean(legacy.get("target_uniprot_id"))
        symbol = clean(legacy.get("approved_symbol")).upper()
        if old_id in reviewed_by_accession:
            canonical = old_id
            id_status = "current_reviewed_accession"
        else:
            candidates = sorted(symbol_index.get(symbol, set()))
            canonical = candidates[0] if len(candidates) == 1 else ";".join(candidates)
            id_status = (
                "mapped_by_approved_symbol"
                if len(candidates) == 1
                else "ambiguous_or_unresolved_identifier"
            )
        final = final_by_accession.get(canonical, {})
        if final:
            scope_status = (
                "included_core_or_extended"
                if final.get("release_tier_v5") in {"A", "B", "C"}
                else "external_support_review"
            )
        else:
            scope_status = "not_supported_as_membrane_in_v5"
        legacy_resolution.append(
            {
                **legacy,
                "canonical_uniprot_id_v5": canonical,
                "identifier_resolution_status_v5": id_status,
                "membrane_scope_resolution_v5": scope_status,
                "release_tier_v5": clean(final.get("release_tier_v5")),
                "independent_membrane_sources_v5": clean(
                    final.get("independent_membrane_sources_v5")
                ),
                "review_reason_v5": clean(final.get("review_reason_v5")),
            }
        )

    write_tsv(
        RELEASE / "human_membrane_protein_master_v5.tsv",
        master_rows,
        V4_FIELDS + V5_FIELDS,
    )
    write_tsv(
        RELEASE / "human_integral_membrane_view_v5.tsv",
        [r for r in master_rows if r["release_tier_v5"] == "A"],
        V4_FIELDS + V5_FIELDS,
    )
    write_tsv(
        RELEASE / "human_core_membrane_view_v5.tsv",
        [r for r in master_rows if r["release_tier_v5"] in {"A", "B"}],
        V4_FIELDS + V5_FIELDS,
    )
    write_tsv(
        RELEASE / "human_extended_membrane_view_v5.tsv",
        [r for r in master_rows if r["release_tier_v5"] in {"A", "B", "C"}],
        V4_FIELDS + V5_FIELDS,
    )
    write_tsv(
        REVIEW / "human_membrane_review_queue_v5.tsv",
        [r for r in master_rows if r["release_tier_v5"] == "R"],
        V4_FIELDS + V5_FIELDS,
    )
    write_tsv(CROSSWALKS / "source_evidence_matrix_v5.tsv", evidence_matrix)
    write_tsv(
        CROSSWALKS / "unreviewed_canonical_resolution_v5.tsv",
        unreviewed_resolution,
    )
    write_tsv(
        REVIEW / "unreviewed_82_no_reviewed_gene_symbol_match_v5.tsv",
        [
            row
            for row in unreviewed_resolution
            if row["v4_gene_symbol_resolution_status"]
            == "no_reviewed_gene_symbol_match"
        ],
    )
    write_tsv(
        CROSSWALKS / "legacy_377_resolution_v5.tsv", legacy_resolution
    )
    write_tsv(CROSSWALKS / "single_pass_resolution_v5.tsv", topology_rows)
    write_tsv(CROSSWALKS / "classification_resolution_v5.tsv", classification_rows)
    write_tsv(REVIEW / "cross_source_conflict_review_v5.tsv", conflicts)

    stats = {
        "release_date": TODAY,
        "reviewed_reference_proteome_rows": len(reviewed_raw),
        "v4_master_rows": len(v4_rows),
        "v5_union_rows": len(master_rows),
        "v5_new_reviewed_external_union_rows": sum(
            r["record_status_v5"] == "new_reviewed_external_union" for r in master_rows
        ),
        "release_tiers": dict(Counter(r["release_tier_v5"] for r in master_rows)),
        "external_source_coverage": {
            "HPA_predicted_membrane": sum(
                r["hpa_predicted_membrane_v5"] == "1" for r in master_rows
            ),
            "HPA_plasma_membrane_location": sum(
                r["hpa_plasma_membrane_location_v5"] == "1" for r in master_rows
            ),
            "UniTmp_HTP": sum(r["htp_present_v5"] == "1" for r in master_rows),
            "Membranome": sum(
                r["membranome_present_v5"] == "1" for r in master_rows
            ),
            "OPM": sum(r["opm_present_v5"] == "1" for r in master_rows),
            "PDBTM": sum(r["pdbtm_present_v5"] == "1" for r in master_rows),
        },
        "unreviewed_candidates": len(unreviewed_resolution),
        "unreviewed_resolution_status": dict(
            Counter(r["canonical_resolution_status"] for r in unreviewed_resolution)
        ),
        "unreviewed_without_single_canonical": sum(
            r["canonical_resolution_status"] != "absorbed_to_reviewed_canonical"
            for r in unreviewed_resolution
        ),
        "legacy_rows": len(legacy_resolution),
        "legacy_identifier_status": dict(
            Counter(r["identifier_resolution_status_v5"] for r in legacy_resolution)
        ),
        "legacy_membrane_scope_status": dict(
            Counter(r["membrane_scope_resolution_v5"] for r in legacy_resolution)
        ),
        "single_pass_original_unresolved": len(topology_rows),
        "single_pass_resolution": dict(
            Counter(r["resolved_topology_v5"] for r in topology_rows)
        ),
        "v4_unclassified": sum(
            clean(r.get("functional_primary_class")) == "unclassified"
            or clean(r.get("functional_subclass")) == "unclassified"
            for r in v4_rows
        ),
        "v5_unclassified_all_union": sum(
            r["functional_primary_class_v5"] == "unclassified" for r in master_rows
        ),
        "v5_unclassified_retained_v4": sum(
            r["record_status_v5"] == "retained_from_v4"
            and r["functional_primary_class_v5"] == "unclassified"
            for r in master_rows
        ),
        "classification_primary_counts": dict(
            Counter(r["functional_primary_class_v5"] for r in master_rows)
        ),
        "conflict_review_rows": len(conflicts),
        "normalized_source_rows": {
            "HPA_membrane_related": len(hpa_rows),
            "UniTmp_HTP": len(htp_rows),
            "Membranome_human": len(membranome_rows),
            "OPM_human_structures": len(opm_rows),
            "PDBTM_structures": len(pdbtm_rows),
            "TCDB_human": len(tcdb_rows),
            "GPCRdb_human": len(gpcrdb_rows),
            "GtoPdb_human_mappings": len(gtopdb_rows),
        },
    }
    (REPORTS / "V5_BUILD_STATS.json").write_text(
        json.dumps(stats, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(stats, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
